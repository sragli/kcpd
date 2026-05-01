defmodule Kcpd do
  @moduledoc """
  Kernel Change Point Detection (KCPD).

  Identifies structural breaks in a time series using kernel-based cost
  functions minimised via exact dynamic programming (DYNP).

  The cost of a segment `S = [a, b)` is derived from the Maximum Mean
  Discrepancy (MMD):

      C(a, b) = Σᵢ∈S k(xᵢ, xᵢ) − (1/|S|) Σᵢ∈S Σⱼ∈S k(xᵢ, xⱼ)

  ## Complexity

  Precomputation of the kernel matrix and 2-D prefix sums: `O(n²)`.
  DYNP search: `O(n² · K)` time, `O(n²)` space.

  ## Supported kernels

    * `:rbf`       — Radial Basis Function: `exp(-‖xᵢ − xⱼ‖² / 2σ²)`
    * `:linear`    — Linear:               `xᵢᵀ xⱼ`
    * `:laplacian` — Laplacian:            `exp(-‖xᵢ − xⱼ‖₁ / σ)`
    * Any 2-arity anonymous function `fn xi, xj -> ... end`

  ## Quick start

      iex> Kcpd.detect([0, 0, 0, 5, 5, 5], 1)
      [3, 6]

  """

  @doc """
  Detect `n_bkps` change points in `signal`.

  ## Parameters

    - `signal`  – list of observations. Each element may be a number
                  (univariate) or a list of numbers (multivariate).
    - `n_bkps`  – number of change points to detect (non-negative integer).
    - `opts`    – keyword options:

      | Key          | Default  | Description                                                          |
      |--------------|----------|----------------------------------------------------------------------|
      | `:kernel`    | `:rbf`   | `:rbf`, `:linear`, `:laplacian`, or a 2-arity function              |
      | `:bandwidth` | `1.0`    | Bandwidth σ. Use `:auto` for the median pairwise-distance heuristic. |

  ## Returns

  A sorted list of `n_bkps + 1` integers — the exclusive end positions
  (1-based) of each segment. The last element is always `length(signal)`.

  ## Examples

      iex> Kcpd.detect([0, 0, 0, 5, 5, 5], 1)
      [3, 6]

      iex> Kcpd.detect([0.0, 0.0, 1.0, 1.0, 0.0, 0.0], 2)
      [2, 4, 6]

  """
  def detect(signal, n_bkps, opts \\ [])

  def detect(signal, 0, _opts) when is_list(signal), do: [length(signal)]

  def detect(signal, n_bkps, opts)
      when is_list(signal) and is_integer(n_bkps) and n_bkps >= 1 do
    n = length(signal)

    if n_bkps >= n do
      raise ArgumentError,
            "n_bkps (#{n_bkps}) must be less than the signal length (#{n})"
    end

    kernel_opt = Keyword.get(opts, :kernel, :rbf)
    bandwidth_opt = Keyword.get(opts, :bandwidth, 1.0)

    signal_arr = List.to_tuple(signal)
    bandwidth = resolve_bandwidth(bandwidth_opt, signal_arr, n, kernel_opt)
    kernel_fn = build_kernel(kernel_opt, bandwidth)

    {prefix, diag} = precompute(signal_arr, n, kernel_fn)
    cost_fn = fn a, b -> segment_cost(a, b, prefix, diag) end

    dynp(n, n_bkps, cost_fn)
  end

  # ---------------------------------------------------------------------------
  # Kernel construction
  # ---------------------------------------------------------------------------

  defp build_kernel(:rbf, bw) do
    fn xi, xj -> :math.exp(-sq_dist(xi, xj) / (2.0 * bw * bw)) end
  end

  defp build_kernel(:linear, _bw) do
    fn xi, xj -> dot(xi, xj) end
  end

  defp build_kernel(:laplacian, bw) do
    fn xi, xj -> :math.exp(-l1_dist(xi, xj) / bw) end
  end

  defp build_kernel(f, _bw) when is_function(f, 2), do: f

  defp sq_dist(a, b) when is_number(a), do: (a - b) * (a - b)

  defp sq_dist(a, b) when is_list(a) do
    Enum.zip_reduce(a, b, 0.0, fn x, y, acc -> acc + (x - y) * (x - y) end)
  end

  defp l1_dist(a, b) when is_number(a), do: abs(a - b)

  defp l1_dist(a, b) when is_list(a) do
    Enum.zip_reduce(a, b, 0.0, fn x, y, acc -> acc + abs(x - y) end)
  end

  defp dot(a, b) when is_number(a), do: a * b

  defp dot(a, b) when is_list(a) do
    Enum.zip_reduce(a, b, 0.0, fn x, y, acc -> acc + x * y end)
  end

  # ---------------------------------------------------------------------------
  # Bandwidth: median pairwise-distance heuristic
  # ---------------------------------------------------------------------------

  defp resolve_bandwidth(:auto, arr, n, kernel) when kernel in [:rbf, :laplacian] do
    dists =
      for i <- 0..(n - 2), j <- (i + 1)..(n - 1), do: sq_dist(elem(arr, i), elem(arr, j))

    case dists do
      [] ->
        1.0

      _ ->
        sorted = Enum.sort(dists)
        median_sq = Enum.at(sorted, div(length(sorted), 2))
        max(:math.sqrt(max(median_sq, 0.0)), 1.0e-8)
    end
  end

  defp resolve_bandwidth(bw, _arr, _n, _kernel), do: bw

  # ---------------------------------------------------------------------------
  # Kernel matrix precomputation and 2-D prefix sums
  #
  # prefix[{i, j}]  = Σ_{r=0}^{i-1} Σ_{c=0}^{j-1} K[r][c]   (i, j ∈ 0..n)
  # diag[i]         = Σ_{r=0}^{i-1} K[r][r]                   (i ∈ 0..n)
  #
  # These allow O(1) segment cost evaluation via the rectangle-sum identity.
  # ---------------------------------------------------------------------------

  defp precompute(arr, n, kernel_fn) do
    # Upper-triangle (including diagonal) of the symmetric kernel matrix.
    k_upper =
      Map.new(
        for i <- 0..(n - 1), j <- i..(n - 1), do: {{i, j}, kernel_fn.(elem(arr, i), elem(arr, j))}
      )

    k_at = fn i, j ->
      if i <= j, do: k_upper[{i, j}], else: k_upper[{j, i}]
    end

    # Build the (n+1) x (n+1) 2-D prefix sum.
    # Boundary rows/columns are 0; recurrence:
    #   P[i][j] = K[i-1][j-1] + P[i-1][j] + P[i][j-1] - P[i-1][j-1]
    prefix =
      Enum.reduce(0..n, %{}, fn i, acc ->
        Enum.reduce(0..n, acc, fn j, m ->
          val =
            if i == 0 or j == 0 do
              0.0
            else
              k_at.(i - 1, j - 1) +
                Map.get(m, {i - 1, j}, 0.0) +
                Map.get(m, {i, j - 1}, 0.0) -
                Map.get(m, {i - 1, j - 1}, 0.0)
            end

          Map.put(m, {i, j}, val)
        end)
      end)

    # Diagonal prefix sum: diag[i] = Σ_{r<i} K[r][r]
    diag =
      Enum.reduce(1..n, %{0 => 0.0}, fn i, acc ->
        Map.put(acc, i, acc[i - 1] + k_at.(i - 1, i - 1))
      end)

    {prefix, diag}
  end

  # ---------------------------------------------------------------------------
  # Segment cost C(a, b) for the segment [a, b)   (0-indexed, |S| = b − a)
  #
  # The rectangle-sum identity gives the intra-segment kernel sum in O(1):
  #   Σᵢ∈S Σⱼ∈S K[i][j] = P[b][b] - P[a][b] - P[b][a] + P[a][a]
  # ---------------------------------------------------------------------------

  defp segment_cost(a, b, prefix, diag) when b > a do
    sz = b - a
    diag_sum = diag[b] - diag[a]
    k_sum = prefix[{b, b}] - prefix[{a, b}] - prefix[{b, a}] + prefix[{a, a}]
    diag_sum - k_sum / sz
  end

  # ---------------------------------------------------------------------------
  # Dynamic programming (DYNP) — exact globally optimal segmentation
  #
  # dp[{k, t}] = {min_cost, prev_t}
  #   best total cost to cover [0, t) with exactly k segments, and the
  #   start of the last segment (= end of the (k-1)-th segment).
  # ---------------------------------------------------------------------------

  defp dynp(n, n_bkps, cost_fn) do
    # k = 1: each t is a single segment [0, t).
    dp =
      Map.new(1..n, fn t ->
        {{1, t}, {cost_fn.(0, t), 0}}
      end)

    # Fill k = 2 .. n_bkps + 1.
    dp =
      Enum.reduce(2..(n_bkps + 1), dp, fn k, dp_acc ->
        Enum.reduce(k..n, dp_acc, fn t, m ->
          {best_cost, best_s} =
            Enum.reduce((k - 1)..(t - 1), {:infinity, 0}, fn s, {mc, ms} ->
              case m[{k - 1, s}] do
                nil ->
                  {mc, ms}

                {prev_cost, _} ->
                  c = prev_cost + cost_fn.(s, t)
                  if c < mc, do: {c, s}, else: {mc, ms}
              end
            end)

          Map.put(m, {k, t}, {best_cost, best_s})
        end)
      end)

    backtrack(dp, n_bkps + 1, n, [])
  end

  # Recover the breakpoint positions by walking the prev_t pointers.
  defp backtrack(_dp, 1, t, acc), do: Enum.sort([t | acc])

  defp backtrack(dp, k, t, acc) do
    {_, s} = dp[{k, t}]
    backtrack(dp, k - 1, s, [t | acc])
  end
end

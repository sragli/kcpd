# Kernel Change Point Detection

An Elixir implementation of the **Kernel Change Point Detection (KCPD)** algorithm. It locates structural breaks in univariate or multivariate time series by minimising a kernel-based cost function via exact dynamic programming (DYNP).

## How it works

KCPD frames change point detection as a **penalised optimal partitioning** problem. Given a signal of length `n` and a desired number of change points `K`, the algorithm finds the partition into `K + 1` segments that minimises the total cost.

### 1. Kernel matrix

A symmetric `n × n` kernel matrix `K` is computed from the signal, where `K[i][j] = k(xᵢ, xⱼ)`. The kernel `k` is a positive-definite function that measures similarity between observations. Only the upper triangle is stored; the lower triangle is accessed by symmetry.

### 2. Segment cost (MMD-based)

The cost of a segment `S = [a, b)` measures how spread out the observations are *within* that segment. It is derived from the Maximum Mean Discrepancy (MMD):

```
C(a, b) = Σᵢ∈S k(xᵢ, xᵢ) − (1/|S|) Σᵢ∈S Σⱼ∈S k(xᵢ, xⱼ)
```

Intuitively, the first term is the sum of self-similarities and the second is the mean pairwise similarity. A homogeneous segment (all observations drawn from the same distribution) has a low cost; a segment spanning a structural change has a high cost.

### 3. O(1) segment cost via prefix sums

Naively evaluating `C(a, b)` requires O(|S|²) work. Instead, two auxiliary structures are built once from `K`:

- **2-D prefix sum** `P[i][j] = Σ_{r<i} Σ_{c<j} K[r][c]` — allows any rectangular sum over `K` to be retrieved in O(1) using the identity `Σ_{r=a}^{b-1} Σ_{c=a}^{b-1} K[r][c] = P[b][b] − P[a][b] − P[b][a] + P[a][a]`.
- **Diagonal prefix sum** `D[i] = Σ_{r<i} K[r][r]` — gives the self-similarity sum `Σᵢ∈S k(xᵢ, xᵢ) = D[b] − D[a]` in O(1).

Both structures are computed in O(n²) time and space.

### 4. Dynamic programming (DYNP)

With O(1) cost queries available, the globally optimal segmentation is found by the classic DYNP recurrence:

```
dp[k][t] = min_{s < t} ( dp[k-1][s] + C(s, t) )
```

where `dp[k][t]` is the minimum total cost to cover `[0, t)` with exactly `k` segments. The base case `dp[1][t] = C(0, t)` is filled for all `t`, and then `k` is incremented up to `K + 1`. The optimal breakpoints are recovered by backtracking the stored argmin pointers.

The full search runs in **O(n² · K)** time.

## Installation

Add `:kcpd` to your `mix.exs` dependencies:

```elixir
def deps do
  [
    {:kcpd, "~> 0.1.0"}
  ]
end
```

## Usage

```elixir
# Univariate — detect one change point in a step signal
Kcpd.detect([0, 0, 0, 5, 5, 5], 1)
#=> [3, 6]

# Two change points
Kcpd.detect([0.0, 0.0, 1.0, 1.0, 0.0, 0.0], 2)
#=> [2, 4, 6]

# Multivariate (list of lists)
Kcpd.detect([[0, 0], [0, 0], [0, 0], [5, 5], [5, 5], [5, 5]], 1)
#=> [3, 6]

# Laplacian kernel with auto bandwidth
Kcpd.detect(signal, 2, kernel: :laplacian, bandwidth: :auto)

# Custom kernel function
Kcpd.detect(signal, 1, kernel: fn xi, xj -> :math.exp(-abs(xi - xj)) end)
```

The return value is a sorted list of **exclusive end positions** (1-based). A signal of length `n` with `K` change points produces a list of `K + 1` integers; the last element is always `n`.

## Options

| Option       | Default  | Description                                                             |
|--------------|----------|-------------------------------------------------------------------------|
| `:kernel`    | `:rbf`   | `:rbf`, `:linear`, `:laplacian`, or a 2-arity `fn xi, xj -> ...` function |
| `:bandwidth` | `1.0`    | Bandwidth σ for RBF / Laplacian. Pass `:auto` to use the median pairwise-distance heuristic. |

## Kernels

| Name          | Formula                              |
|---------------|--------------------------------------|
| `:rbf`        | `exp(-‖xᵢ − xⱼ‖² / 2σ²)`           |
| `:linear`     | `xᵢᵀ xⱼ`                            |
| `:laplacian`  | `exp(-‖xᵢ − xⱼ‖₁ / σ)`             |

## Complexity

| Phase              | Time   | Space  |
|--------------------|--------|--------|
| Kernel matrix      | O(n²)  | O(n²)  |
| 2-D prefix sums    | O(n²)  | O(n²)  |
| DYNP (K segments)  | O(n²K) | O(n²)  |

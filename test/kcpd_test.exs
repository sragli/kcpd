defmodule KCPDTest do
  use ExUnit.Case
  doctest KCPD

  describe "detect/3 — zero breakpoints" do
    test "returns [n] immediately" do
      assert KCPD.detect([1, 2, 3, 4], 0) == [4]
    end
  end

  describe "detect/3 — RBF kernel (default)" do
    test "single change point in step signal" do
      assert KCPD.detect([0, 0, 0, 5, 5, 5], 1) == [3, 6]
    end

    test "two change points" do
      assert KCPD.detect([0.0, 0.0, 1.0, 1.0, 0.0, 0.0], 2) == [2, 4, 6]
    end

    test "last element is always the signal length" do
      signal = List.duplicate(0.0, 10) ++ List.duplicate(5.0, 10)
      bkps = KCPD.detect(signal, 1)
      assert List.last(bkps) == 20
    end

    test "single change point in long signal" do
      signal = List.duplicate(0.0, 50) ++ List.duplicate(10.0, 50)
      assert KCPD.detect(signal, 1) == [50, 100]
    end
  end

  describe "detect/3 — alternative kernels" do
    test "linear kernel detects step change" do
      assert KCPD.detect([0.0, 0.0, 0.0, 5.0, 5.0, 5.0], 1, kernel: :linear) == [3, 6]
    end

    test "laplacian kernel detects step change" do
      assert KCPD.detect([0.0, 0.0, 0.0, 5.0, 5.0, 5.0], 1, kernel: :laplacian) == [3, 6]
    end

    test "custom 2-arity kernel function" do
      kernel = fn xi, xj -> :math.exp(-abs(xi - xj)) end
      assert KCPD.detect([0.0, 0.0, 0.0, 5.0, 5.0, 5.0], 1, kernel: kernel) == [3, 6]
    end
  end

  describe "detect/3 — bandwidth option" do
    test "auto bandwidth detects step change" do
      signal = [0.0, 0.0, 0.0, 5.0, 5.0, 5.0]
      assert KCPD.detect(signal, 1, bandwidth: :auto) == [3, 6]
    end

    test "explicit bandwidth" do
      signal = [0.0, 0.0, 0.0, 5.0, 5.0, 5.0]
      assert KCPD.detect(signal, 1, bandwidth: 2.0) == [3, 6]
    end
  end

  describe "detect/3 — multivariate signal" do
    test "detects change in 2-D signal" do
      signal = [[0, 0], [0, 0], [0, 0], [5, 5], [5, 5], [5, 5]]
      assert KCPD.detect(signal, 1) == [3, 6]
    end
  end

  describe "detect/3 — error handling" do
    test "raises when n_bkps >= signal length" do
      assert_raise ArgumentError, fn -> KCPD.detect([1, 2, 3], 3) end
    end
  end
end

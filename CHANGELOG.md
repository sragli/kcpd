# Changelog

## [0.1.0] - 2026-05-01

### Added

- `KCPD.detect/3` — exact change point detection via kernel-based dynamic programming (DYNP).
- Built-in kernels: `:rbf` (Radial Basis Function, default), `:linear`, and `:laplacian`.
- Support for custom 2-arity kernel functions passed as anonymous functions.
- `:bandwidth` option accepting a numeric value or `:auto` (median pairwise-distance heuristic).
- Univariate (number) and multivariate (list of numbers) signal support.
- O(1) segment cost evaluation via 2-D prefix sums of the kernel matrix.

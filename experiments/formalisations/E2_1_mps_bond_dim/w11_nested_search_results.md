# w11_nested_search — results

**Result (S206):** Nested `1 + (1 + 3 + 3 + 3)` BlockTriangular search
over rows in [1, row_pool):

* rows in [1, 11): **0 candidates** (structurally obstructed).
* rows in [1, 22): **200+ candidates**, best with max_row = 17 →
  forces j ≥ 2 since 11^j ≥ 17 ⇔ j ≥ 2.

Top candidate (best max_row = 17):
```
ρ = (0, 1, 2, 8, 16, 3, 4, 10, 5, 11, 17)
σ = (10, 1, 0, 4, 6, 2, 7, 8, 3, 5, 9)
```
detA = detB = detC = 1. Six new prime helpers
`{113, 131, 181, 191, 193, 197}` (count 6).

**Falsification statement:** the conjectured nested
`1 + (1 + 3 + 3 + 3)` decomposition over rows [0, 11) (single-session
viable for j ≥ 1) does not exist. The j ≥ 2 weakening is the closest
single-session viable variant.

See `w11_blocktriangular_search_results.md` for context.

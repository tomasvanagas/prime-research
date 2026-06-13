# W=14 j=1 analysis — results

Results of `w14_j1_analysis.py`. This script analyses the W=14 j=1
inner 6×6 structure for the orthogonal corner case
`mps_bond_dim_W_eq_14_d_eq_j_plus_1`.

**Findings.**
* Rank of 14×14 j=1 slab = 7 (matches E2.1's `phi(14) + 1 = 7`).
* 587 invertible inner 6×6 submatrices (rows in `[1, 14)` × live cols
  `{0, 2, 4, 8, 10, 12}`).
* 209 of those admit a `(5+1)` BT structure where one corner row × col
  is a single nonzero entry. Only row 8 has support size 1 (chiP(113)
  at col 0 in the live-col selection); thus the "1" block is forced
  to be row 8 × col 0 in every valid (5+1) decomposition.
* **0 of those 209 admit a further `(3+2)` or `(2+3)` BT of the inner
  5×5** — the residual 5×5 is shape-atomic for parts ≤ 3 in this
  partition.

**Implication.** The (5+1) → (3, 2) / (2, 3) cascade fails. Need to try
shapes that split row 8's pinning differently — specifically `(2, 3, 1)`
or `(2, 4)` or others where row 8 is paired with a different row group.

**See `w14_j1_shape_search_results.md`** for the full shape sweep that
identified `(2, 3, 1)` as the unique parts-≤-3 shape that works.

The comprehensive findings + Lean closure summary is in
`w14_blocktriangular_search_results.md`.

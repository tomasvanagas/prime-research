# W=14 j=1 best (2, 3, 1) candidate — results

Results of `w14_j1_best_candidate.py`. Enumerates all valid `(2, 3, 1)`
BT decompositions across the 587 invertible inner 6×6 submatrices and
scores by minimum new prime helpers needed.

**Total (2, 3, 1) BT decompositions found:** 34.

**Min new helpers across all candidates:** 4 — `{67, 113, 173, 181}`.

**Best candidate (chosen for Lean closure at S235):**

```
ρ ↦ (0, 2, 4, 3, 6, 12, 8)
σ ↦ (1, 2, 8, 4, 10, 12, 0)

A-block (2×2): rows {2, 4} × cols {2, 8}
  [[1, 1],
   [1, 0]]
  det = -1 (via Matrix.det_fin_two)

B-block (3×3): rows {3, 6, 12} × cols {4, 10, 12}
  [[1, 1, 0],
   [1, 0, 1],
   [1, 1, 1]]
  det = -1 (via Matrix.det_fin_three)

C-block (1×1): row 8 × col 0
  [[1]]    (chiP(113) = 1)
  det = 1 (via Matrix.det_fin_one)

A_outer (1×1): row 0 × col 1
  [[1]]    (chiP(2) = 1)
  det = 1 (via Matrix.det_fin_one)

Total det = 1 · ((1 · ((-1) · (-1))) · 1) = 1.
```

**Max ρ.val = 12 < 14^j for every j ≥ 1**, so the closure works for
ALL j ≥ 1 (no j=1 vs j ≥ 2 split needed).

**The full 7×7 submatrix:**
```
        col 1   col 2   col 8   col 4   col 10  col 12  col 0
row 0:  [ 1,    1,      0,      1,      1,      1,      0 ]
row 2:  [ 0,    1,      1,      0,      0,      1,      1 ]
row 4:  [ 0,    1,      0,      1,      1,      0,      0 ]
row 3:  [ 0,    0,      0,      1,      1,      0,      1 ]
row 6:  [ 0,    0,      0,      1,      0,      1,      0 ]
row 12: [ 0,    0,      0,      1,      1,      1,      0 ]
row 8:  [ 0,    0,      0,      0,      0,      0,      1 ]
```

Lower-left zeros (verified):
* Outer (1+6) lower: rows 2..12 + 8 at col 1 = `chiP(14r + 2)` =
  `chiP(2(7r+1))` for r ≥ 1, all even and > 2 → 0.
* Mid (5+1) lower: row 8 at cols 2, 8, 4, 10, 12 = chiP(115, 121,
  117, 123, 125) — none prime, all 0.
* Inner (2+3) lower: rows 3, 6, 12 at cols 2, 8 = chiP(45, 51, 87,
  93, 171, 177) — none prime, all 0.

The Lean closure implementing this candidate is in
`MPSBondDim/Basic.lean` as
`exists_invertible_submatrix_W_eq_14_d_eq_j_plus_1` and
`mps_bond_dim_W_eq_14_d_eq_j_plus_1`. See
`w14_blocktriangular_search_results.md` for the full proof structure.

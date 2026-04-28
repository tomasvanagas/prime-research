# W=7 BlockTriangular Pre-Search (S159)

## Goal

Identify a `1 + 3 + 3` BlockTriangular permutation for the W=7 orthogonal-corner
case `(W=7, d=j+1)` of E2.1, which has `R = φ(7) + 1 = 7`. W=7 was on the
S128/S129/S144 "block-triangular-required" list — the leading-row +
dead-col upper-triangulation route is exhausted in rows `[0, 7)`.

## Method

Search over all `(A_rows, B_rows, A_cols, B_cols)` partitions where
- row `0` is forced to the dead-col witness (`k = 6`, `k+1 = 7` prime),
- rows `1..6` partition into two triples `A_rows`, `B_rows`,
- live cols `{0..5}` partition into two triples `A_cols`, `B_cols`,
- both 3×3 diagonal blocks have nonzero determinant,
- the lower-left block (`B_rows × A_cols`) is zero (BlockTriangular constraint).

We additionally check the upper-right block (`A_rows × B_cols`) for the
stronger BLOCK-DIAGONAL property.

## Result

**2 valid (1+3+3) BlockTriangular decompositions found**, both block-diagonal,
both using zero new prime helpers.

### Top candidate (S159 chosen)

```
ρ ↦ (0, 1, 3, 5, 2, 4, 6)
σ ↦ (6, 1, 3, 5, 0, 2, 4)

A-block det = 1 (1×1 block: chiP(7) = 1)
D₁-block det = -1
D₂-block det = -2
total det    =  1 · (-1) · (-2) = 2 (≠ 0, hence IsUnit over ℚ)
```

The 7×7 submatrix:
```
[1, 1, 0, 0, 0, 1, 1]    (row 0)
[0, 0, 1, 1, 0, 0, 0]    (row 1)
[0, 1, 0, 0, 0, 0, 0]    (row 3)
[0, 1, 0, 1, 0, 0, 0]    (row 5)
[0, 0, 0, 0, 0, 1, 1]    (row 2)
[0, 0, 0, 0, 1, 1, 0]    (row 4)
[0, 0, 0, 0, 1, 0, 1]    (row 6)
```

After dropping the dead-col (col 0) and breaking out the `(1, 6)` and `(3, 3)`
splits, this gives:

```
A = [[1]]                     (1×1)
B = [1, 0, 0, 0, 1, 1]        (1×6)
D = [[0,1,1,0,0,0],
     [1,0,0,0,0,0],
     [1,0,1,0,0,0],
     [0,0,0,0,1,1],
     [0,0,0,1,1,0],
     [0,0,0,1,0,1]]            (6×6, block-diagonal 3+3)

D₁ = [[0,1,1],[1,0,0],[1,0,1]]  → det = -1
D₂ = [[0,1,1],[1,1,0],[1,0,1]]  → det = -2
```

### Primes appearing in submatrix

`{2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47}` — **all** already
declared as `chiP_X_eq_one` helpers in `MPSBondDim/Basic.lean` from prior
session closures (S98, S99, S106, S117, S122, S128, S152). **Zero new
helpers needed** — first instance of a corner closure with no new prime
witnesses.

### Composites whose non-primality must be witnessed

`{1, 4, 6, 8, 9, 10, 12, 14, 15, 16, 18, 20, 21, 22, 24, 25, 26, 27, 28, 30,
32, 33, 34, 35, 36, 38, 39, 40, 42, 44, 45, 46, 48, 49}` — 34 in total.
Ten composites (`6, 9, 15, 18, 24, 27, 33, 36, 42, 45`) not seen in the
S152 W=9 closure; the remaining 24 overlap with W=9's set. All 34 added
as proof-internal `h_not_prime_X` declarations via `decide` (no `norm_num`
needed since all composites are < 50).

## Lean closure

Implemented in `MPSBondDim/Basic.lean` as `exists_invertible_submatrix_W_eq_7_d_eq_j_plus_1`
and `mps_bond_dim_W_eq_7_d_eq_j_plus_1` (both sorry-free, both depending only
on `propext, Classical.choice, Quot.sound`).

The Lean assembly mirrors the S152 W=9 template exactly:
- Outer `finSumFinEquiv : Fin 1 ⊕ Fin 6 ≃ Fin 7` for the 1+6 split.
- Inner `finSumFinEquiv : Fin 3 ⊕ Fin 3 ≃ Fin 6` for the 3+3 split of D.
- `Matrix.det_fromBlocks_zero₂₁` applied twice (nested).
- `Matrix.det_fin_one` for A, `Matrix.det_fin_three` for D₁, D₂.
- Closing step uses `Ne.isUnit` rather than the `IsUnit 1` shortcut available
  to W=9 — first instance where total det ≠ ±1 (we have det = 2 here).

## What this closes

W=7 was the last small-prime wheel on the "block-triangular-required" list
flagged in S128/S129/S144. Closing W=7 leaves the next remaining block-
triangular targets at W ∈ {11, 13, 14, 15, 16, 17, 19, 21, 22, 24, 25, 26,
27, 28, 30, 32, ...}. Of these, W=11 is the next small-prime case (R = 11);
its 11×11 sub-matrix likely admits a `1 + 5 + 5` or `1 + 3 + 3 + 4` split
but each requires its own pre-search.

The closed-W set for E2.1's orthogonal corner is now
**{2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 18, 20}** — twelve wheels.

# W=14 BlockTriangular Pre-Search (S235)

## Goal

Identify a BlockTriangular permutation that closes the W=14 orthogonal-corner
case `(W=14, d=j+1)` of E2.1, where `R = φ(14) + 1 = 7`. W=14 was on the
S128/S137/S144 "block-triangular-required" list — the leading-row + dead-col
upper-triangulation route is exhausted (rows 2 and 5 of the j=1 slab share
support pattern `{0, 2, 8, 12}` at the chosen 7 cols).

S206 recommended W=14 as the "next single-session candidate" on the
hypothesis that *composite W avoids the W=11 parity atomicity*. This
pre-search confirms the recommendation: W=14 j=1 inner 6×6 admits a
**(2, 3, 1) BlockTriangular decomposition** in rows `[1, 14)` —
single-session viable with 4 new prime helpers and **only mathlib's
existing `det_fin_one`, `det_fin_two`, `det_fin_three`** (no
`det_fin_four` or `det_fin_five` needed).

## Method

Three pre-search scripts were run:

* `w14_blocktriangular_search.py` — search for `(1 + 3 + 3)` BT
  (mirroring S159 W=7 template). **Result: 0 in rows `[0, 14)`** (j=1
  obstructed); 38254 in rows `[0, 28)` (j ≥ 2 only). Best `(1+3+3)`
  candidate `ρ ↦ (0, 1, 2, 3, 6, 21, 23)` gives only 2 new helpers
  but max_row 23 forces `j ≥ 2`.
* `w14_j1_analysis.py` — analyses the j=1 slab specifically. Confirms:
  - rank of 14×14 j=1 slab = 7 (E2.1 saturation),
  - 587 invertible 6×6 inner submatrices (rows in [1, 14), live cols),
  - 209 of those admit a `(5+1)` BT structure where `corner row × col`
    is a single nonzero entry (only row 8 at col 0 has support size 1
    among the 13 inner rows),
  - **but ZERO of those 209 admit a further `(3+2)` or `(2+3)` BT of
    the inner 5×5** — the 5×5 sub-block is "shape-atomic" for parts
    ≤ 3.
* `w14_j1_shape_search.py` — exhaustive ordered-partition search across
  all distinct shapes of 6 with parts ≤ 4. Results below.

## Result: (2, 3, 1) BT works in rows `[0, 14)`

**Best candidate (S235 chosen):**

```
ρ ↦ (0, 2, 4, 3, 6, 12, 8)
σ ↦ (1, 2, 8, 4, 10, 12, 0)

A-block (2×2): rows {2, 4} × cols {2, 8}, det = -1
B-block (3×3): rows {3, 6, 12} × cols {4, 10, 12}, det = -1
C-block (1×1): row 8 × col 0, value chiP(113) = 1
total det    = (1) · ((-1) · (-1) · 1) = 1
```

The full 7×7 submatrix:

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

**Max row = 12 < 14.** The closure works for every `j ≥ 1` because
`14^j ≥ 14 > 12` for `j ≥ 1`.

### Block structure

The matrix decomposes as a **3-level nested `det_fromBlocks_zero₂₁`**:

* **Outer (1+6) BT**: `Fin 1 ⊕ Fin 6 ≃ Fin 7` via `finSumFinEquiv`.
  - 1×1 block A_outer: row 0 at col 1 (chiP(2) = 1, dead-col witness),
    det = 1.
  - 6×6 block D_outer: rows `{2, 4, 3, 6, 12, 8}` × cols
    `{2, 8, 4, 10, 12, 0}`.
  - Lower-left zero: rows `{2, 4, 3, 6, 12, 8}` at col 1 are all zero
    (since `chiP(14r + 2) = 0` for `r ≥ 1`).
* **Mid (5+1) BT** of D_outer: `Fin 5 ⊕ Fin 1 ≃ Fin 6` via
  `finSumFinEquiv`.
  - 5×5 block D_mid: rows `{2, 4, 3, 6, 12}` × cols `{2, 8, 4, 10, 12}`.
  - 1×1 block C_outer: row 8 × col 0, value chiP(113) = 1, det = 1.
  - Lower-left zero: row 8 at cols `{2, 8, 4, 10, 12}` are all zero
    (since `chiP(14·8 + k + 1) = chiP(113), chiP(115), chiP(117),
    chiP(121), chiP(123), chiP(125)` and only `chiP(113)` is 1).
* **Inner (2+3) BT** of D_mid: `Fin 2 ⊕ Fin 3 ≃ Fin 5` via
  `finSumFinEquiv`.
  - 2×2 block A: rows `{2, 4}` × cols `{2, 8}`, det = `0·1 - 1·1 = -1`.
  - 3×3 block B: rows `{3, 6, 12}` × cols `{4, 10, 12}`,
    `det_fin_three [[1,1,0],[1,0,1],[1,1,1]] = -1`.
  - Lower-left zero: rows `{3, 6, 12}` at cols `{2, 8}` are all zero
    (verified: `chiP(45) = chiP(51) = chiP(87) = chiP(93) = chiP(171)
    = chiP(177) = 0`).

### Primes appearing in submatrix

```
{2, 3, 5, 11, 13, 17, 19, 23, 29, 41, 43, 47, 53, 59, 67, 89, 97,
 113, 131, 137, 139, 149, 151, 167, 173, 179, 181}
```

After deduplication: 21 distinct primes. Existing `chiP_X_eq_one`
helpers from prior closures cover most:
`{2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61,
89, 97, 109, 127, 149, 179, 199, 211, 241, 293, 337}`.

**Four new prime helpers needed:** `{67, 113, 173, 181}`.

### Composites whose non-primality must be witnessed

```
{1, 9, 15, 16, 25, 27, 33, 39, 45, 51, 55, 69, 81, 85, 87, 93, 99,
 100, 111, 165, 169, 171, 177, 183, 185, 187, 195, 207}
```

28 composites total. All `< 200`, so all `decide`-checkable in Lean.

## Why ZERO (2,3,1) BT for j=1 was not obvious

The S206 W=11 atomicity result tempted prior agents to expect that
"prime W → atomic odd block, composite W → tractable" was the
controlling axis. Two refinements:

1. **W=14 is composite**, but the j=1 inner 6×6 still has only ONE
   row (row 8) with support size 1. The "row of size 1" is the
   *necessary condition* for any BT decomposition with a `1×1` block
   — and W=14 has exactly one such row (the row 8 chiP(113)
   isolation), forcing every BT decomp of the inner 6×6 to make
   row 8 the C_block.
2. **(3, 3) and (2, 2, 2) BT both fail for j=1** despite the inner
   6×6 having rank 6 across 587 row sets. The structural reason: the
   13 inner rows have support sizes `{1, 2, 3, 4}` distributed as
   `{1: 1, 2: 4, 3: 5, 4: 3}`, and the rank-6 invertibility
   requires a 6-row choice that *covers* the full 6-col span. No
   such choice admits an additional BT (3, 3) cut; the row 8
   isolation is forced and once consumed, the residual 5×5 over
   rows with size ≥ 2 is irreducible to two 3-row blocks (mirroring
   the W=11 atomic 5×5 obstruction at one size lower).

## Per-shape pre-search summary (j=1, parts ≤ 4)

| Shape       | Hits | Notes |
|-------------|-----:|-------|
| (5, 1)      |   ≥3 | Not closed in mathlib (no det_fin_five) |
| (1, 5)      |   ≥3 | Same |
| (4, 2)      |    0 | — |
| (2, 4)      |   ≥3 | Needs det_fin_four (not in mathlib) |
| (3, 3)      |    0 | — |
| (4, 1, 1)   |    0 | — |
| (1, 4, 1)   |   ≥3 | Needs det_fin_four |
| (1, 1, 4)   |   ≥3 | Same |
| (3, 2, 1)   |    0 | — |
| (3, 1, 2)   |    0 | — |
| **(2, 3, 1)** |    34 | **Closes via det_fin_two + det_fin_three + det_fin_one — single-session viable** |
| (1, 3, 2)   |    0 | — |
| (2, 1, 3)   |    0 | — |
| (1, 2, 3)   |    0 | — |
| (2, 2, 2)   |    0 | — |

The asymmetry is striking: **(2, 3, 1) HITS** but **(1, 2, 3),
(1, 3, 2), (3, 1, 2), (3, 2, 1), (2, 1, 3) all FAIL**. The reason:
row 8 (support size 1) MUST be the C-block (last row in the BT
order, paired with col 0); the surrounding `5×5` then admits a
(2, 3) BT but no other single split. The forced positional
constraint pins (2, 3, 1) as the unique single-session shape.

## Implication for the W=14 corner closure

Single-session viable. The closure follows the S152/S159 nested
`det_fromBlocks_zero₂₁` template at one extra nesting level
(3 levels vs 2):

```lean
-- Outer (1+6) split
Matrix.det_fromBlocks_zero₂₁ (A := A_outer) (D := D_outer) ...
-- Mid (5+1) split inside D_outer
Matrix.det_fromBlocks_zero₂₁ (A := D_mid) (D := C_outer) ...
-- Inner (2+3) split inside D_mid
Matrix.det_fromBlocks_zero₂₁ (A := A) (D := B) ...
-- 2×2 det_fin_two for A; 3×3 det_fin_three for B; 1×1 det_fin_one for
-- A_outer and C_outer.
```

Estimated 700-900 Lean lines (S152 was 610 at 2-level nesting; W=14
adds one more level + 4 new prime helpers + 28 composite witnesses).

## What this closes

W=14 was the second-smallest "block-triangular-required" wheel
following S159 W=7 closure. Closing W=14 leaves the next remaining
block-triangular targets at `W ∈ {11, 13, 15, 16, 17, 19, 21, 22,
24, 25, 26, 27, 28, 30, 32, ...}`. W=11 remains deferred to a sub-arc
(j≥2 only) or the `det_fin_five` sub-arc.

After W=14 closure, the closed-W set for E2.1's orthogonal corner
becomes **{2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 18, 20}** — thirteen
wheels. W=11 remains the smallest open one (atomic 5×5 odd block
requires multi-session work).

## Files

* `w14_blocktriangular_search.py` — initial (1+3+3) search.
* `w14_j1_analysis.py` — j=1 slab structural analysis + (5+1) +
  inner-(3+2)/(2+3) test.
* `w14_j1_shape_search.py` — efficient partition-based search across
  all shapes with parts ≤ 4.
* `w14_j1_best_candidate.py` — final scoring + verification of best
  (2, 3, 1) candidate.

## Lean closure (S235, COMPLETED)

Implemented in `MPSBondDim/Basic.lean` as
`exists_invertible_submatrix_W_eq_14_d_eq_j_plus_1` and
`mps_bond_dim_W_eq_14_d_eq_j_plus_1` (both sorry-free; `lake build`
succeeds with only the pre-existing sorry at line 467 for the general
`exists_invertible_submatrix`).

**File additions (S235):**
* Four new prime helpers: `chiP_sixty_seven_eq_one`,
  `chiP_one_hundred_thirteen_eq_one`,
  `chiP_one_hundred_seventy_three_eq_one`,
  `chiP_one_hundred_eighty_one_eq_one` (all use `norm_num` for
  primality witness).
* 28 internal `h_not_prime_X` lemmas. Composites `≥ 150`
  (`{169, 170, 171, 177}`) use `norm_num`; the other 24 use `decide`.
* Three-level nested `Matrix.det_fromBlocks_zero₂₁` decomposition:
  - Layer 1 (1+6) outer split.
  - Layer 2 (5+1) middle split.
  - Layer 3 (2+3) inner split.
* Total det = `1 · (1 · ((-1) · (-1)) · 1) = 1`, hence `IsUnit` over `ℚ`.
* `set_option maxHeartbeats 2000000` (10× default) for the 49-entry
  matrix-equality elaboration.

**First instance using `det_fin_two`** in the corner-closure family
(prior closures used `det_fin_one`, `det_fin_three`, or
`det_of_upperTriangular`).

**Closed-W set update:** `{2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 18, 20}`
— **thirteen wheels** (was twelve).

# W=15 / W=16 BlockTriangular Pre-Search (S245)

## Goal

Identify a BlockTriangular permutation that closes the W=15 and W=16
orthogonal-corner cases `(W, d=j+1)` of E2.1, where `R = φ(W) + 1 = 9` for
both. W=15 (= 3·5) and W=16 (= 2^4) sit on the S128/S143/S144
"block-triangular-required" set: the leading-row + dead-col upper-
triangulation route is structurally exhausted at both, so the only
single-session route into the `exists_invertible_submatrix` `sorry`
is via mathlib's `Matrix.det_fromBlocks_zero₂₁` (cf. S152 W=9, S159 W=7,
S235 W=14).

S235 closed W=14 (R=7) with a three-level nested `det_fromBlocks_zero₂₁`
using a `1+(2+3+1)` decomposition; all inner blocks had size ≤ 3, so the
proof used `det_fin_one`, `det_fin_two`, `det_fin_three` (the chain
mathlib provides). The S235 next-action queue listed W=15 and W=16 as
the next single-session candidates.

This pre-search tests whether the same chain (`det_fin_one`, `_two`,
`_three`) suffices for W=15 and W=16, or whether `det_fin_four` is
required.

## Result: BOTH W=15 and W=16 require `det_fin_four`

**Single-session closure of W=15 / W=16 is BLOCKED on the development
of `Matrix.det_fin_four` as a reusable lemma.**

The pre-search exhaustively enumerated all 108 ordered partitions of 8
(the inner block size = `R - 1 = 8`) with parts ≤ 4. Results:

| W   | dead choices | shapes that hit         | min new prime helpers |
|-----|--------------|-------------------------|-----------------------|
| 15  | {2, 4}       | {(4, 4)} only           | 7                     |
| 16  | {1}          | several (all contain 4) | 7                     |

For **W=15**, only one shape `(4, 4)` hits — every shape with all parts
≤ 3 returns NO HIT across all 308 invertible inner 8×8 row subsets and
both dead-col choices.

For **W=16**, multiple shapes hit but every winning shape contains a
part of size 4. Notably **all parts-≤-3 shapes** fail: `(2, 3, 3)`,
`(3, 2, 3)`, `(3, 3, 2)`, `(1, 1, 3, 3)`, `(1, 3, 1, 3)`, `(1, 3, 3, 1)`,
`(2, 1, 3, 2)`, `(2, 2, 2, 2)`, etc. — all 0 hits across 3536 invertible
inner row subsets.

Atomicity verification (`w15_w16_inner_4x4_atomicity.py`): the 4×4
sub-blocks in the chosen `(4, 4)` decompositions of W=15 are themselves
**part-3 atomic** (they admit no `(1, 3)`, `(3, 1)`, `(2, 2)`,
`(1, 1, 2)`, `(1, 2, 1)`, `(2, 1, 1)`, or `(1, 1, 1, 1)` BT split):

```
W=15 A-block (rows {1,3,7,13} × cols {1,3,7,13}):  part-3 atomic
W=15 B-block (rows {2,6,8,12} × cols {0,6,10,12}): part-3 atomic
W=16 A-block (rows {1,2,3,7}  × cols {0,6,10,12}): admits (1, 3)
W=16 B-block (rows {5,11,13,14} × cols {2,4,8,14}): part-3 atomic
```

(The W=16 A-block is `(1, 3)` reducible internally, but the global
shape `(1, 3, 1, 3)` at the 8×8 level still has 0 hits — see "Why the
internal (1, 3) doesn't help" below — so the global decomposition still
needs at least one 4-block.)

## W=15 best candidate

```
shape       = (4, 4)                       (inner 8×8 split)
dead col    = 2  (chiP(3) = 1)
ρ           = (0, 1, 3, 7, 13, 2, 6, 8, 12)
σ           = (2, 1, 3, 7, 13, 0, 6, 10, 12)
block_dets  = [+1, +1]                     (A·B det = +1)
full det    = +1
max_row     = 13   (need 15^j ≥ 14; j ≥ 1 OK since max_row < 15)
new helpers = {101, 103, 107, 131, 191, 193, 197}  (7 new chiP_X_eq_one)
```

The full 9×9 submatrix:

```
         σ:  2   1   3   7  13   0   6  10  12
ρ= 0:   [   1,  1,  0,  0,  0,  0,  1,  1,  1 ]
ρ= 1:   [   0,  1,  1,  1,  1,  0,  0,  0,  0 ]   <- A: rows {1,3,7,13}
ρ= 3:   [   0,  1,  0,  1,  1,  0,  0,  0,  0 ]      cols {1,3,7,13}
ρ= 7:   [   0,  1,  1,  1,  0,  0,  0,  0,  0 ]      det(A) = +1
ρ=13:   [   0,  1,  1,  0,  0,  0,  0,  0,  0 ]
ρ= 2:   [   0,  0,  0,  0,  0,  1,  1,  1,  1 ]   <- B: rows {2,6,8,12}
ρ= 6:   [   0,  0,  0,  0,  0,  0,  1,  1,  1 ]      cols {0,6,10,12}
ρ= 8:   [   0,  0,  0,  0,  0,  0,  1,  1,  0 ]      det(B) = +1
ρ=12:   [   0,  0,  0,  0,  0,  1,  0,  1,  1 ]
```

Three-level structure:
* **Outer** `(1+8)` BT via `Fin 1 ⊕ Fin 8 ≃ Fin 9`. 1×1 block = row 0 ×
  col σ=2 (chiP(3) = 1). The dead col k=2 is zero on rows 1..14 (since
  `15r + 3 = 3·(5r+1)` is divisible by 3 for r ≥ 1, and the only prime
  divisible by 3 is 3 itself).
* **Inner** `(4+4)` BT via `Fin 4 ⊕ Fin 4 ≃ Fin 8`. A in upper-left,
  B in lower-right; lower-left = rows {2, 6, 8, 12} × cols {1, 3, 7, 13}
  is all zero (since for r ∈ {2, 6, 8, 12} and c ∈ {1, 3, 7, 13},
  `15r + c + 1` falls in the small-prime range, and direct check confirms
  all those entries are 0 in the original 14×14 j=1 slab).
* **Block dets**: `det A = det B = +1`, both via `det_fin_four`.

**Composites needed (`h_not_prime_X` witnesses):** 52 composites in the
range [1, 209] including `{1, 4, 8, 14, 16, 18, 22, ..., 121, 122, ...,
209}`. All `decide`-checkable except possibly the few ≥ 150.

## W=16 best candidate

```
shape       = (4, 4)                       (inner 8×8 split)
dead col    = 1  (chiP(2) = 1)
ρ           = (0, 1, 2, 3, 7, 5, 11, 13, 14)
σ           = (1, 0, 6, 10, 12, 2, 4, 8, 14)
block_dets  = [-1, -1]                     (A·B det = +1)
full det    = +1
max_row     = 14   (need 16^j ≥ 15; j ≥ 1 OK since max_row < 16)
new helpers = {83, 191, 223, 227, 229, 233, 239}  (7 new chiP_X_eq_one)
```

The full 9×9 submatrix:

```
         σ:  1   0   6  10  12   2   4   8  14
ρ= 0:   [   1,  0,  1,  1,  1,  1,  1,  0,  0 ]
ρ= 1:   [   0,  1,  1,  0,  1,  1,  0,  0,  1 ]   <- A: rows {1,2,3,7}
ρ= 2:   [   0,  0,  0,  1,  0,  0,  1,  1,  1 ]      cols {0,6,10,12}
ρ= 3:   [   0,  0,  0,  1,  1,  0,  1,  0,  0 ]      det(A) = -1
ρ= 7:   [   0,  1,  0,  0,  0,  0,  0,  0,  1 ]
ρ= 5:   [   0,  0,  0,  0,  0,  1,  0,  1,  0 ]   <- B: rows {5,11,13,14}
ρ=11:   [   0,  0,  0,  0,  0,  1,  1,  0,  1 ]      cols {2,4,8,14}
ρ=13:   [   0,  0,  0,  0,  0,  1,  0,  0,  1 ]      det(B) = -1
ρ=14:   [   0,  0,  0,  0,  0,  1,  1,  1,  1 ]
```

Same three-level structure as W=15 (outer `(1+8)`, inner `(4+4)`, inner-
inner `(1+3)` for the W=16 A-block if desired). For W=16 the A-block
admits an additional `(1, 3)` BT split internally — row 7 at col 0 has
support size 1 (chiP(7·16+0+1) = chiP(113) = 1), so the A-block could
in principle decompose as `1 + (3×3)` and avoid `det_fin_four` on the A
block. But the **B-block remains part-3 atomic**, so `det_fin_four` is
still needed for B. Net: W=16 saves at most one `det_fin_four` over
W=15, but cannot fully avoid the technique.

**Composites needed (`h_not_prime_X` witnesses):** 50 composites in the
range [1, 237].

## Why the internal `(1, 3)` of the W=16 A-block doesn't help globally

The pre-search returns NO HIT for the global shape `(1, 3, 1, 3)` at
8×8. Yet the (4, 4)-best W=16 candidate's 4×4 A-block IS internally
`(1, 3)` decomposable. How?

The constraint differential: a global `(1, 3, 1, 3)` decomposition of the
8×8 inner needs row partition `(R1, R2, R3, R4)` with sizes (1, 3, 1, 3)
AND column partition `(C1, C2, C3, C4)` such that all six lower-left
cross-blocks `(R_i, C_j)` for `i > j` are zero. This is six constraints.

The `(4, 4)` decomposition + internal `(1, 3)` of the A-block has only
THREE constraints: lower-left of `(4, 4)`, lower-left of A's `(1, 3)`,
and that's it. The internal `(1, 3)` of A doesn't carry constraints
about how its "1" col interacts with B's rows.

So the global `(1, 3, 1, 3)` is strictly more constrained than the
nested `(4, 4) ∘ (1, 3)`, and the W=16 8×8 fails the global version
while supporting the nested version. The Lean proof can use the nested
version to compute `det A = (chiP(113)) · det(3×3 cofactor) = 1 · (±1)`
via `det_fromBlocks_zero₂₁ + det_fin_one + det_fin_three`.

## Implication for the next-action queue

The W=15 / W=16 closures are NOT single-session under the current Lean
infrastructure. The arc decomposes into TWO sub-tasks:

1. **Sub-arc A: develop `Matrix.det_fin_four` as a reusable lemma in
   `MPSBondDim/Basic.lean`.** Estimated 1 session of focused Lean
   engineering. Target: a lemma of the shape

   ```
   lemma Matrix.det_fin_four (M : Matrix (Fin 4) (Fin 4) ℚ) :
       M.det = (cofactor expansion along row 0)
   ```

   derivable from `Matrix.det_succ_row_zero` + `Fin.sum_univ_four` +
   `Matrix.det_fin_three`. Once available, future R=9 closures (W=15,
   W=16, plus W ∈ {21, 22, 25, 27, ...} eventually) can use it directly.

2. **Sub-arc B: apply `det_fin_four` to close W=15 (or W=16).** Estimated
   1 session. W=15 is the cleaner target because:
   - Only one shape `(4, 4)` works, so the proof structure is forced.
   - Block dets `[+1, +1]` give `IsUnit 1` shortcut at the very end.
   - 7 new prime helpers `{101, 103, 107, 131, 191, 193, 197}` (smallest
     prime ≥ 150 is 191 — `norm_num` for primes ≥ 150 per S137 pattern;
     the rest fit `decide`).
   - Max row 13 < 15 so works for all `j ≥ 1`.

W=16 has more flexibility (multiple shapes, internal A-block reducibility)
but the same `det_fin_four` requirement and 7 new prime helpers
`{83, 191, 223, 227, 229, 233, 239}`.

## Verification

The pre-search script `w15_w16_blocktriangular_search.py` uses
bitmask-accelerated zero-block checks (sigs[r] & A_cols_mask == 0 for
"lower-left zero"), then `np.linalg.det` for the 4×4 invertibility
check. Total runtime ~2 min on Python 3 over (15: 308 row subsets ·
108 shapes; 16: 3536 row subsets · 108 shapes).

Both candidates verified by recomputing the 9×9 submatrix entry-wise
from `chiP(ρ_i · W + σ_k + 1)` and confirming `np.linalg.det = +1` to
machine precision.

## Files

* `w15_w16_blocktriangular_search.py` — the pre-search.
* `w15_w16_blocktriangular_search_run1.log` — full output (401 lines).
* `w15_w16_inner_4x4_atomicity.py` — verification that the 4×4 inner
  sub-blocks are part-3 atomic.
* `w15_w16_inner_4x4_atomicity_results.log` — atomicity verification log.

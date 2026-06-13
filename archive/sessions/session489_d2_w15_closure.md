# Session 489 — Arc 2 / Sub-arc D-2: W=15 orthogonal-corner closure

**Date:** 2026-05-01.
**Mode:** ARC CONTINUATION (Arc 2 — Lean Formalisation Track).
**Run #** 491 (per `.run_state` post-session).
**Self-graded:** **B-grade** (substantive refinement of an existing
edge: another `mps_bond_dim` orthogonal-corner instance closed
unconditionally, sorry-free, axiom-pure. First corner closure to use the
S423 `det_fin_four` primitive — landing the planned Sub-arc D-2 from
the S245 / S423 next-action queue. Not A-grade because no new
mathematical content beyond a 14th specialisation of the same
identity, with the technique (nested `det_fromBlocks_zero₂₁`) already
established at S152/S159/S235.)

## Mission

Per RESEARCH_AGENDA.md Arc 2 next-action queue (post-S423): execute the
planned Sub-arc D-2 — close the orthogonal corner case `(W = 15, d = j +
1)` of E2.1 unconditionally, using the S245 pre-search candidate and
the S423 `det_fin_four` primitive. Statement target:

```
mps_bond_dim_W_eq_15_d_eq_j_plus_1 : ∀ j ≥ 1, (unfolding 15 (j+1) j).rank = 9
```

## What was produced

> **Theorem (S489, sorry-free, axiom-pure).**
> For every `j ≥ 1`, `(unfolding 15 (j+1) j).rank = 9`.
>
> `#print axioms` shows only `propext, Classical.choice, Quot.sound`.

**Location:** `MPSBondDim/MPSBondDim/Basic.lean` lines ~5933–6717, just
before `end E2_1`. Eight new declarations:

* 7 prime helpers: `chiP_one_hundred_one_eq_one`,
  `chiP_one_hundred_three_eq_one`, `chiP_one_hundred_seven_eq_one`,
  `chiP_one_hundred_thirty_one_eq_one`,
  `chiP_one_hundred_ninety_one_eq_one`,
  `chiP_one_hundred_ninety_three_eq_one`,
  `chiP_one_hundred_ninety_seven_eq_one` (all via `norm_num`).
* `exists_invertible_submatrix_W_eq_15_d_eq_j_plus_1`: the prime exhibit.
* `mps_bond_dim_W_eq_15_d_eq_j_plus_1`: the rank equality.

## Proof structure

* **Permutation** (from S245 pre-search):
  - `ρ : Fin 9 → Fin (15^j)` mapping `(0..8) ↦ (0, 1, 3, 7, 13, 2, 6, 8, 12)`.
  - `σ : Fin 9 → Fin (15^((j+1)-j))` mapping `(0..8) ↦ (2, 1, 3, 7, 13, 0, 6, 10, 12)`.
  - Max `ρ.val = 13 < 15 ≤ 15^j` for `j ≥ 1`. Max `σ.val = 13 < 15`.

* **Matrix** (9×9 `chiP` submatrix, equals to the explicit `Mexp`):
  ```
            σ:  2   1   3   7  13   0   6  10  12
  ρ= 0:   [    1,  1,  0,  0,  0,  0,  1,  1,  1 ]
  ρ= 1:   [    0,  1,  1,  1,  1,  0,  0,  0,  0 ]   <- A: rows {1,3,7,13}
  ρ= 3:   [    0,  1,  0,  1,  1,  0,  0,  0,  0 ]      cols {1,3,7,13}
  ρ= 7:   [    0,  1,  1,  1,  0,  0,  0,  0,  0 ]      det(A) = +1
  ρ=13:   [    0,  1,  1,  0,  0,  0,  0,  0,  0 ]
  ρ= 2:   [    0,  0,  0,  0,  0,  1,  1,  1,  1 ]   <- D: rows {2,6,8,12}
  ρ= 6:   [    0,  0,  0,  0,  0,  0,  1,  1,  1 ]      cols {0,6,10,12}
  ρ= 8:   [    0,  0,  0,  0,  0,  0,  1,  1,  0 ]      det(D) = +1
  ρ=12:   [    0,  0,  0,  0,  0,  1,  0,  1,  1 ]
  ```

* **Two-level nested `det_fromBlocks_zero₂₁`**:
  - **Layer 1 (outer, 1+8)** via `finSumFinEquiv : Fin 1 ⊕ Fin 8 ≃ Fin 9`:
    splits off row 0 + dead-col witness `chiP 3 = 1`. Lower-left
    `1..8 × {col 2}` is zero because `15r + 3 = 3·(5r+1) ≠ 3` is
    composite for `r ≥ 1`. `A_outer.det = 1` via `det_fin_one`.
  - **Layer 2 (inner, 4+4)** via `finSumFinEquiv : Fin 4 ⊕ Fin 4 ≃ Fin 8`:
    splits the 8×8 inner D_outer block into a top-left 4×4 (rows
    {1,3,7,13} × cols {1,3,7,13}) with all small primes, and a
    bottom-right 4×4 (rows {2,6,8,12} × cols {0,6,10,12}) with mid-range
    primes. Lower-left 4×4 is zero (the underlying chiP values for
    rows {2,6,8,12} × cols {1,3,7,13} are all composite — verified via
    52 non-primality witnesses).
  - **Inner block determinants** (computed via `simp [det_fin_four]`):
    - `A_in.det = 1` for `!![1, 1, 1, 1; 1, 0, 1, 1; 1, 1, 1, 0; 1, 1, 0, 0]`.
    - `D_in.det = 1` for `!![1, 1, 1, 1; 0, 1, 1, 1; 0, 1, 1, 0; 1, 0, 1, 1]`.
  - **Total determinant**: `1 · (1 · 1) = 1`, hence `IsUnit`.

* **Upper-bound subtlety**: the matrix has `15^((j+1)-j) = 15` columns,
  so `rank_le_width` gives only `rank ≤ 15`. We instead cite the general
  `upper_bound`, which evaluates to `φ(15) · 15^0 + 1 = 8 · 1 + 1 = 9`.

## Key debugging milestones

1. **First build attempt (`maxHeartbeats 2000000`)**: timeout at
   `isDefEq` on the `h_fromBlocks_outer` step (line 6660,
   `all_goals fin_cases i <;> fin_cases j <;> rfl`). The 1×8 + 8×1 +
   8×8 = 81 sub-cases each requiring `rfl` exceeds the 10× budget.
   Bumped to `4000000` (20×) — succeeded.
2. **Second build attempt (`maxHeartbeats 4000000`)**: errors at
   `det_fin_four` applications for `A_in` and `D_in`. The pattern
   `rw [det_fin_four]; norm_num` left unsolved goals because `norm_num`
   doesn't unfold `(!![...]) i j` to specific entries. Replaced with
   `simp [det_fin_four]` (parallel to W=14's `simp [Matrix.det_fin_two]`
   pattern) — succeeded.
3. **Third build (`maxHeartbeats 4000000`, `simp [det_fin_four]`)**:
   build succeeded with only warnings (unused simp args from the matrix
   equality, which are inherited from the W=14 template). `#print
   axioms` confirms axiom-pure.

## Why not A-grade

This closure adds one more wheel `W = 15` to the closed-W set
`{2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 18, 20}` — now `{2, 3, 4, 5, 6, 7,
8, 9, 10, 12, 14, **15**, 18, 20}`. The technique (nested
`det_fromBlocks_zero₂₁`) is from S152/S159/S235 verbatim; the
determinant primitive (`det_fin_four`) is from S423; the (ρ, σ)
candidate is from S245. The session executed a planned single-session
target with a known shape — substantive rigor work but not novel
mathematical content. **Honest grading: B**, not A.

A-grade would have required either:
(a) a closure with a *new* technique (e.g., `det_fin_five` for W=11),
    or
(b) a new structural insight beyond the closure itself (e.g., a unified
    framework spanning multiple wheels), or
(c) closing a previously-open mathematical question rather than a
    pre-planned formalisation step.

## Self-evaluation (per CLAUDE.md)

1. **What was produced that wasn't in the project before?**
   - The W=15 mps_bond_dim instance is new (Lean theorem,
     `mps_bond_dim_W_eq_15_d_eq_j_plus_1`).
   - 7 new prime helpers `chiP_X_eq_one` for X ∈ {101, 103, 107, 131,
     191, 193, 197}.
   - First Lean corner closure to use the S423 `det_fin_four` primitive.
   - Empirical fact: `set_option maxHeartbeats 4000000` is the
     necessary budget for `R = 9` matrix equality at fromBlocks.
2. **What edges did the work compose / cite?**
   - E2.1 (the MPS bond-dim identity itself). The closure is the W=15
     specialisation of that edge.
   - Sub-arc D-1 (S423): cites `det_fin_four` directly.
   - S152/S159/S235: cites the nested `det_fromBlocks_zero₂₁` template.
   - S245: cites the pre-search permutation and shape uniqueness.
3. **Why was the session not duplicate?** It was a planned closure of a
   pre-staged target. Sub-arc D-2 was the explicit next-action listed
   in RESEARCH_AGENDA.md after S423.
4. **Next-action for the next agent (single-session)**:
   - **Option A: Sub-arc D-3 — close W=16.** Permutation (from S245):
     `ρ = (0, 1, 2, 3, 7, 5, 11, 13, 14)`,
     `σ = (1, 0, 6, 10, 12, 2, 4, 8, 14)`, dead col 1, block dets
     `[-1, -1]`, full det `+1`. 6 new prime helpers
     `{83, 223, 227, 229, 233, 239}` (191 already declared at S489).
     Estimated ~860 Lean lines (similar to W=15). The W=16 A-block
     can use inner `(1, 3)` fromBlocks instead of `det_fin_four`,
     saving one `det_fin_four` invocation, but B-block still needs it.
   - **Option B: Sub-arc D-4 — develop `det_fin_five`** for W=11.
     The S423 D-1 proof is a working template (cofactor expansion via
     `det_succ_row_zero` + per-cofactor `det_fin_four` + 16
     `decide`-checked `Fin.succAbove` resolutions for the
     `q ∈ {2, 3}` regime + `ring` on 120 monomials). Estimated 1
     session for `det_fin_five`, then 1 more for W=11 closure
     (Path B from S206).

## File touch summary

* `experiments/formalisations/E2_1_mps_bond_dim/MPSBondDim/MPSBondDim/Basic.lean`:
  +860 lines (W=15 closure + 7 prime helpers).
* `RESEARCH_AGENDA.md`: Arc 2 status updated to "IN PROGRESS — Run #491 /
  Session 489 (Sub-arc D-2: closed W=15 ...)"; W=15 milestone added
  with full proof description; closed-W set incremented from 13 to 14;
  next-action queue updated to point at Sub-arc D-3 (W=16) or Sub-arc
  D-4 (`det_fin_five` for W=11).
* `experiments/formalisations/E2_1_mps_bond_dim/mps_bond_dim_notes.md`:
  S489 paragraph appended after the S423 paragraph.

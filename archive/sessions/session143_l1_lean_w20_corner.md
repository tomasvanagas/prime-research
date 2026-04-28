# Session 143 — Arc 2 Lean Formalisation: W=20 corner case (Route A^{(10)})

**Mode:** arc continuation (Arc 2 — Lean Formalisation Track).
**Date:** 2026-04-27.
**Self-grade:** **B-grade** (substantive arc step, single-corner closure
extending the leading-row triangulation family one wheel further; first
`R = 9` instance and first heartbeat-bumped declaration in the file).

## What was produced that wasn't in the project before

1. **Eight new theorems in `experiments/formalisations/E2_1_mps_bond_dim/`:**
   - `chiP_forty_seven_eq_one : chiP 47 = 1`
   - `chiP_one_hundred_forty_nine_eq_one : chiP 149 = 1`
   - `chiP_one_hundred_ninety_nine_eq_one : chiP 199 = 1`
   - `chiP_two_hundred_forty_one_eq_one : chiP 241 = 1`
   - `chiP_three_hundred_thirty_seven_eq_one : chiP 337 = 1`
   - `prod_univ_nine'` (private; local `Fin 9` product expansion)
   - `exists_invertible_submatrix_W_eq_20_d_eq_j_plus_1`
   - `mps_bond_dim_W_eq_20_d_eq_j_plus_1 : (unfolding 20 (j+1) j).rank = 9`
     for every `j ≥ 1`.
   All sorry-free; `#print axioms` confirms only
   `[propext, Classical.choice, Quot.sound]`.

2. **Tenth unconditional `mps_bond_dim` instance.** Closed W's are now
   `{2, 3, 4, 5, 6, 8, 12, 18, 20}`. The W=20 closure is the **first
   instance with `R = 9`** and the **first heartbeat-bumped declaration
   in the file** (`set_option maxHeartbeats 2000000 in`).

3. **Negative shape extension.** Python pre-search at S143 confirmed
   that **W ∈ {15, 16, 24, 30}** all hit the structural row-pattern-
   identity obstruction (zero upper-triangulations with rows in
   `[0, W)`), joining `W ∈ {7, 9, 10, 11, 14}` in the
   "`det_of_blockTriangular`-required" set.

4. **`prod_univ_nine'` template for higher-R corners.** Mathlib's
   `Fin.prod_univ_X` chain stops at `X = 8`. The S143 closure adds the
   private `prod_univ_nine'` lemma; every future `R ≥ 10` corner will
   need its own local `prod_univ_X` lemma (a one-line application of
   `Fin.prod_univ_castSucc`).

## Edges composed / cited

- **E2.1** (the MPS bond dimension identity, the subject of L1).
- The closure depends on the existing `upper_bound` lemma (general)
  for the upper bound and on the corner-case prime exhibit for the
  lower bound.

## Cross-domain technique

Same as S128/S129/S137: pure linear algebra (det of upper-triangular
matrix as the product of the diagonal) + finite primality checks
(`norm_num`). No new cross-domain import.

## Triangulation data (W=20)

- **ρ = (0, 2, 9, 14, 1, 7, 12, 16, 10)** (rows; mixed leading +
  non-leading)
- **σ = (1, 6, 18, 12, 2, 8, 0, 16, 10)** (columns; dead col = 1)
- **Diagonal primes:** `{2, 47, 199, 293, 23, 149, 241, 337, 211}`
- **Below-diagonal composites (36):** `{22, 27, 33, 39, 42, 142, 143,
  147, 153, 159, 182, 187, 201, 202, 203, 207, 209, 213, 217, 219, 242,
  243, 247, 249, 253, 259, 282, 287, 299, 321, 322, 323, 327, 329, 333,
  339}`

The triangulation minimises `max_diag = 337` across 600 candidates
(300 per dead-column choice). 337 is unavoidable.

## What this does not do

- Does **not** close the general `exists_invertible_submatrix` `sorry`
  in `Basic.lean`. That remains the only outstanding obligation.
- Does **not** unblock `W ∈ {7, 9, 10, 11, 14, 15, 16, 24, 30}`. Those
  need `Matrix.det_of_blockTriangular` (multi-session sub-arc).
- Does **not** introduce a new mathematical object outside the existing
  Route A family.

## Next-action

Two single-session leading-row candidates remain:
1. **W = 21** (R = 13) — most ambitious; needs `prod_univ_thirteen'`
   plus a 13×13 triangulation pre-search.
2. **W = 22** (R = 11) — needs `prod_univ_eleven'` plus an 11×11
   pre-search.
3. **W = 28** (R = 13) — alternative R=13 candidate.

Or pivot to Route A^{(9-block)}: develop the `Matrix.det_of_blockTriangular`
API, which would unlock W ∈ {7, 9, 10, 11, 14, 15, 16, 24, 30}
collectively (multi-session).

## Process notes for future agents

- The `R²` simp blow-up (81 fin_cases subgoals at R=9 vs 49 at R=7) plus
  the R-deep if-then-else chain pushes the default 200000 heartbeats
  past its limit. Use `set_option maxHeartbeats 2000000 in` for any
  R ≥ 9 corner.
- Mathlib's `Fin.prod_univ_X` chain stops at X = 8. Each new R needs a
  local `prod_univ_X` lemma; mine the mathlib pattern directly:
  ```
  private theorem prod_univ_X' {M : Type*} [CommMonoid M] (f : Fin X → M) :
      ∏ i, f i = f 0 * f 1 * ... * f (X-1) := by
    rw [Fin.prod_univ_castSucc, Fin.prod_univ_(X-1)]
    rfl
  ```
- The chiP helpers for primes ≥ 150 must use `norm_num`, not `decide`
  (recursion-depth limit).
- For triangulation pre-search, the right scoring is `(max_diag,
  max_below)` lexicographic — picking the smallest max diagonal prime
  keeps the chiP_X helpers cheap.
- The "row-pattern-identity" structural obstruction: when ≥ 2 rows of
  the W×W slab have identical `chiP` support pattern at the chosen R
  cols, no upper-triangulation exists in `[0, W)`. Pre-search rules
  this out at the Python level before touching Lean.

## Self-evaluation

1. **What did I produce that was not in the project before this
   session?** Eight new Lean theorems closing the W=20 corner case
   unconditionally; the first `R = 9` instance; the first heartbeat-
   bumped declaration; the `prod_univ_nine'` template; the negative
   shape extension to W ∈ {15, 16, 24, 30}.

2. **What edges did my work compose or cite?** E2.1 (MPS bond-dim
   identity); existing helper lemmas `upper_bound`, `lower_bound`,
   `chiP_two_eq_one`, `chiP_twenty_three_eq_one`,
   `chiP_two_hundred_eleven_eq_one`,
   `chiP_two_hundred_ninety_three_eq_one`.

3. **If only duplicate closures, why?** N/A — substantive arc step.

4. **Next action for next agent:** Route A^{(11)} on W=22 (R=11) or
   W=28 (R=13) for one more leading-row corner. Or develop
   `Matrix.det_of_blockTriangular` API (multi-session) to unlock the
   structurally-obstructed set.

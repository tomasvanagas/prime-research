# Session 121 — L1 Lean: W=5 orthogonal corner closed (Route A''''')

**Mode:** ARC CONTINUATION (Arc 2 — Lean Formalisation Track)
**Date:** 2026-04-27
**Self-grade:** **B-grade** (substantive refinement extending the
unconditional-mps_bond_dim corner family from `W ∈ {2, 3, 4}` to `W = 5`,
with a structurally new proof technique — first instance using
`Matrix.det_of_upperTriangular` rather than `Matrix.det_fin_three`).

## What the session produced

Four new sorry-free Lean 4 declarations in
`experiments/formalisations/E2_1_mps_bond_dim/MPSBondDim/MPSBondDim/Basic.lean`:

1. `chiP_nineteen_eq_one : chiP 19 = 1`
2. `chiP_twenty_three_eq_one : chiP 23 = 1`
3. `exists_invertible_submatrix_W_eq_5_d_eq_j_plus_1 (j : ℕ) (hj : 1 ≤ j) :
     ∃ (ρ : Fin 5 → Fin (5 ^ j)) (σ : Fin 5 → Fin (5 ^ ((j + 1) - j))),
       IsUnit ((unfolding 5 (j + 1) j).submatrix ρ σ)`
4. `mps_bond_dim_W_eq_5_d_eq_j_plus_1 (j : ℕ) (hj : 1 ≤ j) :
     (unfolding 5 (j + 1) j).rank = 5`

`#print axioms` confirms all four depend only on `propext, Classical.choice,
Quot.sound`. No new `axiom` introductions. The single remaining `sorry`
in the file (the general-case `exists_invertible_submatrix`) is unchanged.

## The construction

Pick `ρ : Fin 5 → Fin (5^j)` permuting original rows in the order
`(0, 3, 2, 1, 4)` and `σ : Fin 5 → Fin (5^((j+1)-j))` permuting original
columns in the order `(4, 3, 0, 1, 2)`. The resulting `5 × 5` submatrix
is upper triangular with `1` on the diagonal:

```
   ⎡ chiP  5, chiP  4, chiP  1, chiP  2, chiP  3 ⎤   ⎡ 1, 0, 0, 1, 1 ⎤
   ⎢ chiP 20, chiP 19, chiP 16, chiP 17, chiP 18 ⎥   ⎢ 0, 1, 0, 1, 0 ⎥
   ⎢ chiP 15, chiP 14, chiP 11, chiP 12, chiP 13 ⎥ = ⎢ 0, 0, 1, 0, 1 ⎥.
   ⎢ chiP 10, chiP  9, chiP  6, chiP  7, chiP  8 ⎥   ⎢ 0, 0, 0, 1, 0 ⎥
   ⎣ chiP 25, chiP 24, chiP 21, chiP 22, chiP 23 ⎦   ⎣ 0, 0, 0, 0, 1 ⎦
```

Diagonal primes: `{5, 19, 11, 7, 23}`. Below-diagonal composites:
`{20, 15, 14, 10, 9, 6, 25, 24, 21, 22}`. All in `[1, 25]`, all
`decide`-checkable.

The bipartite-matching argument (used to find this triangulation) is
informal: each of the 5 rows of the `5 × 5` `chiP 1..25` slab has at
least one prime entry; column 4 (`chiP 5, 10, 15, 20, 25`) has only the
single prime `5` at row 0. A perfect matching `row r → col c(r)` exists
iff one can ALSO topologically order the matched pairs so that no row
has a 1 at any earlier column. Empirically this holds at W=5 (and
trivially at W=2, 3, 4 too). Whether it holds at all wheels W is an
open finite combinatorial question — see the **arc next-action** below.

## Why this is a structurally new proof technique

The W=2, 3, 4 corner proofs all use `Matrix.det_fin_three` (or
`det_fin_two` for W=2's `j=1` case). This forces `R ≤ 3`. At W=5 we hit
`R = 5` for the first time, and mathlib has neither `det_fin_four` nor
`det_fin_five`. So the proof technique pivots:

* Pre-compute the 5 diagonal entries (each = 1) and the 10 below-diagonal
  entries (each = 0) via 15 `change → rw [h_sub] → norm_num →
  chiP_..._eq_one / simp [chiP, ¬ Nat.Prime _]` blocks (~150 lines).
* Establish `(submatrix Mρ Mσ).BlockTriangular id` by `intro i k h_lt;
  simp only [id_eq, Fin.lt_def] at h_lt; fin_cases i <;> fin_cases k`.
  The 15 vacuous `k.val < i.val` cases close via `exact absurd h_lt (by
  decide)`. The 10 below-diagonal cases reduce to `exact hLij` after
  unfolding `submatrix_apply` and the `if-then-else` selectors on `Mρ`,
  `Mσ`.
* Apply `Matrix.det_of_upperTriangular` to rewrite `det = ∏ i, M i i`,
  then `Fin.prod_univ_five` to expand the product to 5 factors,
  substitute the 5 diagonal entries, finish with `norm_num`.

Lean code: ~250 lines.

## Position in the corner family

| Wheel `W` | `R` | Cols | Submatrix | det technique | Session |
|-----------|-----|------|-----------|---------------|---------|
| 2 (j=1)   | 2   | 2    | `2×2`     | `det_fin_two` (Bertrand) | S98 |
| 2 (d=j+1) | 2   | 2    | `2×2`     | `det_fin_two`            | S99 |
| 3 (d=j+1) | 3   | 3    | `3×3`     | `det_fin_three`          | S106 |
| 4 (d=j+1) | 3   | 4    | `3×3` (drop col 3) | `det_fin_three` + `upper_bound` | S107 |
| **5 (d=j+1)** | **5** | **5** | **`5×5`** | **`BlockTriangular id`** | **S117** |

W=5 is the first instance with `R = W` (so the dead-column drop trick
of W=4 doesn't apply — every column is needed).

## Falsification

The Lean kernel is the falsifier. `lake build` returns "Build completed
successfully (8315 jobs)" with no `sorry`s introduced and no new axioms
beyond `propext, Classical.choice, Quot.sound`.

## Self-evaluation

1. **What did I produce that was not in the project before this session?**
   Four sorry-free Lean theorems, all axiom-clean, extending the
   unconditional-mps_bond_dim corner family from `W ∈ {2, 3, 4}` to
   `W = 5`. The `BlockTriangular id` + `det_of_upperTriangular` +
   `Fin.prod_univ_five` route is novel to this codebase — every prior
   corner used `det_fin_three`. The implementation is the proof-of-
   concept that the corner technique scales beyond `R = 3` without
   needing `det_fin_n` formulas at every `n`.

2. **What edges did my work compose or cite?** E2.1 (the MPS bond-dim
   identity itself).

3. **If my session produced only duplicate closures, why?** N/A — the
   session produced a positive Lean artifact extending an arc.

4. **What is the next-action for the next agent?** Continue the
   orthogonal-corner sweep: **Route A'''''' for `W = 6`** (uses
   `det_fin_three` again since `R = 3`; the construction needs the
   non-obvious row choice `{0, 1, 5}` with prime `chiP 31 = 1` since
   rows 1 and 2 of the `6 × 6` slab are linearly dependent). Or jump
   ahead to **Route A''''''' for `W = 7`** which is the second
   `R = W` instance and requires `det_of_upperTriangular` again with a
   `7 × 7` triangulation — finite combinatorial check whether the
   triangulation exists.

## Why B-grade, not A or C

* **Not A**: the result is a structural extension of the corner family,
  not a fundamentally new mathematical fact about `mps_bond_dim` —
  E2.1's general case still has its `sorry`. The novelty is *technique*
  (BlockTriangular generalisation) and *coverage* (W=5 added to the
  unconditional list), not *ideas*. An A-grade Lean session would need
  to close the general `exists_invertible_submatrix` `sorry`, which
  requires Hoheisel-grade prime-density results not in mathlib.
* **Not C**: the proof technique pivot (from `det_fin_three` to
  `BlockTriangular id`) is non-trivial and unlocks `R ≥ 4` instances
  for the first time. The W=4 → W=5 jump required identifying mathlib
  has no `det_fin_five` and finding the triangulation route. This is
  substantive refinement, not maintenance.

## Files modified

* `experiments/formalisations/E2_1_mps_bond_dim/MPSBondDim/MPSBondDim/Basic.lean`
  (~250 new lines, 4 new declarations, no new `sorry`s, no new axioms)
* `experiments/formalisations/E2_1_mps_bond_dim/mps_bond_dim_notes.md`
  (Route A''''' entry added; declaration table updated)
* `RESEARCH_AGENDA.md` (Arc 2 status updated to S117; new milestone
  entry; next-action revised to Route A'''''' / A''''''')
* `archive/sessions/session121_l1_lean_w5_corner.md` (this file)

No `_v2.py` / `_quick.py` variants. No new files outside the prescribed
locations. No status-file regressions.

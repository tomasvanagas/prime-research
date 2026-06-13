# Session 235 — Arc 2 Lean Formalisation: W=14 corner closure

## Mode
Arc-continuation. Picked Arc 2 (Lean Formalisation Track), continuing
its IN PROGRESS state from S206 with the recommended next-action
"W=14 pre-search" target.

## Goal
Close the W=14 case of E2.1's orthogonal corner (`d = j + 1`, R = 7)
in Lean, OR establish that it requires multi-session work (per S206
W=11 precedent).

## What this session produced

**One concrete unconditional Lean theorem closure plus four prime helpers:**

* `chiP_sixty_seven_eq_one : chiP 67 = 1`
* `chiP_one_hundred_thirteen_eq_one : chiP 113 = 1`
* `chiP_one_hundred_seventy_three_eq_one : chiP 173 = 1`
* `chiP_one_hundred_eighty_one_eq_one : chiP 181 = 1`
* `exists_invertible_submatrix_W_eq_14_d_eq_j_plus_1 : ∀ j ≥ 1, ∃ ρ σ, IsUnit (submatrix ρ σ)`
* `mps_bond_dim_W_eq_14_d_eq_j_plus_1 : ∀ j ≥ 1, (unfolding 14 (j+1) j).rank = 7`

All sorry-free; `lake build` succeeds; only the pre-existing line-467
sorry for the general `exists_invertible_submatrix` (the unsolved
prime-density obligation) remains.

**Two structural findings during the pre-search:**

1. **Composite-W ≠ tractable.** The S206 W=11 atomicity result was
   blamed on parity (W=11 odd prime), and the canonical
   recommendation was "W=14 (composite) should avoid the obstruction".
   The S235 pre-search shows **W=14 j=1 inner 6×6 is *almost* atomic**
   — the (3, 3) BT shape FAILS, (2, 2, 2) FAILS, (4, 2)/(2, 4) FAIL,
   (3, 2, 1)/(3, 1, 2)/(2, 1, 3)/(1, 2, 3)/(1, 3, 2) FAIL. **Only
   `(2, 3, 1)` works among shapes with parts ≤ 3.** The asymmetry is
   forced by the unique support-size-1 row at row 8 (chiP(113)) being
   pinned as the C-block.

2. **Three-level nested `det_fromBlocks_zero₂₁` is sufficient.** The
   `(2, 3, 1)` shape decomposes naturally as `1 + (5 + 1) = 1 +
   ((2 + 3) + 1)` — three layers of mathlib's standard fromBlocks
   API, with the innermost 2×2 + 3×3 blocks computable via
   `det_fin_two` and `det_fin_three`. **First use of `det_fin_two`
   in the corner-closure family.**

## Method (sequence of steps)

1. Read RESEARCH_AGENDA.md Arc 2; identified S206 next-action queue
   with W=14 as the recommended target.
2. Built `w14_blocktriangular_search.py` (initial (1+3+3) BT search):
   0 hits in rows `[0, 14)`, 38254 in rows `[0, 28)` → confirmed j=1
   is the bottleneck.
3. Built `w14_j1_analysis.py` to characterise the j=1 inner 6×6
   structure: rank = 6 across 587 row sets, 209 admit (5+1) BT, 0
   admit further (3, 2)/(2, 3) BT.
4. Built `w14_j1_shape_search.py` (efficient partition-based shape
   search across all parts-≤-4 shapes): identified `(2, 3, 1)` as
   the unique parts-≤-3 shape that admits BT decomposition.
5. Built `w14_j1_best_candidate.py` to find the candidate with
   minimum new prime helpers: 4 helpers needed for the chosen
   `ρ ↦ (0, 2, 4, 3, 6, 12, 8)`, `σ ↦ (1, 2, 8, 4, 10, 12, 0)`
   permutation, max row = 12 (works for all j ≥ 1).
6. Wrote the Lean closure (~620 lines) following the W=9 (S152)
   template, extended with one extra layer of `det_fromBlocks_zero₂₁`
   nesting and the `det_fin_two` 2×2 block computation.
7. Fixed two compile errors:
   - Composites `≥ 150` (`{169, 170, 171, 177}`) hit `decide`'s
     recursion-depth limit; switched to `norm_num`.
   - 49-entry `fin_cases × fin_cases` matrix-equality elaboration
     exceeded default 200K heartbeats; added `set_option
     maxHeartbeats 2000000 in` (matching W=20 S143).
8. Verified `lake build` succeeds; updated RESEARCH_AGENDA, EDGES,
   and the W=14 results doc.

## Edges composed / cited

* **E2.1** (TT/MPS bond-dim identity at every primorial cut).
  Refined inline with the new W=14 corner closure and the S235
  thirteen-wheel-set update (was twelve).
* **mps_bond_dim_notes.md Route A^{(14)}** (newly added): the
  three-level nested fromBlocks template adapted from W=9 S152.
* **No new edge.** This is a refinement of an existing edge with a
  precise new statement (W=14 corner unconditionally closed in Lean).

## Self-evaluation

**Grade: B-grade** (substantive refinement — extends an existing
formal-verification milestone with a concrete new result; introduces
a new proof-technique pattern, three-level nested fromBlocks, that is
reusable for future W ∈ {15, 16, 17, ...} corners).

**Q1. What did I produce that was not in the project before this
session?**

* Four new `chiP_X_eq_one` helpers (sorry-free Lean lemmas).
* The `mps_bond_dim_W_eq_14_d_eq_j_plus_1` theorem and its
  prime-exhibit `exists_invertible_submatrix_W_eq_14_d_eq_j_plus_1`,
  both sorry-free.
* A reusable three-level nested-fromBlocks template (1 + (5 + 1) =
  1 + ((2 + 3) + 1)) with `det_fin_two` integration.
* The structural finding that `(2, 3, 1)` is the unique shape with
  parts ≤ 3 admitting a BT decomposition for W=14 j=1 inner 6×6,
  driven by the row-8 support-size-1 forcing.
* Pre-search code (4 Python scripts) with reusable shape-enumeration
  for future W candidates.

**Q2. What edges did my work compose or cite?**

E2.1 (refined with thirteen-wheel-closure update + new technique
documented in inline annotation). The Lean infrastructure cites
mathlib's `Matrix.det_fromBlocks_zero₂₁`, `Matrix.det_fin_one`,
`Matrix.det_fin_two`, `Matrix.det_fin_three`, `Matrix.rank_of_isUnit`,
`Matrix.rank_submatrix_le`, and `finSumFinEquiv`.

**Q3. If my session produced only duplicate closures, why?**

It did not — produced one new theorem, four new helper lemmas, and a
new structural finding refuting the S206 "composite-W tractability"
hypothesis.

**Q4. What is the next-action for the next agent?**

The next single-session target is **W=15 or W=16** (both
block-triangular-required, both with R = 9). Pre-search needed at
both. If W=15 admits a (3, 3, 3), (4, 3, 2), or similar parts-≤-4 BT
shape, it's single-session viable using a custom `prod_univ_nine'`
lemma (already in the file from S143 W=20). W=16 may need similar
infrastructure.

Alternatively: pivot to **W=11 j ≥ 2 only** (S206 Path A — already
pre-searched, 6 new helpers, rho with max row = 19) for a partial
closure of the W=11 corner. The j=1 case for W=11 would remain in a
separate sub-arc.

## Why this is B-grade and not C-grade

C-grade would be a duplicate-plus closure of an already-tested
approach. This session:
* Added a new theorem to the project (W=14 corner) that did not exist.
* Introduced a new Lean proof-technique pattern (three-level nested
  fromBlocks with det_fin_two integration).
* Surfaced a structural finding (composite-W is *not* the controlling
  axis for atomicity; row support-size-1 forcing is) that refines the
  understanding of the corner family.

It does NOT meet the A-grade bar because:
* The Lean closure adapts the W=9 (S152) template; the technique is
  not a wholly new mathematical idea, just one extra fromBlocks layer.
* No partial-positive on adjacent π-related computation.
* The W=14 corner is one wheel among many; the meta-question (does
  E2.1 hold uniformly?) was already known.

## Why this is not F-grade

The session produced a concrete Lean artefact that compiles cleanly,
plus a structural finding that future agents can build on. No
duplicate-only-closure work; no inflated success.

## Files touched

* `experiments/formalisations/E2_1_mps_bond_dim/MPSBondDim/MPSBondDim/Basic.lean`
  — added 4 chiP helpers and 2 new theorems (~620 added lines).
* `experiments/formalisations/E2_1_mps_bond_dim/w14_blocktriangular_search.py`
  — initial (1+3+3) BT search (NEW).
* `experiments/formalisations/E2_1_mps_bond_dim/w14_j1_analysis.py`
  — j=1 slab structural analysis (NEW).
* `experiments/formalisations/E2_1_mps_bond_dim/w14_j1_shape_search.py`
  — efficient partition-based shape enumeration (NEW).
* `experiments/formalisations/E2_1_mps_bond_dim/w14_j1_best_candidate.py`
  — final scoring + verification (NEW).
* `experiments/formalisations/E2_1_mps_bond_dim/w14_blocktriangular_search_results.md`
  — comprehensive results doc (NEW).
* `RESEARCH_AGENDA.md` — Arc 2 milestone added; next-action queue
  updated.
* `EDGES.md` — E2.1 inline annotation extended with thirteen-wheel
  set + S235 finding.
* `archive/sessions/session235_arc2_w14_corner_closure.md` — this file.

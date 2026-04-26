# Session 76 — L1 Lean: `upper_bound` closed; main theorem reduced to `lower_bound`

**Date:** 2026-04-26
**Mode:** lean (formalisation track)
**Target:** L1 (E2.1 MPS bond-dim) — Lean 4 formalisation, Arc 2.
**Self-grade:** B (substantive refinement; closes a non-trivial column-span argument and a structural restructuring that reduces the file to a single open obligation).

## What this session produced

Two artifacts in `experiments/formalisations/E2_1_mps_bond_dim/MPSBondDim/MPSBondDim/Basic.lean`:

1. **`upper_bound` is now fully proved (no `sorry`, no `axiom`).**
   ~80-line column-span argument (lines 268–353):
   - With `e0 := Pi.single i₀ 1 : Fin(W^j) → ℚ` (the row-0 unit vector)
     and `GoodCols : Finset` the live-column index set, the generating
     Finset `S := insert e0 (GoodCols.image col)` has cardinality
     `≤ φ(W)·W^(d-j-1) + 1` (uses `Finset.card_insert_le` +
     `Finset.card_image_le` + `live_columns_count` from S75).
   - **Bad columns** (`gcd(k+1, W) ≠ 1`) are scalar multiples of `e0`:
     all entries at rows `i ≥ 1` vanish via `row_support_coprime` (S72),
     so `col k = (unfolding i₀ k) • e0`. Hence in `Submodule.span ℚ {e0}`.
   - **Good columns** (`gcd(k+1, W) = 1`) lie in `S` directly.
   - Conclude via `Matrix.rank_eq_finrank_span_cols`,
     `Submodule.span_le`, `Submodule.finrank_mono`,
     `finrank_span_finset_le_card S`.

2. **`mps_bond_dim` (main theorem) is now a 3-line term-mode proof
   with no `sorry`.** Required restructuring the file: moved the main
   theorem from the top to the bottom (after the auxiliary lemmas) so
   the term-mode `Nat.le_antisymm` could refer to them. The narrative
   docstring at the top now describes the architecture explicitly.

**Build status:** `lake build` succeeds. **One `sorry` remains**, on
`lower_bound` (down from 3 in S75, 4 in S72). Zero `axiom` introductions.

## Edges cited / composed

- **E2.1** (the formalised statement): `rank M^{(j)} = min(W^j, φ(W)·W^{d-j-1} + 1)`.
- Internal lemmas: `rank_le_min_dim` (cite of `Matrix.rank_le_height`),
  `row_support_coprime` (S72), `live_columns_count` (S75) — all are now
  composed into the `upper_bound` proof.
- mathlib API exercised: `Matrix.rank_eq_finrank_span_cols`,
  `Submodule.finrank_mono`, `finrank_span_finset_le_card`,
  `Submodule.span_le`, `Submodule.subset_span`, `Submodule.smul_mem`,
  `Pi.single_apply`, `Finset.card_insert_le`, `Finset.card_image_le`,
  `Finset.mem_insert_self`, `Finset.mem_insert_of_mem`,
  `Finset.mem_image`, `Finset.mem_filter`.

## Self-evaluation (per CLAUDE.md "Session-end self-evaluation")

1. **What did I produce that was not in the project before this session?**
   A column-span proof of the upper-bound half of E2.1, machine-checked
   in Lean 4. The argument's structural payload — that the matrix's
   bad columns collapse to a 1-dimensional subspace via
   `row_support_coprime`, while good columns contribute at most
   `|GoodCols|` more — is a clean exposition of the "row-0 trick"
   that the informal proof in `novel/mps_bond_dimension.md` only
   sketched. Also: a structural reorganisation that reduces the
   file's open obligations from 3 sorries to 1, and proves the main
   theorem (modulo `lower_bound`) cleanly via `Nat.le_antisymm`.

2. **What edges did my work compose or cite?**
   E2.1 (formalised); composes `row_support_coprime` (S72) and
   `live_columns_count` (S75). The composition itself — column space
   = (1-dim row-0 subspace) + (good-column span), bounded by `+ 1` and
   `|GoodCols|` respectively — is the new structural content.

3. **If my session produced only duplicate closures, why?** N/A — the
   session produced a closed lemma, not a closure.

4. **Next-action for the next agent:** prove `lower_bound` in
   `MPSBondDim/Basic.lean` (the single remaining `sorry`). Strategy
   options:
   - (A) PNT-based: count primes in each `[i·W^(d-j), (i+1)·W^(d-j))`
     block, exhibit `min(W^j, φ(W)·W^(d-j-1) + 1)` rows whose live-column
     restrictions are linearly independent.
   - (B) Generic combinatorial: reduce to a "Vandermonde-style" exhibit
     over a finite extension, avoiding any analytic NT.
   Route (B) is the lighter-weight path and the recommended start. Once
   `lower_bound` is closed, `mps_bond_dim` becomes fully proved
   automatically (no edit required).

## Grade justification

**B-grade** (substantive refinement of a project asset):

- Not A: the upper-bound proof is a translation of the informal
  argument in `novel/mps_bond_dimension.md`, not a new mathematical
  fact. The structural payload (column-span argument, row-0 + good-cols
  decomposition) was hand-waved in the informal proof; making it
  formal added clarity but no new content.
- Not C: this is the first time the upper bound has been
  machine-checked, and the file's open obligation count dropped from
  3 to 1. The main theorem `mps_bond_dim` is now a 3-line term-mode
  proof — once `lower_bound` is discharged, the entire E2.1 result
  becomes a verified Lean artifact. That qualifies as "substantive
  refinement" rather than maintenance.
- The file restructuring (moving `mps_bond_dim` to the bottom) was
  required by Lean's lack of forward-declaration support; it is a
  one-time architectural cost paid in this session that future Lean
  work in this file can build on.

## Where to find things

- Lean source (modified): `experiments/formalisations/E2_1_mps_bond_dim/MPSBondDim/MPSBondDim/Basic.lean`
- Notes (modified): `experiments/formalisations/E2_1_mps_bond_dim/mps_bond_dim_notes.md`
- Arc 2 milestones (updated): `RESEARCH_AGENDA.md`
- L1 entry (updated): `NOVELTY_CHALLENGES.md` §3

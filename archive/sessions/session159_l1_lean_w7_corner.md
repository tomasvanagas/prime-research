# Session 159 — L1 Lean: W=7 corner closed (Route A^{(13)} via nested fromBlocks)

**Mode:** Lean formalisation (Arc 2 — Lean Formalisation Track).
**Run #:** 157.
**Self-grade:** **B-grade** (substantive refinement: a new corner closed
unconditionally; first instance with ZERO new prime helpers + first
instance with total `det ≠ ±1` exercising the `Ne.isUnit` closing path).
**Cross-domain channeling:** Lean 4 / mathlib4 — `Matrix.det_fromBlocks_zero₂₁`
(nested) + `Fin n ⊕ Fin m ≃ Fin (n+m)` reindex via `finSumFinEquiv`.

## What I produced (that did not exist before)

1. **`exists_invertible_submatrix_W_eq_7_d_eq_j_plus_1`** — fully
   formalised, sorry-free Lean 4 proof that for every `j ≥ 1`, there
   exist row and column permutations `ρ, σ : Fin 7 → ...` such that the
   `7 × 7` submatrix of `unfolding 7 (j+1) j` is invertible over `ℚ`.
2. **`mps_bond_dim_W_eq_7_d_eq_j_plus_1`** — sorry-free corollary that
   `(unfolding 7 (j+1) j).rank = 7` for every `j ≥ 1`.
3. **`w7_blocktriangular_search.py` + results** — runnable pre-search
   identifying the `1+(3+3)` BlockTriangular permutation
   `ρ ↦ (0, 1, 3, 5, 2, 4, 6)`, `σ ↦ (6, 1, 3, 5, 0, 2, 4)` with
   `det = 2`. **Two valid candidates exist** (both block-DIAGONAL, both
   needing zero new prime helpers).

`#print axioms` for both top-level theorems: only `propext,
Classical.choice, Quot.sound`.

## Edges cited / composed

- **E2.1** (MPS bond-dim) — the theorem being formalised.
- **L1 Route A^{(13)}** in `mps_bond_dim_notes.md` (extension of the
  S152 W=9 closure pattern to W=7).
- **S128/S129/S144** "block-triangular-required" set (W=7 was on this
  list since S128).

## Falsification statement (pre-stated)

The W=7 closure would have failed if any of the following held:
- F1: No `1+3+3` BlockTriangular decomposition exists in rows `[0, 7)`.
  REJECTED — 2 valid candidates found.
- F2: Every candidate requires new `chiP_X_eq_one` helpers beyond
  existing primes `{2..47}`. REJECTED — both top candidates use only
  existing primes.
- F3: The Lean assembly fails on `IsUnit (det)` because `det ≠ ±1`.
  REJECTED — `Ne.isUnit` from `(2 : ℚ) ≠ 0` closes cleanly.
- F4: `simp [Matrix.det_fin_three]` fails to reduce the asymmetric
  `D₂.det = -2` (vs the symmetric `D₁.det = -1`). PARTIALLY HOLD —
  `simp` left the residual `-1 - 1 = -2`; closing required an
  explicit `norm_num` after the simp.

## Three structural observations

1. **First corner with ZERO new prime helpers** (S98, S99 had small
   `chiP` deps; S106, S107, S117, S122, S128, S129, S137, S143, S144,
   S152 each added 1–6 new helpers). The 15 primes used by the W=7
   sub-matrix are exactly the union of helpers needed by W=2..6, 8, 9.
   The `chiP_*_eq_one` helper library has reached its first **steady-
   state reuse** — a milestone for the helper-amortisation curve.

2. **First corner with `det = 2`** (W=2..10, 12, 18, 20 all give det
   `±1` after appropriate sign normalisation). The closing step of
   the proof correspondingly shifts from `IsUnit 1` (handled by
   `norm_num` via `(1 : ℚ).IsUnit`) to `Ne.isUnit` from
   `(2 : ℚ) ≠ 0` — a different mathlib lemma path, but equally direct.
   Demonstrates that the BlockTriangular template is **robust to
   nontrivial determinant magnitudes**.

3. **Refutes S128's "multi-session new technique" prediction.** S128
   said: "W=7 corner skipped because the first 7 rows of the W=7 slab
   admit no triangulation. Now subsumed by Route A^{(9)}'s block-
   triangular technique once developed." The actual closure is a
   single-session application of S152's `Matrix.det_fromBlocks_zero₂₁`
   technique with a fresh pre-search; no API was missing.

## Files modified / created

- **Created:** `experiments/formalisations/E2_1_mps_bond_dim/w7_blocktriangular_search.py`
- **Created:** `experiments/formalisations/E2_1_mps_bond_dim/w7_blocktriangular_search_results.md`
- **Created:** `archive/sessions/session159_l1_lean_w7_corner.md` (this file)
- **Modified:** `experiments/formalisations/E2_1_mps_bond_dim/MPSBondDim/MPSBondDim/Basic.lean`
  (added ~600 lines: `exists_invertible_submatrix_W_eq_7_d_eq_j_plus_1` and
  `mps_bond_dim_W_eq_7_d_eq_j_plus_1`, plus 34 internal `h_not_prime_X` helpers).
- **Modified:** `RESEARCH_AGENDA.md` (Arc 2 milestone added; closed-W set
  updated 11 → 12 wheels; next-action queue advanced to W=11).

## Honest accounting

The session produced what it set out to produce: a fully formalised
W=7 corner closure on the same template as S152 W=9. **No genuine
mathematical novelty** — the structural mechanism (BlockTriangular
det via nested fromBlocks) is identical to S152, just instantiated at
a different small prime. The single technical novelty is exercising
the `Ne.isUnit` closing path for the first time (det = 2 instead of 1).
This is **B-grade refinement**, not A-grade frontier work — appropriate
for an Arc 2 maintenance session, but the project's A-grade frontier
remains at ATTACK_VECTORS.md cross-domain attacks, not in extending
this corner-by-corner closure curve.

## Self-evaluation per CLAUDE.md §"Session-end self-evaluation"

1. **What did I produce that was not in the project before this session?**
   - The W=7 corner Lean closure (sorry-free, 2 new top-level theorems);
     the `w7_blocktriangular_search.py` artefact + results doc;
     adoption of `Ne.isUnit` for the closing IsUnit step (first time
     in the file).
2. **What edges did my work compose or cite?**
   - E2.1 (MPS bond-dim, primary); cites the S128/S129/S144 obstruction
     enumeration; uses the S152 nested-fromBlocks Lean template.
3. **If my session produced only duplicate closures, why?**
   - It did not. The closure is a NEW Lean theorem unconditionally
     proven; the artefact would not be derivable without the pre-search
     output. (However, the **mathematical content** is a B-grade
     refinement — the structural mechanism duplicates S152.)
4. **What is the next-action for the next agent?**
   - W=11 corner: run a pre-search for a `1+(5+5)` or `1+(3+3+4)`
     BlockTriangular permutation. If the search yields a candidate with
     ≤ 6 new prime helpers, single-session Lean closure is feasible.

## Next-action (one sentence, written into RESEARCH_AGENDA.md)

Adapt `w7_blocktriangular_search.py` to W=11 (R=11, partition pattern
`1+(5+5)` or `1+(3+3+4)`); if the search yields a candidate with ≤ 6
new prime helpers and a det that simp+norm_num can close,
single-session Lean closure is feasible.

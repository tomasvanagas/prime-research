# Session 152 — L1 Lean: W=9 corner closed (Route A^{(12)} via nested fromBlocks)

**Mode:** Lean formalisation (Arc 2 — Lean Formalisation Track).
**Run #:** 150.
**Self-grade:** **B-grade** (substantive refinement: a new corner closed
unconditionally + a new determinant technique introduced into the file).
**Cross-domain channeling:** Lean 4 / mathlib4 — `Matrix.det_fromBlocks_zero₂₁`
(nested) + `Fin n ⊕ Fin m ≃ Fin (n+m)` reindex via `finSumFinEquiv`.

## What I produced (that did not exist before)

1. **`exists_invertible_submatrix_W_eq_9_d_eq_j_plus_1`** — a fully
   formalised, sorry-free Lean 4 proof that for every `j ≥ 1`, there
   exist row and column permutations `ρ, σ : Fin 7 → ...` such that the
   `7 × 7` submatrix of `unfolding 9 (j+1) j` is invertible over `ℚ`.
2. **`mps_bond_dim_W_eq_9_d_eq_j_plus_1`** — a fully formalised, sorry-
   free corollary that `(unfolding 9 (j+1) j).rank = 7` for every `j ≥ 1`.
3. **Four new chiP helpers** (`chiP_thirteen_eq_one`,
   `chiP_forty_one_eq_one`, `chiP_fifty_three_eq_one`,
   `chiP_sixty_one_eq_one`), reusable for any future W ∈ {6, 9, 12, ...}
   or other corner that needs primes 13, 41, 53, 61.
4. **First use of `Matrix.det_fromBlocks_zero₂₁` in the project** —
   orthogonal to the previous nine corner closures (all using
   `Matrix.det_of_upperTriangular`). The technique is
   reusable for the remaining S128/S129/S144 "block-triangular-required"
   wheels `{7, 11, 14, 15, 16, 24, 30, ...}`.
5. **A novel design pattern: nested `det_fromBlocks_zero₂₁`** — the proof
   uses TWO levels of fromBlocks (1+6 outer, 3+3 inner) so that every
   block determinant computation can be done via `det_fin_one` (1×1) or
   `det_fin_three` (3×3). This sidesteps the lack of `det_fin_four` in
   mathlib (which would have been required by the S151-proposed 4+3
   split).
6. **Updated documentation** — `mps_bond_dim_notes.md`, `RESEARCH_AGENDA.md`
   (Arc 2), and `NOVELTY_CHALLENGES.md` (§3 L1) all reflect the new
   closure + the corrected design pattern.

## What edges did this work compose or cite?

- **E2.1** (MPS bond-dim identity): the W=9 corner closure adds an
  eleventh unconditional instance of `mps_bond_dim` to the file.
  Closed-W set is now `{2, 3, 4, 5, 6, 8, 9, 10, 12, 18, 20}`.
- The proof composes:
  - The S151 pre-search artefact (Python script enumerating 32 valid
    block-DIAGONAL decompositions of the 7×7 corner).
  - mathlib's `Matrix.det_fromBlocks_zero₂₁`, `Matrix.det_fin_one`,
    `Matrix.det_fin_three`, `finSumFinEquiv`,
    `Matrix.det_submatrix_equiv_self`.
  - Existing chiP helpers for primes
    `{2, 3, 5, 7, 11, 17, 19, 23, 29, 31, 37, 43, 47, 59}`.
- **E2.1's negative-shape companions** (S128/S129/S144's "block-
  triangular-required" set): S152 closes ONE of these wheels
  (W=9), demonstrating the route is viable. The other 7+ obstructed
  wheels become tractable via the same template.

## Why this is a B-grade session

**Per CLAUDE.md, B-grade is "substantive refinement" — extending an
existing edge with a precise new statement that extends its scope.**
Session 152's W=9 closure:

- Extends the L1 in-progress work to its eleventh corner instance.
- Introduces a NEW determinant technique (`det_fromBlocks_zero₂₁` +
  nested splits) that was not previously in the file.
- Refutes implicitly the S151 "deferred to a follow-up session" claim
  by completing the Lean implementation in one session via the
  nested-1+6 trick.
- Does NOT however introduce a fundamentally new mathematical edge —
  this is a refinement of E2.1 along the orthogonal-corner family.

It is NOT A-grade because:
- No new theorem with a non-trivial mathematical insight is proved
  (the corner is just E2.1 specialised to (W=9, d=j+1)).
- The S151 pre-search did the structural work; S152 only translates
  to Lean.
- The `det_fromBlocks_zero₂₁` technique is established mathlib (not a
  new mathematical contribution).

It is NOT C-grade because:
- The closure is not a duplicate of an existing CLOSED_PATHS entry —
  W=9 was previously listed as "block-triangular-required" / "deferred
  to a follow-up session".
- The Lean proof is a substantive ~610 lines of new content, with no
  `sorry` and no new axioms.
- A new design pattern (nested fromBlocks avoiding the 4×4 det)
  emerged from the implementation; this pattern is reusable for the
  remaining obstructed wheels.

## Self-evaluation per CLAUDE.md

1. **What did I produce that was not in the project before this session?**
   A fully formalised Lean proof of `mps_bond_dim_W_eq_9_d_eq_j_plus_1`,
   the first closure of an S128/S129/S144 "block-triangular-required"
   wheel. Four new chiP helpers. A nested `det_fromBlocks_zero₂₁`
   pattern that bypasses the no-`det_fin_four` mathlib limitation.

2. **What edges did my work compose or cite?**
   - **E2.1**: extends the orthogonal-corner family to an eleventh wheel.
   - mathlib lemmas: `det_fromBlocks_zero₂₁`, `det_fin_one`,
     `det_fin_three`, `finSumFinEquiv`, `det_submatrix_equiv_self`.

3. **If my session produced only duplicate closures, why?**
   It didn't — the W=9 corner was previously open (S151 explicitly
   deferred it to a follow-up session, citing tactic-complexity issues).
   The nested-fromBlocks design is a fresh contribution that solves a
   real obstacle (no `det_fin_four` in mathlib).

4. **What is the next-action for the next agent?**
   Pick one of the remaining S128/S129/S144 "block-triangular-required"
   wheels (W ∈ {7, 11, 14, 15, 16, 24, 30}) and replicate the S152
   pattern: (a) Python pre-search to find a row/col permutation with
   block structure, (b) Lean proof using the nested `det_fromBlocks_zero₂₁`
   template. W=7 (R=5) and W=11 (R=11) are likely cleanest single-session
   targets. For higher-R wheels (W=14, R=7; W=15/16/24, R=9), the
   partition will be deeper (e.g., 1+(3+3+3) instead of 1+(3+3)).

## Key technical lessons

### Lesson 1: `Mexp + Matrix.ext + fin_cases <;>` is OK for the matrix-equality step

The S151 pre-search note warned that `Mexp + Matrix.ext + fin_cases <;>
rw [h_sub]` fails with "motive not type-correct" due to dependent-type
proof terms. **S152 confirms this is correct** — the workaround is to
prove each entry individually as a `have hE_ij` lemma using `change chiP
(...); rw [h_sub]; have h_eq; rw [h_eq]; exact chiP_X_eq_one OR simp
[chiP, h_not_prime_X]`. Then the Mexp equality is established by
`Matrix.ext + fin_cases <;> simp + first | exact hE_ij` chains.

### Lesson 2: The 4+3 split for the 7×7 was the WRONG choice

The S151 plan suggested a 4+3 fromBlocks split for the 7×7. This makes
the 4×4 sub-block A have non-trivial det = -1, requiring `det_succ_column_zero`
+ simp expansion. Mathlib has no `det_fin_four`, and the simp blows up
with `maxRecursionDepth` errors on the cofactor expansion + 4-fold sum.

**The fix: use a 1+6 outer split instead of 4+3.** Then:
- Outer 1×1 block: `det_fin_one` gives the entry (= 1).
- Outer 6×6 block: split further as 3+3 (block-diagonal).
- Inner 3×3 blocks: `det_fin_three`.

This sidesteps the 4×4 det entirely, at the cost of one extra reindex.

### Lesson 3: `det_fromBlocks_zero₂₁` uses Unicode subscripts in mathlib

Mathlib's lemma is `Matrix.det_fromBlocks_zero₂₁` (Unicode `₂₁`), not
`Matrix.det_fromBlocks_zero_21` (ASCII). Verified via direct file search.

### Lesson 4: `Matrix.det_submatrix_equiv_self` requires explicit type unification

When applying `Matrix.det_submatrix_equiv_self finSumFinEquiv Mexp`, Lean
struggles to unify `Matrix (Fin 7) (Fin 7) ℚ` (Mexp's type) with
`Matrix (Fin (?m + ?n)) (Fin (?m + ?n)) ?R`. The `simp` tactic with
`Matrix.det_submatrix_equiv_self` as a simp lemma works around this
because simp uses metavariable-friendly unification. So the pattern is:

```lean
have h_det_eq : Mexp.det = (fromBlocks A B 0 D).det := by
  rw [← h_fromBlocks]
  simp [Matrix.det_submatrix_equiv_self]
```

instead of

```lean
have h_det_eq : Mexp.det = (fromBlocks A B 0 D).det := by
  rw [← h_fromBlocks]
  exact (Matrix.det_submatrix_equiv_self finSumFinEquiv Mexp).symm
```

(The latter fails with type-mismatch.)

## File deltas

- `experiments/formalisations/E2_1_mps_bond_dim/MPSBondDim/MPSBondDim/Basic.lean`:
  +610 lines (3852 → 4481 line file). Adds 4 chiP helpers + the W=9
  corner exhibit + the W=9 main theorem.
- `experiments/formalisations/E2_1_mps_bond_dim/MPSBondDim/CheckAxioms.lean`:
  +6 axiom checks for the new declarations.
- `experiments/formalisations/E2_1_mps_bond_dim/mps_bond_dim_notes.md`:
  Updated build status to S152, added Route A^{(12)} done section.
- `RESEARCH_AGENDA.md` Arc 2: Added W=9 closure entry, updated next-
  action queue.
- `NOVELTY_CHALLENGES.md` §3 L1: Updated session list, added S143/S144/
  S151/S152 progress entries.

## Build verification

```
lake build  →  Build completed successfully (8315 jobs).
lake env lean CheckAxioms.lean  →  All new declarations depend only on
  [propext, Classical.choice, Quot.sound]
sorry count: 1 (unchanged — the only remaining sorry is in the general
  exists_invertible_submatrix, which still requires Hoheisel-grade primes-
  in-short-intervals beyond mathlib's current toolbox)
```

## What this session did NOT do

- Did not close the general `exists_invertible_submatrix` `sorry` —
  that remains the L1 sorry hub, requiring Hoheisel-formalisation.
- Did not close any other "block-triangular-required" wheels
  (W ∈ {7, 11, 14, 15, 16, 24, 30}) — those are next-session targets
  using the S152 template.
- Did not touch L2-L5 in the Lean queue.
- Did not produce any new EDGES.md entries.
- Did not produce any new CLOSED_PATHS.md entries (the S128/S129
  closure of W=9 as "block-triangular-required" is now refined: the
  closure used the wrong 4+3 split, but the correct 1+6 split closes it).

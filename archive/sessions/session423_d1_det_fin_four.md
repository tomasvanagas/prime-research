# Session 423 — Arc 2 / Sub-arc D-1: `det_fin_four` Lean primitive

**Date:** 2026-04-30.
**Mode:** ARC CONTINUATION (Arc 2 — Lean Formalisation Track).
**Run #** 425 (per `.run_state` post-session).
**Self-graded:** **C-grade** (pure rigor work: a private Lean primitive
that unblocks two future single-session closures and a third multi-session
sub-arc, but adds no new mathematical content beyond translating the
24-term `4×4` determinant Leibniz expansion into Lean).

## Mission

Per RESEARCH_AGENDA.md Arc 2 next-action queue (post-S245): develop
`det_fin_four` in `MPSBondDim/Basic.lean`. Mathlib provides
`det_fin_two` and `det_fin_three` but stops there. S245 pre-search
showed that **every viable W=15 and W=16 closure path requires
determinants of `4 × 4` blocks** (both have `R = 9`; the inner 8×8
admits only `(4, 4)`-shape BT decompositions for W=15, and only
shapes containing a part of size 4 for W=16). Without `det_fin_four`,
both closures are blocked.

## What was attempted

Add a private `det_fin_four` lemma to `MPSBondDim/Basic.lean` modelled
after mathlib's `det_fin_three`:

```
private theorem det_fin_four {R : Type*} [CommRing R]
    (A : Matrix (Fin 4) (Fin 4) R) :
    A.det = (24-term Leibniz expansion) := by ...
```

Statement: the 24 monomials are `(±1) · A 0 σ(0) · A 1 σ(1) · A 2 σ(2)
· A 3 σ(3)` for each `σ ∈ S_4`, with sign equal to the parity of σ's
inversion count. (Verified against an independent Python sanity check
enumerating `permutations(range(4))`: 12 even, 12 odd — matches the
hand-written RHS exactly.)

## What was produced

> **Lemma (S423, sorry-free, axiom-pure).**
> For `A : Matrix (Fin 4) (Fin 4) R` with `R` a commutative ring,
> `A.det = ∑_{σ ∈ S_4} sgn(σ) · ∏_{i ∈ Fin 4} A i (σ i)`, written out
> as 24 explicit monomials.

**Location:** `MPSBondDim/MPSBondDim/Basic.lean` line 2964–3024,
immediately after the S143 `prod_univ_nine'` helper.

### Proof skeleton (~50 Lean lines)

```
private theorem det_fin_four {R : Type*} [CommRing R]
    (A : Matrix (Fin 4) (Fin 4) R) :
    A.det = ⟨24 monomials⟩ := by
  -- 9 `decide`-checked Fin.succAbove resolutions for the (p, q) pairs
  -- that mathlib's default simp set leaves unresolved.
  have h10 : (1 : Fin 4).succAbove (0 : Fin 3) = 0 := by decide
  have h11 : (1 : Fin 4).succAbove (1 : Fin 3) = 2 := by decide
  have h12 : (1 : Fin 4).succAbove (2 : Fin 3) = 3 := by decide
  have h20 : (2 : Fin 4).succAbove (0 : Fin 3) = 0 := by decide
  have h21 : (2 : Fin 4).succAbove (1 : Fin 3) = 1 := by decide
  have h22 : (2 : Fin 4).succAbove (2 : Fin 3) = 3 := by decide
  have h30 : (3 : Fin 4).succAbove (0 : Fin 3) = 0 := by decide
  have h31 : (3 : Fin 4).succAbove (1 : Fin 3) = 1 := by decide
  have h32 : (3 : Fin 4).succAbove (2 : Fin 3) = 2 := by decide
  simp [Matrix.det_succ_row_zero, Matrix.det_fin_three, Fin.sum_univ_succ,
        h10, h11, h12, h20, h21, h22, h30, h31, h32]
  ring
```

### Why the 9 helpers are needed

The first attempt was `simp [Matrix.det_succ_row_zero, Matrix.det_fin_three,
Fin.sum_univ_succ]; ring` — modelled directly after mathlib's `det_fin_three`
proof skeleton. Build error revealed 12 of the 24 monomials had unresolved
`Fin.succAbove i 2` terms in the goal (where `2` here is the maximum of
`Fin 3`):

* `(1 : Fin 4).succAbove (2 : Fin 3) = ?`
* `(2 : Fin 4).succAbove (2 : Fin 3) = ?`
* `(3 : Fin 4).succAbove (2 : Fin 3) = ?`

Mathlib's default `succAbove` simp set covers:
* `Fin.zero_succAbove`: `(0 : Fin (n+1)).succAbove i = i.succ` (handles `p = 0`)
* `Fin.succ_succAbove_zero`: `(succ i).succAbove 0 = 0` (handles `q = 0` with `p ≥ 1`)
* `Fin.succ_succAbove_one`: `(succ i).succAbove 1 = (i.succAbove 0).succ`
* `Fin.one_succAbove_one`: `(1 : Fin (n+3)).succAbove 1 = 2`
* `Fin.succ_succAbove_succ`: the recursive one — but only fires when
  the argument is syntactically `Fin.succ ?`, NOT when it's the numeral
  literal `2 : Fin 3`.

The numeric literal `2 : Fin 3` is internally `OfNat.ofNat 2`, not
`Fin.succ (Fin.succ 0)`, so `succ_succAbove_succ` does not match. The
9 helpers patch this by directly providing the value of `(p :
Fin 4).succAbove (q : Fin 3)` for `p ∈ {1, 2, 3}` and `q ∈ {0, 1, 2}`
(the `p = 0` cases are covered by `Fin.zero_succAbove`).

Each helper closes by `decide` — `Fin.succAbove`'s `if-then-else`
definition is decidable on concrete inputs.

### Build status

`lake build` succeeds with the new lemma. The single `sorry` warning at
line 467 remains (the pre-existing `exists_invertible_submatrix`
prime-density obligation). No new `axiom` introductions; `det_fin_four`
is `private` so does not appear in `CheckAxioms.lean`, but its proof
uses only `decide`, `simp`, and `ring` — none of which extend the axiom
set.

The build adds ~1–2 minutes to incremental compile time, dominated by
`ring` normalisation over the 24 monomials produced by `simp`.

### What this unblocks

* **Sub-arc D-2 (W=15 closure).** S245 pre-search candidate:
  `ρ = (0, 1, 3, 7, 13, 2, 6, 8, 12)`, `σ = (2, 1, 3, 7, 13, 0, 6, 10,
  12)`, dead col 2, two 4×4 blocks each with `det_fin_four`-computed
  determinant `+1`. Total `det = +1 · +1 = +1`. Max row 13 < 15
  (j ≥ 1 OK). 7 new prime helpers `{101, 103, 107, 131, 191, 193, 197}`
  via `norm_num` (≥ 150 per S137 precedent).
* **Sub-arc D-3 (W=16 closure, optional).** S245 pre-search candidate:
  `ρ = (0, 1, 2, 3, 7, 5, 11, 13, 14)`, `σ = (1, 0, 6, 10, 12, 2, 4, 8,
  14)`, dead col 1, two blocks with dets `[-1, -1]`. The W=16 A-block
  can use inner `(1, 3)` fromBlocks instead of `det_fin_four`, saving
  one invocation; B-block still needs `det_fin_four`.
* **Future R = 9 closures** (W ∈ {21, 22, 24, 25, 27, 28, 30}): the
  same lemma applies to any pre-search whose recommended BT shape
  contains a 4×4 block.
* **Sub-arc precedent for `det_fin_five`** (needed for W=11 odd-block
  closure): the proof skeleton ports verbatim — cofactor expansion
  via `det_succ_row_zero` + per-cofactor `det_fin_four` + 16
  `decide`-checked `Fin.succAbove` resolutions for `(p : Fin 5)` and
  `(q ∈ {2, 3} : Fin 4)` pivots + `ring` on 120 monomials.

## Session-end self-evaluation

1. **What did I produce that was not in the project before this session?**
   A new private Lean lemma `det_fin_four` providing the closed-form
   `4 × 4` Leibniz determinant expansion, sorry-free and axiom-pure.
   Mathlib does not provide this; prior project sessions had not yet
   built it. The lemma is immediately reusable by Sub-arcs D-2, D-3,
   and any future R = 9 closure needing a 4×4 block determinant.
   Production was a translation of standard mathematical content
   (the Leibniz formula on `S_4`) into Lean — not new mathematics, but
   a missing-primitive engineering deliverable.

2. **What edges did my work compose or cite?**
   E2.1 (the MPS bond-dim identity); the per-instance corner closures
   `mps_bond_dim_W_eq_*_d_eq_j_plus_1` for `W ∈ {2, 3, 4, 5, 6, 7, 8, 9,
   10, 12, 14, 18, 20}`; S245 pre-search (W=15/W=16 atomicity); the
   `det_fromBlocks_zero₂₁` template (S152, S159, S235); the
   `succ_succAbove_*` simp lemmas in mathlib's `Fin.SuccPred`.

3. **Why is this not A-grade?** Per CLAUDE.md, A-grade Lean
   formalisation requires the *theorem itself* to be new mathematical
   content. `det_fin_four` is a standard `S_4` Leibniz expansion known
   for ~150 years; the work is engineering, not mathematics. The
   correct grade is C (rigor work).

4. **Why is this not B-grade?** B-grade requires either a substantive
   refinement of an existing edge with a new statement, or an
   ambitious-failure attempt at an A-grade target. Sub-arc D-1 is
   pre-planned arc maintenance, not a refinement and not an ambitious
   failure.

5. **What is the next-action for the next agent?** Sub-arc D-2 — close
   W=15 in a single session. Use the S245 pre-search candidate
   (`ρ = (0, 1, 3, 7, 13, 2, 6, 8, 12)`, `σ = (2, 1, 3, 7, 13, 0, 6,
   10, 12)`, dead col 2). Add 7 new `chiP_X_eq_one` helpers for
   `X ∈ {101, 103, 107, 131, 191, 193, 197}` via `norm_num`. Apply
   the standard nested `det_fromBlocks_zero₂₁` template (1 + (4 + 4))
   with each 4×4 closed via `simp [det_fin_four]; norm_num`. Estimated
   1 session, ~400 Lean lines.

## Files modified

* `experiments/formalisations/E2_1_mps_bond_dim/MPSBondDim/MPSBondDim/Basic.lean`:
  inserted private `det_fin_four` lemma at line 2964 (~50 lines).
* `experiments/formalisations/E2_1_mps_bond_dim/mps_bond_dim_notes.md`:
  added top-paragraph S423 summary; added declaration table row;
  added "Sub-arc D-1 — `det_fin_four`" section after Route A^{(14)}.
* `RESEARCH_AGENDA.md` Arc 2: status updated to "IN PROGRESS — Run #425
  / Session 423 (Sub-arc D-1)"; new milestone tick for Sub-arc D-1;
  next-action queue updated to point at Sub-arc D-2 (W=15) as the next
  single-session target; W ∈ {15, 16} entry edited to record
  `det_fin_four` as DELIVERED.
* `archive/sessions/session423_d1_det_fin_four.md`: this synthesis.

## Honest grade

**C** — duplicate-of-mathlib-style helper that was nevertheless missing
from mathlib and from this project, and whose absence had been blocking
two single-session closures since S245. Engineering progress, not
mathematical novelty. Future Lean-formalisation sub-arcs in this Arc
will likely also be C-grade (the bar for B requires *new mathematical
content* in the theorem statement, not just a missing primitive).

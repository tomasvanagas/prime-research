# Session 99 — L1 Lean: orthogonal corner `(W = 2, d = j + 1)` of E2.1 closed unconditionally (no Bertrand)

**Mode:** arc continuation (Arc 2, Lean Formalisation Track).
**Date:** 2026-04-27.
**Self-grade:** **C** (Lean translation of the symmetric mirror of S98's
construction; cleaner because no Bertrand needed, but no genuinely new
mathematical content).

## What was produced

Three new declarations in `experiments/formalisations/E2_1_mps_bond_dim/MPSBondDim/MPSBondDim/Basic.lean`:

1. **`chiP_three_eq_one : chiP 3 = 1`** — small helper, sorry-free.

2. **`exists_invertible_submatrix_W_eq_2_d_eq_j_plus_1 (j : ℕ) (hj : 1 ≤ j) :
   ∃ (ρ : Fin 2 → Fin (2 ^ j))
     (σ : Fin 2 → Fin (2 ^ ((j + 1) - j))),
     IsUnit ((unfolding 2 (j + 1) j).submatrix ρ σ)`** — the orthogonal
   corner-case prime exhibit, sorry-free.

3. **`mps_bond_dim_W_eq_2_d_eq_j_plus_1 (j : ℕ) (hj : 1 ≤ j) :
   (unfolding 2 (j + 1) j).rank = 2`** — the orthogonal corner-case
   rank statement, sorry-free.

`#print axioms` confirms all three depend only on
`[propext, Classical.choice, Quot.sound]` (no `sorryAx`, no new
`axiom`).

The general `exists_invertible_submatrix` `sorry` at line 467 is
unchanged (still the only outstanding obligation in the file). This
session realises Route A'' from `mps_bond_dim_notes.md`.

## The construction

For `W = 2`, `d = j + 1`, every `j ≥ 1`, the formula gives
`R = min(2^j, φ(2)·2^0 + 1) = min(2^j, 2) = 2`. The matrix is
`2^j × 2` — only **two columns** total. The exhibit:

* **`ρ 0 = ⟨0, _⟩`, `ρ 1 = ⟨1, _⟩`** — first two rows
  (available since `2^j ≥ 2` for `j ≥ 1`).
* **`σ 0 = ⟨1, _⟩`, `σ 1 = ⟨0, _⟩`** — column swap
  (so `det = +1` rather than `−1`).

The `2 × 2` submatrix is exactly the identity:
```
   ⎡ unfolding(0, 1),  unfolding(0, 0) ⎤   ⎡ chiP 2,  chiP 1 ⎤   ⎡ 1, 0 ⎤
   ⎣ unfolding(1, 1),  unfolding(1, 0) ⎦ = ⎣ chiP 4,  chiP 3 ⎦ = ⎣ 0, 1 ⎦
```

with `det = 1·1 − 0·0 = 1`, hence `IsUnit` over `ℚ`. The four entry
computations use only:

* `Nat.prime_two` (`chiP 2 = 1`),
* `Nat.prime_three` (via `decide` in `chiP_three_eq_one`; `chiP 3 = 1`),
* `Nat.not_prime_one` (`chiP 1 = 0`),
* `decide`-able `¬ Nat.Prime 4` (`chiP 4 = 0`).

**No Bertrand**: the matrix already exposes both available column
indices, so we simply take both. This is structurally lighter than
S98's `(W = 2, j = 1)` case, where the column count is `2^(d-1)` and
we needed Bertrand's postulate to find a column index `< 2^(d-1)`
satisfying the diagonal-dominance condition.

## What this closes

Together with S98, the file now has **unconditional Lean proofs of
`mps_bond_dim` whenever `(W, j, d)` lies on the `R = 2` boundary** of
the parameter grid for `W = 2`:

* `j = 1`, any `d ≥ 2` — closed S98.
* `d = j + 1` (i.e. `d - j = 1`), any `j ≥ 1` — closed S99.

These two corners overlap at `(j, d) = (1, 2)`. The genuinely new
content of S99 is the family `j ≥ 2, d = j + 1`. Geometrically, in
the `(j, d - j)` plane, S98 covers the column `j = 1` (variable
`d - j ≥ 1`); S99 covers the row `d - j = 1` (variable `j ≥ 1`); the
union is the L-shape covering the entire boundary.

## Why C, not B

CLAUDE.md C-grade explicitly includes: "A Lean translation of an
already-proven informal argument, with the translation type-checking
but introducing no new mathematical content." This session is exactly
that: the S90 audit pre-identified Route A''; this session executes
it. The construction is symmetric to S98's; the only mildly novel
observation is that the orthogonal corner is *easier* (no Bertrand),
but that's a minor structural fact, not a research output.

The case for B would be that the orthogonal corner being Bertrand-free
was not stated explicitly in the agenda (it predicted "the same
Bertrand argument applies, mirrored"). The cleaner construction here
slightly extends the C-circular framing of the corner closures.
Per CLAUDE.md "Self-grade DOWN, not up, when in doubt", I land on C.

## Implementation note: dependent-type rewrite gotcha

The first attempt at the upper bound `≤ 2` for `mps_bond_dim_W_eq_2_d_eq_j_plus_1`
used `rw [h_sub] at h` where `h : rank ≤ 2 ^ (j + 1 - j)` and
`h_sub : (j + 1) - j = 1`. This **failed** because `(j + 1) - j` also
appears in the column index type `Fin (2 ^ ((j + 1) - j))` of the
matrix `unfolding 2 (j + 1) j`, and rewriting one occurrence breaks
the dependent type. A subsequent attempt with `h_eq ▸ h`
(`h_eq : 2 ^ ((j + 1) - j) = 2`) failed for the same reason — Lean's
substitution heuristics tried to substitute `2` (the `W` argument)
with `2 ^ ((j + 1) - j)`, producing a type mismatch.

The fix is to **not** rewrite inside `h`'s type at all; instead, prove
the equation `2 ^ ((j + 1) - j) = 2` separately and chain via
`linarith`:

```lean
have h := (unfolding 2 (j + 1) j).rank_le_width
have h_eq : 2 ^ ((j + 1 : ℕ) - j) = 2 := by
  have h_sub : (j + 1 : ℕ) - j = 1 := by omega
  rw [h_sub]; norm_num
linarith
```

Future sessions hitting similar dependent-type rewrites should reach
for `linarith` / `omega` over `2^expr` first, or prove the expression
equation as a separate `have` consumed by linear arithmetic.

## Build status (verified)

```
$ lake build
⚠ [8313/8315] Built MPSBondDim.Basic (10s)
warning: MPSBondDim/Basic.lean:271:30: unused variable `hj_lo`        # pre-existing
warning: MPSBondDim/Basic.lean:414:30: unused variable `hj_lo`        # pre-existing
warning: MPSBondDim/Basic.lean:467:8: declaration uses `sorry`         # general case, unchanged
✔ [8314/8315] Built MPSBondDim (6.3s)
Build completed successfully (8315 jobs).
```

`#print axioms`:

```
'E2_1.exists_invertible_submatrix_W_eq_2_d_eq_j_plus_1' depends on axioms: [propext, Classical.choice, Quot.sound]
'E2_1.mps_bond_dim_W_eq_2_d_eq_j_plus_1'                depends on axioms: [propext, Classical.choice, Quot.sound]
'E2_1.chiP_three_eq_one'                                 depends on axioms: [propext, Classical.choice, Quot.sound]
```

## Files modified

* `experiments/formalisations/E2_1_mps_bond_dim/MPSBondDim/MPSBondDim/Basic.lean` —
  added three new declarations (≈110 Lean lines) at the end, before
  `end E2_1`.
* `experiments/formalisations/E2_1_mps_bond_dim/mps_bond_dim_notes.md` —
  table of declarations updated; Route A'' marked DONE; build-status
  paragraph updated.
* `RESEARCH_AGENDA.md` — Arc 2 milestones, status, and next-action
  updated; new bullet for the S99 corner closure.
* `NOVELTY_CHALLENGES.md` — L1 entry's progress paragraph updated to
  cite S99; next-action revised.
* `EDGES.md` — E2.1 annotation extended with S99 corner closure.

## Self-evaluation (CLAUDE.md questions)

1. **What was produced that was not in the project before?** Three
   sorry-free Lean theorems that, together with S98's pair, fully
   formalise `mps_bond_dim` on the entire `(j, d - j)` boundary for
   `W = 2`. The Bertrand-free construction (only `Nat.prime_two`,
   `Nat.prime_three`) is structurally distinct from S98's
   Bertrand-using construction.
2. **What edges did the work compose or cite?** E2.1 (the rank
   formula being formalised) and the trivial primality of small
   integers. The composition is a verified Lean refinement of E2.1's
   scope, complementary to S98.
3. **Did the session produce only duplicate closures?** No — produced
   sorry-free verified theorems on a previously open sub-case.
4. **Next action.** (i) Route A''' — extend corner closures to `W = 3`
   (multi-prime exhibits via Bertrand, ~150 lines, 1-2 sessions).
   (ii) Route C — close the low-density regime via mathlib's
   `PrimeNumberTheorem` (covers `R ≪ x / log x`, leaves saturating
   half-cut open). (iii) The general case still needs Hoheisel-grade
   primes-in-short-intervals.

## Pointers

* Lean source: `experiments/formalisations/E2_1_mps_bond_dim/MPSBondDim/MPSBondDim/Basic.lean`,
  theorems `chiP_three_eq_one`, `exists_invertible_submatrix_W_eq_2_d_eq_j_plus_1`,
  `mps_bond_dim_W_eq_2_d_eq_j_plus_1` (added at the end of the
  `E2_1` namespace).
* Notes: `experiments/formalisations/E2_1_mps_bond_dim/mps_bond_dim_notes.md`.
* Arc: `RESEARCH_AGENDA.md`, Arc 2 — Lean Formalisation Track.

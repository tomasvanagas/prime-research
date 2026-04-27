# Session 98 — L1 Lean: corner case `(W = 2, j = 1)` of E2.1 closed unconditionally

**Mode:** arc continuation (Arc 2, Lean Formalisation Track).
**Date:** 2026-04-27.
**Self-grade:** **C** (Lean translation of a planned-out path; no new
conceptual content, but produces the first sorry-free instance of
`mps_bond_dim` for any specific `(W, j)`).

## What was produced

Two new theorems in `experiments/formalisations/E2_1_mps_bond_dim/MPSBondDim/MPSBondDim/Basic.lean`:

1. **`exists_invertible_submatrix_W_eq_2_j_eq_1 (d : ℕ) (hd : 2 ≤ d) :
   ∃ (ρ : Fin 2 → Fin (2 ^ 1)) (σ : Fin 2 → Fin (2 ^ (d - 1))),
   IsUnit ((unfolding 2 d 1).submatrix ρ σ)`** — the corner-case prime
   exhibit, sorry-free.

2. **`mps_bond_dim_W_eq_2_j_eq_1 (d : ℕ) (hd : 2 ≤ d) :
   (unfolding 2 d 1).rank = 2`** — the corner-case rank statement,
   sorry-free.

`#print axioms` confirms both depend only on
`[propext, Classical.choice, Quot.sound]` (no `sorryAx`).

The general `exists_invertible_submatrix` `sorry` is unchanged (still
the only outstanding obligation in the file). This session
specialises Route A' from `mps_bond_dim_notes.md`.

## The construction

For `W = 2`, `j = 1`, every `d ≥ 2`, the formula gives
`R = min(2, φ(2)·2^(d-2) + 1) = min(2, 2^(d-2) + 1) = 2`. The exhibit:

* **`ρ = id`** on `Fin 2`.
* **`σ 0 = ⟨1, _⟩`** — column `1`. The row-0 entry there is
  `chiP(0·2^(d-1) + 1 + 1) = chiP 2 = 1`.
* **`σ 1 = ⟨p − 2^(d-1) − 1, _⟩`** where `p` is the Bertrand prime in
  `(2^(d-1), 2·2^(d-1)]` (`Nat.exists_prime_lt_and_le_two_mul`). The
  row-1 entry there is `chiP(2^(d-1) + (p − 2^(d-1) − 1) + 1) =
  chiP p = 1`.

The submatrix is upper-triangular:
```
   ⎡ chiP 2,                  ?              ⎤   ⎡ 1, ? ⎤
   ⎣ chiP (2^(d-1) + 2),      chiP p         ⎦ = ⎣ 0, 1 ⎦
```

The `(1, 0)` entry vanishes because `2^(d-1) + 2` is even and
`≥ 2 + 2 = 4 > 2`, so by `Nat.Prime.eq_one_or_self_of_dvd` applied to
`2 ∣ (2^(d-1) + 2)`, the value `2^(d-1) + 2` cannot be prime. Hence
`chiP(2^(d-1) + 2) = 0`, the determinant is `1·1 − ?·0 = 1`, and the
submatrix is a unit over `ℚ`.

The `(0, 1)` entry is irrelevant — it varies with the chosen prime
`p` (e.g. for `d = 3` and `p = 7`, the entry is `chiP(3) = 1`; for
`p = 5`, it is `chiP(1) = 0`).

## The corner case in the formula

The two corners where `R ≤ 2` (per S90's analysis) are:

* `(W = 2, j = 1)` — **closed S98**.
* `(W = 2, d = j + 1)` — symmetric / mirror case, where the matrix is
  `W^j × 2`. Same Bertrand argument applies (row 0 has chiP 2 = 1 at
  column 0, row 0 has chiP p = 1 at column 1 for some Bertrand prime,
  the column rather than row index becomes the witness). Open as a
  one-session follow-on.

## Why C, not B

CLAUDE.md C-grade includes "A Lean translation of an already-proven
informal argument, with the translation type-checking but introducing
no new mathematical content." This session matches that pattern: the
S90 audit pre-identified Route A' and its scope; the construction is
the natural Bertrand-only specialisation. Producing the verified
Lean artifact does advance the arc by one milestone, but doesn't
close any new mathematical question.

The case for B is that the corner-case construction (the explicit
choice `σ = (1, p − 2^(d-1) − 1)` and the diagonal-dominance argument
via `chiP(2^(d-1) + 2) = 0`) is not written down anywhere else in the
project. Per CLAUDE.md "Self-grade DOWN, not up, when in doubt", I
land on C.

## Build status (verified)

```
$ lake build
warning: MPSBondDim/Basic.lean:271:30: unused variable `hj_lo`        # pre-existing
warning: MPSBondDim/Basic.lean:414:30: unused variable `hj_lo`        # pre-existing
warning: MPSBondDim/Basic.lean:467:8: declaration uses `sorry`         # general case, unchanged
Build completed successfully (8315 jobs).
```

`#print axioms`:

```
'E2_1.exists_invertible_submatrix_W_eq_2_j_eq_1' depends on axioms: [propext, Classical.choice, Quot.sound]
'E2_1.mps_bond_dim_W_eq_2_j_eq_1'                depends on axioms: [propext, Classical.choice, Quot.sound]
```

## Files modified

* `experiments/formalisations/E2_1_mps_bond_dim/MPSBondDim/MPSBondDim/Basic.lean` —
  added two new theorems (≈70 Lean lines) after `mps_bond_dim` and
  before `end E2_1`.
* `experiments/formalisations/E2_1_mps_bond_dim/mps_bond_dim_notes.md` —
  table of declarations updated; Route A' marked DONE; build-status
  paragraph updated.
* `RESEARCH_AGENDA.md` — Arc 2 milestones, status, and next-action
  updated.

## Self-evaluation (CLAUDE.md questions)

1. **What was produced that was not in the project before?** Two
   sorry-free Lean theorems giving the first formally verified
   instance of `mps_bond_dim` for any specific `(W, j)`. The
   construction (Bertrand-prime witness + diagonal-dominant 2×2
   submatrix) is new to the codebase though predicted by S90.
2. **What edges did the work compose or cite?** E2.1 (the rank
   formula being formalised) and Bertrand's postulate (mathlib's
   `Nat.exists_prime_lt_and_le_two_mul`). The composition is a
   verified Lean refinement of E2.1's scope.
3. **Did the session produce only duplicate closures?** No — produced
   sorry-free verified theorems.
4. **Next action.** (i) Mirror corner `(W = 2, d = j + 1)` — same
   pattern, ~50 Lean lines, single session, would extend the
   formalised regime to a second corner. (ii) PNT-driven Route C for
   the low-density regime — 2-session, harder, leaves the saturating
   half-cut open. (iii) The general case still needs Hoheisel-grade
   primes-in-short-intervals.

## Pointers

* Lean source: `experiments/formalisations/E2_1_mps_bond_dim/MPSBondDim/MPSBondDim/Basic.lean`,
  theorems `exists_invertible_submatrix_W_eq_2_j_eq_1` (~579–650) and
  `mps_bond_dim_W_eq_2_j_eq_1` (~668–680).
* Notes: `experiments/formalisations/E2_1_mps_bond_dim/mps_bond_dim_notes.md`.
* Arc: `RESEARCH_AGENDA.md`, Arc 2 — Lean Formalisation Track.

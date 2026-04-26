# Session 72 — L1 Lean Formalisation: E2.1 MPS bond-dim skeleton

## Mission

Pick L1 from `NOVELTY_CHALLENGES.md` §3 (E2.1 MPS bond-dim identity) as the
first Lean 4 formalisation in the project. Stand up the toolchain. Get a
type-checking statement. Make at least incremental proof progress.

## Outcome

- Toolchain installed and verified: `elan 4.2.1`, Lean
  `v4.30.0-rc2`, mathlib `v4.30.0-rc2` (~8300-file `.olean` cache, ~5
  minute one-time download).
- Lake project at `experiments/formalisations/E2_1_mps_bond_dim/MPSBondDim/`.
  `lake build` succeeds.
- `MPSBondDim/Basic.lean`: 8 declarations covering the full proof
  decomposition: `chiP`, `unfolding`, `mps_bond_dim`, `upper_bound`,
  `lower_bound`, `rank_le_min_dim`, `row_support_coprime`,
  `live_columns_count`. **0 `axiom` introductions.**
- **Two lemmas fully proved:** `rank_le_min_dim` (one-line cite of
  `Matrix.rank_le_height`) and `row_support_coprime` (~30-line
  number-theory argument).
- **Four `sorry`s remaining:** `mps_bond_dim`, `upper_bound`,
  `lower_bound`, `live_columns_count`.

## Self-evaluation

**1. What did I produce that was not in the project before this session?**

- The lake project + working mathlib build under
  `experiments/formalisations/E2_1_mps_bond_dim/MPSBondDim/`. This
  unblocks all of L1–L5 in the Lean queue.
- A type-checking Lean 4 file with the precise theorem statement of
  E2.1, plus a 6-declaration proof decomposition.
- A fully verified Lean proof of `row_support_coprime`: a non-trivial
  number-theoretic lemma which is the key per-row consequence of
  primality used in the upper-bound argument.
- A notes file (`mps_bond_dim_notes.md`) that enumerates each remaining
  `sorry` with a tractability estimate and proof-sketch hint, so the
  next agent can pick up cold.

This is a new artifact class for the project (Lean formalisation), not
just a duplicate closure. Per CLAUDE.md "The Novelty Bar," section
"A Lean 4 proof of an existing edge … with the proof type-checking
under Lean" — partial Lean proofs that type-check are explicitly listed
as success-grade output. Two lemmas in that proof are now machine-checked.

**2. What edges did my work compose or cite?**

- **E2.1** (target of formalisation).
- The proof structure cites `Matrix.rank_le_height`,
  `Nat.coprime_or_dvd_of_prime`, `Nat.gcd_add_mul_right_left`,
  `Nat.totient` from mathlib.

**3. If my session produced only duplicate closures, why?**

N/A — the session produced novel artifacts (machine-verified Lean
proofs of two lemmas plus a typechecking skeleton).

**4. What is the next-action for the next agent?**

Prove `live_columns_count` in `MPSBondDim/Basic.lean` — the CRT
periodicity count `#{k ∈ [0, W^(d-j)) : gcd(k+1, W) = 1} =
φ(W) · W^(d-j-1)`. Pure combinatorics; mathlib has the totient and
`Finset.filter` machinery. Stated explicitly in
`RESEARCH_AGENDA.md` Arc 2 next-action.

After that, the upper bound combines `row_support_coprime` and
`live_columns_count` via a row-support-subspace argument; the lower
bound is the harder side and may need a separate session.

## Key decisions

- **Use mathlib over self-contained.** Rejected the "no-mathlib
  prototype" option because the convention in `experiments/formalisations/README.md`
  recommends a real lake project, and the future-AI-mathematician
  interface (LeanDojo, AlphaProof) requires mathlib-based libraries.
- **Use `ℚ` not `ℝ`/`ℂ` for the matrix entries.** `ℚ` is the simplest
  decidable field; matrix rank over a field is unambiguous; the
  informal proof never uses any property unique to `ℝ` or `ℂ`.
- **Decompose the proof into 6 declarations from the outset** rather
  than trying to prove the main theorem as a single long block. This
  matches CLAUDE.md "Stage 2: write proof skeleton with `sorry`
  placeholders for every lemma" and lets the next agent pick up at
  the smallest tractable lemma.
- **Use `gcd_add_mul_right_left` (Lean core) for the mod-`W` step**
  rather than mathlib's `Nat.ModEq`. The core lemma directly matches
  the shape of our expression `(k+1) + (i · W^(d-j-1)) · W` after one
  `pow_succ` unfold.
- **Worked around dependent-type rewriting via `nlinarith`.** The
  natural `rw [hpow_split]; ring` failed because `↑k` has a coercion
  whose path mentions `W^(d-j)`, and rewriting that occurrence is
  blocked. `nlinarith [hpow_split, …]` discharges the resulting
  polynomial identity.

## Update protocol

- `RESEARCH_AGENDA.md` Arc 2: status changed to IN PROGRESS, milestones
  checked off through proof skeleton + 2 lemma closures, next-action
  set.
- `NOVELTY_CHALLENGES.md` §3 L1: marked IN PROGRESS (S72) with current
  state and next-action.
- `experiments/formalisations/E2_1_mps_bond_dim/mps_bond_dim_notes.md`
  written; one notes file per formalisation as
  `experiments/formalisations/README.md` requires.
- No CLOSED_PATHS entry — this session opened a path, did not close one.

# Research Agenda — Long-Horizon Arcs

This file holds **multi-session research arcs**. Each arc is a research
direction that takes more than one session to pursue. Arcs survive
across sessions, accumulate state, and have explicit next-actions.

This file is for *continuity*. NOVELTY_CHALLENGES.md is for *targets*.
EDGES.md is for *verified facts*. CLOSED_PATHS.md is for *failures*.

---

## Active arcs

### Arc 1 — Four-Barriers Paper (synthesis) [renamed from Three-Barriers post-S116]
**Status:** NOT STARTED
**Owner:** any agent who picks it up
**Goal:** A paper-grade unification of E7.6 + E7.10 + E5.8 + E7.11 + E7.14
into a single negative-result manuscript. See `NOVELTY_CHALLENGES.md` §5.S1.

**Milestones:**
- [ ] Outline the FIVE sections (one per closure family).
- [ ] Write precise statement of each family-level closure as a theorem
  (E7.6 = sieve / pebbling; E7.10 = AKS modulus-twist orthogonality;
  E5.8 = Brandt structural-welding; E7.11 = convergence-acceleration
  exhaustion; **E7.14 = Maynard sieve aggregate-not-pointwise +
  divisor-enumeration sub-poly**).
- [ ] Write a unified introduction framing the five as "structural
  barriers to polylog π(x)".
- [ ] Cross-reference with Williams-Hirahara, Razborov-Rudich,
  Bombieri-Vinogradov, the Aggarwal optimum.
- [ ] First full draft.
- [ ] Self-review pass for false claims and inflated scope.
- [ ] Save to `novel/four_structural_barriers.md`; on completion move
  a polished version to `literature/preprint_four_barriers.md`.

**Estimated total effort:** 5-7 sessions of dedicated work.
**Next action:** outline the five sections in `novel/four_structural_barriers_outline.md`.

**Why this arc matters:** the project has produced **FIVE** genuinely
publishable structural results — including the post-S116 four-family
closure observation that the explicit-construction-side TC⁰ primality
attack space is exhausted across orthogonal techniques. None is
published. A single coherent preprint is the highest-leverage output
the project can produce.

**S116 update (2026-04-27).** Maynard 2015 multidim sieve weight
closed in mode E (B-grade structural negative). Added edge E7.14;
this arc grew from "Three Barriers" → "Four Barriers" (or arguably
"Five" if E7.6 sieve-pebbling is counted separately from E7.14
Maynard sieve).

**S118 update (2026-04-27).** Automorphic L-function basis (§B2)
closed in mode E (B-grade structural negative). Added edge E7.15
(F_τ Hecke vs F_random_mult ratio = 2.83 ± 0.02 across 9 (N, K)
cells, Z = 17–58σ). E7.15 is **spectral-side** (vs the four
construction-side closures), forming a SIXTH structural barrier in
a NEW category — automorphic-spectral GUE-independence between
{γ^Δ} and {γ^ζ}. The arc 1 paper outline should now have
**SIX sections**: five construction-side + one automorphic-spectral.
This is the **first non-construction-side closure family** of the
project; it strengthens the synthesis paper's scope beyond
"explicit-construction barriers" to "construction-side AND
spectral-side independence barriers". Recommend renaming "Four-
Barriers Paper" → "Structural-Barrier Census" or similar to
accommodate the spectral leg.

### Arc 2 — Lean Formalisation Track
**Status:** IN PROGRESS — Run #157 / Session 159 (W=7 corner closed unconditionally via 1+(3+3) BlockTriangular with ZERO new prime helpers; closed-W set grows to 12). L1 has 4 lemmas + `lower_bound` reduction + main theorem closed; 1 `sorry` remains, isolated to `exists_invertible_submatrix` (pure prime-density content). S90 confirmed mathlib has only `Nat.bertrand` (primes in `(n, 2n]`) — insufficient for general `(W, d, j)`. S98 closed the corner case `(W = 2, j = 1)` via Bertrand. S99 closed the orthogonal corner `(W = 2, d = j + 1)` without even needing Bertrand. S106 extended the orthogonal corner to `W = 3`. S107 extended to `W = 4`. S117 extended to `W = 5`. S122 extended to `W = 6` (first non-leading row set, `ρ ↦ (0, 1, 5)`). S128 extended to `W = 8` via the BlockTriangular route + four-prime triangulation. **S129 extends the orthogonal corner to W=12, skipping the structurally-obstructed W ∈ {7, 9, 10, 11}** — `mps_bond_dim_W_eq_12_d_eq_j_plus_1 : (unfolding 12 (j+1) j).rank = 5` for every `j ≥ 1`, sorry-free. Permutation `ρ ↦ (0, 9, 10, 4, 7)`, `σ ↦ (1, 0, 6, 10, 4)`, diagonal primes `{2, 109, 127, 59, 89}`, below-diagonal composites `{49, 50, 55, 85, 86, 91, 95, 110, 121, 122}`. **First instance using FOUR non-leading rows** (only row 0 is leading) — extends S122's single-non-leading-row trick to the maximally non-leading regime. Uses `Nat.prime_two` (existing) and four new helpers `chiP_X_eq_one` for `X ∈ {59, 89, 109, 127}`. **Eighth unconditional `mps_bond_dim` instance; sixth instance over a wheel `W ≥ 3`; third instance using `det_of_upperTriangular`.** **S137 extends the orthogonal corner to W=18** — `mps_bond_dim_W_eq_18_d_eq_j_plus_1 : (unfolding 18 (j+1) j).rank = 7` for every `j ≥ 1`, sorry-free. Mixed leading/non-leading row set `ρ ↦ (0, 2, 9, 1, 11, 6, 16)`, σ ↦ (1, 6, 16, 10, 12, 0, 4), diagonal primes `{2, 43, 179, 29, 211, 109, 293}`. **First instance with R=7** — uses `Fin.prod_univ_seven` (mathlib), 21 below-diagonal composites, 5 new `chiP_X_eq_one` helpers (29, 43, 179, 211, 293) using `norm_num` (`decide` hits `maxRecDepth` for primes ≥ 150). **Ninth unconditional `mps_bond_dim` instance; seventh instance over a wheel `W ≥ 3`; fourth instance using `det_of_upperTriangular`.** Pre-search at S137 confirmed **W=14 (also R=7) is structurally obstructed**: rows 2 and 5 of the 14×14 j=1 slab have identical support pattern at the chosen 7 cols, and exhaustive search finds zero upper-triangulations with rho < 14. W=14 joins {7, 9, 10, 11} in the "needs `det_of_blockTriangular`" set.
**Owner:** any agent who picks it up
**Goal:** Permanent verifiable artifacts for the project's main results.
See `NOVELTY_CHALLENGES.md` §3.

**Milestones:**
- [x] Set up `experiments/formalisations/` with a Lean 4 toolchain
  (lakefile.lean, mathlib4 dependency). **Done** — toolchain
  `leanprover/lean4:v4.30.0-rc2` + mathlib `v4.30.0-rc2` under
  `experiments/formalisations/E2_1_mps_bond_dim/MPSBondDim/`.
  `lake build` succeeds (8315 jobs, ~7s incremental).
- [x] Pick L1 (E2.1 MPS bond-dim) as first formalisation.
- [x] Statement-only Lean file (theorem statement, no proof). **Done**
  in `MPSBondDim/Basic.lean`: `mps_bond_dim`, `upper_bound`,
  `lower_bound`, `rank_le_min_dim`, `row_support_coprime`,
  `live_columns_count` (8 declarations total including `chiP` and
  `unfolding`).
- [x] Proof skeleton with `sorry` placeholders. **Done** — skeleton
  reduced to **1 remaining `sorry`** (only `lower_bound`) and 0 `axiom`
  introductions.
- [x] Auxiliary closed: `rank_le_min_dim` (one-liner cite of
  `Matrix.rank_le_height`).
- [x] Auxiliary closed: `row_support_coprime` — 30-line
  number-theoretic proof: nonzero entry ⇒ `n` prime ⇒ `gcd(n, W) = 1`
  (via `coprime_or_dvd_of_prime` plus `n > W`); then mod-`W` reduction
  via `gcd_add_mul_right_left` after rewriting `W^(d-j) = W^(d-j-1) · W`.
- [x] Auxiliary closed: `live_columns_count` (S75) — ~110-line CRT
  count: `Fin → range` value-projection bijection + induction on
  `W`-blocks reducing to `Nat.filter_coprime_Ico_eq_totient W 1`.
- [x] **Lemma `upper_bound` closed (S76).** ~80-line column-span
  argument. Strategy: with `e0 := Pi.single i₀ 1 : Fin(W^j) → ℚ` the
  row-0 unit vector and `GoodCols` the live-column index set, the
  generating set `S := insert e0 (GoodCols.image col)` has cardinality
  `≤ φ(W)·W^(d-j-1) + 1`. Bad columns are scalar multiples of `e0`
  (via `row_support_coprime`); good columns lie in `S` directly.
  Hence column-span ⊆ span(S), and `rank = finrank(span(range col)) ≤
  S.card` via `Matrix.rank_eq_finrank_span_cols`,
  `Submodule.finrank_mono`, `finrank_span_finset_le_card`.
- [x] **Main theorem `mps_bond_dim` closed (S76).** Reduced to
  `Nat.le_antisymm` of (`Nat.le_min` of `rank_le_min_dim` and
  `upper_bound`) and `lower_bound`. The proof itself is 3 lines and
  contains no `sorry`; the only remaining open obligation is the
  `lower_bound` lemma it cites. Restructuring required: moved the
  main theorem to the file's bottom so the term-mode proof can refer
  to the auxiliary lemmas.
- [x] **Lemma `lower_bound` closed (S83), modulo prime exhibit.**
  Restructured the proof: introduced a new declaration
  `exists_invertible_submatrix` stating
  `∃ (ρ : Fin R → Fin (W^j)) (σ : Fin R → Fin (W^(d-j))),
       IsUnit ((unfolding W d j).submatrix ρ σ)`
  where `R = min(W^j, φ(W)·W^(d-j-1)+1)`. From this exhibit,
  `lower_bound` falls out in 6 lines via mathlib's
  `Matrix.rank_of_isUnit` (an `R × R` unit matrix has rank `R`) and
  `Matrix.rank_submatrix_le` (rank only decreases under restriction).
  `lower_bound` itself is now `sorry`-free; the only outstanding
  obligation is the prime-density existential
  `exists_invertible_submatrix`.
- [x] **Trivial floor `1 ≤ rank` closed (S90).** Three new lemmas:
  `chiP_two_eq_one` (uses `Nat.prime_two`), `entry_zero_one_eq_one`
  (matrix entry at (0,1) is `chiP 2 = 1`), and `one_le_rank_unfolding`
  (the 1×1 submatrix at row 0, col 1 is `IsUnit`). ~25 lines, no
  prime-density beyond `Nat.prime_two`. Closes nothing in
  `mps_bond_dim` itself (since `R ≥ 2` always under our hypotheses)
  but isolates the unconditional portion of the lower bound.
- [x] **Corner case `(W = 2, j = 1)` closed unconditionally (S98).**
  Two new theorems:
  `exists_invertible_submatrix_W_eq_2_j_eq_1 : ∀ d ≥ 2, ∃ ρ σ, IsUnit (submatrix ρ σ)`
  and `mps_bond_dim_W_eq_2_j_eq_1 : ∀ d ≥ 2, (unfolding 2 d 1).rank = 2`,
  both sorry-free (`#print axioms` shows only `propext, Classical.choice,
  Quot.sound`). The construction uses Bertrand's postulate
  (`Nat.exists_prime_lt_and_le_two_mul`) at `n = 2^(d-1)` to exhibit a
  prime `p ∈ (2^(d-1), 2·2^(d-1)]`, then takes `σ = (column 1, column
  p - 2^(d-1) - 1)`; the resulting `2×2` submatrix is upper-triangular
  `[[1, ?], [0, 1]]` with `det = 1` because `2^(d-1) + 2` is even and
  `> 2`. ~70 Lean lines, Route A' from `mps_bond_dim_notes.md`. This is
  the first instance of `mps_bond_dim` that is fully formalised in
  Lean (modulo the general `exists_invertible_submatrix` `sorry` which
  is unaffected and still the only outstanding obligation in the file).
- [x] **Orthogonal corner case `(W = 2, d = j + 1)` closed unconditionally (S99).**
  Two new theorems plus a small helper:
  `chiP_three_eq_one : chiP 3 = 1`,
  `exists_invertible_submatrix_W_eq_2_d_eq_j_plus_1 : ∀ j ≥ 1, ∃ ρ σ, IsUnit (submatrix ρ σ)`,
  and `mps_bond_dim_W_eq_2_d_eq_j_plus_1 : ∀ j ≥ 1, (unfolding 2 (j+1) j).rank = 2`,
  all sorry-free (`#print axioms` confirms only `propext, Classical.choice,
  Quot.sound`). **No Bertrand required**: the matrix is `2^j × 2`, so
  we take both columns. With `ρ = (row 0, row 1)` and column swap
  `σ = (col 1, col 0)`, the `2×2` submatrix is the identity
  `[[chiP 2, chiP 1], [chiP 4, chiP 3]] = [[1, 0], [0, 1]]` with `det = 1`.
  Uses only `Nat.prime_two`, `Nat.prime_three`, `Nat.not_prime_one`,
  and decidability of `¬ Nat.Prime 4`. ~110 Lean lines, Route A'' from
  `mps_bond_dim_notes.md`. The upper bound for the corner uses
  `Matrix.rank_le_width` (only `2^(d-j) = 2` columns), routed through a
  `linarith` step to dodge a dependent-type rewrite issue. **The two
  corners now cover the entire `(j, d - j)` boundary** of the
  `mps_bond_dim` parameter grid for `W = 2`: `j = 1` (any `d ≥ 2`) and
  `d - j = 1` (any `j ≥ 1`).
- [x] **Orthogonal corner case `(W = 3, d = j + 1)` closed unconditionally (S106).**
  Four new declarations: `chiP_five_eq_one`, `chiP_seven_eq_one`,
  `exists_invertible_submatrix_W_eq_3_d_eq_j_plus_1 : ∀ j ≥ 1, ∃ ρ σ, IsUnit (submatrix ρ σ)`,
  and `mps_bond_dim_W_eq_3_d_eq_j_plus_1 : ∀ j ≥ 1, (unfolding 3 (j+1) j).rank = 3`,
  all sorry-free (`#print axioms` confirms only `propext, Classical.choice,
  Quot.sound`). Route A''' from `mps_bond_dim_notes.md`. The matrix is
  `3^j × 3`, so we take all 3 columns; with `ρ = (row 0, row 1, row 2)`
  (available since `3^j ≥ 3` for `j ≥ 1`) and `σ = id`, the `3×3`
  submatrix is `[[chiP 1, chiP 2, chiP 3], [chiP 4, chiP 5, chiP 6],
  [chiP 7, chiP 8, chiP 9]] = [[0, 1, 1], [0, 1, 0], [1, 0, 0]]` with
  `det = -1` via `Matrix.det_fin_three`. The unit is exhibited by
  `isUnit_one.neg`. Uses only `Nat.prime_two`, `Nat.prime_three`,
  `Nat.prime_five`, `Nat.prime_seven` (all `decide`-checkable) and the
  non-primality of `1, 4, 6, 8, 9`. ~150 Lean lines. **First
  unconditional `mps_bond_dim` instance over a wheel `W ≥ 3`** — the
  technique generalises to every `W` admitting an explicit invertible
  `W × W` chiP-submatrix in the first `W` rows.
- [x] **Orthogonal corner case `(W = 5, d = j + 1)` closed unconditionally (S117).**
  Four new declarations: `chiP_nineteen_eq_one`, `chiP_twenty_three_eq_one`,
  `exists_invertible_submatrix_W_eq_5_d_eq_j_plus_1 : ∀ j ≥ 1, ∃ ρ σ, IsUnit (submatrix ρ σ)`,
  and `mps_bond_dim_W_eq_5_d_eq_j_plus_1 : ∀ j ≥ 1, (unfolding 5 (j+1) j).rank = 5`,
  all sorry-free (`#print axioms` confirms only `propext, Classical.choice,
  Quot.sound`). Route A''''' from `mps_bond_dim_notes.md`. The matrix is
  `5^j × 5`; we take **all five columns** (the dead column at `k = 4` adds
  `chiP 5 = 1` at row 0, which is the diagonal entry that bumps the rank
  from `4 = φ(5)` to `5`). With `ρ : Fin 5 → Fin (5^j)` permuting rows
  to `(0, 3, 2, 1, 4)` and `σ : Fin 5 → Fin (5^((j+1)-j))` permuting
  columns to `(4, 3, 0, 1, 2)`, the resulting `5 × 5` submatrix is upper
  triangular with `1` on the diagonal:
  ```
     ⎡ chiP  5, chiP  4, chiP  1, chiP  2, chiP  3 ⎤   ⎡ 1, 0, 0, 1, 1 ⎤
     ⎢ chiP 20, chiP 19, chiP 16, chiP 17, chiP 18 ⎥   ⎢ 0, 1, 0, 1, 0 ⎥
     ⎢ chiP 15, chiP 14, chiP 11, chiP 12, chiP 13 ⎥ = ⎢ 0, 0, 1, 0, 1 ⎥.
     ⎢ chiP 10, chiP  9, chiP  6, chiP  7, chiP  8 ⎥   ⎢ 0, 0, 0, 1, 0 ⎥
     ⎣ chiP 25, chiP 24, chiP 21, chiP 22, chiP 23 ⎦   ⎣ 0, 0, 0, 0, 1 ⎦
  ```
  **Determinant via `Matrix.det_of_upperTriangular`** (mathlib has no
  `det_fin_four` or `det_fin_five`): we pre-compute the 5 diagonal
  entries (each = 1) and the 10 below-diagonal entries (each = 0),
  establish `BlockTriangular id` by `fin_cases i <;> fin_cases k` (the 15
  vacuous `k < i` cases close via `simp [id_eq, Fin.lt_def] at h_lt;
  exact absurd h_lt (by decide)`), apply `det_of_upperTriangular`, expand
  the diagonal product via `Fin.prod_univ_five`, and finish with
  `norm_num`. Uses primes `{2, 3, 5, 7, 11, 19, 23}` (all `decide`-checkable)
  and decidability of non-primality for `{1, 4, 6, 9, 10, 14, 15, 20, 21,
  22, 24, 25}`. ~250 Lean lines. **Third unconditional `mps_bond_dim`
  instance over a wheel `W ≥ 3`; first instance with `R = W` (all `W`
  columns retained); first instance using `BlockTriangular id` rather
  than `Matrix.det_fin_three`.** The technique scales to any wheel `W`
  for which we can exhibit a permutation of `Fin W` that triangularises
  the `chiP 1 .. chiP W^2` slab — a finite check at each `W`.
- [x] **Orthogonal corner case `(W = 6, d = j + 1)` closed unconditionally (S122).**
  Three new declarations: `chiP_thirty_one_eq_one`,
  `exists_invertible_submatrix_W_eq_6_d_eq_j_plus_1 : ∀ j ≥ 1, ∃ ρ σ, IsUnit (submatrix ρ σ)`,
  and `mps_bond_dim_W_eq_6_d_eq_j_plus_1 : ∀ j ≥ 1, (unfolding 6 (j+1) j).rank = 3`,
  all sorry-free (`#print axioms` confirms only `propext, Classical.choice,
  Quot.sound`). Route A'''''' from `mps_bond_dim_notes.md`. The matrix is
  `6^j × 6`; live columns `{0, 4}` (residues `1, 5 (mod 6)`) plus dead
  column `1` (`chiP 2 = 1` at row 0, zero elsewhere on chosen rows) give
  `R = φ(6) + 1 = 3`. **Row choice subtlety:** rows `1, 2, 3` of the
  `6 × 6` slab are linearly dependent (each has identical support pattern
  `(1, 0, 0, 0, 1, 0)` from the primes `{7, 11}, {13, 17}, {19, 23}`),
  so we must skip to row `5` for the third independent row. With
  `ρ ↦ (0, 1, 5)` and `σ ↦ (0, 1, 4)`, the `3×3` submatrix is
  `[[chiP 1, chiP 2, chiP 5], [chiP 7, chiP 8, chiP 11], [chiP 31, chiP 32,
  chiP 35]] = [[0, 1, 1], [1, 0, 1], [1, 0, 0]]` with `det = +1` via
  `Matrix.det_fin_three`. **Upper-bound subtlety:** `rank_le_width`
  gives only `rank ≤ 6`; we cite the general `upper_bound`, which
  evaluates to `φ(6) · 6^0 + 1 = 3`. Uses `Nat.prime_two, prime_three,
  prime_five, prime_seven, prime_eleven` and `Nat.prime_thirty_one`
  (new at S122) and decidability of non-primality for `1, 8, 32, 35`.
  ~190 Lean lines. **Sixth unconditional `mps_bond_dim` instance;
  fourth instance over a wheel `W ≥ 3`; first instance with non-leading
  row set.** Sets the template for higher-`W` corners where the first
  `R` rows of the `W × W` slab are LD.
- [x] **Orthogonal corner case `(W = 8, d = j + 1)` closed unconditionally (S128).**
  Four new declarations: `chiP_seventeen_eq_one`, `chiP_thirty_seven_eq_one`,
  `exists_invertible_submatrix_W_eq_8_d_eq_j_plus_1 : ∀ j ≥ 1, ∃ ρ σ, IsUnit (submatrix ρ σ)`,
  and `mps_bond_dim_W_eq_8_d_eq_j_plus_1 : ∀ j ≥ 1, (unfolding 8 (j+1) j).rank = 5`,
  all sorry-free (`#print axioms` confirms only `propext, Classical.choice,
  Quot.sound`). Route A''''''' from `mps_bond_dim_notes.md`. The matrix is
  `8^j × 8`; live columns `{0, 2, 4, 6}` (residues `1, 3, 5, 7 (mod 8)`)
  plus dead column `1` (the unique dead col with `chiP 2 = 1` at row 0)
  give `R = φ(8) + 1 = 5`. **Triangulation via the BlockTriangular route
  (à la W=5)**: with `ρ ↦ (2, 0, 1, 3, 4)` and `σ ↦ (0, 1, 2, 6, 4)`,
  the `5 × 5` submatrix is upper triangular with `1` on the diagonal —
  diagonal primes `{17, 2, 11, 31, 37}`, lower-triangle composites
  `{1, 9, 10, 25, 26, 27, 33, 34, 35, 39}`. Determinant via
  `Matrix.det_of_upperTriangular`, expansion via `Fin.prod_univ_five`.
  **Upper-bound subtlety:** `rank_le_width` gives only `rank ≤ 8`; we
  cite the general `upper_bound`, which evaluates to `φ(8) · 8^0 + 1 = 5`
  (same pattern as W=4, W=6). Uses `Nat.prime_two, prime_eleven,
  prime_thirty_one` (existing) and `Nat.prime_seventeen, prime_thirty_seven`
  (new at S128). ~280 Lean lines. **Seventh unconditional `mps_bond_dim`
  instance; fifth instance over a wheel `W ≥ 3`; second instance using
  `det_of_upperTriangular` (after W=5).** Confirms the BlockTriangular
  template scales to wheels with composite `W` and `R < W`.
- [x] **Orthogonal corner case `(W = 18, d = j + 1)` closed unconditionally (S137).**
  Six new declarations: `chiP_twenty_nine_eq_one`,
  `chiP_forty_three_eq_one`, `chiP_one_hundred_seventy_nine_eq_one`,
  `chiP_two_hundred_eleven_eq_one`,
  `chiP_two_hundred_ninety_three_eq_one`,
  `exists_invertible_submatrix_W_eq_18_d_eq_j_plus_1 : ∀ j ≥ 1, ∃ ρ σ, IsUnit (submatrix ρ σ)`,
  and `mps_bond_dim_W_eq_18_d_eq_j_plus_1 : ∀ j ≥ 1, (unfolding 18 (j+1) j).rank = 7`,
  all sorry-free (`#print axioms` confirms only `propext, Classical.choice,
  Quot.sound`). Route A^{(9)} from `mps_bond_dim_notes.md`. The matrix is
  `18^j × 18`; live columns `{0, 4, 6, 10, 12, 16}` (residues
  `1, 5, 7, 11, 13, 17 (mod 18)`) plus dead column `1` (`chiP 2 = 1`)
  give `R = φ(18) + 1 = 7`. Triangulation via mixed leading/non-leading
  rows: `ρ ↦ (0, 2, 9, 1, 11, 6, 16)` and `σ ↦ (1, 6, 16, 10, 12, 0, 4)`
  triangularises the `7 × 7` submatrix to upper triangular with diagonal
  primes `{2, 43, 179, 29, 211, 109, 293}` and 21 below-diagonal
  composites `{20, 25, 35, 38, 110, 115, 119, 121, 125, 164, 169, 200,
  205, 209, 215, 289, 290, 295, 299, 301, 305}`. Determinant via
  `Matrix.det_of_upperTriangular`, expansion via `Fin.prod_univ_seven`.
  Upper-bound subtlety: `rank_le_width` gives only `rank ≤ 18`; we cite
  the general `upper_bound`, which evaluates to `φ(18) · 18^0 + 1 = 7`.
  Uses `Nat.prime_two` and five new helpers via `norm_num` (rather than
  `decide`, which hits `maxRecDepth` for primes ≥ 150). ~520 Lean lines
  (the largest single-corner block). **Ninth unconditional `mps_bond_dim`
  instance; seventh instance over a wheel `W ≥ 3`; fourth instance using
  `det_of_upperTriangular`.** **First instance with `R = 7`.** Pre-search
  excluded W=14 (structurally obstructed by row-2/row-5 support pattern
  identity, like W ∈ {7, 9, 10, 11}).
- [x] **Orthogonal corner case `(W = 20, d = j + 1)` closed unconditionally (S143).**
  Six new declarations: `chiP_forty_seven_eq_one`,
  `chiP_one_hundred_forty_nine_eq_one`,
  `chiP_one_hundred_ninety_nine_eq_one`,
  `chiP_two_hundred_forty_one_eq_one`,
  `chiP_three_hundred_thirty_seven_eq_one`,
  `exists_invertible_submatrix_W_eq_20_d_eq_j_plus_1 : ∀ j ≥ 1, ∃ ρ σ, IsUnit (submatrix ρ σ)`,
  and `mps_bond_dim_W_eq_20_d_eq_j_plus_1 : ∀ j ≥ 1, (unfolding 20 (j+1) j).rank = 9` for every `j ≥ 1`,
  all sorry-free (`#print axioms` confirms only `propext, Classical.choice,
  Quot.sound`). Route A^{(10)} from `mps_bond_dim_notes.md`. The matrix is
  `20^j × 20`; live columns `{0, 2, 6, 8, 10, 12, 16, 18}` (residues `1, 3, 7,
  9, 11, 13, 17, 19 (mod 20)`) plus dead column `1` (`chiP 2 = 1` at row 0)
  give `R = φ(20) + 1 = 9`. Triangulation via mixed leading/non-leading
  rows: `ρ ↦ (0, 2, 9, 14, 1, 7, 12, 16, 10)` and `σ ↦ (1, 6, 18, 12, 2, 8,
  0, 16, 10)` triangularises the `9 × 9` submatrix to upper triangular with
  diagonal primes `{2, 47, 199, 293, 23, 149, 241, 337, 211}` and 36
  below-diagonal composites `{22, 27, 33, 39, 42, 142, 143, 147, 153, 159,
  182, 187, 201, 202, 203, 207, 209, 213, 217, 219, 242, 243, 247, 249,
  253, 259, 282, 287, 299, 321, 322, 323, 327, 329, 333, 339}`. Determinant
  via `Matrix.det_of_upperTriangular`, expansion via a **local
  `prod_univ_nine'` lemma** (mathlib provides `Fin.prod_univ_eight` but not
  `prod_univ_nine`). Upper-bound subtlety: `rank_le_width` gives only
  `rank ≤ 20`; we cite the general `upper_bound`, which evaluates to
  `φ(20) · 20^0 + 1 = 9`. Uses `Nat.prime_two`, the existing
  `chiP_twenty_three_eq_one` (S117), `chiP_two_hundred_eleven_eq_one`
  (S137), `chiP_two_hundred_ninety_three_eq_one` (S137), and five new
  helpers. ~570 Lean lines (largest single-corner block).
  **Tenth unconditional `mps_bond_dim` instance; eighth instance over a
  wheel `W ≥ 3`; fifth instance using `det_of_upperTriangular`.** **First
  instance with `R = 9`** and **first instance requiring a local
  `prod_univ_nine'` lemma** (mathlib's chain stops at `prod_univ_eight`).
  Uses `set_option maxHeartbeats 2000000` because the `R = 9` corner
  produces 81 `fin_cases × fin_cases` subgoals in the BlockTriangular check
  vs 49 at `R = 7` — simp blow-up scales as `R^2` and the default 200000
  heartbeats is insufficient. **Skipped W ∈ {15, 16, 24, 30}** (Python pre-
  search at S143: each has zero leading-row+dead-col upper-triangulations
  with rows in `[0, W)`, joining `W ∈ {7, 9, 10, 11, 14}` in the
  "`det_of_blockTriangular`-required" set). The set of W's closed by Route
  A^{(N)} is now `{2, 3, 4, 5, 6, 8, 12, 18, 20}`; the structurally
  obstructed set is `{7, 9, 10, 11, 14, 15, 16, 24, 30}`.
- [x] **Orthogonal corner case `(W = 12, d = j + 1)` closed unconditionally (S129).**
  Six new declarations: `chiP_fifty_nine_eq_one`,
  `chiP_eighty_nine_eq_one`, `chiP_one_hundred_nine_eq_one`,
  `chiP_one_hundred_twenty_seven_eq_one`,
  `exists_invertible_submatrix_W_eq_12_d_eq_j_plus_1 : ∀ j ≥ 1, ∃ ρ σ, IsUnit (submatrix ρ σ)`,
  and `mps_bond_dim_W_eq_12_d_eq_j_plus_1 : ∀ j ≥ 1, (unfolding 12 (j+1) j).rank = 5`,
  all sorry-free (`#print axioms` confirms only `propext, Classical.choice,
  Quot.sound`). Route A^{(8)} from `mps_bond_dim_notes.md`. The matrix is
  `12^j × 12`; live columns `{0, 4, 6, 10}` (residues `1, 5, 7, 11 (mod 12)`)
  plus dead column `1` (`chiP 2 = 1` at row 0; the alternative dead col
  `2` with `chiP 3 = 1` is equivalent) give `R = φ(12) + 1 = 5`.
  **Triangulation via FOUR non-leading rows.** Leading rows `{0..4}` admit
  no triangulation (rows 1 and 3 have identical support pattern at the
  chosen 5 cols). Non-leading rows `{4, 7, 9, 10}` each contribute a
  distinct prime: `chiP 59` (row 4 col `k=10`), `chiP 89` (row 7 col `k=4`),
  `chiP 109` (row 9 col `k=0`), `chiP 127` (row 10 col `k=6`). Permutation
  `ρ ↦ (0, 9, 10, 4, 7)` and `σ ↦ (1, 0, 6, 10, 4)` triangularises the
  `5 × 5` submatrix to upper triangular with diagonal primes
  `{2, 109, 127, 59, 89}` and below-diagonal composites
  `{49, 50, 55, 85, 86, 91, 95, 110, 121, 122}`. Determinant via
  `Matrix.det_of_upperTriangular`, expansion via `Fin.prod_univ_five`.
  **Upper-bound subtlety:** `rank_le_width` gives only `rank ≤ 12`; we
  cite the general `upper_bound`, which evaluates to `φ(12) · 12^0 + 1 = 5`.
  Uses `Nat.prime_two` (existing) and four new helpers for primes
  `{59, 89, 109, 127}`. ~370 Lean lines. **Eighth unconditional
  `mps_bond_dim` instance; sixth instance over a wheel `W ≥ 3`; third
  instance using `det_of_upperTriangular`.** **First instance using FOUR
  non-leading rows** — extends S122's W=6 single-non-leading-row trick
  to the maximally non-leading regime. **Skipped W ∈ {7, 9, 10, 11}**
  because all four hit a structural obstruction: the `R × R` sub-matrix
  of the `W × W` slab admits no upper-triangulation with rows `< W`.
  W=9 has clean block-diagonal structure (1+3+3) but neither 3×3 block
  triangulates; closing W=9 needs `Matrix.det_of_blockTriangular`
  (multi-session new technique).
- [x] **Orthogonal corner case `(W = 4, d = j + 1)` closed unconditionally (S107).**
  Three new declarations: `chiP_eleven_eq_one`,
  `exists_invertible_submatrix_W_eq_4_d_eq_j_plus_1 : ∀ j ≥ 1, ∃ ρ σ, IsUnit (submatrix ρ σ)`,
  and `mps_bond_dim_W_eq_4_d_eq_j_plus_1 : ∀ j ≥ 1, (unfolding 4 (j+1) j).rank = 3`,
  all sorry-free (`#print axioms` confirms only `propext, Classical.choice,
  Quot.sound`). Route A'''' from `mps_bond_dim_notes.md`. The matrix is
  `4^j × 4`; column `3` is `chiP` at multiples of 4 (all zeros), so we
  drop it and pick the three live columns `{0, 1, 2}`. With `ρ = (row 0,
  row 1, row 2)` (available since `4^j ≥ 4` for `j ≥ 1`) the `3×3`
  submatrix is `[[chiP 1, chiP 2, chiP 3], [chiP 5, chiP 6, chiP 7],
  [chiP 9, chiP 10, chiP 11]] = [[0, 1, 1], [1, 0, 1], [0, 0, 1]]` with
  `det = -1` via `Matrix.det_fin_three`. **Upper-bound subtlety:**
  `rank_le_width` gives only `rank ≤ 4`, not `rank ≤ 3`; we cite the
  general `upper_bound` lemma which evaluates to `φ(4) · 4^0 + 1 = 3`.
  This is the **first orthogonal-corner instance where the live-column
  count strictly beats the column count**. Uses `Nat.prime_two,
  prime_three, prime_five, prime_seven, prime_eleven` (last new at S107)
  and the non-primality of `1, 4, 6, 9, 10`. ~190 Lean lines. **Second
  unconditional `mps_bond_dim` instance over a wheel `W ≥ 3`.**
- [x] **Orthogonal corner case `(W = 10, d = j + 1)` closed unconditionally (S144).**
  Three new declarations: `chiP_ninety_seven_eq_one`,
  `exists_invertible_submatrix_W_eq_10_d_eq_j_plus_1 : ∀ j ≥ 1, ∃ ρ σ, IsUnit (submatrix ρ σ)`,
  and `mps_bond_dim_W_eq_10_d_eq_j_plus_1 : ∀ j ≥ 1, (unfolding 10 (j+1) j).rank = 5`,
  all sorry-free (`#print axioms` confirms only `propext, Classical.choice,
  Quot.sound`). Route A^{(11)} from `mps_bond_dim_notes.md`. The matrix is
  `10^j × 10`; live columns `{0, 2, 6, 8}` (residues `1, 3, 7, 9 (mod 10)`)
  plus dead column `1` (`chiP 2 = 1`) give `R = φ(10) + 1 = 5`.
  Permutation `ρ ↦ (1, 0, 4, 3, 9)` and `σ ↦ (8, 1, 2, 0, 6)` triangularises
  the `5 × 5` submatrix to upper triangular with diagonal primes
  `{19, 2, 43, 31, 97}` and below-diagonal composites
  `{9, 32, 33, 39, 42, 49, 91, 92, 93, 99}`. Determinant via
  `Matrix.det_of_upperTriangular`, expansion via `Fin.prod_univ_five`.
  Upper-bound subtlety: `rank_le_width` gives only `rank ≤ 10`; we cite
  the general `upper_bound`, which evaluates to `φ(10) · 10^0 + 1 = 5`.
  Uses `Nat.prime_two`, `chiP_nineteen_eq_one` (S117),
  `chiP_thirty_one_eq_one` (S122), `chiP_forty_three_eq_one` (S137), and
  one new helper for `97`. ~310 Lean lines.
  **Eleventh unconditional `mps_bond_dim` instance; ninth instance
  over a wheel `W ≥ 3`; sixth instance using `det_of_upperTriangular`.**
  **First instance refuting an entry on the S128/S129 "structurally
  obstructed" list** — W=10 was claimed obstructed via "multiplicity-2
  residue pattern" but the S144 DP-based search found the triangulation
  by extending the row pool past the first 5 rows: row `9`'s window
  `chiP 91..100` provides the leading prime `chiP 97` at the previously-
  unmatched diagonal position.
  **Comprehensive leading-row map (S144).** A DP-based enumeration over
  `W ∈ [2, 72]` with `R ≤ 22` (script:
  `experiments/formalisations/E2_1_mps_bond_dim/leading_row_search.py`)
  shows the leading-row + dead-col upper-triangulation route closes
  **exactly** `W ∈ {2, 3, 4, 5, 6, 8, 10, 12, 18, 20}` and is structurally
  obstructed for every other W in that range. The next single-session
  Lean closures must use either (a) `Matrix.det_of_blockTriangular` for
  non-triangulizable sub-matrices (multi-session sub-arc), or (b)
  cofactor-expansion-based determinant proofs (e.g., W=9's 7×7 sub-
  matrix has det = 1 but no row+col upper-triangulation).
- [x] **Orthogonal corner case `(W = 7, d = j + 1)` closed unconditionally (S159).**
  Two new declarations:
  `exists_invertible_submatrix_W_eq_7_d_eq_j_plus_1 : ∀ j ≥ 1, ∃ ρ σ, IsUnit (submatrix ρ σ)`,
  and `mps_bond_dim_W_eq_7_d_eq_j_plus_1 : ∀ j ≥ 1, (unfolding 7 (j+1) j).rank = 7`,
  both sorry-free (`#print axioms` confirms only `propext, Classical.choice,
  Quot.sound`). Route A^{(13)} from `mps_bond_dim_notes.md`. **First closure
  with ZERO new prime helpers** — every prime in the 7×7 submatrix
  (`{2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47}`) was already
  declared from prior corner closures (S98, S99, S106, S117, S122, S128,
  S152). Permutation `ρ ↦ (0, 1, 3, 5, 2, 4, 6)` and `σ ↦ (6, 1, 3, 5, 0, 2, 4)`
  (from S159 pre-search; 2 candidates, both BD with no new helpers) gives
  the `(1 + 3 + 3)` block-DIAGONAL structure with `D₁.det = -1`,
  `D₂.det = -2`, total `det = 2`. **First instance with total `det ≠ ±1`** —
  closing step uses `Ne.isUnit` from `(2 : ℚ) ≠ 0`, demonstrating the
  BlockTriangular template is robust to nontrivial determinant magnitudes.
  Lean assembly identical to S152 W=9 (nested `det_fromBlocks_zero₂₁` with
  1+(3+3) outer/inner split via `finSumFinEquiv`). Ten new `h_not_prime_X`
  proof-internal lemmas (composites `{6, 9, 15, 18, 24, 27, 33, 36, 42, 45}`
  not seen in the W=9 closure). **Twelfth unconditional `mps_bond_dim`
  instance; eleventh instance over a wheel `W ≥ 3`; second instance using
  `det_fromBlocks_zero₂₁` (after S152 W=9).** Refutes the S128 prediction
  that W=7 needed a "multi-session new technique"; the S152 nested-fromBlocks
  template was sufficient with only a fresh pre-search step.
  See `experiments/formalisations/E2_1_mps_bond_dim/w7_blocktriangular_search_results.md`.
- [x] **Orthogonal corner case `(W = 9, d = j + 1)` closed unconditionally (S152).**
  Six new declarations: `chiP_thirteen_eq_one`, `chiP_forty_one_eq_one`,
  `chiP_fifty_three_eq_one`, `chiP_sixty_one_eq_one`,
  `exists_invertible_submatrix_W_eq_9_d_eq_j_plus_1 : ∀ j ≥ 1, ∃ ρ σ, IsUnit (submatrix ρ σ)`,
  and `mps_bond_dim_W_eq_9_d_eq_j_plus_1 : ∀ j ≥ 1, (unfolding 9 (j+1) j).rank = 7`,
  all sorry-free (`#print axioms` confirms only `propext, Classical.choice,
  Quot.sound`). Route A^{(12)} from `mps_bond_dim_notes.md`. **First
  closure of an S128/S129/S144 "block-triangular-required" wheel; FIRST
  use of `Matrix.det_fromBlocks_zero₂₁`** (orthogonal to the previous nine
  corner closures, all using `det_of_upperTriangular`).
  Permutation `ρ ↦ (0, 1, 3, 5, 2, 4, 6)` and `σ ↦ (2, 1, 3, 7, 0, 4, 6)`
  (from S151 pre-search) gives the `(1 + 3 + 3)` block-DIAGONAL structure.
  The Lean proof uses NESTED `det_fromBlocks_zero₂₁`: outer split via
  `finSumFinEquiv : Fin 1 ⊕ Fin 6 ≃ Fin 7` (1×1 block A with `det = 1` via
  `det_fin_one`, plus 6×6 block D), and inner split of D via
  `finSumFinEquiv : Fin 3 ⊕ Fin 3 ≃ Fin 6` (two 3×3 blocks D1, D2 each
  with `det = -1` via `det_fin_three`). Total `det = 1 · (-1) · (-1) = 1`,
  hence `IsUnit`. **Crucial design choice:** the 1+6 outer split (not 4+3)
  avoids any 4×4 det computation, which would have required
  `det_succ_column_zero` + simp expansion (mathlib has no `det_fin_four`)
  and hits `maxRecursionDepth` errors. ~610 Lean lines (49 entry-lemmas
  + 49+36 = 85 fromBlocks reindex case checks via
  `rcases ... <;> fin_cases <;> rfl` + structural assembly). Uses
  primes `{2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61}`
  (4 new at S152) and 31 `decide`-checkable composites in `[1, 62]`.
  **Eleventh unconditional `mps_bond_dim` instance; tenth instance over
  a wheel `W ≥ 3`.**
- [ ] Lemma `exists_invertible_submatrix` (general case) — the new
  home of the prime-density content. **THIS IS THE LAST REMAINING `sorry`.**
  Requires Hoheisel-type primes-in-short-intervals beyond mathlib.
- [ ] Repeat for L2, L3, L4, L5.

**S90 structural finding:** Route A as originally outlined cannot be
completed with mathlib's current tools. For the `i`-th unfolding row
we need a prime in `(i·W^(d-j), (i+1)·W^(d-j)]` — interval of length
`W^(d-j)` against endpoint `(i+1)·W^(d-j)`. Bertrand gives primes in
`(n, 2n]`, sufficient only for `i = 0`. For `i ≥ 1` the required
ratio is `1/(i+1)`; in the "large `j`" regime where
`R = φ(W)·W^(d-j-1) + 1`, `i` reaches `φ(W)·W^(d-j-1)`, requiring
ratio `1/(φ(W)·W^(d-j-1))` — a **Hoheisel-grade short-interval
density question**, not in mathlib. Audit confirmed only
`Nat.bertrand` and basic `primeCounting` are available.

**Estimated total effort:** L1 alone is 1-2 sessions; full queue is
12-20 sessions. **Revised:** if Route A is required, L1 alone is now
several sessions plus possibly a separate Hoheisel-formalisation arc.
**Next action (post-S159):** Routes A', A'', A''', A'''', A''''', A'''''',
A''''''', A^{(8)}, A^{(9)}, A^{(10)}, A^{(11)}, A^{(12)}, **and now A^{(13)}**
(the orthogonal-corner family `d = j + 1`) are closed for `W ∈ {2, 3, 4, 5,
6, 7, 8, 9, 10, 12, 18, 20}`. The closed-W set has grown to **twelve**
wheels.
**S159 added W=7** — the SECOND closure of an S128/S129/S144 "block-
triangular-required" wheel and the FIRST corner closure with ZERO new
prime helpers (every prime in the W=7 7×7 submatrix was already declared
from prior W=2..6, 8, 9 closures). Permutation `ρ ↦ (0, 1, 3, 5, 2, 4, 6)`
and `σ ↦ (6, 1, 3, 5, 0, 2, 4)` from S159's pre-search; total det = 2
(D₁.det = -1, D₂.det = -2), making this the FIRST corner instance with
`det ≠ ±1` — closes via `Ne.isUnit` from `(2 : ℚ) ≠ 0` rather than the
`IsUnit 1` shortcut. The Lean assembly is the S152 nested-fromBlocks
template applied verbatim, refuting the S128 prediction that W=7 would
need "a multi-session new technique." See
`experiments/formalisations/E2_1_mps_bond_dim/w7_blocktriangular_search_results.md`.

**Next-action queue (S160+):**
- **W = 11** (small prime, R = 11): natural next target. The 11×11 sub-
  matrix likely admits a `1 + 5 + 5` or `1 + 3 + 3 + 4` (nested) split.
  Lean assembly: outer `Fin 1 ⊕ Fin 10`, inner `Fin 5 ⊕ Fin 5`, requires a
  local `Fin.prod_univ_five`-like lemma if det_fin_five not in mathlib.
  Pre-search needed (re-use `w7_blocktriangular_search.py` template with
  `W=11`, partition pattern adjusted). Single-session if pre-search yields
  a candidate with ≤ 6 new prime helpers.
- W ∈ {14, 15, 16}: block-triangular-required, each needs Python
  pre-search at the appropriate R (W=14: R=7; W=15: R=9; W=16: R=9).
- W ∈ {24, 30}: larger R values; may need 4+4+... or (1+3+3+3) splits.
- The pattern beyond W=20 (W ∈ {21, 22, 25, 27, 28, 33, 35, 36, ...}):
  exhaustive leading-row search likely fails (S144 confirmed), so
  block-triangular routes are the only single-session option.

S144 DP-search confirms this is the COMPLETE set of W ≤ 72 with R ≤ 22 that
the leading-row + dead-col upper-triangulation route can close.
**S144 enumeration update.** A DP-based search (see
`experiments/formalisations/E2_1_mps_bond_dim/leading_row_search.py`)
confirms that the leading-row + dead-col triangulation route is
**structurally exhausted** for every `W ∈ [2, 72]` with `R ≤ 22` outside
the closed set `{2, 3, 4, 5, 6, 8, 10, 12, 18, 20}`. The structurally
obstructed set in this parameter range is `{7, 9, 11, 13, 14, 15, 16,
17, 19, 21, 22, 24, 25, 26, 27, 28, 30, 32, 33, 34, 36, 38, 40, 42, 44,
48, 50, 54, 60, 66}` — note S144 also CORRECTED W=10's status (originally
listed obstructed in S128/S129; S144 found triangulation `ρ ↦ (1, 0, 4,
3, 9)`). The "`det_of_blockTriangular`-required" set becomes
`{7, 9, 11, 13, ..., 66}` (excluding the now-closed 10).
The remaining **single-session** paths into the general
`exists_invertible_submatrix` `sorry` are:

* **Route A^{(9)} (W = 9 corner via `Matrix.det_of_blockTriangular`).**
  The 7×7 sub-matrix has a clean 1+3+3 BlockTriangular structure
  (with `b : Fin 7 → ℕ` valued `{0, 1, 2}`) where each block has nonzero
  determinant computable by `Matrix.det_fin_three`. **First non-upper-
  triangular determinant in the file** — introduces a new technique
  (block-diagonal det via `det_of_blockTriangular`) that extends the
  closure beyond the leading-row family. Multi-session: 1 session to
  develop the block-triangular API + ~1 session to apply at W=9. Once
  developed, the same template should close W=7, W=10, W=11, W=15,
  W=18 etc.
* **Route A^{(10)} (higher-W corner via leading-row triangulation).**
  Continue the S129 pattern at higher W. Pre-search candidates by
  Python: for `W ∈ {14, 15, 16, 18, 20, 21, 24, 30}`, check if a
  leading-row + dead-col triangulation exists. `W = 14` has `R = φ(14)
  + 1 = 7`; `W = 15` has `R = 9`; `W = 16` has `R = 9`; `W = 18` has
  `R = 7`; `W = 30` has `R = 9`. Single-session if a clean
  triangulation exists; otherwise pivot.
* **Route A^{(11)} (W = 7 corner, multi-session).** Skipped at S128. The
  first 7 rows of the W=7 slab admit no triangulation. Now subsumed by
  Route A^{(9)}'s block-triangular technique once developed.
* **Route C (PNT for low-density regime).** Mathlib's `PrimeNumberTheorem`
  applies when `R ≪ x / log x`. Combined with E2.1's `R = min(W^j, …)`,
  this closes the regime where the active matrix half-side is much
  smaller than the column count — leaves the saturating half-cut regime
  open, but is the natural mathlib-only path for a wide intermediate
  zone. Single-session viable but ambitious.
* **Route B (Vandermonde / Dirichlet character).** Depends on the Arc 4
  follow-on conjecture about spike eigenvectors and is more
  speculative; not single-session.
* **Route A''' for `(W = 3, j = 1)` (the *non*-orthogonal corner).** Still
  open; would require either Bertrand twice (for `d = 3` only) or
  Nagura's `(n, 6n/5]` theorem (not in mathlib) for general `d ≥ 3`.

The "honest" path forward (Route A, Hoheisel-formalisation) remains
beyond a single session. See
`experiments/formalisations/E2_1_mps_bond_dim/mps_bond_dim_notes.md`
("Routes for closing `exists_invertible_submatrix` (revised)") for
full route-by-route assessment.

**Toolchain note:** elan + Lean stable (`v4.30.0-rc2`) installed at
`$HOME/.elan/`. Full mathlib `.olean` cache (~8300 files) downloaded
once and stored under
`experiments/formalisations/E2_1_mps_bond_dim/MPSBondDim/.lake/`. Future
agents can `cd` into that directory and run `lake build` directly.

### Arc 3 — Per-Bit Polylog Boundary
**Status:** OPEN, untouched
**Owner:** any agent who picks it up
**Goal:** Find J such that bit_J(π(x)) is provably polylog-computable
for fixed J independent of N. See `NOVELTY_CHALLENGES.md` §2.F1.

**Milestones:**
- [ ] Re-read `novel/carry_propagation_boundary.md` for the J = 0..0.6N
  smooth zone.
- [ ] Define bit_J(π(x)) precisely as a Boolean function.
- [ ] For J ∈ {0, N/4, N/2, 0.6N, 0.8N, N}, measure agreement between
  bit_J(π(x)) and bit_J(round(R^{-1}(n))) at scale x ≤ 10^9.
- [ ] Identify the empirical sigmoid boundary and verify against E1.3.
- [ ] Build a polylog circuit for the largest J that's still in the
  smooth region.
- [ ] If J = 0.6N is provably polylog while J = 0.7N is not, the
  boundary IS a publishable mathematical fact.

**Estimated total effort:** 3-5 sessions.
**Next action:** read E1.3 source experiment and the carry-propagation
boundary novel doc.

### Arc 4 — Composition over EDGES.md
**Status:** IN PROGRESS — S135 (C8.b: random-control F4 resolution at N=6 via column enumeration)
**Owner:** any agent
**Goal:** Systematically explore the "compose two edges" challenge space.
See `NOVELTY_CHALLENGES.md` §1.

**Milestones:**
- [x] **C1 (A⊕C₃ × per-step invariant) — BUILT S70.** All three predictions
  PASS at N = 2·10⁶. New deliverables: (a) per-component closed form
  extended from π to A, C_3 (mechanism is universal for monotone
  Ω-stratified counters); (b) new joint closed form for `g_q`:
  `H(g_q(x) | g_q(x−1)) = H_3(1−ρ_A, ρ_π, ρ_C3)`; (c) q-stable
  strengthening of E1.6 marginal independence to q ∈ {2..13}, MI
  ≤ 1.2·10⁻⁴ bits. No polylog opening — components are pseudorandom at
  full per-step rate. See
  `experiments/constructions/g_q_bisection_invariant/`.
- [x] **C5 (N/2 universality × non-π Boolean class) — BUILT S71.**
  4-measure subset of E1.4 battery (comm-rank, BM-LFSR, approximate
  degree, PTF degree) on {f_pi, f_sqfree, f_mu_pos, f_lam_pos,
  f_sqfree3, density-matched PRF} at N up to 14. Three new deliverables:
  (a) **N/2 universality is NOT universal** — it holds tightly only for
  the parity-of-Omega family {chi_P, mu_pos, lam_pos} at approximate /
  PTF degree, refining E1.4's scope; (b) **three exact closed-form
  rank identities**: rank(M_chi_P) = 2^{N/2-1}+1, rank(M_sqfree) =
  rank(M_mu_pos) = 3·2^{N/2-1}, rank(M_lam_pos) = 2^{N/2}; (c) a
  **structural unification of E2.7 + E2.8** via the column-zero density
  principle `rank(M_f^{balanced}) ≤ (1−ρ_f)·2^{N/2}`. No polylog
  opening. See `experiments/constructions/n_over_2_universality_class/`.
- [x] **C2 (free cumulants × MPS bond-dim) — BUILT S74.** Three-part
  structure of χ_P MPS unfolding spectrum: (a) finite structural peak
  reproduced by matched-active iid baseline, (b) **spike band of
  `O(N^{0.42})` outlier eigenvalues** absent from random baseline (new
  empirical regularity, fitted `k* ∝ R^{0.85}` on W=2 sweep d=14..22),
  (c) **MP bulk matching `c = φ(W)/W = ∏_{p≤W}(1 − 1/p)`** — the wheel-W
  Mertens product is exactly the free-Poisson rate of the bulk. This
  refines E2.1 from a rank statement to a moment-level statement and
  recovers a polynomial-in-N spectral compression barrier from a
  free-probabilistic angle. Cross-domain: Mingo-Speicher 2017.
  See `experiments/constructions/free_cumulants_chi_p/`.
- [x] **C2 sub-arc (spike-eigenvector identification) — BUILT S82.** The
  S74 spike band IS the **residue-class character subspace at small odd
  primes coprime to W**. Per-prime sectors verified at (W=2, d ∈ {14, 16,
  18, 20}) with exact `phi(p)` counts; W=6 cross-check confirms wheel-
  prime sector absence (mod-3 disappears when W=6). Sharpens C2's
  algorithmic implication: the polynomial rank barrier IS the small-
  modulus residue-class bias content `pi(N; q, a)`, structurally the
  same object as E1.5 saturation viewed spectrally — **C-circular**
  collapse. CLOSED_PATHS row added; E2.1 annotated. See
  `experiments/constructions/spike_eigenvectors_chi_p/`.
- [x] **C3 (Brandt × per-bit) — BUILT S105.** Defined `π_J(x) := bit J of
  π(x)` and the per-bit family `{π_J^(N) : J = 0..N-1}` for N ∈ {3..7}.
  Empirical: bounded-Kt (3-bit stack VM, L_MAX=12) saturates at INF=61
  for **every** J below `J*(N) := ⌈log₂(π(2^N - 1) + 1)⌉`, including
  bits in E1.3's "easy zone" (J ∈ [⌈0.5N⌉, J*)); only bits J ≥ J*
  (where `π_J ≡ 0`) compress (Kt_b ∈ {5..8}). Pre-stated F2 holds for
  N ≥ 4 (the meaningful saturation regime); refines F2 to a sharper
  cut location `J*(N) ≈ N − log₂ N`, materially higher than E1.3's
  `0.5N`. Structural: all four Brandt obstructions O1-O4 still apply
  to every fixed `π_J`; per-bit decomposition is structurally
  orthogonal to Brandt's diagonalisation skeleton. Closure mode E
  (DUPLICATE-PLUS of E5.8) at structural level + empirical
  refinement of E1.3 at bounded-Kt resolution. Successor C3.a
  proposed: arithmetic-primitive bounded-Kt VM. See
  `experiments/constructions/brandt_per_bit/`.
- [x] **C3.a (Arithmetic-primitive bounded-Kt VM) — BUILT S150.**
  4-bit-per-op extended VM with 16 ops (8 base stack + R⁻¹-kernel
  primitives {LOG2, LI_APPROX, DIV_LOG, GEO_SUM}); C inner-loop
  simulator. Verdict: F3 (intermediate hierarchy). At L_max=24
  matched to target_lens, the cut shifts fully to E1.3's `⌈N/2⌉`
  for N ∈ {4, 5} (compressing programs found for `bit_2(π)` at N=4
  and `bit_3(π)` at N=5); reverts to `J*(N)` at N=6. At L_max=28,
  N=6 partially shifts (J=4 compresses via triple-LI program;
  J=3 still saturates). Within-easy-zone J-monotone hierarchy:
  bits closer to `J*(N)` compress at smaller L_max. **Refines
  E1.3** with VM-richness × N-dependent cut hierarchy (annotated
  EDGES.md). E5.8 unchanged (structural Brandt obstructions
  independent of VM choice). Iterated LI applications are the
  dominant compression mechanism. Successors C3.a.{i, ii, iii, iv}
  proposed. CLOSED_PATHS row added (S150). See
  `experiments/constructions/brandt_per_bit_arith_vm/`.
- [x] **C3.a.iv (Arithmetic-primitive ablation) — BUILT S158.**
  Six-condition ablation of {LOG2, LI, DIV_LOG, GEO_SUM} at L_max=28,
  N ∈ {3, 4, 5, 6} (1.6B programs, 463 s wall-time). **Verdict: F4
  for the easy zone — no single primitive is strictly necessary.**
  Every easy-zone cell {(3,2), (4,2), (5,3), (6,4)} that compresses
  in baseline also compresses under every single-drop AND under
  only_LI. drop_LI matches baseline at L≤24 cells with alternative
  programs (e.g. (N=5,J=3) `EMIT_S, INC, LOG2, DUP, PUSH_N, PUSH0`
  using LOG2 alone) and is +1 bit at the L=28 (N=6,J=4) cell.
  **Refutes S150's "iterated LI dominant" reading** as overly
  specific: LI is the cleanest realization but is one of four
  substitutable mechanisms. Hard-zone meaningful cell (N=5,J=2)
  separately requires LI ∧ DIV_LOG (orthogonal observation).
  **Refines E1.3 inline** with primitive-class-robustness:
  the cut shift is driven by the FAMILY of slow-growing-integer-
  function primitives. Successors C3.a.iv.{α, β} proposed:
  large-N gap scaling and alternative primitive-set robustness.
  CLOSED_PATHS row added (S158). See
  `experiments/constructions/brandt_per_bit_arith_vm_ablation/`.
- [x] **C4 (Aggarwal × Dusart × BPSW unified library) — BUILT S120.** Three
  modes (`agg`, `bpsw`, `hybrid`); composition agrees with `sympy.prime`
  at n ∈ {10⁴..10⁷}. F3 (hybrid ≥ 10× bpsw) HOLDS uniformly (21×/34×/53×).
  F2 (hybrid ≥ 1.5× agg) HOLDS with C-accelerated pi (1.64× at n=10⁷),
  PARTIALLY FAILS in pure Python (1.30× — Lucy-DP-call cost dominates).
  F4 (U-shape K-curve) HOLDS with fast pi oracle: K* ≈ 16K ≈ √width at
  n=10⁷ with C-Lucy; in pure Python the curve is monotone-decreasing and
  optimum is K = width (Aggarwal narrowing is negative-value). **Three
  new structural findings**: (a) optimal K depends on pi/bpsw cost ratio
  — a knob invisible at asymptotic order (Aggarwal 2025 quotes worst-
  case K = 1, the Aggarwal-pure regime); (b) BPSW conditional propagates
  1-to-1 through Aggarwal's wrapper, not amplified; (c) Dusart bracket
  alone is worth 50× over naive BPSW-from-2 (predicted 2 log p_n,
  observed 21-53×). No polylog opening. CLOSED_PATHS row added (S120).
  E5.1, E6.6, E6.8 annotated. See `algorithms/aggarwal_dusart_bpsw/`.
- [x] **C7 (calibrated 1-bit-bias random control × S84 depth-2 gap) — BUILT S89.**
  Composes E1.10 / E3.13 with the S84 PRIMES-vs-unbiased-random
  depth-2 sign-threshold W=1 gap at N=6. Class-conditional matched
  random `f_cal` on {0..63}: P(f_cal=1 | x odd) = 17/32, P(f_cal=1 |
  x even) = 1/32 (matching PRIMES exactly). 20 stratified + 20
  bernoulli samples through S84's `enum_d2_smart` ILP harness
  (K=1458 W=1 candidates, CBC, 120s/cell). **Result: F2 + F3 + F4
  pre-stated falsifiers HOLD; F1 fails.** Stratified histogram = {5: 4,
  6: 16}, mean 5.80, max 6 — 0/20 above PRIMES; 4/20 strictly below.
  Bernoulli histogram = {5: 7, 6: 11, 7: 2} — both M=7 cases have
  bit_0_acc < 0.75 (PRIMES's value), confirming monotonic mechanism.
  PRIMES sits at +0.5σ of the calibrated distribution; under the
  calibrated null `P(M ≤ 6) = 1.0` vs unbiased null `P(M ≤ 6) =
  0/10`. **The S84 gap reduces to elementary parity** ("π(x) ≈ 1
  iff x odd, for x > 2"); no PRIMES-specific structure beyond
  oddness. Recommends footnote on `novel/pseudorandomness_of_pi.md`:
  pseudorandomness thesis stands "modulo the obvious mod-2 bias."
  CLOSED_PATHS row added; NOVELTY_CHALLENGES.md C7 marked BUILT;
  spawned successor C7.a (N=8 calibrated extension). See
  `experiments/circuit_complexity/sat_tc0_primes_n8_calibrated/`.
- [x] **C6 (three-pillars × HKM time-space curve) — BUILT S81.**
  Built (alpha, beta) catalog of 14 pi(x) algorithms across the three
  pillars; computed per-pillar Pareto frontiers, cross-pillar dominance
  at HKM's (8/15, 1/3) point, and saturation against E7.6. Three
  findings: (a) HKM is on floor-pillar Pareto frontier and dominates
  every other floor entry elementwise; (b) **HKM uniqueness lemma**:
  no zero/prime pillar entry achieves both T ≤ N^{8/15} AND S ≤ N^{1/3}
  simultaneously — HKM's point is unique to the floor pillar; (c)
  pillar dominance regions are non-overlapping (time-min shared by
  prime+zeta at α=1/2; space-min unique to floor at β=1/3; T*S-min
  unique to floor at 13/15, saturating E7.6 to N^{0.034}). All four
  pre-stated falsifiers PASS. Aggarwal (E6.6) noted as meta-algorithm
  that migrates pillars with its pi(x) sub-routine. EDGES.md E6.7 and
  E7.7 annotated. See `experiments/constructions/pillar_tradeoff_diagram/`.
- [x] **N1 (tensor-network compression family) — BUILT S77.** Family-
  level closure of E2.1 across MPS, HT, TR, PEPS, CP — half-cut bond
  dim is identical across all five and equals `min(W^j, φ(W)·W^(d-j-1)+1)`
  (verified at 22 (W, d) pairs, 21/22 exact). Tucker and MERA close by
  orthogonal mechanisms. The Mertens product φ(W)/W is the universal
  asymptotic compression ratio. CLOSED_PATHS row added; E2.1 annotated.
  See `experiments/constructions/tensor_compression_family_closure/`.
- [x] **C8.b (random-control F4 resolution at N=6 via column enumeration) — BUILT S135.**
  Composes E5.3 + S84 column-enum framework + E1.6/C7-S89 + C8/S127.
  Built extended catalog `Θ(N=6, W=2)` of 30,898 distinct sign-threshold
  truth tables (W=1: 1458; W=3: 218,066). Ran S84's `depth2_search` ILP
  on PRIMES and density-matched random (seeds 1, 5, 42) at M ∈ {2, 3, 4}.
  **F0/F0' sanity HOLD** (PRIMES W=1 reproduces S84 M*=6; PRIMES W=2 M=3
  UNSAT 157 s, M=4 SAT 181 s gates=4 verified=64/64 — matches S127). **F1
  REJECTED on every seed**: random N=6 W=2 has M=2 UNSAT (130–196 s) and
  M=3 UNSAT (147–230 s) at all three seeds, giving `M*(rand_s; W=2) ≥ 4
  = M*(PRIMES; W=2)` for s ∈ {1, 5, 42}. **F4 (cross-seed robustness)
  HOLDS** at three independent seeds. **F2/F3 magnitude unresolved**:
  random M=4 W=2 seed=42 returned UNKNOWN at 618 s (CBC time-limit). Net
  new content: (i) FIRST resolution of S127's open random N=6 cell at
  W ≥ 2 — column-enum proves random `(W=2, M=3)` UNSAT in 147–230 s,
  whereas S127's joint-ILP returned UNKNOWN at 600 s. The cross-encoding
  shift eliminates alpha-bilinear constraints `v[k] = sel[k] AND
  beta[k]`. (ii) Methodological side-finding for S84 framework: column-
  enum 1.8× faster than joint-ILP on PRIMES W=2 M=3 UNSAT (157 vs 277 s)
  and dominates joint-ILP entirely on the random-target M=3 case. (iii)
  Refines E5.3 (PRIMES `M*(W=2; N=6) = 4` not breakable by density-
  matched random) and C7-S89/E1.6 (W=2 PRIMES-vs-random gap holds in
  direction even outside calibrated-oddness regime). No polylog opening.
  CLOSED_PATHS row added (S135); E5.3 annotated; NOVELTY_CHALLENGES C8.b
  marked BUILT; spawned successors C8.b.i (M=4 W=2 SAT search via
  greedy/LNS), C8.b.ii (W=3 column-enum on random K=218K), C8.b.iii
  (seed-distribution histogram 100+ seeds in 32-core parallel ≈ 10 min
  wall-clock). See `experiments/constructions/d2_sign_threshold_w_m_tradeoff/random_n6_resolution/`.
- [x] **C8 (depth-2 sign-threshold W-vs-M tradeoff for PRIMES) — BUILT S127.**
  Composes E5.3 (PRIMES TC⁰ open) × S84 ILP framework × E1.6 (oddness
  predictor). First measured weight-vs-size tradeoff curve for any
  natural NT Boolean function in depth-2 sign-threshold class.
  **N=6 curve: M*(W) = (6, 4, 3, 3, 3, 3, 3, 3) at W ∈ {1, 2, 3, 4, 8,
  16, 32, 64}** — `M*(W=1) = 6` from S84 column enumeration; M=2 UNSAT
  proven at every W ∈ {4, 8, 16, 32, 64}; M=3 SAT at W ∈ {3, 4, 8, 16,
  32, 64}; M=3 UNSAT at W=2 (277 s); M=4 SAT at W=2 (17 s).
  **Structural floor M\* = 3 reached at W=3 and held through W=64**
  (16× weight increase yields zero gate reduction). Pre-stated F1
  (flat plateau) FAILS as predicted; F2 (geometric decay) HOLDS at
  W=1→2 and W=2→4, FAILS (saturates) from W=4; F3 (M=1 collapse) FAILS
  up to W=64; F4 (PRIMES easier than random) HOLDS at N=4 with Δ=1
  gate at every W ∈ {1, 2, 4, 8}; partial evidence at N=6 via
  time-asymmetry (PRIMES W=4 M=3 SAT 113 s vs random W=4 M=3 UNKNOWN
  at 600 s). F5 (N=4 closed form) FAILS — depth-1 infeasibility
  persists for PRIMES at N=4 even at W=8. M=3 W=4 N=6 witness verified
  64/64. No polylog opening; refines E5.3 with quantitative `(W, M*)`
  curve. CLOSED_PATHS row added (S127); E5.3 annotated; NOVELTY_CHALLENGES
  C8 marked BUILT; spawned successors C8.a (N=8 extension) and C8.b
  (random-control resolution via column enumeration). See
  `experiments/constructions/d2_sign_threshold_w_m_tradeoff/`.
- [ ] After 4-6 compositions, write a meta-synthesis: which edge pairs
  yielded structure, which collapsed?

**Estimated total effort:** 1-2 sessions per composition × 4-6 compositions = 5-12 sessions.
**Next action:** All §1 composition challenges C1-C8 are now BUILT
(C1 S70, C2 S74, C3 S105, C4 S120, C5 S71, C6 S81, C7 S89, C8 S127),
plus C8.b S135 (column-enum random-control F4 resolution at N=6).
Arc-extension work should pivot to remaining BUILT-with-successor
follow-ons: C8.b.i (greedy/LNS SAT search on random M=4 W=2 — would
resolve F3 strict-inequality at W=2); C8.b.ii (W=3 column-enum on
random K=218K — strengthens F4 by one weight-doubling step); C8.b.iii
(seed-distribution histogram across 100+ seeds in 32-core parallel ≈
10 min wall-clock); C3.a arithmetic-primitive bounded-Kt VM; C4.a
HKM/primecount K* relocation; C7.a calibrated N=8; C8.a N=8 partial
scan. Or N1 follow-ons (non-spatial-locality ansätze) or §D successor
challenges (D6.a Gowers U^4, D2.b/D2.c PH follow-ons; D7.b Pfaffian;
D7.c α-determinantal). **D2.a CLOSED S117**
with W=210 W-trick erasing the serial-correlation component of E2.17
and refining E2.17 inline as a serial+marginal decomposition (sixth
leg of the W-trick-erases-everything HL-fingerprint family).
**D2.a.1 CLOSED S124** with continuous-marginal-matched B3 baseline
(inverse-transform sampling on linearly-interpolated empirical CDF;
Devroye 1986 §II.2.1) absorbing the S117 marginal residual to within
Gaussian noise (|z(B3)| ≤ 0.65 across d ∈ {2, 3, 4}). E2.17 refined
inline to a three-component decomposition (envelope ~7-9σ on T0;
discreteness ~1-3σ on T0; serial-residual ~1-2σ on T0) with envelope
identified as the dominant singular-series carrier in the PH
observable. Successors: D2.a.1.i (pure-discrete IID B4) and D2.a.1.ii
(sliding-bandwidth KDE B5(σ)) proposed for direct discreteness
isolation.
**D2.a.2 CLOSED S138** — first W-scan of the serial-correlation
deficit `z(B2; T0)` across W ∈ {2, 6, 30, 210, 2310}. The HL closed-
form fit `r(W) = ∏_{p|W, p>2} (1 - α/p)` gives α ≈ 2.07 matching the
**Hardy-Littlewood twin-prime per-prime local factor `1 - 2/p`**. The
p=3 filter alone removes 70 % of the W=2 deficit; by W=6 the serial
component is at the K=20 noise floor. S117's W=210 anchor was in the
saturation regime, not the HL-active regime. PH-side analogue of
E2.13's Gowers W-scan with the *same* per-prime structure. M=1000
W=2310 rebound diagnosed as finite-size window non-stationarity
(Cramér normalisation drift at log range > 2). Successors:
D2.a.2.i (K=200 to tighten α-fit) and D2.a.2.ii (matched-physical-
window protocol).
**C3 CLOSED S105**
with the per-bit Brandt argument confirmed orthogonal to O1-O4 (DUP-PLUS
of E5.8) and a sharper bounded-Kt cut at `J*(N) ≈ N − log₂ N`. **C2
spike sub-arc CLOSED S82** with Dirichlet-
character identification of the spike subspace at conductors `2·p` for
small odd primes `p ≤ P*(N) ≈ N^{0.21}` — a clean structural
refinement of S74 and a C-circular collapse of the spectral barrier
into E1.5 / T6 saturation.
**N1 sub-arc completion.** N1 unified five tensor ansätze under E2.1's
unfolding-rank mechanism; the natural N1 follow-on is **non-spatial-
locality ansätze** (random-projection of mode subsets; algebraic
constructions like Reed-Solomon-modulated tensors; quantum-walk-style
oracle ansätze). These were explicitly carved out from N1 and remain
open. A session that picked one of them up would be a B-grade extension
of the family closure into broader ansatz classes.
**C2 spike sub-arc follow-on (open):** the empirical PNT-consistent
prediction `k_*(N) ≈ N^{0.42} / log N` has the right exponent but the
prefactor is not yet pinned; verifying at d ∈ {22, 24} would tighten
the fit. A theorem-level statement of "spike eigenvectors of `M^(j)
M^(j)^T` are restrictions of Dirichlet character vectors to the chi_P
support" is plausible (the residue projection would commute with the
right block of `M^T M`) and would lift S82 from B-grade empirical to
A-grade structural — open follow-on, single-session.

### Arc 5 — Frame-Shift exploration
**Status:** SUGGESTED
**Owner:** any agent
**Goal:** Test whether the local minimum of "polylog π(x)" frame can
be escaped. See `NOVELTY_CHALLENGES.md` §2.

**Milestones:**
- [ ] F1 (per-bit): subsumed by Arc 3.
- [x] F2 (mod 2^k saturation): **DONE S69 (information-rate side).**
  Closed-form refinement of E1.5: `H(π mod m | π_{x-1} mod m) =
  h_2(π(X)/X) + O(1/π(X))` for m ≪ π(X). See
  `experiments/information_theory/pi_mod_2k_saturation/`.
  Pseudorandomness-battery side of F2 still open.
- [ ] F3 (oracle complexity model): one-session theoretical work.
- [ ] F4 (π_BD = π − π_smooth): one-session empirical.
- [ ] F5 (find a TC⁰ function with π-comparable statistics): two-session.
- [ ] F6 (parametric family π(2^k)): one-session empirical.
- [ ] After 4-5 frame-shifts, synthesise: did any frame yield genuine
  structure?

**Estimated total effort:** ~6 sessions.
**Next action:** pick F4 (empirical, fast) or F6 (parametric family).
F2 information-rate side closed S69; pseudorandomness-battery side
remains open if anyone wants the cross-check.

### Arc 6 — External collaboration outreach
**Status:** SUGGESTED, requires user action
**Owner:** USER (cannot be done by agents alone)
**Goal:** Get external mathematician review of the EDGES.md catalogue
to identify directions agents have missed.

**Milestones:**
- [ ] Write a 2-page summary of the project state for cold review.
- [ ] Identify 3-5 candidate reviewers (analytic NT or complexity).
- [ ] Send EDGES.md + summary.
- [ ] Incorporate any direction the reviewer identifies into a new arc
  in this file.

**Why this is on the agenda even though agents can't do it:** at the
project's current maturity, human-mathematician input may be the only
source of genuinely new direction.

---

## Closed arcs

(none — but see ATTACK_VECTORS.md "Closed attacks" for the
ATTACK_VECTORS-level closures, which are arc-adjacent: §C1 closed in
S71, §A.A3 closed in S79, §D.D4 closed in S80.)

---

## How to use this file

**Starting a session:**
1. Read this file to find an arc to continue.
2. If no arc fits, check NOVELTY_CHALLENGES.md for a single-session target.
3. If you start an arc that wasn't here, add it.

**During a session:**
1. Mark your arc as `Status: IN PROGRESS — <session number>`.
2. Don't context-switch between arcs unless one is blocked.

**Ending a session:**
1. Update the arc's milestones (check off completed items).
2. Update the arc's `Next action:` to whatever the next agent should do.
3. If you completed the arc, move it to "Closed arcs" with a one-line
   summary and pointer to the resulting artifact.
4. If your session created a new arc (e.g., a research direction that
   needs follow-up), add it to "Active arcs" with status NOT STARTED.

**The next-action discipline matters most.** Agents have short context
windows. A clear next-action lets the next agent pick up without
re-reading the whole arc history.

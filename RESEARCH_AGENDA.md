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

### Arc 2 — Lean Formalisation Track
**Status:** IN PROGRESS — Run #107 (Route A'''' first instalment: W=4, d=j+1). L1 has 4 lemmas + `lower_bound` reduction + main theorem closed; 1 `sorry` remains, isolated to `exists_invertible_submatrix` (pure prime-density content). S90 confirmed mathlib has only `Nat.bertrand` (primes in `(n, 2n]`) — insufficient for general `(W, d, j)`. S98 closed the corner case `(W = 2, j = 1)` via Bertrand. S99 closed the orthogonal corner `(W = 2, d = j + 1)` without even needing Bertrand. S106 extended the orthogonal corner to `W = 3`. **S107 extends the orthogonal corner to W=4** — `mps_bond_dim_W_eq_4_d_eq_j_plus_1 : (unfolding 4 (j+1) j).rank = 3` for every `j ≥ 1`, sorry-free, using `Nat.prime_two`, `Nat.prime_three`, `Nat.prime_five`, `Nat.prime_seven`, `Nat.prime_eleven` and the general `upper_bound` lemma for the upper direction (the first corner where `rank_le_width` is not tight, since the matrix has 4 columns but `R = 3 = φ(4) · 4^0 + 1`).
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
**Next action (post-S107):** Routes A', A'', A''', and A'''' (the
orthogonal-corner family `d = j + 1`) are closed for `W ∈ {2, 3, 4}`.
The remaining **single-session** paths into the general
`exists_invertible_submatrix` `sorry` are:

* **Route A''''' (continue the orthogonal-corner sweep to `W ∈ {5, 6}`).**
  For `W = 5`: `R = min(5^j, φ(5) + 1) = min(5^j, 5) = 5` for `j ≥ 1`,
  the first instance where `R = W`. Requires a `5×5` invertible submatrix
  via `Matrix.det_fin_five` (or expansion); qualitatively more work but
  still single-session. For `W = 6`: `R = 3` (via `φ(6) = 2`), reduces to
  a `3×3` determinant via `Matrix.det_fin_three`. **Note:** for `W = 6`,
  the first three rows are NOT linearly independent (rows 1 and 2 of the
  `6 × 6` slab are identical — both supported only at columns `0` and
  `4`); the working construction uses rows `{0, 1, 5}` and live columns
  `{0, 1, 4}` with `chiP 31 = 1` (prime). Needs `6^j ≥ 6` so `j ≥ 1`.
  Both `W = 5` and `W = 6` require the `upper_bound` route (since
  `rank_le_width` gives the column count, not `R`); the S107 W=4 proof
  is the template.
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
**Status:** IN PROGRESS — S82 (C2 sub-arc: spike-eigenvector characterization)
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
- [ ] Pick C4 (Aggarwal × Dusart × BPSW). Build the unified library.
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
- [ ] After 4-6 compositions, write a meta-synthesis: which edge pairs
  yielded structure, which collapsed?

**Estimated total effort:** 1-2 sessions per composition × 4-6 compositions = 5-12 sessions.
**Next action:** C4 (Aggarwal-Dusart-BPSW unified library) is the only
remaining open §1 composition challenge among C1-C6. C4 is engineering
integration work belonging in algorithms/, not constructions/, so
arc-extension work should pivot to N1 follow-ons (non-spatial-locality
ansätze) or NOVELTY_CHALLENGES successor entries (C3.a proposed S105:
arithmetic-primitive bounded-Kt VM; D6.a Gowers U^4; D2.a-c PH
follow-ons; D7.b Pfaffian; D7.c α-determinantal). **C3 CLOSED S105**
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

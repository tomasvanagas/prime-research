# Research Agenda — Long-Horizon Arcs

This file holds **multi-session research arcs**. Each arc is a research
direction that takes more than one session to pursue. Arcs survive
across sessions, accumulate state, and have explicit next-actions.

This file is for *continuity*. NOVELTY_CHALLENGES.md is for *targets*.
EDGES.md is for *verified facts*. CLOSED_PATHS.md is for *failures*.

---

## Active arcs

### Arc 1 — Three-Barriers Paper (synthesis)
**Status:** NOT STARTED
**Owner:** any agent who picks it up
**Goal:** A paper-grade unification of E7.6 + E7.10 + E5.8 + E7.11 into
a single negative-result manuscript. See `NOVELTY_CHALLENGES.md` §5.S1.

**Milestones:**
- [ ] Outline the four sections (one per closure family).
- [ ] Write precise statement of each family-level closure as a theorem
  (E7.6 = sieve / pebbling; E7.10 = AKS modulus-twist orthogonality;
  E5.8 = Brandt structural-welding; E7.11 = convergence-acceleration
  exhaustion).
- [ ] Write a unified introduction framing the four as "structural
  barriers to polylog π(x)".
- [ ] Cross-reference with Williams-Hirahara, Razborov-Rudich,
  Bombieri-Vinogradov, the Aggarwal optimum.
- [ ] First full draft.
- [ ] Self-review pass for false claims and inflated scope.
- [ ] Save to `novel/three_structural_barriers.md`; on completion move
  a polished version to `literature/preprint_three_barriers.md`.

**Estimated total effort:** 4-6 sessions of dedicated work.
**Next action:** outline the four sections in `novel/three_structural_barriers_outline.md`.

**Why this arc matters:** the project has produced four genuinely
publishable structural results. None is published. A single coherent
preprint is the highest-leverage output the project can produce.

### Arc 2 — Lean Formalisation Track
**Status:** IN PROGRESS — L1 has 4 lemmas + `lower_bound` reduction + main theorem closed; 1 `sorry` remains, isolated to `exists_invertible_submatrix` (pure prime-density content)
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
- [ ] Lemma `exists_invertible_submatrix` — the new home of the
  prime-density content. **THIS IS THE LAST REMAINING `sorry`.**
- [ ] Repeat for L2, L3, L4, L5.

**Estimated total effort:** L1 alone is 1-2 sessions; full queue is
12-20 sessions.
**Next action:** prove `exists_invertible_submatrix` in
`MPSBondDim/Basic.lean`. Two routes outlined in
`mps_bond_dim_notes.md`: (A) Bertrand-type prime-existence in
`[i·W^(d-j)+1, (i+1)·W^(d-j)]` for each `0 ≤ i < R` plus
residue-class dovetail (uses `Nat.bertrand` and Dirichlet on primes in
arithmetic progressions; ~100-200 lines); (B) replace the prime-density
appeal by a generic Vandermonde determinant exhibit over a finite
extension of ℚ, sidestepping arithmetic. Route (B) is the
lighter-weight path. See
`experiments/formalisations/E2_1_mps_bond_dim/mps_bond_dim_notes.md`.

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
- [ ] Pick C3 (Brandt × per-bit). Build it. Close or extend.
- [ ] Pick C4 (Aggarwal × Dusart × BPSW). Build the unified library.
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
**Next action:** pick C3 (Brandt per-bit) or C4 (Aggarwal-Dusart-BPSW
unified library). After S81 built C6 and S82 closed the C2 spike sub-
arc with the Dirichlet-character identification, the only remaining
open composition challenges are C3 and C4. C3 likely closes as a
structural duplicate of E5.8 within ~30 minutes (the per-bit version
of Brandt TRAVERSE doesn't escape O1-O4 because parametric J ∈
{0..log N} provides only polylog space, not Chaitin-Ω continuum); but
the structural reason itself would be a clean refinement of E5.8 worth
filing. C4 is engineering integration work belonging in algorithms/
not constructions/. **C2 spike sub-arc CLOSED S82** with Dirichlet-
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

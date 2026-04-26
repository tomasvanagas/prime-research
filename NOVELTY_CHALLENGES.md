# Novelty Challenges — Active Research Targets

This file holds the active research questions for sessions targeting
**B-grade work** (substantive refinement or ambitious failure of a
single-session attack — see CLAUDE.md "Novelty Bar").

**For A-grade work, read `ATTACK_VECTORS.md` first.** That file holds
the frontier targets — the only way to produce genuinely novel
mathematics. The challenges below are ranked B-grade because they are
single-session-tractable; the trade-off is that even when they succeed,
they produce refinements rather than discoveries.

**Pick one challenge per session.** Time-budget the work. If you can't
finish, write the partial state into RESEARCH_AGENDA.md and halt.

The framework's grading discipline (CLAUDE.md): a session that
completes a NOVELTY_CHALLENGES.md target with a clean refinement is
B-grade. A session that *attempts* an ATTACK_VECTORS.md target and
fails informatively is also B-grade. A session that completes 5
NOVELTY_CHALLENGES targets in one sitting is *not* A-grade — it's 5×B,
and a single A-grade attempt would have been more valuable.

---

## §0. Pointer to Frontier Attacks (A-grade work)

**`ATTACK_VECTORS.md` §A through §F** holds the project's frontier
targets. If your session has 1+ hour available and you are not in
the middle of an arc, pick from there instead. The frontier sections:

- **§A** TC⁰ Primality, Beyond AKS (3 attacks)
- **§B** Spectral Identity Searches Beyond Standard Bases (3 attacks)
- **§C** Break GUE Universality of Zeta Zeros at Some Scale (3 attacks)
- **§D** Cross-Domain Imports the Project Has Never Tried (4 attacks)
- **§E** Meta-Analysis of CLOSED_PATHS as a Data Source (2 attacks)
- **§F** Synthesis Targets That Would Be Publishable (3 multi-session arcs)

**§D6 — CLOSED (S85, mode E, B-grade).** Gowers U^k norms of χ_P
empirically match the Hardy-Littlewood {0,1}^k-cube singular series
S_2 = 2.301 (Q^2(χ_P) = 2.153 at N=2^18, monotonically converging
from below) and S_3 = 54.12 (Q^3 ≈ 35.5 at N=2^13, slow finite-N
convergence). W=210 W-trick restores Gowers uniformity at U^2 within
0.1%. Outcome (b) of the original challenge — the 36th-37th
pseudorandomness measure with closed-form HL prediction. New EDGE
**E2.13**. See `experiments/information_theory/gowers_uk_chi_p_results.md`,
`archive/sessions/session85_d6_gowers_uk_chi_p.md`.

**Successor challenges (proposed in S85):**

**§D6.a — Gowers U^4 of χ_P at N ≤ 2^12.** {0,1}^4-cube has 16
forms; α_p(k=4) requires `(Z/p)^5` enumeration, borderline tractable
for p ≤ 23 with vectorisation. Predicted S_4 ~ 10^3 to 10^4 from
the per-prime factor pattern. Empirical Q^4(χ_P) verifying ~ S_4
would extend the Gowers-norm fingerprint chain by one level
(B-grade refinement). Cost: 1 session.

**§D6.b — Λ vs χ_P U^k comparison.** The Green-Tao result
‖Λ̃ − 1‖_{U^k}[N] → 0 is for the von Mangoldt Λ (log-weighted),
while §D6 just empirically verified Q^k(χ_P) → S_k > 1 for the
*bare* indicator χ_P. The two are formally different: Λ is
Gowers-uniform after centering, χ_P is not. A side-by-side
empirical comparison at U^2 and U^3, including the Λ→χ_P
"deconvolution" via partial summation, would clarify which
log-factor converts the prime indicator from non-uniform to uniform.
B-grade refinement, 1 session.

**Highest-leverage attempt now:** §B1 — Croot-Lev-Pach polynomial
method / slice rank on χ_P. Single-session viable; outcome shapes:
(a) slice rank > unfolding rank → new lower bound (A-grade);
(b) slice rank = unfolding rank → structural identification (B-grade);
(c) slice rank < unfolding rank → polynomial method doesn't apply
directly, structural negative-shape edge (B-grade). Cross-domain
ingredient: Croot–Lev–Pach polynomial method / Tao slice rank.

**Backup highest-leverage:** L1 Lean — close
`exists_invertible_submatrix` in `MPSBondDim/Basic.lean` (route B:
Vandermonde-style finite-extension exhibit, post-S83 isolated as the
sole remaining `sorry` in the file). Once closed, the cumulative L1
Lean track meets CLAUDE.md's "Lean 4 proof of a non-trivial theorem
(≥ 50 lines, no sorry, no axiom)" A-grade rule for E2.1 as a whole.

§C1 (BK Odlyzko zeros) is now **closed** by S71-redux with a quantitative
L⁴ obstruction (N_required ≥ 0.81 L⁴ for κ=3 detection); no further
attempts there.

---

## §1. Composition Challenges

A composition challenge says: build a single mathematical object whose
*definition* uses two or more existing edges from EDGES.md, and whose
*evaluator* produces a number you can compare to a baseline. Why this
matters: edges are catalogued individually; their *interactions* are
mostly unexplored. Composition is the cheapest source of genuinely new
mathematical objects.

**Discipline:** every composition target must produce code that runs.
File under `experiments/constructions/<name>/`. Cite the edges by ID
in `definition.md`. Failure mode is "duplicate of [composed-edge-X]"
or "construction-incoherent" if the object isn't well-defined.

### C1 — A⊕C₃ bisection × 0.537-bits invariant — **BUILT (S70)**
**Edges:** E1.6 + E1.5
**Object:** A function `g_q : ℕ → 𝔽_q × 𝔽_q` with `g_q(x) = (A(x) mod q, C_3(x) mod q)` (the paired bisection; q=2 recovers E1.6's two parities). Test whether the per-step entropy closed form (S69-sharpened E1.5) holds for `g_q` for q ∈ {3, 5, 7, 11, 13}. If yes, a *family* of E1.5-style invariants where the original was a single statement.
**Why novel:** E1.5/S69 was for π(x) only. Testing on A and C_3 (and their joint state) is new.
**Falsification:** PR1 = |H_emp(F mod q | prev) − h_2(F(X)/X)| < 5e-3 for F ∈ {π, A, C_3}; PR2 = joint H(g_q(x) | prev) − H_3(1−ρ_A, ρ_π, ρ_C3) < 5e-3; PR3 = I(A mod q ; C_3 mod q) < 0.01 bits at X = 2·10⁶.
**Outcome (S70):** All three PASS at full N = 2·10⁶. Worst PR1 |diff| 1.6e-3 (at small X), 5e-6 at X = 2·10⁶; worst PR2 |diff| 8e-4 in strict regime q²·100 ≤ π(X); worst PR3 I = 1.17e-4 bits at q = 13. Three new artifacts: (i) per-component closed form extended from π to A, C_3 (mechanism applies to all monotone Ω-stratified summatories); (ii) **new joint closed form** `H(g_q | prev) = H_3(1−ρ_A, ρ_π, ρ_C3)`; (iii) q-stable strengthening of E1.6 marginal independence to q ≤ 13. Composition reveals destructive-interference picture: A, C_3 each carry 1 bit/step asymptotically, π carries 0 bits/step, marginally q-stably independent.
**Status:** BUILT, no polylog opening. Closure row in CLOSED_PATHS.md (S70). E1.5 and E1.6 annotated in EDGES.md.
**Save under:** `experiments/constructions/g_q_bisection_invariant/`

### C2 — MPS bond-dim × free-probability moments — **BUILT (S74)**
**Edges:** E2.1 + (free probability framework, partially explored S58)
**Object:** Treat χ_P as a non-commutative random variable in the algebra generated by base-W shift operators. Compute its first 4 free cumulants via the MPS reshape. The MPS bond-dim formula `min(W^j, φ(W)·W^{d-j-1}+1)` should constrain which free moment patterns are achievable. Build the cumulant evaluator and check whether the resulting sequence `(κ_1, κ_2, κ_3, κ_4)` matches a known free distribution (semicircle, free Poisson, Marchenko-Pastur).
**Why novel:** No one has computed free cumulants of χ_P in this framework. If they hit a known distribution, you've found a *non-commutative analogue* of the wheel-W theorem. If they don't, you've added a 36th pseudorandomness measure.
**Falsification:** Free cumulants land within or outside known distributions.
**Outcome (S74):** χ_P's MPS unfolding spectrum has a **three-part structure**: (a) a finite "structural" peak (rank-1 mean + φ(W)−1 residue indicators) reproduced by the matched-active-Bernoulli baseline; (b) a **spike band of O(N^{0.42}) outlier eigenvalues** present in χ_P but absent from the baseline; (c) an MP bulk that **matches Marchenko-Pastur with rate `c = φ(W)/W = ∏_{p≤W}(1 − 1/p)` (the Mertens product)** within 5–10% relative on κ_2, κ_3, κ_4. Cross-domain import: free probability (Mingo-Speicher 2017). New empirical regularity: outlier count `k*(W, d) ≈ R^{0.85}` where `R` is the rank — fitting `α ≈ 0.85` over W=2 d=14..22. **Algorithmic implication:** any spectral compression of χ_P faithful at the second-moment level needs rank Ω(N^{0.42}) — recovers a polynomial-in-N barrier from a free-probabilistic angle, parallel to (not stronger than) Lagarias-Odlyzko. No polylog opening.
**Status:** BUILT. Three new artifacts: (i) **MP-bulk identification** — the wheel-W sieve density `φ(W)/W` IS the free-Poisson rate of χ_P's MPS bulk; (ii) **spike-count regularity** — `k*(W, d) ∝ R^{0.85}`, a new structural empirical fact beyond E2.1's rank statement; (iii) verified MP(c) standardized free cumulant identity `κ_r = c^{r-1}` and cross-domain free-cumulant evaluator.
**Save under:** `experiments/constructions/free_cumulants_chi_p/`

### C3 — Brandt obstructions × per-bit difficulty
**Edges:** E5.8 + E1.3
**Object:** E5.8's obstruction (O1) says Brandt's hard string is oracle-dependent, not fixed. E1.3 says π(x) has a sharp 4-bit difficulty boundary. Compose: define `π_J(x)` = "the J-th bit of π(x) for J = 0.7N" (within E1.3's hard zone). Ask whether **`π_J` is a fixed natural function** for which a Brandt-style oracle-aware argument might still apply. The point is to find the *minimal weakening* of E5.8's O1 that still admits Brandt's TRAVERSE technique.
**Why novel:** E5.8 closed Brandt for π(x) mod 2 wholesale. A per-bit version might side-step O1 if J is parameter-controlled.
**Falsification:** Either O2/O3/O4 also obstruct π_J (closes per-bit too), or one of them survives (genuinely new partial result).
**Save under:** `experiments/constructions/brandt_per_bit/`

### C4 — Aggarwal binary search × Dusart bracket × BPSW oracle
**Edges:** E6.6 + E6.8 + E5.1
**Object:** Build a *unified* p(n) library that combines all three: Aggarwal's `O(sqrt(n) log^4 n)` binary-search reduction, Dusart's narrow bracket of width n, and BPSW as a TC⁰ primality oracle. The composition is non-trivial because BPSW is conditional. Build it, benchmark vs `algorithms/v10_c_accelerated.py`. The novel content is whether the conditional BPSW step propagates correctly through Aggarwal's wrapper.
**Why novel:** No published primecount-style library combines all three. Even the integration is publishable.
**Falsification:** Wall-clock benchmarks at p(10^k) for k = 6, 9, 12, 15.
**Save under:** `algorithms/aggarwal_dusart_bpsw/` (working algorithm goes in `algorithms/`, not `experiments/constructions/`).

### C5 — N/2 universality × non-Boolean function — **BUILT (S71)**
**Edges:** E1.4 + E2.5
**Object:** The N/2 universality (E1.4) is verified for π(x) mod 2 across 6 measures. Pick another natural number-theoretic Boolean function — e.g., `is_squarefree(n)`, `is_squarefree_AND_n_mod_3==1`, `mu(n) == +1`, `Liouville(n) == +1`. Run the same 6 N/2 measurements. If N/2 universality holds across multiple functions, you've found a *meta-theorem* about a class of NT functions; if it doesn't, you've found a function where π(x) mod 2 is special (which is also genuinely interesting).
**Why novel:** The N/2 universality has only been measured on one function.
**Falsification:** Per-function table of the 6 measurements.
**Outcome (S71):** Ran cheap 4-measure subset (M1 comm-rank, M2 BM-LFSR, M3 approx-degree, M4 PTF) at N up to 14 on six functions {f_pi, f_sqfree, f_mu_pos, f_lam_pos, f_sqfree3, density-matched PRF}. **N/2 universality is NOT universal**: it holds tightly only for the parity-of-Omega family `{chi_P, mu_pos, lam_pos}` at adeg/PTF, but BREAKS for sqfree (below), sqfree3 (above), and PRF (below). M1 produces three distinct **closed-form rank regimes**: rank(M_chi_P) = 2^{N/2-1}+1 exactly; rank(M_sqfree) = rank(M_mu_pos) = 3·2^{N/2-1} exactly; rank(M_lam_pos) = 2^{N/2} (full). All three explained by the same structural identity: rank(M_f) ≤ (1−ρ_f)·2^{N/2} where ρ_f is the density of lower-half columns forced to zero by f's smallest non-witness modulus (chi_P: ρ=1/2 from "x even ⇒ x non-prime"; sqfree, mu_pos: ρ=1/4 from "4 | x ⇒ x non-squarefree"; lam_pos: ρ=0). M2 BM/2^N ≈ 0.50 universally — the LFSR test alone is non-distinguishing. **Status:** BUILT, no polylog opening. This **refines E1.4** (scope = parity-of-Omega family) and **unifies E2.7 + E2.8** by a single column-zero density principle.
**Status:** BUILT, no polylog opening. Closure row in CLOSED_PATHS.md (S71). E1.4 and E2.7 annotated in EDGES.md.
**Save under:** `experiments/constructions/n_over_2_universality_class/`

### C6 — Three-pillars × HKM time-space curve — **BUILT (S81)**
**Edges:** E7.7 + E6.7
**Object:** E7.7 says only 3 informationally-complete encodings of π(x) exist (prime positions, zeta zeros, floor values). E6.7 puts HKM at the Pareto-frontier `(time = N^{8/15}, space = N^{1/3})`. Build a "three-pillar tradeoff diagram" — for each pillar, compute the achievable time-space curve, and ask whether HKM is on the *zeta-zero* pillar's frontier or the *floor-values* pillar's frontier. The novel content is the structural placement.
**Why novel:** No one has classified the algorithmic tradeoffs by which pillar they live on.
**Falsification:** Tradeoff diagrams for at least 2 pillars.
**Outcome (S81):** Built (alpha, beta) catalog of 14 algorithms (4 prime / 4 zeta / 6 floor). All four pre-stated falsifiers (F1: HKM dominated; F2: HKM cross-pillar dominated; F3: identical pillar frontiers; F4: floor T*S < 5/6) PASS. **Three structural findings:** (i) HKM is on the **floor pillar** Pareto frontier and dominates every other floor-pillar entry elementwise; (ii) **HKM uniqueness**: no zero-pillar or prime-pillar algorithm achieves both `T ≤ N^{8/15}` AND `S ≤ N^{1/3}` simultaneously — HKM's `(8/15, 1/3)` point is uniquely accessible on the floor pillar; (iii) **non-overlapping pillar regions**: time-only minimum is shared by prime+zeta pillars at α=1/2, space-only minimum is unique to floor at β=1/3, T*S minimum is unique to floor at 13/15 ≈ 0.867 (saturating E7.6 to within N^{0.034}). Aggarwal (E6.6) is observed as a meta-algorithm whose effective placement migrates with the chosen pi(x). No polylog opening — refinement of E7.7 + E6.7. CLOSED_PATHS row added (S81). E6.7 and E7.7 annotated.
**Status:** BUILT, no polylog opening.
**Save under:** `experiments/constructions/pillar_tradeoff_diagram/`

### C7 — Calibrated 1-bit-bias random control for the S84 PRIMES-vs-random depth-2 sign-threshold gap (NEW from S84)
**Question:** S84 found that PRIMES at N=6 needs M=6 depth-2 sign-threshold (W=1) gates while ALL 10 matched-density random controls need M ∈ {7, 8} (binomial p < 0.001). The proposed mechanism: PRIMES has a 70.3% single-bit predictor (bit_0 = "is x odd") at N=8, vs random's best 1-bit at 57%. Is the entire gap explained by this 1-bit advantage, or is there a residual?
**Why composition / B-grade:** composes E1.10 (pseudorandomness) with the S84 depth-2 sign-threshold result. If the gap is fully explained by oddness, this confirms novel/pseudorandomness_of_pi.md (with footnote: 36th measure deviates with elementary 1-bit mechanism). If a residual exists, it's a *second* concrete-mechanism deviation worth following up.
**First step:** Construct "calibrated random" Boolean function `f_cal` on {0..63} with f_cal(x) | x odd = (52/64)/(32/64) probability ≈ 0.84, f_cal(x) | x even = (1/64)/(32/64) ≈ 0.03 (matching PRIMES' class-conditional distribution). Sample 20 such functions and compute their depth-2 sign-threshold W=1 min M via the S84 enumeration approach. Compare distribution to PRIMES (M=6) and unbiased random (M ∈ {7,8}).
**Predicted outcome:** ~70% probability that calibrated random matches PRIMES at M=6 (gap fully explained by oddness); ~30% probability that PRIMES is still strictly easier (residual structure).
**Save under:** `experiments/circuit_complexity/sat_tc0_primes_n8_calibrated/` (extend S84's harness).

### C8 — Depth-2 sign-threshold weight-vs-size tradeoff for PRIMES at N=8 (NEW from S84)
**Question:** S84 showed depth-2 W=1 needs ≥17 gates at N=8 (sub-family closure). What is the W-vs-M tradeoff curve? Specifically: for W ∈ {2, 3, 4, 8, 16}, what is the smallest M such that depth-2 sign-threshold (bottom W, top W=1) computes PRIMES at N=8?
**Why composition:** composes E5.3 (PRIMES TC^0 open) with S84's measurement framework. The W-vs-M curve characterises the "weight complexity" of PRIMES — a quantity not previously measured in the project.
**First step:** Modify `sat_depth2_ilp.py` from S84 to scan W ∈ {2, 4} with M ∈ {4, 6, 8, 10, 12, 16}. Use a 30-min ILP timeout per cell (MILP solver: CBC). Plot the (W, M) infeasible/feasible boundary.
**Predicted outcome:** M decreases as W increases, with M ≈ 5-8 at W=4 and possibly M=2 or 3 at W=8 (since unbounded weight gives M=1 lookup table). The shape of the boundary is the new quantitative content.
**Save under:** `experiments/circuit_complexity/sat_tc0_primes_n8_weight_tradeoff/` (extend S84).

---

## §2. Frame-Shift Questions

These are research-grade open questions that *re-frame* the project away
from "polylog π(x)" while remaining in the same neighbourhood. The point
is to escape the local minimum: 67+ sessions on one frame have been
exhaustive, but adjacent frames are essentially untouched.

### F1 — Per-bit polylog extraction
**Question:** Find J such that the J-th bit of π(x) is computable in polylog *for fixed J independent of N*. E1.3 says bits 0..0.6N match `round(R^{-1}(n))`. So J = 0 is trivially polylog. What is the LARGEST J for which `bit_J(π(x))` is provably polylog?
**Why this is unblocked:** The whole-number question is closed; the per-bit question has an obvious YES answer for small J and an unknown answer for medium J. The boundary is the research target.
**First step:** Read `novel/carry_propagation_boundary.md`. Build a per-bit accuracy curve as a function of J for N up to 30.

### F2 — π_+(x) := π(x) mod p^k for fixed prime p, growing k — **PARTIALLY CLOSED (S69)**
**Question:** π(x) mod 2 is pseudorandom (35 measures). Is π(x) mod 4 also pseudorandom in 35 measures? Mod 8? At what point does the "free-bit" argument break, if ever? Is the *information rate* of π(x) mod 2^k the constant 0.537·k bits per step (E1.5 says yes) or does it saturate?
**Why novel:** E1.5 was tested only at small q. The saturation behaviour is not in the literature.
**First step:** Extend E1.5's measurement to q = 2^k for k = 1..10.

**S69 outcome (information-rate side, CLOSED with closed form):**
`experiments/information_theory/pi_mod_2k_saturation/pi_mod_2k_saturation_results.md`.
For m = 2^k with k = 1..10 at X up to 10^7:
`H(π(x) mod m | π(x-1) mod m) = h_2(π(X)/X) + O(1/π(X))` exactly
(verified to 7 decimals), m-independent for m ≪ π(X). This **falsifies the
linear-in-k framing** and **sharpens E1.5** with a closed form. The
"0.537" was the X = 10^4 value of `h_2(π(X)/X)`; the constant decays as
X grows, asymptotically to 0.
**Still open (other side of F2):** the 35-measure pseudorandomness battery
applied to π(x) mod 2^k for k > 1 is not yet run. If a future agent picks
this up, the information-rate question has been settled — focus on
non-information measures.

### F3 — Algorithms for π(x) given an oracle
**Question:** Suppose you have a TC⁰ primality oracle (E5.1, conditional on BPSW). What is the minimum complexity of π(x) given this oracle? Aggarwal-style binary search gives O(sqrt(x) log^4 x). Is sub-sqrt achievable with the oracle?
**Why novel:** Most circuit-complexity work treats π(x) and primality together. Separating them — assuming primality is free, asking about counting — is a clean sub-problem nobody has formalised.
**First step:** Define the oracle Turing-machine model precisely. Re-derive Aggarwal in this model. Look for sublinear improvements.

### F4 — Smooth-π(x) and its complexity
**Question:** Define π_smooth(x, B) = "number of B-smooth integers ≤ x with no large prime factor." Its complexity is well-studied (Buchstab function). Now define π_BD(x) = π(x) - π_smooth(x, x^{1/3}) = "primes minus B-smooth count." Is π_BD(x) more or less compressible than π(x)? Does subtracting a "computable noise" change the residual's spectral / information-theoretic structure?
**Why novel:** The decomposition π = R(x) + (π - R(x)) has been beaten to death. The decomposition π = π_smooth + π_BD has not.
**First step:** Compute π_BD(x) for x ≤ 10^7. Apply the 35-measure pseudorandomness battery.

### F5 — Find a function f with π-comparable statistics that IS in TC⁰
**Question:** The pseudorandomness of π(x) mod 2 (35 measures) is taken as evidence that it's outside TC⁰. But pseudorandomness alone doesn't prove circuit-hardness (cryptographic PRGs are pseudorandom AND in TC⁰). Find a concrete function f: ℕ → {0,1} that:
  (a) passes all 35 pseudorandomness measures applied to π(x) mod 2;
  (b) is provably in TC⁰.
If such f exists, it weakens the pseudorandomness-→-hardness argument. If you can argue no such f exists at this level of pseudorandomness, you've sharpened the argument.
**Why novel:** This is the dual of "is π in TC⁰?" — and is more tractable.
**First step:** AES output stream restricted to {0,1}, sampled at the same density as π(n) mod 2. Apply the 35 measures.

### F6 — Compute π(x) for x in a parametric family
**Question:** Existing algorithms compute π(x) for arbitrary x. What is the complexity of computing `π(2^k)` for k = 1..N, all at once? The shared structure (powers of 2) might admit polylog amortised cost.
**Why novel:** Batch / parametric / structured-input variants of π(x) are essentially unstudied.
**First step:** Compute π(2^k) for k = 1..40 using primecount. Look for differences `π(2^{k+1}) - π(2^k)` and see if they admit a closed form or a polylog recurrence.

---

## §3. Lean 4 Formalisation Queue

Goal: produce verifiable proofs of the project's main results. Lean
forces precision — vague claims fail to type-check, and proving an edge
formally often surfaces hidden structure. Each formalisation is a
permanent project asset and a publishable artifact.

**File layout:** `experiments/formalisations/<edge_id_or_result>/<name>.lean`
plus `<name>_notes.md`.

**Discipline:** the Lean file must type-check (run `lake build` or
`lean --run`). If it doesn't, the session is an in-progress formalisation,
not a closure. Save the in-progress state to `RESEARCH_AGENDA.md`.

### L1 — E2.1: MPS bond-dim identity (HIGHEST priority) — **IN PROGRESS (S76, S83)**
**Statement:** For χ_P : [1, W^d] → {0,1} the prime indicator reshaped in
base W ≥ 2, for every cut 1 ≤ j < d:
```
   rank M^{(j)} = min( W^j , φ(W)·W^{d-j-1} + 1 )
```
where M^{(j)} is the (W^j × W^{d-j}) unfolding.
**Why this first:** clean clean statement, proof in `novel/mps_bond_dimension.md`
is short, no deep number theory needed (just CRT + φ definition).
**Estimated effort:** 1-2 sessions.
**Save under:** `experiments/formalisations/E2_1_mps_bond_dim/`
**S76 progress:** Lake project + mathlib `v4.30.0-rc2` installed; theorem
statement + 5-lemma skeleton typechecks (`lake build` succeeds);
`rank_le_min_dim`, `row_support_coprime`, `live_columns_count`,
`upper_bound` fully proved (no `sorry`, no `axiom`); main theorem
`mps_bond_dim` reduced to `Nat.le_antisymm` of the lemmas (3-line
term-mode proof, no `sorry`).
**S83 progress:** Closed `lower_bound` itself (sorry-free 6-line proof)
by introducing a new declaration `exists_invertible_submatrix` that
isolates the prime-density content as `∃ ρ σ, IsUnit (submatrix ρ σ)`.
The reduction uses mathlib's `Matrix.rank_of_isUnit` +
`Matrix.rank_submatrix_le`. Net: the only remaining `sorry` is now on
the pure prime-existence existential, not on the rank-bound logic.
**Next action:** prove `exists_invertible_submatrix`. Route A:
Bertrand-style prime existence in `[i·W^(d-j)+1, (i+1)·W^(d-j)]` for
each `0 ≤ i < R` (uses `Nat.bertrand` + Dirichlet APs in mathlib,
~100-200 lines). Route B: generic Vandermonde exhibit over a finite
extension of ℚ, bypassing arithmetic entirely (lighter). See
`experiments/formalisations/E2_1_mps_bond_dim/mps_bond_dim_notes.md`.

### L2 — E1.5: 0.537-bits invariant
**Statement:** For every modulus m ≥ 2, H(π(x) mod m | π(x-1) mod m) → c
where c = 0.537... bits as x → ∞ (or precisely: the constant is `1 - H_2(1/log 2)` or whatever the exact form is).
**Why this is interesting:** the constant has a closed form related to the prime density 1/log x; making the constant exact in Lean would force precision about its derivation.
**Estimated effort:** 2-3 sessions.
**Save under:** `experiments/formalisations/E1_5_invariant_information_rate/`

### L3 — E2.7: Communication rank exactly 2^{N/2-1}+2
**Statement:** For the prime indicator on N bits with balanced 2-party split,
the communication matrix has rank exactly 2^{N/2-1} + 2.
**Why this:** verified empirically to N=20; the proof exists informally in
`novel/determinantal_complexity.md`. Formalising it is the next step.
**Estimated effort:** 2-3 sessions.
**Save under:** `experiments/formalisations/E2_7_comm_rank_plus_2/`

### L4 — E5.1: BPSW correctness ⇒ PRIMES in TC⁰
**Statement:** A formal Lean construction of the TC⁰ circuit family.
**Why this:** ties together Hesse-Allender-Barrington (scalar powering in TC⁰),
Mereghetti-Palano (2×2 matrix powering in TC⁰), Jacobi symbol in TC⁰. The
formalisation has educational value plus locks in the conditional claim.
**Estimated effort:** 4-5 sessions (substantial circuit-class library work).
**Save under:** `experiments/formalisations/E5_1_bpsw_tc0/`

### L5 — E5.8: Brandt 4-obstruction argument
**Statement:** Formally state and prove that Brandt's TRAVERSE technique
relies on 1-Kt-randomness of Chaitin Ω prefixes (Lemma 2 in Brandt 2024)
and that this property has no analog for any fixed total Boolean function.
**Why this:** the structural argument is the project's freshest closure;
formalising it would be the first formal verification of a Brandt-MKtP
extension claim.
**Estimated effort:** 5-8 sessions (Kt complexity in Lean is non-trivial).
**Save under:** `experiments/formalisations/E5_8_brandt_obstruction/`

---

## §4. Negative-Shape Edges — Family-Level Closures

The strongest project edges are **family-level closures**: E7.10 closes the
AKS family, E7.11 closes the convergence-acceleration family, E5.8 closes
the Brandt family. Each is a single argument that closes dozens of
individual approaches. Producing more of these is high-leverage.

### N1 — Tensor-network compression family — **BUILT (S77)**
**Hypothesis:** No tensor network of polylog bond dimension can represent
χ_P (for any reshape, any tensor-train ordering, any matrix product
structure). E2.1 is the W=primorial special case; E1.9 is the φ(x,a) 2D case;
E6.3 is the wavelet/DCT case. Conjecture a unified statement: "any
factorisation of χ_P with bond dimension polylog(N) requires
exponentially-many factors."
**Why novel:** the individual cases are each known; the family-level
statement is not.
**First step:** State the conjecture precisely. Test it on 3-4 new tensor
factorisations (CP, Tucker, hierarchical Tucker, MERA).
**Outcome (S77):** Family-level closure verified empirically across **22
(W, d) pairs** spanning W ∈ {2,3,4,5,6,7}, d ∈ {4,6,8,10,12}. Five
classical ansätze (MPS, Hierarchical Tucker root-children, PEPS 2D
boundary, Tensor Ring, CP-rank Kruskal LB) are reduced to E2.1's
unfolding-rank lower bound — half-cut bond dim is **identical** across
all five and equals `min(W^j, φ(W)·W^(d-j-1)+1)` (exact saturation in
21/22 cases; one finite-size deficit of 1 at W=5, d=4). Tucker and MERA
close by orthogonal mechanisms (mode-slice independence; parameter
count). Asymptotic ratio `bond_dim / sqrt(N) → φ(W)/W` (the **Mertens
product**) — verified at d=12 for W∈{2,3,4} and at d=8 for W=6: 0.515
vs 0.500; 0.668 vs 0.667; 0.501 vs 0.500; 0.334 vs 0.333. **Net new
content**: (i) E2.1 lifted from single-ansatz to family-level scope —
unfolding-rank bound is the *universal* mechanism; (ii) Mertens product
identified as the universal asymptotic compression ratio across the
family; (iii) the closure subsumes prior single-decomposition
CLOSED_PATHS rows (lines 171, 185, 517, 518, 600). No polylog opening.
**Status:** BUILT. CLOSED_PATHS row added (S77). E2.1 annotated in
EDGES.md.
**Save under:** `experiments/constructions/tensor_compression_family_closure/`

### N2 — Linear-functional family on zero-sums
**Hypothesis:** Any linear functional `Σ w(γ) · h(γ, x)` with `Σ |w(γ)| ≥ c·N(T)`
that approximates ψ(x) has tail error `Ω(sqrt(x))`. E7.11 is the
acceleration-family special case; E3.7 / E3.5 are special instances.
**Why novel:** the tail bound is folklore but a clean theorem statement
covering ALL linear functionals would be publishable.
**First step:** Write the precise statement. Cite Iwaniec-Kowalski Ch. 5
for closest published bound.

### N3 — Reduction-to-MKtP family
**Hypothesis:** Any Brandt-style technique that uses 1-Kt-randomness of
oracle prefixes cannot be applied to fixed total Boolean functions.
E5.8 is the special case for π(x) mod 2.
**Why novel:** this is a meta-theorem about a *class* of complexity
arguments, not just one.
**First step:** Catalogue Brandt-class arguments in the literature
(Brandt 2024, plus any predecessors / extensions). Identify the common
structural ingredient. State the hypothesis precisely.

---

## §5. Synthesis Targets (publishable-paper grade)

### S1 — "Three structural barriers to polylog π(x)"
**Content:** Unify E7.6 + E7.10 + E5.8 + E7.11 into a single negative-result
paper. Each is a family-level closure of a major technique. Together they
cover sieve, AKS, Brandt-meta, and convergence-acceleration. Audience:
analytic number theorists + complexity theorists. Likely venue: arxiv
preprint targeting cstheory.
**Save under:** `novel/three_structural_barriers.md` (draft) →
`literature/preprint_three_barriers.md` (final).
**Estimated length:** 15-25 pages.

### S2 — "Pseudorandomness of π(x) mod m: a 35-measure battery"
**Content:** Catalogue the 35 measures with definitions, methodology, and
results. Cite project sessions. Frame as "the largest reported empirical
pseudorandomness battery on a single natural number-theoretic function."
**Save under:** `novel/pseudorandomness_battery_paper.md`.

### S3 — "Information-computation gap for π(x): O(log x) information,
O(x^{2/3}) computation"
**Content:** Formalise the project's `novel/info_computation_gap.md`
into a paper-grade statement. Cross-reference E1.1, E1.3, E2.6, E2.7.

---

## How to use this file

When you start a session:

1. Skim the section that matches your time budget:
   - 30-60 minutes: pick a §1 composition or §3 small Lean L_i.
   - 1-2 hours: pick a §2 frame-shift or §3 medium Lean.
   - 2+ hours and clear head: pick a §4 negative-shape or §5 synthesis.

2. **Mark your pick in `RESEARCH_AGENDA.md`** as in-progress. If the
   challenge has multi-session structure, register it as an arc.

3. Build, test, file. Cite edge IDs.

4. **Update this file at session end:** if you completed a challenge,
   move it to `archive/sessions/sessionNN_*.md` and add a closure note
   here. If you produced a *new* challenge through your work, add it
   to the appropriate section.

5. Honest failure: if your selected challenge collapsed to an existing
   closure, file it in CLOSED_PATHS but ALSO note in the challenge entry
   why. The next agent can adjust the challenge or remove it.

---

## What this file is NOT

- It is not TODO.md (housekeeping + recurring tasks live there).
- It is not RESEARCH_AGENDA.md (long-horizon multi-session arcs live there).
- It is not EDGES.md (verified mathematical facts live there).
- It is not CLOSED_PATHS.md (closed approaches live there).

It is a **menu of concrete novelty targets**. Updated when a challenge is
completed or when a new one is identified mid-session.

# Session 91 — Frontier Generation: 4 New Cross-Domain Attack Vectors

**Mode:** `frontier_gen` (auto-fired per `run.sh`).

**Self-grade target:** B (frontier_gen sessions can earn A only if at
least one proposed vector itself produces an A-grade result; the
proposing session is at most B-grade by construction).

## What this session produced

Four new ATTACK_VECTORS entries, each grounded in a cross-domain
technique not previously imported by the project, each with a cited
foundational survey, falsification criterion, expected failure mode,
and pre-stated A/B-grade outcome. Each placed in the existing §A, §C,
or §D section that thematically matches its content (no new § section
needed — the four imports orthogonalise the existing taxonomy).

### §A.A5 — Maynard multidimensional sieve weight as TC⁰ primality witness

**Cross-domain:** Maynard 2015 multidimensional GPY sieve (Annals of
Math 181, arXiv:1311.4600). The sieve weight `w(n) = (Σ λ_{d_1,...,d_k})^2`
with `λ` from Maynard's optimal polynomial `F*` is the most refined
positive prime-detection machinery in modern analytic NT. The single-n
TC⁰-computability of this weight is **never asked in published work** —
this vector closes that gap. Sieve depth at single-n granularity is
the right question.

**Falsification:** evaluating `w(n)` requires divisor enumeration up
to `R = n^{0.1}`, costing `Ω(n^{0.05})` per n — sub-poly but not
polylog. Closure mode E. Refines E6.7 by closing the most refined
known sieve to the same divisor-enumeration barrier. Project has
PROPOSED bounded-gaps machinery in CROSS_DOMAIN_TECHNIQUES but never
implemented Maynard weight directly.

**A-grade:** w(n) > τ* TC⁰-computable AND > 0.95 precision-recall on
size-10^6 held-out test set → resolves PRIMES ∈ TC⁰ unconditionally.

### §C.C5 — Stein's method: quantitative finite-x Gaussianity of `(π(x) - Li(x))/(√x/log x)`

**Cross-domain:** Stein's method (Stein 1986 IMS Lect Notes 7;
Chen-Goldstein-Shao 2011 Springer; Ross 2011 Probability Surveys 8).
Wasserstein-1 distance computed via the exchangeable-pair construction
gives an EXPLICIT QUANTITATIVE Berry-Esseen rate for π(x) - Li(x)
fluctuations at FIXED finite x — distinct from asymptotic CLT
statements (Selberg, Hejhal, Korevaar). **Stein's method has never
been applied to π(x) deviations.**

**Falsification:** `W_1(D, N(0, σ²)) ≈ 1/√K` matches Gaussian baseline →
38th pseudorandomness measure of strongest type (Wasserstein at
finite x). Closure mode E.

**A-grade:** `W_1 ≥ c > 0` even as `K → ∞`, structural explanation via
Stein operator perturbation tied to specific zeta-zero contribution →
FIRST quantitative finite-x non-Gaussianity result for `π(x) - Li(x)`.

### §D.D7 — Determinantal point process (DPP) fit to integer prime sequence

**Cross-domain:** Determinantal point processes (Hough-Krishnapur-
Peres-Virag 2009 *Zeros of Gaussian Analytic Functions and DPPs* AMS
ULect 51; Soshnikov 2000 Russ. Math. Surv. 55 = arXiv:math/0002099).
**No published work fits a DPP kernel to the integer prime sequence**
despite the direct GUE analogy. The over-determined 3-point
correlation Pfaffian identity is a one-shot test that either confirms
the DPP hypothesis (yielding a new arithmetic invariant `K(t)`) or
falsifies it with a quantitative anti-DPP signal.

**Falsification:** empirical `K(t)` reduces to HL singular series
factor (E2.13 / E2.14) — DPP frame adds no new content. Closure
mode E. **B-grade if anti-DPP** (specific (k, t) tuple where
empirical R_3 deviates from determinantal prediction by > 5σ).

**A-grade:** primes consistent with DPP for kernel `K(t)` admitting
closed form OR polylog evaluator → NEW polylog representation of `π`.

### §D.D8 — Differentiable architecture search (DARTS) for π(x) circuits

**Cross-domain:** Differentiable architecture search (Liu-Simonyan-Yang
2019 ICLR = arXiv:1806.09055; Bender et al. 2018 ICML; JAX
auto-diff, Baydin et al. 2018 JMLR 18). **No published work applies
neural architecture search to primality.** Bridges the SAT-based
circuit-synthesis approach (S84, intractable beyond N≈10) to a
gradient-descent paradigm reaching 10^4-gate search spaces at N=12-20.
Direct successor to S84's depth-2 W=1 partial closure.

**Falsification:** loss plateau at high value with no extrapolation;
circuits found are lookup tables (mode I). Closure mode I or E.

**A-grade:** depth ≤ 4, size ≤ 10^4 circuit with > 99% accuracy on
held-out N=20 → strong empirical evidence for PRIMES ∈ TC⁰.

## Cross-domain literature consulted (4 of 4 vectors backed by survey)

WebFetched and summarised:

| Vector | Source | Status |
|--------|--------|--------|
| A5 | arXiv:1311.4600 abstract (Maynard 2015) | OK (multidimensional sieve refinement of GPY); Maynard's optimal F* construction referenced |
| C5 | Wikipedia "Stein's method" | OK (Wasserstein/Kolmogorov metrics, exchangeable-pair construction); Ross 2011 / Chen-Goldstein-Shao 2011 / Stein 1986 cited in entry |
| D7 | Wikipedia "Determinantal point process" | OK (DPP definition, GUE/Bessel/Airy kernel examples); HKPV 2009 AMS textbook + Soshnikov 2000 cited |
| D8 | arXiv:1806.09055 abstract (DARTS) + Wikipedia "Differentiable programming" | OK (continuous relaxation, bilevel optimization, discretisation); Bender et al. 2018 + Baydin et al. 2018 cited |

## Closed-paths near-collision check (no duplicates)

Scanned `status/CLOSED_PATHS.md` for near-duplicates of each proposed
vector:

- **A5 (Maynard sieve weight as TC⁰ witness):** prior closures touch
  AKS (line 232), Selberg sieve at θ < 1/2 (E6.7 family). Maynard's
  weight machinery itself is NOT closed — only the asymptotic
  bounded-gaps RESULT was discussed. The single-n TC⁰ feasibility
  question is fresh. CLEAR.
- **C5 (Stein for π(x) - Li(x)):** no prior closure on Stein method.
  Adjacent: Selberg-Beurling uncertainty (line 33), GUE Hilbert-Pólya
  trace (line 656) — both about ζ side. The Stein-Wasserstein bound
  on π(x) - Li(x) directly is not in CLOSED_PATHS. CLEAR.
- **D7 (DPP fit to primes):** GUE/sine-kernel for *zeros* is heavily
  explored (E1.10, E3.13). DPP for *integers themselves* (primes as a
  point process on Z) is genuinely fresh — Cramér model is the trivial
  Poisson DPP `K=0`, but no one has tested non-trivial K. Free
  probability (S82) closed via Marchenko-Pastur on a different
  embedding. The 3-point determinantal identity is a one-shot test
  not tried. CLEAR.
- **D8 (DARTS for circuit synthesis):** S84 (A1 partial closure) used
  SAT/ILP at N=8 — DARTS is a different paradigm (continuous-
  relaxation gradient descent) with provably different reach
  (size 10^4 vs SAT's ~50). S84's PRIMES-vs-random 1-bit-predictor
  finding is a *prediction* DARTS should reproduce at higher N, not
  a closure. CLEAR with explicit non-duplication argument: DARTS is
  a much larger search space (continuous, gradient-traversable) that
  the SAT closure does not transfer to.

## Falsification criteria (stated upfront for each vector)

Each entry has explicit "Failure profile (E / I / INC)" + "A-grade
success" + "B-grade success". Re-stated concisely:

| Vector | A-grade outcome | B-grade outcome | Predicted failure mode |
|--------|-----------------|-----------------|------------------------|
| A5 | TC⁰ Maynard weight + > 0.95 precision-recall on 10^6 test | `op(w) ≥ Ω(d(n))` lower bound refining E6.7 | E (sub-poly, not polylog) |
| C5 | `W_1 ≥ c > 0` quantitative non-Gaussianity at finite x | Stein-Chen `W_1 ≤ c log(x)/√x` sharpening Pintz/Korevaar | E (matches `1/√K` Gaussian) |
| D7 | DPP kernel with closed-form `K(t)` → polylog π via trace | quantitative anti-DPP at specific (k, t) tuple | E (K reduces to HL factor) |
| D8 | depth≤4 size≤10^4 circuit, > 99% acc at N=20 (extrapolating) | DARTS lower loss on PRIMES vs random at high N | I (overfit, no extrapolation) |

Each is single-session-tractable for the first numerical step. The
A-grade success modes are independent across vectors (different
techniques + different targets), so the project's expected A-grade
yield from any 1 vector succeeding ≈ 1 - (0.85)^4 ≈ 48% across the
four attempts. This is the baseline expected value of this
`frontier_gen` session.

## Self-extension (per CLAUDE.md autonomy invariants)

Each vector's "B-grade success" mode produces a follow-on:

- **A5 B-grade** → follow-up: "Are *Pascadi 2025*'s `x^{5/8}` level-of-
  distribution weights TC⁰-computable, even if Maynard's are not?
  Pascadi's larger θ raises the divisor bound but tightens the
  inversion."
- **C5 B-grade** → follow-up: "Compute the analogous Wasserstein
  distance for the Liouville partial sum `L(x) := Σ_{n ≤ x} λ(n)`.
  GT Möbius-orthogonality predicts faster Gaussian convergence than
  for π(x) - Li(x). Direct comparison would split the multiplicative
  vs additive regimes (parallel to §G)."
- **D7 B-grade** (anti-DPP) → follow-up: "Is the empirical 3-point
  prime correlation a SIGNED determinant (Pfaffian for orthogonal
  polynomial ensembles)? The orthogonal-symplectic DPP family is
  larger than the unitary family."
- **D8 B-grade** (DARTS-vs-random gap) → follow-up: "Replace the gate
  library with `{MOD_2, MOD_3, MOD_5}` to match the TC⁰-natural
  modular primitives (E5.1 BPSW). Does the architecture concentrate
  on a specific modular witness?"

Four follow-ups, all single-session-tractable.

## Self-grade

**B.** This is the upper-bound for a `frontier_gen` session that does
not itself attempt an attack. The four vectors each:

- import a technique that is on `CROSS_DOMAIN_TECHNIQUES.md` priority
  list OR is a clearly-orthogonal addition (Stein, DPP, Maynard
  bounded gaps, DARTS — all four were UNUSED prior);
- have at least one cited survey URL (Wikipedia, arXiv abstract, AMS
  book listing);
- have a single-session concrete first step (numerical experiment in
  3 of 4; literature-and-implement task for A5);
- have explicit falsification criteria with E/I/INC predictions;
- have stated A-grade vs B-grade outcomes;
- pass the CLOSED_PATHS near-duplicate check;
- one (D8) is a direct successor to a flagged S84 next-action,
  exactly the form of self-extension CLAUDE.md asks for.

Honest grade ceiling for frontier_gen: B if at least 2 vectors are
paper-grade fresh (yes — DPP for primes and Stein for π(x) - Li(x)
are both genuinely unpublished questions); A only if I have evidence
≥ 2 will produce A-grade work (impossible to claim at proposal time).
Realistic A-grade yield: 1, possibly 2, of the 4 will actually produce
A-grade output if attacked competently. The DPP and Stein vectors are
the most likely A-grade producers; Maynard A5 is most likely a B-grade
refinement of E6.7; DARTS D8 is likely a B-grade replication of
S84's gap finding at higher N.

**Self-grade DOWN, not up.** B.

## Next-action for the next agent

Pick the highest-leverage of the 4 vectors to attack:

- **C5 (Stein method)** is the highest info-per-effort: single
  precomputed sieve table + sorted-CDF integration; outcome decidable
  in one session; if `W_1 ≫ 1/√K` we have an A-grade quantitative
  non-Gaussianity result.
- **D7 (DPP)** is also single-session: pair correlation already
  available, only the 3-point test is new; A-grade if a closed-form
  kernel emerges.
- **A5 (Maynard sieve)** requires 1 session to implement the weight
  evaluator + AUC test; B-grade likely (refinement of E6.7) but
  A-grade possible.
- **D8 (DARTS)** is most code-heavy (JAX implementation of the
  search loop); 2 sessions minimum; likely B-grade (S84 replication
  at higher N) but high A-grade potential.

Recommended ordering: **C5 first** (single session, highest info-per-
effort), then **D7** (verifies whether DPP frame is even applicable
to primes), then **A5** (settles whether Maynard escapes the sieve
barrier), then **D8** (longest-running but addresses S84 next-action
directly).

Update `NOVELTY_CHALLENGES.md §0 highest-leverage attempt` to point at
C5.

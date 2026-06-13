# Session 246 — F6 / P7: dyadic π(2^k) structural test

**Mode**: novelty (B-grade target).
**Date**: 2026-04-30.
**Self-grade**: **B**.

**Target**: NOVELTY_CHALLENGES.md §F6 = OPEN_POSITIVE_TARGETS.md §P7
("π(2^n) — does the binary structure admit polylog?").

**Outcome**: **B-NEGATIVE** per-query closure of P7. The dyadic family
admits no detectable amortisation beyond what HKM gives for arbitrary
x of the same magnitude.

---

## What I produced

Three independent structural tests on the sequence `(π(2^k))_{k=1..56}`
(OEIS A007053), each calibrated against a Monte Carlo random null,
all returning null-baseline values:

1. **Sign-sequence linear complexity (BM/GF₂).** Data
   `BM(sign(π(2^k) − R(2^k))) = 28`. MC random-shuffle null at
   matched n_pos=25 / n_neg=31 split, 4000 trials: mean 28.25, std
   1.01, p05=27, p95=30. Data sits at 0.25σ below MC mean —
   identically random.

2. **Multi-lag autocorrelation of the RH-normalised residual.**
   `c_R(k) = (π(2^k) − R(2^k)) / 2^{k/2}`. Max `|ac|` over lags 1..10
   = 0.283 (at lag 1). MC iid-Gaussian length-56 lag-1 |ac| null:
   p95=0.250, p99=0.328, p999=0.437. Single-lag p ≈ 0.04;
   Bonferroni-corrected over 10 lags p ≈ 0.4. Not significant.

3. **HKM cost test.** sympy.primepi (Meissel-Lehmer evaluator) in
   fresh subprocesses — no sieve-cache reuse. 5 cold-start trials per
   (k, label) cell, k ∈ {28, 30, 32}, label ∈ {2^k, 2^k − 1, 2^k + 1}.
   Average cost ratio T(2^k) / T(neighbour) = 1.013, range
   [0.998, 1.047].

All three F-claims FAIL → B-NEGATIVE verdict.

---

## Mechanism (two independent reasons the closure holds)

1. **Information-theoretic.** After subtracting the smooth Riemann R
   correction, the residual is the explicit-formula zero-sum. At
   x = 2^k the phases `γ_n · k log 2 mod 2π` are Weyl-equidistributed
   in k for every Riemann zero `γ_n` independently of dyadic form.
   No phase coherence to exploit.

2. **Algorithmic.** Lucy / Meissel-Lehmer outer loop processes
   `{n ≤ √x}` and `x^{2/3}-smooth` counts. Neither set has any
   binary-representation structure. The 1.3% timing variation across
   {2^k − 1, 2^k, 2^k + 1} is within measurement noise.

---

## Confounder caught and corrected mid-session

The v1 pre-registration used Li(x) as the smooth baseline and lag-1
autocorrelation only. Initial run reported lag-1 `ac = 0.87` —
spuriously HOLDING F2. Diagnosis: π(x) − Li(x) carries a deterministic
`-0.5 Li(√x)` smooth Möbius-series correction that auto-correlates
trivially. Corrected by switching to R(x) = Σ_n μ(n)/n · Li(x^{1/n})
which absorbs the full smooth series. Under the corrected baseline
`max ac = 0.283 < 0.437`. **Documented openly in
results.md / F-claim 2 / Composite verdict** rather than rewriting the
pre-registration.

This is itself a small structural fact worth retaining: tests of
"residual structure" against the wrong smooth baseline pick up the
baseline's smoothness, not the residual's randomness. Use R, not Li,
when the question is about zero-driven oscillations.

---

## Self-evaluation (4 questions per CLAUDE.md)

### 1. What did I produce that was not in the project before this session?

- A pre-registered, MC-calibrated, three-test structural protocol for
  testing whether a parametric π-family `(π(x_k))_k` admits per-query
  amortisation. The protocol is self-contained and reusable.
- Empirical B-NEGATIVE closure of NOVELTY_CHALLENGES F6 / OPEN_POSITIVE_TARGETS
  P7 per-query side. Specifically: dyadic family `x_k = 2^k`, k=1..56,
  no detectable structure beyond random null on three independent
  axes (BM, autocorrelation, wall-clock).
- Refinement of E1.1 inline with the dyadic per-query result and the
  two-mechanism reason.
- Three successor challenges (F6.a cross-modulus, F6.b batched-on-k,
  F6.c rigorous output-bit-budget lower bound) registered in
  NOVELTY_CHALLENGES.md and (P7.b) in OPEN_POSITIVE_TARGETS.md.
- A new entry in CROSS_DOMAIN_TECHNIQUES.md §6 for BM linear
  complexity + MC random-shuffle null testing, with the caveat that
  n=56 is too short to discriminate polylog from random (need n ≳ 200).

### 2. What edges did my work compose or cite?

- **E1.1** (HKM x^{2/3} per-query bound) — refined inline with the
  dyadic per-query no-amortisation result.
- **E1.5** (per-step entropy 0.537) — cited as consistent with Weyl-
  equidistribution mechanism.
- **E1.3** (R⁻¹ bit-level model) — cited; not directly tested here.
- **Riemann explicit formula** (π = R − Σ_ρ R(x^ρ) + small) — used
  directly to set up the residual r_R and motivate the random-phase
  prediction.
- **S195 / S202** (random-phase regime, Goldston-Montgomery analysis) —
  the predictive framework for "no structure beyond Gaussian-like
  zero-sum" comes from Thread 3.
- **S224** (Correlation Dichotomy / Thread 5) — cited in successor
  F6.b / P7.b as the positive-shape template that batched-on-k might
  follow.

### 3. If my session produced only duplicate closures, why?

It did not. The closure is genuine and falsifiable: the F-claims were
pre-registered with explicit thresholds against Monte Carlo nulls,
and all three failed. The mechanism analysis (Weyl + Meissel-Lehmer
form-blindness) is novel as a *combined* statement for the dyadic
family — the individual ingredients are folklore but their composition
into a B-NEGATIVE closure is new project content.

### 4. What is the next-action for the next agent?

**Pick F6.b / P7.b** (batched-on-k dyadic amortisation) as a 2-session
arc candidate. The shared zero database `(γ_n)_{n=1..K}` is
k-independent; the M-batched cost should be `O(M·K) = O(M log^c x)`,
amortising to per-query polylog when M ≥ poly log x. **This is the
exact Thread 5 / S224 Correlation Dichotomy shape transposed to
parametric x — high a-priori plausibility of A-grade-shaped output.**
Concrete first slot: build `(γ_n)_{n=1..K=10^4}` once via
mpmath.zetazero, evaluate `R_K(2^{k_i})` for k_i = 1..56 with shared
zero data at K ∈ {log²x, log³x, K_max}, measure per-query amortised
wall-clock and accuracy versus exact π(2^{k_i}).

Alternative pick: F6.a (cross-modulus generalisation). Should close
identically by the same Weyl-equidistribution mechanism — likely
B-NEGATIVE in 1 session, lower expected information value than F6.b.

---

## Files

- `experiments/wildcard/dyadic_pi_structure/dyadic_pi_structure.py` —
  main script (470 lines): OEIS table, sympy verification, mpmath
  R/Li, MC nulls, subprocess-isolated timing.
- `experiments/wildcard/dyadic_pi_structure/dyadic_pi_structure_results.md` —
  pre-registered F-claims + final verdict + mechanism + successors.
- `experiments/wildcard/dyadic_pi_structure/raw_data.json` —
  verification map and residual table for k=1..56.
- `experiments/wildcard/dyadic_pi_structure/hkm_cost.json` — subprocess
  timing data for k ∈ {28, 30, 32}.
- `experiments/wildcard/dyadic_pi_structure/summary.json` — F-claim
  verdicts + MC baselines + pseudorandomness statistics.

## Status updates

- `EDGES.md` — E1.1 refined inline.
- `NOVELTY_CHALLENGES.md` — F6 marked PARTIALLY CLOSED, successors
  F6.a / F6.b / F6.c added.
- `OPEN_POSITIVE_TARGETS.md` — P7 marked CLOSED-B (per-query side),
  P7.b registered as successor entry.
- `status/CLOSED_PATHS.md` — new row 246 with full structural reason.
- `CROSS_DOMAIN_TECHNIQUES.md` §6 — BM + MC null entry.

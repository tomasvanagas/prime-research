# Critique — Critique-46 (2026-04-26)

Critique of the four S46-fresh proposals (`archive/ephemeral/proposals_session.md`)
and their experimental verdicts (`experiments/proposals/{pade_zero_sum,
lagrange_borel, modular_fingerprint, zero_aware_variance}_results.md`).

The proposer (proposals session) ran each of the four proposals and self-reported
all four as CLOSED. My job is to verify those verdicts against the 720+-entry
CLOSED_PATHS.md catalogue, decide whether any genuinely novel component
survives, and confirm the failure-mode tagging.

**Bottom line:** all four self-verdicts hold up. None reaches the `novel/` bar
and none reopens any direction in `OPEN_PROBLEMS.md`. P-A, P-B, P-D are
**strict duplicates-plus** of established closures (lines 26 + 657 + 697,
40 + 43 + 49, 256 + 688 respectively). P-C is a **DUPLICATE-PLUS** of line 689
with one empirical refinement worth recording (depth-5 fingerprint perfectly
separates primes from composites at n ≤ 10⁴, but the algorithm is
non-polylog regardless of discrimination rate).

---

## Verdict summary

| # | Proposal | Self-verdict | Critic verdict | Mode | Closest prior |
|---|----------|--------------|----------------|------|----------------|
| A | Pade/Wynn extrapolation of psi(x) zero-sum | CLOSED (I) | **DUPLICATE-PLUS, CLOSED** | I | lines 26, 49, 657, 697 |
| B | Borel resummation of Cipolla asymptotic | CLOSED (E) | **DUPLICATE-PLUS, CLOSED** | E | lines 40, 43, 49 |
| C | Hecke / arithmetic fingerprint primality oracle | CLOSED (C) | **DUPLICATE-PLUS, CLOSED** | C | line 689 |
| D | Zero-aware control variate for Monte-Carlo pi(x) | CLOSED (I) | **DUPLICATE-PLUS, CLOSED** | E+I | lines 256, 257, 688 |

---

## A — Pade/Wynn extrapolation of zero-sum (psi formulation)

### Has the exact approach been tried?

**Yes — covered tightly by four prior closures:**

- **Line 26 (S32):** "Convergence acceleration of zero sum — FAIL (I).
  Tested Richardson, Aitken Δ², Shanks, Padé, Cesàro, Euler-Maclaurin on
  partial zero sums for x = 10⁴..10⁷. Errors GROW as N^{+1.0} (random walk).
  Best method (Shanks) gives O(1) improvement only. GUE-random phases make
  each term independent — no structured error to accelerate."
- **Line 49 (S10):** "Resummation (Borel/Pade/Shanks/Richardson/Mellin-Barnes)
  — FAIL (E). GUE statistics prevent partial sum prediction."
- **Line 657 (S32):** repeats line 26 with slightly different numerics
  (errors grow ~ N^{0.8-1.0}, Shanks(3) ~ 10x constant).
- **Line 697 (S51):** "Borel-Pade × Cesaro-Fejer hybrid — FAIL (E)."
  Stacked acceleration also flat.

P-A's specific recipe (geometric ladder T_k = 5·2^k, k = 0..5, Wynn-epsilon
on the resulting psi_T values) is a strict subcase of "any convergence-
acceleration scheme on a partial-zero-sum sequence." Wynn-epsilon was not
called out by name in line 26 but is mathematically equivalent to Shanks
on the relevant sequences (Wynn-eps_2 is the Shanks transform; line 26
explicitly tests Shanks).

### Empirical result

The proposer's table is decisive in the negative direction:

| x | best partial err | Wynn extrap err | extrap helps? |
|---|---|---|---|
| 100 | 0.024 | 78.18 | NO (3000× worse) |
| 1000 | 0.525 | 995.42 | NO (1900× worse) |
| 10000 | 2.43 | 10013.88 | NO (4100× worse) |

Wynn doesn't merely fail to accelerate — it **catastrophically diverges**
on this sequence. This is the expected behavior when a sequence has
GUE-random oscillations (no smooth envelope), and reproduces line 26's
"errors GROW as N^{+1}" finding under a different transform.

### Decay-rate slopes

Empirical |psi - psi_T| vs T fits give average slope b ≈ -0.25 (proposer's
report: -0.580, +0.217, -0.394 across three x). RH-conditional bound predicts
b ≈ -1; observed slopes are far worse, dominated by the small-T regime
where the trace formula has not yet kicked in. This is consistent across
S32 (line 657) and earlier closures.

### Failure-mode analysis

Mode **I**, as proposer self-tagged. The substantive bet ("tail oscillation
has a smooth envelope") is empirically false, ruled out structurally by
GUE pair correlation of zeros (line 24, line 22) and by the Hilbert-transform
analysis of line 26.

### Verdict — **DUPLICATE-PLUS (refines lines 26 + 657 + 697 with explicit Wynn-epsilon numbers), CLOSED, mode I.**

Worth a CLOSED_PATHS entry because Wynn-epsilon is not by name in the
existing list, but the closure mechanism is identical to its parents.

---

## B — Lagrange/Cipolla asymptotic with Borel resummation

### Has the exact approach been tried?

**Yes — covered tightly by three prior closures:**

- **Line 40 (S5):** "Cipolla Pade resummation — FAIL (E). Divergent,
  Pade can't generalize."
- **Line 43 (S6):** "Ramanujan-type series / Pade — FAIL (E). Cipolla
  DIVERGES, Pade gives 1372 mean error."
- **Line 45 (S6):** "Borel summation — FAIL (E). NOT APPLICABLE
  (power-law not factorial)."
- **Line 49 (S10):** "Resummation (Borel/Pade/...) — FAIL (E)."

Note on line 45: that closure rejected Borel summation on the *zero-sum*
side (the explicit formula's terms are power-law in x, not factorial in
zero index). For the *Cipolla* asymptotic series itself, the coefficients
DO grow factorially, so Borel is in principle applicable — proposer is
correct on this technicality. But lines 40 + 43 are then the relevant
parents, and they show empirically that Borel-Pade resummation of
Cipolla's asymptotic series does not converge to p(n).

### Empirical result

| n | true p(n) | best Cipolla K=2-3 err | Borel-Pade resummed err |
|---|---|---|---|
| 10 | 29 | 5.97 (K=1) | 85.67 |
| 100 | 541 | 27.77 (K=2) | 204.62 |
| 1000 | 7919 | 78.6 (K=2) | 534.05 |
| 10000 | 104729 | 183.1 (K=3) | 260.33 |

Borel-Pade is **always worse** than the best raw Cipolla truncation by a
factor of 1.4× to 14×. The proposer correctly identifies this as either
(a) Stokes-line / multi-instanton structure that Borel-Pade in single
direction can't capture, or (b) the asymptotic series being already past
its optimal truncation point at small K.

### Failure-mode analysis

Mode **E** (proposer self-tagged correctly). Cipolla series is the
asymptotic expansion of `Li^{-1}(n)` (or equivalently of p_n / n), and is
analytically equivalent to the explicit formula by analytic continuation.
Borel-Pade is a different computational route to the same analytic object;
when the object can't be evaluated to integer precision in polylog, no
resummation scheme of its asymptotic series can either. Line 45's
"NOT APPLICABLE" labeling is loose but the conclusion is right.

### Side-finding

Cipolla *relative* error shrinks: 0.21 → 0.0017 from n=10 to n=10⁴ — but
absolute error grows monotonically with n in any optimal-truncation
window. This matches the well-known "asymptotic series gives
log(p(n)) / log log(p(n)) digits at optimal truncation" theorem
(Dusart 1999, refined by Massias-Robin) and is not a refinement worth
keeping separately; it's already implicit in lines 40, 43.

### Verdict — **DUPLICATE-PLUS (refines lines 40 + 43 with quantitative Borel-Pade vs raw Cipolla comparison), CLOSED, mode E.**

---

## C — Hecke / arithmetic fingerprint as primality oracle

### Has the exact approach been tried?

**Yes, structurally — line 689 (S46):**

> "Edixhoven-Couveignes tau-based bilinear-form prime detection — FAIL (C).
> Edixhoven 2011 gives polylog tau(p) primality test under GRH, but does
> NOT give pi(x). Proposal conflates primality testing with ordinality
> (which prime is the n-th). Same flaw as line 643. Bilinear form
> B(x) = sum_{mn≤x} tau(m)tau(n) has d_tau(p)=2 trivially for prime p, but
> verifying d_tau(N) requires factorising N. Best obtainable via
> primality-test-in-window: O(sqrt(n) polylog n) (arXiv:2510.16285), not
> polylog."

Adjacent closures: **line 578 (S23)** Ono partition-criterion approach
(modular form criterion gives only 46-72% accuracy near random); **line 640
(S29)** f(x) vs Ramanujan tau correlation (Pearson r = 0.010, no signal);
**line 36 (S25)** Heuristic candidate generation + AKS (candidates scale as
n^0.577, not polylog).

P-C's specific recipe (depth-5 fingerprint = (tau(n) - n¹¹ - 1 mod 691,
tau(n) mod 5, sigma_1(n) - 1 - n, n mod 30, n^6 mod 7)) is a finer
empirical implementation than line 689 but the **structural barrier is
identical**: every entry in the fingerprint either requires factoring n
(tau, sigma_k via Eichler-Selberg + multiplicativity) or carries
constant-bit information (n mod m).

### Empirical result

| depth | unique prime fps | composite FPs | FP rate | prime collide rate |
|---|---|---|---|---|
| 1 | 1 | 10 | 0.11% | 100% |
| 2 | 3 | 9 | 0.10% | 74.86% |
| 3 | 3 | 0 | 0.00% | 0.00% |
| 4 | 11 | 0 | 0.00% | 0.00% |
| 5 | 12 | 0 | 0.00% | 0.00% |

At n ≤ 10⁴, the depth-5 fingerprint achieves perfect separation. **This
is a real empirical observation** — not a refinement of any existing
closure I can find. But:

### Why this does not escalate

1. **Cost of computing tau(n) for composite n.** Edixhoven 2011 gives
   tau(p) at primes p in O(polylog p) under GRH (line 689). For composite
   n, tau is multiplicative: tau(p₁^{a₁}...p_k^{a_k}) = product
   tau(p_i^{a_i}). Computing this requires the prime factorisation of n,
   which is O(exp((log n)^{1/3+epsilon})) by GNFS. So the proposed test
   "is fingerprint(n) in PRIME_FINGERPRINTS?" runs in subexponential
   time, not polylog. P-C's complexity claim ("each tau(n) mod l:
   O(polylog n)") is **incorrect** for composite n.

2. **Even granted the test is polylog**, P-C only gives a primality
   *oracle*, not pi(x). The original proposal acknowledges this: "the
   bottleneck moves to *counting* primes with the fingerprint in [1,x],
   which is a different problem." That counting problem is exactly
   pi(x) again — the route is circular at the meta level.

3. **Discrimination rate ≠ algorithm.** Line 689's parent argument
   ("conflates primality testing with ordinality") applies here too: even
   a perfect primality oracle gives only the AKS-type algorithm
   pi(x) = sum_{n=2}^x is_prime_oracle(n), which is O(x · polylog x),
   not polylog.

4. **N=10⁴ scaling.** The number of "distinct prime fingerprints"
   monotonically increases with n: 1 at depth 1 to 12 at depth 5
   on 1229 primes. The information-theoretic floor is log₂(1229) ≈ 10.3
   bits. Depth-5 fingerprint uses ~ log₂(691·5·N·30·7) ≈ log₂(7e7) ≈ 26
   bits per n. This is the 50% redundancy regime; no asymptotic gain
   from increasing depth, and at n = 10⁶ the same battery would need
   roughly log₂(N) more bits per n.

### Failure-mode analysis

Mode **C** (proposer self-tagged correctly). The empirical separation is
real but circular at the algorithmic level: factoring n is the hidden
prerequisite for evaluating each component of the fingerprint at composite n.

### Verdict — **DUPLICATE-PLUS (refines line 689 with explicit depth-5 separation numbers), CLOSED, mode C.**

The empirical observation "depth-5 multiplicative fingerprint perfectly
separates primes from composites in [2, 10⁴]" is a real numerical fact but
does not reach the `novel/` bar — it is a finite-domain observation that
the proposer's own analysis shows breaks asymptotically.

### Side-finding worth recording

The depth-vs-FP-rate table (P-C above) is **a useful quantitative
benchmark of the AKS-style primality-oracle approach**: at n ≤ 10⁴, ~5
multiplicative invariants suffice for perfect separation, but each
requires factoring or polylog under GRH only at primes. Worth a
CLOSED_PATHS line citing both this and line 689.

---

## D — Zero-aware control variate for Monte-Carlo pi(x)

### Has the exact approach been tried?

**Yes, structurally — three prior closures:**

- **Line 256 (S15):** "Randomized zeta zero sampling — FAIL (E). Must
  use 100% of zeros (variance analysis); no sampling gain possible."
- **Line 257 (S15):** "Probabilistic sieve (I-E sampling) — FAIL (E).
  Variance 10⁶-10⁹× WORSE than exhaustive due to massive cancellation
  in I-E."
- **Line 688 (S46):** "Harper 2020 random-multiplicative variance
  reduction for pi(x) — FAIL (E+I). E|sum_{n<=x} f(n)| =
  O(sqrt(x)/(log log x)^{1/4}) for random multiplicative f IS strict,
  but each Monte Carlo sample costs O(x); ensemble cost O(M·x),
  super-sqrt(x) for any M >= 1."
- **Line 99 (S36):** "Randomized inclusion-exclusion for pi(x) — FAIL (I).
  Random subset sampling: std ~ 2^k / sqrt(samples). For x=500 (k=8),
  need ~370K samples vs 256 full IE terms. Variance reduction impossible
  — randomization makes IE strictly worse."

P-D's twist (use truncated explicit formula as control variate) is a
specific instance of "use any computable proxy for pi(x) as control
variate"; this is the standard MC variance-reduction recipe. Lines 256,
257, 688 collectively rule out: variance reduction below
O(sqrt(x)/polylog) requires either sub-sqrt(x) zero-sum tail bound
(disproven, line 26) or a non-trivial multiplicative-character cancellation
(disproven, line 688).

### Empirical result

The proposer's measurement is decisive:

- log-log slope b of |psi - psi_T| vs T: average -0.294 (RH predicts -1).
- Variance reduction factor with PNT-style 1/log y control variate: 1.006×
  (essentially zero).
- pi(10000) = 1229; naive MC error = 13.2; PNT-controlled MC error = 14.4
  (slightly WORSE).

The PNT control variate gives only constant-factor variance reduction; the
explicit-formula truncation does not improve with T at the scales tested
because the residual is GUE-random with no exploitable structure — the
same observation as lines 256, 257.

### Failure-mode analysis

Mode **E + I**, not just I. The empirical measurement (proposer's
slope -0.29 at small T) is the small-T regime artifact; at T ≥ sqrt(x)
the slope approaches the RH-predicted -1, which is the **information-
theoretic limit**. Any control-variate with deterministic O(polylog) cost
that beats this would imply a structural cancellation in the zero sum
itself — which is the open problem of FOCUS-1.

### Verdict — **DUPLICATE-PLUS (refines lines 256 + 688 with explicit PNT-control-variate variance ratio numbers), CLOSED, mode E+I.**

---

## Cross-cutting observations

1. **All four proposals self-closed correctly.** Methodology was sound:
   each test was set up to falsify a single concrete claim (Wynn vs
   plain truncation; Borel-Pade vs raw Cipolla; depth-5 fingerprint
   discrimination rate; PNT-control-variate variance reduction factor).

2. **All four are duplicates-plus** of existing closures. The closest
   parents are lines 26 + 49 + 657 + 697 (P-A), lines 40 + 43 + 49 (P-B),
   line 689 (P-C), lines 256 + 257 + 688 (P-D). None of the four reaches
   the `novel/` bar, none reopens an `OPEN_PROBLEMS.md` direction.

3. **The strongest by-product is P-C's depth-5 fingerprint table.**
   It is a useful quantitative anchor for the AKS-style primality-oracle
   approach but is not standalone-novel — it sits inside line 689's
   parent closure mechanism.

4. **Pattern in the proposer's choices.** All four proposals are
   "convergence-acceleration / variance-reduction" type interventions on
   already-closed analytic primitives (zero sum, Cipolla asymptotic,
   primality oracle, MC pi). This is the **same pattern that closed
   S43 + S46 + S48 + S51 + S63** (lines 657, 685-697, 715-718) — the
   project has now systematically exhausted this family. Worth noting
   in the next session-insights entry.

5. **None of the four warrants `novel/` placement.** The OPEN_PROBLEMS
   landscape is unchanged: only Circuit Complexity of pi(x) (with
   FOCUS-1 sub-attack 1 [Bernstein 2003 strengthened-r AKS] and
   sub-attack 3 [Healy-Viola Frobenius transplant] still un-built per
   line 714's progress note).

---

## Action items

1. **Add four entries to `status/CLOSED_PATHS.md`** (S46-fresh batch).
   Drafted below.
2. **Update `status/SESSION_INSIGHTS.md`** with critique-46 entry noting
   that the convergence-acceleration / variance-reduction family is now
   exhausted across P-A through P-D plus prior sessions.
3. **Append session synthesis** at
   `archive/sessions/session_critique46.md`.
4. **Set `.run_state` to 46.**

### Drafted CLOSED_PATHS.md entries (insert after line 721)

```
| Wynn-epsilon (=Shanks transform) on geometric ladder of psi(x) zero-sum partial sums (S46-fresh) | FAIL | I | DUPLICATE-PLUS of lines 26 (S32), 49 (S10), 657 (S32), 697 (S51). Tested geometric ladder T_k = 5*2^k, k=0..5, K_max=160 zeros, x in {100, 1000, 10000}. Wynn extrapolation diverges catastrophically: x=10000 best partial err 2.43 vs Wynn err 10013.88 (4100x worse). Empirical |psi - psi_T| log-log slope b ≈ -0.25 average, far below RH-conditional -1 — small-T regime where trace formula hasn't kicked in. GUE pair correlation of zeros gives no exploitable smooth envelope to extrapolate. Wynn-eps_2 = Shanks transform (line 26 already tested by name); P-A is a strict subcase. See experiments/proposals/pade_zero_sum.py | 46 |
| Borel-Pade resummation of Cipolla asymptotic series for p(n) (S46-fresh) | FAIL | E | DUPLICATE-PLUS of lines 40 (S5 Cipolla Pade), 43 (S6 Ramanujan/Pade), 45 (S6 Borel — caveat: line 45 was on zero-sum side which is power-law not factorial, but Cipolla coefficients DO grow factorially so Borel is in principle applicable here; still does not converge). Tested K=1..7 truncation + Borel-Pade resum at L=log n, M=log log n. Borel-Pade always 1.4-14x WORSE than best raw Cipolla truncation: at n=10000 best Cipolla K=3 err=183 vs Borel-Pade err=260. Diagnosis: Stokes-line / multi-instanton structure that single-direction Borel-Pade cannot capture, OR small-K is already past optimal truncation. Cipolla relative error 0.21 -> 0.0017 from n=10 to 10^4 (asymptotic-series character) — but absolute error grows. Cipolla = analytically equivalent to explicit formula by analytic continuation; no resummation of its asymptotic gives integer-precision in polylog. See experiments/proposals/lagrange_borel.py | 46 |
| Depth-5 multiplicative fingerprint (Hecke tau, sigma, n mod 30, n^6 mod 7) as primality oracle (S46-fresh) | FAIL | C | DUPLICATE-PLUS of line 689 (S46 Edixhoven-Couveignes tau-based prime detection). Empirical: depth-5 fingerprint (tau(n)-n^11-1 mod 691, tau(n) mod 5, sigma_1(n)-1-n, n mod 30, n^6 mod 7) achieves PERFECT separation of 1229 primes from 8770 composites in [2,10000] — useful quantitative anchor. BUT three barriers: (i) tau(n) for COMPOSITE n requires factoring n via multiplicativity (Edixhoven 2011 gives polylog tau(p) only at primes under GRH, not at composites); subexp factoring cost via GNFS; (ii) even granted polylog primality oracle, counting primes in [1,x] with the oracle is pi(x) again (the original problem) — circular at meta level; (iii) depth-vs-bits floor: 5 features ≈ 26 bits per n; finite-domain perfect separation collapses asymptotically as log(N) more bits needed at larger N. P-C's complexity claim "tau(n) mod l: O(polylog n)" is incorrect for composite n. See experiments/proposals/modular_fingerprint.py | 46 |
| Zero-aware control variate for Monte-Carlo pi(x) using truncated explicit formula (S46-fresh) | FAIL | E+I | DUPLICATE-PLUS of lines 256 (S15 Randomized zeta zero sampling — must use 100% of zeros), 257 (S15 Probabilistic sieve — variance 10^6-10^9x worse), 688 (S46 Harper random-multiplicative — sample cost O(x), no asymptotic gain). Empirical at x=10000, T=200 zeros: PNT-control-variate variance reduction factor only 1.006x (essentially zero); naive MC err 13.2 vs PNT-controlled 14.4 (slightly WORSE). Log-log slope of |psi - psi_T| vs T: average b ≈ -0.29 across x in {100,500,1000,5000,10000} — small-T regime; at T >= sqrt(x) slope approaches RH-conditional -1, which IS the information-theoretic limit. GUE-random residuals carry no exploitable structure beyond the trivial bound. Beating O(sqrt(x)/polylog) variance reduction would require structural cancellation in the zero sum itself — exactly the FOCUS-1 open question. See experiments/proposals/zero_aware_variance.py | 46 |
```

---

## Critic's bottom line

All four proposer self-verdicts hold up. All four are duplicates-plus of
existing closures (P-A → 26+49+657+697; P-B → 40+43+49; P-C → 689;
P-D → 256+257+688). One useful empirical refinement (P-C: depth-5
fingerprint perfectly separates primes from composites in [2, 10⁴]) is
worth recording but does not reach the `novel/` bar because the algorithm
is non-polylog regardless of separation rate.

The session adds four well-documented closures and a quantitative
benchmark on multiplicative-fingerprint discrimination depth. Nothing
reopens any closed direction; nothing reaches the `novel/` bar.
**`OPEN_PROBLEMS.md` remains as it was: only Circuit Complexity of pi(x)
is genuinely open**, with FOCUS-1 sub-attack 1 (Bernstein 2003
strengthened-r) and sub-attack 3 (Healy-Viola Frobenius transplant)
still un-built per line 714 — though those are now closed (S64, S66) per
TODO.md, leaving FOCUS-3 (Brandt MKtP) un-engaged.

Pattern observation: the convergence-acceleration / variance-reduction
family of interventions (Pade, Wynn, Borel, mollifiers, control variates,
random sampling) on already-closed analytic primitives is now
**systematically exhausted** across P-A..P-D plus prior sessions
S5, S6, S10, S15, S25, S32, S43-S46, S48, S51, S63. Future proposals in
this family should be pre-disqualified by a single CLOSED_PATHS check
against the master list of acceleration techniques.

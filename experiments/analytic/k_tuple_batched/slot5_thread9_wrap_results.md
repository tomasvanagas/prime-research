# Thread 9 (P4) slot 5 (FINAL) — structural decomposition + Thread 9 wrap

**Session:** S433
**Date:** 2026-04-30
**Code:** `slot5_thread9_wrap.py`
**Outputs:**
  - `slot5_extended_t.csv`       28 cells: T(w), Δ(w), log w / log log w, alpha_emp, coef_log_loglog
  - `slot5_decomposition.csv`    full / multi-K decomposition rows
  - `slot5_gm_residual.csv`      GM-zero residual analysis cells
  - `slot5_run.log`
**Self-grade:** **B** — slot 5 (FINAL) produces (a) extension of slot 4's
T(w)/Δ(w) numerics from K=200K → K=350K (1.75×); (b) **REFUTATION of
slot 4's "−12 C₂²" identification as a *unique fit*** by demonstrating
that Model B (Δ/w = A·log w·log log w + B·log w + C) is **fit-range
unstable** when extended to slot 5's range — A varies −5.53 (K_min=50K)
to −8.88 (K_min=200K), so slot 4's reported A = −5.36 (matching −12 C₂²
to 2.4%) was a 16-cell local minimum, not the asymptotic; (c) full
single-prime decomposition S_1(K=100K) = −2.04 (25.3% of empirical
S(K) = −8.03), partial cross-prime S_2(K=100K, primes ≤ 200) = −1.30
(16.2%), remainder 58.4% from cross-prime tail (primes > 200) +
higher-order. Not A because the closed-form derivation goal is
*unreachable* in this slot, AND the slot 4 structural identification
is empirically *unjustified* on extended data. Not C because the
slot 4 conjecture is genuinely **REFUTED** as a unique structural
identification at the 35% confidence level, and the partial-positive
characterization of Thread 9 is now precisely calibrated.

## Mission (slot 5 FINAL, from `.commit_state` recommended_next_action)

> "Slot 5 of Thread 9 (FINAL slot). Highest-yield: rigorous derivation
> of the −12 C₂² log·log log coefficient via Goldston-Yıldırım 2007
> partial-sum machinery on HL 4-tuple singular series. Alternative slot
> 5: (b) Goldston-Montgomery zero-residual analysis of F_pred − F_emp
> ≈ ±0.003; (c) wrap synthesizing the 5-slot Thread 9 arc as partial-
> positive Correlation-Dichotomy-shape result."

**Outcome:** Option (a) is **unreachable** — the asymptotic Δ(w) ∼ −12
C₂² · w · log w · log log w fails empirically when slot 4's fit range
is extended (slot 4's −5.36 was a finite-K artefact of insufficient
range, just as slot 4 had refuted slot 3's −5.72 at smaller range).
Option (b) is **undertested** — slot 5 only has slot 4's single
corrected F_pred cell. Option (c) is the **realised work**: Thread 9
wrap as partial-positive with explicit falsification of the slot 4
structural identification.

## Slot context (Thread 9 progression)

- **Slot 1 (S429, B):** 5-8× sieve-shared batched-x speedup at M=64;
  HL approximation HL_H(x; w) = C(H)·w/log^k x matches empirical
  ≤ 0.34 σ_Pois at 6 cells.
- **Slot 2 (S430, B):** F_HL_kt = σ_eff/σ_pois = 0.87 ± 0.03 measured
  for k=2 across 3 decades, 18 cells (σ_eff/σ_pois decade-stable to
  ±0.045 wide).
- **Slot 3 (S431, B):** prime-by-prime cancellation identity verified
  to bit precision (max dev 4×10⁻¹⁶, 168 primes); Theorem 1
  (windowed-twin Gallagher-Poisson under HL); F_pred match F_emp
  <0.5% wide regime; empirical fit Δ ≈ −5.72·w·log(w) + 24 C₂²·w.
- **Slot 4 (S432, B):** identified slot 3 software bug (off-by-one
  in singular-series tail); REFUTED slot 3's α = 5.72 as
  finite-w artefact (rolling-band α grows 6.43 → 9.73); Model B fit
  on K up to 200K gave Δ/w = −5.36·log w·log log w + 9.30·log w − 22.37
  (matched −12 C₂² = −5.230 to 2.4%); structural candidate
  Δ ∼ −12 C₂²·w·log w·log log w.
- **Slot 5 (THIS, S433, B):** EXTENSION + DECOMPOSITION + WRAP.

## What slot 5 produced (NEW content)

### F16. Slot 4's −12 C₂² structural identification REFUTED on extended data

**Setup.** Slot 4 fitted Model B (Δ/w = A·log w·log log w + B·log w
+ C) over K = 100..200,000 (22 cells), reporting A = −5.36, B = 9.30,
C = −22.37 with RMS rel err 1.5%. Slot 4 then identified A ≈ −12 C₂²
= −5.230 (2.4% rel diff) as the structural asymptotic candidate,
suggesting the conjecture Δ(w) ∼ −12 C₂²·w·log w·log log w as w → ∞.

**Slot 5 extension** computed T(w), Δ(w) at K = 250K, 300K, 350K
(extending log w from 14.0 to 14.6, an extra ≈0.6 in log w). The
extended cells are:

| K | log w | Δ(w) | Δ/w | −Δ/(w·log w) | −Δ/(w·log w·log log w) |
|---|------:|-----:|----:|------:|------:|
| 200,000 | 13.998 | −1.0961e8 | −91.34 | 6.526 | 2.473 |
| 250,000 | 14.221 | −1.4064e8 | −93.76 | 6.593 | 2.484 |
| 300,000 | 14.403 | −1.7234e8 | −95.74 | 6.647 | 2.492 |
| 350,000 | 14.558 | −2.0465e8 | −97.45 | 6.694 | 2.500 |

**Slot 5 Model B fits at multiple K_min:**

| K_min | n cells | A | B | C | RMS | \|A − (−12 C₂²)\|/\|...\| |
|------:|--------:|------:|------:|------:|-----:|------:|
| 1,000 | 19 | −7.252 | +15.55 | −41.06 | 0.068 | **38.7%** |
| 10,000 | 13 | −7.174 | +15.27 | −40.04 | 0.037 | **37.2%** |
| 50,000 | 8 | −5.527 | +9.31 | −17.51 | 0.008 | 5.7% |
| 100,000 | 6 | −6.184 | +11.70 | −26.65 | 0.006 | 18.2% |
| 200,000 | 4 | −8.884 | +21.60 | −65.48 | 0.005 | **70.0%** |

**Conclusion.** Model B coefficient A is **NOT stable** under K_min
variation when fit range is extended to K = 350,000. A varies from
−5.53 (K_min=50K, 5.7% from −12 C₂²) to −8.88 (K_min=200K, 70% off).
Slot 4's reported A = −5.36 was based on K_min=1000 fit on K up to
200,000 — but my K_min=1000 fit on K up to 350,000 gives A = −7.25,
35% off −12 C₂². Adding 3 high-K cells shifts the fit by 35%.

This means: **Model B is NOT a unique asymptotic identification.**
The 3-parameter fit y = A·log w·log log w + B·log w + C is highly
degenerate (log w and log w·log log w differ by only the slowly-
varying log log w factor, leading to fit collinearity over the
empirical log w ∈ [6, 15] range).

Slot 4's structural candidate "Δ(w) ∼ −12 C₂²·w·log w·log log w" is
**REFUTED as a unique closed-form identification**.

**The empirical residual −Δ/(w·log w·log log w) does grow toward −5.23
direction**: 2.24 (K=100) → 2.40 (K=10K) → 2.50 (K=350K). But the rate
is so slow that at K = 350,000 we are still 52% below −12 C₂². To
reach within 5% of −5.23, extrapolation suggests w ≳ 10⁵⁰ — far
outside any feasible numerical range.

### F17. Single-prime + cross-prime decomposition at K = 100,000

**Setup.** Use the Fourier expansion of S_4(m)/⟨S_4⟩ = ∏_{p≥5}(1 +
ε_p(m mod p)), where ε_p has zero mean over m mod p. Expanding the
product:

```
G(m) − 1 = ∑_p ε_p(m) + ∑_{p<q} ε_p(m) ε_q(m) + ...
```

Sum over k = 1..K:

```
S(K) = sum_p T_p(K)  +  sum_{p<q} T_{pq}(K)  +  higher
```

T_p(K) = ∑_{k=1}^K ε_p(6k mod p) is the **single-prime partial-cycle
sum** at prime p.

T_{pq}(K) = ∑_{k=1}^K ε_p(6k mod p) · ε_q(6k mod q) is the **cross-
prime second-order** at primes p, q.

**Single-prime S_1(K):** computed exactly for ALL primes p ∈ [5, 6K+10]
(splitting into p ≤ K via partial-cycle iteration, and K < p ≤ cap via
analytical n-spike formula). Values:

```
K = 10,000:  S_1 = -1.920
K = 30,000:  S_1 = -1.927
K = 50,000:  S_1 = -2.134
K = 80,000:  S_1 = -2.300
K = 100,000: S_1 = -2.035
```

S_1 is bounded around −2.0 and grows slowly (consistent with single-
prime ≈ −∑_{p≤K} 1/(p − 2) ≈ −log log K + Mertens-correction =
−log log 100,000 − 0.13 = −2.51 − 0.13 = −2.38; the empirical
−2.035 is slightly less negative, suggesting partial-cycle leftover
fluctuations and a partial cancellation against expectations).

**Cross-prime S_2(K, primes ≤ 200):** computed exactly by direct
iteration over the partial cycle of length K mod (pq) for each pair
(p, q). 946 pairs × average K mod pq ≈ 19000 = ≈18M iterations,
0.4 seconds.

```
K = 10,000:  S_2(p,q≤200) = -2.415
K = 30,000:  S_2 = -1.275
K = 50,000:  S_2 = -2.058
K = 80,000:  S_2 = -2.288
K = 100,000: S_2 = -1.303
```

S_2 fluctuates around −1.5..−2.4 with K (variance from specific K mod
pq distributions across all 946 (p, q) pairs).

**Empirical S(K) and decomposition:**

| K | S_emp | S_1 | S_2(p,q≤200) | remainder | S_1% | S_2% | rem% |
|---|------:|----:|----:|----:|----:|----:|----:|
| 10,000 | −6.85 | −1.92 | −2.41 | −2.51 | 28.0% | 35.3% | 36.7% |
| 30,000 | −6.66 | −1.93 | −1.28 | −3.46 | 28.9% | 19.2% | 51.9% |
| 50,000 | −8.01 | −2.13 | −2.06 | −3.82 | 26.6% | 25.7% | 47.7% |
| 80,000 | −9.10 | −2.30 | −2.29 | −4.51 | 25.3% | 25.1% | 49.6% |
| 100,000 | −8.03 | −2.04 | −1.30 | −4.70 | 25.3% | 16.2% | 58.4% |

**Conclusion.**

- **Single-prime captures 25–29%** (full S_1 over all primes p ≤ K
  plus boundary tail). This refines slot 4's reported "32%" (slot 4
  used a heuristic Mertens estimate, not exact computation).
- **Cross-prime partial (primes ≤ 200) captures 16–35%** with
  fluctuations.
- **Remainder = 37–58%** from cross-prime tail (primes > 200) plus
  3rd-order and higher correlations.

The decomposition does NOT close the gap. The cross-prime *tail*
(primes p > 200, q > 200) and higher-order terms together account
for half of S(K). A complete derivation would require:

1. The cross-prime second-order sum extended to all primes (computa-
   tional cost ≈ 10⁹ iterations at K=100K, manageable but not done in
   slot 5).
2. Third-order ∑_{p<q<r} T_{pqr}(K) to bound higher-order contribu-
   tions.

Neither is in slot 5's scope. The decomposition empirically validates
that the leading-order is NOT pure single-prime, but does not produce
a closed form.

### F18. Goldston-Montgomery zero residual analysis — UNDERTESTED

Slot 3 / slot 4 found the systematic residual F_pred − F_emp ≈ ±0.003
at the cleanest cell (10⁸ wide, w=4071). Slot 5's GM-zero hypothesis
test:

| x_anchor | w | F_emp | F_pred (slot 4 corr.) | F_diff | log w / log x |
|------:|----:|----:|----:|----:|----:|
| 10⁸ | 4071 | 0.9113 | 0.9080 | +0.0033 | 0.451 |
| 10⁷ | 3117 | 0.8774 | (not corrected) | (n/a) | 0.499 |
| 10⁶ | 2290 | 0.8215 | (not corrected) | (n/a) | 0.560 |

Slot 5 has only ONE clean cell (10⁸ wide). The hypothesis "F_pred −
F_emp scales as log w / log x" cannot be tested without re-running
slot 4's F_pred at 10⁶ and 10⁷ wide (slot 4 only corrected the 10⁸
wide cell explicitly). **Slot 5 declines this option** as undertested.

### F19. Thread 9 PARTIAL-POSITIVE characterization

Thread 9 closes with the following calibrated partial-positive
characterization (analogous to Thread 5 S224 Correlation Dichotomy,
Thread 7 S260 Approximate-π named ε, Thread 8 P2 Conditional batched
π_h):

**Achievement (positive):**
- (a) Slot 1 measurable speedup: **5–8× per-x throughput improvement**
  for batched correlated narrow-window k-tuple counts at M = 64,
  growing with √x · log log x as predicted (slot 1 magnitude smaller
  than Thread 5's 33× because per-query baseline is segmented sieve
  not Lucy DP).
- (b) Slot 2 measured constant: **F_HL_kt = 0.87 ± 0.03 across 3
  decades** for k=2 windowed, decade-stable to ±0.045 (5× tighter
  than Thread 8's F_HL ∈ [0.36, 0.70]).
- (c) Slot 3 conditional theorem: **Theorem 1 (windowed-twin
  Gallagher-Poisson under HL):** Var[N₂(x; w)] / E[N₂] → 1 as
  x → ∞ with w/x → 0. F_pred matches F_emp to <0.5% at the
  cleanest cell (10⁸ wide, slot 3 narrative preserved under slot 4
  bug fix).

**Limitations (open / negative):**
- (d) Slot 4 structural candidate −12 C₂² · log w · log log w
  REFUTED at the 35% level by slot 5 fit-instability analysis. The
  closed-form structural identification of Δ(w)'s asymptotic is
  **OPEN**.
- (e) Single-prime + small-cross-prime decomposition captures only
  41% of S(K). Cross-prime tail (primes > 200) plus higher-order
  contributions account for the remainder. Goldston-Yıldırım 2007
  partial-sum machinery would be needed to close, beyond slot 5
  scope.
- (f) Goldston-Montgomery zero residual analysis (slot 4's slot-5
  alternative target) UNDERTESTED with only 1 corrected F_pred cell.

**Conditional theorem statement (Thread 9 wrap):**

```
Conditional on:
  (HL-4) Hardy-Littlewood 4-tuple conjecture for N_2(x; w),
  (CD)   Slot 5 conjecture: Δ(w) = T(w) - 2 C_2² w² has asymptotic
         Δ(w) = -A_∞(w) · w · log w · log log w · (1 + o(1))
         with A_∞(w) → A_* ∈ [5.0, 9.0] as w → ∞ (interval REFINED
         from slot 4's point estimate -12 C_2² = 5.23 to slot 5's
         empirical 5-9 range based on Model B fit instability),

For x → ∞ at any fixed w(x) = polylog(x) (e.g. w = log² x · const):
  F²(x; w) := Var[N_2(x; w)] / E[N_2(x; w)]
            = 1 + Δ(w) / (C_2 · w · log² x) + ε_GM(x, w)
where ε_GM is the unaccounted Goldston-Montgomery zero contribution
(slot 3/4 systematic ±0.003 residual).

For w/x → 0 with x fixed and increasing: F²(x; w) → 1.
```

**Falsifiability (Thread 9 wrap):**

- F1: empirical α coefficient should converge to a definite value
  in [5.0, 9.0] as K → ∞. Slot 5: K = 350K shows coefficient
  −Δ/(w·log w·log log w) = 2.50, growing slowly toward unknown
  asymptote. Verification at K = 10⁹ would test convergence.
- F2: A complete cross-prime decomposition (all primes) should
  account for >90% of S(K). Slot 5: only 41% accounted with primes
  ≤ 200. Computing primes ≤ 10K needs ~10⁹ iterations.
- F3: GM zero contribution should produce F_pred − F_emp ≈ const ·
  log w / log x. Slot 5: only 1 cell tested.

## Why B-grade (CLAUDE.md grading)

**Not A** because:
- Closed-form derivation (slot 5 mission target option (a)) is
  unreachable. The Goldston-Yıldırım 2007 machinery would require
  Selberg-Delange-style partial-sum analysis of the HL 4-tuple
  singular series, beyond slot 5 scope.
- The asymptotic Δ(w) ∼ −12 C₂²·w·log w·log log w is **REFUTED as a
  unique fit** but not REPLACED by a different definitive structural
  shape. We have a refined empirical interval [5.0, 9.0] for the
  coefficient, but no closed-form.
- The decomposition (single-prime + small-prime cross-prime) accounts
  for 41% of S(K), not the >90% needed for a closure.
- Goldston-Montgomery zero analysis is undertested.

**Not C** because:
- Slot 4's structural identification with −12 C₂² is now **REFUTED at
  the 35% level** as a unique fit. This is a substantive negative
  result that revises the framework's understanding of Thread 9's
  asymptotic.
- The single-prime decomposition (S_1 = full single-prime sum,
  computed analytically for K < p ≤ cap) is the *first exact*
  computation in the project (slot 4 had only Mertens-style estimate).
- Cross-prime second-order contribution explicitly computed (16% at
  K=100K with primes ≤ 200) — first quantitative cross-prime measure-
  ment in the project.
- Thread 9 partial-positive characterization is precisely calibrated
  in slot 5 in a way that allows future framework agents to know the
  exact scope (positive: slots 1+2+3 algorithmic and theorem-level;
  open: closed-form Δ(w) asymptotic).

This is the project's **first multi-slot demonstration that empirical
fit to a structural constant (−12 C₂² in slot 4) was finite-K-fitted
artefact**, paralleling slot 4's REFUTATION of slot 3's −5.72 fit.
The pattern (each slot REFUTES the previous slot's structural
conjecture at a higher K range) suggests the asymptotic shape is
genuinely complex and may require analytic-NT techniques beyond the
project's current toolkit.

## Edges composed / cited

- **E1.5** (information density of π) — σ_pois prediction for
  Poisson baseline.
- **HL 1923** (k-tuple conjecture) — provides S_4.
- **Hardy-Littlewood prime quadruplet constant** S_4(0,2,6,8) =
  4.15118 — slot 4 bug detection reference.
- **Mertens 1874** (∑_{p ≤ x} 1/p ∼ log log x) — single-prime
  heuristic basis (slot 4 + slot 5).
- **Gallagher 1976** ("Distribution of primes in short intervals")
  — slot 3 Theorem 1 (windowed-twin Gallagher-Poisson under HL).
- **Goldston-Yıldırım 2007** — would close the cross-prime tail
  if applied; not done in slot 5.
- **Goldston-Montgomery 1987** — slot 5 GM zero-residual hypothesis
  (undertested).
- **Bombieri-Davenport 1966** — slot 3 cancellation identity context.
- **S429 / S430 / S431 / S432** — Thread 9 prior slots.

## Cross-domain ingredient (status update)

`CROSS_DOMAIN_TECHNIQUES.md`:

- **Pair-correlation analysis applied to HL k-tuple singular series**
  (slot 3 entry, mode I). Slot 5 result: **REFINES** to mode I-with-
  caveat — empirical fit (slot 4 Model B) is NOT structurally unique
  on extended data; slot 5 falsifies the −12 C₂² identification.
- **Selberg-Delange partial-sum machinery** (slot 4 PROPOSED) status
  remains PROPOSED. Slot 5 does NOT apply this technique (unreachable).
- **Fourier expansion of singular-series local factors**: slot 5 used
  this technique to compute S_2 (cross-prime second-order) directly.
  ε_p Fourier modes c_p(ξ) = O(1/p²) bounded; this is a **NEW
  technique-USE in the project** (mode I — implementation match).

## Falsifiability (slot 5 hypothesis status)

- **H15 (slot 4 −12 C₂² as unique structural fit):** **REFUTED** by
  Model B fit instability under K_min variation on extended K = 350K
  data. Coefficient varies −5.5 to −8.9 across reasonable fit ranges.
- **H16 (single-prime captures > 50% of S(K)):** **REFUTED** by full
  S_1 computation: captures 25–29% at K = 10K..100K.
- **H17 (cross-prime over primes ≤ 200 captures > 50%):** **REFUTED**
  by direct computation: captures 16–35%.
- **H18 (combined S_1 + S_2(p,q≤200) captures > 80%):** **REFUTED**
  by computation: captures 41–61%, leaving ≥ 39% unexplained.
- **H19 (GM zero residual ε_GM(x, w) ~ const · log w / log x):**
  **UNDERTESTED** — slot 5 has only 1 corrected F_pred cell.
- **H20 (Thread 9 partial-positive characterization at slots 1+2+3
  achievement):** **CONFIRMED** — slot 5 wrap explicitly delineates
  the positive (slots 1, 2, 3 algorithmic, empirical, conditional-
  theorem-level) from the open (slot 4 / slot 5 sub-leading
  asymptotic).

## Self-extension (per CLAUDE.md autonomy invariants)

This session BUILT a slot-5 wrap-mode synthesis. It REFUTES slot 4's
structural identification (substantive negative result), and produces
the Thread 9 partial-positive characterization. No new
ATTACK_VECTORS entry needed (Thread 9 itself was the attack; slot 5
closes it as DONE_PARTIAL_POSITIVE_CONDITIONAL).

`OPEN_POSITIVE_TARGETS.md` §P4 status: ACTIVE → CLOSED-CONDITIONAL.

**Successor entries** in `OPEN_POSITIVE_TARGETS.md`:
- **P4-extension-a (PROPOSED):** Apply Goldston-Yıldırım 2007 partial-
  sum machinery rigorously to the HL 4-tuple singular series Cesaro
  sum ∑_m S_4(0,2,m,m+2) · (w − m). If derived, would resolve slot 5's
  open structural question (Δ asymptotic shape).
- **P4-extension-b (PROPOSED):** GM zero-residual analysis at slot 3
  cells (re-run F_pred at 10⁶ wide, 10⁷ wide with corrected S_4) to
  test the slot 4 / slot 5 GM hypothesis.

## Self-evaluation (CLAUDE.md 4-question)

1. **What did I produce that was not in the project before this
   session?**
   - REFUTATION of slot 4's −12 C₂² identification as a unique fit:
     Model B coefficient varies −5.5 to −8.9 on extended data.
   - Full single-prime decomposition S_1(K) computed exactly via
     analytic n-spike formula for K < p ≤ cap. First exact computation
     in the project.
   - Cross-prime second-order S_2(K, primes ≤ 200) computed by direct
     Fourier-style iteration. First explicit cross-prime contribution
     measurement in the project.
   - Quantitative decomposition table at K = 10K..100K: S_1 25-29%,
     S_2 16-35%, remainder 37-58%.
   - Thread 9 partial-positive characterization with conditional
     theorem statement and 6 explicit falsifiers.

2. **What edges did my work compose or cite?** E1.5, HL 1923, HL
   prime quadruplet constant, Mertens 1874, Gallagher 1976,
   Goldston-Yıldırım 2007, Goldston-Montgomery 1987, Bombieri-
   Davenport 1966, S429-S432.

3. **If my session produced only duplicate closures, why?** It did
   not. Slot 5 produced (i) a substantive REFUTATION of slot 4's
   structural identification (extending fit range exposes
   instability), (ii) the first exact decomposition of S(K) into
   single-prime + cross-prime + remainder at K=100K, and (iii) the
   Thread 9 partial-positive wrap characterization.

4. **What is the next-action for the next agent?** Thread 9 is DONE
   (5 slots used). `.commit_state` should advance to next thread.
   Per CLAUDE.md, after Thread 5 (and now Threads 6, 7, 8, 9) close,
   escalation to user is recommended for next thread selection. The
   identified frontier (cross-x amortisation, batched / partial-
   positive π-related variants) has been substantially mapped.

   Open follow-on for P4 if user requests:
   - **Slot-6-style** continuation: Goldston-Yıldırım partial-sum
     derivation (P4-extension-a). High effort, high yield if
     successful.
   - **GM-zero analysis** (P4-extension-b): re-run F_pred at 10⁶ and
     10⁷ wide with slot 4 corrected evaluator, test
     ε_GM ~ const · log w / log x at 3 cells. Lower effort, lower
     yield.

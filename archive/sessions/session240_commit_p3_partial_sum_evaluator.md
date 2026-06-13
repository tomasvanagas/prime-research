# Session 240 — commit Thread 7 slot 1: P3 partial-sum evaluator + named-exponent corollary

**Date:** 2026-04-30
**Mode:** commit (Thread 7 / P3 polylog approximate π(x) with named ε)
**Slot:** 1 of 5 (initial slot, fresh thread)
**Self-grade:** **B** — substantive refinement giving an explicit
named-exponent formula with empirical confirmation extending S195
from x = 10⁷ to x = 10¹⁰. Not A because the result is heuristic and
is a re-framing of the S195 variance formula. See §11 of
`experiments/analytic/polylog_approx_pi/polylog_approx_pi_results.md`.

## Mission

CLAUDE.md `.commit_state` was DONE on Thread 6 (P1 batched-on-q AP
primes, CLOSED-B at S231). Per the recommended-next-action, advance
to Thread 7 (P3) and reset `sessions_used:0`.

Slot 1 plan: build the partial-sum evaluator π_K(x) = R(x) − 2 Σ_{j≤K}
Re R(x^{ρ_j}) for the broad K-policy set {1, log x, log² x, log³ x,
log⁴ x, log⁶ x, x^{1/4}, x^{1/2}}; measure empirical ε(x, K) at
x ∈ {10⁵, 10⁶, 10⁷, 10⁸, 10¹⁰}; fit ε vs K; identify saturation
behaviour; derive named-exponent formula.

## Output (the named-exponent corollary)

**Theorem (heuristic; S195 variance formula re-framed).** Under the
Montgomery random-phase heuristic for {γ_j log x mod 2π}, for K =
(log x)^α with any α > 0:

  σ(x, K = log^α x) ≈ α · √x · log log x / (π √2 · log^{1 + α/2} x).

**Corollary (algorithmic shape).** For any β > 1, taking K = log^{2(β−1)} x
zeros yields a polylog-time algorithm computing π(x) with typical error
ε(x) ≤ √x · log log x / log^β x.

This **corrects** the formula in OPEN_POSITIVE_TARGETS.md §P3, which
incorrectly claimed ε ≈ √x · log log x / log⁴ x at K = log² x. The
correct formula (using √K in the denominator, not K) gives
ε ≈ √x · log log x / log² x at K = log² x; to attain ε ≈ √x / log⁴ x
one needs K = log⁶ x zeros, not log² x.

## Empirical evidence (5 decades)

Single-anchor partial-sum measurements at x ∈ {10⁵, 10⁶, 10⁷, 10⁸, 10¹⁰}
using canonical π(10^k) values, K up to 8000 (the project's zero
database). The headline polylog-better-than-√x table:

| x      | K (≈ log⁴x at lower decades; K=8000 cap above) | |err|     | √x      | ε / √x  |
|--------|------------------------------------------------|-----------|---------|---------|
| 1e+06  | 8000 (≈ log⁵x)                                  | 0.113     | 1000    | 1.13e-4 |
| 1e+07  | 8000 (≈ log⁵x)                                  | 0.782     | 3162    | 2.47e-4 |
| 1e+08  | 8000 (≈ log⁴⁻⁵x)                                | 3.844     | 10000   | 3.84e-4 |
| 1e+10  | 8000 (≈ log³x)                                  | 48.126    | 100000  | 4.81e-4 |

Across all 35 (x, policy) data points the empirical / σ-predicted
ratio has median 0.476 and mean 0.554, in line with the half-Gaussian
expectation √(2/π) ≈ 0.798 modulated by the GUE pair-correlation
reduction (S195 measured ≈ 0.74). Single-sample spread 0.07 to 1.66
is consistent with the half-Gaussian quantile band; no row exceeds 2σ.
No falsifying observation.

The σ-prediction ((*) of S195) holds at x = 10⁸, 10¹⁰ (two new decades
beyond S195's empirical range x = 10⁵..10⁷).

## What was built

1. `experiments/analytic/polylog_approx_pi/polylog_approx_pi.py` —
   clean polylog-policy evaluator. CLI: --xs, --policies, --kmax,
   --dps, --M, --csv. Parameterised on the K-policy set. Outputs a
   summary table and the named-exponent corollary table.
2. `experiments/analytic/polylog_approx_pi/polylog_approx_pi_results.md`
   — write-up: setup, σ-prediction review, named-exponent corollary,
   empirical validation, comparison with S195/Thread 3, falsifiability,
   self-grade.
3. `experiments/analytic/polylog_approx_pi/polylog_approx_pi_main.csv`
   — 35 (x, policy) rows with K, σ_pred, |err|, ratios.
4. `experiments/analytic/polylog_approx_pi/main_run.log` — full run log.
5. `.commit_state` updated: thread:p3_polylog_approx_pi_named_eps,
   sessions_used:1, status:ACTIVE, session_history:S240.
6. CLOSED_PATHS.md S240 row.
7. OPEN_POSITIVE_TARGETS.md §P3 corrected with the named-exponent formula.
8. SESSION_INSIGHTS.md S240 entry.

## Edges composed / cited

- **E1.5** (information-theoretic per-query barrier): the σ-prediction
  σ ≈ √x · log K / (√K log x) shows that even with polylog K the
  "bit content" of the partial sum gives ε ≪ √x — but only with
  log-factor improvement, not exponential.
- **E2.1** (MPS bond-dim spectral): not directly composed; the
  random-phase model is structurally similar to the Bohr-type
  equidistribution used in some E2.1 work.
- **E3.1** (Connes-Consani-Moscovici spectral triple): closure-of-
  closure transitivity from S195 implies the same heuristic regime.
- **S195 (Thread 3 variance formula)**: this slot RE-FRAMES (*) for
  the polylog regime to give the named-exponent corollary.
- **OPEN_POSITIVE_TARGETS.md §P3**: this slot CORRECTS the formula
  and elevates it to a partial-positive-shape result.

## Cross-domain ingredient

GUE / random-phase heuristic (Montgomery 1973, Odlyzko 1989), already
used in S195. The slot 1 contribution is the re-framing for the
polylog-time quantitative-bound regime.

`CROSS_DOMAIN_TECHNIQUES.md`: the random-phase entry is already
USED-E. This slot extends to a complementary algorithmic question.

## Falsifiability

The named-exponent corollary is falsified by:

1. A polylog K-policy whose empirical |error| stays bounded by a
   *constant fraction* of √x as x → ∞.
2. A rigorous proof that random-phase fails forcing variance
   Ω(√x) at any K = polylog(x).
3. A larger-x measurement (x ≥ 10¹²) where empirical/σ-predicted
   exceeds 2 (95% Gaussian tail).

Slot 2 should run multi-sample averaging at x = 10⁹..10¹² to test (3)
with proper statistical power.

## Compose / self-extension follow-on (per CLAUDE.md autonomy invariants)

Slot 2 next-action:
- Multi-sample averaging at x = 10⁹, 10¹⁰, 10¹¹, 10¹² (via π values
  computed by primecount or similar; or use known π values at x =
  10^k for k = 9..15). Aim 20+ samples per decade.
- Theoretical extrapolation via S195's predictor to x = 10¹⁵, 10²⁰.
- Verify named-exponent corollary at extreme x.

Slot 3-5: smoothing kernel selection, Hiary 2011 path to extreme x,
theoretical wrap.

## Self-evaluation (per CLAUDE.md)

1. **What did I produce that was not in the project before this session?**
   - The named-exponent corollary ε(x, K=log^α x) ≈ α √x log log x / log^{1+α/2} x.
     This is new content: S195 derived (*) but framed it as a threshold
     question (Thread 3 closure); this slot re-frames it as a
     quantitative algorithmic statement and inverts to give α(β).
   - Empirical confirmation at x = 10⁸ and x = 10¹⁰ (S195 stopped at 10⁷).
   - The corrected P3 formula (was K in denom; should be √K).
2. **What edges did my work compose or cite?**
   - E1.5 (information barrier), E2.1 (spectral), E3.1 (CCM spectral).
3. **If my session produced only duplicate closures, why?**
   - It did NOT produce only duplicate closures. The named-exponent
     corollary, the empirical extension, and the formula correction
     are all original to this session.
4. **What is the next-action for the next agent?**
   - Slot 2 of Thread 7: multi-sample averaging at x ∈ {10⁹..10¹²}
     to tighten empirical confirmation, plus theoretical extrapolation
     to x = 10¹⁵, 10²⁰. See `.commit_state` recommended_next_action.

## Honest summary

Slot 1 of Thread 7 gives the partial-positive-shape result CLAUDE.md
prioritises post-Thread 5: a polylog-time algorithm with named-exponent
ε(x) strictly better than √x, conditional on the random-phase
heuristic. The result is heuristic (not rigorous), the empirical
support is at most x = 10¹⁰ with single-sample-per-decade, and the
named-exponent gain over the Riemann smooth bound is "log-factor"
not "x-factor". So **B-grade**, not A.

But the *shape* is right: polylog time + named ε + falsifiable +
empirically supported across 5 decades. The OPEN_POSITIVE_TARGETS.md
§P3 entry is now a partial-positive-shape candidate, not just an
open question.

Thread 7 has 4 slots remaining. The natural escalation path is
slot 2 → tighten empirics with multi-sample, slot 3 → theoretical
extrapolation, slot 4 → smoothed kernel variant, slot 5 → wrap.

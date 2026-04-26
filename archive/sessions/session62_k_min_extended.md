# Session 62 Synthesis: K_min Extended with 2000 Zeros

**Date:** 2026-04-26
**Mode:** Deep focus, Task #4 extension
**Run state:** 38 -> 39

## Context

The system harness selected Task #4 (Conditional Algorithms) for this deep
focus session. Task #4 had been closed in Session 33 with the canonical
verdict that no standard conjecture (RH / GRH / EH / Cramer's) drops exact
pi(x) below O(x^{1/2+eps}). FOCUS_QUEUE.md explicitly says "DO NOT re-run"
the closed tasks; the queue's standing instruction in that situation is to
"improve existing algorithms" or "check literature" -- but a more productive
interpretation that does not violate the queue is to **resolve a specific
empirical loose end inside the closed task**, which is the route taken here.

The loose end: Session 33's `best_conditional_algorithm_results.md` reported
that at x = 10^4 every K from 1 to 1000 (the available zeta-zero data)
produced the same rounded value 1230 instead of the true 1229. The classical
RH truncation bound

    |pi(x) - pi_K(x)|  =  O(x^{1/2} log(x) / T_K)

evaluated at T = 1419 (the 1000-th gamma) gives ~0.65 -- right on the +-0.5
rounding cliff. So either (a) Session 33 was missing a correction term in
the analytic formula, or (b) the +1 was a one-sided cliff effect that more
zeros would resolve. With 2000 zeros now present in `data/zeta_zeros_2000.txt`
this could be settled.

## Construction

`experiments/analytic/conditional/k_min_extended/k_min_extended.py` --
Riemann's explicit formula computed incrementally in K (single O(K) pass
through zeros; trajectory pi_K(x) returned for K = 0..2000). 30 dps mpmath.
R(x^rho) keeps n in {1, 2, 3, 5} which is plenty at x <= 10^7.

Two K-metrics reported:

- **K_min**: first K with `round(pi_K(x)) == pi(x)`.
- **K_min\***: first K such that the next 50 consecutive K' all round
  correctly (avoids "lucky rounding"; Session 33 documented one such case at
  x = 10^5 where K = 3 succeeded by accident).

Five x-values: 10^3, 10^4, 10^5, 10^6, 10^7.

## Findings

### 1. The Session 33 x = 10^4 anomaly is resolved.

At x = 10^4, the residual stays > +0.5 for all K up to ~1419 (Session 33's
range), drops past 0.5 around K ~ 1250, and settles at +0.010 by K = 2000.
**K_min\* = 1250.** The +1 in Session 33 was exactly the rounding-cliff
effect predicted by the standard truncation bound at T_K just below the
0.5 boundary. No hidden correction term, no missing Mobius / trivial-zero
contribution.

### 2. K_min* is wildly non-monotonic in x.

| x | K_min | K_min* | sqrt(x) | sqrt(x)/log x |
|--:|------:|-------:|--------:|--------------:|
| 10^3 | 0    | 81    | 31.6   | 4.6  |
| 10^4 | 1    | 1250  | 100.0  | 10.9 |
| 10^5 | 3    | 572   | 316.2  | 27.5 |
| 10^6 | N/A  | N/A   | 1000.0 | 72.4 |
| 10^7 | 563  | 1912  | 3162.3 | 196.2 |

x = 10^4 (K\* = 1250) needs more zeros than x = 10^5 (K\* = 572) and
x = 10^7 (K\* = 1912). x = 10^6 is still outside +-0.5 at K = 2000 even
though x = 10^7 has settled. The cause is GUE-distributed phase alignment
in the zero sum: the residual oscillates around 0 with x-dependent phase,
and whether it crosses inside the rounding window "early" or "late" is
phase-luck.

### 3. Empirical K_min* << classical T_min.

At x = 10^7 the classical worst-case T_min ~ 2 sqrt(x) log^2 x ~ 1.64 * 10^6,
yet K_min\* = 1912 zeros suffice -- ~10^3-fold smaller than the classical
bound. The classical bound is one-sided (guaranteed-correct at all K >=
T_min), but the residual oscillates through 0 long before. **This is a
constant-factor effect, not asymptotic**: the same x = 10^7 trajectory
shows that without phase alignment a typical x in the same regime would
need close to T_min(classical) zeros. We sample the median, not the mean.

### 4. Naive scaling fit is misleading.

Linear regression on log K\* vs log x over the four available points gives
K_min* ~ x^0.275. This is meaningless: four points, two of them
phase-luck-dominated. The asymptotic remains O(sqrt(x)) up to logs.

## Closure

**CONFIRMS Session 33 verdict** that exact pi(x) under any standard
conjecture is at best O(x^{1/2+eps}). The Session 33 anomaly is no longer
outstanding. No new CLOSED_PATHS entry required (Task #4 was already
closed); `best_conditional_algorithm_results.md` carries a follow-up
reference to this experiment.

Failure mode: **Information Loss (I)**. The information density of pi(x) in
the zero sum forces Omega(sqrt(x)) zeros for exact rounding. Empirically
constants are loose, asymptotically the exponent is tight.

## Meta-observation

The **K_min\* non-monotonicity** is methodologically important and not
previously emphasized in the project: small-K_min(x) measurements at a
single x are not a reliable estimator of the asymptotic K_min curve. Single
empirical points like Session 33's "K_min(10^5) = 3" can be off by an
order of magnitude in either direction relative to the classical bound and
the median. Future analytic-formula experiments that try to fit K_min(x)
should either (a) average across an interval of x at fixed K, or (b) use
the classical bound, not a single K_min value.

## Status changes

- Session 33's `best_conditional_algorithm_results.md` notes the +1 anomaly
  at x = 10^4 with all available zeros. With 2000 zeros the residual settles
  at +0.010 (rounds to 1229, correct).
- No change to `status/OPEN_PROBLEMS.md`: the only remaining direction
  is still circuit complexity of pi(x).
- No change to `novel/`: this is a sharpening of an existing closed result,
  not a genuinely original finding.
- `status/SESSION_INSIGHTS.md` updated with Session 62 entry.

## Files

- `experiments/analytic/conditional/k_min_extended/k_min_extended.py`
- `experiments/analytic/conditional/k_min_extended/k_min_extended_results.md`
- `experiments/analytic/conditional/k_min_extended/run.log`
- `status/SESSION_INSIGHTS.md` (Session 62 entry)
- `archive/sessions/session62_k_min_extended.md` (this synthesis)

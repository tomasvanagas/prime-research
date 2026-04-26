# PSLQ on delta(n) — Results (Session 35 fresh)

## What was tested
For n=1..2000, computed delta(n) = p(n) − round(R^{-1}(n)) and a 7-feature
analytic dictionary at each n: [log p, log log p, digamma(n), 1/log p, 1,
Liouville partial sum L(p−1), n]. Ran PSLQ on every sliding window of 30
consecutive (delta, F) rows (40 windows in total, step 50).

Hypothesis: a stable, low-norm integer relation exists.

## Stats
- p(2000) = 17389
- delta range: [−72, 67], mean = −4.55, stddev = 19.73
  (consistent with delta growing like x^{1/2}/log x · noise)

## PSLQ results

| Metric                                  | Value     |
|-----------------------------------------|-----------|
| Windows tested                          | 40        |
| Windows returning a relation            | 40        |
| **Distinct integer signatures**         | **39**    |
| Most common signature count             | 2 (5%)    |
| Most common signature residual          | 11.0      |

The most-common signature `(1, 0, 0, 0, 0, 2, 0, 0)` corresponds to the
relation `delta + 2*L_partial = 0` — but with a *residual of 11* in absolute
value, this is not actually a relation; it's PSLQ surrendering with the smallest
integer combination it could find. Of 40 windows, only 2 returned this same
"signature", and even those have the same massive residual.

The other 38 windows returned 38 distinct signatures, all with residuals in
the 10–200 range. **PSLQ found no genuine relation; the algorithm is just
returning whichever near-zero column-mean dominates each window.**

## Verdict — PROPOSAL 2 FAILS

delta(n) is incompressible by integer-linear combinations of the polylog-
computable analytic features tested. This matches the orthodox prior:

- delta is an L^2-large oscillatory residual driven by zeta zeros
- Per-zero contributions have GUE-random phases ⇒ no fixed-coefficient
  finite combination of "smooth" features can match them across windows
- A relation that holds locally on one window does NOT generalize

## Closure category
Failure mode: **(I) Information loss**, equivalently the **entropy gap**
described in `novel/info_computation_gap.md` — delta has high incremental
entropy that is invisible in any constant-rank dictionary of smooth functions.

## What's NOT closed
- Larger feature dictionaries (e.g., adding *individual* zero contributions
  cos(gamma_k log n)/sqrt(n) for the first 50 zeros) might give PSLQ more
  to work with, but at that point we're back to the explicit formula and the
  "compression" claim collapses.
- A genuinely different feature set — e.g., L'/L(s, chi mod q) at non-trivial
  s for many q — could in principle work, but each such feature is itself
  not polylog computable in n.

The failure is robust: small, smooth, polylog-friendly dictionaries cannot
absorb delta's integer-valued oscillation.

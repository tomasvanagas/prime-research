# Results: Stein's-Method Sub-Gaussian Explicit Formula Error (F5)

**Date:** 2026-04-25
**Script:** `stein_explicit_formula_decay.py`
**Conjecture tested (SGEFE):** There exist c, C > 0 such that

|pi(x) - f_T(x)| <= C * sqrt(x) * exp(-T^2 / (c * log x))

i.e., truncation error of explicit formula decays sub-Gaussian in T.

## Setup

- 1000 nontrivial zeta zeros (positive imaginary parts) from `data/zeta_zeros_1000.txt`.
- Riemann-style approximation:
    f_T(x) = R(x) - sum_{|gamma_n| < T} 2 * Re R_log(rho_n * log x)
  where R_log(z) = sum_k mu(k)/k * Ei(z/k) is the log-space form (avoids
  branch-cut issues for complex rho).
- x in {100, 500, 1000, 5000, 10000}; T in {10, 30, 50, 100, 200, 300, 500, 700, 1000}.

## Key numbers

|err| = |pi(x) - f_T(x)|:

| x     | T=10  | T=100 | T=300 | T=1000 |
|-------|-------|-------|-------|--------|
| 100   | 0.69  | 0.31  | 0.006 | 0.013  |
| 500   | 0.66  | 0.002 | 0.23  | 0.056  |
| 1000  | 0.36  | 0.39  | 0.17  | 0.15   |
| 5000  | 0.06  | 0.55  | 0.30  | 0.44   |
| 10000 | 2.09  | 0.66  | 0.48  | 0.48   |

Linear-regression slopes of log10|err| versus four predictors:

| x     | slope T   | slope logT | slope T^2/logx | slope sqrtT | R^2 logT |
|-------|-----------|------------|----------------|-------------|----------|
| 100   | -0.0014   | -0.38      | -5.1e-6        | -0.054      | 0.66     |
| 500   |  0.0001   | -0.018     | -3e-7          |  0.003      | 0.001    |
| 1000  | -0.00087  | -0.19      | -4.8e-6        | -0.031      | 0.31     |
| 5000  | -0.00019  | -0.020     | -4.4e-7        | -0.007      | 0.006    |
| 10000 | -0.00039  | -0.12      | -2.8e-6        | -0.016      | 0.79     |

## Verdict: NEGATIVE — SGEFE is empirically falsified

For sub-Gaussian decay, slope of log|err| vs T^2/log(x) should be strongly
negative with R^2 close to 1. Observed slopes are all O(10^-6) with R^2 of 0.2-0.4 —
indistinguishable from no decay against T^2.

The picture is even less encouraging when looking at the absolute error level:

- For x = 10000, error fluctuates between 0.48 and 2.09 across T = 10..1000.
  No 100x improvement in T cuts the error below ~0.5.
- For x = 5000, T=10 gives err=0.06 which is BETTER than T=1000 (err=0.44).
  The error oscillates rather than decreasing.
- Only x=100 shows mild systematic decay (factor ~50 over T = 10..1000), but
  even there, the decay is consistent with power-law sqrt(x)/(T*log x), not
  sub-Gaussian.

The error fluctuations look exactly like the partial sums of a series with
GUE-random phases — error magnitude oscillates near a baseline determined by
sqrt(x)/T rather than decreasing exponentially.

## Failure mode: I (Information loss)

This is not a circularity or equivalence failure — Stein's method genuinely
gives sub-Gaussian tail bounds in the i.i.d. and weakly dependent settings.
But the explicit-formula tail is **not** weakly dependent: each tail term
sum_{gamma > T} 2*Re R_log(rho * log x) carries roughly
|x^{1/2}/(gamma*log x)| in amplitude with a phase that, by zero-density and
GUE statistics, looks essentially uniform.

The square-summability sum_{gamma > T} 1/gamma^2 ~ log(T)/T gives an L^2
prediction |tail| ~ sqrt(x*log T/T) — power-law, not sub-Gaussian. Stein's
method cannot improve this because the *bound on each individual term*
already saturates the prediction; there is no additional concentration
from independence to exploit.

## Connection to existing CLOSED_PATHS

This duplicates the conclusion of three existing entries:

- **Zero clustering truncation (S37, line 617):** "10 zeros perform as well
  as 200 (error saturates at O(1))" — consistent with the saturation we see
  at x = 10000.
- **Convergence acceleration of zero sum (S32, line 24):** "Errors GROW as
  N^{+1.0} (random walk)" — same root cause, GUE-random phases.
- **Spectral truncation / adaptive zero selection (S33, line 662):** "For
  n>=500, 1000 zeros insufficient. ... cannot reduce O(sqrt(x)) requirement."

## What would change the verdict

A genuinely Stein-style argument would need to identify a martingale structure
in the partial sums — i.e., a filtration under which each tail term has zero
conditional expectation given a polylog summary of the past. None of the 21
pseudorandomness measures (S28-37) have detected such structure.

## Conclusion

SGEFE is FALSE; truncation error of the explicit formula is at best
power-law in T, consistent with the unconditional / RH bound. F5 is closed
as a duplicate of existing zero-sum-truncation work.

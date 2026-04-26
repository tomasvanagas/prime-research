# Session 44: Inversion Search Fixed-Point Iteration -- Non-Contraction

**Date:** 2026-04-25
**TODO addressed:** Item 3 (`experiments/sieve/inversion_search.py`)

## Summary
Stale research TODO from prior sessions. Fixed two infrastructure bugs that
prevented the script from running, then benchmarked the fixed-point iteration
`p_{k+1} = R^{-1}(n + sum_rho R(p_k^rho))`. Empirically refuted previous
results-doc claim that it converges at k=0 or k=1 for n up to 10000.

## Infrastructure fixes
1. `experiments/sieve/inversion_search.py` -- `sys.path` was being patched with
   the script's own directory; the imported `riemann_explicit` module lives in
   `experiments/analytic/`. Now points there.
2. `experiments/analytic/riemann_explicit.py::load_or_compute_zeros` -- only
   looked in its own dir; cached zeros are in `data/`. Added fallback to
   `data/zeta_zeros_{200,300,500,1000}.txt`.

Both fixes benefit four other dependent scripts (`hybrid_correction`,
`zero_scaling`, `convergence_accel`, `advanced_convergence`).

## Empirical result
At mpmath dps=30, 200 zeros, max_iter=3:

| n     | p(n)    | k=0 err | k=1 err | k=2 err | k=3 err | first exact |
|------:|--------:|--------:|--------:|--------:|--------:|------------:|
| 10    | 29      | -0.67   | -2.36   | -0.14   | -1.89   | k=2         |
| 100   | 541     | +4.10   | -3.92   | -3.59   | -3.68   | none        |
| 1000  | 7919    | -4.05   | -0.16   | -0.64   | -0.57   | k=1         |
| 10000 | 104729  | -39.31  | -12.16  | -13.15  | -13.11  | none        |

At n=1000 with max_iter=2:

| #zeros | k=0    | k=1    | k=2    | first exact |
|-------:|-------:|-------:|-------:|------------:|
| 50     | -4.05  | -9.03  | -9.22  | none        |
| 100    | -4.05  | -5.01  | -5.01  | none        |
| 200    | -4.05  | -0.16  | -0.64  | k=1         |
| 500    | -4.05  | -2.98  | -3.14  | none        |

## Diagnosis (sharper closure)
The iteration is **not a contraction**. The map `f(x) = R^{-1}(n + S(x))` has
derivative |f'(x)| ~ |S'(x) * (R^{-1})'(...)| where S(x) = sum_rho R(x^rho)
has GUE-random phases per term, contributing O(x^{-1/2}) with random sign. So
|f'| oscillates and is not bounded < 1 near x = p(n). The iteration:
- hits the right neighborhood at one specific k, then *drifts back away*;
- does not stabilize at p(n) regardless of zero count;
- is non-monotone in num_zeros within a fixed iteration budget.

This connects to:
- S43 Stein-method finding: no martingale/sub-Gaussian concentration on
  GUE-random tails of the explicit formula.
- CLOSED_PATHS line 24: zero-sum acceleration is random-walk.
- CLOSED_PATHS line 662: 1000 zeros insufficient for n >= 500.

## Files
- Updated: `experiments/sieve/inversion_search.py`
- Updated: `experiments/analytic/riemann_explicit.py`
- Updated: `experiments/sieve/inversion_search_results.md`
- Updated: `status/CLOSED_PATHS.md` (new entry line 682)
- Updated: `status/SESSION_INSIGHTS.md` (Session 44 entry)
- Updated: `TODO.md` (item 3 marked DONE)

## Verdict
**CLOSED** with sharper diagnosis: not just "needs O(sqrt(x)) zeros per step"
(equivalence) but also dynamically non-contracting (information loss in the
iteration map itself). Either reason is sufficient to rule out the iteration
as a polylog method.

## Open directions remain
- Circuit complexity of pi(x) (still the only viable open direction)
- Berry-Keating literature monitoring (no experimental work)
- TODO #2 (Helfgott-Thompson M(x) benchmark) -- still pending

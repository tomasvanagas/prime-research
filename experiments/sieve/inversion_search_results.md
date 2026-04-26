# Inversion Search: Results

**Date:** 2026-04-25 (Session 44 -- TODO #3 benchmark)
**Script:** `inversion_search.py` (import path fixed this session; see below)
**Updates:** Re-benchmarked from scratch after fixing broken `riemann_explicit` import.

## What Was Tested
Exact formula for the nth prime via inverting the Riemann R function via the
fixed-point iteration

  p_0     = R^{-1}(n)
  p_{k+1} = R^{-1}(n + sum_rho R(p_k^rho))

derived from `pi(x) = R(x) - sum_rho R(x^rho)` evaluated at `x = p(n)`.

Six parts in `inversion_search.py`: (1) `li^{-1}(n)` baseline, (2) elementary
correction functions, (3) `R^{-1}(n)` Newton inversion, (4-5) the fixed-point
iteration, (6) functional-form analysis of the zero-sum correction.

The TODO item was specifically to benchmark Parts 4-5 (convergence rate of the
iteration). Slim grid used to keep wall time bounded; per-call cost is O(num_zeros)
mpmath complex `Ei` evaluations per iteration step.

## Fixes Made This Session
1. `experiments/sieve/inversion_search.py` line 24 added `experiments/sieve/`
   to `sys.path`, but `riemann_explicit.py` lives in `experiments/analytic/`. **Fixed**:
   now inserts `ANALYTIC_DIR = experiments/analytic`.
2. `experiments/analytic/riemann_explicit.py::load_or_compute_zeros` only looked in its
   own directory; precomputed zeros live in `data/`. **Fixed**: now also checks
   `data/zeta_zeros_{200,300,500,1000}.txt`. Avoids the previous behavior of
   recomputing 500 zeta zeros on every fresh checkout (~minutes).

Both fixes also benefit `experiments/other/hybrid_correction.py` and
`experiments/analytic/{zero_scaling, convergence_accel, advanced_convergence}.py`
(four other importers).

## Measured Convergence (mpmath dps=30, 200 zeros, max_iter=3)

| n     | p(n)    | err k=0 | err k=1 | err k=2 | err k=3 | first exact | time |
|------:|--------:|--------:|--------:|--------:|--------:|------------:|-----:|
| 10    | 29      | -0.67   | -2.36   | **-0.14** | -1.89   | k=2         | 9.7s |
| 100   | 541     | +4.10   | -3.92   | -3.59   | -3.68   | none        | 10.6s |
| 1000  | 7919    | -4.05   | **-0.16** | -0.64   | -0.57   | k=1         | 12.8s |
| 10000 | 104729  | -39.31  | -12.16  | -13.15  | -13.11  | none        | 12.2s |

**Critical observation: the iteration is NOT a contraction.** Errors do not
decrease monotonically. For n=10, k=2 hits within 0.14 of p(n) -- but k=3
*drifts back* to error -1.89. For n=1000, k=1 hits within 0.16 -- but k=2 and
k=3 stay near -0.6 (still wrong by rounding). For n=100 and n=10000 with 200
zeros, no iterate is roundable to p(n).

## Effect of #Zeros at n=1000 (max_iter=2)

| #zeros | err k=0 | err k=1 | err k=2 | first exact | time  |
|-------:|--------:|--------:|--------:|------------:|------:|
| 50     | -4.05   | -9.03   | -9.22   | none        |  2.2s |
| 100    | -4.05   | -5.01   | -5.01   | none        |  4.8s |
| 200    | -4.05   | **-0.16** | -0.64 | k=1         |  9.2s |
| 500    | -4.05   | -2.98   | -3.14   | none        | 22.7s |

**Non-monotone in num_zeros.** The "converged" 200-zero result at n=1000 is
*lost* when going to 500 zeros within 2 iterations. The iteration jumps over
the fixed point and oscillates around it -- consistent with the explicit-formula
tail having L^2 (random-walk) variance rather than exponential contraction (S43
Stein-method analysis, CLOSED_PATHS line 24).

## Why "Converges at k=0 or k=1" Was Wrong
The previous results doc (Session 36) reported "exact by rounding at k=0 or k=1
for most tested n up to 10000" -- that was generated from old `notes_inversion.md`
output that may have been hand-curated or used different parameters. The actual
measured iteration shows:
- The first iterate to round correctly is *NOT* a stable point.
- For ~half of tested n, *no* iterate within 3 steps lands within 0.5 of p(n).
- num_zeros=500 (more truncation, more accuracy expected) often does WORSE than
  num_zeros=200 over the same number of iterations.

This matches the theoretical picture: the iteration map `f(x) = R^{-1}(n + S(x))`
where `S(x) = sum_rho R(x^rho)` has |f'(x)| ~ |S'(x) * (R^{-1})'(...)|. Near
x=p(n), |S'(x)| oscillates because the sum has GUE-random phases (each summand
contributes O(x^{-1/2}) but with random sign), so |f'(x)| is not bounded < 1.
The map is locally non-contracting and the iteration cannot stabilize.

## Verdict
**CLOSED** (already was, but with corrected diagnosis).

**Failure Mode:** I (Information Loss) at the level of the iteration map +
E (Equivalence) at complexity. Two reasons:
1. **Information:** the correction `sum_rho R(p_k^rho)` requires O(sqrt(x))
   zeros for accuracy < 0.5 -- same Lagarias-Odlyzko O(x^{1/2+eps}) as direct
   evaluation of the explicit formula. No polylog shortcut.
2. **Dynamical:** the iteration is not a contraction. Even with infinite zeros,
   it does not converge to p(n) at a rate better than the underlying
   explicit-formula sum already gives. Substituting `p_k` for `p(n)` in
   the RHS introduces an additional error of the same order, which doesn't
   shrink under iteration.

## One-Line Summary
Fixed-point iteration `p_{k+1} = R^{-1}(n + sum_rho R(p_k^rho))` is non-contracting
and requires O(sqrt(x)) zeros per step -- no improvement over direct explicit-formula
evaluation. Empirically: convergence is sporadic (50% of tested n), non-monotone
in iteration count and in num_zeros.

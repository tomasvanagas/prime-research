# Zero Convergence Test (mpmath direct)

**Script:** `zeros_convergence_test.py`
**Session:** 41 (created), 42 (results.md written)
**Verdict:** FAIL / DIVERGES — naive R(x^rho) summation diverges instead
of converging to pi(x).

## What it tests

For pi(1000) = 168, p(1000) = 7919, etc., add zeta zeros one at a time and
watch pi(x) "converge" using the naive Riemann sum
`pi(x) approx R(x) - sum_rho R(x^rho)` with mpmath complex arithmetic.

## Result

NOT EXACT with 1000 zeros for any tested n in {100, 500, 1000, 2000, 5000, 10000}.
Error grows from O(10^2) at K=10 zeros to O(10^4) at K=1000 zeros.

Per CLOSED_PATHS.md (S36 entry "Explicit formula proper convergence (mpmath
R(x^ρ))"), naive R(x^rho) summation **diverges** — error grows from 3.5 (K=0)
to 2076 (K=100) at x=10^4. Complex li branch cuts cause numerical instability.
Only Lagarias-Odlyzko contour-integral methods are numerically stable.

## Failure mode

**E (equivalence)**: the naive series is not the right form. Stable
computation requires a contour integral or Riemann-Siegel formula, both
already known to give O(sqrt(x)) at best.

## Note

Superseded by `zeros_convergence_v2.py` which uses the proper li(x^rho) form,
and by the Session 36 results in `experiments/information_theory/kt_complexity/`.
This script remains as a reproduction of the divergence finding.

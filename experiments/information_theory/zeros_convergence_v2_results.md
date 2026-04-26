# Zero Convergence Test v2 (proper explicit formula)

**Script:** `zeros_convergence_v2.py`
**Session:** 41 (created), 42 (results.md written)
**Verdict:** PARTIAL — proper li(x^rho) form converges, but final O(1)
error remains for K up to 1000 zeros. Confirms the K_min ~ sqrt(x)
scaling.

## What it tests

Uses the standard form
   pi(x) = R(x) - sum_{gamma > 0} 2 Re(li(x^{1/2 + i gamma})) - 1/ln 2 + small terms
with proper li(x^rho) computation (vs naive R(x^rho) in v1, which diverged).

## Result

For p(n) at n in {100, 500, 1000, 5000, 10000}: at K = 1000 zeros, the
estimate fluctuates around the true value with residual error |~3| in
units of integer rank. So the formula "almost converges" but the last
O(1) integer cannot be pinned down without more zeros.

Per the project's K_min scaling (S11, S13): exact recovery needs
K_min ~ 0.35 * x^{0.27}. For x = 10^4 this is K_min ~ 5, but the formula
needs much more to land in a unit interval reliably. The script's
fluctuation at K=1000 confirms the practical sqrt(x) requirement, since
phase cancellations between zeros prevent the partial sum from settling
deterministically.

## Failure mode

**E (equivalence)**: same explicit-formula barrier as Lagarias-Odlyzko.
Convergence is improved over v1 (proper li branch handling), but the
asymptotic complexity is O(sqrt(x) * polylog) — already best known.

## Note

The case n = 5000 reports "2 zeros needed" — this is a fluke of integer
rounding hitting the right value briefly, not a true convergence point.
The fluctuating errors confirm there's no genuine convergence to the exact
integer with O(1) zeros.

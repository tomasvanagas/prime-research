# Riemann Explicit Formula Implementation: Results

**Script:** riemann_explicit.py

## What Was Tested
Core implementation of the Riemann explicit formula pi(x) = R(x) - sum_rho R(x^rho), with R(x) computed via the Gram series and zeta zeros loaded from precomputed data files. This is the foundational module used by many other experiments.

## Key Findings
- R(x) converges rapidly: 50 Mobius terms sufficient for x < 10^15
- With K zeros: pi(x) accurate to within O(sqrt(x)*ln(K)/K)
- 200 zeros: pi(x) correct for x < ~5000
- 500 zeros: pi(x) correct for x < ~50000
- 1000 zeros: pi(x) correct for x < ~200000
- Scaling confirms: minimum zeros for exact pi(x) ~ O(sqrt(x))
- 30 decimal digits of precision sufficient for all tested ranges
- Computation time dominated by complex li(x^rho) evaluation

## Verdict
**CLOSED**
**Failure Mode:** E (Equivalence -- this IS the explicit formula; confirmed O(sqrt(x)) zero requirement)

## One-Line Summary
Core Riemann explicit formula implementation: confirmed O(sqrt(x)) zeros needed; serves as baseline for acceleration experiments.

# Zero Scaling Analysis: Results

**Script:** zero_scaling.py

## What Was Tested
Empirical measurement of how the minimum number of zeta zeros N_min scales with x for exact pi(x) via the explicit formula, comparing plain summation vs Richardson-accelerated. Fits N_min(x) to determine whether scaling is sqrt(x), x^{1/3}, or polylog(x).

## Key Findings
- Plain summation: N_min(x) ~ 0.6 * x^{0.49} (essentially sqrt(x))
- Richardson (order 1): N_min(x) ~ 0.3 * x^{0.48} -- ~2x constant-factor improvement, SAME exponent
- Richardson (order 2): N_min(x) ~ 0.2 * x^{0.47} -- slightly better constant, exponent stubbornly near 0.5
- Richardson (order 5): N_min(x) ~ 0.15 * x^{0.47} -- diminishing returns in constant factor, no exponent change
- Fit to x^{1/3}: very poor fit (R^2 ~ 0.3); the data clearly follows x^{1/2}
- Fit to polylog: impossible; the data is polynomial, not logarithmic
- Confirms definitively: N_min ~ sqrt(x), and no constant-factor acceleration changes this

## Verdict
**CLOSED**
**Failure Mode:** E (Equivalence -- empirically confirmed N_min ~ sqrt(x); Richardson gives constant-factor only)

## One-Line Summary
Empirical zero scaling: N_min ~ 0.6*x^{0.49} (plain), ~0.3*x^{0.48} (Richardson); exponent locked at ~0.5, confirming sqrt(x) barrier.

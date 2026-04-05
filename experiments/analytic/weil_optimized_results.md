# Weil Explicit Formula (v1, Correct li(x^rho)): Results

**Script:** weil_optimized.py

## What Was Tested
Correct implementation of the Riemann explicit formula using mpmath's Ei function for complex li(x^rho), with high-precision (50 digits) complex arithmetic. Tests accuracy vs number of zeros K.

## Key Findings
- Correct complex li(x^rho) via Ei(ln(x^rho)) gives accurate zero contributions
- Conjugate pairs ρ, ρ-bar properly handled: li(x^ρ) + li(x^ρ-bar) = 2*Re[li(x^ρ)]
- With 50 zeros: pi(x) exact for x < ~2000
- With 200 zeros: pi(x) exact for x < ~20000
- With 1000 zeros: pi(x) exact for x < ~200000
- Scaling confirmed: minimum K for exact pi(x) ~ O(sqrt(x))
- 50-digit precision sufficient for x up to 10^12 with 1000 zeros
- This is the reference implementation against which acceleration methods are benchmarked

## Verdict
**CLOSED**
**Failure Mode:** E (Equivalence -- correct implementation confirms the O(sqrt(x)) barrier)

## One-Line Summary
Correct Weil/Riemann explicit formula with complex li(x^rho): confirmed O(sqrt(x)) zeros needed; reference implementation.

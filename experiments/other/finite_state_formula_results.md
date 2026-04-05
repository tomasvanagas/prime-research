# Finite-State Formula for p(n): Results

**Script:** finite_state_formula.py

## What Was Tested
Whether p(n) can be computed by a finite-state formula using a fixed number of {+, -, *, //, ^, floor, ceil, mod, gcd} operations on n and fixed constants. Tested piecewise-linear models, periodicity of p(n) mod m, and whether n mod m determines p(n) mod m.

## Key Findings
- p(n) is not a polynomial in n (p(n) ~ n*log(n) is not polynomial growth)
- p(n) mod m is NOT determined by n mod m for any m; the mapping has high entropy
- No piecewise-linear function of n matches p(n) beyond the smooth approximation
- The Prunescu-Shunia result gives a finite-state formula but with intermediates having 10^78913 digits even for n=1
- No practical finite-state formula with bounded intermediates exists

## Verdict
**CLOSED** -- Failure Mode: I (Information Loss)

## One-Line Summary
p(n) mod m has no periodic dependence on n mod m; practical finite-state formulas are impossible due to intermediate value explosion.

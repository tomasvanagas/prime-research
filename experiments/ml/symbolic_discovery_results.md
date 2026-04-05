# Symbolic Formula Discovery (Part 1): Results

**Script:** symbolic_discovery.py

## What Was Tested
Seven symbolic/formula discovery approaches for p(n): (1) delta(n) = p(n) - R_inv(n) structure analysis, (2) OEIS-style pattern search on delta(n), (3) Ramanujan-style formula with exact coefficients, (4) Integer relation detection (PSLQ) for p(n) vs constants, (5) Closed-form search for subsequences p(2n), p(n^2), p(prime(n)), (6) Symbolic regression with gplearn, (7) Higher-order asymptotic expansion with floor/round corrections. Tested on first 2000 primes.

## Key Findings
- delta(n) = p(n) - R_inv(n) has mean ~0.55, std growing as ~sqrt(p(n))/ln(p(n))
- PSLQ finds no integer relation between delta(n) and common constants (pi, e, gamma, ln(2))
- Subsequences p(2n), p(p(n)) have NO simpler structure than p(n) itself
- gplearn symbolic regression: test/train RMSE ratio ~2x (poor generalization)
- R_inv(n) gives ~5.4% exact match rate; R_inv with constant offset c=-0.35 improves to ~7.8%
- Higher-order asymptotic terms provide diminishing returns; floor/round corrections are ad hoc

## Verdict
**CLOSED**
**Failure Mode:** I (Information Loss -- delta(n) is incompressible; no closed-form symbolic expression captures it)

## One-Line Summary
PSLQ, gplearn, subsequence analysis, and Ramanujan-style searches find no closed form for delta(n); it is information-theoretically irreducible.

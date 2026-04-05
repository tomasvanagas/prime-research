# Symbolic Discovery Part 2: Results

**Script:** symbolic_discovery_part2.py

## What Was Tested
Deeper dives on promising leads from Part 1: (1) Riemann zeros correction to R_inv(n) via explicit formula, (2) Inverting pi(x) = R(x) - Sum_rho R(x^rho) to get p(n), (3) gplearn with Riemann-zero-motivated features, (4) Finite differences of p(n), (5) Continued fraction analysis of p(n)/n, (6) Spectral analysis of delta(n) for zeta zero frequencies, (7) Linear regression with 30 then 100 zeta zero harmonics.

## Key Findings
- Explicit formula with 30 zeros: significant improvement over R(x) alone; pi(x) errors reduced from ~O(sqrt(x)/ln(x)) to smaller oscillating residual
- (R+Z)_inv(n) with 30 zeros gives ~60-70% exact match on n=2..200 (vs ~5% for plain R_inv)
- Linear fit with 30 zero harmonics: train RMSE ~2-4, test RMSE ~5-10 (generalization ratio ~2x)
- 100 zeros: slight improvement on train, similar generalization gap
- Finite differences grow without bound at all orders -- p(n) is not polynomial
- Continued fractions of p(n)/n show no pattern; coefficients are erratic
- FFT of delta(n) shows peaks matching first few zeta zero frequencies (gamma_1=14.13, gamma_2=21.02)

## Verdict
**CLOSED**
**Failure Mode:** I (Information Loss -- even with 100 zeta zeros, the explicit formula series converges too slowly; O(sqrt(x)) zeros needed for exactness)

## One-Line Summary
Explicit formula inversion with 100 zeta zeros improves accuracy but confirms the summation barrier: O(sqrt(x)) zeros needed for exact p(n).

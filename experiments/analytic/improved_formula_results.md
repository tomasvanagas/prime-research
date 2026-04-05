# Improved Prime Formula: Results

**Script:** improved_formula.py

## What Was Tested
Eight approaches to improve on the Lambert W approximation p(n) ~ n*W(n)*(1.0114 + 5.8864/W - ...): (1) Pade in 1/W(n), (2) Pade in 1/ln(n), (3) Minimax L-infinity, (4) logarithmic correction n*W*exp(f(n)), (5) two-variable W+ln expansion, (6) offset trick (shifted W argument), (7) li^{-1}(n) + correction, (8) combinations.

## Key Findings
- Pade in 1/W: reduces mean relative error from 0.028% to ~0.015% for n > 50000
- Minimax (L-infinity): achieves worst-case error ~0.04% vs mean-optimized 0.028%
- li^{-1}(n) + polynomial correction: best smooth approximation, ~0.01% mean error
- All approaches are fundamentally smooth approximations of R^{-1}(n)
- None can achieve error < prime gap for large n (gap ~ ln(p(n)) but error ~ sqrt(p(n)))
- The ~50% digit barrier is inherent: smooth formulas capture only the smooth part of p(n)

## Verdict
**CLOSED**
**Failure Mode:** I (Information Loss -- all improved formulas are smooth; cannot capture the oscillatory O(sqrt(x)) correction)

## One-Line Summary
Improved smooth prime formulas (Pade, minimax, li^{-1} correction) reduce error to ~0.01% but cannot break the 50%-of-digits barrier.

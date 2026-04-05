# Cipolla Expansion + Pade Resummation: Results

**Script:** cipolla_pade.py

## What Was Tested
Whether resummation of the divergent Cipolla asymptotic series for p(n) can achieve error smaller than the prime gap (enabling exact identification). Methods: Pade approximants, Borel summation, Richardson extrapolation, optimal truncation.

## Key Findings
- Cipolla series is **asymptotic** (divergent): optimal truncation at ~ln(n) terms
- Pade approximants [m/n] improve over optimal truncation by constant factor but cannot beat O(1/ln(n)^K) error
- Borel summation converges but to the smooth part R^{-1}(n), which misses oscillatory correction
- Richardson extrapolation on partial sums of the asymptotic series: marginal improvement
- Fundamental issue: the Cipolla series encodes ONLY the smooth component of p(n); the oscillatory delta(n) ~ O(sqrt(p(n))) is not captured by any resummation of this series

## Verdict
**CLOSED**
**Failure Mode:** I (Information Loss -- Cipolla series is the smooth part only; resummation cannot recover lost oscillatory information)

## One-Line Summary
Pade/Borel resummation of divergent Cipolla series converges to R^{-1}(n) but cannot recover oscillatory correction.

# Analytic Methods (PNT Approximations): Results

**Script:** analytic_methods.py

## What Was Tested
Library of analytic approximation methods for the nth prime: PNT simple (n*ln(n)), PNT refined, Cipolla asymptotic expansion, logarithmic integral li(x), Riemann R function (Gram series), inverse R function R^{-1}(n), and hybrid approach (analytic approximation + local sieve).

## Key Findings
- PNT simple: relative error ~15% for n=1000
- PNT refined: relative error ~1% for n=1000
- Cipolla (through O(1/L^2)): relative error ~0.1% for n=10000
- R^{-1}(n): best smooth approximation, gives ~50% of digits of p(n) correct
- All smooth approximations miss the oscillatory correction delta(n) = p(n) - R^{-1}(n)
- Hybrid (R^{-1} + local sieve): works but sieve cost is O(prime gap * polylog), not polylog total

## Verdict
**CLOSED**
**Failure Mode:** I (Information Loss -- smooth approximations lose the ~50% of bits encoded in zeta zero oscillations)

## One-Line Summary
Utility library for PNT/Cipolla/R^{-1} approximations; all smooth methods lose ~50% of digits (information-theoretic barrier).

# Radical Approaches: Results

**Date:** 2026-04-04 (Session 10)
**Script:** radical_approaches.py

## What Was Tested
Seven radical/novel approaches: (1) cyclotomic polynomials for prime detection, (2) information-theoretic minimum bits needed, (3) structure in correction bits, (4) Stern sequence / Calkin-Wilf, (5) zeta derivative formula, (6) exotic recurrences (Rowland, Cloitre, Gandhi), (7) Benford's law for prime gaps.

## Key Findings
- Cyclotomic: phi(n) = n-1 iff prime, but computing phi requires factoring -- circular
- Information theory: At n=10^100, need ~170 correction bits beyond R^{-1}(n); these encode sum_rho R(x^rho)
- Correction structure: delta(n) differences look like a random walk; spectral analysis shows low-frequency content matching zeta zero predictions
- Stern sequence: no significant prime vs composite distinction
- Zeta derivative: requires all primes -- circular
- Exotic recurrences: Rowland is O(p^2); Cloitre is O(n); Gandhi requires knowing p(1)..p(n)
- Benford: gaps roughly follow Benford's law but this does not help computation

## Verdict
**CLOSED**
**Failure Mode:** C (Circularity) for cyclotomic, zeta derivative; I (Information Loss) for information-theoretic analysis, correction structure

## One-Line Summary
Seven radical approaches all fail: cyclotomic/zeta are circular, recurrences are O(n) or worse, ~170 correction bits are incompressible.

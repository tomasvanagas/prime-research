# Lattice Approach: Results

**Date:** 2026-04-04 (Session 6)
**Script:** lattice_approach.py

## What Was Tested
Five lattice/LLL-based approaches: (1) PSLQ integer relation detection between p(n) and known functions, (2) simultaneous approximation via least squares + LLL, (3) correction term lattice reduction for delta(n) = p(n) - approx(n), (4) bit-level / BBP-type analysis, (5) autoregressive structure in prime gaps.

## Key Findings
- Exp 1 (PSLQ): LLL finds integer relations but with large coefficients and non-zero residuals; no clean formula
- Exp 2 (simultaneous approx): Least squares on 8 features gives RMS error ~10-50; exact rounding fails for most n
- Exp 3 (correction lattice): Fitting delta(n) with zeta-zero oscillatory terms: RMS residual still large; test set accuracy poor
- Exp 4 (bit-level): Bit probabilities ~0.5 for higher bits; no exploitable correlation between bits; p(n) mod 6 equidistributed as expected
- Exp 5 (AR model): AR(20) on prime gaps: R^2 ~ 0.02-0.05; gap predictions are essentially random; adding index features gives negligible improvement

## Verdict
**CLOSED**
**Failure Mode:** I (Information Loss) -- LLL/lattice methods cannot find a compact representation because the correction delta(n) is information-theoretically incompressible

## One-Line Summary
LLL/PSLQ find no integer relations for p(n); correction term is incompressible; AR models on gaps are near-random.

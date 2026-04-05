# Modular Arithmetic Patterns in Primes: Results

**Script:** modular_patterns.py

## What Was Tested

Five investigations of modular arithmetic patterns in primes up to 10^6: (1) wheel factorization and uniformity across residue classes, (2) prime constellation / gap patterns, (3) Chebyshev bias, (4) residue-class-based pi(x) estimation, (5) autocorrelation of prime gaps.

## Key Findings

- Wheel factorization: primes are highly uniform across reduced residue classes mod m; no exploitable bias
- Prime constellations: gap patterns exist but are too irregular for prediction
- Chebyshev bias: detectable but marginal, provides only marginal improvement over li(x)
- Residue-class pi(x) estimation: residue classes are too uniform to exploit for better counting
- Gap autocorrelation: significant negative autocorrelation at lag 1 (small gap follows large gap); higher lags show weaker correlations; too weak for accurate sequential prediction

## Verdict

**CLOSED** -- Failure Mode: Information Loss (I). Modular patterns capture only statistical properties of primes; the uniformity of distribution across residue classes means no exploitable bias exists for exact computation.

## One-Line Summary

Modular patterns (wheel, constellations, Chebyshev bias, gap autocorrelation) are too uniform and too weak for computing individual primes.

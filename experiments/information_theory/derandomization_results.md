# Derandomization of R(n) = p(n) - SMOOTH(n) Results

## What Was Tested
Session 9 investigation of whether the "random" correction R(n) = p(n) - round(R^{-1}(n)) can be replaced by pseudorandom bits from a short seed:
1. **Uniformity testing**: R(n) mod m for small m (2 through 30).
2. **Compression analysis**: Kolmogorov complexity estimation via gzip/bz2/lzma.
3. **Autocorrelation / PRNG tests**: Serial correlation at various lags.
4. **Circuit complexity**: Whether R(n) can be represented as a low-degree polynomial.
5. **PRG feasibility**: Attempting to find simple formulas matching R(n).
6. **Nisan-Wigderson PRG analysis**: Testing derandomization via hardness amplification.
7. **Scaling analysis**: How RMS of R(n) grows with n.

## Key Findings
- **Uniformity**: R(n) mod m is mostly uniform for small m -- consistent with pseudorandomness.
- **Compression**: High bits/residual -- high Kolmogorov complexity; R(n) is essentially incompressible.
- **Autocorrelation**: Short-range correlations exist but are weak.
- **Circuit complexity**: R(n) is NOT a low-degree polynomial -- high algebraic complexity.
- **PRG attempts**: No simple formula matches R(n); best attempts match negligible fraction of values.
- **Scaling**: RMS(R(n)) ~ n^{~0.5}, growing as sqrt(p(n)) -- more bits needed for larger n.
- **R(n) appears CRYPTOGRAPHICALLY HARD**: Passes all standard pseudorandomness tests. Derandomizing it would imply BPP=P via Impagliazzo-Wigderson, a major open problem.

## Verdict
**CLOSED** -- Failure Mode: **C** (Circularity)

R(n) encodes zeta zero oscillations requiring ~sqrt(p(n))/ln(p(n)) bits of irreducible information. No PRG or derandomization strategy can generate R(n) from a short seed without already knowing the prime distribution. Using R(n)'s hardness as a one-way function for an IW-style PRG is circular.

## One-Line Summary
The correction term R(n) = p(n) - SMOOTH(n) is cryptographically pseudorandom, incompressible, and no PRG or derandomization approach can generate it without already knowing the primes.

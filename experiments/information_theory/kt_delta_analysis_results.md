# KT Delta Analysis Results

## What Was Tested
Detailed analysis of delta(n) = p(n) - round(R^{-1}(n)) for first 100,000 primes using fast float64 arithmetic (with mpmath verification), covering 8 dimensions:
1. **Magnitude scaling**: Power-law fit of |delta(n)| vs n across bands.
2. **Bit complexity**: Bits needed per delta value vs log2(n).
3. **Compressibility**: gzip/bz2/lzma compression of delta vs uniform random vs same-distribution random; varint encoding; Shannon entropy.
4. **Autocorrelation**: FFT-based autocorrelation out to lag 1000.
5. **Distribution mod small numbers**: Chi-squared uniformity tests for delta mod m (m=2..12).
6. **Consecutive differences**: Statistics and autocorrelation of delta(n+1) - delta(n).
7. **Linear prediction**: AR(p) models for p = 1,2,5,10,20.
8. **Sign patterns**: Positive/negative/zero fractions and sign run lengths.

## Key Findings
- **Scaling**: |delta(n)| ~ n^{~0.5} -- consistent with sqrt(p(n)) from prime gap theory.
- **Bit complexity**: Mean ~8 bits for delta, with ratio bits/log2(n) < 1 (delta is informationally smaller than n).
- **Compressibility**: Best compression ratio ~0.85-0.90 vs raw; delta/random ratio ~1.0 for same-distribution random. Shannon entropy ~7-8 bits/symbol. Nearly incompressible.
- **Autocorrelation**: Very few significant autocorrelations beyond white noise expectation. Mean |autocorrelation| close to 1/sqrt(N).
- **Mod structure**: Some non-uniformity for mod 2 (reflecting even/odd prime structure) but most moduli show uniformity.
- **Consecutive differences**: dd = gap(n) - Delta_R_inv(n); essentially unpredictable. Most common values are small even numbers.
- **Linear prediction**: AR(1) through AR(20) all achieve R^2 near 0 -- negligible variance reduction. Delta is NOT linearly predictable from its past.
- **Sign patterns**: ~50/50 positive/negative split; mean run length ~2.0 (consistent with iid).

## Verdict
**CLOSED** -- Failure Mode: **I** (Information Loss)

delta(n) appears pseudorandom with magnitude ~n^{0.5}. No exploitable linear or short-range structure detected at N=100,000. The correction term encodes contributions of O(sqrt(p(n))) Riemann zeros with random phases, consistent with the information-theoretic barrier.

## One-Line Summary
At N=100,000, delta(n) is pseudorandom (incompressible, uncorrelated, linearly unpredictable) with |delta| ~ n^{0.5}, confirming the zeta-zero oscillation barrier to exact prime computation.

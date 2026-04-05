# Sublinear Correction: Results

**Script:** sublinear_correction.py

## What Was Tested
Whether partial compressibility of delta(n) = pi(n) - R(n) (DCT captures 99% energy in ~10% of coefficients) can be exploited for a hybrid algorithm: compute R(x) in O(polylog), compute "spectrally dominant" correction in O(x^alpha) for alpha < 2/3, round to get pi(x).

## Key Findings
- DCT of delta(n): 99% energy in ~10.4% of coefficients for N = 10^4 -- confirms partial compressibility
- However, identifying WHICH coefficients matter requires computing delta(n) first -- circular
- The dominant DCT coefficients change location as N grows -- no fixed sparse support
- Combining DCT compression with Lagarias-Odlyzko analytic method: the analytic method already achieves O(x^{1/2+epsilon}), and DCT doesn't reduce this
- The Deléglise-Rivat O(x^{2/3}/ln^2(x)) combinatorial method is not improved by spectral structure
- Conclusion: compressibility is a post-hoc observation, not a computational shortcut -- you need the answer to find the compression

## Verdict
**CLOSED**
**Failure Mode:** Circularity (C) -- identifying the dominant spectral components of delta(n) requires computing delta(n) first; compression is retrospective, not predictive.

## One-Line Summary
Sublinear correction via DCT compression: 10% energy concentration exists but identifying dominant coefficients requires delta(n) first -- circular.

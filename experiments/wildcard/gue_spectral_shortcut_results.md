# GUE Spectral Shortcut: Results

**Script:** gue_spectral_shortcut.py

## What Was Tested
Whether the GUE (Gaussian Unitary Ensemble) determinantal structure of zeta zero spacings enables computing the oscillatory sum S(x) = sum_j x^{rho_j}/rho_j without evaluating each term individually. Tests GUE CLT (Central Limit Theorem) predictions against actual zeta zero sums.

## Key Findings
- GUE CLT predicts variance O(log N) for linear statistics of N eigenvalues -- but this applies to RANDOM GUE matrices, not the fixed zeta zero configuration
- For actual zeta zeros: the oscillatory sum has specific values determined by the exact zero locations, not just their statistical distribution
- The CLT gives the correct ORDER of fluctuations but not the exact value -- error is O(sqrt(log x)), far from the < 0.5 needed for exact rounding
- Determinantal structure (sine kernel) describes correlations but does not provide a computational shortcut for summing individual terms
- The GUE model is descriptive (explains statistics) not prescriptive (cannot replace exact computation)

## Verdict
**CLOSED**
**Failure Mode:** Information Loss (I) -- GUE statistics capture the distribution of the zero sum but not its exact value; O(sqrt(log x)) error remains.

## One-Line Summary
GUE spectral shortcut: CLT gives correct fluctuation order O(sqrt(log x)) but not exact zero-sum value -- statistics != computation.

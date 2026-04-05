# SVD Spectral Decay of pi(x) Communication Matrix: Results

**Script:** `svd_spectral_decay.py`
**Session:** 17

## What Was Tested
Distribution of singular values of the pi(x) communication matrix beyond the top-2 (which capture >99.99% of energy). Whether the oscillatory singular values decay geometrically (enabling low-rank approximation) or are flat (random-matrix-like). Also compared eigenvalue distribution of residual with GUE predictions.

## Key Findings
- Top-2 singular values capture the smooth R(x) contribution (>99.99% of Frobenius norm)
- Remaining singular values are flat (not geometrically decaying) -- random-matrix-like
- Eigenvalue distribution of the oscillatory residual matches GUE semicircle law
- The flat spectrum means no finite-rank approximation captures the oscillatory part
- This confirms the oscillatory contribution requires O(2^{N/2}) singular values

## Verdict
**CLOSED**
**Failure Mode:** Information loss (oscillatory singular values are flat/GUE; no low-rank approximation possible)

## One-Line Summary
Oscillatory singular values of pi(x) communication matrix are flat (GUE-like); no low-rank approximation beyond the smooth part.

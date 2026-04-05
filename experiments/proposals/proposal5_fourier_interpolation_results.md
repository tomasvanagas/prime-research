# Proposal 5: Fourier Interpolation of the Oscillatory Part -- Results

**Script:** proposal5_fourier_interpolation.py

## What Was Tested
View the oscillatory correction S(x) = pi(x) - R(x) as a quasi-periodic function in log-space with "frequencies" gamma_k (zeta zero imaginary parts). Sample S at known points (costing O(x^{2/3}) each), fit Fourier coefficients, then predict S at new points in O(K) time.

## Key Findings
- Sampling S(x) at training points works: can compute S = pi(x) - R(x) exactly for small x via primepi.
- Fitting K coefficients to K+1 samples gives a linear system that can be solved.
- Prediction on nearby test points: moderate accuracy for small x (100-10000 range).
- Critical failure: the number of significant zeros K is NOT polylog(x). For error < 1 at height x, need K ~ sqrt(x) * polylog(x) zeros.
- Even with "sparse Fourier" tricks, the signal is NOT effectively sparse -- GUE statistics mean all zeros contribute comparably, no dominant subset.
- Preprocessing cost: O(K * x_max^{2/3}) per sample point, and K ~ sqrt(x) zeros needed, so total is worse than direct computation.

## Verdict
**CLOSED**
**Failure Mode:** Information Loss (I) -- the oscillatory sum requires sqrt(x) * polylog zeros for unit accuracy; not effectively sparse.

## One-Line Summary
Fourier interpolation of S(x): requires sqrt(x) zero-frequencies, not polylog; GUE statistics prevent effective sparsity.

# Proposal 9: Sparse Fourier Transform on the Explicit Formula -- Results

**Script:** proposal9_sparse_fourier_explicit.py

## What Was Tested
Apply sparse Fourier transform (Hassanieh et al. 2012) to the zero-sum oscillatory part of the explicit formula. If the "effective number of zeros" contributing at a given x is polylog(x), sparse recovery could identify dominant contributions in O(polylog) time.

## Key Findings
- The explicit formula oscillatory sum S(x) = sum_rho x^rho / (rho * ln(x)) is a trigonometric sum with frequencies gamma_k * ln(x) / (2*pi).
- Sparse Fourier recovery requires the signal to be k-sparse with k = O(polylog). But the zero contributions have comparable magnitudes: each ~ x^{1/2} / |gamma_k|, summing to sqrt(x) overall.
- The GUE-random spacing of zeros means no dominant subset -- ALL zeros with gamma up to ~sqrt(x) contribute non-negligibly.
- Tested with 200 precomputed zeta zeros: truncation at K zeros leaves residual error that only decays as ~sqrt(x)/K, not exponentially.
- For unit accuracy, need K ~ sqrt(x) zeros -- the signal is NOT effectively sparse.
- Sparse Fourier techniques assume sparsity; when the sparsity condition fails, they give no advantage over direct computation.

## Verdict
**CLOSED**
**Failure Mode:** Information Loss (I) -- the zero-sum is not effectively sparse; all sqrt(x) zeros contribute comparably.

## One-Line Summary
Sparse Fourier on explicit formula: zero-sum is not k-sparse for k = polylog; all sqrt(x) zeros contribute, no recovery shortcut.

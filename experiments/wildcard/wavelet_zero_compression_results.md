# Wavelet Compression of Zeta Zero Sum — Results

**Experiment:** Is the oscillatory correction C(x) = pi(x) - li(x) compressible in wavelet or Fourier basis?

**Verdict:** CLOSED — Partially compressible (N^0.75 coefficients), far from polylog.

## Key Findings

1. **Raw signal appears ultra-sparse** — 99% of energy in 1 coefficient. But this is an artifact: the DC (mean) component dominates because |mean(C)| >> std(C).

2. **Detrended oscillatory part is NOT sparse**:
   - After removing mean: 99.9% energy requires ~47% of coefficients
   - After removing linear trend: ~50% of coefficients
   - After removing degree-5 polynomial: ~53-70% of coefficients
   - The more smooth structure you remove, the MORE random the residual looks

3. **Sparsity scaling (the critical test)**: For the detrended correction at x ~ 100000:
   - Wavelet 99.9% count ~ N^0.758
   - Fourier 99.9% count ~ N^0.734
   - Both are sublinear but FAR from O(polylog) = O((log N)^c)

4. **Power spectrum consistently f^{-1.7}** across all scales (1000 to 500000). This corresponds to Hurst exponent H ≈ 0.35 (anti-persistent process), consistent with oscillatory prime correction.

## Numbers
| Scale | N | Raw 99% | Detrended 99.9% wavelet | Detrended 99.9% fourier |
|---|---|---|---|---|
| 10000 | 2048 | 2 (0.1%) | 877 (42.8%) | 907 (44.3%) |
| 100000 | 2048 | 1 (0.0%) | 973 (47.5%) | 1109 (54.2%) |
| 500000 | 2048 | 1 (0.0%) | 699 (34.1%) | 866 (42.3%) |

## Scaling Law
99.9% energy coefficient count ~ N^{0.75} (average of wavelet and Fourier)

For x = 10^100 with N ~ x^{1/2} = 10^50 correction terms:
- Would need ~ (10^50)^{0.75} = 10^{37.5} wavelet coefficients
- Still enormous, though better than 10^50 (25% compression)

## Failure Mode
**INFORMATION**: The oscillatory correction encodes zeta zero phases which are information-theoretically dense. The 1/f^{1.7} spectrum means partial compressibility exists but not enough for polylog. The N^{0.75} scaling matches theoretical expectations: the correction has Kolmogorov complexity O(N^{1-epsilon}) for small epsilon.

**Session:** 40 (Fresh Perspective)

# Proposal 22: Compressed Sensing Recovery of delta(n) — Results

## Idea
Model delta(n) = p(n) - R^{-1}(n) as a signal, decompose it in a "zeta zero basis" {cos(gamma_j * log(R^{-1}(n)))}, and use compressed sensing (L1 minimization) to recover a sparse representation from subsampled measurements.

## Key Results

### Sparsity analysis
delta(n) is NOT sparse in the zeta zero basis:
| Basis size | Residual std | Significant coeffs | Sparsity ratio |
|-----------|-------------|-------------------|---------------|
| 10 funcs | 20.69 | 13 (of 20) | 65% dense |
| 20 funcs | 18.42 | 24 (of 40) | 60% dense |
| 50 funcs | 13.96 | 69 (of 100) | 69% dense |
| 100 funcs | 11.47 | 94 (of 200) | 47% dense |

delta(n) std = 33.29, so even 100 basis functions only capture ~66% of the variance.

### CS recovery from subsampling
Subsampling degrades gracefully — even 10% of samples gives similar train error. But **test error is always ~30**, meaning NO generalization beyond the training range.

### Extrapolation failure
Training on n=10..5000, testing on n=5001..10000:
- Correct prime identification rate: **11-16%** (essentially random)
- The learned representation does NOT extrapolate

### Fourier sparsity
delta(n) in standard Fourier basis:
- Top 5 coefficients: 38.8% energy
- Top 50 coefficients: 87.2% energy  
- Top 100 coefficients: 92.8% energy

This is moderate concentration, not sparse enough for CS guarantees.

## Verdict: CLOSED
delta(n) is not sparse in any tested basis (zeta zeros, Fourier). It requires O(N) coefficients to represent accurately, consistent with the information-theoretic barrier. CS cannot bypass the fundamental incompressibility of the oscillatory part.

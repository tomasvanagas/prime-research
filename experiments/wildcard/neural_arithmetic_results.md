# Neural Arithmetic: Results

**Script:** neural_arithmetic.py

## What Was Tested
Whether the correction delta(n) = pi(n) - R(n) has a compact neural network representation: (1) small network generalization on delta(n), (2) random Fourier feature projection for effective dimension, (3) NTK eigenvalue alignment, plus Kolmogorov complexity estimation via zlib compression.

## Key Findings
- Neural network (numpy-only, small MLP): fits training data but does NOT generalize to held-out n -- overfits completely
- Random Fourier features: need O(N) features to capture delta(n) -- no low effective dimension detected
- Compression test: zlib compression ratio of delta(n) sequence is near 1.0 (incompressible), confirming high Kolmogorov complexity
- This is consistent with delta(n) encoding ~50% of the bits of p(n), which are information-theoretically tied to zeta zeros
- TC^0 circuit representation would require polylog width, but delta has no detected threshold-circuit-friendly structure

## Verdict
**CLOSED**
**Failure Mode:** Information Loss (I) -- delta(n) is incompressible (zlib ratio ~1.0), does not generalize under neural approximation, and has no low-dimensional structure.

## One-Line Summary
Neural arithmetic on delta(n): no generalization, incompressible (zlib ~1.0), O(N) random features needed -- no compact representation.

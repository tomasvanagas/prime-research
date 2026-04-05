# Neural Architecture Experiments: Results

**Script:** neural_prime.py

## What Was Tested
Five neural architectures implemented from scratch in numpy (no PyTorch/TF) to learn delta(n) = p(n) - round(R_inv(n)): (1) Transformer on base-10 digit encoding, (2) Number-theoretic features + MLP, (3) Fourier Neural Operator (FNO), (4) Symbolic regression via genetic programming, (5) Kolmogorov-Arnold Network (KAN).

## Key Findings
- Transformer on digits: learns digit distribution statistics but cannot predict individual corrections; test accuracy ~random
- Number-theoretic MLP: features like n mod q, Omega(n), mu(n) have near-zero correlation with delta(n)
- FNO: captures spectral structure of delta(n) on training range but fails to extrapolate (generalization ratio >3x)
- GP symbolic regression: best formulas are O(sqrt(n)/ln(n)) scaling with random-looking residual
- KAN: approximates smooth component well but delta(n) correction is not in any finite KAN basis
- All architectures converge to the same ~50% digit accuracy as R_inv(n) -- the smooth approximation ceiling

## Verdict
**CLOSED**
**Failure Mode:** I (Information Loss -- all neural architectures learn only the smooth PNT approximation; the oscillatory correction from zeta zeros is unlearnable from finite training data)

## One-Line Summary
Five numpy-based neural architectures (Transformer, MLP, FNO, GP, KAN) all hit the smooth approximation ceiling; delta(n) is unlearnable.

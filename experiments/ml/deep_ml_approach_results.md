# Deep ML Approach: Results

**Script:** deep_ml_approach.py

## What Was Tested
Six ML methods for predicting p(n) using sklearn/numpy (no PyTorch): MLP on log-space, residue network (p(n) mod q patterns), autoregressive gap prediction, custom symbolic regression, neural ODE approximation via scipy, and correction model for delta(n) = p(n) - R_inv(n). Trained on first 100K-200K primes.

## Key Findings
- MLP regressor on log-space features achieves low relative error (~0.01%) but near-zero exact matches on test set
- Residue patterns (p(n) mod small q) are predictable only to the extent of Dirichlet density -- no exploitable structure beyond PNT
- Gap prediction with autoregressive features shows no generalization beyond mean gap ~ln(n)
- Custom symbolic regression finds formulas equivalent to known asymptotic expansions (n*ln(n) + n*ln(ln(n)) - n + ...)
- Neural ODE approximation converges to the smooth PNT trajectory, missing the oscillatory correction
- Correction model for delta(n) overfits on train, fails on test -- delta is effectively random

## Verdict
**CLOSED**
**Failure Mode:** I (Information Loss -- ML models cannot capture the ~50% of digits encoded in zeta zero oscillations; they learn only the smooth part)

## One-Line Summary
MLP/GBR/RF/symbolic regression on prime features learn only the PNT smooth approximation; exact correction delta(n) is unlearnable (random-like).

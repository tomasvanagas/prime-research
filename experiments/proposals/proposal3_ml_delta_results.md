# Proposal 3: ML Surrogate Model for delta(n) -- Results

**Script:** proposal3_ml_delta.py

## What Was Tested
Train a machine learning model on features of n (residues mod small primes, fractional parts of n * gamma_k for zeta zeros, scale features) to predict delta(n) = p(n) - round(R^{-1}(n)). Test whether delta(n) has learnable structure that generalizes.

## Key Findings
- Features include: log(n), sqrt(n), n^{1/3}, residues mod first 10 primes, fractional parts of n * gamma_k for 20 zeta zeros.
- Linear models capture essentially zero variance of delta(n) -- the residual is noise-like.
- The zeta-zero-based features (fractional parts of n * gamma_k) are the most informative, but require knowing the zeros themselves.
- delta(n) distribution is approximately Gaussian (CLT over many zero contributions), consistent with incompressible randomness.
- Even partial prediction (to within +/- C) would require capturing the cumulative effect of ~sqrt(x) zeta zeros -- the information content scales with sqrt(x), not polylog(x).

## Verdict
**CLOSED**
**Failure Mode:** Information Loss (I) -- delta(n) encodes GUE-random zeta zero phases; ML cannot compress what is information-theoretically incompressible.

## One-Line Summary
ML surrogate for delta(n): no learnable structure beyond trivial features; residual is GUE-random noise.

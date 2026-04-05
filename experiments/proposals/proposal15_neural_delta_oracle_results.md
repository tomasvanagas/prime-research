# Proposal 15: Neural Delta Oracle + Certification -- Results

**Script:** proposal15_neural_delta_oracle.py

## What Was Tested
Train a neural network to predict delta(n) = p(n) - R^{-1}(n), then verify the candidate via primality testing and local prime counting. Also explores the "bootstrapping" idea: if p(n-1) is known, find p(n) by searching the next O(log^2 x) numbers.

## Key Findings
- Features: n, log(n), log(log(n)), fractional parts of n * gamma_k for zeta zeros.
- Neural prediction accuracy: for small n (training range), moderate fit. For extrapolation to larger n, prediction degrades rapidly -- delta(n) is non-stationary and its statistics change with n.
- The verification step is the bottleneck: even if the neural net predicts delta(n) to within the prime gap, verifying pi(candidate) = n requires counting primes, which is O(x^{2/3}).
- Bootstrapping from p(n-1): finding next prime costs O(log^2(x) * polylog) via AKS, but bootstrapping from p(1) to p(n) costs O(n) steps -- not polylog in n.
- For bootstrapping to help, need a "starting point" p(n_0) with n_0 close to target n. Computing any single p(n_0) is the original hard problem.
- The catch-22: prediction requires generalization that delta(n) doesn't support; verification requires pi(x) computation that costs O(x^{2/3}).

## Verdict
**CLOSED**
**Failure Mode:** Circularity (C) -- prediction doesn't generalize; verification requires the hard computation we're trying to avoid.

## One-Line Summary
Neural delta oracle: prediction doesn't generalize to new n; verification via pi(x) costs O(x^{2/3}), defeating the purpose.

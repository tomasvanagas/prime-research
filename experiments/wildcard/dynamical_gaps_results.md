# Dynamical Gaps: Results

**Script:** dynamical_gaps.py

## What Was Tested
Deep analysis of prime gaps as a dynamical system: AR(k) models on raw and log-normalized gaps, residue-conditional maps g(n+1) = f(g(n), p(n) mod m), hidden Markov models on quantized gaps, entropy/conditional entropy/mutual information, modular structure, cumulative sum complexity, and substitution/morphism approximations.

## Key Findings
- AR(k) for k up to 50: R^2 < 0.01, gaps are linearly unpredictable from history
- MI(g_n; g_{n+1}) ~ 0.38 bits out of ~3.7 bits entropy -- only 10.3% mutual information
- Residue-conditional maps: conditioning on p(n) mod m provides negligible additional predictability
- HMM on quantized gaps: no hidden states improve prediction beyond marginal distribution
- Gap sequence is only ~9% more compressible than i.i.d. with same marginal -- near Cramer random model
- Substitution approximations fail: block complexity grows linearly with block length

## Verdict
**CLOSED**
**Failure Mode:** Information Loss (I) -- prime gaps are near-random with only 10% mutual information between consecutive gaps; no dynamical structure to exploit.

## One-Line Summary
Dynamical gap analysis (AR, HMM, MI, morphisms) on 10^5 primes: MI ~0.38 bits (10%), 9% compressible -- near Cramer random model.

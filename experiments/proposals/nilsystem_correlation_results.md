# Green-Tao Nilsystem Correlation for Prime Prediction -- Results

**Script:** nilsystem_correlation.py

## What Was Tested
Fit nilsequences (1-step: periodic/Fourier functions; 2-step: bracket quadratics) to the prime indicator 1_P(n) and measure residual structure after removing nilsequence correlations. Key question: can nilsequence-based prediction beat R^{-1}(n)?

## Key Findings
- 1-step nilsequences (Fourier modes): best individual frequencies capture the mod-q biases (primes avoid 0 mod q). Total variance explained is small -- density 1/ln(n) accounts for most, residual is noise-like.
- 2-step nilsequences (bracket quadratics n -> F(alpha*n, beta*n^2)): capture quadratic patterns in prime distribution, but explanatory power is marginal beyond 1-step.
- After removing nilsequence correlations, the residual has the same information-theoretic complexity as before -- the nilsequences only capture the smooth/structured part that R(x) already captures.
- Green-Tao theory guarantees primes CORRELATE with nilsequences (inverse theorem for Gowers norms), but the correlation is with the DENSITY, not with individual prime positions.
- Nilsequence prediction accuracy does NOT beat R^{-1}(n) for locating specific primes. Both capture the smooth part; neither captures the oscillatory/random part.
- The residual after all nilsequence corrections still has O(log n) bits of entropy per prime.

## Verdict
**CLOSED**
**Failure Mode:** Equivalence (E) -- nilsequences capture the same smooth structure as R(x); residual is the same oscillatory/random component.

## One-Line Summary
Nilsystem correlations: capture smooth density (same as R(x)); residual after removal has full O(log n) bits of entropy, no prediction gain.

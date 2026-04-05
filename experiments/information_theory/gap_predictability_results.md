# Gap-Based Predictability of pi(x) Results

## What Was Tested
Analysis of whether prime gap sequences have exploitable predictive structure, using primes up to 200,000:
1. **Gap statistics**: Mean, std, CV, max gap, distribution mod 6.
2. **Autoregressive prediction**: AR(p) models for p = 1,2,3,5,10,20,50 on normalized gaps g/ln(p), with train/test split.
3. **Autocorrelation of gaps**: At lags 1-100, with significance testing.
4. **Conditional gap distribution**: P(g_{n+1} | g_n) for specific previous gap values.
5. **Information-theoretic analysis**: gzip/bz2/lzma compression of gap bytes vs iid random with same distribution; Shannon entropy; conditional entropy H(g_{n+1}|g_n); mutual information I(g_n; g_{n+1}).
6. **Cramer model comparison**: KS test of normalized gaps vs exponential distribution.

## Key Findings
- **AR prediction**: Improvement over baseline is less than 5% at ALL orders tested. Gaps are essentially unpredictable from their history.
- **Autocorrelation**: Weak but statistically significant at small lags (Hardy-Littlewood effect). No useful long-range structure.
- **Conditional distribution**: After gap=2, the next gap has slightly different mean, but the effect is tiny.
- **Compression**: Gaps are only marginally more compressible than iid random with same distribution (~few percent).
- **Shannon entropy**: Near maximum for the number of distinct gap values.
- **Mutual information**: Knowing g_n reduces uncertainty about g_{n+1} by only ~2-5% -- far too small for sublinear counting.
- **Cramer model**: Normalized gaps pass KS test reasonably well -- consistent with Cramer random model.
- **Conditional entropy**: H(g_{n+1}|g_n) is close to marginal entropy H(g) -- minimal information gain from conditioning.

## Verdict
**CLOSED** -- Failure Mode: **I** (Information Loss)

Prime gaps have mild short-range correlations (Hardy-Littlewood) but no useful long-range structure. The mutual information between consecutive gaps is far too small (~2-5%) to enable any sublinear counting shortcut. The gap sequence is information-theoretically close to the Cramer random model.

## One-Line Summary
AR models, autocorrelation, compression, and mutual information analysis all confirm prime gaps are near-iid with only ~2-5% conditional information gain -- far insufficient for a counting shortcut.

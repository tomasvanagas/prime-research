# Proposal 21: Zero Clustering Truncation — Results

## Idea
Exploit GUE statistics of zeta zero spacings to group zeros into clusters, representing each cluster by a single term with a coherence correction factor. Reduces the number of evaluations from O(sqrt(x)) individual zeros to fewer clusters.

## Key Results

### Accuracy vs number of zeros
Even 5 zeros bring the explicit formula correction to within ~1 of the true pi(x), compared to ~14 mean error for R^{-1} alone. But more zeros don't consistently improve — the 200-zero correction has HIGHER mean error (1.479) than 10 zeros (1.262).

This suggests the explicit formula with only 200 zeros is NOT converging — it oscillates around the true value, consistent with the O(sqrt(x)) zeros needed theory.

### Clustering results
- Radius 1.0: 175 clusters from 200 zeros — minimal compression
- Radius 2.0: 117 clusters — errors comparable to individual
- Radius 5.0: 61 clusters — errors sometimes WORSE (coherence loss)
- **Radius 10.0: 35 clusters — surprisingly BETTER errors** (0.025, 0.054, 0.009)

The radius-10 result is anomalous — likely due to cancellations within clusters happening to align well for these specific test points. Not generalizable.

### Tail prediction
The tail (contribution of zeros beyond a cutoff) is NOT predictable from low zeros. The ratio tail/partial varies wildly (-3.4 to +3.3), confirming that each zero carries independent information.

### Error statistics (n=10..5000)
| Method | Mean Error | Median | Max | Std |
|--------|-----------|--------|-----|-----|
| R^{-1} only | 14.020 | 11.018 | 66.376 | 11.790 |
| 10 zeros | 1.262 | 1.007 | 3.753 | 0.934 |
| 50 zeros | 1.411 | 1.131 | 5.448 | 1.142 |
| 200 zeros | 1.479 | 1.128 | 5.670 | 1.218 |

## Verdict: CLOSED
Clustering reduces constants but not asymptotics. The tail contribution from unseen zeros is unpredictable, confirming that O(sqrt(x)) zeros are fundamentally needed. The approach cannot achieve O(polylog).

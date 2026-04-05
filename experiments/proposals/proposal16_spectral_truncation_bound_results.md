# Proposal 16: Spectral Truncation with Adaptive Zero Selection — Results

## Idea
Instead of summing over ALL zeta zeros in the explicit formula, select a small
O(polylog(n)) subset of "resonant" zeros (those whose gamma * log(x) is near
a multiple of pi), then compute only their contributions.

## Results

### R^{-1}(n) quality (with fixed Gram series)
| n | p(n) | R_inv(n) | delta(n) |
|---|------|----------|----------|
| 10 | 29 | 29.34 | -0.34 |
| 100 | 541 | 536.48 | 4.52 |
| 1000 | 7919 | 7922.57 | -3.57 |
| 5000 | 48611 | 48554.95 | 56.05 |
| 10000 | 104729 | 104767.79 | -38.79 |

### Zeros needed for exact p(n) via adaptive selection
| n | Zeros needed | log(n) | log^2(n) |
|---|-------------|--------|----------|
| 10 | 1 | 2.30 | 5.30 |
| 50 | 2 | 3.91 | 15.30 |
| 100 | 2 | 4.61 | 21.21 |
| 500 | >1000 | 6.21 | 38.62 |
| 1000 | >1000 | 6.91 | 47.72 |
| 2000 | 24 | 7.60 | 57.77 |
| 5000 | >1000 | 8.52 | 72.54 |
| 10000 | >1000 | 9.21 | 84.83 |

### Key findings
- Adaptive selection wins over sequential in ~41% of cases (not consistently better)
- For small n (< 100), 1-2 zeros can give exact answers
- For n >= 500, 1000 zeros are typically insufficient
- The leading-term approximation for zero contributions is too crude
- Even with correct multi-term expansions, the number of zeros needed grows
  as O(sqrt(p(n))) because the sum over zeros converges as inverse-square-root

## Verdict: CLOSED
The spectral truncation idea fails because:
1. The leading-term zero contribution formula is too inaccurate for large x
2. "Resonant" zeros don't form a predictable small subset — which zeros matter
   depends on x in an essentially unpredictable way
3. The number of significant zeros grows as O(x^{1/2}) by the explicit formula's
   error term, matching the Lagarias-Odlyzko bound
4. Adaptive selection provides marginal improvement over sequential but cannot
   reduce O(x^{1/2}) zeros to O(polylog)

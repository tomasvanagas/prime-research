# Proposal 21: Zero Clustering Truncation -- Results

## What was tested

Whether grouping nearby Riemann zeta zeros into clusters (exploiting GUE spacing statistics) can reduce the number of explicit-formula terms from O(sqrt(x)) to O(polylog(x)). The idea: represent each cluster by a single representative zero weighted by cluster size, with a sinc-like coherence correction for phase spread. Four sub-experiments were run using 200 precomputed zeta zeros.

## Key numerical results

### Test 1: Accuracy vs number of zeros (corrections to pi(x) estimate)

| n | p(n) | R_inv | 5 zeros | 10 zeros | 20 zeros | 50 zeros | 100 zeros | 200 zeros |
|---|------|-------|---------|----------|----------|----------|-----------|-----------|
| 10 | 29 | 29.13 | -0.09 | -0.09 | -0.06 | 0.06 | 0.15 | 0.34 |
| 50 | 229 | 230.32 | -0.60 | -0.63 | -0.60 | -0.24 | -0.25 | -0.18 |
| 100 | 541 | 536.40 | -0.54 | -0.37 | -0.53 | -0.44 | -0.75 | -0.72 |
| 500 | 3571 | 3575.05 | -0.16 | -0.25 | -0.35 | 0.31 | 0.06 | 0.18 |
| 1000 | 7919 | 7922.70 | -1.11 | -0.50 | 0.15 | -0.26 | 0.14 | 0.64 |

Adding more zeros does NOT monotonically improve the correction. All corrections remain O(1) regardless of zero count.

### Test 2: Clustered vs individual zeros

| Radius | Clusters | n=100 err (ind/cls) | n=500 err (ind/cls) | n=1000 err (ind/cls) |
|--------|----------|---------------------|---------------------|----------------------|
| 1.0 | 175 | 0.717 / 0.594 | 0.180 / 0.255 | 0.640 / 0.678 |
| 2.0 | 117 | 0.717 / 0.657 | 0.180 / 0.295 | 0.640 / 0.056 |
| 5.0 | 61 | 0.717 / 0.305 | 0.180 / 0.672 | 0.640 / 1.515 |
| 10.0 | 35 | 0.717 / 0.025 | 0.180 / 0.054 | 0.640 / 0.009 |

Clustering sometimes produces lower errors than individual zeros (especially radius=10), but this is accidental cancellation, not systematic improvement. The errors are erratic and non-monotonic across n values.

### Test 3: Tail prediction via GUE statistics

The tail-to-partial ratio is wildly unstable:
- At n=100: ratios range from 0.36 to 0.94
- At n=500: ratios range from -1.73 to -0.42
- At n=1000: ratios range from -3.42 to 3.33

The tail (contribution of high zeros) is **not predictable** from the low zeros. Sign flips and magnitude swings indicate the tail contributions carry independent information.

### Test 4: Error distribution (n=10..5000)

| Method | Mean err | Median err | Max err | Std |
|--------|----------|------------|---------|-----|
| R^{-1} only | 14.020 | 11.018 | 66.376 | 11.790 |
| 10 zeros | 1.262 | 1.007 | 3.753 | 0.934 |
| 50 zeros | 1.411 | 1.131 | 5.448 | 1.142 |
| 200 zeros | 1.479 | 1.128 | 5.670 | 1.218 |

Critical finding: 10 zeros performs as well as 200 zeros (mean error 1.26 vs 1.48). Adding zeros beyond ~10 actually increases error slightly due to accumulated floating-point oscillations. The correction saturates at O(1) error with very few zeros.

## Verdict: FAIL

## Failure mode: I -- Information Loss

The clustering approach reduces the number of terms but does NOT reduce the residual error below O(1). The fundamental problem: with only 200 zeros, the explicit formula already cannot push below O(1) error for the range tested. Clustering preserves roughly the same O(1) error with fewer terms, but this is irrelevant -- the bottleneck is not the number of low-lying zeros but the tail of the zero sum (the ~10^48 zeros needed at large x). The tail contributions are informationally independent and unpredictable from low zeros, as Test 3 demonstrates with erratic, sign-flipping tail/partial ratios.

## Conclusion

Zero clustering can reduce constants in the explicit formula computation but cannot change the asymptotic zero count required. The experiment confirms that even a modest O(1) correction from zeros saturates quickly (~10 zeros suffice for the tested range), and the remaining error is dominated by the infinite tail of higher zeros whose contributions are mutually incoherent and carry independent phase information. At scale (x ~ 10^100), this tail requires O(x^{1/2}) zeros -- clustering into O(polylog) representatives would discard the phase information encoded in individual zero positions, which is exactly the information needed to resolve pi(x) to integer accuracy. This is another manifestation of the core barrier: the oscillatory part of pi(x) is information-theoretically incompressible.

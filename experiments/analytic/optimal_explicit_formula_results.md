# Optimal Explicit Formula Configuration: Results

**Script:** optimal_explicit_formula.py

## What Was Tested
Systematic optimization of the explicit formula pi(x) = R(x) - sum_rho R(x^rho): testing all combinations of (1) number of zeros K, (2) smoothing kernels (none, Gaussian, Riesz, Cesaro, Fejer), (3) precision per zero, (4) conjugate pair handling.

## Key Findings
- Raw summation: converges as K^{-0.01} (barely converges) due to oscillatory terms
- Cesaro kernel (1-k/K): error ~1/K, moderate improvement
- Riesz kernel (1-(k/K)^2): slightly better than Cesaro
- Fejer kernel (1-k/K)^2: best among polynomial kernels, error ~1/K^2
- Gaussian kernel exp(-(k/K)^2 * sigma): best overall with sigma~2, error drops exponentially for low K then plateaus
- Optimal configuration: Gaussian kernel with sigma=2 gives ~10x improvement over raw sum
- BUT: all kernels converge to the same asymptotic regime O(sqrt(x)/K) for large K
- Precision per zero: 30 digits sufficient for x < 10^12; higher precision needed for larger x
- Conjugate pairs must be summed together to maintain real arithmetic

## Verdict
**CLOSED**
**Failure Mode:** E (Equivalence -- optimal kernel gives constant-factor improvement; O(sqrt(x)) zeros still required for exactness)

## One-Line Summary
Optimal explicit formula: Gaussian kernel with sigma~2 gives ~10x improvement; asymptotic O(sqrt(x)) zero scaling unchanged.

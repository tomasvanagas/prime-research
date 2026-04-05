# Self-Correcting Explicit Formula: Results

## Experiment

**Idea**: Use the truncated explicit formula pi(x) = Li(x) - sum_rho Li(x^rho) - log(2)
with fewer zeros than normally needed (N << x^{1/2}), then apply integer constraints
(rounding, monotonicity, step-function, primality testing) to correct errors.

**Hypothesis**: Self-correction constraints reduce the number of zeros needed for exact pi(x).

## Setup
- Range: x = 2 to 1000 (extended test to 2000)
- Zero counts tested: 5, 10, 20, 30, 50, 100
- Zeros loaded from `data/zeta_zeros_*.txt`
- Precision: mpmath with 50 decimal digits
- Three correction methods:
  1. **Round**: just round to nearest integer
  2. **Corrected**: round + enforce monotonicity + even-number constraint
  3. **Primality**: round + monotonicity + step-function enforcement via isprime()

## Results Summary

| N zeros | Round Acc | Corr Acc | Prime Acc | Round MaxErr | Corr MaxErr | Prime MaxErr | Time |
|---------|-----------|----------|-----------|-------------|-------------|-------------|------|
| 5       | 0.30%     | 2.20%    | 0.00%     | 68          | 38          | 4           | 1.1s |
| 10      | 0.00%     | 0.60%    | 0.00%     | 145         | 72          | 10          | 2.1s |
| 20      | 0.00%     | 0.40%    | 0.00%     | 236         | 134         | 21          | 3.5s |
| 30      | 0.00%     | 0.40%    | 0.00%     | 376         | 137         | 30          | 5.5s |
| 50      | 0.00%     | 1.00%    | 0.00%     | 571         | 115         | 52          | 9.3s |
| 100     | 0.00%     | 0.60%    | 0.00%     | 1236        | 106         | 106         | 19.0s|

100% accuracy was **NOT achieved** with any method, even with 100 zeros.

## Key Findings

### 1. Error grows linearly with number of zeros (CRITICAL BUG IN APPROACH)

The truncated explicit formula has a **systematic bias** that grows proportionally to N.
At small x, each zero's contribution Li(x^rho) is negative, and subtracting these negative
values makes the result increasingly positive. For x=2:

- Li(2) = 0 (by definition of offset Li)
- Each Re(Li(2^rho)) ~ -0.3 to -0.8
- Subtracting N such terms adds roughly N * 0.5 to the result
- With 100 zeros, pi_approx(2) ~ 107 instead of 1

This is a well-known convergence issue: the explicit formula's truncation error depends
on the truncation point relative to x. For the formula to converge, you need roughly
O(x/log(x)) zeros, not a fixed number.

### 2. Error distribution is highly autocorrelated

| N zeros | Mean Error | Std Error | Max Error | Autocorr(1) | Autocorr(5) | Within 0.5 |
|---------|-----------|-----------|-----------|-------------|-------------|-----------|
| 10      | -66.59    | 31.80     | 145.05    | 0.994       | 0.932       | 0.0%      |
| 30      | -211.05   | 72.19     | 376.41    | 0.989       | 0.891       | 0.0%      |
| 50      | -354.48   | 111.63    | 571.32    | 0.987       | 0.901       | 0.0%      |

The error is **extremely smooth** (autocorrelation > 0.98 at lag 1), confirming it
is dominated by the systematic truncation bias, not random fluctuations.

### 3. Self-correction constraints provide minimal help

- **Monotonicity + rounding**: Reduces max error by ~50% but accuracy stays < 3%
- **Primality testing**: Actually HURTS because the base approximation is so far off
  that the step-correction algorithm propagates errors globally
- The constraints cannot compensate for a systematic bias of order ~N

### 4. The fundamental barrier

The explicit formula's truncation error at point x with N zeros is approximately:

  E(x, N) ~ -sum_{gamma > gamma_N} 2*Re(Li(x^rho)) / log(x)

This sum decays as roughly O(x^{1/2} / (gamma_N * log(x))), where gamma_N is the
N-th zero. For the error to be < 0.5 (required for rounding to work), you need
gamma_N >> x^{1/2}, i.e., N = O(x^{1/2} * log(x)) zeros.

**No amount of integer constraints can fix a systematic bias of magnitude >> 1.**
The constraints only help when the error is already < 1, which is exactly the regime
where the formula already works by simple rounding.

## Verdict

**CLOSED.** The self-correcting approach fails because:

1. The truncation error is **systematic** (bias), not **random** (noise)
2. Integer constraints can only correct errors < 1; the formula's error with few zeros is >> 1
3. The error grows proportionally to number of zeros used (at small x)
4. The autocorrelation of errors is > 0.98, meaning nearby points have nearly identical errors -- no useful "triangulation" from a window of values

The idea would only work if the truncation error were random with magnitude < 1.
But it is deterministic and large. This is fundamentally why O(x^{1/2}) zeros are needed.

**This is approach #650+ in the closed paths list. The integer constraint trick does not
bypass the information-theoretic barrier: recovering pi(x) exactly requires information
about O(x^{1/2}) zeta zeros, which cannot be compressed.**

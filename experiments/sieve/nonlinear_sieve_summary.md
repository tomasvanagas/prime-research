# Nonlinear Sieve Experiments: Summary of Findings

**Date:** 2026-04-04 (Session 14)
**Experiments:** 6 scripts, 25+ sub-experiments

## Executive Summary

**ALL nonlinear sieve approaches tested FAIL to reduce complexity below O(x^{2/3}).**

Nonlinear operations (products, bitwise, thresholds) CAN break the Selberg parity barrier
(they distinguish primes from semiprimes) but do NOT reduce computational cost. The barrier
is not parity -- it is INFORMATION: O(sqrt(x)) distinct floor values exist, and no nonlinear
combination of fewer values can reconstruct pi(x).

## Detailed Results

### 1. Products of Floor Values (`nonlinear_floor_products.py`)

- **Correction term** delta(a,b,x) = floor(x/a)*floor(x/b) - floor(x/(a*b)) does NOT
  separate primes from composites. Value overlap is nearly complete.
- **Polynomial fit** (degree 2 in K floor values):
  - K=10, quadratic: max_err=1.34, 131/198 exact on training set
  - Generalization: max_err ~10, near 0% exact on test set
  - **No generalization.** The polynomial overfits the training range.
- **Correction term sums:** K=10 correction terms give max_err=2.85. Poor.
- **Information content:** x=10000 has only 199 distinct floor values (~2*sqrt(x)).
  Products expand this to ~6882, but still O(x) information is needed.

### 2. Comparisons and Thresholding (`nonlinear_comparisons.py`)

- **Parity of floor(x/k):** Correlation with primality is ~0 (0.009 at x=100).
  No signal in parity patterns.
- **Gap function** g(k,x) = floor(x/k) - floor(x/(k+1)):
  - Only O(sqrt(x)) values of k have nonzero gap
  - Most primes (96%+) are INVISIBLE to the gap function (gap=0 at prime k)
  - Only ~4-7% of primes above sqrt(x) appear as "jump points" floor(x/m)
- **Modular floor values** floor(x/k) mod m: Correlation with primality <0.12.
- **Bitwise operations:** XOR/AND/OR of floor values carry no useful prime signal.
  popcount_sum/x ~ 1.33 while pi(x)/x ~ 0.25 -- no connection.

### 3. Multiplicative Structure (`nonlinear_multiplicative_structure.py`)

- **GCD structure:** GCD(floor(x/k), k) does not distinguish primes from composites.
- **Floor value collisions:** Only ~30% of collision groups contain primes. Not a
  useful discriminator.
- **Quadratic sieve analog:** sum floor(x/p)^2 / sum floor(x/k)^2 ~ 0.71 across
  all tested x. Ratio is stable but encodes prime DENSITY not pi(x).
- **Nonlinear Mobius:** pi(x) as quadratic in M(x/k) values: max_err=14.6 on
  training, 44.4 on test. Complete failure.

### 4. Bitwise Operations / TC^0 (`nonlinear_bitwise_tc0.py`)

- **Bit decomposition:** Linear on 150 bit-features (10 bits x 15 floor values):
  max_err=2.17, 312/498 exact. Adding AND pairs: max_err=1.53, 428/498 exact.
  Still ~15% error rate and no generalization guarantee.
- **Polynomial degree vs floor values needed:**
  - K=2, deg=3: train max_err=1.58, test max_err=25.31 (CATASTROPHIC overfit)
  - K=12, deg=2: train max_err=1.76, test max_err=5.74 (poor generalization)
  - K=20, deg=2: train 285/290 exact, test max_err=10.49 (TOTAL overfit)
  - **Pattern: increasing K or degree improves training, worsens testing.**
- **Mutual information:**
  - I(floor(x/1); pi(x)) = 7.95 bits (floor(x/1) = x determines pi(x))
  - I(floor(x/2); pi(x)) = 7.65 bits (95% of info in x/2 alone)
  - PAIRS of floor values barely add info beyond the individual (joint info =
    max of individual infos). Products/XOR sometimes LOSE info.
  - **Critical insight:** floor(x/k) for small k are highly correlated.
    Diversity of information requires k up to O(sqrt(x)).

### 5. Selberg Parity Barrier (`nonlinear_selberg_parity.py`)

- **Parity barrier confirmed:** Legendre sieve (linear) correctly finds all primes
  above sqrt(x) but cannot COUNT them efficiently (2^{pi(sqrt(x))} terms).
- **Split sieve:** Running sieves on disjoint prime sets and multiplying outputs
  = conjunction = full Legendre sieve. No improvement.
- **Prime vs semiprime distinguishing:** Nonlinear ops CAN distinguish (semiprimes
  always have a small divisor, primes don't). But this requires O(sqrt(n))
  divisibility tests per number n, giving O(x^{3/2}) total -- WORSE.
- **Key conclusion:** Parity barrier is about CORRECTNESS of linear sieves.
  The EFFICIENCY barrier is different: it comes from the O(sqrt(x)) distinct
  floor values that parameterize the Dirichlet hyperbola.

### 6. Floor Value Identities (`nonlinear_floor_identities.py`)

- **Recursive identities:** pi(x) ~ c*pi(x/2) + d*pi(x/3) + linear floor terms:
  max_err ~2.6 with 3-5 recursive terms. The "correction" (residual error)
  is ~2-3, comparable to prime gaps -- cannot be made zero with few terms.
- **Integer coefficient search:** No simple rational formula exists. Best-fit
  coefficients are irrational with no pattern.
- **Higher-order differences:** Delta_1 floor(x/k) = [k|x] (divisor indicator).
  Products of differences give lcm-divisibility, not primality.
- **Key: overfitting scales with degree.**
  - K=20, deg=2: 231 features, train 285/290 exact, test max_err=10.49
  - The model memorizes the training range via coefficients but fails outside.

## Theoretical Analysis

### Why Nonlinear Sieves Cannot Help

1. **Floor value bottleneck:** There are exactly 2*floor(sqrt(x)) - 1 distinct
   values in {floor(x/k) : 1 <= k <= x}. ALL information about primes up to x
   is encoded in these O(sqrt(x)) values plus their multiplicities.

2. **No nonlinear compression:** A degree-d polynomial in K floor values has
   O(K^d) terms. For this to equal pi(x) EXACTLY for all x, we need K^d >= sqrt(x),
   giving K >= x^{1/(2d)}. Even degree 10 requires K >= x^{0.05}, which is
   NOT polylog.

3. **Generalization failure is structural:** The polynomial fits work on the training
   range because they have enough free parameters. But pi(x) is essentially a
   step function (jumps at primes), and polynomials in floor values are piecewise
   polynomial -- the pieces don't align outside the training range.

4. **Parity barrier vs efficiency barrier:** Nonlinear operations solve the
   parity problem (CAN distinguish primes from semiprimes). But the EFFICIENCY
   barrier is that evaluating the Legendre/Meissel-Lehmer formula requires
   visiting O(x^{2/3}) of the O(sqrt(x)) floor-value groups. Nonlinear operations
   don't reduce the number of groups that must be visited.

## Verdict

**CLOSED.** Nonlinear sieve approaches:
- Break parity barrier: YES (in principle)
- Reduce complexity below O(x^{2/3}): NO
- Mode: **E** (Equivalence) -- nonlinear combinations of floor values still
  require O(sqrt(x)) distinct values, equivalent to Meissel-Lehmer.

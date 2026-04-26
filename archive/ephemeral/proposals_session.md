# Proposal Session — 2026-04-26

Four concrete proposals for computing p(n) in O(polylog n) time, each with
a runnable test in `experiments/proposals/`.

---

## Proposal A: Pade–Borel Acceleration of the Explicit-Formula Zero Sum

### Idea
The Riemann explicit formula gives
```
pi(x) = R(x) - sum_rho R(x^rho) - log(2) + integral...
```
Truncating the zero sum at height T leaves error of order x/sqrt(T) * log(x),
so we apparently need T = sqrt(x) zeros for exact integer recovery.

**Hypothesis.** The partial sums S_T(x) = sum_{|rho|<=T} R(x^rho) form a
sequence whose tail has structured oscillation (driven by GUE pair correlation
of zeros). If a Pade or Wynn-epsilon resummation of S_T can be applied to a
short geometric ladder of partial sums (T_k = T_0 * 2^k for k=0..log log x),
we may recover S_inf from O(polylog) zeros instead of O(sqrt x).

The substantive bet: the *correlations* between zeros (Montgomery–Odlyzko
GUE) carry redundant information about the tail. A few zeros plus their
known correlation structure may suffice.

### Pseudocode
```
def pi_pade(x, num_zeros=20):
    rhos = load_zeros(num_zeros)
    Ts   = [rhos[2**k - 1].imag for k in range(log2(num_zeros))]
    S    = [partial_zero_sum(x, T) for T in Ts]
    S_inf = wynn_epsilon_extrapolate(S, Ts)
    return round(R(x) - S_inf - log(2))
```

### Complexity
- O(polylog(x)) zeros; each R(x^rho) is O(log x) with Riemann–Siegel arithmetic.
- Wynn-epsilon: O((log x)^3).

### Key assumption
The truncation tail S_inf - S_T(x) is *oscillatory in T with a smooth envelope*,
not just bounded by O(x/T^{1/2}). Pade/Wynn need the smooth envelope.

### Test design
For x = 100, 1000, 10000:
- Compute S_T for T at heights of zeros 5, 10, 20, 40, 80.
- Wynn-epsilon to extrapolate.
- Compare R(x) - S_extrap to true pi(x).

Code: `experiments/proposals/pade_zero_sum.py`

---

## Proposal B: Lagrange–Bürmann Inversion with Borel Resummation

### Idea
The Cipolla / asymptotic expansion of p(n):
```
p(n) = n*ln(n) + n*(ln ln n - 1) + n*(ln ln n - 2)/ln n + O(n*(ln ln n)^2/ln^2 n)
```
is a divergent (asymptotic) series. The standard prescription truncates at
the optimal index, getting roughly half the digits.

**Fresh angle.** Apply Borel summation to the *full formal series* obtained
by Lagrange–Bürmann inversion of pi(x) ~ Li(x) (not pi(x) ~ x/log x).
The formal series for Li^{-1}(n) has cleaner combinatorial structure: its
coefficients are signed Stirling-like numbers.

If the Borel transform of this formal series has finite radius of convergence
(unproven), the resummed series gives p(n) to arbitrary precision in
polylog(n) operations.

### Pseudocode
```
def p_lagrange_borel(n, K=30):
    coeffs = lagrange_inversion_coeffs_of_li(K)   # K coeffs of formal series
    borel  = [coeffs[k] / factorial(k) for k in range(K)]
    return integrate_borel(borel, log(n))
```

### Complexity
- Coefficient generation: O(K^2) symbolic arithmetic.
- Borel integration: O(K log K).
- Total O((log n)^c) if K = O(polylog n) suffices.

### Key assumption
The Borel transform converges for the inversion series of Li. If it does,
the analytic continuation gives p(n) without a zero-sum correction.
(If not, the experiment will show divergence and we close cleanly.)

### Test design
For n in {10, 100, 1000, 10000}:
- Generate first 30 Lagrange coefficients of Li^{-1}.
- Compute partial sums and Borel-resummed values.
- Plot |truncated value - true p(n)| vs K. If error decays exponentially in K,
  this is polylog. If it bottoms out at sqrt(p(n)), it's the standard barrier.

Code: `experiments/proposals/lagrange_borel.py`

---

## Proposal C: Modular-Form Fingerprint as a Polylog Primality Oracle

### Idea
Hecke eigenvalues of weight-k cusp forms encode arithmetic data of primes.
For Ramanujan's Delta, tau(p) is computable in O(polylog(p)) via Eichler–
Selberg + class numbers. The Ramanujan congruence tau(p) ≡ p^11 + 1 (mod 691)
holds for *all* primes p.

**Question.** Does there exist a finite collection of modular forms whose
joint Hecke fingerprint identifies primes among integers, in the sense that
*for every prime p and every composite n*, at least one form distinguishes
them?

If yes, we have a polylog-time membership test for primes (faster than the
AKS log^c bound, conjecturally). Combined with binary search, we get p(n)
in polylog. The bottleneck moves to *counting* primes with the fingerprint
in [1,x], which is a different problem.

This proposal is mainly a falsifiable conjecture: measure the discrimination
rate of small fingerprints on n <= 10000.

### Pseudocode
```
def fingerprint(n):
    return (tau_mod(n, 691), tau_mod(n, 5), j_invariant_mod(n, 11), ...)

def is_prime_oracle(n):
    return fingerprint(n) in PRIME_FINGERPRINTS

# Test: do all primes p <= 10000 have a unique fingerprint pattern? Do composites?
```

### Complexity
- Each tau(n) mod l: O(polylog n) via the standard recursion plus modular
  inverse on Hecke eigenvalues at prime power levels.
- Fingerprint comparison: O(polylog n) per check.

### Key assumption
A polylog-size collection of Hecke eigenvalues (mod small primes) separates
primes from composites among n <= N. This is a strengthening of known
modular form / Galois representation properties.

### Test design
For n = 1..10000:
- Compute fingerprint via tau(n) mod 691 (using sympy).
- Augment with sigma_k(n) mod l for several small l, k.
- Measure: how many composites collide with primes? Build a confusion matrix.
- If discrimination rate exceeds 99% with ~10 features, this is promising.

Code: `experiments/proposals/modular_fingerprint.py`

---

## Proposal D: Weighted Random Sampling with Zero-Aware Variance Reduction

### Idea
A naive Monte-Carlo estimator pi(x) ~ x * mean(is_prime(U_k)) has variance
~ x. Importance-weighted by 1/log k drops variance by a log factor.

**Fresh angle.** Use the explicit formula as a *control variate*: the
truncated zero sum sum_{|rho|<T} R(x^rho) is computable in O(polylog) and
is correlated with the random fluctuations of pi up to height T. Subtracting
this control variate from the Monte-Carlo estimator should drop variance by
factor T.

If T = polylog(x) zeros suffice to make residual variance polylog, then
O(polylog) random samples suffice and the algorithm is polylog.

This is empirically falsifiable: measure variance of the residual
(Monte-Carlo estimate of pi(x)) - (Riemann explicit formula truncated at
polylog zeros) for x growing.

### Pseudocode
```
def pi_control_variate(x, num_zeros, num_samples):
    rhos = load_zeros(num_zeros)
    expl = R(x) - sum(R(x**rho) for rho in rhos) - log(2)
    # Monte-Carlo of the residual fluctuation
    samples = [is_prime(uniform_int(1,x)) - 1/log(k) for k in random_indices]
    correction = x * mean(samples)
    return round(expl + correction)
```

### Complexity
- Computing the explicit-formula truncation: O(num_zeros * log x).
- Variance of residual after control variate: ?
- Conjectured polylog if zero correlations capture most variance.

### Key assumption
Var(pi(x) - explicit_formula_truncated_at_T) = O(x/T^{1+epsilon}) for some
epsilon > 0. (Standard bounds give epsilon=0.) The substantive claim is
that GUE statistics enable variance reduction beyond the trivial bound.

### Test design
For x = 1000, 10000, 100000:
- Compute explicit formula truncation at T = 10, 50, 200 zeros.
- Compute residual = pi(x) - truncation.
- Plot residual^2 vs T on log-log scale.
- If slope < -1, the bet wins; if = -1 (the predicted classical bound),
  it's exactly the sqrt-x barrier we already know.

Code: `experiments/proposals/zero_aware_variance.py`

---

## Summary

| Proposal | Polylog if... | Falsifiable? |
|----------|---------------|--------------|
| A. Pade–Borel zero sum | tail oscillation has smooth envelope | yes, n <= 10^4 |
| B. Lagrange–Borel inversion | Borel transform of Li^{-1} converges | yes, n <= 10^4 |
| C. Modular fingerprint | small Hecke fingerprint separates primes | yes, n <= 10^4 |
| D. Zero-aware control variate | variance reduction beats sqrt(x) | yes, x <= 10^5 |

Each test below n = 10000 returns either "barrier confirmed" or
"polylog plausible — escalate." All four were chosen so that a single small
experiment cleanly distinguishes the two outcomes.

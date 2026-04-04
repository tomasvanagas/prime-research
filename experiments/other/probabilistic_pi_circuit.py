#!/usr/bin/env python3
"""
Session 13: Probabilistic circuit for pi(x)

Key insight: Even if DETERMINISTIC circuits for pi(x) need size 2^{Theta(N)},
what about RANDOMIZED circuits (BPP-type)?

If PRIMES is in TC^0 (say via BPSW), then each primality test is O(poly(N))
size and constant depth. A randomized circuit could:

1. Uniformly sample M random numbers n_1, ..., n_M from [1, x]
2. Test each for primality (TC^0, cost poly(N) each)
3. Count how many are prime: C = #{i : n_i is prime}
4. Estimate pi(x) ≈ x * C / M

For EXACT pi(x), we need error < 0.5. By Chernoff:
P(|C/M - pi(x)/x| > 0.5/x) < 2*exp(-2*M*(0.5/x)^2)

For this to be < delta: M > x^2 * ln(2/delta) / 2

So M = Theta(x^2) random samples → circuit size Theta(x^2 * poly(N)).
This is EXPONENTIAL in N. WORSE than deterministic.

BUT: What about CLEVERER random estimation?

Idea: Binary search with randomized pi(x) estimation at each step.
If we can estimate pi(x) to within ±0.5 with O(polylog) circuit...

The issue is variance: Var(estimator) = pi(x)/x * (1-pi(x)/x) / M ≈ 1/(M*ln(x)).
For error < 0.5: M > 4*x / (pi(x)*ln(x)) ≈ 4*ln(x). Only O(log x) samples!

WAIT. Let me re-derive this more carefully.

We want: |hat{pi}(x) - pi(x)| < 0.5 with high probability.
hat{pi}(x) = x * (1/M) * sum_i 1_{n_i prime}.
E[hat{pi}(x)] = pi(x).
Var[hat{pi}(x)] = x^2 * (p*(1-p))/M where p = pi(x)/x ≈ 1/ln(x).
StdDev = x * sqrt(p*(1-p)/M).

For |error| < 0.5 w.h.p.: x * sqrt(p/M) < 0.5
→ M > 4 * x^2 * p = 4 * x^2 / ln(x) ≈ 4 * x * pi(x).

So M = Theta(x * pi(x)) = Theta(x^2 / ln(x)). Still exponential in N.

The fundamental issue: each sample gives ~1 bit of information,
and we need ~log(pi(x)) ≈ log(x/ln(x)) ≈ N bits of information.
But the VARIANCE of each sample is high (~p = 1/ln(x)), so many
samples are needed to reduce the variance.

ALTERNATIVE: Use importance sampling with proposal distribution.
If we could sample from a distribution biased toward primes...

P(n is prime) ≈ 1/ln(n). So the "natural" importance weight is ln(n).
Under importance sampling with g(n) ∝ 1/(n*ln(n)):
hat_pi = (1/M) * sum_i 1_{n_i prime} / g(n_i)

But sampling from g(n) requires knowing which n are prime... circular!

What about using the SMOOTH part as a guide? We know R(x) ≈ pi(x).
Can we use R to define a non-uniform sampling distribution that
reduces variance?

Test this numerically.
"""

import numpy as np
from sympy import primepi, isprime
from scipy.special import expi
import time

def li(x):
    if x <= 1: return 0.0
    return float(expi(np.log(x)))

def R_density(n, x):
    """Approximate density of primes near n, from R(x) theory."""
    if n < 2: return 0.0
    return 1.0 / np.log(max(n, 2.01))

def experiment_sampling_variance():
    """How many samples do we need for exact pi(x)?"""
    print("="*70)
    print("EXPERIMENT: Random Sampling Circuit for pi(x)")
    print("="*70)

    for x in [100, 1000, 10000, 100000]:
        N = int(np.ceil(np.log2(x)))
        pi_x = int(primepi(x))
        p = pi_x / x

        print(f"\nx = {x} (N={N} bits), pi(x) = {pi_x}, p = {p:.4f}")

        # Theoretical M needed
        # For Chebyshev: P(|hat - pi| > 0.5) < Var/(0.5)^2 = 4*x^2*p*(1-p)/M
        # For this < 0.01: M > 400 * x^2 * p * (1-p)
        M_theory = int(400 * x**2 * p * (1-p))
        print(f"  Theoretical M for 99% confidence: {M_theory:,}")
        print(f"  Circuit size: {M_theory} * poly({N}) ≈ {M_theory * N**2:,}")
        print(f"  vs deterministic: x^{2/3} ≈ {int(x**(2/3)):,}")

        # Empirical test
        trials = 100
        for M in [10, 100, 1000, 10000]:
            if M > x * 100: break
            successes = 0
            for trial in range(trials):
                samples = np.random.randint(1, x+1, size=M)
                count = sum(1 for n in samples if isprime(int(n)))
                estimate = x * count / M
                if abs(estimate - pi_x) < 0.5:
                    successes += 1
            print(f"  M={M:6d}: exact in {successes}/{trials} trials "
                  f"({100*successes/trials:.0f}%)")

def experiment_importance_sampling():
    """Can importance sampling reduce the number of samples?"""
    print("\n" + "="*70)
    print("EXPERIMENT: Importance Sampling with PNT density")
    print("="*70)

    for x in [1000, 10000]:
        pi_x = int(primepi(x))

        print(f"\nx = {x}, pi(x) = {pi_x}")

        trials = 100
        for M in [50, 200, 1000, 5000]:
            if M > x * 10: break

            # Uniform sampling
            uniform_exact = 0
            # PNT-weighted sampling (sample more from regions with higher prime density)
            pnt_exact = 0

            for trial in range(trials):
                # Uniform
                samples = np.random.randint(2, x+1, size=M)
                count = sum(1 for n in samples if isprime(int(n)))
                estimate = (x-1) * count / M + 1  # +1 for prime 2
                if abs(estimate - pi_x) < 0.5:
                    uniform_exact += 1

                # PNT-weighted: sample more from small numbers (higher prime density)
                # Use 1/ln(n) weighting via rejection sampling
                pnt_samples = []
                attempts = 0
                while len(pnt_samples) < M and attempts < M * 20:
                    n = np.random.randint(2, x+1)
                    # Accept with probability proportional to 1/ln(n)
                    if np.random.random() < 1.0 / np.log(n):
                        pnt_samples.append(n)
                    attempts += 1

                if len(pnt_samples) > 0:
                    # Importance weight: w_i = ln(n_i) (correct for biased sampling)
                    count_w = sum(np.log(n) for n in pnt_samples if isprime(int(n)))
                    total_w = sum(np.log(n) for n in pnt_samples)
                    if total_w > 0:
                        # Estimate p = pi(x)/(x-1)
                        p_est = count_w / total_w
                        est = p_est * (x-1) + 1
                        if abs(est - pi_x) < 0.5:
                            pnt_exact += 1

            print(f"  M={M:5d}: uniform {uniform_exact}/{trials} ({100*uniform_exact/trials:.0f}%), "
                  f"PNT-weighted {pnt_exact}/{trials} ({100*pnt_exact/trials:.0f}%)")

def experiment_hash_counting():
    """
    Alternative: Can we compute pi(x) by hashing?

    If h is a random hash function h: [x] -> [K], then
    #{n ≤ x prime : h(n) = 0} ≈ pi(x)/K.

    By computing this count for K = 1, 2, 4, 8, ..., we can
    estimate pi(x) at different scales.

    But the counting still requires testing all n... circular!

    Unless we can compute #{n ≤ x : n prime AND h(n)=0} without
    testing each n individually. This would require:
    - A hash function computable in TC^0 (yes, many are)
    - A way to count n with h(n)=0 AND n prime without enumeration
    - This is EXACTLY the same problem as pi(x) but restricted to
      the hash bucket. Size of bucket ≈ x/K.
    - No advantage unless K is very large, which reduces accuracy.

    Verdict: Hash-based counting doesn't help.
    """
    print("\n" + "="*70)
    print("EXPERIMENT: Hash-Based Counting (theoretical analysis)")
    print("="*70)
    print("  Hash function h: [x] -> [K]")
    print("  Count primes in bucket 0: C_0 = #{n ≤ x : prime(n) AND h(n)=0}")
    print("  Estimate: pi(x) ≈ K * C_0")
    print()
    print("  Problem: computing C_0 requires testing all x/K numbers in bucket")
    print("  Cost: x/K primality tests + counting = O(x/K * poly(N))")
    print("  For exact result: K=1, cost = O(x * poly(N)). WORSE than sieve.")
    print("  For approximate: K large → C_0 has high variance.")
    print()
    print("  The hash doesn't REDUCE computation — it just PARTITIONS it.")
    print("  Each partition requires the same total work as the full problem.")
    print()
    print("  VERDICT: Hash counting is EQUIVALENT to direct counting. CLOSED.")

if __name__ == "__main__":
    experiment_sampling_variance()
    experiment_importance_sampling()
    experiment_hash_counting()

    print("\n" + "="*70)
    print("CONCLUSION")
    print("="*70)
    print("""
Random sampling for exact pi(x) requires M = Theta(x^2/ln(x)) samples.
This is EXPONENTIAL in N = log(x). Even importance sampling with PNT
density as proposal distribution doesn't help — the variance of each
sample is Theta(1/ln(x)), requiring Theta(x/ln(x)) samples to achieve
error < 0.5 in the ABSOLUTE count.

Probabilistic circuits cannot beat deterministic circuits for EXACT pi(x).
The reason: exact counting requires information about EVERY prime up to x,
and each sample only reveals O(1) bits.

This is the COUNTING barrier: even if individual primality tests are cheap
(TC^0), counting how many pass requires looking at (nearly) all candidates.
""")

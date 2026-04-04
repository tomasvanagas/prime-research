"""
Session 8: Multiplicative function inversion approach.

The prime indicator 1_P(n) is NOT multiplicative. But we can express it via:
1_P(n) = -Σ_{d|n} μ(d) · [n/d > 1 and n/d = 1] ... that's just checking divisors.

Novel idea: Express the prime indicator as a Dirichlet convolution:
  1_P = f * g  where f and g are "simpler" functions.

We know: Λ(n) = -Σ_{d|n} μ(d) ln(d) (von Mangoldt)
And: 1_P(n) = Σ_{k: p^k=n} 1 for some prime p... that's circular.

Actually: Ω(n) = Σ_{p^k||n} k (total prime factors with multiplicity)
And: ω(n) = Σ_{p|n} 1 (distinct prime factors)

For prime p: ω(p) = 1, Ω(p) = 1
For prime power p^k: ω(p^k) = 1, Ω(p^k) = k
For composite n with ≥2 distinct factors: ω(n) ≥ 2

So: n is prime iff ω(n) = 1 and Ω(n) = 1
Equivalently: n is prime iff n > 1 and Ω(n) = 1

Can we compute ω(n) or Ω(n) without factoring?

Another angle: Liouville's function λ(n) = (-1)^{Ω(n)}
Σ_{n≤x} λ(n) = O(√x) under RH. But computing this sum is O(x^{2/3}).

Novel: What if we compute the "prime count" via MULTIPLICATIVE functions?
π(x) = Σ_{n≤x} 1_P(n)
     = Σ_{n≤x} Σ_{d|n} f(d) for some f (if 1_P can be expressed as Σ f(d))

Since 1_P is not multiplicative, its Dirichlet series doesn't factor over primes.
P(s) = Σ p^{-s} = Σ_{k≥1} μ(k)/k · log ζ(ks)

Can we evaluate P(s) at specific s to extract individual primes?
P(s) is analytic for Re(s) > 1 with natural boundary at Re(s) = 0.

Key idea: P(s) - 1/p_n^s represents the sum WITHOUT the n-th prime.
If we know P(s) and P_without_n(s), then p_n = (P(s) - P_without_n(s))^{-1/s}... circular.

Let me try something computational: Dirichlet series evaluation.
"""

import numpy as np
from sympy import primerange, mobius, factorint
from mpmath import mp, mpf, log as mplog, zeta, fsum, power
import time

mp.dps = 30

def prime_zeta(s, N=10000):
    """Compute P(s) = Σ_{p≤N} p^{-s} via Möbius inversion of log ζ."""
    # P(s) = Σ_{k=1}^∞ μ(k)/k · log ζ(ks)
    # Converges for Re(s) > 1/2 (but slowly)
    s = mpf(s)
    result = mpf(0)
    for k in range(1, 30):
        mu_k = int(mobius(k))
        if mu_k == 0:
            continue
        ks = k * s
        if ks > 50:
            break  # log ζ(ks) ≈ 0
        log_zeta_ks = mplog(zeta(ks))
        result += mpf(mu_k) / k * log_zeta_ks
    return result

def prime_zeta_direct(s, N=10000):
    """Direct computation of P(s) for comparison."""
    primes = list(primerange(2, N))
    return fsum(power(p, -mpf(s)) for p in primes)

print("=" * 60)
print("MULTIPLICATIVE FUNCTION INVERSION")
print("=" * 60)

# Test: P(s) computation via Möbius vs direct
print("\n--- Test 1: Prime zeta function P(s) ---")
for s in [2, 3, 4, 5]:
    p_mobius = prime_zeta(s)
    p_direct = prime_zeta_direct(s)
    print(f"  P({s}): Möbius={float(p_mobius):.10f}, Direct={float(p_direct):.10f}, diff={float(abs(p_mobius-p_direct)):.2e}")

# Novel idea: Can we extract p(n) from P(s) evaluated at multiple s values?
# P(s) = 2^{-s} + 3^{-s} + 5^{-s} + ...
# If we evaluate P(s) at s=2,3,4,...,K, we get K equations in infinitely many unknowns.
# But the contributions from large primes are exponentially small.
# For s=2: P(2) ≈ 0.4522 (dominated by 2^{-2} = 0.25)
# For large s: P(s) → 2^{-s} (dominated by smallest prime)

print("\n--- Test 2: Extracting primes from P(s) values ---")
# Newton's identities: from power sums S_k = Σ p_i^k we can recover p_i
# But we have P(s) = Σ p^{-s}, not power sums.
# S_k = Σ p^k → ∞ (diverges). So we can't use Newton's identities directly.

# Instead: from P(s) = Σ p^{-s} for various s, we can use Prony's method
# to extract the "frequencies" (which are the primes)

# Prony's method: Given signal y(s) = Σ a_k · r_k^s, recover r_k
# Here: y(s) = P(s) = Σ p^{-s} = Σ (1/p)^s, so r_k = 1/p_k, a_k = 1

# Compute P(s) at integer points s = 2, 3, ..., 2M+1
M = 10  # Try to recover 10 primes
samples = []
for s in range(2, 2*M + 2):
    ps = float(prime_zeta_direct(s, N=200))
    samples.append(ps)

print(f"  P(s) samples for s=2..{2*M+1}: {[f'{x:.6f}' for x in samples[:5]]}")

# Form Hankel matrix for Prony
H = np.array([[samples[i+j] for j in range(M)] for i in range(M)])
h = np.array(samples[M:2*M])

# Solve Hc = -h for polynomial coefficients
try:
    c = np.linalg.solve(H, -h)
    # Roots of z^M + c[M-1]z^{M-1} + ... + c[0]
    poly = np.concatenate([np.array([1]), c[::-1]])
    roots = np.roots(poly)

    # r_k = 1/p_k, so p_k = 1/r_k
    recovered_primes = sorted([1/r.real for r in roots if abs(r.imag) < 0.1 and r.real > 0])

    actual_primes = list(primerange(2, 50))
    print(f"\n  Recovered 'primes' via Prony: {[f'{p:.2f}' for p in recovered_primes[:10]]}")
    print(f"  Actual first 10 primes: {actual_primes[:10]}")

    # How many are close?
    matched = 0
    for rp in recovered_primes:
        for ap in actual_primes:
            if abs(rp - ap) < 0.5:
                matched += 1
                break
    print(f"  Primes recovered correctly: {matched}/{len(recovered_primes)}")
except np.linalg.LinAlgError:
    print("  Hankel matrix singular — Prony method fails")

# Test 3: Can we use Prony more carefully with high precision?
print("\n--- Test 3: High-precision Prony extraction ---")

M = 5  # Just try to extract first 5 primes
samples_hp = []
for s in range(2, 2*M + 2):
    ps = float(prime_zeta_direct(mpf(s), N=200))
    samples_hp.append(ps)

H_hp = np.array([[samples_hp[i+j] for j in range(M)] for i in range(M)])
h_hp = np.array(samples_hp[M:2*M])

try:
    c_hp = np.linalg.solve(H_hp, -h_hp)
    poly_hp = np.concatenate([np.array([1]), c_hp[::-1]])
    roots_hp = np.roots(poly_hp)
    primes_hp = sorted([1/r.real for r in roots_hp if abs(r.imag) < 0.1 and r.real > 0])
    print(f"  Recovered primes (M=5): {[f'{p:.4f}' for p in primes_hp]}")
    print(f"  Actual: [2, 3, 5, 7, 11]")
except:
    print("  Failed")

# Test 4: Novel — use DERIVATIVES of P(s)
print("\n--- Test 4: Derivatives of P(s) ---")
# P'(s) = -Σ p^{-s} ln(p)
# P''(s) = Σ p^{-s} (ln p)^2
# At s → ∞: P(s) → 2^{-s}, so P'(s)/P(s) → -ln(2)
# At finite s: the ratio P'(s)/P(s) is a weighted average of -ln(p_k)

# P'(s)/P(s) = -Σ w_k(s) ln(p_k) where w_k = p_k^{-s}/P(s)

for s in [2, 5, 10, 20]:
    primes = list(primerange(2, 200))
    ps_val = sum(p ** (-s) for p in primes)
    ps_deriv = sum(-(p ** (-s)) * np.log(p) for p in primes)
    ratio = ps_deriv / ps_val
    effective_prime = np.exp(-ratio)
    weights = [(p ** (-s)) / ps_val for p in primes[:5]]
    print(f"  s={s}: P'(s)/P(s) = {ratio:.4f}, 'effective prime' = {effective_prime:.2f}, "
          f"top weights = {[f'{w:.3f}' for w in weights]}")

print(f"\n  At large s, the 'effective prime' → 2 (smallest prime dominates)")
print(f"  At s=2, it's ~2.6 (mix of small primes)")
print(f"  This gives a smooth function of s, not individual primes")

# Test 5: Sieve-like Dirichlet convolution
print("\n--- Test 5: Can we compute π(x) via Dirichlet convolution? ---")

# π(x) = Σ_{n≤x} 1_P(n)
# Where 1_P(n) = Σ_{d|n} μ(d) · [ω(n/d) ≥ 1] · [n/d is prime power]
# More precisely: 1_P = Λ/(ln) as Dirichlet series (sum only k=1 terms)
# P(s) = Σ μ(k)/k · log ζ(ks)

# The Möbius inversion: 1_P(n) = Σ_{d|n} μ(d) · f(n/d)
# where f(m) = 1 if m is a prime power, 0 otherwise... wait

# Actually: Λ(n) = Σ_{d|n} μ(d) · (-ln(n/d)) = -Σ_{d|n} μ(d)ln(n/d)
# And prime indicator: 1_P(n) = Σ_{k≥1} [n is k-th power of a prime] · ... hmm

# The simplest Dirichlet relation:
# Ω(n) = Σ_{d|n} Λ(d)/ln(d) where we define Λ(1)/ln(1) = 0
# But this needs all prime power divisors

# There's no Dirichlet convolution that avoids the barrier
print("  π(x) via Dirichlet convolution requires knowing divisor structure → O(x^{2/3})")

print("\n" + "=" * 60)
print("CONCLUSIONS (Multiplicative Inversion)")
print("=" * 60)
print("""
1. P(s) = Σ p^{-s}: computable via Möbius inversion of log ζ, but requires ζ values
2. Prony method on P(s) samples: Recovers first few primes only (~2 correct)
   - Hankel matrix ill-conditioned, cannot extract beyond p≈7
3. P(s) derivatives: give weighted averages, not individual primes
4. Dirichlet convolution: no shortcut for π(x)

The multiplicative structure of ζ(s) = Π(1-p^{-s})^{-1} encodes primes,
but EXTRACTING individual primes requires either:
  - Computing P(s) to high precision (needs ζ → needs O(T^{1/2}) terms)
  - Solving an exponential system (Prony → ill-conditioned for >5 terms)

VERDICT: Multiplicative inversion does NOT bypass the O(x^{2/3}) barrier.
""")

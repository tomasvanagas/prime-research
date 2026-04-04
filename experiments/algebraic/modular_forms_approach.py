"""
Session 5: MODULAR FORMS approach to prime counting

NOVEL IDEA: The Ramanujan tau function τ(n) is the coefficient of q^n
in Δ(q) = q * prod_{n>=1} (1-q^n)^24, a weight 12 modular form.

By Deligne (proving the Ramanujan conjecture): |τ(p)| < 2*p^{11/2}
for primes p. More importantly, τ(p) has specific properties:
- τ(p) ≡ 1 + p^11 (mod 691)  [Ramanujan congruence]
- |τ(p)| < 2*p^{11/2}
- τ(p) determines the angle θ_p via τ(p) = 2*p^{11/2}*cos(θ_p)
  where θ_p follows the Sato-Tate distribution

QUESTION: Can the Sato-Tate distribution or Ramanujan congruences
help predict which integers are prime?

Also: The Hecke eigenform gives a_p for each prime p.
Reversing: if we know a_p, can we determine p?

APPROACH: Use modular form coefficients as a "prime detector"
a(n) has special values at primes vs composites.
"""

import numpy as np
import sympy
from sympy import prime as sympy_prime, isprime, primepi
from collections import Counter
import time


def compute_tau(n_max):
    """
    Compute Ramanujan's tau function tau(n) for n = 1..n_max
    using the product formula: Delta(q) = q * prod_{k>=1} (1-q^k)^24

    We work with coefficients of the q-expansion.
    """
    # Initialize coefficients array
    coeffs = np.zeros(n_max + 1, dtype=np.int64)
    coeffs[0] = 1  # Start with 1

    # Multiply by (1-q^k)^24 for each k
    for k in range(1, n_max + 1):
        # (1-q^k)^24 = sum_{j=0}^{24} C(24,j) * (-1)^j * q^{kj}
        # Apply this to the existing coefficients
        new_coeffs = np.zeros(n_max + 1, dtype=np.int64)

        for j in range(25):  # j = 0 to 24
            if k * j > n_max:
                break
            sign = (-1) ** j
            binom = 1
            for i in range(j):
                binom = binom * (24 - i) // (i + 1)

            coeff_j = sign * binom
            for m in range(n_max + 1):
                if m + k * j <= n_max:
                    new_coeffs[m + k * j] += coeff_j * coeffs[m]

        coeffs = new_coeffs

    # Delta(q) = q * prod, so tau(n) = coeffs[n-1]
    tau = np.zeros(n_max + 1, dtype=np.int64)
    for n in range(1, n_max + 1):
        tau[n] = coeffs[n - 1]

    return tau


def compute_tau_fast(n_max):
    """
    Faster computation using Ramanujan's formula:
    tau(n) can be computed using sigma functions.

    Ramanujan's identity: tau(n) = (65/756)*sigma_11(n) - (691/756)*sigma_5(n)
    ... wait, that's not right. Let me use a different approach.

    Actually, use the recurrence:
    n*tau(n) = tau(1)*sigma_11(n) - sum_{k=1}^{n-1} sigma_11(k)*tau(n-k) ... not standard.

    Simplest: use the partition-like expansion via Euler products.
    Or just use the convolution directly.
    """
    # Use the Jacobi theta function approach or direct q-expansion
    # For speed, compute via the Eisenstein series representation:
    # Delta = (E_4^3 - E_6^2) / 1728

    # E_4 = 1 + 240 * sum_{n>=1} sigma_3(n) * q^n
    # E_6 = 1 - 504 * sum_{n>=1} sigma_5(n) * q^n

    # Compute sigma_k(n) = sum_{d|n} d^k
    def sigma_k(n, k):
        result = 0
        for d in range(1, int(n**0.5) + 1):
            if n % d == 0:
                result += d**k
                if d != n // d:
                    result += (n // d)**k
        return result

    # Compute E_4 and E_6 coefficients
    e4 = np.zeros(n_max + 1, dtype=object)
    e6 = np.zeros(n_max + 1, dtype=object)
    e4[0] = 1
    e6[0] = 1

    for n in range(1, n_max + 1):
        e4[n] = 240 * sigma_k(n, 3)
        e6[n] = -504 * sigma_k(n, 5)

    # E_4^3 by convolution
    e4_sq = np.zeros(n_max + 1, dtype=object)
    for n in range(n_max + 1):
        s = 0
        for k in range(n + 1):
            s += e4[k] * e4[n - k]
        e4_sq[n] = s

    e4_cube = np.zeros(n_max + 1, dtype=object)
    for n in range(n_max + 1):
        s = 0
        for k in range(n + 1):
            s += e4_sq[k] * e4[n - k]
        e4_cube[n] = s

    # E_6^2 by convolution
    e6_sq = np.zeros(n_max + 1, dtype=object)
    for n in range(n_max + 1):
        s = 0
        for k in range(n + 1):
            s += e6[k] * e6[n - k]
        e6_sq[n] = s

    # Delta = (E_4^3 - E_6^2) / 1728
    # tau(n) = coefficient of q^n in Delta
    # Delta = sum_{n>=1} tau(n) * q^n
    # So Delta[0] = 0, Delta[n] = (e4_cube[n] - e6_sq[n]) / 1728

    tau = np.zeros(n_max + 1, dtype=object)
    for n in range(1, n_max + 1):
        val = e4_cube[n] - e6_sq[n]
        assert val % 1728 == 0, f"Not divisible: {val} at n={n}"
        tau[n] = val // 1728

    return tau


def analyze_tau_primes():
    """
    Analyze the Ramanujan tau function at primes vs composites.
    """
    print("=" * 80)
    print("RAMANUJAN TAU FUNCTION: PRIME vs COMPOSITE ANALYSIS")
    print("=" * 80)

    n_max = 500
    print(f"\nComputing tau(n) for n=1..{n_max}...")
    t0 = time.time()
    tau = compute_tau_fast(n_max)
    t1 = time.time()
    print(f"Computed in {t1-t0:.2f}s")

    # Verify: tau(1) = 1, tau(2) = -24, tau(3) = 252, tau(4) = -1472
    known = {1: 1, 2: -24, 3: 252, 4: -1472, 5: 4830, 6: -6048,
             7: -16744, 8: 84480, 9: -113643, 10: -115920}

    print("\nVerification:")
    for n, expected in known.items():
        actual = int(tau[n])
        match = "OK" if actual == expected else f"FAIL (got {actual})"
        print(f"  tau({n}) = {actual} (expected {expected}) — {match}")

    # Analyze tau at primes
    prime_tau = []
    composite_tau = []
    for n in range(2, n_max + 1):
        if isprime(n):
            prime_tau.append((n, int(tau[n])))
        else:
            composite_tau.append((n, int(tau[n])))

    print(f"\nPrime tau values (first 20):")
    for p, t in prime_tau[:20]:
        # Normalized: tau(p) / (2 * p^{11/2})
        norm = t / (2 * p**5.5)
        print(f"  tau({p:3d}) = {t:>15d}, normalized = {norm:>8.4f}")

    # Check Ramanujan congruence: tau(p) ≡ 1 + p^11 (mod 691)
    print(f"\nRamanujan congruence tau(p) ≡ 1 + p^11 (mod 691):")
    for p, t in prime_tau[:15]:
        expected_mod = (1 + pow(p, 11, 691)) % 691
        actual_mod = t % 691
        if actual_mod < 0:
            actual_mod += 691
        match = "OK" if actual_mod == expected_mod else "FAIL"
        print(f"  p={p:3d}: tau(p) mod 691 = {actual_mod:3d}, 1+p^11 mod 691 = {expected_mod:3d} — {match}")

    # Can tau values distinguish primes from composites?
    print(f"\nStatistics:")
    p_vals = [abs(t) for _, t in prime_tau]
    c_vals = [abs(t) for _, t in composite_tau if abs(t) > 0]

    p_norms = [abs(t) / n**5.5 for n, t in prime_tau]
    c_norms = [abs(t) / n**5.5 for n, t in composite_tau if abs(t) > 0]

    print(f"  |tau(p)| / p^5.5: mean={np.mean(p_norms):.4f}, std={np.std(p_norms):.4f}")
    print(f"  |tau(c)| / c^5.5: mean={np.mean(c_norms):.4f}, std={np.std(c_norms):.4f}")

    # Key: for composites n = p*q: tau(n) = tau(p)*tau(q) - p^11 * tau(n/p^2)
    # This multiplicativity is the key property.
    # For primes: tau(p) is "random" in [-2*p^{11/2}, 2*p^{11/2}]
    # For composites: tau is determined by factorization

    # Can we use tau to DETECT primes?
    # Yes, but computing tau(n) requires knowing the factorization of n,
    # OR computing the modular form coefficient directly (O(n^2) by convolution).

    print(f"\n  VERDICT: tau(n) can distinguish primes from composites")
    print(f"  (primes have |tau(p)/p^5.5| in Sato-Tate distribution)")
    print(f"  BUT: computing tau(n) costs O(n^2) — WORSE than direct primality testing")

    # Could Sato-Tate help PREDICT primes?
    # The angle θ_p = arccos(tau(p)/(2*p^{11/2})) follows sin^2(θ) distribution.
    # But this tells us NOTHING about which integers ARE prime.
    # It only describes the DISTRIBUTION of tau at prime inputs.

    print(f"\n  Sato-Tate distribution check:")
    angles = [np.arccos(np.clip(t / (2 * p**5.5), -1, 1)) for p, t in prime_tau if p > 2]
    if angles:
        # Should follow density (2/π) * sin²(θ) on [0, π]
        hist, bins = np.histogram(angles, bins=20, range=(0, np.pi))
        print(f"  Angle distribution (should be sin²-shaped):")
        for i in range(len(hist)):
            mid = (bins[i] + bins[i+1]) / 2
            expected = len(angles) / 20 * (2/np.pi) * np.sin(mid)**2 * np.pi
            bar = "#" * int(hist[i] / max(hist) * 40) if max(hist) > 0 else ""
            print(f"    θ={mid:.2f}: {hist[i]:3d} (expected ~{expected:.1f}) {bar}")


def analyze_modular_form_detection():
    """
    Can ANY modular form efficiently detect primes?

    Key insight: For a Hecke eigenform f of weight k and level N,
    the Fourier coefficients a_f(p) at primes satisfy:
    |a_f(p)| ≤ 2 * p^{(k-1)/2}  (Ramanujan-Petersson conjecture, proved by Deligne)

    For composite n = p*q: a_f(n) = a_f(p)*a_f(q) (multiplicativity)

    So: a_f(n) tells us about the FACTORIZATION of n.
    If n is prime: a_f(n) is "random" in its range.
    If n is composite: a_f(n) is determined by its factors.

    But computing a_f(n) for a single n is as hard as factoring n!
    (Unless we compute the entire sequence up to n.)

    CONCLUSION: Modular forms ENCODE prime information but COMPUTING
    their coefficients is at least as hard as detecting primes.
    """
    print("\n" + "=" * 80)
    print("MODULAR FORM PRIME DETECTION: FEASIBILITY ANALYSIS")
    print("=" * 80)

    print("""
  ANALYSIS:

  For a Hecke eigenform f of weight k, level N:

  1. COMPUTING a_f(n) for a single n:
     - If n is prime: O(N * k * n^{1/2}) via modular symbols
     - If n has known factorization: O(log(n) * k) via multiplicativity
     - Directly from q-expansion: O(n^2) by convolution

  2. COMPUTING a_f(n) for ALL n ≤ x:
     - O(x^2) by convolution
     - O(x * log(x)) via Hecke operators (with precomputation)

  3. USING a_f to count primes:
     - a_f is multiplicative → a_f(n) = prod_{p^k || n} (something depending on a_f(p))
     - Detecting primes via a_f: n is prime iff a_f(n) is NOT determined by
       any proper factorization. But checking this IS factoring.

  4. The Sato-Tate distribution:
     - a_f(p) / (2*p^{(k-1)/2}) ~ cos(θ_p) where θ_p follows sin²(θ)*2/π
     - This is STATISTICAL — gives distribution, not individual values
     - Cannot predict whether a specific n is prime

  5. L-functions of modular forms:
     - L(f, s) = sum_n a_f(n) / n^s
     - The analytic continuation and functional equation are known
     - The zeros of L(f, s) are related to prime distribution via explicit formula
     - BUT this is the SAME barrier as Riemann zeta — just a different L-function
     - Langlands program: all L-functions face the same fundamental limitation

  VERDICT: Modular forms provide BEAUTIFUL theory but NO computational shortcut.
  The arithmetic of modular forms is governed by the same spectral theory
  (zeta zeros / L-function zeros) that creates the barrier in the first place.
""")


if __name__ == "__main__":
    analyze_tau_primes()
    analyze_modular_form_detection()

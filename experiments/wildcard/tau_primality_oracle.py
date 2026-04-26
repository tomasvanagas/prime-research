"""
Ramanujan tau function as a polylog primality oracle.

Idea: tau(n) is the n-th Fourier coefficient of the modular form Delta
of weight 12. For prime p, Lehmer-Ramanujan congruences hold:
    tau(p) ≡ 1 + p^11 (mod 691)
    tau(p) ≡ p + p^10 (mod 25)
    tau(p) ≡ p + p^4 (mod 7)
    etc.

For composite n = a*b: tau(n) = tau(a)tau(b) (when gcd=1, multiplicativity).

KEY QUESTION: Is there a tau-based test that detects primality in O(polylog n)?
If tau(n) computed for arbitrary n in polylog (via Schoof-Pila style on the
ell-adic Galois representation rho_ell of Delta) AND a congruence cleanly
separates primes from composites, we'd have a TC^0-like primality oracle.

Schoof's algorithm computes tau(p) in polylog(p), assuming p is prime.
For arbitrary n, computing tau(n) requires knowing the factorization OR a
direct formula. Direct formulas exist (e.g. Niebur's recursion) but in
exp(O(log n)) time, not polylog.

This experiment:
1. Verify that tau(p) ≡ 1 + p^11 (mod 691) holds for primes.
2. Check if it FAILS for composites - giving a separator.
3. Estimate the computational complexity of evaluating tau(n) for general n.
"""

from sympy import isprime, primerange
from sympy.ntheory.modular import crt
import math


def sigma_k(n, k):
    """sum of d^k over divisors d of n."""
    s = 0
    d = 1
    while d * d <= n:
        if n % d == 0:
            s += d ** k
            if d != n // d:
                s += (n // d) ** k
        d += 1
    return s


def tau_recursion(N):
    """Compute tau(1..N) via Niebur's recurrence:
    n^2 * tau(n) = sigma_1(n)*(35*sigma_3(n) - 250*sigma_3(n)*?) ...
    For simplicity, use direct convolution:
        Delta(q) = q * prod_{n>=1}(1-q^n)^24
    Compute via Euler-style product expansion to N terms.
    """
    # Compute eta(q)^24 = prod (1-q^n)^24 truncated to degree N
    coeffs = [0] * (N + 1)
    coeffs[0] = 1
    for n in range(1, N + 1):
        # multiply by (1-q^n)^24 = sum_{k=0..24} C(24,k) (-q^n)^k
        new = [0] * (N + 1)
        for k in range(0, 25):
            if k * n > N:
                break
            sign = (-1) ** k
            from math import comb
            c = sign * comb(24, k)
            for j in range(N + 1 - k * n):
                new[j + k * n] += c * coeffs[j]
        coeffs = new
    # Delta = q * eta^24, so tau(n) = coeffs[n-1]
    tau = [0] * (N + 1)
    for n in range(1, N + 1):
        tau[n] = coeffs[n - 1]
    return tau


def test_lehmer_congruences(N=200):
    """For all primes p <= N, check tau(p) mod 691 == (1 + p^11) mod 691."""
    tau = tau_recursion(N)
    print(f"Computed tau(1..{N})")
    fails = []
    for p in primerange(2, N + 1):
        lhs = tau[p] % 691
        rhs = (1 + pow(p, 11, 691)) % 691
        if lhs != rhs:
            fails.append((p, lhs, rhs))
    print(f"Lehmer mod 691: {len(fails)} failures among primes up to {N}")

    # Now test composites: does tau(n) mod 691 != (1 + n^11) mod 691?
    print("\nFor composites n in [4..N], check if tau(n) mod 691 deviates from (1+n^11) mod 691:")
    composite_match = 0
    composite_fail = 0
    for n in range(4, N + 1):
        if isprime(n):
            continue
        lhs = tau[n] % 691
        rhs = (1 + pow(n, 11, 691)) % 691
        if lhs == rhs:
            composite_match += 1
        else:
            composite_fail += 1
    print(f"  Composites where tau(n) ≡ 1+n^11 (mod 691): {composite_match}")
    print(f"  Composites where tau(n) ≢ 1+n^11 (mod 691): {composite_fail}")
    print(f"  -> separator strength: {composite_fail}/{composite_match + composite_fail}")
    return composite_match, composite_fail


def test_multi_modulus_separator(N=500):
    """
    Combine multiple congruences. For prime p:
      tau(p) ≡ 1 + p^11 (mod 691)
      tau(p) ≡ 1 + p (mod 5)         [actually mod 25, simplified]
      tau(p) ≡ 1 + p^4 (mod 7)
    Use product congruence and check separation strength on composites.
    """
    tau = tau_recursion(N)

    moduli = [
        (691, lambda p: (1 + pow(p, 11, 691)) % 691),
        (7, lambda p: (1 + pow(p, 4, 7)) % 7),
        (5, lambda p: (1 + pow(p, 7, 5)) % 5),
    ]

    prime_pass = 0
    prime_total = 0
    for p in primerange(2, N + 1):
        prime_total += 1
        passes = all((tau[p] % m) == fn(p) for m, fn in moduli)
        if passes:
            prime_pass += 1

    composite_pass = 0
    composite_total = 0
    for n in range(4, N + 1):
        if isprime(n):
            continue
        composite_total += 1
        passes = all((tau[n] % m) == fn(n) for m, fn in moduli)
        if passes:
            composite_pass += 1

    print(f"\nMulti-congruence test (mod 691, 7, 5):")
    print(f"  Primes passing: {prime_pass}/{prime_total} (should be 100%)")
    print(f"  Composites passing: {composite_pass}/{composite_total}")
    print(f"  -> false-positive rate on composites: {composite_pass / max(1, composite_total):.4f}")
    print(f"  -> if this rate were 0, we'd have a primality test from tau alone.")


def main():
    print("=" * 60)
    print("Ramanujan tau as primality oracle")
    print("=" * 60)
    test_lehmer_congruences(200)
    test_multi_modulus_separator(500)
    print("\nVerdict context:")
    print("- tau(n) for general n: O(n) classical, O(sqrt(n)) via Schoof-style if n is prime.")
    print("- For composite n, no known polylog algorithm to compute tau(n) directly.")
    print("- So even if the congruences perfectly separated primes, the oracle would")
    print("  not be polylog-time. Multi-congruence false-positive rate quantifies the gap.")


if __name__ == "__main__":
    main()

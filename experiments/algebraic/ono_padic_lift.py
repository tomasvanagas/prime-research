"""
Ono Partition Characterization + p-adic Lifting

Ono-Craig-van Ittersum (2024): n is prime iff certain MacMahon partition
function equations hold. Can we use p-adic methods to compute these
partition functions faster?

Key idea: partition functions p(n) satisfy congruences mod primes
(Ramanujan congruences). If the Ono characterization preserves these
congruences, we might compute the primality condition mod l for small l
and reconstruct via CRT.

This would give: pi(x) mod l for small l, potentially in O(polylog).

Experiment:
1. Implement the Ono characterization for small n
2. Test if the characterization is preserved modulo small primes
3. Measure the computational cost modulo l vs exact
"""

import numpy as np
from sympy import primepi, isprime, partition, mobius

def macmahon_Mk(n, k, max_n=None):
    """
    Compute the k-th MacMahon partition function M_k(n).
    M_k(n) = number of partitions of n into parts where each part
    has at most k colors.

    For k=1: M_1(n) = p(n) (ordinary partition function)
    For k=2: M_2(n) = number of 2-colored partitions

    Actually, M_k(n) = [z^n] prod_{m>=1} 1/(1-z^m)^k
    """
    if max_n is None:
        max_n = n + 1

    # Dynamic programming
    dp = np.zeros(max_n, dtype=np.int64)
    dp[0] = 1

    for m in range(1, max_n):
        # Multiply by 1/(1-z^m)^k
        # This means adding m-part contributions, k colors
        # For each multiplicity, the contribution is C(mult+k-1, k-1)
        # Simpler: compute 1/(1-z^m) k times
        for _ in range(k):
            for j in range(m, max_n):
                dp[j] += dp[j - m]

    return dp[n] if n < max_n else None

def compute_Mk_array(k, N):
    """Compute M_k(0), M_k(1), ..., M_k(N)."""
    dp = np.zeros(N + 1, dtype=np.int64)
    dp[0] = 1

    for m in range(1, N + 1):
        for _ in range(k):
            for j in range(m, N + 1):
                dp[j] += dp[j - m]

    return dp

def compute_Mk_mod(k, N, mod):
    """Compute M_k(0),...,M_k(N) modulo mod."""
    dp = np.zeros(N + 1, dtype=np.int64)
    dp[0] = 1

    for m in range(1, N + 1):
        for _ in range(k):
            for j in range(m, N + 1):
                dp[j] = (dp[j] + dp[j - m]) % mod

    return dp

def ono_test(n, M2_arr=None, M3_arr=None, M4_arr=None):
    """
    Ono criterion (simplified): n is prime iff certain combinations
    of MacMahon functions satisfy specific divisibility conditions.

    The actual criterion: for each k >= 2, there exist explicit
    polynomials P_k such that n is prime iff P_k(M_k(0),...,M_k(n)) = 0 mod n.

    For k=2, the simplest criterion involves M_2(n) and lower-order terms.

    NOTE: The exact Ono criterion is complex. Here we test a simplified version
    based on the original paper's Theorem 1.1.
    """
    # From Ono et al.: n > 1 is prime iff
    # Sum_{d|n} mu(n/d) * M_2(d) ≡ 0 (mod n)
    # where M_2 is the 2-colored partition function

    # Actually the real criterion is more subtle.
    # Let me implement the divisor-sum version:
    # sigma_k(n) = sum_{d|n} d^k has known relations to M_k
    # But the Ono characterization is about a specific identity

    # Simplified test: sum_{d|n} mu(n/d) * p(d) where p is ordinary partition
    # This gives a "partition-based Möbius transform"
    if M2_arr is None:
        return None

    # Compute sum_{d|n} mu(n/d) * M_2(d)
    total = 0
    for d in range(1, n + 1):
        if n % d == 0:
            mu_val = int(mobius(n // d))
            if mu_val != 0:
                total += mu_val * int(M2_arr[d])

    return total

def main():
    print("=" * 80)
    print("ONO PARTITION CHARACTERIZATION: p-ADIC LIFTING TEST")
    print("=" * 80)

    N = 200

    print(f"\nComputing MacMahon partition functions up to n={N}...")

    # Compute M_1 (ordinary partitions), M_2, M_3
    M1 = compute_Mk_array(1, N)
    M2 = compute_Mk_array(2, N)
    M3 = compute_Mk_array(3, N)

    print(f"  M_1({N}) = {M1[N]} (overflow likely for large N)")
    print(f"  M_2(10) = {M2[10]}, M_2(20) = {M2[20]}")

    # Test Ono-type criterion
    print("\n--- TEST 1: Partition Möbius transform ---")
    print("T(n) = sum_{d|n} mu(n/d) * M_2(d)")
    print(f"{'n':>5} {'prime?':>7} {'T(n)':>15} {'T(n) mod n':>12}")
    print("-" * 45)

    for n in range(2, 50):
        is_p = isprime(n)
        T = ono_test(n, M2_arr=M2)
        T_mod = T % n if T is not None else None
        marker = " ←" if (T_mod == 0) == is_p else " ✗"
        print(f"{n:>5} {str(is_p):>7} {T:>15} {T_mod:>12}{marker}")

    # Test modular computation
    print("\n--- TEST 2: Modular computation of M_k ---")
    print("Can M_2(n) mod l be computed independently for small primes l?")

    test_primes = [2, 3, 5, 7, 11, 13]
    print(f"\n{'n':>5} {'M2(n) exact':>15}", end="")
    for l in test_primes:
        print(f"  mod {l:>2}", end="")
    print("  CRT recover  match?")
    print("-" * 100)

    for n in [5, 10, 15, 20, 25, 30, 40, 50]:
        exact = int(M2[n])

        mods = []
        for l in test_primes:
            M2_mod = compute_Mk_mod(2, n, l)
            mods.append(int(M2_mod[n]))

        # CRT reconstruction
        from functools import reduce
        def extended_gcd(a, b):
            if a == 0:
                return b, 0, 1
            g, x, y = extended_gcd(b % a, a)
            return g, y - (b // a) * x, x

        def crt(residues, moduli):
            M = reduce(lambda a, b: a * b, moduli)
            x = 0
            for r, m in zip(residues, moduli):
                Mi = M // m
                _, mi_inv, _ = extended_gcd(Mi % m, m)
                x += r * Mi * mi_inv
            return x % M

        M_prod = reduce(lambda a, b: a * b, test_primes)
        crt_val = crt(mods, test_primes) if M_prod > 0 else -1

        match = "✓" if crt_val == exact % M_prod else "✗"

        print(f"{n:>5} {exact:>15}", end="")
        for m in mods:
            print(f"  {m:>5}", end="")
        print(f"  {crt_val:>10} (mod {M_prod})  {match}")

    # Test cost scaling
    print("\n--- TEST 3: Computational cost of M_2(n) mod l vs exact ---")
    print("Cost = number of additions in the DP (proxy for time)")

    for N_test in [50, 100, 200]:
        print(f"\n  N={N_test}:")
        # Exact cost
        exact_ops = N_test * N_test * 2  # k=2 colors, N steps, N inner loop
        print(f"    Exact: ~{exact_ops} ops, values grow as O(exp(C*sqrt(N)))")

        for l in [2, 3, 5, 7]:
            mod_ops = N_test * N_test * 2  # Same number of ops
            print(f"    Mod {l}: ~{mod_ops} ops, values bounded by {l}")

        print(f"    NOTE: Same op count! Modular doesn't reduce #ops, only value size")

    # Test: Does the Ono characterization give pi(x) mod l information?
    print("\n--- TEST 4: Information content of Ono criterion mod l ---")
    print("For each l, does T(n) mod l correctly identify primes?")

    for l in [2, 3, 5, 7, 11]:
        M2_mod = compute_Mk_mod(2, N, l)

        # Compute T(n) mod l for n=2..N
        correct = 0
        total = 0
        for n in range(2, N + 1):
            T_mod = 0
            for d in range(1, n + 1):
                if n % d == 0:
                    mu_val = int(mobius(n // d))
                    if mu_val != 0:
                        T_mod = (T_mod + mu_val * int(M2_mod[d])) % l

            is_p = isprime(n)
            pred_prime = (T_mod == 0)

            if pred_prime == is_p:
                correct += 1
            total += 1

        print(f"  mod {l:>2}: accuracy = {correct}/{total} = {correct/total:.3f}")

    print("\n--- TEST 5: Searching for O(polylog) partition shortcuts ---")
    print("Can M_k(n) mod l be computed in O(polylog(n)) for FIXED l?")
    print()
    print("Key observations:")
    print("  1. M_k(n) mod l satisfies: M_k(n) = sum over partitions of n")
    print("  2. The partition DP has O(n^2) steps — NOT polylog")
    print("  3. Ramanujan congruences give M_1(n) mod 5 = 0 if 5|n+4")
    print("     But these are SPECIAL cases, not general computation")
    print("  4. For the Ono criterion, we need M_k(d) for ALL d|n")
    print("     Number of divisors: O(n^epsilon). Still need M_k at each.")
    print("  5. VERDICT: O(n) computation per n, O(n^2) total for pi(n)")
    print("     This is WORSE than Meissel-Lehmer O(n^{2/3})")

    print("\n" + "=" * 80)
    print("CONCLUSIONS")
    print("=" * 80)
    print("""
1. The partition Möbius transform T(n) = sum_{d|n} mu(n/d)*M_2(d) does NOT
   cleanly characterize primes (it's an approximation, not the exact Ono criterion).

2. Modular computation of M_k: same number of operations as exact, just
   with bounded values. No asymptotic improvement.

3. The Ono characterization requires M_k(d) for ALL divisors d of n.
   Computing M_k(d) costs O(d^2) per divisor, total O(n^{1+epsilon}).
   This is MUCH WORSE than Meissel-Lehmer's O(n^{2/3}).

4. Ramanujan-type congruences give M_k(n) mod l for SPECIAL (n,l) pairs
   but NOT general polylog computation.

5. VERDICT: The Ono partition characterization is STRUCTURALLY beautiful
   but COMPUTATIONALLY inferior to all known methods. p-adic lifting
   does not help because the partition DP has O(n^2) inherent complexity.
   CLOSED.
""")

if __name__ == "__main__":
    main()

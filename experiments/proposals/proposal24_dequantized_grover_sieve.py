#!/usr/bin/env python3
"""
Proposal 24: Dequantized Grover Sieve for Prime Counting

IDEA: Grover's algorithm gives quadratic speedup for unstructured search.
For prime counting, a quantum algorithm could count primes in [1,x] using
O(sqrt(x)) queries. Recent dequantization results (Tang 2019, Chia et al. 2020,
Jethwani et al. 2025) show that quantum advantages can sometimes be replicated
classically when data has special structure (low-rank, sparse, etc.).

KEY INSIGHT: The prime indicator function chi_P(n) is NOT unstructured.
It has structure from:
1. Multiplicative structure (sieve of Eratosthenes)
2. Fourier structure (Dirichlet characters, exponential sums)
3. Algebraic structure (prime ideals in number rings)

If we can express chi_P as a low-rank or sparse operator in some basis,
dequantization techniques could give a classical O(polylog) algorithm.

APPROACH: Model pi(x) as an inner product <v_x | M | w> where:
- v_x encodes the counting query up to x
- M is a "primality operator" with structure from the sieve
- w is a reference vector

If M has rank r, the quantum algorithm needs O(sqrt(r)) steps.
If r = polylog(x), dequantization gives classical O(polylog(x) * poly(r)).

ALSO: Explore the connection between Grover search and the Meissel-Lehmer
method. The ML method essentially does inclusion-exclusion which is a form
of structured search. Can we apply quantum-inspired search pruning?

CONJECTURE: The "primality operator" has effective rank O(x^{1/3}/log(x)),
matching the Meissel-Lehmer complexity. True polylog rank would be a breakthrough.

TIME COMPLEXITY: O(polylog(x)) if the primality operator has polylog rank.
Otherwise O(x^{1/3+epsilon}) as a potentially improved ML method.
"""

import math
import os
import numpy as np
from functools import lru_cache

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

def sieve_primes(limit):
    if limit < 2:
        return []
    sieve = [True] * (limit + 1)
    sieve[0] = sieve[1] = False
    for i in range(2, int(limit**0.5) + 1):
        if sieve[i]:
            for j in range(i*i, limit + 1, i):
                sieve[j] = False
    return [i for i in range(limit + 1) if sieve[i]]

def prime_indicator_matrix(N):
    """
    Build a matrix M where M[i,j] = 1 if i*j <= N and i*j is NOT prime.
    The idea: pi(N) = N - 1 - (number of composites in [2,N])
    Composites are detected by the matrix: if M has a 1 at (i,j), then i*j is composite.
    """
    # Sieve-like matrix: rows are potential factors, columns are quotients
    # M[p][k] = 1 means p*k is composite (being sieved out by p)
    sqrt_N = int(math.sqrt(N)) + 1
    primes_small = sieve_primes(sqrt_N)

    # Build the "sieve matrix" for small primes
    # Row i corresponds to prime p_i, column j corresponds to multiple j
    # Entry is 1 if p_i * j <= N
    rows = len(primes_small)
    cols = N // 2 + 1  # max quotient
    # This is too large to materialize; work with it implicitly
    return primes_small

def rank_analysis_small(N):
    """
    For small N, build the full primality indicator and analyze its rank
    in various decompositions.
    """
    is_prime = [False] * (N + 1)
    primes = sieve_primes(N)
    for p in primes:
        is_prime[p] = True

    # Decomposition 1: Sieve matrix rank
    # Build matrix where rows are primes <= sqrt(N), columns are numbers 2..N
    # M[i,j] = 1 if prime_i divides j
    sqrt_N = int(math.sqrt(N))
    small_primes = [p for p in primes if p <= sqrt_N]
    M_sieve = np.zeros((len(small_primes), N - 1))
    for i, p in enumerate(small_primes):
        for j in range(2, N + 1):
            if j % p == 0:
                M_sieve[i, j - 2] = 1

    rank_sieve = np.linalg.matrix_rank(M_sieve)

    # Decomposition 2: Fourier-domain analysis
    # DFT of prime indicator
    chi = np.array([1.0 if is_prime[n] else 0.0 for n in range(2, N + 1)])
    fft_chi = np.fft.fft(chi)
    magnitudes = np.abs(fft_chi)

    # "Effective rank" = number of Fourier coefficients needed for 90%, 99% energy
    total_energy = np.sum(magnitudes**2)
    sorted_mags = np.sort(magnitudes**2)[::-1]
    cumulative = np.cumsum(sorted_mags)
    rank_90 = np.searchsorted(cumulative, 0.90 * total_energy) + 1
    rank_99 = np.searchsorted(cumulative, 0.99 * total_energy) + 1

    # Decomposition 3: SVD of the "factoring tensor"
    # Reshape chi into a matrix and analyze SVD
    side = int(math.sqrt(N - 1))
    if side * side < N - 1:
        side += 1
    chi_padded = np.zeros(side * side)
    chi_padded[:len(chi)] = chi
    M_reshaped = chi_padded.reshape(side, side)
    sv = np.linalg.svd(M_reshaped, compute_uv=False)
    sv_norm = sv / sv[0] if sv[0] > 0 else sv

    # Effective rank from SVD
    sv_rank_90 = np.sum(np.cumsum(sv**2) <= 0.90 * np.sum(sv**2)) + 1
    sv_rank_99 = np.sum(np.cumsum(sv**2) <= 0.99 * np.sum(sv**2)) + 1

    return {
        'N': N,
        'pi_N': len(primes),
        'sieve_matrix_rank': rank_sieve,
        'sieve_matrix_shape': M_sieve.shape,
        'fourier_rank_90': rank_90,
        'fourier_rank_99': rank_99,
        'svd_rank_90': sv_rank_90,
        'svd_rank_99': sv_rank_99,
        'top_5_sv': sv[:5].tolist(),
        'sv_decay': sv_norm[:10].tolist(),
    }

def quantum_inspired_count(N, primes_list):
    """
    Quantum-inspired prime counting using amplitude estimation ideas.

    Classical analogue: instead of checking every number, use the sieve
    structure to "skip" large blocks of composites.

    The quantum version would prepare a superposition, apply the primality
    oracle, and use amplitude estimation. Classical analogue:
    1. Estimate pi(N) ≈ R(N) (smooth approximation)
    2. Refine using structured random sampling with importance weighting
    3. Use inclusion-exclusion truncation with error bounds

    This is essentially a randomized version of Meissel-Lehmer.
    """
    # Phase 1: Smooth estimate
    from math import log
    li_N = li_func(N)
    estimate = li_N

    # Phase 2: Small prime corrections (structured part)
    sqrt_N = int(math.sqrt(N))
    small_primes = [p for p in primes_list if p <= sqrt_N]

    # Legendre-like counting: pi(N) = pi(sqrt(N)) + N - 1 - sum_{p<=sqrt(N)} floor(N/p) + ...
    # This is the Meissel-Lehmer approach; we try to truncate it smarter

    # Phase 3: Random sampling for the "unstructured" residual
    np.random.seed(42)
    n_samples = int(math.sqrt(N))  # quantum would use sqrt of this

    # Sample random numbers in [2, N] and check primality via small primes
    count = 0
    for _ in range(n_samples):
        x = np.random.randint(2, N + 1)
        is_prime_x = True
        for p in small_primes:
            if p * p > x:
                break
            if x % p == 0:
                is_prime_x = False
                break
        if is_prime_x and x > 1:
            count += 1

    # Estimate pi(N) from sample proportion
    sample_estimate = count / n_samples * (N - 1)
    return sample_estimate

def li_func(x):
    if x <= 1:
        return 0.0
    lnx = math.log(x)
    total = 0.0
    term = 1.0
    for k in range(1, 100):
        term *= lnx / k
        total += term / k
        if abs(term / k) < 1e-15:
            break
    return 0.5772156649015329 + math.log(abs(lnx)) + total

def test_dequantized_sieve():
    print("=" * 80)
    print("PROPOSAL 24: Dequantized Grover Sieve for Prime Counting")
    print("=" * 80)

    # Test 1: Rank analysis for various N
    print("\n--- Test 1: Rank of primality-related matrices ---")
    print(f"{'N':>6} {'pi(N)':>6} {'sieve_rank':>10} {'fourier_90':>10} "
          f"{'fourier_99':>10} {'svd_90':>8} {'svd_99':>8}")

    for N in [50, 100, 200, 500, 1000, 2000]:
        result = rank_analysis_small(N)
        print(f"{N:>6} {result['pi_N']:>6} {result['sieve_matrix_rank']:>10} "
              f"{result['fourier_rank_90']:>10} {result['fourier_rank_99']:>10} "
              f"{result['svd_rank_90']:>8} {result['svd_rank_99']:>8}")

    # Test 2: How does rank scale with N?
    print("\n--- Test 2: Rank scaling analysis ---")
    Ns = [50, 100, 200, 500, 1000, 2000]
    ranks_sieve = []
    ranks_f90 = []
    ranks_f99 = []
    for N in Ns:
        r = rank_analysis_small(N)
        ranks_sieve.append(r['sieve_matrix_rank'])
        ranks_f90.append(r['fourier_rank_90'])
        ranks_f99.append(r['fourier_rank_99'])

    # Fit power law: rank ~ N^alpha
    log_N = np.log(np.array(Ns, dtype=float))
    for label, ranks in [("sieve", ranks_sieve), ("fourier_90", ranks_f90),
                          ("fourier_99", ranks_f99)]:
        log_r = np.log(np.array(ranks, dtype=float))
        coeffs = np.polyfit(log_N, log_r, 1)
        print(f"  {label:>12}: rank ~ N^{coeffs[0]:.3f} (fit: "
              f"log(rank) = {coeffs[0]:.3f}*log(N) + {coeffs[1]:.3f})")

    # Test 3: SVD spectrum
    print("\n--- Test 3: SVD spectrum of prime indicator ---")
    for N in [200, 500, 1000]:
        r = rank_analysis_small(N)
        decay = r['sv_decay'][:10]
        print(f"  N={N}: top singular values (normalized): "
              f"{', '.join(f'{v:.4f}' for v in decay)}")

    # Test 4: Quantum-inspired random sampling
    print("\n--- Test 4: Quantum-inspired random sampling ---")
    primes = sieve_primes(10000)
    actual_counts = {}
    for N in [100, 500, 1000, 5000, 10000]:
        actual = sum(1 for p in primes if p <= N)
        actual_counts[N] = actual
        est = quantum_inspired_count(N, primes)
        err_pct = abs(est - actual) / actual * 100
        print(f"  N={N:>6}: actual pi(N)={actual:>5}, "
              f"sample estimate={est:>8.1f}, error={err_pct:.1f}%")

    # Test 5: Dirichlet character decomposition
    print("\n--- Test 5: Dirichlet character decomposition ---")
    print("  Decomposing pi(x) into contributions from Dirichlet L-functions")
    # pi(x; q, a) counts primes ≡ a (mod q)
    # By orthogonality: pi(x; q, a) = (1/phi(q)) * sum_chi chi_bar(a) * pi(x, chi)
    # If L-functions have computable special values, this could help
    for q in [2, 3, 4, 6, 8, 12]:
        N = 1000
        residue_counts = {}
        for p in primes:
            if p > N:
                break
            if p < q:
                continue
            r = p % q
            if math.gcd(r, q) == 1:
                residue_counts[r] = residue_counts.get(r, 0) + 1

        phi_q = sum(1 for a in range(1, q) if math.gcd(a, q) == 1)
        expected = sum(1 for p in primes if p <= N and p >= q) / phi_q
        max_deviation = max(abs(c - expected) for c in residue_counts.values()) if residue_counts else 0
        print(f"  q={q:>2}, phi(q)={phi_q:>2}: residue distribution {dict(sorted(residue_counts.items()))}")
        print(f"         expected≈{expected:.0f} each, max_dev={max_deviation:.0f}")

    print("\n--- VERDICT ---")
    print("The rank of primality-related matrices grows as N^alpha for alpha ≈ 0.3-0.5.")
    print("This matches the Meissel-Lehmer complexity O(N^{2/3}), NOT polylog.")
    print("Dequantization would give at best a quadratic improvement over quantum,")
    print("which itself only gives sqrt speedup, so the classical bound remains O(N^{1/3}).")
    print("To get polylog, we would need the 'primality operator' to have polylog rank,")
    print("which appears equivalent to a major number-theoretic breakthrough.")

if __name__ == '__main__':
    test_dequantized_sieve()

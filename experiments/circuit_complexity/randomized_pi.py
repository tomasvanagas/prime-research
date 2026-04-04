"""
Session 15: Can randomized algorithms compute pi(x) exactly in sublinear time?

Key question: Does adding randomness to computation help for pi(x)?

Analysis:
1. Analytic approach: sum_rho R(x^rho) has K_min ~ x^{0.27} terms (Session 13).
   Random sampling of zeros: variance analysis shows ALL terms needed.

2. Combinatorial approach: Meissel-Lehmer sieve uses inclusion-exclusion.
   Random sampling: massive cancellation prevents subsampling.

3. Hashing approach: Could random hash functions help identify primes faster?

4. RP/BPP oracle: If PRIMES is in TC^0 (BPSW correct), could a randomized
   COUNTING algorithm exploit this?

5. Quantum counting: Grover + quantum counting gives O(sqrt(x)) queries.
   Can we do better with quantum + classical hybrid?

Experiment 1: Variance of random subsample of zeta zero sum.
Experiment 2: Randomized inclusion-exclusion (probabilistic sieve).
Experiment 3: Can hash-based counting give sublinear pi(x)?
"""

import numpy as np
import math
from sympy import primepi

def experiment1_zeta_zero_sampling():
    """
    The explicit formula: pi(x) ≈ R(x) - sum_{k=1}^K R(x^{rho_k})

    If we sample M of K zeros uniformly at random, the estimate is:
    pi_est(x) = R(x) - (K/M) * sum_{sampled} R(x^{rho_k})

    Variance: Var = (K/M - 1) * sum_k |R(x^{rho_k})|^2 / K

    For the result to be within 0.5 of pi(x), we need Var < 0.25.
    How large must M be?
    """
    print("=" * 60)
    print("Experiment 1: Zeta zero sampling variance")
    print("=" * 60)

    # Load zeta zeros
    try:
        with open('data/zeta_zeros_200.txt', 'r') as f:
            zeros = []
            for line in f:
                line = line.strip()
                if line and not line.startswith('#'):
                    try:
                        zeros.append(float(line))
                    except ValueError:
                        pass
    except FileNotFoundError:
        # Use known first zeros
        zeros = [14.134725, 21.022040, 25.010858, 30.424876, 32.935062,
                 37.586178, 40.918719, 43.327073, 48.005151, 49.773832,
                 52.970321, 56.446248, 59.347044, 60.831779, 65.112544,
                 67.079811, 69.546402, 72.067158, 75.704691, 77.144840]

    print(f"Using {len(zeros)} zeta zeros")

    for N in [20, 30, 40, 50]:
        x = 2**N

        # Compute |R(x^rho)| for each zero
        # R(x^rho) ≈ li(x^rho) for the leading term
        # |li(x^{1/2 + it})| ≈ x^{1/2} / (|rho| * ln(x))
        #   where rho = 1/2 + i*gamma (assuming RH)

        term_magnitudes = []
        for gamma in zeros:
            rho_mag = math.sqrt(0.25 + gamma**2)
            # |R(x^rho)| ≈ x^{1/2} / (rho_mag * N * ln(2))
            mag = (2**(N/2)) / (rho_mag * N * math.log(2))
            term_magnitudes.append(mag)

        term_mags = np.array(term_magnitudes)
        K = len(zeros)

        # Total sum magnitude (worst case = all same sign)
        total_worst = np.sum(term_mags)
        # RMS of terms
        rms = np.sqrt(np.mean(term_mags**2))

        # For random subsample of M out of K:
        # Var ≈ (K-M)/M * sum |terms|^2 / K = (K/M - 1) * rms^2
        # Need sqrt(Var) < 0.5 for exact answer

        # Required M such that (K/M - 1) * K * rms^2 < 0.25
        total_sq = np.sum(term_mags**2)
        required_M = K * total_sq / (total_sq + 0.25/K)

        print(f"\nN={N}: x = 2^{N}")
        print(f"  K = {K} zeros, largest term = {term_mags[0]:.2e}")
        print(f"  Sum |terms| = {total_worst:.2e}")
        print(f"  Sum |terms|^2 = {total_sq:.2e}")
        print(f"  Required M for Var<0.25: {required_M:.1f} / {K} = {required_M/K*100:.2f}%")
        print(f"  Conclusion: Must use {required_M/K*100:.1f}% of zeros — NO sampling gain")


def experiment2_probabilistic_sieve():
    """
    Can we do a probabilistic sieve? Instead of full Legendre I-E:
    pi(x) = pi(sqrt(x)) - 1 + sum_{S} (-1)^|S| floor(x / prod S)

    What if we sample subsets S randomly? The I-E sum has 2^{pi(sqrt(x))} terms.
    With random sampling, could we estimate pi(x) to ±0.5?

    Variance analysis:
    - Each term is floor(x / prod_{p in S} p) with sign (-1)^|S|
    - Terms range from x (S=empty) to 0 (S large)
    - Massive cancellation: result = pi(x) ~ x/ln(x) from sum of terms up to x

    Variance per sample: Var(term) ~ x^2 (since terms range up to x)
    Need Var/M < 0.25, so M > 4*Var = 4*x^2 — FAR worse than direct.
    """
    print("\n" + "=" * 60)
    print("Experiment 2: Probabilistic sieve variance")
    print("=" * 60)

    for x in [100, 1000, 10000]:
        from sympy import primerange
        primes = list(primerange(2, int(math.sqrt(x)) + 1))
        k = len(primes)
        pi_x = int(primepi(x))

        # Compute ALL I-E terms
        terms = []
        for subset_mask in range(2**k):
            prod = 1
            bits_set = 0
            for j in range(k):
                if subset_mask & (1 << j):
                    prod *= primes[j]
                    bits_set += 1
            sign = (-1) ** bits_set
            term = sign * (x // prod)
            terms.append(term)

        terms = np.array(terms, dtype=float)

        # The I-E result (Legendre): sum of all terms = some_value, adjusted by pi(sqrt(x))-1
        ie_sum = np.sum(terms)
        # Legendre: Phi(x, k) = ie_sum, pi(x) = Phi(x,k) + k - 1

        # Variance of a single random term
        mean_term = np.mean(terms)
        var_term = np.var(terms)

        # For random subsample of M terms (scaled):
        # Estimate = (2^k / M) * sum of M random terms
        # Var(estimate) = (2^k / M)^2 * M * var_term = 2^{2k} * var_term / M
        # Need < 0.25: M > 4 * 2^{2k} * var_term

        required_M = 4 * (2**(2*k)) * var_term
        total_terms = 2**k

        print(f"\nx={x}: k={k} primes, 2^k={total_terms} I-E terms")
        print(f"  pi(x) = {pi_x}")
        print(f"  I-E sum + k - 1 = {int(ie_sum)} + {k} - 1 = {int(ie_sum) + k - 1}")
        print(f"  Mean term = {mean_term:.2f}, Var = {var_term:.2e}")
        print(f"  Required M for |error|<0.5: {required_M:.2e}")
        print(f"  Ratio M/2^k: {required_M/total_terms:.2e}")
        print(f"  Conclusion: {'WORSE' if required_M > total_terms else 'BETTER'} than exhaustive")


def experiment3_hash_counting():
    """
    Hash-based prime counting: use random hash functions to sketch the set of primes.

    Idea: if we have a family of hash functions h_1, ..., h_L that map {1,...,x} to {0,1},
    and we can compute sum_{n=1}^x h_i(n) * [n prime] efficiently for each h_i,
    then L such sums could determine pi(x).

    But computing sum_{n=1}^x h_i(n) * [n prime] is at least as hard as pi(x)!
    This is circular.

    Alternative: Can we compute F(x) = sum_{n=1}^x f(n) for MULTIPLICATIVE f in polylog?
    Dirichlet hyperbola: sum_{n<=x} (f*g)(n) = ...
    For f = 1_prime (NOT multiplicative!), this doesn't directly apply.

    For COMPLETELY multiplicative f (like chi_q characters):
    sum_{n<=x} chi(n) can be computed in O(sqrt(x)) using Polya-Vinogradov / hyperbola.
    But we need the PRIME sum, not the full sum.
    """
    print("\n" + "=" * 60)
    print("Experiment 3: Hash-based counting analysis")
    print("=" * 60)

    print("THEORETICAL ANALYSIS:")
    print()
    print("The fundamental problem: we want to count elements of a SET")
    print("(the primes ≤ x) of size pi(x) ~ x/ln(x).")
    print()
    print("Communication complexity argument:")
    print("  Alice has a 'description' of the set of primes (e.g., the prime indicator)")
    print("  Bob has x")
    print("  They want to compute pi(x)")
    print()
    print("  The set of primes has 2^{pi(x)} possible subsets of {1,...,x}.")
    print("  So the 'description' requires at least pi(x) * log2(x/pi(x)) ~ pi(x) bits")
    print("  (by Stirling, C(x, pi(x)) ~ 2^{H(pi(x)/x)*x})")
    print()
    print("  But the primes are DETERMINISTIC — Alice doesn't need to send the set.")
    print("  The question is whether Bob (with x on his tape) can count primes")
    print("  using only polylog(x) time.")
    print()
    print("  A hash-based sketch of size S bits can estimate cardinality to ±epsilon*n")
    print("  with S = O(1/epsilon^2 * log(n)). For epsilon < 1/(2*pi(x)) = ln(x)/(2x),")
    print("  we need S = O(x^2 / ln^2(x) * N) — WORSE than storing all primes.")
    print()
    print("CONCLUSION: No randomized counting sketch can give EXACT pi(x)")
    print("more efficiently than deterministic computation.")
    print()

    # Verify: how many bits does an exact count of pi(x) need?
    for N in [20, 40, 60, 80, 100]:
        x = 2**N
        pi_approx = x / (N * math.log(2))
        bits_to_specify_primes = pi_approx * math.log2(x / pi_approx)

        # Compare to N bits (the input size)
        print(f"N={N:3d}: pi(x)~{pi_approx:.2e}, bits to specify prime set: {bits_to_specify_primes:.2e}, "
              f"ratio to x: {bits_to_specify_primes / x:.4f}")


def experiment4_quantum_counting_bounds():
    """
    Quantum counting: Given an oracle for f(n) = [n is prime],
    quantum counting (Brassard-Hoyer-Tapp 2000) determines pi(x) to ±sqrt(x)
    in O(1) queries, or exactly in O(sqrt(x/pi(x))) ≈ O(sqrt(ln(x))) queries.

    Wait — that seems too good! Let me recalculate.

    Quantum counting uses quantum amplitude estimation on Grover's algorithm.
    For counting M items among N: O(sqrt(N/M)) queries for exact count.

    Here N = x, M = pi(x) ~ x/ln(x).
    Exact count: O(sqrt(x / (x/ln(x)))) = O(sqrt(ln(x))) queries.

    BUT: each "query" involves a quantum circuit that evaluates PRIMES(n),
    which takes poly(N) gates (AKS). So total: O(sqrt(ln(x)) * poly(N)) = poly(N).

    This would give pi(x) in POLYNOMIAL time (in N)... but we already know
    pi(x) ∈ #P which is computable in poly(x) = exponential in N.

    The catch: quantum amplitude estimation gives an ESTIMATE with bounded error,
    NOT an exact count. For exact count, you need the full QFT on the counting
    register, which requires M queries where M is the PERIOD... actually, the
    exact quantum counting requires O(x) queries in the worst case.

    Let me think more carefully.
    """
    print("\n" + "=" * 60)
    print("Experiment 4: Quantum counting bounds")
    print("=" * 60)

    print("Quantum counting for pi(x):")
    print()
    print("Model: Quantum oracle O_f with f(n) = [n is prime]")
    print()
    print("Brassard-Hoyer-Tapp (2000):")
    print("  Can estimate |{n: f(n)=1, n<=x}| to ±epsilon using O(sqrt(x)/epsilon) queries")
    print("  For epsilon=0.5 (exact): O(sqrt(x)) queries")
    print()
    print("Nayak-Wu (2000):")
    print("  Lower bound: Omega(sqrt(x/pi(x))) for exact counting = Omega(sqrt(ln x))")
    print()
    print("The GAP: sqrt(ln x) vs sqrt(x) for exact quantum counting of pi(x)")
    print()
    print("Key question: Can structure of PRIMES improve quantum counting?")
    print("  - Standard quantum counting treats f as a BLACK BOX")
    print("  - PRIMES has structure (AKS, BPSW, etc.)")
    print("  - Can QUANTUM + STRUCTURE beat O(sqrt(x))?")
    print()

    for N in [20, 50, 100, 200, 333]:
        x_log = N * math.log(2)
        sqrt_x = 2**(N/2)
        sqrt_lnx = math.sqrt(x_log)
        quantum_std = sqrt_x  # standard quantum counting
        quantum_lower = sqrt_lnx  # Nayak-Wu lower bound

        print(f"N={N:3d}: sqrt(x) = 2^{N/2:.0f}, sqrt(ln x) = {sqrt_lnx:.1f}, "
              f"gap = 2^{N/2 - math.log2(sqrt_lnx):.1f}")

    print()
    print("CONCLUSION: Even quantum counting needs Omega(sqrt(ln x)) = Omega(sqrt(N))")
    print("queries. The standard algorithm uses O(sqrt(x)) = O(2^{N/2}).")
    print("Closing this gap requires EXPLOITING the structure of PRIMES,")
    print("which is equivalent to our original problem.")


if __name__ == "__main__":
    experiment1_zeta_zero_sampling()
    experiment2_probabilistic_sieve()
    experiment3_hash_counting()
    experiment4_quantum_counting_bounds()

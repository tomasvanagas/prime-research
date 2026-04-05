"""
Approximate polynomial degree of the prime indicator over GF(p).

The Razborov-Smolensky method: if f has approximate degree omega(polylog N)
over GF(p) for ALL primes p, then f ∉ ACC^0.

If isPrime has approximate degree O(polylog N) over some GF(p),
it could potentially be in ACC^0, and pi(x) might then be in #ACC^0.

Approximate degree over GF(p): minimum d such that there exists a polynomial
g over GF(p) of degree d with g(x) = f(x) for at least (1 - 1/(2p)) fraction
of inputs.

We compute this for the prime indicator function on {0,1}^N for N = 4..14.
"""

import numpy as np
from sympy import isprime
from itertools import combinations
import sys

def get_truth_table(N):
    """Get the prime indicator for all N-bit integers."""
    return np.array([1 if isprime(x) else 0 for x in range(2**N)], dtype=np.int8)

def gf_p_multilinear_coeffs(f, N, p):
    """
    Compute the unique multilinear polynomial over GF(p) that equals f
    on all of {0,1}^N.

    Over GF(p), any function {0,1}^N -> GF(p) has a unique multilinear
    polynomial representation (the ANF over GF(p)).

    Returns dict: frozenset of variable indices -> coefficient mod p
    """
    n_points = 2**N
    # Use Moebius inversion on the Boolean lattice
    # ANF coefficient of monomial S = sum over T subset of S of (-1)^{|S|-|T|} f(T)
    # Over GF(p), this is sum over T subset of S of (-1)^{|S|-|T|} f(T) mod p

    coeffs = {}
    f_mod = f % p

    # For each subset S of [N], compute the coefficient
    # We'll limit to subsets of manageable size
    max_degree = N  # compute all

    # Zeta-Moebius transform approach
    # Transform f to get ANF coefficients
    g = f_mod.copy().astype(np.int64)

    for i in range(N):
        step = 1 << i
        for j in range(n_points):
            if j & step:
                g[j] = (g[j] - g[j ^ step]) % p

    # g[S] is now the ANF coefficient for subset S (where S is the bitmask)
    return g

def compute_degree_distribution(anf, N, p):
    """Compute the distribution of ANF coefficients by degree."""
    degree_counts = {}
    degree_nonzero = {}

    for S in range(2**N):
        deg = bin(S).count('1')
        if deg not in degree_counts:
            degree_counts[deg] = 0
            degree_nonzero[deg] = 0
        degree_counts[deg] += 1
        if anf[S] % p != 0:
            degree_nonzero[deg] += 1

    return degree_counts, degree_nonzero

def compute_approx_degree(f, N, p, epsilon):
    """
    Find the minimum degree d such that there exists a degree-d polynomial
    over GF(p) that agrees with f on at least (1-epsilon) fraction of {0,1}^N.

    Method: for each degree d from 0 to N, find the best degree-d polynomial
    approximation to f over GF(p) by solving a system or using LP.

    For small N, we can use exhaustive search over all possible sets of monomials.
    For larger N, we use a greedy approach.
    """
    n_points = 2**N
    f_mod = f % p

    # Build the monomial evaluation matrix
    # Each column corresponds to a monomial (subset of variables)
    # Each row corresponds to an input point

    # For degree d: include all monomials of degree <= d
    # Number of such monomials: sum_{k=0}^{d} C(N, k)

    results = {}

    for d in range(N + 1):
        # Build monomial evaluation matrix for degree <= d
        monomial_indices = []
        for S in range(2**N):
            if bin(S).count('1') <= d:
                monomial_indices.append(S)

        n_monomials = len(monomial_indices)
        if n_monomials > 5000:
            # Too many monomials for direct solve, skip
            results[d] = None
            continue

        # Build matrix A where A[x, j] = prod_{i in S_j} x_i for x in {0,1}^N
        A = np.zeros((n_points, n_monomials), dtype=np.int64)
        for j, S in enumerate(monomial_indices):
            for x in range(n_points):
                # Monomial S evaluated at point x: 1 if all bits of S are set in x
                A[x, j] = 1 if (x & S) == S else 0

        # Over GF(p), find c minimizing number of disagreements
        # This is: minimize |{x : (A*c)(x) != f(x) mod p}|
        # For GF(2), this is a syndrome decoding problem (NP-hard in general)
        # But for small problems, we can use least squares over reals and round

        # Use the GF(p) ANF: restrict to degree-d terms
        anf = gf_p_multilinear_coeffs(f, N, p)

        # Compute residual when truncating to degree d
        reconstructed = np.zeros(n_points, dtype=np.int64)
        for j, S in enumerate(monomial_indices):
            coeff = int(anf[S]) % p
            if coeff != 0:
                for x in range(n_points):
                    if (x & S) == S:
                        reconstructed[x] = (reconstructed[x] + coeff) % p

        # Count agreements
        agreements = np.sum((reconstructed % p) == (f_mod % p))
        agreement_frac = agreements / n_points
        disagreement_frac = 1 - agreement_frac

        results[d] = {
            'n_monomials': n_monomials,
            'agreements': int(agreements),
            'agreement_frac': agreement_frac,
            'disagreement_frac': disagreement_frac,
            'sufficient': disagreement_frac <= epsilon
        }

        if disagreement_frac <= epsilon:
            return d, results

    return N, results

def main():
    print("=" * 80)
    print("APPROXIMATE POLYNOMIAL DEGREE OF PRIME INDICATOR")
    print("=" * 80)

    primes_for_test = [2, 3, 5, 7]

    for N in [4, 6, 8, 10, 12, 14]:
        if N > 14:
            print(f"\nN={N}: skipped (too large)")
            continue

        print(f"\n{'='*60}")
        print(f"N = {N} (x ∈ {{0, ..., {2**N - 1}}})")
        print(f"{'='*60}")

        f = get_truth_table(N)
        n_primes = np.sum(f)
        density = n_primes / 2**N
        print(f"  Primes: {n_primes}/{2**N} (density {density:.4f})")

        for p in primes_for_test:
            print(f"\n  --- GF({p}) ---")

            # Compute ANF
            anf = gf_p_multilinear_coeffs(f, N, p)

            # Degree distribution
            _, degree_nonzero = compute_degree_distribution(anf, N, p)

            max_deg = max(d for d, c in degree_nonzero.items() if c > 0)
            total_nonzero = sum(degree_nonzero.values())

            print(f"  Exact degree: {max_deg}")
            print(f"  Nonzero coefficients: {total_nonzero} / {2**N}")

            # Show degree distribution
            print(f"  {'deg':>5} {'nonzero':>8} {'total':>8} {'frac':>8}")
            for d in sorted(degree_nonzero.keys()):
                total_at_d = 1
                for j in range(d):
                    total_at_d = total_at_d * (N - j) // (j + 1)
                frac = degree_nonzero[d] / total_at_d if total_at_d > 0 else 0
                if degree_nonzero[d] > 0:
                    print(f"  {d:>5} {degree_nonzero[d]:>8} {total_at_d:>8} {frac:>8.3f}")

            # Compute approximate degree for epsilon = 1/(2p)
            epsilon = 1.0 / (2 * p)
            print(f"\n  Approx degree (eps={epsilon:.3f}):")

            approx_d, details = compute_approx_degree(f, N, p, epsilon)

            for d in sorted(details.keys()):
                if details[d] is not None:
                    info = details[d]
                    marker = " ← SUFFICIENT" if info['sufficient'] else ""
                    print(f"    deg {d:>3}: agreement={info['agreement_frac']:.4f}, "
                          f"disagree={info['disagreement_frac']:.4f}{marker}")
                    if info['sufficient']:
                        break

            print(f"  Approximate degree over GF({p}): {approx_d}")

    print("\n" + "=" * 80)
    print("SCALING ANALYSIS")
    print("=" * 80)

    print("\nApprox degree vs N for each GF(p):")
    print(f"{'N':>5}", end="")
    for p in primes_for_test:
        print(f"  GF({p})", end="")
    print("  exact")
    print("-" * 50)

    for N in [4, 6, 8, 10, 12]:
        f = get_truth_table(N)
        print(f"{N:>5}", end="")

        for p in primes_for_test:
            epsilon = 1.0 / (2 * p)
            ad, _ = compute_approx_degree(f, N, p, epsilon)
            print(f"  {ad:>5}", end="")

        # Exact degree
        anf2 = gf_p_multilinear_coeffs(f, N, 2)
        _, dnz2 = compute_degree_distribution(anf2, N, 2)
        ed = max(d for d, c in dnz2.items() if c > 0)
        print(f"  {ed:>5}")

    print("\nIf approx degree grows as O(N), isPrime ∉ ACC^0 (conditionally)")
    print("If approx degree grows as O(polylog N), isPrime might be in ACC^0")

    print("\n" + "=" * 80)
    print("CONCLUSIONS")
    print("=" * 80)

if __name__ == "__main__":
    main()

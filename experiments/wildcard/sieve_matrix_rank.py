#!/usr/bin/env python3
"""
Experiment: Rank structure of the sieve matrix.

Key insight: The sieve of Eratosthenes can be viewed as a binary matrix problem.
Define M[n][j] = 1 if prime_j divides n, for n=1..x, primes p_j <= y.

Then: pi(x) - pi(y) + 1 = |{n <= x : row n is all zeros}|
     = x - rank stuff...

Questions:
1. What is the rank of M over GF(2)? Over R? Over GF(p)?
2. Does the rank grow much slower than min(x, pi(y))?
3. Can we compute |{zero rows}| from the SVD/rank structure without examining all rows?
4. Is there a polynomial-size "sketch" of M that preserves the zero-row count?
5. What's the singular value spectrum? If it decays fast, low-rank approx works.

Also: The "transfer matrix" approach - computing the sieve as a product of
per-prime projection matrices.
"""

import numpy as np
from sympy import primerange, isprime, factorint
from collections import Counter
import time

def sieve_matrix(x, y):
    """Build binary matrix M[n][j] = 1 if prime_j | n, for n=1..x, primes <= y."""
    primes = list(primerange(2, y+1))
    M = np.zeros((x, len(primes)), dtype=np.int8)
    for j, p in enumerate(primes):
        for n in range(p, x+1, p):
            M[n-1][j] = 1
    return M, primes

def count_zero_rows(M):
    """Count rows that are all zeros (these correspond to 1 and numbers coprime to all sieving primes)."""
    return np.sum(np.all(M == 0, axis=1))

def analyze_rank_structure(x_values, y_func):
    """Analyze how rank of sieve matrix grows with x."""
    print("=" * 70)
    print("SIEVE MATRIX RANK ANALYSIS")
    print("=" * 70)

    for x in x_values:
        y = y_func(x)
        M, primes = sieve_matrix(x, y)

        # Rank over R
        rank_R = np.linalg.matrix_rank(M.astype(float))

        # Rank over GF(2) - approximate via reduced row echelon
        rank_GF2 = gf2_rank(M)

        # SVD spectrum
        if min(M.shape) > 1:
            sv = np.linalg.svd(M.astype(float), compute_uv=False)
            sv_normalized = sv / sv[0] if sv[0] > 0 else sv
            # How many singular values capture 99% of energy?
            energy = np.cumsum(sv**2) / np.sum(sv**2) if np.sum(sv**2) > 0 else np.zeros_like(sv)
            k99 = np.searchsorted(energy, 0.99) + 1
        else:
            sv_normalized = []
            k99 = 0

        # Zero rows = numbers coprime to primorial(y)
        n_zero = count_zero_rows(M)
        expected_zero = x * np.prod([1 - 1/p for p in primes])  # Euler product

        print(f"\nx={x}, y={y}, primes={len(primes)}")
        print(f"  Matrix shape: {M.shape}")
        print(f"  Rank over R:    {rank_R}")
        print(f"  Rank over GF(2): {rank_GF2}")
        print(f"  Max possible rank: {min(M.shape)}")
        print(f"  Rank ratio (R):  {rank_R/min(M.shape):.3f}")
        print(f"  Zero rows: {n_zero} (expected by Euler: {expected_zero:.1f})")
        print(f"  k for 99% energy: {k99}/{len(primes)}")
        if len(sv_normalized) > 0:
            print(f"  Top 5 singular values (normalized): {sv_normalized[:5]}")
            print(f"  Bottom 5 singular values: {sv[max(0,len(sv)-5):]}")

def gf2_rank(M):
    """Compute rank of binary matrix over GF(2)."""
    A = M.copy() % 2
    rows, cols = A.shape
    rank = 0
    for col in range(cols):
        # Find pivot
        pivot = None
        for row in range(rank, rows):
            if A[row, col] == 1:
                pivot = row
                break
        if pivot is None:
            continue
        # Swap
        A[[rank, pivot]] = A[[pivot, rank]]
        # Eliminate
        for row in range(rows):
            if row != rank and A[row, col] == 1:
                A[row] = (A[row] + A[rank]) % 2
        rank += 1
    return rank

def transfer_matrix_experiment(x_max, y):
    """
    The sieve as a product of transfer matrices.

    For each prime p, define the projection: remove multiples of p.
    The sieve applies these projections sequentially.

    Key question: does the "transfer matrix" product have low rank?
    """
    print("\n" + "=" * 70)
    print("TRANSFER MATRIX EXPERIMENT")
    print("=" * 70)

    primes = list(primerange(2, y+1))

    # State: indicator vector of length x_max
    # Each prime p acts as: v[n] *= (1 - [p|n])

    # Instead of full vectors, represent the sieve as a product of diagonal matrices
    # D_p = diag(1 - [p|1], 1 - [p|2], ..., 1 - [p|x_max])
    # Sieve result = D_{p_k} ... D_{p_2} D_{p_1} * ones

    # The product of diagonal matrices is diagonal, so rank = number of survivors
    # This doesn't help directly.

    # Better: work modulo some base B.
    # For each prime p, the pattern of divisibility by p in [1..B*p] repeats with period p.
    # If B = lcm of all primes <= y (= primorial(y)), the pattern repeats with period B.
    # But B = primorial(y) is huge.

    # Alternative: use CRT. Work mod p for each sieving prime.
    # The sieve state mod p is determined by n mod p.
    # This gives a "tensor product" structure: state in Z/p1 x Z/p2 x ... x Z/pk

    # The number of coprime residues = phi(primorial(y)) = primorial(y) * prod(1-1/p)
    # And pi(x) - pi(y) + 1 ≈ x * phi(primorial(y)) / primorial(y)

    # The key insight: via CRT, the "state space" is Z/primorial(y).
    # A number n is coprime to primorial(y) iff (n mod p1, n mod p2, ..., n mod pk)
    # avoids (0, *, ..., *), (*, 0, ..., *), etc.

    # This is an inclusion-exclusion on the tensor factors.
    # The "rank" of this tensor is the number of surviving residue classes = phi(primorial)

    # But can we compute x * phi(primorial(y)) / primorial(y) + correction more efficiently?

    # The correction to count integers in [1,x] coprime to primorial(y):
    # S(x, y) = sum_{d | primorial(y)} mu(d) * floor(x/d)

    # This has 2^pi(y) terms. Can we group them by size of d?

    # Group by omega(d) (number of prime factors of d):
    # S(x,y) = sum_{k=0}^{pi(y)} (-1)^k sum_{S subset primes, |S|=k} floor(x / prod(S))

    # For large x, floor(x/d) ≈ x/d, so the sum ≈ x * prod(1-1/p) = known.
    # The fractional parts {x/d} are the hard part.

    # Idea: can we compute sum of {x/d} mu(d) efficiently?

    x = x_max
    print(f"\nAnalyzing sieve sum structure for x={x}, y={y}")
    print(f"Sieving primes: {primes}")
    print(f"Number of Mobius terms: 2^{len(primes)} = {2**len(primes)}")

    # Compute the Mobius sum exactly
    from itertools import combinations

    main_term = 0
    fractional_sum = 0
    terms_by_size = {}

    for k in range(len(primes) + 1):
        sign = (-1) ** k
        term_sum = 0
        frac_sum = 0
        for S in combinations(primes, k):
            d = 1
            for p in S:
                d *= p
            floor_xd = x // d
            frac_xd = x / d - floor_xd
            term_sum += sign * floor_xd
            frac_sum += sign * frac_xd
        terms_by_size[k] = (term_sum, frac_sum)
        main_term += term_sum
        fractional_sum += frac_sum

    euler_approx = x * np.prod([1 - 1/p for p in primes])
    actual = sum(1 for n in range(1, x+1) if all(n % p != 0 for p in primes))

    print(f"\nExact count (brute): {actual}")
    print(f"Mobius sum: {main_term}")
    print(f"Euler product approx: {euler_approx:.2f}")
    print(f"Fractional correction: {fractional_sum:.6f}")
    print(f"\nTerms by subset size k:")
    for k, (ts, fs) in terms_by_size.items():
        print(f"  k={k}: integer_part_sum={ts}, fractional_sum={fs:.6f}")

    # KEY QUESTION: Does the fractional sum have structure?
    # If it's always small and predictable, we can get exact answers cheaply.
    print(f"\n  Total fractional sum magnitude: {abs(fractional_sum):.6f}")
    print(f"  This must be close to an integer (it's {round(fractional_sum)} + {fractional_sum - round(fractional_sum):.6f})")

def fractional_part_analysis(x_range, y):
    """
    Deep analysis of the fractional parts in the Mobius sum.

    S(x,y) = sum_d mu(d) floor(x/d) = sum_d mu(d) (x/d - {x/d})
           = x * prod(1-1/p) - sum_d mu(d) {x/d}

    The fractional correction F(x) = sum_d mu(d) {x/d} must satisfy:
    S(x,y) = round(x * prod(1-1/p) - F(x))

    If F(x) is smooth/predictable, this gives a fast algorithm!
    """
    print("\n" + "=" * 70)
    print("FRACTIONAL PART ANALYSIS")
    print("=" * 70)

    from itertools import combinations
    primes = list(primerange(2, y+1))

    # Generate all squarefree divisors of primorial(y) with their Mobius values
    divisors = []
    for k in range(len(primes) + 1):
        for S in combinations(primes, k):
            d = 1
            for p in S:
                d *= p
            divisors.append((d, (-1)**k))

    euler_prod = np.prod([1 - 1/p for p in primes])

    F_values = []
    S_values = []

    for x in x_range:
        F = sum(mu * (x/d - x//d) for d, mu in divisors)
        S = sum(mu * (x // d) for d, mu in divisors)
        F_values.append(F)
        S_values.append(S)

    F_values = np.array(F_values)
    x_arr = np.array(x_range)

    print(f"\nPrimes sieved: {primes}")
    print(f"Euler product: {euler_prod:.6f}")
    print(f"Number of divisors: {len(divisors)}")

    # Statistics of F(x)
    print(f"\nFractional correction F(x) statistics:")
    print(f"  Mean:   {np.mean(F_values):.4f}")
    print(f"  Std:    {np.std(F_values):.4f}")
    print(f"  Min:    {np.min(F_values):.4f}")
    print(f"  Max:    {np.max(F_values):.4f}")
    print(f"  Range:  {np.max(F_values) - np.min(F_values):.4f}")

    # Autocorrelation of F(x)
    F_centered = F_values - np.mean(F_values)
    if np.std(F_centered) > 0:
        autocorr = np.correlate(F_centered, F_centered, mode='full')
        autocorr = autocorr[len(autocorr)//2:]
        autocorr /= autocorr[0]

        print(f"\n  Autocorrelation at lags 1-10:")
        for lag in range(1, min(11, len(autocorr))):
            print(f"    lag {lag}: {autocorr[lag]:.4f}")

    # Periodicity: F(x) should have period primorial(y) up to boundary effects
    prim_y = 1
    for p in primes:
        prim_y *= p
    print(f"\n  Primorial({y}) = {prim_y}")
    print(f"  Expected periodicity in F(x): period {prim_y}")

    # Check: is F(x) = F(x + primorial) exactly?
    if len(x_range) > prim_y + 10:
        diffs = []
        for i in range(min(20, len(x_range) - prim_y)):
            x1 = x_range[i]
            x2 = x1 + prim_y
            if x2 in x_range:
                idx2 = x_range.index(x2)
                diffs.append(F_values[idx2] - F_values[i])
        if diffs:
            print(f"  F(x+prim) - F(x) for first cases: {diffs[:5]}")

    # KEY TEST: Can F(x) be predicted from a short formula?
    # Try polynomial fit
    for deg in [1, 2, 3, 5]:
        coeffs = np.polyfit(x_arr, F_values, deg)
        pred = np.polyval(coeffs, x_arr)
        residual = np.max(np.abs(F_values - pred))
        rmse = np.sqrt(np.mean((F_values - pred)**2))
        print(f"\n  Polynomial degree {deg}: max_residual={residual:.4f}, rmse={rmse:.4f}")


def rank_growth_scaling():
    """How does the effective rank scale?"""
    print("\n" + "=" * 70)
    print("RANK GROWTH SCALING")
    print("=" * 70)

    # Test: for y = x^{1/3}, how does rank of sieve matrix scale with x?
    results = []
    for x in [30, 60, 100, 200, 500, 1000, 2000]:
        y = int(x ** (1/3)) + 1
        if y < 2:
            y = 2
        M, primes = sieve_matrix(x, y)
        if M.shape[1] == 0:
            continue

        rank_R = np.linalg.matrix_rank(M.astype(float))

        sv = np.linalg.svd(M.astype(float), compute_uv=False)
        if sv[0] > 0:
            energy = np.cumsum(sv**2) / np.sum(sv**2)
            k90 = np.searchsorted(energy, 0.90) + 1
            k99 = np.searchsorted(energy, 0.99) + 1
        else:
            k90 = k99 = 0

        n_primes = len(primes)
        results.append((x, y, n_primes, rank_R, k90, k99))
        print(f"x={x:5d}, y={y:2d}, #primes={n_primes}, rank={rank_R}, k90={k90}, k99={k99}")

    if len(results) > 2:
        xs = [r[0] for r in results]
        ranks = [r[3] for r in results]
        # Fit rank ~ x^alpha
        log_xs = np.log(xs)
        log_ranks = np.log([max(r, 1) for r in ranks])
        alpha = np.polyfit(log_xs, log_ranks, 1)[0]
        print(f"\nRank scaling: rank ~ x^{alpha:.3f}")


def nullspace_structure():
    """
    Analyze the nullspace of the sieve matrix.
    If the nullspace has nice structure, we can count zero rows efficiently.
    """
    print("\n" + "=" * 70)
    print("NULLSPACE STRUCTURE")
    print("=" * 70)

    for x in [100, 500, 1000]:
        y = int(x ** 0.5) + 1
        M, primes = sieve_matrix(x, y)

        # Compute SVD
        U, S, Vt = np.linalg.svd(M.astype(float), full_matrices=True)

        # Null space of M^T: left singular vectors with zero singular values
        # These correspond to "hidden symmetries" in the sieve

        nullity = M.shape[0] - np.linalg.matrix_rank(M.astype(float))
        n_zero_rows = count_zero_rows(M)
        n_survivors = n_zero_rows  # Numbers coprime to primorial(y), including 1
        actual_pi = sum(1 for n in range(2, x+1) if isprime(n))

        print(f"\nx={x}, y={y}")
        print(f"  Matrix: {M.shape}, Rank: {M.shape[0]-nullity}, Nullity: {nullity}")
        print(f"  Zero rows (coprime to P(y)): {n_zero_rows}")
        print(f"  pi(x)={actual_pi}, pi(y)={sum(1 for p in primes)}")
        print(f"  pi(x)-pi(y)+1 should equal coprime count for y=sqrt(x): {actual_pi - len(primes) + 1}")
        print(f"  Note: coprime count includes non-primes > y that are coprime to P(y)")

        # The zero rows of M are exactly the numbers n where gcd(n, primorial(y)) = 1
        # This includes 1, primes > y, and products of primes > y (semiprimes etc.)


def recursive_sieve_idea():
    """
    THE BIG IDEA: Recursive sieve via hierarchical matrix factorization.

    Standard sieve: apply projection for each prime sequentially.
    Meissel-Lehmer: split into cases based on prime size, recurse.

    New idea: Can we factorize the sieve into a LOGARITHMIC depth tree?

    Level 0: Full range [1, x]
    Level 1: Split primes into two groups, sieve each independently, combine
    Level 2: Further split, etc.

    The combination step is the key. If we sieve by primes in set A and set B
    independently, the combined count is:

    S(x, A ∪ B) = S(x, A) + S(x, B) - x + (overcounting correction)

    Wait, that's inclusion-exclusion again. The correction is exactly the hard part.

    But what if the correction is LOW-RANK?

    S(x, A ∪ B) = |{n ≤ x : gcd(n, prod(A)*prod(B)) = 1}|
                 = Σ_{d | prod(A)} Σ_{e | prod(B)} μ(d)μ(e) floor(x/(de))

    = Σ_{d | prod(A)} μ(d) * Σ_{e | prod(B)} μ(e) floor(x/(de))
    = Σ_{d | prod(A)} μ(d) * S(x/d, B)

    This IS the Meissel-Lehmer recursion! The key question: can S(x/d, B)
    be represented compactly for all d simultaneously?

    If S(·, B) can be represented by a polynomial (or low-degree function),
    then we can evaluate Σ μ(d) S(x/d, B) efficiently.
    """
    print("\n" + "=" * 70)
    print("RECURSIVE SIEVE IDEA: Low-rank factorization")
    print("=" * 70)

    # Test: for fixed set B of primes, how well can S(x, B) be approximated
    # by a polynomial in x?

    from itertools import combinations

    for B in [[2,3], [2,3,5], [2,3,5,7], [2,3,5,7,11], [2,3,5,7,11,13]]:
        primorial_B = 1
        for p in B:
            primorial_B *= p
        euler_B = np.prod([1 - 1/p for p in B])

        # Compute S(x, B) for x = 1 to 1000
        x_vals = np.arange(1, 1001)
        S_vals = np.zeros(len(x_vals))

        divisors_B = []
        for k in range(len(B) + 1):
            for S in combinations(B, k):
                d = 1
                for p in S:
                    d *= p
                divisors_B.append((d, (-1)**k))

        for i, x in enumerate(x_vals):
            S_vals[i] = sum(mu * (x // d) for d, mu in divisors_B)

        # S(x, B) = x * euler_B + oscillatory_correction
        correction = S_vals - x_vals * euler_B

        # Is the correction periodic with period primorial(B)?
        period = primorial_B
        if period <= 500:
            periodic_residuals = []
            for i in range(period, min(len(x_vals), 2*period)):
                periodic_residuals.append(abs(correction[i] - correction[i - period]))
            max_periodic_error = max(periodic_residuals) if periodic_residuals else float('inf')
        else:
            max_periodic_error = float('inf')

        print(f"\nB = {B}")
        print(f"  primorial = {primorial_B}, euler = {euler_B:.6f}")
        print(f"  Correction range: [{np.min(correction):.2f}, {np.max(correction):.2f}]")
        print(f"  Correction is periodic with period {period}: max_error = {max_periodic_error}")

        # Since correction is periodic, S(x, B) = x * euler_B + f(x mod primorial_B)
        # where f is a known function on {0, 1, ..., primorial_B - 1}

        # This means: to compute S(x, B), we only need x mod primorial_B!
        # And the Euler product gives the linear term.

        # KEY INSIGHT: S(x, A∪B) = Σ_{d|prod(A)} μ(d) S(x/d, B)
        #            = Σ_{d|prod(A)} μ(d) [x/d * euler_B + f((x/d) mod primorial_B)]
        #            = euler_B * S(x, A) + Σ_{d|prod(A)} μ(d) f(floor(x/d) mod primorial_B)

        # The second sum involves floor(x/d) mod primorial_B for various d.
        # floor(x/d) mod m can be computed in O(1) per d.
        # So the total cost is O(2^|A|) for the sum over d.

        # If we split primes into groups of size log(log(x)), each group A has
        # 2^|A| = log(x) terms. And we have pi(y)/log(log(x)) groups.
        # Total: pi(y) * log(x) / log(log(x)) ... still not polylog for y = x^{1/3}.

        print(f"  → S(x,B) = {euler_B:.4f}*x + f(x mod {primorial_B}) [EXACTLY PERIODIC]")


if __name__ == "__main__":
    print("SIEVE MATRIX RANK EXPERIMENT")
    print("Testing whether the sieve has exploitable low-rank structure\n")

    # Experiment 1: Basic rank analysis
    analyze_rank_structure(
        x_values=[50, 100, 200, 500, 1000],
        y_func=lambda x: int(x**0.5) + 1
    )

    # Experiment 2: Transfer matrix / Mobius sum analysis
    transfer_matrix_experiment(1000, 7)  # Sieve by {2,3,5,7}

    # Experiment 3: Fractional part analysis (the key to exactness)
    fractional_part_analysis(list(range(1, 501)), 7)

    # Experiment 4: Rank growth scaling
    rank_growth_scaling()

    # Experiment 5: Nullspace structure
    nullspace_structure()

    # Experiment 6: Recursive sieve idea
    recursive_sieve_idea()

    print("\n" + "=" * 70)
    print("SUMMARY OF KEY FINDINGS")
    print("=" * 70)
    print("""
    1. Sieve matrix rank: Is it much less than min(x, pi(y))?
    2. Fractional correction: Is it smooth/periodic/predictable?
    3. Recursive factorization: Can the combination step be made cheap?

    The critical question: Is there a way to avoid the 2^pi(y) blowup
    in the Mobius inclusion-exclusion that doesn't reduce to Meissel-Lehmer?
    """)

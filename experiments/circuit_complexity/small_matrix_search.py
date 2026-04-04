"""
Small Matrix Search for det(M) = pi(x)
Session 14 - April 2026

Goal: For small x, search for integer matrices M of size k x k (k << x)
whose determinant equals pi(x). If such matrices exist with k = O(N) or
O(N^2) where N = log2(x), this would prove pi(x) in GapL.

Approaches:
1. Random integer matrix search (baseline)
2. Structured search: companion matrices, triangular matrices
3. Optimization: gradient-free search over integer matrices
4. Constructive: try to build M from number-theoretic data

For x=10..100, pi(x) ranges from 4 to 25, so we need det(M) to hit
specific small integers. This is easy for 1x1 (trivially M = [[pi(x)]]),
but the question is whether there's a SYSTEMATIC construction that
scales, not just ad-hoc solutions.
"""

import numpy as np
from itertools import product
import math
import sympy
from collections import defaultdict
import time

def pi(x):
    """Prime counting function."""
    return int(sympy.primepi(x))

def search_2x2_det(target):
    """
    Find all 2x2 integer matrices with entries in [-B, B] whose det = target.
    det([[a,b],[c,d]]) = ad - bc = target
    """
    solutions = []
    B = 10
    for a in range(-B, B+1):
        for d in range(-B, B+1):
            # ad - bc = target => bc = ad - target
            bc = a*d - target
            # Find all factorizations of bc
            if bc == 0:
                # b=0 or c=0
                for b in range(-B, B+1):
                    if b == 0:
                        for c in range(-B, B+1):
                            solutions.append(((a,0),(c,d)))
                    else:
                        solutions.append(((a,b),(0,d)))
            else:
                for b in range(-B, B+1):
                    if b != 0 and bc % b == 0:
                        c = bc // b
                        if -B <= c <= B:
                            solutions.append(((a,b),(c,d)))
    return solutions

def search_systematic_construction(x_max=100):
    """
    Try to find a SYSTEMATIC matrix construction M(x) such that det(M(x)) = pi(x).

    Idea 1: Use the Redheffer-like construction but compressed.
    The Redheffer matrix R[i][j] = 1 if j|i or j=1, and det(R) = M(x) (Mertens).
    For pi(x), we need something different.

    Idea 2: Use inclusion-exclusion. pi(x) = (x-1) - sum_{p<=sqrt(x)} |multiples of p in [2,x]|
    + sum_{p<q} |multiples of pq in [2,x]| - ...
    This is an alternating sum that CAN be expressed as a determinant (via Cauchy-Binet or LGV).

    Idea 3: Build a matrix from the Legendre sieve.
    """
    print("=== Systematic Matrix Construction Search ===\n")

    # Idea 2: Inclusion-exclusion as determinant
    # The Legendre formula: pi(x) - pi(sqrt(x)) + 1 = sum_{S subset primes<=sqrt(x)} (-1)^|S| floor(x / prod(S))
    # This is exactly the permanent of a certain matrix... or can be a determinant?

    # Actually, inclusion-exclusion sum = det of a matrix by the following trick:
    # For primes p1, ..., pk (primes <= sqrt(x)):
    # sum_{S} (-1)^|S| f(S) = det(I - A) where A encodes the contributions

    # Let's try: define a matrix M of size (k+1) x (k+1) where k = pi(sqrt(x))
    # Such that det(M) relates to pi(x)

    for x in [10, 20, 30, 50, 100]:
        N = math.ceil(math.log2(x)) if x > 1 else 1
        target = pi(x)
        sqrtx = int(math.isqrt(x))
        primes_below_sqrt = [p for p in range(2, sqrtx+1) if sympy.isprime(p)]
        k = len(primes_below_sqrt)

        print(f"x={x}, N={N}, pi(x)={target}, sqrt(x)={sqrtx}, "
              f"primes<={sqrtx}: {primes_below_sqrt} (k={k})")

        # Legendre formula:
        # pi(x) = pi(sqrt(x)) - 1 + sum_{S subset primes<=sqrt(x)} (-1)^|S| floor(x / prod(S))
        # where S ranges over all subsets

        legendre_sum = 0
        for mask in range(1 << k):
            prod_S = 1
            bits = 0
            for i in range(k):
                if mask & (1 << i):
                    prod_S *= primes_below_sqrt[i]
                    bits += 1
            legendre_sum += ((-1)**bits) * (x // prod_S)

        # legendre_sum should equal pi(x) - pi(sqrt(x)) + 1
        # because it counts integers in [1,x] with no prime factor <= sqrt(x)
        # (i.e., 1 and primes in (sqrt(x), x])

        expected = target - pi(sqrtx) + 1  # 1 + primes in (sqrt(x), x]
        print(f"  Legendre I-E sum = {legendre_sum}, expected = {expected}, match = {legendre_sum == expected}")

        # Now: can this inclusion-exclusion sum be expressed as a SMALL determinant?
        # The I-E sum over 2^k subsets...
        # Key insight: the matrix-tree theorem and LGV lemma express certain
        # alternating sums as determinants.

        # Attempt: construct a (k+1) x (k+1) matrix where:
        # Row 0 represents the "total" count x
        # Row i represents "sieving by prime p_i"
        # And the determinant gives the I-E sum

        # A well-known trick: det(I + D) where D is diagonal with entries d_i
        # = prod(1 + d_i) = sum_{S} prod_{i in S} d_i
        # But we need ALTERNATING signs and floor(x/prod(S)), not prod of independent terms.

        # Another approach:
        # det of [[1, a1, a2, ..., ak],
        #         [1, 1+b1, 0, ..., 0],
        #         [1, 0, 1+b2, ..., 0],
        #         ...
        #         [1, 0, 0, ..., 1+bk]]
        # This doesn't directly give I-E.

        # Direct construction attempt:
        # We want det(M) = sum_{S subset [k]} (-1)^|S| floor(x / prod_{i in S} p_i)
        # Expand det by Leibniz formula: det(M) = sum_{sigma in S_n} sgn(sigma) prod M[i][sigma(i)]

        # For k=1 (one prime p): sum = floor(x/1) - floor(x/p) = x - floor(x/p)
        # Need 2x2 matrix with det = x - floor(x/p)
        # [[x, floor(x/p)], [1, 1]] has det = x - floor(x/p). YES!
        # Or [[1, 1], [floor(x/p), x]] has det = x - floor(x/p). Also works.

        # For k=2 (primes p, q):
        # sum = x - floor(x/p) - floor(x/q) + floor(x/(pq))
        # Need 3x3 matrix with this determinant?

        # Try upper triangular + correction:
        # For the I-E formula with k primes, this is the PERMANENT of a 0-1 matrix
        # evaluated at specific points... actually no, it's simpler.

        # The I-E sum is: sum_{S} (-1)^|S| f(prod(S))
        # where f(d) = floor(x/d)

        # This equals prod_{i=1}^{k} (1 - T_{p_i}) applied to f(1) = x
        # where T_p is the operator f(d) -> f(dp)

        # In multiplicative number theory, this is the Legendre sieve.
        # The question: can prod(I - T_{p_i}) be represented as a small matrix?

        print(f"  Sieve as operator product: need product of {k} operators (I - T_p)")
        print()

def search_random_matrices(target, size, n_trials=100000, entry_bound=5):
    """Search for size x size integer matrices with det = target."""
    count = 0
    for _ in range(n_trials):
        M = np.random.randint(-entry_bound, entry_bound+1, size=(size, size))
        d = int(round(np.linalg.det(M)))
        if d == target:
            count += 1
            if count <= 3:
                return M, count
    return None, count

def explore_triangular_constructions():
    """
    For triangular matrices, det = product of diagonal.
    So det(M) = pi(x) iff product of diagonal entries = pi(x).
    This is trivial (diagonal = [pi(x), 1, 1, ..., 1]).

    The INTERESTING question is: can the entries be logspace-computable
    functions of x? For a 1x1 matrix, M = [[pi(x)]] requires computing
    pi(x) to fill the entry, which is circular.

    For larger matrices, we need entries that are SIMPLE functions of x
    but whose determinant encodes the complex computation of pi(x).
    """
    print("=== Triangular / Structured Matrix Constructions ===\n")

    # The Redheffer matrix for Mertens: M(x) = det(R_x)
    # R_x[i][j] = 1 if j|i or j=1, size x x x
    # This is NOT logspace-sized, but the entries ARE simple (divisibility check)

    # Can we build a SMALL analog for pi(x)?

    # Idea: Define a matrix indexed by (d, d') where d, d' are in a small set S
    # S should have size poly(log x)

    # The floor-value set {floor(x/k)} has size O(sqrt(x)) - too large
    # But what if we only use a subset?

    # The Lucy DP uses floor values and prime sieving steps.
    # If we could "compress" the floor-value set to poly(log x) values...

    # Key observation: pi(x) depends on ALL floor values (proven in Session 12,
    # the linear transformation is full-rank on ~80% of values).
    # So we CAN'T just drop most floor values.

    # But maybe a DIFFERENT basis? Instead of floor values, use some other set of
    # O(poly(N)) intermediate quantities?

    print("Observation: any poly(N)-size matrix construction for pi(x) must")
    print("avoid computing individual floor values (there are O(sqrt(x)) of them).")
    print("It must use a FUNDAMENTALLY DIFFERENT set of intermediate values.\n")

    # What if the matrix entries involve ANALYTIC quantities?
    # E.g., values of zeta(s) at specific points?
    # det of a matrix of zeta values?

    # Selberg trace formula: relates spectral data to geometric data via a determinant
    # Fredholm determinant of an operator = product over eigenvalues
    # det(I - K) where K is a trace-class operator on L^2

    # For the zeta function: det(I - p^{-s}) over primes = 1/zeta(s)
    # So zeta(s) = 1 / prod_p (1 - p^{-s})
    # And pi(x) = sum_{p <= x} 1

    # But computing zeta(s) at enough precision requires O(sqrt(x)) terms anyway.

    print("Analytic approach also requires O(sqrt(x)) zeta zeros or Dirichlet terms.")
    print("No obvious poly(N)-size matrix construction.\n")

    # Explicit small matrices for small x
    print("--- Explicit small matrix constructions for small x ---")
    for x in [10, 20, 30, 50, 100]:
        target = pi(x)
        N = math.ceil(math.log2(x))
        print(f"\nx={x}, N={N}, pi(x)={target}")

        # 2x2: det = ad - bc = target. Many solutions.
        # Interesting: can we find a,b,c,d as SIMPLE functions of x?

        # E.g., [[x, x-target], [1, 1]] has det = x - (x-target) = target. Trivial.
        # This is circular because we need to know target=pi(x) to build the matrix.

        # What about [[x, floor(x/2)], [2, 1]]?
        # det = x*1 - floor(x/2)*2 = x - 2*floor(x/2) = x mod 2
        # That gives 0 or 1, not pi(x).

        # What about a matrix built from Euler products?
        # For 2x2: can we express pi(x) as a determinant of quantities involving
        # partial Euler products?

        # Partial Euler product: P(s, x) = prod_{p <= x} (1 - p^{-s})
        # log P(s, x) = sum_{p <= x} log(1 - p^{-s})
        # d/ds log P(s, x) = sum_{p <= x} p^{-s} log(p) / (1 - p^{-s})
        # At s -> infinity: P -> 1, d/ds P -> 0
        # At s = 0: P(0, x) = prod_{p<=x} (1-1) = 0 if x >= 2

        # Not directly useful for extracting pi(x).

        # Counting interpretation: pi(x) = #{n in [2,x] : n is prime}
        # = #{n in [2,x]} - #{n in [2,x] : n is composite}
        # = (x-1) - #{composites in [2,x]}
        # = (x-1) - sum_{n=2}^{x} (1 - chi_prime(n))

        # This is just rewriting, not helpful.

        print(f"  Trivial 1x1: M = [[{target}]], but this is circular.")
        print(f"  Trivial 2x2: M = [[{x}, {x-target}], [1, 1]], also circular.")

def attempt_compressed_redheffer(x):
    """
    The Redheffer matrix R_x is x x x with R[i][j] = 1 if j|i or j=1.
    det(R_x) = M(x) = Mertens function.

    To get pi(x), note: pi(x) = sum_{n=1}^{x} |mu(n)| * chi_prime(n)
    Alternatively: consider modifying R to count primes instead of Mertens.

    Key idea: analyze the RANK and structure of R_x.
    If rank(R_x) << x, we might compress it.
    """
    print(f"\n=== Compressed Redheffer Analysis for x = {x} ===")

    # Build Redheffer matrix
    R = np.zeros((x, x), dtype=float)
    for i in range(1, x+1):
        for j in range(1, x+1):
            if j == 1 or (i % j == 0):
                R[i-1][j-1] = 1

    # Compute determinant (should be M(x))
    det_R = int(round(np.linalg.det(R)))
    M_x = sum(sympy.mobius(n) for n in range(1, x+1))
    print(f"  det(R_{x}) = {det_R}, M({x}) = {M_x}, match = {det_R == M_x}")

    # Analyze rank and singular values
    U, sigma, Vt = np.linalg.svd(R)
    rank = np.sum(sigma > 1e-10)

    print(f"  Size: {x} x {x}")
    print(f"  Rank: {rank}")
    print(f"  Top 10 singular values: {sigma[:10].round(3)}")
    print(f"  Bottom 10 singular values: {sigma[-10:].round(6)}")

    # How many singular values are needed to capture most of the "information"?
    total_energy = np.sum(sigma**2)
    cumulative = np.cumsum(sigma**2) / total_energy
    for threshold in [0.90, 0.95, 0.99, 0.999]:
        k = np.searchsorted(cumulative, threshold) + 1
        print(f"  SVs for {threshold:.1%} energy: {k} (of {x})")

    # The key: rank is always x (full rank) because det != 0 (when M(x) != 0)
    # But maybe there's low-rank structure AFTER removing the trivial parts?

    # Decompose R = J + D where J has all 1s in column 1, D has divisibility
    J = np.zeros((x, x))
    J[:, 0] = 1
    D = R - J  # D[i][j] = 1 if j|i and j > 1, else 0

    U_D, sigma_D, _ = np.linalg.svd(D)
    rank_D = np.sum(sigma_D > 1e-10)
    print(f"\n  D = R - J (divisibility without column 1):")
    print(f"  Rank of D: {rank_D}")
    print(f"  Top 10 SVs of D: {sigma_D[:10].round(3)}")

    return det_R, rank, sigma

def build_prime_counting_matrix(x):
    """
    Try to build a matrix whose determinant equals pi(x) directly.

    Approach: modify the Redheffer construction.

    The Redheffer det = M(x) = sum mu(n).
    We want pi(x) = sum chi_prime(n).

    Note: chi_prime(n) = mu(n)^2 * chi_squarefree_prime(n)
    Actually chi_prime(n) = 1 iff n has exactly one prime factor (itself).

    Alternative: pi(x) = sum_{n=2}^{x} sum_{d|n, d<n} mu(d)  ... no, that's different.

    Actually, there IS a matrix for pi(x):
    Define Q_x as x x x matrix with Q[i][j] = mu(j) if j|i, else 0.
    Then (Q * 1)[i] = sum_{j|i} mu(j) = [i == 1] (Mobius inversion)
    So Q * ones_vector = e_1.

    For pi(x), we need sum_{n=2}^{x} [n is prime]
    = sum_{n=2}^{x} [Omega(n) == 1 and n is squarefree]

    Hmm, this doesn't simplify to a nice matrix.

    Let's try a different approach: the Smith determinant.
    det of (gcd(i,j)) matrix = prod_{k=1}^{n} phi(k)
    det of (lcm(i,j)) matrix doesn't have a nice form

    What about: det of matrix A where A[i][j] = [j | i] * w(j)?
    det(A) = prod_{k=1}^{n} w(k) (since A is lower triangular with diagonal w(k))

    So if we want det = pi(x), we need prod w(k) = pi(x).
    One way: w(k) = 1 for all k except w(1) = pi(x). Circular again.
    """
    print(f"\n=== Prime-Counting Matrix Construction for x = {x} ===")

    target = pi(x)
    print(f"  Target: pi({x}) = {target}")

    # Smith-type: det(gcd(i,j))_{i,j=1}^{n} = prod phi(k)
    n = x
    GCD = np.array([[math.gcd(i+1, j+1) for j in range(n)] for i in range(n)], dtype=float)
    det_gcd = np.linalg.det(GCD)
    prod_phi = 1
    for k in range(1, n+1):
        prod_phi *= sympy.totient(k)
    print(f"  det(GCD matrix) = {det_gcd:.0f}, prod phi(k) = {prod_phi}")

    # LCM matrix
    LCM = np.array([[math.lcm(i+1, j+1) for j in range(n)] for i in range(n)], dtype=float)
    det_lcm = np.linalg.det(LCM)
    print(f"  det(LCM matrix) = {det_lcm:.6e}")

    # Divisor matrix: D[i][j] = 1 if j|i
    DIV = np.array([[1 if (i+1) % (j+1) == 0 else 0 for j in range(n)] for i in range(n)], dtype=float)
    det_div = np.linalg.det(DIV)
    print(f"  det(Divisor matrix) = {det_div:.0f} (should be 1, it's lower triangular with 1s on diag)")

    # Weighted divisor: D[i][j] = mu(j) if j|(i+1), but only if relevant
    # What if we weight rows/columns to make det = pi(x)?

    # Another approach: characteristic polynomial
    # For any matrix A, det(A - lambda*I) = char poly, and the constant term = det(A)
    # But this doesn't help construct A.

    # Direct approach: build a block matrix
    # [[I_k, B], [C, D]] where k is small and det = det(D - C*B) if I_k is identity
    # = det(D) - trace adjustments...
    # This is Schur complement: det = det(I_k) * det(D - C * I_k^{-1} * B) = det(D - CB)

    return target

def experiment_lucy_as_matrix_product(x):
    """
    Express the Lucy DP as a product of affine transformations,
    then analyze the resulting matrix.

    S_{step+1} = A_p * S_{step} + b_p
    where A_p = I - E_p (E_p is the sieve-step matrix)
    and b_p adds back pi(p-1) for affected entries.

    The full computation: S_final = A_last * ... * A_2 * S_0 + accumulated_affine

    The product M = A_last * ... * A_2 is a |V| x |V| matrix.
    Question: does M have special structure? Low rank? Sparse?
    """
    print(f"\n=== Lucy DP as Matrix Product for x = {x} ===")

    V = sorted(set(x // k for k in range(1, x+1)) | set(range(1, int(math.isqrt(x))+1)))
    V = sorted(set(V))
    n = len(V)
    v_to_idx = {v: i for i, v in enumerate(V)}

    sqrtx = int(math.isqrt(x))
    primes = [p for p in range(2, sqrtx+1) if sympy.isprime(p)]

    print(f"  |V| = {n}, primes up to sqrt({x}) = {primes}")

    # Build each A_p matrix
    # A_p = I, except for rows corresponding to v >= p^2:
    # A_p[v_idx][v_idx] = 1 (diagonal)
    # A_p[v_idx][floor(v/p)_idx] -= 1
    # A_p[v_idx][(p-1)_idx] += 1  (wait, this is the affine term from S[p-1])

    # Actually, let's be precise about the affine decomposition.
    # S_new[v] = S_old[v] - (S_old[floor(v/p)] - S_old[p-1])  for v >= p^2
    # S_new[v] = S_old[v] - S_old[floor(v/p)] + S_old[p-1]
    #
    # BUT S_old[p-1] at step p equals pi(p-1), which is a CONSTANT at that step.
    # So the linear part is: S_new[v] = S_old[v] - S_old[floor(v/p)]
    # And the affine part is: + pi(p-1) for each v >= p^2
    #
    # Wait, let me reconsider. If we treat S as a vector over V, then:
    # For v >= p^2: S_new[v] = S_old[v] - S_old[floor(v/p)] + S_old[p-1]
    # This is LINEAR in S_old (since p-1 is in V and S_old[p-1] is a component).
    # So A_p is truly linear, not just affine!
    #
    # A_p[v, v] = 1  (always)
    # A_p[v, floor(v/p)] -= 1  (if v >= p^2)
    # A_p[v, p-1] += 1  (if v >= p^2)

    # Build the product M = A_last * ... * A_2
    M = np.eye(n)

    for p in primes:
        A_p = np.eye(n)
        for v in V:
            if v >= p * p:
                vi = v_to_idx[v]
                vpi = v_to_idx[v // p]
                pm1i = v_to_idx[p - 1]
                A_p[vi][vpi] -= 1
                A_p[vi][pm1i] += 1
        M = A_p @ M

    # Now S_final = M @ S_initial
    # S_initial[v] = v - 1 for all v
    S_initial = np.array([v - 1 for v in V], dtype=float)
    S_final = M @ S_initial

    # pi(x) = S_final[x]
    xi = v_to_idx[x]
    computed_pi = int(round(S_final[xi]))
    actual_pi = pi(x)
    print(f"  pi({x}) via matrix product: {computed_pi}, actual: {actual_pi}, match: {computed_pi == actual_pi}")

    # Analyze M
    print(f"  M is {n}x{n}")
    rank_M = np.linalg.matrix_rank(M)
    print(f"  Rank of M: {rank_M}")

    U, sigma, Vt = np.linalg.svd(M)
    print(f"  Top 10 singular values: {sigma[:10].round(4)}")
    if n > 10:
        print(f"  Bottom 10 singular values: {sigma[-10:].round(6)}")

    det_M = np.linalg.det(M)
    print(f"  det(M) = {det_M:.6f}")

    # Sparsity
    nonzero = np.sum(np.abs(M) > 1e-10)
    print(f"  Nonzero entries: {nonzero} / {n*n} ({nonzero/(n*n):.1%})")

    # Key question: does M have low-rank structure?
    total_energy = np.sum(sigma**2)
    cumulative = np.cumsum(sigma**2) / total_energy
    for threshold in [0.90, 0.95, 0.99, 0.999]:
        k = np.searchsorted(cumulative, threshold) + 1
        print(f"  SVs for {threshold:.1%} energy: {k} / {n}")

    # The CRUCIAL row: the one for x
    # pi(x) = M[xi] . S_initial
    # So pi(x) = sum_v M[xi][v_idx] * (v - 1)
    row_x = M[xi]
    nonzero_in_row = np.sum(np.abs(row_x) > 1e-10)
    print(f"\n  Row for x={x}: {nonzero_in_row} nonzero entries out of {n}")
    print(f"  Row entries: {row_x.round(4)}")

    # This row gives the EXACT linear combination of initial S-values
    # that yields pi(x). If most entries are zero, we might not need all floor values.

    # Actually, the individual A_p matrices are sparse, but their product M
    # may be dense. Let's check how the density grows.

    return M, V

def experiment_legendre_ie_as_det():
    """
    Inclusion-exclusion formula for the Legendre sieve.

    pi(x) - pi(sqrt(x)) + 1 = sum_{S subset P} (-1)^|S| floor(x / prod(S))
    where P = {primes <= sqrt(x)}.

    This is 2^|P| terms. Can we express it as a determinant of an |P| x |P| matrix?

    Key insight: the matrix permanent is sum_{sigma} prod a_{i,sigma(i)}
    and the determinant is sum_{sigma} sgn(sigma) prod a_{i,sigma(i)}.

    Our sum is over SUBSETS, not permutations. But there's a classic trick:

    sum_{S subset [k]} (-1)^|S| f(S) = det(M) where M is constructed from f.

    Specifically, if f(S) = g(prod_{i in S} p_i) for a multiplicative input,
    and g(d) = floor(x/d), then:

    Consider the matrix M[i][j] for i,j in {0,1,...,k}:
    Actually, this is related to the "Bonferroni" / "sieve" matrix approach.
    """
    print("\n=== Inclusion-Exclusion as Determinant ===\n")

    for x in [30, 50, 100, 200]:
        target = pi(x)
        sqrtx = int(math.isqrt(x))
        P = [p for p in range(2, sqrtx+1) if sympy.isprime(p)]
        k = len(P)

        # Compute the I-E sum directly
        ie_sum = 0
        for mask in range(1 << k):
            prod_S = 1
            bits = bin(mask).count('1')
            for i in range(k):
                if mask & (1 << i):
                    prod_S *= P[i]
            ie_sum += ((-1)**bits) * (x // prod_S)

        # ie_sum = 1 + #{primes in (sqrt(x), x]} = 1 + pi(x) - pi(sqrt(x))
        pi_from_ie = ie_sum - 1 + pi(sqrtx)

        print(f"x={x}, primes<={sqrtx}: {P}, k={k}")
        print(f"  I-E sum = {ie_sum}, pi(x) = {target}, reconstructed = {pi_from_ie}, match = {pi_from_ie == target}")

        # Now try to express this as a determinant
        # Method: the multilinear extension
        # Define f(t_1, ..., t_k) = sum_{S} prod_{i in S} t_i * prod_{i not in S} (1-t_i) * floor(x/prod_{i in S} p_i)
        # At t_i = -1: f(-1,...,-1) = sum_S (-1)^|S| * (-2)^{k-|S|} * floor(x/prod_S)
        # Hmm, not quite right.

        # Simpler: define the polynomial
        # Q(t) = prod_{i=1}^{k} (1 - t/p_i)  evaluated at t=... something
        # prod(1 - t/p_i) = sum_{S} (-1)^|S| t^|S| / prod_{i in S} p_i
        # Not directly floor(x/prod_S).

        # Direct matrix attempt for small k:
        if k <= 4:
            # Build all possible (k+1)x(k+1) matrices with entries from {floor(x/d) : d | prod(P)} + {0, 1, -1}
            # and check if det = ie_sum
            # This is too large for brute force. Try structured approaches.

            # Approach: LGV lemma
            # Build a DAG with k+1 source nodes and k+1 sink nodes
            # Path weights are floor(x/d) for various d
            # Non-intersecting path systems give determinants

            # For k=2, P = {p, q}:
            # ie_sum = x - floor(x/p) - floor(x/q) + floor(x/(pq))
            # = (x - floor(x/p)) - (floor(x/q) - floor(x/(pq)))
            # = (x - floor(x/q)) - (floor(x/p) - floor(x/(pq)))

            # As 2x2 det:
            # det([[x, floor(x/q)], [floor(x/p)/floor(x/(pq)), 1]])
            # = x*1 - floor(x/q) * floor(x/p)/floor(x/(pq))
            # Nope, entries must be integers.

            # Try: det([[a, b], [c, d]]) = ad - bc = ie_sum
            # With a, b, c, d from {x, floor(x/p), floor(x/q), floor(x/(pq))}

            if k >= 2:
                p, q = P[0], P[1]
                vals = {
                    'x': x,
                    'x/p': x // p,
                    'x/q': x // q,
                    'x/pq': x // (p*q),
                    '1': 1,
                    '0': 0,
                }

                print(f"  Floor values: {vals}")
                print(f"  ie_sum = {ie_sum}")

                # Search 2x2 determinants
                found_2x2 = False
                for (n1, v1), (n2, v2), (n3, v3), (n4, v4) in product(vals.items(), repeat=4):
                    if v1*v4 - v2*v3 == ie_sum:
                        if not found_2x2:
                            print(f"  Found 2x2: det([[{n1},{n2}],[{n3},{n4}]]) = {v1}*{v4} - {v2}*{v3} = {ie_sum}")
                            found_2x2 = True
                            break

                if not found_2x2:
                    print(f"  No 2x2 determinant found with these entries.")

                    # Try with additional entries
                    extra = {'x-1': x-1, 'p': p, 'q': q, 'p-1': p-1, 'q-1': q-1}
                    vals.update(extra)
                    for (n1, v1), (n2, v2), (n3, v3), (n4, v4) in product(vals.items(), repeat=4):
                        if v1*v4 - v2*v3 == ie_sum:
                            print(f"  Found 2x2 (extended): det([[{n1},{n2}],[{n3},{n4}]]) = {v1}*{v4} - {v2}*{v3} = {ie_sum}")
                            break

def main():
    print("=" * 70)
    print("SMALL MATRIX SEARCH: det(M) = pi(x)")
    print("Session 14 - Determinant approach to prime counting")
    print("=" * 70)

    # 1. Systematic construction attempts
    search_systematic_construction()

    # 2. Triangular / structured constructions
    explore_triangular_constructions()

    # 3. Compressed Redheffer for small x
    for x in [10, 20, 30, 50]:
        attempt_compressed_redheffer(x)

    # 4. Lucy DP as matrix product
    for x in [30, 50, 100]:
        experiment_lucy_as_matrix_product(x)

    # 5. Inclusion-exclusion as determinant
    experiment_legendre_ie_as_det()

    # 6. Build prime-counting matrix
    for x in [10, 15, 20]:
        build_prime_counting_matrix(x)

    print("\n" + "=" * 70)
    print("CONCLUSIONS")
    print("=" * 70)
    print("""
Summary of findings:

1. TRIVIAL constructions (1x1 or 2x2 with pi(x) as entry) are CIRCULAR.
   They require knowing pi(x) to build the matrix.

2. The Redheffer matrix gives M(x) = Mertens, but is x*x (exponential).
   Its rank is FULL (= x when M(x) != 0). SVD shows no obvious low-rank
   structure that could compress it to poly(log x) size.

3. The Lucy DP matrix product M = A_last * ... * A_2 is |V| x |V| where
   |V| = O(sqrt(x)). The individual A_p are sparse, but the product M
   becomes DENSE. The row for x has O(sqrt(x)) nonzero entries = ALL
   floor values contribute. This confirms Session 12 finding.

4. The inclusion-exclusion (Legendre sieve) has 2^k terms where k = pi(sqrt(x)).
   For small k, 2x2 determinants CAN sometimes express the I-E sum using
   floor values, but this doesn't scale (k grows as sqrt(x)/ln(x)).

5. No poly(N)-size matrix construction found. The fundamental obstacle:
   pi(x) depends on O(sqrt(x)) independent pieces of information
   (floor values or zeta zeros), and no known compression exists.

KEY INSIGHT: A GapL algorithm would need to avoid BOTH floor-value sets
AND zeta zero sums. It would need a fundamentally new set of O(poly(N))
intermediate quantities whose determinant yields pi(x).
""")

if __name__ == "__main__":
    main()

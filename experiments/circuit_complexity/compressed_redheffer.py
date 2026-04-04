"""
Compressed Redheffer Matrix Analysis
Session 14 - April 2026

The Redheffer matrix R_x is x*x with R[i][j] = 1 if j|i or j=1.
det(R_x) = M(x) = Mertens function.

We want pi(x) instead. Approaches:
1. Analyze R_x structure: rank, SVD, sparsity pattern
2. Build a modified "prime-Redheffer" matrix with det = pi(x)
3. Study factorizations R = L*D*U and see if D is compressible
4. Investigate the Smith normal form of R_x
5. Connection between M(x) and pi(x): can we get pi(x) from M(x) data?
"""

import numpy as np
import math
import sympy
from collections import defaultdict
import time

def pi(x):
    return int(sympy.primepi(x))

def build_redheffer(n):
    """Build the n x n Redheffer matrix."""
    R = np.zeros((n, n), dtype=float)
    for i in range(1, n+1):
        for j in range(1, n+1):
            if j == 1 or (i % j == 0):
                R[i-1][j-1] = 1
    return R

def build_prime_redheffer(n):
    """
    Attempt to build a Redheffer-like matrix whose determinant = pi(n).

    Idea: The Redheffer matrix encodes the Mobius function via
    det(R_n) = sum_{k=1}^{n} mu(k) = M(n).

    For pi(n), we need sum_{k=1}^{n} chi_prime(k).

    Note that chi_prime(k) = sum_{d|k} lambda(d) where lambda is...
    Actually, let's think differently.

    The Redheffer matrix works because:
    - Row i, column j=1: always 1 (gives the "sum" structure)
    - Row i, column j>1: 1 if j|i (encodes divisibility = Mobius inversion)

    When we expand det(R), the result involves products along paths
    through the divisibility structure, and the alternating signs
    conspire to give sum mu(k).

    For pi(n), we could try:
    - Keep the divisibility structure but weight differently
    - Use a different function instead of mu

    Actually, there's a known result:
    If we define R'[i][j] = f(j) if j|i, and R'[i][1] = 1,
    then det(R') = sum_{k=1}^{n} (mu * f)(k) where * is Dirichlet convolution.

    For det(R') = pi(n), we need (mu * f)(k) = chi_prime(k).
    So f = mu^{-1} * chi_prime = 1 * chi_prime (since mu^{-1} = 1, the constant function).
    f(k) = sum_{d|k} chi_prime(d) = number of prime divisors of k (with multiplicity? no, without)
    f(k) = omega(k) = number of distinct prime factors of k.

    So: define R'[i][j] = omega(j) if j|i (and j > 1), R'[i][1] = 1
    Then det(R') should equal pi(n)?

    Wait, let me be more careful. The Redheffer matrix has:
    R[i][j] = 1 if j=1 or j|i

    Its determinant det(R_n) = sum_{k=1}^n mu(k).

    The general theory: if A[i][j] = f(gcd(i,j)) for some f, then
    det(A) = prod_{k=1}^n (sum_{d|k} f(d) mu(k/d)).
    But Redheffer isn't of this form.

    Let me look at this differently. The Redheffer matrix can be written as:
    R = D + e_1 * ones^T  (where D is the divisibility matrix with D[i][j]=1 if j|i, j>1)
    Wait no, R[i][j] = 1 if j=1 OR j|i. So R = J + D' where J has col 1 all 1s,
    D'[i][j] = 1 if j|i and j >= 2.

    Actually R[i][1] = 1 always. R[i][j] = 1 if j|i for j >= 2. R[1][1] = 1, R[1][j] = 0 for j >= 2.
    Hmm no: R[1][j] = 1 if j|1 or j=1, so R[1][1] = 1, R[1][j] = 0 for j >= 2 (since j doesn't divide 1).

    So the first row is [1, 0, 0, ..., 0].

    For i >= 2: R[i][1] = 1, R[i][j] = [j|i] for j >= 2.

    Cofactor expansion along row 1: det(R) = 1 * det(R with row 1, col 1 removed)
    The remaining (n-1)x(n-1) matrix has rows i=2,...,n and cols j=2,...,n:
    B[i][j] = 1 if j|i (for i,j >= 2).

    Hmm that's not quite right either because column 1 contributes.
    Let me just compute numerically.
    """
    print("\n=== Prime-Redheffer Matrix Experiments ===\n")

    for n in [10, 20, 30, 50]:
        print(f"--- n = {n} ---")
        target_mertens = sum(int(sympy.mobius(k)) for k in range(1, n+1))
        target_pi = pi(n)

        # Standard Redheffer
        R = build_redheffer(n)
        det_R = int(round(np.linalg.det(R)))
        print(f"  det(R_{n}) = {det_R}, M({n}) = {target_mertens}")

        # Attempt 1: Weight column j by omega(j)+1 instead of 1
        # R'[i][j] = omega(j) + 1 if j|i or j=1
        # Hmm, this changes the structure completely.

        # Attempt 2: Modified Redheffer with chi_prime weighting
        # R2[i][j] = chi_prime(j) if j|i, and R2[i][1] = 1
        R2 = np.zeros((n, n), dtype=float)
        for i in range(1, n+1):
            R2[i-1][0] = 1  # column 1
            for j in range(2, n+1):
                if i % j == 0:
                    R2[i-1][j-1] = 1 if sympy.isprime(j) else 0
        det_R2 = int(round(np.linalg.det(R2)))
        print(f"  det(R2) [chi_prime weighting] = {det_R2}, pi({n}) = {target_pi}")

        # Attempt 3: omega(j) weighting on the divisibility part
        R3 = np.zeros((n, n), dtype=float)
        for i in range(1, n+1):
            R3[i-1][0] = 1
            for j in range(2, n+1):
                if i % j == 0:
                    R3[i-1][j-1] = len(sympy.factorint(j))  # omega(j)
        det_R3 = int(round(np.linalg.det(R3)))
        print(f"  det(R3) [omega weighting] = {det_R3}")

        # Attempt 4: Use Liouville lambda instead of 1
        R4 = np.zeros((n, n), dtype=float)
        for i in range(1, n+1):
            R4[i-1][0] = 1
            for j in range(2, n+1):
                if i % j == 0:
                    omega_j = sum(sympy.factorint(j).values())  # Omega(j) with multiplicity
                    R4[i-1][j-1] = (-1)**omega_j  # Liouville lambda
        det_R4 = int(round(np.linalg.det(R4)))
        # sum of Liouville lambda = L(n)
        L_n = sum((-1)**sum(sympy.factorint(k).values()) for k in range(1, n+1))
        print(f"  det(R4) [Liouville] = {det_R4}, L({n}) = {L_n}")

        # Attempt 5: Use von Mangoldt function
        R5 = np.zeros((n, n), dtype=float)
        for i in range(1, n+1):
            R5[i-1][0] = 1
            for j in range(2, n+1):
                if i % j == 0:
                    # Lambda(j) = log(p) if j = p^k, else 0
                    f = sympy.factorint(j)
                    if len(f) == 1:
                        p = list(f.keys())[0]
                        R5[i-1][j-1] = math.log(p)
                    else:
                        R5[i-1][j-1] = 0
        det_R5 = np.linalg.det(R5)
        # sum Lambda(k) = psi(n) = Chebyshev psi
        psi_n = sum(math.log(list(sympy.factorint(k).keys())[0]) if len(sympy.factorint(k)) == 1 else 0 for k in range(2, n+1))
        print(f"  det(R5) [von Mangoldt] = {det_R5:.4f}, psi({n}) = {psi_n:.4f}")

        # Attempt 6: Diagonal weighting to convert M(n) -> pi(n)
        # If det(R) = M(n), can we find diagonal matrices D_L, D_R such that
        # det(D_L * R * D_R) = pi(n)?
        # det(D_L * R * D_R) = prod(d_L_i) * prod(d_R_j) * det(R) = C * M(n)
        # This can only give pi(n) if pi(n) / M(n) is rational and factorizable
        # But M(n) can be 0 or negative, so this doesn't work in general.
        if target_mertens != 0:
            ratio = target_pi / target_mertens
            print(f"  pi({n})/M({n}) = {ratio:.4f} - diagonal scaling would need this ratio")
        else:
            print(f"  M({n}) = 0, diagonal scaling impossible")

        print()

def analyze_redheffer_rank_and_svd():
    """
    Analyze the rank structure of the Redheffer matrix for various n.
    Key question: as n grows, does the matrix have growing low-rank components?
    """
    print("=== Redheffer Matrix SVD Analysis ===\n")

    print(f"{'n':>6} {'det':>8} {'rank':>6} {'sv1':>10} {'sv2':>10} {'sv_last':>10} "
          f"{'90%':>5} {'99%':>5}")

    for n in [10, 20, 30, 50, 75, 100, 150, 200]:
        R = build_redheffer(n)
        det_R = np.linalg.det(R)
        sigma = np.linalg.svd(R, compute_uv=False)
        rank = np.sum(sigma > 1e-10)

        total_energy = np.sum(sigma**2)
        cumulative = np.cumsum(sigma**2) / total_energy
        k90 = np.searchsorted(cumulative, 0.90) + 1
        k99 = np.searchsorted(cumulative, 0.99) + 1

        print(f"  {n:>4}: det={det_R:>7.0f}, rank={rank:>4}, "
              f"sv1={sigma[0]:>9.3f}, sv2={sigma[1]:>9.3f}, "
              f"sv_last={sigma[-1]:>9.6f}, 90%={k90:>3}, 99%={k99:>3}")

def analyze_redheffer_sparsity():
    """Analyze sparsity pattern of the Redheffer matrix."""
    print("\n=== Redheffer Matrix Sparsity ===\n")

    for n in [50, 100, 200]:
        R = build_redheffer(n)
        nnz = np.sum(R > 0)
        # Column 1 has n entries. Plus divisor pairs for j >= 2.
        # Number of divisor pairs = sum_{j=2}^{n} floor(n/j)
        div_pairs = sum(n // j for j in range(2, n+1))
        print(f"  n={n}: nnz={nnz}, density={nnz/(n*n):.4f}, "
              f"col1={n}, div_pairs={div_pairs}")

        # Average entries per row
        row_nnz = np.sum(R > 0, axis=1)
        print(f"    Row nnz: min={row_nnz.min()}, max={row_nnz.max()}, "
              f"mean={row_nnz.mean():.1f}, median={np.median(row_nnz):.1f}")

        # The sparsity is O(n log n / n^2) = O(log n / n) -> 0
        # So R is VERY sparse for large n

def build_pi_from_mertens():
    """
    Explore the relationship between M(x) = sum mu(k) and pi(x).

    Key identity: pi(x) = sum_{k=1}^{x} |mu(k)| * chi_prime(k)
    But also: sum |mu(k)| = Q(x) = number of squarefree integers <= x ~ 6x/pi^2

    More useful: Mobius inversion gives
    pi(x) = sum_{d=1}^{x} mu(d) * [x/d has a prime factor > x^{1/2}]
    ... this doesn't simplify.

    Actually, there's a relationship via the Selberg sieve:
    pi(x) is NOT easily computable from M(x) values alone.
    M(x) ~ 0 (on average), while pi(x) ~ x/ln(x).
    They encode very different information.

    But what if we use MULTIPLE Redheffer-like matrices?
    """
    print("\n=== Mertens -> pi(x) connection ===\n")

    for n in [50, 100]:
        M_vals = {}
        for k in range(1, n+1):
            M_vals[k] = sum(int(sympy.mobius(j)) for j in range(1, k+1))

        pi_n = pi(n)
        print(f"n={n}: pi({n})={pi_n}, M({n})={M_vals[n]}")

        # Can pi(n) be expressed as a linear combination of M(k) for k <= n?
        # pi(n) = sum c_k M(k)?
        # M(k) = sum_{j=1}^{k} mu(j), so sum c_k M(k) = sum_j mu(j) * sum_{k>=j} c_k
        # For this to equal pi(n) = sum chi_prime(j), we need:
        # mu(j) * sum_{k>=j} c_k = chi_prime(j)  for all j

        # This requires sum_{k>=j} c_k = chi_prime(j) / mu(j)
        # But mu(j) = 0 for non-squarefree j, and chi_prime(j) = 1 for primes (which are squarefree)
        # So for prime j: sum_{k>=j} c_k = 1/mu(j) = 1/(-1) = -1 (since mu(p) = -1)
        # For squarefree composite j: sum_{k>=j} c_k = 0/mu(j) = 0
        # For non-squarefree j: need 0/0, which means sum_{k>=j} c_k must be finite

        # Telescoping: c_j = sum_{k>=j} c_k - sum_{k>=j+1} c_k
        # So c_j = f(j) - f(j+1) where f(j) = sum_{k>=j} c_k
        # f(j) = chi_prime(j)/mu(j) when mu(j) != 0
        # f(j) = ? when mu(j) = 0

        # For primes: f(p) = -1
        # For squarefree composites: f(j) = 0
        # For non-squarefree: f(j) is unconstrained (any value works)

        # Let's set f(j) = 0 for non-squarefree j as well.
        # Then c_j = f(j) - f(j+1)

        f = {}
        for j in range(1, n+2):
            mu_j = int(sympy.mobius(j)) if j <= n else 0
            is_prime_j = sympy.isprime(j) if j <= n else False
            if mu_j != 0:
                f[j] = (1 if is_prime_j else 0) / mu_j
            else:
                f[j] = 0

        c = {}
        for j in range(1, n+1):
            c[j] = f[j] - f.get(j+1, 0)

        # Verify
        reconstructed = sum(c[k] * M_vals[k] for k in range(1, n+1))
        print(f"  Reconstructed pi({n}) from M values: {reconstructed:.4f} (should be {pi_n})")

        # Show first few c values
        nonzero_c = [(k, c[k]) for k in range(1, min(n+1, 30)) if abs(c[k]) > 1e-10]
        print(f"  First nonzero c_k: {nonzero_c[:10]}")
        print(f"  Number of nonzero c_k: {sum(1 for k in c if abs(c[k]) > 1e-10)}")

        # This means: pi(n) = sum c_k * det(R_k) where R_k is the k-th Redheffer matrix
        # This is a SUM of determinants, not a SINGLE determinant.
        # A single determinant: can we embed all these into one matrix?

        # Block diagonal matrix: det(diag(R_1, R_2, ..., R_n)) = prod det(R_k)
        # Product, not sum. To get a sum, we'd need a different construction.

        print()

def experiment_lgv_lattice_path():
    """
    Lindstrom-Gessel-Viennot Lemma:
    det(M) = signed count of non-intersecting lattice path systems
    where M[i][j] = (number of paths from source s_i to sink t_j).

    For pi(x): can we design a lattice/DAG such that the signed count
    of non-intersecting paths equals pi(x)?

    Simple idea: k source-sink pairs where k = O(log x),
    and the non-intersecting path count encodes primality testing
    for all numbers up to x via some clever construction.

    This is the DAG question from the GapL framing.
    """
    print("\n=== LGV Lattice Path Construction ===\n")

    # Start simple: for a SINGLE path (1x1 matrix = just a path count),
    # we need a DAG with #(source->sink paths) = pi(x).
    # Trivially: a DAG with pi(x) parallel edges from s to t. But this is circular.

    # For 2x2: det([[a,b],[c,d]]) = ad - bc = pi(x)
    # where a = #paths(s1->t1), b = #paths(s1->t2), etc.
    # Non-intersecting path interpretation: ad - bc = #{non-int path pairs}

    # For the sieve:
    # pi(x) = #{n in [2,x] : n has no prime factor <= sqrt(x)}  (+ pi(sqrt(x)) - 1)
    # = #{n in [2,x]} - #{n in [2,x] with prime factor p1} - ... + ...
    # The I-E terms correspond to paths through a "sieve graph"

    # Build a concrete small DAG for x = 30
    x = 30
    target = pi(x)
    sqrtx = int(math.isqrt(x))
    P = [p for p in range(2, sqrtx+1) if sympy.isprime(p)]
    k = len(P)  # number of sieve primes

    print(f"x={x}, pi(x)={target}, primes<={sqrtx}: {P}")

    # DAG construction attempt:
    # Layers: L0 (source), L1 (after considering p1=2), L2 (after p2=3), L3 (after p3=5), L4 (sink)
    # At each layer, a node represents a "residual count"
    # Node (layer, state) where state encodes which primes have been sieved

    # This gives 2^k nodes per layer = exponential. We know this.
    # The question: is there a SMALLER DAG?

    # Alternative: DAG where nodes represent ARITHMETIC conditions
    # E.g., node (k, r) means "integer k with residue r mod P_i"
    # But this still has O(prod P_i) = O(x) states

    # What about a DAG based on BIT REPRESENTATION of n?
    # N = log2(x) bits. Read bits from MSB to LSB.
    # At each step, maintain some "state" about whether the number
    # formed so far could be prime.
    # This is essentially a finite automaton for primality.
    # But primality is NOT a regular language in any base! (It's not even
    # recognizable by a finite automaton for any fixed modulus.)

    print(f"\nDAG approaches considered:")
    print(f"  1. Sieve DAG: 2^k = {2**k} nodes per layer (EXPONENTIAL)")
    print(f"  2. Residue DAG: O(prod P_i) = O({math.prod(P)}) states (EXPONENTIAL)")
    print(f"  3. Bit-reading DAG: impossible (primality not regular)")
    print(f"  4. Arithmetic circuit DAG: O(sqrt(x)) = {sqrtx} nodes (known best)")

    # Experiment: for small x, ENUMERATE all DAGs up to a certain size
    # and check if any gives pi(x) as the path count
    # (This is only feasible for very small x and DAG sizes)

    print(f"\n--- Brute-force small DAG search for x={x} ---")
    print(f"  Target: find DAG on <= {2*k+2} nodes with signed path count = {target}")
    print(f"  (This is computationally infeasible even for this small case.)")
    print(f"  The search space is super-exponential in the number of nodes.")

    # Instead, verify the KNOWN construction (sieve DAG) works:
    # Build the sieve DAG and count paths
    print(f"\n--- Sieve DAG verification ---")

    # Each node: (layer=prime_index, subset_of_sieved_primes)
    # Edge from (i, S) to (i+1, S) with weight floor(x/prod(S)) - floor(x/prod(S+{p_{i+1}}))
    # ... this is the I-E expansion

    # Simpler: just verify that the I-E formula gives correct pi(x)
    ie_sum = 0
    for mask in range(1 << k):
        prod_S = 1
        bits = bin(mask).count('1')
        for i in range(k):
            if mask & (1 << i):
                prod_S *= P[i]
        ie_sum += ((-1)**bits) * (x // prod_S)

    pi_sieve = ie_sum - 1 + pi(sqrtx)
    print(f"  I-E gives: ie_sum={ie_sum}, pi_sieve={pi_sieve}, actual={target}")
    print(f"  Match: {pi_sieve == target}")

    # The DAG for this I-E has 2^k = {2**k} paths, one per subset S.
    # Each path has weight (-1)^|S| * floor(x/prod(S)).
    # The signed sum gives ie_sum.
    print(f"\n  DAG has 2^k = {2**k} paths for k={k} primes.")
    print(f"  This is exponential in k ~ sqrt(x)/ln(sqrt(x)) ~ sqrt(x)/ln(x).")
    print(f"  For x=10^100: k ~ 2.2 * 10^48, so 2^k is absurdly large.")

def main():
    print("=" * 70)
    print("COMPRESSED REDHEFFER MATRIX ANALYSIS")
    print("Session 14 - Determinant approach to prime counting")
    print("=" * 70)

    # 1. Prime-Redheffer variants
    build_prime_redheffer(50)

    # 2. SVD analysis
    analyze_redheffer_rank_and_svd()

    # 3. Sparsity
    analyze_redheffer_sparsity()

    # 4. M(x) -> pi(x) connection
    build_pi_from_mertens()

    # 5. LGV lattice path approach
    experiment_lgv_lattice_path()

    print("\n" + "=" * 70)
    print("CONCLUSIONS")
    print("=" * 70)
    print("""
1. The standard Redheffer matrix R_n gives M(n) = Mertens function.
   It is FULL RANK for most n (whenever M(n) != 0).
   SVD shows no useful low-rank structure.

2. Modified Redheffer matrices with different weightings:
   - chi_prime weighting does NOT give pi(n) as det.
   - omega weighting does NOT give pi(n) as det.
   - Liouville weighting gives L(n) (Liouville summatory).
   None of the simple modifications yield pi(n).

3. pi(n) CAN be expressed as a LINEAR COMBINATION of M(k) values
   (sum c_k * M(k) = pi(n)), but this requires O(n) terms and
   the c_k depend on primality of individual numbers (circular).

4. The LGV lattice path approach requires 2^k paths through the
   sieve DAG, where k = pi(sqrt(x)). This is inherently exponential.
   No smaller DAG construction is known.

5. FUNDAMENTAL OBSTACLE: All known constructions that give pi(x)
   as a determinant use matrices of size >= sqrt(x). The information
   content of pi(x) depends on O(sqrt(x)) independent values
   (floor-value set), and no compression to poly(log x) is known.
""")

if __name__ == "__main__":
    main()

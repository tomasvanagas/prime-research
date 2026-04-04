"""
Session 15: Can pi(x) be expressed as det/perm of a poly(N)-size matrix?

GapL question: pi(x) = det(A(x)) where A is m×m with m = poly(N), N = log2(x)?

Known:
- Redheffer matrix: det(R_n) = M(n), but R is n×n = 2^N × 2^N (exponential)
- Session 14: floor-value-based matrices need exponential size
- Session 14: Lucy DP matrices are unipotent, full-rank product, no structure

NEW ANGLES:
1. Can we express pi(x) as det of a matrix whose entries are BITS of x?
2. Can we find a poly(N)-size matrix using modular arithmetic?
3. Smith normal form / invariant factors approach?
4. Can we search for small matrices empirically?

Experiment 1: For small N, exhaustively search for matrices with det = pi(x).
Experiment 2: Analyze what constraints the matrix entries must satisfy.
Experiment 3: Test if pi(x) has a poly(N)-size PERMANENT representation.
"""

import numpy as np
from sympy import primepi, Matrix, isprime
from itertools import product as iproduct
import math

def experiment1_small_det_search():
    """
    For x = 2^N with small N, search for m×m matrices A with entries
    that are simple functions of the bits of x, such that det(A) = pi(x).

    Start with m=2 (2x2 matrix, 4 entries), N up to 8.
    Each entry is a linear function of the bits of x.
    """
    print("=" * 60)
    print("Experiment 1: Small determinant search")
    print("=" * 60)

    # For m=2: det = a*d - b*c
    # Need: a*d - b*c = pi(x) for all x in some range
    # Each of a,b,c,d is a linear function of bits of x

    # For N bits, each entry = c0 + c1*x1 + c2*x2 + ... + cN*xN
    # where xi are the bits

    # This is a system of 2^N equations (one per x) in 4*(N+1) unknowns
    # For it to have a solution, 4*(N+1) >= 2^N... at N=4, 20 >= 16. Possible!

    for N in range(3, 9):
        x_max = 2**N
        xs = list(range(1, x_max + 1))
        pi_vals = [int(primepi(x)) for x in xs]

        # Bits of x: x1 (LSB), x2, ..., xN
        def bits(x, N):
            return [(x >> i) & 1 for i in range(N)]

        # For 2x2 det: ad - bc = pi(x)
        # Let each of a,b,c,d be affine in bits: f(x) = w0 + sum wi*xi
        # Then ad - bc is degree 2 in bits — a quadratic form

        # Approach: pi(x) must be expressible as a quadratic form in bits
        # Build the quadratic model: pi(x) = sum_{i<=j} q_{ij} * xi * xj + sum_i l_i * xi + c

        bit_matrix = np.array([bits(x, N) for x in xs])  # (2^N) x N

        # Build quadratic features
        features = [np.ones(len(xs))]  # constant
        for i in range(N):
            features.append(bit_matrix[:, i])  # linear
        for i in range(N):
            for j in range(i, N):
                features.append(bit_matrix[:, i] * bit_matrix[:, j])  # quadratic

        X = np.column_stack(features)
        y = np.array(pi_vals, dtype=float)

        # Solve least squares
        result = np.linalg.lstsq(X, y, rcond=None)
        coeffs = result[0]
        residual = y - X @ coeffs
        max_error = np.max(np.abs(residual))

        # Check if ANY quadratic form in bits exactly represents pi(x)
        exact = max_error < 0.01

        n_features = X.shape[1]
        print(f"N={N}: x_max={x_max:4d}, pi(x_max)={pi_vals[-1]:4d}, "
              f"features={n_features:4d}, max_residual={max_error:.4f}, "
              f"exact_quadratic={'YES' if exact else 'NO'}")

        if not exact and N <= 5:
            # Try cubic
            cubic_features = list(features)
            for i in range(N):
                for j in range(i, N):
                    for k in range(j, N):
                        cubic_features.append(bit_matrix[:, i] * bit_matrix[:, j] * bit_matrix[:, k])

            X3 = np.column_stack(cubic_features)
            result3 = np.linalg.lstsq(X3, y, rcond=None)
            residual3 = y - X3 @ result3[0]
            max_error3 = np.max(np.abs(residual3))
            exact3 = max_error3 < 0.01
            print(f"  -> cubic: features={X3.shape[1]:4d}, max_residual={max_error3:.4f}, "
                  f"exact={'YES' if exact3 else 'NO'}")


def experiment2_det_rank_analysis():
    """
    The matrix A with det(A) = pi(x) must satisfy:
    - det is a polynomial of degree m in the entries
    - entries are functions of x (or bits of x)
    - So det(A(x)) is a polynomial of degree m*d in bits (if entries have degree d)

    We know: prime indicator has ANF degree Theta(N) over GF(2).
    pi(x) = sum of prime indicators = a polynomial in bits of x.

    What is the degree of pi(x) as a REAL polynomial in bits of x?
    If degree = d, then we need m*entry_degree >= d for a det representation.
    """
    print("\n" + "=" * 60)
    print("Experiment 2: Degree of pi(x) as polynomial in bits")
    print("=" * 60)

    for N in range(3, 13):
        x_max = 2**N
        xs = list(range(1, x_max + 1))
        pi_vals = np.array([int(primepi(x)) for x in xs], dtype=float)

        def bits(x, N):
            return [(x >> i) & 1 for i in range(N)]

        bit_matrix = np.array([bits(x, N) for x in xs])

        # Try polynomials of increasing degree
        best_degree = None
        for degree in range(1, min(N+1, 16)):
            # Build monomial basis up to given degree
            # For efficiency, use itertools for small N
            if N > 10 and degree > 3:
                break  # too many features

            from itertools import combinations_with_replacement
            monomials = []
            for d in range(degree + 1):
                for combo in combinations_with_replacement(range(N), d):
                    mono = np.ones(len(xs))
                    for idx in combo:
                        mono *= bit_matrix[:, idx]
                    monomials.append(mono)

            X = np.column_stack(monomials)

            if X.shape[1] > X.shape[0]:
                # More features than samples — always fits
                result = np.linalg.lstsq(X, pi_vals, rcond=None)
                residual = pi_vals - X @ result[0]
                max_err = np.max(np.abs(residual))
                if max_err < 0.1:
                    best_degree = degree
                    break
                continue

            result = np.linalg.lstsq(X, pi_vals, rcond=None)
            residual = pi_vals - X @ result[0]
            max_err = np.max(np.abs(residual))

            if max_err < 0.1:
                best_degree = degree
                break

        if best_degree:
            print(f"N={N:2d}: pi(x) is degree-{best_degree} polynomial in bits "
                  f"(over R, {len(monomials)} monomials)")
        else:
            print(f"N={N:2d}: pi(x) requires degree > {degree} in bits")


def experiment3_permanent_representation():
    """
    Can pi(x) = perm(B) for some small 0/1 matrix B(x)?

    The permanent of an m×m 0/1 matrix counts perfect matchings.
    Range: 0 to m! (since it's the permanent of a 0/1 matrix).

    For pi(x) ~ x/ln(x), we need m such that m! >= pi(x), i.e., m >= O(N).
    But can we actually CONSTRUCT such a matrix?

    Approach: For each small x, find the smallest 0/1 matrix with perm = pi(x).
    """
    print("\n" + "=" * 60)
    print("Experiment 3: Permanent representation")
    print("=" * 60)

    from sympy.combinatorics.permutations import Permutation

    def permanent(M):
        """Compute permanent of a small matrix by Ryser formula"""
        n = M.shape[0]
        if n == 0:
            return 1
        if n == 1:
            return M[0, 0]

        # Ryser's formula
        result = 0
        for subset in range(1, 2**n):
            col_sum = np.zeros(n)
            sign = 1
            bits_set = 0
            for j in range(n):
                if subset & (1 << j):
                    col_sum += M[:, j]
                    bits_set += 1

            prod = 1
            for i in range(n):
                prod *= col_sum[i]

            if bits_set % 2 == 0:
                result -= prod
            else:
                result += prod

        if n % 2 == 0:
            result = -result
        return int(round(result))

    # For small values of pi(x), find minimum matrix size
    test_values = [(x, int(primepi(x))) for x in [8, 16, 32, 64, 128, 256]]

    for x, target in test_values:
        N = int(math.log2(x))
        # Minimum m: m! >= target, so m >= ...
        m = 1
        while math.factorial(m) < target:
            m += 1

        print(f"x={x:4d} (N={N:2d}): pi(x)={target:4d}, min m for m!≥pi(x): {m}")

        # For very small targets, try to find actual 0/1 matrices
        if target <= 20:
            found = False
            for size in range(2, 6):
                # Random search for 0/1 matrix with perm = target
                for trial in range(10000):
                    M = np.random.randint(0, 2, (size, size))
                    if permanent(M) == target:
                        print(f"  Found {size}x{size} matrix with perm={target}")
                        found = True
                        break
                if found:
                    break
            if not found:
                print(f"  No small (≤5×5) 0/1 matrix found with perm={target}")


def experiment4_det_from_bits():
    """
    Direct search: for N=4,5,6, can we find a 2x2 or 3x3 matrix whose entries
    are polynomials in bits of x, with det = pi(x)?

    For 2x2: det = ad - bc. Need to express pi(x) as a product minus product.
    For 3x3: det = a(ei-fh) - b(di-fg) + c(dh-eg).
    """
    print("\n" + "=" * 60)
    print("Experiment 4: Searching for det = pi(x) via bit polynomials")
    print("=" * 60)

    for N in range(4, 8):
        x_max = 2**N
        xs = list(range(1, x_max + 1))
        pi_vals = np.array([int(primepi(x)) for x in xs])

        def bits(x, N):
            return [(x >> i) & 1 for i in range(N)]

        bit_matrix = np.array([bits(x, N) for x in xs])

        # For 2x2 det: need pi(x) = f1(bits)*f2(bits) - f3(bits)*f4(bits)
        # where each fi is a low-degree polynomial in bits

        # Equivalently: pi(x) = P(bits) - Q(bits) where P,Q are products of
        # two linear forms each.

        # Try: a = sum ai*xi + a0, d = sum di*xi + d0, etc.
        # det = (a0 + sum ai xi)(d0 + sum di xi) - (b0 + sum bi xi)(c0 + sum ci xi)

        # This is a rank-2 quadratic form (difference of two rank-1 forms)
        # Can pi(x) be expressed as a rank-2 quadratic form in bits?

        # Build all degree-2 monomials
        from itertools import combinations_with_replacement
        monos_1 = [np.ones(len(xs))]  # constant
        for i in range(N):
            monos_1.append(bit_matrix[:, i])

        monos_2 = []
        for i in range(N):
            for j in range(i, N):
                monos_2.append(bit_matrix[:, i] * bit_matrix[:, j])

        # pi(x) as quadratic: pi(x) = c + sum li*xi + sum qij*xi*xj
        X_full = np.column_stack(monos_1 + monos_2)
        result = np.linalg.lstsq(X_full, pi_vals.astype(float), rcond=None)
        fitted = X_full @ result[0]
        max_err = np.max(np.abs(fitted - pi_vals))

        # Check if the quadratic form has rank 2
        # Extract the quadratic coefficient matrix Q
        n_lin = N + 1  # constant + N linear
        q_coeffs = result[0][n_lin:]

        Q = np.zeros((N, N))
        idx = 0
        for i in range(N):
            for j in range(i, N):
                if i == j:
                    Q[i, j] = q_coeffs[idx]
                else:
                    Q[i, j] = q_coeffs[idx] / 2
                    Q[j, i] = q_coeffs[idx] / 2
                idx += 1

        rank_Q = np.linalg.matrix_rank(Q, tol=0.01)
        eigenvalues = np.linalg.eigvalsh(Q)
        pos_eig = np.sum(eigenvalues > 0.01)
        neg_eig = np.sum(eigenvalues < -0.01)

        print(f"N={N}: quadratic max_err={max_err:.4f}, "
              f"Q rank={rank_Q}, eigenvalues: {pos_eig}+/{neg_eig}- "
              f"({'EXACT' if max_err < 0.01 else 'NOT exact'})")
        if max_err < 0.01 and rank_Q <= 4:
            print(f"  -> pi(x) IS a rank-{rank_Q} quadratic form in bits!")
            print(f"  -> This gives a {rank_Q}x{rank_Q} determinant representation!")


def experiment5_characteristic_polynomial():
    """
    Alternative approach: instead of det(A) = pi(x) directly,
    can we find a matrix A(x) whose CHARACTERISTIC POLYNOMIAL
    has pi(x) as a coefficient or root?

    Or: can pi(x) be a TRACE of a matrix power? tr(A^k) = sum eigenvalues^k.
    Session 14 showed pi(x) is NOT an LRS, but maybe a TRACE of a
    VARIABLE matrix (A depends on x)?
    """
    print("\n" + "=" * 60)
    print("Experiment 5: Trace/characteristic polynomial approaches")
    print("=" * 60)

    # For small x, check: is there a small matrix A(x) with tr(A) = pi(x)?
    # Trivially yes (1x1 matrix [pi(x)]). But can entries be simple functions of x?

    # More interesting: tr(A^x) = pi(x)?
    # This would make pi(x) = sum lambda_i^x, a linear recurrence.
    # Session 14 proved this is NOT possible. Confirmed.

    # What about tr(A(x)) where A(x) is a structured matrix depending on x?
    # E.g., A(x) = sum_{k} x_k * B_k for basis matrices B_k?
    # Then tr(A(x)) = sum_k x_k * tr(B_k) = linear in bits = NOT pi(x)
    # (pi(x) is not linear in bits)

    # What about det(xI - A) for a FIXED matrix A?
    # This is the characteristic polynomial of A. Evaluating at x gives a polynomial
    # of degree m (the matrix size). pi(x) is not a polynomial in x (it's a step function).
    # So this can't work directly.

    # Key insight: pi(x) is a STEP FUNCTION of x. It increases by 1 at each prime.
    # As a function of the INTEGER x, it can be expressed as:
    # pi(x) = sum_{p prime, p <= x} 1 = sum_{n=2}^x [n is prime]

    # This is NOT a polynomial in x. But it IS a polynomial in the BITS of x
    # (since any Boolean function of N bits is a multilinear polynomial over R).

    print("Key insight: pi(x) is a multilinear polynomial of degree N in bits of x.")
    print("Any Boolean function on N bits has a unique multilinear representation.")
    print("The degree of this polynomial = N (from GF(2) ANF degree = N, Session 13).")
    print()
    print("For det representation: det of m×m matrix with degree-d entries")
    print("gives degree m*d polynomial. Need m*d >= N.")
    print("If entries are LINEAR in bits (d=1): need m >= N.")
    print("If entries are QUADRATIC (d=2): need m >= N/2.")
    print("If entries are CONSTANT (d=0): impossible (det is constant).")
    print()
    print("The question: does the SPECIFIC polynomial pi(bits) have a")
    print("small determinantal representation?")
    print()

    # Determinantal complexity: the smallest m such that pi(x) = det(L)
    # where L is m×m with affine-linear entries in bits of x.
    # This is a well-studied question in algebraic complexity!

    # Known: for the permanent polynomial of n×n matrix,
    #   det complexity >= n^2/2 (Mignon-Ressayre 2004)
    #   det complexity = 2^{Theta(n)} conjectured

    # For pi(x): we have a polynomial of degree N in N variables.
    # What is its determinantal complexity?

    print("DETERMINANTAL COMPLEXITY of pi(x):")
    print("This is exactly the algebraic complexity question!")
    print("If det_complexity(pi) = poly(N), then pi(x) is in GapL.")
    print("If det_complexity(pi) = 2^{Omega(N)}, then pi(x) is NOT in GapL.")


if __name__ == "__main__":
    experiment1_small_det_search()
    experiment2_det_rank_analysis()
    experiment3_permanent_representation()
    experiment4_det_from_bits()
    experiment5_characteristic_polynomial()

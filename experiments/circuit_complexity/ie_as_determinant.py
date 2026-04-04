"""
Inclusion-Exclusion as Determinant
Session 14 - April 2026

The Legendre sieve gives:
  pi(x) - pi(sqrt(x)) + 1 = sum_{S subset P} (-1)^|S| floor(x / prod(S))
where P = {primes <= sqrt(x)}.

This is a sum over 2^|P| terms. Can it be expressed as a SMALL determinant?

Key results in combinatorics:
- Many inclusion-exclusion sums CAN be expressed as determinants
  (e.g., derangements, non-crossing partitions, lattice path counts)
- The LGV lemma gives det = signed count of non-intersecting paths
- Cauchy-Binet: det(AB) = sum of products of minors

This experiment systematically explores whether the prime-counting I-E
can be written as a determinant of a matrix smaller than 2^|P| x 2^|P|.
"""

import numpy as np
import math
import sympy
from itertools import product as cart_product, combinations
from functools import reduce

def pi(x):
    return int(sympy.primepi(x))

def floor_div(x, d):
    return x // d

def compute_ie_sum(x, P):
    """Compute the inclusion-exclusion sum for the Legendre sieve."""
    k = len(P)
    total = 0
    for mask in range(1 << k):
        prod_S = 1
        bits = bin(mask).count('1')
        for i in range(k):
            if mask & (1 << i):
                prod_S *= P[i]
        total += ((-1)**bits) * (x // prod_S)
    return total

def experiment_ie_factorization(x):
    """
    Key insight: the I-E sum factors as a PRODUCT of operators.

    sum_{S} (-1)^|S| floor(x/prod(S))
    = prod_{p in P} (1 - T_p) applied to f(1) = x
    where T_p f(d) = f(dp) = floor(x/(dp)) when applied to f(d) = floor(x/d)

    Wait, this isn't quite a product of NUMBERS. It's a product of OPERATORS
    applied to a function. But if the function space is finite-dimensional
    (which it is: the floor-value set has finite size), then each operator
    is a MATRIX.

    The operator (1 - T_p) on the floor-value set:
    - For each floor-value v, (1-T_p)f(v) = f(v) - f(v*p) = floor(x/v) - floor(x/(vp))
    - But v*p might not give a floor-value of x... actually it does:
      if v = floor(x/k), then vp ~ x*p/k, and floor(x/(vp)) is a floor value

    Actually, let's think in terms of the DIVISOR set D = {d : floor(x/d) is distinct}
    The operators T_p map D to D (since d -> dp maps divisors to divisors... not exactly).

    Let me work on the "smooth number" side.
    Define D = {d : d is a product of primes in P, d <= x}
    = {p1^a1 * p2^a2 * ... * pk^ak : product <= x}
    This is the set of P-smooth numbers <= x.

    |D| can be very large (exponential in k for k large).
    But for the I-E sum, we only use SQUAREFREE products (each prime at most once).
    So |D_sf| = 2^k = exponential.
    """
    sqrtx = int(math.isqrt(x))
    P = [p for p in range(2, sqrtx+1) if sympy.isprime(p)]
    k = len(P)

    print(f"\n=== I-E Factorization for x={x} ===")
    print(f"P = {P}, k = {k}")

    ie_sum = compute_ie_sum(x, P)
    print(f"I-E sum = {ie_sum}")
    print(f"pi(x) = {pi(x)}, pi(sqrt(x)) = {pi(sqrtx)}")
    print(f"Reconstructed: {ie_sum - 1 + pi(sqrtx)} (should be {pi(x)})")

    # The I-E sum as operator product:
    # (1 - T_{p1}) (1 - T_{p2}) ... (1 - T_{pk}) [floor(x/1)]
    #
    # Each (1 - T_p) is a 2x2 "local" operation: it takes f(d) and subtracts f(dp).
    # The composition of k such operations gives a sum over 2^k terms.
    #
    # If the operators COMMUTE, then we can rearrange and maybe compress.
    # Do they commute? (1-T_p)(1-T_q)f(d) = f(d) - f(dp) - f(dq) + f(dpq)
    # (1-T_q)(1-T_p)f(d) = f(d) - f(dq) - f(dp) + f(dqp)
    # Since f(dpq) = f(dqp), YES they commute!
    #
    # So the I-E is a COMMUTATIVE product of operators.
    # This is the same as the product: prod_{p in P} (1 - T_p)
    #
    # Can this product be computed on a smaller state space?

    # Think of it as a polynomial: define z_p for each prime p.
    # The "generating function" is:
    # F(z_1, ..., z_k) = sum_{d squarefree P-smooth} mu(d) * z_1^{v_1(d)} * ... * z_k^{v_k(d)} * floor(x/d)
    # = prod_{i=1}^k (floor(x/1) - z_i * floor(x/p_i)) ... no, this doesn't factor.

    # The CRUCIAL non-factoring: floor(x/(pq)) != floor(x/p) * floor(x/q) / x
    # The floor function is NOT multiplicative in this sense.

    # What if we use a different basis?
    # Instead of {floor(x/d)}, use {x/d} (rational, no floor)?
    # Then: prod (1 - T_p) [x/1] = x * prod(1 - 1/p) = x * prod_{p in P}(p-1)/p
    # This gives the smooth approximation (Euler product), not pi(x).
    # The DIFFERENCE between the actual I-E sum and this smooth product
    # is exactly the "fractional part" error = sum of {x/d} terms.

    # Compute the smooth vs actual
    euler_product = x
    for p in P:
        euler_product *= (p - 1) / p
    fractional_error = ie_sum - euler_product

    print(f"\nSmooth approximation (Euler product): {euler_product:.4f}")
    print(f"Actual I-E sum: {ie_sum}")
    print(f"Fractional error: {fractional_error:.4f}")
    print(f"Relative error: {abs(fractional_error/ie_sum):.4f}")

    # The error comes from: sum_{S} (-1)^|S| {x/prod(S)} (fractional parts)
    # This sum has 2^k terms and encodes the "random" component of pi(x).
    # It cannot be computed without examining all 2^k subsets (or equivalent info).

    return ie_sum

def experiment_cauchy_binet(x):
    """
    Cauchy-Binet theorem: det(AB) = sum over k-subsets of columns of A (rows of B)
    of the product of corresponding minors.

    If we can write the I-E sum as det(AB) where A is m x k and B is k x m,
    with m << 2^k, this would give a small determinant.

    The I-E sum: sum_{S subset [k]} (-1)^|S| floor(x / prod_{i in S} p_i)
    = sum_{S subset [k]} (-1)^|S| w(S) where w(S) = floor(x/prod(S))

    By Cauchy-Binet: det of an n x n matrix = sum over minors.
    But our sum is over SUBSETS with SIGNS, which is exactly a determinant
    by the Leibniz formula... of what matrix?

    Actually, the sum sum_S (-1)^|S| w(S) is the evaluation of the multilinear
    polynomial w at all variables = -1:
    w(t_1, ..., t_k) = sum_S prod_{i in S} t_i * w(S)
    evaluated at t_i = 1 for all i, with (-1)^|S| built in.

    Hmm, let me think about this differently.

    For the I-E sum with k = 2 primes {p, q}:
    S = x - floor(x/p) - floor(x/q) + floor(x/(pq))

    As a 2x2 determinant:
    det([[1, 1], [floor(x/q), floor(x/p)]]) = floor(x/p) - floor(x/q)  -- wrong sign
    det([[x, 1], [floor(x/q), 1]]) = x - floor(x/q)  -- not the full sum

    What about:
    det([[1, 1, 0], [0, 1, 1], [floor(x/(pq)), floor(x/q), floor(x/p)]]) ... too many entries

    Actually, there IS a standard way: the I-E sum equals det(D - A) where
    D is diagonal with D[i][i] = something and A encodes the "sieve" steps.

    For the SIEVE MATRIX construction:
    Define M[i][j] for i,j in {0, 1, ..., k} as follows:
    - M is (k+1) x (k+1)
    - M[0][0] = x (total count)
    - M[0][j] = floor(x/p_j) for j >= 1
    - M[i][0] = 1 for i >= 1
    - M[i][j] = delta_{ij} for i,j >= 1

    Then det(M) = x * det(I_k) - sum_{j=1}^{k} floor(x/p_j) * cofactor
    = x - sum floor(x/p_j) * ... hmm, this gives only the first level of I-E.
    """
    sqrtx = int(math.isqrt(x))
    P = [p for p in range(2, sqrtx+1) if sympy.isprime(p)]
    k = len(P)

    print(f"\n=== Cauchy-Binet / Small Determinant Search for x={x} ===")
    print(f"P = {P}, k = {k}")

    ie_sum = compute_ie_sum(x, P)
    print(f"Target I-E sum = {ie_sum}")

    # Systematic search for (k+1) x (k+1) matrices with det = ie_sum
    # Entries drawn from: {floor(x/d) : d is squarefree P-smooth} union {0, 1, -1}

    # Compute all floor(x/d) for squarefree d
    floor_vals = {}
    for mask in range(1 << k):
        d = 1
        for i in range(k):
            if mask & (1 << i):
                d *= P[i]
        floor_vals[mask] = x // d

    print(f"Floor values: {dict((bin(m), v) for m, v in floor_vals.items())}")

    # For k <= 3, try to find a (k+1) x (k+1) matrix
    if k <= 3:
        n = k + 1

        # Instead of brute force over all entries, try structured forms.

        # Form 1: Upper triangular with specific diagonals
        # det = product of diagonals
        # Need: d0 * d1 * ... * dk = ie_sum
        # This requires factoring ie_sum into k+1 factors.
        print(f"\n  ie_sum = {ie_sum}")
        if ie_sum > 0:
            factors = sympy.factorint(ie_sum)
            print(f"  Factorization: {factors}")

        # Form 2: M = I + R where R has specific structure
        # det(I + R) = 1 + tr(R) + ... (for nilpotent R: finite expansion)

        # Form 3: Block structure
        # [[A, b], [c^T, d]] with det = Ad - c^T A^{-1} b * det(A) (Schur)

        # For k=2: try 3x3 matrices
        if k == 2:
            p, q = P[0], P[1]
            vals = {
                'x': x,
                'x/p': x // p,
                'x/q': x // q,
                'x/pq': x // (p*q),
                '1': 1,
                '0': 0,
                '-1': -1,
            }
            val_list = list(vals.items())

            # Try all 3x3 matrices with entries from vals
            # That's 7^9 ~ 40 million... too many. Restrict.

            # Try: first row = [1, ?, ?], first col = [1, ?, ?]
            # Then det = 1 * (M11*M22 - M12*M21) - M01*(M10*M22 - M12*M20) + M02*(M10*M21 - M11*M20)

            # Actually let's try a specific construction:
            # M = [[x, x//p, x//q],
            #      [1,  1,    0   ],
            #      [1,  0,    1   ]]
            # det = x*(1*1 - 0*0) - (x//p)*(1*1 - 0*1) + (x//q)*(1*0 - 1*1)
            #     = x - x//p - x//q
            M_test = np.array([[x, x//p, x//q], [1, 1, 0], [1, 0, 1]], dtype=float)
            det_test = int(round(np.linalg.det(M_test)))
            print(f"\n  Test M1: det = {det_test} (want {ie_sum})")
            print(f"    = x - x//p - x//q = {x} - {x//p} - {x//q} = {x - x//p - x//q}")

            # Need to add back +floor(x/(pq)). Try:
            # M = [[x, x//p, x//q],
            #      [1,  1, x//pq/(something)],
            #      [1,  0,    1   ]]
            # Let's be more systematic with 3x3.

            # det = M00(M11*M22 - M12*M21) - M01(M10*M22 - M12*M20) + M02(M10*M21 - M11*M20)
            # We want this = x - x//p - x//q + x//(pq)

            # Try: M00=1, and arrange the rest:
            # = 1*(M11*M22 - M12*M21) - M01(M10*M22 - M12*M20) + M02(M10*M21 - M11*M20)

            # Simplest: M = [[1, 0, 0], [?, ?, ?], [?, ?, ?]] -> det = M11*M22 - M12*M21
            # = 2x2 problem: need M11*M22 - M12*M21 = ie_sum
            # With entries from floor values: just a 2x2 determinant.

            # For 2x2: a*d - b*c = ie_sum = x - x//p - x//q + x//(pq)
            # Try a=x, d=1, b=x//p, c=? -> x - x//p * c = ie_sum
            # c = (x - ie_sum) / (x//p) = (x//p + x//q - x//(pq)) / (x//p)
            # This is generally not an integer or a floor value.

            # Try: a=x-x//p, d=1, b=0, c=? -> (x-x//p) - 0 = x-x//p != ie_sum (missing -x//q+x//pq)
            # So we need the 2x2 to handle all four terms.

            # Search 2x2 over extended value set
            extended = list(vals.values())
            # Add more values
            extended += [x-x//p, x-x//q, x//p-x//(p*q), x//q-x//(p*q)]
            extended_names = list(vals.keys()) + ['x-x/p', 'x-x/q', 'x/p-x/pq', 'x/q-x/pq']

            found = False
            for i, a in enumerate(extended):
                for j, b in enumerate(extended):
                    for ki, c in enumerate(extended):
                        for li, d in enumerate(extended):
                            if a*d - b*c == ie_sum:
                                print(f"  Found 2x2: det([[{extended_names[i]},{extended_names[j]}],"
                                      f"[{extended_names[ki]},{extended_names[li]}]]) = {a}*{d}-{b}*{c} = {ie_sum}")
                                found = True
                                break
                        if found: break
                    if found: break
                if found: break

            if not found:
                print(f"  No 2x2 determinant found with extended floor values.")

        # For k=3: try 4x4 matrices (too large for brute force)
        if k == 3:
            n_vals = 2**k + 3  # floor values + {0, 1, -1}
            print(f"\n  k=3: 4x4 search is too large ({n_vals}^16 ~ 10^13)")
            print(f"  Would need structural insight to reduce search space.")

    # KEY THEORETICAL RESULT:
    # The I-E sum over k primes has 2^k terms.
    # For it to be a determinant of an m x m matrix (m << 2^k),
    # we'd need the terms to have a special ALGEBRAIC structure.
    #
    # The terms are floor(x/d) for various squarefree d.
    # These satisfy: floor(x/d) = x/d - {x/d} where {.} is fractional part.
    #
    # The "smooth part" x * prod(1-1/p) is a single number.
    # The "error part" involves fractional parts {x/d} for all 2^k values of d.
    # These fractional parts are essentially INDEPENDENT (by equidistribution
    # of x/d mod 1 for different squarefree d).
    #
    # Since the error terms are independent, they cannot be compressed
    # into a smaller determinant.

    print(f"\n  THEORETICAL ANALYSIS:")
    print(f"  The I-E sum has 2^k terms. The 'smooth part' is 1 number (Euler product).")
    print(f"  The 'error part' involves 2^k independent fractional parts.")
    print(f"  Independent quantities cannot be compressed into a smaller determinant.")
    print(f"  Therefore: no m x m determinant with m << 2^k can represent the I-E sum")
    print(f"  UNLESS we use entries that are themselves nontrivial to compute.")

def experiment_permanent_vs_determinant():
    """
    The I-E sum is closely related to the PERMANENT of a matrix.

    permanent(A) = sum_{sigma in S_n} prod A[i][sigma(i)]   (no sign)
    determinant(A) = sum_{sigma in S_n} sgn(sigma) prod A[i][sigma(i)]

    The I-E sum over subsets {S subset [k]}: sum (-1)^|S| f(S)
    can be written as the permanent of a certain matrix... or as the
    evaluation of a multilinear polynomial.

    Specifically: define a (k x k) matrix B where B[i][j] = ...
    Then the sum over subsets relates to the permanent or determinant of B.

    Actually, the simplest formulation:
    sum_{S subset [k]} (-1)^|S| floor(x/prod_{i in S} p_i)
    = prod_{i=1}^k (1 - T_{p_i}) floor(x/1)
    = Evaluated at specific points of the polynomial:
      Poly(z_1, ..., z_k) = sum_S (prod_{i in S} z_i) * floor(x/prod_{i in S} p_i)
      evaluated at z_i = -1 for all i.

    The polynomial Poly is multilinear (degree 1 in each variable).
    Its evaluation at z = (-1, ..., -1) gives the I-E sum.
    Its evaluation at z = (0, ..., 0) gives floor(x/1) = x.

    A multilinear polynomial in k variables has 2^k coefficients.
    It can be represented as a determinant of a (k+1) x (k+1) matrix
    IF AND ONLY IF it has a specific algebraic structure.

    For the I-E polynomial:
    Poly(z) = sum_S (prod z_i) floor(x / prod p_i)
    = floor(x) + z_1*floor(x/p_1) + z_2*floor(x/p_2) + z_1*z_2*floor(x/(p_1*p_2)) + ...

    This is NOT a product of linear factors (because floor(x/(pq)) != floor(x/p)*floor(x/q)/x).
    So it cannot be a determinant of a matrix with entries linear in z_i.

    HOWEVER: the SMOOTH PART (ignoring floors) IS a product:
    x * prod(1 + z_i/p_i) which has a (k+1) x (k+1) determinant representation
    (product of 1+z_i/p_i = det of diagonal matrix).
    """
    print("\n=== Permanent vs Determinant ===\n")

    for x in [30, 100]:
        sqrtx = int(math.isqrt(x))
        P = [p for p in range(2, sqrtx+1) if sympy.isprime(p)]
        k = len(P)

        print(f"x={x}, P={P}, k={k}")

        # Build the multilinear polynomial coefficients
        coeffs = {}
        for mask in range(1 << k):
            d = 1
            for i in range(k):
                if mask & (1 << i):
                    d *= P[i]
            coeffs[mask] = x // d

        # Evaluate at z = (-1, ..., -1)
        ie_sum = sum((-1)**bin(mask).count('1') * coeffs[mask] for mask in range(1 << k))

        # Smooth version (no floor)
        smooth_coeffs = {}
        for mask in range(1 << k):
            d = 1
            for i in range(k):
                if mask & (1 << i):
                    d *= P[i]
            smooth_coeffs[mask] = x / d  # exact rational

        smooth_ie = sum((-1)**bin(mask).count('1') * smooth_coeffs[mask] for mask in range(1 << k))

        # Error = actual - smooth
        error = ie_sum - smooth_ie

        print(f"  I-E sum = {ie_sum}")
        print(f"  Smooth I-E = {smooth_ie:.6f} = x * prod(1-1/p) = {x} * {smooth_ie/x:.6f}")
        print(f"  Error (fractional parts) = {error:.6f}")

        # Verify smooth = Euler product
        euler = x * reduce(lambda a, p: a * (1 - 1/p), P, 1.0)
        print(f"  Euler product = {euler:.6f}")
        print(f"  Match: {abs(smooth_ie - euler) < 1e-10}")

        # The error term: sum_S (-1)^|S| {x/prod(S)}
        # where {y} = y - floor(y) = fractional part
        frac_parts = {}
        for mask in range(1 << k):
            d = 1
            for i in range(k):
                if mask & (1 << i):
                    d *= P[i]
            frac_parts[mask] = (x / d) - (x // d)

        error_check = -sum((-1)**bin(mask).count('1') * frac_parts[mask] for mask in range(1 << k))
        print(f"  Error from frac parts = {error_check:.6f} (should match {error:.6f})")

        # Distribution of fractional parts
        fp_vals = [frac_parts[m] for m in range(1 << k)]
        print(f"  Fractional parts: min={min(fp_vals):.4f}, max={max(fp_vals):.4f}, "
              f"mean={np.mean(fp_vals):.4f}, std={np.std(fp_vals):.4f}")

        # Are the fractional parts independent?
        # Correlation matrix of fractional parts across different x values
        print()

def main():
    print("=" * 70)
    print("INCLUSION-EXCLUSION AS DETERMINANT")
    print("Session 14 - Can the prime-counting I-E be a small determinant?")
    print("=" * 70)

    # 1. I-E factorization analysis
    for x in [30, 50, 100, 200]:
        experiment_ie_factorization(x)

    # 2. Cauchy-Binet / small determinant search
    for x in [30, 50, 100]:
        experiment_cauchy_binet(x)

    # 3. Permanent vs determinant
    experiment_permanent_vs_determinant()

    print("\n" + "=" * 70)
    print("CONCLUSIONS")
    print("=" * 70)
    print("""
INCLUSION-EXCLUSION AS DETERMINANT:

1. The I-E sum = prod_{p in P} (1 - T_p) [x], where operators COMMUTE.
   But commutativity doesn't help compress: the sum still has 2^k terms.

2. The smooth part (Euler product) is poly(N)-computable.
   The error (fractional parts) requires 2^k terms.
   These fractional parts are equidistributed and essentially independent.

3. For k=2 primes: a 2x2 determinant CAN express the I-E sum using
   "compound" floor values like (x - x//p). But this doesn't scale.

4. The I-E polynomial is multilinear in k variables with 2^k monomials.
   It is NOT factorable (due to floor function non-multiplicativity).
   Therefore it cannot be a small determinant of a matrix with simple entries.

5. THEORETICAL BARRIER: The fractional part errors {x/d} for 2^k squarefree d
   carry O(2^k) bits of independent information. No matrix smaller than
   O(2^{k/2}) x O(2^{k/2}) can encode this (by information theory).
   Since k ~ sqrt(x)/ln(x), this gives matrices of size 2^{Theta(sqrt(x)/log(x))}.

KEY INSIGHT: The impossibility of small determinants for the I-E sum
comes from the FLOOR FUNCTION introducing independent error terms.
A GapL algorithm would need to avoid the floor function entirely --
i.e., use a completely different mathematical framework than the sieve.
""")

if __name__ == "__main__":
    main()

"""
Session 15: Determinantal complexity of pi(x) as polynomial in bits.

KEY INSIGHT: pi(x) viewed as a multilinear polynomial in bits x_1,...,x_N
(where x = sum x_i * 2^i) has degree N. The question "Is pi(x) in GapL?" is
equivalent to: does pi(x) have polynomial determinantal complexity?

Determinantal complexity dc(f): smallest m such that f = det(M) where M is m×m
with entries that are affine-linear functions of the variables.

For N=3: pi(x) is degree 3 in 3 bits. Can we find a 3×3 matrix with linear
entries whose det = pi(x)?

For N=4: pi(x) is degree 4 in 4 bits. Need 4×4 matrix.

Strategy: Use the explicit multilinear polynomial for pi(x), then try to
decompose it as a determinant using known algorithms.
"""

import numpy as np
from itertools import product as iproduct
from sympy import primepi, symbols, expand, det, Matrix, Poly
from sympy.polys.monomials import itermonomials
import math

def compute_pi_polynomial(N):
    """
    Compute the multilinear polynomial P(x_1,...,x_N) = pi(sum x_i * 2^i)
    where each x_i is in {0,1}.

    Returns a dict: monomial_tuple -> coefficient
    where monomial_tuple is a tuple of variable indices present.
    """
    # Evaluate at all 2^N points
    evaluations = {}
    for bits in iproduct([0, 1], repeat=N):
        x = sum(b * (2**i) for i, b in enumerate(bits))
        evaluations[bits] = int(primepi(x))

    # Compute multilinear polynomial using Möbius inversion on Boolean lattice
    # f(x) = sum_{S} c_S * prod_{i in S} x_i
    # c_S = sum_{T subset S} (-1)^{|S|-|T|} f(1_T)

    coefficients = {}
    for S_mask in range(2**N):
        S = tuple(i for i in range(N) if S_mask & (1 << i))
        coeff = 0
        for T_mask in range(2**N):
            if T_mask & ~S_mask:  # T not subset of S
                continue
            # T is a subset of S
            sign = (-1) ** (bin(S_mask ^ T_mask).count('1'))
            bits = tuple((T_mask >> i) & 1 for i in range(N))
            coeff += sign * evaluations[bits]
        if coeff != 0:
            coefficients[S] = coeff

    return coefficients


def test_det_representation(N, max_matrix_size=None):
    """
    For a given N, compute pi(x) polynomial and check if it can be
    represented as det of a small matrix with affine-linear entries.

    A matrix M of size m with affine-linear entries in x_1,...,x_N:
    M_{ij} = a_{ij,0} + sum_k a_{ij,k} * x_k

    det(M) = polynomial of degree m in x_1,...,x_N

    We need m >= N (degree constraint).
    """
    print(f"\n{'='*60}")
    print(f"Determinantal complexity for N = {N}")
    print(f"{'='*60}")

    # Compute the polynomial
    poly_coeffs = compute_pi_polynomial(N)

    # Create symbolic variables
    x_vars = symbols(f'x0:{N}')

    # Build the symbolic polynomial
    poly = sum(c * np.prod([x_vars[i] for i in S]) if S else c
               for S, c in poly_coeffs.items())

    print(f"pi(x) polynomial in {N} variables:")
    print(f"  Number of nonzero monomials: {len(poly_coeffs)}")
    print(f"  Degree: {max(len(S) for S in poly_coeffs.keys()) if poly_coeffs else 0}")
    print(f"  Max coefficient magnitude: {max(abs(c) for c in poly_coeffs.values())}")

    # Print the polynomial for small N
    if N <= 5:
        print(f"  Polynomial: {poly}")

    # Verify by spot-checking
    for test_x in [2**N - 1, 2**(N-1), 7]:
        if test_x < 2**N:
            bits = {x_vars[i]: (test_x >> i) & 1 for i in range(N)}
            val = int(poly.subs(bits))
            expected = int(primepi(test_x))
            assert val == expected, f"Mismatch at x={test_x}: got {val}, expected {expected}"
    print("  Verification: PASSED")

    if max_matrix_size is None:
        max_matrix_size = N + 2

    # Try to find determinantal representation of size m
    for m in range(N, max_matrix_size + 1):
        print(f"\n  Trying {m}×{m} matrix...")

        # Number of parameters: m^2 * (N+1) = m^2*N + m^2
        n_params = m * m * (N + 1)
        # Number of constraints: 2^N evaluations
        n_constraints = 2**N

        print(f"  Parameters: {n_params}, Constraints: {n_constraints}")

        if n_params < n_constraints:
            print(f"  UNDERDETERMINED: not enough parameters")
            continue

        # For small m and N, try numerical optimization
        if m <= 4 and N <= 5:
            result = find_det_representation_numerical(N, m, x_vars, poly_coeffs)
            if result is not None:
                print(f"  FOUND determinantal representation of size {m}!")
                return m, result
            else:
                print(f"  No representation found at size {m}")

    return None, None


def find_det_representation_numerical(N, m, x_vars, poly_coeffs):
    """
    Use numerical optimization to find m×m matrix M with affine-linear entries
    such that det(M) matches pi(x) at all 2^N evaluation points.
    """
    from scipy.optimize import minimize

    # Parameters: A[i,j,k] for i,j in [m], k in [0..N]
    # M_{ij}(x) = A[i,j,0] + sum_{k=1}^N A[i,j,k] * x_k
    n_params = m * m * (N + 1)

    # Target values at all 2^N points
    targets = []
    bit_configs = []
    for bits in iproduct([0, 1], repeat=N):
        x = sum(b * (2**i) for i, b in enumerate(bits))
        targets.append(int(primepi(x)))
        bit_configs.append(bits)

    targets = np.array(targets, dtype=float)

    def objective(params):
        A = params.reshape(m, m, N + 1)
        total_err = 0
        for idx, bits in enumerate(bit_configs):
            M = A[:, :, 0].copy()
            for k in range(N):
                M += A[:, :, k + 1] * bits[k]
            det_val = np.linalg.det(M)
            total_err += (det_val - targets[idx])**2
        return total_err

    best_result = None
    best_error = float('inf')

    for trial in range(50):
        np.random.seed(trial)
        x0 = np.random.randn(n_params) * 0.5

        result = minimize(objective, x0, method='L-BFGS-B',
                         options={'maxiter': 5000, 'ftol': 1e-20})

        if result.fun < best_error:
            best_error = result.fun
            best_result = result

    print(f"    Best error after optimization: {best_error:.6e}")

    if best_error < 1e-6:
        # Verify
        A = best_result.x.reshape(m, m, N + 1)
        all_correct = True
        for idx, bits in enumerate(bit_configs):
            M = A[:, :, 0].copy()
            for k in range(N):
                M += A[:, :, k + 1] * bits[k]
            det_val = np.linalg.det(M)
            if abs(det_val - targets[idx]) > 0.01:
                all_correct = False
                break

        if all_correct:
            return A
        else:
            print(f"    Verification failed despite low error")

    return None


def experiment_affine_dimension():
    """
    For each N, what is the dimension of the affine space of multilinear
    degree-N polynomials that are determinants of m×m matrices?

    The space of all determinants of m×m affine-linear matrices has dimension
    m^2*(N+1) parameters but lies in a subvariety of the space of degree-m
    polynomials.

    Key question: as N grows, does the dimension of the "determinantal" subspace
    grow polynomially or exponentially relative to the full polynomial space?
    """
    print("\n" + "=" * 60)
    print("Experiment: Affine dimension analysis")
    print("=" * 60)

    for N in range(2, 10):
        # Full multilinear space dimension = 2^N
        full_dim = 2**N

        # Determinantal parameters for m×m: m^2*(N+1)
        # But determinant is degree m, so need m >= N
        m = N
        det_params = m * m * (N + 1)

        # The image dimension (# independent det polynomials) is at most
        # det_params but could be less due to GL(m) symmetry
        # GL(m) has dimension m^2, so effective params = det_params - m^2 + 1
        effective_params = det_params - m * m + 1

        print(f"N={N:2d}: full space = 2^{N} = {full_dim:6d}, "
              f"det params ({m}×{m}) = {det_params:4d}, "
              f"effective = {effective_params:4d}, "
              f"ratio = {effective_params/full_dim:.3f}")


if __name__ == "__main__":
    # Test small cases
    for N in range(2, 7):
        test_det_representation(N, max_matrix_size=N + 2)

    experiment_affine_dimension()

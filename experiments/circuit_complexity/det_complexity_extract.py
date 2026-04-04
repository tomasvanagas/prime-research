"""
Session 15: Extract and analyze determinantal representations found for pi(x).
Try harder optimization for N=5,6.
"""

import numpy as np
from itertools import product as iproduct
from sympy import primepi
from scipy.optimize import minimize, differential_evolution
import math

def compute_targets(N):
    """Compute pi(x) for all 2^N values of x."""
    targets = []
    bit_configs = []
    for bits in iproduct([0, 1], repeat=N):
        x = sum(b * (2**i) for i, b in enumerate(bits))
        targets.append(int(primepi(x)))
        bit_configs.append(bits)
    return np.array(targets, dtype=float), bit_configs


def build_matrix(A, bits, m, N):
    """Build m×m matrix from parameters A and bit configuration."""
    M = A[:, :, 0].copy()
    for k in range(N):
        M += A[:, :, k + 1] * bits[k]
    return M


def objective(params, m, N, targets, bit_configs):
    A = params.reshape(m, m, N + 1)
    total_err = 0
    for idx, bits in enumerate(bit_configs):
        M = build_matrix(A, bits, m, N)
        det_val = np.linalg.det(M)
        total_err += (det_val - targets[idx])**2
    return total_err


def find_and_extract(N, m, n_trials=200):
    """Find det representation and extract matrix structure."""
    print(f"\n{'='*60}")
    print(f"N={N}, searching for {m}×{m} determinantal representation")
    print(f"{'='*60}")

    targets, bit_configs = compute_targets(N)
    n_params = m * m * (N + 1)

    best_A = None
    best_error = float('inf')

    for trial in range(n_trials):
        np.random.seed(trial * 7 + 13)
        x0 = np.random.randn(n_params) * (1.0 if trial < n_trials//2 else 3.0)

        result = minimize(
            objective, x0, args=(m, N, targets, bit_configs),
            method='L-BFGS-B',
            options={'maxiter': 10000, 'ftol': 1e-25}
        )

        if result.fun < best_error:
            best_error = result.fun
            best_A = result.x.reshape(m, m, N + 1)

            if best_error < 1e-8:
                break

    print(f"Best error: {best_error:.2e}")

    if best_error < 1e-4:
        # Verify and round
        all_correct = True
        for idx, bits in enumerate(bit_configs):
            M = build_matrix(best_A, bits, m, N)
            det_val = np.linalg.det(M)
            if abs(det_val - targets[idx]) > 0.1:
                all_correct = False
                print(f"  FAIL at bits={bits}: det={det_val:.4f}, target={targets[idx]}")
                break

        if all_correct:
            print("FOUND! Extracting matrix structure...")

            # Print the matrix
            var_names = [f"x{i}" for i in range(N)]
            print(f"\nM = M_0 + sum_k x_k * M_k")

            print(f"\nM_0 (constant part):")
            print(np.round(best_A[:, :, 0], 3))

            for k in range(N):
                if np.max(np.abs(best_A[:, :, k+1])) > 0.001:
                    print(f"\nM_{k} (coefficient of x{k}):")
                    print(np.round(best_A[:, :, k+1], 3))

            # Check sparsity
            total_entries = m * m * (N + 1)
            nonzero = sum(np.abs(best_A[:, :, k]) > 0.01
                         for k in range(N+1)).sum()
            print(f"\nSparsity: {nonzero}/{total_entries} entries nonzero")

            # Try to find integer approximation
            # Scale to make entries close to integers
            for scale in np.arange(0.5, 10, 0.5):
                A_scaled = best_A * scale
                A_int = np.round(A_scaled)
                if np.max(np.abs(A_scaled - A_int)) < 0.1:
                    # Verify integer matrix
                    success = True
                    for idx, bits in enumerate(bit_configs):
                        M = build_matrix(A_int, bits, m, N)
                        det_val = round(np.linalg.det(M))
                        target_scaled = round(targets[idx] * scale**m)
                        if det_val != target_scaled:
                            success = False
                            break
                    if success:
                        print(f"\nInteger representation at scale {scale}:")
                        print(f"det(M) = pi(x) * {scale}^{m}")
                        for k in range(N+1):
                            if np.max(np.abs(A_int[:, :, k])) > 0.01:
                                label = f"M_0" if k == 0 else f"M_{k-1} (x{k-1})"
                                print(f"\n{label}:")
                                print(A_int[:, :, k].astype(int))
                        break

            return best_A
        else:
            print("Verification FAILED")
    else:
        print(f"Best error {best_error:.2e} too large")

    return None


def try_structured_construction(N, m):
    """
    Try specific matrix structures:
    1. Companion matrix style
    2. Tridiagonal (like elementary symmetric polynomials)
    3. Toeplitz
    4. Circulant
    """
    print(f"\n{'='*60}")
    print(f"Trying structured constructions for N={N}, m={m}")
    print(f"{'='*60}")

    targets, bit_configs = compute_targets(N)

    # Construction 1: Tridiagonal (Grenet-style for elementary symmetric)
    # For e_k, the matrix is tridiagonal with specific entries
    # Adapt: M = I + diag(coeffs) + superdiag(1s)
    # This gives det = polynomial in the diagonal entries

    # For m=N, a tridiagonal matrix has 3*N-2 + N*(N+1) parameters
    # (diagonal: N entries affine in x, super/sub-diagonal: N-1 entries affine in x)
    # That's (3*N - 2) * (N+1) parameters for general tridiagonal

    print(f"\nStructured approach: tridiagonal matrix")
    n_tri_params = (3 * m - 2) * (N + 1)

    # Parameterize: diagonal d_i = a_i0 + sum_k a_ik * x_k
    #               superdiag s_i = b_i0 + sum_k b_ik * x_k
    #               subdiag l_i = c_i0 + sum_k c_ik * x_k

    def objective_tri(params):
        d_params = params[:m * (N+1)].reshape(m, N+1)
        s_params = params[m*(N+1):m*(N+1)+(m-1)*(N+1)].reshape(m-1, N+1)
        l_params = params[m*(N+1)+(m-1)*(N+1):].reshape(m-1, N+1)

        total_err = 0
        for idx, bits in enumerate(bit_configs):
            M = np.zeros((m, m))
            for i in range(m):
                M[i, i] = d_params[i, 0] + sum(d_params[i, k+1] * bits[k] for k in range(N))
            for i in range(m-1):
                M[i, i+1] = s_params[i, 0] + sum(s_params[i, k+1] * bits[k] for k in range(N))
                M[i+1, i] = l_params[i, 0] + sum(l_params[i, k+1] * bits[k] for k in range(N))

            det_val = np.linalg.det(M)
            total_err += (det_val - targets[idx])**2
        return total_err

    best_err_tri = float('inf')
    for trial in range(100):
        np.random.seed(trial * 11 + 7)
        x0 = np.random.randn(n_tri_params) * 2.0
        result = minimize(objective_tri, x0, method='L-BFGS-B',
                         options={'maxiter': 5000, 'ftol': 1e-20})
        if result.fun < best_err_tri:
            best_err_tri = result.fun
            if best_err_tri < 1e-8:
                break

    print(f"  Best tridiagonal error: {best_err_tri:.2e}")
    if best_err_tri < 1e-4:
        print("  SUCCESS with tridiagonal structure!")
    else:
        print("  Failed — pi(x) does NOT have a small tridiagonal det representation")

    # Construction 2: Upper triangular + rank-1 perturbation
    # det(I + u*v^T + D) where D is diagonal
    print(f"\nStructured approach: I + rank-1 + diagonal")
    n_r1_params = (m + m + m) * (N + 1)  # u, v, diagonal

    def objective_r1(params):
        idx_end = m * (N+1)
        u_params = params[:idx_end].reshape(m, N+1)
        v_params = params[idx_end:2*idx_end].reshape(m, N+1)
        d_params = params[2*idx_end:3*idx_end].reshape(m, N+1)

        total_err = 0
        for idx, bits in enumerate(bit_configs):
            u = np.array([u_params[i, 0] + sum(u_params[i, k+1] * bits[k] for k in range(N)) for i in range(m)])
            v = np.array([v_params[i, 0] + sum(v_params[i, k+1] * bits[k] for k in range(N)) for i in range(m)])
            d = np.array([d_params[i, 0] + sum(d_params[i, k+1] * bits[k] for k in range(N)) for i in range(m)])

            M = np.diag(d) + np.outer(u, v)
            det_val = np.linalg.det(M)
            total_err += (det_val - targets[idx])**2
        return total_err

    best_err_r1 = float('inf')
    for trial in range(100):
        np.random.seed(trial * 13 + 3)
        x0 = np.random.randn(n_r1_params)
        result = minimize(objective_r1, x0, method='L-BFGS-B',
                         options={'maxiter': 5000, 'ftol': 1e-20})
        if result.fun < best_err_r1:
            best_err_r1 = result.fun
            if best_err_r1 < 1e-8:
                break

    print(f"  Best rank-1 error: {best_err_r1:.2e}")
    if best_err_r1 < 1e-4:
        print("  SUCCESS with rank-1 structure!")
    else:
        print("  Failed — pi(x) does NOT have a small rank-1 det representation")


if __name__ == "__main__":
    # Extract representations for N=2,3,4
    for N in [2, 3, 4]:
        find_and_extract(N, N, n_trials=200)

    # Try harder for N=5 (degree 4, so maybe 4x4 works?)
    find_and_extract(5, 4, n_trials=200)  # degree is 4, not 5
    find_and_extract(5, 5, n_trials=200)

    # Try N=6
    find_and_extract(6, 6, n_trials=200)

    # Try structured constructions
    for N in [3, 4, 5]:
        try_structured_construction(N, N)

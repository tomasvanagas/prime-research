#!/usr/bin/env python3
"""
Test whether top oscillatory singular vectors of pi(x) communication matrix
correspond to zeta zero contributions from the explicit formula.

The explicit formula says:
  pi(x) ≈ R(x) - sum_rho R(x^rho)

where rho = 1/2 + i*gamma_k. The R(x^rho) terms contain oscillations at
frequencies gamma_k in log-space:
  R(x^rho) ~ x^{1/2} * cos(gamma_k * ln(x)) / (gamma_k * ln(x))

If the SVD decomposes these contributions, we should see:
- Top osc SVs correlating with low-frequency zeta zeros (small gamma_k)
- Right singular vectors resembling cos(gamma_k * ln(b))

This would confirm the SVD is just re-discovering the explicit formula
decomposition, meaning no NEW structure exists.
"""

import numpy as np
import os

def sieve(n):
    if n < 2:
        return []
    is_prime = [True] * (n + 1)
    is_prime[0] = is_prime[1] = False
    for i in range(2, int(n**0.5) + 1):
        if is_prime[i]:
            for j in range(i*i, n + 1, i):
                is_prime[j] = False
    return is_prime

def pi_func(is_prime):
    pi = [0] * len(is_prime)
    count = 0
    for i in range(len(is_prime)):
        if is_prime[i]:
            count += 1
        pi[i] = count
    return pi

# First 50 zeta zero imaginary parts
ZETA_ZEROS = [
    14.134725, 21.022040, 25.010858, 30.424876, 32.935062,
    37.586178, 40.918719, 43.327073, 48.005151, 49.773832,
    52.970321, 56.446248, 59.347044, 60.831779, 65.112544,
    67.079811, 69.546402, 72.067158, 75.704691, 77.144840,
    79.337375, 82.910381, 84.735493, 87.425275, 88.809111,
    92.491899, 94.651344, 95.870634, 98.831194, 101.317851,
    103.725538, 105.446623, 107.168611, 111.029536, 111.874659,
    114.320220, 116.226680, 118.790783, 121.370126, 122.946829,
    124.256819, 127.516684, 129.578704, 131.087689, 133.497737,
    134.756510, 138.116043, 139.736209, 141.123707, 143.111846,
]


def build_comm_matrix(N, pi_vals):
    half = N // 2
    rows = 2 ** (N - half)
    cols = 2 ** half
    M = np.zeros((rows, cols), dtype=np.float64)
    for a in range(rows):
        for b in range(cols):
            x = a * (2 ** half) + b
            if x < len(pi_vals):
                M[a, b] = pi_vals[x]
    return M


def test_zeta_connection(N, pi_vals):
    """Test if oscillatory SVs correspond to zeta zero contributions."""
    half = N // 2
    M = build_comm_matrix(N, pi_vals)

    # SVD
    U, S, Vt = np.linalg.svd(M, full_matrices=False)

    # Remove smooth part (top 2 SVs)
    M_smooth = S[0] * np.outer(U[:, 0], Vt[0, :]) + S[1] * np.outer(U[:, 1], Vt[1, :])
    M_osc = M - M_smooth
    U_osc, S_osc, Vt_osc = np.linalg.svd(M_osc, full_matrices=False)

    # Build zeta-zero basis vectors for Bob's coordinates
    cols = 2 ** half
    b_vals = np.arange(cols)

    # For each zeta zero gamma_k, create vectors cos(gamma_k * ln(a*2^half + b))
    # and sin(gamma_k * ln(a*2^half + b)) evaluated over Bob's inputs
    # For a FIXED Alice value a, the pattern over b is:
    # cos(gamma_k * ln(a*2^half + b))
    # We need to check if the RIGHT singular vectors correlate with these

    print(f"\n{'='*60}")
    print(f"ZETA ZERO ↔ SVD CONNECTION for N={N}")
    print(f"{'='*60}")

    # For each top oscillatory SV, compute correlation with zeta-zero oscillations
    n_osc = min(20, len(S_osc[S_osc > 1e-10]))
    n_zeros = min(30, len(ZETA_ZEROS))

    print(f"\nCorrelation of top {n_osc} osc right-SVs with cos/sin(gamma_k * ln(x))")
    print(f"Using Bob's coordinates b=0..{cols-1} with Alice=1 (to avoid x=0)")

    # Evaluate basis functions at Bob's x values (with Alice=1)
    alice_val = 1
    x_vals = alice_val * (2**half) + b_vals
    ln_x = np.log(x_vals.astype(float) + 1e-10)

    zeta_cos = np.zeros((n_zeros, cols))
    zeta_sin = np.zeros((n_zeros, cols))
    for k in range(n_zeros):
        zeta_cos[k] = np.cos(ZETA_ZEROS[k] * ln_x)
        zeta_sin[k] = np.sin(ZETA_ZEROS[k] * ln_x)

    # Normalize
    for k in range(n_zeros):
        norm_c = np.linalg.norm(zeta_cos[k])
        norm_s = np.linalg.norm(zeta_sin[k])
        if norm_c > 0:
            zeta_cos[k] /= norm_c
        if norm_s > 0:
            zeta_sin[k] /= norm_s

    # For each oscillatory SV, find best matching zeta zero
    print(f"\n{'SV_idx':>6} {'SV_val':>10} {'best_zero':>10} {'gamma_k':>10} {'corr':>10} {'type':>6}")
    print("-" * 60)

    best_matches = []
    for i in range(n_osc):
        v = Vt_osc[i, :]
        v_norm = v / (np.linalg.norm(v) + 1e-15)

        best_corr = 0
        best_k = 0
        best_type = 'cos'

        for k in range(n_zeros):
            corr_c = abs(np.dot(v_norm, zeta_cos[k]))
            corr_s = abs(np.dot(v_norm, zeta_sin[k]))
            if corr_c > best_corr:
                best_corr = corr_c
                best_k = k
                best_type = 'cos'
            if corr_s > best_corr:
                best_corr = corr_s
                best_k = k
                best_type = 'sin'

        print(f"{i:>6} {S_osc[i]:>10.4f} {best_k:>10} {ZETA_ZEROS[best_k]:>10.4f} {best_corr:>10.4f} {best_type:>6}")
        best_matches.append((i, best_k, best_corr, best_type))

    # Overall: what fraction of oscillatory variance is explained by zeta-zero basis?
    # Project oscillatory matrix onto zeta-zero subspace
    # For the full matrix, we need to consider all Alice values
    rows = 2 ** (N - half)

    # Build full zeta-zero basis for the matrix
    # Each zero contributes a rank-2 component (cos and sin)
    explained_var = 0
    total_osc_var = np.sum(S_osc[S_osc > 1e-10]**2)

    # For each Alice value a, build cos/sin vectors and project
    full_basis = []
    for k in range(n_zeros):
        # Each zero gives two basis matrices (cos and sin components)
        M_cos = np.zeros_like(M_osc)
        M_sin = np.zeros_like(M_osc)
        for a in range(rows):
            x_a = a * (2**half) + b_vals
            ln_xa = np.log(x_a.astype(float) + 1e-10)
            # Weight by x^{-1/2} / (gamma_k * ln(x)) as explicit formula suggests
            weight = np.where(x_a > 1, x_a.astype(float)**(-0.5) / (ZETA_ZEROS[k] * np.maximum(ln_xa, 1e-10)), 0)
            M_cos[a, :] = weight * np.cos(ZETA_ZEROS[k] * ln_xa)
            M_sin[a, :] = weight * np.sin(ZETA_ZEROS[k] * ln_xa)

        norm_c = np.linalg.norm(M_cos)
        norm_s = np.linalg.norm(M_sin)
        if norm_c > 0:
            full_basis.append(M_cos / norm_c)
        if norm_s > 0:
            full_basis.append(M_sin / norm_s)

    if full_basis:
        # Orthogonalize and project
        # Stack as vectors
        M_osc_vec = M_osc.flatten()
        basis_vecs = np.array([b.flatten() for b in full_basis])

        # Gram-Schmidt orthogonalization
        Q = np.zeros_like(basis_vecs)
        for i in range(len(basis_vecs)):
            q = basis_vecs[i].copy()
            for j in range(i):
                q -= np.dot(q, Q[j]) * Q[j]
            norm = np.linalg.norm(q)
            if norm > 1e-10:
                Q[i] = q / norm
            else:
                Q[i] = 0

        # Project M_osc onto this subspace
        coeffs = Q @ M_osc_vec
        explained_var = np.sum(coeffs**2)

        print(f"\nZeta-zero subspace analysis ({n_zeros} zeros, {2*n_zeros} basis vectors):")
        print(f"  Total oscillatory variance: {total_osc_var:.4f}")
        print(f"  Explained by zeta basis:    {explained_var:.4f}")
        print(f"  Fraction explained:         {100*explained_var/total_osc_var:.2f}%")
        print(f"  Residual variance:          {total_osc_var - explained_var:.4f}")

    return best_matches


def main():
    max_x = 2**20 + 100
    print(f"Sieving primes up to {max_x}...")
    is_prime = sieve(max_x)
    pi_vals = pi_func(is_prime)

    for N in [10, 12, 14, 16, 18, 20]:
        if 2**N > max_x:
            break
        test_zeta_connection(N, pi_vals)


if __name__ == '__main__':
    main()

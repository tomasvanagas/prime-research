"""
Session 14: Algebraic structure analysis of Lucy DP transition matrices.

The Lucy DP computes pi(x) via a recurrence on the floor-value set V = {floor(x/k)}.
Each sieve step at prime p applies a LINEAR transformation to the state vector.
We analyze whether these matrices have exploitable algebraic structure:
1. Displacement rank (Toeplitz-like, Cauchy-like, Vandermonde-like)
2. Spectral properties (eigenvalue distribution, spectral gap)
3. Sparsity patterns
4. Kronecker/tensor structure
"""

import numpy as np
from math import isqrt
from sympy import primerange, isprime
import sys

def floor_value_set(x):
    """Return sorted floor-value set {floor(x/k) : 1 <= k <= x}."""
    vals = set()
    k = 1
    while k <= x:
        v = x // k
        vals.add(v)
        k = x // v + 1 if v > 0 else x + 1
    return sorted(vals)

def lucy_dp_matrix(x, p, V):
    """
    Construct the transition matrix M_p for sieve step at prime p.
    State: S(v) for v in V.
    Update: S(v) -= S(floor(v/p)) - (j-1) for v >= p^2.

    The affine update S_new = M_p @ S_old + b_p.
    Returns (M_p, b_p).
    """
    n = len(V)
    v_to_idx = {v: i for i, v in enumerate(V)}

    M = np.eye(n)
    b = np.zeros(n)

    for i, v in enumerate(V):
        if v >= p * p:
            vp = v // p
            if vp in v_to_idx:
                j = v_to_idx[vp]
                M[i, j] -= 1  # S(v) -= S(v//p)
                # The "+= (j-1)" part: j-1 = number of primes already sieved
                # This is a constant added to b
                # Actually, in the recurrence, it's S(v,j) = S(v,j-1) - [S(v//p,j-1) - (j-1)]
                # where j is the prime INDEX. We handle this in b.

    return M, b

def displacement_rank(M, shift_type='toeplitz'):
    """
    Compute displacement rank of M.
    Toeplitz-like: rank(Z_0 M - M Z_1) where Z is the shift matrix.
    Cauchy-like: rank(D_0 M - M D_1) where D is diagonal.
    """
    n = M.shape[0]

    if shift_type == 'toeplitz':
        Z0 = np.zeros((n, n))
        Z0[1:, :-1] = np.eye(n-1)
        Z1 = Z0.copy()
        disp = Z0 @ M - M @ Z1
    elif shift_type == 'cauchy':
        # Use floor values as diagonal entries
        d0 = np.arange(n, dtype=float)
        d1 = np.arange(n, dtype=float) + 0.5  # Shifted to avoid degeneracy
        D0 = np.diag(d0)
        D1 = np.diag(d1)
        disp = D0 @ M - M @ D1

    # Rank via SVD
    sv = np.linalg.svd(disp, compute_uv=False)
    tol = max(sv) * max(n, n) * np.finfo(float).eps * 10
    rank = np.sum(sv > tol)
    return rank, sv

def analyze_matrix(M, label):
    """Comprehensive analysis of a single transition matrix."""
    n = M.shape[0]
    print(f"\n{'='*60}")
    print(f"Matrix: {label} (size {n}x{n})")
    print(f"{'='*60}")

    # Sparsity
    nnz = np.count_nonzero(M - np.eye(n))  # Non-identity entries
    print(f"Non-identity entries: {nnz} / {n*n} ({100*nnz/n/n:.1f}%)")

    # Rank
    sv = np.linalg.svd(M, compute_uv=False)
    tol = max(sv) * n * np.finfo(float).eps * 10
    rank = np.sum(sv > tol)
    print(f"Rank: {rank} / {n}")
    print(f"Condition number: {sv[0]/sv[-1] if sv[-1] > 0 else 'inf':.2e}")

    # Eigenvalues
    eigvals = np.linalg.eigvals(M)
    unique_eigvals = len(set(np.round(eigvals.real, 6)))
    print(f"Distinct eigenvalues (approx): {unique_eigvals}")
    print(f"Eigenvalue 1 multiplicity: {np.sum(np.abs(eigvals - 1) < 1e-6)}")
    print(f"Max |eigenvalue|: {max(abs(eigvals)):.6f}")
    print(f"Min |eigenvalue|: {min(abs(eigvals)):.6f}")

    # Displacement ranks
    dr_t, sv_t = displacement_rank(M, 'toeplitz')
    dr_c, sv_c = displacement_rank(M, 'cauchy')
    print(f"Displacement rank (Toeplitz-like): {dr_t}")
    print(f"Displacement rank (Cauchy-like): {dr_c}")

    # Is M-I low rank?
    M_minus_I = M - np.eye(n)
    sv_diff = np.linalg.svd(M_minus_I, compute_uv=False)
    tol_diff = max(sv_diff) * n * np.finfo(float).eps * 10 if max(sv_diff) > 0 else 1e-10
    rank_diff = np.sum(sv_diff > tol_diff)
    print(f"Rank of M-I: {rank_diff}")

    return {
        'rank': rank, 'displacement_toeplitz': dr_t, 'displacement_cauchy': dr_c,
        'rank_M_minus_I': rank_diff, 'eigvals': eigvals
    }

def analyze_composed_matrix(matrices, labels, V):
    """Analyze the composed product of all transition matrices."""
    n = matrices[0].shape[0]
    product = np.eye(n)
    for M in matrices:
        product = M @ product

    print(f"\n{'='*60}")
    print(f"COMPOSED PRODUCT of {len(matrices)} matrices (size {n}x{n})")
    print(f"{'='*60}")

    sv = np.linalg.svd(product, compute_uv=False)
    tol = max(sv) * n * np.finfo(float).eps * 10
    rank = np.sum(sv > tol)
    print(f"Rank: {rank} / {n}")

    P_minus_I = product - np.eye(n)
    sv_diff = np.linalg.svd(P_minus_I, compute_uv=False)
    tol_diff = max(sv_diff) * n * np.finfo(float).eps * 10 if max(sv_diff) > 0 else 1e-10
    rank_diff = np.sum(sv_diff > tol_diff)
    print(f"Rank of Product-I: {rank_diff}")

    # Singular value decay
    print(f"Singular values (top 10): {sv[:10].round(4)}")
    if len(sv) > 10:
        print(f"Singular values (bottom 5): {sv[-5:].round(6)}")

    # Low-rank approximation error
    for r in [1, 2, 5, 10, min(20, n)]:
        if r < len(sv):
            err = np.sqrt(sum(sv[r:]**2)) / np.sqrt(sum(sv**2))
            print(f"Rank-{r} approx relative error: {err:.6f}")

    # Check if product has Toeplitz/circulant structure
    dr_t, _ = displacement_rank(product, 'toeplitz')
    print(f"Displacement rank (Toeplitz-like): {dr_t}")

    return product

def check_kronecker_structure(M, V, x, p):
    """
    Check if M_p has Kronecker/tensor product structure.
    If M_p = A ⊗ B for small A, B, this would enable compression.
    """
    n = M.shape[0]
    M_diff = M - np.eye(n)

    # Check rank of M-I (this is the "sieve correction")
    sv = np.linalg.svd(M_diff, compute_uv=False)
    tol = max(sv) * n * np.finfo(float).eps * 10 if max(sv) > 0 else 1e-10
    rank = np.sum(sv > tol)

    # M-I should have rank related to the number of v >= p^2
    affected = sum(1 for v in V if v >= p*p)
    print(f"  Rank(M_p - I) = {rank}, affected values = {affected}")

def main():
    for x in [100, 500, 1000, 5000]:
        print(f"\n{'#'*70}")
        print(f"# x = {x}")
        print(f"{'#'*70}")

        V = floor_value_set(x)
        n = len(V)
        primes = list(primerange(2, isqrt(x) + 1))

        print(f"Floor-value set size: {n}")
        print(f"Primes to sieve: {len(primes)} ({primes})")

        matrices = []
        for p in primes:
            M, b = lucy_dp_matrix(x, p, V)
            matrices.append(M)

        # Analyze individual matrices
        for p, M in zip(primes, matrices):
            info = analyze_matrix(M, f"M_{p} (x={x})")
            check_kronecker_structure(M, V, x, p)

        # Analyze composed product
        if len(matrices) > 0:
            product = analyze_composed_matrix(matrices, [f"M_{p}" for p in primes], V)

            # Verify: does the product give the right answer?
            # Initial state: S(v, 0) = v - 1 for each v in V
            S0 = np.array([v - 1 for v in V], dtype=float)

            # Apply affine transformations
            S = S0.copy()
            for j, (p, M) in enumerate(zip(primes, matrices)):
                # S(v) -= S(v//p) - j  for v >= p^2
                for i, v in enumerate(V):
                    if v >= p * p:
                        vp = v // p
                        v_to_idx = {vv: ii for ii, vv in enumerate(V)}
                        if vp in v_to_idx:
                            S[i] -= S[v_to_idx[vp]] - j

            pi_x = int(S[-1])  # V is sorted, so last element is x
            from sympy import primepi
            true_pi = primepi(x)
            print(f"\nComputed pi({x}) = {pi_x}, true = {true_pi}, match = {pi_x == true_pi}")

if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""
Robust tensor rank computation for chi_P using gradient descent only.
Avoids ALS numerical issues. Computes for N=6,9,12.
"""

import numpy as np
from sympy import isprime
import time
from itertools import product as iprod

np.random.seed(42)

def build_prime_tensor(N):
    k = N // 3
    d = 2 ** k
    T = np.zeros((d, d, d), dtype=np.float64)
    count = 0
    for a in range(d):
        for b in range(d):
            for c in range(d):
                x = a + b * d + c * (d * d)
                if isprime(x):
                    T[a, b, c] = 1.0
                    count += 1
    return T, count

def reconstruct(A, B, C):
    """Reconstruct tensor from CP factors A[d,r], B[d,r], C[d,r]."""
    return np.einsum('ir,jr,kr->ijk', A, B, C)

def gd_tensor_rank(T, rank, n_restarts=60, max_iter=3000, lr=0.005):
    """Find rank-r approximation via gradient descent with Adam."""
    d = T.shape[0]
    best_resid = float('inf')

    for restart in range(n_restarts):
        rng = np.random.default_rng(restart * 7 + 13)
        # Initialize with small random values
        A = rng.standard_normal((d, rank)) * 0.3
        B = rng.standard_normal((d, rank)) * 0.3
        C = rng.standard_normal((d, rank)) * 0.3

        # Adam optimizer
        mA, mB, mC = np.zeros_like(A), np.zeros_like(B), np.zeros_like(C)
        vA, vB, vC = np.zeros_like(A), np.zeros_like(B), np.zeros_like(C)
        beta1, beta2, eps = 0.9, 0.999, 1e-8

        for it in range(max_iter):
            approx = reconstruct(A, B, C)
            resid_tensor = approx - T
            loss = np.sum(resid_tensor**2)

            if loss < 1e-20:
                return np.sqrt(loss), restart

            # Gradients
            gA = np.einsum('ijk,jr,kr->ir', resid_tensor, B, C) * 2
            gB = np.einsum('ijk,ir,kr->jr', resid_tensor, A, C) * 2
            gC = np.einsum('ijk,ir,jr->kr', resid_tensor, A, B) * 2

            t = it + 1
            for param, grad, m, v in [(A, gA, mA, vA), (B, gB, mB, vB), (C, gC, mC, vC)]:
                m[:] = beta1 * m + (1 - beta1) * grad
                v[:] = beta2 * v + (1 - beta2) * grad**2
                m_hat = m / (1 - beta1**t)
                v_hat = v / (1 - beta2**t)
                param -= lr * m_hat / (np.sqrt(v_hat) + eps)

        final_resid = np.sqrt(np.sum((reconstruct(A, B, C) - T)**2))
        if final_resid < best_resid:
            best_resid = final_resid

    return best_resid, -1

def unfolding_ranks(T):
    d = T.shape[0]
    r1 = np.linalg.matrix_rank(T.reshape(d, d*d))
    r2 = np.linalg.matrix_rank(T.transpose(1,0,2).reshape(d, d*d))
    r3 = np.linalg.matrix_rank(T.transpose(2,0,1).reshape(d, d*d))
    return r1, r2, r3

def unfolding_ranks_f2(T_int):
    """Rank over F_2."""
    d = T_int.shape[0]
    results = []
    for perm in [(0,1,2), (1,0,2), (2,0,1)]:
        M = T_int.transpose(perm).reshape(d, d*d).copy()
        # Gaussian elimination over F_2
        rank = 0
        nrows, ncols = M.shape
        for col in range(ncols):
            pivot = None
            for row in range(rank, nrows):
                if M[row, col] % 2 == 1:
                    pivot = row
                    break
            if pivot is None:
                continue
            M[[rank, pivot]] = M[[pivot, rank]]
            for row in range(nrows):
                if row != rank and M[row, col] % 2 == 1:
                    M[row] = (M[row] + M[rank]) % 2
            rank += 1
        results.append(rank)
    return tuple(results)

def compute_slice_rank(T):
    """Upper bound on slice rank via greedy subtraction of rank-1 slices."""
    d = T.shape[0]
    R = T.copy()
    slices = 0

    for _ in range(d * d):  # max possible
        if np.max(np.abs(R)) < 1e-10:
            break

        # Find best mode-1 slice
        best_score = 0
        best_mode = -1
        best_idx = -1
        best_matrix = None

        for mode in range(3):
            R_unf = np.moveaxis(R, mode, 0)
            for i in range(d):
                s = R_unf[i]
                score = np.sum(s**2)
                if score > best_score:
                    best_score = score
                    best_mode = mode
                    best_idx = i
                    best_matrix = s.copy()

        if best_score < 1e-10:
            break

        # Subtract best slice
        if best_mode == 0:
            unit = np.zeros(d)
            unit[best_idx] = 1.0
            R -= np.einsum('i,jk->ijk', unit, best_matrix)
        elif best_mode == 1:
            unit = np.zeros(d)
            unit[best_idx] = 1.0
            R -= np.einsum('ik,j->ijk', best_matrix, unit)
        else:
            unit = np.zeros(d)
            unit[best_idx] = 1.0
            R -= np.einsum('ij,k->ijk', best_matrix, unit)

        slices += 1

    return slices

def main():
    print("=" * 72)
    print("ROBUST TENSOR RANK COMPUTATION FOR chi_P")
    print("=" * 72)

    for N in [6, 9, 12]:
        k = N // 3
        d = 2 ** k

        print(f"\n{'='*72}")
        print(f"N = {N}, tensor {d}x{d}x{d}, x ∈ [0, {2**N})")
        print(f"{'='*72}")

        T, n_primes = build_prime_tensor(N)
        T_int = T.astype(int)
        density = n_primes / (2**N)
        print(f"Primes: {n_primes}/{2**N} (density={density:.4f})")

        # Unfolding ranks
        ur_R = unfolding_ranks(T)
        ur_F2 = unfolding_ranks_f2(T_int)
        print(f"\nMode unfolding ranks (R):  {ur_R}")
        print(f"Mode unfolding ranks (F2): {ur_F2}")
        max_unf = max(ur_R)

        # True tensor rank via GD
        print(f"\n--- Tensor Rank Search (GD with Adam) ---")
        n_restarts = 80 if d <= 8 else (30 if d <= 16 else 15)
        max_iter = 4000 if d <= 8 else 2000

        found_rank = None
        for r in range(max(ur_R), d*d + 1):
            t0 = time.time()
            resid, restart = gd_tensor_rank(T, r, n_restarts=n_restarts, max_iter=max_iter)
            elapsed = time.time() - t0

            status = "EXACT!" if resid < 1e-6 else f"resid={resid:.6e}"
            print(f"  Rank {r:3d}: {status} ({elapsed:.1f}s)")

            if resid < 1e-6:
                found_rank = r
                break

            # Stop searching if we've gone far enough
            if r > max_unf * 3:
                print(f"  Stopping search at rank {r} (3x max unfolding)")
                found_rank = r  # upper bound
                break

        print(f"  ==> tensor_rank_R(chi_P, N={N}) {'=' if found_rank and resid < 1e-6 else '<='} {found_rank}")

        # Slice rank
        print(f"\n--- Slice Rank ---")
        sr = compute_slice_rank(T)
        print(f"  Slice rank (R, greedy): {sr}")

        # Random comparison
        print(f"\n--- Random Comparison (5 samples) ---")
        rng = np.random.default_rng(999)
        rand_ranks = []
        rand_slices = []
        for trial in range(5):
            T_rand = (rng.random((d, d, d)) < density).astype(np.float64)

            # Find rank
            for r in range(max(ur_R), d*d + 1):
                resid, _ = gd_tensor_rank(T_rand, r, n_restarts=max(n_restarts//3, 15), max_iter=max_iter)
                if resid < 1e-6:
                    rand_ranks.append(r)
                    break
                if r > max_unf * 3:
                    rand_ranks.append(r)
                    break

            sr_rand = compute_slice_rank(T_rand)
            rand_slices.append(sr_rand)
            print(f"  Random {trial}: rank<={rand_ranks[-1]}, slice={sr_rand}")

        avg_rand = np.mean(rand_ranks) if rand_ranks else 0
        print(f"\n  chi_P tensor rank: {'≤' if resid >= 1e-6 else ''}{found_rank}")
        print(f"  Random mean rank:  {avg_rand:.1f} ({rand_ranks})")
        print(f"  chi_P slice rank:  {sr}")
        print(f"  Random mean slice: {np.mean(rand_slices):.1f} ({rand_slices})")
        if avg_rand > 0:
            print(f"  Ratio (prime/random rank): {found_rank/avg_rand:.3f}")

if __name__ == "__main__":
    main()

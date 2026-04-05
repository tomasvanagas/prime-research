#!/usr/bin/env python3
"""
TRUE TENSOR RANK of the prime indicator function chi_P(x).

For N-bit numbers, partition bits into 3 equal groups:
  x = a + b * 2^{N/3} + c * 2^{2N/3}
  T[a][b][c] = 1 if x is prime, 0 otherwise.

Computes:
  1. Mode unfolding ranks (sanity check)
  2. True tensor rank over R via gradient descent (Adam) with many restarts
  3. Tensor rank over F_2 via greedy search
  4. Slice rank over R and F_2
  5. Comparison with random {0,1} tensors of same density

For N = 6, 9, 12  (tensor sizes 4x4x4, 8x8x8, 16x16x16).
"""

import numpy as np
from sympy import isprime
import time
import warnings
warnings.filterwarnings('ignore')

np.random.seed(42)

# ============================================================
# 1. BUILD TENSORS
# ============================================================

def build_prime_tensor(N):
    """Build T[a][b][c] = chi_P(a + b*2^{N/3} + c*2^{2N/3})."""
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


def build_random_tensor_same_density(shape, density, rng):
    """Build a random {0,1} tensor with approximately the given density."""
    return (rng.random(shape) < density).astype(np.float64)


# ============================================================
# 2. MODE UNFOLDING RANKS
# ============================================================

def unfolding_ranks(T):
    """Compute mode-1, mode-2, mode-3 unfolding ranks over R."""
    d1, d2, d3 = T.shape
    M1 = T.reshape(d1, -1)
    r1 = np.linalg.matrix_rank(M1)
    M2 = np.transpose(T, (1, 0, 2)).reshape(d2, -1)
    r2 = np.linalg.matrix_rank(M2)
    M3 = np.transpose(T, (2, 0, 1)).reshape(d3, -1)
    r3 = np.linalg.matrix_rank(M3)
    return r1, r2, r3


def rank_F2(M):
    """Rank of binary matrix over F_2 via Gaussian elimination."""
    M = M.copy().astype(int)
    rows, cols = M.shape
    rank = 0
    for col in range(cols):
        pivot = None
        for row in range(rank, rows):
            if M[row, col] == 1:
                pivot = row
                break
        if pivot is None:
            continue
        M[[rank, pivot]] = M[[pivot, rank]]
        for row in range(rows):
            if row != rank and M[row, col] == 1:
                M[row] = (M[row] + M[rank]) % 2
        rank += 1
    return rank


def unfolding_ranks_F2(T_int):
    """Compute mode unfolding ranks over F_2."""
    d1, d2, d3 = T_int.shape
    r1 = rank_F2(T_int.reshape(d1, -1))
    r2 = rank_F2(np.transpose(T_int, (1, 0, 2)).reshape(d2, -1))
    r3 = rank_F2(np.transpose(T_int, (2, 0, 1)).reshape(d3, -1))
    return r1, r2, r3


# ============================================================
# 3. TRUE TENSOR RANK OVER R (Adam gradient descent)
# ============================================================

def cp_decompose_adam(T, rank, max_iter=5000, lr=0.02, seed=None):
    """
    Find CP decomposition of given rank using Adam optimizer.
    Returns residual.
    """
    if seed is not None:
        rng = np.random.RandomState(seed)
    else:
        rng = np.random.RandomState()

    d1, d2, d3 = T.shape
    T_norm = np.linalg.norm(T)
    if T_norm == 0:
        return 0.0

    scale = (T_norm / rank) ** (1.0 / 3.0)
    A = rng.randn(d1, rank) * scale * 0.5
    B = rng.randn(d2, rank) * scale * 0.5
    C = rng.randn(d3, rank) * scale * 0.5

    beta1, beta2, eps = 0.9, 0.999, 1e-8
    mA = np.zeros_like(A); vA = np.zeros_like(A)
    mB = np.zeros_like(B); vB = np.zeros_like(B)
    mC = np.zeros_like(C); vC = np.zeros_like(C)

    best_resid = np.inf
    patience = 500
    no_improve = 0

    for it in range(1, max_iter + 1):
        T_hat = np.einsum('ir,jr,kr->ijk', A, B, C)
        diff = T_hat - T
        resid = np.linalg.norm(diff)

        if resid < best_resid:
            best_resid = resid
            no_improve = 0
            if resid < 1e-9:
                return resid
        else:
            no_improve += 1
            if no_improve > patience:
                break

        # Check for NaN
        if np.isnan(resid):
            return np.inf

        gA = np.einsum('ijk,jr,kr->ir', diff, B, C)
        gB = np.einsum('ijk,ir,kr->jr', diff, A, C)
        gC = np.einsum('ijk,ir,jr->kr', diff, A, B)

        # Gradient clipping
        for g in [gA, gB, gC]:
            gnorm = np.linalg.norm(g)
            if gnorm > 100:
                g *= 100.0 / gnorm

        t = it
        for param, grad, m, v in [(A, gA, mA, vA), (B, gB, mB, vB), (C, gC, mC, vC)]:
            m[:] = beta1 * m + (1 - beta1) * grad
            v[:] = beta2 * v + (1 - beta2) * grad**2
            m_hat = m / (1 - beta1**t)
            v_hat = v / (1 - beta2**t)
            param -= lr * m_hat / (np.sqrt(v_hat) + eps)

    return best_resid


def find_tensor_rank_R(T, max_rank=None, n_restarts=50, label="", time_limit=None):
    """Find minimum rank over R by trying Adam at increasing ranks."""
    d1, d2, d3 = T.shape
    if max_rank is None:
        max_rank = max(d1, d2, d3) ** 2

    r1, r2, r3 = unfolding_ranks(T)
    lower = max(r1, r2, r3)
    nnz = int(np.sum(T != 0))

    print(f"  [{label}] Shape: {T.shape}, nonzeros: {nnz}")
    print(f"  [{label}] Unfolding ranks (R): ({r1}, {r2}, {r3}), lower bound = {lower}")
    print(f"  [{label}] Searching tensor rank from {lower} to {max_rank}...")

    threshold = 1e-6
    start_time = time.time()

    for rank_try in range(lower, max_rank + 1):
        if time_limit and (time.time() - start_time) > time_limit:
            print(f"  [{label}] TIME LIMIT reached at rank {rank_try}")
            return rank_try, "TIMEOUT"

        best_resid = np.inf
        for restart in range(n_restarts):
            lr = [0.02, 0.01, 0.005][restart % 3]
            resid = cp_decompose_adam(T, rank_try, max_iter=6000, lr=lr,
                                      seed=restart * 1000 + rank_try)
            best_resid = min(best_resid, resid)
            if resid < threshold:
                break

            if time_limit and (time.time() - start_time) > time_limit:
                break

        if best_resid < threshold:
            print(f"  [{label}] Rank {rank_try}: FOUND (resid={best_resid:.2e})")
            return rank_try, "FOUND"
        else:
            print(f"  [{label}] Rank {rank_try}: best resid = {best_resid:.6e}")

    print(f"  [{label}] WARNING: reached max_rank={max_rank}")
    return max_rank, "UPPER_BOUND"


# ============================================================
# 4. TENSOR RANK OVER F_2 (Greedy)
# ============================================================

def find_tensor_rank_F2(T_int, max_rank=None):
    """Greedy upper bound on tensor rank over F_2."""
    d1, d2, d3 = T_int.shape
    if max_rank is None:
        max_rank = max(d1, d2, d3) ** 2

    r1, r2, r3 = unfolding_ranks_F2(T_int)
    lower = max(r1, r2, r3)
    print(f"  [F2] Unfolding ranks over F_2: ({r1}, {r2}, {r3}), lower bound = {lower}")

    def all_nonzero_vecs(d):
        vecs = []
        for v_int in range(1, 2**d):
            v = np.array([(v_int >> i) & 1 for i in range(d)], dtype=np.int8)
            vecs.append(v)
        return vecs

    vecs1 = all_nonzero_vecs(d1)
    vecs2 = all_nonzero_vecs(d2)
    vecs3 = all_nonzero_vecs(d3)

    residual = T_int.copy().astype(np.int8)
    rank_used = 0
    total = len(vecs1) * len(vecs2) * len(vecs3)
    exhaustive = total <= 500000

    while np.any(residual):
        best_improvement = -np.inf
        best_abc = None
        current_weight = int(np.sum(residual))

        if exhaustive:
            for a in vecs1:
                for b in vecs2:
                    ab = np.outer(a, b)
                    for c in vecs3:
                        R1 = ab[:, :, np.newaxis] * c[np.newaxis, np.newaxis, :]
                        new_w = int(np.sum((residual - R1) % 2))
                        imp = current_weight - new_w
                        if imp > best_improvement:
                            best_improvement = imp
                            best_abc = (a.copy(), b.copy(), c.copy())
        else:
            n_samples = 200000
            for _ in range(n_samples):
                a = vecs1[np.random.randint(len(vecs1))]
                b = vecs2[np.random.randint(len(vecs2))]
                c = vecs3[np.random.randint(len(vecs3))]
                ab = np.outer(a, b)
                R1 = ab[:, :, np.newaxis] * c[np.newaxis, np.newaxis, :]
                new_w = int(np.sum((residual - R1) % 2))
                imp = current_weight - new_w
                if imp > best_improvement:
                    best_improvement = imp
                    best_abc = (a.copy(), b.copy(), c.copy())

        if best_abc is None or best_improvement <= 0:
            idx = tuple(np.argwhere(residual != 0)[0])
            a = np.zeros(d1, dtype=np.int8); a[idx[0]] = 1
            b = np.zeros(d2, dtype=np.int8); b[idx[1]] = 1
            c = np.zeros(d3, dtype=np.int8); c[idx[2]] = 1
            best_abc = (a, b, c)

        a, b, c = best_abc
        R1 = np.einsum('i,j,k->ijk', a, b, c).astype(np.int8)
        residual = (residual - R1) % 2
        rank_used += 1

        remaining = int(np.sum(residual))
        if rank_used <= 3 or rank_used % 5 == 0 or remaining == 0:
            print(f"    [F2 greedy] step {rank_used}: remaining 1s = {remaining} (was {current_weight})")

        if rank_used > max_rank:
            break

    print(f"  [F2] Greedy F_2 rank: {rank_used}")
    return rank_used, lower


# ============================================================
# 5. SLICE RANK
# ============================================================

def compute_slice_rank_R(T):
    """Slice rank over R (greedy)."""
    d1, d2, d3 = T.shape
    residual = T.copy()
    n_slices = 0

    while np.linalg.norm(residual) > 1e-10:
        best_norm = 0
        best_info = None
        for mode in range(3):
            dims = [d1, d2, d3]
            for j in range(dims[mode]):
                if mode == 0: sl = residual[j, :, :]
                elif mode == 1: sl = residual[:, j, :]
                else: sl = residual[:, :, j]
                norm = np.linalg.norm(sl)
                if norm > best_norm:
                    best_norm = norm
                    best_info = (mode, j, sl.copy())
        if best_info is None or best_norm < 1e-12:
            break
        mode, j, sl = best_info
        if mode == 0: residual[j, :, :] -= sl
        elif mode == 1: residual[:, j, :] -= sl
        else: residual[:, :, j] -= sl
        n_slices += 1
    return n_slices


def compute_slice_rank_F2(T_int):
    """Slice rank over F_2 (greedy)."""
    d1, d2, d3 = T_int.shape
    residual = T_int.copy().astype(np.int8)
    n_slices = 0

    while np.any(residual):
        best_weight = 0
        best_info = None
        for mode in range(3):
            dims = [d1, d2, d3]
            for j in range(dims[mode]):
                if mode == 0: sl = residual[j, :, :]
                elif mode == 1: sl = residual[:, j, :]
                else: sl = residual[:, :, j]
                w = int(np.sum(sl != 0))
                if w > best_weight:
                    best_weight = w
                    best_info = (mode, j, sl.copy())
        if best_info is None or best_weight == 0:
            break
        mode, j, sl = best_info
        if mode == 0: residual[j, :, :] = (residual[j, :, :] - sl) % 2
        elif mode == 1: residual[:, j, :] = (residual[:, j, :] - sl) % 2
        else: residual[:, :, j] = (residual[:, :, j] - sl) % 2
        n_slices += 1
    return n_slices


# ============================================================
# MAIN
# ============================================================

def main():
    print("=" * 72)
    print("TRUE TENSOR RANK OF PRIME INDICATOR chi_P")
    print("=" * 72)
    print("Method: Adam gradient descent with multiple restarts")
    print("Threshold for 'exact': residual < 1e-6 (for {0,1} tensors)")
    print()

    results = {}

    # Configuration per N
    configs = {
        6:  {'n_restart': 80, 'n_restart_rand': 40, 'n_random': 5,
             'time_limit': None, 'time_limit_rand': 60},
        9:  {'n_restart': 50, 'n_restart_rand': 25, 'n_random': 3,
             'time_limit': 600, 'time_limit_rand': 300},
        12: {'n_restart': 20, 'n_restart_rand': 10, 'n_random': 3,
             'time_limit': 600, 'time_limit_rand': 300},
    }

    for N in [6, 9, 12]:
        cfg = configs[N]
        k = N // 3
        d = 2 ** k
        print(f"\n{'='*72}")
        print(f"N = {N} bits, tensor size {d}x{d}x{d}, x in [0, {2**N})")
        print(f"{'='*72}")

        t0 = time.time()
        T, n_primes = build_prime_tensor(N)
        T_int = T.astype(np.int8)
        density = n_primes / (2 ** N)
        print(f"\nPrimes in range: {n_primes}/{2**N} (density={density:.4f})")
        print(f"Tensor built in {time.time()-t0:.2f}s")

        # --- Unfolding ranks ---
        print(f"\n--- Mode Unfolding Ranks ---")
        ur_R = unfolding_ranks(T)
        ur_F2 = unfolding_ranks_F2(T_int)
        print(f"  Over R:   mode-1={ur_R[0]}, mode-2={ur_R[1]}, mode-3={ur_R[2]}")
        print(f"  Over F_2: mode-1={ur_F2[0]}, mode-2={ur_F2[1]}, mode-3={ur_F2[2]}")

        # --- True tensor rank over R ---
        print(f"\n--- True Tensor Rank over R ---")
        max_r = min(d * d, n_primes + 2)
        t0 = time.time()
        tr_R, status = find_tensor_rank_R(
            T, max_rank=max_r, n_restarts=cfg['n_restart'],
            label="chi_P", time_limit=cfg['time_limit']
        )
        time_R = time.time() - t0
        print(f"  >>> tensor_rank_R(chi_P, N={N}) {'<=' if status != 'FOUND' else '='} {tr_R}  [{time_R:.1f}s] ({status})")

        # --- Tensor rank over F_2 ---
        print(f"\n--- Tensor Rank over F_2 ---")
        t0 = time.time()
        tr_F2, lb_F2 = find_tensor_rank_F2(T_int, max_rank=max_r)
        time_F2 = time.time() - t0
        print(f"  >>> tensor_rank_F2(chi_P, N={N}) <= {tr_F2} (greedy UB)  [{time_F2:.1f}s]")
        print(f"  >>> F_2 lower bound (unfolding): {lb_F2}")

        # --- Slice rank ---
        print(f"\n--- Slice Rank ---")
        sr_R = compute_slice_rank_R(T)
        sr_F2 = compute_slice_rank_F2(T_int)
        print(f"  Slice rank over R:   {sr_R}")
        print(f"  Slice rank over F_2: {sr_F2}")

        # --- Random comparison ---
        n_random = cfg['n_random']
        print(f"\n--- Random {{0,1}} Tensor Comparison ({n_random} samples, density={density:.4f}) ---")
        rng = np.random.default_rng(123 + N)
        random_ranks_R = []
        random_ranks_F2_ub = []
        random_slice_R = []
        random_slice_F2 = []

        for trial in range(n_random):
            T_rand = build_random_tensor_same_density((d, d, d), density, rng)
            T_rand_int = T_rand.astype(np.int8)

            rr_R, rr_status = find_tensor_rank_R(
                T_rand, max_rank=max_r,
                n_restarts=cfg['n_restart_rand'],
                label=f"rand{trial}",
                time_limit=cfg['time_limit_rand']
            )
            random_ranks_R.append(rr_R)

            rr_F2, _ = find_tensor_rank_F2(T_rand_int, max_rank=max_r)
            random_ranks_F2_ub.append(rr_F2)

            rsr_R = compute_slice_rank_R(T_rand)
            rsr_F2 = compute_slice_rank_F2(T_rand_int)
            random_slice_R.append(rsr_R)
            random_slice_F2.append(rsr_F2)

        print(f"\n  Random tensor ranks (R):         {random_ranks_R}")
        print(f"  Random tensor ranks (F_2, UB):   {random_ranks_F2_ub}")
        print(f"  Random slice ranks (R):          {random_slice_R}")
        print(f"  Random slice ranks (F_2):        {random_slice_F2}")

        results[N] = {
            'shape': (d, d, d),
            'n_primes': n_primes,
            'density': density,
            'unfolding_R': ur_R,
            'unfolding_F2': ur_F2,
            'tensor_rank_R': tr_R,
            'tensor_rank_R_status': status,
            'tensor_rank_F2_ub': tr_F2,
            'tensor_rank_F2_lb': lb_F2,
            'slice_rank_R': sr_R,
            'slice_rank_F2': sr_F2,
            'random_ranks_R': random_ranks_R,
            'random_ranks_F2_ub': random_ranks_F2_ub,
            'random_slice_R': random_slice_R,
            'random_slice_F2': random_slice_F2,
        }

    # ============================================================
    # SUMMARY
    # ============================================================
    print(f"\n\n{'='*72}")
    print("SUMMARY TABLE")
    print(f"{'='*72}")
    hdr = (f"{'N':>3} | {'Shape':>10} | {'#P':>4} | {'Dens':>5} | "
           f"{'UnfR':>10} | {'UnfF2':>10} | "
           f"{'TR_R':>5} | {'TR_F2':>6} | {'SR_R':>5} | {'SR_F2':>6} | "
           f"{'RandR':>10} | {'RandF2':>10}")
    print(hdr)
    print("-" * len(hdr))
    for N, r in results.items():
        avg_rr = np.mean(r['random_ranks_R']) if r['random_ranks_R'] else 0
        avg_rf2 = np.mean(r['random_ranks_F2_ub']) if r['random_ranks_F2_ub'] else 0
        stat = "" if r['tensor_rank_R_status'] == 'FOUND' else "<="
        print(f"{N:>3} | {str(r['shape']):>10} | {r['n_primes']:>4} | {r['density']:>5.3f} | "
              f"{str(r['unfolding_R']):>10} | {str(r['unfolding_F2']):>10} | "
              f"{stat}{r['tensor_rank_R']:>4} | {'<=' + str(r['tensor_rank_F2_ub']):>6} | "
              f"{r['slice_rank_R']:>5} | {r['slice_rank_F2']:>6} | "
              f"{avg_rr:>10.1f} | {avg_rf2:>10.1f}")

    print(f"\n{'='*72}")
    print("KEY OBSERVATIONS")
    print(f"{'='*72}")
    for N, r in results.items():
        d = r['shape'][0]
        max_unf = max(r['unfolding_R'])
        avg_rr = np.mean(r['random_ranks_R']) if r['random_ranks_R'] else 0
        avg_rf2 = np.mean(r['random_ranks_F2_ub']) if r['random_ranks_F2_ub'] else 0
        print(f"\nN={N} (d={d}, tensor {d}x{d}x{d}):")
        print(f"  tensor_rank_R = {r['tensor_rank_R']} ({r['tensor_rank_R_status']})")
        print(f"  max_unfolding_rank = {max_unf}")
        ratio = r['tensor_rank_R'] / max_unf if max_unf > 0 else float('inf')
        print(f"  tensor_rank / unfolding_rank = {ratio:.2f}")
        if avg_rr > 0:
            print(f"  PRIME rank_R / avg RANDOM rank_R = {r['tensor_rank_R']}/{avg_rr:.1f} = {r['tensor_rank_R']/avg_rr:.2f}")
        print(f"  slice_rank_R = {r['slice_rank_R']}, slice_rank_F2 = {r['slice_rank_F2']}")
        print(f"  F_2 tensor rank: [{r['tensor_rank_F2_lb']}, {r['tensor_rank_F2_ub']}]")
        print(f"  Random comparison:")
        print(f"    R-ranks:  {r['random_ranks_R']}")
        print(f"    F2-ranks: {r['random_ranks_F2_ub']}")

    print(f"\n{'='*72}")
    print("INTERPRETATION")
    print(f"{'='*72}")
    print("""
True tensor rank measures the minimum number of rank-1 (outer product) terms
needed to exactly represent chi_P as a 3-way tensor. This is strictly harder
than the unfolding (matrix) rank.

- If tensor_rank(prime) < tensor_rank(random): primes have LESS multilinear
  complexity than random, meaning they have exploitable structure. This could
  potentially be leveraged for faster algorithms.
- If tensor_rank(prime) ~ tensor_rank(random): primes look generic at this
  scale (consistent with pseudorandomness).
- If tensor_rank(prime) > tensor_rank(random): primes require MORE multilinear
  terms, making them harder to decompose.
- Slice rank bounds partition number / multiparty communication complexity.

Note: R-rank from optimization is an UPPER BOUND (may not be the global min).
F_2 rank from greedy is also an UPPER BOUND.
""")

    print("Done.")


if __name__ == "__main__":
    main()

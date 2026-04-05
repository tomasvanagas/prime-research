#!/usr/bin/env python3
"""
k-party NOF communication complexity of pi(x) -- v2 (refined analysis).

Key insight from v1: ALL mode-unfolding ranks are FULL (= dimension of that mode).
This means the multilinear rank is trivially maximal and doesn't capture
the NOF structure well.

Better measures for k-party NOF complexity:
1. The DISCREPANCY method: for k>=3 NOF, the complexity relates to the
   largest "cylinder intersection" needed, not just unfolding rank.
2. Compute rank MODULO small primes (helps detect algebraic structure).
3. Compute the rank after subtracting the "smooth" part (R^{-1}(n)).
4. For the critical NOF measure: compute the CUT RANK for the partition
   {player i} vs {all other players}. In k-party NOF, this is the
   mode-i unfolding rank. But the NOF complexity is O(log(cut_rank))
   for the BEST player -- this IS what we computed.

Actually, re-examining the v1 results more carefully:

The key formula emerging is:
  max_rank = 2^(ceil(N/k))    when k divides N or ceil(N/k) >= 2
  max_rank = 2                 when k = N (each player has 1 bit)

But the NOF complexity lower bound is log2(max_rank) = ceil(N/k).

For k = O(log N) parties: ceil(N/k) ~ N/log(N), which still grows.
For k = N/c parties: ceil(N/k) = c, which is CONSTANT.

The critical question is whether the NOF complexity is truly Theta(N/k) or
whether it could be O(polylog(N)) for some fixed k.

This script focuses on:
1. Verifying the exact formula for max_rank vs (N,k)
2. Computing rank over GF(2) and GF(p) to detect algebraic structure
3. The normalized "information per party" = log2(max_rank) / (N/k)
4. Attempting to find if rank drops below full for the LARGEST group
5. Testing N=18 (feasible for small k)
"""

import numpy as np
from itertools import product
import time
import sys

def sieve(n):
    if n < 2:
        return [False, False]
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

def balanced_partition(N, k):
    base = N // k
    extra = N % k
    sizes = [base + 1] * extra + [base] * (k - extra)
    return sizes

def build_tensor_flat(N, k, pi_vals):
    """Build function table and reshape into k-dimensional tensor."""
    sizes = balanced_partition(N, k)
    dims = [2**s for s in sizes]
    total = 2**N
    T = np.zeros(total, dtype=np.float64)
    for x in range(min(total, len(pi_vals))):
        T[x] = pi_vals[x]
    T = T.reshape(dims)
    return T, sizes

def mode_unfolding_rank(T, mode):
    """Rank of mode-i unfolding (rows = mode i, cols = all others)."""
    k = T.ndim
    axes = [mode] + [j for j in range(k) if j != mode]
    T_perm = np.transpose(T, axes)
    nrows = T.shape[mode]
    ncols = int(np.prod([T.shape[j] for j in range(k) if j != mode]))
    M = T_perm.reshape(nrows, ncols)
    return np.linalg.matrix_rank(M, tol=1e-10)

def rank_mod_p(T, mode, p):
    """Rank of mode-i unfolding computed modulo p."""
    k = T.ndim
    axes = [mode] + [j for j in range(k) if j != mode]
    T_int = T.astype(np.int64)
    T_perm = np.transpose(T_int, axes)
    nrows = T_int.shape[mode]
    ncols = int(np.prod([T_int.shape[j] for j in range(k) if j != mode]))
    M = T_perm.reshape(nrows, ncols) % p
    # Gaussian elimination over GF(p)
    M = M.astype(np.int64) % p
    rows, cols = M.shape
    rank = 0
    used_cols = []
    for col in range(cols):
        # Find pivot
        pivot = -1
        for row in range(rank, rows):
            if M[row, col] % p != 0:
                pivot = row
                break
        if pivot == -1:
            continue
        # Swap
        M[[rank, pivot]] = M[[pivot, rank]]
        # Normalize
        inv = pow(int(M[rank, col]), p-2, p)  # Fermat's little theorem
        M[rank] = (M[rank] * inv) % p
        # Eliminate
        for row in range(rows):
            if row != rank and M[row, col] % p != 0:
                M[row] = (M[row] - M[row, col] * M[rank]) % p
        rank += 1
    return rank


def compute_nof_discrepancy(T, sizes):
    """
    Compute a discrepancy-based measure.

    For k-party NOF, the DISCREPANCY of f is:
      disc_k(f) = max over cylinder intersections C: |sum_{x in C} (-1)^{f(x)}| / 2^N

    For real-valued f (like pi(x)), we instead measure the
    approximation error: how well can f be approximated by
    a function with low NOF complexity?

    We approximate this by computing the singular value decay of the
    mode-i unfoldings and checking if the spectrum concentrates.
    """
    k = len(sizes)
    results = {}

    for mode in range(k):
        axes = [mode] + [j for j in range(k) if j != mode]
        T_perm = np.transpose(T, axes)
        nrows = T.shape[mode]
        ncols = int(np.prod([T.shape[j] for j in range(k) if j != mode]))
        M = T_perm.reshape(nrows, ncols)

        # SVD
        U, S, Vh = np.linalg.svd(M, full_matrices=False)

        # Normalized singular values
        S_norm = S / S[0] if S[0] > 0 else S

        # Effective rank at various thresholds
        eff_ranks = {}
        for thresh in [0.01, 0.001, 0.0001]:
            eff_ranks[thresh] = np.sum(S_norm > thresh)

        # Fraction of variance in top-k singular values
        total_var = np.sum(S**2)
        cumvar = np.cumsum(S**2) / total_var if total_var > 0 else np.zeros_like(S)

        results[mode] = {
            'singular_values': S[:min(10, len(S))],
            'normalized_sv': S_norm[:min(10, len(S))],
            'effective_ranks': eff_ranks,
            'cumulative_variance': cumvar[:min(10, len(cumvar))],
            'full_rank': len(S),
            'numerical_rank': np.sum(S > 1e-10),
        }

    return results


def analyze_pi_minus_smooth(N, k, pi_vals):
    """
    Subtract the smooth approximation R^{-1}(n) ~ n*log(n) and analyze
    the residual's tensor structure.

    The smooth part is easily computable; the hard part is the oscillatory residual.
    If the residual has LOW tensor rank, that's evidence pi(x) could be in lower classes.
    """
    sizes = balanced_partition(N, k)
    dims = [2**s for s in sizes]
    total = 2**N

    # Build pi(x) tensor
    T_pi = np.zeros(total, dtype=np.float64)
    for x in range(min(total, len(pi_vals))):
        T_pi[x] = pi_vals[x]

    # Build smooth approximation: pi(x) ~ x/ln(x) for x >= 2
    T_smooth = np.zeros(total, dtype=np.float64)
    for x in range(2, min(total, len(pi_vals))):
        T_smooth[x] = x / np.log(x)

    # Residual
    T_resid = T_pi - T_smooth

    T_pi = T_pi.reshape(dims)
    T_smooth = T_smooth.reshape(dims)
    T_resid = T_resid.reshape(dims)

    # Compute ranks for all three
    results = {'pi': {}, 'smooth': {}, 'residual': {}}
    for name, T in [('pi', T_pi), ('smooth', T_smooth), ('residual', T_resid)]:
        max_rank = 0
        for mode in range(k):
            r = mode_unfolding_rank(T, mode)
            results[name][mode] = r
            if r > max_rank:
                max_rank = r
        results[name]['max'] = max_rank

    return results


def verify_rank_formula(all_results):
    """
    From v1 data, the observed pattern is:
      max_rank = 2^(size of largest group)

    because the mode-unfolding along the largest group has dimensions
    (2^max_size) x (2^(N - max_size)), and the rank equals the smaller dimension.

    For balanced partitions: max group size = ceil(N/k), so max_rank = 2^ceil(N/k).

    But wait -- in v1, all modes had FULL rank. That means:
      mode-i rank = 2^(sizes[i]) for ALL i.
    So max_rank = 2^(max(sizes)) = 2^ceil(N/k).

    This is NOT meaningful for NOF lower bounds because the NOF complexity
    is O(log(rank)) = O(ceil(N/k)), and the question is whether this is
    the TRUE complexity or just an artifact.

    The mode-unfolding rank being full means every mode's fibers are
    linearly independent -- pi(x) cannot be written as a low-rank
    approximation in ANY single mode. But this doesn't rule out
    efficient NOF protocols.
    """
    print("\n" + "="*70)
    print("VERIFICATION: Is max_rank always 2^ceil(N/k)?")
    print("="*70)

    all_match = True
    for (N, k), result in sorted(all_results.items()):
        sizes = balanced_partition(N, k)
        predicted = 2**max(sizes)
        actual = result['max_rank']
        match = "OK" if actual == predicted else "MISMATCH"
        if actual != predicted:
            all_match = False
        print(f"  N={N:>2}, k={k:>2}: sizes={str(sizes):>30}, "
              f"predicted=2^{max(sizes)}={predicted:>6}, actual={actual:>6} [{match}]")

    print(f"\n  All match: {all_match}")
    return all_match


def main():
    max_N = 18
    max_x = 2**max_N + 100
    print(f"Sieving primes up to {max_x}...")
    is_prime = sieve(max_x)
    pi_vals = pi_func(is_prime)
    print(f"  pi({max_x}) = {pi_vals[-1]}")

    # ============================================================
    # PART 1: Collect all (N,k) results and verify formula
    # ============================================================
    print("\n" + "="*70)
    print("PART 1: SYSTEMATIC SCAN & FORMULA VERIFICATION")
    print("="*70)

    all_results = {}
    test_Ns = [6, 8, 10, 12, 14, 16, 18]
    test_ks = [2, 3, 4, 5, 6, 7, 8]

    for N in test_Ns:
        if 2**N > len(pi_vals):
            continue
        for k in test_ks:
            if k > N or k < 2:
                continue
            sizes = balanced_partition(N, k)
            if min(sizes) < 1:
                continue
            # Memory check
            if 2**N * 8 > 2e9:
                continue

            T, sizes = build_tensor_flat(N, k, pi_vals)
            max_rank = 0
            mode_ranks = {}
            for mode in range(k):
                r = mode_unfolding_rank(T, mode)
                mode_ranks[mode] = r
                if r > max_rank:
                    max_rank = r

            all_results[(N, k)] = {
                'N': N, 'k': k, 'sizes': sizes,
                'dims': [2**s for s in sizes],
                'mode_ranks': mode_ranks,
                'max_rank': max_rank,
                'nof_lb': np.log2(max_rank),
            }

    verify_rank_formula(all_results)

    # ============================================================
    # PART 2: Rank modulo small primes
    # ============================================================
    print("\n" + "="*70)
    print("PART 2: RANK MODULO SMALL PRIMES")
    print("="*70)
    print("  If rank_mod_p < full rank, there's algebraic structure over GF(p)")

    print(f"\n{'N':>4} {'k':>3} {'R_rank':>8} {'GF(2)':>8} {'GF(3)':>8} {'GF(5)':>8} {'GF(7)':>8}")
    print("-" * 55)

    for N in [8, 10, 12, 14]:
        for k in [2, 3, 4]:
            if (N, k) not in all_results:
                continue
            sizes = balanced_partition(N, k)
            T, _ = build_tensor_flat(N, k, pi_vals)

            # Find the mode with largest group (gives most info)
            max_mode = sizes.index(max(sizes))

            r_real = all_results[(N, k)]['mode_ranks'][max_mode]

            r_gf = {}
            for p in [2, 3, 5, 7]:
                r_gf[p] = rank_mod_p(T, max_mode, p)

            print(f"{N:>4} {k:>3} {r_real:>8} {r_gf[2]:>8} {r_gf[3]:>8} {r_gf[5]:>8} {r_gf[7]:>8}")

    # ============================================================
    # PART 3: SVD spectrum / discrepancy analysis
    # ============================================================
    print("\n" + "="*70)
    print("PART 3: SVD SPECTRUM OF MODE UNFOLDINGS")
    print("="*70)
    print("  Fast spectral decay => low effective rank => easier NOF protocols")

    for N in [12, 14, 16]:
        for k in [2, 3, 4, 5]:
            if (N, k) not in all_results:
                continue
            T, sizes = build_tensor_flat(N, k, pi_vals)
            disc = compute_nof_discrepancy(T, sizes)

            # Report for mode with largest group
            max_mode = sizes.index(max(sizes))
            info = disc[max_mode]

            print(f"\n  N={N}, k={k}, mode={max_mode} (hides {sizes[max_mode]} bits):")
            print(f"    Full rank: {info['full_rank']}, Numerical rank: {info['numerical_rank']}")
            sv_str = ', '.join(f'{s:.2f}' for s in info['normalized_sv'])
            print(f"    Top normalized SVs: [{sv_str}]")
            cv_str = ', '.join(f'{c:.4f}' for c in info['cumulative_variance'])
            print(f"    Cumulative variance: [{cv_str}]")
            print(f"    Effective ranks: {info['effective_ranks']}")

    # ============================================================
    # PART 4: Residual after smooth part
    # ============================================================
    print("\n" + "="*70)
    print("PART 4: RANK OF pi(x), x/ln(x), AND RESIDUAL")
    print("="*70)
    print("  If residual has lower rank, NOF complexity may be reducible")

    print(f"\n{'N':>4} {'k':>3} {'pi_rank':>10} {'smooth_rank':>12} {'resid_rank':>12} {'resid/pi':>10}")
    print("-" * 65)

    for N in [8, 10, 12, 14, 16]:
        for k in [2, 3, 4, 5]:
            if (N, k) not in all_results:
                continue
            res = analyze_pi_minus_smooth(N, k, pi_vals)
            pr = res['pi']['max']
            sr = res['smooth']['max']
            rr = res['residual']['max']
            ratio = rr / pr if pr > 0 else 0
            print(f"{N:>4} {k:>3} {pr:>10} {sr:>12} {rr:>12} {ratio:>10.4f}")

    # ============================================================
    # PART 5: The critical scaling question
    # ============================================================
    print("\n" + "="*70)
    print("PART 5: THE CRITICAL SCALING QUESTION")
    print("="*70)

    print("""
For k-party NOF with balanced partition:
  - Mode-i unfolding rank = 2^(sizes[i]) = FULL RANK for ALL modes
  - This gives NOF lower bound = max(sizes) = ceil(N/k)

This means:
  - For k=2: NOF >= N/2  (exponential in N)
  - For k=3: NOF >= N/3
  - For k=O(1): NOF >= N/k = Theta(N) -- still linear in N
  - For k=N/c: NOF >= c = O(1) -- constant!
  - For k=O(log N): NOF >= N/log(N) -- still superlogarithmic

The mode-unfolding rank measure CANNOT distinguish between:
  (a) pi(x) has inherently high NOF complexity (not in ACC^0)
  (b) pi(x) is in TC^0 but the mode-unfolding is too coarse to see it

The fact that ALL unfoldings are FULL RANK is actually expected for ANY
non-degenerate function (most random functions have this property).
The mode-unfolding rank being full is a NECESSARY but NOT SUFFICIENT
condition for high NOF complexity.
""")

    # ============================================================
    # PART 6: Better lower bound: rank over GF(2) structure
    # ============================================================
    print("="*70)
    print("PART 6: DETAILED GF(2) ANALYSIS -- PARITY STRUCTURE")
    print("="*70)
    print("  pi(x) mod 2 = parity of number of primes <= x")
    print("  If this has low NOF complexity, pi(x) mod 2^k could be computed")
    print("  hierarchically.\n")

    for N in [8, 10, 12, 14, 16]:
        if 2**N > len(pi_vals):
            continue
        for k in [2, 3, 4, 5, 6, 7, 8]:
            if k > N:
                continue
            sizes = balanced_partition(N, k)
            if min(sizes) < 1:
                continue

            T, _ = build_tensor_flat(N, k, pi_vals)

            # Compute rank of pi(x) mod 2 tensor
            T_mod2 = (T.astype(np.int64) % 2).astype(np.float64)

            max_rank_mod2 = 0
            for mode in range(k):
                r = mode_unfolding_rank(T_mod2, mode)
                if r > max_rank_mod2:
                    max_rank_mod2 = r

            # Also pi(x) mod 4
            T_mod4 = (T.astype(np.int64) % 4).astype(np.float64)
            max_rank_mod4 = 0
            for mode in range(k):
                r = mode_unfolding_rank(T_mod4, mode)
                if r > max_rank_mod4:
                    max_rank_mod4 = r

            expected = 2**max(sizes)
            print(f"  N={N:>2}, k={k:>2}: pi rank={expected:>5}, "
                  f"pi%2 rank={max_rank_mod2:>5}, "
                  f"pi%4 rank={max_rank_mod4:>5}, "
                  f"ratio(mod2)={max_rank_mod2/expected:.4f}")

    # ============================================================
    # PART 7: Summary for N=18 (largest feasible)
    # ============================================================
    print("\n" + "="*70)
    print("PART 7: N=18 ANALYSIS")
    print("="*70)

    N = 18
    if 2**N <= len(pi_vals):
        for k in [2, 3, 4, 5, 6, 8]:
            sizes = balanced_partition(N, k)
            if min(sizes) < 1:
                continue

            t0 = time.time()
            T, _ = build_tensor_flat(N, k, pi_vals)
            t_build = time.time() - t0

            t0 = time.time()
            max_rank = 0
            mode_ranks = {}
            for mode in range(k):
                r = mode_unfolding_rank(T, mode)
                mode_ranks[mode] = r
                if r > max_rank:
                    max_rank = r
            t_rank = time.time() - t0

            expected = 2**max(sizes)
            full = (max_rank == expected)
            print(f"  k={k:>2}: sizes={str(sizes):>25}, max_rank={max_rank:>6}, "
                  f"expected={expected:>6}, full={'YES' if full else 'NO':>3}, "
                  f"time={t_build+t_rank:.2f}s")

    # ============================================================
    # FINAL COMPREHENSIVE TABLE
    # ============================================================
    print("\n" + "="*70)
    print("FINAL COMPREHENSIVE TABLE")
    print("="*70)

    print(f"\n{'N':>3} {'k':>3} {'ceil(N/k)':>9} {'max_rank':>9} {'=2^ceil?':>8} "
          f"{'log2(r)':>8} {'NOF_lb':>7} {'for ACC0':>8}")
    print("-" * 70)

    for (N, k) in sorted(all_results.keys()):
        r = all_results[(N, k)]
        ceil_nk = max(r['sizes'])
        predicted = 2**ceil_nk
        match = "YES" if r['max_rank'] == predicted else "NO"
        nof_lb = ceil_nk  # log2(max_rank) = ceil(N/k)
        # For ACC^0, need NOF complexity O(polylog) for k = O(1)
        acc0 = "poly" if nof_lb <= 4 else "exp"
        print(f"{N:>3} {k:>3} {ceil_nk:>9} {r['max_rank']:>9} {match:>8} "
              f"{r['nof_lb']:>8.2f} {nof_lb:>7} {acc0:>8}")


if __name__ == '__main__':
    main()

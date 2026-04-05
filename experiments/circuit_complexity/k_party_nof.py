#!/usr/bin/env python3
"""
k-party NOF (Number-On-Forehead) communication complexity of pi(x).

Extends the 3-party analysis to k=2,3,4,5,6,7,8 parties.

For k-party NOF with balanced partition of N bits into k groups:
- Build k-dimensional tensor T[g1, g2, ..., gk] = pi(concat(g1,...,gk))
- Compute multilinear rank = tuple of ranks of mode-i unfoldings
- Compute maximum slice rank (rank of flattening one mode vs all others)
- The NOF complexity for player i is lower-bounded by log2(mode-i unfolding rank)

Key question: does max unfolding rank scale as 2^{N/k} or drop to poly(N)?
If it drops for k = O(log N), this would suggest pi(x) is in ACC^0 / TC^0.

Results from Session 19:
- 2-party: rank = 2^{N/2-1} + 2 for balanced partition
- 3-party: balanced cut rank ~ 2^{N/3}
"""

import numpy as np
from itertools import product
import time
import sys

def sieve(n):
    """Sieve of Eratosthenes."""
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
    """Cumulative prime counting function."""
    pi = [0] * len(is_prime)
    count = 0
    for i in range(len(is_prime)):
        if is_prime[i]:
            count += 1
        pi[i] = count
    return pi


def balanced_partition(N, k):
    """Partition N bits into k groups as evenly as possible.
    Returns list of group sizes summing to N.
    E.g., N=10, k=3 -> [4, 3, 3] or [3, 4, 3] etc.
    We use: first (N % k) groups get ceil(N/k), rest get floor(N/k).
    """
    base = N // k
    extra = N % k
    sizes = [base + 1] * extra + [base] * (k - extra)
    return sizes


def build_k_tensor(N, k, pi_vals):
    """Build k-dimensional tensor T[g1,...,gk] = pi(x) where
    x = g1 * 2^(s2+s3+...+sk) + g2 * 2^(s3+...+sk) + ... + gk.

    sizes[i] = number of bits in group i.
    """
    sizes = balanced_partition(N, k)
    dims = [2**s for s in sizes]

    total = 1
    for d in dims:
        total *= d
    assert total == 2**N, f"Dimension mismatch: product={total}, expected={2**N}"

    # Build the tensor as flat array then reshape
    T = np.zeros(total, dtype=np.float64)
    for x in range(min(total, len(pi_vals))):
        T[x] = pi_vals[x]

    T = T.reshape(dims)
    return T, sizes


def compute_multilinear_rank(T, sizes):
    """Compute the multilinear rank: for each mode i, unfold T into a matrix
    (mode-i fibers vs all other modes) and compute its rank.

    Mode-i unfolding: rows indexed by mode i (dim = 2^sizes[i]),
                      cols indexed by all other modes (dim = product of other dims).

    Returns dict mapping mode index to rank.
    """
    k = len(sizes)
    dims = T.shape
    results = {}

    for mode in range(k):
        # Move target mode to front, then flatten
        axes = [mode] + [j for j in range(k) if j != mode]
        T_perm = np.transpose(T, axes)
        nrows = dims[mode]
        ncols = int(np.prod([dims[j] for j in range(k) if j != mode]))
        M = T_perm.reshape(nrows, ncols)
        rank = np.linalg.matrix_rank(M, tol=1e-10)
        results[mode] = {
            'rank': rank,
            'nrows': nrows,
            'ncols': ncols,
            'max_possible': min(nrows, ncols),
            'bits_hidden': sizes[mode],
            'bits_visible': N - sizes[mode] if 'N' in dir() else sum(sizes) - sizes[mode],
        }

    return results


def compute_slice_ranks(T, sizes):
    """Compute slice ranks along each mode.

    Slice rank along mode i = rank of the matrix obtained by taking each
    "slice" along mode i as a vector and computing the rank of the resulting matrix.
    This is equivalent to the mode-i unfolding rank.

    Also compute: for each pair of modes (i,j), the rank of the
    matrix obtained by grouping modes {i,j} vs rest.
    """
    k = len(sizes)
    dims = T.shape

    # Single-mode unfoldings (same as multilinear rank)
    single_mode = compute_multilinear_rank(T, sizes)

    # Pairwise groupings: group two modes together vs rest
    pair_results = {}
    if k >= 3:
        for i in range(k):
            for j in range(i+1, k):
                # Group modes i,j vs all others
                group_ij = [i, j]
                group_rest = [m for m in range(k) if m not in group_ij]

                axes = group_ij + group_rest
                T_perm = np.transpose(T, axes)
                nrows = dims[i] * dims[j]
                ncols = int(np.prod([dims[m] for m in group_rest]))
                M = T_perm.reshape(nrows, ncols)
                rank = np.linalg.matrix_rank(M, tol=1e-10)
                pair_results[(i,j)] = {
                    'rank': rank,
                    'nrows': nrows,
                    'ncols': ncols,
                    'max_possible': min(nrows, ncols),
                    'bits_grouped': sizes[i] + sizes[j],
                }

    return single_mode, pair_results


def analyze_k_party(N, k, pi_vals, verbose=True):
    """Full k-party NOF analysis for N-bit inputs with k parties."""
    sizes = balanced_partition(N, k)
    dims = [2**s for s in sizes]

    # Check feasibility
    total_entries = 2**N
    tensor_size_mb = total_entries * 8 / (1024**2)  # float64
    if tensor_size_mb > 2000:  # 2GB limit
        if verbose:
            print(f"  SKIP: N={N}, k={k} would need {tensor_size_mb:.0f} MB")
        return None

    t0 = time.time()
    T, sizes = build_k_tensor(N, k, pi_vals)
    t_build = time.time() - t0

    t0 = time.time()
    ml_rank = compute_multilinear_rank(T, sizes)
    t_rank = time.time() - t0

    # Fix the bits_visible field
    total_N = sum(sizes)
    for mode in ml_rank:
        ml_rank[mode]['bits_visible'] = total_N - sizes[mode]

    max_rank = max(ml_rank[m]['rank'] for m in ml_rank)
    min_rank = min(ml_rank[m]['rank'] for m in ml_rank)

    if verbose:
        print(f"\n{'='*70}")
        print(f"{k}-PARTY NOF ANALYSIS: N={N}, partition={sizes}, dims={dims}")
        print(f"Tensor shape: {T.shape}, build={t_build:.3f}s, rank={t_rank:.3f}s")
        print(f"{'='*70}")

        for mode in range(k):
            r = ml_rank[mode]
            print(f"  Mode {mode} (player hides {r['bits_hidden']} bits, sees {r['bits_visible']}): "
                  f"rank = {r['rank']} / {r['max_possible']} "
                  f"(ratio={r['rank']/r['max_possible']:.4f})")

        print(f"\n  Max mode-unfolding rank: {max_rank}")
        print(f"  Min mode-unfolding rank: {min_rank}")
        print(f"  NOF complexity lower bound: log2(max_rank) = {np.log2(max_rank):.2f}")

        # Expected scaling
        expected = 2**(min(sizes))  # 2^{N/k} for balanced
        print(f"  Expected 2^(N/k) = 2^({min(sizes)}) = {expected}")
        print(f"  Actual max rank / expected: {max_rank / expected:.4f}")

    return {
        'N': total_N,
        'k': k,
        'sizes': sizes,
        'dims': dims,
        'ml_rank': ml_rank,
        'max_rank': max_rank,
        'min_rank': min_rank,
        'nof_lb': np.log2(max_rank),
    }


def run_pairwise_analysis(N, k, pi_vals):
    """For k>=3, also compute pairwise grouping ranks."""
    sizes = balanced_partition(N, k)
    dims = [2**s for s in sizes]
    total_entries = 2**N
    tensor_size_mb = total_entries * 8 / (1024**2)
    if tensor_size_mb > 2000:
        return None

    T, sizes = build_k_tensor(N, k, pi_vals)
    single, pairs = compute_slice_ranks(T, sizes)

    print(f"\n  Pairwise grouping ranks for k={k}, N={N}:")
    for (i,j), info in sorted(pairs.items()):
        print(f"    Modes ({i},{j}) [{info['bits_grouped']} bits] vs rest: "
              f"rank = {info['rank']} / {info['max_possible']}")

    return single, pairs


def main():
    # Precompute pi(x) for all x up to 2^18
    max_N = 18
    max_x = 2**max_N + 100
    print(f"Sieving primes up to {max_x}...")
    is_prime = sieve(max_x)
    pi_vals = pi_func(is_prime)
    print(f"  pi({max_x}) = {pi_vals[-1]}")

    # ============================================================
    # PHASE 1: Systematic scan over (N, k) pairs
    # ============================================================
    print("\n" + "="*70)
    print("PHASE 1: SYSTEMATIC k-PARTY NOF SCAN")
    print("="*70)

    # Results table: (N, k) -> result dict
    all_results = {}

    # For each N, try k = 2, 3, 4, ..., min(8, N)
    test_Ns = [6, 8, 10, 12, 14, 16]
    test_ks = [2, 3, 4, 5, 6, 7, 8]

    for N in test_Ns:
        for k in test_ks:
            if k > N:
                continue
            sizes = balanced_partition(N, k)
            if min(sizes) < 1:
                continue
            # Skip if too large
            if 2**N > len(pi_vals):
                continue

            result = analyze_k_party(N, k, pi_vals, verbose=True)
            if result is not None:
                all_results[(N, k)] = result

    # ============================================================
    # PHASE 2: Pairwise grouping analysis for k>=3
    # ============================================================
    print("\n" + "="*70)
    print("PHASE 2: PAIRWISE GROUPING ANALYSIS")
    print("="*70)

    for N in [8, 10, 12]:
        for k in [3, 4, 5, 6]:
            if k > N:
                continue
            sizes = balanced_partition(N, k)
            if min(sizes) < 1:
                continue
            if 2**N > len(pi_vals):
                continue
            run_pairwise_analysis(N, k, pi_vals)

    # ============================================================
    # PHASE 3: Scaling analysis
    # ============================================================
    print("\n" + "="*70)
    print("PHASE 3: SCALING ANALYSIS")
    print("="*70)

    # For each k, how does max_rank scale with N?
    print(f"\n{'N':>4} {'k':>3} {'partition':>20} {'max_rank':>10} {'log2(rank)':>10} "
          f"{'2^(N/k)':>10} {'ratio':>8} {'rank/rows':>10}")
    print("-" * 90)

    for N in sorted(set(r['N'] for r in all_results.values())):
        for k in sorted(set(r['k'] for r in all_results.values())):
            if (N, k) not in all_results:
                continue
            r = all_results[(N, k)]
            expected = 2**(min(r['sizes']))
            ratio = r['max_rank'] / expected
            # rank as fraction of smallest dimension
            min_dim = min(r['dims'])
            rank_frac = r['max_rank'] / min_dim
            print(f"{N:>4} {k:>3} {str(r['sizes']):>20} {r['max_rank']:>10} "
                  f"{r['nof_lb']:>10.2f} {expected:>10} {ratio:>8.4f} {rank_frac:>10.4f}")

    # ============================================================
    # PHASE 4: Key question - does complexity drop for k ~ log(N)?
    # ============================================================
    print("\n" + "="*70)
    print("PHASE 4: COMPLEXITY vs NUMBER OF PARTIES")
    print("="*70)

    for N in sorted(set(r['N'] for r in all_results.values())):
        print(f"\nN = {N}:")
        print(f"  {'k':>3} {'N/k':>5} {'max_rank':>10} {'log2(rank)':>10} {'2^(N/k)':>10} {'rank/2^(N/k)':>12}")
        print(f"  " + "-" * 60)

        for k in sorted(set(r['k'] for r in all_results.values())):
            if (N, k) not in all_results:
                continue
            r = all_results[(N, k)]
            nk = N / k
            expected = 2**(N // k)  # floor
            ratio = r['max_rank'] / expected
            print(f"  {k:>3} {nk:>5.1f} {r['max_rank']:>10} {r['nof_lb']:>10.2f} "
                  f"{expected:>10} {ratio:>12.4f}")

    # ============================================================
    # PHASE 5: Fit the scaling exponent
    # ============================================================
    print("\n" + "="*70)
    print("PHASE 5: SCALING EXPONENT FITS")
    print("="*70)

    # For fixed k, fit log2(max_rank) = alpha * N + beta
    for k in sorted(set(r['k'] for r in all_results.values())):
        Ns_k = []
        ranks_k = []
        for N in sorted(set(r['N'] for r in all_results.values())):
            if (N, k) in all_results:
                Ns_k.append(N)
                ranks_k.append(all_results[(N, k)]['max_rank'])

        if len(Ns_k) >= 3:
            log_ranks = np.log2(np.array(ranks_k, dtype=float))
            coeffs = np.polyfit(Ns_k, log_ranks, 1)
            alpha = coeffs[0]
            beta = coeffs[1]
            print(f"\n  k={k}: log2(max_rank) ~ {alpha:.4f} * N + {beta:.2f}")
            print(f"         => max_rank ~ 2^({alpha:.4f} * N)")
            print(f"         Expected for 2^(N/k): alpha = {1/k:.4f}")
            print(f"         Actual alpha / expected: {alpha * k:.4f}")
            print(f"         Data: {list(zip(Ns_k, ranks_k))}")

            if alpha * k > 1.1:
                print(f"         => Rank grows FASTER than 2^(N/k) -- strong lower bound")
            elif alpha * k > 0.9:
                print(f"         => Rank grows ~ 2^(N/k) -- consistent with high NOF complexity")
            elif alpha * k > 0.5:
                print(f"         => Rank grows slower than 2^(N/k) -- some structure exploitable")
            else:
                print(f"         => Rank grows much slower -- possible ACC^0/TC^0 evidence!")

    # ============================================================
    # PHASE 6: Critical test - k = N/2 and k = N
    # ============================================================
    print("\n" + "="*70)
    print("PHASE 6: EXTREME PARTITIONS (k near N)")
    print("="*70)

    for N in [8, 10, 12, 14, 16]:
        if 2**N > len(pi_vals):
            continue
        for k in [N//2, N-1, N]:
            if k < 2 or k > N:
                continue
            sizes = balanced_partition(N, k)
            if min(sizes) < 1:
                continue
            result = analyze_k_party(N, k, pi_vals, verbose=True)
            if result is not None:
                all_results[(N, k)] = result

    # ============================================================
    # FINAL SUMMARY
    # ============================================================
    print("\n" + "="*70)
    print("FINAL SUMMARY TABLE")
    print("="*70)

    print(f"\n{'N':>4} {'k':>3} {'partition':>25} {'max_rank':>10} {'log2(r)':>8} "
          f"{'min_dim':>8} {'rank/dim':>9} {'2^(N/k)':>8} {'r/2^(N/k)':>10}")
    print("-" * 100)

    for (N, k) in sorted(all_results.keys()):
        r = all_results[(N, k)]
        min_dim = min(r['dims'])
        expected = 2**(min(r['sizes']))
        ratio = r['max_rank'] / expected
        rank_frac = r['max_rank'] / min_dim
        print(f"{N:>4} {k:>3} {str(r['sizes']):>25} {r['max_rank']:>10} "
              f"{r['nof_lb']:>8.2f} {min_dim:>8} {rank_frac:>9.4f} "
              f"{expected:>8} {ratio:>10.4f}")

    print("\n\nKey interpretation:")
    print("  rank/dim close to 1.0 => rank is maximal (exponential NOF complexity)")
    print("  rank/dim << 1         => structure exists (lower NOF complexity)")
    print("  r/2^(N/k) ~ 1        => rank scales as expected 2^(N/k)")
    print("  r/2^(N/k) << 1       => rank drops faster (good for ACC^0/TC^0)")

    return all_results


if __name__ == '__main__':
    results = main()

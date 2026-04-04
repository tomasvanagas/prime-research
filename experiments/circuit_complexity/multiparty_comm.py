#!/usr/bin/env python3
"""
Multi-party communication complexity of pi(x).

In the Number-On-Forehead (NOF) model with k players:
- Player i sees ALL bits EXCEPT the i-th group
- Low NOF complexity implies small ACC circuits

For pi(x), we split N bits into 3 groups and build a 3D tensor.
The tensor rank gives lower bounds.

Key question: Does the 3-party complexity grow with N, or is it bounded?
If bounded (or slowly growing), this suggests TC^0 circuits might exist.

Also: test CUT RANK (bipartite rank for 3-party via cylinder intersections).
"""

import numpy as np
from itertools import product

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


def three_party_tensor(N, pi_vals):
    """Build 3D tensor T[a,b,c] = pi(a*2^{2N/3} + b*2^{N/3} + c) for balanced 3-way split."""
    if N % 3 != 0:
        # Use unbalanced split: (N-2k, k, k) where k = N//3
        k = N // 3
        k2 = N - 2*k
    else:
        k = N // 3
        k2 = k

    dim_a = 2**k2
    dim_b = 2**k
    dim_c = 2**k

    T = np.zeros((dim_a, dim_b, dim_c), dtype=np.float64)
    for a in range(dim_a):
        for b in range(dim_b):
            for c in range(dim_c):
                x = a * (2**(2*k)) + b * (2**k) + c
                if x < len(pi_vals):
                    T[a, b, c] = pi_vals[x]
    return T, (k2, k, k)


def nof_cut_ranks(T, dims):
    """Compute cut ranks for all 3 partitions in the NOF model.

    In NOF, player i's complexity is bounded by log(cut_rank_i).
    cut_rank_i = rank of the matrix obtained by unfolding T along dimension i.
    """
    k2, k, k3 = dims
    dim_a, dim_b, dim_c = T.shape

    results = {}

    # Cut rank for player A (A sees B,C on forehead):
    # Unfold: rows = a, cols = (b,c) → matrix of shape (dim_a, dim_b*dim_c)
    M_a = T.reshape(dim_a, dim_b * dim_c)
    rank_a = np.linalg.matrix_rank(M_a, tol=1e-10)
    results['A'] = {'rank': rank_a, 'max_rank': min(dim_a, dim_b * dim_c),
                     'bits': k2, 'other_bits': k + k3}

    # Cut rank for player B:
    # Unfold: rows = b, cols = (a,c) → matrix of shape (dim_b, dim_a*dim_c)
    M_b = np.transpose(T, (1, 0, 2)).reshape(dim_b, dim_a * dim_c)
    rank_b = np.linalg.matrix_rank(M_b, tol=1e-10)
    results['B'] = {'rank': rank_b, 'max_rank': min(dim_b, dim_a * dim_c),
                     'bits': k, 'other_bits': k2 + k3}

    # Cut rank for player C:
    # Unfold: rows = c, cols = (a,b) → matrix of shape (dim_c, dim_a*dim_b)
    M_c = np.transpose(T, (2, 0, 1)).reshape(dim_c, dim_a * dim_b)
    rank_c = np.linalg.matrix_rank(M_c, tol=1e-10)
    results['C'] = {'rank': rank_c, 'max_rank': min(dim_c, dim_a * dim_b),
                     'bits': k3, 'other_bits': k2 + k}

    return results


def analyze_3party(N, pi_vals):
    """Full 3-party analysis for N-bit inputs."""
    T, dims = three_party_tensor(N, pi_vals)
    k2, k, k3 = dims

    print(f"\n{'='*60}")
    print(f"3-PARTY NOF ANALYSIS for N={N}")
    print(f"Split: ({k2}, {k}, {k3}) bits → tensor shape {T.shape}")
    print(f"{'='*60}")

    # Cut ranks
    cuts = nof_cut_ranks(T, dims)
    for player in ['A', 'B', 'C']:
        r = cuts[player]
        print(f"  Player {player}: {r['bits']} bits hidden, sees {r['other_bits']} bits")
        print(f"    Cut rank = {r['rank']} / {r['max_rank']} (ratio = {r['rank']/r['max_rank']:.4f})")

    # NOF complexity lower bound
    max_cut_rank = max(cuts[p]['rank'] for p in ['A', 'B', 'C'])
    nof_lb = np.log2(max_cut_rank)
    print(f"\n  NOF complexity lower bound: {nof_lb:.2f} bits")
    print(f"  2-party balanced rank for comparison: {2**(N//2 - 1) + 2}")

    # Compare with 2-party balanced
    two_party_rank = 2**(N//2 - 1) + 2
    print(f"\n  Key ratio: 3-party max_cut_rank / 2-party_rank = {max_cut_rank / two_party_rank:.4f}")

    return cuts


def try_different_3splits(N, pi_vals):
    """Try all possible 3-way splits to find minimum cut rank."""
    print(f"\n{'='*60}")
    print(f"ALL 3-WAY SPLITS for N={N}")
    print(f"{'='*60}")

    best_max_cut = float('inf')
    best_split = None

    for k1 in range(1, N-1):
        for k2 in range(1, N-k1):
            k3 = N - k1 - k2
            if k3 < 1:
                continue
            # Only consider k1 <= k2 <= k3 to avoid duplicates
            sorted_k = tuple(sorted([k1, k2, k3]))
            if (k1, k2, k3) != sorted_k:
                continue

            dim1, dim2, dim3 = 2**k1, 2**k2, 2**k3

            # Build tensor
            T = np.zeros((dim1, dim2, dim3), dtype=np.float64)
            for a in range(dim1):
                for b in range(dim2):
                    for c in range(dim3):
                        x = a * (2**(k2+k3)) + b * (2**k3) + c
                        if x < len(pi_vals):
                            T[a, b, c] = pi_vals[x]

            # Cut ranks (unfold along each mode)
            ranks = []
            for mode in range(3):
                axes = list(range(3))
                axes.remove(mode)
                M = np.transpose(T, [mode] + axes)
                s = M.shape
                M = M.reshape(s[0], s[1] * s[2])
                r = np.linalg.matrix_rank(M, tol=1e-10)
                ranks.append(r)

            max_cut = max(ranks)
            if max_cut < best_max_cut:
                best_max_cut = max_cut
                best_split = (k1, k2, k3)

            print(f"  ({k1},{k2},{k3}): cut_ranks = ({ranks[0]},{ranks[1]},{ranks[2]}), max = {max_cut}")

    print(f"\n  BEST split: {best_split}, max cut rank = {best_max_cut}")
    print(f"  NOF lower bound: {np.log2(best_max_cut):.2f} bits")

    return best_split, best_max_cut


def scaling_analysis(results_by_N):
    """How does 3-party complexity scale with N?"""
    print(f"\n{'='*60}")
    print("SCALING ANALYSIS: 3-PARTY NOF vs 2-PARTY")
    print(f"{'='*60}")

    print(f"\n{'N':>4} {'best_3p_rank':>14} {'2p_rank':>10} {'3p/2p':>8} {'log2(3p)':>10}")
    print("-" * 55)

    for N, best_rank in sorted(results_by_N.items()):
        two_p = 2**(N//2 - 1) + 2
        ratio = best_rank / two_p
        lg = np.log2(best_rank)
        print(f"{N:>4} {best_rank:>14} {two_p:>10} {ratio:>8.4f} {lg:>10.2f}")

    Ns = sorted(results_by_N.keys())
    ranks = [results_by_N[N] for N in Ns]

    if len(Ns) >= 3:
        # Fit log2(rank) vs N
        log_ranks = np.log2(np.array(ranks, dtype=float))
        coeffs = np.polyfit(Ns, log_ranks, 1)
        print(f"\n  Fit: log2(best_3p_rank) ≈ {coeffs[1]:.2f} + {coeffs[0]:.4f} * N")
        print(f"  => rank ≈ 2^({coeffs[0]:.4f} * N)")
        print(f"  Exponent: {coeffs[0]:.4f} (compare: 2-party = 0.5)")

        if coeffs[0] < 0.4:
            print(f"  => 3-party complexity grows SLOWER than 2-party!")
        elif coeffs[0] > 0.45:
            print(f"  => 3-party complexity is similar to 2-party (exponential)")
        else:
            print(f"  => Borderline; need more data points")


def main():
    max_x = 2**18 + 100  # Limit for 3-party (tensor gets big)
    print(f"Sieving primes up to {max_x}...")
    is_prime = sieve(max_x)
    pi_vals = pi_func(is_prime)

    results_by_N = {}

    # Balanced 3-way split for small N
    for N in [6, 9, 12, 15, 18]:
        if 2**N > max_x:
            break
        analyze_3party(N, pi_vals)

    # Try all 3-way splits for small N
    for N in [6, 8, 9, 10, 12, 14, 15, 16, 18]:
        if 2**N > max_x:
            break
        # Skip if tensor would be too large
        if 2**N > 2**16:
            print(f"\nSkipping N={N} for exhaustive split search (too large)")
            continue
        best_split, best_rank = try_different_3splits(N, pi_vals)
        results_by_N[N] = best_rank

    scaling_analysis(results_by_N)


if __name__ == '__main__':
    main()

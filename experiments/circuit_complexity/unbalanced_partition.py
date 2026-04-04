#!/usr/bin/env python3
"""
Unbalanced communication complexity partitions for pi(x).

For N-bit numbers, Alice has the top k bits, Bob has the bottom (N-k) bits.
We form the communication matrix M[a][b] = pi(a * 2^(N-k) + b)
and compute its rank over the rationals (via numpy float64).

Session 17 found: for balanced partition (k = N/2), rank = 2^{N/2-1} + 2.
Question: Does rank depend on which partition we use? Is there a "good" split?
"""

import numpy as np
import time
import sys

def sieve_pi(max_val):
    """Compute pi(x) for all x in [0, max_val] using sieve of Eratosthenes."""
    if max_val < 2:
        return np.zeros(max_val + 1, dtype=np.int64)
    is_prime = np.ones(max_val + 1, dtype=bool)
    is_prime[0] = is_prime[1] = False
    for i in range(2, int(max_val**0.5) + 1):
        if is_prime[i]:
            is_prime[i*i::i] = False
    return np.cumsum(is_prime).astype(np.int64)


def compute_ranks_for_N(N):
    """Compute communication matrix rank for all partitions k=1..N-1."""
    max_val = (1 << N) - 1
    pi = sieve_pi(max_val)

    results = {}
    for k in range(1, N):
        n_bob = N - k
        n_alice_vals = 1 << k
        n_bob_vals = 1 << n_bob

        # Build communication matrix M[a][b] = pi(a * 2^(N-k) + b)
        # a in [0, 2^k - 1], b in [0, 2^(N-k) - 1]
        M = np.zeros((n_alice_vals, n_bob_vals), dtype=np.float64)
        for a in range(n_alice_vals):
            base = a * n_bob_vals
            M[a, :] = pi[base:base + n_bob_vals]

        r = np.linalg.matrix_rank(M)
        results[k] = {
            'rank': r,
            'alice_bits': k,
            'bob_bits': n_bob,
            'matrix_shape': (n_alice_vals, n_bob_vals),
            'max_possible_rank': min(n_alice_vals, n_bob_vals),
        }

    return results


def main():
    N_values = [8, 10, 12, 14, 16]

    all_results = {}

    for N in N_values:
        t0 = time.time()
        print(f"Computing N={N} (2^N = {1<<N})...", flush=True)
        results = compute_ranks_for_N(N)
        elapsed = time.time() - t0
        all_results[N] = results
        print(f"  Done in {elapsed:.2f}s", flush=True)

    # Print results
    print("\n" + "="*80)
    print("COMMUNICATION MATRIX RANK FOR pi(x) — UNBALANCED PARTITIONS")
    print("="*80)

    for N in N_values:
        results = all_results[N]
        balanced_k = N // 2
        balanced_rank = results[balanced_k]['rank']
        formula_rank = (1 << (N // 2 - 1)) + 2  # 2^{N/2-1} + 2

        print(f"\nN = {N}  (x in [0, {(1<<N)-1}], pi(max) = {sieve_pi((1<<N)-1)[-1]})")
        print(f"  Balanced partition (k={balanced_k}): rank = {balanced_rank}, "
              f"formula 2^{{N/2-1}}+2 = {formula_rank}")
        print(f"  {'k':>3s}  {'Alice':>6s}  {'Bob':>6s}  {'Matrix':>12s}  "
              f"{'Rank':>6s}  {'MaxRank':>8s}  {'Rank/Max':>8s}  {'Rank/Balanced':>14s}")
        print(f"  {'-'*3}  {'-'*6}  {'-'*6}  {'-'*12}  "
              f"{'-'*6}  {'-'*8}  {'-'*8}  {'-'*14}")

        for k in sorted(results.keys()):
            r = results[k]
            rank = r['rank']
            max_rank = r['max_possible_rank']
            ratio_max = rank / max_rank if max_rank > 0 else 0
            ratio_balanced = rank / balanced_rank if balanced_rank > 0 else 0
            print(f"  {k:3d}  {r['alice_bits']:6d}  {r['bob_bits']:6d}  "
                  f"{str(r['matrix_shape']):>12s}  {rank:6d}  {max_rank:8d}  "
                  f"{ratio_max:8.4f}  {ratio_balanced:14.4f}")

    # Summary: minimum rank across all partitions for each N
    print("\n" + "="*80)
    print("SUMMARY: MINIMUM RANK ACROSS ALL PARTITIONS")
    print("="*80)
    print(f"{'N':>4s}  {'MinRank':>8s}  {'AtK':>4s}  {'BalancedRank':>13s}  "
          f"{'Formula':>8s}  {'Min/Balanced':>13s}")
    print(f"{'-'*4}  {'-'*8}  {'-'*4}  {'-'*13}  {'-'*8}  {'-'*13}")

    for N in N_values:
        results = all_results[N]
        balanced_k = N // 2
        balanced_rank = results[balanced_k]['rank']
        formula_rank = (1 << (N // 2 - 1)) + 2

        min_rank = min(r['rank'] for r in results.values())
        min_k = min(k for k, r in results.items() if r['rank'] == min_rank)

        ratio = min_rank / balanced_rank if balanced_rank > 0 else 0
        print(f"{N:4d}  {min_rank:8d}  {min_k:4d}  {balanced_rank:13d}  "
              f"{formula_rank:8d}  {ratio:13.4f}")

    # Check symmetry: rank(k) vs rank(N-k)
    print("\n" + "="*80)
    print("SYMMETRY CHECK: rank(k) vs rank(N-k)")
    print("="*80)
    for N in N_values:
        results = all_results[N]
        print(f"\nN = {N}:")
        for k in range(1, N // 2 + 1):
            r1 = results[k]['rank']
            r2 = results[N - k]['rank']
            sym = "==" if r1 == r2 else "!="
            print(f"  rank(k={k}) = {r1},  rank(k={N-k}) = {r2}  {sym}")

    # Key analysis: rank as function of min(k, N-k)
    print("\n" + "="*80)
    print("RANK vs min(k, N-k) — LOOKING FOR PATTERN")
    print("="*80)
    for N in N_values:
        results = all_results[N]
        print(f"\nN = {N}:")
        print(f"  {'min(k,N-k)':>10s}  {'Rank':>8s}  {'2^min-1+2':>10s}  {'Match?':>7s}")
        seen = set()
        for k in range(1, N):
            m = min(k, N - k)
            if m in seen:
                continue
            seen.add(m)
            rank = results[k]['rank']
            predicted = (1 << (m - 1)) + 2 if m >= 1 else 1
            match = "YES" if rank == predicted else "NO"
            print(f"  {m:10d}  {rank:8d}  {predicted:10d}  {match:>7s}")


if __name__ == '__main__':
    main()

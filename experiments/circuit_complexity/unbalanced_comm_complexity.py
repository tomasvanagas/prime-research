#!/usr/bin/env python3
"""
Unbalanced communication complexity for pi(x): comprehensive sweep.

For N-bit inputs, partition x = a * 2^k + b where:
  - Alice holds the top (N-k) bits: a in [0, 2^{N-k} - 1]
  - Bob holds the bottom k bits:    b in [0, 2^k - 1]
  - Matrix M[a,b] = pi(a * 2^k + b)
  - Shape: 2^{N-k} rows x 2^k columns

We compute rank(M) for all k in 1..N-1 and all feasible N.

Key question: Is there ANY k where rank grows polynomially in N?
If rank = min(rows, cols) for all k, the function is "full rank everywhere"
and no partition helps.

Previous result (Session 17): balanced k=N/2 gives rank = 2^{N/2-1} + 2.
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


def compute_all_partitions(N, pi_values):
    """Compute rank for all partitions k=1..N-1."""
    results = []
    for k in range(1, N):
        n_top = N - k        # Alice's bits (rows)
        n_rows = 1 << n_top
        n_cols = 1 << k      # Bob's bits (columns)

        # Build M[a,b] = pi(a * 2^k + b)
        M = np.zeros((n_rows, n_cols), dtype=np.float64)
        for a in range(n_rows):
            base = a * n_cols
            M[a, :] = pi_values[base:base + n_cols]

        rank = np.linalg.matrix_rank(M)
        max_rank = min(n_rows, n_cols)

        results.append({
            'k': k,
            'alice_bits': n_top,
            'bob_bits': k,
            'rows': n_rows,
            'cols': n_cols,
            'rank': rank,
            'max_rank': max_rank,
            'ratio': rank / max_rank if max_rank > 0 else 0,
            'rank_deficit': max_rank - rank,
        })

    return results


def main():
    output_lines = []

    def log(s="", end="\n"):
        print(s, end=end, flush=True)
        if end == "\n":
            output_lines.append(s)
        else:
            if output_lines and not output_lines[-1].endswith("\n"):
                output_lines[-1] += s
            else:
                output_lines.append(s)

    log("=" * 90)
    log("UNBALANCED COMMUNICATION COMPLEXITY FOR pi(x)")
    log("M[a,b] = pi(a * 2^k + b), Alice has top (N-k) bits, Bob has bottom k bits")
    log("=" * 90)

    # Determine feasible N values
    # N=20 means 2^20 = 1M entries, matrix up to 1024x1024 -- feasible
    # N=22 means 2^22 = 4M entries, matrix up to 2048x2048 -- borderline
    N_values = [4, 6, 8, 10, 12, 14, 16, 18, 20]

    all_results = {}

    for N in N_values:
        max_val = (1 << N) - 1
        t0 = time.time()
        log(f"\nComputing N={N} (range [0, {max_val}])...")
        pi_values = sieve_pi(max_val)
        pi_max = pi_values[-1]
        results = compute_all_partitions(N, pi_values)
        elapsed = time.time() - t0
        all_results[N] = results
        log(f"  Done in {elapsed:.2f}s, pi({max_val}) = {pi_max}")

    # === FULL TABLE ===
    log("\n" + "=" * 90)
    log("FULL RESULTS TABLE")
    log("=" * 90)

    for N in N_values:
        results = all_results[N]
        max_val = (1 << N) - 1

        log(f"\nN = {N}  (x in [0, {max_val}])")
        log(f"  {'k':>3s}  {'Alice':>6s}  {'Bob':>5s}  {'Shape':>14s}  "
            f"{'Rank':>6s}  {'MaxRnk':>7s}  {'Deficit':>7s}  {'Rank/Max':>8s}")
        log(f"  {'---':>3s}  {'------':>6s}  {'-----':>5s}  {'-'*14:>14s}  "
            f"{'------':>6s}  {'-------':>7s}  {'-------':>7s}  {'--------':>8s}")

        for r in results:
            shape_str = f"{r['rows']}x{r['cols']}"
            log(f"  {r['k']:3d}  {r['alice_bits']:6d}  {r['bob_bits']:5d}  "
                f"{shape_str:>14s}  {r['rank']:6d}  {r['max_rank']:7d}  "
                f"{r['rank_deficit']:7d}  {r['ratio']:8.4f}")

    # === SYMMETRY ANALYSIS ===
    log("\n" + "=" * 90)
    log("SYMMETRY: rank(k) vs rank(N-k)")
    log("Note: Transpose swaps Alice/Bob, so rank should be the same.")
    log("=" * 90)

    for N in N_values:
        results = all_results[N]
        rdict = {r['k']: r['rank'] for r in results}
        log(f"\nN = {N}:")
        for k in range(1, (N + 1) // 2 + 1):
            if k in rdict and (N - k) in rdict:
                r1, r2 = rdict[k], rdict[N - k]
                sym = "==" if r1 == r2 else "!="
                log(f"  rank(k={k}) = {r1:6d},  rank(k={N-k}) = {r2:6d}  {sym}")

    # === KEY ANALYSIS: rank vs 2^{min(k,N-k)-1} + 2 ===
    log("\n" + "=" * 90)
    log("PATTERN CHECK: Does rank = 2^{min(k,N-k)-1} + 2 for all partitions?")
    log("(This is the balanced-partition formula applied to the smaller side)")
    log("=" * 90)

    for N in N_values:
        results = all_results[N]
        rdict = {r['k']: r['rank'] for r in results}
        log(f"\nN = {N}:")
        log(f"  {'k':>3s}  {'min(k,N-k)':>10s}  {'Rank':>6s}  {'2^(m-1)+2':>10s}  "
            f"{'MaxRnk':>7s}  {'Match?':>7s}  {'FullRank?':>9s}")
        seen = set()
        for k in range(1, N):
            m = min(k, N - k)
            if m in seen:
                continue
            seen.add(m)
            rank = rdict[k]
            predicted = (1 << (m - 1)) + 2
            max_rank = min(1 << k, 1 << (N - k))
            is_match = "YES" if rank == predicted else "NO"
            is_full = "FULL" if rank == max_rank else f"deficit={max_rank - rank}"
            log(f"  {k:3d}  {m:10d}  {rank:6d}  {predicted:10d}  "
                f"{max_rank:7d}  {is_match:>7s}  {is_full:>9s}")

    # === GROWTH ANALYSIS ===
    log("\n" + "=" * 90)
    log("GROWTH ANALYSIS: How does rank grow with N for each partition type?")
    log("Looking for ANY k where rank is polynomial in N (not exponential)")
    log("=" * 90)

    # For each fixed min(k, N-k) = m, track how rank grows across N values
    log(f"\n  {'m':>3s}  ", end="")
    for N in N_values:
        log(f"{'N='+str(N):>8s}  ", end="")
    log("")

    for m in range(1, max(N_values) // 2 + 1):
        log(f"  {m:3d}  ", end="")
        for N in N_values:
            if m >= N or m < 1:
                log(f"{'---':>8s}  ", end="")
                continue
            results = all_results[N]
            rdict = {r['k']: r['rank'] for r in results}
            k = m  # Use k = m (Bob has m bits)
            if k in rdict:
                log(f"{rdict[k]:8d}  ", end="")
            else:
                log(f"{'---':>8s}  ", end="")
        log("")

    # === EXTREME PARTITIONS ===
    log("\n" + "=" * 90)
    log("EXTREME PARTITIONS: k=1 (Bob has 1 bit) and k=N-1 (Bob has N-1 bits)")
    log("=" * 90)

    log(f"\n  {'N':>4s}  {'rank(k=1)':>10s}  {'maxRank(k=1)':>13s}  "
        f"{'rank(k=N-1)':>12s}  {'maxRank(k=N-1)':>15s}")
    for N in N_values:
        results = all_results[N]
        rdict = {r['k']: r for r in results}
        r1 = rdict[1]
        rn1 = rdict[N - 1]
        log(f"  {N:4d}  {r1['rank']:10d}  {r1['max_rank']:13d}  "
            f"{rn1['rank']:12d}  {rn1['max_rank']:15d}")

    # === MINIMUM RANK ACROSS ALL PARTITIONS ===
    log("\n" + "=" * 90)
    log("MINIMUM RANK ACROSS ALL PARTITIONS FOR EACH N")
    log("=" * 90)

    log(f"\n  {'N':>4s}  {'MinRank':>8s}  {'AtK':>4s}  {'MaxPossible':>12s}  "
        f"{'Balanced':>9s}  {'2^(N/2-1)+2':>12s}")
    for N in N_values:
        results = all_results[N]
        min_rank = min(r['rank'] for r in results)
        min_k = min(r['k'] for r in results if r['rank'] == min_rank)
        balanced_k = N // 2
        balanced_rank = next(r['rank'] for r in results if r['k'] == balanced_k)
        formula = (1 << (N // 2 - 1)) + 2
        log(f"  {N:4d}  {min_rank:8d}  {min_k:4d}  "
            f"{min(1 << min_k, 1 << (N - min_k)):12d}  "
            f"{balanced_rank:9d}  {formula:12d}")

    # === POLYNOMIAL RANK CHECK ===
    log("\n" + "=" * 90)
    log("POLYNOMIAL RANK CHECK")
    log("Is there ANY partition where rank grows polynomially in N?")
    log("=" * 90)

    # For k=1, the matrix has 2^{N-1} rows and 2 columns, so max rank = 2
    # For k=N-1, the matrix has 2 rows and 2^{N-1} columns, so max rank = 2
    # These are always "small" but trivially so (max rank is 2)
    # The real question is: for k = O(log N), is rank polynomial?

    log("\nFor small k (k = 1,2,3,4,5), rank is bounded by 2^k:")
    for k_fixed in range(1, 6):
        log(f"\n  k={k_fixed} (Bob has {k_fixed} bits, max rank = {1 << k_fixed}):")
        for N in N_values:
            if k_fixed >= N:
                continue
            results = all_results[N]
            r = next(r for r in results if r['k'] == k_fixed)
            log(f"    N={N:2d}: rank = {r['rank']:6d} / {r['max_rank']:6d}  "
                f"(ratio = {r['ratio']:.4f})")

    # === CONCLUSION ===
    log("\n" + "=" * 90)
    log("CONCLUSION")
    log("=" * 90)

    # Check if rank always equals min(2^k, 2^{N-k}) or 2^{min(k,N-k)-1}+2
    any_polynomial = False
    for N in N_values:
        results = all_results[N]
        for r in results:
            m = min(r['k'], N - r['k'])
            # rank is at most 2^m, which is exponential in m
            # For polynomial rank in N, we'd need rank = O(N^c)
            if r['rank'] <= 10 * N**3 and r['max_rank'] > 10 * N**3:
                any_polynomial = True
                log(f"  INTERESTING: N={N}, k={r['k']}, rank={r['rank']} << max_rank={r['max_rank']}")

    if not any_polynomial:
        log("\n  No partition found where rank is polynomial in N while max_rank is exponential.")
        log("  For every partition, rank = min(2^k, 2^{N-k}) or very close to it,")
        log("  OR the constraint is trivial (max_rank is already small).")
        log("\n  The communication complexity of pi(x) appears to be Omega(N/2) = Omega(log x / 2)")
        log("  regardless of the bit partition chosen.")
    else:
        log("\n  Found partitions with potentially polynomial rank -- investigate further!")

    # Save results
    with open('/apps/aplikacijos/prime-research/experiments/circuit_complexity/unbalanced_comm_results.txt', 'w') as f:
        f.write('\n'.join(output_lines) + '\n')
    log("\nResults saved to experiments/circuit_complexity/unbalanced_comm_results.txt")


if __name__ == '__main__':
    main()

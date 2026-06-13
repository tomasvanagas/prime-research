"""
Thread 11 / Slot 1 — bounded-direction line cover (numpy).

For large N (e.g. 10^6), the all-pairs O(pi(N)^2) line enumeration is
infeasible in pure Python. This script restricts the line search to
canonical directions (a, b) with gcd(|a|,|b|)=1 and max(|a|,|b|) <= K_max.
Each such direction defines a partition of the primes into lines (by
intercept). We then greedily cover.

This gives an UPPER BOUND on the true greedy line cover (since we
allow fewer lines). Empirically the dominant lines on the Ulam spiral
are the slope-1 diagonals — captured at K_max = 1 — so K_max = 30 is
generous.

CLI:
  python p11_ulam_bounded.py --N 1000000 --K 30 --baseline-trials 3
"""
from __future__ import annotations

import argparse
import csv
import math
import os
import random
import sys
import time
from collections import defaultdict
from math import gcd

import numpy as np


def sieve_primes(N: int) -> np.ndarray:
    if N < 2:
        return np.array([], dtype=np.int64)
    sieve = np.ones(N + 1, dtype=bool)
    sieve[0] = sieve[1] = False
    for i in range(2, int(N ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i :: i] = False
    return np.nonzero(sieve)[0].astype(np.int64)


def ulam_coords_array(N: int) -> tuple[np.ndarray, np.ndarray]:
    """Compute Ulam coords for n=1..N, return arrays X, Y of size N+1.
    X[0], Y[0] are unused (set to 0).
    """
    X = np.zeros(N + 1, dtype=np.int32)
    Y = np.zeros(N + 1, dtype=np.int32)
    if N < 1:
        return X, Y
    # n=1 at (0, 0). Already initialised.
    if N == 1:
        return X, Y
    directions = [(1, 0), (0, 1), (-1, 0), (0, -1)]
    dir_idx = 0
    leg_length = 1
    n = 1
    x, y = 0, 0
    while n < N:
        for _ in range(2):
            dx, dy = directions[dir_idx]
            for _ in range(leg_length):
                x += dx
                y += dy
                n += 1
                X[n] = x
                Y[n] = y
                if n >= N:
                    return X, Y
            dir_idx = (dir_idx + 1) % 4
        leg_length += 1
    return X, Y


def canonical_directions(K: int) -> list[tuple[int, int]]:
    """All (a, b) with gcd(|a|,|b|)=1 and max(|a|,|b|) <= K, canonical
    (positive a, or a=0 and positive b)."""
    dirs = [(0, 1)]  # vertical line
    for a in range(1, K + 1):
        for b in range(-K, K + 1):
            if gcd(a, abs(b)) == 1:
                dirs.append((a, b))
    return dirs


def bounded_greedy_line_cover(prime_X: np.ndarray, prime_Y: np.ndarray, K: int, verbose: bool = False) -> tuple[int, list]:
    """Bounded-direction greedy line cover using numpy.

    For each direction (a, b) in canonical directions, compute intercept
    c = b*x - a*y for every prime. Bucket primes by (a, b, c). The set
    of buckets across all directions is the candidate line set. Greedy:
    pick the largest bucket, mark its primes covered, update other
    buckets containing those primes. Repeat.

    Returns (L, chosen_lines).
    """
    import heapq
    M = len(prime_X)
    dirs = canonical_directions(K)
    if verbose:
        print(f"#   bounded-greedy: M={M} primes, {len(dirs)} canonical directions, K_max={K}")
        sys.stdout.flush()

    # For each direction, compute intercepts for all primes and bucket.
    # line_indices: line_key (a, b, c) -> list of prime indices on it.
    line_indices: dict[tuple[int, int, int], list[int]] = {}
    point_to_lines: list[list[tuple[int, int, int]]] = [[] for _ in range(M)]
    for (a, b) in dirs:
        intercepts = b * prime_X.astype(np.int64) - a * prime_Y.astype(np.int64)
        # Group prime indices by intercept
        order = np.argsort(intercepts, kind="stable")
        sorted_ints = intercepts[order]
        # Find run boundaries
        # Manually iterate: faster for moderate M
        if M == 0:
            continue
        i = 0
        while i < M:
            j = i
            while j < M and sorted_ints[j] == sorted_ints[i]:
                j += 1
            if j - i >= 2:
                key = (a, b, int(sorted_ints[i]))
                idxs = order[i:j].tolist()
                line_indices[key] = idxs
                for ii in idxs:
                    point_to_lines[ii].append(key)
            i = j

    if verbose:
        print(f"#   distinct multi-prime lines (bounded): {len(line_indices)}")
        sys.stdout.flush()

    # Greedy with heap.
    line_remaining: dict[tuple[int, int, int], set] = {k: set(v) for k, v in line_indices.items()}
    heap: list = []
    for k, members in line_remaining.items():
        heapq.heappush(heap, (-len(members), k))
    covered = set()
    chosen = []
    while heap:
        neg_sz, key = heapq.heappop(heap)
        if key not in line_remaining:
            continue
        actual = len(line_remaining[key])
        if actual != -neg_sz:
            if actual >= 2:
                heapq.heappush(heap, (-actual, key))
            continue
        if actual < 2:
            continue
        members = line_remaining.pop(key)
        chosen.append((key, sorted(members)))
        covered.update(members)
        for i in members:
            for other in point_to_lines[i]:
                if other == key:
                    continue
                if other in line_remaining:
                    line_remaining[other].discard(i)
                    new_sz = len(line_remaining[other])
                    if new_sz < 2:
                        del line_remaining[other]
                    else:
                        heapq.heappush(heap, (-new_sz, other))
    L = len(chosen) + (M - len(covered))
    return L, chosen


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--N", type=int, default=1_000_000)
    p.add_argument("--K", type=int, default=30, help="max direction component for canonical directions")
    p.add_argument("--baseline-trials", type=int, default=3)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--top-k-lines", type=int, default=20)
    p.add_argument("--out-dir", type=str, default=".")
    args = p.parse_args()

    rng = np.random.default_rng(args.seed)
    pyrng = random.Random(args.seed)
    out_dir = args.out_dir
    os.makedirs(out_dir, exist_ok=True)

    print(f"# Thread 11 / Slot 1 — Ulam-spiral bounded-direction line cover, N={args.N}, K={args.K}")
    sys.stdout.flush()

    t0 = time.time()
    X, Y = ulam_coords_array(args.N)
    primes = sieve_primes(args.N)
    pi_N = len(primes)
    print(f"# pi(N)={pi_N}, ulam built in {time.time()-t0:.2f}s")
    sys.stdout.flush()

    # Step 1: prime line cover
    prime_X = X[primes]
    prime_Y = Y[primes]
    t0 = time.time()
    L_prime, chosen_prime = bounded_greedy_line_cover(prime_X, prime_Y, args.K, verbose=True)
    print(f"# L_Ulam_bounded(primes, N={args.N}, K={args.K}) = {L_prime}, time {time.time()-t0:.2f}s")
    print(f"#   chosen multi-lines: {len(chosen_prime)}, singletons: {L_prime - len(chosen_prime)}")
    sys.stdout.flush()

    chosen_sorted = sorted(chosen_prime, key=lambda x: -len(x[1]))[: args.top_k_lines]
    print(f"\n# Top-{args.top_k_lines} densest greedy-chosen lines on primes (Ulam, bounded):")
    print(f"# {'rank':>4} {'a':>4} {'b':>4} {'intercept':>10} {'count':>6}")
    for r, ((a, b, c), members) in enumerate(chosen_sorted, 1):
        print(f"# {r:>4} {a:>4} {b:>4} {c:>10} {len(members):>6}")

    # Step 2: matched-baseline
    print(f"\n# Matched-baseline (random integers in [1,N])")
    L_random_list = []
    for trial in range(args.baseline_trials):
        t0 = time.time()
        sample_idxs = rng.choice(np.arange(1, args.N + 1), size=pi_N, replace=False)
        rx = X[sample_idxs]
        ry = Y[sample_idxs]
        L_r, _ = bounded_greedy_line_cover(rx, ry, args.K, verbose=False)
        L_random_list.append(L_r)
        print(f"#   trial {trial+1}: L_random = {L_r}, time {time.time()-t0:.2f}s")
        sys.stdout.flush()
    L_random_mean = sum(L_random_list) / len(L_random_list)
    L_random_var = sum((x - L_random_mean) ** 2 for x in L_random_list) / max(1, len(L_random_list) - 1)
    L_random_std = math.sqrt(L_random_var) if len(L_random_list) > 1 else 0.0
    print(f"#   L_random_mean = {L_random_mean:.1f} +/- {L_random_std:.1f}")
    if L_random_std > 0:
        z = (L_random_mean - L_prime) / L_random_std
        print(f"#   z-score (random - prime) / std = {z:+.2f}sigma  (positive => primes have FEWER lines)")
    else:
        z = float("nan")

    # Save summary CSV
    summary_path = os.path.join(out_dir, f"summary_bounded_N{args.N}_K{args.K}.csv")
    with open(summary_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["N", "K", "pi_N", "L_primes", "L_random_mean", "L_random_std",
                    "n_baseline_trials", "z_score"])
        w.writerow([args.N, args.K, pi_N, L_prime, f"{L_random_mean:.2f}", f"{L_random_std:.2f}",
                    args.baseline_trials, f"{z:.3f}" if not math.isnan(z) else "NA"])

    toplines_path = os.path.join(out_dir, f"top_lines_bounded_N{args.N}_K{args.K}.csv")
    with open(toplines_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["rank", "a", "b", "intercept", "prime_count"])
        for r, ((a, b, c), members) in enumerate(chosen_sorted, 1):
            w.writerow([r, a, b, c, len(members)])

    print(f"\n# wrote {summary_path}, {toplines_path}")
    print(f"# DONE bounded N={args.N}")


if __name__ == "__main__":
    main()

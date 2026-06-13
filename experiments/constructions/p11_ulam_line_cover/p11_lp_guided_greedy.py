"""
Thread 11 / Slot 4 — LP-guided greedy variant for line cover.

Rather than thresholding the LP solution at fixed tau, sort active lines
by x_l * |line| (expected coverage in LP) and greedily pick the line of
maximum *uncovered* size, restricted to LP-active support.

Compare with: pure greedy (no LP), iterated rounding (tau=0.5), and
LP-guided greedy (this script).

CLI:
  python p11_lp_guided_greedy.py --embedding ulam --N 10000 --K 5
  python p11_lp_guided_greedy.py --embedding ulam --N 100000 --K 5
"""
from __future__ import annotations

import argparse
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p11_ulam_bounded import (
    sieve_primes,
    ulam_coords_array,
    bounded_greedy_line_cover,
)
from p11_alt_embeddings import residue_class_coords, polynomial_image_coords
from p11_lp_relaxation import enumerate_lines
from p11_iterated_rounding import solve_set_cover_lp


def lp_guided_greedy(line_indices: dict, M: int, lp_eps: float = 1e-6,
                      verbose: bool = True) -> int:
    """LP-guided greedy: solve LP once, then greedy among LP-active lines.

    1. Solve LP. Collect active lines L_act = {l : x_l > eps}.
    2. Greedy on L_act: pick line covering most uncovered primes.
    3. Repeat until done. Cover remaining with singletons.
    """
    t0 = time.time()
    lp_val, x, line_keys = solve_set_cover_lp(set(range(M)), line_indices)
    t_lp = time.time() - t0
    if verbose:
        print(f"# LP solved in {t_lp:.2f}s, LP value = {lp_val:.4f}")

    # Active lines
    active_idx = [i for i in range(len(x)) if x[i] > lp_eps]
    if verbose:
        print(f"# active LP lines: {len(active_idx)} / {len(line_keys)}")

    # Build line membership for active lines
    active_lines = []
    for i in active_idx:
        key = line_keys[i]
        members = set(line_indices[key])
        active_lines.append((key, members, x[i]))

    # Greedy
    covered = set()
    chosen = []
    while True:
        best_size = 0
        best_idx = -1
        for j, (key, members, w) in enumerate(active_lines):
            uncov_in_line = len(members - covered)
            if uncov_in_line > best_size:
                best_size = uncov_in_line
                best_idx = j
        if best_size < 2 or best_idx < 0:
            break
        key, members, w = active_lines[best_idx]
        chosen.append(key)
        covered |= members
    L_total = len(chosen) + (M - len(covered))
    if verbose:
        print(f"# LP-guided greedy: {len(chosen)} multi-lines + "
              f"{M - len(covered)} singletons = {L_total}")
    return L_total


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--embedding", type=str, default="ulam",
                   choices=["residue", "poly", "ulam"])
    p.add_argument("--q", type=int, default=210)
    p.add_argument("--N", type=int, default=10_000)
    p.add_argument("--K", type=int, default=5)
    p.add_argument("--baseline-seed", type=int, default=None)
    p.add_argument("--out-dir", type=str, default=".")
    args = p.parse_args()

    primes = sieve_primes(args.N)
    pi_N = len(primes)
    print(f"# Thread 11 / Slot 4 — LP-guided greedy, "
          f"embedding={args.embedding}, N={args.N}, K={args.K}")
    print(f"# pi(N)={pi_N}")
    sys.stdout.flush()

    if args.embedding == "ulam":
        X32, Y32 = ulam_coords_array(args.N)
        X = X32.astype(np.int64); Y = Y32.astype(np.int64)
        name = "ulam"
    elif args.embedding == "residue":
        X, Y = residue_class_coords(args.N, args.q)
        name = f"residue_q{args.q}"
    else:
        X, Y = polynomial_image_coords(args.N, args.q)
        name = f"poly_q{args.q}"

    if args.baseline_seed is not None:
        rng = np.random.default_rng(args.baseline_seed)
        sample_idxs = rng.choice(np.arange(1, args.N + 1), size=pi_N, replace=False)
        prime_X = X[sample_idxs]; prime_Y = Y[sample_idxs]
        label = f"{name}_baseline_seed{args.baseline_seed}"
    else:
        prime_X = X[primes]; prime_Y = Y[primes]
        label = f"{name}_primes"

    M = pi_N
    t0 = time.time()
    L_greedy, _ = bounded_greedy_line_cover(prime_X, prime_Y, args.K, verbose=False)
    print(f"# bounded greedy L = {L_greedy} (in {time.time()-t0:.2f}s)")

    line_indices = enumerate_lines(prime_X, prime_Y, args.K)
    print(f"# {len(line_indices)} multi-prime candidate lines")

    L_guided = lp_guided_greedy(line_indices, M, verbose=True)

    out_path = os.path.join(args.out_dir,
                             f"lp_guided_{label}_N{args.N}_K{args.K}.txt")
    with open(out_path, "w") as f:
        f.write(f"label: {label}\nN: {args.N}\nK: {args.K}\nM: {M}\n")
        f.write(f"L_greedy: {L_greedy}\n")
        f.write(f"L_lp_guided: {L_guided}\n")
        f.write(f"diff: {L_guided - L_greedy}\n")
    print(f"# wrote {out_path}")


if __name__ == "__main__":
    main()

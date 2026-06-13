"""
Thread 11 / Slot 3 — extract LP solution lines and dual weights.

Solves the LP relaxation, returns (a, b, c) -> x_l > eps for inspection,
and matches against bounded-greedy chosen lines. Also extracts LP dual
y_p (per-prime fractional weight) to identify "expensive" primes.

CLI:
  python p11_lp_dual_inspect.py --embedding ulam --N 10000 --K 5
"""
from __future__ import annotations

import argparse
import os
import sys
import time

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import csr_matrix

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p11_ulam_bounded import (
    sieve_primes,
    ulam_coords_array,
    bounded_greedy_line_cover,
)
from p11_alt_embeddings import residue_class_coords, polynomial_image_coords
from p11_lp_relaxation import enumerate_lines


def solve_lp_with_dual(M: int, line_indices: dict):
    on_a_line = [False] * M
    line_keys = list(line_indices.keys())
    L = len(line_keys)
    rows_idx, cols_idx, data = [], [], []
    for col, key in enumerate(line_keys):
        for prime_idx in line_indices[key]:
            on_a_line[prime_idx] = True
            rows_idx.append(prime_idx)
            cols_idx.append(col)
            data.append(1)
    n_singletons = 0
    singleton_keys = []
    for p in range(M):
        if not on_a_line[p]:
            rows_idx.append(p)
            cols_idx.append(L + n_singletons)
            data.append(1)
            singleton_keys.append(p)
            n_singletons += 1
    n_vars = L + n_singletons
    A = csr_matrix((data, (rows_idx, cols_idx)), shape=(M, n_vars))
    A_ub = -A
    b_ub = -np.ones(M)
    c = np.ones(n_vars)
    bounds = [(0.0, None)] * n_vars
    res = linprog(c, A_ub=A_ub, b_ub=b_ub, bounds=bounds, method="highs",
                  options={"presolve": True})
    if not res.success:
        raise RuntimeError(res.message)

    x = res.x
    # Dual: linprog returns marginals (Lagrange multipliers).
    dual = None
    if hasattr(res, "ineqlin") and res.ineqlin is not None:
        dual = -res.ineqlin.marginals  # marginals are for A_ub @ x <= b_ub; we want for A @ x >= 1.

    line_weights = []
    for col, key in enumerate(line_keys):
        if x[col] > 1e-7:
            line_weights.append((key, x[col], len(line_indices[key])))
    line_weights.sort(key=lambda t: -t[1])

    n_singletons_active = sum(1 for col in range(L, n_vars) if x[col] > 1e-7)
    return res.fun, x, dual, line_weights, n_singletons_active, line_keys


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--embedding", type=str, default="ulam",
                   choices=["residue", "poly", "ulam"])
    p.add_argument("--q", type=int, default=210)
    p.add_argument("--N", type=int, default=10_000)
    p.add_argument("--K", type=int, default=5)
    p.add_argument("--top", type=int, default=20)
    p.add_argument("--out-dir", type=str, default=".")
    p.add_argument("--baseline-seed", type=int, default=None,
                   help="if set, run on matched-baseline random points")
    args = p.parse_args()

    primes = sieve_primes(args.N)
    pi_N = len(primes)
    print(f"# pi(N)={pi_N}")

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

    print(f"# {label} N={args.N} K={args.K}")
    line_indices = enumerate_lines(prime_X, prime_Y, args.K)
    M = pi_N
    print(f"# {len(line_indices)} multi-prime lines")
    sys.stdout.flush()

    L_lp, x, dual, line_weights, n_sing_active, line_keys = solve_lp_with_dual(M, line_indices)
    print(f"# LP optimum = {L_lp:.4f}")
    print(f"# active lines (x > 0): {len(line_weights)}, active singletons: {n_sing_active}")

    # Top lines by weight
    print(f"\n## Top-{args.top} active LP lines (highest weight)")
    print(f"# {'rank':>4} {'a':>4} {'b':>4} {'intercept':>12} "
          f"{'#prime':>8} {'x_l':>10}")
    integer_count = 0
    for r, (key, w, sz) in enumerate(line_weights[:args.top], 1):
        a, b, c_int = key
        is_int = abs(w - round(w)) < 1e-6
        if is_int and abs(round(w) - 1.0) < 1e-6:
            integer_count += 1
        print(f"# {r:>4} {a:>4} {b:>4} {c_int:>12} {sz:>8} {w:>10.6f}")

    # Lines at exactly x=1.0 (the integer-feasible part of LP solution)
    n_x1 = sum(1 for (_, w, _) in line_weights if abs(w - 1.0) < 1e-6)
    n_frac = sum(1 for (_, w, _) in line_weights if w < 1 - 1e-6 and w > 1e-6)
    print(f"\n# lines at x=1.0: {n_x1}; fractional lines (0<x<1): {n_frac}")

    # Total weight by direction (a, b)
    direction_weight = {}
    for (a, b, c_int), w, sz in line_weights:
        direction_weight[(a, b)] = direction_weight.get((a, b), 0.0) + w
    print(f"\n## Top directions by total weight (sum_l x_l per direction)")
    sorted_dirs = sorted(direction_weight.items(), key=lambda t: -t[1])[:15]
    for (a, b), w in sorted_dirs:
        print(f"#  ({a:>3}, {b:>3}): total_weight = {w:.4f}")


if __name__ == "__main__":
    main()

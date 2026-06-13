"""
Thread 11 / Slot 4 — true integer optimum via MILP (HiGHS branch-and-bound).

Test if the integrality gap LP < OPT < greedy is real (i.e. the structural
gap is 1.0) or a heuristic artefact (i.e. true OPT = LP). Use scipy.optimize.milp
on the set-cover IP at small N (10^3, 5000, 10^4 only — IP is NP-hard).

Compare:
  LP      lower bound (slot 3)
  greedy  upper bound (slot 1)
  iter    iterated rounding upper bound (slot 4)
  OPT     true integer optimum (this script, MILP)
"""
from __future__ import annotations

import argparse
import os
import sys
import time

import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import csr_matrix

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p11_ulam_bounded import (
    sieve_primes,
    ulam_coords_array,
    bounded_greedy_line_cover,
)
from p11_lp_relaxation import enumerate_lines


def build_set_cover_matrix(M: int, line_indices: dict):
    on_a_line = [False] * M
    line_keys = list(line_indices.keys())
    L = len(line_keys)
    rows, cols, data = [], [], []
    for col, key in enumerate(line_keys):
        for prime_idx in line_indices[key]:
            on_a_line[prime_idx] = True
            rows.append(prime_idx); cols.append(col); data.append(1)
    n_singletons = 0
    for p in range(M):
        if not on_a_line[p]:
            rows.append(p); cols.append(L + n_singletons); data.append(1)
            n_singletons += 1
    n_vars = L + n_singletons
    A = csr_matrix((data, (rows, cols)), shape=(M, n_vars))
    return A, line_keys, n_singletons


def solve_ip(M: int, line_indices: dict, time_limit: float = 600.0):
    A, line_keys, n_singletons = build_set_cover_matrix(M, line_indices)
    n_vars = A.shape[1]
    c = np.ones(n_vars)
    constraints = LinearConstraint(A, lb=np.ones(M), ub=np.full(M, np.inf))
    integrality = np.ones(n_vars)  # all binary 0/1 effectively
    bounds = Bounds(lb=0.0, ub=1.0)
    options = {"time_limit": time_limit, "disp": True}
    t0 = time.time()
    res = milp(c, constraints=constraints, integrality=integrality,
               bounds=bounds, options=options)
    t = time.time() - t0
    if not res.success:
        print(f"# MILP did not finish: {res.message}, status={res.status}")
        if hasattr(res, "fun") and res.fun is not None:
            print(f"# best feasible found = {res.fun}")
        return None, t, res
    return float(res.fun), t, res


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--N", type=int, default=10_000)
    p.add_argument("--K", type=int, default=5)
    p.add_argument("--time-limit", type=float, default=600.0)
    p.add_argument("--out-dir", type=str, default=".")
    args = p.parse_args()

    primes = sieve_primes(args.N)
    pi_N = len(primes)
    print(f"# Thread 11 / Slot 4 — MILP integer optimum, N={args.N}, K={args.K}")
    print(f"# pi(N)={pi_N}")
    sys.stdout.flush()

    X32, Y32 = ulam_coords_array(args.N)
    X = X32.astype(np.int64); Y = Y32.astype(np.int64)
    prime_X = X[primes]; prime_Y = Y[primes]
    M = pi_N

    t0 = time.time()
    L_greedy, _ = bounded_greedy_line_cover(prime_X, prime_Y, args.K, verbose=False)
    print(f"# greedy L = {L_greedy} (in {time.time()-t0:.2f}s)")

    line_indices = enumerate_lines(prime_X, prime_Y, args.K)
    print(f"# {len(line_indices)} multi-prime lines")
    sys.stdout.flush()

    print(f"# solving MILP, time_limit={args.time_limit}s ...")
    sys.stdout.flush()
    L_ip, t_ip, res = solve_ip(M, line_indices, args.time_limit)
    if L_ip is not None:
        print(f"# MILP integer optimum = {L_ip:.4f} in {t_ip:.2f}s")
    else:
        print(f"# MILP did not converge in {t_ip:.2f}s")

    out_path = os.path.join(args.out_dir, f"milp_ulam_N{args.N}_K{args.K}.txt")
    with open(out_path, "w") as f:
        f.write(f"N: {args.N}\nK: {args.K}\nM: {M}\n")
        f.write(f"L_greedy: {L_greedy}\n")
        if L_ip is not None:
            f.write(f"L_ip: {L_ip}\n")
            f.write(f"t_ip: {t_ip:.2f}\n")
        else:
            f.write(f"L_ip: timeout\n")
            f.write(f"t_ip: {t_ip:.2f}\n")
    print(f"# wrote {out_path}")


if __name__ == "__main__":
    main()

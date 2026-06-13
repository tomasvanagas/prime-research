"""
Thread 11 / Slot 3 — LP-relaxation lower bound on minimum line cover.

Slots 1-2 (S484, S485) measured greedy line-cover sizes for primes vs
matched-baseline random points across two embeddings (Ulam spiral and
residue-class / polynomial-image grids). Greedy is an upper bound on
the integer optimum; the LP relaxation is a *lower* bound. Slot 3
asks: does the LP relaxation give a strict improvement over greedy
(integrality gap > 1), or does the wheel-sieve / Ulam-greedy already
hit the LP floor?

Cross-domain ingredient: Stanley 1989 *Adv. Math.* matroid-theoretic
line-cover-LP bound; LP duality on the set-cover hypergraph.

Setup:
- universe U = {primes p <= N}, embedded under Phi: N -> Z^2.
- candidate lines L = lines through >=2 prime points. Restrict to
  bounded directions (a, b) with gcd=1, max <= K_max for tractability.
- LP: minimise sum_{l in L} x_l subject to sum_{l: p in l} x_l >= 1
  for all p, x_l >= 0.

Two concrete tests (per the slot-3 plan in .commit_state):
  (a) residue_q=210 at N=10^5 with greedy L_p=49 — does LP push <49?
  (b) Ulam at N=10^4 with greedy L_p=91 — does LP give same?

If LP = greedy in both cases, the wheel-sieve / Ulam-greedy is
LP-tight: no fractional cover does better. This is the structural
"close" outcome forecast in slot-2.

If LP < greedy, there's slack — interesting integrality gap.

CLI:
  python p11_lp_relaxation.py --embedding residue --q 210 --N 100000
  python p11_lp_relaxation.py --embedding ulam --N 10000 --K 5
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
    canonical_directions,
    bounded_greedy_line_cover,
)
from p11_alt_embeddings import residue_class_coords, polynomial_image_coords


def enumerate_lines(prime_X: np.ndarray, prime_Y: np.ndarray, K: int):
    """Return dict: (a, b, c) -> sorted list of prime indices on the line.
    Using the same bounded-direction parameterisation as the greedy.
    """
    M = len(prime_X)
    dirs = canonical_directions(K)
    line_indices: dict[tuple[int, int, int], list[int]] = {}
    for (a, b) in dirs:
        intercepts = b * prime_X.astype(np.int64) - a * prime_Y.astype(np.int64)
        order = np.argsort(intercepts, kind="stable")
        sorted_ints = intercepts[order]
        i = 0
        while i < M:
            j = i
            while j < M and sorted_ints[j] == sorted_ints[i]:
                j += 1
            if j - i >= 2:
                key = (a, b, int(sorted_ints[i]))
                idxs = order[i:j].tolist()
                line_indices[key] = idxs
            i = j
    return line_indices


def solve_lp(M: int, line_indices: dict, verbose: bool = False) -> tuple[float, float, dict]:
    """Solve the set-cover LP relaxation.

    Variables: one per multi-prime line + one per singleton prime
    (unique-line variable, cost 1).
    Constraints: each prime covered to >= 1 in fractional weight.
    Objective: min sum x_l.

    Returns (lp_value, runtime_seconds, info).
    """
    # Identify singletons (primes not on any line).
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

    # Add a singleton column for each prime not on any multi-line.
    n_singletons = 0
    for p in range(M):
        if not on_a_line[p]:
            rows_idx.append(p)
            cols_idx.append(L + n_singletons)
            data.append(1)
            n_singletons += 1

    n_vars = L + n_singletons
    A = csr_matrix((data, (rows_idx, cols_idx)), shape=(M, n_vars))
    # linprog uses A_ub @ x <= b_ub; rewrite "A @ x >= 1" as "(-A) @ x <= -1".
    A_ub = -A
    b_ub = -np.ones(M)
    c = np.ones(n_vars)
    bounds = [(0.0, None)] * n_vars
    if verbose:
        print(f"#   LP: {M} primes, {L} multi-lines, {n_singletons} singletons, "
              f"{n_vars} vars, {len(data)} nonzeros")
        sys.stdout.flush()
    t0 = time.time()
    res = linprog(c, A_ub=A_ub, b_ub=b_ub, bounds=bounds,
                  method="highs", options={"presolve": True})
    t = time.time() - t0
    if not res.success:
        raise RuntimeError(f"LP did not solve: {res.message}")
    info = {
        "status": res.status,
        "n_lines": L,
        "n_singletons": n_singletons,
        "n_vars": n_vars,
        "n_nonzeros": len(data),
        "n_iter": getattr(res, "nit", None),
    }
    return float(res.fun), t, info


def run_embedding(name: str, X_full: np.ndarray, Y_full: np.ndarray,
                  primes: np.ndarray, K: int) -> dict:
    M = len(primes)
    prime_X = X_full[primes]
    prime_Y = Y_full[primes]

    print(f"\n## {name}: M={M} primes, K_max={K}")
    sys.stdout.flush()

    # Greedy
    t0 = time.time()
    L_greedy, _ = bounded_greedy_line_cover(prime_X, prime_Y, K, verbose=False)
    t_greedy = time.time() - t0
    print(f"   greedy L = {L_greedy} (in {t_greedy:.2f}s)")
    sys.stdout.flush()

    # Enumerate lines
    t0 = time.time()
    line_indices = enumerate_lines(prime_X, prime_Y, K)
    t_enum = time.time() - t0
    print(f"   {len(line_indices)} multi-prime candidate lines (in {t_enum:.2f}s)")
    sys.stdout.flush()

    # LP
    L_lp, t_lp, lp_info = solve_lp(M, line_indices, verbose=True)
    print(f"   LP relaxation = {L_lp:.4f} (in {t_lp:.2f}s)")
    print(f"   integrality_gap = greedy / LP = {L_greedy / max(L_lp, 1e-9):.4f}")
    sys.stdout.flush()

    return {
        "name": name,
        "M": M,
        "K": K,
        "L_greedy": L_greedy,
        "L_lp": L_lp,
        "n_lines": lp_info["n_lines"],
        "n_singletons": lp_info["n_singletons"],
        "n_vars": lp_info["n_vars"],
        "n_nonzeros": lp_info["n_nonzeros"],
        "t_greedy": t_greedy,
        "t_enum": t_enum,
        "t_lp": t_lp,
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--embedding", type=str, default="residue",
                   choices=["residue", "poly", "ulam"])
    p.add_argument("--q", type=int, default=210, help="modulus for residue/poly embeddings")
    p.add_argument("--N", type=int, default=100_000)
    p.add_argument("--K", type=int, default=5)
    p.add_argument("--out-dir", type=str, default=".")
    p.add_argument("--baseline-trials", type=int, default=0,
                   help="if >0, also LP-solve matched-baseline random points")
    p.add_argument("--seed", type=int, default=42)
    args = p.parse_args()

    print(f"# Thread 11 / Slot 3 — LP relaxation, embedding={args.embedding}, "
          f"q={args.q}, N={args.N}, K={args.K}")
    sys.stdout.flush()

    primes = sieve_primes(args.N)
    pi_N = len(primes)
    print(f"# pi(N)={pi_N}")

    if args.embedding == "residue":
        X, Y = residue_class_coords(args.N, args.q)
        name = f"residue_q{args.q}_N{args.N}"
    elif args.embedding == "poly":
        X, Y = polynomial_image_coords(args.N, args.q)
        name = f"poly_q{args.q}_N{args.N}"
    elif args.embedding == "ulam":
        X32, Y32 = ulam_coords_array(args.N)
        X = X32.astype(np.int64)
        Y = Y32.astype(np.int64)
        name = f"ulam_N{args.N}"
    else:
        raise ValueError(args.embedding)

    r = run_embedding(name, X, Y, primes, args.K)

    baseline_results = []
    if args.baseline_trials > 0:
        rng = np.random.default_rng(args.seed)
        for trial in range(args.baseline_trials):
            sample_idxs = rng.choice(np.arange(1, args.N + 1),
                                      size=pi_N, replace=False)
            br = run_embedding(f"{name}_baseline_t{trial}",
                               X, Y, sample_idxs, args.K)
            baseline_results.append(br)

    out_path = os.path.join(args.out_dir,
                             f"lp_summary_{r['name']}_K{args.K}.txt")
    with open(out_path, "w") as f:
        for k, v in r.items():
            f.write(f"{k}: {v}\n")
        for i, br in enumerate(baseline_results):
            f.write(f"\n# baseline trial {i}\n")
            for k, v in br.items():
                f.write(f"{k}: {v}\n")

    if baseline_results:
        L_lp_b = [b["L_lp"] for b in baseline_results]
        L_g_b = [b["L_greedy"] for b in baseline_results]
        mean_lp = sum(L_lp_b) / len(L_lp_b)
        mean_g = sum(L_g_b) / len(L_g_b)
        print(f"\n# baseline LP mean: {mean_lp:.2f}, greedy mean: {mean_g:.2f}")
        print(f"# primes vs baseline: LP {r['L_lp']:.2f} vs {mean_lp:.2f}, "
              f"greedy {r['L_greedy']} vs {mean_g:.2f}")

    print(f"# wrote {out_path}")
    print(f"# DONE LP {r['name']} K={args.K}")


if __name__ == "__main__":
    main()

"""
Thread 11 / Slot 4 — iterated LP rounding for integer line cover.

Slot 3 (S486) found a stable LP-vs-greedy integrality gap for primes on
the Ulam spiral: greedy is ~1.22-1.26 above LP at N in {10^4, 10^5}.
Greedy(primes) / LP(primes) ~ 1.22 means greedy is leaving 22% on the
table. The LP optimum LP_p ~ 0.78 sqrt(N) is achievable by integer
solutions only if a rounding scheme exists.

Cross-domain ingredient: iterated LP rounding (Singh-Lau, Lavi-Swamy,
Raghavan-Thompson 1987 randomised rounding for set cover). Standard
guarantee: O(log n)-approximation matches greedy for unstructured set
cover (Chvatal 1979); the structural hope here is that the prime LP
puts heavy weight on a *small* set of HL-rich diagonals, so iterated
threshold rounding may achieve much better than O(log) factor.

Algorithm (iterated rounding with threshold tau and residual LP):
  1. Solve LP. Get fractional solution x*.
  2. For all lines l with x*_l >= tau: round to 1 (select l in cover).
  3. Mark all primes covered by these lines.
  4. Re-solve LP on the residual problem (uncovered primes).
  5. Repeat until residual LP integral or all primes covered.
  6. Final fallback: cover remaining primes with singleton lines.

Compare: pure greedy vs iterated rounding, in terms of integer cover
size produced.

CLI:
  python p11_iterated_rounding.py --embedding ulam --N 10000 --K 5 --tau 0.5
  python p11_iterated_rounding.py --embedding ulam --N 100000 --K 5 --tau 0.3
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


def solve_set_cover_lp(uncovered_set: set, line_indices: dict) -> tuple[float, np.ndarray, list]:
    """Solve LP on residual: cover all primes in uncovered_set with
    fractional weighting on lines (intersected with uncovered).
    Returns (lp_value, x, line_keys).
    """
    if not uncovered_set:
        return 0.0, np.array([]), []
    uncov_list = sorted(uncovered_set)
    uncov_index = {p: i for i, p in enumerate(uncov_list)}
    M = len(uncov_list)

    # Restrict each line to uncovered members
    line_keys: list = []
    line_members: list = []
    for key, members in line_indices.items():
        restricted = [uncov_index[p] for p in members if p in uncov_index]
        if len(restricted) >= 2:
            line_keys.append(key)
            line_members.append(restricted)

    rows, cols, data = [], [], []
    for col, members in enumerate(line_members):
        for r in members:
            rows.append(r); cols.append(col); data.append(1)
    L = len(line_keys)
    # Singletons for primes not on any restricted multi-line
    on_line = np.zeros(M, dtype=bool)
    for r in rows:
        on_line[r] = True
    n_sing = 0
    for r in range(M):
        if not on_line[r]:
            rows.append(r); cols.append(L + n_sing); data.append(1)
            n_sing += 1
    n_vars = L + n_sing
    A = csr_matrix((data, (rows, cols)), shape=(M, n_vars))
    A_ub = -A
    b_ub = -np.ones(M)
    c = np.ones(n_vars)
    bounds = [(0.0, None)] * n_vars
    res = linprog(c, A_ub=A_ub, b_ub=b_ub, bounds=bounds, method="highs",
                  options={"presolve": True})
    if not res.success:
        raise RuntimeError(res.message)
    return float(res.fun), res.x[:L], line_keys


def iterated_rounding(line_indices: dict, M: int, tau: float = 0.5,
                       max_rounds: int = 20, verbose: bool = True) -> dict:
    """Iterated threshold rounding."""
    uncovered = set(range(M))
    chosen: list = []
    history: list = []

    for rnd in range(max_rounds):
        if not uncovered:
            break
        lp_val, x, line_keys = solve_set_cover_lp(uncovered, line_indices)
        if not line_keys:
            # Only singletons feasible
            history.append({"round": rnd, "uncovered_in": len(uncovered),
                            "lp_value": lp_val, "rounded": 0,
                            "covered_by_round": 0})
            break
        # Round lines with x >= tau to 1
        rounded_idx = [i for i in range(len(x)) if x[i] >= tau - 1e-9]
        rounded_keys = [line_keys[i] for i in rounded_idx]

        # Apply rounded lines (using ORIGINAL line_indices membership;
        # cover everything they cover, not just restricted)
        newly_covered = set()
        for key in rounded_keys:
            for p in line_indices[key]:
                if p in uncovered:
                    newly_covered.add(p)

        history.append({
            "round": rnd,
            "uncovered_in": len(uncovered),
            "lp_value": lp_val,
            "rounded": len(rounded_keys),
            "covered_by_round": len(newly_covered),
        })
        if verbose:
            print(f"#   round {rnd}: uncov_in={len(uncovered)}, "
                  f"lp={lp_val:.2f}, rounded={len(rounded_keys)}, "
                  f"newly_covered={len(newly_covered)}")
            sys.stdout.flush()

        if not rounded_keys or not newly_covered:
            # No progress with this tau; lower tau adaptively
            new_tau = max(tau * 0.7, 0.05)
            if abs(new_tau - tau) < 1e-3:
                # Final fallback: round the highest-weight line
                if len(x) == 0:
                    break
                best = int(np.argmax(x))
                rounded_keys = [line_keys[best]]
                newly_covered = set()
                for p in line_indices[line_keys[best]]:
                    if p in uncovered:
                        newly_covered.add(p)
                if not newly_covered:
                    break
            else:
                if verbose:
                    print(f"#   no progress at tau={tau}, lowering to {new_tau}")
                tau = new_tau
                continue

        chosen.extend(rounded_keys)
        uncovered -= newly_covered

    # Cover remaining with singletons (one line per prime)
    n_singletons = len(uncovered)
    L_total = len(chosen) + n_singletons
    return {
        "L_iterated": L_total,
        "n_multi_lines": len(chosen),
        "n_singletons": n_singletons,
        "n_rounds": len(history),
        "history": history,
        "chosen": chosen,
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--embedding", type=str, default="ulam",
                   choices=["residue", "poly", "ulam"])
    p.add_argument("--q", type=int, default=210)
    p.add_argument("--N", type=int, default=10_000)
    p.add_argument("--K", type=int, default=5)
    p.add_argument("--tau", type=float, default=0.5,
                   help="initial rounding threshold; adapts down on no-progress")
    p.add_argument("--baseline-seed", type=int, default=None)
    p.add_argument("--out-dir", type=str, default=".")
    args = p.parse_args()

    primes = sieve_primes(args.N)
    pi_N = len(primes)
    print(f"# Thread 11 / Slot 4 — iterated rounding, "
          f"embedding={args.embedding}, N={args.N}, K={args.K}, tau={args.tau}")
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

    # Reference: greedy
    t0 = time.time()
    L_greedy, _ = bounded_greedy_line_cover(prime_X, prime_Y, args.K, verbose=False)
    t_greedy = time.time() - t0
    print(f"# greedy L = {L_greedy} (in {t_greedy:.2f}s)")

    # Reference: LP relaxation lower bound (one-shot, no rounding)
    line_indices = enumerate_lines(prime_X, prime_Y, args.K)
    print(f"# {len(line_indices)} multi-prime candidate lines")
    t0 = time.time()
    lp_init, _, _ = solve_set_cover_lp(set(range(M)), line_indices)
    t_lp_init = time.time() - t0
    print(f"# LP one-shot lower bound = {lp_init:.2f} (in {t_lp_init:.2f}s)")

    # Iterated rounding
    t0 = time.time()
    res = iterated_rounding(line_indices, M, tau=args.tau, max_rounds=30, verbose=True)
    t_iter = time.time() - t0

    print(f"\n# iterated_rounding L = {res['L_iterated']} "
          f"(multi={res['n_multi_lines']}, sing={res['n_singletons']}) "
          f"in {t_iter:.2f}s, {res['n_rounds']} rounds")
    print(f"# integrality gap (greedy / LP)  = {L_greedy / max(lp_init, 1e-9):.4f}")
    print(f"# integrality gap (iterated / LP) = {res['L_iterated'] / max(lp_init, 1e-9):.4f}")
    print(f"# improvement greedy -> iterated:  {L_greedy - res['L_iterated']:+d} "
          f"({100*(res['L_iterated']/L_greedy - 1):+.2f}% relative)")

    out_path = os.path.join(args.out_dir,
                             f"iter_round_{label}_N{args.N}_K{args.K}_tau{args.tau}.txt")
    with open(out_path, "w") as f:
        f.write(f"label: {label}\n")
        f.write(f"N: {args.N}\nK: {args.K}\ntau: {args.tau}\nM: {M}\n")
        f.write(f"L_greedy: {L_greedy}\n")
        f.write(f"LP_one_shot: {lp_init:.4f}\n")
        f.write(f"L_iterated: {res['L_iterated']}\n")
        f.write(f"n_multi_lines: {res['n_multi_lines']}\n")
        f.write(f"n_singletons: {res['n_singletons']}\n")
        f.write(f"n_rounds: {res['n_rounds']}\n")
        f.write(f"t_greedy: {t_greedy:.2f}\n")
        f.write(f"t_lp_init: {t_lp_init:.2f}\n")
        f.write(f"t_iter: {t_iter:.2f}\n")
        f.write("# history:\n")
        for h in res["history"]:
            f.write(f"#   {h}\n")
    print(f"# wrote {out_path}")


if __name__ == "__main__":
    main()

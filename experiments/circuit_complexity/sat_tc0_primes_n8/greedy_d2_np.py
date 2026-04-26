"""
Numpy-vectorised greedy depth-2 sign-threshold heuristic.
"""
from __future__ import annotations
import argparse, json, time
import numpy as np

import sys; sys.path.insert(0, '.')
from enum_d2_smart import enumerate_w1_thresholds_pruned, primes_table, random_table

def greedy_d2_np(target_arr: np.ndarray, theta_mat: np.ndarray, max_M: int = 100):
    """
    target_arr: (2^N,) ∈ {0, 1}
    theta_mat:  (K, 2^N) ∈ {0, 1}
    Returns dict with selected gates and signs.
    """
    K, NX = theta_mat.shape
    pos = (target_arr == 1)
    neg = ~pos
    s_top = np.zeros(NX, dtype=np.int32)
    selected = []
    for step in range(max_M):
        # For each candidate k and sign a in {-1, +1}, compute new s_top.
        # gap = min_{x: target=1} (s_top + a*theta[k]) - max_{x: target=0} (s_top + a*theta[k])
        # We want to maximize this gap. If > 0, feasible.
        # Vectorise over (k, a)
        # new_s_top[k, x] = s_top[x] + a * theta_mat[k, x]
        # For a=+1:
        ns_pos = s_top[None, :] + theta_mat  # (K, NX)
        ns_neg = s_top[None, :] - theta_mat
        # min over pos inputs, max over neg inputs (per gate)
        # use masked min/max via where
        big = 10**9
        ns_pos_pos = np.where(pos[None, :], ns_pos, +big)
        ns_pos_neg = np.where(neg[None, :], ns_pos, -big)
        ns_neg_pos = np.where(pos[None, :], ns_neg, +big)
        ns_neg_neg = np.where(neg[None, :], ns_neg, -big)
        gap_pos = ns_pos_pos.min(axis=1) - ns_pos_neg.max(axis=1)
        gap_neg = ns_neg_pos.min(axis=1) - ns_neg_neg.max(axis=1)
        # Find best (k, a)
        idx_p = int(gap_pos.argmax()); best_p = int(gap_pos[idx_p])
        idx_n = int(gap_neg.argmax()); best_n = int(gap_neg[idx_n])
        if best_p >= best_n:
            best_k, best_a, best_g = idx_p, +1, best_p
        else:
            best_k, best_a, best_g = idx_n, -1, best_n
        s_top = s_top + best_a * theta_mat[best_k]
        selected.append({"k": best_k, "alpha": best_a, "gap": best_g})
        if best_g >= 1:  # feasible: s_min_pos - s_max_neg >= 1
            s_max_neg = int(s_top[neg].max(initial=-10**9))
            T_top = s_max_neg + 1
            return {"M": step + 1, "T_top": T_top, "selected": selected, "verified": True}
    return {"M": max_M, "verified": False, "selected": selected}

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--N", type=int, default=8)
    parser.add_argument("--k-max", type=int, default=5)
    parser.add_argument("--max-M", type=int, default=80)
    parser.add_argument("--target", choices=["primes", "random", "both"], default="both")
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()
    cands = enumerate_w1_thresholds_pruned(args.N, k_max=args.k_max)
    print(f"N={args.N}, k_max={args.k_max}: K = {len(cands)} candidates")
    theta_mat = np.array([list(c["tt"]) for c in cands], dtype=np.int32)
    out = {"N": args.N, "k_max": args.k_max, "K": len(cands)}
    if args.target in ("primes", "both"):
        p = primes_table(args.N)
        target_arr = np.array(p, dtype=np.int32)
        t0 = time.time()
        r = greedy_d2_np(target_arr, theta_mat, args.max_M)
        r["elapsed_s"] = time.time() - t0
        out["primes"] = r
        print(f"PRIMES: greedy M = {r['M']} (verified={r['verified']})  [{r['elapsed_s']:.1f}s]")
    if args.target in ("random", "both"):
        rt = random_table(args.N, sum(primes_table(args.N)), seed=args.seed)
        target_arr = np.array(rt, dtype=np.int32)
        t0 = time.time()
        r = greedy_d2_np(target_arr, theta_mat, args.max_M)
        r["elapsed_s"] = time.time() - t0
        out["random"] = r
        print(f"RANDOM: greedy M = {r['M']} (verified={r['verified']})  [{r['elapsed_s']:.1f}s]")
    def cleanup(o):
        if isinstance(o, dict): return {k: cleanup(v) for k, v in o.items()}
        if isinstance(o, list): return [cleanup(v) for v in o]
        if isinstance(o, np.integer): return int(o)
        if isinstance(o, np.ndarray): return o.tolist()
        try: json.dumps(o); return o
        except (TypeError, ValueError): return str(o)
    with open(f"greedy_d2_np_n{args.N}_k{args.k_max}.json", "w") as f:
        json.dump(cleanup(out), f, indent=2)
    print(f"Wrote greedy_d2_np_n{args.N}_k{args.k_max}.json")

if __name__ == "__main__":
    main()

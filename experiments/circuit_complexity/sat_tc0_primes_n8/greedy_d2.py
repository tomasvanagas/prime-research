"""
Greedy / column-generation heuristic for depth-2 sign-threshold.

Quickly find an UPPER BOUND on the min-M depth-2 circuit by:
1. Enumerate W=1 sign-threshold candidates.
2. Iteratively pick the candidate (and sign) that maximally improves the
   current LP relaxation of "is target a sign-weighted threshold of selected gates".
3. Stop when feasibility achieved.

This gives a quick upper bound on the depth-2 sign-threshold complexity,
complementing the exact ILP lower bound.
"""
from __future__ import annotations
import argparse, json, random, time
from typing import List, Tuple, Dict
import pulp

import sys; sys.path.insert(0, '.')
from enum_d2_smart import enumerate_w1_thresholds_pruned, primes_table, random_table

def feasibility_with_gates(target: List[int], N: int, gate_tts: List[Tuple[int,...]],
                           gate_alphas: List[int]) -> Tuple[bool, int]:
    """Given selected gates and signs, can we choose T_top integer to match target?"""
    # s_top(x) = sum_j alpha_j * theta_j(x)
    domain = list(range(2**N))
    s_top = [sum(a * tt[x] for a, tt in zip(gate_alphas, gate_tts)) for x in domain]
    # Need T_top such that target=1 ⇒ s_top >= T_top, target=0 ⇒ s_top <= T_top - 1
    s_max_neg = max((s_top[x] for x in domain if target[x] == 0), default=-10**9)
    s_min_pos = min((s_top[x] for x in domain if target[x] == 1), default=+10**9)
    if s_min_pos > s_max_neg:
        return True, s_max_neg + 1  # T_top = s_max_neg + 1 works (min T satisfying both)
    return False, 0

def greedy_d2(target: List[int], N: int, candidates: List[Dict],
              max_M: int = 100) -> Dict:
    """Greedy: at each step, pick the candidate (and sign) that minimizes
    the number of misclassified inputs."""
    domain = list(range(2**N))
    selected_tts = []
    selected_alphas = []
    selected_idxs = []
    for step in range(max_M):
        # Try each candidate with each sign
        best_improvement = -1
        best_choice = None
        for k, c in enumerate(candidates):
            for a in [-1, +1]:
                new_tts = selected_tts + [c["tt"]]
                new_alphas = selected_alphas + [a]
                feas, T_top = feasibility_with_gates(target, N, new_tts, new_alphas)
                if feas:
                    return {"M": step + 1, "T_top": T_top,
                            "selected_idxs": selected_idxs + [k],
                            "selected_alphas": new_alphas,
                            "verified": True}
                # Score: maximize the gap (s_min_pos - s_max_neg)
                s_top = [sum(aa * tt[x] for aa, tt in zip(new_alphas, new_tts)) for x in domain]
                s_max_neg = max((s_top[x] for x in domain if target[x] == 0), default=-10**9)
                s_min_pos = min((s_top[x] for x in domain if target[x] == 1), default=+10**9)
                gap = s_min_pos - s_max_neg  # >0 means feasible
                if gap > best_improvement:
                    best_improvement = gap
                    best_choice = (k, a, c)
        if best_choice is None:
            break
        k, a, c = best_choice
        selected_tts.append(c["tt"])
        selected_alphas.append(a)
        selected_idxs.append(k)
        # Status
        s_top = [sum(aa * tt[x] for aa, tt in zip(selected_alphas, selected_tts)) for x in domain]
        s_max_neg = max((s_top[x] for x in domain if target[x] == 0), default=-10**9)
        s_min_pos = min((s_top[x] for x in domain if target[x] == 1), default=+10**9)
        print(f"  step {step+1}: picked gate {k} alpha={a:+d}, gap = {s_min_pos - s_max_neg}")
    return {"M": max_M, "verified": False, "selected_idxs": selected_idxs,
            "selected_alphas": selected_alphas}

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--N", type=int, default=8)
    parser.add_argument("--k-max", type=int, default=5)
    parser.add_argument("--max-M", type=int, default=80)
    parser.add_argument("--target", choices=["primes", "random", "both"], default="both")
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()
    print(f"N={args.N}, k_max={args.k_max}")
    cands = enumerate_w1_thresholds_pruned(args.N, k_max=args.k_max)
    print(f"K = {len(cands)} candidates")
    out = {"N": args.N, "k_max": args.k_max, "K": len(cands)}
    if args.target in ("primes", "both"):
        print("\n== PRIMES ==")
        p = primes_table(args.N)
        t0 = time.time()
        out["primes"] = greedy_d2(p, args.N, cands, args.max_M)
        out["primes"]["elapsed_s"] = time.time() - t0
        print(f"PRIMES greedy upper bound: M = {out['primes']['M']}  ({out['primes']['elapsed_s']:.1f}s)")
    if args.target in ("random", "both"):
        print("\n== RANDOM ==")
        rt = random_table(args.N, sum(primes_table(args.N)), seed=args.seed)
        t0 = time.time()
        out["random"] = greedy_d2(rt, args.N, cands, args.max_M)
        out["random"]["elapsed_s"] = time.time() - t0
        print(f"RANDOM greedy upper bound: M = {out['random']['M']}  ({out['random']['elapsed_s']:.1f}s)")
    def cleanup(o):
        if isinstance(o, dict): return {k: cleanup(v) for k, v in o.items()}
        if isinstance(o, list): return [cleanup(v) for v in o]
        if isinstance(o, tuple): return list(o)
        try: json.dumps(o); return o
        except (TypeError, ValueError): return str(o)
    with open(f"greedy_d2_n{args.N}_k{args.k_max}.json", "w") as f:
        json.dump(cleanup(out), f, indent=2)
    print(f"\nWrote greedy_d2_n{args.N}_k{args.k_max}.json")

if __name__ == "__main__":
    main()

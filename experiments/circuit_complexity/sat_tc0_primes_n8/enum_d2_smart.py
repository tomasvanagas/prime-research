"""
Enumerated depth-2 sign-threshold search with smarter candidate generation.

Key idea: prune candidate bottom-layer threshold functions by:
- Limiting the number of nonzero weights (k_max).
- Removing constants and complements (alpha=±1 covers complements).
- For N=8 this brings K from 35K to 1-5K depending on k_max.

Provides:
- enumerate_w1_thresholds_pruned(N, k_max): enumerate W=1 sign-threshold
  with at most k_max nonzero weights.
- depth2_search(target, N, K_thetas, M_list): for each M, ILP-solve.
- Bottom-layer interpretability: per-gate decoding into "is x_a XOR x_b XOR ... ≥ T?"
"""
from __future__ import annotations
import argparse, json, random, time
from itertools import product, combinations
from typing import List, Tuple, Dict, Optional, Set
import pulp

def sieve(n):
    is_p = [False, False] + [True]*max(0, n-1)
    for i in range(2, int(n**0.5)+1):
        if is_p[i]:
            for j in range(i*i, n+1, i): is_p[j] = False
    return is_p

def primes_table(N): return [1 if sieve(2**N - 1)[k] else 0 for k in range(2**N)]

def random_table(N, w, seed):
    rng = random.Random(seed)
    t = [0]*(2**N); idx = list(range(2**N)); rng.shuffle(idx)
    for k in idx[:w]: t[k] = 1
    return t

def bits(x, N): return [(x>>i)&1 for i in range(N)]

def enumerate_w1_thresholds_pruned(N: int, k_max: int = None,
                                    drop_complements: bool = True) -> List[Dict]:
    """
    Returns list of dicts: {'weights': tuple, 'T': int, 'tt': tuple_truth_table}
    with at most k_max nonzero weights.
    """
    if k_max is None:
        k_max = N
    seen: Set[Tuple[int, ...]] = set()
    out = []
    domain_bits = [bits(x, N) for x in range(2**N)]
    for k in range(1, k_max + 1):
        for support in combinations(range(N), k):
            for sign_pattern in product([-1, 1], repeat=k):
                w = [0]*N
                for ii, idx in enumerate(support):
                    w[idx] = sign_pattern[ii]
                w = tuple(w)
                sums = [sum(w[i] * b[i] for i in range(N)) for b in domain_bits]
                s_min, s_max = min(sums), max(sums)
                for T in range(s_min, s_max + 2):
                    tt = tuple(1 if s >= T else 0 for s in sums)
                    if all(v == 0 for v in tt) or all(v == 1 for v in tt):
                        continue
                    if drop_complements:
                        # Canonical form: prefer the form with fewer 1s, or lex-smaller
                        ctt = tuple(1 - v for v in tt)
                        canon = min(tt, ctt)
                        if canon in seen:
                            continue
                        seen.add(canon)
                    else:
                        if tt in seen:
                            continue
                        seen.add(tt)
                    out.append({"w": w, "T": T, "tt": tt})
    return out

def depth2_search(target: List[int], N: int, candidates: List[Dict],
                  M_list: List[int], time_limit: int = 600,
                  greedy_warm_start: bool = False) -> Dict:
    """
    For each M in M_list, find smallest depth-2 sign-threshold representation
    using ≤ M of the candidate bottom-layer threshold functions.
    """
    K = len(candidates)
    domain = list(range(2**N))
    Theta = [c["tt"] for c in candidates]
    print(f"  Search over K={K} candidates, target weight = {sum(target)}")
    results = {"K": K}

    for M in M_list:
        print(f"  M={M:3d}: ", end="", flush=True)
        prob = pulp.LpProblem(f"d2_M{M}", pulp.LpMinimize)
        sel = [pulp.LpVariable(f"sel_{k}", cat='Binary') for k in range(K)]
        beta = [pulp.LpVariable(f"beta_{k}", cat='Binary') for k in range(K)]
        v = [pulp.LpVariable(f"v_{k}", cat='Binary') for k in range(K)]
        T_top = pulp.LpVariable("T_top", lowBound=-(M+2), upBound=M+2, cat='Integer')
        prob += pulp.lpSum(sel)  # minimize selected count (≤ M)
        prob += pulp.lpSum(sel) <= M
        for k in range(K):
            prob += v[k] <= sel[k]
            prob += v[k] <= beta[k]
            prob += v[k] >= sel[k] + beta[k] - 1
        for xi, x in enumerate(domain):
            ones_k = [k for k in range(K) if Theta[k][xi] == 1]
            s_top = 2 * pulp.lpSum(v[k] for k in ones_k) - pulp.lpSum(sel[k] for k in ones_k)
            if target[x] == 1:
                prob += s_top >= T_top, f"top1_{xi}"
            else:
                prob += s_top <= T_top - 1, f"top0_{xi}"
        solver = pulp.PULP_CBC_CMD(msg=0, timeLimit=time_limit, threads=4)
        t0 = time.time()
        status = prob.solve(solver)
        elapsed = time.time() - t0
        sname = pulp.LpStatus[status]
        if sname == "Optimal":
            chosen = [k for k in range(K) if sel[k].value() > 0.5]
            sols = []
            for k in chosen:
                a = 2*int(round(beta[k].value())) - 1
                sols.append({"theta_idx": k, "alpha": a,
                             "weights": list(candidates[k]["w"]),
                             "threshold": candidates[k]["T"]})
            T_top_v = int(round(T_top.value()))
            ok = True; n_correct = 0
            for xi in range(2**N):
                s = sum(s["alpha"] * Theta[s["theta_idx"]][xi] for s in sols)
                out = 1 if s >= T_top_v else 0
                if out == target[xi]: n_correct += 1
                else: ok = False
            results[M] = {"status": "sat", "elapsed_s": elapsed,
                         "chosen_count": len(chosen), "T_top": T_top_v,
                         "verified": ok, "correct": n_correct,
                         "sols": sols}
            print(f"sat ({elapsed:.1f}s) gates={len(chosen)} verified={ok}")
            if ok:
                results["min_M"] = len(chosen)
                results["sols_min"] = sols
                results["T_top_min"] = T_top_v
                return results
        elif sname == "Infeasible":
            results[M] = {"status": "unsat", "elapsed_s": elapsed}
            print(f"unsat ({elapsed:.1f}s)")
        else:
            results[M] = {"status": "unknown", "elapsed_s": elapsed, "lp_status": sname}
            print(f"unknown ({elapsed:.1f}s)")
    return results

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--N", type=int, default=8)
    parser.add_argument("--Mlist", type=str, default="6,8,10,12,16,20")
    parser.add_argument("--k-max", type=int, default=4)
    parser.add_argument("--time-limit", type=int, default=600)
    parser.add_argument("--target", choices=["primes", "random", "both"], default="both")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--out", type=str, default="enum_d2_smart_results.json")
    args = parser.parse_args()
    M_list = [int(x) for x in args.Mlist.split(",")]
    print(f"N={args.N}, k_max={args.k_max}")
    t0 = time.time()
    cands = enumerate_w1_thresholds_pruned(args.N, k_max=args.k_max)
    print(f"Enumerated {len(cands)} candidates in {time.time()-t0:.1f}s "
          f"(k_max={args.k_max}, dedup with complements)")
    out = {"N": args.N, "k_max": args.k_max, "K": len(cands), "M_list": M_list}
    p = primes_table(args.N); n_p = sum(p)
    print(f"PRIMES weight: {n_p}/{2**args.N}")
    if args.target in ("primes","both"):
        print("== PRIMES ==")
        out["primes"] = depth2_search(p, args.N, cands, M_list, args.time_limit)
    if args.target in ("random","both"):
        rt = random_table(args.N, n_p, seed=args.seed)
        print("== RANDOM ==")
        out["random"] = depth2_search(rt, args.N, cands, M_list, args.time_limit)
    def cleanup(o):
        if isinstance(o, dict): return {k: cleanup(v) for k, v in o.items()}
        if isinstance(o, list): return [cleanup(v) for v in o]
        if isinstance(o, tuple): return list(o)
        try: json.dumps(o); return o
        except (TypeError, ValueError): return str(o)
    with open(args.out, "w") as f:
        json.dump(cleanup(out), f, indent=2)
    print(f"\nWrote {args.out}")

if __name__ == "__main__":
    main()

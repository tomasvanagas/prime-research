"""
Pre-enumerate all distinct W=1 sign-threshold functions on N variables,
then ask: smallest M such that PRIMES is a sign-weighted threshold of
M of them?

For N=8, |W=1 threshold functions| ~ 50K which is feasible to enumerate.

Once enumerated as a 2^N x K Boolean matrix Theta, the depth-2 problem
becomes:

  min_{S, alpha, T_top}  |S|
  s.t.  for each x in [0..2^N-1]:
          (sum_{j in S} alpha_j * Theta[x,j] >= T_top)  iff  PRIMES(x) = 1
  with alpha_j in {-1, +1}, T_top integer.

We encode this as ILP and ask for feasibility at increasing M.

This is faster than the from-scratch search because the bottom-layer
parameters are eliminated; only column selection + top weights remain.
"""
from __future__ import annotations
import argparse, json, random, time
from itertools import product
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

def enumerate_w1_thresholds(N: int) -> List[Tuple[int, ...]]:
    """
    Enumerate all distinct Boolean functions f(x) = sign(sum_i w_i x_i - T)
    with w_i ∈ {-1, 0, +1}, T ∈ {-(N+1), ..., N+1}.

    Returns list of truth tables (as tuples of 0/1 length 2^N).

    De-duplicates: identical truth tables collapsed.
    Excludes constants (all-0, all-1).
    Includes both f and 1-f for completeness (we'll get both via signing).
    """
    seen: Set[Tuple[int, ...]] = set()
    domain_bits = [bits(x, N) for x in range(2**N)]
    for w in product([-1, 0, 1], repeat=N):
        if all(wi == 0 for wi in w):
            continue  # constant
        # Possible sums: for each x, sum(w_i * x_i)
        sums = [sum(w[i] * b[i] for i in range(N)) for b in domain_bits]
        s_min, s_max = min(sums), max(sums)
        # Thresholds that produce distinct functions: T from s_min to s_max + 1
        for T in range(s_min, s_max + 2):
            tt = tuple(1 if s >= T else 0 for s in sums)
            if not (all(v == 0 for v in tt) or all(v == 1 for v in tt)):
                seen.add(tt)
    return sorted(seen)

def depth2_min_M_via_enumeration(target: List[int], N: int, M_list: List[int],
                                  thetas: List[Tuple[int,...]], time_limit: int = 300,
                                  alpha_signed: bool = True) -> Dict:
    """For each M in M_list, determine if there's a depth-2 sign-threshold
    expression using ≤ M of the enumerated bottom gates."""
    K = len(thetas)
    print(f"  Enumerated K={K} W=1 sign-threshold functions on N={N}")
    domain = list(range(2**N))
    # Theta[x][k] = thetas[k][x]
    Theta = thetas  # list of tuples, Theta[k][x]
    results = {"K": K}
    for M in M_list:
        print(f"  M={M:3d}: ", end="", flush=True)
        prob = pulp.LpProblem(f"d2enum_M{M}", pulp.LpMinimize)
        # Use binary "gate selected" variable
        sel = [pulp.LpVariable(f"sel_{k}", cat='Binary') for k in range(K)]
        if alpha_signed:
            # alpha_k ∈ {-1, +1} when sel_k=1, else don't matter
            beta = [pulp.LpVariable(f"beta_{k}", cat='Binary') for k in range(K)]
            # alpha_k = 2*beta_k - 1
        T_top = pulp.LpVariable("T_top", lowBound=-(M+1), upBound=M+1, cat='Integer')
        # Limit total gates
        prob += pulp.lpSum(sel) <= M, "max_gates"
        # Objective: minimize total gates (just gives a feasible solution)
        prob += pulp.lpSum(sel)

        # For each input, contribution from gate k is:
        #   sel_k * alpha_k * Theta[k][x]
        # = sel_k * (2*beta_k - 1) * Theta[k][x]
        # = Theta[k][x] * (2 * (sel_k AND beta_k) - sel_k)
        # = 2*u_kx - sel_k * Theta[k][x]
        # where u_kx = sel_k AND beta_k AND Theta[k][x]  (Theta[k][x] is constant)
        # If Theta[k][x] = 0, contribution = 0; if Theta[k][x] = 1, contribution = 2*v_k - sel_k
        # where v_k = sel_k AND beta_k.

        # v_k = sel AND beta (binary):
        v = [pulp.LpVariable(f"v_{k}", cat='Binary') for k in range(K)]
        for k in range(K):
            prob += v[k] <= sel[k]
            prob += v[k] <= beta[k]
            prob += v[k] >= sel[k] + beta[k] - 1

        # s_top[x] = 2 * sum_{k: Theta[k][x]=1} v_k - sum_{k: Theta[k][x]=1} sel_k
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
                sols.append({"theta_idx": k, "alpha": a, "theta": list(Theta[k])})
            T_top_v = int(round(T_top.value()))
            # verify
            ok = True
            n_correct = 0
            for xi in range(2**N):
                s = sum(s["alpha"] * s["theta"][xi] for s in sols)
                out = 1 if s >= T_top_v else 0
                if out == target[xi]:
                    n_correct += 1
                else:
                    ok = False
            results[M] = {"status": "sat", "elapsed_s": elapsed,
                         "chosen_count": len(chosen), "T_top": T_top_v,
                         "verified": ok, "correct": n_correct,
                         "sample_thetas": sols[:5]}
            print(f"sat ({elapsed:.1f}s) gates={len(chosen)} verified={ok} correct={n_correct}/{2**N}")
            if ok:
                results["min_M"] = len(chosen)
                results["full_sols"] = sols
                results["full_T_top"] = T_top_v
                return results
        elif sname == "Infeasible":
            results[M] = {"status": "unsat", "elapsed_s": elapsed}
            print(f"unsat ({elapsed:.1f}s)")
        else:
            results[M] = {"status": "unknown", "elapsed_s": elapsed, "lp_status": sname}
            print(f"unknown:{sname} ({elapsed:.1f}s)")
    return results

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--N", type=int, default=8)
    parser.add_argument("--Mlist", type=str, default="4,6,8,10,12,16,20,24,32")
    parser.add_argument("--time-limit", type=int, default=300)
    parser.add_argument("--target", choices=["primes", "random", "both"], default="both")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--out", type=str, default="enumerated_depth2_results.json")
    args = parser.parse_args()
    M_list = [int(x) for x in args.Mlist.split(",")]
    print(f"Enumerating W=1 sign-threshold functions on N={args.N}")
    t0 = time.time()
    thetas = enumerate_w1_thresholds(args.N)
    print(f"Found {len(thetas)} distinct functions [{time.time()-t0:.1f}s]")
    out = {"N": args.N, "K": len(thetas), "M_list": M_list}
    p = primes_table(args.N); n_p = sum(p)
    print(f"PRIMES weight: {n_p}/{2**args.N}")
    if args.target in ("primes","both"):
        print("== PRIMES ==")
        out["primes"] = depth2_min_M_via_enumeration(p, args.N, M_list, thetas, args.time_limit)
    if args.target in ("random","both"):
        rt = random_table(args.N, n_p, seed=args.seed)
        print("== RANDOM ==")
        out["random"] = depth2_min_M_via_enumeration(rt, args.N, M_list, thetas, args.time_limit)

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

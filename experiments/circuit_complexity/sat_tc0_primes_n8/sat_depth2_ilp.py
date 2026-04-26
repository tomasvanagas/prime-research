"""
ILP-based depth-2 sign-threshold circuit search for PRIMES at N=8.

Cleaner / more scalable than the Z3 If-then-else encoding.  Uses PuLP +
CBC.  Variables and Big-M linearisation:

For each bottom gate j ∈ [M]:
  w[j][i] ∈ {-W..W} integer (input weight)
  T[j]    ∈ {-NW..NW} integer (threshold)
  theta[j][x] ∈ {0,1} for each input x ∈ [2^N]   (gate output on x)

Big-M linking:
  s[j][x] := sum_i w[j][i]*x_i           (just an expression, not var)
  s[j][x] - T[j] >= -B*(1 - theta[j][x])     # if theta=1 then s >= T
  s[j][x] - T[j] <= -1 + B*theta[j][x]       # if theta=0 then s <= T-1
  with B = 2*N*W + 2  (loose enough)

Top layer: alpha[j] ∈ {-1,+1} via beta[j] ∈ {0,1} (alpha = 2*beta - 1).
  v[j][x] := beta[j] AND theta[j][x]
  s_top[x] = 2 * sum_j v[j][x] - sum_j theta[j][x]
  Constraint: s_top[x] >= T_top  iff f(x)=1
  Implemented via:
    f(x)=1 ⇒ s_top[x] >= T_top
    f(x)=0 ⇒ s_top[x] <= T_top - 1

AND linearisation: v[j][x] = beta[j] AND theta[j][x]
  v <= beta;  v <= theta;  v >= beta + theta - 1.

Symmetry-breaking:
  T[0] <= T[1] <= ... <= T[M-1]  (canonical ordering)

This is intentionally a tighter formulation than the Z3 version — it
avoids If-then-else by encoding theta as Boolean directly.
"""
from __future__ import annotations
import argparse, json, random, time
from typing import List, Optional, Tuple, Dict
import pulp

def sieve(n: int) -> List[bool]:
    is_p = [False, False] + [True] * max(0, n - 1)
    for i in range(2, int(n**0.5) + 1):
        if is_p[i]:
            for j in range(i*i, n + 1, i):
                is_p[j] = False
    return is_p

def primes_table(N: int) -> List[int]:
    s = sieve(2**N - 1)
    return [1 if s[k] else 0 for k in range(2**N)]

def random_table(N: int, weight: int, seed: int) -> List[int]:
    rng = random.Random(seed)
    table = [0] * (2**N)
    idx = list(range(2**N))
    rng.shuffle(idx)
    for k in idx[:weight]:
        table[k] = 1
    return table

def bits(x: int, N: int) -> List[int]:
    return [(x >> i) & 1 for i in range(N)]

def encode_depth2(target: List[int], N: int, M: int, W: int = 1,
                  symmetry_break: bool = True) -> pulp.LpProblem:
    """Build the depth-2 sign-threshold ILP (alpha in {-1,+1})."""
    domain = list(range(2**N))
    bb = [bits(x, N) for x in domain]
    B = 2 * N * W + 2  # Big-M
    NW = N * W

    prob = pulp.LpProblem(f"d2_M{M}_W{W}", pulp.LpMinimize)
    prob += 0  # feasibility, no objective

    # Bottom-layer weights
    w = [[pulp.LpVariable(f"w_{j}_{i}", lowBound=-W, upBound=W, cat='Integer')
          for i in range(N)] for j in range(M)]
    T = [pulp.LpVariable(f"T_{j}", lowBound=-NW, upBound=NW + 1, cat='Integer')
         for j in range(M)]

    # Bottom-layer outputs (Boolean per input)
    theta = [[pulp.LpVariable(f"theta_{j}_{x}", cat='Binary')
              for x in domain] for j in range(M)]

    # Top-layer alpha encoded as beta ∈ {0,1}, alpha = 2*beta - 1
    beta = [pulp.LpVariable(f"beta_{j}", cat='Binary') for j in range(M)]
    # Top-layer threshold
    T_top = pulp.LpVariable("T_top", lowBound=-(M+1), upBound=M+1, cat='Integer')

    # AND linearisation
    v = [[pulp.LpVariable(f"v_{j}_{x}", cat='Binary')
          for x in domain] for j in range(M)]

    # Big-M linking for theta
    for j in range(M):
        for xi, x in enumerate(domain):
            s_jx = pulp.lpSum(w[j][i] * bb[xi][i] for i in range(N))
            # if theta=1 then s_jx >= T[j]
            prob += s_jx - T[j] >= -B * (1 - theta[j][xi]), f"th1_{j}_{xi}"
            # if theta=0 then s_jx <= T[j] - 1
            prob += s_jx - T[j] <= -1 + B * theta[j][xi], f"th0_{j}_{xi}"

    # AND constraints: v[j][x] = beta[j] AND theta[j][x]
    for j in range(M):
        for xi, x in enumerate(domain):
            prob += v[j][xi] <= beta[j], f"v1_{j}_{xi}"
            prob += v[j][xi] <= theta[j][xi], f"v2_{j}_{xi}"
            prob += v[j][xi] >= beta[j] + theta[j][xi] - 1, f"v3_{j}_{xi}"

    # Top-layer output: s_top[x] = 2*sum_j v[j][x] - sum_j theta[j][x]
    for xi, x in enumerate(domain):
        s_top_x = 2 * pulp.lpSum(v[j][xi] for j in range(M)) \
                  - pulp.lpSum(theta[j][xi] for j in range(M))
        if target[x] == 1:
            prob += s_top_x >= T_top, f"top1_{xi}"
        else:
            prob += s_top_x <= T_top - 1, f"top0_{xi}"

    # Symmetry breaking
    if symmetry_break:
        for j in range(M - 1):
            prob += T[j] <= T[j+1], f"sb_T_{j}"

    return prob, {"w": w, "T": T, "theta": theta, "beta": beta, "v": v, "T_top": T_top}

def solve_one(target: List[int], N: int, M: int, W: int = 1,
              time_limit: int = 300) -> Tuple[str, Optional[Dict]]:
    prob, vars_ = encode_depth2(target, N, M, W)
    solver = pulp.PULP_CBC_CMD(msg=0, timeLimit=time_limit, threads=4)
    t0 = time.time()
    status_int = prob.solve(solver)
    elapsed = time.time() - t0
    status_name = pulp.LpStatus[status_int]
    out = {"elapsed_s": elapsed, "status_name": status_name}

    if status_name == "Optimal":
        # Feasibility found
        try:
            T_top_v = int(round(vars_["T_top"].value()))
            T_v = [int(round(vars_["T"][j].value())) for j in range(M)]
            w_v = [[int(round(vars_["w"][j][i].value())) for i in range(N)] for j in range(M)]
            beta_v = [int(round(vars_["beta"][j].value())) for j in range(M)]
            alpha_v = [2*b - 1 for b in beta_v]
            out.update({"T_top": T_top_v, "T": T_v, "w": w_v, "alpha": alpha_v})
            return "sat", out
        except Exception as e:
            out["err"] = str(e)
            return "unknown", out
    elif status_name == "Infeasible":
        return "unsat", out
    elif status_name == "Not Solved":
        return "unknown", out
    else:
        return "unknown", out

def verify(model: Dict, target: List[int], N: int) -> Tuple[bool, int]:
    w = model["w"]; T = model["T"]; alpha = model["alpha"]; T_top = model["T_top"]
    M = len(T)
    correct = 0
    for x in range(2**N):
        bb = bits(x, N)
        theta = [int(sum(w[j][i] * bb[i] for i in range(N)) >= T[j]) for j in range(M)]
        s_top = sum(alpha[j] * theta[j] for j in range(M))
        out = 1 if s_top >= T_top else 0
        if out == target[x]:
            correct += 1
    return correct == 2**N, correct

def scan(target: List[int], N: int, M_list: List[int], W: int = 1,
         time_limit: int = 300, label: str = "primes") -> Dict:
    results = {}
    for M in M_list:
        print(f"  [{label}] N={N} W={W} M={M:3d}: ", end="", flush=True)
        status, model = solve_one(target, N, M, W=W, time_limit=time_limit)
        t = model["elapsed_s"]
        results[M] = {"status": status, "elapsed_s": t}
        if status == "sat":
            ok, ncorr = verify(model, target, N)
            results[M]["verified"] = ok
            results[M]["model"] = model
            print(f"{status:8s} ({t:6.1f}s) verified={ok} correct={ncorr}/{2**N}")
            if ok:
                # Found a solution; could keep going but typically stop
                pass
        else:
            print(f"{status:8s} ({t:6.1f}s)")
    return results

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--N", type=int, default=8)
    parser.add_argument("--Mlist", type=str, default="4,8,12,16,24,32,48")
    parser.add_argument("--W", type=int, default=1)
    parser.add_argument("--time-limit", type=int, default=300)
    parser.add_argument("--target", choices=["primes", "random", "both"], default="both")
    parser.add_argument("--out", type=str, default="sat_depth2_ilp_results.json")
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()
    M_list = [int(x) for x in args.Mlist.split(",")]
    print(f"N={args.N}  M_list={M_list}  W={args.W}  time_limit={args.time_limit}s")
    p = primes_table(args.N)
    n_p = sum(p)
    print(f"PRIMES weight: {n_p}/{2**args.N}")
    out = {"N": args.N, "W": args.W, "M_list": M_list, "n_primes": n_p}
    if args.target in ("primes", "both"):
        out["primes"] = scan(p, args.N, M_list, W=args.W, time_limit=args.time_limit, label="primes")
    if args.target in ("random", "both"):
        rt = random_table(args.N, n_p, seed=args.seed)
        out["random"] = scan(rt, args.N, M_list, W=args.W, time_limit=args.time_limit, label="random")

    def cleanup(o):
        if isinstance(o, dict):
            return {k: cleanup(v) for k, v in o.items()}
        if isinstance(o, list):
            return [cleanup(v) for v in o]
        try:
            json.dumps(o); return o
        except (TypeError, ValueError):
            return str(o)
    with open(args.out, "w") as f:
        json.dump(cleanup(out), f, indent=2)
    print(f"\nWrote {args.out}")

if __name__ == "__main__":
    main()

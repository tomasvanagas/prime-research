"""
SAT/ILP-based search for small TC^0 circuits computing PRIMES at N=8.

Frontier attack §A1 from ATTACK_VECTORS.md: find a primality witness using
only NC^1-known scalar primitives, by enumerating small TC^0 circuits via
SAT and seeing what falls out.

The proper TC^0 basis is THRESHOLD GATES (linear thresholds over Boolean
inputs), not {AND, OR, XOR}.  Prior work in this project (S20, S28,
sat_circuit_minimization.py) measured min Boolean-circuit size in the
{AND, OR, XOR} basis and found pi(x) mod 2 ~ 2^(0.73N) BDD-complex.  This
script instead asks: what is the SHALLOWEST and SMALLEST threshold-gate
circuit that EXACTLY computes the prime indicator on all 256 inputs?

Approach
--------
1. PTF degree.  Compute the exact polynomial-threshold-function degree
   of PRIMES at N=8 via LP feasibility on the monomial basis.  This is a
   classical lower-bound technique (Aspnes-Beigel-Furst-Rudich, Sherstov)
   on threshold circuit complexity.

2. Depth-1 sign-threshold.  LP test whether PRIMES is itself a (signed,
   integer-weight, possibly large-weight) linear threshold function.

3. Depth-2 sign-threshold via Z3.  Encode "exists M-bottom-gate, depth-2
   sign-threshold circuit matching PRIMES on N=8" with bounded weights W,
   for increasing M.  Find smallest feasible M for each W.

4. Comparison.  Run same protocol on (a) random Boolean functions with
   matched density and (b) is_prime (no carry / no pi(x) cumulative
   structure).  Establish whether PRIMES is harder than random of the
   same weight.

Why this might surprise us:
- If small depth-2 circuits exist with bounded weights, that's a
  partial positive A-grade result: PRIMES has compact non-trivial
  TC^0 representation at small N.  Doesn't immediately give polylog
  PRIMES, but contradicts the implicit working assumption.
- If even huge bounded-weight depth-2 circuits are infeasible, that's
  a quantitative negative bound complementing the prior measurements.

Cross-domain ingredient: SAT/ILP solver + formal threshold-circuit
synthesis, as in Niemetz et al. and the SAT-based logic synthesis
community.  Has not been attempted in this project on threshold
gates (only on {AND, OR, XOR}).

Failure profile:
- Z3 times out at small M; we report bounds + scaling, B-grade.
- All M up to budget infeasible at W=1, but trivially feasible at
  large W (Muroga); we report W-vs-M tradeoff curve, B-grade.
- We find a small circuit (M ~ 10, W ~ 1); that would be A-grade.

Usage:
    python sat_tc0_primes_n8.py [--phase ptf|depth1|depth2|all] [--N 8]
                                 [--Mmax 32] [--Wmax 4] [--timeout 600]
"""
from __future__ import annotations
import argparse, json, math, random, time, sys
from typing import List, Dict, Optional, Tuple, Callable

import z3
import pulp

# -- target functions -----------------------------------------------------

def sieve(n: int) -> List[bool]:
    is_p = [False, False] + [True] * (n - 1)
    for i in range(2, int(n**0.5) + 1):
        if is_p[i]:
            for j in range(i*i, n + 1, i):
                is_p[j] = False
    return is_p

def is_prime_table(N: int) -> List[int]:
    s = sieve(2**N - 1 if 2**N >= 2 else 2)
    return [1 if (k <= 2**N - 1 and s[k]) else 0 for k in range(2**N)]

def pi_mod2_table(N: int) -> List[int]:
    """parity of #{p prime : p <= x}"""
    s = sieve(2**N - 1)
    out = [0] * (2**N)
    cnt = 0
    for k in range(2**N):
        if k <= 2**N - 1 and s[k]:
            cnt += 1
        out[k] = cnt & 1
    return out

def random_table(N: int, weight: int, seed: int = 12345) -> List[int]:
    rng = random.Random(seed)
    table = [0] * (2**N)
    indices = list(range(2**N))
    rng.shuffle(indices)
    for k in indices[:weight]:
        table[k] = 1
    return table

def bits_of(x: int, N: int) -> List[int]:
    return [(x >> i) & 1 for i in range(N)]


# -- 1. PTF degree via LP -------------------------------------------------

def ptf_degree_lp(table: List[int], N: int, max_degree: int = 8) -> int:
    """
    Find smallest d such that there's a degree-d polynomial p with
    sign(p(x)) = 2*table[x] - 1 (i.e., p(x) > 0 iff table[x]=1).

    Uses LP feasibility: minimize a slack epsilon, want epsilon > 0.
    p(x) = sum_{|S|<=d} c_S * prod_{i in S} x_i

    Returns smallest d for which LP is feasible (modulo numerical eps);
    returns -1 if even d=N is infeasible (shouldn't happen).
    """
    from itertools import combinations
    domain = list(range(2**N))
    bits = [bits_of(x, N) for x in domain]
    targets = [2*table[x] - 1 for x in domain]  # +/- 1

    for d in range(0, max_degree + 1):
        # monomials of degree <= d
        monos = []
        for k in range(d + 1):
            for S in combinations(range(N), k):
                monos.append(S)
        # LP: maximize epsilon s.t.
        #   targets[x] * sum_S c_S * prod_{i in S} bits[x][i] >= epsilon for all x
        #   |c_S| <= 1 (to bound)
        prob = pulp.LpProblem(f"ptf_d{d}", pulp.LpMaximize)
        eps = pulp.LpVariable("eps", lowBound=-100, upBound=100)
        c = {S: pulp.LpVariable(f"c_{'_'.join(map(str, S))}", lowBound=-1, upBound=1)
             for S in monos}
        prob += eps
        for xi, x in enumerate(domain):
            mono_vals = []
            for S in monos:
                v = 1
                for i in S:
                    v &= bits[xi][i]
                mono_vals.append(v)
            poly = pulp.lpSum(c[S] * mv for S, mv in zip(monos, mono_vals))
            prob += targets[xi] * poly >= eps, f"con_{xi}"
        solver = pulp.PULP_CBC_CMD(msg=0, timeLimit=120)
        status = prob.solve(solver)
        if status == 1 and eps.value() > 1e-6:
            return d
    return -1


# -- 2. Depth-1 sign-threshold via LP ------------------------------------

def is_threshold_function(table: List[int], N: int, weight_bound: Optional[int] = None) -> Tuple[bool, Optional[Dict]]:
    """
    Test if table[x] = 1 iff sum w_i x_i >= T, for some real (or integer
    if weight_bound provided) weights w_1..w_N and threshold T.

    Returns (feasible, {weights, threshold}).
    """
    from itertools import product
    domain = list(range(2**N))
    bits = [bits_of(x, N) for x in domain]
    prob = pulp.LpProblem("threshold", pulp.LpMaximize)
    if weight_bound is None:
        w = [pulp.LpVariable(f"w_{i}", lowBound=-100, upBound=100) for i in range(N)]
        T = pulp.LpVariable("T", lowBound=-100, upBound=100)
    else:
        w = [pulp.LpVariable(f"w_{i}", lowBound=-weight_bound, upBound=weight_bound, cat='Integer') for i in range(N)]
        T = pulp.LpVariable("T", lowBound=-N*weight_bound, upBound=N*weight_bound, cat='Integer')
    eps = pulp.LpVariable("eps", lowBound=-10, upBound=10)
    prob += eps
    for xi, x in enumerate(domain):
        s = pulp.lpSum(w[i] * bits[xi][i] for i in range(N)) - T
        if table[x] == 1:
            prob += s >= eps, f"pos_{xi}"
        else:
            prob += -s >= eps, f"neg_{xi}"
    solver = pulp.PULP_CBC_CMD(msg=0, timeLimit=60)
    status = prob.solve(solver)
    if status == 1 and eps.value() > 1e-6:
        return True, {"w": [w[i].value() for i in range(N)], "T": T.value()}
    return False, None


# -- 3. Depth-2 sign-threshold via Z3 ------------------------------------

def search_depth2(table: List[int], N: int, M: int, W: int = 1, W_top: int = 1,
                  timeout_s: int = 600) -> Tuple[str, Optional[Dict]]:
    """
    Search for a depth-2 sign-threshold circuit with:
      - M bottom-layer gates, each with weights in [-W..W] and integer threshold
      - 1 top-layer gate, weights in [-W_top..W_top], integer threshold.

    Returns (status, model_dict) where status in {'sat', 'unsat', 'unknown'}.
    """
    s = z3.Solver()
    s.set("timeout", timeout_s * 1000)

    # Variables
    w = [[z3.Int(f"w_{j}_{i}") for i in range(N)] for j in range(M)]
    T = [z3.Int(f"T_{j}") for j in range(M)]
    alpha = [z3.Int(f"alpha_{j}") for j in range(M)]
    T_top = z3.Int("T_top")

    # Domain
    for j in range(M):
        for i in range(N):
            s.add(w[j][i] >= -W, w[j][i] <= W)
        s.add(T[j] >= -N*W, T[j] <= N*W + 1)
        s.add(alpha[j] >= -W_top, alpha[j] <= W_top)
    s.add(T_top >= -M*W_top, T_top <= M*W_top + 1)

    # Symmetry breaking: order bottom gates lexicographically (cheap canonical form)
    # Sort by lex order of (T_j, w_j_0, w_j_1, ...).
    for j in range(M - 1):
        # Strict lex tie-break: enforce T_j <= T_{j+1}
        s.add(T[j] <= T[j+1])

    # For each input, one Boolean var per bottom gate output
    for x_int in range(2**N):
        bits = bits_of(x_int, N)
        sums = [z3.Sum([w[j][i] * bits[i] for i in range(N)]) for j in range(M)]
        # Bottom-layer outputs as If-then-else
        theta = [z3.If(sums[j] >= T[j], 1, 0) for j in range(M)]
        # Top-layer sum
        s_top = z3.Sum([alpha[j] * theta[j] for j in range(M)])
        if table[x_int] == 1:
            s.add(s_top >= T_top)
        else:
            s.add(s_top <= T_top - 1)

    t0 = time.time()
    res = s.check()
    elapsed = time.time() - t0

    if res == z3.sat:
        m = s.model()
        return "sat", {
            "w": [[m[w[j][i]].as_long() for i in range(N)] for j in range(M)],
            "T": [m[T[j]].as_long() for j in range(M)],
            "alpha": [m[alpha[j]].as_long() for j in range(M)],
            "T_top": m[T_top].as_long(),
            "elapsed_s": elapsed,
        }
    elif res == z3.unsat:
        return "unsat", {"elapsed_s": elapsed}
    else:
        return "unknown", {"elapsed_s": elapsed, "reason": str(s.reason_unknown())}


def verify_depth2(model: Dict, table: List[int], N: int) -> Tuple[bool, int]:
    """Verify a depth-2 solution against the truth table.  Returns (correct, n_correct)."""
    w = model["w"]; T = model["T"]; alpha = model["alpha"]; T_top = model["T_top"]
    M = len(T)
    correct = 0
    bad = []
    for x in range(2**N):
        bits = bits_of(x, N)
        theta = [int(sum(w[j][i] * bits[i] for i in range(N)) >= T[j]) for j in range(M)]
        s_top = sum(alpha[j] * theta[j] for j in range(M))
        out = 1 if s_top >= T_top else 0
        if out == table[x]:
            correct += 1
        else:
            bad.append((x, out, table[x]))
    return correct == 2**N, correct


# -- 4. Driver ------------------------------------------------------------

def run_all(N: int, Mmax: int, Wmax: int, timeout_s: int):
    print(f"\n{'='*68}\nSAT TC^0 search for PRIMES at N={N}\n{'='*68}\n")
    primes_table = is_prime_table(N)
    n_primes = sum(primes_table)
    print(f"# primes in [0, {2**N - 1}] = {n_primes}")

    # Phase 1: PTF degree
    print(f"\n[1] PTF degree of PRIMES at N={N}")
    d = ptf_degree_lp(primes_table, N, max_degree=N)
    print(f"    PTF degree (real-weight, monomial basis) = {d}")

    # Random and is_prime control: PTF degrees
    print(f"\n    Control: PTF degree of random function with same weight")
    rand = random_table(N, n_primes)
    d_rand = ptf_degree_lp(rand, N, max_degree=N)
    print(f"    PTF degree (random, weight={n_primes}) = {d_rand}")

    # Phase 2: Depth-1 sign-threshold
    print(f"\n[2] Depth-1 (single threshold gate) for PRIMES at N={N}")
    feas, sol = is_threshold_function(primes_table, N, weight_bound=None)
    print(f"    Real-weight feasible: {feas}")
    if feas:
        print(f"    weights: {sol['w']}, T: {sol['T']}")
    feas_int, _ = is_threshold_function(primes_table, N, weight_bound=2**N)
    print(f"    Integer-weight (|w|<=2^N) feasible: {feas_int}")

    # Phase 3: Depth-2 sign-threshold scan
    print(f"\n[3] Depth-2 sign-threshold scan: PRIMES at N={N}")
    print(f"    Vary M (bottom gates) in 1..{Mmax} and W in 1..{Wmax}")
    results = {"primes": {}, "random": {}}

    target_funcs = [("primes", primes_table)]
    if N <= 8:
        rand_table_local = random_table(N, n_primes, seed=42)
        target_funcs.append(("random", rand_table_local))

    for tname, ttable in target_funcs:
        results[tname] = {}
        for W in range(1, Wmax + 1):
            results[tname][W] = {}
            min_M_sat = None
            for M in [1, 2, 3, 4, 6, 8, 12, 16, 20, 24, 32, 48, 64][:None]:
                if M > Mmax:
                    break
                print(f"    [{tname}] W={W} M={M:3d}: ", end="", flush=True)
                status, model = search_depth2(ttable, N, M, W=W, W_top=1, timeout_s=timeout_s)
                t = model.get("elapsed_s", 0) if model else 0
                print(f"{status:8s} ({t:6.1f}s)")
                results[tname][W][M] = {"status": status, "elapsed_s": t}
                if status == "sat":
                    ok, ncorr = verify_depth2(model, ttable, N)
                    results[tname][W][M]["verified"] = ok
                    results[tname][W][M]["model"] = model
                    if ok:
                        min_M_sat = M
                        break
                if status == "unknown" and t > 0.9 * timeout_s:
                    # giving up at this W
                    break
            results[tname][W]["min_M_sat"] = min_M_sat
            print(f"    [{tname}] W={W}: min M (depth-2 sat) = {min_M_sat}")

    return {
        "N": N,
        "n_primes": n_primes,
        "ptf_degree_primes": d,
        "ptf_degree_random": d_rand,
        "depth1_real_feasible": feas,
        "depth1_int_feasible": feas_int,
        "depth2": results,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--N", type=int, default=8)
    parser.add_argument("--Mmax", type=int, default=20)
    parser.add_argument("--Wmax", type=int, default=2)
    parser.add_argument("--timeout", type=int, default=300)
    parser.add_argument("--phase", choices=["ptf", "depth1", "depth2", "all"], default="all")
    parser.add_argument("--out", type=str, default="sat_tc0_primes_n8_results.json")
    args = parser.parse_args()

    if args.phase == "all":
        out = run_all(args.N, args.Mmax, args.Wmax, args.timeout)
    elif args.phase == "ptf":
        primes = is_prime_table(args.N)
        d = ptf_degree_lp(primes, args.N, max_degree=args.N)
        out = {"N": args.N, "ptf_degree": d}
        print(f"PTF degree of PRIMES at N={args.N} = {d}")
    elif args.phase == "depth1":
        primes = is_prime_table(args.N)
        feas, sol = is_threshold_function(primes, args.N)
        out = {"N": args.N, "depth1_feasible": feas, "solution": sol}
        print(f"Depth-1 feasible: {feas}, sol = {sol}")
    elif args.phase == "depth2":
        primes = is_prime_table(args.N)
        for M in [1, 2, 4, 8, 16]:
            if M > args.Mmax: break
            status, model = search_depth2(primes, args.N, M, W=1, W_top=1, timeout_s=args.timeout)
            print(f"M={M}: {status}")
            if status == "sat":
                ok, ncorr = verify_depth2(model, primes, args.N)
                print(f"  verified: {ok}, correct: {ncorr}/{2**args.N}")
                break
        out = {"N": args.N, "result": status}

    with open(args.out, "w") as f:
        # Filter out non-serialisable model objects
        def cleanup(o):
            if isinstance(o, dict):
                return {k: cleanup(v) for k, v in o.items()}
            if isinstance(o, list):
                return [cleanup(v) for v in o]
            try:
                json.dumps(o)
                return o
            except (TypeError, ValueError):
                return str(o)
        json.dump(cleanup(out), f, indent=2)
    print(f"\nWrote results to {args.out}")


if __name__ == "__main__":
    main()

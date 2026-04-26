"""
Test depth-1 threshold representability of PRIMES at N=4..8 with various
weight bounds.

Returns: depth-1 sign-threshold complexity of PRIMES is INFINITY (it's
not a threshold function at all) — i.e. infeasible even with unbounded
real weights.
"""
from __future__ import annotations
import json, time
from typing import List, Optional, Tuple, Dict
import pulp, random

def sieve(n: int) -> List[bool]:
    is_p = [False, False] + [True] * max(0, n - 1)
    for i in range(2, int(n**0.5) + 1):
        if is_p[i]:
            for j in range(i*i, n+1, i):
                is_p[j] = False
    return is_p

def primes_table(N: int) -> List[int]:
    s = sieve(2**N - 1)
    return [1 if s[k] else 0 for k in range(2**N)]

def random_table(N: int, w: int, seed: int) -> List[int]:
    rng = random.Random(seed)
    t = [0]*(2**N); idx = list(range(2**N)); rng.shuffle(idx)
    for k in idx[:w]: t[k] = 1
    return t

def bits(x, N): return [(x>>i)&1 for i in range(N)]

def threshold_lp(table: List[int], N: int) -> Tuple[bool, Optional[Dict], float]:
    """LP for real-weight depth-1 threshold."""
    domain = list(range(2**N))
    bb = [bits(x,N) for x in domain]
    prob = pulp.LpProblem("th", pulp.LpMaximize)
    w = [pulp.LpVariable(f"w_{i}", lowBound=-1000, upBound=1000) for i in range(N)]
    T = pulp.LpVariable("T", lowBound=-1000, upBound=1000)
    eps = pulp.LpVariable("eps", lowBound=-10, upBound=10)
    prob += eps
    for xi, x in enumerate(domain):
        s = pulp.lpSum(w[i] * bb[xi][i] for i in range(N)) - T
        if table[x] == 1:
            prob += s >= eps
        else:
            prob += -s >= eps + 1e-6
    solver = pulp.PULP_CBC_CMD(msg=0, timeLimit=60)
    t0 = time.time()
    status = prob.solve(solver)
    elapsed = time.time() - t0
    feas = (status == 1) and (eps.value() is not None) and (eps.value() > 1e-6)
    sol = None
    if feas:
        sol = {"w": [w[i].value() for i in range(N)], "T": T.value(), "eps": eps.value()}
    return feas, sol, elapsed

def main():
    out = {}
    for N in [4,5,6,7,8]:
        p = primes_table(N)
        feas, sol, t = threshold_lp(p, N)
        wp = sum(p)
        print(f"N={N}: PRIMES (weight {wp}/{2**N}) depth-1 threshold = {feas}  [{t:.2f}s]")
        out[N] = {"primes_d1_feasible": feas, "n_primes": wp, "time_s": t}
        # Also try random matched
        rd1 = []
        for seed in [11,22,33,44,55]:
            r = random_table(N, wp, seed)
            f, _, _ = threshold_lp(r, N)
            rd1.append(f)
        out[N]["random_d1_feasibility"] = rd1
        print(f"        random matched (5 seeds): {rd1.count(True)}/5 feasible")
    with open("depth1_threshold_results.json", "w") as f:
        json.dump(out, f, indent=2)
    print("\nWrote depth1_threshold_results.json")

if __name__ == "__main__":
    main()

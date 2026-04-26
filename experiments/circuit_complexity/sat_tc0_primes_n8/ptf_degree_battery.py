"""
PTF degree of PRIMES vs random controls at N=4..8.

PTF degree d_pm(f) is a known LOWER BOUND on the SIZE of any depth-2
threshold circuit computing f, via Aspnes-Beigel-Furst-Rudich (1994):
size(f) >= 2^(Omega(d_pm(f))) for unweighted threshold-of-thresholds.

Sherstov 2008 sharpened this for sign-rank vs PTF degree.

If PTF(PRIMES at N) = PTF(random with same weight at N), that's strong
evidence PRIMES has no special TC^0 structure beyond what density forces.

If PTF(PRIMES) < PTF(random), then PRIMES has special structure
(unlikely given the project's pseudorandomness measurements).

If PTF(PRIMES) > PTF(random), genuine HARDNESS evidence.
"""
from __future__ import annotations
import json, random, time
from itertools import combinations
from typing import List, Optional
import pulp

def sieve(n: int) -> List[bool]:
    is_p = [False, False] + [True] * max(0, n - 1)
    for i in range(2, int(n**0.5) + 1):
        if is_p[i]:
            for j in range(i*i, n + 1, i):
                is_p[j] = False
    return is_p

def is_prime_table(N: int) -> List[int]:
    s = sieve(2**N - 1 if 2**N >= 2 else 2)
    return [1 if (k <= 2**N - 1 and s[k]) else 0 for k in range(2**N)]

def pi_mod2_table(N: int) -> List[int]:
    s = sieve(2**N - 1)
    out = [0] * (2**N)
    cnt = 0
    for k in range(2**N):
        if k <= 2**N - 1 and s[k]:
            cnt += 1
        out[k] = cnt & 1
    return out

def random_table(N: int, weight: int, seed: int) -> List[int]:
    rng = random.Random(seed)
    table = [0] * (2**N)
    indices = list(range(2**N))
    rng.shuffle(indices)
    for k in indices[:weight]:
        table[k] = 1
    return table

def bits_of(x: int, N: int) -> List[int]:
    return [(x >> i) & 1 for i in range(N)]

def ptf_feasible(table: List[int], N: int, d: int, time_limit: int = 120) -> Optional[bool]:
    """Test if there's a degree-d real-coeff polynomial p s.t. sign(p(x)) = 2*table[x]-1."""
    domain = list(range(2**N))
    bits = [bits_of(x, N) for x in domain]
    targets = [2*table[x] - 1 for x in domain]
    monos = []
    for k in range(d + 1):
        for S in combinations(range(N), k):
            monos.append(S)
    prob = pulp.LpProblem("ptf", pulp.LpMaximize)
    eps = pulp.LpVariable("eps", lowBound=-100, upBound=100)
    c = {S: pulp.LpVariable(f"c_{'_'.join(map(str, S)) or 'const'}", lowBound=-1, upBound=1)
         for S in monos}
    prob += eps
    # precompute monomial values per (x, S)
    for xi in range(len(domain)):
        terms = []
        for S in monos:
            v = 1
            for i in S:
                v &= bits[xi][i]
            if v:
                terms.append(c[S])
        poly = pulp.lpSum(terms)
        prob += targets[xi] * poly >= eps
    solver = pulp.PULP_CBC_CMD(msg=0, timeLimit=time_limit)
    status = prob.solve(solver)
    if status != 1:
        return None  # solver failed / timeout
    return eps.value() is not None and eps.value() > 1e-6

def find_ptf_degree(table: List[int], N: int, max_d: Optional[int] = None) -> int:
    if max_d is None:
        max_d = N
    for d in range(0, max_d + 1):
        feas = ptf_feasible(table, N, d)
        if feas is True:
            return d
        if feas is None:
            return -1  # solver timed out
    return -1

def main():
    out = {}
    for N in [4, 5, 6, 7, 8]:
        print(f"\n=== N={N} ===")
        t0 = time.time()
        primes = is_prime_table(N)
        wp = sum(primes)
        d_p = find_ptf_degree(primes, N)
        print(f"  PRIMES (weight {wp}/{2**N}): PTF degree = {d_p}  [{time.time()-t0:.1f}s]")

        t0 = time.time()
        pi2 = pi_mod2_table(N)
        wpi2 = sum(pi2)
        d_pi2 = find_ptf_degree(pi2, N)
        print(f"  pi(x) mod 2 (weight {wpi2}/{2**N}): PTF degree = {d_pi2}  [{time.time()-t0:.1f}s]")

        # 5 random controls
        d_rands = []
        for seed in [11, 22, 33, 44, 55]:
            t0 = time.time()
            r = random_table(N, wp, seed)
            d_r = find_ptf_degree(r, N)
            d_rands.append(d_r)
            print(f"  RANDOM (weight {wp}, seed {seed}): PTF degree = {d_r}  [{time.time()-t0:.1f}s]")

        out[N] = {
            "n_primes": wp,
            "ptf_degree_primes": d_p,
            "ptf_degree_pi_mod2": d_pi2,
            "ptf_degree_random": d_rands,
            "n_inputs": 2**N,
        }

    with open("ptf_degree_results.json", "w") as f:
        json.dump(out, f, indent=2)
    print("\nWrote ptf_degree_results.json")
    print("\nSummary:")
    print(f"{'N':>3} {'PRIMES':>8} {'pi mod 2':>10} {'random med':>12}")
    for N in sorted(out):
        rs = sorted(out[N]['ptf_degree_random'])
        med = rs[len(rs)//2]
        print(f"{N:>3} {out[N]['ptf_degree_primes']:>8} {out[N]['ptf_degree_pi_mod2']:>10} {med:>12}")

if __name__ == "__main__":
    main()

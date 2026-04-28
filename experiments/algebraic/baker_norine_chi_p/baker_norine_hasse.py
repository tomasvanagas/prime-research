"""
Test Baker-Norine RR on the Hasse diagram of divisibility — sparser graph
where primes are NOT all leaves.

Hasse(N) := ([1, N], edges (a, b) where b = p*a for some prime p, a < b).

This graph has |E| = sum_{p prime} floor(N/p) ~ N log log N (Mertens).
"""

import json
import time
import argparse
import random
import statistics
from collections import deque
from sympy import isprime, primerange

from baker_norine_chi_p import (
    graph_stats, is_connected, q_reduce, is_winnable,
    prime_divisor, random_matched_divisor, fire_vertex
)


def hasse_graph(N):
    """V = [1, N], edges (a, p*a) for prime p with p*a <= N."""
    V = list(range(1, N + 1))
    adj = {v: [] for v in V}
    primes = list(primerange(2, N + 1))
    for a in range(1, N + 1):
        for p in primes:
            if p * a > N:
                break
            adj[a].append(p * a)
            adj[p * a].append(a)
    return V, adj


def ratio_graph(N, k_max=4):
    """Generalised Hasse: edges (a, b) with b/a in {2, ..., k_max} (any integer)."""
    V = list(range(1, N + 1))
    adj = {v: [] for v in V}
    for a in range(1, N + 1):
        for k in range(2, k_max + 1):
            b = k * a
            if b > N:
                break
            adj[a].append(b)
            adj[b].append(a)
    return V, adj


def winnable_subtraction_set(D, q, V, adj):
    return [v for v in V if is_winnable({u: D[u] - (1 if u == v else 0) for u in V}, q, V, adj)]


def reduced_signature(D, q, V, adj):
    Dr = q_reduce(D, q, V, adj)
    return {
        "D_red_at_q": Dr[q],
        "support": [v for v in V if v != q and Dr[v] > 0],
        "max_chip": max(Dr[v] for v in V if v != q),
        "total_nonq": sum(Dr[v] for v in V if v != q),
        "Dr": dict(Dr),
    }


def run(N, graph_fn, name, n_random, seed):
    print(f"\n=== N={N}, graph={name} ===", flush=True)
    V, adj = graph_fn(N)
    n, m, g = graph_stats(V, adj)
    if not is_connected(V, adj):
        print(f"  DISCONNECTED — skip")
        return None
    q = V[0]
    primes = [v for v in V if isprime(v)]
    deg = len(primes)
    print(f"  |V|={n}, |E|={m}, g={g}, deg(D_P)={deg}, q={q}", flush=True)

    # Prime divisor
    DP = prime_divisor(V)
    t0 = time.time()
    sig_P = reduced_signature(DP, q, V, adj)
    W_P = winnable_subtraction_set(DP, q, V, adj)
    t_P = time.time() - t0
    print(f"  D_P: D'(q)={sig_P['D_red_at_q']}, |supp|={len(sig_P['support'])}, "
          f"max={sig_P['max_chip']}, total_nonq={sig_P['total_nonq']}, |W|={len(W_P)} ({t_P:.1f}s)", flush=True)
    print(f"    support: {sig_P['support'][:30]}", flush=True)

    # Random matched controls
    rng = random.Random(seed)
    rand_data = []
    for i in range(n_random):
        DR = random_matched_divisor(V, deg, rng)
        sig_R = reduced_signature(DR, q, V, adj)
        W_R = winnable_subtraction_set(DR, q, V, adj)
        rand_data.append({
            "D_red_at_q": sig_R["D_red_at_q"],
            "supp_size": len(sig_R["support"]),
            "max_chip": sig_R["max_chip"],
            "W_size": len(W_R),
        })

    # Compute z-scores
    def z(P_val, R_vals):
        m = statistics.mean(R_vals)
        s = statistics.stdev(R_vals) if len(R_vals) > 1 else 0
        return (P_val - m) / s if s > 0 else (0.0 if P_val == m else float('inf')), m, s

    zq, mq, sq = z(sig_P["D_red_at_q"], [d["D_red_at_q"] for d in rand_data])
    zsupp, msupp, ssupp = z(len(sig_P["support"]), [d["supp_size"] for d in rand_data])
    zmax, mmax, smax = z(sig_P["max_chip"], [d["max_chip"] for d in rand_data])
    zW, mW, sW = z(len(W_P), [d["W_size"] for d in rand_data])

    print(f"  Z-scores (D_P vs {n_random} random):")
    print(f"    D'(q):    P={sig_P['D_red_at_q']:5d}  R={mq:6.2f}±{sq:.2f}  z={zq}")
    print(f"    |supp|:   P={len(sig_P['support']):5d}  R={msupp:6.2f}±{ssupp:.2f}  z={zsupp}")
    print(f"    max chip: P={sig_P['max_chip']:5d}  R={mmax:6.2f}±{smax:.2f}  z={zmax}")
    print(f"    |W|:      P={len(W_P):5d}  R={mW:6.2f}±{sW:.2f}  z={zW}")

    return {
        "N": N,
        "graph": name,
        "n": n, "m": m, "g": g,
        "deg": deg,
        "sig_P": sig_P,
        "W_P": W_P,
        "rand_data": rand_data,
        "z_at_q": zq,
        "z_supp": zsupp,
        "z_max": zmax,
        "z_W": zW,
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--N", type=int, nargs="+", default=[32, 64, 128])
    p.add_argument("--n_random", type=int, default=20)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--out", default="hasse_results.json")
    args = p.parse_args()

    out = []
    for N in args.N:
        for name, fn in [("hasse", hasse_graph),
                         ("ratio2", lambda N: ratio_graph(N, 2)),
                         ("ratio3", lambda N: ratio_graph(N, 3))]:
            res = run(N, fn, name, args.n_random, args.seed)
            if res is not None:
                out.append(res)
            with open(args.out, "w") as f:
                json.dump(out, f, indent=2, default=str)


if __name__ == "__main__":
    main()

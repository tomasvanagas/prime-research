"""
D17 — Discrete Morse theory on the divisibility-poset Hasse diagram.

We compute Forman's discrete-Morse critical-cell counts on the 1-skeleton
of the divisibility-poset order complex on [1, N], for N in {64, 256, 1024, 4096}.

Greedy random elementary collapse:
  - At each step, pick a uniformly random degree-1 vertex (a "free face").
  - Pair it with its unique incident edge (collapse).
  - Repeat until no degree-1 vertex remains.
  - Output critical-cell counts m_0 = #(remaining vertices),
    m_1 = #(remaining edges).

Run R seeds per N; report (min, mean, max).
Compare to Erdős–Rényi G(N, |E_div|) baseline (same vertex count and edge count
but uniform-random edges).

Forman's Morse-Euler formula: m_0 - m_1 = chi(graph) = N - |E|.

Cross-domain refs:
  Forman 2002 "User's guide to discrete Morse theory" Sém. Lothar. Combin. 48
  Benedetti-Lutz 2014 Exp. Math. 23, 66 (random discrete Morse)
  Joswig-Pfetsch 2006 Discrete Comput. Geom. (optimal-matching NP-hardness)
"""

import argparse
import json
import random
import time
from math import log
from sympy import primerange, primepi


def build_divisibility_hasse(N):
    """Hasse diagram of (Z_{>=1} ∩ [1, N], |). Edges = covering relations.

    n covers m  iff  m | n and n/m is prime.
    """
    primes = list(primerange(2, N + 1))
    vertices = list(range(1, N + 1))
    edges = []
    for m in range(1, N + 1):
        for p in primes:
            if m * p > N:
                break
            edges.append((m, m * p))
    return vertices, edges


def random_morse_collapse(num_vertices, edges, seed):
    """Greedy random Morse collapse on the 1-skeleton.

    Returns (m_0, m_1, n_steps).
    """
    rng = random.Random(seed)
    # Build adjacency
    adj = [set() for _ in range(num_vertices + 2)]  # 1-indexed
    for (u, v) in edges:
        adj[u].add(v)
        adj[v].add(u)

    n_alive_vertices = num_vertices
    n_alive_edges = len(edges)
    deg1 = [v for v in range(1, num_vertices + 1) if len(adj[v]) == 1]

    while deg1:
        # Pick random
        idx = rng.randint(0, len(deg1) - 1)
        v = deg1[idx]
        deg1[idx] = deg1[-1]
        deg1.pop()
        if len(adj[v]) != 1:
            continue
        u = next(iter(adj[v]))
        adj[v].remove(u)
        adj[u].remove(v)
        n_alive_vertices -= 1
        n_alive_edges -= 1
        if len(adj[u]) == 1:
            deg1.append(u)

    return n_alive_vertices, n_alive_edges


def er_baseline(num_vertices, num_edges, seed):
    """Erdős–Rényi G(n, m) on n=num_vertices vertices with m=num_edges
    random simple edges (no self-loops, no multi-edges)."""
    rng = random.Random(seed)
    edges = set()
    while len(edges) < num_edges:
        u = rng.randint(1, num_vertices)
        v = rng.randint(1, num_vertices)
        if u == v:
            continue
        if u > v:
            u, v = v, u
        edges.add((u, v))
    return list(edges)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N-list", type=int, nargs="+", default=[64, 256, 1024, 4096])
    ap.add_argument("--seeds", type=int, default=20)
    ap.add_argument("--baselines", type=int, default=20)
    ap.add_argument("--out", type=str,
                    default="d17_discrete_morse_divisibility_data.json")
    args = ap.parse_args()

    results = {"N_list": args.N_list, "seeds": args.seeds,
               "baselines": args.baselines, "by_N": {}}

    for N in args.N_list:
        t0 = time.time()
        vertices, edges = build_divisibility_hasse(N)
        n_v = len(vertices)
        n_e = len(edges)
        chi_div = n_v - n_e
        pi_N = int(primepi(N))

        # Run R seeds of greedy random collapse on the divisibility Hasse diagram
        div_runs = []
        for s in range(args.seeds):
            m0, m1 = random_morse_collapse(n_v, edges, seed=s)
            div_runs.append((m0, m1))
            assert m0 - m1 == chi_div, (m0, m1, chi_div)

        # Run R baseline ER(n, m) graphs and apply same collapse
        er_runs = []
        for s in range(args.baselines):
            er_edges = er_baseline(n_v, n_e, seed=10_000 + s)
            m0, m1 = random_morse_collapse(n_v, er_edges, seed=20_000 + s)
            er_runs.append((m0, m1))
            assert m0 - m1 == chi_div, (m0, m1, chi_div, n_v, n_e)

        div_m0 = [r[0] for r in div_runs]
        div_m1 = [r[1] for r in div_runs]
        er_m0 = [r[0] for r in er_runs]
        er_m1 = [r[1] for r in er_runs]

        elapsed = time.time() - t0
        rec = {
            "N": N,
            "n_vertices": n_v,
            "n_edges": n_e,
            "chi": chi_div,
            "pi_N": pi_N,
            "log2_N": log(N, 2),
            "log_N": log(N),
            "div_m0_min": min(div_m0),
            "div_m0_mean": sum(div_m0) / len(div_m0),
            "div_m0_max": max(div_m0),
            "div_m1_min": min(div_m1),
            "div_m1_mean": sum(div_m1) / len(div_m1),
            "div_m1_max": max(div_m1),
            "er_m0_min": min(er_m0),
            "er_m0_mean": sum(er_m0) / len(er_m0),
            "er_m0_max": max(er_m0),
            "er_m1_min": min(er_m1),
            "er_m1_mean": sum(er_m1) / len(er_m1),
            "er_m1_max": max(er_m1),
            "elapsed_s": round(elapsed, 2),
        }
        results["by_N"][N] = rec
        print(f"N={N}: |V|={n_v} |E|={n_e} pi(N)={pi_N} chi={chi_div}")
        print(f"  div  m0 min/mean/max = {rec['div_m0_min']}/{rec['div_m0_mean']:.1f}/{rec['div_m0_max']}")
        print(f"  div  m1 min/mean/max = {rec['div_m1_min']}/{rec['div_m1_mean']:.1f}/{rec['div_m1_max']}")
        print(f"  ER   m0 min/mean/max = {rec['er_m0_min']}/{rec['er_m0_mean']:.1f}/{rec['er_m0_max']}")
        print(f"  ER   m1 min/mean/max = {rec['er_m1_min']}/{rec['er_m1_mean']:.1f}/{rec['er_m1_max']}")
        print(f"  elapsed: {elapsed:.2f}s")

    with open(args.out, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nWrote {args.out}")


if __name__ == "__main__":
    main()

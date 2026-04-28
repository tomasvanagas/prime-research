"""
D45 — Baker-Norine graph Riemann-Roch / chip-firing rank of the prime divisor.

Computes the divisor rank r_G(D_P^N) of the prime divisor on:
  (A) divisibility graph Γ_N on [1, N] (1 included to ensure connectivity)
  (B) coprimality graph Γ'_N on [2, N]

Compares against matched-degree random divisors. Tests for Brill-Noether
specialness: is r(D) > max{-1, deg(D) - g}?

Cross-domain:
- Baker-Norine 2007 Adv. Math. 215, 766 = arXiv:math/0608360
- Dhar 1990 Phys. Rev. Lett. 64, 1613 (chip-firing burning)
"""

import sys
import json
import time
import argparse
import random
from math import gcd
from collections import deque
from sympy import isprime, primerange


# -----------------------------------------------------------------------------
# Graph construction
# -----------------------------------------------------------------------------

def divisibility_graph(N):
    """V = [1, N], edges between (a, b) with a | b, a != b. Connected because
    1 divides everything."""
    V = list(range(1, N + 1))
    adj = {v: [] for v in V}
    for a in range(1, N + 1):
        for b in range(2 * a, N + 1, a):
            adj[a].append(b)
            adj[b].append(a)
    return V, adj


def coprimality_graph(N):
    """V = [2, N], edges between (a, b) with gcd(a, b) = 1."""
    V = list(range(2, N + 1))
    adj = {v: [] for v in V}
    for i, a in enumerate(V):
        for b in V[i + 1:]:
            if gcd(a, b) == 1:
                adj[a].append(b)
                adj[b].append(a)
    return V, adj


def graph_stats(V, adj):
    n = len(V)
    m = sum(len(adj[v]) for v in V) // 2
    g = m - n + 1
    return n, m, g


def is_connected(V, adj):
    if not V:
        return True
    seen = {V[0]}
    queue = deque([V[0]])
    while queue:
        u = queue.popleft()
        for w in adj[u]:
            if w not in seen:
                seen.add(w)
                queue.append(w)
    return len(seen) == len(V)


# -----------------------------------------------------------------------------
# Chip firing primitives
# -----------------------------------------------------------------------------

def fire_vertex(D, v, adj):
    """Fire vertex v: D[v] -= deg(v), D[u] += 1 for each neighbor u."""
    D[v] -= len(adj[v])
    for u in adj[v]:
        D[u] += 1


def fire_set(D, S, adj):
    """Fire all vertices in S simultaneously."""
    Sset = set(S)
    for v in S:
        D[v] -= len(adj[v])
        for u in adj[v]:
            D[u] += 1


def dhars_burn(D, q, V, adj):
    """Run Dhar's burning algorithm starting at q. Return None if D is
    superstable on V \\ {q}, else return the unburnt set (which is a
    fireable subset)."""
    burnt = {q}
    while True:
        added = False
        for v in V:
            if v in burnt:
                continue
            burnt_count = 0
            for u in adj[v]:
                if u in burnt:
                    burnt_count += 1
            if burnt_count > D[v]:
                burnt.add(v)
                added = True
        if not added:
            break
    if len(burnt) == len(V):
        return None
    return [v for v in V if v not in burnt]


def q_reduce(D_in, q, V, adj, max_iter=200000):
    """Compute the q-reduced form of D.

    Phase 1: chip-lending. Repeatedly: pick vertex v != q with D[v] < 0,
    BFS shortest path from q to v in the graph, fire all vertices on the
    path except v. This propagates one chip from q toward v.

    Phase 2: superstabilize. Run Dhar's burning; while violators exist,
    fire them.
    """
    D = dict(D_in)
    iters = 0

    # Phase 1: chip lending
    while iters < max_iter:
        bad = [v for v in V if v != q and D[v] < 0]
        if not bad:
            break
        # Pick the most-deficit vertex
        v = min(bad, key=lambda w: D[w])
        # BFS from q to v
        parent = {q: None}
        queue = deque([q])
        while queue:
            u = queue.popleft()
            if u == v:
                break
            for w in adj[u]:
                if w not in parent:
                    parent[w] = u
                    queue.append(w)
        # Reconstruct path
        path = []
        u = v
        while u is not None:
            path.append(u)
            u = parent[u]
        path.reverse()
        # Fire each vertex on the path except v (the target)
        for u in path[:-1]:
            fire_vertex(D, u, adj)
        iters += 1

    if iters >= max_iter:
        raise RuntimeError(f"Phase 1 did not converge after {max_iter} iters")

    # Phase 2: superstabilize via Dhar's burn
    while iters < max_iter:
        violators = dhars_burn(D, q, V, adj)
        if violators is None:
            break
        fire_set(D, violators, adj)
        iters += 1

    if iters >= max_iter:
        raise RuntimeError(f"Phase 2 did not converge after {max_iter} iters")

    return D


# -----------------------------------------------------------------------------
# Rank computation
# -----------------------------------------------------------------------------

def is_winnable(D, q, V, adj):
    """Test if D is linearly equivalent to an effective divisor.
    Equivalent: q-reduced form has D'(q) >= 0."""
    D_red = q_reduce(D, q, V, adj)
    return D_red[q] >= 0


def rank(D, q, V, adj, max_k=3, time_limit=120.0):
    """Compute Baker-Norine rank up to max_k. Returns (rank, info_dict).

    Algorithm: r >= k iff for every effective E of degree k, |D - E| != ∅.
    Test by enumerating multisets of size k.
    """
    from itertools import combinations_with_replacement

    info = {}
    # k=0: D winnable?
    if not is_winnable(D, q, V, adj):
        return -1, {}
    info[0] = "winnable"

    for k in range(1, max_k + 1):
        t0 = time.time()
        all_pass = True
        first_fail = None
        n_tested = 0
        for E_support in combinations_with_replacement(V, k):
            if time.time() - t0 > time_limit:
                info[k] = f"TIMEOUT after {n_tested} multisets"
                return k - 1, info  # conservative: return last confirmed rank
            D_minus_E = dict(D)
            for v in E_support:
                D_minus_E[v] -= 1
            if not is_winnable(D_minus_E, q, V, adj):
                all_pass = False
                first_fail = E_support
                break
            n_tested += 1
        if all_pass:
            info[k] = f"r >= {k} (tested {n_tested} multisets)"
        else:
            info[k] = f"FAILED at E={first_fail} (tested {n_tested} multisets)"
            return k - 1, info
    return max_k, info


# -----------------------------------------------------------------------------
# Divisor builders
# -----------------------------------------------------------------------------

def prime_divisor(V):
    """D_P(v) = 1 if v prime, 0 otherwise."""
    return {v: (1 if isprime(v) else 0) for v in V}


def random_matched_divisor(V, deg, rng, exclude_leaves=False, adj=None):
    """Place `deg` chips uniformly at random on |V| vertices (no two on
    same vertex, like the prime divisor). If exclude_leaves, skip vertices
    with degree 1."""
    candidates = [v for v in V if not (exclude_leaves and adj is not None and len(adj[v]) <= 1)]
    if deg > len(candidates):
        raise ValueError(f"deg {deg} > |candidates| {len(candidates)}")
    chosen = rng.sample(candidates, deg)
    return {v: (1 if v in set(chosen) else 0) for v in V}


# -----------------------------------------------------------------------------
# Main experiment driver
# -----------------------------------------------------------------------------

def run_experiment(N, graph_type, n_random, max_k, q_choice, time_limit, seed):
    """Run for one (N, graph_type) cell."""
    print(f"\n=== N={N}, graph={graph_type}, q={q_choice}, max_k={max_k} ===", flush=True)
    if graph_type == "div":
        V, adj = divisibility_graph(N)
    elif graph_type == "coprime":
        V, adj = coprimality_graph(N)
    else:
        raise ValueError(graph_type)
    n, m, g = graph_stats(V, adj)
    connected = is_connected(V, adj)
    print(f"  |V|={n}, |E|={m}, genus g={g}, connected={connected}", flush=True)

    if not connected:
        print(f"  SKIP: graph disconnected", flush=True)
        return {"N": N, "graph_type": graph_type, "skipped": "disconnected"}

    # Pick q
    if q_choice == "low":
        q = V[0]
    elif q_choice == "max_deg":
        q = max(V, key=lambda v: len(adj[v]))
    else:
        q = q_choice  # explicit

    primes = [v for v in V if isprime(v)]
    deg_DP = len(primes)
    bn_bound = max(-1, deg_DP - g)
    print(f"  deg(D_P)={deg_DP}, BN expected r = max(-1, deg-g) = {bn_bound}", flush=True)

    # PRIMES rank
    D_P = prime_divisor(V)
    t0 = time.time()
    r_P, info_P = rank(D_P, q, V, adj, max_k=max_k, time_limit=time_limit)
    t_P = time.time() - t0
    print(f"  r(D_P) = {r_P} ({t_P:.1f}s), info={info_P}", flush=True)

    # Random controls
    rng = random.Random(seed)
    r_randoms = []
    for i in range(n_random):
        D_R = random_matched_divisor(V, deg_DP, rng)
        t0 = time.time()
        r_R, info_R = rank(D_R, q, V, adj, max_k=max_k, time_limit=time_limit)
        t_R = time.time() - t0
        r_randoms.append(r_R)
        print(f"    rand {i+1}/{n_random}: r={r_R} ({t_R:.1f}s)", flush=True)

    # z-score
    if len(r_randoms) > 1:
        import statistics
        mean_R = statistics.mean(r_randoms)
        sd_R = statistics.stdev(r_randoms)
        z = (r_P - mean_R) / sd_R if sd_R > 0 else float('inf') if r_P != mean_R else 0.0
    else:
        mean_R = r_randoms[0] if r_randoms else None
        sd_R = 0
        z = None

    return {
        "N": N,
        "graph_type": graph_type,
        "q": q,
        "n_vertices": n,
        "n_edges": m,
        "genus": g,
        "deg_DP": deg_DP,
        "BN_bound": bn_bound,
        "r_DP": r_P,
        "info_DP": info_P,
        "time_DP": t_P,
        "n_random": len(r_randoms),
        "r_random": r_randoms,
        "mean_random": mean_R,
        "sd_random": sd_R,
        "z_score": z,
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--N", type=int, nargs="+", default=[16, 32, 64, 128])
    p.add_argument("--graph", choices=["div", "coprime", "both"], default="both")
    p.add_argument("--n_random", type=int, default=20)
    p.add_argument("--max_k", type=int, default=2)
    p.add_argument("--q", default="low")
    p.add_argument("--time_limit", type=float, default=60.0)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--out", default="results.json")
    args = p.parse_args()

    results = []
    graphs = ["div", "coprime"] if args.graph == "both" else [args.graph]
    for N in args.N:
        for gt in graphs:
            res = run_experiment(N, gt, args.n_random, args.max_k, args.q,
                                 args.time_limit, args.seed)
            results.append(res)
            with open(args.out, "w") as f:
                json.dump(results, f, indent=2, default=str)

    print("\n=== SUMMARY ===")
    for r in results:
        if "skipped" in r:
            print(f"  N={r['N']:3d} {r['graph_type']:>7s}: SKIPPED ({r['skipped']})")
            continue
        print(f"  N={r['N']:3d} {r['graph_type']:>7s}: r(D_P)={r['r_DP']:3d}  "
              f"deg={r['deg_DP']:3d}  g={r['genus']:5d}  BN={r['BN_bound']:3d}  "
              f"rand_mean={r['mean_random']:.2f} sd={r['sd_random']:.2f} z={r['z_score']}")


if __name__ == "__main__":
    main()

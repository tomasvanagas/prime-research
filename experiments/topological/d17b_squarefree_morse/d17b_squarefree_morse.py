"""
C11 / D17.b — Discrete Morse on the SQUAREFREE-only divisibility Hasse diagram.

Vertices V_sqfree(N) = {n in [1, N] : mu(n) != 0}.
Covering edges: (m, mp) where m, mp both squarefree and p prime
                (equivalently m | mp, mp/m = p, p does not divide m).

This is exactly the Boolean lattice on primes <= N truncated to subsets S
with prod(S) <= N. The full Boolean lattice's *order complex* is shellable
with one critical cell, but the 1-skeleton (Hasse diagram) is the n-cube
graph — which has minimum degree n and is not graph-collapsible.

Run greedy random elementary Morse collapse (Forman 2002; Benedetti-Lutz
2014), report m_0, m_1, compare to ER(|V|, |E|) baseline. Also decompose
collapse cascade into waves to look for a structural identity analogous to
the D17 closed form `collapses(N) = (π(N) - π(N/2)) + Π_pow(N) + 1`.

Pre-stated falsifiers (see d17b_squarefree_morse_results.md):

  F1 (A-grade): m_0(H_N^sqfree) = O(polylog N).
  F2 (B-grade ER match): |m_0(div) - m_0(ER)| / |V_sqfree(N)| <= 0.02
                           uniformly in N.
  F3 (B-grade closed form): m_0 admits an exact arithmetic
                            decomposition into terms involving π and ω.
  F4 (Morse-rigidity): random Morse output deterministic across seeds.
"""

import argparse
import json
import random
import time
from collections import Counter
from math import log
from sympy import primerange, primepi, mobius, factorint


def squarefree_in_range(N):
    """Return list of squarefree integers in [1, N]."""
    is_sqfree = [True] * (N + 1)
    is_sqfree[0] = False
    p = 2
    while p * p <= N:
        sq = p * p
        for k in range(sq, N + 1, sq):
            is_sqfree[k] = False
        p += 1
    return [n for n in range(1, N + 1) if is_sqfree[n]]


def build_squarefree_hasse(N):
    """Hasse diagram of the squarefree divisibility poset on [1, N].

    Vertices: squarefree integers in [1, N].
    Edges: (m, mp) for m squarefree, p prime, p does not divide m, mp <= N.
    """
    sqfree = squarefree_in_range(N)
    primes = list(primerange(2, N + 1))
    sqfree_set = set(sqfree)

    edges = []
    for m in sqfree:
        for p in primes:
            mp = m * p
            if mp > N:
                break
            # squarefree condition: p does not divide m
            if m % p == 0:
                continue
            edges.append((m, mp))
    return sqfree, edges


def omega(n):
    """Number of distinct prime factors of n."""
    if n == 1:
        return 0
    return len(factorint(n))


def random_morse_collapse(vertices, edges, seed):
    """Greedy random elementary Morse collapse on the 1-skeleton.

    `vertices` is the *list* of vertex labels (not necessarily {1..n}).
    Returns (m_0, m_1) = (#remaining vertices, #remaining edges).
    """
    rng = random.Random(seed)
    label_to_idx = {v: i for i, v in enumerate(vertices)}
    n = len(vertices)
    adj = [set() for _ in range(n)]
    for (u, v) in edges:
        iu, iv = label_to_idx[u], label_to_idx[v]
        adj[iu].add(iv)
        adj[iv].add(iu)

    n_alive_vertices = n
    n_alive_edges = len(edges)
    deg1 = [i for i in range(n) if len(adj[i]) == 1]

    while deg1:
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
    """ER G(n, m) on n=num_vertices vertices with m=num_edges
    distinct random simple edges (no self-loops, no multi-edges).
    Returns edges as pairs of integer indices in [0, num_vertices)."""
    rng = random.Random(seed)
    edges = set()
    while len(edges) < num_edges:
        u = rng.randint(0, num_vertices - 1)
        v = rng.randint(0, num_vertices - 1)
        if u == v:
            continue
        if u > v:
            u, v = v, u
        edges.add((u, v))
    return list(edges)


def er_collapse(num_vertices, edges, seed):
    """Run greedy collapse on an ER-shape graph with vertices labeled 0..n-1."""
    rng = random.Random(seed)
    adj = [set() for _ in range(num_vertices)]
    for (u, v) in edges:
        adj[u].add(v)
        adj[v].add(u)

    n_alive_vertices = num_vertices
    n_alive_edges = len(edges)
    deg1 = [v for v in range(num_vertices) if len(adj[v]) == 1]

    while deg1:
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


def cascade_breakdown(N):
    """Run a single deterministic-order collapse and report wave structure.

    Returns counts of vertices peeled at each "wave". Wave 0 = initial
    leaves; wave k+1 = vertices that became leaves only after the wave-k
    peeling cascade.
    """
    vertices, edges = build_squarefree_hasse(N)
    label_to_idx = {v: i for i, v in enumerate(vertices)}
    n = len(vertices)
    adj = [set() for _ in range(n)]
    for (u, v) in edges:
        iu, iv = label_to_idx[u], label_to_idx[v]
        adj[iu].add(iv)
        adj[iv].add(iu)

    deg1 = [i for i in range(n) if len(adj[i]) == 1]
    initial_leaves = [vertices[i] for i in deg1]

    wave_count = []
    waves_label_record = {}  # vertex_label -> wave number
    current_wave_set = set(deg1)
    wave_id = 0
    while deg1:
        # Process one full wave at a time.
        next_wave = []
        wave_size = 0
        # We need to peel everyone currently in deg1 (plus everyone that
        # becomes degree-1 *as a consequence* of this wave's peeling — those
        # belong to wave_id + 1).
        # To do this cleanly: snapshot deg1, then peel them in a fixed order.
        snapshot = list(set(deg1))
        deg1 = []
        for v in snapshot:
            if len(adj[v]) != 1:
                continue
            u = next(iter(adj[v]))
            adj[v].remove(u)
            adj[u].remove(v)
            waves_label_record[vertices[v]] = wave_id
            wave_size += 1
            if len(adj[u]) == 1:
                deg1.append(u)
        wave_count.append(wave_size)
        wave_id += 1
        current_wave_set = set(deg1)

    # Remaining vertices and edges: critical
    n_alive = sum(1 for i in range(n) if len(adj[i]) > 0 or i not in waves_label_record and True)
    # Better: m_0 = vertices not collapsed = those without an entry in waves_label_record? No,
    # vertices with adj == 0 may still be present (isolated). m_0 = vertices not paired in collapse.
    # Each collapse pairs (vertex, edge). So m_0 = n - sum(wave_count).
    m_0 = n - sum(wave_count)
    m_1 = len(edges) - sum(wave_count)
    return {
        "N": N,
        "n_vertices": n,
        "n_edges": len(edges),
        "m_0": m_0,
        "m_1": m_1,
        "wave_sizes": wave_count,
        "n_initial_leaves": len(initial_leaves),
        "n_waves": len(wave_count),
    }


def omega_distribution(N):
    """For squarefree n in [1, N], count by ω(n)."""
    sqfree = squarefree_in_range(N)
    counts = Counter()
    for n in sqfree:
        counts[omega(n)] += 1
    return dict(counts)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N-list", type=int, nargs="+",
                    default=[64, 128, 256, 512, 1024, 2048, 4096])
    ap.add_argument("--seeds", type=int, default=20)
    ap.add_argument("--baselines", type=int, default=20)
    ap.add_argument("--determinism-N", type=int, nargs="+",
                    default=[64, 256, 1024])
    ap.add_argument("--determinism-seeds", type=int, default=200)
    ap.add_argument("--out", type=str,
                    default="d17b_squarefree_morse_data.json")
    args = ap.parse_args()

    results = {
        "N_list": args.N_list,
        "seeds": args.seeds,
        "baselines": args.baselines,
        "by_N": {},
        "determinism": {},
        "cascade": {},
    }

    for N in args.N_list:
        t0 = time.time()
        vertices, edges = build_squarefree_hasse(N)
        n_v = len(vertices)
        n_e = len(edges)
        chi = n_v - n_e
        pi_N = int(primepi(N))
        pi_Nhalf = int(primepi(N // 2))
        omega_dist = omega_distribution(N)

        # divisibility-Hasse runs
        div_runs = []
        for s in range(args.seeds):
            m0, m1 = random_morse_collapse(vertices, edges, seed=s)
            div_runs.append((m0, m1))
            assert m0 - m1 == chi, (m0, m1, chi)

        # ER baseline runs
        er_runs = []
        for s in range(args.baselines):
            er_e = er_baseline(n_v, n_e, seed=10_000 + s)
            m0, m1 = er_collapse(n_v, er_e, seed=20_000 + s)
            er_runs.append((m0, m1))
            assert m0 - m1 == chi, (m0, m1, chi, n_v, n_e)

        div_m0 = [r[0] for r in div_runs]
        div_m1 = [r[1] for r in div_runs]
        er_m0 = [r[0] for r in er_runs]
        er_m1 = [r[1] for r in er_runs]

        elapsed = time.time() - t0
        rec = {
            "N": N,
            "n_vertices_sqfree": n_v,
            "n_edges": n_e,
            "chi": chi,
            "pi_N": pi_N,
            "pi_N_half": pi_Nhalf,
            "primes_in_(N/2,N]": pi_N - pi_Nhalf,
            "omega_distribution": omega_dist,
            "log2_N": log(N, 2),
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
        print(f"N={N}: |V_sqf|={n_v} |E|={n_e} chi={chi} pi(N)={pi_N} "
              f"pi(N)-pi(N/2)={pi_N - pi_Nhalf}")
        print(f"  div  m0 min/mean/max = "
              f"{rec['div_m0_min']}/{rec['div_m0_mean']:.1f}/{rec['div_m0_max']}; "
              f"m0/|V| = {rec['div_m0_mean']/n_v:.4f}")
        print(f"  ER   m0 min/mean/max = "
              f"{rec['er_m0_min']}/{rec['er_m0_mean']:.1f}/{rec['er_m0_max']}; "
              f"m0/|V| = {rec['er_m0_mean']/n_v:.4f}")
        print(f"  Δ(div - ER mean) = {rec['div_m0_mean'] - rec['er_m0_mean']:+.2f}; "
              f"elapsed: {elapsed:.2f}s")

    # determinism test (different seeds, same N)
    print("\n=== Determinism test ===")
    for N in args.determinism_N:
        vertices, edges = build_squarefree_hasse(N)
        outputs = Counter()
        for s in range(args.determinism_seeds):
            m0, m1 = random_morse_collapse(vertices, edges, seed=s * 991 + 17)
            outputs[(m0, m1)] += 1
        results["determinism"][N] = {str(k): v for k, v in outputs.items()}
        print(f"N={N}: distinct (m0, m1) outputs (out of {args.determinism_seeds}): "
              f"{len(outputs)}; top: {outputs.most_common(3)}")

    # cascade breakdown for closed-form hunting
    print("\n=== Cascade wave breakdown ===")
    for N in args.N_list:
        rec = cascade_breakdown(N)
        results["cascade"][N] = rec
        print(f"N={N}: m_0={rec['m_0']} m_1={rec['m_1']} "
              f"n_waves={rec['n_waves']} initial_leaves={rec['n_initial_leaves']} "
              f"wave_sizes={rec['wave_sizes'][:8]}...")

    with open(args.out, "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nWrote {args.out}")


if __name__ == "__main__":
    main()

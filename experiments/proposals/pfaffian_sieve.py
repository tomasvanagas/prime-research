#!/usr/bin/env python3
"""
PROPOSAL: Planar Graph Matching / Pfaffian Approach to Prime Counting

Key idea: The inclusion-exclusion sieve resembles a permanent/partition-function.
For planar graphs, the FKT algorithm computes such quantities in polynomial time.
Question: does the sieve computation have bounded treewidth?

Two encodings tested:
  1. Sieve bipartite graph: rows=integers 1..x, cols=primes<=sqrt(x),
     edge (i,p) if p | i.  Check planarity and treewidth.
  2. Lucy DP graph: vertices = distinct floor(x/k) values,
     edges from the DP recurrence. Check treewidth/pathwidth.

If treewidth is bounded (or grows slowly), FPT algorithms could help.
If treewidth grows as Omega(sqrt(x)/log(x)), the approach is closed.

Author: Claude (Session 26)
Date: 2026-04-05
"""

import time
import math
import sys

try:
    import networkx as nx
    from networkx.algorithms.approximation import treewidth_min_degree
except ImportError:
    print("ERROR: networkx required. Install with: pip install networkx")
    sys.exit(1)

import sympy


# ---------------------------------------------------------------------------
# Utility: list of primes up to n
# ---------------------------------------------------------------------------
def primes_up_to(n):
    """Simple sieve of Eratosthenes."""
    if n < 2:
        return []
    sieve = [True] * (n + 1)
    sieve[0] = sieve[1] = False
    for i in range(2, int(n**0.5) + 1):
        if sieve[i]:
            for j in range(i * i, n + 1, i):
                sieve[j] = False
    return [i for i in range(2, n + 1) if sieve[i]]


# ===========================================================================
# EXPERIMENT 1: Sieve Bipartite Graph
# ===========================================================================
def build_sieve_bipartite_graph(x):
    """
    Build bipartite graph G = (Integers, Primes, E)
    where integers are 1..x, primes are primes <= sqrt(x),
    and edge (i, p) exists iff p divides i.

    Returns the graph and basic stats.
    """
    sqrt_x = int(math.isqrt(x))
    ps = primes_up_to(sqrt_x)

    G = nx.Graph()

    # Add integer nodes (prefix 'i') and prime nodes (prefix 'p')
    for i in range(1, x + 1):
        G.add_node(f"i{i}", bipartite=0)
    for p in ps:
        G.add_node(f"p{p}", bipartite=1)

    edge_count = 0
    for p in ps:
        for mult in range(p, x + 1, p):
            G.add_edge(f"i{mult}", f"p{p}")
            edge_count += 1

    return G, len(ps), edge_count


def test_sieve_graph_planarity(x_values):
    """Test planarity and treewidth of the sieve bipartite graph."""
    print("=" * 70)
    print("EXPERIMENT 1: Sieve Bipartite Graph")
    print("=" * 70)
    print(f"{'x':>8} | {'#primes':>8} | {'#nodes':>8} | {'#edges':>8} | "
          f"{'planar':>8} | {'tw_approx':>10} | {'time(s)':>8}")
    print("-" * 70)

    results = []
    for x in x_values:
        t0 = time.time()
        G, n_primes, n_edges = build_sieve_bipartite_graph(x)
        n_nodes = G.number_of_nodes()

        # Check planarity
        is_planar = nx.check_planarity(G)[0]

        # Approximate treewidth (min-degree heuristic)
        tw, _ = treewidth_min_degree(G)

        elapsed = time.time() - t0

        print(f"{x:>8} | {n_primes:>8} | {n_nodes:>8} | {n_edges:>8} | "
              f"{str(is_planar):>8} | {tw:>10} | {elapsed:>8.3f}")

        results.append({
            'x': x, 'n_primes': n_primes, 'n_nodes': n_nodes,
            'n_edges': n_edges, 'planar': is_planar, 'treewidth': tw
        })

    return results


# ===========================================================================
# EXPERIMENT 2: Lucy DP Graph
# ===========================================================================
def lucy_dp_values(x):
    """
    Compute the set of distinct floor(x/k) values (the 'Lucy DP' set).
    These are the subproblem states in the Meissel-Lehmer method.
    Returns sorted list of distinct values.
    """
    sqrt_x = int(math.isqrt(x))
    vals = set()
    for k in range(1, sqrt_x + 1):
        vals.add(x // k)
        vals.add(k)
    return sorted(vals)


def build_lucy_dp_graph(x):
    """
    Build a graph where:
    - Vertices = distinct floor(x/k) values
    - Edges connect v to v' if the Lucy DP recurrence links them:
      S(v, p) depends on S(v, p-1) and S(v // p, p-1)

    For each value v and each prime p <= sqrt(v):
        edge(v, floor(v/p))
    """
    vals = lucy_dp_values(x)
    val_set = set(vals)
    sqrt_x = int(math.isqrt(x))
    ps = primes_up_to(sqrt_x)

    G = nx.Graph()
    for v in vals:
        G.add_node(v)

    edge_count = 0
    for v in vals:
        for p in ps:
            if p * p > v:
                break
            target = v // p
            if target in val_set:
                if not G.has_edge(v, target):
                    G.add_edge(v, target)
                    edge_count += 1

    return G, len(vals), edge_count


def test_lucy_dp_graph(x_values):
    """Test treewidth of the Lucy DP graph."""
    print()
    print("=" * 70)
    print("EXPERIMENT 2: Lucy DP Recurrence Graph")
    print("=" * 70)
    print(f"{'x':>8} | {'#values':>8} | {'#edges':>8} | "
          f"{'planar':>8} | {'tw_approx':>10} | {'time(s)':>8}")
    print("-" * 70)

    results = []
    for x in x_values:
        t0 = time.time()
        G, n_vals, n_edges = build_lucy_dp_graph(x)

        is_planar = nx.check_planarity(G)[0]
        tw, _ = treewidth_min_degree(G)

        elapsed = time.time() - t0

        print(f"{x:>8} | {n_vals:>8} | {n_edges:>8} | "
              f"{str(is_planar):>8} | {tw:>10} | {elapsed:>8.3f}")

        results.append({
            'x': x, 'n_vals': n_vals, 'n_edges': n_edges,
            'planar': is_planar, 'treewidth': tw
        })

    return results


# ===========================================================================
# EXPERIMENT 3: Treewidth Growth Rate Analysis
# ===========================================================================
def analyze_growth(results, label, x_key='x', tw_key='treewidth'):
    """Fit treewidth growth to power law: tw ~ x^alpha."""
    print(f"\n--- {label}: Treewidth Growth Analysis ---")

    if len(results) < 3:
        print("Not enough data points for fitting.")
        return None

    # Filter out zero treewidths
    filtered = [(r[x_key], r[tw_key]) for r in results if r[tw_key] > 0]
    if len(filtered) < 2:
        print("Not enough nonzero treewidth values.")
        return None

    xs = [math.log(f[0]) for f in filtered]
    tws = [math.log(f[1]) for f in filtered]

    # Least squares fit: log(tw) = alpha * log(x) + beta
    n = len(xs)
    sum_x = sum(xs)
    sum_y = sum(tws)
    sum_xy = sum(a * b for a, b in zip(xs, tws))
    sum_x2 = sum(a * a for a in xs)

    denom = n * sum_x2 - sum_x * sum_x
    if abs(denom) < 1e-12:
        print("Degenerate fit.")
        return None

    alpha = (n * sum_xy - sum_x * sum_y) / denom
    beta = (sum_y - alpha * sum_x) / n
    C = math.exp(beta)

    print(f"  Fit: treewidth ~ {C:.4f} * x^{alpha:.4f}")
    print(f"  Exponent alpha = {alpha:.4f}")

    # Compare with theoretical expectations
    if alpha > 0.4:
        print(f"  --> Treewidth grows as ~x^{alpha:.2f}: UNBOUNDED, approach fails.")
    elif alpha > 0.1:
        print(f"  --> Treewidth grows moderately: probably unbounded for large x.")
    else:
        print(f"  --> Treewidth grows slowly: WORTH INVESTIGATING further!")

    return alpha


# ===========================================================================
# EXPERIMENT 4: Reduced sieve graph (only composite structure)
# ===========================================================================
def build_prime_coprimality_graph(x):
    """
    Build a graph where vertices = primes up to sqrt(x),
    edge (p, q) if there exists an integer <= x divisible by both p and q.
    (This always exists if p*q <= x, so the graph is nearly complete
    for small primes, but the interesting structure is for larger ones.)

    More useful: the 'sieve dependency graph' where vertices = primes,
    edge (p,q) weighted by |{n <= x : p|n and q|n}| = floor(x/(p*q)).
    """
    sqrt_x = int(math.isqrt(x))
    ps = primes_up_to(sqrt_x)

    G = nx.Graph()
    for p in ps:
        G.add_node(p)

    for i, p in enumerate(ps):
        for j in range(i + 1, len(ps)):
            q = ps[j]
            if p * q > x:
                break
            # Number of integers <= x divisible by both p and q
            count = x // (p * q)
            if count > 0:
                G.add_edge(p, q, weight=count)

    return G


def test_prime_dependency_graph(x_values):
    """Test treewidth of the prime-prime dependency graph."""
    print()
    print("=" * 70)
    print("EXPERIMENT 3: Prime Dependency Graph (primes as vertices)")
    print("=" * 70)
    print(f"{'x':>8} | {'#primes':>8} | {'#edges':>8} | "
          f"{'planar':>8} | {'tw_approx':>10} | {'density':>8} | {'time(s)':>8}")
    print("-" * 70)

    results = []
    for x in x_values:
        t0 = time.time()
        G = build_prime_coprimality_graph(x)
        n_primes = G.number_of_nodes()
        n_edges = G.number_of_edges()
        max_edges = n_primes * (n_primes - 1) // 2 if n_primes > 1 else 1
        density = n_edges / max_edges if max_edges > 0 else 0

        is_planar = nx.check_planarity(G)[0]
        tw = 0
        if n_primes > 1:
            tw, _ = treewidth_min_degree(G)

        elapsed = time.time() - t0

        print(f"{x:>8} | {n_primes:>8} | {n_edges:>8} | "
              f"{str(is_planar):>8} | {tw:>10} | {density:>8.4f} | {elapsed:>8.3f}")

        results.append({
            'x': x, 'n_primes': n_primes, 'n_edges': n_edges,
            'planar': is_planar, 'treewidth': tw, 'density': density
        })

    return results


# ===========================================================================
# MAIN
# ===========================================================================
def main():
    print("Pfaffian / Planar Graph Approach to Prime Counting")
    print("=" * 70)
    print()

    # Range of x values to test
    # Sieve bipartite graph: small x only (graphs get huge)
    sieve_x_vals = [20, 30, 50, 80, 100, 150, 200, 300, 500]

    # Lucy DP graph: can go larger since there are only O(sqrt(x)) nodes
    lucy_x_vals = [100, 500, 1000, 5000, 10000, 50000, 100000, 500000, 1000000]

    # Prime dependency graph
    dep_x_vals = [100, 500, 1000, 5000, 10000, 50000, 100000, 500000, 1000000]

    # --- Experiment 1: Sieve bipartite graph ---
    sieve_results = test_sieve_graph_planarity(sieve_x_vals)

    # --- Experiment 2: Lucy DP graph ---
    lucy_results = test_lucy_dp_graph(lucy_x_vals)

    # --- Experiment 3: Prime dependency graph ---
    dep_results = test_prime_dependency_graph(dep_x_vals)

    # --- Growth analysis ---
    print()
    print("=" * 70)
    print("TREEWIDTH GROWTH ANALYSIS")
    print("=" * 70)

    alpha_sieve = analyze_growth(sieve_results, "Sieve Bipartite Graph")
    alpha_lucy = analyze_growth(lucy_results, "Lucy DP Graph")
    alpha_dep = analyze_growth(dep_results, "Prime Dependency Graph")

    # --- Verdict ---
    print()
    print("=" * 70)
    print("VERDICT")
    print("=" * 70)

    all_planar_sieve = all(r['planar'] for r in sieve_results)
    all_planar_lucy = all(r['planar'] for r in lucy_results)

    print(f"  Sieve graph planar for all tested x? {all_planar_sieve}")
    print(f"  Lucy DP graph planar for all tested x? {all_planar_lucy}")

    if alpha_sieve is not None:
        print(f"  Sieve treewidth exponent: {alpha_sieve:.4f}")
    if alpha_lucy is not None:
        print(f"  Lucy DP treewidth exponent: {alpha_lucy:.4f}")
    if alpha_dep is not None:
        print(f"  Prime dep. treewidth exponent: {alpha_dep:.4f}")

    # Check if any encoding gives bounded/slowly-growing treewidth
    bounded = False
    if alpha_lucy is not None and alpha_lucy < 0.1:
        bounded = True
        print("\n  ** Lucy DP graph has slowly growing treewidth -- INVESTIGATE **")
    if alpha_dep is not None and alpha_dep < 0.1:
        bounded = True
        print("\n  ** Prime dependency graph has slowly growing treewidth -- INVESTIGATE **")

    if not bounded:
        print("\n  All graph encodings show unbounded treewidth growth.")
        print("  The sieve computation is NOT close to planar.")
        print("  Pfaffian/FKT approach: CLOSED.")
        print("  --> The inclusion-exclusion structure has high treewidth,")
        print("      consistent with the information-theoretic barrier.")
    else:
        print("\n  Some encodings show promise -- further investigation needed.")

    print()
    print("Done.")


if __name__ == "__main__":
    main()

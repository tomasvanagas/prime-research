"""Follow-up analysis: verify deterministic-output for divisibility-Hasse greedy
collapse + characterise the structural identity collapses(N) ≈ π(N) − π(N/2)."""

import random
from sympy import primerange, primepi
from collections import Counter

from d17_discrete_morse_divisibility import (
    build_divisibility_hasse, random_morse_collapse
)


def test_determinism(N, n_seeds=200):
    vertices, edges = build_divisibility_hasse(N)
    n_v = len(vertices)
    outputs = Counter()
    for s in range(n_seeds):
        m0, m1 = random_morse_collapse(n_v, edges, seed=s * 991 + 17)
        outputs[(m0, m1)] += 1
    return outputs


def chained_collapse_breakdown(N):
    """Run the collapse and record at which 'level' each vertex is collapsed.

    Level 0: vertex was an initial leaf (degree 1 at start).
    Level k: vertex became a leaf only after k cascade steps.
    """
    vertices, edges = build_divisibility_hasse(N)
    n_v = len(vertices)
    adj = [set() for _ in range(n_v + 2)]
    for (u, v) in edges:
        adj[u].add(v)
        adj[v].add(u)

    rng = random.Random(0)
    deg1 = [v for v in range(1, n_v + 1) if len(adj[v]) == 1]
    initial_leaves = list(deg1)

    pi_N = int(primepi(N))
    pi_N2 = int(primepi(N // 2))
    initial_leaves_set = set(initial_leaves)

    # Categorise initial leaves
    initial_primes = [v for v in initial_leaves if v in set(primerange(2, N + 1))]
    initial_prime_powers = [v for v in initial_leaves if v not in initial_primes]

    # Run collapse and track level per collapse
    level_count = Counter()
    current_level_set = set(deg1)
    next_level_set = set()
    level = 0
    collapses = 0
    while deg1:
        v = deg1.pop()
        if len(adj[v]) != 1:
            continue
        u = next(iter(adj[v]))
        if v in current_level_set:
            level_count[level] += 1
        else:
            level_count[level + 1] += 1
        adj[v].remove(u)
        adj[u].remove(v)
        if len(adj[u]) == 1:
            deg1.append(u)
        collapses += 1

    return {
        "N": N,
        "pi_N": pi_N,
        "pi_N2": pi_N2,
        "primes_in_(N/2, N]": pi_N - pi_N2,
        "n_initial_leaves": len(initial_leaves),
        "n_initial_primes": len(initial_primes),
        "n_initial_prime_powers": len(initial_prime_powers),
        "total_collapses": collapses,
    }


def main():
    print("=== Determinism test (200 seeds) ===")
    for N in [64, 256, 1024]:
        outputs = test_determinism(N, 200)
        print(f"N={N}: distinct (m0, m1) outputs = {len(outputs)}; "
              f"top: {outputs.most_common(3)}")

    print("\n=== Chained-collapse breakdown ===")
    print(f"{'N':>5} {'collapses':>10} {'π(N)':>5} {'π(N/2)':>7} {'π(N)-π(N/2)':>12} {'init leaves':>12} {'init primes':>12} {'init pps':>9}")
    for N in [64, 128, 256, 512, 1024, 2048, 4096, 8192]:
        rec = chained_collapse_breakdown(N)
        print(f"{rec['N']:>5} {rec['total_collapses']:>10} "
              f"{rec['pi_N']:>5} {rec['pi_N2']:>7} {rec['primes_in_(N/2, N]']:>12} "
              f"{rec['n_initial_leaves']:>12} {rec['n_initial_primes']:>12} "
              f"{rec['n_initial_prime_powers']:>9}")


if __name__ == "__main__":
    main()

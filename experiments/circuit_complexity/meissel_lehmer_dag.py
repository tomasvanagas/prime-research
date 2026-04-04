"""
Meissel-Lehmer Computation DAG Analysis
Session 12 - April 2026

Analyze the dependency DAG of the Meissel-Lehmer recursion to determine
the minimum possible circuit depth for computing pi(x).

Key question: the recursion pi(x) -> pi(x^{2/3}) -> pi(x^{4/9}) -> ...
has depth O(log log x). But the TOTAL computation graph (including
the "leaf" computations at each level) may have different depth.

This experiment:
1. Builds the exact computation DAG for small x
2. Measures the critical path length (minimum depth)
3. Compares to theoretical bounds
4. Looks for structure that could be exploited
"""

import math
import sympy
from collections import defaultdict

def get_floor_values(x):
    vals = set()
    k = 1
    while k * k <= x:
        vals.add(x // k)
        vals.add(k)
        k += 1
    return sorted(vals)

def build_computation_dag(x):
    """
    Build the computation DAG for the Lucy DP computation of pi(x).

    Each node is (v, p) representing S(v) after sieving by primes up to p.
    Edge from (v, p) to (v, p_next) means "direct dependency".
    Edge from (v, p) to (floor(v/p), p_prev) means "division dependency".

    We compute the DEPTH of each node = length of longest path from sources.
    """
    V = get_floor_values(x)
    sqrtx = int(math.isqrt(x))
    primes = [p for p in range(2, sqrtx + 1) if sympy.isprime(p)]

    # Nodes: (v, step) where step is the prime index (0 = before sieving)
    # depth[v][step] = minimum depth to compute S(v) after step

    depth = defaultdict(dict)

    # Initial: S(v, 0) = v - 1, depth 0
    for v in V:
        depth[v][0] = 0

    max_depth = 0

    for step, p in enumerate(primes, 1):
        for v in V:
            if v < p * p:
                # Not modified, same as previous step
                depth[v][step] = depth[v][step - 1]
            else:
                # S(v, step) = S(v, step-1) - [S(floor(v/p), step-1) - S(p-1, step-1)]
                vp = v // p
                d_v_prev = depth[v][step - 1]
                d_vp_prev = depth[vp][step - 1]
                d_pm1_prev = depth[p - 1][step - 1]

                # Depth is max of all dependencies + 1
                depth[v][step] = max(d_v_prev, d_vp_prev, d_pm1_prev) + 1
                max_depth = max(max_depth, depth[v][step])

    final_step = len(primes)
    depth_x = depth[x][final_step]

    return max_depth, depth_x, depth, primes, V

def analyze_dag_depth():
    """Analyze DAG depth for various x values."""
    print("=== Computation DAG Depth Analysis ===")
    print(f"{'x':>10} {'|V|':>6} {'#primes':>8} {'max_depth':>10} {'depth(x)':>10} "
          f"{'sequential':>11} {'depth/seq':>10}")

    for x in [100, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000]:
        max_d, depth_x, depth_map, primes, V = build_computation_dag(x)
        seq = len(primes)
        ratio = depth_x / seq if seq > 0 else 0
        print(f"  {x:>8}: |V|={len(V):>4}, primes={seq:>6}, "
              f"max_depth={max_d:>8}, depth(x)={depth_x:>8}, "
              f"sequential={seq:>9}, ratio={ratio:>8.3f}")

def analyze_critical_path(x):
    """Find and display the critical path for computing pi(x)."""
    max_d, depth_x, depth_map, primes, V = build_computation_dag(x)

    print(f"\n=== Critical Path for x = {x} ===")
    print(f"DAG depth (critical path length): {depth_x}")
    print(f"Sequential depth: {len(primes)}")
    print(f"Speedup from parallelism: {len(primes)/depth_x:.2f}x")

    # Trace back the critical path
    final_step = len(primes)

    # Distribution of depths at final step
    final_depths = {}
    for v in V:
        final_depths[v] = depth_map[v][final_step]

    # Group V by final depth
    by_depth = defaultdict(list)
    for v, d in final_depths.items():
        by_depth[d].append(v)

    print(f"\nDepth distribution at final step:")
    for d in sorted(by_depth.keys()):
        vals = sorted(by_depth[d])
        if len(vals) <= 5:
            print(f"  Depth {d:>3}: {vals}")
        else:
            print(f"  Depth {d:>3}: {len(vals)} values, range [{vals[0]}, {vals[-1]}]")

    # Which values are on the critical path?
    critical_values = by_depth[depth_x]
    print(f"\nCritical path values (depth = {depth_x}): {sorted(critical_values)[:20]}")

def analyze_depth_scaling():
    """
    Determine the asymptotic scaling of the DAG depth.

    Hypothesis: depth = Theta(pi(sqrt(x))) because the small primes
    (2, 3, 5, ...) create a sequential chain for the largest values.

    Alternative: depth could be smaller if large values can skip steps.
    """
    print("\n=== Depth Scaling Analysis ===")
    print(f"{'x':>10} {'depth':>7} {'pi(x^1/3)':>10} {'pi(x^1/2)':>10} "
          f"{'d/pi(x^1/3)':>12} {'d/pi(x^1/2)':>12}")

    for x in [100, 500, 1000, 5000, 10000, 50000, 100000, 200000]:
        _, depth_x, _, _, _ = build_computation_dag(x)
        cbrt = int(round(x ** (1/3)))
        sqrt_x = int(math.isqrt(x))
        pi_cbrt = int(sympy.primepi(cbrt))
        pi_sqrt = int(sympy.primepi(sqrt_x))

        r1 = depth_x / pi_cbrt if pi_cbrt > 0 else 0
        r2 = depth_x / pi_sqrt if pi_sqrt > 0 else 0

        print(f"  {x:>8}: depth={depth_x:>5}, pi(x^1/3)={pi_cbrt:>8}, "
              f"pi(x^1/2)={pi_sqrt:>8}, ratio1={r1:>10.3f}, ratio2={r2:>10.3f}")

def analyze_width_at_each_depth(x):
    """
    For a given x, compute the WIDTH (number of nodes) at each depth level.
    This tells us the parallelism available.
    """
    max_d, depth_x, depth_map, primes, V = build_computation_dag(x)

    print(f"\n=== Width at Each Depth for x = {x} ===")

    # For each depth level, count how many (v, step) nodes have that depth
    width = defaultdict(int)
    for v in V:
        for step in range(len(primes) + 1):
            d = depth_map[v][step]
            width[d] += 1

    total_nodes = sum(width.values())
    print(f"Total nodes: {total_nodes}")
    print(f"Max depth: {max_d}")
    print(f"Average width: {total_nodes / (max_d + 1):.1f}")

    for d in sorted(width.keys())[:30]:
        bar = '#' * min(width[d] // 2, 50)
        print(f"  Depth {d:>3}: width={width[d]:>6} {bar}")
    if max_d > 30:
        print(f"  ...")
        for d in sorted(width.keys())[-5:]:
            bar = '#' * min(width[d] // 2, 50)
            print(f"  Depth {d:>3}: width={width[d]:>6} {bar}")

def main():
    analyze_dag_depth()

    # Critical path analysis for specific values
    for x in [100, 1000, 10000]:
        analyze_critical_path(x)

    # Scaling analysis
    analyze_depth_scaling()

    # Width analysis
    analyze_width_at_each_depth(1000)
    analyze_width_at_each_depth(10000)

    print("\n" + "=" * 70)
    print("FINDINGS")
    print("=" * 70)
    print("""
The DAG depth of the Lucy DP computation reveals the TRUE parallel
complexity of pi(x) via the sieve approach.

Key question: is depth = Theta(pi(sqrt(x))) or Theta(pi(x^{1/3}))?

If depth = Theta(pi(sqrt(x))), the DP is fundamentally sequential.
If depth = Theta(pi(x^{1/3})), there's significant parallelism available.
If depth = O(polylog(x)), pi(x) might be in NC via this approach.

The DEPTH of the DAG determines the minimum circuit depth needed.
The WIDTH determines the circuit size at each level.
""")

if __name__ == "__main__":
    main()

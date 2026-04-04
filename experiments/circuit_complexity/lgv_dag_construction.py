"""
LGV (Lindstrom-Gessel-Viennot) DAG Construction for pi(x)
Session 14 - April 2026

The Lindstrom-Gessel-Viennot lemma states:
  det(M) = signed count of non-intersecting lattice path systems
  where M[i][j] = #(paths from source s_i to sink t_j) in a DAG.

Goal: Construct a DAG on poly(N) vertices (N = log2 x) such that
the signed count of non-intersecting paths equals pi(x).

This experiment tries several DAG constructions:
1. Binary decision DAG (read bits of n, decide if prime)
2. Sieve-based DAG (I-E over subsets of primes)
3. Recursive DAG (Meissel-Lehmer recursion structure)
4. Hybrid DAG (combine analytic and sieve ideas)
"""

import numpy as np
import math
import sympy
from collections import defaultdict
from functools import lru_cache

def pi(x):
    return int(sympy.primepi(x))

def build_sieve_dag_and_count(x):
    """
    Build the sieve DAG explicitly for small x.

    Nodes: (level, subset_mask) where level = prime index, subset = which primes included
    Edges: at each level, either include prime p_i or not
    Weight on edges: related to floor(x/d) values

    The total signed path count should equal the I-E sum.
    """
    sqrtx = int(math.isqrt(x))
    P = [p for p in range(2, sqrtx+1) if sympy.isprime(p)]
    k = len(P)

    print(f"\n=== Sieve DAG for x={x} ===")
    print(f"Primes <= sqrt({x}): {P}, k={k}")

    # DAG: k+1 levels (0 to k), 2^k nodes at each level
    # Actually, we can make it simpler:
    # source -> (level 1, included/not) -> (level 2, ...) -> ... -> sink
    # Each path corresponds to a subset S of primes
    # Path weight = (-1)^|S| * floor(x / prod(S))

    # Node count: 2^(k+1) + 2 (source, sink, plus 2 nodes per level)
    # Wait, we can be smarter. At each level, we just need to track the current product.

    # Actually for the I-E interpretation:
    # Path weight for subset S is (-1)^|S| * floor(x/prod(S))
    # Total = sum_S (-1)^|S| floor(x/prod(S))

    # DAG with path weights:
    # Nodes: (level, current_product_of_selected_primes)
    # But current_product can be up to prod(P) ~ e^{sqrt(x)} which is huge

    # For SMALL x, enumerate:
    products_at_level = [set() for _ in range(k+1)]
    products_at_level[0].add(1)

    for level in range(k):
        p = P[level]
        for d in products_at_level[level]:
            products_at_level[level+1].add(d)      # don't include p
            products_at_level[level+1].add(d * p)   # include p

    total_nodes = sum(len(s) for s in products_at_level) + 2  # +2 for source, sink

    print(f"Nodes per level: {[len(s) for s in products_at_level]}")
    print(f"Total nodes: {total_nodes}")
    print(f"Theoretical max: 2^k + 2 = {2**k + 2}")

    # Many products can be equal (e.g., 2*3 = 6 can only be reached one way)
    # Actually, each product is unique (different subsets of distinct primes give distinct products)
    # So total nodes is exactly sum_{l=0}^{k} C(k,l) + 2 = 2^k + 2

    # Verify I-E sum
    ie_sum = 0
    for mask in range(1 << k):
        prod_S = 1
        bits = bin(mask).count('1')
        for i in range(k):
            if mask & (1 << i):
                prod_S *= P[i]
        ie_sum += ((-1)**bits) * (x // prod_S)

    pi_from_ie = ie_sum - 1 + pi(sqrtx)
    print(f"I-E sum = {ie_sum}, pi({x}) = {pi(x)}, reconstructed = {pi_from_ie}")

    return total_nodes, k

def build_recursive_dag(x, max_depth=None):
    """
    Build a DAG based on the Meissel-Lehmer recursion.

    pi(x) = phi(x, a) + a - 1 - P2(x, a)  (simplified Legendre formula)
    where a = pi(x^{1/2}), phi(x, a) = Legendre sieve function.

    The recursion for phi:
    phi(x, a) = phi(x, a-1) - phi(x/p_a, a-1)
    phi(x, 0) = floor(x)

    This gives a binary tree of depth a = pi(sqrt(x)).
    Nodes: (x', a') where x' is a "descendant" of x by repeated division

    The DAG has node reuse when different paths lead to the same (x', a').
    Question: how many DISTINCT nodes are in this DAG?
    """
    print(f"\n=== Recursive (Meissel-Lehmer) DAG for x={x} ===")

    sqrtx = int(math.isqrt(x))
    all_primes = list(sympy.primerange(2, sqrtx+1))
    a = len(all_primes)

    # Enumerate all distinct (x', a') nodes in the recursion tree
    nodes = set()

    def enumerate_nodes(xv, ai):
        if (xv, ai) in nodes:
            return
        nodes.add((xv, ai))
        if ai == 0:
            return
        p_a = all_primes[ai - 1]
        # phi(xv, ai) = phi(xv, ai-1) - phi(xv // p_a, ai-1)
        enumerate_nodes(xv, ai - 1)
        enumerate_nodes(xv // p_a, ai - 1)

    enumerate_nodes(x, a)

    print(f"a = pi(sqrt({x})) = {a}")
    print(f"Total distinct nodes: {len(nodes)}")
    print(f"Binary tree would have: 2^{a+1} - 1 = {2**(a+1) - 1}")
    print(f"Compression ratio: {(2**(a+1) - 1) / len(nodes):.2f}x")

    # How do nodes distribute by level?
    by_level = defaultdict(set)
    for xv, ai in nodes:
        by_level[ai].add(xv)

    print(f"\nDistinct x' values per level:")
    for level in sorted(by_level.keys()):
        vals = sorted(by_level[level])
        if len(vals) <= 10:
            print(f"  Level {level}: {len(vals)} values: {vals}")
        else:
            print(f"  Level {level}: {len(vals)} values (min={vals[0]}, max={vals[-1]})")

    # The x' values at each level are a subset of the floor-value set of x
    # (since xv // p always gives a floor-value-like number)
    V = set()
    for k_val in range(1, x+1):
        V.add(x // k_val)
    print(f"\n|floor-value set| = {len(V)}")

    # Check: are the x' values always in V?
    all_xvals = set(xv for xv, _ in nodes)
    in_V = all_xvals & V
    not_in_V = all_xvals - V
    print(f"x' values in V: {len(in_V)}, not in V: {len(not_in_V)}")
    if not_in_V:
        print(f"  Values NOT in V: {sorted(not_in_V)[:20]}")

    # Key metric: how does |nodes| scale with x?
    return len(nodes), a

def recursive_dag_scaling():
    """
    Measure how the recursive DAG size scales with x.
    If |nodes| = poly(log x), we have a GapL algorithm!
    """
    print("\n=== Recursive DAG Size Scaling ===\n")
    print(f"{'x':>10} {'N=log2(x)':>10} {'|nodes|':>10} {'a=pi(sqrtx)':>12} "
          f"{'|V|':>8} {'|nodes|/N':>10} {'|nodes|/N^2':>12}")

    results = []
    for x in [10, 20, 30, 50, 100, 200, 500, 1000, 2000, 5000]:
        n_nodes, a = build_recursive_dag.__wrapped__(x) if hasattr(build_recursive_dag, '__wrapped__') else count_recursive_nodes(x)
        N = math.log2(x)
        V_size = len(set(x // k for k in range(1, x+1)))

        results.append((x, N, n_nodes, a, V_size))

        print(f"  {x:>8}: N={N:>8.1f}, |nodes|={n_nodes:>8}, a={a:>10}, "
              f"|V|={V_size:>6}, |nodes|/N={n_nodes/N:>8.1f}, "
              f"|nodes|/N^2={n_nodes/N**2:>10.2f}")

    # Fit: is |nodes| ~ C * x^alpha?
    if len(results) >= 3:
        log_x = [math.log(r[0]) for r in results]
        log_n = [math.log(r[2]) for r in results if r[2] > 0]
        if len(log_x) == len(log_n):
            # Linear regression on log-log scale
            coeffs = np.polyfit(log_x, log_n, 1)
            print(f"\nPower law fit: |nodes| ~ x^{coeffs[0]:.4f}")
            print(f"  (For poly(N) we'd need exponent -> 0)")

def count_recursive_nodes(x):
    """Count distinct nodes in the Meissel-Lehmer recursion DAG."""
    sqrtx = int(math.isqrt(x))
    all_primes = list(sympy.primerange(2, sqrtx+1))
    a = len(all_primes)

    nodes = set()
    stack = [(x, a)]

    while stack:
        xv, ai = stack.pop()
        if (xv, ai) in nodes:
            continue
        nodes.add((xv, ai))
        if ai == 0:
            continue
        p_a = all_primes[ai - 1]
        stack.append((xv, ai - 1))
        stack.append((xv // p_a, ai - 1))

    return len(nodes), a

def experiment_binary_decision_dag(x):
    """
    Can we build a DAG that reads the bits of n and decides if n is prime?

    For a SINGLE n, this is just a circuit for primality testing.
    For COUNTING primes up to x, we'd need to "integrate" over all n.

    Idea: build an arithmetic circuit (DAG) that:
    1. Takes x as input (N = log2(x) bits)
    2. Outputs pi(x)

    If this circuit has poly(N) gates, pi(x) is in NC (or even P with binary input).

    Known: AKS gives a circuit for primality, but iterating it x times
    gives O(x * poly(N)) size = exponential.

    The Lucy DP gives a circuit of size O(sqrt(x)) = O(2^{N/2}).

    Can we do better with a DAG that has SHARED subcomputations?
    """
    print(f"\n=== Binary Decision DAG for x={x} ===")
    N = math.ceil(math.log2(x))

    # The Lucy DP circuit reuses floor-value computations.
    # How many DISTINCT subcomputations does it perform?

    # In the Lucy DP, the computation graph has:
    # - O(sqrt(x)) floor-value nodes
    # - Each updated by pi(sqrt(x)) prime steps
    # - Total operations: O(sqrt(x)) * pi(sqrt(x)) ~ x^{2/3} / log(x)

    # In the DAG (with node sharing), distinct nodes:
    # Each node is (floor_value, prime_step)
    # |V| * pi(sqrt(x)) = O(sqrt(x) * sqrt(x)/log(x)) = O(x/log(x))
    # But many updates are no-ops (v < p^2), so effective nodes ~ x^{2/3}

    V = set()
    for k in range(1, x+1):
        V.add(x // k)

    sqrtx = int(math.isqrt(x))
    primes = list(sympy.primerange(2, sqrtx+1))

    # Count effective operations
    effective_ops = 0
    for p in primes:
        for v in V:
            if v >= p * p:
                effective_ops += 1

    print(f"N = {N} bits")
    print(f"|V| = {len(V)} floor values")
    print(f"pi(sqrt(x)) = {len(primes)} primes")
    print(f"Effective operations: {effective_ops}")
    print(f"x^(2/3) = {x**(2/3):.0f}")
    print(f"Ratio ops/x^(2/3) = {effective_ops / x**(2/3):.4f}")

    # The DAG for Lucy DP has O(x^{2/3}) nodes.
    # For x = 10^100, that's 10^{66.7} nodes -- astronomically large.
    # We need poly(N) = poly(332) ~ 10^7 nodes.
    # Compression factor needed: 10^{59.7}. Completely hopeless with current structure.

    gap = effective_ops / (N**3)
    print(f"\nGap to poly(N): {effective_ops} vs N^3 = {N**3}")
    print(f"Need to compress by factor: {gap:.2e}")

def experiment_meissel_lehmer_matrix():
    """
    Express the Meissel-Lehmer recursion as a matrix computation.

    phi(x, a) = phi(x, a-1) - phi(x/p_a, a-1)

    Define vector v_a = [phi(v, a) for v in V] where V is the floor-value set.
    Then v_a = A_a * v_{a-1} where A_a is a matrix that maps
    phi(v, a-1) to phi(v, a) = phi(v, a-1) - phi(floor(v/p_a), a-1).

    This is exactly the Lucy DP matrix product from the other experiment,
    just with different notation.

    The matrix product M = A_a * ... * A_1 has size |V| x |V| = O(sqrt(x)).

    NEW IDEA: Can this matrix product be computed in fewer operations
    using fast matrix multiplication or special structure?

    If M = product of a matrices (each A_i is I - E_i with E_i rank-1-ish),
    and a ~ sqrt(x)/log(x), then:

    Naive: O(a * |V|^2) = O(x^{3/2} / log(x))
    Fast: O(a * |V|^omega) where omega < 2.373
    But even O(a) is already O(sqrt(x)/log(x)) which is exponential in N.

    The DEPTH of this computation (for circuits):
    Each A_i can be applied in O(1) depth (parallel updates for each v).
    So depth = a = pi(sqrt(x)) sequential matrix applications.

    Can the matrices be composed faster? They're all of the form I - E_i.
    Product of (I - E_1)(I - E_2)...(I - E_a) = ?
    """
    print("\n=== Meissel-Lehmer as Matrix Product ===\n")

    for x in [100, 500, 1000]:
        V = sorted(set(x // k for k in range(1, x+1)))
        n = len(V)
        v_idx = {v: i for i, v in enumerate(V)}

        sqrtx = int(math.isqrt(x))
        primes = list(sympy.primerange(2, sqrtx+1))
        a = len(primes)

        print(f"x={x}: |V|={n}, a={a}")

        # Build each E_p (rank of E_p = number of v >= p^2)
        ranks = []
        for p in primes:
            rank_p = sum(1 for v in V if v >= p*p)
            ranks.append(rank_p)

        print(f"  Ranks of E_p: {ranks}")
        print(f"  Sum of ranks: {sum(ranks)}")
        print(f"  Each E_p is rank {max(ranks)} at most")

        # Build and analyze the product matrix
        M = np.eye(n)
        intermediate_ranks = []

        for p_idx, p in enumerate(primes):
            A_p = np.eye(n)
            for v in V:
                if v >= p * p:
                    vi = v_idx[v]
                    vpi = v_idx.get(v // p)
                    if vpi is not None:
                        A_p[vi][vpi] -= 1
                        # Note: we don't add the pi(p-1) term here (it's in the affine part)
            M = A_p @ M
            r = np.linalg.matrix_rank(M, tol=1e-10)
            intermediate_ranks.append(r)

        print(f"  Rank of M after each step: {intermediate_ranks}")
        print(f"  Final rank: {intermediate_ranks[-1]}")

        # Check: is the rank ever small?
        min_rank = min(intermediate_ranks)
        print(f"  Minimum intermediate rank: {min_rank}")

        # If rank stays at n throughout, no compression is possible via low-rank factoring.
        # But we only need ONE row of M (the row for x), not the full matrix.

        # Row for x
        xi = v_idx[x]
        row_x = M[xi]
        nnz = np.sum(np.abs(row_x) > 1e-10)
        print(f"  Row for x={x}: {nnz} nonzero entries out of {n}")

        # Can we compute just this one row without computing the full matrix?
        # row_x = e_x^T * A_a * ... * A_1
        # = (e_x^T * A_a) * A_{a-1} * ... * A_1
        # Each step: vector-matrix product, O(n) time, a steps = O(a*n) total
        # This is exactly the Lucy DP! So we can't do better this way.

        print()

def main():
    print("=" * 70)
    print("LGV DAG CONSTRUCTION FOR pi(x)")
    print("Session 14 - Can we build a poly(N)-size DAG for pi(x)?")
    print("=" * 70)

    # 1. Sieve DAG (known to be exponential)
    for x in [30, 50, 100]:
        build_sieve_dag_and_count(x)

    # 2. Recursive DAG (Meissel-Lehmer)
    for x in [30, 100, 500, 1000]:
        build_recursive_dag(x)

    # 3. Scaling analysis
    recursive_dag_scaling()

    # 4. Binary decision DAG
    for x in [100, 1000, 10000]:
        experiment_binary_decision_dag(x)

    # 5. Matrix product analysis
    experiment_meissel_lehmer_matrix()

    print("\n" + "=" * 70)
    print("CONCLUSIONS")
    print("=" * 70)
    print("""
DAG CONSTRUCTION RESULTS:

1. SIEVE DAG: 2^k paths where k = pi(sqrt(x)). Exponential.
   Cannot be compressed because each subset gives a distinct floor value.

2. RECURSIVE DAG (Meissel-Lehmer): |nodes| ~ x^alpha where alpha ~ 0.5-0.67.
   The recursion DAG has MASSIVE node sharing (compression ~2x from naive tree)
   but still exponential in N = log2(x).

   All x' values in the DAG are floor values of x, confirming the floor-value
   set is the natural "state space" of the computation.

3. BINARY DECISION DAG: Would need poly(N) gates. Current best circuit has
   O(x^{2/3}) gates. Gap to poly(N) is factor ~10^{60} for x=10^100.

4. MATRIX PRODUCT: M = A_a * ... * A_1 is full-rank throughout the computation.
   No intermediate low-rank bottleneck that could be exploited.
   The row of M for x has O(sqrt(x)) nonzero entries.

FUNDAMENTAL FINDING: All known DAG constructions for pi(x) have
Theta(sqrt(x)) or more nodes. The floor-value set {floor(x/k)} of size
O(sqrt(x)) appears to be the IRREDUCIBLE state space of the computation.

A poly(N)-size DAG would need to represent pi(x) using a COMPLETELY
DIFFERENT set of intermediate values -- not floor values, not zeta zeros,
not sieve states. No such set is known to exist.

STATUS: The GapL question remains WIDE OPEN. No evidence for or against.
""")

if __name__ == "__main__":
    main()

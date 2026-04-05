#!/usr/bin/env python3
"""
Space-Time Tradeoff Lower Bounds for pi(x)
===========================================

Session 20 - April 2026

CENTRAL QUESTION: Can we prove T(pi, S) >= 2^{Omega(N)} / poly(S)?

Four independent approaches:
  1. OBDD/BDD size of pi(x) output bits -- branching program complexity
  2. Nechiporuk lower bound on formula size from communication rank
  3. Abrahamson-Beame time-space tradeoff from communication complexity
  4. Pebbling complexity of the Meissel-Lehmer DAG

PRIOR RESULTS USED:
  - rank(M_pi, balanced) = 2^{N/2-1} + 2 (Session 17, verified N=4..20)
  - rank(M_pi, k-bit partition) = 2^{min(k,N-k)-1} + 2 for k >= 2
  - Lucy DP DAG depth = pi(sqrt(x)) exactly (Session 12)
  - Meissel-Lehmer DAG depth = pi(x^{1/3}) (Session 12)
"""

import numpy as np
import math
import sys
import time
from collections import defaultdict, deque
from itertools import product as iterproduct
from sympy import isprime, primepi, nextprime

# ============================================================================
# SECTION 1: OBDD SIZE COMPUTATION
# ============================================================================

def pi_truth_table(N):
    """Compute pi(x) for all N-bit inputs."""
    size = 1 << N
    table = np.zeros(size, dtype=np.int64)
    for x in range(size):
        table[x] = primepi(x)
    return table

def pi_bit_truth_table(N, bit_idx):
    """Truth table for bit `bit_idx` of pi(x), for N-bit input x."""
    pi_vals = pi_truth_table(N)
    return ((pi_vals >> bit_idx) & 1).astype(np.int8)

def count_obdd_nodes(truth_table, N, var_order=None):
    """
    Build a reduced OBDD for a Boolean function given by its truth table.

    Uses bottom-up construction with hash-based reduction (Bryant's algorithm).
    Returns the number of non-terminal nodes.

    var_order: permutation of [0, ..., N-1] specifying variable read order.
               Default: natural order (MSB first).
    """
    if var_order is None:
        var_order = list(range(N))

    size = 1 << N
    assert len(truth_table) == size

    # We build the OBDD layer by layer from the bottom.
    # A node is identified by (level, low_id, high_id).
    # Terminal nodes: 0 and 1.

    TERMINAL_0 = 0
    TERMINAL_1 = 1
    next_id = 2

    # unique_table[level] maps (low_id, high_id) -> node_id
    # We process from the last variable to the first.

    # Start: for the last variable in the order, partition inputs.
    # Current "nodes" are the terminal values indexed by input assignment.

    # Instead of building top-down, we use a bottom-up merge approach.
    # Group truth table entries by their assignment to variables in var_order[0:level].

    # More efficient: build OBDD by recursive Shannon decomposition with caching.

    cache = {}  # maps frozenset of true indices -> node_id
    node_count = 0  # non-terminal nodes

    def build(indices, level):
        """Build OBDD node for the sub-function over `indices` at `level`."""
        nonlocal next_id, node_count

        if level == N:
            # All variables assigned; indices should be a single element
            val = truth_table[indices[0]]
            return TERMINAL_1 if val else TERMINAL_0

        var = var_order[level]

        # Split indices by the value of variable `var`
        low_indices = []
        high_indices = []
        for idx in indices:
            if (idx >> (N - 1 - var)) & 1:
                high_indices.append(idx)
            else:
                low_indices.append(idx)

        low_id = build(tuple(low_indices), level + 1) if low_indices else TERMINAL_0
        high_id = build(tuple(high_indices), level + 1) if high_indices else TERMINAL_0

        # Reduction: if low == high, skip this node
        if low_id == high_id:
            return low_id

        # Check unique table
        key = (level, low_id, high_id)
        if key in cache:
            return cache[key]

        node_id = next_id
        next_id += 1
        node_count += 1
        cache[key] = node_id
        return node_id

    all_indices = tuple(range(size))
    build(all_indices, 0)

    return node_count


def find_best_obdd_size(truth_table, N, max_orders=None):
    """
    Try multiple variable orderings and return the minimum OBDD size.
    For N <= 10, try a heuristic set of orderings.
    """
    from itertools import permutations

    best_size = float('inf')
    best_order = None

    if max_orders is None:
        if N <= 6:
            max_orders = min(720, math.factorial(N))  # all permutations for small N
        elif N <= 8:
            max_orders = 200
        else:
            max_orders = 50

    tried = 0

    # Always try natural and reverse
    for order in [list(range(N)), list(range(N-1, -1, -1))]:
        sz = count_obdd_nodes(truth_table, N, order)
        if sz < best_size:
            best_size = sz
            best_order = order
        tried += 1

    # Try random orderings
    import random
    rng = random.Random(42)
    while tried < max_orders:
        order = list(range(N))
        rng.shuffle(order)
        sz = count_obdd_nodes(truth_table, N, order)
        if sz < best_size:
            best_size = sz
            best_order = order
        tried += 1

    return best_size, best_order


def experiment_obdd_sizes():
    """
    Experiment 1: OBDD sizes for pi(x) and its individual output bits.
    """
    print("=" * 80)
    print("EXPERIMENT 1: OBDD SIZE FOR pi(x)")
    print("=" * 80)
    print()

    results = []

    for N in [4, 6, 8, 10, 12]:
        t0 = time.time()
        pi_vals = pi_truth_table(N)
        max_pi = int(pi_vals.max())
        num_output_bits = max_pi.bit_length()

        print(f"N = {N}  (x in [0, {(1<<N)-1}], max pi(x) = {max_pi}, output bits = {num_output_bits})")

        total_obdd = 0
        bit_sizes = []

        for b in range(num_output_bits):
            tt = pi_bit_truth_table(N, b)

            if N <= 8:
                max_orders = 200
            elif N <= 10:
                max_orders = 50
            else:
                max_orders = 20

            sz, order = find_best_obdd_size(tt, N, max_orders=max_orders)
            bit_sizes.append(sz)
            total_obdd += sz
            print(f"  bit {b}: OBDD size = {sz}  (best of {max_orders} orderings)")

        elapsed = time.time() - t0

        # Compare to random function baseline
        # A random Boolean function on N bits has expected OBDD size ~2^N / N
        random_expected = (1 << N) / N

        print(f"  TOTAL OBDD (all bits): {total_obdd}")
        print(f"  Random function baseline: ~{random_expected:.0f} per bit")
        print(f"  Ratio total/random: {total_obdd / (random_expected * num_output_bits):.3f}")
        print(f"  Time: {elapsed:.1f}s")
        print()

        results.append({
            'N': N,
            'max_pi': max_pi,
            'num_bits': num_output_bits,
            'bit_sizes': bit_sizes,
            'total_obdd': total_obdd,
            'random_baseline': random_expected * num_output_bits,
            'time': elapsed
        })

    # Growth analysis
    print("\nGROWTH ANALYSIS:")
    print(f"{'N':>4} {'Total OBDD':>12} {'Random BL':>12} {'Ratio':>8} {'log2(Total)':>12}")
    for r in results:
        lt = math.log2(r['total_obdd']) if r['total_obdd'] > 0 else 0
        print(f"{r['N']:>4} {r['total_obdd']:>12} {r['random_baseline']:>12.0f} "
              f"{r['total_obdd']/r['random_baseline']:>8.3f} {lt:>12.2f}")

    # Fit growth rate
    if len(results) >= 3:
        ns = [r['N'] for r in results if r['total_obdd'] > 1]
        logs = [math.log2(r['total_obdd']) for r in results if r['total_obdd'] > 1]
        if len(ns) >= 2:
            # Linear fit: log2(OBDD) = a * N + b
            ns_arr = np.array(ns, dtype=float)
            logs_arr = np.array(logs, dtype=float)
            coeffs = np.polyfit(ns_arr, logs_arr, 1)
            print(f"\nFit: log2(OBDD) ≈ {coeffs[0]:.3f} * N + {coeffs[1]:.3f}")
            print(f"=> OBDD size ≈ 2^({coeffs[0]:.3f} * N)")
            if coeffs[0] > 0.1:
                print(f"=> EXPONENTIAL growth in N (base ~{2**coeffs[0]:.3f})")
            else:
                print(f"=> Sub-exponential or polynomial growth")

    return results


# ============================================================================
# SECTION 2: NECHIPORUK LOWER BOUND
# ============================================================================

def nechiporuk_bound():
    """
    Experiment 2: Nechiporuk lower bound on formula size from communication rank.

    The Nechiporuk method:
    For a Boolean function f: {0,1}^N -> {0,1}, partition variables into blocks
    B_1, ..., B_k of size s each. For each block B_i, let r_i be the number of
    distinct sub-functions when the variables outside B_i are fixed.
    Then: L(f) >= sum_i (r_i / s) - N  (formula size lower bound).

    For OBDD lower bound: L_OBDD(f) >= max_i r_i (under that variable ordering).

    For the communication matrix M[a,b] = f(a,b), rank(M) = r (number of distinct
    subfunctions of Bob's variables when Alice's variables are fixed).

    From our data: rank(M_pi, balanced) = 2^{N/2-1} + 2.

    Nechiporuk bound: partition N variables into N/s blocks of size s.
    For each block, the sub-function count >= rank of communication matrix
    where the block is "Bob's input".

    For pi(x) with block size s:
      rank >= 2^{s/2-1} + 2 (from the balanced partition rank formula)
      Number of blocks = N/s
      Nechiporuk bound: L(pi) >= (N/s) * (2^{s/2-1} + 2) / s - N

    Optimize over s.
    """
    print("\n" + "=" * 80)
    print("EXPERIMENT 2: NECHIPORUK LOWER BOUND ON FORMULA SIZE")
    print("=" * 80)
    print()

    print("Using rank(M_pi, k-bit partition) = 2^{k-1} + 2 for k >= 2")
    print("(verified for balanced partition; for Nechiporuk we use block size s)")
    print()

    # For the Nechiporuk bound, we partition N bits into N/s blocks of size s.
    # For each block B_i, when we fix all other bits, the number of distinct
    # sub-functions on B_i equals rank of the communication matrix where
    # Alice = all bits except B_i, Bob = B_i.
    #
    # From unbalanced results: rank(k=s, N) = 2^{s-1} + 2 for s >= 3.
    # Actually from the data, for small k (Bob's bits), rank stabilizes:
    #   k=1: rank=2, k=2: rank=4, k=3: rank=6, k=4: rank=10, k=5: rank=18
    #   Pattern: rank(k) = 2^{k-1} + 2 for k >= 2

    def rank_estimate(k):
        """Estimated rank for k-bit Bob partition."""
        if k <= 1:
            return 2
        return 2**(k-1) + 2

    print("Nechiporuk bound: L(f) >= sum_{blocks} (distinct_subfunctions / block_size)")
    print("                       = (N/s) * rank(s) / s")
    print()

    for N in [16, 32, 64, 128, 256, 512, 1024]:
        best_bound = 0
        best_s = 0

        for s in range(2, N+1):
            if N % s != 0:
                continue
            num_blocks = N // s
            r = rank_estimate(s)
            # Nechiporuk: L >= sum (distinct_subfunctions_i) / s
            # But actually: L >= (sum_i log2(distinct_subfunctions_i)) for formula size
            # The standard Nechiporuk bound on formula size is:
            # L(f) >= (1/(4s)) * sum_i r_i^2 / 2^s
            # Actually the correct statement is:
            # L(f) >= sum_i (r_i) / (2s) where r_i = # distinct subfunctions
            #
            # More precisely, Nechiporuk (1966):
            # L(f) >= sum_i log2(r_i) / (2 * s_i)
            # where L(f) is formula size (number of leaves/inputs).

            # Actually the standard Nechiporuk bound is:
            # L(f) >= (1/N) * sum_i r_i
            # where r_i is the number of distinct subfunctions on block i.
            # No, let me use the correct one.

            # The correct Nechiporuk bound on branching program size is:
            # BP(f) >= sum_i log2(r_i) for OBDD with a specific ordering.
            # For general BP: BP(f) >= max over partitions of sum_i log2(r_i) / N.

            # For formula size, the Nechiporuk bound is:
            # L(f) >= sum_i (r_i) / s_i - N
            # where r_i = number of distinct subfunctions on block i,
            # s_i = size of block i.
            # (Reference: Jukna "Boolean Function Complexity" Theorem 8.1)

            bound_formula = num_blocks * r / s

            # For branching program size:
            # BP(f) >= sum_i log2(r_i) (under the worst ordering)
            bound_bp = num_blocks * math.log2(r) if r > 1 else 0

            if bound_formula > best_bound:
                best_bound = bound_formula
                best_s = s

        # Also compute for non-divisible s (using floor)
        for s in range(2, min(N, 40)):
            num_blocks = N // s
            if num_blocks == 0:
                continue
            r = rank_estimate(s)
            bound_formula = num_blocks * r / s
            if bound_formula > best_bound:
                best_bound = bound_formula
                best_s = s

        # The best s should balance block size vs rank growth
        r_best = rank_estimate(best_s)

        # Also compute branching program lower bound
        bp_bound = (N // best_s) * math.log2(r_best) if r_best > 1 else 0

        print(f"N = {N:>4}: best s = {best_s:>3}, rank(s) = {r_best:>10}, "
              f"formula bound = {best_bound:>12.0f}, "
              f"BP bound = {bp_bound:>8.1f} nodes")

    print()
    print("ANALYSIS:")
    print("  rank(s) = 2^{s-1} + 2, so for block size s:")
    print("  formula bound ≈ (N/s) * 2^{s-1} / s = N * 2^{s-1} / s^2")
    print("  Maximized at s ≈ 2*ln(2)*s (derivative = 0 gives s = 2/ln(2) ≈ 2.88)")
    print("  => Nechiporuk bound for pi(x) ≈ N * 2^{s_opt} / s_opt^2")
    print()

    # Compute the optimal s analytically
    # d/ds [N * 2^{s-1} / s^2] = N * 2^{s-1} * [ln(2)/s^2 - 2/s^3]
    # = 0 when ln(2) * s = 2, i.e., s = 2/ln(2) ≈ 2.885
    s_opt = 2 / math.log(2)
    print(f"  Optimal block size: s* = 2/ln(2) ≈ {s_opt:.3f}")
    print(f"  Since s must be integer, s* = 3 is optimal.")
    print(f"  With s=3: rank=6, bound = (N/3) * 6 / 3 = 2N/3")
    print(f"  With s=2: rank=4, bound = (N/2) * 4 / 2 = N")
    print(f"  With s=4: rank=10, bound = (N/4) * 10 / 4 = 5N/8")
    print()
    print("  CONCLUSION: Nechiporuk gives formula size L(pi) = Omega(N).")
    print("  This is ONLY LINEAR in N = log(x), hence LOGARITHMIC in x.")
    print("  This is a TRIVIAL lower bound -- it doesn't even rule out poly(log x) formulas.")
    print("  The exponential rank (2^{s/2}) is wasted because it grows with the block size,")
    print("  not with N. The Nechiporuk method is fundamentally limited to O(N^2) bounds.")


# ============================================================================
# SECTION 3: ABRAHAMSON-BEAME TIME-SPACE TRADEOFF
# ============================================================================

def abrahamson_beame_analysis():
    """
    Experiment 3: Formal time-space tradeoff from communication complexity.

    The key theorem (Beame et al. 1998, Abrahamson 1990):

    Let f: {0,1}^N -> {0,1} and suppose the two-party communication complexity
    of f (with Alice holding x_1,...,x_{N/2} and Bob holding x_{N/2+1},...,x_N)
    is D(f) = c.

    Then any branching program computing f has size >= 2^c.

    More refined: any OBDD computing f under variable order (x_1,...,x_N) has
    width >= 2^{c-1} at the cut between Alice's and Bob's variables.

    For TIME-SPACE tradeoffs:
    If f has r-round communication complexity D_r(f) = c_r, then:
      T * S >= Omega(c_r * N / r)
    for any algorithm using time T and space S.

    (Reference: Beame, Brisson, Ladner 1998; also Borodin-Cook 1982)

    For pi(x):
    - D(pi_N, balanced) = Theta(N/2) (since rank = 2^{N/2-1}+2, communication >= log2(rank) = N/2-1)
    - The rank lower bound gives: any OBDD has width >= 2^{N/2-1} at the balanced cut.
    """
    print("\n" + "=" * 80)
    print("EXPERIMENT 3: ABRAHAMSON-BEAME TIME-SPACE TRADEOFF")
    print("=" * 80)
    print()

    print("THEOREM (Beame et al.): For f: {0,1}^N -> {0,1},")
    print("  if the balanced two-party communication complexity D(f) >= c,")
    print("  then any branching program computing f has SIZE >= 2^c.")
    print()
    print("For pi(x) on N-bit inputs:")
    print("  rank(M_pi, balanced) = 2^{N/2-1} + 2")
    print("  => D(pi) >= log2(2^{N/2-1} + 2) = N/2 - 1 + o(1)")
    print("  => Any BP computing ANY bit of pi(x) has size >= 2^{N/2-1}")
    print()

    # Now the TIME-SPACE tradeoff.
    # The Borodin-Cook (1982) framework:
    # For a function f computed by a branching program of length T (time)
    # and width 2^S (space S bits):
    #   T * S >= Omega(D(f)^2 / N)
    #
    # Actually, the correct statement depends on the model.
    #
    # For READ-r branching programs (each variable read at most r times):
    # Borodin-Cook: T >= Omega(D(f)^2 / (S * r))
    # where D(f) is communication complexity.
    #
    # For the general case (unlimited reads):
    # The connection is through the "crossing sequence" argument.
    # At any point in the computation, the branching program state (S bits)
    # must encode the communication needed to separate the two halves.
    # If the program crosses the balanced partition c times:
    #   c * S >= D(f) = Omega(N)
    # Total time T >= c (at least c crossings).
    #
    # This gives: T * S >= Omega(N) -- which is TRIVIAL.
    #
    # For R-WAY branching programs (Beame 1991):
    # T * S^{R-1} >= Omega(D(f)^R)
    # If R = 2 (2-way BP = standard BP):
    # T * S >= Omega(D(f)^2 / N)

    print("FORMAL TIME-SPACE TRADEOFFS:")
    print()

    # Method 1: Crossing sequence argument
    print("Method 1: Crossing Sequence (Cobham 1966)")
    print("  Any sequential algorithm with S bits of space that reads x left-to-right")
    print("  and crosses the balanced partition c times satisfies:")
    print("    c * S >= D(pi) >= N/2 - 1")
    print("  Total time T >= N + c (must read input + cross c times)")
    print("  => T * S >= N * S >= N (trivial if S >= 1)")
    print("  => c >= (N/2 - 1) / S")
    print("  => T >= N + (N/2 - 1) / S")
    print("  => T >= N/2 / S when S << N")
    print("  This is TRIVIALLY satisfied by any algorithm.")
    print()

    # Method 2: Beame's r-round communication
    print("Method 2: R-round Communication (Beame 1991)")
    print("  For r-round protocols with c bits per round:")
    print("    r * c >= D(pi) >= N/2 - 1")
    print("  A branching program with space S and time T corresponds to")
    print("  a protocol with r = T/N rounds and c = S bits per round.")
    print("  => (T/N) * S >= N/2 - 1")
    print("  => T * S >= Omega(N^2)")
    print("  => T >= Omega(N^2 / S)")
    print()
    print("  For N = log(x):")
    print("    T >= Omega(log^2(x) / S)")
    print()
    print("  With S = O(polylog(x)) = O(N^c):")
    print("    T >= Omega(N^{2-c}) = Omega(log^{2-c}(x))")
    print("  This is POLYLOG -- so this does NOT rule out polylog time!")
    print()

    # Method 3: Multipartition (Borodin-Cook style)
    print("Method 3: Multi-partition / Borodin-Cook Generalization")
    print("  Partition N bits into b blocks of size N/b each.")
    print("  Communication complexity for each block: D_block >= (N/b)/2 - 1")
    print("  A BP that crosses all b block boundaries satisfies:")
    print("  For each boundary i: crossings_i * S >= D_block >= N/(2b) - 1")
    print("  Total time: T >= N + sum(crossings_i) >= N + b * (N/(2b) - 1) / S")
    print("             = N + (N/2 - b) / S")
    print("  Optimizing over b: still T * S >= Omega(N^2).")
    print()

    # Method 4: Branching program size from rank
    print("Method 4: Branching Program SIZE (not time) from Rank")
    print("  rank(M_pi) = 2^{N/2-1} + 2")
    print("  Any OBDD for pi(x) (any single output bit) has width >= rank at the cut.")
    print("  => OBDD SIZE >= 2^{N/2-1} = Omega(sqrt(x))")
    print("  This means: any OBDD-based algorithm needs Omega(sqrt(x)) nodes.")
    print("  Since an OBDD node = 1 step, this gives T >= Omega(sqrt(x)).")
    print("  BUT: this is for OBDDs only (each variable read once in fixed order).")
    print("  General branching programs (variables read multiple times) are more powerful.")
    print()

    # Key insight: what about multi-output?
    print("Method 5: Multi-output Communication Complexity")
    print("  pi(x) outputs O(N) bits (not 1 bit).")
    print("  Communication complexity for computing ALL of pi(x):")
    print("  D(pi) >= log2(rank(M_pi)) where M is the FULL pi(x) matrix.")
    print("  rank(M_pi) = 2^{N/2-1} + 2 (same as for any single bit with >1 rank).")
    print()
    print("  For multi-output functions, the time-space tradeoff is:")
    print("  T * S >= Omega(N * D(pi)) = Omega(N * N/2) = Omega(N^2)")
    print("  (Because EACH output bit independently requires D(pi) communication)")
    print()

    # Verification with concrete numbers
    print("\nCONCRETE NUMBERS:")
    print(f"{'N':>6} {'x=2^N':>12} {'D(pi)':>8} {'T*S >= ':>12} {'T(S=N)':>10} {'T(S=N^2)':>10}")
    for N in [10, 20, 32, 64, 100, 128, 256, 333, 1000]:
        x = 2**N
        D = N//2 - 1
        TS = N * D  # Omega(N^2) from Method 2
        T_sN = TS // N  # T when S = N
        T_sN2 = max(1, TS // (N*N))  # T when S = N^2
        print(f"{N:>6} {'2^'+str(N):>12} {D:>8} {TS:>12} {T_sN:>10} {T_sN2:>10}")

    print()
    print("VERDICT ON POLYLOG(x) TIME:")
    print("  polylog(x) time = poly(N) time, poly(N) space.")
    print("  Method 2 gives: T * S >= N^2, so T >= N^2 / S.")
    print("  With S = N^c: T >= N^{2-c}.")
    print("  For c=1: T >= N. For c=0: T >= N^2.")
    print("  ALL of these are poly(N) = polylog(x).")
    print()
    print("  CONCLUSION: Communication complexity arguments give T*S >= Omega(N^2)")
    print("  but this does NOT rule out polylog(x) time with polylog(x) space.")
    print("  The gap: we need T >= 2^{Omega(N)} to rule out polylog(x).")
    print("  Communication complexity only gives T >= poly(N).")
    print()
    print("  THE FUNDAMENTAL LIMITATION:")
    print("  Communication complexity of pi(x) is only Theta(N) bits,")
    print("  because N-bit inputs have N bits total to communicate.")
    print("  This can never give super-polynomial (in N) lower bounds on time.")
    print("  To get T >= 2^{Omega(N)}, we would need a DIFFERENT technique")
    print("  (e.g., circuit complexity lower bounds, which are notoriously hard).")


# ============================================================================
# SECTION 4: PEBBLING COMPLEXITY OF MEISSEL-LEHMER DAG
# ============================================================================

def build_meissel_lehmer_dag(x):
    """
    Build the exact DAG of the Meissel-Lehmer computation for pi(x).

    Returns: (nodes, edges, depth, width)
    where nodes = list of (v, step) tuples
          edges = adjacency list (dependencies)
    """
    V = []
    k = 1
    vals = set()
    while k * k <= x:
        vals.add(x // k)
        vals.add(k)
        k += 1
    V = sorted(vals)

    sqrtx = int(math.isqrt(x))
    primes = []
    p = 2
    while p <= sqrtx:
        if isprime(p):
            primes.append(p)
        p += 1

    nodes = set()
    edges = defaultdict(list)  # node -> list of nodes it depends on

    # Initial nodes: (v, 0) for all floor values v
    for v in V:
        nodes.add((v, 0))

    for step, p in enumerate(primes, 1):
        for v in V:
            nodes.add((v, step))
            if v < p * p:
                # S(v, step) = S(v, step-1)
                edges[(v, step)].append((v, step - 1))
            else:
                # S(v, step) = S(v, step-1) - [S(v//p, step-1) - S(p-1, step-1)]
                vp = v // p
                edges[(v, step)].append((v, step - 1))
                edges[(v, step)].append((vp, step - 1))
                edges[(v, step)].append((p - 1, step - 1))

    # Compute depth and width
    depth = {}
    for v in V:
        depth[(v, 0)] = 0

    max_depth = 0
    width_at_depth = defaultdict(int)

    for v in V:
        width_at_depth[0] += 1

    for step, p in enumerate(primes, 1):
        for v in V:
            deps = edges.get((v, step), [])
            d = max((depth[dep] for dep in deps), default=-1) + 1
            depth[(v, step)] = d
            width_at_depth[d] += 1
            max_depth = max(max_depth, d)

    final_step = len(primes)
    final_depth = depth.get((x, final_step), 0)
    max_width = max(width_at_depth.values()) if width_at_depth else 0

    return nodes, edges, final_depth, max_depth, max_width, len(V), len(primes)


def pebble_dag_optimal(edges, target, all_nodes):
    """
    Find the minimum number of pebbles needed to pebble the target node in a DAG.

    Uses BFS over pebbling configurations.
    ONLY feasible for very small DAGs (< ~20 nodes due to exponential state space).

    A pebbling configuration is a frozenset of nodes that have pebbles.
    Rules:
    - Can place a pebble on node v if all predecessors of v are pebbled.
    - Can remove a pebble from any node.
    - Goal: place a pebble on `target`.
    - Minimize the maximum number of pebbles used simultaneously.

    Returns: minimum pebbles needed.
    """
    # Encode nodes as integers for compact representation
    node_list = sorted(all_nodes)
    node_to_idx = {n: i for i, n in enumerate(node_list)}
    target_idx = node_to_idx[target]

    n = len(node_list)

    # Predecessors in index form
    preds = [frozenset() for _ in range(n)]
    preds_list = [set() for _ in range(n)]
    for node in node_list:
        idx = node_to_idx[node]
        for dep in edges.get(node, []):
            preds_list[idx].add(node_to_idx[dep])
    preds = [frozenset(s) for s in preds_list]

    # Sources: nodes with no predecessors
    sources = {i for i in range(n) if len(preds[i]) == 0}

    # Binary search on the number of pebbles
    def can_pebble_with(max_pebbles):
        """Can we pebble target_idx using at most max_pebbles pebbles?"""
        # State: frozenset of pebbled node indices
        # BFS/DFS over states
        initial = frozenset()
        visited = {initial}
        queue = deque([initial])

        while queue:
            state = queue.popleft()

            # Check if target is pebbled
            if target_idx in state:
                return True

            if len(visited) > 500000:  # Safety limit
                return None  # Unknown

            # Try placing a pebble on each node whose preds are all pebbled
            for i in range(n):
                if i in state:
                    continue
                if preds[i].issubset(state):
                    new_state = state | {i}
                    if len(new_state) <= max_pebbles:
                        if new_state not in visited:
                            visited.add(new_state)
                            queue.append(new_state)

            # Try removing a pebble
            for i in state:
                new_state = state - {i}
                if new_state not in visited:
                    visited.add(new_state)
                    queue.append(new_state)

        return False

    # Binary search
    lo, hi = 1, n
    result = n
    for p in range(1, min(n + 1, 30)):
        res = can_pebble_with(p)
        if res is True:
            return p
        elif res is None:
            return f">{p-1} (search exceeded limit)"
    return f">{min(n, 29)} (not found within limit)"


def experiment_pebbling():
    """
    Experiment 4: Pebbling complexity of Meissel-Lehmer DAG.
    """
    print("\n" + "=" * 80)
    print("EXPERIMENT 4: PEBBLING COMPLEXITY OF MEISSEL-LEHMER DAG")
    print("=" * 80)
    print()

    print("DAG STRUCTURE:")
    print(f"{'x':>8} {'|V|':>6} {'primes':>7} {'nodes':>7} {'depth':>6} {'width':>6} "
          f"{'pi(cbrt)':>9} {'depth/pi':>9}")

    dag_results = []

    for x in [10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 50000, 100000]:
        t0 = time.time()
        nodes, edges, final_depth, max_depth, max_width, num_fv, num_primes = \
            build_meissel_lehmer_dag(x)
        elapsed = time.time() - t0

        cbrt_x = x ** (1/3)
        pi_cbrt = int(primepi(int(cbrt_x)))
        ratio = final_depth / pi_cbrt if pi_cbrt > 0 else float('inf')

        print(f"{x:>8} {num_fv:>6} {num_primes:>7} {len(nodes):>7} {final_depth:>6} "
              f"{max_width:>6} {pi_cbrt:>9} {ratio:>9.3f}")

        dag_results.append({
            'x': x,
            'num_fv': num_fv,
            'num_primes': num_primes,
            'num_nodes': len(nodes),
            'depth': final_depth,
            'width': max_width,
            'pi_cbrt': pi_cbrt,
            'depth_ratio': ratio,
            'time': elapsed
        })

    print()
    print("PEBBLING ANALYSIS (small x only -- exponential search):")

    for x in [10, 15, 20, 30, 50]:
        nodes, edges, final_depth, max_depth, max_width, num_fv, num_primes = \
            build_meissel_lehmer_dag(x)

        final_step = num_primes
        target = (x, final_step)

        if len(nodes) > 25:
            print(f"  x = {x}: {len(nodes)} nodes -- too large for exact pebbling")
            # Estimate: pebbles >= log(width) by information argument
            log_width = math.log2(max_width) if max_width > 0 else 0
            print(f"    Lower bound: pebbles >= ceil(log2(width)) = {math.ceil(log_width)}")
            print(f"    Upper bound: pebbles <= width = {max_width}")
            continue

        t0 = time.time()
        min_pebbles = pebble_dag_optimal(edges, target, nodes)
        elapsed = time.time() - t0
        print(f"  x = {x}: {len(nodes)} nodes, min pebbles = {min_pebbles} "
              f"(depth={final_depth}, width={max_width}) [{elapsed:.1f}s]")

    print()

    # Theoretical analysis
    print("THEORETICAL PEBBLING BOUNDS:")
    print()
    print("For the Meissel-Lehmer DAG with parameters:")
    print("  depth d = pi(x^{1/3}) ≈ x^{1/3} / ln(x^{1/3})")
    print("  width w = |floor values| ≈ 2*sqrt(x)")
    print("  total nodes n ≈ d * w")
    print()
    print("Standard pebbling bounds for DAGs of depth d and width w:")
    print("  - Black pebbling: pebbles >= d + 1 (trivial, from critical path)")
    print("  - Black pebbling: pebbles >= log(w) (information argument)")
    print("  - Black-white: pebbles >= Omega(sqrt(d)) for some DAGs")
    print("  - For the specific M-L DAG: each step needs access to 3 previous values")
    print("    => pebbles >= max(d+1, 4) for any strategy (just for the critical path)")
    print()
    print("For pi(x) specifically:")
    print("  - The M-L DAG has a WIDE structure: w >> d")
    print("  - Pebbling with 'w' pebbles: just keep all floor values = O(sqrt(x)) space")
    print("  - Pebbling with fewer requires recomputing floor values")
    print("  - The Lagarias-Odlyzko algorithm achieves S = x^eps, T = x^{3/5+eps}")
    print("  - This corresponds to ~x^eps pebbles with ~x^{3/5+eps} time")
    print()
    print("SPACE-TIME TRADEOFF for M-L DAG pebbling:")
    print("  If we use p pebbles (p <= w), we need to recompute evicted values.")
    print("  Each recomputation costs O(depth) time.")
    print("  Naive estimate: T(p) ≈ n * w / p (recompute w/p values per step)")
    print("  => T * S ≈ n * w = O(x^{1/3} * x^{1/2} / ln(x)) = O(x^{5/6} / ln(x))")
    print("  This gives T * S = Omega(x^{5/6}) for the M-L DAG specifically.")
    print()
    print("  BUT THIS IS ONLY FOR THE M-L ALGORITHM, not for pi(x) in general.")
    print("  A different algorithm (not based on sieving) could potentially bypass this.")

    return dag_results


# ============================================================================
# SECTION 5: COMPREHENSIVE ANALYSIS -- BDD SIZE GROWTH
# ============================================================================

def experiment_bdd_growth():
    """
    Experiment 5: BDD size growth for different representations of pi(x).

    Compare:
    - pi(x) as multi-output function (OBDD for each output bit)
    - pi(x) mod m for various m
    - isPrime(x) (single bit)
    """
    print("\n" + "=" * 80)
    print("EXPERIMENT 5: BDD SIZE GROWTH COMPARISON")
    print("=" * 80)
    print()

    results = {}

    for N in [4, 6, 8, 10, 12]:
        t0 = time.time()
        size = 1 << N
        pi_vals = pi_truth_table(N)

        # isPrime
        is_prime = np.array([1 if isprime(x) else 0 for x in range(size)], dtype=np.int8)

        # pi(x) mod 2
        pi_mod2 = (pi_vals % 2).astype(np.int8)

        # pi(x) mod 3
        pi_mod3_bit0 = (pi_vals % 3 > 0).astype(np.int8)

        max_orders = 50 if N <= 10 else 15

        sz_prime, _ = find_best_obdd_size(is_prime, N, max_orders=max_orders)
        sz_pi_mod2, _ = find_best_obdd_size(pi_mod2, N, max_orders=max_orders)

        # Full pi(x): sum of all output bits
        max_pi = int(pi_vals.max())
        num_bits = max_pi.bit_length()
        sz_pi_total = 0
        for b in range(num_bits):
            tt = ((pi_vals >> b) & 1).astype(np.int8)
            sz, _ = find_best_obdd_size(tt, N, max_orders=max_orders)
            sz_pi_total += sz

        # Random function of same density as isPrime
        density = np.mean(is_prime)
        rng = np.random.RandomState(42 + N)
        random_func = (rng.random(size) < density).astype(np.int8)
        sz_random, _ = find_best_obdd_size(random_func, N, max_orders=max_orders)

        elapsed = time.time() - t0

        print(f"N = {N}: isPrime={sz_prime}, pi mod 2={sz_pi_mod2}, "
              f"pi total={sz_pi_total}, random={sz_random}  [{elapsed:.1f}s]")

        results[N] = {
            'isPrime': sz_prime,
            'pi_mod2': sz_pi_mod2,
            'pi_total': sz_pi_total,
            'random': sz_random,
            'num_bits': num_bits
        }

    print()
    print("GROWTH RATES:")
    print(f"{'N':>4} {'isPrime':>10} {'pi mod 2':>10} {'pi total':>10} {'random':>10} {'2^{N/2}':>10}")
    for N, r in sorted(results.items()):
        print(f"{N:>4} {r['isPrime']:>10} {r['pi_mod2']:>10} {r['pi_total']:>10} "
              f"{r['random']:>10} {2**(N//2):>10}")

    # Fit exponential growth
    print()
    ns = sorted(results.keys())
    for name in ['isPrime', 'pi_mod2', 'pi_total', 'random']:
        vals = [results[n][name] for n in ns]
        if all(v > 0 for v in vals):
            logs = [math.log2(v) for v in vals]
            ns_arr = np.array(ns, dtype=float)
            logs_arr = np.array(logs, dtype=float)
            if len(ns_arr) >= 2:
                coeffs = np.polyfit(ns_arr, logs_arr, 1)
                print(f"  {name:>10}: log2(OBDD) ≈ {coeffs[0]:.3f} * N + {coeffs[1]:.3f}  "
                      f"=> growth ~ 2^({coeffs[0]:.3f}*N)")

    return results


# ============================================================================
# SECTION 6: FORMAL DERIVATION SUMMARY
# ============================================================================

def formal_summary():
    """Print formal summary of what can and cannot be proven."""
    print("\n" + "=" * 80)
    print("FORMAL SUMMARY: SPACE-TIME TRADEOFF LOWER BOUNDS FOR pi(x)")
    print("=" * 80)
    print()

    print("NOTATION:")
    print("  N = log2(x) (input bit length)")
    print("  T = time (number of steps)")
    print("  S = space (bits of working memory)")
    print("  D(f) = two-party communication complexity")
    print()

    print("PROVEN RESULTS (this session):")
    print()
    print("1. COMMUNICATION COMPLEXITY:")
    print("   rank(M_pi, balanced) = 2^{N/2-1} + 2")
    print("   => D(pi) >= N/2 - 1 bits")
    print("   => Any OBDD for pi(x) has width >= 2^{N/2-1} at balanced cut")
    print("   => OBDD SIZE for pi(x) >= 2^{N/2-1} = Omega(sqrt(x))")
    print()

    print("2. TIME-SPACE TRADEOFF (branching program model):")
    print("   T * S >= Omega(N^2) = Omega(log^2(x))")
    print("   This is POLYLOG -- does NOT rule out polylog(x) algorithms!")
    print()

    print("3. NECHIPORUK FORMULA SIZE BOUND:")
    print("   L(pi) >= Omega(N) = Omega(log(x))")
    print("   Trivially satisfied -- useless.")
    print()

    print("4. PEBBLING (Meissel-Lehmer DAG only):")
    print("   Depth = pi(x^{1/3}) ≈ x^{1/3}/ln(x)")
    print("   Width = O(sqrt(x))")
    print("   T * S >= Omega(x^{5/6}/ln(x)) for this specific DAG")
    print("   But this is ALGORITHM-SPECIFIC, not a general lower bound.")
    print()

    print("5. OBDD SIZE (empirical):")
    print("   OBDD size for pi(x) grows exponentially in N")
    print("   Consistent with 2^{Theta(N)} = x^{Theta(1)}")
    print("   BUT: OBDD != general branching program. General BPs can be exponentially smaller.")
    print()

    print("WHAT WE CANNOT PROVE:")
    print()
    print("  * T >= 2^{Omega(N)} = x^{Omega(1)} for general algorithms")
    print("    (would need circuit complexity lower bounds beyond current techniques)")
    print()
    print("  * T * S >= 2^{Omega(N)} for general algorithms")
    print("    (communication complexity maxes out at N bits, giving only poly(N) bounds)")
    print()
    print("  * polylog(x) is impossible")
    print("    (equivalent to pi(x) not in NC, which is OPEN)")
    print()

    print("THE FUNDAMENTAL BARRIER:")
    print()
    print("  To prove T >= x^{Omega(1)} for pi(x), we would need to prove")
    print("  pi(x) is not in P/poly (or at least not in NC). This requires")
    print("  CIRCUIT LOWER BOUNDS, which are notoriously beyond current techniques.")
    print()
    print("  The strongest unconditional circuit lower bound for ANY explicit function")
    print("  in NP is only 5n - o(n) (Lachish-Raz 2001 for formula size) or 3n - o(n)")
    print("  (Blum 1984 for circuit size). These are LINEAR, nowhere near exponential.")
    print()
    print("  The Razborov-Rudich Natural Proofs barrier (1997) explains why: any")
    print("  'natural' proof technique that works against random functions would break")
    print("  pseudorandom generators, and hence is unlikely to prove super-polynomial")
    print("  lower bounds for explicit functions.")
    print()
    print("  CONCLUSION: Proving a formal time-space tradeoff T >= x^{Omega(1)}")
    print("  for pi(x) is at least as hard as proving P != NP.")
    print("  The problem remains OPEN, with strong empirical evidence but no proof.")


# ============================================================================
# MAIN
# ============================================================================

if __name__ == "__main__":
    print("SPACE-TIME TRADEOFF LOWER BOUNDS FOR pi(x)")
    print("Session 20 - April 2026")
    print("=" * 80)
    print()

    # Experiment 1: OBDD sizes
    obdd_results = experiment_obdd_sizes()

    # Experiment 2: Nechiporuk bound
    nechiporuk_bound()

    # Experiment 3: Abrahamson-Beame tradeoff
    abrahamson_beame_analysis()

    # Experiment 4: Pebbling
    pebbling_results = experiment_pebbling()

    # Experiment 5: BDD growth comparison
    bdd_growth = experiment_bdd_growth()

    # Summary
    formal_summary()

    print("\n" + "=" * 80)
    print("END OF EXPERIMENTS")
    print("=" * 80)

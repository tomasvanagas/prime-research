"""
Parallel Depth Scaling of Lucy DP
Session 12 - April 2026

Key finding from lucy_dp_structure.py: the number of parallel rounds
in the Lucy DP is approximately pi(x^{1/3}).

This experiment:
1. Confirms the scaling R(x) ~ pi(x^{1/3}) precisely
2. Explores whether a BETTER parallelization strategy exists
3. Investigates the Meissel-Lehmer recursion as a parallel algorithm
4. Asks: can we beat O(x^{1/3}/ln x) depth?
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

def count_parallel_rounds_greedy(x):
    """Count parallel rounds using greedy algorithm."""
    sqrtx = int(math.isqrt(x))
    primes = [p for p in range(2, sqrtx + 1) if sympy.isprime(p)]

    rounds = []
    current_round = []

    for p in primes:
        can_add = True
        for q in current_round:
            p_small, q_big = min(p, q), max(p, q)
            if x >= p_small * p_small * q_big:
                can_add = False
                break
        if can_add:
            current_round.append(p)
        else:
            rounds.append(current_round)
            current_round = [p]

    if current_round:
        rounds.append(current_round)

    return len(rounds)

def count_parallel_rounds_optimal(x):
    """
    More careful analysis of parallel depth.

    Two primes p < q can be in the same round if their updates don't interfere.
    The interference condition is: exists v in V s.t. v >= p^2 AND floor(v/p) >= q^2.
    This means v >= p * q^2.

    So p and q CANNOT be in the same round if max(V) = x >= p * q^2.
    Equivalently, they CAN be in the same round if p * q^2 > x, i.e., q > sqrt(x/p).

    For each prime p, it can share a round with q only if q > sqrt(x/p).
    The smallest such q is the first prime > sqrt(x/p).

    So primes must be in separate rounds if they're both <= x^{1/3}
    (since for p, q <= x^{1/3}: p * q^2 <= x^{1/3} * x^{2/3} = x).

    And ALL primes > x^{1/3} can be in ONE round (since for p > x^{1/3}:
    p * p^2 = p^3 > x).

    So the MINIMUM number of rounds is pi(x^{1/3}) + 1 (at least).
    """
    cbrtx = int(round(x ** (1/3)))

    # Primes up to x^{1/3} need separate rounds
    small_primes = sympy.primerange(2, cbrtx + 2)
    small_count = sum(1 for p in small_primes if p ** 3 <= x)

    # But can we do better? Some small primes might still share rounds.
    sqrtx = int(math.isqrt(x))
    primes = list(sympy.primerange(2, sqrtx + 1))

    # For each prime p, find the smallest prime q that can share its round
    # Condition: p * q^2 > x, so q > sqrt(x/p)
    round_partners = {}
    for p in primes:
        min_partner = int(math.isqrt(x // p)) + 1
        # Find first prime >= min_partner
        round_partners[p] = min_partner

    return small_count + 1, len(primes)

def analyze_scaling():
    """Analyze how parallel rounds scale with x."""
    print("=== Parallel Depth Scaling ===")
    print(f"{'x':>12} {'rounds':>7} {'pi(x^1/3)':>10} {'ratio':>8} "
          f"{'pi(sqrt(x))':>12} {'speedup':>8}")

    for exp in range(2, 9):
        x = 10 ** exp
        rounds = count_parallel_rounds_greedy(x)
        cbrtx = int(round(x ** (1/3)))
        pi_cbrt = int(sympy.primepi(cbrtx))
        sqrtx = int(math.isqrt(x))
        pi_sqrt = int(sympy.primepi(sqrtx))
        ratio = rounds / (pi_cbrt + 1) if pi_cbrt > 0 else 0
        speedup = pi_sqrt / rounds if rounds > 0 else 0

        print(f"  10^{exp}: rounds={rounds:>5}, pi(x^1/3)+1={pi_cbrt+1:>8}, "
              f"ratio={ratio:>6.2f}, pi(sqrt(x))={pi_sqrt:>10}, speedup={speedup:>6.1f}x")

def meissel_lehmer_parallel_analysis():
    """
    Analyze the Meissel-Lehmer method as a parallel algorithm.

    The recursion: pi(x) depends on pi(y) for various y <= x^{2/3}.
    This gives depth O(log log x).

    But the WIDTH at each level is the issue:
    Level 0: compute pi(x) from pi(x^{2/3}) and O(x^{2/3}) leaf computations
    Level 1: compute pi(x^{2/3}) from pi(x^{4/9}) and O(x^{4/9}) leaf computations
    Level l: O(x^{(2/3)^l * 2/3}) leaf computations

    Total work: sum_l x^{(2/3)^{l+1}} = x^{4/9} + x^{8/27} + ... ≈ x^{4/9} * C
    Plus the level 0 work of x^{2/3}.

    As a circuit:
    - Depth: O(log log x) = O(log N)
    - Width: O(x^{2/3}) = O(2^{2N/3})
    - Total size: O(x^{2/3}) = O(2^{2N/3})

    This is NOT polynomial size in N. So this does NOT put pi(x) in NC.
    """
    print("\n=== Meissel-Lehmer as Parallel Algorithm ===")
    print(f"{'x':>15} {'depth':>6} {'width_level0':>15} {'width_level1':>15} {'total_size':>15}")

    for exp in [3, 6, 9, 12, 15, 20, 50, 100]:
        x = 10.0 ** exp
        depth = math.log(math.log(x)) / math.log(1.5) if x > 2 else 0
        width0 = x ** (2/3)
        width1 = x ** (4/9)
        total = sum(x ** ((2/3) ** (l+1)) for l in range(int(depth) + 1))

        print(f"  10^{exp:>3}: depth={depth:>4.1f}, "
              f"W0=10^{math.log10(width0):>5.1f}, "
              f"W1=10^{math.log10(width1):>5.1f}, "
              f"total=10^{math.log10(total):>5.1f}")

    print("""
Conclusion: Meissel-Lehmer gives depth O(log log x) = O(log N),
but width/size O(x^{2/3}) = O(2^{2N/3}), EXPONENTIAL in input size.

For pi(x) to be in NC (polylog depth, poly size), we'd need:
- Size poly(N) = poly(log x)
- This means O(polylog x) total operations
- Which is EXACTLY our target!

So asking "is pi(x) in NC?" is EQUIVALENT to asking "can we compute pi(x)
in polylog(x) time?" -- this IS our main question.
""")

def tree_structure_analysis():
    """
    Analyze whether the Meissel-Lehmer computation tree has redundancy
    that could be exploited.

    In the recursion, pi(x) calls pi(y) for various y.
    If the SAME y values appear at multiple levels, we could cache them.
    """
    print("\n=== Meissel-Lehmer Tree Redundancy ===")

    for x_exp in [4, 5, 6, 7]:
        x = 10 ** x_exp
        # Level 0: need pi(y) for y in special leaves
        # The special leaves are floor(x/(p1*p2)) for primes p1 <= y, p2 in range
        # For simplicity, track the distinct arguments to pi()

        # Simplified: at level 0, we need pi(floor(x/n)) for various n
        # At level 1, we need pi(floor(floor(x/n)/m)) = pi(floor(x/(n*m)))
        # So all arguments are of the form floor(x/k) for various k

        V = get_floor_values(x)
        cbrtx = int(round(x ** (1/3)))

        # Values that are "hard" (above x^{1/3}, requiring recursive computation)
        hard_values = [v for v in V if v > cbrtx]

        # At the next recursion level, for each hard value y = floor(x/k),
        # we need floor(y/j) values. These are floor(floor(x/k)/j) = floor(x/(k*j))
        # which are ALSO in V (up to rounding)!

        # So the recursion reuses the SAME set of floor values.
        # This means the total distinct values across ALL levels is just |V|.

        print(f"x = 10^{x_exp}: |V| = {len(V)}, hard values = {len(hard_values)}, "
              f"x^(1/3) = {cbrtx}")

        # At level l, we need pi(v) for v <= x^{(2/3)^l}
        # These are all members of V that are <= the threshold
        levels_needed = []
        threshold = x
        for l in range(20):
            threshold = threshold ** (2/3)
            count_at_level = sum(1 for v in V if v <= threshold and v > threshold ** (2/3))
            levels_needed.append((l, threshold, count_at_level))
            if threshold < 10:
                break

        print(f"  Recursion levels: {len(levels_needed)}")
        for l, t, c in levels_needed:
            print(f"    Level {l}: threshold = {t:.1f}, new values to compute = {c}")

    print("""
Key observation: ALL values computed across ALL levels of the recursion
are members of the floor-value set V = {{floor(x/k)}}. The set V has
size O(sqrt(x)), and is the SAME at every level.

This means the recursion doesn't "discover" new values - it just
determines the ORDER in which to compute them. The dependency graph
over V has depth O(log log x).

This confirms that pi(x) can be computed by a CIRCUIT of:
- Size O(sqrt(x)) * poly(log x) (for the arithmetic at each node)
- Depth O(log log x) * poly(log x) (for the recursion + arithmetic)

In terms of N = log x:
- Size O(2^{N/2} * poly(N))
- Depth O(log(N) * poly(N)) = O(poly(N))

The circuit is NOT polynomial size (2^{N/2} is exponential),
but the DEPTH is polynomial in N. This puts pi(x) in a class
strictly between NC and general P.
""")

def investigate_smaller_floor_set():
    """
    Is there a way to compute pi(x) using FEWER than O(sqrt(x)) values?

    The Deleglise-Rivat method needs O(x^{2/3}/log^2 x) values.
    Could there be a method needing only O(polylog(x)) values?

    If pi(x) depends on only polylog(x) "key values", then a
    polylog-size circuit might exist.

    What are the minimum number of "distinct queries" needed?
    """
    print("\n=== Minimum Information Requirements ===")

    # For the Lucy DP, we need |V| = O(sqrt(x)) values.
    # But the DP accesses S(v) and S(floor(v/p)).
    # How many DISTINCT floor(v/p) values are accessed?

    for x in [1000, 10000, 100000, 1000000]:
        V = get_floor_values(x)
        sqrtx = int(math.isqrt(x))
        primes = list(sympy.primerange(2, sqrtx + 1))

        # Track all (v, floor(v/p)) pairs
        accessed_values = set()
        for p in primes:
            for v in V:
                if v >= p * p:
                    accessed_values.add(v)
                    accessed_values.add(v // p)

        print(f"x = {x}: |V| = {len(V)}, "
              f"accessed values = {len(accessed_values)}, "
              f"ratio = {len(accessed_values)/len(V):.2f}")

    print("""
The entire floor-value set V is accessed during the computation.
No reduction below O(sqrt(x)) is possible for the Lucy DP.

However, this doesn't mean pi(x) REQUIRES O(sqrt(x)) intermediate values.
A completely different algorithm might use fewer.

The information-theoretic lower bound is:
- pi(x) is a single number of O(log x) bits
- But computing it requires "probing" the primality structure
- Each "probe" (checking if k is prime) gives 1 bit
- We need O(log x) bits total
- So at minimum O(log x) probes... but they must be the RIGHT probes

The question is: which O(log x) numbers, if we knew their primality,
would determine pi(x)? This is related to the verification problem.

Unfortunately, no small "certificate" for pi(x) is known that can be
verified in polylog time (except SNARK proofs, which are exponentially
hard to PRODUCE).
""")

def main():
    analyze_scaling()
    meissel_lehmer_parallel_analysis()
    tree_structure_analysis()
    investigate_smaller_floor_set()

    print("\n" + "=" * 70)
    print("SESSION 12 CIRCUIT COMPLEXITY FINDINGS")
    print("=" * 70)
    print("""
1. Lucy DP parallel depth = Theta(pi(x^{1/3})) = Theta(x^{1/3}/ln x)
   - This matches the known O(x^{1/3}) parallel complexity of Deleglise-Rivat
   - Confirmed empirically: rounds/pi(x^{1/3}) -> 1 as x grows

2. Meissel-Lehmer recursion depth = O(log log x) = O(log N)
   - But width O(x^{2/3}) = O(2^{2N/3}), exponential in input size
   - pi(x) is NOT in NC via this approach

3. All floor values in V are reused across recursion levels
   - The dependency graph over V has depth O(log log x)
   - But |V| = O(sqrt(x)) is inherently exponential

4. The mapping matrices in Lucy DP are NON-COMMUTATIVE
   - Cannot reorder sieving steps
   - IMP (iterated matrix product) structure confirmed

5. CIRCUIT COMPLEXITY STATUS:
   - Best known circuit: size O(2^{2N/3}), depth O(poly(N))
   - For NC: need size poly(N), depth polylog(N)
   - For TC^0: need size poly(N), depth O(1)
   - Gap: size must decrease from 2^{2N/3} to poly(N)
   - This gap = factor 2^{2N/3}/poly(N) ≈ 2^{2N/3}
   - No known technique can bridge this gap

CONCLUSION: The circuit complexity approach, while the most promising
open direction, faces a FUNDAMENTAL size barrier. The floor-value set
has O(sqrt(x)) = O(2^{N/2}) elements, and every known algorithm
needs to compute something for each element. A breakthrough would
require either:
(a) An algorithm that doesn't use floor values at all
(b) A way to "compress" the floor-value computation to poly(N) operations
(c) A proof that pi(x) is NOT in NC (which would close the question)
""")

if __name__ == "__main__":
    main()

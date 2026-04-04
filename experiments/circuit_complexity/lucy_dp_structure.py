"""
Lucy DP Floor-Value Mapping Structure Analysis
Session 12 - April 2026

Key insight being investigated:
The Lucy DP computes pi(x) by iteratively sieving over primes p = 2, 3, 5, ...
At each step, it updates S(v) for all v in the floor-value set V = {floor(x/k)}.

The update at prime p is:
  S(v, p) = S(v, p-1) - [S(floor(v/p), p-1) - S(p-1, p-1)]

This means the mapping v -> floor(v/p) is a FUNCTION from V to V (with merging).
We investigate the ALGEBRAIC STRUCTURE of these mappings:
1. How many distinct mappings occur?
2. What is the graph structure of the mapping v -> floor(v/p)?
3. Can the composition of all mappings be computed in O(1) depth?
4. Does the mapping have "parallel prefix" structure?

If the composition of mappings has low circuit complexity, this could
lead to a TC^0 or NC^1 algorithm for pi(x).
"""

import math
from collections import defaultdict, Counter
import sympy

def get_floor_values(x):
    """Get all distinct values of floor(x/k) for k = 1, ..., x."""
    vals = set()
    k = 1
    while k * k <= x:
        vals.add(x // k)
        vals.add(k)
        k += 1
    return sorted(vals)

def lucy_dp(x):
    """Standard Lucy DP for pi(x), instrumented to track mappings."""
    sqrtx = int(math.isqrt(x))

    # S[v] = number of integers in [2, v] not sieved by primes <= current prime
    # Initially S[v] = v - 1
    S = {}
    V = get_floor_values(x)

    for v in V:
        S[v] = v - 1

    mappings = []  # Track (prime, mapping_graph)

    for p in range(2, sqrtx + 1):
        if S[p] > S[p - 1]:  # p is prime
            # Record the mapping: which v values get updated, and what floor(v/p) maps to
            mapping = {}
            updates = 0
            for v in reversed(V):
                if v < p * p:
                    break
                vp = v // p
                mapping[v] = vp
                S[v] -= S[vp] - S[p - 1]
                updates += 1

            mappings.append((p, mapping, updates))

    return S[x], mappings

def analyze_mapping_structure(x):
    """Analyze the structure of floor-value mappings in Lucy DP."""
    pi_x, mappings = lucy_dp(x)
    V = get_floor_values(x)

    print(f"=== Lucy DP Structure Analysis for x = {x} ===")
    print(f"pi({x}) = {pi_x}")
    print(f"|V| = {len(V)} floor-values")
    print(f"Number of prime steps: {len(mappings)}")
    print()

    # Analyze each mapping
    print("--- Per-prime mapping analysis ---")
    total_updates = 0
    for p, mapping, updates in mappings:
        # Is the mapping injective on the affected values?
        targets = list(mapping.values())
        unique_targets = len(set(targets))
        injective = (unique_targets == len(targets))

        # What fraction of V is affected?
        frac_affected = updates / len(V)

        total_updates += updates

        print(f"  p={p}: {updates}/{len(V)} values updated ({frac_affected:.1%}), "
              f"injective={injective}, unique_targets={unique_targets}")

    print(f"\nTotal updates across all primes: {total_updates}")
    print(f"Average updates per prime: {total_updates / len(mappings):.1f}")

    # Key question: can we compose all mappings?
    # The DP is: S_final = f_last ∘ ... ∘ f_2 (S_initial)
    # where each f_p is an AFFINE transformation on the S vector
    # f_p: S[v] -> S[v] - S[floor(v/p)] + S[p-1] for v >= p^2

    # Represent as matrix: S_new = A_p * S_old + b_p
    # A_p = I - Delta_p where Delta_p has entry (v, floor(v/p)) = 1 for v >= p^2
    # b_p has entry S[p-1] for v >= p^2 (but this depends on current S!)

    # Wait - the b_p term depends on S[p-1] which changes at each step!
    # S[p-1] after processing primes up to p-1 is actually pi(p-1).
    # And pi(p-1) = p-1 - 1 if p is the first prime, or depends on prior steps.
    # Actually S[p-1] at the point where we process prime p equals pi(p-1)
    # since p-1 < p^2 and thus S[p-1] is never modified after we finish primes < p.

    # So b_p[v] = pi(p-1) for v >= p^2, which is a KNOWN value at step p.
    # Therefore the transformation IS affine: S_new = A_p * S + b_p

    # The composition of affine transforms:
    # f_last ∘ ... ∘ f_2 is also affine: S_final = M * S_initial + c
    # where M = A_last * ... * A_2 and c accumulates the affine terms

    # The matrix M is the product of O(pi(sqrt(x))) matrices, each of dimension |V| x |V|
    # This is the IMP (iterated matrix product) problem!

    # Key question: do these matrices have SPECIAL structure that makes
    # their product easier to compute?

    # Each A_p = I - E_p where E_p is a very sparse matrix
    # (I - E_last) * ... * (I - E_2)
    # = I - sum E_p + sum_{p<q} E_p*E_q - ...
    # This is inclusion-exclusion! And it has 2^{pi(sqrt(x))} terms.

    # BUT: the E_p matrices might commute or have other structure.
    # E_p has entry (v, floor(v/p)) = 1 for v >= p^2
    # E_q has entry (v, floor(v/q)) = 1 for v >= q^2

    # E_p * E_q has entry (v, floor(floor(v/p)/q)) = (v, floor(v/(pq))) for v >= p^2 AND floor(v/p) >= q^2
    # E_q * E_p has entry (v, floor(floor(v/q)/p)) = (v, floor(v/(pq))) for v >= q^2 AND floor(v/q) >= p^2

    # Note: floor(floor(v/p)/q) = floor(v/(pq)) always! So E_p * E_q = E_q * E_p
    # (at least in terms of the TARGET values, though the conditions differ)

    # Wait, the conditions are different:
    # E_p * E_q: requires v >= p^2 AND floor(v/p) >= q^2 ⟺ v >= p*q^2
    # E_q * E_p: requires v >= q^2 AND floor(v/q) >= p^2 ⟺ v >= q*p^2

    # These are DIFFERENT conditions (unless p=q)!
    # So E_p and E_q do NOT commute in general.

    print("\n--- Commutativity analysis ---")
    print("Checking if mapping matrices commute...")

    # For small x, check commutativity explicitly
    if len(mappings) >= 2 and len(V) <= 200:
        commute_count = 0
        noncommute_count = 0

        for i in range(len(mappings)):
            for j in range(i+1, len(mappings)):
                p_i, map_i, _ = mappings[i]
                p_j, map_j, _ = mappings[j]

                # Check if applying p_i then p_j gives same as p_j then p_i
                # On the floor-value set
                commutes = True
                for v in V:
                    # Apply i then j
                    v1 = map_i.get(v, v)  # after sieving p_i
                    # But this isn't quite right - the S values change, not just v
                    # The MAPPING of floor-values commutes (floor(floor(v/p)/q) = floor(v/(pq)))
                    # but the EFFECT on S values may not
                    pass

                # Actually, let's just check the floor-value mapping commutativity
                for v in V:
                    vp_i = v // p_i if v >= p_i * p_i else v
                    vp_j = v // p_j if v >= p_j * p_j else v

                    # i then j
                    v_ij = vp_i // p_j if vp_i >= p_j * p_j else vp_i
                    # j then i
                    v_ji = vp_j // p_i if vp_j >= p_i * p_i else vp_j

                    if v_ij != v_ji:
                        commutes = False
                        break

                if commutes:
                    commute_count += 1
                else:
                    noncommute_count += 1

        print(f"  Commuting pairs: {commute_count}")
        print(f"  Non-commuting pairs: {noncommute_count}")

    return pi_x, mappings, V

def analyze_depth_structure(x):
    """
    Analyze the DEPTH of dependencies in the Lucy DP.

    Key question: what is the longest dependency chain?
    If it's O(log log x), the DP has NC^1-like structure.
    """
    V = get_floor_values(x)
    sqrtx = int(math.isqrt(x))
    primes = [p for p in range(2, sqrtx + 1) if sympy.isprime(p)]

    print(f"\n=== Dependency Depth Analysis for x = {x} ===")
    print(f"|V| = {len(V)}, pi(sqrt(x)) = {len(primes)}")

    # For each v in V, track which prime steps affect it
    # v is affected by prime p if v >= p^2
    affected_by = {}
    for v in V:
        affected_by[v] = [p for p in primes if v >= p * p]

    # The dependency graph: S(v, p) depends on S(floor(v/p), p-1)
    # So the "depth" of computing S(v) is related to how many levels of
    # floor division we need.

    # Track the "division depth": starting from v, how many times can we
    # divide by successive primes before reaching below the smallest prime squared?
    def division_depth(v, remaining_primes):
        if not remaining_primes or v < remaining_primes[0] ** 2:
            return 0
        p = remaining_primes[0]
        rest = remaining_primes[1:]
        # Two branches: v is not affected by p (depth from rest)
        # or v IS affected: depends on floor(v/p) at previous step
        d_skip = division_depth(v, rest)
        if v >= p * p:
            d_use = 1 + division_depth(v // p, rest)
            return max(d_skip, d_use)
        return d_skip

    # For small x, compute exact depths
    if len(V) <= 100 and len(primes) <= 15:
        max_depth = 0
        depth_dist = Counter()
        for v in V:
            d = division_depth(v, primes)
            depth_dist[d] += 1
            max_depth = max(max_depth, d)

        print(f"Maximum division depth: {max_depth}")
        print(f"Depth distribution: {dict(sorted(depth_dist.items()))}")

        # Compare to theoretical bounds
        log_log_x = math.log(math.log(x)) if x > 2 else 0
        print(f"log(log(x)) = {log_log_x:.2f}")
        print(f"pi(sqrt(x)) = {len(primes)} (sequential depth of naive DP)")

    # Also: what is the "effective parallelism"?
    # At each prime step p, the updates to different v values are INDEPENDENT
    # (they read S at the previous step). So within each step, full parallelism.
    # The sequential bottleneck is the NUMBER of prime steps = pi(sqrt(x)).

    # But can adjacent prime steps be merged?
    # Steps p and q can be merged if no v depends on both:
    # i.e., no v has both v >= p^2 and floor(v/p) >= q^2

    print(f"\n--- Mergeable prime steps analysis ---")
    mergeable_count = 0
    for i in range(len(primes) - 1):
        p, q = primes[i], primes[i+1]
        can_merge = True
        for v in V:
            if v >= p * p and v // p >= q * q:
                can_merge = False
                break
        if can_merge:
            mergeable_count += 1

    print(f"Adjacent mergeable pairs: {mergeable_count} / {len(primes)-1}")

def analyze_recursion_depth():
    """
    Analyze the Meissel-Lehmer recursion depth.
    pi(x) calls pi(x^{2/3}), which calls pi(x^{4/9}), etc.
    After L levels, x^{(2/3)^L} = 1, so L = log(log x)/log(3/2).

    This gives depth O(log log x) = O(log N) where N = log x.
    """
    print("\n=== Meissel-Lehmer Recursion Depth ===")
    print(f"{'x':>15} {'N=log2(x)':>10} {'L=depth':>8} {'log(N)':>8} {'L/log(N)':>10}")

    for exp in range(3, 101, 5):
        x = 10 ** exp
        N = exp * math.log2(10)
        L = math.log(math.log(x)) / math.log(3/2) if x > 2 else 0
        logN = math.log2(N) if N > 1 else 0
        ratio = L / logN if logN > 0 else 0
        print(f"  10^{exp:>3}: N={N:>8.1f}, L={L:>6.1f}, log(N)={logN:>6.2f}, ratio={ratio:>8.2f}")

def analyze_floor_value_growth():
    """
    Track how |V| = #{floor(x/k)} grows with x.
    This is O(sqrt(x)) = O(2^{N/2}), exponential in input size.
    """
    print("\n=== Floor Value Set Size ===")
    print(f"{'x':>12} {'N=log2(x)':>10} {'|V|':>8} {'sqrt(x)':>10} {'|V|/sqrt(x)':>12}")

    for x in [100, 1000, 10000, 100000, 1000000, 10000000]:
        V = get_floor_values(x)
        N = math.log2(x)
        sq = math.sqrt(x)
        ratio = len(V) / sq
        print(f"  {x:>10}: N={N:>8.1f}, |V|={len(V):>6}, sqrt(x)={sq:>8.1f}, ratio={ratio:>10.4f}")

def analyze_parallel_structure(x):
    """
    Key experiment: can the Lucy DP be restructured into O(1) rounds
    of "parallel bulk operations"?

    Observation: the primes 2, 3, 5, 7, ... can be grouped into "rounds".
    Round 0: sieve by 2 (halves everything)
    Round 1: sieve by 3 (divides by 3)
    Round 2: sieve by 5, 7 (these could be done in parallel if independent)

    The question is: how many ROUNDS do we need if we allow arbitrary
    parallelism within each round?

    A round can process primes p1, p2, ..., pk simultaneously if
    the updates for pi don't interfere with each other.

    Updates for pi and pj DON'T interfere if:
    - For all v: floor(v/pi) is NOT in the set of v values being updated by pj
    - i.e., for all v >= pi^2: floor(v/pi) < pj^2
    - i.e., v/pi < pj^2, so v < pi * pj^2

    This means primes pi and pj can be processed in parallel for ALL v
    only if max(v) < pi * pj^2 AND max(v) < pj * pi^2.
    Since max(v) = x, this requires x < min(pi * pj^2, pj * pi^2).
    For pi < pj: x < pi * pj^2 is the binding constraint.

    For x = 10^6, primes up to 1000:
    Can p=997 and p=991 be parallel? Need 10^6 < 991 * 997^2 ≈ 10^9. YES!
    Can p=2 and p=3 be parallel? Need 10^6 < 2 * 9 = 18. NO!
    """
    V = get_floor_values(x)
    sqrtx = int(math.isqrt(x))
    primes = [p for p in range(2, sqrtx + 1) if sympy.isprime(p)]

    print(f"\n=== Parallel Round Structure for x = {x} ===")

    # Greedy: assign primes to rounds. Two primes p < q can be in same round
    # if for all v: processing p doesn't affect q's inputs and vice versa.
    # Sufficient condition: x < p * q^2 AND x < q * p^2
    # ⟺ x < min(p * q^2, q * p^2) = p * q * min(q, p) = p * q * p (since p < q)
    # ⟺ x < p^2 * q

    # Greedy algorithm: process primes in order, group into rounds
    rounds = []
    current_round = []

    for p in primes:
        can_add = True
        for q in current_round:
            # Check if p and q can be processed simultaneously
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

    print(f"Number of parallel rounds (greedy): {len(rounds)}")
    print(f"Sequential steps (naive): {len(primes)}")
    print(f"Speedup: {len(primes) / len(rounds):.2f}x")
    print(f"log(log(x)): {math.log(math.log(x)):.2f}" if x > 2 else "")

    for i, r in enumerate(rounds[:10]):
        print(f"  Round {i}: primes {r}")
    if len(rounds) > 10:
        print(f"  ... ({len(rounds) - 10} more rounds)")
        print(f"  Last round: {rounds[-1]}")

    return rounds

def main():
    print("=" * 70)
    print("LUCY DP FLOOR-VALUE MAPPING STRUCTURE ANALYSIS")
    print("Session 12 - Investigating circuit complexity of pi(x)")
    print("=" * 70)

    # Small case analysis
    for x in [100, 1000, 10000]:
        analyze_mapping_structure(x)
        print()

    # Depth analysis
    for x in [100, 1000, 10000]:
        analyze_depth_structure(x)

    # Recursion depth scaling
    analyze_recursion_depth()

    # Floor value set growth
    analyze_floor_value_growth()

    # Parallel structure
    for x in [100, 1000, 10000, 100000, 1000000]:
        analyze_parallel_structure(x)

    print("\n" + "=" * 70)
    print("CONCLUSIONS")
    print("=" * 70)
    print("""
Key findings:
1. The Lucy DP mapping matrices do NOT commute (floor-value conditions differ)
2. The sequential depth is pi(sqrt(x)), but the EFFECTIVE parallel depth
   may be much smaller if primes can be grouped into rounds
3. The Meissel-Lehmer recursion has depth O(log log x) = O(log N)
4. |V| = O(sqrt(x)) = O(2^{N/2}), exponential in input size
5. The parallel round structure determines the actual circuit depth needed

If the number of rounds is O(log log x) or O(1), this suggests
TC^0 or NC^1 feasibility. If it's Omega(sqrt(x)/log(x)), then
the DP is inherently sequential.
""")

if __name__ == "__main__":
    main()

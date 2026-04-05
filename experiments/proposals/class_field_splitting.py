#!/usr/bin/env python3
"""
Class Field Tower Prime Reconstruction -- Experimental Test

IDEA: A prime p splits in Q(sqrt(d)) iff Legendre(d, p) = 1.
The splitting pattern across m quadratic fields determines p mod M
for some large M. If we can compute pi_split(x, K) faster than pi(x),
we might reconstruct p(n) in polylog time.

TESTS:
1. Enumerate fundamental discriminants d in [-100, 100]
2. For primes p(n) with n <= 5000, compute splitting patterns
3. Measure information content: how many bits per field?
4. Given m fields, how many prime candidates survive?
5. Key question: is pi_split(x, K) any easier than pi(x)?

References:
- Chebotarev density theorem: density of primes splitting in K = 1/[K:Q]
- For quadratic fields, this is exactly 1/2 (up to finitely many exceptions)
- Closed path: "Chebotarev density theorem" (FAIL, C, session 4)
- Closed path: "Number field algebraic decomposition" (FAIL, E, session 7)
"""

import time
import math
from collections import defaultdict
from sympy import isprime, nextprime, primerange, legendre_symbol, jacobi_symbol
from sympy import factorint

# ============================================================
# Part 0: Utility functions
# ============================================================

def fundamental_discriminants(bound):
    """Return fundamental discriminants d with |d| <= bound."""
    fund_discs = []
    for d in range(-bound, bound + 1):
        if d == 0 or d == 1:
            continue
        # d is a fundamental discriminant if:
        # d = 1 mod 4 and d is squarefree, OR
        # d = 0 mod 4, d/4 is squarefree, and d/4 = 2 or 3 mod 4
        if d % 4 == 1 or (d < 0 and (-d) % 4 == 3):
            # Check if d is squarefree
            abs_d = abs(d)
            if abs_d == 1:
                continue
            factors = factorint(abs_d)
            if all(e == 1 for e in factors.values()):
                fund_discs.append(d)
        elif d % 4 == 0:
            q = d // 4
            abs_q = abs(q)
            if abs_q == 0:
                continue
            factors = factorint(abs_q)
            if all(e == 1 for e in factors.values()):
                if q % 4 == 2 or q % 4 == 3 or q % 4 == -1 or q % 4 == -2:
                    fund_discs.append(d)
    return sorted(fund_discs)


def splitting_pattern(p, discriminants):
    """
    For a prime p and list of fundamental discriminants,
    return the splitting pattern: +1 (splits), -1 (inert), 0 (ramifies).
    Uses the Kronecker symbol (d/p).
    """
    pattern = []
    for d in discriminants:
        if p == 2:
            # Special handling for p=2
            if d % 2 == 0:
                pattern.append(0)  # ramifies
            elif d % 8 in (1, 7):
                pattern.append(1)  # splits
            else:
                pattern.append(-1)  # inert
        else:
            if d % p == 0:
                pattern.append(0)  # ramifies
            else:
                pattern.append(legendre_symbol(d, p))
    return tuple(pattern)


# ============================================================
# Part 1: Enumerate fundamental discriminants
# ============================================================

print("=" * 70)
print("PART 1: Fundamental discriminants in [-100, 100]")
print("=" * 70)

fund_discs = fundamental_discriminants(100)
print(f"Number of fundamental discriminants: {len(fund_discs)}")
print(f"First 20: {fund_discs[:20]}")
print(f"Last 20: {fund_discs[-20:]}")

# For the experiment, use positive fundamental discriminants (real quadratic fields)
# and a few negative ones
pos_discs = [d for d in fund_discs if d > 0]
neg_discs = [d for d in fund_discs if d < 0]
print(f"\nPositive fundamental discriminants: {len(pos_discs)}")
print(f"Negative fundamental discriminants: {len(neg_discs)}")

# ============================================================
# Part 2: Compute splitting patterns for p(n), n <= 5000
# ============================================================

print("\n" + "=" * 70)
print("PART 2: Splitting patterns for first 5000 primes")
print("=" * 70)

primes_list = list(primerange(2, 50000))[:5000]
print(f"Number of primes: {len(primes_list)}")
print(f"p(1)={primes_list[0]}, p(5000)={primes_list[-1]}")

# Use a subset of discriminants for tractability
# Pick small discriminants that give independent information
test_discs = [d for d in fund_discs if abs(d) <= 60]
print(f"\nUsing {len(test_discs)} discriminants with |d| <= 60")

t0 = time.time()
patterns = {}
for i, p in enumerate(primes_list):
    patterns[p] = splitting_pattern(p, test_discs)
t1 = time.time()
print(f"Computed all patterns in {t1 - t0:.3f}s")

# ============================================================
# Part 3: Information content -- bits per field
# ============================================================

print("\n" + "=" * 70)
print("PART 3: Information content per quadratic field")
print("=" * 70)

# For each discriminant d, count how many primes split vs are inert
# Chebotarev says it should be ~50/50
for d in test_discs[:10]:
    splits = sum(1 for p in primes_list if patterns[p][test_discs.index(d)] == 1)
    inerts = sum(1 for p in primes_list if patterns[p][test_discs.index(d)] == -1)
    ramifs = sum(1 for p in primes_list if patterns[p][test_discs.index(d)] == 0)
    frac = splits / (splits + inerts) if (splits + inerts) > 0 else 0
    print(f"  d={d:4d}: splits={splits}, inert={inerts}, ramifies={ramifs}, "
          f"split_frac={frac:.4f}")

print("\n  (Chebotarev predicts split_frac -> 0.5000)")

# Measure: using first m discriminants, how many distinct patterns are there?
print("\nDistinct patterns vs number of fields used:")
for m in [1, 2, 5, 10, 15, 20, 30, 40, 50, len(test_discs)]:
    if m > len(test_discs):
        break
    sub_patterns = set()
    for p in primes_list:
        sub_patterns.add(patterns[p][:m])
    bits = math.log2(len(sub_patterns)) if len(sub_patterns) > 1 else 0
    print(f"  m={m:3d} fields: {len(sub_patterns):6d} distinct patterns "
          f"({bits:.1f} bits, theoretical max = {m:.1f} bits)")

# ============================================================
# Part 4: Prime candidate filtering
# ============================================================

print("\n" + "=" * 70)
print("PART 4: How many prime candidates survive the splitting filter?")
print("=" * 70)

# For a specific target prime, how many other primes share its pattern?
# This tells us how much the splitting data narrows down the answer.
test_indices = [100, 500, 1000, 2500, 5000]
for idx in test_indices:
    if idx > len(primes_list):
        break
    target = primes_list[idx - 1]
    target_pattern = patterns[target]

    # Count primes with the same pattern
    matches = [p for p in primes_list if patterns[p] == target_pattern]
    print(f"\n  p({idx}) = {target}")
    print(f"    Primes sharing full pattern ({len(test_discs)} fields): "
          f"{len(matches)} out of {len(primes_list)}")
    if len(matches) <= 20:
        print(f"    Matching primes: {matches}")

    # Also show for fewer fields
    for m in [5, 10, 20, 30]:
        if m > len(test_discs):
            break
        sub_target = target_pattern[:m]
        sub_matches = sum(1 for p in primes_list if patterns[p][:m] == sub_target)
        print(f"    With {m:2d} fields: {sub_matches} matches")

# ============================================================
# Part 5: CRT reconstruction -- what modulus can we determine?
# ============================================================

print("\n" + "=" * 70)
print("PART 5: CRT reconstruction from splitting data")
print("=" * 70)

# The splitting pattern of p in Q(sqrt(d)) determines (d/p),
# which is a function of p mod |d| (for odd primes not dividing d).
# So splitting in m fields determines p mod lcm(d_1, ..., d_m).

# But we only get 1 bit per field (split or inert), NOT the full residue.
# Quadratic reciprocity: (d/p) depends on p mod 4|d|.
# Actually, for fundamental discriminant d, (d/p) depends on p mod |d|.

# What modulus M can the first m discriminants determine?
small_pos = [abs(d) for d in test_discs if d > 0][:30]
small_neg = [abs(d) for d in test_discs if d < 0][:30]
all_abs = [abs(d) for d in test_discs]

print("\nModulus achievable from discriminants:")
running_lcm = 1
for i, d in enumerate(all_abs[:40]):
    running_lcm = math.lcm(running_lcm, d)
    if (i + 1) in [1, 2, 3, 5, 10, 15, 20, 25, 30, 35, 40]:
        log_M = math.log2(running_lcm)
        print(f"  First {i+1:2d} discs: lcm = {running_lcm} "
              f"({log_M:.1f} bits)")

print(f"\n  Full lcm of all {len(all_abs)} discriminants: "
      f"{math.log2(math.lcm(*all_abs)):.1f} bits")
print(f"  For p(10^100) ~ 10^102, we need ~339 bits")

# ============================================================
# Part 6: THE KEY QUESTION -- Is pi_split easier than pi?
# ============================================================

print("\n" + "=" * 70)
print("PART 6: Is pi_split(x, K) easier to compute than pi(x)?")
print("=" * 70)

print("""
ANALYSIS: For K = Q(sqrt(d)), we have:
  pi_split(x, d) = #{p <= x : (d/p) = 1}

By orthogonality of characters:
  pi_split(x, d) = (1/2) * [pi(x) + sum_{chi mod |d|, chi != chi_0} chi(d) * pi(x, chi)]

where pi(x, chi) = sum_{p<=x} chi(p).

KEY INSIGHT: pi_split(x, d) is a LINEAR COMBINATION of pi(x) and
Dirichlet character sums. Computing pi_split is AT LEAST as hard as
computing pi(x), because:

  pi(x) = 2 * pi_split(x, d) - (correction from pi(x, chi))

The character sums pi(x, chi) are computed via L(s, chi) which has
its OWN zeros -- the generalized Riemann zeros. These are just as
hard to handle as the zeros of zeta(s).
""")

# Empirical verification: pi_split(x, d) vs pi(x)/2
print("Empirical pi_split vs pi(x)/2:")
for d in [5, -3, 8, -7, 12, 13, -4]:
    if d not in test_discs:
        continue
    d_idx = test_discs.index(d)
    N = len(primes_list)
    pi_split = sum(1 for p in primes_list if patterns[p][d_idx] == 1)
    deviation = pi_split - N / 2
    rel_dev = deviation / math.sqrt(N)
    print(f"  d={d:4d}: pi_split={pi_split}, pi(x)/2={N/2:.0f}, "
          f"deviation={deviation:+.0f}, "
          f"deviation/sqrt(N)={rel_dev:+.2f}")

print("""
The deviations are O(sqrt(N)) -- consistent with the prime number
theorem for arithmetic progressions. This deviation is EXACTLY the
information encoded in the generalized Riemann zeros of L(s, chi_d).
""")

# ============================================================
# Part 7: Correlation between splitting patterns
# ============================================================

print("=" * 70)
print("PART 7: Independence of splitting patterns across fields")
print("=" * 70)

# Are different fields giving independent bits?
# Check pairwise correlation
import itertools

print("Pairwise correlations between first 10 discriminants:")
disc_subset = test_discs[:10]
for i, j in itertools.combinations(range(len(disc_subset)), 2):
    d1, d2 = disc_subset[i], disc_subset[j]
    corr_sum = 0
    count = 0
    for p in primes_list:
        s1 = patterns[p][i]
        s2 = patterns[p][j]
        if s1 != 0 and s2 != 0:  # exclude ramified
            corr_sum += s1 * s2
            count += 1
    corr = corr_sum / count if count > 0 else 0
    if abs(corr) > 0.05:
        print(f"  ({d1:4d}, {d2:4d}): correlation = {corr:+.4f}  <-- NOTABLE")

print("\n(Only showing correlations > 0.05 in absolute value)")
print("Most pairs should have near-zero correlation (independence).")

# ============================================================
# Part 8: Can we avoid computing individual primes?
# ============================================================

print("\n" + "=" * 70)
print("PART 8: Circularity analysis")
print("=" * 70)

print("""
To reconstruct p(n) from splitting data, we need:
  1. The COUNTS pi_split(x, K_i) for each field K_i
  2. These counts at the right value of x (which IS p(n))

CIRCULARITY CHECK:
  - To compute pi_split(x, d) = #{p <= x : (d/p) = 1}, we need to
    enumerate primes up to x. This is exactly as hard as computing pi(x).

  - Alternative: use the explicit formula for pi(x, chi):
    pi(x, chi) = li(x) - sum_{rho} li(x^rho) + ...
    where rho ranges over zeros of L(s, chi_d).

  - This requires the zeros of L(s, chi_d), which are ADDITIONAL
    quantities beyond the zeros of zeta(s). More zeros = more work.

  - Even if we KNEW all the splitting counts, CRT reconstruction
    gives p(n) mod M. To get p(n) exactly, we need M > p(n),
    which means lcm of discriminants must be > p(n) ~ n ln n.
    The number of fundamental discriminants d with |d| <= D is O(D).
    We need lcm(d_1,...,d_m) > n ln n.
    By prime number theorem for discriminants, lcm grows as e^D,
    so we need D ~ log(n ln n) ~ log n discriminants.
    That's O(log n) fields -- fine for polylog.

  - BUT each field requires computing its OWN set of L-function zeros.
    Total zero computation: O(log n) fields * O(x^{1/2+eps}) per field
    = O(x^{1/2+eps} * log n). WORSE than just using zeta alone.

VERDICT: The splitting approach INCREASES total computation because
we need zeros of MULTIPLE L-functions instead of just zeta(s).
The information is there, but accessing it is HARDER, not easier.
""")

# ============================================================
# Part 9: Summary
# ============================================================

print("=" * 70)
print("SUMMARY")
print("=" * 70)

print(f"""
EXPERIMENT: Class Field Tower Prime Reconstruction
DISCRIMINANTS TESTED: {len(test_discs)} fundamental discriminants, |d| <= 60
PRIMES TESTED: {len(primes_list)} (up to p(5000) = {primes_list[-1]})

FINDINGS:
1. INFORMATION: Each quadratic field gives ~1 bit about p.
   With {len(test_discs)} fields, we get {math.log2(len(set(patterns[p] for p in primes_list))):.1f}
   effective bits -- enough to distinguish all 5000 primes.

2. CRT MODULUS: lcm of {len(all_abs)} discriminants gives
   {math.log2(math.lcm(*all_abs)):.1f} bits -- sufficient for small primes,
   but we need ~339 bits for p(10^100).

3. INDEPENDENCE: Splitting patterns across different fields are
   approximately independent (low pairwise correlation), confirming
   ~1 bit per field.

4. CIRCULARITY (FATAL): Computing pi_split(x, d) requires either:
   (a) Enumerating primes up to x (circular, same cost as pi(x))
   (b) Explicit formula with zeros of L(s, chi_d) -- HARDER than
       using zeta(s) alone, since we need zeros of MULTIPLE L-functions.

5. NET EFFECT: The splitting approach INCREASES computation cost.
   Instead of O(x^{{1/2+eps}}) for one L-function, we need
   O(x^{{1/2+eps}} * log n) for O(log n) L-functions.

CLASSIFICATION: FAIL (C + E)
  C = Circularity: need primes to compute splitting counts
  E = Equivalence: explicit formula route is strictly harder

This confirms the closed-path entries:
  - "Chebotarev density theorem" (session 4)
  - "Number field algebraic decomposition" (session 7)
""")

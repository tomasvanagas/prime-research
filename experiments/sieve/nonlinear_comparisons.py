"""
Experiment: Comparisons and thresholding of floor values for prime detection.

Key idea: Nonlinear operations like mod 2, sign, thresholding on floor(x/k)
might detect primes more directly than linear combinations.

The gap g(k,x) = floor(x/k) - floor(x/(k+1)) counts how many multiples of k
are in (x/(k+1), x/k]. For k=1 this is roughly x/2. For prime p <= x,
floor(x/p) - floor(x/(p+1)) >= 1 always.

Key nonlinear operations to test:
1. floor(x/k) mod m for various m
2. sign(floor(x/k) - floor(x/(k+1)) - threshold)
3. Parity of floor(x/k)
4. Indicator: floor(x/k) != floor(x/(k+1)) (nonzero gap)
"""

import numpy as np
from sympy import primepi, isprime, primerange
import math

# ============================================================
# Experiment 1: Parity of floor values
# ============================================================
print("=" * 70)
print("EXPERIMENT 1: Parity patterns in floor(x/k)")
print("=" * 70)

# For fixed x, look at the binary string b_k = floor(x/k) mod 2 for k=1..x
# Does this binary string encode pi(x)?

for x in [30, 50, 100, 200]:
    parity_str = ""
    for k in range(1, x + 1):
        parity_str += str((x // k) % 2)

    # Count 1s
    ones = parity_str.count('1')
    pi_x = primepi(x)

    # Check parity at prime positions
    prime_parities = []
    composite_parities = []
    for k in range(2, x + 1):
        if isprime(k):
            prime_parities.append((x // k) % 2)
        else:
            composite_parities.append((x // k) % 2)

    print(f"\nx={x}: pi(x)={pi_x}")
    print(f"  Total 1s in parity string: {ones}/{x}")
    print(f"  Parity at prime k: mean={np.mean(prime_parities):.3f}")
    print(f"  Parity at composite k: mean={np.mean(composite_parities):.3f}")
    print(f"  Difference: {abs(np.mean(prime_parities) - np.mean(composite_parities)):.3f}")

# ============================================================
# Experiment 2: Gap function g(k,x) = floor(x/k) - floor(x/(k+1))
# ============================================================
print("\n" + "=" * 70)
print("EXPERIMENT 2: Gap function and prime detection")
print("=" * 70)

# g(k,x) = floor(x/k) - floor(x/(k+1))
# For k=p prime: floor(x/p) counts multiples of p up to x
# g(p,x) = #{multiples of p in (x/(p+1), x/p]} -- at least 1 if p <= x

for x in [100, 500, 1000]:
    gaps = {}
    for k in range(1, x):
        gaps[k] = (x // k) - (x // (k + 1))

    # How often is g(k,x) = 0?
    zero_gaps = [k for k, g in gaps.items() if g == 0]
    nonzero_prime = sum(1 for k in range(2, x) if gaps[k] > 0 and isprime(k))
    nonzero_composite = sum(1 for k in range(2, x) if gaps[k] > 0 and not isprime(k))
    zero_prime = sum(1 for k in range(2, x) if gaps[k] == 0 and isprime(k))
    zero_composite = sum(1 for k in range(2, x) if gaps[k] == 0 and not isprime(k))

    pi_x = primepi(x)
    print(f"\nx={x}: pi(x)={pi_x}")
    print(f"  Nonzero gaps: {len(gaps) - len(zero_gaps)} (prime: {nonzero_prime}, composite: {nonzero_composite})")
    print(f"  Zero gaps: {len(zero_gaps)} (prime: {zero_prime}, composite: {zero_composite})")
    print(f"  Distinct floor(x/k) values: {len(set(x // k for k in range(1, x + 1)))}")

    # The number of distinct floor(x/k) values = number of nonzero gaps + 1
    # This is O(sqrt(x))

    # g(k,x) = 1 filter: how many primes have gap exactly 1?
    gap1_prime = sum(1 for k in range(2, x) if gaps[k] == 1 and isprime(k))
    gap1_composite = sum(1 for k in range(2, x) if gaps[k] == 1 and not isprime(k))
    print(f"  Gap=1: {gap1_prime} primes, {gap1_composite} composites")

    # For large k (k > sqrt(x)), g(k,x) is 0 or 1
    sqrt_x = int(math.sqrt(x))
    large_k_prime = sum(1 for k in range(sqrt_x, x) if gaps.get(k, 0) == 1 and isprime(k))
    large_k_composite = sum(1 for k in range(sqrt_x, x) if gaps.get(k, 0) == 1 and not isprime(k))
    print(f"  Gap=1 for k>sqrt(x)={sqrt_x}: {large_k_prime} primes, {large_k_composite} composites")

# ============================================================
# Experiment 3: Nonlinear combination -- products of gap indicators
# ============================================================
print("\n" + "=" * 70)
print("EXPERIMENT 3: Products of gap indicators for prime counting")
print("=" * 70)

# Idea: For k > sqrt(x), floor(x/k) is distinct iff gap(k) >= 1
# The indicator I(k) = min(1, gap(k)) = 1 iff floor(x/k) != floor(x/(k+1))
# For k in (sqrt(x), x]: I(k) = 1 iff there's at least one multiple of some
# number in the range. But this doesn't directly count primes.

# Better idea: Mobius-like inclusion-exclusion with nonlinear terms
# pi(x) = sum_{k=2}^{x} [k is prime]
# [k is prime] = [k not divisible by any p < k]
# = product_{p < k, p prime} (1 - [p | k])
# But this is circular. However:
# [p | k] = [floor(k/p)*p == k] = [k - floor(k/p)*p == 0]
# which involves a product (floor * p) and a comparison.

# A TC^0-friendly version:
# [p | k] can be computed as: (k mod p == 0) which is a comparison after modular reduction
# The issue is we need to iterate over primes p, which is circular.

# Instead, Eratosthenes without knowing primes:
# composite(k) = OR_{j=2}^{sqrt(k)} [j | k]
# prime(k) = NOT composite(k) AND k >= 2
# This uses O(sqrt(x)) divisibility checks per number, O(x * sqrt(x)) total

# But with floor values of x:
# For each k, [k | x] = (x mod k == 0) = (x - k * floor(x/k) == 0)
# This counts divisors of x, not primes up to x.

# Key question: can we use floor(x/k) values to count primes WITHOUT
# testing each number individually?

# Let's try: sum of nonlinear functions of floor values
x = 100
pi_x = primepi(x)
print(f"\nx={x}, pi(x)={pi_x}")

# Attempt: sum_{k=2}^{x} f(floor(x/k), floor(x/(k-1)), k) = pi(x)?
# where f is some nonlinear function

# Observation: floor(x/k) = floor(x/(k+1)) + gap(k)
# For k > sqrt(x): gap(k) in {0, 1}
# gap(k) = 1 iff there exists m such that x/(k+1) < m <= x/k
# i.e., k <= x/m < k+1, i.e., m = floor(x/k) and k = floor(x/m)

# The set of k where gap(k)=1 for k > sqrt(x) has size O(sqrt(x))
# These are exactly the k = floor(x/m) for m = 1, 2, ..., sqrt(x)

# Can we extract pi(x) from knowing WHICH k values have gap=1?
# gap(k)=1 for k>sqrt(x) iff k = floor(x/m) for some m <= sqrt(x)
# These are the "jump points" of the Dirichlet hyperbola

# How many of these jump points are prime?
sqrt_x = int(math.sqrt(x))
jump_primes = 0
jump_total = 0
for m in range(1, sqrt_x + 1):
    k = x // m
    if k > sqrt_x:
        jump_total += 1
        if isprime(k):
            jump_primes += 1

# Primes > sqrt(x) and <= x that are NOT jump points
primes_above_sqrt = list(primerange(sqrt_x + 1, x + 1))
non_jump_primes = len(primes_above_sqrt) - jump_primes
print(f"  sqrt(x)={sqrt_x}")
print(f"  Jump points above sqrt(x): {jump_total}")
print(f"  Jump points that are prime: {jump_primes}")
print(f"  Primes above sqrt(x): {len(primes_above_sqrt)}")
print(f"  Primes NOT at jump points: {non_jump_primes}")
print(f"  (These primes are 'invisible' to the gap function)")

# ============================================================
# Experiment 4: Modular floor values -- floor(x/k) mod m
# ============================================================
print("\n" + "=" * 70)
print("EXPERIMENT 4: floor(x/k) mod m patterns")
print("=" * 70)

for x in [100, 500]:
    pi_x = primepi(x)
    for m in [2, 3, 6]:
        # For each k, compute floor(x/k) mod m
        mod_vals = {k: (x // k) % m for k in range(1, x + 1)}

        # Sum of floor(x/k) mod m
        total = sum(mod_vals.values())

        # Sum only at prime k
        prime_sum = sum(mod_vals[k] for k in range(2, x + 1) if isprime(k))

        # Sum only at composite k
        comp_sum = sum(mod_vals[k] for k in range(2, x + 1) if not isprime(k) and k > 1)

        print(f"\nx={x}, mod {m}: total={total}, prime_sum={prime_sum}, comp_sum={comp_sum}")

        # Correlation with primality
        prime_indicator = np.array([1 if isprime(k) else 0 for k in range(2, x + 1)])
        mod_vector = np.array([(x // k) % m for k in range(2, x + 1)])
        corr = np.corrcoef(prime_indicator, mod_vector)[0, 1]
        print(f"  Correlation of (floor(x/k) mod {m}) with primality: {corr:.4f}")

# ============================================================
# Experiment 5: XOR / AND of floor values
# ============================================================
print("\n" + "=" * 70)
print("EXPERIMENT 5: Bitwise operations on floor values")
print("=" * 70)

for x in [50, 100, 200]:
    pi_x = primepi(x)

    # XOR of all floor(x/k) for k = 1..x
    xor_all = 0
    for k in range(1, x + 1):
        xor_all ^= (x // k)

    # XOR of floor(x/p) for primes p
    xor_primes = 0
    for p in primerange(2, x + 1):
        xor_primes ^= (x // p)

    # AND of all floor(x/k)
    and_all = (x // 1)
    for k in range(2, x + 1):
        and_all &= (x // k)

    # Popcount of floor(x/k) summed
    popcount_sum = sum(bin(x // k).count('1') for k in range(1, x + 1))

    print(f"\nx={x}: pi(x)={pi_x}")
    print(f"  XOR of all floor(x/k): {xor_all}")
    print(f"  XOR of floor(x/p) for primes: {xor_primes}")
    print(f"  AND of all floor(x/k): {and_all}")
    print(f"  Sum of popcount(floor(x/k)): {popcount_sum}")

    # Can popcount_sum predict pi(x)?
    # Ratio
    print(f"  popcount_sum / x = {popcount_sum / x:.4f}")
    print(f"  pi(x) / x = {pi_x / x:.4f}")

# ============================================================
# Experiment 6: Threshold detection
# ============================================================
print("\n" + "=" * 70)
print("EXPERIMENT 6: Threshold-based prime indicators")
print("=" * 70)

# For k in range, define T(k,x) = sign(floor(x/k) - floor(x/(k+1)) - 1)
# T(k,x) = 1 if gap > 1, 0 if gap = 1, -1 if gap = 0
# For small k: gap is large (roughly x/k^2), so T=1
# For k ~ sqrt(x): gap transitions from >1 to 0 or 1

for x in [100, 500, 1000]:
    pi_x = primepi(x)
    sqrt_x = int(math.sqrt(x))

    # Count primes using threshold on gaps
    # For k > sqrt(x): gap in {0,1}
    # sum of gaps for k > sqrt(x) = floor(x/sqrt(x)) - 1 (telescoping)

    gap_sum_large_k = (x // (sqrt_x + 1))  # floor(x/(sqrt(x)+1)) - floor(x/x) but not quite

    # Actually: sum_{k=sqrt(x)+1}^{x} gap(k) = floor(x/(sqrt(x)+1)) - 0 = floor(x/(sqrt(x)+1))
    # because floor(x/x)=1 and the sum telescopes
    actual_gap_sum = sum((x // k) - (x // (k + 1)) for k in range(sqrt_x + 1, x))

    # This sum counts the number of distinct floor values in the range (sqrt(x), x]
    # which is floor(x/(sqrt(x)+1))

    print(f"\nx={x}, sqrt(x)={sqrt_x}, pi(x)={pi_x}")
    print(f"  Sum of gap(k) for k>sqrt(x): {actual_gap_sum}")
    print(f"  floor(x/(sqrt(x)+1)) = {x // (sqrt_x + 1)}")
    print(f"  Number of k>sqrt(x) with gap=1: {sum(1 for k in range(sqrt_x+1, x) if (x//k)-(x//(k+1))==1)}")

    # The gap=1 positions above sqrt(x) are floor(x/m) for m=1..sqrt(x)
    # So there are sqrt(x) of them.
    # Of these, how many are prime? That's pi(x) - pi(sqrt(x)) approximately
    # No wait, not all primes > sqrt(x) appear as floor(x/m).

    # Let's check: which primes > sqrt(x) are of the form floor(x/m)?
    floor_set = set(x // m for m in range(1, sqrt_x + 1))
    primes_in_set = sum(1 for v in floor_set if v > sqrt_x and isprime(v))
    primes_above = len(list(primerange(sqrt_x + 1, x + 1)))

    print(f"  floor(x/m) values above sqrt(x): {len([v for v in floor_set if v > sqrt_x])}")
    print(f"  Of these, primes: {primes_in_set}")
    print(f"  Total primes in (sqrt(x), x]: {primes_above}")
    print(f"  Fraction of primes captured: {primes_in_set/primes_above:.3f}")

print("\nDone.")

#!/usr/bin/env python3
"""
Session 10: Parallel Complexity of π(x) / p(n)

KEY INSIGHT: π(x) requires O(x^{2/3}) WORK but how much DEPTH?
If the Meissel-Lehmer formula can be parallelized to polylog depth,
then π(x) ∈ NC and p(n) can be computed in polylog time with
polynomial processors!

The Meissel-Lehmer formula has structure:
  π(x) = φ(x, a) + a - 1 - P₂(x, a)
where a = π(x^{1/3})

The key operations are:
1. φ(x, a) = sum over smooth numbers (O(x^{2/3}) terms, parallelizable?)
2. P₂(x, a) = sum over pairs of primes
3. Recursive calls: π(x/p) for various p

QUESTION: What's the circuit depth of computing π(x)?
"""

import math
from sympy import prime, primepi, isprime
import time

# ============================================================
# PART 1: Meissel-Lehmer structure analysis
# ============================================================
print("=" * 60)
print("PART 1: Meissel-Lehmer Depth Analysis")
print("=" * 60)

# The recursive structure of Meissel-Lehmer:
# π(x) depends on π(x^{1/3}), which depends on π(x^{1/9}), etc.
# Recursion depth = log₃(log x)

for log_x in [10, 20, 50, 100, 200, 1000]:
    x = 10**log_x
    # Recursion depth: how many times do we take cube root before x < 10?
    depth = 0
    val = log_x  # We track log₁₀(x) for convenience
    while val > 1:  # Stop when x < 10 (can enumerate directly)
        val = val / 3
        depth += 1

    # At each level, the number of subproblems
    # Level 0: 1 problem of size x
    # Level 1: ~x^{2/3}/log(x) subproblems of size ~x^{1/3}
    # Level 2: ~x^{2/9}/log(x^{1/3}) subproblems of size ~x^{1/9}
    # etc.

    level_sizes = []
    y = log_x
    for d in range(depth):
        # Number of subproblems at this level ≈ x^{2/3 * (1/3)^d} / log terms
        exponent = (2/3) * (1/3)**d
        num_probs = 10**(log_x * exponent) if log_x * exponent < 20 else float('inf')
        level_sizes.append((d, y, exponent, num_probs))
        y = y / 3

    print(f"\nx = 10^{log_x}:")
    print(f"  Recursion depth: {depth}")
    print(f"  Level structure:")
    for d, size_log, exp, num in level_sizes[:6]:
        if num < 1e20:
            print(f"    Level {d}: ~10^{size_log:.1f} per subproblem, "
                  f"~10^{log_x * exp:.1f} subproblems")
        else:
            print(f"    Level {d}: ~10^{size_log:.1f} per subproblem, "
                  f"~10^{log_x * exp:.0f} subproblems")

    # Total WORK = sum of subproblems × cost per subproblem
    # Total DEPTH = recursion depth × cost per level
    # If each level can be parallelized (sum reduction):
    # Depth per level = O(log(num_subproblems)) = O(log x)
    # Total depth = recursion_depth × O(log x) = O(log x × log log x)

    parallel_depth = depth * log_x * math.log(10) / math.log(2)  # in bits
    print(f"  Estimated parallel depth: O({depth} × {log_x:.0f} × log2(10)) ≈ {parallel_depth:.0f} binary ops")
    print(f"  This is O(log²(x) × log(log(x))) — POLYLOGARITHMIC!")

# ============================================================
# PART 2: Can the Meissel-Lehmer sum be parallelized?
# ============================================================
print("\n" + "=" * 60)
print("PART 2: Parallelizability of Meissel-Lehmer")
print("=" * 60)

# The main sum in Meissel-Lehmer is:
# φ(x, a) = #{n ≤ x : gcd(n, P(a)) = 1} where P(a) = p₁p₂...p_a
#
# This can be computed via inclusion-exclusion:
# φ(x, a) = x - Σᵢ ⌊x/pᵢ⌋ + Σᵢ<ⱼ ⌊x/(pᵢpⱼ)⌋ - ...
#
# The number of terms in the full inclusion-exclusion is 2^a (exponential in a).
# But with the recursive formula:
# φ(x, a) = φ(x, a-1) - φ(x/p_a, a-1)
#
# This recursion has depth a = π(x^{1/3}).
# For x = 10^100: a = π(10^{33}) ≈ 10^33/ln(10^33) ≈ 1.3 × 10^31
#
# So the recursion depth of φ is ~10^31, which is NOT polylog!
#
# BUT: the Deleglise-Rivat optimization truncates the recursion and
# uses direct formulas for the "ordinary" and "special" leaves.
# The truncation reduces depth to O(x^{1/3}/log x) ≈ 10^31 still...
#
# However, the PARALLEL depth of computing a sum of N independent terms
# is O(log N). The question is whether the terms are independent.

print("φ(x, a) recursion analysis:")
print("  φ(x, a) = φ(x, a-1) - φ(x/p_a, a-1)")
print("  This recursion has depth a = π(x^{1/3})")
print()

for log_x in [10, 20, 50, 100]:
    x_cube_root = 10**(log_x/3)
    a_approx = x_cube_root / (log_x/3 * math.log(10))
    print(f"  x = 10^{log_x}: a ≈ {a_approx:.2e}, depth of φ recursion = {a_approx:.2e}")

print()
print("CRITICAL: The φ recursion has depth O(x^{1/3}) — NOT polylog!")
print("This is the SEQUENTIAL bottleneck of Meissel-Lehmer.")
print()

# BUT: can we avoid the φ recursion entirely?
# The Lucy_Hedgehog method computes π(x) using dynamic programming
# on the set {⌊x/k⌋ : k = 1,...,x}. This set has O(√x) elements.
# The DP fills the table in O(x^{2/3}/log x) time.
#
# Is this DP parallelizable?
# The recurrence is: S(v, p) = S(v, p-1) - (S(v/p, p-1) - S(p-1, p-1))
# Each level p depends on the PREVIOUS level — sequential in p!
# The number of levels = π(√x), which is ~√x / ln(√x).
# For x = 10^100: ~10^50 / 115 ≈ 10^48 sequential steps!
#
# Even within each level, the ~√x updates can be done in parallel.
# But the DEPTH is O(π(√x)) = O(√x / log x) — NOT polylog.

print("Lucy_Hedgehog DP analysis:")
print("  Recurrence: S(v,p) = S(v,p-1) - (S(v/p,p-1) - S(p-1,p-1))")
print("  Sequential depth = π(√x)")
for log_x in [10, 20, 50, 100]:
    depth = 10**(log_x/2) / (log_x/2 * math.log(10))
    print(f"  x = 10^{log_x}: depth = {depth:.2e}")

print()
print("BOTH Meissel-Lehmer AND Lucy_Hedgehog have non-polylog depth!")
print("The sequential bottleneck is inherent in known sieve-based methods.")

# ============================================================
# PART 3: Is there a PARALLEL sieve?
# ============================================================
print("\n" + "=" * 60)
print("PART 3: Parallel Sieve Complexity")
print("=" * 60)

# The sieve of Eratosthenes can be parallelized:
# For each prime p ≤ √x, mark multiples p², p²+p, ..., x
# These markings are INDEPENDENT across different primes!
# So all primes can sieve in PARALLEL.
#
# Depth: O(√x / log x) primes to sieve with, but they're independent
# → depth O(x/p) for the largest prime p ~ √x → O(√x) per prime
# But since primes sieve independently → depth O(√x)
#
# Actually, the parallel sieve works like this:
# 1. Find all primes ≤ √x (recursively, depth D(√x))
# 2. For each prime, mark composites (parallel, depth O(1) per mark with O(x) processors)
# 3. Count unmarked (parallel reduction, depth O(log x))
#
# Step 1 requires D(√x) depth recursively.
# D(x) = D(√x) + O(log x)
# D(x) = O(log²x) by master theorem!
#
# BUT: step 1 requires O(√x) processors at the top level,
# and step 2 requires O(x) processors... which is O(x) total.
# The CIRCUIT SIZE is O(x), which is NOT polynomial in log(x).

print("Parallel sieve of Eratosthenes:")
print("  Step 1: Find primes ≤ √x recursively")
print("  Step 2: For each prime, mark multiples (independent, parallel)")
print("  Step 3: Count unmarked (parallel reduction)")
print()
print("  Depth of step 1: D(x) = D(√x) + O(log x) = O(log²x)")
print("  BUT: circuit SIZE (processors) = O(x) — polynomial in x, not polylog!")
print()
print("  For p(n): x = p(n) ≈ n·ln(n). Circuit size = O(n·ln(n)).")
print("  This is polynomial in n — so p(n) is in NC with polynomial-in-n processors!")
print("  BUT: we want polylog(n) TIME on a SEQUENTIAL machine,")
print("  which would require polylog(n) depth AND polylog(n) processors.")
print()

# The key question: is π(x) in NC with polylog(log x) processors?
# This seems unlikely because the sieve requires O(x) work total.
# NC = polylog depth + polynomial processors.
# For π(x), polynomial means polynomial in x (the INPUT).
# But our input is just n (log(x) bits), so x = 2^{input_size}.
# Processors polynomial in INPUT = polynomial in log(x) = polylog(x).

print("COMPLEXITY CLASSIFICATION:")
print()
print("  Input: n (log₂(n) bits)")
print("  Output: p(n) (log₂(p(n)) ≈ log₂(n) + log₂(ln(n)) bits)")
print()
print("  In terms of INPUT SIZE s = log₂(n):")
print("  - Sieve: O(2^s) work — EXPONENTIAL in input size!")
print("  - Meissel-Lehmer: O(2^{2s/3}) work — still exponential")
print("  - Even Lucy_Hedgehog: O(2^{2s/3}) — exponential")
print()
print("  ALL known methods are EXPONENTIAL in the input size!")
print("  This means p(n) is NOT known to be in P (polynomial in input)!")
print()
print("  More precisely:")
print("    p(n) ∈ EXPTIME (computable in exponential time)")
print("    p(n) might or might not be in P")
print("    p(n) is NOT known to be in NC")
print()
print("  This is EXACTLY the gap the information theory agent identified!")
print("  K(p(n)|n) = O(s) bits, but computing those bits takes 2^{Θ(s)} time.")

# ============================================================
# PART 4: What would polylog p(n) mean?
# ============================================================
print("\n" + "=" * 60)
print("PART 4: Implications of polylog p(n)")
print("=" * 60)

print("""
If p(n) could be computed in polylog(n) time, it would mean:

1. π(x) ∈ P (polynomial in log x) — a MAJOR complexity result
2. The prime-counting function has a "shortcut" that avoids
   both the explicit formula AND the combinatorial sieve
3. This would likely imply new results about the Riemann zeta function

What's known:
  - PRIMES ∈ P (AKS, 2002): Testing if a single number is prime
  - π(x) UNKNOWN complexity: Counting primes up to x
  - No proof that π(x) is NOT in P (no lower bound beyond Ω(log x))
  - No proof that π(x) IS in P either

The problem is OPEN. Unlike the "provably impossible" claim from
previous sessions, the computational complexity of π(x) and p(n)
is genuinely UNKNOWN.

WHAT WE KNOW FOR CERTAIN:
  - The explicit formula approach needs Ω(√x) work (proven, GUE barrier)
  - The combinatorial sieve approach needs Ω(x^{1/3}) work (conjectured)
  - Both are exponential in the input size
  - No approach avoiding these two frameworks is known
  - But NO unconditional lower bound prevents a polylog algorithm

OPEN QUESTION: Does p(n) have polynomial-in-input-size circuit complexity?
This is equivalent to asking: is there a Boolean circuit of size polylog(n)
that computes p(n)?
""")

# ============================================================
# PART 5: Practical test — can we beat known methods even slightly?
# ============================================================
print("=" * 60)
print("PART 5: Attempting Slightly Better Than Known Methods")
print("=" * 60)

# The best we can do is R^{-1}(n) + local search using primality testing.
# R^{-1}(n) has error ~√p(n)/log(p(n)). We need to search within this window.
# Primality testing each candidate: O(log^6(p)) per test (AKS) or O(log^2(p)) (Miller-Rabin under GRH)
# Window size: ~2√p/log(p)
# Total candidates in window: ~2√p/(log p)²
#
# Total work: O(√p/log²p × log²p) = O(√p) — same as analytic method!
#
# So even the "brute force within the window" approach is O(√p) = O(√(n ln n)).

print("R^{-1}(n) + local primality testing:")
for log_n in [10, 20, 50, 100]:
    n = 10**log_n
    p_approx = n * log_n * math.log(10)
    window = 2 * math.sqrt(p_approx) / math.log(p_approx)
    candidates = window / math.log(p_approx)
    test_cost = (log_n * math.log(10))**2  # per test (Miller-Rabin under GRH)
    total = candidates * test_cost
    print(f"  n=10^{log_n}: window={window:.2e}, candidates={candidates:.2e}, "
          f"total work={total:.2e}")

print()
print("This is O(√(n·ln(n))) = O(n^{1/2+ε}) — exponential in input size.")
print("No improvement over direct methods.")

print("\n" + "=" * 60)
print("CONCLUSION")
print("=" * 60)
print("""
The computational complexity of p(n) is:
  - UPPER BOUND: O(2^{2s/3}) where s = log₂(n) [Meissel-Lehmer]
  - LOWER BOUND: Ω(s) [trivial — must read input]
  - Gap: EXPONENTIAL

This gap is the fundamental OPEN PROBLEM.
No approach in 375+ attempts has narrowed it.
But equally, no proof has shown it can't be narrowed.

The problem of finding p(n) in polylog(n) time remains OPEN.
""")

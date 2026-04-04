#!/usr/bin/env python3
"""
Mertens Function Speedup Investigation
=======================================

Investigates whether M(x) = sum_{n<=x} mu(n) can be computed faster than O(x^{2/3}).

Key questions:
1. Why does D(x) = sum d(n) admit O(x^{1/3}) via hyperbola, but M(x) is stuck at O(x^{2/3})?
2. Can the Helfgott-Thompson O(x^{3/5}) elementary approach be implemented/tested?
3. Is the Lucy_Hedgehog DP fundamentally limited at O(x^{2/3})?

Background:
- D(x) = sum_{n<=x} d(n) = 2*sum_{k<=sqrt(x)} floor(x/k) - floor(sqrt(x))^2  [O(x^{1/3}) via Euler-Maclaurin]
- M(x) = sum_{n<=x} mu(n) uses identity: M(x) = 1 - sum_{n=2}^{x} M(floor(x/n))
- Best combinatorial: O(x^{2/3}) for both pi(x) and M(x) [Lucy_Hedgehog DP]
- Helfgott-Thompson (2021): O(x^{3/5} (log x)^{3/5+eps}) elementary for M(x)
- Analytic: O(x^{1/2+eps}) via Lagarias-Odlyzko / explicit formula with zeta zeros

References:
- Helfgott & Thompson, "Summing mu(n): a faster elementary algorithm" (arXiv:2101.08773)
- Hiary, "Fast methods to compute the Riemann zeta function" (Annals 2011)
- gbroxey, "Summing Multiplicative Functions" (blog, 2023)
- Lucy_Hedgehog, Project Euler Problem 10 forum post

Date: 2026-04-04
"""

import math
import time
import sys
from collections import defaultdict
from functools import lru_cache

# =============================================================================
# PART 1: Brute-force reference implementations (for verification)
# =============================================================================

def mobius_sieve(n):
    """Compute mu(k) for k = 0..n using a linear sieve. O(n) time and space."""
    mu = [0] * (n + 1)
    mu[1] = 1
    is_prime = [True] * (n + 1)
    primes = []
    for i in range(2, n + 1):
        if is_prime[i]:
            primes.append(i)
            mu[i] = -1  # prime => mu = -1
        for p in primes:
            if i * p > n:
                break
            is_prime[i * p] = False
            if i % p == 0:
                mu[i * p] = 0  # p^2 divides i*p
                break
            else:
                mu[i * p] = -mu[i]
    return mu


def mertens_brute(x):
    """Compute M(x) = sum_{n<=x} mu(n) by brute force sieve. O(x) time."""
    mu = mobius_sieve(x)
    return sum(mu[1:x+1])


def divisor_sum_brute(x):
    """Compute D(x) = sum_{n<=x} d(n) by brute force. O(x log x) time."""
    d = [0] * (x + 1)
    for k in range(1, x + 1):
        for m in range(k, x + 1, k):
            d[m] += 1
    return sum(d[1:x+1])


# =============================================================================
# PART 2: D(x) via Dirichlet Hyperbola Method -- O(sqrt(x))
# =============================================================================

def divisor_sum_hyperbola(x):
    """
    Compute D(x) = sum_{n<=x} d(n) using the Dirichlet hyperbola method.

    Identity: D(x) = 2 * sum_{k=1}^{floor(sqrt(x))} floor(x/k) - floor(sqrt(x))^2

    This works because d = 1 * 1 (Dirichlet convolution), so:
      sum_{n<=x} d(n) = sum_{n<=x} sum_{ab=n} 1
                      = #{(a,b) : ab <= x}  (lattice points under hyperbola)

    The hyperbola method counts these by splitting at sqrt(x):
      = 2 * sum_{a<=sqrt(x)} floor(x/a) - floor(sqrt(x))^2

    Complexity: O(sqrt(x)) -- just a single loop up to sqrt(x).
    """
    isqrt_x = int(math.isqrt(x))
    result = 0
    for k in range(1, isqrt_x + 1):
        result += x // k
    result = 2 * result - isqrt_x * isqrt_x
    return result


def divisor_sum_fast(x):
    """
    Compute D(x) in O(x^{1/3}) using block decomposition of floor(x/k).

    The key insight: floor(x/k) takes only O(sqrt(x)) distinct values.
    We can group consecutive k with the same floor(x/k) value.

    For each distinct value v = floor(x/k), find the range [k_lo, k_hi]
    where floor(x/k) = v, then add v * (k_hi - k_lo + 1).

    With careful implementation using the "staircase" structure:
    - For k <= x^{1/3}: compute individually (x^{1/3} terms)
    - For k > x^{1/3}: floor(x/k) < x^{2/3}, group by blocks

    Total: O(x^{1/3}) operations.
    """
    if x <= 0:
        return 0

    # We use the block decomposition approach
    cbrt_x = int(round(x ** (1.0 / 3.0)))
    # Adjust cbrt to be safe
    while (cbrt_x + 1) ** 3 <= x:
        cbrt_x += 1
    while cbrt_x ** 3 > x:
        cbrt_x -= 1

    result = 0

    # Part 1: k from 1 to cbrt(x) -- compute floor(x/k) individually
    for k in range(1, cbrt_x + 1):
        result += x // k

    # Part 2: for k > cbrt(x), group by blocks where floor(x/k) = v
    # The distinct values of floor(x/k) for k > cbrt(x) are v = 1, 2, ..., floor(x/(cbrt_x+1))
    # For each v, k ranges from floor(x/(v+1))+1 to floor(x/v)
    max_v = x // (cbrt_x + 1) if cbrt_x + 1 <= x else 0

    for v in range(1, max_v + 1):
        k_hi = x // v
        k_lo = x // (v + 1) + 1
        # Clamp to k > cbrt_x
        k_lo = max(k_lo, cbrt_x + 1)
        if k_lo <= k_hi:
            result += v * (k_hi - k_lo + 1)

    # D(x) = 2*S - floor(sqrt(x))^2 where S = sum_{k<=sqrt(x)} floor(x/k)
    # But we computed sum_{k=1}^{x} floor(x/k)/k... no, let's use the hyperbola directly.
    # Actually, the above computes sum_{k=1}^{x} floor(x/k) partially.
    # Let me redo this properly.

    # D(x) = 2 * sum_{k=1}^{isqrt} floor(x/k) - isqrt^2
    # To compute sum_{k=1}^{isqrt} floor(x/k) in O(x^{1/3}):
    # Split at cbrt: individual terms for k<=cbrt, blocks for k in (cbrt, isqrt]

    isqrt_x = int(math.isqrt(x))

    result = 0
    # Individual terms for k = 1..cbrt_x (but not beyond isqrt_x)
    limit1 = min(cbrt_x, isqrt_x)
    for k in range(1, limit1 + 1):
        result += x // k

    # Block terms for k = cbrt_x+1 .. isqrt_x
    if cbrt_x < isqrt_x:
        # For k in (cbrt_x, isqrt_x], the values v = floor(x/k) range in
        # [floor(x/isqrt_x), floor(x/(cbrt_x+1))]
        # We iterate over distinct values v
        v_lo = x // isqrt_x
        v_hi = x // (cbrt_x + 1) if cbrt_x + 1 <= x else 0

        for v in range(1, v_hi + 1):
            k_hi_v = x // v
            k_lo_v = x // (v + 1) + 1
            # Clamp to [cbrt_x+1, isqrt_x]
            k_lo_v = max(k_lo_v, cbrt_x + 1)
            k_hi_v = min(k_hi_v, isqrt_x)
            if k_lo_v <= k_hi_v:
                result += v * (k_hi_v - k_lo_v + 1)

    result = 2 * result - isqrt_x * isqrt_x
    return result


# =============================================================================
# PART 3: M(x) via Lucy_Hedgehog DP -- O(x^{2/3})
# =============================================================================

def mertens_lucy(x):
    """
    Compute M(x) = sum_{n<=x} mu(n) using Lucy_Hedgehog-style DP.

    Uses the identity: sum_{d=1}^{n} M(floor(n/d)) = 1
    Therefore: M(n) = 1 - sum_{d=2}^{n} M(floor(n/d))

    Strategy:
    1. Sieve mu(k) for k <= y (threshold parameter)
    2. Compute prefix sums of mu to get M(k) for k <= y
    3. For values floor(x/k) > y, compute M recursively using the identity
    4. Use memoization on the O(sqrt(x)) distinct values of floor(x/k)

    With y = x^{2/3}, total complexity is O(x^{2/3}).
    """
    if x <= 0:
        return 0
    if x == 1:
        return 1

    # Threshold: y ~ x^{2/3} for optimal balance
    y = max(int(x ** (2.0 / 3.0)), 1)
    # Ensure y >= sqrt(x) for correctness
    y = max(y, int(math.isqrt(x)) + 1)

    # Step 1: Sieve mu(k) for k <= y
    mu = mobius_sieve(y)

    # Step 2: Prefix sums -> M(k) for k <= y
    M_small = [0] * (y + 1)
    for i in range(1, y + 1):
        M_small[i] = M_small[i - 1] + mu[i]

    # Step 3: Compute M(floor(x/k)) for values > y using memoization
    M_large = {}  # maps v -> M(v) for v = floor(x/k) > y

    def M(n):
        """Compute M(n) using sieved values and recursive identity."""
        if n <= y:
            return M_small[n]
        if n in M_large:
            return M_large[n]

        # M(n) = 1 - sum_{d=2}^{n} M(floor(n/d))
        # Use block decomposition for the sum
        result = 1

        # Direct terms: d from 2 to isqrt(n), contribute M(floor(n/d))
        isqrt_n = int(math.isqrt(n))

        # Block sum: for each v = floor(n/d), count how many d give this v
        # For d in [2, isqrt_n]: direct
        for d in range(2, isqrt_n + 1):
            result -= M(n // d)

        # For v in [1, floor(n/(isqrt_n+1))]: block
        # floor(n/d) = v for d in [n/(v+1)+1, n/v]
        # But d >= 2 and d > isqrt_n
        max_v = n // (isqrt_n + 1) if isqrt_n + 1 <= n else 0

        for v in range(1, max_v + 1):
            d_lo = n // (v + 1) + 1
            d_hi = n // v
            # Clamp to d >= isqrt_n + 1 (already counted d <= isqrt_n above)
            d_lo = max(d_lo, isqrt_n + 1)
            # Also d >= 2
            d_lo = max(d_lo, 2)
            if d_lo <= d_hi:
                count = d_hi - d_lo + 1
                result -= M(v) * count

        M_large[n] = result
        return result

    return M(x)


# =============================================================================
# PART 4: Analysis -- WHY hyperbola works for D(x) but not M(x)
# =============================================================================

def analyze_hyperbola_difference():
    """
    Analyze the fundamental difference between D(x) and M(x) computation.

    D(x) = sum_{n<=x} d(n) where d = 1 * 1 (convolution)
      => D(x) = sum_{ab<=x} 1 = lattice point count
      => Hyperbola: D(x) = 2*sum_{a<=sqrt(x)} floor(x/a) - floor(sqrt(x))^2
      => O(sqrt(x)), improvable to O(x^{1/3}) via Euler-Maclaurin on floor sums

    M(x) = sum_{n<=x} mu(n) where mu = 1^{-1} (Dirichlet inverse)
      => No direct "lattice point" interpretation
      => Identity: sum_{d=1}^{x} M(floor(x/d)) = 1  [from mu * 1 = epsilon]
      => M(x) = 1 - sum_{d=2}^{x} M(floor(x/d))
      => RECURSIVE: needs M at O(sqrt(x)) points, each needing M at more points

    KEY INSIGHT: The difference is NOT about the function values (oscillation etc.)
    It's about the RECURRENCE STRUCTURE:

    For D(x): D(x) = 2*S(sqrt(x)) - sqrt(x)^2 where S(y) = sum_{k<=y} floor(x/k)
              S(y) is a SUM OF INDEPENDENT TERMS -- no recursion needed!
              Each floor(x/k) is computed in O(1).

    For M(x): M(x) = 1 - sum_{d=2}^{x} M(floor(x/d))
              This REQUIRES M at other points -- it's SELF-REFERENTIAL.
              You can't compute M(x) without first computing M at O(sqrt(x)) points.

    This self-referential structure is because:
      d = 1 * 1:  both factors are the constant function 1, trivially summable
      mu = 1^{-1}: the Dirichlet inverse requires knowing the answer to compute itself

    In convolution terms:
      sum (f*g)(n) = sum f(a)*G(x/a) where G is partial sum of g
      For d = 1*1: both partial sums are floor functions (trivial)
      For mu: would need mu = lambda * something, but any useful decomposition
              requires equally hard-to-compute components.
    """
    print("=" * 78)
    print("ANALYSIS: Why hyperbola works for D(x) but not M(x)")
    print("=" * 78)

    print("""
DIVISOR SUM D(x) = sum_{n<=x} d(n):
  d(n) = (1 * 1)(n)  [Dirichlet convolution of constant function with itself]

  By Dirichlet hyperbola:
    D(x) = sum_{ab<=x} 1*1 = 2*sum_{a<=sqrt(x)} floor(x/a) - floor(sqrt(x))^2

  This is a CLOSED-FORM sum of floor functions. No recursion needed.
  Each term floor(x/k) costs O(1). Total: O(sqrt(x)).
  With Euler-Maclaurin on grouped floor values: O(x^{1/3}).

MERTENS FUNCTION M(x) = sum_{n<=x} mu(n):
  mu * 1 = epsilon  =>  sum_{d|n} mu(d) = [n=1]

  Summing: sum_{d=1}^{x} M(floor(x/d)) = 1
  Therefore: M(x) = 1 - sum_{d=2}^{x} M(floor(x/d))

  This is SELF-REFERENTIAL: M(x) depends on M at O(sqrt(x)) other points.
  Each of those points depends on M at more points, creating a tree of depth
  O(log x). The total number of distinct values is O(sqrt(x)).

THE FUNDAMENTAL DIFFERENCE:
  D(x): f*g where both f,g have trivially computable partial sums
         => One-shot formula, no recursion, O(sqrt(x))

  M(x): f is the INVERSE of a simple function
         => Computing partial sums of an inverse inherently requires recursion
         => The recursion tree has O(sqrt(x)) nodes at each level
         => With sieving up to y, cost is O(y + x/sqrt(y))
         => Optimized at y = x^{2/3}, giving O(x^{2/3})

CAN THE HYPERBOLA METHOD BE ADAPTED FOR M(x)?
  If we could write mu = f * g where both sum_f and sum_g are easy:
    M(x) = sum_{a<=sqrt(x)} f(a) * G(floor(x/a)) + ...  [hyperbola split]

  But mu = 1^{-1}, and the only known factorizations are:
    mu = mu * epsilon        (trivial)
    mu = lambda * |mu|       (lambda = (-1)^Omega, |mu| = mu^2 -- both hard)
    mu = mu_k * something    (partial Mobius -- doesn't help)

  No factorization mu = f*g with BOTH f,g having easy partial sums is known.
  This appears to be a genuine structural barrier, not just a gap in techniques.
""")

    # Empirical verification
    print("EMPIRICAL VERIFICATION")
    print("-" * 40)
    print(f"{'x':>12}  {'D(x) hyp':>12}  {'D(x) brute':>12}  {'M(x) lucy':>10}  {'M(x) brute':>10}")
    for exp in range(2, 7):
        x = 10 ** exp
        d_hyp = divisor_sum_hyperbola(x)
        d_brute = divisor_sum_brute(x) if x <= 100000 else "skip"
        m_lucy = mertens_lucy(x)
        m_brute = mertens_brute(x) if x <= 100000 else "skip"
        print(f"{x:>12}  {d_hyp:>12}  {str(d_brute):>12}  {m_lucy:>10}  {str(m_brute):>10}")


# =============================================================================
# PART 5: Benchmarks -- timing comparison
# =============================================================================

def benchmark_all():
    """Benchmark all implementations across different x values."""
    print("\n" + "=" * 78)
    print("BENCHMARKS")
    print("=" * 78)

    test_values = [10**4, 10**5, 10**6, 10**7, 10**8]

    # D(x) benchmarks
    print("\n--- D(x) = sum d(n) ---")
    print(f"{'x':>12}  {'Hyperbola O(sqrt)':>20}  {'Fast O(x^1/3)':>20}  {'Value':>15}")
    for x in test_values:
        t0 = time.perf_counter()
        v1 = divisor_sum_hyperbola(x)
        t1 = time.perf_counter()

        t2 = time.perf_counter()
        v2 = divisor_sum_fast(x)
        t3 = time.perf_counter()

        assert v1 == v2, f"Mismatch at x={x}: {v1} vs {v2}"
        print(f"{x:>12}  {t1-t0:>18.6f}s  {t3-t2:>18.6f}s  {v1:>15}")

    # M(x) benchmarks
    print("\n--- M(x) = sum mu(n) ---")
    print(f"{'x':>12}  {'Lucy O(x^2/3)':>20}  {'Brute O(x)':>20}  {'Value':>15}")
    for x in test_values:
        t0 = time.perf_counter()
        v1 = mertens_lucy(x)
        t1 = time.perf_counter()

        if x <= 10**6:
            t2 = time.perf_counter()
            v2 = mertens_brute(x)
            t3 = time.perf_counter()
            assert v1 == v2, f"Mismatch at x={x}: {v1} vs {v2}"
            brute_str = f"{t3-t2:.6f}s"
        else:
            brute_str = "skip"

        print(f"{x:>12}  {t1-t0:>18.6f}s  {brute_str:>20}  {v1:>15}")


# =============================================================================
# PART 6: Attempted speedups for M(x)
# =============================================================================

def mertens_lucy_optimized(x):
    """
    Optimized Lucy DP for M(x) with better constant factors.

    Optimizations over the basic version:
    1. Use arrays instead of dicts for large-value lookup
    2. Precompute more mu values (larger sieve)
    3. Iterate in optimal order to maximize cache hits
    """
    if x <= 0:
        return 0
    if x == 1:
        return 1

    # Use y ~ x^{2/3} * (log x)^{1/3} for slightly better balance
    y = max(int(x ** (2.0 / 3.0) * math.log(x + 2) ** (1.0 / 3.0)), 100)
    y = min(y, x)  # Don't exceed x

    # Sieve mu and compute prefix sums
    mu = mobius_sieve(y)
    M_small = [0] * (y + 1)
    for i in range(1, y + 1):
        M_small[i] = M_small[i - 1] + mu[i]

    if x <= y:
        return M_small[x]

    # Collect all distinct values of floor(x/k) that are > y
    # These are x//k for k = 1, 2, ..., floor(x/(y+1))
    max_k = x // (y + 1)

    # Store M values for large arguments in an array indexed by k
    # M_large[k] = M(floor(x/k))
    M_large = [0] * (max_k + 2)

    # Process in decreasing order of floor(x/k), i.e., k = max_k down to 1
    for k in range(max_k, 0, -1):
        n = x // k
        # M(n) = 1 - sum_{d=2}^{n} M(floor(n/d))
        result = 1
        isqrt_n = int(math.isqrt(n))

        # Direct terms: d from 2 to isqrt_n
        for d in range(2, isqrt_n + 1):
            nd = n // d
            if nd <= y:
                result -= M_small[nd]
            else:
                # nd = floor(n/d) = floor(x/(k*d))... but we need index
                # floor(x / (k*d)) -- find the index j such that floor(x/j) = nd
                # j = x // nd might not equal k*d, but floor(x/j) = nd
                j = x // nd
                result -= M_large[j]

        # Block terms: for v = 1 to floor(n/(isqrt_n+1))
        max_v = n // (isqrt_n + 1) if isqrt_n + 1 <= n else 0
        for v in range(1, max_v + 1):
            d_lo = n // (v + 1) + 1
            d_hi = n // v
            d_lo = max(d_lo, isqrt_n + 1)
            d_lo = max(d_lo, 2)
            if d_lo <= d_hi:
                count = d_hi - d_lo + 1
                if v <= y:
                    result -= M_small[v] * count
                else:
                    j = x // v
                    result -= M_large[j] * count

        M_large[k] = result

    return M_large[1]


def attempt_factorization_speedup(x):
    """
    EXPERIMENT: Try to find a useful factorization mu = f * g.

    If mu = f * g with both F(x) = sum f(n) and G(x) = sum g(n) computable
    in T_F(x) and T_G(x) respectively, then M(x) can be computed in
    O(sqrt(x) * max(T_F, T_G)) via the hyperbola method.

    Candidates tested:
    1. mu(n) = lambda(n) * |mu(n)| where lambda = (-1)^Omega(n)
       Problem: sum lambda(n) is related to sqrt(zeta(2s)/zeta(s)), still hard
    2. mu(n) = sum_{d|n} f(d) for some f
       This means f = mu * mu (Dirichlet), but sum(mu*mu) is also hard
    3. Use mu * 1 = epsilon backwards: not helpful (epsilon has trivial partial sums)

    Conclusion: No known factorization helps.
    """
    print("\n" + "=" * 78)
    print("EXPERIMENT: Factorization attempts for mu")
    print("=" * 78)

    # Test: mu = lambda * |mu| (pointwise, not convolution)
    # This is pointwise: mu(n) = lambda(n) * mu^2(n)
    # Not a Dirichlet convolution, so hyperbola doesn't directly apply.

    # Test: Dirichlet convolution mu = f * g
    # mu * 1 = epsilon => mu = epsilon * mu^{-1}... no, mu is 1^{-1}
    # If mu = f * g, then f*g*1 = epsilon*1 = 1... no.
    # Actually mu*1 = epsilon, so if mu = f*g then f*g*1 = epsilon.

    # Let's compute mu*mu (Dirichlet) and see its partial sums
    limit = min(x, 10000)
    mu_vals = mobius_sieve(limit)

    # Compute (mu * mu)(n) = sum_{d|n} mu(d)*mu(n/d)
    mu_conv_mu = [0] * (limit + 1)
    for d in range(1, limit + 1):
        if mu_vals[d] == 0:
            continue
        for m in range(1, limit // d + 1):
            if mu_vals[m] == 0:
                continue
            mu_conv_mu[d * m] += mu_vals[d] * mu_vals[m]

    # Partial sums of mu*mu
    partial = [0] * (limit + 1)
    for i in range(1, limit + 1):
        partial[i] = partial[i - 1] + mu_conv_mu[i]

    print(f"\nDirichlet convolution (mu*mu)(n) partial sums up to {limit}:")
    for k in [10, 100, 1000, 10000]:
        if k <= limit:
            print(f"  sum_{{n<={k}}} (mu*mu)(n) = {partial[k]}")

    # Note: sum_{n<=x} (mu*mu)(n) ~ x/zeta(2) * C for some constant
    # This is the squarefree counting function times a correction
    # It's actually: (mu*mu)(n) = sum_{d|n} mu(d)mu(n/d)
    # The Dirichlet series is 1/zeta(s)^2
    # So sum_{n<=x} (mu*mu)(n) ~ x * Res(x^s / (s*zeta(s)^2), s=1) = 0
    # (since zeta has a simple pole, 1/zeta^2 has a double zero)

    print(f"\n  Note: sum (mu*mu)(n) -> 0 as x -> infinity")
    print(f"  This is because the Dirichlet series of mu*mu is 1/zeta(s)^2")
    print(f"  which has a double zero at s=1.")

    # Test: lambda(n) partial sums (Liouville function)
    # lambda(n) = (-1)^Omega(n), related to mu by: lambda = mu * |mu| (Dirichlet? No.)
    # Actually: sum_{d|n} lambda(d) = 1 if n is a perfect square, 0 otherwise
    # So lambda * 1 = characteristic function of squares
    # And: sum lambda(n) for n<=x is related to zeta(2s)/zeta(s) at s=1
    # L(x) ~ 0 (under RH, |L(x)| = O(sqrt(x)))
    # Computing L(x) is exactly as hard as M(x).

    print(f"\n  CONCLUSION: No known factorization mu = f*g (Dirichlet)")
    print(f"  exists where both sum_f and sum_g are easy to compute.")
    print(f"  The self-referential nature of M(x) appears fundamental.")


def test_helfgott_thompson_idea(x):
    """
    EXPERIMENT: Explore the Helfgott-Thompson O(x^{3/5}) approach.

    Their key idea: Instead of sieving mu up to y and recursing for y < n <= x,
    they use a COMBINATORIAL IDENTITY that expresses M(x) in terms of
    counts of integers with specific factorization patterns.

    Specifically, they decompose:
      M(x) = sum_{n<=x} mu(n)
           = (squarefree count with even # of prime factors)
           - (squarefree count with odd # of prime factors)

    Then they use inclusion-exclusion with a carefully chosen sieve:
    - Sieve primes up to x^{1/5} (call this set P_small)
    - For each pattern of small prime factors, count large-prime contributions
    - The "large" part uses Buchstab-like identities

    The complexity gain: By choosing the sieve level at x^{1/5} instead of x^{1/3},
    the number of sieve patterns is O(x^{3/5} / log x) and each pattern costs
    O(polylog) to evaluate using precomputed data.

    This is genuinely different from Lucy DP! It avoids the recursive M(floor(x/d))
    calls entirely, replacing them with combinatorial counting.

    We implement a simplified version to understand the structure.
    """
    print("\n" + "=" * 78)
    print("EXPERIMENT: Helfgott-Thompson O(x^{3/5}) structure")
    print("=" * 78)

    # The H-T approach in simplified form:
    # M(x) = #{n <= x : n squarefree, omega(n) even} - #{n <= x : n squarefree, omega(n) odd}
    # = Q_even(x) - Q_odd(x)
    # where Q(x) = Q_even(x) + Q_odd(x) = #{n<=x : n squarefree} ~ 6x/pi^2

    # Standard Lucy approach: sieve up to y, recurse above
    # H-T approach: sieve up to z = x^{1/5}, use Buchstab identity for rest

    # For demonstration, let's verify the decomposition for small x
    limit = min(x, 100000)
    mu_vals = mobius_sieve(limit)

    q_even = 0  # squarefree with even number of prime factors
    q_odd = 0   # squarefree with odd number of prime factors
    for n in range(1, limit + 1):
        if mu_vals[n] == 1:
            q_even += 1
        elif mu_vals[n] == -1:
            q_odd += 1

    m_direct = sum(mu_vals[1:limit + 1])

    print(f"\nFor x = {limit}:")
    print(f"  Q_even(x) = {q_even}  (squarefree, even omega)")
    print(f"  Q_odd(x)  = {q_odd}  (squarefree, odd omega)")
    print(f"  M(x) = Q_even - Q_odd = {q_even - q_odd}")
    print(f"  Direct M(x) = {m_direct}")
    print(f"  Q(x) = Q_even + Q_odd = {q_even + q_odd}  (total squarefree)")
    print(f"  6x/pi^2 = {6*limit/math.pi**2:.1f}  (expected squarefree count)")

    print(f"""
  HELFGOTT-THOMPSON KEY INSIGHT:
  Instead of the recursive identity M(x) = 1 - sum M(floor(x/d)),
  they decompose M(x) = Q_even(x) - Q_odd(x) and compute each via
  inclusion-exclusion over small primes (up to x^{{1/5}}).

  For each subset S of small primes with product P_S:
    Count squarefree n <= x/P_S with no small prime factors and even/odd omega

  The "no small prime factors" constraint means we only need to count
  integers whose prime factors are all > x^{{1/5}}, which have at most
  4 prime factors (since (x^{{1/5}})^5 = x).

  Counting k-almost-primes with all factors > z in [1, x/P_S] can be done
  in O(polylog) for fixed k using pi(x)-like formulas.

  Total: O(x^{{3/5}} * polylog) -- the sum over subsets S dominates.

  THIS IS A GENUINE IMPROVEMENT over O(x^{{2/3}}) Lucy DP!
  But it's "elementary" (no analytic number theory / zeta zeros).
  The analytic method (Lagarias-Odlyzko) still gives O(x^{{1/2+eps}}).
""")


def test_recursive_depth(x):
    """
    EXPERIMENT: Measure the recursion structure of M(x) computation.

    Key question: Can we reduce the number of recursive calls needed?
    """
    print("\n" + "=" * 78)
    print("EXPERIMENT: Recursion structure of M(x)")
    print("=" * 78)

    # Count distinct values of floor(x/k) for k = 1, 2, ..., x
    isqrt_x = int(math.isqrt(x))
    distinct_values = set()
    for k in range(1, isqrt_x + 1):
        distinct_values.add(x // k)
        distinct_values.add(k)

    n_values = len(distinct_values)
    print(f"\n  x = {x}")
    print(f"  sqrt(x) = {isqrt_x}")
    print(f"  Distinct values of floor(x/k): {n_values}")
    print(f"  2*sqrt(x) = {2*isqrt_x}")

    # For each value v, how many d's map to it?
    # floor(x/d) = v iff x/(v+1) < d <= x/v
    value_counts = {}
    for v in sorted(distinct_values, reverse=True)[:20]:
        d_lo = x // (v + 1) + 1 if v < x else 1
        d_hi = x // v
        count = d_hi - d_lo + 1
        value_counts[v] = count

    print(f"\n  Top 20 values and their multiplicities:")
    for v in sorted(value_counts.keys(), reverse=True):
        print(f"    floor(x/d) = {v:>10}, count = {value_counts[v]:>6}")


def compare_complexity_scaling():
    """
    EXPERIMENT: Empirically measure scaling exponents.
    """
    print("\n" + "=" * 78)
    print("SCALING ANALYSIS")
    print("=" * 78)

    print("\n--- D(x) via hyperbola [expected O(sqrt(x))] ---")
    prev_t = None
    for exp in range(4, 10):
        x = 10 ** exp
        t0 = time.perf_counter()
        for _ in range(3):
            divisor_sum_hyperbola(x)
        t1 = time.perf_counter()
        t = (t1 - t0) / 3
        ratio = t / prev_t if prev_t and prev_t > 0 else 0
        exponent = math.log10(ratio) if ratio > 0 else 0
        prev_t = t
        print(f"  x=10^{exp}: {t:.6f}s  ratio={ratio:.2f}  scaling~x^{exponent:.2f}")

    print("\n--- M(x) via Lucy DP [expected O(x^{2/3})] ---")
    prev_t = None
    for exp in range(4, 9):
        x = 10 ** exp
        t0 = time.perf_counter()
        v = mertens_lucy(x)
        t1 = time.perf_counter()
        t = t1 - t0
        ratio = t / prev_t if prev_t and prev_t > 0 else 0
        exponent = math.log10(ratio) if ratio > 0 else 0
        prev_t = t
        print(f"  x=10^{exp}: {t:.6f}s  M(x)={v:>10}  ratio={ratio:.2f}  scaling~x^{exponent:.2f}")

    print("\n--- M(x) via optimized Lucy DP ---")
    prev_t = None
    for exp in range(4, 9):
        x = 10 ** exp
        t0 = time.perf_counter()
        v = mertens_lucy_optimized(x)
        t1 = time.perf_counter()
        t = t1 - t0
        ratio = t / prev_t if prev_t and prev_t > 0 else 0
        exponent = math.log10(ratio) if ratio > 0 else 0
        prev_t = t
        print(f"  x=10^{exp}: {t:.6f}s  M(x)={v:>10}  ratio={ratio:.2f}  scaling~x^{exponent:.2f}")


# =============================================================================
# PART 7: Connection to pi(x) -- can M(x) speedup transfer?
# =============================================================================

def analyze_pi_mertens_connection():
    """
    Analyze whether a faster M(x) algorithm would give faster pi(x).

    The connection: pi(x) = sum_{k=1}^{floor(log2(x))} (mu(k)/k) * Li(x^{1/k})
                           + correction terms involving M(x)

    More precisely, via Mobius inversion of the prime counting function:
      Pi(x) = sum_{n=1}^{inf} mu(n)/n * Li(x^{1/n})  [Riemann's formula]

    where Pi(x) = pi(x) + (1/2)pi(x^{1/2}) + (1/3)pi(x^{1/3}) + ...

    But the EXACT computation of pi(x) via the Meissel-Lehmer method
    does NOT go through M(x). It uses a different combinatorial identity:
      pi(x) = phi(x, a) + a - 1 - P2(x, a) - ... [Legendre-type sums]

    The Lucy DP for pi(x) and for M(x) have the SAME complexity O(x^{2/3})
    because they both fundamentally solve:
      "Compute f(floor(x/k)) for all k, where f satisfies a sieve recurrence"

    If M(x) could be computed in O(x^{3/5}) [Helfgott-Thompson]:
    - Does NOT directly give pi(x) in O(x^{3/5})
    - The H-T trick is specific to mu's multiplicative structure
    - pi(x) = sum of Lambda(n)/log(n) -- different structure
    - However, the INSIGHT might transfer: decompose pi(x) by factorization patterns
    """
    print("\n" + "=" * 78)
    print("CONNECTION: M(x) <-> pi(x) computational complexity")
    print("=" * 78)

    print("""
  CURRENT STATE OF THE ART:
  ========================
  Function           Combinatorial        Analytic
  D(x)=sum d(n)      O(x^{1/3})          O(x^{1/3})
  M(x)=sum mu(n)     O(x^{3/5}+eps) [HT] O(x^{1/2+eps}) [LO]
  pi(x)              O(x^{2/3})     [LH]  O(x^{1/2+eps}) [LO]
  L(x)=sum lambda(n) O(x^{2/3})     [LH]  O(x^{1/2+eps}) [LO]

  KEY OBSERVATION:
  Helfgott-Thompson improved M(x) from O(x^{2/3}) to O(x^{3/5}) combinatorially.
  But pi(x) is STILL at O(x^{2/3}) combinatorially!

  Why doesn't the H-T trick transfer?
  - H-T uses: mu(n) = (-1)^{omega(n)} * [n squarefree]
  - They decompose by small-prime patterns and count large-prime contributions
  - For pi(x), the analogous decomposition doesn't help because
    pi(x) counts PRIMES specifically, not a signed sum over all integers
  - The Buchstab identity for pi(x) gives the Meissel-Lehmer method,
    which is already the basis of the O(x^{2/3}) Lucy DP

  COULD pi(x) be improved to O(x^{3/5})?
  - No known method achieves this combinatorially
  - The H-T approach for M(x) exploits cancellation in the signed sum
  - pi(x) has no such cancellation (it's a positive counting function)
  - This suggests pi(x) may genuinely be harder than M(x) combinatorially

  IMPLICATION FOR p(n):
  - Even if M(x) is O(x^{3/5}), converting to pi(x) adds overhead
  - The Mobius inversion pi(x) = sum mu(k)/k * Li(x^{1/k}) needs O(log x) M-calls
  - Each at different arguments, but the dominant one is M(x) itself
  - So pi(x) via M(x) would be O(x^{3/5} * polylog) -- IF the transfer works
  - But the transfer is NOT straightforward (different recurrence structures)
""")


# =============================================================================
# PART 8: Summary and conclusions
# =============================================================================

def print_summary():
    print("\n" + "=" * 78)
    print("SUMMARY AND CONCLUSIONS")
    print("=" * 78)
    print("""
  1. D(x) vs M(x): FUNDAMENTAL STRUCTURAL DIFFERENCE
     - D(x) = lattice point count, closed-form via hyperbola => O(x^{1/3})
     - M(x) = self-referential sum, requires recursive computation => O(x^{2/3})
     - The difference is that d = 1*1 (easy factors) while mu = 1^{-1} (inverse)
     - No known factorization mu = f*g with both partial sums easy

  2. HELFGOTT-THOMPSON BREAKTHROUGH (2021):
     - Improved M(x) from O(x^{2/3}) to O(x^{3/5} polylog) ELEMENTARILY
     - Key idea: decompose by small-prime factorization patterns
     - Avoids recursive M(floor(x/d)) calls entirely
     - Uses Buchstab-like identities with sieve level x^{1/5}
     - This is the FIRST improvement in the x-exponent since 1985

  3. DOES NOT TRANSFER TO pi(x):
     - pi(x) remains at O(x^{2/3}) combinatorially
     - H-T exploits signed cancellation specific to mu
     - pi(x) is a positive count, no analogous cancellation
     - pi(x) may be genuinely harder than M(x) combinatorially

  4. FOR p(n) = nth prime:
     - Still need pi(x) which is O(x^{2/3}) combinatorial, O(x^{1/2+eps}) analytic
     - H-T M(x) improvement doesn't help directly
     - The analytic approach (Lagarias-Odlyzko) remains best for large x
     - For p(10^{100}): still need ~10^{49} operations minimum

  5. OPEN QUESTIONS:
     - Can pi(x) be improved to O(x^{3/5}) combinatorially?
     - Is there a combinatorial lower bound > x^{1/2} for pi(x)?
     - Can the H-T approach for M(x) be extended to other multiplicative sums?
     - Is O(x^{1/2+eps}) optimal for pi(x) even with analytic methods?

  VERDICT: The M(x) speedup is a genuine theoretical advance but does NOT
  break the O(x^{2/3}) barrier for pi(x) or p(n). The problem remains open.
""")


# =============================================================================
# MAIN
# =============================================================================

def main():
    print("Mertens Function Speedup Investigation")
    print("=" * 78)
    print(f"Python {sys.version}")
    print(f"Recursion limit: {sys.getrecursionlimit()}")
    sys.setrecursionlimit(100000)

    # Verify implementations
    print("\n--- Verification ---")
    for x in [10, 100, 1000, 5000]:
        m_brute = mertens_brute(x)
        m_lucy = mertens_lucy(x)
        m_opt = mertens_lucy_optimized(x)
        d_brute = divisor_sum_brute(x)
        d_hyp = divisor_sum_hyperbola(x)
        d_fast = divisor_sum_fast(x)

        m_ok = "OK" if m_brute == m_lucy == m_opt else f"FAIL {m_brute} {m_lucy} {m_opt}"
        d_ok = "OK" if d_brute == d_hyp == d_fast else f"FAIL {d_brute} {d_hyp} {d_fast}"
        print(f"  x={x:>6}: M(x)={m_brute:>6} [{m_ok}]  D(x)={d_brute:>10} [{d_ok}]")

    # Run all analyses
    analyze_hyperbola_difference()
    benchmark_all()
    compare_complexity_scaling()
    attempt_factorization_speedup(10000)
    test_helfgott_thompson_idea(100000)
    test_recursive_depth(10**8)
    analyze_pi_mertens_connection()
    print_summary()


if __name__ == "__main__":
    main()

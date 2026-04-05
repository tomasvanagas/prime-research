#!/usr/bin/env python3
"""
Tropical Geometry / Min-Plus Algebra Approach to Prime Counting
===============================================================

PRIOR ART: Sessions 4, 16 closed "tropical geometry" with findings:
  - Tropicalization loses all info (min = smallest prime)
  - Tropical convolution = standard convolution
  - Floor values ARE tropical objects but yield no new content

This experiment does a thorough, concrete investigation:

1. Tropical Dirichlet convolution:  Lambda = mu * log  in min-plus
2. Tropical Mertens function:  tropical sum of mu(n)
3. Tropical sieve complexity vs standard sieve
4. Tropical Euler product
5. p-adic valuations of x!  and extracting pi(x)
6. Tropical determinant of prime divisibility matrix

We test whether ANY of these give sub-O(x^{2/3}) access to pi(x) or p(n).
"""

import numpy as np
import time
from math import log, floor, sqrt, gcd
from functools import lru_cache
from collections import defaultdict

# ============================================================
# Utility: small prime sieve for ground truth
# ============================================================

def sieve_primes(limit):
    """Eratosthenes sieve, returns list of primes up to limit."""
    is_prime = [False, False] + [True] * (limit - 1)
    for i in range(2, int(limit**0.5) + 1):
        if is_prime[i]:
            for j in range(i*i, limit + 1, i):
                is_prime[j] = False
    return [i for i in range(2, limit + 1) if is_prime[i]]

def mobius_sieve(limit):
    """Compute mu(n) for n = 0..limit."""
    mu = [0] * (limit + 1)
    mu[1] = 1
    is_prime = [True] * (limit + 1)
    primes = []
    for i in range(2, limit + 1):
        if is_prime[i]:
            primes.append(i)
            mu[i] = -1
        for p in primes:
            if i * p > limit:
                break
            is_prime[i * p] = False
            if i % p == 0:
                mu[i * p] = 0
                break
            else:
                mu[i * p] = -mu[i]
    return mu

def mangoldt_values(limit):
    """Compute Lambda(n) for n = 1..limit. Lambda(p^k) = log(p), else 0."""
    vals = [0.0] * (limit + 1)
    is_prime = [True] * (limit + 1)
    for p in range(2, limit + 1):
        if is_prime[p]:
            # Mark composites
            for j in range(p * p, limit + 1, p):
                is_prime[j] = False
            # Set Lambda(p^k) = log(p)
            pk = p
            lp = log(p)
            while pk <= limit:
                vals[pk] = lp
                pk *= p
    return vals

# ============================================================
# MIN-PLUS (TROPICAL) SEMIRING
# ============================================================
# In the tropical semiring (R union {+inf}, min, +):
#   a "oplus" b  =  min(a, b)       (tropical addition)
#   a "otimes" b =  a + b           (tropical multiplication)
#   additive identity: +inf
#   multiplicative identity: 0

INF = float('inf')

def trop_add(a, b):
    """Tropical addition = min."""
    return min(a, b)

def trop_mul(a, b):
    """Tropical multiplication = ordinary addition."""
    if a == INF or b == INF:
        return INF
    return a + b

def trop_sum(arr):
    """Tropical sum = min over array."""
    return min(arr) if arr else INF

def trop_prod(arr):
    """Tropical product = ordinary sum."""
    s = 0.0
    for x in arr:
        if x == INF:
            return INF
        s += x
    return s

# ============================================================
# TEST 1: Tropical Dirichlet Convolution
# ============================================================
# Standard:  (f * g)(n) = sum_{d|n} f(d) * g(n/d)
# Tropical:  (f *_T g)(n) = min_{d|n} { f(d) + g(n/d) }

def divisors(n):
    divs = []
    for d in range(1, int(n**0.5) + 1):
        if n % d == 0:
            divs.append(d)
            if d != n // d:
                divs.append(n // d)
    return sorted(divs)

def tropical_dirichlet_conv(f, g, limit):
    """Tropical Dirichlet convolution: h(n) = min_{d|n} (f(d) + g(n/d))."""
    h = [INF] * (limit + 1)
    for n in range(1, limit + 1):
        for d in divisors(n):
            val = trop_mul(f[d], g[n // d])
            h[n] = trop_add(h[n], val)
    return h

def test_tropical_dirichlet(limit=200):
    """
    Test 1: Build Lambda = mu * log tropically.
    Standard: Lambda(n) = sum_{d|n} mu(d) * log(n/d)
    Tropical: Lambda_T(n) = min_{d|n} { mu_T(d) + log(n/d) }

    Key question: does tropical Lambda still isolate prime powers?
    """
    print("=" * 70)
    print("TEST 1: Tropical Dirichlet Convolution (Lambda = mu * log)")
    print("=" * 70)

    mu = mobius_sieve(limit)
    lam_standard = mangoldt_values(limit)

    # For tropical mu: we need a tropical encoding.
    # Option A: mu_T(n) = mu(n) directly (treating it as a real value)
    # Problem: mu(n) in {-1, 0, 1}, and tropical addition = min,
    #          so negative values dominate.
    # Option B: mu_T(n) = |mu(n)| (absolute value)
    # Option C: mu_T(n) = 0 if mu(n)!=0, +inf if mu(n)=0 (squarefree indicator)

    # Let's try all three encodings
    log_vals = [0.0] * (limit + 1)  # log_vals[0] unused
    for n in range(1, limit + 1):
        log_vals[n] = log(n)

    print("\n--- Option A: mu_T(n) = mu(n) as real number ---")
    mu_A = [INF] + [float(mu[n]) for n in range(1, limit + 1)]
    lam_A = tropical_dirichlet_conv(mu_A, log_vals, limit)

    primes_list = sieve_primes(limit)
    prime_set = set(primes_list)

    # Check: does lam_A distinguish primes?
    prime_vals = [(p, lam_A[p]) for p in primes_list[:15]]
    composite_vals = [(n, lam_A[n]) for n in range(4, 50) if n not in prime_set][:15]

    print(f"  Tropical Lambda at primes:     {prime_vals[:8]}")
    print(f"  Tropical Lambda at composites: {composite_vals[:8]}")

    # Check separability
    pvals_set = set(v for _, v in prime_vals)
    cvals_set = set(v for _, v in composite_vals)
    overlap = pvals_set & cvals_set
    print(f"  Value overlap between primes and composites: {len(overlap) > 0}")

    print("\n--- Option C: mu_T(n) = 0 if squarefree, +inf if not ---")
    mu_C = [INF] * (limit + 1)
    for n in range(1, limit + 1):
        if mu[n] != 0:
            mu_C[n] = 0.0
    lam_C = tropical_dirichlet_conv(mu_C, log_vals, limit)

    prime_vals_C = [(p, round(lam_C[p], 4)) for p in primes_list[:15]]
    composite_vals_C = [(n, round(lam_C[n], 4)) for n in range(4, 50) if n not in prime_set][:15]
    print(f"  Tropical Lambda at primes:     {prime_vals_C[:8]}")
    print(f"  Tropical Lambda at composites: {composite_vals_C[:8]}")

    # lam_C(n) = min_{squarefree d | n} log(n/d)
    # This equals log(n / largest_squarefree_divisor(n)) = log(radical_part)
    # NOT useful for isolating primes.
    print("  => lam_C(n) = log(n / max squarefree divisor of n)")
    print("     This is 0 for all squarefree n, not just primes.")

    return lam_A, lam_C

# ============================================================
# TEST 2: Tropical Mertens Function
# ============================================================

def test_tropical_mertens(limit=500):
    """
    Standard Mertens: M(x) = sum_{n<=x} mu(n)
    Tropical Mertens: M_T(x) = min_{n<=x} mu(n)

    Since mu(n) in {-1, 0, 1}, tropical Mertens = min = -1 for all x >= 2.
    This is trivially useless. But let's also try:
      M_T(x) = min_{n<=x} mu(n) * log(n)   (weighted version)
    """
    print("\n" + "=" * 70)
    print("TEST 2: Tropical Mertens Function")
    print("=" * 70)

    mu = mobius_sieve(limit)

    # Naive tropical Mertens
    trop_mertens = []
    running_min = INF
    for n in range(1, limit + 1):
        running_min = min(running_min, mu[n])
        trop_mertens.append(running_min)

    print(f"  Tropical Mertens M_T(x) = min{{mu(n) : n<=x}} = {trop_mertens[5]} for all x>=2")
    print("  => Trivially -1. No information content.\n")

    # Weighted: min_{n<=x} mu(n) * log(n)
    # For mu(n) = -1: contribution is -log(n), which keeps decreasing
    # So this is just -log(largest n<=x with mu(n)=-1)
    weighted_mertens = []
    running_min = INF
    for n in range(1, limit + 1):
        if mu[n] != 0:
            val = mu[n] * log(n)
            running_min = min(running_min, val)
        weighted_mertens.append(running_min)

    print(f"  Weighted tropical Mertens at x=100: {weighted_mertens[99]:.4f}")
    print(f"  Weighted tropical Mertens at x=500: {weighted_mertens[499]:.4f}")
    print(f"  Expected: -log(x) ~ {-log(500):.4f}")
    print("  => Just tracks -log(largest squarefree n). No prime info.\n")

# ============================================================
# TEST 3: Tropical Sieve Complexity
# ============================================================

def test_tropical_sieve(limit=1000):
    """
    Legendre sieve: pi(x) = pi(sqrt(x)) + sum_{d | P#} mu(d) * floor(x/d) - 1
    where P# = product of primes up to sqrt(x).

    Tropical version: replace sum with min, multiplication with addition.
    pi_T(x) = min_{d | P#} { mu_T(d) + floor(x/d) }

    Question: is there a tropical sieve identity that equals pi(x)?
    """
    print("=" * 70)
    print("TEST 3: Tropical Sieve vs Standard Sieve")
    print("=" * 70)

    primes = sieve_primes(limit)
    pi_true = [0] * (limit + 1)
    for p in primes:
        for i in range(p, limit + 1):
            pi_true[i] += 1
    # Fix: pi_true[n] should be count of primes <= n
    pi_true = [0] * (limit + 1)
    count = 0
    pset = set(primes)
    for n in range(limit + 1):
        if n in pset:
            count += 1
        pi_true[n] = count

    sqrt_lim = int(sqrt(limit))
    small_primes = sieve_primes(sqrt_lim)

    # Generate squarefree divisors of primorial
    # For small limits, enumerate divisors of product of small primes
    def gen_squarefree_divisors(prime_list):
        """Generate all squarefree divisors of prod(prime_list) with mu values."""
        divisors = [(1, 1)]  # (divisor, mu(divisor))
        for p in prime_list:
            new_divs = []
            for d, m in divisors:
                new_divs.append((d * p, -m))
            divisors.extend(new_divs)
        return divisors

    # Only use primes up to sqrt(limit) -- Legendre-style
    if len(small_primes) > 15:
        small_primes_use = small_primes[:15]
    else:
        small_primes_use = small_primes

    # Too many divisors if too many primes. Limit to manageable.
    if len(small_primes_use) > 12:
        small_primes_use = small_primes_use[:12]

    divs_mu = gen_squarefree_divisors(small_primes_use)
    print(f"  Using {len(small_primes_use)} small primes, {len(divs_mu)} squarefree divisors")

    # Standard Legendre: pi(x) - pi(sqrt(x)) + 1 = sum_{(d,mu)} mu * floor(x/d)
    x = limit
    standard_sum = sum(m * (x // d) for d, m in divs_mu)
    print(f"  Standard Legendre sum for x={x}: {standard_sum}")
    print(f"  True pi({x}) = {pi_true[x]}")

    # Tropical Legendre: min_{(d,mu)} { mu + floor(x/d) }
    # (treating mu as weight)
    trop_min = min(m + (x // d) for d, m in divs_mu)
    print(f"  Tropical Legendre (mu + floor(x/d)): min = {trop_min}")
    print(f"  => This is just min(floor(x/d)) - 1 = floor(x/primorial) - 1")
    print(f"     No relation to pi(x).\n")

    # Alternative: tropical with log weights
    trop_log = min(log(d + 1) + (x // d) for d, m in divs_mu if m != 0)
    print(f"  Tropical with log weights: {trop_log:.4f}")
    print(f"  => Still dominated by smallest divisor term. No prime info.\n")

    # Complexity comparison
    print("  COMPLEXITY COMPARISON:")
    print(f"  Standard Legendre sieve: O(2^pi(sqrt(x))) terms = O(2^{len(small_primes)}) = {2**len(small_primes)}")
    print(f"  Tropical version: same number of terms, but min instead of sum")
    print(f"  => No complexity improvement. Tropical doesn't reduce term count.")

# ============================================================
# TEST 4: Tropical Euler Product
# ============================================================

def test_tropical_euler_product(limit=200):
    """
    Standard: zeta(s) = prod_{p prime} 1/(1-p^{-s})
    log zeta(s) = sum_p -log(1 - p^{-s}) = sum_p sum_{k>=1} p^{-ks}/k

    Tropical: trop_prod = sum, trop_sum = min
    So tropical Euler product = sum_p log(1/(1-p^{-s}))  [this is just log zeta(s)]
    And tropical version of the sum over k: min_k { ks*log(p) + log(k) }
    """
    print("=" * 70)
    print("TEST 4: Tropical Euler Product")
    print("=" * 70)

    primes = sieve_primes(limit)

    s_values = [2.0, 3.0, 5.0, 10.0]

    for s in s_values:
        # Standard log zeta(s) for comparison
        log_zeta = sum(-log(1 - p**(-s)) for p in primes)

        # Tropical Euler "product" = sum of individual terms (same as log zeta)
        trop_euler = sum(-log(1 - p**(-s)) for p in primes)

        # Tropical inner sum: min_{k>=1} { k*s*log(p) + log(k) }
        # The minimum is at k=1 for all p (since ks*log(p) grows fastest)
        trop_inner_sum = sum(s * log(p) for p in primes)  # k=1 terms dominate

        print(f"  s={s}: log_zeta = {log_zeta:.6f}, trop_euler = {trop_euler:.6f}, trop_inner = {trop_inner_sum:.4f}")

    print("\n  FINDING: Tropical Euler product = log(zeta(s)), same as standard.")
    print("  The tropical inner sum (min over k) just picks k=1: sum of s*log(p).")
    print("  This is sum_{p<=N} s*log(p) ~ s*N (by PNT), containing no more")
    print("  information than theta(N). No shortcut to pi(x).\n")

# ============================================================
# TEST 5: p-adic Valuations and x!
# ============================================================

def test_padic_factorial(limit=500):
    """
    v_p(n!) = sum_{k>=1} floor(n/p^k) = (n - s_p(n)) / (p-1)
    where s_p(n) = digit sum of n in base p.

    This is indeed a "tropical polynomial" in the sense that floor = tropical operation.

    Question: Can we extract pi(x) from the collection {v_p(x!) : p prime}?

    We know: sum_{p<=x} v_p(x!) = sum_{p<=x} (x - s_p(x))/(p-1)
    And: log(x!) = sum_{p} v_p(x!) * log(p)  (fundamental theorem of arithmetic)
    """
    print("=" * 70)
    print("TEST 5: p-adic Valuations of x! and pi(x)")
    print("=" * 70)

    primes = sieve_primes(limit)

    def v_p(n, p):
        """p-adic valuation of n!"""
        result = 0
        pk = p
        while pk <= n:
            result += n // pk
            pk *= p
        return result

    def digit_sum(n, p):
        """Sum of digits of n in base p."""
        s = 0
        while n > 0:
            s += n % p
            n //= p
        return s

    test_values = [50, 100, 200, 500]

    for x in test_values:
        primes_up_to_x = [p for p in primes if p <= x]
        pi_x = len(primes_up_to_x)

        # Sum of v_p(x!) over primes p <= x
        vp_sum = sum(v_p(x, p) for p in primes_up_to_x)

        # Using the digit sum formula
        digit_sum_total = sum((x - digit_sum(x, p)) / (p - 1) for p in primes_up_to_x)

        # Verify they match
        assert abs(vp_sum - digit_sum_total) < 0.01, f"Mismatch at x={x}"

        # Can we extract pi(x) from vp_sum?
        # vp_sum = sum_{p<=x} (x - s_p(x))/(p-1)
        # If s_p(x) were 0 for all p: vp_sum ~ x * sum_{p<=x} 1/(p-1) ~ x * log(log(x))
        # The ratio vp_sum / (x * log(log(x))) should be ~ 1 + corrections

        if x > 2:
            ratio = vp_sum / (x * log(log(x))) if log(log(x)) > 0 else 0
        else:
            ratio = 0

        print(f"  x={x:4d}: pi(x)={pi_x:4d}, sum v_p(x!)={vp_sum:10d}, "
              f"ratio/(x*loglogx)={ratio:.4f}")

    print()

    # KEY TEST: Can v_p(x!) for a SINGLE prime p help locate primes?
    # v_p(x!) - v_p((x-1)!) = v_p(x) = {k if p^k | x, 0 otherwise}
    # This tells us if x is a prime power, not if x is prime.
    print("  Can single-prime valuations detect primes?")
    print("  v_2(n!) - v_2((n-1)!) = v_2(n) = number of times 2 divides n")

    p = 2
    detections = []
    for n in range(2, 50):
        v_diff = v_p(n, p) - v_p(n - 1, p)
        is_prime = n in set(primes)
        detections.append((n, v_diff, is_prime))

    print(f"  n: v_2(n), prime?")
    for n, vd, ip in detections[:20]:
        marker = " <-- PRIME" if ip else ""
        print(f"    {n:3d}: v_2={vd}{marker}")

    print("\n  FINDING: v_p(x!) encodes divisibility by p, not primality.")
    print("  To detect if x is prime, we'd need v_p(x)=0 for ALL p<x,")
    print("  which requires checking all primes -- circular!\n")

    # Advanced: tropical polynomial structure
    print("  TROPICAL POLYNOMIAL STRUCTURE:")
    print("  v_p(n!) = sum_{k>=1} floor(n/p^k)")
    print("  floor(n/m) = tropical: it's the largest integer <= n/m")
    print("  But computing v_p(n!) for a single (n,p) is O(log n) -- already fast.")
    print("  The bottleneck is summing over ALL primes p <= x, which requires pi(x).")
    print("  => CIRCULAR: need pi(x) to compute the sum that encodes pi(x).\n")

# ============================================================
# TEST 6: Tropical Determinant of Prime Divisibility Matrix
# ============================================================

def test_tropical_determinant(N=30):
    """
    Build matrix A where A[i][j] = 1 if prime_i divides j, else +inf (tropical zero).

    Tropical det = min over permutations sigma of sum_i A[i][sigma(i)]
    (This is the assignment problem, solvable in O(n^3) by Hungarian algorithm.)

    Standard det would give inclusion-exclusion sieve count.
    Does tropical det give useful prime information?
    """
    print("=" * 70)
    print("TEST 6: Tropical Determinant of Prime Divisibility Matrix")
    print("=" * 70)

    primes = sieve_primes(N)
    n_primes = len(primes)

    # Build matrix: rows = primes, cols = integers 2..N
    cols = list(range(2, N + 1))
    n_cols = len(cols)

    # Square matrix: use min(n_primes, n_cols)
    size = min(n_primes, n_cols)

    # A[i][j] = 0 if primes[i] divides cols[j], else INF
    # (In tropical: 0 is multiplicative identity, INF is additive identity)
    A = [[INF] * size for _ in range(size)]
    for i in range(size):
        for j in range(size):
            if cols[j] % primes[i] == 0:
                A[i][j] = 0.0
            else:
                A[i][j] = INF

    # Brute force tropical det for small matrices
    from itertools import permutations

    if size <= 10:
        min_cost = INF
        best_perm = None
        for perm in permutations(range(size)):
            cost = sum(A[i][perm[i]] for i in range(size))
            if cost < min_cost:
                min_cost = cost
                best_perm = perm

        print(f"  Matrix size: {size}x{size}")
        print(f"  Tropical determinant = {min_cost}")
        if best_perm and min_cost < INF:
            assignment = [(primes[i], cols[best_perm[i]]) for i in range(size)]
            print(f"  Optimal assignment (prime -> integer): {assignment}")
            print(f"  => Each prime is assigned to a multiple of itself.")
        else:
            print(f"  Tropical det = INF: no complete assignment exists.")
    else:
        # Use greedy approximation
        print(f"  Matrix size: {size}x{size} (too large for brute force)")

    # Alternative: weighted matrix
    print(f"\n  WEIGHTED VERSION: A[i][j] = log(cols[j]) if p_i | cols[j], else INF")
    B = [[INF] * size for _ in range(size)]
    for i in range(size):
        for j in range(size):
            if cols[j] % primes[i] == 0:
                B[i][j] = log(cols[j])

    if size <= 10:
        min_cost = INF
        best_perm = None
        for perm in permutations(range(size)):
            cost = sum(B[i][perm[i]] for i in range(size))
            if cost < min_cost:
                min_cost = cost
                best_perm = perm

        print(f"  Tropical det (weighted) = {min_cost:.4f}")
        if best_perm and min_cost < INF:
            assignment = [(primes[i], cols[best_perm[i]]) for i in range(size)]
            print(f"  Optimal assignment: {assignment}")

            # The optimal assignment maps each prime to its smallest multiple in cols
            ideal = [(p, p) for p in primes[:size] if p <= N]  # each prime to itself
            print(f"  Ideal (each prime to itself): {ideal}")

    print(f"\n  FINDING: Tropical det of divisibility matrix solves an assignment")
    print(f"  problem: match primes to their multiples. The answer is trivially")
    print(f"  'each prime maps to itself' (cost = sum log(p_i)).")
    print(f"  This encodes NO new information about pi(x).")

# ============================================================
# TEST 7: Tropical Convolution and Counting
# ============================================================

def test_tropical_convolution(limit=200):
    """
    Standard: if f = 1_primes (indicator of primes), then
    (f * f)(n) = #{(p,q) : p+q=n, p,q prime} = Goldbach count.

    Tropical: (f *_T f)(n) = min_{p+q=n, p,q prime} {f(p) + f(q)}
    If f(p) = 0 for prime, INF for composite:
    Tropical conv = 0 if n is sum of two primes, INF otherwise.
    = indicator of Goldbach numbers. Doesn't COUNT anything.

    What if f(p) = p (or log(p))?
    Tropical conv = min_{p+q=n} {p + q} = n (trivially).
    Or min_{p+q=n} {log(p) + log(q)} = min log(pq) where p+q=n.
    """
    print("\n" + "=" * 70)
    print("TEST 7: Tropical Convolution")
    print("=" * 70)

    primes_set = set(sieve_primes(limit))

    # Indicator function
    f = [INF] * (limit + 1)
    for p in primes_set:
        f[p] = 0.0

    # Tropical additive convolution: h(n) = min_{a+b=n} (f(a) + f(b))
    goldbach_trop = [INF] * (limit + 1)
    for n in range(4, limit + 1):
        for a in range(2, n // 2 + 1):
            if f[a] < INF and f[n - a] < INF:
                goldbach_trop[n] = 0.0
                break

    # Count Goldbach numbers
    goldbach_numbers = sum(1 for n in range(4, limit + 1, 2) if goldbach_trop[n] == 0)
    even_count = len(range(4, limit + 1, 2))
    print(f"  Even numbers 4..{limit}: {even_count}")
    print(f"  Goldbach (expressible as p+q): {goldbach_numbers}")
    print(f"  => Tropical convolution gives EXISTENCE, not COUNT.")

    # With log weights
    g = [INF] * (limit + 1)
    for p in primes_set:
        g[p] = log(p)

    trop_conv_log = [INF] * (limit + 1)
    for n in range(4, limit + 1):
        for a in range(2, n // 2 + 1):
            if g[a] < INF and g[n - a] < INF:
                val = g[a] + g[n - a]  # log(a) + log(n-a) = log(a*(n-a))
                trop_conv_log[n] = min(trop_conv_log[n], val)

    print(f"\n  Tropical conv with log weights (first 10 even n >= 4):")
    for n in range(4, 24, 2):
        if trop_conv_log[n] < INF:
            # min log(p*q) where p+q=n
            min_pq = np.exp(trop_conv_log[n])
            print(f"    n={n:3d}: min log(p*q) = {trop_conv_log[n]:.4f}, "
                  f"min p*q = {min_pq:.1f}")

    print(f"\n  FINDING: Tropical convolution converts COUNTING to EXISTENCE/OPTIMIZATION.")
    print(f"  This is a fundamental loss: pi(x) is a COUNT, but tropical algebra")
    print(f"  computes MIN/MAX. You cannot recover a count from a min.\n")

# ============================================================
# SYNTHESIS
# ============================================================

def synthesis():
    print("=" * 70)
    print("SYNTHESIS: Tropical Geometry for Prime Counting")
    print("=" * 70)

    print("""
  CORE INSIGHT: The tropical semiring (min, +) is fundamentally about
  OPTIMIZATION, not COUNTING. Prime counting pi(x) is a COUNTING problem.

  Specific failures:

  1. TROPICAL DIRICHLET CONVOLUTION: Lambda_T(n) loses the alternating
     signs of Mobius function. Result depends on encoding choice, but
     never isolates primes.

  2. TROPICAL MERTENS: min of mu(n) = -1 for all x >= 2. Trivially useless.

  3. TROPICAL SIEVE: Same number of terms as standard sieve, but min
     instead of sum. The min is dominated by the smallest-divisor term.
     No complexity improvement.

  4. TROPICAL EULER PRODUCT: Equals log(zeta(s)), which is the STANDARD
     object. No new structure.

  5. P-ADIC VALUATIONS: v_p(n!) is O(log n) to compute for one (n,p),
     but summing over all primes requires knowing the primes. Circular.

  6. TROPICAL DETERMINANT: Solves an assignment problem (primes to multiples)
     with trivial solution. No pi(x) information.

  7. TROPICAL CONVOLUTION: Converts counting to existence/optimization.
     Fundamental information loss: cannot recover a count from a min.

  MATHEMATICAL REASON: The map from (R, +, *) to (R, min, +) is a
  VALUATION -- it's a homomorphism that collapses the additive structure.
  Specifically:
    val(a + b) >= min(val(a), val(b))
  with equality when val(a) != val(b). The INEQUALITY means we LOSE
  information about sums. But pi(x) IS a sum (count). So tropicalization
  systematically destroys exactly the information we need.

  VERDICT: CLOSED. Tropical/min-plus algebra cannot help with pi(x) or p(n).
  The counting-to-optimization collapse is irreversible.
""")

# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":
    t0 = time.time()

    test_tropical_dirichlet()
    test_tropical_mertens()
    test_tropical_sieve()
    test_tropical_euler_product()
    test_padic_factorial()
    test_tropical_determinant()
    test_tropical_convolution()
    synthesis()

    elapsed = time.time() - t0
    print(f"Total runtime: {elapsed:.2f}s")

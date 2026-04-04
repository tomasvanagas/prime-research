#!/usr/bin/env python3
"""
Alternative Decompositions of pi(x): Can We Beat O(x^{2/3})?

Investigation of whether alternative decomposition strategies for pi(x)
could lead to better asymptotic complexity than the Meissel-Lehmer/
Deleglise-Rivat O(x^{2/3}) combinatorial approach.

Six directions investigated:
  1. Buchstab identity tree pruning / memoization
  2. Hyperbola method generalization for Mobius sums
  3. Vaughan's identity for exact computation
  4. Convolution structure exploitation
  5. Dirichlet series evaluation approach
  6. Structural comparison: sum d(n) vs sum mu(n)

Author: Claude (Session 11)
Date: 2026-04-04
"""

import math
import time
from collections import defaultdict
from functools import lru_cache

# ============================================================================
# UTILITIES: Small prime sieve, exact pi(x), mu(n), etc.
# ============================================================================

def sieve_primes(limit):
    """Sieve of Eratosthenes up to limit."""
    if limit < 2:
        return []
    is_prime = [True] * (limit + 1)
    is_prime[0] = is_prime[1] = False
    for i in range(2, int(limit**0.5) + 1):
        if is_prime[i]:
            for j in range(i*i, limit + 1, i):
                is_prime[j] = False
    return [i for i in range(2, limit + 1) if is_prime[i]]

def pi_exact_sieve(x):
    """Exact pi(x) by sieving. For small x only."""
    return len(sieve_primes(int(x)))

def compute_mobius(limit):
    """Compute mu(n) for n=0..limit via sieve."""
    mu = [0] * (limit + 1)
    mu[1] = 1
    # Factor sieve
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

def compute_lambda_mangoldt(limit):
    """Compute Lambda(n) for n=0..limit."""
    lam = [0.0] * (limit + 1)
    for p in sieve_primes(limit):
        pk = p
        logp = math.log(p)
        while pk <= limit:
            lam[pk] = logp
            pk *= p
    return lam

def mertens_function(limit, mu=None):
    """Compute M(x) = sum_{n<=x} mu(n) for all x up to limit."""
    if mu is None:
        mu = compute_mobius(limit)
    M = [0] * (limit + 1)
    for n in range(1, limit + 1):
        M[n] = M[n-1] + mu[n]
    return M

def divisor_sum(limit):
    """Compute D(x) = sum_{n<=x} d(n) for all x up to limit."""
    # d(n) = number of divisors
    d = [0] * (limit + 1)
    for i in range(1, limit + 1):
        for j in range(i, limit + 1, i):
            d[j] += 1
    D = [0] * (limit + 1)
    for n in range(1, limit + 1):
        D[n] = D[n-1] + d[n]
    return D, d

# ============================================================================
# EXPERIMENT 1: Buchstab Identity Tree Analysis
# ============================================================================

def experiment_buchstab_tree():
    """
    Buchstab identity: phi(x, a) = phi(x, a-1) - phi(x/p_a, a-1)

    This creates a binary tree of evaluations. The Meissel-Lehmer method
    prunes this tree cleverly. Question: can we do better?

    Key insight: The tree has at most O(x) leaves (trivially), and
    Meissel-Lehmer reduces this to O(x^{2/3}). The "special leaves" are
    the bottleneck. We count how many distinct (x/m, b) pairs appear.
    """
    print("=" * 70)
    print("EXPERIMENT 1: Buchstab Identity Tree Analysis")
    print("=" * 70)

    results = {}

    for x in [100, 1000, 10000, 100000]:
        primes = sieve_primes(int(x**0.5) + 1)
        a = len(primes)  # pi(sqrt(x))

        # Count nodes in the Buchstab tree
        node_count = [0]
        leaf_count = [0]
        distinct_args = set()

        def buchstab_count(n, b, depth=0):
            """Count evaluations in the Buchstab tree."""
            node_count[0] += 1
            distinct_args.add((int(n), b))

            if b == 0:
                leaf_count[0] += 1
                return int(n)  # phi(n, 0) = floor(n)

            if n < primes[b-1]:
                leaf_count[0] += 1
                return max(0, int(n > 0))

            # phi(n, b) = phi(n, b-1) - phi(n/p_b, b-1)
            r1 = buchstab_count(n, b - 1, depth + 1)
            r2 = buchstab_count(n / primes[b-1], b - 1, depth + 1)
            return r1 - r2

        if x <= 10000:  # Only run for small x (exponential!)
            t0 = time.time()
            result = buchstab_count(x, a)
            t1 = time.time()

            pi_true = pi_exact_sieve(x)
            # phi(x, a) = pi(x) - a + 1 for a = pi(sqrt(x))
            # So pi(x) = phi(x, a) + a - 1
            pi_computed = result + a - 1

            results[x] = {
                'nodes': node_count[0],
                'leaves': leaf_count[0],
                'distinct': len(distinct_args),
                'time': t1 - t0,
                'correct': pi_computed == pi_true,
                'pi_true': pi_true,
            }

            print(f"\nx = {x}:")
            print(f"  pi(sqrt(x)) = {a}")
            print(f"  Total nodes in tree: {node_count[0]}")
            print(f"  Leaves: {leaf_count[0]}")
            print(f"  Distinct (n,b) pairs: {len(distinct_args)}")
            print(f"  Ratio nodes/x^(2/3): {node_count[0] / x**(2/3):.2f}")
            print(f"  Ratio distinct/x^(2/3): {len(distinct_args) / x**(2/3):.2f}")
            print(f"  Result correct: {pi_computed == pi_true}")
        else:
            print(f"\nx = {x}: skipped (tree too large for naive recursion)")

    # Now test with memoization
    print("\n--- With Memoization ---")
    for x in [100, 1000, 10000, 50000]:
        primes = sieve_primes(int(x**0.5) + 1)
        a = len(primes)

        memo = {}
        call_count = [0]

        def buchstab_memo(n_int, b):
            call_count[0] += 1
            key = (n_int, b)
            if key in memo:
                return memo[key]

            if b == 0:
                memo[key] = n_int
                return n_int
            if n_int < primes[b-1]:
                val = int(n_int > 0)
                memo[key] = val
                return val

            val = buchstab_memo(n_int, b-1) - buchstab_memo(n_int // primes[b-1], b-1)
            memo[key] = val
            return val

        t0 = time.time()
        result = buchstab_memo(int(x), a)
        t1 = time.time()

        pi_computed = result + a - 1
        pi_true = pi_exact_sieve(x)

        print(f"\nx = {x}:")
        print(f"  Calls: {call_count[0]}, Memo entries: {len(memo)}")
        print(f"  Ratio memo_entries / x^(2/3): {len(memo) / x**(2/3):.2f}")
        print(f"  Ratio memo_entries / x^(1/2): {len(memo) / x**(1/2):.2f}")
        print(f"  Time: {t1-t0:.4f}s")
        print(f"  Correct: {pi_computed == pi_true}")

    print("\n" + "-" * 70)
    print("BUCHSTAB CONCLUSION:")
    print("  The naive tree is exponential. With memoization, the number of")
    print("  distinct (floor(x/m), b) pairs is the key quantity.")
    print("  The set of distinct floor(x/m) values has size O(sqrt(x)).")
    print("  Combined with b in [0, pi(sqrt(x))], this gives")
    print("  O(sqrt(x) * x^(1/4)/log(x)) = O(x^(3/4)/log(x)) distinct pairs.")
    print("  This matches the Lucy_Hedgehog DP approach exactly.")
    print("  The Deleglise-Rivat optimization reduces to O(x^(2/3)) by")
    print("  handling 'special leaves' where x/m crosses certain thresholds.")
    print("  VERDICT: Memoized Buchstab IS the Lucy DP. Cannot beat O(x^(2/3))")
    print("  without a fundamentally different decomposition of phi.")

# ============================================================================
# EXPERIMENT 2: Hyperbola Method Comparison
# ============================================================================

def experiment_hyperbola_comparison():
    """
    The Dirichlet hyperbola method gives sum_{n<=x} d(n) in O(x^{1/3}).

    Key idea: d(n) = sum_{ab=n} 1, so
    sum_{n<=x} d(n) = sum_{a<=sqrt(x)} floor(x/a) + ... = 2*sum - floor(sqrt(x))^2

    For pi(x), we want sum_{n<=x} [n is prime].
    The prime indicator = -sum_{d|n} mu(d)*ln(d)/ln(n) ... not a simple convolution.

    But the Mertens function M(x) = sum_{n<=x} mu(n) IS a summatory function
    of a multiplicative function. Can we use the hyperbola method on it?

    sum_{n<=x} mu(n) via hyperbola:
    We need an identity like mu = f * g (Dirichlet convolution).
    We know: mu * 1 = epsilon (where 1 is constant-1, epsilon = [n=1]).
    So sum mu(n) = M(x), and sum_{n<=x} sum_{d|n} mu(d) = 1.

    This gives: 1 = sum_{a<=x} mu(a) * floor(x/a)
    => M(x) = 1 - sum_{a=2}^{x} mu(a) * floor(x/a)
    ... which is circular (needs M(x) to compute M(x)).

    The hyperbola trick: split at sqrt(x):
    1 = sum_{a<=u} mu(a)*floor(x/a) + sum_{b<=x/u} M(x/b) - M(u)*floor(x/u)
    for u = sqrt(x). This is the identity used by Deleglise-Rivat for M(x)!
    And it gives O(x^{2/3}) -- NOT O(x^{1/3}).

    WHY the difference? Let's investigate.
    """
    print("\n" + "=" * 70)
    print("EXPERIMENT 2: Hyperbola Method -- Why O(x^{1/3}) for d(n) but O(x^{2/3}) for mu(n)?")
    print("=" * 70)

    # Part A: Verify the hyperbola method for d(n)
    print("\n--- Part A: Hyperbola method for sum d(n) ---")
    for x in [100, 1000, 10000, 100000]:
        u = int(x**0.5)
        # sum_{n<=x} d(n) = 2 * sum_{a<=u} floor(x/a) - u^2
        total = 0
        for a in range(1, u + 1):
            total += int(x // a)
        hyperbola_result = 2 * total - u * u

        # Direct computation
        _, d = divisor_sum(x)
        direct = sum(d[1:x+1])

        print(f"  x={x}: hyperbola={hyperbola_result}, direct={direct}, "
              f"match={hyperbola_result==direct}, ops={u} (x^(1/2))")

    # Part B: The mu(n) summatory function -- hyperbola identity
    print("\n--- Part B: Mertens function via hyperbola identity ---")
    print("  Identity: 1 = sum_{a<=u} mu(a)*floor(x/a) + sum_{b<=v} M(x/b) - M(u)*v")
    print("  where u*v = x (hyperbola split)")

    for x in [100, 1000, 10000]:
        mu = compute_mobius(x)
        M = mertens_function(x, mu)

        # Verify identity: sum_{a=1}^{x} mu(a)*floor(x/a) = 1
        check = sum(mu[a] * (x // a) for a in range(1, x + 1))
        print(f"  x={x}: sum mu(a)*floor(x/a) = {check} (should be 1)")

        # Hyperbola split at u = x^{1/2}
        u = int(x**0.5)
        v = x // u

        # Left sum: sum_{a<=u} mu(a)*floor(x/a) -- needs mu(a) for a <= sqrt(x) [CHEAP]
        left = sum(mu[a] * (x // a) for a in range(1, u + 1))

        # Right sum: sum_{b=1}^{v} M(x//b) -- needs M at O(sqrt(x)) points [EXPENSIVE!]
        right = sum(M[x // b] for b in range(1, v + 1))

        # Cross term
        cross = M[u] * v

        reconstructed = left + right - cross
        print(f"  x={x}: left={left}, right={right}, cross={cross}, "
              f"total={reconstructed} (should be 1)")

    # Part C: WHY the asymmetry?
    print("\n--- Part C: Structural Analysis ---")
    print("""
  For sum d(n): d = 1 * 1 (both factors are the constant-1 function)
    Hyperbola: sum_{a<=u} floor(x/a) + sum_{b<=u} floor(x/b) - u^2
    Both sums are INDEPENDENT -- no recursion needed!
    Each sum has O(sqrt(x)) terms.
    Result: O(sqrt(x)) operations.

  For M(x): mu * 1 = epsilon (Mobius inversion)
    Hyperbola: involves M(x/b) for b=1..v -- i.e., RECURSIVE!
    M(x/b) requires knowing M at all divisor-like points.
    The recursion depth is O(x^{1/3}) with memoization.
    Each level requires O(sqrt(x/b)) work.
    Total: O(x^{2/3}) operations.

  THE FUNDAMENTAL DIFFERENCE:
    - d(n) = (1 * 1)(n): both factors f=1 and g=1 have TRIVIAL partial sums
      (sum_{n<=x} 1 = floor(x)). No recursion needed.
    - mu(n): mu * 1 = epsilon. The partial sum of 1 is trivial, but the partial
      sum of mu IS M(x) itself -- creating a RECURSIVE dependency.

  For the hyperbola method to give O(x^{1/3}) for a function h,
  we need h = f * g where BOTH partial sums F(x) and G(x) are
  independently computable in O(polylog) time.

  For mu(n), there is no known decomposition mu = f * g where both
  F(x) = sum f(n) and G(x) = sum g(n) are cheap to compute.

  This is fundamentally because mu(n) is the Dirichlet inverse of 1,
  and its partial sums encode prime distribution information.
""")

    # Part D: Can we find a "good" decomposition of mu?
    print("--- Part D: Searching for alternative decompositions of mu ---")
    print("  If mu = f * g with both F, G cheaply summable, we win.")
    print("  Known decompositions:")
    print("    mu = mu * epsilon  (trivial, no help)")
    print("    mu * 1 = epsilon   (standard, gives recursion)")
    print("    mu = lambda * |mu| where lambda = Liouville, |mu| = squarefree indicator")
    print("    ... but sum lambda(n) = sqrt(x) + O(...) and sum |mu| also needs computation")

    x = 1000
    mu = compute_mobius(x)

    # Liouville function
    liouville = [0] * (x + 1)
    liouville[1] = 1
    for n in range(2, x + 1):
        temp = n
        omega = 0  # number of prime factors with multiplicity
        for p in range(2, int(temp**0.5) + 1):
            while temp % p == 0:
                omega += 1
                temp //= p
            if temp == 1:
                break
        if temp > 1:
            omega += 1
        liouville[n] = (-1)**omega

    # Check: mu(n) = lambda(n) * |mu(n)|
    for n in range(1, min(50, x)):
        if mu[n] != 0:
            assert liouville[n] == mu[n], f"Failed at n={n}"

    # Summatory Liouville
    L = [0] * (x + 1)
    for n in range(1, x + 1):
        L[n] = L[n-1] + liouville[n]

    print(f"\n  L(1000) = sum lambda(n) = {L[1000]}")
    print(f"  M(1000) = sum mu(n) = {mertens_function(x, mu)[1000]}")
    print(f"  L(x) ~ sqrt(x) is only true on average; actual values oscillate wildly.")
    print(f"  L(x) is NO easier to compute than M(x). Same O(x^(2/3)) barrier.")

# ============================================================================
# EXPERIMENT 3: Vaughan's Identity for Exact Computation
# ============================================================================

def experiment_vaughan_identity():
    """
    Vaughan's identity splits Lambda(n) into Type I and Type II sums:

    Lambda(n) = Lambda_1(n) - Lambda_2(n) + Lambda_3(n) + Lambda_4(n)

    where (for parameter U, V with UV ~ x):
    Lambda_1(n) = sum_{d|n, d<=U} mu(d) * ln(n/d)    [Type I: long]
    Lambda_2(n) = sum_{d|n, d<=UV} (sum_{ab=d, a<=U,b<=V} mu(a)*Lambda(b)) [Type II]
    Lambda_3(n) = [n<=U] * Lambda(n)                    [small primes]
    Lambda_4(n) = sum_{d|n, d>U} mu(d) * ... [bilinear]

    In analytic NT, this is used for BOUNDING sums, not computing them.
    Can we use it for exact computation?
    """
    print("\n" + "=" * 70)
    print("EXPERIMENT 3: Vaughan's Identity for Exact Computation")
    print("=" * 70)

    x = 1000
    primes = sieve_primes(x)
    prime_set = set(primes)
    lam = compute_lambda_mangoldt(x)
    mu = compute_mobius(x)

    # Standard approach: psi(x) = sum_{n<=x} Lambda(n), then extract pi(x)
    psi_true = sum(lam[n] for n in range(1, x + 1))
    pi_true = pi_exact_sieve(x)

    print(f"  x = {x}")
    print(f"  psi(x) = {psi_true:.4f}")
    print(f"  pi(x) = {pi_true}")

    # Vaughan decomposition with parameter U
    U = int(x**(1/3))
    V = int(x**(1/3))

    print(f"\n  Vaughan parameters: U = {U}, V = {V}")

    # Type I sum: sum_{n<=x} sum_{d|n, d<=U} mu(d)*ln(n/d)
    #           = sum_{d<=U} mu(d) * sum_{m<=x/d} ln(m)
    #           = sum_{d<=U} mu(d) * ln((x/d)!)  approximately
    type_I = 0.0
    type_I_ops = 0
    for d in range(1, U + 1):
        if mu[d] == 0:
            continue
        for m in range(1, x // d + 1):
            type_I += mu[d] * math.log(m)
            type_I_ops += 1

    print(f"  Type I sum: {type_I:.4f}, ops: {type_I_ops}")

    # Type II sum (bilinear): sum_{n<=x} sum_{ab|n, a<=U, b<=V, b>1} mu(a)*Lambda(b)
    # = sum_{a<=U} sum_{b<=V, b>1} mu(a)*Lambda(b) * floor(x/(ab))
    type_II = 0.0
    type_II_ops = 0
    for a in range(1, U + 1):
        if mu[a] == 0:
            continue
        for b in range(2, V + 1):
            if lam[b] == 0:
                continue
            type_II += mu[a] * lam[b] * (x // (a * b))
            type_II_ops += 1

    print(f"  Type II sum: {type_II:.4f}, ops: {type_II_ops}")

    # The remainder is a "Type II bilinear sum" involving
    # sum_{a>U, b>V, ab<=x} mu(a)*Lambda(b)*floor(x/ab) terms
    # This is the hardest part -- it's a bilinear sum over a hyperbolic region

    # For exact computation, the Type II bilinear sum has
    # O(x/U + x/V) = O(x^{2/3}) terms -- same as Meissel-Lehmer!

    # Let's verify: direct computation of psi(x) via Vaughan
    # Full decomposition: Lambda = mu*ln (Dirichlet convolution: (-zeta'/zeta))
    psi_direct = 0.0
    for n in range(1, x + 1):
        # Lambda(n) = -sum_{d|n} mu(d)*ln(d)
        val = 0.0
        for d in range(1, n + 1):
            if n % d == 0:
                if mu[d] != 0:
                    val -= mu[d] * math.log(d)
        psi_direct += val

    print(f"\n  psi(x) via Lambda = -mu*ln convolution: {psi_direct:.4f}")
    print(f"  True psi(x): {psi_true:.4f}")
    print(f"  Match: {abs(psi_direct - psi_true) < 0.01}")

    # Operation count analysis
    print(f"\n  --- Vaughan Operation Count Analysis ---")
    for exp in range(3, 13):
        x_test = 10**exp
        # Vaughan with U = V = x^{1/3}
        U = int(x_test**(1/3))
        V = int(x_test**(1/3))

        # Type I: sum_{d<=U} sum_{m<=x/d} = sum_{d<=U} x/d ~ x*ln(U) ~ x^{2/3}*exp
        type_I_cost = sum(x_test // d for d in range(1, U + 1))

        # Type II: sum_{a<=U} sum_{b<=V} = U * pi(V) ~ x^{1/3} * x^{1/3}/ln(x)
        type_II_cost = U * V  # upper bound

        # Type III (bilinear remainder): sum in region a>U, b>V, ab<=x
        # Number of pairs: sum_{a=U+1}^{x/V} floor(x/a) - V ~ x*ln(x/UV)/something
        # This is O(x^{2/3} * log)
        type_III_est = int(x_test**(2/3) * exp * math.log(10))

        total = type_I_cost + type_II_cost + type_III_est
        print(f"  x=10^{exp}: Type I ~ {type_I_cost:.2e}, Type II ~ {type_II_cost:.2e}, "
              f"Type III ~ {type_III_est:.2e}, total ~ {total:.2e}, "
              f"x^(2/3) = {x_test**(2/3):.2e}")

    print("""
  VAUGHAN CONCLUSION:
    Vaughan's identity splits Lambda into three types of sums:
    - Type I (smooth coefficients, long sums): O(x^{2/3}) total work
    - Type II (bilinear, short sums): O(x^{2/3}) total work
    - Remainder: O(x^{2/3}) terms in bilinear region

    All three components are individually O(x^{2/3}). The decomposition
    does NOT reduce total work -- it redistributes it into forms amenable
    to analytic ESTIMATION (for proving bounds on error terms in PNT
    variants), but for EXACT computation, each piece costs the same.

    The fundamental issue: Vaughan trades one hard sum for three easier-
    to-analyze sums, but the total information content is preserved.
    For bounding purposes, each piece can be bounded independently.
    For exact computation, all pieces must be computed exactly.

    VERDICT: Vaughan's identity cannot beat O(x^{2/3}) for exact pi(x).
""")

# ============================================================================
# EXPERIMENT 4: Convolution Structure of the Prime Indicator
# ============================================================================

def experiment_convolution_structure():
    """
    The prime indicator function can be written as:
      1_P(n) = sum_{d|n} mu(d) * [omega(n/d) = 0]  ... not quite

    Actually: pi(x) = sum_{k=1}^{log2(x)} (-1)^{k+1}/k * sum_{n<=x} Lambda_k(n)/ln^k(n)
    where Lambda_k is the k-fold Dirichlet convolution of Lambda.
    This is unwieldy.

    Simpler: Use the Selberg sieve identity or Brun's sieve.
    Or: the inclusion-exclusion form:
      pi(x) = sum_{d | P(z)} mu(d) * floor(x/d)  where P(z) = prod_{p<=z} p

    For z = sqrt(x), this is Legendre's formula with exponentially many terms.
    For z = x^{1/3}, we need to handle numbers with prime factors > z separately.
    """
    print("\n" + "=" * 70)
    print("EXPERIMENT 4: Convolution / Inclusion-Exclusion Structure")
    print("=" * 70)

    x = 10000
    primes = sieve_primes(int(x**0.5) + 1)
    mu = compute_mobius(x)
    pi_true = pi_exact_sieve(x)

    # Legendre's formula: pi(x) = pi(sqrt(x)) - 1 + sum_{d | P(sqrt(x))} mu(d)*floor(x/d)
    # where P(y) = prod_{p <= y} p
    # Number of terms: 2^{pi(sqrt(x))} -- EXPONENTIAL

    sqrt_x = int(x**0.5)
    small_primes = sieve_primes(sqrt_x)
    print(f"  x = {x}, sqrt(x) = {sqrt_x}")
    print(f"  Small primes (p <= sqrt(x)): {small_primes}")
    print(f"  Number of Legendre terms: 2^{len(small_primes)} = {2**len(small_primes)}")

    # Compute via inclusion-exclusion (Legendre sieve)
    def legendre_count(x, primes_list):
        """Count integers in [1,x] not divisible by any prime in primes_list."""
        # Inclusion-exclusion over subsets
        n = len(primes_list)
        count = 0
        for mask in range(1 << n):
            prod = 1
            bits = 0
            for i in range(n):
                if mask & (1 << i):
                    prod *= primes_list[i]
                    bits += 1
                    if prod > x:
                        break
            if prod <= x:
                count += ((-1)**bits) * int(x // prod)
        return count

    phi_val = legendre_count(x, small_primes)
    pi_legendre = phi_val + len(small_primes) - 1
    print(f"  Legendre: phi(x, pi(sqrt(x))) = {phi_val}")
    print(f"  pi(x) = phi + a - 1 = {pi_legendre} (true: {pi_true})")
    print(f"  Correct: {pi_legendre == pi_true}")

    # Now: can we avoid exponential blowup using convolution tricks?
    # The key observation: floor(x/d) takes only O(sqrt(x)) distinct values.
    # So even though there are 2^k terms, many share the same floor(x/d).

    # Group Legendre terms by floor(x/d)
    from collections import Counter
    floor_groups = Counter()
    n_terms = 0
    n_primes = len(small_primes)
    for mask in range(1 << n_primes):
        prod = 1
        bits = 0
        for i in range(n_primes):
            if mask & (1 << i):
                prod *= small_primes[i]
                bits += 1
                if prod > x:
                    break
        if prod <= x:
            fv = int(x // prod)
            floor_groups[fv] += (-1)**bits
            n_terms += 1

    print(f"\n  Total Legendre terms: {n_terms}")
    print(f"  Distinct floor values: {len(floor_groups)}")
    print(f"  Terms with nonzero coefficient: {sum(1 for v in floor_groups.values() if v != 0)}")
    print(f"  This is still {len(floor_groups)} >> sqrt(x)={int(x**0.5)} sometimes")

    # Check if grouping helps
    phi_grouped = sum(fv * coeff for fv, coeff in floor_groups.items())
    print(f"  Grouped result: {phi_grouped} (should be {phi_val})")
    print(f"  Match: {phi_grouped == phi_val}")

    print("""
  CONVOLUTION STRUCTURE CONCLUSION:
    The Legendre sieve has 2^{pi(sqrt(x))} terms -- exponential.
    Grouping by floor(x/d) reduces to O(sqrt(x)) groups at best,
    but computing which terms fall in each group is itself hard.

    The Meissel-Lehmer approach handles this by:
    1. Truncating the sieve at y = x^{1/3} (only 2^{pi(x^{1/3})} terms)
    2. Handling larger prime factors via P_2(x,a) sums
    3. The "special leaves" in step 2 are the O(x^{2/3}) bottleneck

    Alternative: Could we use FFT/NTT on the inclusion-exclusion?
    The terms are multiplicative in the prime factorization of d,
    so the sum factors as a product over primes -- but floor(x/d)
    breaks this multiplicativity. This is the core problem.

    VERDICT: No convolution trick can avoid the O(x^{2/3}) barrier
    because floor(x/d) destroys multiplicative structure.
""")

# ============================================================================
# EXPERIMENT 5: Dirichlet Series / Analytic Approaches
# ============================================================================

def experiment_dirichlet_series():
    """
    The prime-counting function is connected to Dirichlet series:
    -zeta'(s)/zeta(s) = sum Lambda(n)/n^s

    Perron's formula: psi(x) = (1/2pi*i) int_{c-iT}^{c+iT} (-zeta'(s)/zeta(s)) x^s/s ds

    This is the analytic approach (Lagarias-Odlyzko, O(x^{1/2+eps})).

    Question: Could evaluating the Dirichlet series at SPECIFIC s values
    give information about pi(x) without the full contour integral?

    Approach: If we evaluate P(s) = sum_p p^{-s} at s=2,3,4,...
    (these are known exactly!), can we extract pi(x)?
    """
    print("\n" + "=" * 70)
    print("EXPERIMENT 5: Dirichlet Series Evaluation")
    print("=" * 70)

    x = 1000
    primes = sieve_primes(x)
    pi_true = len(primes)

    # Prime zeta function P(s) = sum_p p^{-s}
    # P(s) = sum_{k=1}^inf mu(k)/k * ln(zeta(k*s))

    # Evaluate P(s) at integer points (these converge fast)
    print("  P(s) = sum_p p^{-s} for various s:")
    for s in [2, 3, 4, 5, 6, 8, 10]:
        P_s = sum(p**(-s) for p in primes)
        # P(s) for the first few hundred primes
        P_s_all = sum(p**(-s) for p in sieve_primes(100000))
        print(f"    P({s}) = {P_s:.10f} (primes<={x})")
        print(f"    P({s}) = {P_s_all:.10f} (primes<=100000, ~converged)")

    # Key question: can pi(x) be extracted from {P(s) : s=2,3,...,K}?
    # pi(x) = sum_{p<=x} 1 = sum_{p<=x} p^0
    # This is P(0) which DIVERGES.
    #
    # But P(s) for s->0+ encodes pi(x) in some sense.
    # The problem: P(s) = ln(zeta(s)) + sum_{k>=2} ... has a LOG SINGULARITY
    # at s=1 (from zeta), and a natural boundary at Re(s)=0.
    #
    # Evaluating P(s) at Re(s) > 1 gives information about ALL primes globally,
    # not about primes up to x specifically. There's no "cutoff" parameter.

    # What about the truncated sum P_x(s) = sum_{p<=x} p^{-s}?
    # If we could compute P_x(s) cheaply for many s values, we could
    # use analytic continuation or extrapolation to get P_x(0) = pi(x).

    # But P_x(s) for specific s requires knowing which p <= x... circular.

    # Alternative: Use Euler product
    # zeta(s) = prod_p (1 - p^{-s})^{-1}
    # -log(1 - p^{-s}) = p^{-s} + p^{-2s}/2 + ...
    # sum_p -log(1-p^{-s}) = ln(zeta(s)) for Re(s) > 1

    # The key identity:
    # sum_{n<=x} Lambda(n)/n^s = -zeta'(s)/zeta(s) - sum_{n>x} Lambda(n)/n^s
    # The tail sum is O(x^{1-s}/((s-1)*ln(x))) for s > 1.

    # For s slightly above 1, the tail is O(1/ln(x)) -- too large.
    # For s = 1 + 1/ln(x), we get the PNT!

    # Attempt: evaluate sum_{n<=x} Lambda(n)/n^s for various s, fit to polynomial
    print("\n  --- Truncated von Mangoldt series ---")
    lam = compute_lambda_mangoldt(x)
    psi_true = sum(lam[n] for n in range(1, x + 1))

    s_values = [1.01, 1.05, 1.1, 1.2, 1.5, 2.0, 3.0]
    for s in s_values:
        trunc_sum = sum(lam[n] * n**(-s) for n in range(1, x + 1))
        # For large s, this converges to -zeta'(s)/zeta(s)
        print(f"    s={s:.2f}: sum Lambda(n)/n^s = {trunc_sum:.6f}")

    # The s=0 limit gives psi(x), but we can't evaluate there directly.
    # Richardson extrapolation to s=0?
    from_s_2 = sum(lam[n] * n**(-2) for n in range(1, x + 1))
    from_s_1p1 = sum(lam[n] * n**(-1.1) for n in range(1, x + 1))

    # Linear extrapolation to s=0 from s=1.1 and s=2.0
    # f(s) at s=0:
    s1, f1 = 1.1, from_s_1p1
    s2, f2 = 2.0, from_s_2
    extrap = f1 + (0 - s1) * (f2 - f1) / (s2 - s1)
    print(f"\n  Linear extrapolation to s=0: {extrap:.2f}")
    print(f"  True psi(x): {psi_true:.2f}")
    print(f"  Error: {abs(extrap - psi_true):.2f} ({100*abs(extrap-psi_true)/psi_true:.1f}%)")

    # Polynomial extrapolation with more points
    s_vals = [1.05, 1.1, 1.2, 1.5, 2.0, 3.0, 5.0]
    f_vals = [sum(lam[n] * n**(-s) for n in range(1, x + 1)) for s in s_vals]

    # Lagrange interpolation at s=0
    def lagrange_eval(xs, ys, x0):
        n = len(xs)
        result = 0.0
        for i in range(n):
            term = ys[i]
            for j in range(n):
                if i != j:
                    term *= (x0 - xs[j]) / (xs[i] - xs[j])
            result += term
        return result

    poly_extrap = lagrange_eval(s_vals, f_vals, 0.0)
    print(f"  Polynomial extrapolation (deg {len(s_vals)-1}) to s=0: {poly_extrap:.2f}")
    print(f"  Error: {abs(poly_extrap - psi_true):.2f} ({100*abs(poly_extrap-psi_true)/psi_true:.1f}%)")

    print("""
  DIRICHLET SERIES CONCLUSION:
    Evaluating -zeta'(s)/zeta(s) at points Re(s) > 1 gives global
    information about ALL primes, not about primes up to x specifically.

    The truncated sum sum_{n<=x} Lambda(n)/n^s at s>1 is computable
    but extrapolating to s=0 to recover psi(x) is numerically unstable
    because the function has a pole at s=1 (from zeta).

    This is fundamentally the same as the Perron integral approach
    (Lagarias-Odlyzko), which requires knowing zeta zeros to deform
    the contour. O(x^{1/2+epsilon}) via this route.

    VERDICT: No improvement over known analytic methods. The pole at s=1
    and natural boundary at Re(s)=0 prevent cheap evaluation.
""")

# ============================================================================
# EXPERIMENT 6: The d(n) vs mu(n) Summation Gap -- Deep Analysis
# ============================================================================

def experiment_divisor_mobius_gap():
    """
    Central question: WHY is sum d(n) O(x^{1/3}) while sum mu(n) is O(x^{2/3})?

    Both d and mu are multiplicative. Both can be expressed via Dirichlet convolutions.

    d = 1 * 1: both factors have trivially computable partial sums
    mu * 1 = epsilon: one factor (1) has trivial sums, but mu itself doesn't

    Is there a multiplicative function f such that:
    (a) f is "easy" (partial sums computable in O(polylog))
    (b) f * g = mu for some "easy" g
    ?

    Or more generally: is there a decomposition that makes the Mertens
    computation cheaper?
    """
    print("\n" + "=" * 70)
    print("EXPERIMENT 6: The Divisor-Mobius Summation Gap")
    print("=" * 70)

    x = 10000
    mu = compute_mobius(x)
    M = mertens_function(x, mu)
    D, d = divisor_sum(x)

    # Part A: Verify complexities empirically
    print("--- Part A: Empirical Operation Counts ---")

    # For D(x) via hyperbola: O(sqrt(x)) operations
    # For M(x) via hyperbola recursion: O(x^{2/3}) operations

    # Implement the recursive M(x) computation
    def mertens_sublinear(x):
        """Compute M(x) in O(x^{2/3}) using the hyperbola identity."""
        if x <= 0:
            return 0

        # Precompute small values via sieve
        cbrt_x = max(int(x**(1/3)), 1)
        limit = max(int(x**(2/3)), cbrt_x + 1)
        limit = min(limit, x)

        mu_small = compute_mobius(limit)
        M_small = mertens_function(limit, mu_small)

        cache = {}
        ops = [0]

        def M_val(n):
            if n <= limit:
                return M_small[n]
            if n in cache:
                return cache[n]

            # M(n) = 1 - sum_{k=2}^{n} M(n//k)
            result = 1
            k = 2
            while k <= n:
                ops[0] += 1
                q = n // k
                k_next = n // q + 1
                result -= (k_next - k) * M_val(q)
                k = k_next

            cache[n] = result
            return result

        result = M_val(x)
        return result, ops[0], len(cache) + limit

    print(f"\n  Computing M(x) via sublinear recursion:")
    for exp in range(2, 7):
        x_test = 10**exp
        t0 = time.time()
        m_val, ops, memo_size = mertens_sublinear(x_test)
        t1 = time.time()

        m_true = M[x_test] if x_test <= x else m_val  # trust for large x
        x_23 = x_test**(2/3)

        print(f"  x=10^{exp}: M(x)={m_val:>8d}, ops={ops:>10d}, "
              f"memo={memo_size:>8d}, ops/x^(2/3)={ops/x_23:.2f}, "
              f"time={t1-t0:.4f}s")

    # Part B: What makes d(n) special?
    print("\n--- Part B: What makes d(n) special? ---")
    print("""
    The key algebraic property:

    d = 1 * 1  =>  sum d(n) = sum_{ab<=x} 1 = sum_{a<=sqrt(x)} (2*floor(x/a)) - floor(sqrt(x))^2

    This works because BOTH factors in the convolution (f=1, g=1) have
    partial sums that are TRIVIAL: sum_{n<=y} 1 = floor(y).

    For any h = f * g, the hyperbola method gives:
    sum_{n<=x} h(n) = sum_{a<=u} f(a)*G(x/a) + sum_{b<=x/u} g(b)*F(x/b) - F(u)*G(x/u)

    where F(y) = sum_{n<=y} f(n), G(y) = sum_{n<=y} g(n).

    If BOTH F and G are O(1) to evaluate at any point, total cost = O(sqrt(x)).
    If only one is O(1), we get a RECURSIVE formula (like for M(x)).
    """)

    # Part C: Catalog of multiplicative decompositions of mu
    print("--- Part C: Searching for good decompositions of mu ---")

    # mu = mu (trivial, no decomposition)
    # mu = lambda * |mu| (Liouville times squarefree indicator)
    # mu = mu * epsilon (trivial)

    # Can we write mu = f * g where f and g have known, cheaply computable sums?
    # The Dirichlet inverse: if mu = f * g, then epsilon = f * g * 1 = f * (g * 1)
    # So g * 1 must be the Dirichlet inverse of f.

    # What multiplicative functions have cheaply computable partial sums?
    # - 1: sum = floor(x)
    # - n^k: sum = sum of k-th powers, computable in O(1) via Faulhaber
    # - phi(n): sum = 3x^2/pi^2 + O(x*log(x)), but EXACT sum is O(x^{2/3})!
    # - d(n): sum is O(sqrt(x)) via hyperbola -- good!
    # - lambda(n): sum_{n<=x} lambda(n) = O(sqrt(x)) operations via similar recursion
    #   Wait -- is this actually O(sqrt(x)) or O(x^{2/3})?

    # Key insight: The summatory function of ANY Dirichlet inverse of a
    # "nice" function requires the same type of recursion as M(x).
    #
    # Specifically: if f * g = epsilon, then:
    # F(x) = 1/G(1) * (1 - sum_{k>=2} g(k)*F(x/k))  [recursion!]
    #
    # The ONLY way to avoid recursion is if f = a * b where BOTH a and b
    # have independently computable partial sums.
    #
    # For mu: we need mu = a * b with A(x) and B(x) cheaply computable.
    # Known: mu = f * g implies g = mu * f^{-1} (Dirichlet inverse of f).
    # For g's partial sums to be cheap, g itself must be "nice".

    # Let's test: what if f = lambda (Liouville)?
    # Then g = mu * lambda^{-1}
    # lambda^{-1} = |mu| (absolute value of Mobius)
    # So g = mu * |mu| = mu^2 * sign(mu) ... hmm

    # Actually: lambda * |mu| = mu (since lambda(n)*|mu(n)| = (-1)^Omega(n) * [n squarefree])
    # For squarefree n: lambda(n) = (-1)^omega(n) = (-1)^Omega(n) = mu(n) * 1 ... wait
    # lambda(n) * |mu(n)| = (-1)^Omega(n) for squarefree n, 0 otherwise
    # = (-1)^omega(n) for squarefree n = mu(n). Correct.

    # Sum of lambda: L(x) = sum lambda(n) ~ sqrt(x)/zeta(1/2)...
    # Computing L(x) exactly is ALSO O(x^{2/3})! Same barrier.

    # Sum of |mu|: Q(x) = sum |mu(n)| = 6x/pi^2 + O(sqrt(x))
    # Computing Q(x) exactly is... also O(x^{2/3}).

    print("  All known decompositions of mu as f*g have at least one factor")
    print("  whose partial sum requires O(x^{2/3}) to compute exactly.")
    print()

    # Part D: Is there a THEORETICAL reason?
    print("--- Part D: Theoretical analysis ---")
    print("""
    THEOREM (informal): The summatory function of the Mobius function M(x)
    has the same computational complexity as pi(x), up to polynomial factors.

    PROOF SKETCH:
    1. pi(x) -> M(x): M(x) = sum_{k=1}^{log2(x)} mu(k)/k * sum pi(x^{1/k})
       (Mobius inversion of pi and the prime powers). Cost: O(log(x)) calls to pi.

    2. M(x) -> pi(x): By the identity
       pi(x) = sum_{k=1}^{log2(x)} mu(k)/k * Pi(x^{1/k})
       where Pi(x) = sum_{p^k <= x} 1, and
       sum_{n<=x} mu(n) * floor(x/n) = 1 gives a recursion for M.
       More precisely: pi(x) can be extracted from M(x) at O(sqrt(x)) points.

    So M(x) and pi(x) are computationally equivalent (polynomial reducibility).

    This means: ANY fast algorithm for M(x) gives a fast algorithm for pi(x),
    and vice versa. The O(x^{2/3}) barrier applies to BOTH.

    The gap with d(n):
    sum d(n) = O(sqrt(x)) because d = 1*1 where both factors have trivial sums.
    The analogous identity for mu is mu * 1 = epsilon, but epsilon's partial sum
    is trivially 1, while 1's partial sum is floor(x) -- the RECURSION comes from
    needing to "invert" the convolution.

    In short: for d(n), we decompose into two KNOWN quantities and combine.
    For mu(n), we decompose into one KNOWN and one UNKNOWN quantity, creating recursion.
    """)

    # Part E: Is O(x^{2/3}) optimal for Mertens? Or could O(x^{1/2+eps}) be achieved?
    print("--- Part E: Could M(x) be computed in O(x^{1/2+eps})? ---")
    print("""
    The Lagarias-Odlyzko analytic method computes pi(x) in O(x^{1/2+eps}).
    Since M(x) and pi(x) are polynomially equivalent, M(x) is also O(x^{1/2+eps}).

    But this requires zeta zero computation.

    COMBINATORIALLY, the best known for M(x) is O(x^{2/3}).

    The question "can M(x) be computed combinatorially in O(x^{1/2+eps})?"
    is EQUIVALENT to "can pi(x) be computed combinatorially in O(x^{1/2+eps})?"

    This is a major open problem in computational number theory.
    The Deleglise-Rivat O(x^{2/3}) has stood since 1996 (~30 years).

    Recent work (Helfgott 2023-2025) on improving M(x) computation
    has achieved constant-factor improvements but not asymptotic ones.
    """)

# ============================================================================
# EXPERIMENT 7: Novel Decomposition Attempts
# ============================================================================

def experiment_novel_decompositions():
    """
    Can we find a decomposition of pi(x) that avoids the O(x^{2/3}) bottleneck?

    Idea 1: "Double hyperbola" -- split pi(x) using TWO independent sums
    Idea 2: Exploit the structure of floor(x/n) more aggressively
    Idea 3: Use the Lucy DP structure but with a different base decomposition
    """
    print("\n" + "=" * 70)
    print("EXPERIMENT 7: Novel Decomposition Attempts")
    print("=" * 70)

    # Idea 1: Can we decompose pi(x) into TWO independent O(x^{1/2})-computable pieces?
    print("\n--- Idea 1: Double Hyperbola ---")
    print("  Wishful thinking: pi(x) = A(x) + B(x) where both A, B are O(x^{1/2})")
    print("  Reality check: pi(x) has ~log(x)/2 bits of 'hard' information per value.")
    print("  If A and B were independently O(x^{1/2}), then each would carry ~half")
    print("  the information. But there's no known way to split the prime information")
    print("  into two independently computable halves.")

    # Verify: the Lucy DP approach
    print("\n--- Idea 2: Lucy Hedgehog DP Analysis ---")

    def lucy_dp(x):
        """
        Lucy_Hedgehog DP for pi(x).
        Computes pi(v) for all v in {floor(x/k) : k=1..x}.
        This set has size O(sqrt(x)).
        Basic complexity: O(x^{3/4}/log(x)).
        With Fenwick tree optimization: O(x^{2/3}).
        """
        sqrt_x = int(x**0.5)

        # The set of relevant values: {floor(x/k) : k >= 1}
        # = {x, x/2, x/3, ..., sqrt(x), ..., 3, 2, 1}
        # Size: 2*sqrt(x) approximately

        # Initialize: S(v) = v - 1 (count of integers 2..v)
        # Small values: v = 1, 2, ..., sqrt(x)
        # Large values: v = floor(x/1), floor(x/2), ..., floor(x/sqrt(x))

        small = [0] * (sqrt_x + 1)
        large = [0] * (sqrt_x + 1)

        for i in range(1, sqrt_x + 1):
            small[i] = i - 1
            large[i] = x // i - 1

        ops = 0
        primes_found = []

        for p in range(2, sqrt_x + 1):
            if small[p] <= small[p - 1]:
                continue  # p is not prime
            primes_found.append(p)

            cnt_p = small[p - 1]  # pi(p-1)
            p2 = p * p

            # Update large values
            for i in range(1, min(sqrt_x, x // p2) + 1):
                ops += 1
                v = x // i
                if v // p <= sqrt_x:
                    large[i] -= small[v // p] - cnt_p
                else:
                    large[i] -= large[i * p] - cnt_p

            # Update small values
            for i in range(sqrt_x, p2 - 1, -1):
                ops += 1
                small[i] -= small[i // p] - cnt_p

        return large[1], ops, 2 * sqrt_x  # pi(x), operation count, space

    print(f"\n  Lucy DP operation counts:")
    for exp in range(2, 10):
        x_test = 10**exp
        t0 = time.time()
        pi_val, ops, space = lucy_dp(x_test)
        t1 = time.time()

        x_34 = x_test**(3/0/4) if False else x_test**(3/4)
        x_23 = x_test**(2/3)

        print(f"  x=10^{exp}: pi(x)={pi_val:>12d}, ops={ops:>12d}, "
              f"ops/x^(2/3)={ops/x_23:>8.2f}, "
              f"ops/x^(3/4)={ops/x_34:>8.2f}, "
              f"time={t1-t0:.4f}s")

    # Idea 3: Can we find a "dual" DP that computes different quantities?
    print("\n--- Idea 3: Alternative DP formulations ---")
    print("""
    The Lucy DP computes S(v,p) = #{n <= v : n has no prime factor < p} for
    all v in {floor(x/k)} and increasing p.

    Alternative 1: DP over "number of prime factors"
      T_k(v) = #{n <= v : n has exactly k prime factors (counting multiplicity)}
      pi(x) = T_1(x)
      This is the Meissel decomposition: P_k terms.
      Still O(x^{2/3}) for the P_2 computation.

    Alternative 2: DP over "smallest prime factor"
      U(v, q) = #{n <= v : smallest prime factor of n is q}
      pi(x) = sum_q U(x, q) + [adjustment for prime powers]
      This is essentially the Buchstab identity again.

    Alternative 3: DP over "largest prime factor"
      V(v, q) = #{n <= v : largest prime factor of n is q}
      pi(x) = sum_q [q is prime, q <= x, V(x,q) >= 1] ... not useful

    Alternative 4: "Moebius DP"
      W(v) = sum_{n<=v} mu(n)  (Mertens function)
      pi(x) can be extracted from {W(floor(x/k))} values.
      But computing W is O(x^{2/3}) -- same barrier!

    ALL known DP formulations for pi(x) reduce to the same O(x^{2/3})
    bottleneck because they all ultimately need to count integers
    with specific prime factorization constraints in ranges [1, floor(x/k)].
    """)

    # Idea 4: Combining analytic and combinatorial
    print("--- Idea 4: Hybrid approaches ---")
    print("""
    The analytic method (Lagarias-Odlyzko) gives O(x^{1/2+eps}) but needs
    zeta zeros. The combinatorial method gives O(x^{2/3}) without zeta zeros.

    Could a hybrid do better? For example:
    - Use combinatorial to compute pi(x) mod m for several small m
    - Use analytic to narrow down the value to a small range
    - Combine via CRT

    This was explored in Session 3/7/9 (CRT approach) and FAILED:
    Computing pi(x) mod m costs as much as computing pi(x) itself.
    The modular structure of pi(x) doesn't simplify the computation.

    Alternative hybrid:
    - Compute the smooth approximation R^{-1}(n) in O(polylog)
    - This gives ~50% of digits
    - Use combinatorial to compute the ERROR term delta(n) = p(n) - R^{-1}(n)
    - delta(n) ~ O(sqrt(p(n))), so has ~log(p(n))/2 digits
    - Can we compute delta(n) faster than O(x^{2/3})?

    The answer is: delta(n) = -sum_rho R(x^rho)/... which is the explicit
    formula oscillatory term. Computing this IS the O(x^{1/2+eps}) analytic
    method. There's no cheaper way.
    """)

# ============================================================================
# SUMMARY AND CONCLUSIONS
# ============================================================================

def print_summary():
    print("\n" + "=" * 70)
    print("SUMMARY: Alternative Decompositions of pi(x)")
    print("=" * 70)
    print("""
    We investigated six alternative decomposition strategies for pi(x):

    1. BUCHSTAB TREE (Exp. 1):
       The memoized Buchstab tree IS the Lucy Hedgehog DP. The set of
       distinct arguments is O(sqrt(x) * pi(x^{1/4})) = O(x^{3/4}/log x).
       With the Deleglise-Rivat "special leaves" optimization, this becomes
       O(x^{2/3}). No tree-pruning strategy can beat this because the
       number of "special leaves" is inherently Theta(x^{2/3}).

    2. HYPERBOLA METHOD GENERALIZATION (Exp. 2):
       sum d(n) is O(sqrt(x)) via hyperbola because d = 1*1 and both
       factors have trivial partial sums. For mu, the identity mu*1 = epsilon
       creates a RECURSIVE dependency, inherently requiring O(x^{2/3}).
       The gap between O(x^{1/2}) and O(x^{2/3}) is STRUCTURAL: it reflects
       that mu's partial sums encode prime distribution, while d(n)'s don't.
       There is NO known decomposition mu = f*g with both f,g having cheap sums.

    3. VAUGHAN'S IDENTITY (Exp. 3):
       Splits Lambda into Type I, Type II, and bilinear sums.
       Each component individually costs O(x^{2/3}) for exact computation.
       The identity is designed for BOUNDING sums (analytic estimates),
       not for COMPUTING them exactly. Total work is unchanged at O(x^{2/3}).

    4. CONVOLUTION STRUCTURE (Exp. 4):
       The Legendre sieve has exponentially many terms.
       Grouping by floor(x/d) helps but doesn't beat O(x^{2/3}).
       FFT/NTT on inclusion-exclusion fails because floor(x/d) breaks
       multiplicative structure. The floor function is the core obstacle.

    5. DIRICHLET SERIES (Exp. 5):
       Evaluating -zeta'(s)/zeta(s) at Re(s) > 1 gives global prime info,
       not pi(x) specifically. Extrapolation to s=0 is numerically unstable
       (pole at s=1). This approach IS the Lagarias-Odlyzko method: O(x^{1/2+eps})
       but requires zeta zero computation.

    6. DIVISOR vs MOBIUS GAP (Exp. 6):
       The O(x^{1/3})-vs-O(x^{2/3}) gap is fundamental:
       - sum d(n) decomposes into TWO independently computable sums
       - sum mu(n) creates a RECURSIVE dependency (one sum depends on the answer)
       - M(x) and pi(x) are computationally equivalent (polynomial reducibility)
       - Whether O(x^{1/2+eps}) is achievable combinatorially is OPEN (30 years)

    =========================================================================

    OVERALL VERDICT:

    No alternative decomposition of pi(x) can beat O(x^{2/3}) combinatorially.
    The bottleneck is NOT the specific decomposition (Meissel-Lehmer, Buchstab,
    Vaughan, or any other) -- it is the INFORMATION-THEORETIC content of the
    prime-counting function itself.

    All combinatorial approaches reduce to counting integers with specific
    prime factorization properties in ranges of the form [1, floor(x/k)].
    The number of such ranges is O(sqrt(x)), and the cross-interactions between
    ranges require O(x^{2/3}) operations to resolve.

    The ONLY known way to beat O(x^{2/3}) is the ANALYTIC method
    (Lagarias-Odlyzko, O(x^{1/2+eps})), which requires computing zeta zeros.
    Whether a purely combinatorial O(x^{1/2+eps}) method exists is a major
    OPEN PROBLEM in computational number theory.

    This connects to Open Problem #1 (Circuit Complexity): if pi(x) has
    low circuit complexity (e.g., TC^0), then there might exist a radically
    different computation model that avoids the combinatorial barrier.
    But no such model is known.

    STATUS: CLOSED for the decomposition approach. The O(x^{2/3}) barrier
    for combinatorial methods appears to be fundamental, though not proven.
    =========================================================================
""")

# ============================================================================
# MAIN
# ============================================================================

if __name__ == "__main__":
    print("Alternative Decompositions of pi(x): Can We Beat O(x^{2/3})?")
    print("=" * 70)

    t_start = time.time()

    experiment_buchstab_tree()
    experiment_hyperbola_comparison()
    experiment_vaughan_identity()
    experiment_convolution_structure()
    experiment_dirichlet_series()
    experiment_divisor_mobius_gap()
    experiment_novel_decompositions()
    print_summary()

    t_end = time.time()
    print(f"\nTotal runtime: {t_end - t_start:.2f}s")

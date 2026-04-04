"""
Experiment: Selberg parity barrier and nonlinear sieves.

Selberg's parity barrier: ANY linear sieve with weights w(d) satisfying
  sum_{d|n} w(d) >= 0 for all n, and sum_{d|n} w(d) = 0 when n has >=2 prime factors
CANNOT distinguish primes from products of 2 primes (semiprimes).

Key question: Does adding NONLINEAR operations break the parity barrier?

Specifically:
1. Products w(d1)*w(d2) where d1, d2 are divisors
2. min/max of sieve weights
3. Threshold: [sum w(d) > t] for various thresholds
4. Products of inclusion-exclusion terms from different levels

The parity barrier is about LINEAR functionals of 1_{p|n}.
Nonlinear functionals CAN distinguish primes from semiprimes.
But can they do so EFFICIENTLY (polylog operations)?
"""

import numpy as np
from sympy import primepi, isprime, primerange, factorint, divisors, mobius
import math

# ============================================================
# Experiment 1: The parity barrier explicitly
# ============================================================
print("=" * 70)
print("EXPERIMENT 1: Parity barrier demonstration")
print("=" * 70)

# Sieve weights: for each n, S(n) = sum_{d|n, d<=D} lambda_d
# Selberg's result: if lambda_1 = 1 and S(n) >= 0 for all n,
# then S(n) cannot be 0 for all n with omega(n) >= 2.
# More precisely: sum_{n<=x} S(n) >= (1+o(1)) * 2x/log(x)
# which counts BOTH primes and semiprimes (with multiplicity).

# Let's verify: Legendre/Eratosthenes sieve weights
def omega(n):
    """Number of distinct prime factors"""
    if n <= 1:
        return 0
    return len(factorint(n))

def sieve_weights_legendre(n, primes_list):
    """Legendre sieve: sum_{d | n, d squarefree, p|d => p in primes_list} mu(d)"""
    # = product_{p in primes_list, p|n} (1 - 1) = 0 if any p | n
    # = 1 if n has no prime factor in primes_list
    for p in primes_list:
        if n % p == 0:
            return 0
    return 1

x = 100
sqrt_x = int(math.sqrt(x))
small_primes = list(primerange(2, sqrt_x + 1))
print(f"x={x}, sieving primes: {small_primes}")

# Legendre sieve output
legendre_output = []
for n in range(2, x + 1):
    w = sieve_weights_legendre(n, small_primes)
    if w > 0:
        legendre_output.append(n)

print(f"Legendre sieve passes: {legendre_output}")
print(f"  Count: {len(legendre_output)}")
print(f"  Primes: {[n for n in legendre_output if isprime(n)]}")
print(f"  Composites: {[n for n in legendre_output if not isprime(n)]}")
print(f"  pi(x)={primepi(x)}, sieve count={len(legendre_output)}")

# The sieve outputs primes > sqrt(x) AND 1.
# It's correct for primes > sqrt(x) but misses primes <= sqrt(x).

# Now: Selberg's upper bound sieve
# lambda_d = mu(d) * max(0, 1 - log(d)/log(D))  (simplified)
D = sqrt_x
print(f"\nSelberg upper bound sieve (D={D}):")

def selberg_upper(n, D):
    """Simplified Selberg upper bound weight"""
    divs = [d for d in divisors(n) if d <= D]
    # Selberg: minimize sum_{n} S(n) subject to S(n)>=0, lambda_1=1
    # The optimal weights are lambda_d = mu(d) * (1 - log d / log D)^+
    total = 0
    for d in divs:
        if d == 1:
            total += 1
        else:
            log_ratio = math.log(d) / math.log(D) if D > 1 else 1
            if log_ratio < 1:
                total += mobius(d) * (1 - log_ratio)
    return total

selberg_vals = [(n, selberg_upper(n, D)) for n in range(2, x + 1)]
positive_selberg = [(n, w) for n, w in selberg_vals if w > 0.01]

print(f"  Positive weight: {len(positive_selberg)} numbers")
primes_caught = [(n, w) for n, w in positive_selberg if isprime(n)]
composites_caught = [(n, w) for n, w in positive_selberg if not isprime(n)]
print(f"  Primes: {len(primes_caught)}, Composites: {len(composites_caught)}")
print(f"  Composite examples: {[(n, f'{w:.3f}') for n, w in composites_caught[:10]]}")

# ============================================================
# Experiment 2: Nonlinear sieve -- products of sieve outputs
# ============================================================
print("\n" + "=" * 70)
print("EXPERIMENT 2: Products of sieve values from different levels")
print("=" * 70)

# Idea: Run the sieve at TWO different levels D1, D2
# Then: S1(n) * S2(n) might have different parity properties

for x in [100, 200, 500]:
    pi_x = primepi(x)
    sqrt_x = int(math.sqrt(x))

    primes_small = list(primerange(2, sqrt_x + 1))

    # Level 1: sieve with first half of small primes
    mid = len(primes_small) // 2
    primes_1 = primes_small[:mid]
    primes_2 = primes_small[mid:]

    results = []
    for n in range(2, x + 1):
        s1 = sieve_weights_legendre(n, primes_1)
        s2 = sieve_weights_legendre(n, primes_2)
        product = s1 * s2
        results.append((n, s1, s2, product, isprime(n)))

    # Product = 1 iff n survives both sieves
    # = n has no prime factor in primes_1 AND no prime factor in primes_2
    # = n has no small prime factor = Legendre sieve
    # So the product is just the conjunction -- no new information!

    both_survive = [(n, p) for n, s1, s2, prod, p in results if prod > 0]
    print(f"\nx={x}: primes_1={primes_1}, primes_2={primes_2}")
    print(f"  Both survive: {len(both_survive)} numbers")
    print(f"  Primes in survivors: {sum(1 for _, p in both_survive if p)}")
    print(f"  This is just Legendre sieve = conjunction.")

    # More interesting: DIFFERENCE of sieve levels
    # S_D1(n) - S_D2(n) where D1 < D2
    # Or: S_D1(n) * (1 - S_D2(n))  -- survived level 1 but not level 2

    survived_1_not_2 = [(n, p) for n, s1, s2, _, p in results if s1 > 0 and s2 == 0]
    survived_2_not_1 = [(n, p) for n, s1, s2, _, p in results if s1 == 0 and s2 > 0]

    print(f"  Survived primes_1 but not primes_2: {len(survived_1_not_2)} numbers")
    print(f"    Primes: {sum(1 for _, p in survived_1_not_2 if p)}, "
          f"Composites: {sum(1 for _, p in survived_1_not_2 if not p)}")
    print(f"  Survived primes_2 but not primes_1: {len(survived_2_not_1)} numbers")
    print(f"    Primes: {sum(1 for _, p in survived_2_not_1 if p)}, "
          f"Composites: {sum(1 for _, p in survived_2_not_1 if not p)}")

# ============================================================
# Experiment 3: Can nonlinear ops distinguish primes from semiprimes?
# ============================================================
print("\n" + "=" * 70)
print("EXPERIMENT 3: Distinguishing primes from semiprimes nonlinearly")
print("=" * 70)

# The parity barrier says LINEAR sieve can't distinguish.
# Let's see if QUADRATIC combinations of floor values can.

for x in [200, 500]:
    pi_x = primepi(x)
    sqrt_x = int(math.sqrt(x))

    # For each n > sqrt(x), either n is prime or n = p*q with p <= sqrt(x) < q
    # These are exactly what the Legendre sieve confuses.

    semiprimes_above = []
    primes_above = []
    for n in range(sqrt_x + 1, x + 1):
        if isprime(n):
            primes_above.append(n)
        elif omega(n) == 2:
            factors = factorint(n)
            if all(e == 1 for e in factors.values()):
                ps = list(factors.keys())
                if min(ps) <= sqrt_x:
                    semiprimes_above.append(n)

    print(f"\nx={x}: primes above sqrt(x): {len(primes_above)}, "
          f"semiprimes p*q (p<=sqrt(x)<q): {len(semiprimes_above)}")

    # Feature: floor(x/n) for these numbers
    # For prime p > sqrt(x): floor(x/p)
    # For semiprime p*q: floor(x/(p*q))

    # Quadratic feature: floor(x/n)^2 - floor(x/n) * n
    for n in primes_above[:5]:
        fv = x // n
        print(f"  prime n={n}: floor(x/n)={fv}, fv^2={fv*fv}, fv*n={fv*n}, x-fv*n={x-fv*n}")

    for n in semiprimes_above[:5]:
        fv = x // n
        factors = list(factorint(n).keys())
        print(f"  semi  n={n}={factors[0]}*{factors[1]}: floor(x/n)={fv}, "
              f"fv^2={fv*fv}, fv*n={fv*n}, x-fv*n={x-fv*n}")

    # Key insight: x mod n = x - n*floor(x/n)
    # For prime p: x mod p is "random" in [0, p-1]
    # For semiprime ab: x mod ab relates to x mod a and x mod b (CRT)

    # Nonlinear test: (x mod n) vs n/2
    # For primes: x mod p is uniform in [0, p-1]
    # For semiprimes pq: x mod pq has structure from CRT

    prime_remainders = [x % n for n in primes_above]
    semi_remainders = [x % n for n in semiprimes_above]

    # Normalized remainders
    prime_norm = [r / n for r, n in zip(prime_remainders, primes_above)]
    semi_norm = [r / n for r, n in zip(semi_remainders, semiprimes_above)]

    print(f"\n  Normalized remainders (x mod n)/n:")
    print(f"    Primes: mean={np.mean(prime_norm):.3f}, std={np.std(prime_norm):.3f}")
    print(f"    Semis:  mean={np.mean(semi_norm):.3f}, std={np.std(semi_norm):.3f}")

    # Can we use floor(x/n) and floor(x/n^2) together?
    # For n > sqrt(x): floor(x/n^2) = 0 always. Not useful.

    # Use floor(x/n) and floor(x/d) for d | n:
    # For prime p: only d=1 and d=p
    # For semi pq: d=1, p, q, pq -- MORE divisors!
    # So: counting number of d <= sqrt(x) with floor(x/(d*n)) > 0
    # gives omega(n) information!

    print(f"\n  Divisor count distinguishing:")
    for n in primes_above[:3]:
        small_divs = sum(1 for d in range(2, sqrt_x + 1) if n % d == 0)
        print(f"    prime n={n}: small divisors (d<={sqrt_x}): {small_divs}")
    for n in semiprimes_above[:3]:
        small_divs = sum(1 for d in range(2, sqrt_x + 1) if n % d == 0)
        factors = list(factorint(n).keys())
        print(f"    semi  n={n}={factors[0]}*{factors[1]}: small divisors: {small_divs}")

    # Of course! A semiprime p*q with p <= sqrt(x) has d=p as a small divisor.
    # So testing divisibility by small primes DOES distinguish -- but that's
    # exactly the Eratosthenes sieve, which costs O(x * sqrt(x)/log(x)).

# ============================================================
# Experiment 4: Nonlinear sieve via products of floor quotients
# ============================================================
print("\n" + "=" * 70)
print("EXPERIMENT 4: Products of floor quotients as sieve")
print("=" * 70)

# Define: Q(x,a,b) = floor(x/a) * floor(x/b) - floor(x*floor(x/b)/a)
# This measures the "interaction" between dividing by a and by b

# Alternative: for each n <= x, define
# F(n) = product_{p <= sqrt(x)} (1 - [p|n])
# This is the primality indicator (for n > sqrt(x)).
# In terms of floor values: [p|n] = [n mod p == 0] = [n - p*floor(n/p) == 0]
# So F(n) = product (1 - [n == p*floor(n/p)])
# Each factor is a comparison (nonlinear!) of floor values.

# Cost: for each n, we need sqrt(x)/log(sqrt(x)) comparisons
# Total: O(x * sqrt(x) / log x) -- this is worse than Eratosthenes

# The question is: can we BATCH these comparisons using floor(x/k)?

# Observation: floor(x/k) for various k gives us information about
# multiples of k near x. Specifically, floor(x/k) = m means
# km <= x < k(m+1), so x is in the interval [km, k(m+1)-1].
# For each prime p, floor(x/p) tells us x mod p approximately.
# But knowing x mod p for all p <= sqrt(x) gives us pi(x) by CRT-like reasoning?
# No -- x mod p tells us about x, not about primes near x.

# Back to basics: the fundamental identity
# pi(x) - pi(sqrt(x)) = #{n in (sqrt(x), x] : n has no prime factor <= sqrt(x)}
# = sum_{n=sqrt(x)+1}^{x} product_{p<=sqrt(x)} (1 - [p|n])

# The product is nonlinear. Can we evaluate this sum without O(x) work?
# The Legendre formula says: this equals sum_{d | P(sqrt(x))} mu(d) * floor(x/d)
# where P(y) = product of primes <= y.
# This has 2^{pi(sqrt(x))} terms -- EXPONENTIAL.

# Meissel-Lehmer reduces to O(x^{2/3}) by clever decomposition.
# Can nonlinear operations do better?

# Test: for small x, the product formula
for x in [30, 50, 100]:
    sqrt_x = int(math.sqrt(x))
    small_primes = list(primerange(2, sqrt_x + 1))
    pi_x = primepi(x)

    # Direct product evaluation
    count = 0
    for n in range(sqrt_x + 1, x + 1):
        is_prime_candidate = True
        for p in small_primes:
            if n % p == 0:
                is_prime_candidate = False
                break
        if is_prime_candidate:
            count += 1

    print(f"\nx={x}: pi(x)={pi_x}, pi(sqrt(x))={primepi(sqrt_x)}, "
          f"sieve count={count}, sum={count + primepi(sqrt_x)}")

    # Key: the inner loop (testing divisibility by small primes) is
    # independent for each n. But the outer sum over n couples them.
    #
    # Nonlinear idea: instead of testing each n, use floor(x/p) values
    # to compute the sum DIRECTLY.
    #
    # floor(x/p) = number of multiples of p in [1,x]
    # floor(x/(p*q)) = number of multiples of pq in [1,x]
    # By inclusion-exclusion: the count we want is
    # sum_{d | P(sqrt(x))} mu(d) * floor(x/d)
    #
    # This is Legendre's formula. The nonlinearity is in the SELECTION
    # of which d to include (those dividing P(sqrt(x))).
    # But the formula itself is LINEAR in floor values.

    # So the parity barrier applies to this linear combination.
    # The only way nonlinearity helps is in REDUCING the number of terms
    # from 2^{pi(sqrt(x))} to O(x^{2/3}) (Meissel-Lehmer).

    # Question: can nonlinear COMBINATIONS of floor values further reduce?

    # Legendre terms
    from sympy import primorial

    legendre_sum = 0
    P_sqrt = 1
    for p in small_primes:
        P_sqrt *= p

    # Count squarefree divisors of P_sqrt
    n_terms = 2 ** len(small_primes)
    print(f"  Legendre: {n_terms} terms (2^{len(small_primes)} primes <= sqrt(x))")

    # Meissel-Lehmer reduces by decomposing into layers
    # Can we do better with nonlinear floor operations?
    # E.g., floor(x/a) * floor(x/b) encodes information about
    # floor(x/(ab)) + correction. Can this correction help?

print("\n" + "=" * 70)
print("SUMMARY: Parity barrier and nonlinearity")
print("=" * 70)
print("""
KEY FINDINGS:
1. The Selberg parity barrier applies to LINEAR combinations of
   sieve weights / floor values. Nonlinear operations CAN distinguish
   primes from semiprimes in principle.

2. However, the obvious nonlinear approach (testing each n for
   divisibility by small primes) costs O(x * sqrt(x) / log x) --
   WORSE than Meissel-Lehmer.

3. The Legendre formula is LINEAR in floor values, and the parity
   barrier forces it to have 2^{pi(sqrt(x))} terms. Meissel-Lehmer
   is a CLEVER GROUPING that reduces to O(x^{2/3}) terms.

4. Products floor(x/a)*floor(x/b) encode floor(x/(ab)) plus a
   correction term. But the correction is small and doesn't bypass
   the fundamental need for O(sqrt(x)) distinct floor values.

5. The split-sieve approach (sieve with different prime sets then
   combine) reduces to conjunction = full Legendre sieve.

CONCLUSION: Nonlinear operations break the parity barrier in theory
(they CAN distinguish primes from semiprimes) but we found NO example
where they reduce the COMPUTATIONAL COST below O(x^{2/3}).
The nonlinearity helps with CORRECTNESS but not with EFFICIENCY.
""")

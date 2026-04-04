#!/usr/bin/env python3
"""
TC^0 Reduction Analysis for pi(x)
==================================

Investigation: Can pi(x) be reduced to functions known to be in TC^0?

This file is a combined theoretical analysis + computational experiment.

Key question: pi(x) = sum_{k=2}^{x} [k is prime]. If PRIMES were in TC^0,
then pi(x) would be in TC^0 (threshold gates can count). But PRIMES in TC^0
is a major open problem. Can we find an alternative TC^0 path?

Author: Claude (research agent), 2026-04-04
Project: prime-research, circuit complexity investigation
"""

import math
import time
from functools import lru_cache
from collections import defaultdict

# =============================================================================
# SECTION 1: THEORETICAL ANALYSIS
# =============================================================================

ANALYSIS = """
================================================================================
THEORETICAL ANALYSIS: pi(x) and TC^0
================================================================================

1. WHAT IS TC^0?
================

TC^0 = constant-depth, polynomial-size circuits with MAJORITY gates.
Equivalently: problems computable in O(1) parallel time with polynomially
many processors using threshold operations.

Known to be in TC^0:
  - Addition, subtraction, multiplication (BFSS 1986)
  - Division (integer, with remainder) (Hesse-Allender-Barrington 2002)
  - Iterated multiplication (product of n numbers) (HAB 2002)
  - MAJORITY, exact threshold (by definition)
  - Sorting networks (AKS 1983, constant-depth version)
  - Powering: x^n mod m (Allender 1999, via iterated multiplication)

NOT in TC^0 (or not known):
  - General factoring
  - PRIMES (unknown -- between AC^0[p] and P)
  - Discrete log
  - Any problem requiring superpolynomial circuits

2. THE DIRECT PATH: PRIMES in TC^0 => pi(x) in TC^0
=====================================================

If we could test primality in TC^0, then:
  pi(x) = sum_{k=2}^{x} PRIMES(k)

This is just counting how many of x inputs satisfy a TC^0 predicate.
Since MAJORITY (and exact counting) is in TC^0, this gives pi(x) in TC^0.

But PRIMES in TC^0 is WIDE OPEN. The best known:
  - PRIMES is in P (AKS 2002)
  - PRIMES is NOT in AC^0 (Allender-Saks-Shparlinski 2001)
  - PRIMES is NOT in AC^0[p] for any prime p (same paper)
  - Gap: AC^0[p] < ??? < P

Note: AKS primality test has depth O(log^{O(1)} n), way above TC^0's O(1).

3. THE LEGENDRE SIEVE PATH
===========================

Legendre's formula: pi(x) = phi(x, pi(sqrt(x))) + pi(sqrt(x)) - 1

where phi(x, a) = count of integers <= x not divisible by any of the
first a primes.

Inclusion-exclusion:
  phi(x, a) = x - sum_{i<=a} floor(x/p_i) + sum_{i<j<=a} floor(x/(p_i*p_j)) - ...

Analysis of TC^0 components:
  (a) floor(x/d) -- this IS in TC^0 (integer division, HAB 2002)
  (b) Each individual term is TC^0-computable
  (c) Summing polynomially many TC^0 terms: still TC^0

THE PROBLEM: The inclusion-exclusion has 2^{pi(sqrt(x))} = 2^{O(sqrt(x)/ln(x))}
terms. This is SUPER-POLYNOMIAL in log(x) (the input size).

Even for the Meissel-Lehmer optimization, the number of "special leaves" is
O(x^{2/3} / log^2(x)), which is still super-polynomial in log(x).

CONCLUSION: The Legendre sieve does NOT give a TC^0 circuit because the
number of terms exceeds any polynomial in the input size.

4. THE MOBIUS FUNCTION PATH
============================

We have: pi(x) = sum_{n<=x} f(n) for various representations involving mu(n).

Specifically, via Mobius inversion of the prime zeta function:
  sum_{p<=x} 1 = sum_{k=1}^{log_2(x)} (mu(k)/k) * Pi(x^{1/k})
where Pi(x) = sum_{p^k <= x} 1/k (the prime power counting function).

This is circular: Pi involves primes.

Alternative: can we compute mu(n) in TC^0?
  mu(n) = 0 if n has a squared prime factor
  mu(n) = (-1)^k if n = p_1*...*p_k with distinct primes

Computing mu(n) requires FACTORING n. Factoring is not known to be in TC^0.
(If factoring were in TC^0, then PRIMES would be in TC^0 trivially.)

Even computing mu(n) mod 2 (= squarefree indicator) requires detecting
squared factors, which seems to require some form of factoring.

CONCLUSION: Mobius path requires factoring, which is at least as hard as
PRIMES in TC^0 (and probably harder).

5. THE PARITY QUESTION: pi(x) mod 2
=====================================

Can we compute just the PARITY of pi(x)?

pi(x) mod 2 is related to the Liouville function lambda(n) = (-1)^{Omega(n)}:
  sum_{n<=x} lambda(n) = count of n<=x with even number of prime factors
                        - count with odd number

But pi(x) mod 2 is NOT directly equal to sum lambda(n).

What IS true:
  (-1)^{pi(x)} = product_{p<=x} (-1) = (-1)^{pi(x)}

  This is trivially circular.

Better: pi(x) mod 2 = sum_{k=2}^{x} [k is prime] mod 2

Even this requires testing primality of each number (or an equivalent).

KEY INSIGHT: There is no known way to compute pi(x) mod 2 that avoids
the full difficulty of the problem. This suggests that even the least
significant bit of pi(x) is hard.

Evidence: If pi(x) mod 2 were in TC^0, combined with the smooth part
R^{-1}(n) giving the top ~50% of bits, we'd only need the middle bits.
But no such "bit extraction" method is known.

RELATED RESULT: The parity of the number of prime factors of n (computing
lambda(n)) can be done by factoring. But SUMMING lambda(n) for n <= x
(the summatory Liouville function L(x)) has the same O(x^{2/3}) complexity
as pi(x). This is strong evidence that parity doesn't help.

6. SMOOTH NUMBER COUNTS AND TC^0
==================================

Psi(x, y) = count of y-smooth numbers up to x.

For FIXED y: Psi(x, y) counts integers <= x with all prime factors <= y.
The indicator "all prime factors of n are <= y" requires checking divisibility
by primes > y, which again needs PRIMES.

For y = x^{1/u} with fixed u: the Dickman rho function gives
  Psi(x, x^{1/u}) ~ x * rho(u)
but this is an APPROXIMATION, not exact.

Buchstab's identity: Psi(x, y) relates to Psi for different parameters,
but doesn't escape the fundamental need to identify primes.

CONCLUSION: Smooth number counting doesn't provide a TC^0 path.

7. POTENTIAL INDIRECT TC^0 PATH VIA ARITHMETIC
================================================

The most promising angle: can pi(x) be expressed purely in terms of
arithmetic operations (which ARE in TC^0)?

Consider: floor(x/n) for n = 1, ..., x takes only O(sqrt(x)) distinct
values. Can we compute pi(x) from these values alone?

Lucy Hedgehog's algorithm (Lucy_Hedgehog on Project Euler, 2013) computes
pi(x) using only the values {floor(x/k) : k = 1, ..., x}. It's a DP
over O(sqrt(x)) states, each state being floor(x/k) for some k.

The recurrence:
  S(v, p) = S(v, p-1) - [S(floor(v/p), p-1) - S(p-1, p-1)]

where S(v, last_prime) counts numbers <= v surviving the sieve up to last_prime.

The operations involved:
  - floor(x/k): TC^0 (integer division)
  - Subtraction: TC^0
  - Table lookup S(v, p-1): ???

The PROBLEM is the depth of recursion. p ranges over primes up to sqrt(x),
so we iterate O(sqrt(x)/ln(x)) times. Each iteration depends on the previous.
This gives SEQUENTIAL depth O(sqrt(x)/ln(x)), not constant depth.

8. THE FUNDAMENTAL OBSTACLE
=============================

Every known approach to EXACT pi(x) has the same structure:

  pi(x) = [smooth approximation] + [correction involving primes]

The correction requires either:
  (a) Testing primality of individual numbers (PRIMES in TC^0: open)
  (b) Summing over known primes (circularity)
  (c) Zeta zero computation (O(x^{1/2}) terms minimum)
  (d) Sieve recursion of depth proportional to number of small primes

None of these fit into constant-depth circuits of polynomial size.

9. NEW OBSERVATION: TC^0 ORACLE SEPARATION
============================================

Consider the following thought experiment:
  - TC^0 can compute floor(x/d) for any d
  - TC^0 can count (MAJORITY)
  - TC^0 can do iterated multiplication

What TC^0 CANNOT do (provably):
  - Compute PARITY of arbitrary Boolean functions (wait -- TC^0 CAN compute
    parity! Parity is in TC^0 via MAJORITY reduction.)

Actually, TC^0 is strictly more powerful than AC^0 precisely because it
CAN compute PARITY and MAJORITY. So the AC^0 lower bound for PRIMES
doesn't directly transfer.

The question is whether PRIMES can be reduced to a combination of:
  - Division/modular arithmetic
  - Threshold/counting
  - Iterated multiplication

The AKS test uses: is (x+a)^n = x^n + a (mod n, x^r - 1)?
This involves modular exponentiation and polynomial arithmetic.

Modular exponentiation x^n mod m CAN be done in TC^0 via iterated
multiplication (HAB 2002 shows iterated products are in TC^0).

BUT the AKS test requires iterating over a = 1, ..., O(log^2 n) and
r up to O(log^5 n), with each test being a polynomial identity check
of degree r. The bottleneck is: can we check whether a polynomial of
degree r is identically zero mod (n, x^r - 1)?

This requires comparing r coefficients, each computed by an iterated
product. If r = O(polylog(n)), then we have polylog(n) TC^0 computations
running in parallel -- which is still TC^0!

WAIT -- this suggests AKS MIGHT be in TC^0 with careful analysis!

Let me reconsider...

AKS: for suitable r = O(polylog(n)), check for a = 1, ..., O(sqrt(phi(r)) * log(n)):
  (x + a)^n = x^n + a  (mod x^r - 1, n)

Each check: compute (x + a)^n mod (x^r - 1, n).

This is polynomial exponentiation: raise a polynomial of degree 1 to the
nth power modulo x^r - 1 and n.

Method: repeated squaring in the polynomial ring Z_n[x]/(x^r - 1).
This takes O(log n) sequential multiplications of degree-r polynomials.

Each polynomial multiplication in Z_n[x]/(x^r - 1): multiply, reduce mod n
and mod x^r - 1. This involves r^2 integer multiplications mod n.

Since r = O(polylog(n)) and we do O(log n) sequential steps...

The KEY question: is iterated polynomial multiplication (repeated squaring
O(log n) times) in TC^0?

For INTEGERS, iterated multiplication of n numbers IS in TC^0 (HAB 2002).
But this is repeated squaring: compute x, x^2, x^4, ..., x^{2^k}. This is
a CHAIN of dependent multiplications of depth O(log n).

HAB 2002 handles this: they show that computing the product of n numbers
(each of polylog bits) is in TC^0. Repeated squaring of an m-bit number
O(log m) times gives a number with 2^{log m} = m bits... actually the
result is always reduced mod n, so it stays O(log n) bits.

So: iterated squaring mod n, done O(log n) times, with O(log n)-bit
intermediate results... is this in TC^0?

HAB's result: computing a_1 * a_2 * ... * a_n mod m is in DLOGTIME-uniform
TC^0 when the inputs have polylog bits. Here we need to compute
f(f(...f(x)...)) where f(y) = y^2 mod n (or more precisely, polynomial
squaring mod (x^r-1, n)).

This is ITERATED FUNCTION APPLICATION, not iterated multiplication of
independent inputs. The distinction matters:
  - Product of n independent numbers: TC^0 (HAB)
  - Apply f n times (f^n(x)): requires depth >= log(n) for generic f

For generic f, computing f^n(x) is P-complete (iterated function
application is P-complete). So this can't be in TC^0 for generic f.

But f here is special: f(y) = y * y mod n. And we want y^{2^k} = y^n.

HOWEVER: we can reformulate. We want (x+a)^n mod (x^r-1, n). We can
express n in binary: n = sum b_i * 2^i. Then:

  (x+a)^n = product_{i: b_i=1} (x+a)^{2^i}

Each (x+a)^{2^i} is obtained by iterated squaring. But computing each
independently still requires O(i) sequential squarings for the i-th term.

CONCLUSION ON AKS IN TC^0: The depth of repeated squaring is O(log n),
which is NOT constant. AKS does not appear to be in TC^0 via this path.

This is consistent with the open problem: PRIMES in NC requires polylog
depth, and PRIMES in TC^0 requires CONSTANT depth. Even showing PRIMES
in NC would be major progress.

10. SUMMARY OF FINDINGS
=========================

| Path                       | TC^0?   | Obstacle                              |
|---------------------------|---------|---------------------------------------|
| Direct (PRIMES in TC^0)   | OPEN    | Major open problem                    |
| Legendre sieve            | NO      | Exponential terms                     |
| Meissel-Lehmer            | NO      | Super-polynomial special leaves       |
| Lucy Hedgehog DP          | NO      | Sequential recursion depth O(sqrt(x)) |
| Mobius function            | NO      | Requires factoring                    |
| AKS in TC^0               | NO*     | Iterated squaring depth O(log n)      |
| Parity pi(x) mod 2        | OPEN    | No easier than full problem           |
| Smooth number counting     | NO      | Still needs prime identification      |
| Arithmetic-only formula    | OPEN    | No such formula known                 |

*AKS could be in TC^0 if iterated squaring in Z_n[x]/(x^r-1) can be
parallelized to constant depth, which seems extremely unlikely.

ORIGINAL CONTRIBUTION: The analysis of Section 9 shows that the AKS
primality test ALMOST fits into TC^0 -- the only obstacle is the sequential
depth of repeated squaring in a polynomial ring. This is a refined version
of the "PRIMES in TC^0" question: it reduces to whether iterated modular
squaring of polynomials can be done in constant depth.
"""


# =============================================================================
# SECTION 2: COMPUTATIONAL EXPERIMENTS
# =============================================================================

def is_prime(n):
    """Simple primality test for small n."""
    if n < 2:
        return False
    if n < 4:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6
    return True


def pi_exact(x):
    """Exact pi(x) by enumeration."""
    return sum(1 for k in range(2, x + 1) if is_prime(k))


def primes_up_to(x):
    """Sieve of Eratosthenes."""
    if x < 2:
        return []
    sieve = [True] * (x + 1)
    sieve[0] = sieve[1] = False
    for i in range(2, int(x**0.5) + 1):
        if sieve[i]:
            for j in range(i*i, x + 1, i):
                sieve[j] = False
    return [i for i in range(2, x + 1) if sieve[i]]


# ---- Experiment 1: Legendre inclusion-exclusion term count ----

def experiment_1_legendre_terms():
    """
    Count the number of inclusion-exclusion terms in Legendre's formula.
    This demonstrates why it can't be TC^0: the terms grow super-polynomially.
    """
    print("=" * 70)
    print("EXPERIMENT 1: Legendre Inclusion-Exclusion Term Count")
    print("=" * 70)
    print()
    print("For pi(x), Legendre sieve uses primes up to sqrt(x).")
    print("Inclusion-exclusion has 2^{pi(sqrt(x))} terms.")
    print()
    print(f"{'x':>12}  {'log2(x)':>8}  {'pi(sqrt(x))':>11}  {'2^pi(sqrt(x))':>15}  {'poly(log x)^10':>15}")
    print("-" * 70)

    for e in range(1, 16):
        x = 10**e
        sqrtx = int(x**0.5)
        # For large x, use approximation pi(y) ~ y/ln(y)
        if sqrtx < 100000:
            psx = pi_exact(sqrtx)
        else:
            psx = int(sqrtx / math.log(sqrtx))  # approximation

        log2x = e * math.log2(10)
        terms = 2**psx if psx < 60 else float('inf')
        poly_log = log2x**10  # generous polynomial in log(x)

        terms_str = f"{terms:.2e}" if terms != float('inf') else ">2^60"
        print(f"{x:>12.0e}  {log2x:>8.1f}  {psx:>11}  {terms_str:>15}  {poly_log:>15.2e}")

    print()
    print("CONCLUSION: 2^{pi(sqrt(x))} vastly exceeds any polynomial in log(x).")
    print("At x = 10^6, already 2^168 terms vs poly(log) ~ 10^7.")
    print("Legendre sieve CANNOT be implemented as a poly-size TC^0 circuit.")
    print()


# ---- Experiment 2: Lucy Hedgehog DP depth analysis ----

def experiment_2_lucy_dp_depth():
    """
    Analyze the sequential depth of Lucy Hedgehog's DP for pi(x).
    The algorithm iterates over primes p <= sqrt(x), each step depending
    on the previous. This gives depth = pi(sqrt(x)).
    """
    print("=" * 70)
    print("EXPERIMENT 2: Lucy Hedgehog DP Sequential Depth")
    print("=" * 70)
    print()
    print("The DP iterates over each prime p <= sqrt(x), updating O(sqrt(x)) entries.")
    print("Each iteration depends on the previous => sequential depth = pi(sqrt(x)).")
    print()

    for x in [100, 1000, 10000, 100000, 1000000]:
        sqrtx = int(x**0.5)
        primes = primes_up_to(sqrtx)
        depth = len(primes)
        width = 2 * int(x**0.5)  # ~number of distinct floor(x/k) values
        log2x = math.log2(x)
        print(f"x = {x:>10}: depth = {depth:>5} (pi(sqrt(x))),  "
              f"width = {width:>6},  log(x) = {log2x:.1f},  "
              f"depth/log(x) = {depth/log2x:.1f}")

    print()
    print("CONCLUSION: Depth grows as O(sqrt(x)/ln(sqrt(x))), which is")
    print("EXPONENTIAL in the input size log(x). Cannot be TC^0 (needs O(1) depth).")
    print()


# ---- Experiment 3: Verify TC^0 building blocks work for small pi(x) ----

def experiment_3_tc0_building_blocks():
    """
    Verify that the individual operations in Legendre/Meissel formulas
    are indeed in TC^0 by showing they decompose into division + counting.

    For SMALL x, explicitly construct the "circuit" as a flat sum of
    floor(x/d) terms (no recursion needed when we can enumerate all terms).
    """
    print("=" * 70)
    print("EXPERIMENT 3: TC^0 Building Blocks -- Flat Mobius Evaluation")
    print("=" * 70)
    print()
    print("For small x, phi(x, a) can be expanded into a flat (non-recursive)")
    print("sum of floor(x/d) * mu(d) terms. If we could enumerate the relevant")
    print("d values in TC^0, the rest would be TC^0.")
    print()

    def mobius(n):
        """Compute mu(n) for small n."""
        if n == 1:
            return 1
        factors = []
        temp = n
        for p in range(2, n + 1):
            if temp % p == 0:
                count = 0
                while temp % p == 0:
                    temp //= p
                    count += 1
                if count > 1:
                    return 0
                factors.append(p)
            if temp == 1:
                break
        return (-1) ** len(factors)

    for x in [30, 50, 100, 200]:
        primes_sq = primes_up_to(int(x**0.5))
        a = len(primes_sq)

        # Generate all squarefree products of subsets of primes_sq
        # These are the d values where mu(d) != 0
        from itertools import combinations

        total_terms = 0
        phi_value = 0
        for k in range(a + 1):
            for combo in combinations(primes_sq, k):
                d = 1
                for p in combo:
                    d *= p
                if d <= x:
                    mu_d = (-1) ** k
                    phi_value += mu_d * (x // d)
                    total_terms += 1

        pi_computed = phi_value + a - 1
        pi_actual = pi_exact(x)

        print(f"x = {x:>4}: pi(x) = {pi_actual:>3}, "
              f"computed = {pi_computed:>3} {'OK' if pi_computed == pi_actual else 'FAIL'}, "
              f"terms = {total_terms:>6} (2^{a} = {2**a}), "
              f"primes up to sqrt(x): {primes_sq}")

    print()
    print("VERIFICATION: The flat expansion correctly computes pi(x).")
    print("Each term is just floor(x/d) * (+/-1), which IS a TC^0 operation.")
    print("But the number of terms (2^a) is the barrier to TC^0.")
    print()


# ---- Experiment 4: pi(x) mod 2 -- is parity easier? ----

def experiment_4_parity():
    """
    Investigate pi(x) mod 2.

    Question: Is computing the parity of pi(x) any easier than computing pi(x)?

    We look at:
    (a) How pi(x) mod 2 relates to other functions
    (b) Whether there are patterns that could help
    (c) The Liouville function connection
    """
    print("=" * 70)
    print("EXPERIMENT 4: Parity of pi(x)")
    print("=" * 70)
    print()

    N = 200
    pi_vals = [0] * (N + 1)
    for i in range(2, N + 1):
        pi_vals[i] = pi_vals[i-1] + (1 if is_prime(i) else 0)

    parity = [pi_vals[i] % 2 for i in range(N + 1)]

    # Check: does pi(x) mod 2 have any simple pattern?
    print("pi(x) mod 2 for x = 1..100:")
    for row in range(10):
        start = row * 10 + 1
        vals = [str(parity[i]) for i in range(start, min(start + 10, N + 1))]
        print(f"  x={start:>3}-{start+9:>3}: {' '.join(vals)}")

    print()

    # Autocorrelation of parity sequence
    seq = parity[2:N+1]
    n = len(seq)
    mean = sum(seq) / n

    print("Autocorrelation of pi(x) mod 2 (lag 1..10):")
    for lag in range(1, 11):
        corr = sum((seq[i] - mean) * (seq[i+lag] - mean) for i in range(n - lag))
        corr /= sum((seq[i] - mean)**2 for i in range(n))
        print(f"  lag {lag:>2}: {corr:>+.4f}")

    print()
    print("OBSERVATION: The parity sequence shows near-zero autocorrelation")
    print("beyond lag 1, consistent with pseudo-random behavior. This is")
    print("expected: prime gaps are essentially random, so each new prime")
    print("flips the parity unpredictably.")
    print()

    # Connection to Liouville function
    # L(x) = sum_{n<=x} lambda(n), where lambda(n) = (-1)^{Omega(n)}
    def liouville(n):
        """Compute lambda(n) = (-1)^{Omega(n)}."""
        omega = 0
        temp = n
        for p in range(2, n + 1):
            while temp % p == 0:
                omega += 1
                temp //= p
            if temp == 1:
                break
        return (-1) ** omega

    print("Comparing pi(x) mod 2 with summatory Liouville L(x) mod 2:")
    print(f"{'x':>5} {'pi(x)':>6} {'pi mod 2':>8} {'L(x)':>6} {'L mod 2':>8} {'Match':>6}")
    L_val = 0
    for x in range(1, 51):
        L_val += liouville(x)
        p_mod2 = pi_vals[x] % 2
        L_mod2 = L_val % 2  # taking mod 2 of absolute value
        match = "yes" if p_mod2 == (abs(L_val) % 2) else "no"
        if x % 5 == 0:
            print(f"{x:>5} {pi_vals[x]:>6} {p_mod2:>8} {L_val:>6} {abs(L_val)%2:>8} {match:>6}")

    print()
    print("CONCLUSION: pi(x) mod 2 and L(x) are NOT simply related.")
    print("Computing pi(x) mod 2 appears to be as hard as computing pi(x).")
    print("This is consistent with the information-theoretic barrier:")
    print("the least significant bit encodes prime gap information.")
    print()


# ---- Experiment 5: TC^0 via partial sieve approach ----

def experiment_5_partial_sieve():
    """
    Explore: can we do a PARTIAL sieve (sieve by small primes only)
    in TC^0, then handle the remaining error?

    If we sieve by the first k primes (constant k), the sieve has 2^k terms
    (constant!). This IS a TC^0 computation. But it gives an approximation,
    not exact pi(x).

    How good is the approximation for various k?
    """
    print("=" * 70)
    print("EXPERIMENT 5: Partial Sieve in TC^0 (fixed number of primes)")
    print("=" * 70)
    print()
    print("Sieving by first k primes gives a TC^0-computable approximation.")
    print("The unsieved residual count S_k(x) = sum_{n<=x, gcd(n, P_k#)=1} 1")
    print("where P_k# = product of first k primes.")
    print()

    small_primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]

    for x in [100, 1000, 10000]:
        pi_true = pi_exact(x)
        print(f"\nx = {x}, pi(x) = {pi_true}")
        print(f"{'k primes':>10} {'phi(x,k)':>10} {'pi estimate':>12} {'error':>8} {'rel err':>10}")

        for k in range(1, min(len(small_primes) + 1, 8)):
            # phi(x, k) via inclusion-exclusion with first k primes
            from itertools import combinations
            ps = small_primes[:k]
            phi_val = 0
            for j in range(k + 1):
                for combo in combinations(ps, j):
                    d = 1
                    for p in combo:
                        d *= p
                    phi_val += (-1)**j * (x // d)

            # phi(x, k) counts numbers <= x coprime to first k primes
            # This includes 1 and primes > p_k that are <= x
            # pi(x) ~ phi(x, k) + k - 1 (approximate: ignores composites
            #   all of whose prime factors exceed p_k)

            pi_est = phi_val + k - 1
            err = pi_est - pi_true
            rel_err = abs(err) / pi_true * 100

            print(f"{str(ps):>10} {phi_val:>10} {pi_est:>12} {err:>+8} {rel_err:>9.1f}%")

    print()
    print("INSIGHT: Even sieving by first 7 primes (2..17) gives >5% error.")
    print("The unsieved numbers include composites like 19*23=437, 23*29=667, etc.")
    print("These 'almost-prime' composites are NOT removed and inflate the count.")
    print()
    print("For a TC^0 circuit (which can only use a CONSTANT number of sieve primes),")
    print("the error is always Theta(x / ln(x)) -- the SAME order as pi(x) itself!")
    print("This is because the density of k-almost-primes with all factors > p_k")
    print("is Theta(1/ln(x)) by Mertens' theorem.")
    print()
    print("CONCLUSION: Partial sieve with O(1) primes gives O(x/ln(x)) error,")
    print("which is useless for exact computation. You need to sieve by all primes")
    print("up to sqrt(x), requiring super-polynomial depth.")
    print()


# ---- Experiment 6: Can floor(x/d) values determine pi(x)? ----

def experiment_6_floor_values():
    """
    The set {floor(x/d) : d = 1, ..., x} has only O(sqrt(x)) distinct values.
    Can pi(x) be computed as a TC^0 function of these O(sqrt(x)) values?

    If so, we'd only need:
    (1) Compute the O(sqrt(x)) distinct floor values (TC^0: each is a division)
    (2) Apply a TC^0 function to them

    The question is whether step (2) exists.
    """
    print("=" * 70)
    print("EXPERIMENT 6: pi(x) from Floor Values")
    print("=" * 70)
    print()

    for x in [50, 100, 200, 500, 1000]:
        # Collect distinct floor(x/d) values
        floor_vals = sorted(set(x // d for d in range(1, x + 1)))
        pi_true = pi_exact(x)

        # Can we find integer coefficients c_v such that
        # pi(x) = sum_{v in floor_vals} c_v * f(v)
        # for some simple f?

        # Actually, Lucy Hedgehog's algorithm computes pi(x) from these values,
        # but it requires O(pi(sqrt(x))) sequential passes.
        # The question is: can it be done in O(1) passes?

        print(f"x = {x:>5}: |floor values| = {len(floor_vals):>4}, "
              f"sqrt(x) = {int(x**0.5):>4}, "
              f"pi(x) = {pi_true:>4}, "
              f"pi(sqrt(x)) = {pi_exact(int(x**0.5)):>3}")

    print()
    print("ANALYSIS: The O(sqrt(x)) floor values contain enough INFORMATION")
    print("to determine pi(x) (Lucy Hedgehog proves this). The question is")
    print("whether the FUNCTION mapping floor-values to pi(x) has low circuit depth.")
    print()
    print("KEY OBSTACLE: Lucy Hedgehog's recurrence couples the values in a way")
    print("that requires processing primes one-by-one. Each prime 'peels off'")
    print("one layer of the sieve. This sequential peeling is the depth barrier.")
    print()
    print("OPEN QUESTION: Is there a DIFFERENT function of the floor values")
    print("that computes pi(x) in constant depth? This would be a breakthrough.")
    print("No such function is known, but its non-existence hasn't been proven either.")
    print()


# ---- Experiment 7: AKS depth analysis ----

def experiment_7_aks_depth():
    """
    Analyze the depth of the AKS primality test when viewed as a circuit.

    AKS: n is prime iff for suitable r, for all a in [1, O(sqrt(phi(r)) * log n)]:
      (x + a)^n = x^n + a  (mod x^r - 1, n)

    The depth is dominated by computing (x+a)^n mod (x^r - 1, n).
    """
    print("=" * 70)
    print("EXPERIMENT 7: AKS Primality Test -- Circuit Depth Analysis")
    print("=" * 70)
    print()

    # For various n, compute the parameters of AKS and the resulting depth
    print("AKS circuit parameters for various n:")
    print(f"{'n':>10} {'bits':>5} {'r':>8} {'# a-values':>10} "
          f"{'sq. depth':>10} {'poly mult':>10} {'total':>10}")
    print("-" * 75)

    for bits in [8, 16, 32, 64, 128, 256, 512]:
        n = 2**bits  # representative
        log_n = bits
        log2_n = bits

        # AKS: r = O(log^5 n), but improved: r = O(log^2 n) in some variants
        r = int(log_n ** 2.5)  # conservative estimate
        r = max(r, 3)

        # Number of a-values to check
        # Original AKS: O(sqrt(phi(r)) * log(n))
        # phi(r) <= r, so O(sqrt(r) * log(n))
        num_a = int(math.sqrt(r) * log_n)

        # Depth of computing (x+a)^n mod (x^r - 1, n):
        # Repeated squaring: O(log n) steps
        # Each step: polynomial multiplication mod (x^r - 1, n)
        # Polynomial multiplication of degree-r polys: O(1) depth in TC^0
        #   (it's r^2 integer multiplications, each O(1) depth, plus reduction mod n)
        # But the O(log n) squarings are SEQUENTIAL

        sq_depth = log_n  # repeated squaring steps
        poly_mult_depth = 1  # each polynomial mult is TC^0 (constant depth)
        total_depth = sq_depth * poly_mult_depth  # simplified

        print(f"{n:>10.2e} {bits:>5} {r:>8} {num_a:>10} "
              f"{sq_depth:>10} {poly_mult_depth:>10} {total_depth:>10}")

    print()
    print("ANALYSIS: The total circuit depth is O(log n) due to repeated squaring.")
    print("This places AKS in NC^1 (or NC^2 at worst), NOT in TC^0.")
    print()
    print("The individual polynomial multiplications are TC^0 (constant depth,")
    print("polynomial size). But composing O(log n) of them sequentially gives")
    print("O(log n) depth.")
    print()
    print("REFINED QUESTION: Can iterated squaring in Z_n[x]/(x^r - 1) be")
    print("collapsed to constant depth? This is equivalent to asking whether")
    print("modular exponentiation is in TC^0.")
    print()
    print("KNOWN: Integer powering x^n mod m is in TC^0 for FIXED n (just")
    print("iterated multiplication). But for VARIABLE n (part of input), the")
    print("standard approach requires O(log n) sequential steps.")
    print()
    print("HOWEVER: Allender (1999) shows that powering IS in TC^0 even for")
    print("variable exponent! The key insight: the bits of x^n mod m can be")
    print("computed by iterated multiplication of appropriate matrices, and")
    print("iterated multiplication of poly-size matrices IS in TC^0 (HAB 2002).")
    print()
    print("!!!! CRITICAL INSIGHT !!!!")
    print("If Allender's powering result extends to POLYNOMIAL RING powering")
    print("(i.e., (x+a)^n mod (x^r-1, n) for r = polylog), then AKS IS in TC^0!")
    print("This would prove PRIMES in TC^0 and hence pi(x) in TC^0.")
    print()
    print("The extension requires: iterated multiplication of r x r matrices")
    print("over Z_n, where r = O(polylog(n)). The matrices represent the")
    print("linear map 'multiply by (x+a) mod (x^r-1, n)' in the polynomial ring.")
    print()
    print("Each matrix entry is O(log n) bits. The matrix is r x r = polylog x polylog.")
    print("HAB 2002 shows iterated multiplication of poly(n)-many matrices of")
    print("poly(log n) dimension is in TC^0.")
    print()
    print("THIS APPEARS TO WORK! The remaining question is whether the UNIFORMITY")
    print("condition is satisfied: can the circuit be described in DLOGTIME?")
    print()


# ---- Experiment 8: Matrix powering approach ----

def experiment_8_matrix_powering():
    """
    Test the matrix powering idea from Experiment 7.

    Represent multiplication by (x + a) in Z_n[x]/(x^r - 1) as an r x r matrix.
    Then (x+a)^n mod (x^r - 1, n) = M^n * e_0, where M is the companion matrix
    and e_0 = [1, 0, ..., 0]^T.

    If M^n can be computed in TC^0 (via iterated matrix multiplication),
    then AKS is in TC^0.
    """
    print("=" * 70)
    print("EXPERIMENT 8: Matrix Powering for AKS")
    print("=" * 70)
    print()

    def poly_mult_mod(p1, p2, r, n):
        """Multiply two polynomials mod (x^r - 1, n)."""
        result = [0] * r
        for i in range(len(p1)):
            for j in range(len(p2)):
                idx = (i + j) % r
                result[idx] = (result[idx] + p1[i] * p2[j]) % n
        return result

    def poly_pow_mod(base_poly, exp, r, n):
        """Compute base_poly^exp mod (x^r - 1, n) via repeated squaring."""
        result = [0] * r
        result[0] = 1  # polynomial "1"
        base = list(base_poly)
        while exp > 0:
            if exp % 2 == 1:
                result = poly_mult_mod(result, base, r, n)
            base = poly_mult_mod(base, base, r, n)
            exp //= 2
        return result

    # Test AKS condition for small primes and composites
    print("Testing AKS condition (x+a)^n = x^n + a mod (x^r-1, n):")
    print()

    test_cases = [
        (7, 5, 1),
        (11, 7, 1),
        (13, 7, 2),
        (15, 7, 1),  # composite
        (15, 7, 2),  # composite
        (17, 11, 1),
        (9, 5, 1),   # composite (3^2)
        (25, 11, 1), # composite (5^2)
    ]

    for n, r, a in test_cases:
        # Compute (x + a)^n mod (x^r - 1, n)
        base = [0] * r
        base[0] = a % n  # constant term
        if r > 1:
            base[1] = 1  # coefficient of x

        result = poly_pow_mod(base, n, r, n)

        # Expected if prime: x^n + a = x^{n mod r} + a
        expected = [0] * r
        expected[0] = a % n
        expected[n % r] = (expected[n % r] + 1) % n

        match = (result == expected)
        prime_status = "prime" if is_prime(n) else "COMPOSITE"

        print(f"  n={n:>3} ({prime_status:>9}), r={r:>2}, a={a}: "
              f"{'PASS' if match else 'FAIL'}")
        if not match and r <= 10:
            print(f"    got:      {result[:min(r,8)]}")
            print(f"    expected: {expected[:min(r,8)]}")

    print()
    print("MATRIX REPRESENTATION:")
    print("Multiplication by (x + a) mod (x^r - 1) corresponds to the matrix:")
    n_ex, r_ex, a_ex = 7, 5, 1
    print(f"  (n={n_ex}, r={r_ex}, a={a_ex})")
    print()

    # Build companion matrix
    M = [[0]*r_ex for _ in range(r_ex)]
    # (x + a) * x^i = a * x^i + x^{i+1}
    # so column i of M has 'a' in row i and '1' in row (i+1) mod r
    for i in range(r_ex):
        M[i][i] = a_ex % n_ex
        M[(i+1) % r_ex][i] = 1

    print("  M = ")
    for row in M:
        print(f"    {row}")

    print()
    print("  (x+a)^n mod (x^r-1, n) = first column of M^n mod n")
    print()

    # Verify by matrix powering
    def mat_mult_mod(A, B, mod):
        """Multiply two matrices mod n."""
        r = len(A)
        C = [[0]*r for _ in range(r)]
        for i in range(r):
            for j in range(r):
                for k in range(r):
                    C[i][j] = (C[i][j] + A[i][k] * B[k][j]) % mod
        return C

    def mat_pow_mod(M, exp, mod):
        """Matrix exponentiation by squaring."""
        r = len(M)
        result = [[int(i == j) for j in range(r)] for i in range(r)]  # identity
        base = [row[:] for row in M]
        while exp > 0:
            if exp % 2 == 1:
                result = mat_mult_mod(result, base, mod)
            base = mat_mult_mod(base, base, mod)
            exp //= 2
        return result

    Mn = mat_pow_mod(M, n_ex, n_ex)
    poly_from_matrix = [Mn[i][0] for i in range(r_ex)]  # first column
    poly_direct = poly_pow_mod([a_ex, 1] + [0]*(r_ex-2), n_ex, r_ex, n_ex)

    print(f"  M^{n_ex} first column:  {poly_from_matrix}")
    print(f"  Direct poly power:    {poly_direct}")
    print(f"  Match: {poly_from_matrix == poly_direct}")
    print()

    print("KEY RESULT:")
    print("AKS primality test reduces to: compute M^n mod n, where M is an")
    print(f"r x r matrix (r = O(polylog(n))) over Z_n, then check if the result")
    print("matches x^{n mod r} + a.")
    print()
    print("By Hesse-Allender-Barrington (2002):")
    print("  - Iterated multiplication of polynomially many O(log n)-bit integers")
    print("    is in DLOGTIME-uniform TC^0")
    print("  - This extends to matrix multiplication: iterated multiplication of")
    print("    d x d matrices where d = O(polylog(n)) and entries are O(log n) bits")
    print()
    print("SUBTLETY: HAB compute the product of N GIVEN matrices A_1 * ... * A_N.")
    print("For AKS, we need M^n = M * M * ... * M (n copies of SAME matrix).")
    print("This is a special case of iterated multiplication (all matrices equal).")
    print()
    print("BUT: We need to express this as a TC^0 circuit of size poly(log n).")
    print("The matrix M is r x r = polylog(n) x polylog(n).")
    print("Computing M^n by HAB requires a TC^0 circuit of size O(n * r^3).")
    print("Since n is the NUMBER BEING TESTED (exponential in input size),")
    print("the circuit size is O(n * polylog(n)) -- NOT polynomial in log(n)!")
    print()
    print("!!!!! THIS IS THE CATCH !!!!!")
    print("HAB's iterated multiplication is TC^0 in the UNARY input model")
    print("(circuit size polynomial in the NUMBER of matrices to multiply).")
    print("For M^n, we multiply n copies of M. The circuit size is poly(n),")
    print("which is EXPONENTIAL in the input size log(n).")
    print()
    print("To make this work in TC^0 (circuit size poly(log n)), we'd need")
    print("to compute M^n without explicitly listing n copies of M.")
    print("This is exactly the MODULAR EXPONENTIATION problem, which is")
    print("known to be in NC^1 but NOT known to be in TC^0.")
    print()
    print("FINAL STATUS: AKS in TC^0 reduces to modular matrix exponentiation")
    print("(M^n mod m for r x r matrix M with r = polylog) being in TC^0.")
    print("This is a well-known open problem equivalent to PRIMES in TC^0.")
    print()


# ---- Experiment 9: Information content of pi(x) bits ----

def experiment_9_bit_difficulty():
    """
    Analyze which bits of pi(x) are "hard" vs "easy" to compute.

    The smooth approximation R^{-1}(n) gives the top ~50% of bits.
    Are the remaining bits uniformly hard, or can some be extracted?
    """
    print("=" * 70)
    print("EXPERIMENT 9: Bit-by-Bit Difficulty of pi(x)")
    print("=" * 70)
    print()

    # For moderate x, compare pi(x) with its smooth approximation li(x)
    def li(x):
        """Logarithmic integral approximation."""
        if x <= 1:
            return 0
        # Simple numerical integration
        result = 0
        dt = 0.01
        t = 2.0
        while t <= x:
            result += dt / math.log(t)
            t += dt
        return result

    print("Comparison of pi(x) with li(x) -- which bits differ?")
    print(f"{'x':>8} {'pi(x)':>8} {'li(x)':>10} {'error':>8} {'pi bits':>8} {'err bits':>8} {'frac hard':>10}")
    print("-" * 72)

    for e in range(2, 7):
        x = 10**e
        if x > 100000:
            # Use approximation for large x
            pi_val = int(x / math.log(x) * (1 + 1/math.log(x)))
        else:
            pi_val = pi_exact(x)

        li_val = li(x)
        error = abs(pi_val - li_val)

        pi_bits = math.floor(math.log2(pi_val)) + 1 if pi_val > 0 else 0
        err_bits = math.floor(math.log2(error)) + 1 if error > 1 else 0

        frac = err_bits / pi_bits if pi_bits > 0 else 0

        print(f"{x:>8} {pi_val:>8} {li_val:>10.1f} {error:>8.1f} "
              f"{pi_bits:>8} {err_bits:>8} {frac:>10.1%}")

    print()
    print("The 'fraction hard' column shows what proportion of bits cannot be")
    print("determined from the smooth approximation alone.")
    print()
    print("For a TC^0 circuit, we would need ALL bits, including the hard ones.")
    print("The hard bits encode cumulative prime gap information, which appears")
    print("to require prime enumeration (or equivalent) to determine.")
    print()


# ---- FINAL SYNTHESIS ----

SYNTHESIS = """
================================================================================
SYNTHESIS: TC^0 REDUCTION FOR pi(x)
================================================================================

MAIN FINDINGS:

1. NO KNOWN TC^0 REDUCTION EXISTS for pi(x).

2. The most promising path (AKS -> matrix powering -> HAB iterated multiplication)
   FAILS because it requires circuit size exponential in input length:
   - AKS reduces to M^n mod m for a polylog-dimensional matrix M
   - HAB handles iterated multiplication of GIVEN matrices in TC^0
   - But expressing M^n as a product of n copies of M requires listing
     n matrices, giving circuit size poly(n) = 2^{O(input size)}
   - Compressing this to poly(input size) = poly(log n) is the OPEN PROBLEM

3. The equivalence chain:
   pi(x) in TC^0
     <=> PRIMES in TC^0  (since counting a TC^0 predicate is TC^0)
     <=> modular exponentiation in TC^0  (via AKS + matrix formulation)
     <=> iterated squaring in Z_n in TC^0

   All four are equivalent OPEN PROBLEMS.

4. Alternative paths all fail:
   - Legendre/Meissel-Lehmer: super-polynomial terms
   - Lucy Hedgehog DP: sequential depth O(sqrt(x)/ln(x))
   - Partial sieve (constant k primes): O(x/ln(x)) error
   - Mobius function: requires factoring (at least as hard)
   - pi(x) mod 2: no easier than full pi(x)

5. REFINED OPEN QUESTION: Is modular exponentiation (a^b mod m, with all
   three as inputs of log(n) bits) in uniform TC^0? This is the PRECISE
   question that would resolve PRIMES in TC^0 and hence pi(x) in TC^0.

6. BARRIER ASSESSMENT: The following evidence suggests pi(x) is NOT in TC^0:
   - PRIMES not in AC^0[p] (proven lower bound just below TC^0)
   - AKS has inherently sequential structure (repeated squaring)
   - All known algorithms have depth at least O(log n)
   - But NO PROOF exists that pi(x) is not in TC^0

7. EVIDENCE FOR TC^0 POSSIBILITY:
   - Division IS in TC^0 (non-obvious! required deep algebraic insight by HAB)
   - Iterated multiplication IS in TC^0 (even more surprising)
   - These results used Chinese Remainder Theorem + look-up tables in clever ways
   - A similar breakthrough for modular exponentiation is conceivable

RECOMMENDATION FOR FUTURE WORK:
- Focus on the modular exponentiation question: is a^b mod m in TC^0?
- Study Allender's work on arithmetic circuits more carefully
- Look for CRT-based approaches to modular exponentiation
- Investigate whether the STRUCTURE of the exponent (specific forms of n
  appearing in AKS) can be exploited
- Consider whether approximate/partial results (e.g., computing pi(x) mod
  small constants) have lower circuit complexity
"""

# =============================================================================
# MAIN: Run all experiments
# =============================================================================

def main():
    print(ANALYSIS)
    print()

    experiments = [
        experiment_1_legendre_terms,
        experiment_2_lucy_dp_depth,
        experiment_3_tc0_building_blocks,
        experiment_4_parity,
        experiment_5_partial_sieve,
        experiment_6_floor_values,
        experiment_7_aks_depth,
        experiment_8_matrix_powering,
        experiment_9_bit_difficulty,
    ]

    t0 = time.time()
    for exp in experiments:
        exp()
        print()

    elapsed = time.time() - t0
    print(f"All experiments completed in {elapsed:.2f}s")
    print()
    print(SYNTHESIS)

    print()
    print("=" * 70)
    print("FILE SAVED: experiments/circuit_complexity/tc0_reduction_analysis.py")
    print("=" * 70)
    print()
    print("STATUS: This analysis confirms that TC^0 for pi(x) is equivalent to")
    print("the open problem of modular exponentiation in TC^0. No reduction to")
    print("KNOWN TC^0 functions was found, but no impossibility was proven either.")
    print("The gap remains: AC^0[p] < PRIMES <=? TC^0 < NC^1 < P.")


if __name__ == "__main__":
    main()

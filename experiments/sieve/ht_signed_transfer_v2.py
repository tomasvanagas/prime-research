#!/usr/bin/env python3
"""
H-T Signed Cancellation Transfer to pi(x) -- DEEP INVESTIGATION (v2)
=====================================================================
Session 16 experiment.

PREVIOUS RESULT (Session 11, ht_transfer_attempt.py):
  Direct H-T transfer FAILS because pi(x) counts primes (positive, no sign).
  H-T exploits (-1)^{omega(n)} cancellation in mu(n).

THIS EXPERIMENT: Can H-T's technique apply INDIRECTLY?

Five attack vectors:
  1. Explicit formula: R(x^rho) terms oscillate -- apply H-T to the zero sum?
  2. Weighted prime counting: w(n)*Lambda(n) with cancellation, then recover pi(x)?
  3. Buchstab + signed weights: introduce signs, compute in O(x^{3/5}), undo?
  4. M(x) -> pi(x) conversion cost analysis
  5. Concrete numerical experiments for small x

Author: Claude (Session 16)
Date: 2026-04-04
"""

import math
import time
from functools import lru_cache
from collections import defaultdict

# =============================================================================
# PART 0: Utilities
# =============================================================================

def sieve_primes(limit):
    """Sieve of Eratosthenes. Returns boolean array and list of primes."""
    is_prime = [False, False] + [True] * (limit - 1)
    for i in range(2, int(limit**0.5) + 1):
        if is_prime[i]:
            for j in range(i*i, limit + 1, i):
                is_prime[j] = False
    primes = [i for i in range(2, limit + 1) if is_prime[i]]
    return is_prime, primes

def mobius_sieve(limit):
    """Compute mu(n) for n = 0..limit."""
    mu = [0] * (limit + 1)
    mu[1] = 1
    spf = list(range(limit + 1))
    for i in range(2, int(limit**0.5) + 1):
        if spf[i] == i:
            for j in range(i*i, limit + 1, i):
                if spf[j] == j:
                    spf[j] = i
    for n in range(2, limit + 1):
        if spf[n] == n:
            mu[n] = -1
        else:
            p = spf[n]
            if (n // p) % p == 0:
                mu[n] = 0  # p^2 | n
            else:
                mu[n] = -mu[n // p]
    return mu

def von_mangoldt(n):
    """Compute Lambda(n). Returns log(p) if n=p^k, else 0."""
    if n <= 1:
        return 0.0
    d = 2
    while d * d <= n:
        if n % d == 0:
            # d is smallest prime factor
            m = n
            while m % d == 0:
                m //= d
            if m == 1:
                return math.log(d)
            else:
                return 0.0
        d += 1
    return math.log(n)  # n is prime

def big_omega(n):
    """Omega(n) = total number of prime factors with multiplicity."""
    if n <= 1:
        return 0
    count = 0
    d = 2
    while d * d <= n:
        while n % d == 0:
            count += 1
            n //= d
        d += 1
    if n > 1:
        count += 1
    return count

def little_omega(n):
    """omega(n) = number of distinct prime factors."""
    if n <= 1:
        return 0
    count = 0
    d = 2
    while d * d <= n:
        if n % d == 0:
            count += 1
            while n % d == 0:
                n //= d
        d += 1
    if n > 1:
        count += 1
    return count


# =============================================================================
# PART 1: Explicit formula approach -- oscillatory zero sum
# =============================================================================

def analyze_explicit_formula_cancellation():
    """
    APPROACH 1: The explicit formula for pi(x) involves R(x^rho) summed
    over zeta zeros rho = 1/2 + i*gamma. Each term oscillates as
    x^{i*gamma}/rho, creating cancellation. Can H-T apply here?

    The H-T technique works on SUMS OVER INTEGERS (Dirichlet series).
    The zeta zero sum is NOT a sum over integers -- it's a sum over zeros.

    Analysis:
    - H-T exploits: combinatorial structure of squarefree numbers
    - Zero sum has: analytic structure (GUE distribution of spacings)
    - These are FUNDAMENTALLY DIFFERENT types of sums

    The zero sum R(x^rho) ~ li(x^rho) = integral from 0 to x^rho dt/ln(t)
    For rho = 1/2 + i*gamma: x^rho = x^{1/2} * e^{i*gamma*ln(x)}

    The sum is: sum_rho x^{1/2} * e^{i*gamma*ln(x)} / rho
              ~ x^{1/2} * sum_n c_n * e^{i*gamma_n * ln(x)}

    This is an EXPONENTIAL SUM with frequencies gamma_n.
    Cancellation in exponential sums is studied via van der Corput, etc.
    But the cancellation is ALREADY exploited in the explicit formula --
    it gives the O(x^{1/2+eps}) bound. The issue is that we need the
    EXACT value, not just a bound on the tail.

    H-T's contribution for M(x) was reducing O(x^{2/3}) to O(x^{3/5})
    in the COMBINATORIAL (sieve) approach. The analytic approach for M(x)
    is already O(x^{1/2+eps}) using zeta zeros.

    So the question becomes: can H-T's combinatorial trick improve the
    COMBINATORIAL approach for pi(x) from O(x^{2/3}) to O(x^{3/5})?
    Session 11 showed: NO, because pi(x) lacks signed cancellation.

    NEW QUESTION: Can we compute some SIGNED function F(x) in O(x^{3/5})
    that, combined with easily-computed terms, gives pi(x)?
    """
    print("=" * 72)
    print("APPROACH 1: Explicit formula / oscillatory zero sum")
    print("=" * 72)

    print("""
  The explicit formula: pi(x) = R(x) - sum_rho R(x^rho) - corrections

  H-T works on sums over INTEGERS with multiplicative structure.
  The zero sum is over ANALYTIC objects (zeta zeros).
  These are fundamentally different -- H-T does NOT apply to the zero sum.

  However, the Riemann function R(x) = sum_{k=1}^inf mu(k)/k * li(x^{1/k})
  involves mu(k). And M(x) = sum mu(n) can be computed in O(x^{3/5}) by H-T.

  Question: Does computing R(x) benefit from fast M(x)?

  R(x) = li(x) - li(x^{1/2})/2 - li(x^{1/3})/3 + li(x^{1/5})/5 + ...
  The k-th term is mu(k)/k * li(x^{1/k}), which needs mu(k) for k up to log2(x).
  These mu(k) values are trivial to compute (k <= ~300 for x = 10^100).
  So R(x) is already O(polylog) -- no help needed from H-T.

  The bottleneck is the ZERO SUM, not R(x). H-T cannot help with the zero sum.

  VERDICT: FAILS. H-T technique does not apply to sums over zeta zeros.
  The combinatorial and analytic approaches are structurally incompatible.
""")
    return "FAIL"


# =============================================================================
# PART 2: Weighted prime counting with cancellation
# =============================================================================

def analyze_weighted_prime_counting(x_values):
    """
    APPROACH 2: Compute W(x) = sum_{n<=x} w(n) * f(n) where w gives
    cancellation, then extract pi(x) from W(x).

    Several candidate weighting schemes:
    (a) psi(x) = sum Lambda(n): already O(x^{2/3}) via Lucy DP
    (b) sum mu(n)*Lambda(n) = sum over primes (-log p) + corrections
    (c) sum (-1)^Omega(n) * 1_prime(n) = -pi(x) (trivial signs)
    (d) sum mu(n)*ln(n) = -sum Lambda(d)*mu(n/d)... Selberg symmetry formula
    (e) Liouville weighted: sum lambda(n) * (1 if n prime else 0)

    The KEY problem: to get pi(x) from any signed sum, we need to UNDO
    the weighting. If w(n) = mu(n) on primes, we just get -pi(x) back.
    For w(n) to create cancellation, it must sometimes be +1 and sometimes -1
    on primes. But then undoing requires knowing WHICH primes got which sign.
    """
    print("=" * 72)
    print("APPROACH 2: Weighted prime counting with cancellation")
    print("=" * 72)

    for x in x_values:
        is_prime, primes = sieve_primes(x)
        mu = mobius_sieve(x)

        print(f"\n  x = {x}, pi(x) = {len(primes)}")

        # (a) psi(x) = sum Lambda(n)
        psi = sum(von_mangoldt(n) for n in range(1, x + 1))
        print(f"  psi(x) = {psi:.4f}, pi(x)*avg_log ~ {len(primes) * math.log(x/2):.4f}")

        # (b) sum mu(n) * Lambda(n) over n <= x
        # mu(p^k) = -1 if k=1, 0 if k>1. Lambda(p^k) = log(p).
        # So mu(n)*Lambda(n) = -log(p) if n=p (prime), 0 otherwise.
        sum_mu_lambda = sum(mu[n] * von_mangoldt(n) for n in range(1, x + 1))
        minus_sum_log_p = -sum(math.log(p) for p in primes)
        print(f"  sum mu(n)*Lambda(n) = {sum_mu_lambda:.4f}")
        print(f"  -sum log(p) = {minus_sum_log_p:.4f}  (should match)")

        # (c) Selberg's formula: sum_{n<=x} Lambda(n)^2 + sum_{n<=x} Lambda(n)*sum_{d|n} Lambda(d)
        # Related to pair correlation of zeros. But computing this is O(x*log(x)).

        # (d) Liouville-weighted prime count: sum lambda(n) * 1_prime(n) = -pi(x)
        # Because lambda(p) = (-1)^1 = -1 for primes.
        liouville_prime = sum((-1)**big_omega(n) for n in range(2, x+1) if is_prime[n])
        print(f"  sum lambda(p) for p prime = {liouville_prime} = -{len(primes)}")

        # (e) Create a signed function that INCLUDES pi(x) information:
        # F(x) = sum_{n<=x} mu(n) * [n is squarefree with omega(n)=1]
        # This is sum over primes of mu(p) = sum over primes of (-1) = -pi(x). Trivial.

        # (f) Chebyshev-like: sum_{n<=x} (-1)^n * 1_prime(n)
        # This alternates sign on primes. But the sum is NOT multiplicative.
        alt_prime = sum((-1)**p for p in primes)
        print(f"  sum (-1)^p for p prime = {alt_prime}")
        # This is pi_even(x) - pi_odd(x) where pi_even counts even primes (just 2).
        # Equals 1 - (pi(x) - 1) = 2 - pi(x) for x >= 2. Not useful.

        # (g) More interesting: sum_{n<=x} mu(n) * floor(x/n) = 1 (Mertens identity)
        # Can we create a weighted sum involving primes that telescopes?

        # Vaughan's identity: Lambda(n) = Lambda_1(n) - mu_1 * log * 1(n) + ...
        # This decomposes Lambda into Type I and Type II sums.
        # Type I: sum f(d) * g(n/d) where f is "short" (d <= y)
        # Type II: sum a(m) * b(n) where m,n ~ sqrt(x)
        # These are used for BOUNDS on prime sums, not exact computation.

    print("""
  THEORETICAL ANALYSIS:
  ---------------------
  Any weighting w(n) on primes that creates cancellation must satisfy:
    1. w(p) takes both positive and negative values on primes p
    2. sum_{p<=x} w(p) can be computed in O(x^{3/5}) via H-T-like technique
    3. pi(x) can be recovered from {sum_{p<=x} w(p)} efficiently

  Condition 3 is the killer. To recover pi(x) from a signed prime sum,
  we need to know WHICH primes got positive vs negative weight.

  If w(p) depends only on p (not on other primes), then recovering pi(x)
  from sum w(p) is an inversion problem. Possible cases:

  Case w(p) = (-1)^p: sum = 2 - pi(x). Trivially gives pi(x) but the
    sum is NOT a multiplicative function -- can't use H-T.

  Case w(p) = mu(p) = -1: sum = -pi(x). Trivial, gives pi(x) directly.

  Case w(p) = chi(p) for Dirichlet character chi:
    sum chi(p) is computable but does NOT give pi(x) alone.
    Need ALL characters mod q, then sum gives pi(x; q, a) for each residue class.
    Total cost = phi(q) * O(x^{2/3}) per character. WORSE.

  Case w(n) = mu(n) (on all n, not just primes):
    M(x) = sum mu(n). Computable in O(x^{3/5}) by H-T.
    But M(x) does NOT determine pi(x) without additional information.
    The relationship pi(x) = sum_{k>=1} mu(k)/k * pi(x^{1/k}) is circular.

  VERDICT: No weighting scheme creates cancellation AND allows pi(x) recovery.
  The fundamental issue: primes are a SUBSET of integers, and any multiplicative
  weight collapses all primes to the same value (mu(p)=-1, lambda(p)=-1).
""")
    return "FAIL"


# =============================================================================
# PART 3: Buchstab's identity with signed weights
# =============================================================================

def analyze_buchstab_signed(x_values):
    """
    APPROACH 3: Buchstab's identity decomposes pi(x) recursively.
    Can we introduce signs via Mobius inversion to create cancellation?

    Buchstab's identity:
    S(x, y) = #{n <= x : all prime factors of n > y}
    S(x, y) = S(x, z) - sum_{y < p <= z} S(x/p, p)   for y < z

    For pi(x): pi(x) = pi(sqrt(x)) - 1 + S(x, sqrt(x))
    where S(x, sqrt(x)) counts n <= x with smallest prime factor > sqrt(x),
    i.e., primes > sqrt(x).

    Signed Buchstab: What if we compute sum_{n<=x, lpf(n)>y} mu(n)?
    This IS what H-T compute. The signed version of S(x, y) is:
    S_mu(x, y) = sum_{n<=x, lpf(n)>y} mu(n)

    Then M(x) = S_mu(x, 1) and H-T compute this in O(x^{3/5}).

    Can we relate S_mu(x, y) to S(x, y) (the unsigned version)?

    S(x, y) = #{n<=x : lpf(n)>y} = count of y-smooth-free numbers
    S_mu(x, y) = sum_{n<=x, lpf(n)>y} mu(n) = signed count

    The difference is:
    S(x, y) - S_mu(x, y) = 2 * #{n<=x : lpf(n)>y, omega(n) odd, squarefree}
                          + #{n<=x : lpf(n)>y, not squarefree}

    This difference is NOT simpler than S(x, y) itself.
    """
    print("=" * 72)
    print("APPROACH 3: Buchstab's identity with signed weights")
    print("=" * 72)

    for x in x_values:
        is_prime, primes = sieve_primes(x)
        mu = mobius_sieve(x)
        sqrt_x = int(math.sqrt(x))

        print(f"\n  x = {x}, sqrt(x) = {sqrt_x}")

        # Compute S(x, y) and S_mu(x, y) for y = sqrt(x)
        y = sqrt_x

        S_unsigned = 0  # count of n <= x with lpf(n) > y
        S_signed = 0    # sum mu(n) for n <= x with lpf(n) > y
        for n in range(2, x + 1):
            # find smallest prime factor
            d = 2
            spf = n
            while d * d <= n:
                if n % d == 0:
                    spf = d
                    break
                d += 1
            if spf > y:
                S_unsigned += 1
                S_signed += mu[n]

        # S_unsigned should be pi(x) - pi(y) + 1 (the +1 is because n=1 would
        # have "infinite" lpf, but we start from n=2)
        # Actually S_unsigned for n>=2 with lpf>y counts primes in (y, x]
        # (since if n has lpf > y and n <= x, and n composite, then n = p*q
        # with p > y, so n > y^2 = x, contradiction for n <= x).
        # Wait, n could still be <= x if it's prime.
        # Actually if n <= x and lpf(n) > sqrt(x), then n is prime.
        # Because if n = a*b with a,b >= 2 and min(a,b) > sqrt(x), then n > x.
        pi_range = len([p for p in primes if p > y])

        print(f"  S_unsigned(x, sqrt(x)) = {S_unsigned}")
        print(f"  Primes in ({y}, {x}] = {pi_range}")
        print(f"  These should be equal: {S_unsigned == pi_range}")
        print(f"  S_signed(x, sqrt(x)) = {S_signed}")
        print(f"  Difference = {S_unsigned - S_signed}")

        # Now try y = x^{1/5}
        y5 = max(2, int(x ** 0.2))
        S_unsigned_5 = 0
        S_signed_5 = 0
        unsigned_by_omega = defaultdict(int)
        signed_by_omega = defaultdict(int)

        for n in range(2, x + 1):
            d = 2
            spf = n
            while d * d <= n:
                if n % d == 0:
                    spf = d
                    break
                d += 1
            if spf > y5:
                om = little_omega(n)
                S_unsigned_5 += 1
                S_signed_5 += mu[n]
                unsigned_by_omega[om] += 1
                signed_by_omega[om] += mu[n]

        print(f"\n  With y = x^(1/5) = {y5}:")
        print(f"  S_unsigned(x, y) = {S_unsigned_5}")
        print(f"  S_signed(x, y) = {S_signed_5}")
        print(f"  Breakdown by omega:")
        for k in sorted(unsigned_by_omega.keys()):
            print(f"    omega={k}: unsigned={unsigned_by_omega[k]}, "
                  f"signed={signed_by_omega[k]}")

        # The idea: H-T computes S_signed(x, y) for y ~ x^{1/5} in O(x^{3/5}).
        # We want S_unsigned_1 = A_1 = #{primes in (y, x]}.
        # From the breakdown: A_1 = unsigned_by_omega[1]
        # And the signed version: signed_by_omega[1] = -A_1 (since mu(p)=-1)
        # S_signed = -A_1 + A_2 - A_3 + A_4 - A_5
        # S_unsigned = A_1 + A_2 + A_3 + A_4 + A_5

        # From these two equations:
        # A_1 = (S_unsigned - S_signed) / 2 - (A_3 + A_5 - A_3 - A_5 cancels...)
        # Actually: S_unsigned + S_signed = 2*(A_2 + A_4) (even omega terms)
        #           S_unsigned - S_signed = 2*(A_1 + A_3 + A_5) (odd omega terms)

        # So A_1 = (S_unsigned - S_signed)/2 - A_3 - A_5
        # We need A_3 and A_5 separately!

        A_odd = (S_unsigned_5 - S_signed_5) // 2
        A_1 = unsigned_by_omega.get(1, 0)
        A_3 = unsigned_by_omega.get(3, 0)
        A_5 = unsigned_by_omega.get(5, 0)

        print(f"\n  (S_unsigned - S_signed)/2 = {A_odd}")
        print(f"  = A_1 + A_3 + A_5 = {A_1} + {A_3} + {A_5} = {A_1+A_3+A_5}")
        print(f"  So: A_1 = {A_odd} - A_3 - A_5 = {A_odd} - {A_3} - {A_5} = {A_odd - A_3 - A_5}")

        # Computing A_3(x, y) = #{n <= x : n = pqr, p>q>r>y}
        # This is sum_{r > y} sum_{q > r} (pi(x/(qr)) - pi(q))
        # which costs O(x^{2/3}) via Lucy DP.
        # So: even if S_signed is O(x^{3/5}), we still need A_3 at O(x^{2/3}).

        print(f"\n  PROBLEM: A_3(x, y) requires sum_{{r>y}} sum_{{q>r}} pi(x/(qr))")
        print(f"  Each pi evaluation is O(1) with Lucy DP, but there are O(x^{2/3}) of them.")
        print(f"  So recovering A_1 from S_signed costs O(x^{{2/3}}) -- no improvement!")

    print("""
  THEORETICAL ANALYSIS:
  ---------------------
  Buchstab signed approach:
    S_signed(x, y) = sum_{lpf(n)>y} mu(n) -- computable in O(x^{3/5}) by H-T
    S_unsigned(x, y) = sum_{lpf(n)>y} 1 -- this is what we want for pi(x)

  From S_signed and S_unsigned:
    A_1 (primes) = (S_unsigned - S_signed)/2 - A_3 - A_5

  But computing A_3, A_5 costs O(x^{2/3}) anyway (it's Meissel-Lehmer).

  Could we compute BOTH S_signed AND S_unsigned in O(x^{3/5})?
  No! S_unsigned is literally pi(x) minus small primes. If we could
  compute S_unsigned in O(x^{3/5}), we'd have pi(x) in O(x^{3/5}).

  The fundamental asymmetry:
  - M(x) = S_signed(x, 1) has CANCELLATION: the alternating sum A_0-A_1+A_2-...
    allows H-T to bound the contribution of each A_k and show total = O(x^{3/5}).
  - pi(x) = pi(y) + A_1(x, y) requires A_1 ALONE, which is the single largest
    term with NO cancellation partner.

  VERDICT: FAILS. The signed Buchstab version IS what H-T computes (M(x)),
  but undoing the signs to get pi(x) costs O(x^{2/3}) for the A_3, A_5 terms.
""")
    return "FAIL"


# =============================================================================
# PART 4: M(x) -> pi(x) conversion cost
# =============================================================================

def analyze_mx_to_pix_conversion(x_values):
    """
    APPROACH 4: Compute M(x) in O(x^{3/5}) via H-T, then convert to pi(x).

    What identities connect M(x) and pi(x)?

    Identity 1: sum_{d=1}^{x} M(floor(x/d)) = 1
      This is the Mobius inversion of sum_{d|n} mu(d) = [n=1].
      It determines M at all floor(x/d) values from scratch, but doesn't help pi.

    Identity 2: pi(x) = sum_{k=1}^{log2(x)} mu(k)/k * Li(x^{1/k}) + corrections
      The corrections involve the zeta zero sum. Li terms are easy.
      The corrections cost O(x^{1/2+eps}) analytically or O(x^{2/3}) combinatorially.
      Knowing M(x) doesn't help with these corrections.

    Identity 3: pi(x) = R(x) - sum_rho R(x^rho) + integral term
      R(x) uses mu(k) for k <= log2(x) -- trivial, doesn't need M(x).
      The zero sum is the bottleneck.

    Identity 4: psi(x) = x - sum_{rho} x^rho/rho - ln(2pi)
      pi(x) = psi(x)/ln(x) + integral correction (partial summation)
      Getting pi(x) from psi(x) requires partial summation, which is
      itself a sum over prime powers. Not free.

    Identity 5: pi(x) = psi(x)/ln(x) + sum_{k=2}^{log_2(x)} pi(x^{1/k})/k
      Rewriting: pi(x) = psi(x)/ln(x) + pi(x^{1/2})/2 + pi(x^{1/3})/3 + ...
      The pi(x^{1/k}) terms for k>=2 are small (pi(x^{1/2}) ~ 2*sqrt(x)/ln(x))
      and can be computed recursively. But psi(x) itself costs O(x^{2/3}).

    KEY QUESTION: Can M(x) values help compute psi(x)?

    psi(x) = -sum_{n<=x} mu(n) * ln(n) * ... no, that's not right.

    Actually: sum_{n<=x} Lambda(n) = psi(x)
    And: Lambda = -mu * log via Dirichlet convolution: Lambda(n) = -sum_{d|n} mu(d)*log(d)

    So: psi(x) = -sum_{n<=x} sum_{d|n} mu(d)*log(d)
               = -sum_{d<=x} mu(d)*log(d) * floor(x/d)

    If we know mu(d) for all d <= x (from the sieve, O(x) time), this is computable
    in O(x) time. But that's WORSE than O(x^{2/3}).

    Using floor-value compression: the sum -sum mu(d)*log(d)*floor(x/d) has
    O(sqrt(x)) distinct values of floor(x/d). So we need sum mu(d)*log(d)
    for d in blocks. This requires partial sums of mu(d)*log(d), which is
    related to the derivative of 1/zeta(s) at s=1... itself complex.

    ACTUALLY, the key identity for converting M(x) -> psi(x):
    M(x) = sum mu(n) for n<=x
    M_1(x) = sum mu(n)*log(n) for n<=x
    psi(x) = sum Lambda(n) = -M_1(x) + sum_{n<=x} mu(n) * (x/n - floor(x/n) corrections)

    No, let's be precise:
    psi(x) = -sum_{d<=x} mu(d) * log(d) * floor(x/d)
           = -sum_{d<=x} mu(d) * log(d) * (x/d - {x/d})
           = -x * sum mu(d)*log(d)/d + sum mu(d)*log(d)*{x/d}
           ~ x * (1/zeta)'(1) / (1/zeta(1)) ... diverges

    This is getting circular. Let me just compute numerically.
    """
    print("=" * 72)
    print("APPROACH 4: M(x) -> pi(x) conversion cost analysis")
    print("=" * 72)

    for x in x_values:
        is_prime, primes = sieve_primes(x)
        mu = mobius_sieve(x)
        pi_x = len(primes)
        M_x = sum(mu[n] for n in range(1, x + 1))

        print(f"\n  x = {x}: pi(x) = {pi_x}, M(x) = {M_x}")

        # psi(x) via Lambda
        psi_x = sum(von_mangoldt(n) for n in range(1, x + 1))
        print(f"  psi(x) = {psi_x:.4f}")

        # psi(x) via Mobius: psi(x) = -sum mu(d)*log(d)*floor(x/d)
        psi_via_mu = -sum(mu[d] * math.log(d) * (x // d) for d in range(1, x + 1))
        print(f"  psi(x) via -sum mu(d)*log(d)*floor(x/d) = {psi_via_mu:.4f}")

        # Number of distinct floor(x/d) values
        floor_vals = set()
        for d in range(1, x + 1):
            floor_vals.add(x // d)
        print(f"  Distinct floor(x/d) values: {len(floor_vals)} vs sqrt(x)={int(x**0.5)}")

        # M_1(x) = sum mu(n)*log(n)
        M1_x = sum(mu[n] * math.log(n) for n in range(1, x + 1))
        print(f"  M_1(x) = sum mu(n)*log(n) = {M1_x:.4f}")

        # Check: psi(x) ~ -M_1(x) + error terms?
        print(f"  -M_1(x) = {-M1_x:.4f}")
        print(f"  psi(x) + M_1(x) = {psi_x + M1_x:.4f}")
        print(f"  (difference is NOT small -- no simple relation)")

        # The correct identity: psi(x) = sum_{d<=x} mu(d) * sum_{k=1}^{floor(x/d)} log(k)
        # = sum mu(d) * log(floor(x/d)!) -- Stirling gives log(n!) ~ n*log(n) - n
        psi_via_factorial = sum(mu[d] * sum(math.log(k) for k in range(1, x//d + 1))
                               for d in range(1, x + 1))
        print(f"  psi(x) via sum mu(d)*log(floor(x/d)!) = {psi_via_factorial:.4f}")

        # This is the Chebyshev identity: sum_{n<=x} Lambda(n) = sum_{d<=x} mu(d)*log(floor(x/d)!)
        # Since log(n!) = sum_{k<=n} log(k), and sum_{d|m} mu(d) = [m=1].

    print("""
  CONVERSION COST ANALYSIS:
  -------------------------
  To go from M(x) -> pi(x), the cheapest known routes are:

  Route A: M(x) -> (analytic continuation) -> pi(x)
    Cost: O(x^{1/2+eps}) -- the analytic approach using zeta zeros.
    M(x) doesn't help here; both M(x) and pi(x) require the zero sum.

  Route B: M(x) -> psi(x) -> pi(x)
    psi(x) = -sum mu(d)*log(d)*floor(x/d)  -- needs O(sqrt(x)) distinct
    partial sums of mu(d)*log(d), each requiring O(x^{2/3}) work.
    Total: O(x^{2/3} * sqrt(x)) = O(x^{7/6}). WORSE.

    Alternatively, psi(x) = sum mu(d)*log(floor(x/d)!)
    This has O(sqrt(x)) distinct floor values, and log(n!) is O(1) per
    Stirling. But we still need M(x/d) for grouped summation -> O(x^{2/3}).

  Route C: M(x) at MANY points -> pi(x)
    The identity sum_{d=1}^{x} M(floor(x/d)) = 1 is about M itself.
    No known identity gives pi(x) from M values at multiple points
    without additional O(x^{2/3}) cost terms.

  VERDICT: The conversion M(x) -> pi(x) costs at LEAST O(x^{1/2+eps})
  (analytic) or O(x^{2/3}) (combinatorial). Even with M(x) in O(x^{3/5}),
  the conversion is the bottleneck. No improvement for pi(x).
""")
    return "FAIL"


# =============================================================================
# PART 5: Concrete numerical experiments
# =============================================================================

def concrete_experiments(x_values):
    """
    For small x, compute various weighted sums and check:
    (a) Which sums have cancellation (small magnitude relative to length)?
    (b) Which sums determine pi(x)?
    (c) Is there any sum that has BOTH properties?
    """
    print("=" * 72)
    print("APPROACH 5: Concrete numerical experiments")
    print("=" * 72)

    for x in x_values:
        is_prime, primes = sieve_primes(x)
        mu = mobius_sieve(x)
        pi_x = len(primes)

        print(f"\n{'~'*60}")
        print(f"  x = {x}, pi(x) = {pi_x}")
        print(f"{'~'*60}")

        # Collect various summatory functions
        sums = {}

        # 1. M(x) = sum mu(n) -- has cancellation, doesn't determine pi(x) alone
        sums['M(x)'] = sum(mu[n] for n in range(1, x + 1))

        # 2. L(x) = sum lambda(n) -- Liouville, has cancellation
        sums['L(x)'] = sum((-1)**big_omega(n) for n in range(1, x + 1))

        # 3. psi(x) = sum Lambda(n) -- ~x, no cancellation
        sums['psi(x)'] = sum(von_mangoldt(n) for n in range(1, x + 1))

        # 4. sum mu(n)*Lambda(n) = -sum_{p<=x} log(p) -- moderate cancellation
        sums['sum_mu_Lambda'] = sum(mu[n] * von_mangoldt(n) for n in range(1, x + 1))

        # 5. theta(x) = sum_{p<=x} log(p) -- Chebyshev theta, ~x
        sums['theta(x)'] = sum(math.log(p) for p in primes)

        # 6. sum_{n<=x} mu(n) * omega(n) -- signed count weighted by omega
        sums['sum_mu_omega'] = sum(mu[n] * little_omega(n) for n in range(1, x + 1))

        # 7. sum_{n<=x} (-1)^{omega(n)} * omega(n) -- another weighted sum
        sums['sum_alt_omega'] = sum((-1)**little_omega(n) * little_omega(n)
                                    for n in range(1, x + 1))

        # 8. M_1(x) = sum mu(n)*log(n) -- derivative-like, has cancellation
        sums['M_1(x)'] = sum(mu[n] * math.log(n) for n in range(1, x + 1))

        # 9. sum_{n<=x} mu(n) * pi(n) -- convolution of mu and pi
        cumpi = [0] * (x + 1)
        for n in range(2, x + 1):
            cumpi[n] = cumpi[n-1] + (1 if is_prime[n] else 0)
        sums['sum_mu_pi'] = sum(mu[n] * cumpi[n] for n in range(1, x + 1))

        # 10. sum_{d<=x} mu(d) * pi(floor(x/d)) -- Mobius inversion
        sums['sum_mu_pi_xd'] = sum(mu[d] * cumpi[x // d] for d in range(1, x + 1))

        # Identity: sum_{d<=x} mu(d) * F(floor(x/d)) inverts sum_{d<=x} f(d)*floor(x/d) = G(x)
        # If G(x) = sum_{n<=x} 1_{prime}(n) * floor(x/n)... not standard.

        # Actually: sum_{d|n} f(d) = g(n) implies f(n) = sum_{d|n} mu(n/d)*g(d)
        # For Dirichlet series: if G = f*1 then f = G*mu (convolution).
        # pi(x) = sum 1_prime(n) for n<=x.
        # sum_{d<=x} mu(d) * Pi(floor(x/d)) where Pi(y) = sum_{p^k <= y} 1/k
        # inverts Pi(x) = sum_{k>=1} pi(x^{1/k})/k (prime power counting).

        # 11. What about Ramanujan sums? c_q(n) = sum_{(a,q)=1} e^{2pi i an/q}
        # sum c_q(n) over n has structure but doesn't directly count primes.

        print(f"\n  Summatory functions:")
        print(f"  {'Function':<20s} {'Value':>12s} {'|Value|/x':>10s} {'Cancellation':>14s}")
        print(f"  {'-'*56}")
        for name, val in sums.items():
            ratio = abs(val) / x if x > 0 else 0
            canc = "HIGH" if ratio < 0.1 else ("MODERATE" if ratio < 0.5 else "LOW")
            if isinstance(val, float):
                print(f"  {name:<20s} {val:>12.2f} {ratio:>10.4f} {canc:>14s}")
            else:
                print(f"  {name:<20s} {val:>12d} {ratio:>10.4f} {canc:>14s}")

        # KEY TEST: Which determine pi(x)?
        # M(x) alone doesn't determine pi(x) -- different x can give same M(x)
        # psi(x) determines pi(x) via partial summation (but costs O(x^{2/3}))
        # theta(x) determines pi(x) (same)

        # Test: sum_{d<=x} mu(d) * pi(floor(x/d))
        # This should be = 1_{is_prime}(x) + 1_{x=1}? Let's check.
        # Mobius inversion: if g(x) = sum_{d<=x} f(d) then f(x) = sum_{d<=x} mu(d)*g(floor(x/d))
        # With g = pi (= summatory function of 1_prime), f = 1_prime.
        # So sum_{d<=x} mu(d)*pi(floor(x/d)) = 1 if x is prime, 0 otherwise.
        print(f"\n  sum_{{d<=x}} mu(d)*pi(floor(x/d)) = {sums['sum_mu_pi_xd']}")
        print(f"  (should be 1 if x is prime, 0 otherwise)")
        print(f"  x = {x} is {'PRIME' if is_prime[x] else 'NOT PRIME'}")

        # This is a POINTWISE result, not a summatory one. Not useful for pi(x).

        # Test: can we recover pi(x) from M(x) at multiple points?
        # We know M(floor(x/d)) for d = 1, 2, ..., x.
        # From the identity sum M(floor(x/d)) = 1, these are constrained.
        # But there's no identity sum_d alpha_d * M(floor(x/d)) = pi(x).

        # Actually, IS there? Let's search for a linear combination.
        # sum_d c_d * M(floor(x/d)) = pi(x)?
        # Taking Dirichlet series: (sum c_d/d^s) * (1/zeta(s)) = (something related to pi)?
        # 1/zeta(s) has Dirichlet series sum mu(n)/n^s.
        # If sum c_d/d^s = zeta(s) * P(s) where P(s) = sum 1/p^s (prime zeta),
        # then sum c_d * M(floor(x/d)) ~ sum_{n<=x} sum_{d|n} c_d * mu(n/d)
        # = sum_{n<=x} [coefficient of n^{-s} in P(s)]
        # = sum_{p<=x} 1 = pi(x). YES!

        # So c_d = sum_{k: d=pk for some prime p} 1, i.e., c_d = omega(d) if d > 1?
        # Wait, P(s) = sum_p p^{-s}, and zeta(s)*P(s) = sum_n (sum_{p|n} 1)/n^s
        # = sum_n omega(n)/n^s. So c_d = omega(d).

        # Therefore: sum_{d<=x} omega(d) * M(floor(x/d)) = pi(x)?
        # Let's verify!

        lhs = sum(little_omega(d) * sum(mu[k] for k in range(1, x // d + 1))
                  for d in range(1, x + 1))
        # Compute M(floor(x/d)) efficiently
        M_cache = {}
        for d in range(1, x + 1):
            v = x // d
            if v not in M_cache:
                M_cache[v] = sum(mu[k] for k in range(1, v + 1))

        lhs_fast = sum(little_omega(d) * M_cache[x // d] for d in range(1, x + 1))
        print(f"\n  IDENTITY TEST: sum omega(d) * M(floor(x/d)) = {lhs_fast}")
        print(f"  pi(x) = {pi_x}")
        print(f"  Match: {lhs_fast == pi_x}")

        # If this works, can we compute it in O(x^{3/5})?
        # sum_{d<=x} omega(d) * M(floor(x/d))
        # has O(sqrt(x)) distinct values of floor(x/d).
        # For each block [d1, d2] with floor(x/d) = v, we need sum_{d=d1}^{d2} omega(d).
        # This requires partial sums of omega(d), which is sum_{p<=d} floor(d/p).
        # Computing partial sums of omega over a range is O(sqrt(range)) per block? No.
        # sum_{d<=y} omega(d) = sum_{p<=y} floor(y/p) ~ y * log(log(y)) + O(y/log(y))
        # Each such sum is O(sqrt(y)) via hyperbola method.
        # Total: O(sqrt(x)) blocks * O(sqrt(block_size)) work per block.
        # This is O(x^{2/3}) total -- same as Lucy DP.

        # The bottleneck is the same: O(sqrt(x)) distinct floor values,
        # each requiring O(?) work for the omega partial sum.
        # And M(v) for each v requires O(v^{3/5}) if using H-T.
        # Total: sum_{v in floor_vals} v^{3/5} ~ integral x^{3/5} ... still O(x^{2/3}).

    print("""
  KEY DISCOVERY: The identity sum omega(d) * M(floor(x/d)) = pi(x) is CORRECT!
  This follows from the Dirichlet series identity:
    (sum omega(n)/n^s) = zeta(s) * P(s)
  where P(s) = sum_p p^{-s}, and Mobius inversion gives the identity.

  However, computing this sum costs O(x^{2/3}) because:
  1. There are O(sqrt(x)) distinct floor values
  2. Each M(v) value costs O(v^{3/5}) via H-T
  3. Partial sums of omega(d) cost O(sqrt(y)) per range
  4. Total: max(sum v^{3/5}, sqrt(x) * sqrt(x)) = O(x^{2/3})

  Even though M(v) is cheaper at O(v^{3/5}), the OMEGA PARTIAL SUMS
  and the NUMBER OF DISTINCT VALUES prevent improvement below O(x^{2/3}).

  This is a concrete instance of the general barrier:
  Any multiplicative identity that converts M(x) -> pi(x) must sum
  over O(sqrt(x)) floor values, and the secondary terms cost O(x^{2/3}).
""")

    return "FAIL"


# =============================================================================
# PART 6: Novel identity search -- can ANY signed sum give pi(x) cheaply?
# =============================================================================

def novel_identity_search(x_values):
    """
    Systematically search for identities of the form:
      pi(x) = sum_{d<=x} f(d) * G(floor(x/d))
    where:
      - G is a summatory function computable in O(x^{3/5})
      - f(d) partial sums are computable cheaply
      - Total cost < O(x^{2/3})

    The ONLY summatory function known to be O(x^{3/5}) is M(x) (by H-T).

    So we need: pi(x) = sum f(d) * M(floor(x/d)).
    Dirichlet series: F(s) * (1/zeta(s)) = "prime counting Dirichlet series"
    i.e., F(s) = zeta(s) * P(s) = sum omega(n)/n^s.
    So f(d) = omega(d). (Confirmed above.)

    Are there OTHER G functions computable in O(x^{3/5})?
    L(x) = sum lambda(n) is related to M(x) via L(x) = sum M(floor(x/k^2)).
    Computing L(x) from M values costs O(sqrt(x)) evaluations of M.
    So L(x) is also O(x^{3/5}) (dominated by M at largest argument).

    Can we use L(x)? pi(x) = sum g(d) * L(floor(x/d))?
    Dirichlet series: G(s) * (zeta(2s)/zeta(s)) = "prime counting"
    zeta(2s)/zeta(s) = sum lambda(n)/n^s.
    So G(s) = zeta(s)/zeta(2s) * P(s) = (sum |mu(n)|/n^s) * P(s)
    = sum_{n squarefree} omega(n) / n^s.
    So g(d) = omega(d) * |mu(d)| (omega on squarefree part only, 0 on non-squarefree).
    Same cost structure.

    What about sum_{d<=x} mu(d) * f(floor(x/d)) for various f?
    If f is the prime counting function itself... circular.
    """
    print("=" * 72)
    print("APPROACH 6: Novel identity search")
    print("=" * 72)

    for x in x_values:
        is_prime, primes = sieve_primes(x)
        mu = mobius_sieve(x)
        pi_x = len(primes)

        print(f"\n  x = {x}, pi(x) = {pi_x}")

        # Verify: L(x) = sum_{k^2 <= x} M(floor(x/k^2))
        M_vals = {}
        for n in range(1, x + 1):
            M_vals[n] = sum(mu[k] for k in range(1, n + 1))

        L_via_M = sum(M_vals.get(x // (k*k), 0) for k in range(1, int(x**0.5) + 1))
        L_direct = sum((-1)**big_omega(n) for n in range(1, x + 1))
        print(f"  L(x) direct = {L_direct}")
        print(f"  L(x) via M = {L_via_M}")
        print(f"  Match: {L_via_M == L_direct}")

        # Identity: pi(x) = sum omega(d) * M(floor(x/d)) [verified above]
        # Alternative via L(x): pi(x) = sum omega(d)*|mu(d)| * L(floor(x/d))
        L_cache = {}
        for d in range(1, x + 1):
            v = x // d
            if v not in L_cache:
                L_cache[v] = sum((-1)**big_omega(k) for k in range(1, v + 1))

        pi_via_L = sum(little_omega(d) * (1 if mu[d] != 0 else 0) * L_cache[x // d]
                       for d in range(1, x + 1))
        print(f"  pi(x) via sum omega(d)*|mu(d)|*L(floor(x/d)) = {pi_via_L}")
        print(f"  Match: {pi_via_L == pi_x}")

        # Both identities give pi(x) but cost O(x^{2/3}) due to floor value structure.

        # WILD IDEA: What if we could find f such that f partial sums are O(1)?
        # Then total = O(sqrt(x)) * O(M evaluation) = O(sqrt(x)) * O(x^{3/5})?
        # But that's O(x^{3/5} * x^{1/2}) = O(x^{11/10}). WORSE.
        # Wait no. The M evaluations are at floor(x/d) which are all <= x.
        # Using H-T, we precompute M(v) for all O(sqrt(x)) floor values.
        # H-T's algorithm for M(v) at a single v costs O(v^{3/5}).
        # For ALL floor values simultaneously, the cost is dominated by M(x) = O(x^{3/5}).
        # (Because floor values < sqrt(x) can be done by sieving in O(sqrt(x)) total.)
        # So: precompute all M(floor(x/d)) in O(x^{3/5}).
        # Then: compute sum f(d) * M(floor(x/d)) where f = omega.
        # This requires sum_{d in block} omega(d) for each block.
        # partial sums of omega: sum_{d<=y} omega(d) = sum_{p<=y} floor(y/p)
        # This itself costs O(y/log(y)) by PNT, or O(sqrt(y)) by hyperbola.

        # Total cost of the convolution:
        # O(sqrt(x)) blocks, each needing sum of omega(d) over an interval.
        # Computing sum_{d=a}^{b} omega(d) = sum_{p<=b} (floor(b/p) - floor((a-1)/p))
        # = sum over primes p <= b. For a full block, this is O(pi(b)) per block.
        # Total: O(sqrt(x)) * O(pi(sqrt(x))) = O(sqrt(x) * sqrt(x)/log(x)) = O(x/log(x))

        # Wait, that can't be right. Let me reconsider.
        # The blocks are: for each distinct v = floor(x/d), the block is
        # {d : floor(x/d) = v}, which has size ~ x/v^2.
        # For the large blocks (v small), d ~ x/v, and block size ~ 1.
        # For the small blocks (v large, v > sqrt(x)), d is small (d < sqrt(x)),
        # each block has size 1.

        # Actually in Lucy DP / Meissel-Lehmer, the standard approach groups:
        # - "Small" d (d <= sqrt(x)): each d is its own block, O(sqrt(x)) terms
        # - "Large" v (v <= sqrt(x)): these v's are addressed directly

        # For the small-d part: sum_{d=1}^{sqrt(x)} omega(d) * M(floor(x/d))
        # Each omega(d) is O(log(log(d))), precomputable by sieve.
        # M(floor(x/d)) is precomputed. So this part is O(sqrt(x)).

        # For the large-d part: group by v = floor(x/d).
        # sum_{v=1}^{sqrt(x)} M(v) * (sum_{d: floor(x/d)=v} omega(d))
        # The inner sum is over d in [x/(v+1)+1, x/v].
        # sum_{d=a}^{b} omega(d) requires O(b-a) work naively,
        # or O(sqrt(b)) work using hyperbola for each prime p.
        # Total: sum_{v=1}^{sqrt(x)} O(x/v - x/(v+1)) = sum O(x/v^2) ~ O(sqrt(x))
        # using harmonic sum.

        # So total naive cost: O(sqrt(x)) for the convolution.
        # PLUS O(x^{3/5}) for precomputing all M(floor(x/d)).

        # TOTAL: O(x^{3/5}) !!!

        # WAIT. Is this right?? Let me double-check.

        print(f"\n  COST ANALYSIS for sum omega(d) * M(floor(x/d)):")
        print(f"  - Precompute M(v) for all v = floor(x/d): O(x^{{3/5}}) via H-T")
        print(f"  - Precompute omega(d) for d <= sqrt(x): O(sqrt(x)) via sieve")
        print(f"  - Small d part (d <= sqrt(x)): O(sqrt(x)) lookups")
        print(f"  - Large d part (group by v <= sqrt(x)): need partial sums of omega")

        # For large d part, we need:
        # For each v in 1..sqrt(x): sum_{d=x/(v+1)+1}^{x/v} omega(d)
        # = sum_{d<=x/v} omega(d) - sum_{d<=x/(v+1)} omega(d)
        # = sum_{p<=x/v} floor(x/v / p) - sum_{p<=x/(v+1)} floor(x/(v+1) / p)

        # sum_{d<=y} omega(d) = sum_{p<=y} floor(y/p)
        # For y = x/v, this sum has pi(y) ~ y/log(y) terms.
        # Computing it naively costs O(y/log(y)).
        # For v=1: y=x, cost O(x/log(x)). TOO EXPENSIVE!

        # Can we compute sum_{d<=y} omega(d) for all y = x/v faster?
        # sum_{d<=y} omega(d) = sum_{p<=y} floor(y/p) = sum_{p prime} floor(y/p)
        # This is the PRIME HARMONIC SUM, and computing it requires knowing
        # pi(y/k) for all k, which is... pi(x) again!

        # More precisely:
        # sum_{p<=y} floor(y/p) = sum_{k=1}^{y} sum_{p<=y, p|k} 1
        #                       = sum_{k=1}^{y} omega(k)  [tautological]
        # Also = sum_{p<=y} y/p - sum_{p<=y} {y/p}
        #       ~ y * log(log(y)) + O(y/log(y))
        # But we need the EXACT value, not an approximation.

        # sum_{p<=y} floor(y/p) = sum_{j=1}^{y} pi(y/j) - ... no.
        # Actually: sum_{p<=y} floor(y/p) = sum_{j=1}^{y} (pi(min(j*floor(y/j), y)) - pi(...))
        # This is getting circular.

        # Let's reconsider. The standard way to compute sum f(d)*g(floor(x/d)):
        # Split at sqrt(x).
        # Small d (d <= sqrt(x)): compute f(d) directly, g(floor(x/d)) from table.
        # Large v (v <= sqrt(x)): compute F(x/v) - F(x/(v+1)) where F(y) = sum_{d<=y} f(d).
        # Then multiply by g(v).

        # For f = omega: F(y) = sum_{d<=y} omega(d) = sum_{p<=y} floor(y/p).
        # To compute F(y) for all y = x/v, v=1..sqrt(x):
        # Each F(x/v) = sum_{p<=x/v} floor(x/v / p).
        # This requires pi(x/v / k) for various k, which needs... pi at many points.
        # CIRCULAR!

        # Unless we can compute sum_{p<=y} floor(y/p) without knowing the primes.
        # We can write it as sum_{n<=y} Lambda(n)/log(n) * floor(y/n) + corrections for prime powers.
        # But Lambda(n)/log(n) = 1/log(p) on prime powers, not integer-valued.

        # Alternative: compute sum_{d<=y} omega(d) via inclusion-exclusion on primes.
        # This IS the Lucy DP approach and costs O(y^{2/3}).
        # For the largest y = x, cost is O(x^{2/3}).

        print(f"  - PROBLEM: sum_{{d<=y}} omega(d) = sum_{{p<=y}} floor(y/p)")
        print(f"    requires knowing primes up to y, i.e., pi() evaluations.")
        print(f"    Computing this for y = x costs O(x^{{2/3}}) via Lucy DP.")
        print(f"    So the convolution sum is NOT cheaper than O(x^{{2/3}}).")

    print("""
  REFINED ANALYSIS:
  -----------------
  The identity pi(x) = sum omega(d) * M(floor(x/d)) is correct but its
  evaluation cost depends on computing partial sums of omega(d).

  Naive approach:
    F(y) = sum_{d<=y} omega(d) = sum_{p prime, p<=y} floor(y/p)
    This requires pi(y/k) for k=1..y, which is O(y^{2/3}) via Lucy DP.

  So even if M(floor(x/d)) is precomputed in O(x^{3/5}) via H-T, the
  omega partial sums F(y) for the grouping step cost O(x^{2/3}).

  The bottleneck SHIFTS from M(x) to the partial sums of omega(d).
  This is because omega(d) = sum_{p|d} 1, which depends on primes.

  Could we choose a DIFFERENT f(d) whose partial sums are cheap?
  We need: f * mu = 1_prime (Dirichlet convolution), so f = 1_prime * 1 = omega
  (sum over prime divisors). There's no freedom in choosing f.

  Alternative decomposition: pi(x) = sum g(d) * L(floor(x/d))
  where g(d) = omega(d) * |mu(d)|. Same partial sum problem.

  VERDICT: FAILS. The identity exists but the omega partial sums cost O(x^{2/3}),
  which is the SAME bottleneck as Meissel-Lehmer. The H-T speedup for M(x) is
  wasted because the conversion step re-introduces the O(x^{2/3}) cost.
""")
    return "FAIL"


# =============================================================================
# PART 7: Ultimate obstacle analysis
# =============================================================================

def obstacle_analysis():
    """
    Why does H-T's O(x^{3/5}) for M(x) NOT extend to pi(x)?
    A precise, rigorous analysis of the structural barrier.
    """
    print("=" * 72)
    print("FINAL ANALYSIS: Why H-T cannot help pi(x)")
    print("=" * 72)

    print("""
  THE HELFGOTT-THOMPSON TECHNIQUE (summary):
  ==========================================
  H-T compute M(x) = sum_{n<=x} mu(n) in O(x^{3/5+eps}).

  Key ingredients:
  1. The identity M(x) = 1 - sum_{n=2}^{x} M(floor(x/n)) [self-referential]
  2. Sieve: factor n into small (p <= z) and large (p > z) parts
  3. For fixed small-prime signature S, count large-factor contributions
  4. CANCELLATION: mu alternates signs, so A_1 - A_2 + A_3 - ... is small
  5. The key savings: with z = x^{1/5}, almost-primes have <= 4 large factors,
     and the alternating sum has cancellation that makes certain cross-terms vanish
  6. Careful balancing of z gives the O(x^{3/5}) bound

  WHY THIS FAILS FOR pi(x):
  ==========================

  BARRIER 1: No self-referential identity with cancellation
  ---------------------------------------------------------
  pi(x) = pi(sqrt(x)) + S(x, sqrt(x)) - 1 [Legendre]
  = pi(sqrt(x)) + Phi(x, pi(sqrt(x))) - P_2 - P_3 - ... [Meissel-Lehmer]

  The P_k correction terms are ALL positive (they count k-almost-primes).
  There is NO alternating sum. Every term must be computed individually.

  In contrast, M(x) = 1 - sum M(floor(x/n)) has cancellation because
  M values oscillate in sign.

  BARRIER 2: Primes have trivial small-prime structure
  ----------------------------------------------------
  In H-T's sieve, numbers are classified by their small-prime divisors.
  A prime p has NO small prime divisors (unless p <= z itself).

  For M(x): a number n with small-prime signature {2,3,5} and k large prime
  factors contributes (-1)^{3+k} * [squarefree]. The signs create cancellation
  across different k values.

  For pi(x): a prime p > z has empty signature and k=1. There's only ONE
  category to count: primes in (z, x]. No cancellation.

  BARRIER 3: Conversion M(x) -> pi(x) is itself O(x^{2/3})
  ----------------------------------------------------------
  The identity pi(x) = sum omega(d) * M(floor(x/d)) connects M and pi.
  But evaluating this convolution requires partial sums of omega(d),
  which in turn require pi evaluations -- circular at O(x^{2/3}).

  BARRIER 4: The signed Buchstab identity gives M(x), not pi(x)
  ---------------------------------------------------------------
  S_signed(x, y) = sum_{lpf(n)>y} mu(n) is what H-T compute (it gives M(x)).
  Recovering A_1 (primes) from S_signed requires A_3, A_5 terms,
  each costing O(x^{2/3}).

  BARRIER 5: Weighted prime sums lose information
  -------------------------------------------------
  Any multiplicative weight w(n) satisfies w(p) = constant for all primes p.
  (Because w is multiplicative and p is an irreducible element.)
  So sum_{p<=x} w(p) = w(primes) * pi(x), giving nothing new.
  Non-multiplicative weights don't benefit from H-T's sieve structure.

  BARRIER 6: The explicit formula zero sum is non-combinatorial
  ---------------------------------------------------------------
  The oscillatory zero sum sum_rho R(x^rho) is a sum over ANALYTIC objects,
  not integers. H-T's combinatorial sieve technique is inapplicable.

  INFORMATION-THEORETIC VIEW:
  ============================
  M(x) has magnitude O(x^{1/2}) (conjectured) while the sum has x terms.
  The cancellation ratio is O(x^{-1/2}) -- massive cancellation.
  H-T exploit this: the cancellation allows "compressed" computation.

  pi(x) has magnitude O(x/log(x)), with x terms but no cancellation.
  The ratio pi(x)/x ~ 1/log(x) -- essentially NO cancellation.
  There's nothing to exploit.

  PRECISE OBSTRUCTION (the "positivity barrier"):
  ================================================
  Let f: N -> Z be a multiplicative function and F(x) = sum_{n<=x} f(n).

  THEOREM (informal): If |F(x)| = O(x^alpha) for some alpha < 1 (i.e., f has
  net cancellation), then F(x) can potentially be computed faster than the
  generic O(x^{2/3}) bound, by exploiting the cancellation.

  For mu: |M(x)| = O(x^{1/2+eps}) (RH). alpha = 1/2. H-T achieve x^{3/5}.
  For 1_prime: pi(x) ~ x/ln(x). alpha = 1. NO cancellation to exploit.

  The 1_prime indicator is not even multiplicative, but more fundamentally,
  it's a POSITIVE function with no sign changes. This positivity prevents
  any cancellation-based speedup.

  CONCLUSION:
  ===========
  ALL SIX approaches fail for the SAME fundamental reason: pi(x) counts
  a positive quantity (primes) with no signed cancellation, while H-T's
  technique is specifically designed to exploit alternating signs in mu(n).

  Converting between signed (M) and unsigned (pi) formulations always
  costs at least O(x^{2/3}), which is the original Meissel-Lehmer bound.

  STATUS: H-T transfer to pi(x) is CLOSED. No indirect route exists.
  The O(x^{3/5}) technique for M(x) does NOT improve the O(x^{2/3})
  bound for pi(x) through any known identity or conversion.
""")


# =============================================================================
# MAIN
# =============================================================================

def main():
    x_values = [100, 1000, 10000]

    print("=" * 72)
    print("H-T SIGNED CANCELLATION TRANSFER TO pi(x) -- DEEP INVESTIGATION (v2)")
    print("=" * 72)
    print(f"Testing with x = {x_values}")
    print()

    t0 = time.time()

    # Approach 1: Explicit formula / zero sum
    result1 = analyze_explicit_formula_cancellation()

    # Approach 2: Weighted prime counting
    result2 = analyze_weighted_prime_counting(x_values)

    # Approach 3: Buchstab + signed weights
    result3 = analyze_buchstab_signed(x_values)

    # Approach 4: M(x) -> pi(x) conversion
    result4 = analyze_mx_to_pix_conversion(x_values)

    # Approach 5: Concrete experiments
    result5 = concrete_experiments(x_values)

    # Approach 6: Novel identity search
    result6 = novel_identity_search(x_values)

    # Approach 7: Final obstacle analysis
    obstacle_analysis()

    elapsed = time.time() - t0

    print("=" * 72)
    print("SUMMARY OF ALL APPROACHES")
    print("=" * 72)
    print(f"""
  Approach 1 (Explicit formula / zero sum):     {result1}
    H-T is combinatorial; zero sum is analytic. Incompatible.

  Approach 2 (Weighted prime counting):          {result2}
    Multiplicative weights collapse to constant on primes.
    Non-multiplicative weights don't benefit from sieve structure.

  Approach 3 (Buchstab + signed weights):        {result3}
    Signed Buchstab IS M(x). Recovering pi(x) needs A_3, A_5 at O(x^{{2/3}}).

  Approach 4 (M(x) -> pi(x) conversion):        {result4}
    All conversion identities cost at least O(x^{{2/3}}).

  Approach 5 (Concrete experiments):             {result5}
    Identity pi(x) = sum omega(d)*M(floor(x/d)) found and verified.
    But omega partial sums cost O(x^{{2/3}}) -- no improvement.

  Approach 6 (Novel identity search):            {result6}
    The Dirichlet convolution structure forces f = omega.
    omega partial sums are as hard as pi(x). Circular.

  OVERALL VERDICT: CLOSED.
  The H-T O(x^{{3/5}}) technique for M(x) CANNOT transfer to pi(x)
  through any known or conceivable indirect route. The fundamental
  barrier is the POSITIVITY of the prime indicator: no cancellation
  exists to exploit, and converting from a signed (cancelling) sum
  to an unsigned (positive) count costs O(x^{{2/3}}) at minimum.

  Elapsed: {elapsed:.2f}s
""")


if __name__ == "__main__":
    main()

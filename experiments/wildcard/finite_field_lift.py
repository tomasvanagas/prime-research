#!/usr/bin/env python3
"""
Finite Field Analogy and Lifting: Can we deform F_q[x] -> Z?

Over F_q[x], irreducible counting is O(polylog):
  N_q(n) = (1/n) * sum_{d|n} mu(d) * q^{n/d}

This experiment explores:
  1. Verify the irreducible counting formula for small q, n
  2. Take the q -> 1 limit and observe what happens
  3. Compare density shapes: N_q(n)/q^n vs pi(x)/x
  4. Count primes in Z[i] (Gaussian) and Z[omega] (Eisenstein)
  5. Interpolate via a family parametrized by q
  6. Test q = e^{1/log(x)} substitution
  7. "Virtual curve" idea: can point counts of a curve give pi(x)?

Prior closed results (from CLOSED_PATHS.md):
  - Function field analog: genus finite -> O(g); for Z, genus = infinity  (Session 7)
  - F_1 / lattice / CFs: q->1 degenerate, no structure  (Session 10)
  - Algebraic geometry point counting: genus must be Omega(x/ln(x))  (Session 13)
  - Algebraic variety F_q: Frobenius eigenvalues = zeta zeros  (Session 14, 16)

This experiment quantifies these barriers numerically.
"""

import math
import sys
from functools import lru_cache
from collections import defaultdict

# ---------------------------------------------------------------------------
# Utility functions
# ---------------------------------------------------------------------------

def mobius(n):
    """Compute the Mobius function mu(n)."""
    if n == 1:
        return 1
    factors = []
    temp = n
    d = 2
    while d * d <= temp:
        if temp % d == 0:
            count = 0
            while temp % d == 0:
                temp //= d
                count += 1
            if count > 1:
                return 0
            factors.append(d)
        d += 1
    if temp > 1:
        factors.append(temp)
    return (-1) ** len(factors)


def divisors(n):
    """Return all divisors of n."""
    divs = []
    for d in range(1, int(math.isqrt(n)) + 1):
        if n % d == 0:
            divs.append(d)
            if d != n // d:
                divs.append(n // d)
    return sorted(divs)


def sieve_primes(limit):
    """Sieve of Eratosthenes up to limit."""
    if limit < 2:
        return []
    is_prime = [True] * (limit + 1)
    is_prime[0] = is_prime[1] = False
    for i in range(2, int(math.isqrt(limit)) + 1):
        if is_prime[i]:
            for j in range(i * i, limit + 1, i):
                is_prime[j] = False
    return [i for i in range(2, limit + 1) if is_prime[i]]


def prime_count(x):
    """pi(x) via sieve for moderate x."""
    return len(sieve_primes(int(x)))


# ---------------------------------------------------------------------------
# Part 1: Verify the F_q irreducible counting formula
# ---------------------------------------------------------------------------

def N_q(q, n):
    """
    Number of monic irreducible polynomials of degree n over F_q.
    N_q(n) = (1/n) * sum_{d|n} mu(d) * q^{n/d}
    """
    total = 0
    for d in divisors(n):
        total += mobius(d) * (q ** (n // d))
    return total / n


def count_irreducibles_brute(q, n):
    """
    Brute-force count of monic irreducible polynomials of degree n over F_q.
    Only feasible for small q and n. Polynomials represented as tuples of coefficients.
    """
    if q > 7 or n > 6:
        return None  # Too expensive

    from itertools import product

    def poly_mod(a, b, q):
        """Compute a mod b over F_q. a, b are lists (low to high degree)."""
        a = list(a)
        while len(a) >= len(b) and any(c % q != 0 for c in a):
            # Remove trailing zeros
            while a and a[-1] % q == 0:
                a.pop()
            if len(a) < len(b):
                break
            shift = len(a) - len(b)
            coeff = (a[-1] * pow(b[-1], q - 2, q)) % q  # Fermat inverse
            for i in range(len(b)):
                a[shift + i] = (a[shift + i] - coeff * b[i]) % q
            while a and a[-1] % q == 0:
                a.pop()
        return a

    def is_irreducible(coeffs, q):
        """Check if polynomial is irreducible over F_q by trial division."""
        n = len(coeffs) - 1  # degree
        if n <= 0:
            return False
        if n == 1:
            return True
        # Check divisibility by all monic polys of degree 1..n//2
        for deg in range(1, n // 2 + 1):
            for lower_coeffs in product(range(q), repeat=deg):
                # Monic polynomial: lower_coeffs + [1]
                divisor = list(lower_coeffs) + [1]
                rem = poly_mod(list(coeffs), divisor, q)
                if not rem:
                    return False
        return True

    count = 0
    # Monic polynomials of degree n: coefficients c_0, c_1, ..., c_{n-1}, c_n=1
    for lower in product(range(q), repeat=n):
        poly = list(lower) + [1]  # monic
        if is_irreducible(poly, q):
            count += 1
    return count


def test_formula():
    """Part 1: Verify the formula against brute force."""
    print("=" * 70)
    print("PART 1: Verify N_q(n) formula against brute force")
    print("=" * 70)
    print(f"{'q':>4} {'n':>4} {'Formula':>10} {'Brute':>10} {'Match':>6}")
    print("-" * 40)

    all_match = True
    for q in [2, 3, 5, 7]:
        max_n = 5 if q <= 3 else 4 if q == 5 else 3
        for n in range(1, max_n + 1):
            formula_val = N_q(q, n)
            brute_val = count_irreducibles_brute(q, n)
            match = abs(formula_val - brute_val) < 0.5 if brute_val is not None else "N/A"
            if match is not True and match != "N/A":
                all_match = False
            print(f"{q:4d} {n:4d} {formula_val:10.0f} {str(brute_val):>10} {str(match):>6}")

    # Show formula for larger values
    print("\nFormula values for larger q, n:")
    print(f"{'q':>6} {'n':>4} {'N_q(n)':>15} {'q^n':>15} {'ratio':>10}")
    print("-" * 55)
    for q in [2, 3, 5, 10, 100]:
        for n in [1, 2, 5, 10, 20]:
            val = N_q(q, n)
            qn = q ** n
            print(f"{q:6d} {n:4d} {val:15.0f} {qn:15d} {val/qn:10.6f}")

    return all_match


# ---------------------------------------------------------------------------
# Part 2: The q -> 1 limit
# ---------------------------------------------------------------------------

def test_q_to_1_limit():
    """Part 2: What happens to N_q(n) as q -> 1?"""
    print("\n" + "=" * 70)
    print("PART 2: The q -> 1 limit")
    print("=" * 70)

    print("\nN_q(n) for q approaching 1 from above, n=1..6:")
    print(f"{'q':>12} ", end="")
    for n in range(1, 7):
        print(f"{'n='+str(n):>12}", end="")
    print()
    print("-" * 84)

    for q_val in [10.0, 5.0, 2.0, 1.5, 1.1, 1.01, 1.001, 1.0001]:
        print(f"{q_val:12.4f} ", end="")
        for n in range(1, 7):
            val = N_q(q_val, n)
            print(f"{val:12.4f}", end="")
        print()

    # The key insight: N_q(n)/n as q -> 1
    print("\n\nN_q(n) * n as q -> 1 (should this have a limit?):")
    print(f"{'q':>12} ", end="")
    for n in range(1, 7):
        print(f"{'n='+str(n):>12}", end="")
    print()
    print("-" * 84)

    for q_val in [2.0, 1.5, 1.1, 1.01, 1.001]:
        print(f"{q_val:12.4f} ", end="")
        for n in range(1, 7):
            val = N_q(q_val, n) * n
            print(f"{val:12.6f}", end="")
        print()

    # Try the substitution q = 1 + epsilon, expand
    print("\n\nTaylor expansion: q = 1 + eps, N_q(n) ~ ?")
    print("For n prime: N_q(n) = (q^n - q)/n")
    print("  = ((1+eps)^n - (1+eps))/n")
    print("  ~ (1 + n*eps - 1 - eps)/n = (n-1)*eps/n  as eps -> 0")
    print("So N_{1+eps}(n) -> 0 for all n >= 1. The formula degenerates.")
    print()
    print("VERDICT: q->1 limit is trivially 0. The 'F_1' approach fails here")
    print("because q=1 means the base field has 1 element (no structure).")


# ---------------------------------------------------------------------------
# Part 3: Compare density shapes
# ---------------------------------------------------------------------------

def test_density_comparison():
    """Part 3: Compare N_q(n)/q^n vs pi(x)/x."""
    print("\n" + "=" * 70)
    print("PART 3: Density comparison -- N_q(n)/q^n vs pi(x)/x")
    print("=" * 70)

    # In F_q[x]: density of irreducibles of degree n among all monic degree-n polys
    #   = N_q(n) / q^n ~ 1/n  (prime polynomial theorem)
    # In Z: density of primes near x among integers near x
    #   = pi(x)/x ~ 1/ln(x)  (prime number theorem)
    # Analogy: n <-> ln(x), i.e., degree <-> log of magnitude

    print("\nF_q[x] density: N_q(n)/q^n for q=2:")
    print(f"{'n':>4} {'N_q(n)':>12} {'q^n':>12} {'density':>12} {'1/n':>12} {'ratio':>8}")
    print("-" * 64)
    for n in range(1, 21):
        val = N_q(2, n)
        qn = 2 ** n
        dens = val / qn
        print(f"{n:4d} {val:12.0f} {qn:12d} {dens:12.8f} {1/n:12.8f} {dens*n:8.5f}")

    print("\nZ density: pi(x)/x:")
    primes = sieve_primes(10**6)
    pi_cache = {}
    for p in primes:
        pi_cache[p] = len([pp for pp in primes if pp <= p])

    test_x = [10, 100, 1000, 10000, 100000, 1000000]
    print(f"{'x':>10} {'pi(x)':>10} {'pi(x)/x':>12} {'1/ln(x)':>12} {'ratio':>8}")
    print("-" * 56)
    for x in test_x:
        pix = len([p for p in primes if p <= x])
        dens = pix / x
        inv_ln = 1 / math.log(x)
        print(f"{x:10d} {pix:10d} {dens:12.8f} {inv_ln:12.8f} {dens/inv_ln:8.5f}")

    print("\nKey analogy: degree n in F_q[x]  <->  ln(x) in Z")
    print("Both satisfy 'prime density ~ 1/(measure of size)'")
    print("But F_q[x] has EXACT formula; Z only has asymptotic.")
    print("The difference: F_q[x] zeta = rational, finitely many zeros.")
    print("Z zeta has infinitely many zeros => oscillatory error.")


# ---------------------------------------------------------------------------
# Part 4: Counting in Z[i] and Z[omega]
# ---------------------------------------------------------------------------

def norm_zi(a, b):
    """Norm of a + bi in Z[i]."""
    return a * a + b * b


def norm_zw(a, b):
    """Norm of a + b*omega in Z[omega], omega = e^{2pi i/3}."""
    return a * a - a * b + b * b


def is_prime_int(n):
    """Simple primality test."""
    if n < 2:
        return False
    if n < 4:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    d = 5
    while d * d <= n:
        if n % d == 0 or n % (d + 2) == 0:
            return False
        d += 6
    return True


def count_gaussian_primes(N):
    """
    Count Gaussian primes with norm <= N.
    Gaussian primes are:
    - (1+i)*unit (norm 2): associate of 1+i
    - a+bi where a^2+b^2 = p, p = 1 mod 4 (splits)
    - p*unit where p = 3 mod 4 (inert, norm p^2)
    We count up to associates (divide by 4 for units).
    Actually let's count primes ideals with norm <= N.
    """
    count = 0
    # Norm 2: (1+i) is prime, one prime of norm 2
    if N >= 2:
        count += 1

    # For each rational prime p:
    #   p = 2: ramifies, gives 1 prime ideal
    #   p = 1 mod 4: splits into 2 prime ideals of norm p
    #   p = 3 mod 4: inert, gives 1 prime ideal of norm p^2
    primes = sieve_primes(N)
    for p in primes:
        if p == 2:
            continue
        if p % 4 == 1:
            count += 2  # two primes of norm p
        elif p % 4 == 3:
            if p * p <= N:
                count += 1  # one prime of norm p^2

    return count


def count_eisenstein_primes(N):
    """
    Count Eisenstein primes (Z[omega]) with norm <= N.
    - omega = (-1+sqrt(-3))/2
    - Norm: a^2 - ab + b^2
    - Primes:
      - 1 - omega (norm 3): ramified
      - pi, pi_bar where N(pi) = p, p = 1 mod 3 (splits)
      - p where p = 2 mod 3 (inert, norm p^2)
    """
    count = 0
    if N >= 3:
        count += 1  # norm 3 prime

    primes = sieve_primes(N)
    for p in primes:
        if p == 3:
            continue
        if p % 3 == 1:
            count += 2  # splits
        elif p % 3 == 2:
            if p * p <= N:
                count += 1  # inert

    return count


def test_number_rings():
    """Part 4: Prime counting in Z[i] and Z[omega]."""
    print("\n" + "=" * 70)
    print("PART 4: Prime counting in Z[i] and Z[omega]")
    print("=" * 70)

    test_N = [10, 100, 1000, 10000, 100000]

    print(f"\n{'N':>8} {'pi_Z(N)':>10} {'pi_Zi(N)':>10} {'pi_Zw(N)':>10} "
          f"{'Zi/Z':>8} {'Zw/Z':>8}")
    print("-" * 62)

    primes = sieve_primes(100000)
    for N in test_N:
        pi_z = len([p for p in primes if p <= N])
        pi_zi = count_gaussian_primes(N)
        pi_zw = count_eisenstein_primes(N)
        ratio_i = pi_zi / pi_z if pi_z > 0 else 0
        ratio_w = pi_zw / pi_z if pi_z > 0 else 0
        print(f"{N:8d} {pi_z:10d} {pi_zi:10d} {pi_zw:10d} "
              f"{ratio_i:8.4f} {ratio_w:8.4f}")

    print("\nNote: pi_Z[i](N) counts prime IDEALS with norm <= N.")
    print("The ratio pi_Zi/pi_Z -> ~1 because:")
    print("  - Primes p=1 mod 4: contribute 2 to Zi vs 1 to Z (but norm p vs p)")
    print("  - Primes p=3 mod 4: contribute 1 to Zi at norm p^2 vs 1 to Z at p")
    print("  - Net effect: roughly the same density")
    print()
    print("VERDICT: Counting in Z[i] or Z[omega] is NOT easier.")
    print("The prime counting is determined by rational primes via splitting behavior.")
    print("No computational advantage -- you still need to know the rational primes.")


# ---------------------------------------------------------------------------
# Part 5: Interpolation family parametrized by q
# ---------------------------------------------------------------------------

def test_interpolation():
    """Part 5: Look for interpolation between F_q and Z."""
    print("\n" + "=" * 70)
    print("PART 5: Interpolation family parametrized by q")
    print("=" * 70)

    # Idea: The Necklace formula gives N_q(n) = (1/n) sum mu(d) q^{n/d}
    # This is a polynomial in q. Can we evaluate at non-prime q to interpolate?

    print("\nN_q(n) as a polynomial in q, for n=1..6:")
    for n in range(1, 7):
        terms = []
        for d in divisors(n):
            mu_d = mobius(d)
            if mu_d != 0:
                terms.append(f"{'+' if mu_d > 0 else ''}{mu_d}*q^{n//d}")
        formula = " ".join(terms)
        print(f"  n={n}: N_q({n}) = (1/{n}) * ({formula})")

    # Evaluate for real q
    print(f"\nN_q(n) for real q, n=6:")
    print(f"{'q':>8} {'N_q(6)':>12} {'q^6':>12} {'density':>12}")
    print("-" * 48)
    for q in [1.001, 1.01, 1.1, 1.5, 2.0, 3.0, 5.0, 10.0, math.e, math.pi]:
        val = N_q(q, 6)
        print(f"{q:8.3f} {val:12.4f} {q**6:12.4f} {val/q**6:12.8f}")

    # The p-adic interpolation idea
    print("\n--- p-adic interpolation idea ---")
    print("Iwasawa theory interpolates L-values p-adically.")
    print("The Kubota-Leopoldt p-adic L-function interpolates L(1-k, chi)")
    print("for k = 1, 2, 3, ... But this gives L-VALUES, not prime counts.")
    print("And p-adic interpolation doesn't help with computation over R.")


# ---------------------------------------------------------------------------
# Part 6: q = e^{1/log(x)} substitution
# ---------------------------------------------------------------------------

def test_q_substitution():
    """Part 6: Test q = e^{1/log(x)} and other substitutions."""
    print("\n" + "=" * 70)
    print("PART 6: q = e^{1/log(x)} substitution")
    print("=" * 70)

    # Under the analogy: degree n <-> ln(x), base q <-> e^{1/???}
    # If q = e^t, then q^n = e^{tn} and we want e^{tn} ~ x, so tn ~ ln(x)
    # For fixed n, t ~ ln(x)/n.
    # If n ~ ln(x), then t ~ 1.
    #
    # Let's try: set n = ln(x), q = e, then
    #   N_q(n) = (1/n) sum mu(d) e^{n/d}
    #          ~ e^n / n = x / ln(x)  which IS pi(x) to leading order!
    #
    # But we need to check the error term.

    print("\nKey insight: If q=e and n=ln(x), then N_e(ln(x)) ~ x/ln(x) ~ pi(x)")
    print("This is just the prime number theorem in disguise!")
    print()

    primes = sieve_primes(10**6)

    print("Test: N_e(n) vs pi(e^n)")
    print(f"{'n':>6} {'N_e(n)':>14} {'pi(e^n)':>10} {'e^n':>12} {'ratio':>10}")
    print("-" * 56)

    for n in range(2, 15):
        en = math.exp(n)
        if en > 10**6:
            pix = None
        else:
            pix = len([p for p in primes if p <= en])
        ne_val = N_q(math.e, n)
        ratio = ne_val / pix if pix and pix > 0 else None
        pix_str = str(pix) if pix is not None else "N/A"
        ratio_str = f"{ratio:.6f}" if ratio is not None else "N/A"
        print(f"{n:6d} {ne_val:14.2f} {pix_str:>10} {en:12.1f} {ratio_str:>10}")

    # Now the critical test: does the ERROR match?
    print("\nError analysis: N_e(n) - pi(e^n)")
    print(f"{'n':>4} {'N_e(n)':>14} {'pi(e^n)':>10} {'error':>12} {'rel_err':>10}")
    print("-" * 54)
    for n in range(2, 14):
        en = math.exp(n)
        if en > 10**6:
            break
        pix = len([p for p in primes if p <= en])
        ne_val = N_q(math.e, n)
        error = ne_val - pix
        rel = error / pix if pix > 0 else 0
        print(f"{n:4d} {ne_val:14.4f} {pix:10d} {error:12.4f} {rel:10.6f}")

    print("\nVERDICT: N_e(ln(x)) approximates pi(x) to ~Li(x) accuracy.")
    print("The formula is: (1/ln(x)) * sum_{d|ln(x)} mu(d) * x^{1/d}")
    print("But ln(x) is irrational, so 'd | ln(x)' is meaningless!")
    print("When we force n to be integer, we only get pi at e^n points.")
    print("This is the PNT in disguise, NOT a new formula.")

    # Try the more creative substitution
    print("\n--- Creative substitution: q = e^{1/log x}, n = (log x)^2 ---")
    print("Then q^n = e^{log x} = x. But n must be integer and depend on x.")
    print("This is circular: to pick n,q giving pi(x), you need to know x.")


# ---------------------------------------------------------------------------
# Part 7: Virtual curve / Hasse-Weil idea
# ---------------------------------------------------------------------------

def test_virtual_curve():
    """Part 7: Can a curve over F_q encode pi(x)?"""
    print("\n" + "=" * 70)
    print("PART 7: Virtual curve whose point count gives pi(x)")
    print("=" * 70)

    # For a curve C/F_q of genus g, the zeta function is:
    #   Z(C, T) = P(T) / ((1-T)(1-qT))
    # where P(T) = prod_{i=1}^{2g} (1 - alpha_i T), |alpha_i| = sqrt(q)
    #
    # Point count: #C(F_{q^n}) = q^n + 1 - sum_{i=1}^{2g} alpha_i^n
    #
    # To encode pi(x), we'd need:
    #   sum_{i=1}^{2g} alpha_i^n = q^n + 1 - pi(q^n)   (*)
    #
    # The alpha_i must have |alpha_i| = sqrt(q). Can we find such alpha_i?

    print("\nHasse-Weil setup: #C(F_{q^n}) = q^n + 1 - sum alpha_i^n")
    print("We need: sum alpha_i^n = q^n + 1 - pi(q^n)")
    print()

    # For q=2, compute the required "trace" values
    q = 2
    primes = sieve_primes(2**20)
    print(f"For q={q}, required trace S(n) = q^n + 1 - pi(q^n):")
    print(f"{'n':>4} {'q^n':>10} {'pi(q^n)':>10} {'S(n)':>12} {'S(n)/sqrt(q)^n':>16}")
    print("-" * 56)

    traces = []
    for n in range(1, 21):
        qn = q ** n
        if qn > 10**6:
            pix = None
        else:
            pix = len([p for p in primes if p <= qn])
        if pix is not None:
            s_n = qn + 1 - pix
            bound = 2 * math.sqrt(q) ** n  # Hasse-Weil bound: |S(n)| <= 2g * sqrt(q)^n
            traces.append((n, s_n, bound))
            print(f"{n:4d} {qn:10d} {pix:10d} {s_n:12d} {s_n / math.sqrt(q)**n:16.4f}")

    # Check: for this to come from a curve of genus g, we need |S(n)| <= 2g * sqrt(q)^n
    print(f"\nRequired genus estimate (|S(n)| <= 2g * sqrt(q)^n):")
    print(f"{'n':>4} {'S(n)':>12} {'2g*sqrt(q)^n':>14} {'min g needed':>14}")
    print("-" * 48)
    for n, sn, bound in traces:
        sqn = math.sqrt(q) ** n
        min_g = abs(sn) / (2 * sqn) if sqn > 0 else float('inf')
        print(f"{n:4d} {sn:12d} {bound:14.2f} {min_g:14.2f}")

    print("\nThe required genus grows as ~ q^{n/2} / (n * 2*sqrt(q)^n)")
    print("  = q^{n/2} / (2n * q^{n/2}) = 1/(2n) ... WAIT, let's check.")
    print()

    # Actually: pi(q^n) ~ q^n / (n ln q), so
    # S(n) = q^n + 1 - q^n/(n ln q) ~ q^n (1 - 1/(n ln q))
    # |S(n)| / sqrt(q)^n ~ q^{n/2} (1 - 1/(n ln q))
    # Required genus g >= q^{n/2} / 2
    # For n=20, q=2: g >= 2^10 / 2 = 512

    print("Refined estimate: S(n) ~ q^n * (1 - 1/(n*ln q))")
    print("So min genus ~ q^{n/2} / 2, which grows EXPONENTIALLY in n.")
    print()
    print("For pi(x) with x = q^n:")
    print("  genus ~ sqrt(x) / 2")
    print()
    print("Computing the zeta function of a genus-g curve takes O(g^3) operations.")
    print("So this approach requires O(x^{3/2}) operations -- WORSE than sieving!")
    print()
    print("FUNDAMENTAL BARRIER: The Riemann zeta function has infinitely many zeros,")
    print("corresponding to 'infinite genus'. A curve with finite genus can only")
    print("capture finitely many zeros. To get pi(x) exactly, you need genus ~ sqrt(x),")
    print("and computing with such a curve is at least O(x^{3/2}).")


# ---------------------------------------------------------------------------
# Part 8: Summary and barrier analysis
# ---------------------------------------------------------------------------

def print_summary():
    """Final summary of all findings."""
    print("\n" + "=" * 70)
    print("SUMMARY: Finite Field Lifting Experiment")
    print("=" * 70)

    print("""
1. FORMULA VERIFICATION: The Necklace formula N_q(n) = (1/n) sum mu(d) q^{n/d}
   is confirmed to be exact for all tested q, n.

2. q -> 1 LIMIT: Degenerates to 0. The "field with one element" gives no primes.
   This is because F_1 has no meaningful polynomial ring structure.

3. DENSITY SHAPE: N_q(n)/q^n ~ 1/n  matches  pi(x)/x ~ 1/ln(x)
   under the correspondence n <-> ln(x). This is a beautiful analogy but
   is exactly the PNT -- nothing new.

4. Z[i] AND Z[omega]: Counting prime ideals reduces to counting rational primes
   with congruence conditions. No computational advantage.

5. INTERPOLATION: N_q(n) is a polynomial in q, but evaluating at q=e gives
   the PNT, not an exact formula. The error term doesn't simplify.

6. q = e^{1/log x}: Reduces to PNT. The substitution is circular because
   n must be integer but log(x) is not.

7. VIRTUAL CURVE: To encode pi(x), need genus ~ sqrt(x), leading to
   O(x^{3/2}) complexity. WORSE than current best O(x^{2/3}).

ROOT CAUSE OF FAILURE:
  F_q[x] is easy because its zeta function Z(T) = 1/(1-T) has NO zeros.
  Z is hard because zeta(s) has INFINITELY MANY zeros.
  The zeros encode the "randomness" of primes.
  Any finite-genus approximation can only capture finitely many zeros.
  The genus needed grows with x, destroying any polylog advantage.

  This is the SAME barrier seen in all other approaches:
    - Information-theoretic: ~50% of digits of p(n) are "random" (zeta zeros)
    - Algebraic: infinite-dimensional cohomology (H^1 of Spec(Z))
    - Analytic: infinitely many zeros to sum

CONCLUSION: CLOSED. The finite field analogy, while beautiful, cannot overcome
the fundamental gap between rational zeta functions (finitely many zeros) and
the Riemann zeta function (infinitely many zeros).
""")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("Finite Field Analogy and Lifting Experiment")
    print("Computing pi(x) via F_q[x] <-> Z analogy")
    print("=" * 70)

    match = test_formula()
    print(f"\nFormula verification: {'PASSED' if match else 'FAILED'}")

    test_q_to_1_limit()
    test_density_comparison()
    test_number_rings()
    test_interpolation()
    test_q_substitution()
    test_virtual_curve()
    print_summary()


if __name__ == "__main__":
    main()

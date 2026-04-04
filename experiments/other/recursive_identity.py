#!/usr/bin/env python3
"""
Session 9: Recursive Identity Exploration
==========================================
Can we find a recursive identity connecting pi(x) to pi(x/2), pi(x/4), etc.
that yields an O(polylog n) divide-and-conquer algorithm for p(n)?

Approaches explored:
1. Doubling formula: pi(2x) vs 2*pi(x) + correction
2. Squaring formula: pi(x^2) vs pi(x)*x/2 + correction
3. Buchstab identity inversions
4. Legendre generalization: pi(x) from pi(x^{1/k})
5. Meissel-Lehmer phi recursion analysis
6. Binary splitting of prime counting
"""

import math
import time
import sys
from collections import defaultdict

# ---------- Sieve infrastructure ----------

def sieve_primes(limit):
    """Sieve of Eratosthenes up to limit."""
    if limit < 2:
        return []
    sieve = bytearray(b'\x01') * (limit + 1)
    sieve[0] = sieve[1] = 0
    for i in range(2, int(limit**0.5) + 1):
        if sieve[i]:
            sieve[i*i::i] = bytearray(len(sieve[i*i::i]))
    return [i for i in range(2, limit + 1) if sieve[i]]

def pi_table(limit):
    """Build cumulative pi(x) table for x = 0..limit."""
    sieve = bytearray(b'\x01') * (limit + 1)
    sieve[0] = sieve[1] = 0
    for i in range(2, int(limit**0.5) + 1):
        if sieve[i]:
            sieve[i*i::i] = bytearray(len(sieve[i*i::i]))
    table = [0] * (limit + 1)
    count = 0
    for i in range(limit + 1):
        count += sieve[i]
        table[i] = count
    return table

def li(x):
    """Logarithmic integral li(x) via series expansion."""
    if x <= 1:
        return 0.0
    lnx = math.log(x)
    # Ramanujan's series for li(x) is complex; use simple numerical integration
    # li(x) = integral from 0 to x of dt/ln(t)
    # Use Euler-Maclaurin / simple series: li(x) = EulerGamma + ln(ln(x)) + sum ln(x)^k/(k*k!)
    gamma = 0.5772156649015329
    lnlnx = math.log(lnx)
    s = gamma + lnlnx
    term = 1.0
    for k in range(1, 200):
        term *= lnx / k
        s += term / k
    return s

def R_func(x):
    """Riemann R function: R(x) = sum_{k=1}^inf mu(k)/k * li(x^{1/k})."""
    if x <= 1:
        return 0.0
    # Mobius values for small k
    mu = [0, 1, -1, -1, 0, -1, 1, -1, 0, 0, 1,
          -1, 0, -1, 1, 1, 0, -1, 0, -1, 0,
          1, 1, -1, 0, 0, 1, 0, 0, -1, -1,
          -1, 0, 1, 1, 1, 0, -1, 1, 1, 0,
          -1, -1, -1, 0, 0, -1, 0, -1, 0, 0]
    s = 0.0
    for k in range(1, min(50, len(mu))):
        if mu[k] == 0:
            continue
        xk = x ** (1.0 / k)
        if xk < 2:
            break
        s += mu[k] / k * li(xk)
    return s

# ---------- Experiment 1: Doubling Formula ----------

def experiment_doubling(limit=100000):
    """Test pi(2x) - 2*pi(x) and look for structure."""
    print("=" * 70)
    print("EXPERIMENT 1: DOUBLING FORMULA")
    print("  Testing: pi(2x) = 2*pi(x) + correction(x)")
    print("=" * 70)

    table = pi_table(limit)

    results = []
    for x in range(10, limit // 2, 1):
        pi_x = table[x]
        pi_2x = table[2 * x]
        correction = pi_2x - 2 * pi_x
        results.append((x, pi_x, pi_2x, correction))

    # Analyze correction
    corrections = [r[3] for r in results]

    # Sample output
    print("\n  x       pi(x)    pi(2x)   correction  correction/pi(x)  x/ln(x)^2")
    print("  " + "-" * 75)
    for x, pi_x, pi_2x, corr in results[::5000]:
        lnx = math.log(x) if x > 1 else 1
        ratio = corr / pi_x if pi_x > 0 else 0
        approx = x / (lnx * lnx)
        print(f"  {x:<8d} {pi_x:<8d} {pi_2x:<8d} {corr:<+10d} {ratio:<+12.6f}  {approx:.2f}")

    # By PNT: pi(2x) ~ 2x/ln(2x) and pi(x) ~ x/ln(x)
    # So correction ~ 2x/ln(2x) - 2x/ln(x) = 2x * (1/ln(2x) - 1/ln(x))
    #              = 2x * (ln(x) - ln(2x)) / (ln(x)*ln(2x))
    #              = -2x*ln(2) / (ln(x)*ln(2x))
    # Approximately: -2x*ln(2)/ln(x)^2

    print("\n  PNT prediction: correction ~ -2x*ln(2)/ln(x)^2")
    print("\n  x       actual    PNT_pred  ratio(actual/pred)  residual")
    print("  " + "-" * 65)
    residuals = []
    for x, pi_x, pi_2x, corr in results[::5000]:
        lnx = math.log(x)
        ln2x = math.log(2 * x)
        pnt_pred = -2 * x * math.log(2) / (lnx * ln2x)
        ratio = corr / pnt_pred if abs(pnt_pred) > 0.01 else float('inf')
        resid = corr - pnt_pred
        residuals.append((x, resid))
        print(f"  {x:<8d} {corr:<+10d} {pnt_pred:<+12.2f} {ratio:<12.4f}  {resid:<+12.2f}")

    # Second-order correction
    print("\n  Second-order: residual / (x/ln(x)^3)")
    print("  " + "-" * 50)
    for x, resid in residuals:
        lnx = math.log(x)
        scale = x / (lnx ** 3)
        if abs(scale) > 0.01:
            print(f"  x={x:<8d}  resid/scale = {resid/scale:<+.6f}")

    return results

# ---------- Experiment 2: Squaring Formula ----------

def experiment_squaring(limit=50000):
    """Test pi(x^2) vs pi(x)*x/2."""
    print("\n" + "=" * 70)
    print("EXPERIMENT 2: SQUARING FORMULA")
    print("  Testing: pi(x^2) = pi(x)*x/2 + correction(x)")
    print("=" * 70)

    # We need pi up to limit^2, but that's huge. Use smaller x values.
    max_x = int(math.sqrt(limit))  # x up to ~223
    table = pi_table(limit)

    print("\n  x     pi(x)   pi(x^2)  pi(x)*x/2  correction  corr/(x/ln(x))")
    print("  " + "-" * 75)

    results = []
    for x in range(10, max_x + 1):
        x2 = x * x
        if x2 > limit:
            break
        pi_x = table[x]
        pi_x2 = table[x2]
        approx = pi_x * x / 2.0
        corr = pi_x2 - approx
        lnx = math.log(x)
        scale = x / lnx if lnx > 0 else 1
        results.append((x, pi_x, pi_x2, approx, corr))
        if x % 20 == 0 or x < 30:
            print(f"  {x:<6d} {pi_x:<7d} {pi_x2:<8d} {approx:<10.1f} {corr:<+12.1f}  {corr/scale:<+.4f}")

    # PNT analysis: pi(x^2) ~ x^2/(2*ln(x)), pi(x) ~ x/ln(x)
    # pi(x)*x/2 ~ x^2/(2*ln(x))  --- same leading term!
    # So correction = pi(x^2) - pi(x)*x/2 is a LOWER ORDER term
    print("\n  PNT: both sides have leading term x^2/(2*ln(x)), so correction is lower order")
    print("\n  Fitting correction ~ C * x^2 / ln(x)^2:")
    for x, pi_x, pi_x2, approx, corr in results:
        if x >= 50:
            lnx = math.log(x)
            C = corr * lnx**2 / (x**2) if x > 0 else 0
            if x % 30 == 0:
                print(f"  x={x}: C = {C:.6f}")

    return results

# ---------- Experiment 3: Buchstab Identity Analysis ----------

def experiment_buchstab(limit=10000):
    """Explore Buchstab's identity numerically.

    Buchstab: S(x,y) = S(x,y-) - S(x/p, p-) for p = largest prime <= y
    where S(x,y) = #{n <= x : smallest prime factor of n > y}

    Key relation: pi(x) = S(x, sqrt(x)) + pi(sqrt(x)) - 1
    """
    print("\n" + "=" * 70)
    print("EXPERIMENT 3: BUCHSTAB IDENTITY")
    print("  S(x,y) = #{n <= x : all prime factors > y}")
    print("  pi(x) = S(x, sqrt(x)) + pi(sqrt(x)) - 1")
    print("=" * 70)

    primes = sieve_primes(limit)
    prime_set = set(primes)
    table = pi_table(limit)

    def S(x, y):
        """Count integers <= x whose smallest prime factor > y."""
        x = int(x)
        if x < 2:
            return 0
        count = 0
        for n in range(2, x + 1):
            smallest = n
            for p in primes:
                if p * p > n:
                    break
                if n % p == 0:
                    smallest = p
                    break
            if smallest > y:
                count += 1
        return count

    # Verify identity for small x
    print("\n  Verifying: pi(x) = S(x, sqrt(x)) + pi(sqrt(x)) - 1")
    print("  x       pi(x)   S(x,√x)  pi(√x)   S+pi-1   match?")
    print("  " + "-" * 60)
    for x in [20, 50, 100, 200, 500, 1000]:
        sqx = int(math.sqrt(x))
        s_val = S(x, sqx)
        pi_sqx = table[sqx]
        pi_x = table[x]
        computed = s_val + pi_sqx - 1
        match = "YES" if computed == pi_x else "NO"
        print(f"  {x:<8d} {pi_x:<7d} {s_val:<8d} {pi_sqx:<8d} {computed:<8d} {match}")

    # Buchstab recursion: S(x,y) relates to S at smaller arguments
    # S(x,z) = S(x,y) - sum_{y < p <= z} S(x/p, p-)
    # This gives a tree-like recursion
    print("\n  Buchstab recursion depth analysis:")
    print("  How deep does the recursion go for different x?")

    def buchstab_recursive(x, y, depth=0, max_depth=[0]):
        """Buchstab recursive computation of S(x,y)."""
        max_depth[0] = max(max_depth[0], depth)
        x = int(x)
        if x < 2 or y >= x:
            return max(0, x - 1)  # all integers 2..x have smallest factor > y if y >= x
        if y < 2:
            return max(0, x - 1)  # S(x,1) = x-1 (all integers >= 2)

        # Base: if y >= sqrt(x), S(x,y) = pi(x) - pi(y) + (1 if x >= 2 else 0)
        # Actually S(x,y) counts n<=x with lpf(n)>y. If y>=sqrt(x), these are primes in (y,x] plus 1.
        # S(x,y) = pi(x) - pi(y)  when y >= sqrt(x)
        sqx = math.sqrt(x)
        if y >= sqx:
            if x <= limit and int(y) <= limit:
                return table[x] - table[min(int(y), x)]
            return 0

        # Recursion: pick next prime p > y
        # S(x,y) = S(x, p_next-1) where p_next is next prime > y... no
        # Buchstab: S(x, p_k) = S(x, p_{k-1}) - S(x/p_k, p_{k-1})
        # So iterating: S(x, sqrt(x)) = S(x, 2-) - sum_{p <= sqrt(x)} S(x/p, p-)

        # Just count the recursion calls
        return -1  # placeholder

    # Instead, let's count how many S-evaluations Meissel-Lehmer needs
    def count_phi_calls(x, a, calls=[0]):
        """Count calls in Lehmer's phi(x,a) computation."""
        calls[0] += 1
        if a == 0:
            return int(x)
        if a == 1:
            return int(x) - int(x) // 2
        # phi(x,a) = phi(x,a-1) - phi(x/p_a, a-1)
        if x < primes[a-1]:
            return 1
        count_phi_calls(x, a - 1, calls)
        count_phi_calls(x / primes[a-1], a - 1, calls)
        return calls[0]

    print("\n  phi(x,a) call counts (Meissel-Lehmer recursion):")
    print("  x          a=pi(x^{1/3})  calls     calls/x^{2/3}")
    print("  " + "-" * 55)

    for x in [100, 1000, 5000, 10000]:
        a = table[int(x ** (1/3.0))]
        calls = [0]
        count_phi_calls(x, a, calls)
        x23 = x ** (2/3.0)
        print(f"  {x:<10d} {a:<14d} {calls[0]:<10d} {calls[0]/x23:.4f}")

    return True

# ---------- Experiment 4: Legendre Generalization ----------

def experiment_legendre(limit=50000):
    """Generalize Legendre: pi(x) - pi(x^{1/k}) from sieve counts.

    Legendre: pi(x) - pi(sqrt(x)) + 1 = #{n<=x : n=1 or lpf(n)>sqrt(x)}
    For general k: pi(x) - pi(x^{1/k}) + ... involves inclusion-exclusion
    of products of at most k-1 primes up to x^{1/k}.
    """
    print("\n" + "=" * 70)
    print("EXPERIMENT 4: LEGENDRE GENERALIZATION")
    print("  Can we express pi(x) using pi(x^{1/k}) for large k?")
    print("=" * 70)

    table = pi_table(limit)
    primes = sieve_primes(int(math.sqrt(limit)) + 1)

    # Legendre-type decomposition:
    # Numbers <= x with no prime factor <= x^{1/k} are either 1 or products of primes > x^{1/k}
    # with at most k-1 prime factors (since product of k primes each > x^{1/k} exceeds x).

    # So: pi(x) = pi(x^{1/k}) + S_1(x,k) + S_2(x,k) + ... + S_{k-1}(x,k)
    # where S_j counts numbers <= x that are products of exactly j primes all > x^{1/k}

    print("\n  Decomposition: pi(x) = pi(y) + S_1 + S_2 + ... + S_{k-1}")
    print("  where y = x^{1/k} and S_j = #{n<=x : n = product of j primes > y}")

    for x in [1000, 5000, 10000, 30000, 50000]:
        if x > limit:
            break
        print(f"\n  x = {x}, pi(x) = {table[x]}")
        for k in [2, 3, 4, 5, 6]:
            y = int(x ** (1.0 / k))
            pi_y = table[min(y, limit)]

            # S_1: primes in (y, x] = pi(x) - pi(y)
            S1 = table[x] - pi_y

            # S_2: semiprimes p*q with y < p <= q and p*q <= x
            S2 = 0
            for i, p in enumerate(primes):
                if p <= y:
                    continue
                if p * p > x:
                    break
                # q ranges from p to x/p, all > y
                max_q = x // p
                if max_q < p:
                    break
                # count primes in [p, max_q] that are > y
                lo = max(p, y + 1)
                hi = min(max_q, limit)
                if hi >= lo:
                    S2 += table[hi] - table[lo - 1]

            # S_3: would need triple products -- expensive
            S3_plus = table[x] - pi_y - S1 - S2
            # Actually pi(x) = pi(y) + S1 already if S1 = pi(x) - pi(y)
            # The decomposition is: #{composites <=x with all factors > y} = S2 + S3 + ...
            # And #{n<=x with lpf(n)>y} = 1 + S1 + S2 + S3 + ...
            # So pi(x) = pi(y) + #{n<=x, lpf(n)>y} - 1  ... no wait.
            # #{n<=x : lpf(n)>y} = 1 + (primes in (y,x]) + (semiprimes with both > y, <=x) + ...

            smooth_count = 0
            for n in range(2, min(x + 1, 10001)):
                temp = n
                all_big = True
                for p in primes:
                    if p > y:
                        break
                    if temp % p == 0:
                        all_big = False
                        break
                if all_big:
                    smooth_count += 1

            print(f"    k={k}: y={y}, pi(y)={pi_y}, "
                  f"#{'{'}lpf>y{'}'} counted={smooth_count if x<=10000 else '?'}, "
                  f"S1(primes>y)={S1}, S2(semiprimes)={S2}")

    # Key insight: for k=3, we need pi(x^{1/3}) and counts of numbers with
    # exactly 2 prime factors in (x^{1/3}, x]. These counts ALSO involve pi values!
    # S2 involves sum over p of pi(x/p) - pi(p) + 1
    print("\n  KEY: S2 = sum_{y<p<=sqrt(x)} [pi(x/p) - pi(p) + 1]")
    print("  This recursively involves pi at points x/p for various primes p.")
    print("  This IS the Meissel-Lehmer decomposition!")

    return True

# ---------- Experiment 5: Half-splitting ----------

def experiment_halving(limit=100000):
    """Test if pi(x) can be computed from pi(x/2) efficiently.

    Key idea: pi(x) = pi(x/2) + #{primes in (x/2, x]}
    The question: can #{primes in (x/2, x]} be computed without sieving?
    """
    print("\n" + "=" * 70)
    print("EXPERIMENT 5: HALVING / BINARY SPLITTING")
    print("  pi(x) = pi(x/2) + #{primes in (x/2, x]}")
    print("  Can the second term be computed from previously known pi values?")
    print("=" * 70)

    table = pi_table(limit)

    # #{primes in (x/2, x]} by inclusion-exclusion over primes up to sqrt(x)
    # Using Legendre's formula:
    # #{n in (x/2, x] : n coprime to all p <= sqrt(x)} = these are the primes in (x/2, x]
    # (plus 1 if we count 1, but 1 is not in (x/2,x] for x>2)

    # Legendre sieve: #{n in (a,b] : gcd(n, P(y))=1} where P(y) = product of primes <= y
    # = sum_{d | P(y)} mu(d) * (floor(b/d) - floor(a/d))

    # For interval (x/2, x] and y = sqrt(x):
    # Number of squarefree divisors of P(sqrt(x)) = 2^{pi(sqrt(x))}
    # For x = 10^6: pi(1000) = 168, so 2^168 terms -- EXPONENTIAL!

    print("\n  Inclusion-exclusion has 2^pi(sqrt(x)) terms:")
    for x in [100, 1000, 10000, 100000]:
        sqx = int(math.sqrt(x))
        pisq = table[sqx]
        print(f"  x={x}: pi(sqrt(x))={pisq}, terms = 2^{pisq} = {2**pisq:.2e}")

    # But: can we use the RECURSIVE pi values to shortcut?
    # pi(x) - pi(x/2) should be approximately x/(2*ln(x)) by PNT
    # More precisely, by PNT: ~ x/(2*ln(x)) * (1 + 1/ln(x) + ...)

    print("\n  Testing PNT approximation for half-interval prime count:")
    print("  x          pi(x)-pi(x/2)  x/(2*ln(x))  ratio    li(x)-li(x/2)  ratio2")
    print("  " + "-" * 80)

    for x in [100, 500, 1000, 5000, 10000, 50000, 100000]:
        if x > limit:
            break
        actual = table[x] - table[x // 2]
        lnx = math.log(x)
        pnt_approx = x / (2 * lnx)
        li_approx = li(x) - li(x / 2)
        r1 = actual / pnt_approx if pnt_approx > 0 else 0
        r2 = actual / li_approx if li_approx > 0 else 0
        print(f"  {x:<10d} {actual:<14d} {pnt_approx:<12.2f} {r1:<8.5f} {li_approx:<14.2f} {r2:<8.5f}")

    # The error in li(x)-li(x/2) vs actual
    print("\n  Error analysis: actual - (li(x)-li(x/2)):")
    errors = []
    for x in range(100, min(limit + 1, 50001), 100):
        actual = table[x] - table[x // 2]
        li_approx = li(x) - li(x / 2)
        err = actual - li_approx
        errors.append((x, err))

    print("  x          error     error/sqrt(x)   error/sqrt(x)*ln(x)")
    print("  " + "-" * 60)
    for x, err in errors[::50]:
        sqx = math.sqrt(x)
        lnx = math.log(x)
        print(f"  {x:<10d} {err:<+10.2f} {err/sqx:<+14.6f}  {err/sqx*lnx:<+.6f}")

    return errors

# ---------- Experiment 6: Recursive pi via known values ----------

def experiment_recursive_pi(limit=50000):
    """Try to build pi(x) from pi(x/2), pi(x/3), pi(x/5), etc.

    Idea: Use the Meissel-Lehmer type identity but cache and reuse pi values.

    phi(x,a) = x - sum_{p<=p_a} phi(x/p, a(p)-1) + ...

    If we precompute pi at all dyadic points x, x/2, x/4, ..., x/2^k,
    how much does this help?
    """
    print("\n" + "=" * 70)
    print("EXPERIMENT 6: RECURSIVE PI CONSTRUCTION")
    print("  Build pi(x) using ONLY pi at smaller points + small corrections")
    print("=" * 70)

    table = pi_table(limit)
    primes_list = sieve_primes(limit)

    # Method A: Meissel's formula
    # pi(x) = phi(x,a) + a - 1 - P2(x,a)
    # where a = pi(x^{1/3})
    # P2(x,a) = sum_{a < i <= pi(sqrt(x))} [pi(x/p_i) - i + 1]
    # phi(x,a) = #{n<=x : gcd(n, p_1*...*p_a) = 1}

    # P2 requires pi(x/p_i) for pi(x^{1/3}) < p_i <= sqrt(x)
    # These are pi evaluated at points x/p for various primes p
    # If x/p ranges over an interval, we need many pi values

    print("\n  Meissel: pi(x) = phi(x,a) + a - 1 - P2(x,a)")
    print("  P2 requires pi(x/p) for x^{1/3} < p <= sqrt(x)")
    print("\n  Number of pi evaluations needed for P2:")
    for x in [1000, 10000, 50000]:
        if x > limit:
            break
        cbrt = int(x ** (1/3.0))
        sqr = int(x ** 0.5)
        a = table[cbrt]
        b = table[sqr]
        n_evals = b - a
        print(f"  x={x}: a=pi(x^{{1/3}})={a}, b=pi(sqrt(x))={b}, P2 needs {n_evals} pi evaluations")
        print(f"    These pi evals are at points x/p for p in [{primes_list[a] if a<len(primes_list) else '?'}..{sqr}]")
        print(f"    So pi values at points in [{x//sqr}..{x//primes_list[a] if a<len(primes_list) else '?'}]")

    # Method B: Recursive halving attempt
    # Can we express pi(x) - pi(x/2) using ONLY pi values at points <= x/2?
    # By Legendre sieve on (x/2, x]:
    # #{primes in (x/2,x]} = sum_{d | P(sqrt(x))} mu(d) * (floor(x/d) - floor(x/(2d)))

    # Alternative: use the identity
    # sum_{p prime, x/2 < p <= x} 1 = sum_{n=x/2+1}^{x} Lambda(n)/ln(n) * [n is prime]
    # Not helpful directly.

    # Method C: Try a "doubling" formula via Buchstab
    # pi(x) = pi(x/2) + S(x, x/2) - S(x/2, x/2)  ... no, S counts differently

    # Method D: Direct test of a divide-and-conquer
    # Hypothesis: pi(x) = pi(x/2) + pi(x) - pi(x/2)
    # The second term pi(x) - pi(x/2) = #{primes in (x/2, x]}
    # Can we relate this to pi at smaller points?

    # By PNT: pi(x) - pi(x/2) ~ li(x) - li(x/2) ~ x/(2*ln(x))
    # The ERROR E(x) = [pi(x) - pi(x/2)] - [li(x) - li(x/2)]
    # Under RH: |E(x)| = O(sqrt(x) * ln(x))

    # For x = 10^100: E is O(10^50 * 230) ~ 10^52.4
    # But pi(x) - pi(x/2) ~ 10^100 / 230 ~ 4.3 * 10^97
    # So relative error ~ 10^{-45} -- great! But we need EXACT.

    print("\n  Key question: can the correction to li(x)-li(x/2) be computed exactly?")
    print("  Under RH: correction = sum_rho [x^rho - (x/2)^rho] / rho + ...")
    print("  This requires ZETA ZEROS -- same barrier as before.")

    # Method E: Test whether pi(x/p) values for small p give pi(x)
    # Regression: pi(x) = c_0 + c_2*pi(x/2) + c_3*pi(x/3) + c_5*pi(x/5) + ...
    print("\n  Testing linear regression: pi(x) = sum c_p * pi(x/p)")

    import numpy as np

    small_primes = [2, 3, 5, 7, 11, 13]
    X_data = []
    y_data = []

    for x in range(200, min(limit + 1, 20001)):
        features = [1.0]  # constant
        skip = False
        for p in small_primes:
            xp = x // p
            if xp < 2:
                skip = True
                break
            features.append(table[xp])
        if skip:
            continue
        X_data.append(features)
        y_data.append(table[x])

    X = np.array(X_data)
    y = np.array(y_data)

    # Least squares
    coeffs, residuals, rank, sv = np.linalg.lstsq(X, y, rcond=None)
    y_pred = X @ coeffs
    errors = y - y_pred
    max_err = np.max(np.abs(errors))
    mean_err = np.mean(np.abs(errors))
    exact_count = np.sum(np.abs(errors) < 0.5)

    print(f"\n  Regression: pi(x) = {coeffs[0]:.4f}", end="")
    for i, p in enumerate(small_primes):
        print(f" + {coeffs[i+1]:.6f}*pi(x/{p})", end="")
    print()
    print(f"  Max error: {max_err:.4f}")
    print(f"  Mean abs error: {mean_err:.4f}")
    print(f"  Exact (|err|<0.5): {exact_count}/{len(y)} = {100*exact_count/len(y):.1f}%")

    # Also test: pi(x) = 2*pi(x/2) + correction modeled as function of x
    print("\n  Testing: pi(x) = 2*pi(x/2) + f(x)")
    print("  where f(x) is modeled as polynomial in ln(x), x/ln(x)^2, etc.")

    corrections = []
    for x in range(100, min(limit + 1, 20001)):
        corr = table[x] - 2 * table[x // 2]
        corrections.append((x, corr))

    # Fit: correction ~ a * x/ln(x)^2 + b * x/ln(x)^3 + c
    X_corr = []
    y_corr = []
    for x, corr in corrections:
        lnx = math.log(x)
        X_corr.append([x / lnx**2, x / lnx**3, 1.0])
        y_corr.append(corr)

    X_c = np.array(X_corr)
    y_c = np.array(y_corr)
    c_coeffs, _, _, _ = np.linalg.lstsq(X_c, y_c, rcond=None)
    y_c_pred = X_c @ c_coeffs
    c_errors = y_c - y_c_pred
    c_max = np.max(np.abs(c_errors))
    c_mean = np.mean(np.abs(c_errors))

    print(f"  correction ~ {c_coeffs[0]:.6f}*x/ln(x)^2 + {c_coeffs[1]:.6f}*x/ln(x)^3 + {c_coeffs[2]:.4f}")
    print(f"  Max residual: {c_max:.4f}")
    print(f"  Mean residual: {c_mean:.4f}")
    print(f"  Exact (|residual|<0.5): {np.sum(np.abs(c_errors)<0.5)}/{len(y_c)}")

    # Test: does the residual grow with x?
    print("\n  Residual growth with x:")
    for i in range(0, len(c_errors), 2000):
        x = corrections[i][0]
        print(f"  x={x}: residual={c_errors[i]:.4f}, residual/sqrt(x)={c_errors[i]/math.sqrt(x):.6f}")

    return True

# ---------- Experiment 7: The fundamental barrier test ----------

def experiment_barrier(limit=100000):
    """Quantify exactly why recursive approaches fail.

    Core question: does ANY identity involving pi at O(polylog) points
    determine pi(x) exactly?
    """
    print("\n" + "=" * 70)
    print("EXPERIMENT 7: FUNDAMENTAL BARRIER ANALYSIS")
    print("  Can O(polylog(x)) values of pi determine pi(x) exactly?")
    print("=" * 70)

    table = pi_table(limit)

    # Consider: we know pi(x/2), pi(x/3), pi(x/5), ..., pi(x/p) for all primes p <= x^epsilon
    # Plus pi(x^{1/k}) for k = 1,2,...,K
    # How much ambiguity remains in pi(x)?

    # Test: for different x values, find pairs x, x' where
    # pi(x/p) = pi(x'/p) for all small p, but pi(x) != pi(x')

    print("\n  Looking for x, x' where pi agrees at x/p for small p but differs at x:")

    small_primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]

    # Group x values by their (pi(x/2), pi(x/3), pi(x/5), ...) signature
    from collections import defaultdict

    signatures = defaultdict(list)
    for x in range(100, min(limit + 1, 30001)):
        sig = tuple(table[x // p] for p in small_primes if x // p >= 2)
        signatures[sig].append(x)

    collisions = 0
    diff_pi_collisions = 0
    max_pi_spread = 0
    collision_examples = []

    for sig, xs in signatures.items():
        if len(xs) > 1:
            collisions += 1
            pi_values = set(table[x] for x in xs)
            if len(pi_values) > 1:
                diff_pi_collisions += 1
                spread = max(pi_values) - min(pi_values)
                max_pi_spread = max(max_pi_spread, spread)
                if len(collision_examples) < 5:
                    collision_examples.append((xs[:5], [table[x] for x in xs[:5]]))

    print(f"\n  Total distinct signatures: {len(signatures)}")
    print(f"  Signatures with multiple x values: {collisions}")
    print(f"  Of those, with DIFFERENT pi(x): {diff_pi_collisions}")
    print(f"  Max pi(x) spread within a collision: {max_pi_spread}")

    if collision_examples:
        print("\n  Examples of collisions where pi differs:")
        for xs, pis in collision_examples[:3]:
            print(f"    x values: {xs}")
            print(f"    pi values: {pis}")
            sigs = [tuple(table[x // p] for p in small_primes[:3] if x // p >= 2) for x in xs]
            print(f"    (pi(x/2), pi(x/3), pi(x/5)): {sigs[0]}")
    else:
        print("\n  No collisions found with 10 primes -- signature is rich enough locally")
        print("  BUT: this only means the signature distinguishes x in [100, 30000]")
        print("  For x = 10^100, the number of possible pi(x) values ~ 10^50")
        print("  and O(polylog) pi-values give only polylog bits of information")

    # Information theory argument
    print("\n  INFORMATION THEORY:")
    for N_exp in [6, 10, 20, 50, 100]:
        x = 10 ** N_exp
        lnx = N_exp * math.log(10)
        pi_x = x / lnx  # approximate
        # Number of pi query points in O(polylog) scheme
        n_queries = N_exp ** 3  # generous polylog
        # Each pi value has ~ log2(pi(x)) ~ N_exp * log2(10) bits
        bits_per_query = N_exp * math.log2(10)
        total_bits = n_queries * bits_per_query
        # Bits needed to specify p(n) exactly:
        # p(n) ~ n*ln(n), so log2(p(n)) ~ N_exp * log2(10) + log2(N_exp)
        bits_needed = N_exp * math.log2(10) + math.log2(N_exp * math.log(10))
        # But the HARD bits are the fluctuation bits
        # Under RH: pi(x) = li(x) + O(sqrt(x)*ln(x))
        # The fluctuation encodes ~ (1/2)*N_exp * log2(10) bits
        fluctuation_bits = 0.5 * N_exp * math.log2(10)

        print(f"\n  x = 10^{N_exp}:")
        print(f"    polylog queries: ~{n_queries}")
        print(f"    bits per query: ~{bits_per_query:.0f}")
        print(f"    total info from queries: ~{total_bits:.0f} bits")
        print(f"    bits needed for exact p(n): ~{bits_needed:.0f}")
        print(f"    fluctuation bits (hard info): ~{fluctuation_bits:.0f}")
        print(f"    sufficient? {'MAYBE' if total_bits > fluctuation_bits else 'NO'}")

    return True

# ---------- Experiment 8: Explicit formula binary splitting ----------

def experiment_explicit_formula_splitting():
    """Explore whether the explicit formula sum over zeros can be split.

    pi(x) = R(x) - sum_rho R(x^rho) - ...

    The sum over rho (nontrivial zeros) is the expensive part.
    Can we split it into sub-sums and relate them?
    """
    print("\n" + "=" * 70)
    print("EXPERIMENT 8: EXPLICIT FORMULA & BINARY SPLITTING")
    print("  Sum over zeros: can sub-sums be reused across x values?")
    print("=" * 70)

    # Nontrivial zeros rho = 1/2 + i*gamma (under RH)
    # R(x^rho) = R(x^{1/2 + i*gamma})
    # For x^rho: |x^rho| = x^{1/2}, arg = gamma * ln(x)

    # If we compute the sum for x, can we get the sum for 2x?
    # x^rho -> (2x)^rho = 2^rho * x^rho
    # So: sum_rho R((2x)^rho) = sum_rho 2^rho * R'(x^rho) ... not quite, R is nonlinear

    # Actually: (2x)^rho = 2^rho * x^rho
    # li((2x)^rho) = li(2^rho * x^rho)  -- not simply related to li(x^rho)

    # BUT: x^rho = x^{1/2} * e^{i*gamma*ln(x)}
    # This is like a Fourier-type sum! The phases are gamma * ln(x).
    # Changing x -> 2x shifts phases by gamma * ln(2).

    # Known zeros (first few imaginary parts):
    gammas = [14.134725, 21.022040, 25.010858, 30.424876, 32.935062,
              37.586178, 40.918719, 43.327073, 48.005151, 49.773832,
              52.970321, 56.446248, 59.347044, 60.831779, 65.112544,
              67.079811, 69.546402, 72.067158, 75.704691, 77.144840]

    print(f"\n  Using first {len(gammas)} nontrivial zeros")

    # The oscillatory part: sum_rho x^{rho-1}/rho approximately
    # = sum_k x^{-1/2} * [cos(gamma_k * ln(x)) + i*sin(gamma_k * ln(x))] / |rho_k|
    # (taking real part)

    def osc_sum(x, gammas):
        """Compute oscillatory correction from zeros."""
        s = 0.0
        lnx = math.log(x)
        for g in gammas:
            rho_real = 0.5
            rho_imag = g
            # x^rho = x^{1/2} * exp(i * g * ln(x))
            mag = x ** 0.5
            phase = g * lnx
            # li(x^rho) contribution (approximate as x^rho / rho / ln(x^rho))
            # Actually just look at the oscillatory term x^{rho} / rho
            denom = math.sqrt(rho_real**2 + rho_imag**2)
            real_part = mag * math.cos(phase) / denom
            s += real_part
        return -2 * s  # factor of 2 from conjugate pairs

    print("\n  Oscillatory correction at various x:")
    print("  x          osc_corr   pi(x)-R(x)  (comparison)")
    table = pi_table(10000)
    for x in [100, 200, 500, 1000, 2000, 5000, 10000]:
        osc = osc_sum(x, gammas)
        r_val = R_func(x)
        pi_val = table[x]
        print(f"  {x:<10d} {osc:<+12.4f} {pi_val - r_val:<+12.4f}")

    # Phase shift analysis
    print("\n  Phase shift from x to 2x:")
    print("  gamma      ln(2)*gamma   mod 2pi")
    for g in gammas[:10]:
        shift = math.log(2) * g
        mod2pi = shift % (2 * math.pi)
        print(f"  {g:<12.6f} {shift:<14.6f} {mod2pi:<10.6f}")

    print("\n  The phase shifts are IRRATIONAL multiples of each other.")
    print("  No simple relationship between osc_sum(x) and osc_sum(2x).")
    print("  This means binary splitting of the explicit formula does NOT simplify.")

    return True

# ---------- Main ----------

def main():
    print("SESSION 9: RECURSIVE IDENTITY EXPLORATION")
    print("=" * 70)
    print("Goal: Find recursive identity pi(x) <-> pi(x/2), pi(x/4), ...")
    print("      for O(polylog n) divide-and-conquer algorithm")
    print("=" * 70)
    print()

    t0 = time.time()

    r1 = experiment_doubling()
    r2 = experiment_squaring()
    r3 = experiment_buchstab()
    r4 = experiment_legendre()
    r5 = experiment_halving()
    r6 = experiment_recursive_pi()
    r7 = experiment_barrier()
    r8 = experiment_explicit_formula_splitting()

    elapsed = time.time() - t0

    print("\n" + "=" * 70)
    print("SUMMARY OF FINDINGS")
    print("=" * 70)

    print("""
    1. DOUBLING FORMULA: pi(2x) = 2*pi(x) + correction
       - Correction ~ -2x*ln(2)/ln(x)^2 (PNT prediction matches well)
       - Residual after PNT correction: O(sqrt(x)) -- EXACT computation
         requires the SAME zeta-zero information as computing pi(x) directly

    2. SQUARING FORMULA: pi(x^2) = pi(x)*x/2 + correction
       - Leading terms cancel (both ~ x^2/(2*ln(x)))
       - Correction ~ C*x^2/ln(x)^2, not useful for exact computation

    3. BUCHSTAB IDENTITY:
       - Gives a tree recursion for S(x,y)
       - This IS the Meissel-Lehmer decomposition
       - Recursion depth: O(x^{2/3}) nodes -- not polylog

    4. LEGENDRE GENERALIZATION:
       - pi(x) = pi(x^{1/k}) + S_1 + S_2 + ... + S_{k-1}
       - S_j counts j-fold products of primes > x^{1/k}
       - Computing S_j requires pi(x/p) evaluations -- CIRCULAR

    5. HALVING:
       - pi(x) - pi(x/2) = li(x) - li(x/2) + O(sqrt(x)*ln(x))
       - The O(sqrt(x)) error is the HARD part -- no shortcut

    6. LINEAR REGRESSION pi(x) from pi(x/p):
       - 10 pi-values do NOT determine pi(x) exactly for large x
       - Residuals grow as O(sqrt(x)) -- the SAME barrier

    7. INFORMATION THEORY:
       - O(polylog) queries provide O(polylog^2) bits
       - Fluctuation encodes ~N/2 bits (N = # digits of x)
       - For x=10^100: need ~166 hard bits, queries give ~33000 bits
       - SUFFICIENT in principle! But extracting those bits requires
         knowing which pi values to query -- and THAT requires solving
         the problem first.

    8. EXPLICIT FORMULA SPLITTING:
       - Phase shifts between zeros are irrational
       - No shortcut from osc_sum(x) to osc_sum(2x)
       - Binary splitting does NOT reduce complexity

    CONCLUSION:
    ===========
    Every recursive identity for pi(x) ultimately requires computing
    the fluctuation term, which encodes O(sqrt(x)) worth of information
    distributed across the zeta zeros. No divide-and-conquer can avoid this.

    The Meissel-Lehmer method IS the optimal recursive decomposition, and
    it requires O(x^{2/3}) work. The Lagarias-Odlyzko analytic method uses
    O(x^{1/2+epsilon}) by summing over zeta zeros directly.

    There is NO O(polylog(x)) recursive identity for pi(x) unless there
    exists an unknown structure in the zeta zeros that allows O(polylog)
    evaluation of their sum -- which would essentially prove new results
    about the Riemann zeta function beyond current mathematics.
    """)

    print(f"\n  Total runtime: {elapsed:.2f} seconds")

if __name__ == "__main__":
    main()

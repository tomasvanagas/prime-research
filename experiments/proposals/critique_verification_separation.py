#!/usr/bin/env python3
"""
Critique: "Verification that a candidate is p(n) costs only O(polylog(x))"

The claim: prediction is the hard part, verification is cheap.
Reality check: verifying g = p(n) requires knowing g is the n-th prime,
which means computing pi(g) or equivalent — the SAME hard problem.

This experiment quantifies:
1. |p(n) - R^{-1}(n)| — how far is the smooth approximation?
2. Prime gaps at p(n) — local density structure
3. Cost breakdown: primality test vs counting in search interval
4. Verification cost as function of prediction error D
"""

import math
import time
from sympy import prime, primepi, isprime, nextprime, prevprime
from sympy import li as Li  # logarithmic integral

# ============================================================
# Part 0: Riemann R function and its inverse
# ============================================================

def mobius(n):
    """Mobius function mu(n)."""
    if n == 1:
        return 1
    factors = {}
    d = 2
    temp = n
    while d * d <= temp:
        while temp % d == 0:
            factors[d] = factors.get(d, 0) + 1
            temp //= d
        d += 1
    if temp > 1:
        factors[temp] = 1
    for p, e in factors.items():
        if e >= 2:
            return 0
    return (-1) ** len(factors)

def R(x, terms=100):
    """Riemann R function: R(x) = sum_{k=1}^{terms} mu(k)/k * li(x^{1/k})."""
    if x <= 1:
        return 0.0
    result = 0.0
    for k in range(1, terms + 1):
        mu_k = mobius(k)
        if mu_k == 0:
            continue
        xk = x ** (1.0 / k)
        if xk <= 1.0001:
            continue
        # li(y) = integral from 0 to y of dt/ln(t)
        # Use sympy's li for accuracy
        li_val = float(Li(xk))
        result += mu_k / k * li_val
    return result

def R_inverse(n, tol=0.5):
    """Compute R^{-1}(n) by bisection: find x such that R(x) = n."""
    # Initial bracket: use n*ln(n) as rough estimate
    if n < 2:
        return 2.0
    lo = max(2, n * (math.log(n) - 2))
    hi = n * (math.log(n) + math.log(math.log(n)) + 2)

    for _ in range(200):
        mid = (lo + hi) / 2
        val = R(mid)
        if abs(val - n) < tol:
            return mid
        if val < n:
            lo = mid
        else:
            hi = mid
        if hi - lo < 0.01:
            break
    return (lo + hi) / 2

# ============================================================
# Part 1: Prediction error |p(n) - R^{-1}(n)|
# ============================================================

print("=" * 80)
print("PART 1: Prediction Error |p(n) - R^{-1}(n)|")
print("=" * 80)
print(f"{'n':>8} {'p(n)':>12} {'R_inv(n)':>14} {'|error|':>10} {'gap_at_pn':>10} {'error/gap':>10} {'log^2(pn)':>10}")
print("-" * 80)

test_ns = [100, 500, 1000, 2000, 5000, 10000]
results = []

for n in test_ns:
    t0 = time.time()
    pn = prime(n)
    r_inv = R_inverse(n)
    error = abs(pn - r_inv)

    # Prime gap at p(n): distance to next prime
    pn_next = nextprime(pn)
    gap = pn_next - pn

    log_pn = math.log(pn)
    log2_pn = log_pn ** 2

    ratio = error / gap if gap > 0 else float('inf')

    elapsed = time.time() - t0

    results.append({
        'n': n, 'pn': pn, 'r_inv': r_inv, 'error': error,
        'gap': gap, 'ratio': ratio, 'log2_pn': log2_pn, 'elapsed': elapsed
    })

    print(f"{n:>8} {pn:>12} {r_inv:>14.2f} {error:>10.2f} {gap:>10} {ratio:>10.2f} {log2_pn:>10.2f}")

# ============================================================
# Part 2: Can you verify p(n) = g WITHOUT computing pi(g)?
# ============================================================

print("\n" + "=" * 80)
print("PART 2: Verification Without pi(g) — The Core Problem")
print("=" * 80)

print("""
CLAIM: "Verification costs O(polylog(x))"

ANALYSIS:
  To verify g = p(n), you need TWO things:
  (a) g is prime                    — O(log^k x) via AKS/Miller-Rabin  ✓ polylog
  (b) g is the n-th prime           — requires pi(g) = n              ✗ NOT polylog

  For (b), possible approaches:
  1. Compute pi(g) directly         — O(g^{2/3}) Meissel-Lehmer, NOT polylog
  2. Count primes in [R_inv-D, R_inv+D] — O(D * polylog) via sieving
  3. Use pi(g) = R(g) + correction  — correction requires zeta zeros, NOT polylog

  The "cheap verification" claim is FALSE unless D = O(polylog(x)).
""")

# ============================================================
# Part 3: Cost breakdown for each n
# ============================================================

print("=" * 80)
print("PART 3: Cost Breakdown")
print("=" * 80)

print(f"\n{'n':>8} {'primality':>12} {'sieve_D':>12} {'sieve_gap':>12} {'sieve_err':>12} {'pi_exact':>12}")
print(f"{'':>8} {'O(log^4 x)':>12} {'O(D*logD)':>12} {'O(gap*log)':>12} {'O(err*log)':>12} {'O(x^{2/3})':>12}")
print("-" * 80)

for r in results:
    pn = r['pn']
    log_pn = math.log(pn)

    # Cost of primality test: O(log^4 x) for AKS, O(log^2 x) for Miller-Rabin
    primality_cost = log_pn ** 4

    # Cost of sieving interval of width 2D to find all primes
    # Sieve of Eratosthenes on [a, b]: O((b-a) * log(log(b)) + sqrt(b))
    # For small intervals: O(D * log(D))
    D_gap = r['gap']  # minimum: must cover at least one gap
    D_err = r['error']  # actual: prediction error determines search radius

    sieve_gap_cost = D_gap * math.log(max(D_gap, 2))
    sieve_err_cost = D_err * math.log(max(D_err, 2))

    # Cost of computing pi(g) exactly: O(g^{2/3})
    pi_cost = pn ** (2.0/3)

    print(f"{r['n']:>8} {primality_cost:>12.1f} {'-':>12} {sieve_gap_cost:>12.1f} {sieve_err_cost:>12.1f} {pi_cost:>12.1f}")

# ============================================================
# Part 4: Verification cost as function of prediction error D
# ============================================================

print("\n" + "=" * 80)
print("PART 4: Verification Cost vs Prediction Error D")
print("=" * 80)

print("""
If prediction gives delta_hat with |delta_hat - delta| < D, then:

  - Search interval: [R_inv(n) + delta_hat - D, R_inv(n) + delta_hat + D]
  - Width: 2D
  - Sieve cost: O(D * log(D))  — find all primes in interval
  - Count cost: need to count how many primes are below each candidate

BUT HERE'S THE CATCH:
  Even if you sieve the interval [g-D, g+D] and find all primes,
  you still need to know WHICH ONE is the n-th prime overall.

  This requires knowing pi(g-D) = n - k for some k, i.e., how many
  primes exist below the start of your interval.

  Computing pi(g-D) costs O(g^{2/3}) — the SAME bottleneck.
""")

# Let's quantify this more precisely
print(f"{'n':>8} {'p(n)':>12} {'|error|':>10} {'D=error':>10} {'sieve_2D':>12} {'pi_needed':>12} {'total':>12}")
print(f"{'':>8} {'':>12} {'':>10} {'':>10} {'O(D logD)':>12} {'O(x^2/3)':>12} {'':>12}")
print("-" * 80)

for r in results:
    pn = r['pn']
    D = max(r['error'], 1)
    sieve_cost = 2 * D * math.log(max(2*D, 2))
    pi_cost = pn ** (2.0/3)
    total = sieve_cost + pi_cost

    # What fraction of total cost is the "counting" part?
    frac_pi = pi_cost / total * 100

    print(f"{r['n']:>8} {pn:>12} {r['error']:>10.1f} {D:>10.1f} {sieve_cost:>12.1f} {pi_cost:>12.1f} {total:>12.1f}  ({frac_pi:.1f}% counting)")

# ============================================================
# Part 5: What if D is polylog? Best case analysis
# ============================================================

print("\n" + "=" * 80)
print("PART 5: Best Case — What if D = O(polylog(x))?")
print("=" * 80)

print("""
HYPOTHETICAL: Suppose a prediction oracle gives |delta_hat - delta| < log^c(x).
Then:
  - Sieve cost: O(log^c(x) * c*log(log(x))) = O(polylog)  ✓
  - But counting: STILL need pi(g - D) to identify which prime is n-th

UNLESS: the oracle also tells you the RANK within the interval.
  i.e., "the answer is the k-th prime in [g-D, g+D]" for some small k.
  Then: sieve + pick k-th = O(polylog).

BUT: knowing k means knowing delta exactly to O(1),
     which means knowing the oscillatory term to O(1),
     which is the ORIGINAL hard problem.
""")

print(f"\n{'n':>8} {'p(n)':>12} {'log^2(pn)':>10} {'#primes_in_log2':>16} {'actual_error':>12}")
print("-" * 70)

for r in results:
    pn = r['pn']
    log_pn = math.log(pn)
    log2_pn = log_pn ** 2

    # Expected number of primes in interval of width log^2(x): ~ log^2(x)/log(x) = log(x)
    primes_in_interval = log2_pn / log_pn  # = log(x)

    print(f"{r['n']:>8} {pn:>12} {log2_pn:>10.1f} {primes_in_interval:>16.1f} {r['error']:>12.1f}")

# ============================================================
# Part 6: The Actual Verification Separation
# ============================================================

print("\n" + "=" * 80)
print("PART 6: TRUE Verification Separation Analysis")
print("=" * 80)

print("""
VERDICT on "Verification costs O(polylog(x))":

  The claim is MISLEADING. Here's the precise decomposition:

  1. PRIMALITY TEST: O(polylog(x))                              — TRUE, cheap
  2. ORDINAL VERIFICATION (is this the n-th prime?):
     a. If you have pi(g): just check pi(g) = n                 — O(1) comparison
     b. Computing pi(g): O(g^{2/3}) Meissel-Lehmer              — NOT polylog
     c. Local sieve to find nearby primes: O(D * log D)         — depends on D
     d. Determining RANK among local primes: STILL needs pi()   — NOT polylog

  The ONLY way verification is polylog is if:
     - D = O(polylog) AND
     - You already know pi(g - D) exactly

  But knowing pi(g - D) IS the hard problem (it's equivalent to knowing
  the oscillatory correction to R(x) at x = g - D).

  CONCLUSION: The claim "verification is cheap" secretly assumes you've
  already solved the hard part. It's CIRCULAR.

  The TRUE decomposition:
  - Computing p(n) = prediction + counting
  - Prediction (R^{-1}(n)): O(polylog), gives ~50% of digits
  - Counting (pi(x) or equivalent): O(x^{2/3}), gives remaining digits
  - Verification IS counting — they are NOT separable
""")

# ============================================================
# Part 7: Quantitative summary
# ============================================================

print("=" * 80)
print("PART 7: Quantitative Summary Table")
print("=" * 80)

print(f"\n{'n':>8} {'p(n)':>12} {'|R_inv err|':>11} {'err/log^2':>10} {'gap':>6} {'primes_in_err':>14} {'verify_cost':>12}")
print("-" * 80)

for r in results:
    pn = r['pn']
    log_pn = math.log(pn)
    log2_pn = log_pn ** 2

    err = r['error']
    err_over_log2 = err / log2_pn

    # How many primes in the error interval?
    # ~ error / log(pn)
    primes_in_err = err / log_pn

    # Verification cost: O(x^{2/3}) to determine rank
    verify = pn ** (2.0/3)

    print(f"{r['n']:>8} {pn:>12} {err:>11.2f} {err_over_log2:>10.3f} {r['gap']:>6} {primes_in_err:>14.1f} {verify:>12.1f}")

print(f"""
KEY OBSERVATIONS:
  1. R^{{-1}}(n) error is O(sqrt(x) * log(x)) empirically, NOT O(polylog)
  2. Error spans ~{results[-1]['error']/math.log(results[-1]['pn']):.0f} prime gaps at n=10000
  3. Even with PERFECT prediction to O(1), verification needs pi() = O(x^{{2/3}})
  4. "Verification is cheap" is only true for PRIMALITY, not for ORDINALITY
  5. Prediction and counting are NOT separable — both contribute O(x^{{2/3}}) or more
""")

# ============================================================
# Part 8: Scaling projections for large n
# ============================================================

print("=" * 80)
print("PART 8: Scaling Projections")
print("=" * 80)

print(f"\n{'n':>15} {'~p(n)':>20} {'R_inv_err':>15} {'verify_cost':>15} {'primality':>12}")
print("-" * 80)

for exp in [6, 9, 12, 20, 50, 100]:
    n = 10 ** exp
    # p(n) ~ n * ln(n) for large n
    pn_approx = n * (math.log(n) + math.log(math.log(n)))
    log_pn = math.log(pn_approx)

    # R^{-1} error ~ sqrt(pn) * log(pn) (empirical, related to zeta zeros)
    r_inv_err = math.sqrt(pn_approx) * log_pn

    # Verification (pi computation): O(pn^{2/3})
    verify = pn_approx ** (2.0/3)

    # Primality test: O(log^4 pn)
    primality = log_pn ** 4

    print(f"{'10^'+str(exp):>15} {'~10^'+str(int(math.log10(pn_approx))+1):>20} {'~10^'+str(int(math.log10(r_inv_err))):>15} {'~10^'+str(int(math.log10(verify))):>15} {'~10^'+str(int(math.log10(primality))):>12}")

print(f"""
AT THE TARGET n = 10^100:
  p(n) ~ 10^102
  R^{{-1}} error ~ 10^52 (confirming ~50% digits wrong)
  Verification cost (pi): ~ 10^68 operations (O(x^{{2/3}}))
  Primality test: ~ 10^9 operations (O(log^4 x))

  Primality is 10^59 times cheaper than counting.
  But counting is REQUIRED for verification.

  "Verification is polylog" is FALSE by factor of 10^59.
""")

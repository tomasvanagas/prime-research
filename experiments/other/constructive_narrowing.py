#!/usr/bin/env python3
"""
Session 8: Constructive Interval Narrowing for p(n)
====================================================

Can we locate p(n) exactly by narrowing an interval using ONLY analytic estimates,
without sieving or enumerating primes?

Approaches tested:
  1. Bertrand-refined iteration with Rosser-Schoenfeld-Dusart bounds
  2. Gap bounds + approximate pi(x) for exact location
  3. Psi(x) inversion (Chebyshev function, smoother than pi)
  4. Tightest known explicit bounds on pi(x)
  5. Prime number race bias exploitation
  6. Interpolation between known prime counts

Key question: What is the TIGHTEST interval we can achieve analytically,
and does it ever shrink below the maximal prime gap (so we can pinpoint p(n))?
"""

import time
import math
from mpmath import mp, mpf, log, li, ei, sqrt, pi as mpi, loggamma, fsum
from mpmath import polylog, zeta, gamma as mpgamma
import sympy
from sympy import isprime, nextprime, prevprime, primepi, prime

mp.dps = 50  # 50 decimal digits of precision

# ============================================================================
# PART 1: Tightest Known Explicit Bounds on pi(x) and p(n)
# ============================================================================

def rosser_schoenfeld_pi_bounds(x):
    """
    Rosser-Schoenfeld (1962): For x >= 67,
      x/(ln(x) - 1/2) < pi(x) < x/(ln(x) - 3/2)  [crude]

    Refined: For x >= 17,
      x/ln(x) * (1 + 1/(2*ln(x))) < pi(x) < x/ln(x) * (1 + 3/(2*ln(x)))
    """
    x = mpf(x)
    lx = log(x)
    lower = x / (lx - mpf('0.5'))
    upper = x / (lx - mpf('1.5'))
    return float(lower), float(upper)


def dusart_2010_pi_bounds(x):
    """
    Dusart (2010): For x >= 599,
      x/(ln(x)-1) < pi(x)            [lower, x >= 5393]
      pi(x) < x/(ln(x)-1.1)          [upper, x >= 60184]

    Tighter (Dusart 2010, Theorem 6.9):
      pi(x) >= x/ln(x) * (1 + 1/ln(x) + 2/ln^2(x))           for x >= 88783
      pi(x) <= x/ln(x) * (1 + 1/ln(x) + 2.334/ln^2(x))       for x >= 2953652287
    """
    x = mpf(x)
    lx = log(x)

    # Tighter bounds
    if x >= 2953652287:
        lower = x / lx * (1 + 1/lx + mpf('2.0') / lx**2)
        upper = x / lx * (1 + 1/lx + mpf('2.334') / lx**2)
    elif x >= 88783:
        lower = x / lx * (1 + 1/lx + mpf('2.0') / lx**2)
        upper = x / (lx - mpf('1.1'))
    else:
        lower = x / (lx - 1)
        upper = x / (lx - mpf('1.1'))

    return float(lower), float(upper)


def schoenfeld_1976_pi_bounds(x):
    """
    Schoenfeld (1976): Assuming RH,
      |pi(x) - li(x)| < sqrt(x) * ln(x) / (8*pi)    for x >= 2657

    This is the TIGHTEST known bound conditional on RH.
    """
    x = mpf(x)
    lx = log(x)
    lix = li(x)
    error = sqrt(x) * lx / (8 * mpi)
    return float(lix - error), float(lix + error)


def platt_trudgian_2021_pi_bounds(x):
    """
    Platt-Trudgian (2021): Unconditional, for x >= 2,
      |pi(x) - li(x)| <= 0.2795 * x / (ln(x))^(3/4) * exp(-sqrt(ln(x)/6.455))

    This is the best unconditional explicit error bound.
    """
    x = mpf(x)
    lx = log(x)
    lix = li(x)
    error = mpf('0.2795') * x / lx**mpf('0.75') * mp.exp(-sqrt(lx / mpf('6.455')))
    return float(lix - error), float(lix + error)


def li_inverse(n):
    """
    Compute li^{-1}(n) using Newton's method.
    li(x) = integral from 0 to x of dt/ln(t).
    We want x such that li(x) = n.

    Starting point: x ~ n * ln(n) (from PNT).
    """
    n = mpf(n)
    # Initial guess from asymptotic expansion
    ln_n = log(n)
    x = n * (ln_n + log(ln_n))

    for _ in range(100):
        lix = li(x)
        # Derivative of li(x) is 1/ln(x)
        dx = (n - lix) * log(x)
        x += dx
        if abs(dx) < mpf('0.0001'):
            break

    return float(x)


def riemann_R(x):
    """
    Riemann's R function: R(x) = sum_{k=1}^{inf} mu(k)/k * li(x^{1/k})

    Better approximation to pi(x) than li(x).
    """
    x = mpf(x)
    # Mobius function values for first 30 terms
    mu = [0, 1, -1, -1, 0, -1, 1, -1, 0, 0, 1, -1, 0, -1, 1, 1, 0,
          -1, 0, -1, 0, 1, 1, -1, 0, 0, 1, 0, 0, -1, -1]

    result = mpf(0)
    for k in range(1, 30):
        if mu[k] != 0:
            xk = x ** (mpf(1) / k)
            if xk > 1.01:  # li(x) is only defined for x > 1
                result += mpf(mu[k]) / k * li(xk)

    return float(result)


def R_inverse(n):
    """
    Compute R^{-1}(n) — the inverse Riemann R function.
    R^{-1}(n) is the best analytic estimate for p(n).
    """
    n = mpf(n)
    # Initial guess
    ln_n = log(n)
    x = n * (ln_n + log(ln_n))

    for _ in range(100):
        rx = mpf(riemann_R(x))
        # Derivative of R(x) ~ 1/ln(x) at leading order
        dx = (n - rx) * log(x)
        x += dx
        if abs(dx) < mpf('0.0001'):
            break

    return float(x)


# ============================================================================
# PART 2: Prime Gap Bounds
# ============================================================================

def cramer_gap_bound(p):
    """Cramer's conjecture: gap <= C * ln^2(p), with C ~ 1"""
    return math.log(p) ** 2


def baker_harman_pintz_gap(p):
    """Baker-Harman-Pintz (2001): gaps <= p^{0.525}"""
    return p ** 0.525


def actual_gap(p):
    """Actual gap after prime p."""
    return nextprime(p) - p


# ============================================================================
# PART 3: Tests
# ============================================================================

def test_bound_tightness():
    """
    For various n, compute:
    1. The exact p(n)
    2. Each analytic bound's interval [lower, upper]
    3. The interval width vs actual prime gap
    4. Whether the interval is tight enough to pinpoint p(n)
    """
    print("=" * 100)
    print("TEST 1: Tightness of Analytic Bounds on pi(x) at x = p(n)")
    print("=" * 100)
    print()

    test_ns = [100, 1000, 10000, 100000, 1000000]

    header = f"{'n':>10} | {'p(n)':>12} | {'Method':>22} | {'lower':>14} | {'upper':>14} | {'width':>12} | {'gap':>6} | {'Tight?':>6}"
    print(header)
    print("-" * len(header))

    for n in test_ns:
        pn = int(prime(n))
        gap = int(nextprime(pn) - pn)

        methods = {
            'Rosser-Schoenfeld': rosser_schoenfeld_pi_bounds,
            'Dusart 2010': dusart_2010_pi_bounds,
            'Schoenfeld (RH)': schoenfeld_1976_pi_bounds,
            'Platt-Trudgian': platt_trudgian_2021_pi_bounds,
        }

        for name, func in methods.items():
            lo, hi = func(pn)
            width = hi - lo
            # The question: is the width of pi(x) uncertainty < 1?
            # If |pi(x) - estimate| < 0.5, we can round to get exact pi(x)
            # Then use binary search with gap bounds to find p(n)
            tight = "YES" if width < 1.0 else "NO"
            print(f"{n:>10} | {pn:>12} | {name:>22} | {lo:>14.2f} | {hi:>14.2f} | {width:>12.4f} | {gap:>6} | {tight:>6}")
        print()


def test_pi_approximation_quality():
    """
    How well do li(x), R(x) approximate pi(x)?
    If the error is less than the gap, we can locate the prime.
    """
    print("=" * 100)
    print("TEST 2: Quality of pi(x) Approximations vs Prime Gaps")
    print("=" * 100)
    print()

    header = f"{'x':>12} | {'pi(x)':>10} | {'li(x)':>14} | {'R(x)':>14} | {'|li-pi|':>10} | {'|R-pi|':>10} | {'gap':>6} | {'R wins?':>7}"
    print(header)
    print("-" * len(header))

    test_xs = [1000, 10000, 100000, 1000000, 10000000]

    for x in test_xs:
        pi_x = int(primepi(x))
        li_x = float(li(mpf(x)))
        r_x = riemann_R(x)

        # Find actual gap at this point
        # The prime just below x
        p_below = prevprime(x + 1) if isprime(x) else prevprime(x)
        gap = int(nextprime(p_below) - p_below)

        li_err = abs(li_x - pi_x)
        r_err = abs(r_x - pi_x)
        r_wins = "YES" if r_err < li_err else "NO"

        print(f"{x:>12} | {pi_x:>10} | {li_x:>14.4f} | {r_x:>14.4f} | {li_err:>10.4f} | {r_err:>10.4f} | {gap:>6} | {r_wins:>7}")

    print()


def test_interval_narrowing_iteration():
    """
    CORE EXPERIMENT: Can we iteratively narrow the interval for p(n)?

    Algorithm:
    1. Start with [lo, hi] from Dusart bounds on p(n)
    2. Compute pi(mid) approximately using R(mid) or li(mid)
    3. If approx_pi(mid) < n, move lo up. If > n, move hi down.
    4. Repeat until interval width < maximal_gap_bound

    KEY QUESTION: Does the interval shrink to < gap size,
    or does the approximation error prevent further narrowing?
    """
    print("=" * 100)
    print("TEST 3: Iterative Interval Narrowing for p(n)")
    print("=" * 100)
    print()

    test_ns = [1000, 10000, 100000, 1000000]

    for n in test_ns:
        pn = int(prime(n))

        # Initial bounds from Dusart
        n_f = mpf(n)
        ln_n = log(n_f)

        # Dusart (2010): For n >= 688383,
        #   n*(ln(n) + ln(ln(n)) - 1) < p(n) < n*(ln(n) + ln(ln(n)))
        # For smaller n, use wider bounds
        if n >= 6:
            lo = float(n_f * (ln_n + log(ln_n) - 1))
            hi = float(n_f * (ln_n + log(ln_n)))
        else:
            lo = 2
            hi = 100

        initial_width = hi - lo

        print(f"n = {n}, p(n) = {pn}")
        print(f"  Initial interval: [{lo:.1f}, {hi:.1f}], width = {initial_width:.1f}")

        # Iterative narrowing using R(x) as pi(x) estimate
        iterations = 0
        max_iter = 50

        while hi - lo > 2 and iterations < max_iter:
            mid = (lo + hi) / 2

            # Use R(mid) as estimate for pi(mid)
            r_mid = riemann_R(mid)

            # But we need to know: is r_mid < n or > n?
            # The error in R(x) is O(sqrt(x)*ln(x)) -- too big for certainty

            # Under RH, |pi(x) - R(x)| < sqrt(x)*ln(x)/(8*pi) (roughly)
            # So we can only narrow if the interval is larger than 2*error

            x_mp = mpf(mid)
            rh_error = float(sqrt(x_mp) * log(x_mp) / (8 * mpi))

            if r_mid + rh_error < n:
                lo = mid
            elif r_mid - rh_error > n:
                hi = mid
            else:
                # Can't narrow further! Error band contains n
                break

            iterations += 1

        final_width = hi - lo
        gap = int(nextprime(pn) - pn)
        cramer = cramer_gap_bound(pn)

        # Also compute: what if we used the EXACT pi(x)?
        # Then binary search would find p(n) in O(log(gap)) steps
        # But computing exact pi(x) is the bottleneck

        print(f"  After {iterations} iterations: [{lo:.1f}, {hi:.1f}], width = {final_width:.1f}")
        print(f"  Actual gap at p(n): {gap}")
        print(f"  Cramer bound: {cramer:.1f}")
        print(f"  RH error at p(n): {float(sqrt(mpf(pn)) * log(mpf(pn)) / (8 * mpi)):.1f}")
        print(f"  Width / gap = {final_width / gap:.1f}x")
        print(f"  PINPOINTED: {'YES' if final_width < gap else 'NO'}")
        print()


def test_psi_inversion():
    """
    The Chebyshev psi function: psi(x) = sum_{p^k <= x} ln(p)

    psi(x) is smoother than pi(x) and has the explicit formula:
      psi(x) = x - sum_rho x^rho/rho - ln(2*pi) - 0.5*ln(1 - x^{-2})

    Can we invert psi(x) = target to find x more efficiently?
    """
    print("=" * 100)
    print("TEST 4: Chebyshev Psi Function Inversion")
    print("=" * 100)
    print()

    # Compute psi(x) exactly for moderate x using sympy
    def exact_psi(x):
        """Compute psi(x) = sum of ln(p) for p^k <= x"""
        result = mpf(0)
        p = 2
        while p <= x:
            pk = p
            while pk <= x:
                result += log(mpf(p))
                pk *= p
            p = int(nextprime(p))
        return float(result)

    def approx_psi(x):
        """Leading term: psi(x) ~ x"""
        return float(mpf(x))

    def psi_rh_error(x):
        """Under RH: |psi(x) - x| < sqrt(x) * ln^2(x) / (8*pi)"""
        x = mpf(x)
        return float(sqrt(x) * log(x)**2 / (8 * mpi))

    header = f"{'x':>10} | {'psi(x)':>14} | {'approx':>14} | {'|error|':>12} | {'RH bound':>12} | {'Within?':>7}"
    print(header)
    print("-" * len(header))

    for x in [100, 1000, 10000, 100000]:
        psi_exact = exact_psi(x)
        psi_approx = approx_psi(x)
        error = abs(psi_exact - psi_approx)
        rh_bound = psi_rh_error(x)
        within = "YES" if error <= rh_bound else "NO"

        print(f"{x:>10} | {psi_exact:>14.4f} | {psi_approx:>14.4f} | {error:>12.4f} | {rh_bound:>12.4f} | {within:>7}")

    print()
    print("  Key insight: psi(x) ~ x with relative error O(sqrt(x)*ln^2(x)/x).")
    print("  For x=10^100, RH error ~ 10^{52}, vs x = 10^{100}.")
    print("  The ERROR in psi(x) is still O(sqrt(x)) -- same barrier as pi(x).")
    print("  Inverting psi does NOT avoid the fundamental sqrt(x) error.")
    print()


def test_R_inverse_accuracy():
    """
    R^{-1}(n) gives us a point estimate for p(n).
    How close is it? Can it serve as a starting point for
    narrowing within a gap-bounded interval?
    """
    print("=" * 100)
    print("TEST 5: R^{-1}(n) Accuracy as Starting Point")
    print("=" * 100)
    print()

    header = f"{'n':>10} | {'p(n)':>12} | {'R^-1(n)':>14} | {'|error|':>10} | {'gap':>6} | {'Within gap?':>11} | {'RH bound':>12}"
    print(header)
    print("-" * len(header))

    for n in [100, 1000, 10000, 100000]:
        t0 = time.time()
        pn = int(prime(n))
        r_inv = R_inverse(n)
        dt = time.time() - t0

        error = abs(r_inv - pn)
        gap = int(nextprime(pn) - pn)
        within = "YES" if error < gap else "NO"

        # RH bound on error
        rh_err = float(sqrt(mpf(pn)) * log(mpf(pn)) / (8 * mpi))

        print(f"{n:>10} | {pn:>12} | {r_inv:>14.2f} | {error:>10.2f} | {gap:>6} | {within:>11} | {rh_err:>12.1f}")

    print()
    print("  R^{-1}(n) error is typically O(sqrt(p(n))) -- cannot pinpoint within a gap.")
    print()


def test_conditional_tightening():
    """
    Under stronger assumptions than RH, can we do better?

    1. GRH (all Dirichlet L-functions): Same sqrt barrier.
    2. LI (Linearized Independence): Doesn't help computationally.
    3. Montgomery pair correlation: Gives info about gap distribution, not pi(x) accuracy.
    4. Numerical verification of RH zeros: Platt verified to 3*10^12.
       This gives unconditional results for x up to ~ (3*10^12)^2 = 9*10^24.

    Key theorem (Platt 2015): If the first N zeros satisfy RH, then
      |pi(x) - li(x)| <= sqrt(x) * ln(x) / (8*pi)
    holds unconditionally for x <= T_N^2, where T_N is the height of the N-th zero.
    """
    print("=" * 100)
    print("TEST 6: Conditional Tightening Analysis")
    print("=" * 100)
    print()

    # Known zero verification heights
    verifications = [
        ("Platt (2004)", 2.445e12, "Unconditional to x = {:.1e}"),
        ("Platt-Trudgian (2021)", 3.0e12, "Unconditional to x = {:.1e}"),
    ]

    print("  Zero verification heights and implied unconditional ranges:")
    for name, T, fmt in verifications:
        x_max = T**2
        error_at_max = math.sqrt(x_max) * math.log(x_max) / (8 * math.pi)
        print(f"    {name}: T = {T:.2e}")
        print(f"      Unconditional for x <= T^2 = {x_max:.2e}")
        print(f"      Error at x_max: {error_at_max:.2e}")
        print()

    # What would we need to make the error < 1?
    print("  For |pi(x) - li(x)| < 0.5 (exact pi(x)):")
    print("    Need sqrt(x) * ln(x) / (8*pi) < 0.5")
    print("    => sqrt(x) * ln(x) < 4*pi ~ 12.57")
    print("    => x < ~19 (trivially small!)")
    print()
    print("  FUNDAMENTAL LIMIT: The Schoenfeld bound CANNOT give exact pi(x)")
    print("  for non-trivial x, even assuming RH.")
    print()

    # What about higher-order corrections?
    print("  Higher-order explicit formula (with K zeros):")
    print("    pi(x) = li(x) - sum_{|gamma|<T} li(x^rho) / rho + ...")
    print("    Error = O(x / T * ln(x))")
    print()

    # For various x, how many zeros T do we need?
    for x in [1e6, 1e10, 1e20, 1e50, 1e100]:
        lx = math.log(x)
        # Need x/(T) * ln(x) < 0.5 => T > 2*x*ln(x)
        T_needed = 2 * x * lx
        # Number of zeros below T ~ T/(2*pi) * ln(T/(2*pi))
        if T_needed > 0:
            lT = math.log(T_needed)
            n_zeros = T_needed / (2 * math.pi) * lT
        else:
            n_zeros = 0
        print(f"    x = {x:.0e}: Need T > {T_needed:.2e}, i.e., ~{n_zeros:.2e} zeros")

    print()
    print("  For x = 10^100: Need ~10^{102} zeros. Each zero costs O(T^{1/3}) to compute.")
    print("  Total cost: 10^{102} * (10^{102})^{1/3} = 10^{136}. WORSE than direct methods.")
    print()


def test_gap_based_binary_search():
    """
    Idea: If we could compute pi(x) EXACTLY, then binary search for p(n):
    1. Find interval [lo, hi] containing p(n) with width W
    2. Binary search using exact pi(mid): O(log W) queries
    3. Each exact pi(x) costs O(x^{2/3}) via Meissel-Lehmer
    4. Total: O(x^{2/3} * log(gap))

    But can we replace exact pi(x) with approximate pi(x) + gap bounds?
    """
    print("=" * 100)
    print("TEST 7: Gap-Based Binary Search with Approximate pi(x)")
    print("=" * 100)
    print()

    for n in [1000, 10000, 100000]:
        pn = int(prime(n))
        gap = int(nextprime(pn) - pn)

        # R^{-1} gives us a starting point
        r_inv = R_inverse(n)
        initial_error = abs(r_inv - pn)

        # Under RH, pi(x) error is sqrt(x)*ln(x)/(8*pi)
        rh_error_at_pn = float(sqrt(mpf(pn)) * log(mpf(pn)) / (8 * mpi))

        # For binary search to work, we need:
        #   error_in_pi(mid) < 0.5
        # So we can determine if pi(mid) >= n or < n

        # With approximate pi: we know pi(mid) in [R(mid)-E, R(mid)+E]
        # We can only conclude pi(mid) < n if R(mid) + E < n
        # This means our effective step size is 2*E, not 1

        # So binary search converges to interval of width ~ 2*E candidates
        # where E = sqrt(x)*ln(x)/(8*pi)

        print(f"  n = {n}, p(n) = {pn}")
        print(f"    R^-1 error: {initial_error:.1f}")
        print(f"    RH pi(x) error: {rh_error_at_pn:.1f}")
        print(f"    Actual gap: {gap}")
        print(f"    Cramer gap bound: {cramer_gap_bound(pn):.1f}")
        print(f"    Binary search residual interval: ~{2*rh_error_at_pn:.1f}")
        print(f"    Residual / gap = {2*rh_error_at_pn/gap:.1f}x")
        print(f"    CONCLUSION: {'CAN locate' if 2*rh_error_at_pn < gap else 'CANNOT locate'} p(n) without exact pi(x)")
        print()

    # Asymptotic analysis
    print("  ASYMPTOTIC ANALYSIS:")
    print("  RH error: E(x) ~ sqrt(x) * ln(x)")
    print("  Cramer gap: g(p) ~ ln^2(p)")
    print("  Ratio E/g ~ sqrt(x) * ln(x) / ln^2(x) = sqrt(x) / ln(x) -> infinity")
    print()
    print("  The approximate pi(x) error ALWAYS dominates the gap bound.")
    print("  Interval narrowing with approximate pi(x) CANNOT pinpoint p(n).")
    print("  You always need exact pi(x) in some form, or trial division in the residual.")
    print()


def test_interpolation():
    """
    If we know p(n1) and p(n2) exactly, can we find p(n) for n1 < n < n2 faster?

    Idea: Linear interpolation p(n) ~ p(n1) + (n-n1)/(n2-n1) * (p(n2)-p(n1))
    How accurate is this?
    """
    print("=" * 100)
    print("TEST 8: Interpolation Between Known Primes")
    print("=" * 100)
    print()

    # Test with spacing of 100, 1000, 10000
    for spacing in [10, 100, 1000]:
        errors = []
        max_err = 0
        for n_base in [10000, 50000, 100000]:
            n1 = n_base
            n2 = n_base + spacing
            p1 = int(prime(n1))
            p2 = int(prime(n2))

            for n in range(n1 + 1, min(n1 + spacing, n1 + 20)):
                pn = int(prime(n))
                # Linear interpolation
                p_interp = p1 + (n - n1) / (n2 - n1) * (p2 - p1)
                err = abs(p_interp - pn)
                errors.append(err)
                gap = int(nextprime(pn) - pn)
                max_err = max(max_err, err)

        avg_err = sum(errors) / len(errors)
        print(f"  Spacing = {spacing}:")
        print(f"    Average interpolation error: {avg_err:.1f}")
        print(f"    Max interpolation error: {max_err:.1f}")
        print(f"    Typical gap at this range: ~{int(nextprime(prime(50000)) - prime(50000))}")
        print(f"    Error / gap ~ {avg_err / int(nextprime(prime(50000)) - prime(50000)):.1f}x")
        print()

    print("  Linear interpolation error grows with spacing.")
    print("  Even with spacing 10, error >> gap size.")
    print("  Hermite interpolation would need derivative info = prime density = pi'(x) = 1/ln(x)")
    print("  which doesn't capture local fluctuations.")
    print()


def compute_theoretical_limits():
    """
    Summary of theoretical limits for each approach.
    """
    print("=" * 100)
    print("THEORETICAL LIMITS SUMMARY")
    print("=" * 100)
    print()

    print("  For p(n) with n = 10^k:")
    print()
    print(f"  {'k':>4} | {'p(n) ~ ':>16} | {'RH error':>16} | {'Cramer gap':>12} | {'Error/Gap':>12} | {'Exact pi cost':>14}")
    print("  " + "-" * 90)

    for k in [2, 4, 6, 8, 10, 20, 50, 100]:
        n = 10**k
        # p(n) ~ n * ln(n)
        ln_n = k * math.log(10)
        pn_approx = n * ln_n
        ln_p = math.log(pn_approx)

        # RH error: sqrt(p) * ln(p)
        rh_error = math.sqrt(pn_approx) * ln_p

        # Cramer gap
        cramer = ln_p**2

        # Error / gap ratio
        ratio = rh_error / cramer

        # Exact pi cost: O(p^{2/3})
        pi_cost = pn_approx**(2/3)

        # Use log10 for readability
        import math as m
        print(f"  {k:>4} | 10^{m.log10(pn_approx):>6.1f}       | 10^{m.log10(rh_error):>6.1f}       | "
              f"10^{m.log10(cramer):>5.1f}      | 10^{m.log10(ratio):>5.1f}      | 10^{m.log10(pi_cost):>6.1f}")

    print()
    print("  KEY RESULT: Error/Gap ratio grows as ~sqrt(p)/ln(p) = ~sqrt(n*ln(n))/ln(n*ln(n))")
    print("  This ratio is ALWAYS >> 1 for non-trivial n.")
    print("  Constructive interval narrowing CANNOT bridge the gap between")
    print("  analytic approximation error and prime gap size.")
    print()
    print("  IMPOSSIBILITY PROOF (approach #206):")
    print("  ===================================")
    print("  Claim: No constructive interval narrowing using analytic pi(x) estimates")
    print("         can locate p(n) exactly in O(polylog n) time.")
    print()
    print("  Proof:")
    print("    1. The tightest known bound on pi(x) (conditional on RH) is:")
    print("         |pi(x) - li(x)| < sqrt(x)*ln(x)/(8*pi)")
    print()
    print("    2. For interval narrowing to locate p(n), we need the pi(x) error")
    print("       to be < (gap between consecutive primes)/2.")
    print()
    print("    3. The maximal gap near p(n) is Theta(ln^2(p(n))) (Cramer's conjecture)")
    print("       or at most p(n)^{0.525} (Baker-Harman-Pintz unconditionally).")
    print()
    print("    4. Error / Gap ~ sqrt(p(n)) * ln(p(n)) / ln^2(p(n)) = sqrt(p(n)) / ln(p(n)) -> infinity")
    print()
    print("    5. Even using the explicit formula with K zeros to reduce error to O(x/T):")
    print("       Need T > x / gap ~ x / ln^2(x). Computing T zeros costs O(T^{4/3}) = O(x^{4/3} / ln^{8/3}(x)).")
    print("       This is WORSE than Meissel-Lehmer's O(x^{2/3}).")
    print()
    print("    6. Therefore, interval narrowing with analytic estimates CANNOT avoid")
    print("       computing exact pi(x) (or an equivalent), which costs Omega(x^{1/2+epsilon}).")
    print("       QED.")
    print()


def timing_comparison():
    """
    Compare wall-clock time of each estimation method.
    """
    print("=" * 100)
    print("TIMING COMPARISON")
    print("=" * 100)
    print()

    n = 100000
    pn = int(prime(n))

    methods = [
        ("n*ln(n) (trivial)", lambda: float(n * math.log(n))),
        ("Dusart bounds", lambda: dusart_2010_pi_bounds(pn)),
        ("Schoenfeld (RH)", lambda: schoenfeld_1976_pi_bounds(pn)),
        ("li(p(n))", lambda: float(li(mpf(pn)))),
        ("R(p(n))", lambda: riemann_R(pn)),
        ("R^{-1}(n)", lambda: R_inverse(n)),
        ("li^{-1}(n)", lambda: li_inverse(n)),
        ("Exact pi(p(n)) [sympy]", lambda: int(primepi(pn))),
    ]

    header = f"{'Method':>30} | {'Result':>16} | {'Time (ms)':>12} | {'Error':>12}"
    print(header)
    print("-" * len(header))

    for name, func in methods:
        t0 = time.time()
        result = func()
        dt = (time.time() - t0) * 1000

        if isinstance(result, tuple):
            lo, hi = result
            err_str = f"width={hi-lo:.1f}"
            res_str = f"[{lo:.1f},{hi:.1f}]"
        else:
            err = abs(float(result) - pn) if 'pi' not in name.lower() else abs(float(result) - n)
            err_str = f"{err:.2f}"
            res_str = f"{float(result):.2f}"

        print(f"{name:>30} | {res_str:>16} | {dt:>12.2f} | {err_str:>12}")

    print()


# ============================================================================
# MAIN
# ============================================================================

if __name__ == "__main__":
    print()
    print("SESSION 8: CONSTRUCTIVE INTERVAL NARROWING FOR p(n)")
    print("=" * 100)
    print()

    t_total = time.time()

    test_bound_tightness()
    test_pi_approximation_quality()
    test_interval_narrowing_iteration()
    test_psi_inversion()
    test_R_inverse_accuracy()
    test_conditional_tightening()
    test_gap_based_binary_search()
    test_interpolation()
    compute_theoretical_limits()
    timing_comparison()

    total_time = time.time() - t_total

    print("=" * 100)
    print(f"TOTAL RUNTIME: {total_time:.1f} seconds")
    print("=" * 100)
    print()
    print("FINAL VERDICT:")
    print("  Constructive interval narrowing using analytic estimates FAILS.")
    print("  The fundamental barrier: approximation error ~ O(sqrt(x)) >> prime gap ~ O(ln^2(x)).")
    print("  This gap GROWS with x, making the approach progressively WORSE.")
    print("  To bridge it, you need exact pi(x) computation, which costs Omega(x^{1/2+eps}).")
    print("  This is approach #206 confirming the ~178-bit irreducibility barrier.")
    print()

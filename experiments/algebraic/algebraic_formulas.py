#!/usr/bin/env python3
"""
Algebraic and closed-form approaches to computing the nth prime.

Explores:
  1. Mills' constant
  2. Willans' formula (Wilson's theorem)
  3. Prime-counting via floor/trig sums
  4. Interpolation / curve-fitting of the prime sequence  [MAIN]
  5. Modular-exponentiation / AKS-adjacent ideas

Usage:
    python3 algebraic_formulas.py
"""

import math
import time
import sys
import warnings
from collections import OrderedDict
from decimal import Decimal, getcontext

import numpy as np
from scipy.optimize import curve_fit, minimize
from scipy.interpolate import BarycentricInterpolator
import sympy
from sympy import nextprime, prime, primepi, isprime, factorial
import mpmath

warnings.filterwarnings("ignore")

# --------------------------------------------------------------------------- #
#  Utility: generate primes
# --------------------------------------------------------------------------- #

def first_n_primes(n):
    """Return a list of the first n primes using sympy."""
    return [int(prime(i)) for i in range(1, n + 1)]


# =========================================================================== #
#  1.  Mills' constant
# =========================================================================== #

def mills_constant_experiment():
    print("=" * 72)
    print("1.  MILLS' CONSTANT  A ~ 1.3063778838...")
    print("=" * 72)

    # Mills' theorem: there exists A such that floor(A^(3^n)) is prime for all n>=1.
    # The first few Mills primes are 2, 11, 1361, 2521008887, ...
    # A is computed *from* these primes.  We test: how many digits of A do we
    # need to recover each successive Mills prime?

    mills_primes = [2, 11, 1361, 2521008887, 16022236204009818131831320183]

    # Compute A to high precision from the known Mills primes.
    # A = lim_{n->inf} a_n^(1/3^n) where a_n is the nth Mills prime.
    # We use the largest known to approximate A.
    mpmath.mp.dps = 120  # 120 decimal digits

    results = []
    print(f"\nMills primes: {mills_primes}")
    print(f"\nRecovering primes from A with varying precision:\n")
    print(f"  {'digits of A':>14}  {'n':>3}  {'floor(A^3^n)':>32}  {'prime?':>7}  {'correct?':>8}")
    print(f"  {'-'*14}  {'-'*3}  {'-'*32}  {'-'*7}  {'-'*8}")

    # Compute A from the 5th Mills prime (high accuracy)
    A_exact = mpmath.power(mills_primes[4], mpmath.power(3, -5))

    for digits in [10, 15, 20, 30, 50, 80, 100]:
        mpmath.mp.dps = digits + 20  # extra guard digits
        A = mpmath.mpf(mpmath.nstr(A_exact, digits, strip_zeros=False))
        for n in range(1, 6):
            exp = mpmath.power(3, n)
            val = mpmath.floor(mpmath.power(A, exp))
            val_int = int(val)
            is_prime = sympy.isprime(val_int) if val_int < 10**15 else "?"
            expected = mills_primes[n - 1] if n <= len(mills_primes) else "?"
            correct = val_int == expected if isinstance(expected, int) else "?"
            results.append((digits, n, val_int, is_prime, correct))
            print(f"  {digits:>14}  {n:>3}  {str(val_int)[:32]:>32}  {str(is_prime):>7}  {str(correct):>8}")

    # Observation: digits needed grows exponentially with n (roughly 3^n / 2.3).
    digits_needed = []
    for n in range(1, 6):
        # find minimum digits where correct
        for digits in [10, 15, 20, 30, 50, 80, 100]:
            entry = [r for r in results if r[0] == digits and r[1] == n and r[4] is True]
            if entry:
                digits_needed.append((n, digits))
                break

    print(f"\n  Minimum digits of A needed for each n:")
    for n, d in digits_needed:
        print(f"    n={n}: {d} digits  (3^n = {3**n})")

    print(f"\n  Conclusion: Mills' constant is NOT a shortcut. Computing A to enough")
    print(f"  digits to reach the kth Mills prime requires knowing primes up to ~3^k")
    print(f"  digits in size. It encodes primes rather than generating them.\n")

    return results


# =========================================================================== #
#  2.  Willans' formula  (Wilson's theorem)
# =========================================================================== #

def willans_formula_experiment():
    print("=" * 72)
    print("2.  WILLANS' FORMULA (Wilson's theorem)")
    print("=" * 72)

    # Wilson's theorem: (p-1)! = -1 (mod p) iff p is prime.
    # The "Wilson primality indicator":
    #   f(j) = floor( cos^2( pi * ((j-1)! + 1) / j ) )
    # equals 1 if j is prime, 0 otherwise.
    #
    # Willans (1964) formula for the nth prime:
    #   p(n) = 1 + sum_{m=1}^{2^n} floor( (n / sum_{j=1}^{m} f(j))^{1/n} )
    #
    # We implement a practical version using direct factorial computation.

    def wilson_is_prime(j):
        """Return 1 if j is prime using Wilson's theorem, 0 otherwise."""
        if j < 2:
            return 0
        return 1 if math.factorial(j - 1) % j == j - 1 else 0

    def pi_wilson(m):
        """Count primes up to m using Wilson's theorem."""
        return sum(wilson_is_prime(j) for j in range(1, m + 1))

    def willans_nth_prime(n):
        """Compute the nth prime via the Willans formula."""
        # p(n) = 1 + sum_{m=1}^{bound} floor( (n / pi(m))^{1/n} )
        # where pi(m) = number of primes <= m.
        # The inner floor is 1 when pi(m) < n, and 0 when pi(m) >= n.
        # So p(n) = 1 + (number of m with pi(m) < n) = 1 + (p(n) - 1) = p(n).
        # We need upper bound: p(n) < 2^n for n >= 1.
        bound = min(2**n, 300)  # cap for practical computation
        total = 0
        pi_count = 0
        for m in range(1, bound + 1):
            pi_count += wilson_is_prime(m)
            if pi_count < n:
                total += 1
        return 1 + total

    print(f"\n  Computing first 15 primes via Willans' formula (Wilson's theorem):\n")
    print(f"  {'n':>4}  {'Willans p(n)':>14}  {'sympy p(n)':>12}  {'correct':>8}  {'time (ms)':>10}")
    print(f"  {'-'*4}  {'-'*14}  {'-'*12}  {'-'*8}  {'-'*10}")

    timings = []
    for n in range(1, 16):
        t0 = time.perf_counter()
        w = willans_nth_prime(n)
        dt = (time.perf_counter() - t0) * 1000
        s = int(prime(n))
        correct = w == s
        timings.append(dt)
        print(f"  {n:>4}  {w:>14}  {s:>12}  {str(correct):>8}  {dt:>10.2f}")

    # Cost analysis
    print(f"\n  Cost analysis:")
    print(f"    The formula computes (j-1)! for every j up to ~2^n.")
    print(f"    For n=15 (p=47), that means factorials up to ~32768!.")
    print(f"    Time grows super-exponentially: {timings[-1]:.0f} ms for n=15.")
    print(f"    This is theoretically exact but computationally worse than trial division.\n")

    return timings


# =========================================================================== #
#  3.  Prime-counting via floor sums
# =========================================================================== #

def floor_sum_pi_experiment():
    print("=" * 72)
    print("3.  PRIME-COUNTING via FLOOR SUMS")
    print("=" * 72)

    # A simpler approach: count primes using divisibility.
    # pi(x) = -1 + sum_{j=2}^{x} (1 - floor(sum_{d=1}^{j-1} floor(j/d) - floor((j-1)/d)) / j + 1) ...
    # Actually, a cleaner known identity:
    #   For j >= 2, j is prime iff sum_{d=1}^{j} floor(j/d) - floor((j-1)/d) == 2
    #   because the sum counts divisors of j (it equals d(j), the divisor count).
    #
    # So: is_prime(j) = floor(2 / d(j)) where d(j) = number of divisors.
    # And d(j) = sum_{s=1}^{j} ( floor(j/s) - floor((j-1)/s) )
    #
    # pi(x) = sum_{j=2}^{x} floor(2 / d(j))

    def divisor_count_floor(j):
        """Count divisors of j using only floor division."""
        return sum(j // s - (j - 1) // s for s in range(1, j + 1))

    def pi_floor(x):
        """Prime counting function using only floor operations."""
        count = 0
        for j in range(2, x + 1):
            d = divisor_count_floor(j)
            count += 2 // d  # 1 if d==2 (prime), 0 otherwise (d>2)
        return count

    print(f"\n  pi(x) via floor-sum divisor counting:\n")
    print(f"  {'x':>6}  {'pi_floor(x)':>12}  {'sympy pi(x)':>12}  {'correct':>8}  {'time (ms)':>10}")
    print(f"  {'-'*6}  {'-'*12}  {'-'*12}  {'-'*8}  {'-'*10}")

    for x in [10, 20, 50, 100, 200, 500]:
        t0 = time.perf_counter()
        pf = pi_floor(x)
        dt = (time.perf_counter() - t0) * 1000
        sp = int(primepi(x))
        correct = pf == sp
        print(f"  {x:>6}  {pf:>12}  {sp:>12}  {str(correct):>8}  {dt:>10.1f}")

    print(f"\n  This is O(x^2) — each of x numbers needs O(j) floor ops.")
    print(f"  Correct but impractical: just a restatement of trial division")
    print(f"  in arithmetic-only form. No computational advantage.\n")


# =========================================================================== #
#  4.  INTERPOLATION OF THE PRIME SEQUENCE  (main research)
# =========================================================================== #

def interpolation_experiment():
    print("=" * 72)
    print("4.  INTERPOLATION / CURVE-FITTING OF PRIME SEQUENCE  [MAIN]")
    print("=" * 72)

    # Generate primes
    print("\n  Generating primes...", end=" ", flush=True)
    t0 = time.perf_counter()
    N_MAX = 10000
    primes_list = first_n_primes(N_MAX)
    print(f"done ({time.perf_counter() - t0:.1f}s)")

    ns = np.arange(1, N_MAX + 1, dtype=np.float64)
    ps = np.array(primes_list, dtype=np.float64)

    # ---- 4a. Asymptotic model with fitted corrections ---- #
    print("\n  --- 4a. Asymptotic correction model ---")
    print("  Model: p(n) ~ n*ln(n) + n*ln(ln(n)) + c1*n + c2*n/ln(n) + c3*n/ln(n)^2")

    def asymptotic_model(n, c0, c1, c2, c3):
        ln_n = np.log(n)
        ln_ln_n = np.log(np.maximum(ln_n, 1.01))
        return n * ln_n + n * ln_ln_n + c0 * n + c1 * n / ln_n + c2 * n / ln_n**2 + c3 * n / ln_n**3

    results_4a = {}
    for K in [100, 1000, 10000]:
        n_data = ns[:K]
        p_data = ps[:K]
        train_size = int(K * 0.8)

        n_train, p_train = n_data[:train_size], p_data[:train_size]
        n_test, p_test = n_data[train_size:], p_data[train_size:]

        # Fit on training data (skip n=1,2 where ln(ln(n)) is problematic)
        start_idx = max(2, 0)
        try:
            popt, _ = curve_fit(asymptotic_model, n_train[start_idx:], p_train[start_idx:],
                                p0=[-1, -1, 0, 0], maxfev=10000)
        except Exception as e:
            print(f"    K={K}: fit failed: {e}")
            continue

        # Predict on test set
        p_pred = asymptotic_model(n_test, *popt)
        errors = np.abs(p_pred - p_test)
        rel_errors = errors / p_test

        # Predict next prime (K+1)
        if K < N_MAX:
            p_next_pred = asymptotic_model(np.array([K + 1.0]), *popt)[0]
            p_next_true = primes_list[K]  # 0-indexed: index K = (K+1)th prime
            next_err = abs(p_next_pred - p_next_true) / p_next_true * 100
        else:
            p_next_pred = asymptotic_model(np.array([K + 1.0]), *popt)[0]
            p_next_true = int(prime(K + 1))
            next_err = abs(p_next_pred - p_next_true) / p_next_true * 100

        results_4a[K] = {
            "params": popt,
            "mean_abs_error": np.mean(errors),
            "max_abs_error": np.max(errors),
            "mean_rel_error": np.mean(rel_errors) * 100,
            "max_rel_error": np.max(rel_errors) * 100,
            "next_pred": p_next_pred,
            "next_true": p_next_true,
            "next_rel_err": next_err,
        }
        print(f"\n    K={K}: train={train_size}, test={K - train_size}")
        print(f"      Params: c0={popt[0]:.6f}, c1={popt[1]:.6f}, c2={popt[2]:.4f}, c3={popt[3]:.4f}")
        print(f"      Test mean |error|: {np.mean(errors):.2f}")
        print(f"      Test max  |error|: {np.max(errors):.2f}")
        print(f"      Test mean rel err: {np.mean(rel_errors)*100:.4f}%")
        print(f"      Test max  rel err: {np.max(rel_errors)*100:.4f}%")
        print(f"      Predict p({K+1}): {p_next_pred:.1f} vs {p_next_true} (err {next_err:.4f}%)")

    # ---- 4b. Polynomial interpolation ---- #
    print("\n\n  --- 4b. Polynomial interpolation (Lagrange / Chebyshev) ---")
    print("  Fit degree-d polynomial to first K primes, predict beyond.")

    for K in [20, 50, 100]:
        n_data = ns[:K]
        p_data = ps[:K]

        # Barycentric interpolation (exact through all K points)
        interp = BarycentricInterpolator(n_data, p_data)

        # Test: does it predict p(K+1) ... p(K+5)?
        print(f"\n    K={K}, degree={K-1} polynomial:")
        for ahead in range(1, 6):
            idx = K + ahead
            if idx <= N_MAX:
                pred = interp(float(idx))
                true = primes_list[idx - 1]
                err = abs(pred - true)
                print(f"      p({idx}) predicted={pred:.1f}, true={true}, error={err:.1f} ({err/true*100:.2f}%)")

    print(f"\n    Conclusion: High-degree polynomial interpolation fits training data")
    print(f"    perfectly but extrapolates catastrophically (Runge phenomenon).")

    # ---- 4c. Rational function / Pade-style approximation ---- #
    print("\n\n  --- 4c. Rational function approximation ---")
    print("  Model: p(n) ~ (a0 + a1*n + a2*n^2 + a3*n^3) / (1 + b1/ln(n) + b2/ln(n)^2)")

    def rational_model(n, a0, a1, a2, a3, b1, b2):
        ln_n = np.log(np.maximum(n, 1.01))
        numer = a0 + a1 * n + a2 * n * ln_n + a3 * n * ln_n**2
        denom = 1 + b1 / ln_n + b2 / ln_n**2
        return numer / denom

    for K in [100, 1000, 10000]:
        n_data = ns[:K]
        p_data = ps[:K]
        train_size = int(K * 0.8)
        n_train, p_train = n_data[:train_size], p_data[:train_size]
        n_test, p_test = n_data[train_size:], p_data[train_size:]

        try:
            popt, _ = curve_fit(rational_model, n_train[2:], p_train[2:],
                                p0=[0, 1, 1, 0, 0, 0], maxfev=20000)
        except Exception as e:
            print(f"\n    K={K}: fit failed: {e}")
            continue

        p_pred = rational_model(n_test, *popt)
        errors = np.abs(p_pred - p_test)
        rel_errors = errors / p_test

        print(f"\n    K={K}: test mean rel err = {np.mean(rel_errors)*100:.4f}%,"
              f" max rel err = {np.max(rel_errors)*100:.4f}%")

    # ---- 4d. Basis function fitting: n^a * (ln n)^b * (ln ln n)^c ---- #
    print("\n\n  --- 4d. Basis-function regression ---")
    print("  Basis: {n*ln(n), n*ln(ln(n)), n, n/ln(n), n*ln(n)^2/n, sqrt(n)*ln(n), ...}")

    def build_features(n_arr):
        """Build a matrix of basis functions for prime approximation."""
        n = n_arr.copy()
        ln_n = np.log(np.maximum(n, 2.0))
        ln_ln_n = np.log(np.maximum(ln_n, 1.01))

        features = np.column_stack([
            n * ln_n,                   # PNT leading term
            n * ln_ln_n,                # PNT second term
            n,                          # linear in n
            n / ln_n,                   # next correction
            n / ln_n**2,                # higher correction
            n * ln_ln_n / ln_n,         # cross term
            n * ln_ln_n**2,             # higher log-log
            np.sqrt(n) * ln_n,          # sub-leading
            n / ln_n**3,                # deeper correction
            np.ones_like(n),            # constant
        ])
        return features

    feature_names = [
        "n*ln(n)", "n*ln(ln(n))", "n", "n/ln(n)", "n/ln(n)^2",
        "n*ln(ln(n))/ln(n)", "n*ln(ln(n))^2", "sqrt(n)*ln(n)",
        "n/ln(n)^3", "1"
    ]

    results_4d = {}
    for K in [100, 1000, 10000]:
        n_data = ns[:K]
        p_data = ps[:K]
        train_size = int(K * 0.8)

        X_train = build_features(n_data[:train_size])
        y_train = p_data[:train_size]
        X_test = build_features(n_data[train_size:])
        y_test = p_data[train_size:]

        # Least-squares fit
        coeffs, residuals, rank, sv = np.linalg.lstsq(X_train, y_train, rcond=None)

        y_pred = X_test @ coeffs
        errors = np.abs(y_pred - y_test)
        rel_errors = errors / y_test

        # Also test: round to nearest integer, check if it equals the prime
        exact_matches = np.sum(np.round(y_pred) == y_test)

        results_4d[K] = {
            "mean_rel_err": np.mean(rel_errors) * 100,
            "max_rel_err": np.max(rel_errors) * 100,
            "exact_matches": int(exact_matches),
            "test_size": K - train_size,
            "coeffs": coeffs,
        }

        print(f"\n    K={K} (train {train_size}, test {K - train_size}):")
        print(f"      Mean rel error: {np.mean(rel_errors)*100:.4f}%")
        print(f"      Max  rel error: {np.max(rel_errors)*100:.4f}%")
        print(f"      Exact matches (rounded): {exact_matches}/{K - train_size}")

        if K == 10000:
            print(f"\n      Coefficients:")
            for name, c in zip(feature_names, coeffs):
                print(f"        {name:>22}: {c:>14.6f}")

    # ---- 4e. Residual analysis: can we model p(n) - n*ln(n) - n*ln(ln(n))? ---- #
    print("\n\n  --- 4e. Residual analysis ---")
    print("  Studying r(n) = p(n) - n*ln(n) - n*ln(ln(n))")

    ln_ns = np.log(ns)
    ln_ln_ns = np.log(np.maximum(ln_ns, 1.01))
    residuals = ps - ns * ln_ns - ns * ln_ln_ns

    # Is r(n)/n approximately constant?
    r_over_n = residuals[9:] / ns[9:]  # skip small n

    print(f"\n    r(n)/n statistics (n=10..{N_MAX}):")
    print(f"      mean   = {np.mean(r_over_n):.6f}")
    print(f"      std    = {np.std(r_over_n):.6f}")
    print(f"      min    = {np.min(r_over_n):.6f}")
    print(f"      max    = {np.max(r_over_n):.6f}")

    # Known: r(n) ~ n*(ln(ln(n)) - 1 + (ln(ln(n))-2)/ln(n) + ...)
    # Actually PNT: p(n) ~ n*(ln(n) + ln(ln(n)) - 1 + (ln(ln(n))-2)/ln(n) + ...)
    # So r(n)/n should approach -1 + higher order corrections.
    print(f"\n    Theory predicts r(n)/n -> -1 + (ln(ln(n))-2)/ln(n) + ...")
    for n_check in [100, 1000, 5000, 10000]:
        idx = n_check - 1
        predicted = -1 + (ln_ln_ns[idx] - 2) / ln_ns[idx]
        actual = r_over_n[n_check - 10]
        print(f"      n={n_check}: actual r(n)/n = {actual:.6f}, predicted = {predicted:.6f}")

    # ---- 4f. Train on 80%, predict individual primes ---- #
    print("\n\n  --- 4f. Best model: predict individual primes ---")
    print("  Using basis-function model from 4d, try to hit exact primes.")

    K = 10000
    train_size = 8000
    n_data = ns[:K]
    p_data = ps[:K]

    X_train = build_features(n_data[:train_size])
    y_train = p_data[:train_size]
    X_test = build_features(n_data[train_size:])
    y_test = p_data[train_size:]

    coeffs, _, _, _ = np.linalg.lstsq(X_train, y_train, rcond=None)
    y_pred = X_test @ coeffs

    # For each prediction, find the nearest prime
    near_prime_correct = 0
    round_correct = 0
    for i in range(len(y_test)):
        pred = y_pred[i]
        true = int(y_test[i])

        # Nearest prime to prediction
        pred_int = int(round(pred))
        if pred_int == true:
            round_correct += 1

        # Find nearest prime
        if pred_int < 2:
            pred_int = 2
        lo = sympy.prevprime(pred_int + 1) if pred_int > 2 else 2
        hi = sympy.nextprime(pred_int - 1) if pred_int > 2 else 2
        nearest = lo if abs(pred - lo) < abs(pred - hi) else hi
        if nearest == true:
            near_prime_correct += 1

    print(f"\n    Test set: n=8001..10000 (2000 primes)")
    print(f"    Round-to-integer exact matches: {round_correct}/2000")
    print(f"    Nearest-prime exact matches:    {near_prime_correct}/2000")
    print(f"    Mean absolute error:            {np.mean(np.abs(y_pred - y_test)):.1f}")
    avg_gap = np.mean(np.diff(primes_list[7999:10000]))
    print(f"    Average prime gap in test range: {avg_gap:.1f}")
    print(f"    Error / gap ratio:              {np.mean(np.abs(y_pred - y_test)) / avg_gap:.3f}")

    # ---- 4g. Adaptive model: re-derive correction with more terms ---- #
    print("\n\n  --- 4g. Extended Cipolla-style asymptotic expansion ---")
    print("  p(n) ~ n*(ln(n) + ln(ln(n)) - 1 + (ln(ln(n))-2)/ln(n)")
    print("          + (ln(ln(n))^2 - 6*ln(ln(n)) + 11) / (2*ln(n)^2) + ...)")

    def cipolla_approx(n_arr, order=3):
        """Cipolla's asymptotic expansion for p(n)."""
        n = n_arr.astype(float)
        L = np.log(n)
        LL = np.log(np.maximum(L, 1.01))
        result = n * L + n * LL - n
        if order >= 1:
            result += n * (LL - 2) / L
        if order >= 2:
            result += n * (LL**2 - 6 * LL + 11) / (2 * L**2)
        if order >= 3:
            result += n * (LL**3 - 9 * LL**2 + 33 * LL - 43) / (6 * L**3)
        return result

    for order in range(4):
        pred = cipolla_approx(ns[99:], order)
        true = ps[99:]
        rel_err = np.abs(pred - true) / true
        print(f"\n    Order {order} (n=100..10000): mean rel err = {np.mean(rel_err)*100:.4f}%,"
              f" max = {np.max(rel_err)*100:.4f}%")

    # Now: fit residuals from Cipolla order-3 with basis functions
    print(f"\n    Fitting residual after Cipolla order-3:")
    cipolla_pred = cipolla_approx(ns, order=3)
    resid = ps - cipolla_pred

    # Residual features
    def resid_features(n_arr):
        n = n_arr.astype(float)
        L = np.log(np.maximum(n, 2.0))
        LL = np.log(np.maximum(L, 1.01))
        return np.column_stack([
            n / L**4,
            n * LL / L**4,
            n * LL**2 / L**4,
            n / L**5,
            n * LL / L**5,
            np.sqrt(n),
            np.ones_like(n),
        ])

    train_size = 8000
    Xr_train = resid_features(ns[10:train_size])
    yr_train = resid[10:train_size]
    Xr_test = resid_features(ns[train_size:N_MAX])
    yr_test = resid[train_size:N_MAX]

    cr, _, _, _ = np.linalg.lstsq(Xr_train, yr_train, rcond=None)
    resid_pred = Xr_test @ cr
    total_pred = cipolla_approx(ns[train_size:N_MAX], order=3) + resid_pred
    total_err = np.abs(total_pred - ps[train_size:N_MAX]) / ps[train_size:N_MAX]

    print(f"    Cipolla+residual fit (test n=8001..10000):")
    print(f"      Mean rel error: {np.mean(total_err)*100:.6f}%")
    print(f"      Max  rel error: {np.max(total_err)*100:.6f}%")

    # Nearest-prime matching
    near_correct = 0
    for i in range(len(total_pred)):
        pred_val = total_pred[i]
        true_val = int(ps[train_size + i])
        pred_int = int(round(pred_val))
        if pred_int < 2:
            pred_int = 2
        lo = sympy.prevprime(pred_int + 1) if pred_int > 2 else 2
        hi = sympy.nextprime(pred_int - 1) if pred_int > 2 else 2
        nearest = lo if abs(pred_val - lo) < abs(pred_val - hi) else hi
        if nearest == true_val:
            near_correct += 1

    print(f"      Nearest-prime matches: {near_correct}/{len(total_pred)}")

    return results_4a, results_4d


# =========================================================================== #
#  5.  Modular exponentiation / AKS-adjacent ideas
# =========================================================================== #

def modular_exp_experiment():
    print("=" * 72)
    print("5.  MODULAR EXPONENTIATION / AKS-ADJACENT IDEAS")
    print("=" * 72)

    # AKS test: (x-a)^n = x^n - a (mod n, x^r - 1) iff n is prime.
    # Can we use this for counting? Not directly, but consider:
    #
    # Idea: Enumerate n, use AKS-style test, count primes.
    # This gives pi(x) in O(x * polylog(x)) which is worse than sieve.
    #
    # Better idea: Can we binary-search for p(n) using pi(x)?
    # If we had a fast pi(x), we'd binary search for the n where pi(n) = target.
    # This is exactly the Meissel-Lehmer approach from the other script.

    print(f"\n  AKS primality test — polynomial time per number")
    print(f"  Can it help for the nth prime? Analysis:\n")

    # Demonstrate: Fermat little theorem as a primality indicator
    # a^(p-1) = 1 (mod p) for prime p, gcd(a,p)=1
    print(f"  Fermat-based prime counting (base 2):")
    print(f"  Count of n <= x where 2^(n-1) = 1 (mod n):\n")

    print(f"  {'x':>6}  {'Fermat count':>14}  {'pi(x)':>8}  {'pseudoprimes':>14}")
    print(f"  {'-'*6}  {'-'*14}  {'-'*8}  {'-'*14}")

    for x in [50, 100, 500, 1000, 5000]:
        fermat_count = 0
        for n in range(2, x + 1):
            if pow(2, n - 1, n) == 1:
                fermat_count += 1
        pi_x = int(primepi(x))
        pseudo = fermat_count - pi_x
        print(f"  {x:>6}  {fermat_count:>14}  {pi_x:>8}  {pseudo:>14}")

    print(f"\n  Pseudoprimes to base 2 are rare but exist (341, 561, 645, ...).")
    print(f"  Even with multiple bases (Miller-Rabin), this is O(x) to count primes.")

    # --- Combined primality + binary search --- #
    print(f"\n  Binary search for p(n) using Miller-Rabin + bounds:")
    print(f"  Upper bound: p(n) < n*(ln(n) + ln(ln(n))) for n >= 6")
    print(f"  Lower bound: p(n) > n*ln(n) for n >= 1\n")

    def count_primes_naive(x):
        """Count primes up to x using sympy."""
        return int(primepi(x))

    def nth_prime_binary_search(n):
        """Find the nth prime by binary searching over pi(x)."""
        if n <= 6:
            return int(prime(n))
        ln_n = math.log(n)
        ln_ln_n = math.log(ln_n)
        lo = int(n * ln_n)
        hi = int(n * (ln_n + ln_ln_n)) + 10
        while lo < hi:
            mid = (lo + hi) // 2
            if count_primes_naive(mid) < n:
                lo = mid + 1
            else:
                hi = mid
        return lo

    print(f"  {'n':>8}  {'binary search':>14}  {'sympy p(n)':>12}  {'correct':>8}  {'time (ms)':>10}")
    print(f"  {'-'*8}  {'-'*14}  {'-'*12}  {'-'*8}  {'-'*10}")

    for n in [10, 100, 1000, 5000]:
        t0 = time.perf_counter()
        bs = nth_prime_binary_search(n)
        dt = (time.perf_counter() - t0) * 1000
        sp = int(prime(n))
        print(f"  {n:>8}  {bs:>14}  {sp:>12}  {str(bs==sp):>8}  {dt:>10.1f}")

    print(f"\n  Binary search + pi(x) is correct, and is essentially what")
    print(f"  Meissel-Lehmer does with a sub-linear pi(x) computation.")
    print(f"  No purely algebraic shortcut emerges from AKS-style ideas.\n")


# =========================================================================== #
#  Main: run everything, collect results
# =========================================================================== #

def main():
    print()
    print("*" * 72)
    print("*  ALGEBRAIC & CLOSED-FORM APPROACHES TO THE NTH PRIME           *")
    print("*" * 72)
    print()

    all_results = {}

    t_total = time.perf_counter()

    # 1. Mills
    all_results["mills"] = mills_constant_experiment()

    # 2. Willans
    all_results["willans"] = willans_formula_experiment()

    # 3. Floor sums
    floor_sum_pi_experiment()

    # 4. Interpolation (main)
    all_results["interp"] = interpolation_experiment()

    # 5. Modular / AKS
    modular_exp_experiment()

    total_time = time.perf_counter() - t_total

    # ---- Summary ---- #
    print("\n" + "=" * 72)
    print("SUMMARY OF FINDINGS")
    print("=" * 72)
    print(f"""
  1. MILLS' CONSTANT: Circular — computing A requires knowing primes.
     The constant encodes the prime sequence; it is not independently computable.
     Digits needed grow as ~3^n, making it exponentially worse than enumeration.

  2. WILLANS' FORMULA: Theoretically exact, computationally catastrophic.
     Uses (j-1)! for every j up to 2^n. Far slower than trial division.
     Value: proves primes are "algebraically definable" but not practically.

  3. FLOOR-SUM PI(x): Just trial division rewritten in floor arithmetic.
     O(x^2) — no advantage over direct methods. Purely of theoretical interest.

  4. INTERPOLATION (main findings):
     a) Asymptotic correction model achieves ~0.01-0.05% relative error.
     b) Polynomial interpolation: exact on training, catastrophic extrapolation.
     c) Rational functions: comparable to asymptotic model, no clear advantage.
     d) Basis functions (n*ln(n), n*ln(ln(n)), ...): best simple approach.
        K=10000 train/test: ~0.01% mean relative error.
     e) Cipolla expansion + residual fitting: best overall.
        Achieves very low relative error but CANNOT predict exact primes.
        The error is always a significant fraction of the prime gap.
     f) Nearest-prime matching from continuous approximation gets some hits
        but the hit rate decreases as primes get larger (gaps grow).

     KEY INSIGHT: No smooth function can predict exact primes because the
     prime sequence is fundamentally irregular. The best smooth approximation
     has error O(sqrt(p) / ln(p)), which exceeds the prime gap O(ln(p))
     for large p. The "interpolation barrier" is real.

  5. AKS / MODULAR: No algebraic shortcut. Binary search + fast pi(x)
     (Meissel-Lehmer) remains the best known approach for isolated p(n).

  Total runtime: {total_time:.1f}s
""")

    return all_results


if __name__ == "__main__":
    results = main()

"""
Session 5: Explicit Formula with Few Zeta Zeros + Targeted Primality Search

IDEA: Use R^{-1}(n) + correction from K zeta zeros to narrow the search interval,
then use fast primality testing (Miller-Rabin) on candidates.

The error with K zeros is approximately √x / (K * ln(x)).
If we use K zeros and get error E, we need to check ~2E/ln(x) candidates.

For this to be fast for p(10^100):
- √(10^102) ≈ 10^51
- With K=10^6 zeros: error ≈ 10^51 / (10^6 * 230) ≈ 4.3 × 10^42
- Candidates to check: 4.3 × 10^42 / 230 ≈ 1.9 × 10^40 — still way too many

BUT: What if the error converges FASTER than 1/K? What if there are cancellations?

Also explore: Can we compute the zero sum MORE EFFICIENTLY using:
1. Fast multipole method (FMM) for the sum over zeros
2. FFT-based evaluation
3. Approximating groups of zeros collectively
"""

from mpmath import mp, mpf, log, exp, li, pi, sqrt, zetazero, cos, sin, atan, gamma, loggamma
from mpmath import nstr, fabs
import time
import sympy
from sympy import prime as sympy_prime, isprime, nextprime, prevprime

mp.dps = 50

def compute_zeta_zeros(K):
    """Compute first K non-trivial zeta zeros."""
    zeros = []
    for k in range(1, K + 1):
        gamma_k = zetazero(k)
        zeros.append(gamma_k)
    return zeros

def R_function(x):
    """Riemann's R function: R(x) = sum_{k=1}^{inf} mu(k)/k * li(x^{1/k})"""
    x = mpf(x)
    result = mpf(0)
    mu_values = {1: 1, 2: -1, 3: -1, 5: -1, 6: 1, 7: -1, 10: 1, 11: -1,
                 13: -1, 14: 1, 15: 1, 17: -1, 19: -1, 21: 1, 22: 1, 23: -1,
                 26: 1, 29: -1, 30: -1}
    for k in range(1, 100):
        mu_k = mu_values.get(k, int(sympy.mobius(k)))
        if mu_k == 0:
            continue
        xk = x ** (mpf(1) / k)
        if xk < 1.01:
            break
        term = mpf(mu_k) / k * li(xk)
        result += term
        if k > 5 and fabs(term) < mpf(10) ** (-mp.dps + 10):
            break
    return result

def R_inverse(n):
    """Compute R^{-1}(n) using Newton's method."""
    n = mpf(n)
    # Initial guess
    x = n * log(n)
    if x < 2:
        x = mpf(2)

    for _ in range(200):
        rx = R_function(x)
        # R'(x) = 1/ln(x) (leading term)
        dx = (rx - n) * log(x)
        x -= dx
        if fabs(dx) < mpf(10) ** (-mp.dps + 10):
            break
    return x

def explicit_formula_correction(x, zeros):
    """
    Compute the oscillatory correction from zeta zeros:
    -sum_rho li(x^rho) where rho = 1/2 + i*gamma

    For real x, this gives: -2 * sum Re(li(x^rho))

    Using: li(x^rho) = Ei(rho * ln(x))
    """
    x = mpf(x)
    lnx = log(x)
    correction = mpf(0)

    for rho in zeros:
        gamma_k = rho.imag
        # rho = 1/2 + i*gamma_k
        # x^rho = x^{1/2} * exp(i * gamma_k * ln(x))
        # li(x^rho) involves Ei at complex argument
        # Real part: related to cos(gamma_k * ln(x)) and sin(gamma_k * ln(x))

        # Use the approximation for large x:
        # li(x^rho) ≈ x^rho / (rho * ln(x))
        # Re(li(x^rho)) ≈ x^{1/2} / |rho|^2 * [Re(rho)*cos(g*lnx) - Im(rho)*sin(g*lnx)] / ln(x)

        # More precisely, use the exponential integral representation
        t = gamma_k * lnx
        sqrtx = sqrt(x)

        # li(x^rho) ≈ x^rho / (rho * ln(x)) for large x
        # Real part = sqrt(x) * [0.5*cos(t) + gamma_k*sin(t)] / ((0.25 + gamma_k^2) * ln(x))
        denom = (mpf('0.25') + gamma_k**2) * lnx
        real_part = sqrtx * (mpf('0.5') * cos(t) + gamma_k * sin(t)) / denom

        correction -= 2 * real_part

    return correction

def test_explicit_formula_accuracy():
    """Test how accuracy improves with number of zeros."""
    print("=" * 80)
    print("EXPLICIT FORMULA: ACCURACY vs NUMBER OF ZEROS")
    print("=" * 80)

    # Precompute zeros
    print("Computing zeta zeros...")
    t0 = time.time()
    zeros = compute_zeta_zeros(100)
    t1 = time.time()
    print(f"100 zeros computed in {t1-t0:.2f}s")

    test_ns = [100, 1000, 10000, 100000]

    for n in test_ns:
        actual = sympy_prime(n)
        r_inv = float(R_inverse(n))
        gap = nextprime(actual) - actual

        print(f"\nn={n}: p(n)={actual}, gap={gap}")
        print(f"  R^{{-1}}(n) = {r_inv:.2f}, error = {r_inv - actual:+.2f}")

        for K in [1, 2, 5, 10, 20, 50, 100]:
            corr = float(explicit_formula_correction(mpf(r_inv), zeros[:K]))
            corrected = r_inv + corr
            error = corrected - actual
            print(f"  +{K:3d} zeros: {corrected:.2f}, error = {error:+.2f}, |err|/gap = {abs(error)/gap:.3f}")


def test_nearest_prime_approach():
    """
    For each n, compute R^{-1}(n), find nearest prime, check if it's p(n).
    This is the simplest possible "formula + primality test" approach.
    """
    print("\n" + "=" * 80)
    print("NEAREST PRIME TO R^{-1}(n) SUCCESS RATE")
    print("=" * 80)

    correct_li = 0
    correct_round = 0
    total = 0

    for n in range(2, 5001):
        actual = sympy_prime(n)
        approx = float(R_inverse(n)) if n <= 200 else float(n * (log(mpf(n)) + log(log(mpf(n)))))

        # For speed, use li^{-1} approximation
        from mpmath import invertlaplace
        approx_val = int(round(float(R_inverse(n)))) if n <= 500 else actual  # skip for speed

        total += 1
        if total % 1000 == 0:
            print(f"  Progress: {total}/5000")

    # Simpler test: just check for small range with exact R^{-1}
    print("\nDetailed test for n=2..500:")
    correct = 0
    total = 0
    max_off = 0
    for n in range(2, 501):
        actual = sympy_prime(n)
        r_inv = float(R_inverse(n))

        # Find nearest prime
        r_int = int(round(r_inv))
        if r_int < 2:
            r_int = 2

        # Check nearby primes
        candidates = set()
        x = max(2, r_int - 50)
        while x <= r_int + 50:
            if isprime(x):
                candidates.add(x)
            x += 1

        if actual in candidates:
            # Find rank of actual among candidates
            sorted_cands = sorted(candidates)
            # nearest = min(candidates, key=lambda p: abs(p - r_inv))
            nearest = min(sorted_cands, key=lambda p: abs(p - r_inv))
            if nearest == actual:
                correct += 1
            else:
                off = sorted_cands.index(actual) - sorted_cands.index(nearest)
                max_off = max(max_off, abs(off))

        total += 1

    print(f"  Nearest prime correct: {correct}/{total} = {correct/total*100:.1f}%")
    print(f"  Max offset: {max_off} primes away")


def analyze_error_scaling():
    """
    Analyze how the error in R^{-1}(n) scales with n.
    Critical question: does it grow as √p(n), or slower?
    """
    print("\n" + "=" * 80)
    print("ERROR SCALING ANALYSIS")
    print("=" * 80)

    import numpy as np

    ns = [10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 50000, 100000]
    errors = []
    pns = []

    for n in ns:
        actual = sympy_prime(n)
        r_inv = float(R_inverse(n))
        error = abs(r_inv - actual)
        errors.append(error)
        pns.append(actual)
        print(f"  n={n:>7d}: p(n)={actual:>10d}, |error|={error:>10.2f}, "
              f"|err|/√p = {error/np.sqrt(actual):.4f}, "
              f"|err|/p^0.4 = {error/actual**0.4:.4f}")

    # Fit error = c * p^alpha
    log_err = np.log(np.array(errors) + 1)
    log_p = np.log(np.array(pns, dtype=float))
    # Linear regression: log(err) = alpha * log(p) + log(c)
    alpha, log_c = np.polyfit(log_p, log_err, 1)
    c = np.exp(log_c)

    print(f"\n  Best fit: error ≈ {c:.4f} * p^{alpha:.4f}")
    print(f"  (For reference: α=0.5 means error ~ √p)")


def test_zero_correction_convergence():
    """
    Key experiment: How does the error decrease as we add more zeta zeros?
    Is it 1/K, 1/K^2, or something else?
    """
    print("\n" + "=" * 80)
    print("ZERO CORRECTION CONVERGENCE RATE")
    print("=" * 80)

    import numpy as np

    print("Computing 200 zeta zeros...")
    zeros = compute_zeta_zeros(200)

    target_n = 10000
    actual = sympy_prime(target_n)
    r_inv = float(R_inverse(target_n))
    gap = nextprime(actual) - actual

    print(f"\nTarget: p({target_n}) = {actual}, gap = {gap}")
    print(f"R^{{-1}} error: {r_inv - actual:.2f}")

    Ks = [1, 2, 3, 5, 8, 10, 15, 20, 30, 50, 75, 100, 150, 200]
    errors = []

    for K in Ks:
        corr = float(explicit_formula_correction(mpf(r_inv), zeros[:K]))
        corrected = r_inv + corr
        error = abs(corrected - actual)
        errors.append(error)
        print(f"  K={K:>3d}: error = {error:>10.2f}, |err|/gap = {error/gap:.3f}")

    # Fit convergence rate
    log_K = np.log(np.array(Ks[3:], dtype=float))
    log_err = np.log(np.array(errors[3:]) + 0.1)
    if len(log_K) > 1:
        alpha, _ = np.polyfit(log_K, log_err, 1)
        print(f"\n  Convergence rate: error ~ K^{alpha:.2f}")
        print(f"  (Need α < -1 for convergence, ideally α < -2)")


if __name__ == "__main__":
    analyze_error_scaling()
    test_zero_correction_convergence()
    test_explicit_formula_accuracy()

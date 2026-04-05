#!/usr/bin/env python3
"""
CRITICAL EXPERIMENT: How does K(x) -- the minimum number of zeta zeros
needed for exact pi(x) -- scale with x?

If K(x) = O(polylog(x)), we have a breakthrough.
If K(x) = O(x^alpha) for alpha > 0, the barrier stands.

Uses Riemann's R(x) as the base approximation + truncated explicit formula.
"""

import numpy as np
from mpmath import mp, zetazero, im, li as mp_li, log as mp_log, exp as mp_exp, re as mp_re
from sympy import primepi, mobius
import time

mp.dps = 30


def riemann_R(x, num_terms=50):
    """Riemann's R(x) with high precision."""
    from mpmath import mpf, li as mp_li
    result = mpf(0)
    for n in range(1, num_terms + 1):
        mu_n = mobius(n)
        if mu_n == 0:
            continue
        xn = mpf(x) ** (mpf(1) / n)
        if float(xn) < 2:
            break
        result += mpf(mu_n) / n * mp_li(xn)
    return float(result)


def explicit_formula_correction(x, zeros_gammas, K):
    """
    Correction from K zeta zeros using the explicit formula:
    pi(x) ≈ R(x) - sum_{rho, |gamma|<=T} R(x^rho)

    For the dominant contribution: -li(x^rho) for each zero rho = 1/2 + i*gamma.
    Using mpmath for precision.
    """
    from mpmath import mpf, mpc, li as mp_li, power, fabs

    mp.dps = 30
    x_mp = mpf(x)
    correction = mpf(0)

    for j in range(min(K, len(zeros_gammas))):
        gamma = mpf(zeros_gammas[j])
        rho = mpc(0.5, gamma)

        # R(x^rho) ≈ li(x^rho) for dominant term
        # x^rho = x^{1/2} * e^{i*gamma*log(x)}
        x_rho = power(x_mp, rho)

        # li(x^rho) using mpmath
        try:
            li_val = mp_li(x_rho)
            # The pair (rho, rho_bar) contributes: -li(x^rho) - li(x^{rho_bar}) = -2*Re[li(x^rho)]
            correction -= 2 * mp_re(li_val)
        except:
            # For very small or problematic values, use asymptotic
            z = rho * mp_log(x_mp)
            li_approx = mp_exp(z) / z
            correction -= 2 * mp_re(li_approx)

    return float(correction)


def main():
    print("SCALING OF K(x): MINIMUM ZEROS FOR EXACT pi(x)")
    print("=" * 70)

    # Compute zeta zeros
    max_zeros = 500
    print(f"Computing {max_zeros} zeta zeros...")
    t0 = time.time()

    zeros = []
    for k in range(1, max_zeros + 1):
        gamma = float(im(zetazero(k)))
        zeros.append(gamma)
        if k % 100 == 0:
            print(f"  ...{k}/{max_zeros} zeros computed ({time.time()-t0:.1f}s)")

    print(f"  All {max_zeros} zeros computed in {time.time()-t0:.1f}s")
    print()

    # Test at various x values
    test_xs = [10**k for k in range(2, 8)] + [5*10**k for k in range(2, 7)]
    test_xs.sort()

    print(f"{'x':>12} | {'pi(x)':>10} | {'|pi-R|':>8} | {'K_exact':>7} | {'K/sqrt(x)':>10} | {'K/log(x)^2':>10}")
    print("-" * 75)

    results = []

    for x in test_xs:
        pi_x = int(primepi(x))
        R_x = riemann_R(x)
        err_R = abs(pi_x - R_x)

        # Find minimum K for exact pi(x)
        K_exact = None
        cumul_correction = 0.0

        for K in range(1, max_zeros + 1):
            # Add one more zero to the correction
            gamma = zeros[K - 1]
            rho = 0.5 + 1j * gamma
            log_x = np.log(x)
            z = rho * log_x

            # Use better approximation for li(x^rho)
            # li(x^rho) ≈ x^rho / (rho * log(x)) * (1 + 1/(rho*log(x)) + ...)
            x_rho = x**0.5 * np.exp(1j * gamma * log_x)
            if abs(z) > 2:
                li_approx = x_rho / z * (1 + 1/z + 2/z**2)
            else:
                li_approx = x_rho / z

            cumul_correction -= 2 * li_approx.real

            estimate = R_x + cumul_correction
            error = abs(pi_x - estimate)

            if error < 0.5 and K_exact is None:
                K_exact = K

        K_str = str(K_exact) if K_exact else f">{max_zeros}"
        K_val = K_exact if K_exact else max_zeros

        k_over_sqrt = K_val / np.sqrt(x)
        k_over_log2 = K_val / np.log(x)**2

        print(f"{x:>12} | {pi_x:>10} | {err_R:>8.2f} | {K_str:>7} | {k_over_sqrt:>10.6f} | {k_over_log2:>10.4f}")

        results.append((x, pi_x, K_exact))

    print()

    # Now use mpmath for more precise li(x^rho) computation
    print("=" * 70)
    print("HIGH-PRECISION RECOMPUTATION (using mpmath li)")
    print("=" * 70)
    print()

    # Recompute for a subset using high precision
    subset_xs = [100, 1000, 10000, 100000, 1000000]

    print(f"{'x':>12} | {'pi(x)':>10} | {'K_exact (hp)':>12} | {'K/log(x)':>10} | {'K/log(x)^2':>10}")
    print("-" * 65)

    for x in subset_xs:
        pi_x = int(primepi(x))
        R_x = riemann_R(x)

        K_exact = None
        for K in range(1, min(max_zeros + 1, 300)):
            correction = explicit_formula_correction(x, zeros, K)
            estimate = R_x + correction
            error = abs(pi_x - estimate)

            if error < 0.5 and K_exact is None:
                K_exact = K
                break

        K_str = str(K_exact) if K_exact else f">{min(max_zeros, 300)}"
        K_val = K_exact if K_exact else 300

        print(f"{x:>12} | {pi_x:>10} | {K_str:>12} | {K_val/np.log(x):>10.4f} | {K_val/np.log(x)**2:>10.4f}")

    print()
    print("ANALYSIS:")
    print("If K(x) = O(polylog(x)), then K/log(x)^2 should be bounded.")
    print("If K(x) = O(sqrt(x)), then K/sqrt(x) should be bounded.")
    print()

    # Fit K(x) scaling
    valid_results = [(x, K) for x, _, K in results if K is not None and K > 0]
    if len(valid_results) >= 3:
        log_xs = np.log([r[0] for r in valid_results])
        log_Ks = np.log([r[1] for r in valid_results])

        # Fit log(K) = alpha * log(x) + beta
        coeffs = np.polyfit(log_xs, log_Ks, 1)
        alpha = coeffs[0]
        beta = coeffs[1]

        print(f"Power-law fit: K(x) ~ x^{alpha:.4f}")
        print(f"  If alpha ≈ 0: polylog (BREAKTHROUGH)")
        print(f"  If alpha ≈ 0.5: matches sqrt(x) barrier")
        print(f"  If alpha ≈ 0.33: matches x^{1/3}")
        print()

        # Also fit K = a * log(x)^b
        from scipy.optimize import curve_fit
        try:
            def polylog_model(log_x, a, b):
                return a * log_x**b

            popt, pcov = curve_fit(polylog_model, log_xs,
                                   np.array([r[1] for r in valid_results], dtype=float),
                                   p0=[1, 2], maxfev=10000)
            print(f"Polylog fit: K(x) ~ {popt[0]:.4f} * log(x)^{popt[1]:.4f}")
            if popt[1] < 3:
                print(f"  This would be O(log(x)^{popt[1]:.1f}) -- POTENTIALLY POLYLOG!")
            else:
                print(f"  Exponent {popt[1]:.1f} is large -- likely not truly polylog")
        except Exception as e:
            print(f"  Polylog fit failed: {e}")

    print()


if __name__ == "__main__":
    main()

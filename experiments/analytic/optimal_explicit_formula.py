"""
Session 6: Optimal Explicit Formula Configuration

Given that ALL approaches reduce to the explicit formula,
what is the OPTIMAL configuration?

The explicit formula: π(x) = R(x) - Σ_ρ R(x^ρ) + small terms

Key parameters:
- Number of zeros K to use
- Smoothing kernel (none, Gaussian, Riesz, Cesàro)
- Precision of each zero
- Whether to use conjugate pairs

Previous finding: straight summation converges as K^{-0.01}
Maybe a DIFFERENT summation method converges faster?

THIS EXPERIMENT: Systematically test all combinations to find
the absolute best we can do with K zeros.
"""

import numpy as np
from mpmath import mp, mpf, li, log, sqrt, pi, zetazero, im, re, exp, cos, sin
import time
from sympy import primepi

mp.dps = 30

def explicit_formula_with_kernel(x, K, kernel='none'):
    """
    Compute π(x) using explicit formula with K zeros and specified kernel.

    Kernels:
    - 'none': raw sum
    - 'cesaro': (1 - k/K) weighting
    - 'riesz': (1 - (k/K)^2) weighting
    - 'gaussian': exp(-(k/K)^2 * sigma) weighting
    - 'fejer': (1 - k/K)^2 weighting
    """
    x = mpf(x)
    lnx = log(x)

    # R(x) via li terms with Mobius
    mobius = [0, 1, -1, -1, 0, -1, 1, -1, 0, 0, 1, -1, 0, -1, 1, 1, 0, -1, 0, -1, 0]
    R_x = mpf(0)
    for k in range(1, min(20, len(mobius))):
        if mobius[k] == 0:
            continue
        xk = x ** (mpf(1)/k)
        if xk > 1.001:
            R_x += mpf(mobius[k]) / k * li(xk)

    # Zero contributions
    correction = mpf(0)
    for k in range(1, K + 1):
        gamma = float(im(zetazero(k)))
        rho_abs = np.sqrt(0.25 + gamma**2)

        # Amplitude
        amp = 2 * float(x)**0.5 / (rho_abs * float(lnx))
        phase = gamma * float(lnx)

        # Kernel weight
        if kernel == 'none':
            w = 1.0
        elif kernel == 'cesaro':
            w = 1.0 - k / (K + 1)
        elif kernel == 'riesz':
            w = (1.0 - (k / (K + 1))**2)
        elif kernel == 'gaussian':
            sigma = 2.0
            w = np.exp(-(k / (K + 1))**2 * sigma)
        elif kernel == 'fejer':
            w = (1.0 - k / (K + 1))**2
        else:
            w = 1.0

        correction += w * amp * np.cos(phase)

    # Small correction terms
    small_corr = -1/float(lnx)

    return float(R_x) - correction + small_corr

def systematic_test():
    """Systematically test all kernel/K combinations."""
    print("="*70)
    print("SYSTEMATIC EXPLICIT FORMULA OPTIMIZATION")
    print("="*70)

    test_x_values = [100, 1000, 10000, 100000]
    true_pi = {x: int(primepi(x)) for x in test_x_values}

    kernels = ['none', 'cesaro', 'riesz', 'gaussian', 'fejer']
    K_values = [5, 10, 20, 50]

    # Header
    print(f"\n{'K':>3} {'Kernel':>10}", end="")
    for x in test_x_values:
        print(f" | x={x:>6d}(err)", end="")
    print(" | Mean |err|")
    print("-" * 90)

    best_mean_err = float('inf')
    best_config = None

    for kernel in kernels:
        for K in K_values:
            errors = []
            print(f"{K:3d} {kernel:>10}", end="")
            for x in test_x_values:
                est = explicit_formula_with_kernel(x, K, kernel)
                err = est - true_pi[x]
                errors.append(abs(err))
                print(f" | {float(err):>12.2f}", end="")

            mean_err = np.mean(errors)
            print(f" | {float(mean_err):>7.2f}")

            if mean_err < best_mean_err:
                best_mean_err = mean_err
                best_config = (K, kernel)

    print(f"\nBest configuration: K={best_config[0]}, kernel={best_config[1]}")
    print(f"Mean |error|: {best_mean_err:.2f}")

    # Now test the best config at more x values
    print(f"\n--- Detailed test of best config: K={best_config[0]}, kernel={best_config[1]} ---")
    K, kernel = best_config
    for x in [50, 200, 500, 2000, 5000, 20000, 50000, 100000]:
        true = int(primepi(x))
        est = explicit_formula_with_kernel(x, K, kernel)
        err = est - true
        within_1 = abs(err) < 1.0
        print(f"  x={x:>6d}: π(x)={true:>5d}, est={est:.2f}, err={err:+.2f}, "
              f"{'EXACT' if within_1 else ''}")

    # KEY: What K would we need for error < 1?
    print(f"\n--- How many zeros needed for error < 1? ---")
    for x in [100, 1000, 10000]:
        true = int(primepi(x))
        for K in [5, 10, 20, 50, 100, 200]:
            est = explicit_formula_with_kernel(x, K, kernel)
            err = abs(est - true)
            if err < 1.0:
                print(f"  x={x}: Need K={K} zeros for |error| < 1 (actual={err:.4f})")
                break
        else:
            est = explicit_formula_with_kernel(x, 200, kernel)
            err = abs(est - true)
            print(f"  x={x}: Still error={err:.2f} with 200 zeros")

def main():
    print("Session 6: Optimal Explicit Formula")
    print(f"Started: {time.strftime('%H:%M:%S')}\n")
    t0 = time.time()
    systematic_test()
    elapsed = time.time() - t0
    print(f"\n{'='*70}")
    print(f"Total time: {elapsed:.1f}s")

if __name__ == '__main__':
    main()

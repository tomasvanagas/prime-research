#!/usr/bin/env python3
"""
Session 13 Experiment: Optimized linear combinations of li(x^{1/k})

Question: Can we find coefficients c_1, ..., c_K such that
  pi(x) = round(sum_{k=1}^{K} c_k * li(x^{1/k}))
for ALL x in a range?

The Riemann R(x) uses c_k = mu(k)/k. What if DIFFERENT coefficients
converge faster or even give exactness for a wider range?

This probes whether the information barrier is in the BASIS FUNCTIONS
or in the COEFFICIENTS.
"""

import numpy as np
from scipy.special import expi
from sympy import primepi, isprime, mobius, li as sympy_li
import warnings
warnings.filterwarnings('ignore')

def li(x):
    """Logarithmic integral li(x) = integral from 0 to x of dt/ln(t)"""
    if x <= 1:
        return 0.0
    return float(expi(np.log(x)))

def riemann_R(x, K=50):
    """Standard Riemann R(x) = sum_{k=1}^{K} mu(k)/k * li(x^{1/k})"""
    result = 0.0
    for k in range(1, K+1):
        mu_k = int(mobius(k))
        if mu_k != 0:
            result += mu_k / k * li(x ** (1/k))
    return result

def build_basis(x_values, K):
    """Build matrix of li(x^{1/k}) for k=1..K at each x value"""
    n = len(x_values)
    A = np.zeros((n, K))
    for i, x in enumerate(x_values):
        for k in range(K):
            A[i, k] = li(x ** (1.0 / (k+1)))
    return A

def experiment_optimal_coefficients():
    """Find optimal coefficients for li(x^{1/k}) combination"""
    print("="*70)
    print("EXPERIMENT: Optimized li(x^{1/k}) linear combination")
    print("="*70)

    # Training data: pi(x) for x = 100 to 50000
    x_train = list(range(100, 10001, 10))
    pi_train = [int(primepi(x)) for x in x_train]

    # Test data: x = 10001 to 50000
    x_test = list(range(10001, 50001, 50))
    pi_test = [int(primepi(x)) for x in x_test]

    print(f"\nTraining: {len(x_train)} points, x in [100, 10000]")
    print(f"Test: {len(x_test)} points, x in [10001, 50000]")

    for K in [2, 5, 10, 20, 30, 50]:
        print(f"\n--- K = {K} basis functions li(x^{{1/k}}), k=1..{K} ---")

        # Build basis matrices
        A_train = build_basis(x_train, K)
        b_train = np.array(pi_train, dtype=float)

        A_test = build_basis(x_test, K)
        b_test = np.array(pi_test, dtype=float)

        # Solve least squares: min ||A*c - b||^2
        c_opt, residuals, rank, sv = np.linalg.lstsq(A_train, b_train, rcond=None)

        # Compare with Riemann coefficients
        c_riemann = np.zeros(K)
        for k in range(K):
            mu_k = int(mobius(k+1))
            c_riemann[k] = mu_k / (k+1)

        # Predictions
        pred_train_opt = A_train @ c_opt
        pred_test_opt = A_test @ c_opt

        pred_train_rie = A_train @ c_riemann
        pred_test_rie = A_test @ c_riemann

        # Exact matches (after rounding)
        exact_train_opt = sum(1 for i in range(len(x_train))
                             if round(pred_train_opt[i]) == pi_train[i])
        exact_test_opt = sum(1 for i in range(len(x_test))
                            if round(pred_test_opt[i]) == pi_test[i])

        exact_train_rie = sum(1 for i in range(len(x_train))
                             if round(pred_train_rie[i]) == pi_train[i])
        exact_test_rie = sum(1 for i in range(len(x_test))
                            if round(pred_test_rie[i]) == pi_test[i])

        # Max errors
        max_err_train_opt = max(abs(pred_train_opt[i] - pi_train[i]) for i in range(len(x_train)))
        max_err_test_opt = max(abs(pred_test_opt[i] - pi_test[i]) for i in range(len(x_test)))

        max_err_train_rie = max(abs(pred_train_rie[i] - pi_train[i]) for i in range(len(x_train)))
        max_err_test_rie = max(abs(pred_test_rie[i] - pi_test[i]) for i in range(len(x_test)))

        print(f"  Optimal coeffs - Train exact: {exact_train_opt}/{len(x_train)} "
              f"({100*exact_train_opt/len(x_train):.1f}%), max err: {max_err_train_opt:.3f}")
        print(f"  Optimal coeffs - Test exact:  {exact_test_opt}/{len(x_test)} "
              f"({100*exact_test_opt/len(x_test):.1f}%), max err: {max_err_test_opt:.3f}")
        print(f"  Riemann coeffs - Train exact: {exact_train_rie}/{len(x_train)} "
              f"({100*exact_train_rie/len(x_train):.1f}%), max err: {max_err_train_rie:.3f}")
        print(f"  Riemann coeffs - Test exact:  {exact_test_rie}/{len(x_test)} "
              f"({100*exact_test_rie/len(x_test):.1f}%), max err: {max_err_test_rie:.3f}")

        # Show first few optimal coefficients
        print(f"  First 5 optimal coeffs: {c_opt[:5]}")
        print(f"  First 5 Riemann coeffs: {c_riemann[:5]}")

        # Key question: does the max error grow or stay bounded?
        if K >= 10:
            # Check error growth with x for optimal coefficients
            errors_by_x = []
            checkpoints = [1000, 5000, 10000, 20000, 50000]
            for xc in checkpoints:
                if xc <= 10000:
                    idx = (xc - 100) // 10
                    if idx < len(x_train):
                        err = abs(pred_train_opt[idx] - pi_train[idx])
                        errors_by_x.append((xc, err))
                else:
                    idx = (xc - 10001) // 50
                    if idx < len(x_test):
                        err = abs(pred_test_opt[idx] - pi_test[idx])
                        errors_by_x.append((xc, err))

            print(f"  Error growth: ", end="")
            for xc, err in errors_by_x:
                print(f"x={xc}: {err:.3f}  ", end="")
            print()

def experiment_error_structure():
    """Analyze the STRUCTURE of the residual error"""
    print("\n" + "="*70)
    print("EXPERIMENT: Structure of residual error (optimal - exact)")
    print("="*70)

    K = 30
    x_values = list(range(100, 20001, 1))
    pi_values = [int(primepi(x)) for x in x_values]

    A = build_basis(x_values, K)
    b = np.array(pi_values, dtype=float)

    c_opt, _, _, _ = np.linalg.lstsq(A, b, rcond=None)
    pred = A @ c_opt
    errors = pred - b

    print(f"\nUsing K={K} basis functions, {len(x_values)} points")
    print(f"Mean error: {np.mean(errors):.6f}")
    print(f"Std error: {np.std(errors):.6f}")
    print(f"Max |error|: {np.max(np.abs(errors)):.6f}")

    # Is the error correlated with primality?
    prime_mask = np.array([isprime(x) for x in x_values])
    err_at_primes = errors[prime_mask]
    err_at_composites = errors[~prime_mask]

    print(f"\nError at primes: mean={np.mean(err_at_primes):.6f}, std={np.std(err_at_primes):.6f}")
    print(f"Error at composites: mean={np.mean(err_at_composites):.6f}, std={np.std(err_at_composites):.6f}")

    # Fourier analysis of error
    fft_err = np.fft.fft(errors)
    power = np.abs(fft_err)**2
    freqs = np.fft.fftfreq(len(errors))

    # Top 10 frequencies
    top_idx = np.argsort(power[1:len(power)//2])[-10:] + 1  # Skip DC
    print(f"\nTop 10 Fourier modes of residual error:")
    for idx in reversed(top_idx):
        print(f"  freq={freqs[idx]:.6f}, period={1/freqs[idx]:.1f}, power={power[idx]:.2f}")

    # Spectral flatness (1.0 = white noise)
    log_power = np.log(power[1:len(power)//2] + 1e-30)
    spectral_flatness = np.exp(np.mean(log_power)) / np.mean(power[1:len(power)//2])
    print(f"\nSpectral flatness of residual: {spectral_flatness:.4f} (1.0 = white noise)")

    # Autocorrelation
    from numpy.fft import ifft
    acf = np.real(ifft(power))
    acf = acf / acf[0]  # Normalize
    print(f"\nAutocorrelation at lags 1-10: {acf[1:11]}")

    # Key metric: what fraction of error's information is in top-k modes?
    sorted_power = np.sort(power[1:len(power)//2])[::-1]
    total_power = np.sum(sorted_power)
    for k in [1, 5, 10, 50, 100]:
        frac = np.sum(sorted_power[:k]) / total_power
        print(f"  Top {k} modes capture {100*frac:.1f}% of error variance")

def experiment_error_scaling():
    """How does the optimal approximation error scale with x?"""
    print("\n" + "="*70)
    print("EXPERIMENT: Error scaling with x")
    print("="*70)

    K = 20

    # For different ranges of x, train optimal coefficients and measure error
    ranges = [(100, 1000), (1000, 10000), (10000, 50000)]

    for x_lo, x_hi in ranges:
        step = max(1, (x_hi - x_lo) // 1000)
        x_values = list(range(x_lo, x_hi + 1, step))
        pi_values = [int(primepi(x)) for x in x_values]

        A = build_basis(x_values, K)
        b = np.array(pi_values, dtype=float)

        c_opt, _, _, _ = np.linalg.lstsq(A, b, rcond=None)
        pred = A @ c_opt
        errors = np.abs(pred - b)

        exact = sum(1 for i in range(len(x_values)) if round(pred[i]) == pi_values[i])

        print(f"\nRange [{x_lo}, {x_hi}], step={step}, K={K}")
        print(f"  Max |error|: {np.max(errors):.4f}")
        print(f"  Mean |error|: {np.mean(errors):.4f}")
        print(f"  Exact: {exact}/{len(x_values)} ({100*exact/len(x_values):.1f}%)")
        print(f"  Error at x_hi: {errors[-1]:.4f}")

        # Error at specific percentiles of x
        n = len(x_values)
        for pct in [25, 50, 75, 100]:
            idx = min(int(n * pct / 100), n-1)
            print(f"  Error at {pct}th percentile (x={x_values[idx]}): {errors[idx]:.4f}")

if __name__ == "__main__":
    experiment_optimal_coefficients()
    experiment_error_structure()
    experiment_error_scaling()

    print("\n" + "="*70)
    print("CONCLUSION")
    print("="*70)
    print("""
Key questions answered:
1. Do optimized coefficients beat Riemann R(x)?
2. Does the residual error have exploitable structure?
3. How does the max error scale with x?

If max error stays < 0.5 for ALL x with K = O(polylog(x)) terms,
that would imply pi(x) = round(linear combination), giving O(polylog) time.
If max error grows as x^alpha for some alpha > 0, the barrier is confirmed.
""")

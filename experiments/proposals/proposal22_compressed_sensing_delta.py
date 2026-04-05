#!/usr/bin/env python3
"""
Proposal 22: Compressed Sensing Recovery of delta(n) = p(n) - R^{-1}(n)

IDEA: delta(n) has O(log n) bits of information. If delta(n) is sparse in some
transform domain (Fourier, wavelet, or a learned dictionary), we can recover it
from O(k * polylog(n)) random measurements where k is the sparsity.

Key observations:
1. delta(n) = p(n) - R^{-1}(n) oscillates with amplitude ~sqrt(p(n))/log(p(n))
2. The Fourier transform of delta should have peaks at imaginary parts of zeta zeros
3. If delta is k-sparse in the "zeta zero basis", we need O(k log(n/k)) measurements
4. We can cheaply evaluate delta(n) for small n via sieving -- these are our measurements

The approach:
- Build a measurement matrix from delta(n) values for n = 1, ..., M (cheaply computed)
- The "dictionary" is the set of oscillatory functions {cos(gamma_j * log(R^{-1}(n)))}
- Use basis pursuit / L1 minimization to find sparse coefficients
- Extrapolate to large n using the recovered sparse representation

CONJECTURE: delta(n) is O(log(n))-sparse in the zeta zero basis, so
O(log(n)^2) measurements suffice for exact recovery.

TIME COMPLEXITY: O(polylog(n)) for the extrapolation step, but O(n^{2/3}) for
building the measurement set (unless we can subsample cleverly).
"""

import math
import os
import numpy as np
from scipy import optimize

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(SCRIPT_DIR, '..', '..', 'data')

def load_zeros(n=200):
    path = os.path.join(DATA_DIR, f'zeta_zeros_{n}.txt')
    zeros = []
    with open(path) as f:
        for line in f:
            val = line.strip()
            if val:
                zeros.append(float(val))
    return zeros

ZEROS = load_zeros(200)

def li(x):
    if x <= 0 or x <= 1:
        return 0.0
    if abs(x - 1.0) < 0.01:
        return -100.0
    lnx = math.log(x)
    total = 0.0
    term = 1.0
    for k in range(1, 100):
        term *= lnx / k
        total += term / k
        if abs(term / k) < 1e-15:
            break
    return 0.5772156649015329 + math.log(abs(lnx)) + total

def R_function(x):
    if x <= 1:
        return 0.0
    mu = [0, 1, -1, -1, 0, -1, 1, -1, 0, 0, 1,
          -1, 0, -1, 1, 1, 0, -1, 0, -1, 0]
    total = 0.0
    for k in range(1, 21):
        if mu[k] == 0:
            continue
        xk = x ** (1.0 / k)
        if xk <= 1.0001:
            break
        total += mu[k] / k * li(xk)
    return total

def R_inverse(n):
    x = n * math.log(n) + n * math.log(math.log(max(n, 3))) - n
    if x < 2:
        x = 2.0
    for _ in range(50):
        rx = R_function(x)
        if abs(rx - n) < 1e-10:
            break
        # Approximate derivative: 1/log(x)
        deriv = 1.0 / math.log(x)
        x = x - (rx - n) / deriv
        if x < 2:
            x = 2.0
    return x

def sieve_primes(limit):
    if limit < 2:
        return []
    sieve = [True] * (limit + 1)
    sieve[0] = sieve[1] = False
    for i in range(2, int(limit**0.5) + 1):
        if sieve[i]:
            for j in range(i*i, limit + 1, i):
                sieve[j] = False
    return [i for i in range(limit + 1) if sieve[i]]

def build_zeta_basis(ns, num_basis, zeros):
    """
    Build measurement matrix A where A[i,j] = cos(gamma_j * log(R^{-1}(n_i)))
    and a second matrix for sin terms.
    """
    M = len(ns)
    K = min(num_basis, len(zeros))
    # Combine cos and sin into a single matrix (2K columns)
    A = np.zeros((M, 2 * K))
    for i, n in enumerate(ns):
        x = R_inverse(n)
        lnx = math.log(x)
        sqrtx = math.sqrt(x)
        for j in range(K):
            gamma = zeros[j]
            phase = gamma * lnx
            # The actual basis function from the explicit formula:
            # contribution ~ sqrt(x) * cos(gamma*ln(x)) / (1/4 + gamma^2)
            weight = sqrtx / (0.25 + gamma * gamma) / lnx
            A[i, 2*j] = weight * math.cos(phase)
            A[i, 2*j + 1] = weight * math.sin(phase)
    return A

def compressed_sensing_recover(delta_values, A, alpha=0.01):
    """
    Recover sparse coefficients c such that A @ c ≈ delta_values
    using L1-regularized least squares (LASSO).
    """
    M, K = A.shape
    # Normalize columns
    norms = np.linalg.norm(A, axis=0)
    norms[norms == 0] = 1
    A_norm = A / norms

    # L1 minimization via scipy
    def objective(c):
        residual = A_norm @ c - delta_values
        return 0.5 * np.dot(residual, residual) + alpha * np.sum(np.abs(c))

    c0 = np.zeros(K)
    result = optimize.minimize(objective, c0, method='L-BFGS-B',
                              options={'maxiter': 1000, 'ftol': 1e-12})
    # Un-normalize
    c = result.x / norms
    return c

def predict_delta(n, coeffs, zeros):
    """Predict delta(n) using recovered sparse coefficients."""
    x = R_inverse(n)
    lnx = math.log(x)
    sqrtx = math.sqrt(x)
    K = len(coeffs) // 2
    prediction = 0.0
    for j in range(K):
        gamma = zeros[j]
        phase = gamma * lnx
        weight = sqrtx / (0.25 + gamma * gamma) / lnx
        prediction += coeffs[2*j] * weight * math.cos(phase)
        prediction += coeffs[2*j + 1] * weight * math.sin(phase)
    return prediction

def test_compressed_sensing():
    print("=" * 80)
    print("PROPOSAL 22: Compressed Sensing Recovery of delta(n)")
    print("=" * 80)

    primes = sieve_primes(120000)

    # Compute delta(n) = p(n) - R^{-1}(n) for training set
    train_ns = list(range(10, 5001))
    delta_train = np.array([primes[n-1] - R_inverse(n) for n in train_ns])

    print(f"\nTraining on n = 10..5000 ({len(train_ns)} points)")
    print(f"delta(n) stats: mean={np.mean(delta_train):.3f}, "
          f"std={np.std(delta_train):.3f}, "
          f"max_abs={np.max(np.abs(delta_train)):.3f}")

    # Test 1: Sparsity analysis - how sparse is delta in zeta zero basis?
    print("\n--- Test 1: Sparsity in zeta zero basis ---")
    for num_basis in [10, 20, 50, 100]:
        A = build_zeta_basis(train_ns, num_basis, ZEROS)
        # Least squares fit (not sparse)
        coeffs_ls, residuals, _, _ = np.linalg.lstsq(A, delta_train, rcond=None)
        pred_ls = A @ coeffs_ls
        err_ls = np.std(delta_train - pred_ls)

        # Count significant coefficients (> 10% of max)
        max_coeff = np.max(np.abs(coeffs_ls))
        n_significant = np.sum(np.abs(coeffs_ls) > 0.1 * max_coeff)

        print(f"  {num_basis:>3} basis funcs: residual_std={err_ls:.3f}, "
              f"significant_coeffs={n_significant}, "
              f"max_coeff={max_coeff:.3f}")

    # Test 2: Compressed sensing with subsampled measurements
    print("\n--- Test 2: CS recovery from subsampled measurements ---")
    num_basis = 50
    A_full = build_zeta_basis(train_ns, num_basis, ZEROS)

    for subsample_frac in [0.1, 0.2, 0.5, 1.0]:
        M_sub = max(20, int(len(train_ns) * subsample_frac))
        # Random subsample
        np.random.seed(42)
        idx = np.sort(np.random.choice(len(train_ns), M_sub, replace=False))
        A_sub = A_full[idx]
        delta_sub = delta_train[idx]

        coeffs = compressed_sensing_recover(delta_sub, A_sub, alpha=0.1)

        # Predict on held-out test set
        test_ns = list(range(5001, 8001))
        delta_test = np.array([primes[n-1] - R_inverse(n) for n in test_ns])
        A_test = build_zeta_basis(test_ns, num_basis, ZEROS)
        pred_test = A_test @ coeffs

        err_test = np.std(delta_test - pred_test)
        err_train = np.std(delta_train - A_full @ coeffs)

        # Sparsity of recovered vector
        n_nonzero = np.sum(np.abs(coeffs) > 0.01 * np.max(np.abs(coeffs)))

        print(f"  {subsample_frac*100:>5.0f}% samples ({M_sub:>5}): "
              f"train_err={err_train:.3f}, test_err={err_test:.3f}, "
              f"nonzero_coeffs={n_nonzero}")

    # Test 3: Can we use CS to predict delta for large n?
    print("\n--- Test 3: Extrapolation to larger n ---")
    # Train on n=10..5000, predict n=5001..10000
    A_train = build_zeta_basis(train_ns, 50, ZEROS)
    coeffs = compressed_sensing_recover(delta_train, A_train, alpha=0.05)

    test_ranges = [(5001, 6000), (6001, 8000), (8001, 10000)]
    for lo, hi in test_ranges:
        test_ns = list(range(lo, hi + 1))
        delta_test = np.array([primes[n-1] - R_inverse(n) for n in test_ns])
        A_test = build_zeta_basis(test_ns, 50, ZEROS)
        pred_test = A_test @ coeffs
        err = np.std(delta_test - pred_test)
        # How often would rounding give the right prime?
        correct = 0
        for i, n in enumerate(test_ns):
            x_pred = R_inverse(n) + pred_test[i]
            # Find closest prime to x_pred
            idx_prime = min(range(len(primes)),
                          key=lambda j: abs(primes[j] - x_pred))
            if primes[idx_prime] == primes[n-1]:
                correct += 1
        acc = correct / len(test_ns)
        print(f"  n={lo}-{hi}: err_std={err:.3f}, "
              f"correct_prime_rate={acc:.1%}")

    # Test 4: Fourier analysis of delta to check sparsity directly
    print("\n--- Test 4: Fourier analysis of delta(n) ---")
    # Take FFT of delta(n)
    delta_full = np.array([primes[n-1] - R_inverse(n) for n in range(10, 5001)])
    fft = np.fft.rfft(delta_full)
    magnitudes = np.abs(fft)
    sorted_mags = np.sort(magnitudes)[::-1]
    total_energy = np.sum(magnitudes**2)
    for k in [5, 10, 20, 50, 100]:
        top_k_energy = np.sum(sorted_mags[:k]**2)
        print(f"  Top {k:>3} Fourier coeffs capture "
              f"{top_k_energy/total_energy*100:.1f}% of energy")

    print("\n--- VERDICT ---")
    print("Compressed sensing can reduce the CONSTANT factor in the number of")
    print("zeros needed, but the fundamental barrier is that delta(n) encodes")
    print("O(sqrt(n)/log(n)) bits of information from zeta zeros.")
    print("True polylog recovery requires delta to be O(polylog)-sparse,")
    print("which would imply a number-theoretic breakthrough.")

if __name__ == '__main__':
    test_compressed_sensing()

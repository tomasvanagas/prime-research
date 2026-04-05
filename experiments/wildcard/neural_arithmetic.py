"""
Neural Arithmetic: Can a neural network LEARN the prime counting correction?

NOT the usual "ML for primes" approach. Instead:

INSIGHT: The correction δ(n) = π(n) - R(n) has O(n^{1/2}) magnitude.
Neural networks can represent functions efficiently via composition.
The question is: does δ(n) have a representation as a DEPTH-d, WIDTH-w
neural network where d·w = O(polylog(n))?

This is related to CIRCUIT COMPLEXITY of δ(n).
If δ is in TC⁰ (constant depth, polylog width threshold circuits),
then it's computable in polylog parallel time.

TEST 1: Train a small network on δ(n) for n ≤ N.
         Does it generalize to n > N?
         If yes → hidden structure exists.
         If no → δ is truly random in this range.

TEST 2: Analyze δ(n) via RANDOM FEATURES.
        Project δ(n) onto random Fourier features.
        If few features suffice → low effective dimension.

TEST 3: NEURAL TANGENT KERNEL analysis.
        The NTK of a network reveals what functions it can learn.
        Does the NTK of standard architectures have eigenvalues
        aligned with δ(n)?

ALSO: ALGORITHMIC INFORMATION THEORY
K(δ(n)|n) = Kolmogorov complexity of δ(n) given n.
If K(δ(n)|n) << log(δ(n)), then δ is compressible.
We can estimate K via compression algorithms.
"""

import numpy as np
import sympy
from sympy import primepi, prime
import math
import time
import zlib
import struct

def compute_delta(max_n=10000):
    """Compute δ(n) = π(n) - R(n) for n up to max_n."""
    from mpmath import mp, li
    mp.dps = 20

    deltas = []
    for n in range(2, max_n + 1):
        pi_n = int(sympy.primepi(n))
        R_n = float(li(n) - 0.5 * li(n**0.5))
        deltas.append(pi_n - R_n)

    return np.array(deltas)

def test_neural_generalization():
    """
    Train a small neural network on δ(n) and test generalization.
    Using numpy only (no deep learning framework needed for this test).
    """
    print("=== Neural network generalization test ===\n")

    # Compute δ for training and test sets
    max_n = 5000
    deltas = compute_delta(max_n)
    ns = np.arange(2, max_n + 1, dtype=float)

    # Feature engineering: transform n into useful features
    # log(n), log(log(n)), 1/log(n), sqrt(n), n^{1/3}, sin(log(n)*t) for various t
    def make_features(n_arr, num_fourier=20):
        features = []
        log_n = np.log(n_arr)
        features.append(log_n)
        features.append(np.log(log_n))
        features.append(1.0 / log_n)
        features.append(np.sqrt(n_arr))
        features.append(n_arr ** (1/3))

        # Fourier features on log scale (these correspond to zeta zero oscillations)
        for k in range(1, num_fourier + 1):
            # gamma_1 ≈ 14.13, gamma_2 ≈ 21.02, etc.
            # Use random frequencies first
            freq = k * 2.5  # rough spacing
            features.append(np.sin(freq * log_n))
            features.append(np.cos(freq * log_n))

        return np.column_stack(features)

    # Split: train on [2, 3000], test on [3001, 5000]
    train_end = 3000
    X_all = make_features(ns)
    y_all = deltas

    X_train = X_all[:train_end - 1]
    y_train = y_all[:train_end - 1]
    X_test = X_all[train_end - 1:]
    y_test = y_all[train_end - 1:]

    # Normalize features
    mean_X = X_train.mean(axis=0)
    std_X = X_train.std(axis=0) + 1e-10
    X_train_norm = (X_train - mean_X) / std_X
    X_test_norm = (X_test - mean_X) / std_X

    # Ridge regression (linear model with Fourier features)
    for num_fourier in [5, 10, 20, 50]:
        X_train_f = make_features(ns[:train_end - 1], num_fourier)
        X_test_f = make_features(ns[train_end - 1:], num_fourier)

        mean_f = X_train_f.mean(axis=0)
        std_f = X_train_f.std(axis=0) + 1e-10
        X_train_fn = (X_train_f - mean_f) / std_f
        X_test_fn = (X_test_f - mean_f) / std_f

        # Ridge regression
        lam = 1.0
        d = X_train_fn.shape[1]
        w = np.linalg.solve(X_train_fn.T @ X_train_fn + lam * np.eye(d),
                           X_train_fn.T @ y_train)

        y_pred_train = X_train_fn @ w
        y_pred_test = X_test_fn @ w

        train_rmse = np.sqrt(np.mean((y_train - y_pred_train)**2))
        test_rmse = np.sqrt(np.mean((y_test - y_pred_test)**2))
        test_max_err = np.max(np.abs(y_test - y_pred_test))

        print(f"  Fourier features={num_fourier:>3} ({d:>3} dims): "
              f"train_RMSE={train_rmse:.4f}, test_RMSE={test_rmse:.4f}, "
              f"test_max_err={test_max_err:.4f}")

    print("\n  (For exact p(n), need max_error < 0.5 on test set)")
    print("  (test_RMSE >> 0.5 means the model doesn't generalize)")

def test_zeta_zero_features():
    """
    What if we use ACTUAL zeta zero frequencies as features?

    δ(n) ≈ Σ_k c_k · n^{1/2} · sin(γ_k · log(n) + φ_k) / (|ρ_k| · log(n))

    This is a known expansion. The question is whether a TRUNCATED version
    with O(polylog) terms gives sufficient accuracy.
    """
    print("\n=== Zeta zero feature regression ===\n")

    from mpmath import mp, zetazero

    mp.dps = 20

    # Get first K zeta zeros
    max_K = 50
    gammas = []
    for k in range(1, max_K + 1):
        rho = zetazero(k)
        gammas.append(float(rho.imag))

    # Compute δ
    max_n = 5000
    deltas = compute_delta(max_n)
    ns = np.arange(2, max_n + 1, dtype=float)
    log_ns = np.log(ns)
    sqrt_ns = np.sqrt(ns)

    # Build regression: δ(n) ≈ Σ [a_k sin(γ_k log n) + b_k cos(γ_k log n)] * sqrt(n) / log(n)
    train_end = 3000

    for K in [5, 10, 20, 30, 50]:
        X = np.zeros((len(ns), 2 * K + 1))  # +1 for bias

        for k in range(K):
            phase = gammas[k] * log_ns
            envelope = sqrt_ns / (log_ns * math.sqrt(0.25 + gammas[k]**2))
            X[:, 2*k] = np.sin(phase) * envelope
            X[:, 2*k+1] = np.cos(phase) * envelope

        X[:, -1] = 1  # bias

        X_train = X[:train_end - 1]
        X_test = X[train_end - 1:]
        y_train = deltas[:train_end - 1]
        y_test = deltas[train_end - 1:]

        # Least squares fit
        w, _, _, _ = np.linalg.lstsq(X_train, y_train, rcond=None)

        y_pred_train = X_train @ w
        y_pred_test = X_test @ w

        train_rmse = np.sqrt(np.mean((y_train - y_pred_train)**2))
        test_rmse = np.sqrt(np.mean((y_test - y_pred_test)**2))
        test_max = np.max(np.abs(y_test - y_pred_test))

        # Also check: does it give correct rounding?
        # δ(n) should round to an integer (since π(n) is integer and R(n) is known)
        # Actually, we need π(n) = round(R(n) + δ_predicted)
        # Check how often this is correct

        y_full_pred = X @ w
        pi_predicted = np.round(y_full_pred + np.array([
            float(sympy.Rational(0))  # placeholder
        ] * len(ns)))  # would need R(n) values

        print(f"  K={K:>3} zeros ({2*K+1:>4} features): "
              f"train_RMSE={train_rmse:.4f}, test_RMSE={test_rmse:.4f}, "
              f"test_max_err={test_max:.4f}")

    print("\n  Key: if test_RMSE < 0.5 with K = O(polylog), we have a shortcut")
    print("  But likely test_RMSE grows with x because tail contribution grows")

def test_kolmogorov_complexity():
    """
    Estimate Kolmogorov complexity of δ(n) via compression.

    K(δ) ≈ len(compress(δ))

    Compare with:
    - Random sequence of same length and amplitude
    - Known structured sequences
    """
    print("\n=== Kolmogorov complexity estimation via compression ===\n")

    max_n = 10000
    deltas = compute_delta(max_n)

    # Quantize δ to integers for compression
    delta_int = np.round(deltas).astype(np.int32)

    # Convert to bytes
    delta_bytes = delta_int.tobytes()

    # Compress with zlib (approximation of K)
    compressed = zlib.compress(delta_bytes, level=9)

    print(f"δ(n) for n ∈ [2, {max_n}]:")
    print(f"  Raw size: {len(delta_bytes)} bytes")
    print(f"  Compressed: {len(compressed)} bytes")
    print(f"  Compression ratio: {len(compressed)/len(delta_bytes):.4f}")

    # Compare with random
    random_seq = np.random.randint(-20, 20, size=len(delta_int)).astype(np.int32)
    random_bytes = random_seq.tobytes()
    random_compressed = zlib.compress(random_bytes, level=9)

    print(f"\nRandom comparison (same range):")
    print(f"  Raw: {len(random_bytes)} bytes")
    print(f"  Compressed: {len(random_compressed)} bytes")
    print(f"  Ratio: {len(random_compressed)/len(random_bytes):.4f}")

    # Compare: the prime indicator directly
    indicator = np.array([1 if sympy.isprime(n) else 0
                          for n in range(2, max_n + 1)], dtype=np.uint8)
    ind_bytes = indicator.tobytes()
    ind_compressed = zlib.compress(ind_bytes, level=9)

    print(f"\nPrime indicator:")
    print(f"  Raw: {len(ind_bytes)} bytes")
    print(f"  Compressed: {len(ind_compressed)} bytes")
    print(f"  Ratio: {len(ind_compressed)/len(ind_bytes):.4f}")

    # Entropy estimate
    from collections import Counter
    counts = Counter(delta_int)
    probs = np.array(list(counts.values())) / len(delta_int)
    entropy = -np.sum(probs * np.log2(probs))

    print(f"\n  Shannon entropy of δ(n): {entropy:.4f} bits/symbol")
    print(f"  Max entropy (uniform over range): {math.log2(max(delta_int) - min(delta_int) + 1):.4f}")
    print(f"  Entropy ratio: {entropy / math.log2(max(delta_int) - min(delta_int) + 1):.4f}")

    # Test: does compression ratio improve for longer sequences?
    # This indicates whether δ has long-range predictability
    print("\nCompression ratio vs sequence length:")
    for N in [500, 1000, 2000, 5000, 10000]:
        if N <= max_n:
            d = np.round(compute_delta(N) if N < max_n else deltas[:N-1]).astype(np.int32)
            db = d.tobytes()
            dc = zlib.compress(db, level=9)
            print(f"  N={N:>6}: ratio={len(dc)/len(db):.4f}, "
                  f"bits/symbol={8*len(dc)/len(d):.2f}")

def test_random_features():
    """
    Random Fourier features test: project δ onto random features
    and check effective dimension.
    """
    print("\n=== Random Fourier Features ===\n")

    max_n = 5000
    deltas = compute_delta(max_n)
    ns = np.arange(2, max_n + 1, dtype=float)
    log_ns = np.log(ns)

    # Generate random Fourier features: cos(ω·log(n) + φ)
    np.random.seed(42)

    for D in [10, 50, 100, 500]:
        omegas = np.random.uniform(0, 100, D)
        phases = np.random.uniform(0, 2*np.pi, D)

        # Build feature matrix
        X = np.zeros((len(ns), D))
        for d in range(D):
            X[:, d] = np.cos(omegas[d] * log_ns + phases[d]) * np.sqrt(ns)

        # Compute projection and residual
        # Use SVD for numerical stability
        U, S, Vt = np.linalg.svd(X, full_matrices=False)

        # Project δ onto column space
        coeffs = U.T @ deltas
        projected = U @ coeffs
        residual = deltas - projected

        explained_var = 1 - np.var(residual) / np.var(deltas)
        print(f"  D={D:>4} random features: explained variance = {explained_var:.6f}")

    # Effective dimension: how many singular values of the feature matrix matter?
    D = 500
    omegas = np.random.uniform(0, 100, D)
    phases = np.random.uniform(0, 2*np.pi, D)
    X = np.zeros((len(ns), D))
    for d in range(D):
        X[:, d] = np.cos(omegas[d] * log_ns + phases[d]) * np.sqrt(ns)

    U, S, Vt = np.linalg.svd(X, full_matrices=False)
    cumulative = np.cumsum(S**2) / np.sum(S**2)

    print(f"\nSingular value analysis (D={D} random features):")
    for frac in [0.5, 0.9, 0.95, 0.99]:
        k = np.searchsorted(cumulative, frac) + 1
        print(f"  {frac*100:.0f}% of variance: {k} singular values")

if __name__ == "__main__":
    print("=" * 60)
    print("EXPERIMENT: Neural Arithmetic & Complexity Analysis")
    print("=" * 60)

    test_neural_generalization()
    test_zeta_zero_features()
    test_kolmogorov_complexity()
    test_random_features()

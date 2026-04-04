"""
Experiment: Products of floor values and their relation to primes.

Key questions:
1. Does floor(x/a)*floor(x/b) vs floor(x/(a*b)) carry prime information?
2. Can products of floor values detect primality or count primes?
3. Is there a polynomial in floor values that equals pi(x)?

The "correction term" delta(a,b,x) = floor(x/a)*floor(x/b) - floor(x/(a*b))
is always >= 0 and encodes arithmetic information about x, a, b.
"""

import numpy as np
from sympy import primepi, isprime, factorint
from itertools import combinations_with_replacement
import math

def floor_div(x, k):
    return x // k

def correction_term(x, a, b):
    """delta(a,b,x) = floor(x/a)*floor(x/b) - floor(x/(a*b))"""
    return (x // a) * (x // b) - x // (a * b)

# ============================================================
# Experiment 1: Does the correction term detect primes?
# ============================================================
print("=" * 70)
print("EXPERIMENT 1: Correction terms delta(a,b,x) for prime detection")
print("=" * 70)

x_values = range(2, 60)
for a, b in [(2, 3), (2, 5), (3, 5), (2, 7)]:
    deltas_prime = []
    deltas_composite = []
    for x in x_values:
        d = correction_term(x, a, b)
        if isprime(x):
            deltas_prime.append(d)
        else:
            deltas_composite.append(d)
    print(f"\ndelta({a},{b},x):")
    print(f"  Primes mean={np.mean(deltas_prime):.2f}, std={np.std(deltas_prime):.2f}")
    print(f"  Composites mean={np.mean(deltas_composite):.2f}, std={np.std(deltas_composite):.2f}")
    # Check if distributions are separable
    overlap = len(set(deltas_prime) & set(deltas_composite))
    print(f"  Value overlap: {overlap} values appear in both sets")
    print(f"  Prime values: {sorted(set(deltas_prime))[:15]}")
    print(f"  Composite values: {sorted(set(deltas_composite))[:15]}")

# ============================================================
# Experiment 2: Can pi(x) be written as polynomial in floor values?
# ============================================================
print("\n" + "=" * 70)
print("EXPERIMENT 2: Polynomial fit of pi(x) from floor values")
print("=" * 70)

# For each x, the "features" are floor(x/k) for k=1,...,K
# Try to find coefficients c_k such that pi(x) = sum c_k * floor(x/k)
# (linear case first, then quadratic)

def get_floor_features_linear(x, K):
    return [x // k for k in range(1, K + 1)]

def get_floor_features_quadratic(x, K):
    """Linear + all pairwise products of floor(x/k)"""
    linear = [x // k for k in range(1, K + 1)]
    quad = []
    for i in range(K):
        for j in range(i, K):
            quad.append(linear[i] * linear[j])
    return linear + quad

for K in [3, 5, 8, 10]:
    X_range = range(2, 200)

    # Linear features
    A_lin = np.array([get_floor_features_linear(x, K) for x in X_range], dtype=float)
    y = np.array([primepi(x) for x in X_range], dtype=float)

    # Least squares fit
    coeffs_lin, res_lin, _, _ = np.linalg.lstsq(A_lin, y, rcond=None)
    pred_lin = A_lin @ coeffs_lin
    err_lin = np.max(np.abs(pred_lin - y))
    rmse_lin = np.sqrt(np.mean((pred_lin - y) ** 2))

    # Quadratic features
    A_quad = np.array([get_floor_features_quadratic(x, K) for x in X_range], dtype=float)
    coeffs_quad, res_quad, _, _ = np.linalg.lstsq(A_quad, y, rcond=None)
    pred_quad = A_quad @ coeffs_quad
    err_quad = np.max(np.abs(pred_quad - y))
    rmse_quad = np.sqrt(np.mean((pred_quad - y) ** 2))

    # Check exact integer matches
    exact_lin = sum(1 for i, x in enumerate(X_range) if round(pred_lin[i]) == y[i])
    exact_quad = sum(1 for i, x in enumerate(X_range) if round(pred_quad[i]) == y[i])

    n_lin = A_lin.shape[1]
    n_quad = A_quad.shape[1]
    print(f"\nK={K}: {n_lin} linear features, {n_quad} quadratic features, {len(list(X_range))} data points")
    print(f"  Linear:    max_err={err_lin:.4f}, RMSE={rmse_lin:.4f}, exact={exact_lin}/{len(list(X_range))}")
    print(f"  Quadratic: max_err={err_quad:.4f}, RMSE={rmse_quad:.4f}, exact={exact_quad}/{len(list(X_range))}")

# ============================================================
# Experiment 3: Exact polynomial search (small x, enumerate)
# ============================================================
print("\n" + "=" * 70)
print("EXPERIMENT 3: Exact polynomial in floor values for pi(x)")
print("=" * 70)

# Try: can we find INTEGER coefficients c_{i,j} such that
# pi(x) = sum_{i<=j} c_{i,j} * floor(x/i) * floor(x/j) + sum_k d_k * floor(x/k) + e
# for ALL x in a range?

from numpy.linalg import matrix_rank

for K in [4, 6, 8]:
    X_range = list(range(2, 100))

    # Build feature matrix: constant, linear, quadratic
    features = []
    feature_names = ["1"]
    for x in X_range:
        row = [1]  # constant
        fv = [x // k for k in range(1, K + 1)]
        row.extend(fv)
        for i in range(K):
            for j in range(i, K):
                row.append(fv[i] * fv[j])
        features.append(row)

    if not feature_names or len(feature_names) == 1:
        for k in range(1, K + 1):
            feature_names.append(f"f({k})")
        for i in range(K):
            for j in range(i, K):
                feature_names.append(f"f({i+1})*f({j+1})")

    A = np.array(features, dtype=float)
    y = np.array([primepi(x) for x in X_range], dtype=float)

    rank = matrix_rank(A)
    print(f"\nK={K}: {A.shape[1]} features, rank={rank}, {len(X_range)} equations")

    # Solve
    coeffs, res, _, _ = np.linalg.lstsq(A, y, rcond=None)
    pred = A @ coeffs
    max_err = np.max(np.abs(pred - y))
    exact = sum(1 for i in range(len(X_range)) if abs(round(pred[i]) - y[i]) < 0.01)

    print(f"  Max error: {max_err:.6f}")
    print(f"  Exact (after rounding): {exact}/{len(X_range)}")

    # Check generalization on x=100..200
    X_test = list(range(100, 200))
    A_test = []
    for x in X_test:
        row = [1]
        fv = [x // k for k in range(1, K + 1)]
        row.extend(fv)
        for i in range(K):
            for j in range(i, K):
                row.append(fv[i] * fv[j])
        A_test.append(row)
    A_test = np.array(A_test, dtype=float)
    y_test = np.array([primepi(x) for x in X_test], dtype=float)
    pred_test = A_test @ coeffs
    max_err_test = np.max(np.abs(pred_test - y_test))
    exact_test = sum(1 for i in range(len(X_test)) if abs(round(pred_test[i]) - y_test[i]) < 0.01)
    print(f"  Generalization (100-199): max_err={max_err_test:.4f}, exact={exact_test}/{len(X_test)}")

# ============================================================
# Experiment 4: Correction terms summed as potential pi(x) formula
# ============================================================
print("\n" + "=" * 70)
print("EXPERIMENT 4: Sum of correction terms vs pi(x)")
print("=" * 70)

# Consider S(x) = sum_{a,b <= K} w(a,b) * delta(a,b,x) for some weights
# Can we choose weights to make S(x) = pi(x)?

for K in [5, 8, 10]:
    X_range = list(range(2, 150))
    features = []
    for x in X_range:
        row = []
        for a in range(1, K + 1):
            for b in range(a, K + 1):
                row.append(correction_term(x, a, b))
        features.append(row)

    A = np.array(features, dtype=float)
    y = np.array([primepi(x) for x in X_range], dtype=float)

    rank = matrix_rank(A)
    n_feat = A.shape[1]

    coeffs, res, _, _ = np.linalg.lstsq(A, y, rcond=None)
    pred = A @ coeffs
    max_err = np.max(np.abs(pred - y))
    exact = sum(1 for i in range(len(X_range)) if abs(round(pred[i]) - y[i]) < 0.01)

    print(f"\nK={K}: {n_feat} correction terms, rank={rank}")
    print(f"  Max error: {max_err:.6f}, exact={exact}/{len(X_range)}")

print("\n" + "=" * 70)
print("EXPERIMENT 5: How many distinct floor values exist and their info content")
print("=" * 70)

for x in [100, 1000, 10000]:
    floor_vals = set()
    for k in range(1, x + 1):
        floor_vals.add(x // k)
    n_distinct = len(floor_vals)
    pi_x = primepi(x)
    bits_needed = math.log2(x) if x > 0 else 0
    print(f"x={x}: {n_distinct} distinct floor(x/k) values (~2*sqrt({x})={2*int(math.sqrt(x))}), pi(x)={pi_x}")

    # Products of pairs of distinct floor values
    fv_list = sorted(floor_vals)
    product_vals = set()
    for i in range(len(fv_list)):
        for j in range(i, min(i + 50, len(fv_list))):  # sample
            product_vals.add(fv_list[i] * fv_list[j])
    print(f"  Sample pairwise products: {len(product_vals)} distinct values (from {min(50, len(fv_list))} pairs each)")

print("\nDone.")

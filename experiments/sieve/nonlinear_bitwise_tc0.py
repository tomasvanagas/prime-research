"""
Experiment: Bitwise operations on floor values for TC^0 prime counting.

TC^0 can compute: addition, multiplication, comparison, threshold gates.
Also: iterated addition (counting), majority, modular arithmetic.

Key question: Can pi(x) be expressed as a TC^0 circuit of POLYNOMIAL size
using floor values as intermediate results?

The standard Meissel-Lehmer uses O(x^{2/3}) floor values.
If we allow NONLINEAR (quadratic, bitwise) operations on a SMALL set
of floor values, can we reduce the number needed?

Focus: express pi(x) using O(polylog(x)) floor values combined nonlinearly.
"""

import numpy as np
from sympy import primepi, isprime, primerange
import math
from itertools import combinations

# ============================================================
# Experiment 1: Information in individual bits of floor(x/k)
# ============================================================
print("=" * 70)
print("EXPERIMENT 1: Bit decomposition of floor values")
print("=" * 70)

for x in [100, 500, 1000]:
    pi_x = primepi(x)
    N = int(math.log2(x)) + 1  # number of bits in x

    # For each bit position b and each k, extract bit b of floor(x/k)
    # bit(b, floor(x/k)) = (floor(x/k) >> b) & 1

    K = 20  # number of floor values to use

    # Build feature matrix: each row is x', columns are bits of floor(x'/k)
    # We want to predict pi(x') from bits of floor(x'/1), ..., floor(x'/K)

    # But we're fixing x and varying... no, for the circuit model,
    # x is the input and we want to compute pi(x).

    # For the experiment: vary x and see if bit patterns of floor(x/k)
    # for small k determine pi(x).

    print(f"\nx={x}: N={N} bits, pi(x)={pi_x}")

# Let's do the proper experiment: for x in a range, extract bits of floor(x/k)
# for k=1..K, and see if pi(x) is determined.

K = 15
X_range = list(range(2, 500))
N_bits = 10  # use bottom 10 bits

features_bits = []
targets = []
for x in X_range:
    row = []
    for k in range(1, K + 1):
        v = x // k
        for b in range(N_bits):
            row.append((v >> b) & 1)
    features_bits.append(row)
    targets.append(primepi(x))

A = np.array(features_bits, dtype=float)
y = np.array(targets, dtype=float)

# Linear fit using bits
coeffs, _, _, _ = np.linalg.lstsq(A, y, rcond=None)
pred = A @ coeffs
max_err = np.max(np.abs(pred - y))
exact = sum(1 for i in range(len(X_range)) if abs(round(pred[i]) - y[i]) < 0.01)
print(f"\nLinear on bits of floor(x/k), K={K}, {N_bits} bits each:")
print(f"  Features: {A.shape[1]}, max_err={max_err:.4f}, exact={exact}/{len(X_range)}")

# Now try AND/OR/XOR pairs of bits as additional features
print("\nAdding pairwise AND of bit features...")
# Too many pairs for large K*N_bits, so subsample
n_base = K * N_bits
# Select 50 random pairs
np.random.seed(42)
pairs = [(np.random.randint(n_base), np.random.randint(n_base)) for _ in range(100)]
pairs = list(set((min(a, b), max(a, b)) for a, b in pairs if a != b))[:80]

features_and = []
for i, x in enumerate(X_range):
    base = features_bits[i]
    row = list(base)
    for a, b in pairs:
        row.append(base[a] & base[b])  # AND
    features_and.append(row)

A2 = np.array(features_and, dtype=float)
coeffs2, _, _, _ = np.linalg.lstsq(A2, y, rcond=None)
pred2 = A2 @ coeffs2
max_err2 = np.max(np.abs(pred2 - y))
exact2 = sum(1 for i in range(len(X_range)) if abs(round(pred2[i]) - y[i]) < 0.01)
print(f"  Features: {A2.shape[1]}, max_err={max_err2:.4f}, exact={exact2}/{len(X_range)}")

# ============================================================
# Experiment 2: Majority/threshold on floor value bits
# ============================================================
print("\n" + "=" * 70)
print("EXPERIMENT 2: Threshold/majority operations")
print("=" * 70)

# TC^0 can compute MAJORITY and THRESHOLD gates.
# Threshold_t(x_1,...,x_n) = 1 iff at least t of x_i are 1

# For each x, define: T(x,k) = [floor(x/k) > x/(k+1)]
# This is [the gap is positive], i.e., there's at least one integer in (x/(k+1), x/k]

for x_val in [100, 500]:
    pi_x = primepi(x_val)

    # Threshold indicators for various k
    indicators = []
    for k in range(2, x_val):
        # Is floor(x/k) > floor(x/(k+1))?
        ind = 1 if (x_val // k) > (x_val // (k + 1)) else 0
        indicators.append(ind)

    # Sum of indicators = number of distinct floor values - 1
    sum_ind = sum(indicators)

    print(f"\nx={x_val}: pi(x)={pi_x}")
    print(f"  Sum of gap indicators: {sum_ind}")
    print(f"  Distinct floor values: {len(set(x_val // k for k in range(1, x_val + 1)))}")

    # Now: for each k, a more refined indicator
    # I(k) = [floor(x/k) mod k == 0] -- detects when k | x
    divisor_indicators = [1 if x_val % k == 0 else 0 for k in range(2, x_val)]
    n_divisors = sum(divisor_indicators)
    print(f"  Number of divisors of {x_val} in [2,{x_val}): {n_divisors}")

    # Sieve-style: for each k, I(k) = [there exists j < k such that j|k]
    # This is 1 for all composites, 0 for primes
    # Computing this requires testing divisibility for each pair -- O(x * sqrt(x))

    # TC^0 version: for each k, compute k mod j for j=2..sqrt(k)
    # Then take OR (which is a threshold-1 gate)
    # This gives primality testing for each k separately
    # Total: O(x * x^{1/2}) = O(x^{3/2}) operations in TC^0
    # Still too many -- we want O(polylog(x))

# ============================================================
# Experiment 3: Can O(polylog) floor values suffice with nonlinear ops?
# ============================================================
print("\n" + "=" * 70)
print("EXPERIMENT 3: Minimum floor values needed")
print("=" * 70)

# Theoretical bound: floor(x/k) for k=1..K gives us K values
# Each is an N-bit number. Total information: K*N bits.
# pi(x) requires ~N bits to specify.
# So K >= 1 suffices informationally... but can we COMPUTE it?

# Key: floor(x/1) = x (trivially), but pi(x) is a complex function of x.
# The question is whether NONLINEAR operations on {floor(x/k)} for
# k in a SMALL set can compute pi(x).

# Test: for K floor values, what's the best we can do with degree-d polynomials?

for K in [2, 3, 5, 8, 12]:
    X_range = list(range(10, 300))

    for degree in [1, 2, 3]:
        # Features: all monomials of degree <= d in {floor(x/1), ..., floor(x/K)}
        from itertools import combinations_with_replacement

        def get_monomial_features(x, K, degree):
            fv = [x // k for k in range(1, K + 1)]
            features = [1]  # constant
            for d in range(1, degree + 1):
                for combo in combinations_with_replacement(range(K), d):
                    val = 1
                    for idx in combo:
                        val *= fv[idx]
                    features.append(val)
            return features

        # Count features for this degree
        sample_feat = get_monomial_features(100, K, degree)
        n_feat = len(sample_feat)

        if n_feat > len(X_range) - 10:
            continue  # Skip if underdetermined

        A = np.array([get_monomial_features(x, K, degree) for x in X_range], dtype=float)
        y = np.array([primepi(x) for x in X_range], dtype=float)

        coeffs, _, _, _ = np.linalg.lstsq(A, y, rcond=None)
        pred = A @ coeffs
        max_err = np.max(np.abs(pred - y))
        exact = sum(1 for i in range(len(X_range)) if abs(round(pred[i]) - y[i]) < 0.01)

        # Generalization
        X_test = list(range(300, 500))
        A_test = np.array([get_monomial_features(x, K, degree) for x in X_test], dtype=float)
        y_test = np.array([primepi(x) for x in X_test], dtype=float)
        pred_test = A_test @ coeffs
        max_err_test = np.max(np.abs(pred_test - y_test))
        exact_test = sum(1 for i in range(len(X_test)) if abs(round(pred_test[i]) - y_test[i]) < 0.01)

        print(f"K={K:2d}, deg={degree}: {n_feat:4d} features | "
              f"train max_err={max_err:.2f}, exact={exact}/{len(X_range)} | "
              f"test max_err={max_err_test:.2f}, exact={exact_test}/{len(X_test)}")

# ============================================================
# Experiment 4: Information-theoretic analysis
# ============================================================
print("\n" + "=" * 70)
print("EXPERIMENT 4: Mutual information between floor values and pi(x)")
print("=" * 70)

# How much information does each floor(x/k) carry about pi(x)?
# Use binning to estimate mutual information

from collections import Counter

def entropy(values):
    counts = Counter(values)
    total = len(values)
    return -sum((c/total) * math.log2(c/total) for c in counts.values())

def mutual_info(x_vals, y_vals):
    """Estimate mutual information between x and y (both discrete)"""
    h_y = entropy(y_vals)
    # H(Y|X)
    x_groups = {}
    for xi, yi in zip(x_vals, y_vals):
        if xi not in x_groups:
            x_groups[xi] = []
        x_groups[xi].append(yi)

    h_y_given_x = sum(len(g)/len(y_vals) * entropy(g) for g in x_groups.values() if len(g) > 1)
    return h_y - h_y_given_x

X_range = list(range(2, 2000))
pi_vals = [primepi(x) for x in X_range]

print(f"\nH(pi(x)) for x in [2,2000): {entropy(pi_vals):.3f} bits")

for k in [1, 2, 3, 5, 7, 10, 15, 20, 30, 50]:
    fv = [x // k for x in X_range]
    mi = mutual_info(fv, pi_vals)
    print(f"  I(floor(x/{k:2d}); pi(x)) = {mi:.3f} bits")

# Pairs of floor values
print("\nMutual info from PAIRS of floor values:")
for k1, k2 in [(1, 2), (2, 3), (2, 5), (3, 5), (1, 3), (5, 7)]:
    pair_vals = [(x // k1, x // k2) for x in X_range]
    mi = mutual_info(pair_vals, pi_vals)
    print(f"  I((floor(x/{k1}), floor(x/{k2})); pi(x)) = {mi:.3f} bits")

# Products
print("\nMutual info from PRODUCTS of floor values:")
for k1, k2 in [(2, 3), (2, 5), (3, 5)]:
    prod_vals = [(x // k1) * (x // k2) for x in X_range]
    mi = mutual_info(prod_vals, pi_vals)
    print(f"  I(floor(x/{k1})*floor(x/{k2}); pi(x)) = {mi:.3f} bits")

# XOR
print("\nMutual info from XOR of floor values:")
for k1, k2 in [(2, 3), (2, 5), (3, 5)]:
    xor_vals = [(x // k1) ^ (x // k2) for x in X_range]
    mi = mutual_info(xor_vals, pi_vals)
    print(f"  I(floor(x/{k1}) XOR floor(x/{k2}); pi(x)) = {mi:.3f} bits")

print("\nDone.")

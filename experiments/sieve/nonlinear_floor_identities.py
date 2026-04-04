"""
Experiment: Search for novel identities involving nonlinear floor combinations.

The Hermite identity: floor(nx) = sum_{k=0}^{n-1} floor(x + k/n)
can express floor(nx) as a sum of shifted floors.

Question: Is there an identity that expresses pi(x) using a NONLINEAR
combination of O(polylog(x)) floor values?

We search empirically for:
1. Identities of the form pi(x) = F(floor(x/a1), ..., floor(x/ak))
2. Recursive relations: pi(x) = G(pi(x/2), floor(x/3), ...)
3. Novel floor-based formulas that are EXACT for small x
"""

import numpy as np
from sympy import primepi, isprime, primerange
import math
from itertools import product as iproduct

# ============================================================
# Experiment 1: Exhaustive search for simple nonlinear formulas
# ============================================================
print("=" * 70)
print("EXPERIMENT 1: Searching for simple nonlinear pi(x) formulas")
print("=" * 70)

# Template: pi(x) = a*floor(x/2) + b*floor(x/3) + c*floor(x/2)*floor(x/3) + d
# where a,b,c,d are rational constants

# More general: for K divisors and degree D, enumerate formulas and check

def check_formula(formula_func, x_range):
    """Check if formula matches pi(x) for all x in range"""
    errors = []
    for x in x_range:
        try:
            val = formula_func(x)
            target = primepi(x)
            errors.append(val - target)
        except:
            return None, None
    errors = np.array(errors)
    exact = sum(1 for e in errors if abs(e) < 0.01)
    return errors, exact

# Try some hand-crafted formulas based on known identities
x_range = list(range(2, 200))

formulas = {
    "floor(x/1) - floor(x/2) - 1": lambda x: x - x//2 - 1,
    "floor(x/ln(x)) approx": lambda x: int(x / math.log(x)) if x > 1 else 0,
    "x - sum floor(x/p) + sum floor(x/pq) (Legendre-like, p,q<=5)":
        lambda x: x - 1 - (x//2 + x//3 + x//5) + (x//6 + x//10 + x//15) - x//30,
    "floor(x/2)*floor(x/3) - floor(x/6)*x//1":
        lambda x: (x//2)*(x//3) - (x//6)*(x//1),
}

for name, func in formulas.items():
    errors, exact = check_formula(func, x_range)
    if errors is not None:
        print(f"\n{name}:")
        print(f"  Max error: {max(abs(errors)):.1f}, exact: {exact}/{len(x_range)}")
        print(f"  First 20 errors: {list(errors[:20].astype(int))}")

# ============================================================
# Experiment 2: Least squares search for integer-coefficient formulas
# ============================================================
print("\n" + "=" * 70)
print("EXPERIMENT 2: Integer coefficient search")
print("=" * 70)

# Features: 1, floor(x/k) for k=1..K, floor(x/i)*floor(x/j) for i<=j<=K
# Find near-integer coefficients that give exact pi(x)

K = 10
x_range = list(range(2, 500))

features = []
feature_names = ["1"]
for k in range(1, K + 1):
    feature_names.append(f"f{k}")
for i in range(1, K + 1):
    for j in range(i, K + 1):
        feature_names.append(f"f{i}*f{j}")

for x in x_range:
    fv = [x // k for k in range(1, K + 1)]
    row = [1]
    row.extend(fv)
    for i in range(K):
        for j in range(i, K):
            row.append(fv[i] * fv[j])
    features.append(row)

A = np.array(features, dtype=float)
y = np.array([primepi(x) for x in x_range], dtype=float)

# Solve and look for near-integer coefficients
coeffs, _, _, _ = np.linalg.lstsq(A, y, rcond=None)

# Which coefficients are close to simple fractions?
print(f"\nCoefficients (K={K}, {len(feature_names)} features):")
for i, (name, c) in enumerate(zip(feature_names, coeffs)):
    if abs(c) > 1e-6:
        # Check if close to p/q for small q
        best_frac = None
        best_err = abs(c)
        for q in range(1, 13):
            p = round(c * q)
            err = abs(c - p/q)
            if err < best_err:
                best_err = err
                best_frac = (p, q)
        if best_frac and best_err < 0.01:
            print(f"  {name:12s} = {c:12.6f} ~ {best_frac[0]}/{best_frac[1]} (err={best_err:.6f})")
        elif abs(c) > 0.001:
            print(f"  {name:12s} = {c:12.6f}")

pred = A @ coeffs
max_err = np.max(np.abs(pred - y))
exact = sum(1 for i in range(len(x_range)) if abs(round(pred[i]) - y[i]) < 0.01)
print(f"\nFit quality: max_err={max_err:.4f}, exact={exact}/{len(x_range)}")

# ============================================================
# Experiment 3: Recursive floor identities
# ============================================================
print("\n" + "=" * 70)
print("EXPERIMENT 3: Recursive floor-based identities")
print("=" * 70)

# Known: pi(x) = pi(x-1) + [x is prime]
# Also: pi(x) = pi(x/2) + #{primes in (x/2, x]}
# And: #{primes in (y, 2y]} >= 1 for y >= 25 (Bertrand's postulate)

# Can we express pi(x) recursively using FEWER recursive calls?
# pi(x) = alpha * pi(floor(x/2)) + beta * pi(floor(x/3)) + gamma * floor(x/k) + ...

# Fit linear combination of pi(floor(x/k)) values
x_range = list(range(20, 500))

for K_pi in [2, 3, 4, 5]:
    # Features: pi(floor(x/k)) for k=1..K_pi, plus floor(x/k) for k=1..5
    feat = []
    for x in x_range:
        row = [1]
        for k in range(2, K_pi + 1):
            row.append(primepi(x // k))
        for k in range(1, 6):
            row.append(x // k)
        feat.append(row)

    A = np.array(feat, dtype=float)
    y = np.array([primepi(x) for x in x_range], dtype=float)

    coeffs, _, _, _ = np.linalg.lstsq(A, y, rcond=None)
    pred = A @ coeffs
    max_err = np.max(np.abs(pred - y))
    exact = sum(1 for i in range(len(x_range)) if abs(round(pred[i]) - y[i]) < 0.01)

    print(f"\npi(x) ~ sum c_k*pi(floor(x/k)) + sum d_k*floor(x/k), K_pi={K_pi}:")
    print(f"  max_err={max_err:.4f}, exact={exact}/{len(x_range)}")

    # Show coefficients
    names = ["1"] + [f"pi(x/{k})" for k in range(2, K_pi + 1)] + [f"f(x/{k})" for k in range(1, 6)]
    for name, c in zip(names, coeffs):
        if abs(c) > 0.001:
            print(f"    {name:12s}: {c:.6f}")

# ============================================================
# Experiment 4: Floor(x/k) difference sequences
# ============================================================
print("\n" + "=" * 70)
print("EXPERIMENT 4: Higher-order differences of floor(x/k)")
print("=" * 70)

# Define Delta_k f(x) = f(x) - f(x-k) (backward difference)
# For floor values: Delta_1 floor(x/k) = floor(x/k) - floor((x-1)/k)
# = [k | x] (1 if k divides x, 0 otherwise) -- this is well known

# Higher order: Delta_1^2 floor(x/k) = [k|x] - [k|(x-1)]
# = [k|x] - [k|x-1] which is 1 if k|x, -1 if k|x-1, else 0

# Product of differences:
# Delta_1 floor(x/k1) * Delta_1 floor(x/k2) = [k1|x] * [k2|x] = [lcm(k1,k2)|x]
# This is just divisibility testing!

for x_val in [100, 200, 500]:
    print(f"\nx={x_val}, pi(x)={primepi(x_val)}")

    # Sum_{k=2}^{x} Delta_1 floor(x/k) = sum_{k=2}^{x} [k|x] = d(x) - 1
    # where d(x) is the number of divisors of x
    divisor_sum = sum(1 for k in range(2, x_val + 1) if x_val % k == 0)
    print(f"  sum [k|x] for k=2..x = {divisor_sum} = d(x)-1")

    # Now: product_{k=2}^{sqrt(x)} (1 - [k|n]) for each n gives primality
    # But this is the sieve again.

    # Novel approach: look at SECOND differences
    # Delta_1 floor(x/k) at x and at x-1:
    # floor(x/k) - floor((x-1)/k) = [k|x]
    # floor((x-1)/k) - floor((x-2)/k) = [k|(x-1)]
    # Second difference: [k|x] - [k|(x-1)]

    # For prime p: x mod p cycles through 0,1,...,p-1
    # [p|x] = 1 every p steps. So the second difference is nonzero at x=mp and x=mp+1.

    # Cross-product: [k1|x]*[k2|x] = [lcm(k1,k2)|x]
    # For k1=p, k2=q primes: [pq|x] = 1 every pq steps
    # This detects multiples of pq, not primality.

# ============================================================
# Experiment 5: Can degree-2 polynomial in sqrt(x) floor values compute pi(x)?
# ============================================================
print("\n" + "=" * 70)
print("EXPERIMENT 5: Using only O(sqrt(x)) floor values with degree-2 polynomial")
print("=" * 70)

# Meissel-Lehmer uses O(x^{2/3}) floor values linearly.
# Can we use only O(sqrt(x)) floor values if we allow degree 2?

# The distinct floor values {floor(x/k) : k=1,...,x} have size O(sqrt(x))
# Call them v_1 > v_2 > ... > v_M where M ~ 2*sqrt(x)

# Can pi(x) = sum_{i,j} c_{i,j} * v_i * v_j + sum_i d_i * v_i + e ?

for x in [50, 100, 150]:
    pi_x = primepi(x)
    sqrt_x = int(math.sqrt(x))

    # Get distinct floor values
    floor_set = sorted(set(x // k for k in range(1, x + 1)), reverse=True)
    M = len(floor_set)

    # Build quadratic feature matrix for ONE value of x
    # But we need MULTIPLE x values to fit...
    #
    # Actually, the question is: for FIXED x, is pi(x) a function of
    # the floor value SET? The answer is yes -- pi(x) is determined by x.
    # But we need a formula that works for ALL x.

    print(f"\nx={x}: {M} distinct floor values, pi(x)={pi_x}")
    print(f"  Floor values: {floor_set[:20]}...")

# For a formula valid across x values, we need floor(x/k) for SPECIFIC k values
# (not the full set). The question is: how many k values, combined nonlinearly?

# Summary test: for K specific k-values, degree D, what's the error?
print("\n\nSummary: Error vs K (specific divisors) and degree")
print(f"{'K':>3} {'deg':>4} {'features':>8} {'train_err':>10} {'test_err':>10} {'train_exact':>12}")

for K in [5, 10, 15, 20]:
    x_train = list(range(10, 300))
    x_test = list(range(300, 500))

    for deg in [1, 2]:
        from itertools import combinations_with_replacement

        def make_features(x, K, deg):
            fv = [x // k for k in range(1, K + 1)]
            feats = [1]
            for d in range(1, deg + 1):
                for combo in combinations_with_replacement(range(K), d):
                    val = 1
                    for idx in combo:
                        val *= fv[idx]
                    feats.append(val)
            return feats

        n_feat = len(make_features(100, K, deg))
        if n_feat > len(x_train) - 20:
            continue

        A_train = np.array([make_features(x, K, deg) for x in x_train], dtype=float)
        y_train = np.array([primepi(x) for x in x_train], dtype=float)

        c, _, _, _ = np.linalg.lstsq(A_train, y_train, rcond=None)
        p_train = A_train @ c
        train_err = np.max(np.abs(p_train - y_train))
        train_exact = sum(1 for i in range(len(x_train)) if abs(round(p_train[i]) - y_train[i]) < 0.01)

        A_test = np.array([make_features(x, K, deg) for x in x_test], dtype=float)
        y_test = np.array([primepi(x) for x in x_test], dtype=float)
        p_test = A_test @ c
        test_err = np.max(np.abs(p_test - y_test))

        print(f"{K:3d} {deg:4d} {n_feat:8d} {train_err:10.3f} {test_err:10.3f} {train_exact:8d}/{len(x_train)}")

print("\nDone.")

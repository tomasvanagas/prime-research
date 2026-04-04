"""
Session 14: Polynomial identity approach to pi(x).

Core idea: Can pi(x) be expressed as a POLYNOMIAL in a small number of
"easy-to-compute" functions?

Candidate easy functions (all computable in TC^0 or NC^1):
- floor(x/k) for various k
- gcd(x, k) for various k
- x mod k for various k
- floor(log_2(x)), etc.
- Euler's totient phi(k) for small k (fixed, precomputed)

We know pi(x) IS a linear combination of floor values (the Legendre sieve),
but needs exponentially many. Can a POLYNOMIAL (degree > 1) in a small
number of values work?

Also explore: Can pi(x) be expressed as a resultant, discriminant, or
other algebraic invariant of polynomials whose coefficients depend on x?
"""

import numpy as np
from sympy import primepi, factorint, gcd, totient
from math import isqrt, log, floor
from itertools import combinations_with_replacement

def generate_features(x, feature_set='floor'):
    """Generate features (easy-to-compute functions of x) for polynomial fitting."""
    features = {}

    if feature_set in ('floor', 'all'):
        # Floor values for small divisors
        for k in range(1, min(50, x+1)):
            features[f'floor(x/{k})'] = x // k

    if feature_set in ('mod', 'all'):
        # x mod k for small k
        for k in range(2, 31):
            features[f'x mod {k}'] = x % k

    if feature_set in ('log', 'all'):
        # Logarithmic functions
        if x > 0:
            features['floor(log2(x))'] = int(log(x, 2)) if x > 0 else 0
            features['floor(sqrt(x))'] = isqrt(x)
            features['x'] = x

    if feature_set in ('gcd', 'all'):
        # GCD with small numbers (primorial, etc.)
        for k in [2, 6, 30, 210, 2310]:
            features[f'gcd(x,{k})'] = gcd(x, k)

    return features


def polynomial_regression(X_train, y_train, X_test, y_test, degree, feature_names):
    """
    Fit a polynomial of given degree in the features.
    Returns training and test errors.
    """
    from sklearn.preprocessing import PolynomialFeatures
    from sklearn.linear_model import LinearRegression

    poly = PolynomialFeatures(degree=degree, include_bias=True)
    X_train_poly = poly.fit_transform(X_train)
    X_test_poly = poly.transform(X_test)

    reg = LinearRegression(fit_intercept=False)
    reg.fit(X_train_poly, y_train)

    y_train_pred = reg.predict(X_train_poly)
    y_test_pred = reg.predict(X_test_poly)

    train_exact = sum(1 for a, b in zip(y_train, y_train_pred) if abs(a - round(b)) == 0)
    test_exact = sum(1 for a, b in zip(y_test, y_test_pred) if abs(a - round(b)) == 0)

    train_mae = np.mean(np.abs(y_train - y_train_pred))
    test_mae = np.mean(np.abs(y_test - y_test_pred))

    return {
        'train_exact': train_exact, 'test_exact': test_exact,
        'train_mae': train_mae, 'test_mae': test_mae,
        'n_features': X_train_poly.shape[1],
        'coeff_norm': np.linalg.norm(reg.coef_)
    }


def manual_polynomial_features(x_vals, feature_names_list, degree):
    """Generate polynomial features manually (no sklearn dependency)."""
    n = len(x_vals)
    k = len(feature_names_list)

    # For degree 1, just the features
    if degree == 1:
        X = np.array([[x_vals[i][j] for j in range(k)] for i in range(n)])
        names = feature_names_list[:]
        # Add bias
        X = np.column_stack([np.ones(n), X])
        names = ['1'] + names
        return X, names

    # For degree 2, add pairwise products
    if degree == 2:
        X_list = []
        names = ['1']
        X_list.append(np.ones(n))
        for j in range(k):
            X_list.append(np.array([x_vals[i][j] for i in range(n)]))
            names.append(feature_names_list[j])
        for j1 in range(k):
            for j2 in range(j1, k):
                X_list.append(np.array([x_vals[i][j1] * x_vals[i][j2] for i in range(n)]))
                names.append(f'{feature_names_list[j1]}*{feature_names_list[j2]}')
        X = np.column_stack(X_list)
        return X, names

    return None, None


def test_floor_polynomial():
    """Test polynomial combinations of floor values."""
    print("=" * 70)
    print("Testing polynomial combinations of floor(x/k) values")
    print("=" * 70)

    # Generate data
    x_range = range(10, 501)
    x_train = list(range(10, 301))
    x_test = list(range(301, 501))

    for n_floors in [5, 10, 20]:
        floor_ks = list(range(1, n_floors + 1))
        feature_names = [f'fl(x/{k})' for k in floor_ks]

        # Build feature matrices
        X_train_data = [[x // k for k in floor_ks] for x in x_train]
        X_test_data = [[x // k for k in floor_ks] for x in x_test]
        y_train = np.array([int(primepi(x)) for x in x_train], dtype=float)
        y_test = np.array([int(primepi(x)) for x in x_test], dtype=float)

        for degree in [1, 2]:
            X_train, names = manual_polynomial_features(X_train_data, feature_names, degree)
            X_test, _ = manual_polynomial_features(X_test_data, feature_names, degree)

            if X_train is None:
                continue

            # Fit
            try:
                coeffs, residuals, rank, sv = np.linalg.lstsq(X_train, y_train, rcond=None)
                y_train_pred = X_train @ coeffs
                y_test_pred = X_test @ coeffs

                train_exact = sum(1 for a, b in zip(y_train, y_train_pred) if abs(a - round(b)) == 0)
                test_exact = sum(1 for a, b in zip(y_test, y_test_pred) if abs(a - round(b)) == 0)
                train_mae = np.mean(np.abs(y_train - y_train_pred))
                test_mae = np.mean(np.abs(y_test - y_test_pred))

                print(f"\n{n_floors} floor values, degree {degree}: {len(names)} features")
                print(f"  Train: MAE={train_mae:.3f}, exact={train_exact}/{len(x_train)}")
                print(f"  Test:  MAE={test_mae:.3f}, exact={test_exact}/{len(x_test)}")
                print(f"  Coefficient norm: {np.linalg.norm(coeffs):.2e}")
            except Exception as e:
                print(f"\n{n_floors} floor values, degree {degree}: FAILED ({e})")


def test_mixed_features():
    """Test polynomial combinations of mixed features (floor + mod + log)."""
    print("\n" + "=" * 70)
    print("Testing mixed feature polynomial regression")
    print("=" * 70)

    x_train = list(range(10, 301))
    x_test = list(range(301, 501))

    y_train = np.array([int(primepi(x)) for x in x_train], dtype=float)
    y_test = np.array([int(primepi(x)) for x in x_test], dtype=float)

    # Feature sets
    feature_configs = [
        ("floor(x/k) k=1..10", lambda x: [x // k for k in range(1, 11)]),
        ("floor + mod6 + sqrt", lambda x: [x // k for k in range(1, 8)] + [x % 6, isqrt(x)]),
        ("floor + mod30 + log", lambda x: [x // k for k in range(1, 8)] + [x % 30, int(log(x, 2)) if x > 0 else 0, isqrt(x)]),
        ("R^-1(x) + floor(x/k) k=1..5", lambda x: [round(x / log(x) * (1 + 1/log(x) + 2/log(x)**2)) if x > 1 else 0] + [x // k for k in range(1, 6)]),
    ]

    for name, feat_fn in feature_configs:
        n_feat = len(feat_fn(100))
        feat_names = [f'f{i}' for i in range(n_feat)]

        X_train_data = [feat_fn(x) for x in x_train]
        X_test_data = [feat_fn(x) for x in x_test]

        for degree in [1, 2]:
            X_train, names = manual_polynomial_features(X_train_data, feat_names, degree)
            X_test, _ = manual_polynomial_features(X_test_data, feat_names, degree)

            if X_train is None:
                continue

            try:
                coeffs, _, _, _ = np.linalg.lstsq(X_train, y_train, rcond=None)
                y_train_pred = X_train @ coeffs
                y_test_pred = X_test @ coeffs

                train_exact = sum(1 for a, b in zip(y_train, y_train_pred) if abs(a - round(b)) == 0)
                test_exact = sum(1 for a, b in zip(y_test, y_test_pred) if abs(a - round(b)) == 0)
                train_mae = np.mean(np.abs(y_train - y_train_pred))
                test_mae = np.mean(np.abs(y_test - y_test_pred))

                print(f"\n{name}, degree {degree}: {len(names)} features")
                print(f"  Train: MAE={train_mae:.3f}, exact={train_exact}/{len(x_train)}")
                print(f"  Test:  MAE={test_mae:.3f}, exact={test_exact}/{len(x_test)}")
            except Exception as e:
                print(f"\n{name}, degree {degree}: FAILED ({e})")


def test_prime_indicator_polynomial():
    """
    Can the prime indicator 1_prime(n) be expressed as a polynomial (mod 2)
    in a small number of bits/features of n?

    Over GF(2), Session 13 showed ANF degree = Theta(N). But what about
    over Z? The prime indicator as a function {0,1}^N -> {0,1} might have
    a polynomial representation over Z with lower degree.
    """
    print("\n" + "=" * 70)
    print("Prime indicator as polynomial over Z")
    print("=" * 70)

    from sympy import isprime

    for N in [5, 6, 7, 8]:
        max_x = 2**N
        x_vals = list(range(2, max_x))
        y_vals = [1 if isprime(x) else 0 for x in x_vals]

        # Features: the N bits of x
        def bits(x, N):
            return [(x >> i) & 1 for i in range(N)]

        X_data = [bits(x, N) for x in x_vals]
        feat_names = [f'b{i}' for i in range(N)]

        for degree in [1, 2]:
            X, names = manual_polynomial_features(X_data, feat_names, degree)
            if X is None:
                continue

            y = np.array(y_vals, dtype=float)
            coeffs, residuals, rank, sv = np.linalg.lstsq(X, y, rcond=None)
            y_pred = X @ coeffs
            exact = sum(1 for a, b in zip(y, y_pred) if abs(a - round(b)) == 0)
            mae = np.mean(np.abs(y - y_pred))

            print(f"\nN={N} (x up to {max_x}), degree {degree}: {len(names)} features, rank={rank}")
            print(f"  MAE: {mae:.4f}, exact: {exact}/{len(x_vals)} ({100*exact/len(x_vals):.1f}%)")
            print(f"  Condition number: {sv[0]/sv[-1] if sv[-1] > 0 else 'inf':.2e}")


def test_log_floor_identity():
    """
    Test whether there's a POLYNOMIAL IDENTITY involving floor(x/k),
    floor(log(x)), etc. that equals pi(x) EXACTLY.

    Strategy: for x in a range, compute floor values and pi(x), then
    check if the residual after subtracting the best linear combination
    has any structure.
    """
    print("\n" + "=" * 70)
    print("Residual structure analysis")
    print("=" * 70)

    x_range = list(range(10, 500))
    y = np.array([int(primepi(x)) for x in x_range], dtype=float)

    # Best linear combination of floor(x/k) for k=1..20
    floor_feats = np.array([[x // k for k in range(1, 21)] for x in x_range], dtype=float)
    X = np.column_stack([np.ones(len(x_range)), floor_feats])
    coeffs, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
    y_pred = X @ coeffs
    residual = y - y_pred

    print(f"Linear combination of 20 floor values:")
    print(f"  MAE of residual: {np.mean(np.abs(residual)):.4f}")
    print(f"  Max residual: {max(np.abs(residual)):.4f}")
    print(f"  Std of residual: {np.std(residual):.4f}")

    # Analyze residual: is it correlated with any simple function?
    x_arr = np.array(x_range, dtype=float)
    test_funcs = {
        'sqrt(x)': np.sqrt(x_arr),
        'x/ln(x)': x_arr / np.log(x_arr),
        'ln(x)': np.log(x_arr),
        'x mod 6': x_arr % 6,
        'x mod 30': x_arr % 30,
        'isqrt(x)': np.array([isqrt(x) for x in x_range], dtype=float),
        'li(x) - pi(x)': np.array([x/log(x) - int(primepi(x)) for x in x_range]),
    }

    print("\nCorrelation of residual with candidate functions:")
    for name, func in test_funcs.items():
        corr = np.corrcoef(residual, func)[0, 1]
        print(f"  {name:20s}: r = {corr:+.4f}")

    # Autocorrelation of residual
    res_centered = residual - residual.mean()
    if np.std(res_centered) > 1e-10:
        acf = np.correlate(res_centered, res_centered, mode='full')
        acf = acf[len(res_centered)-1:]
        acf /= acf[0]
        print(f"\nAutocorrelation of residual (lags 1-10): {acf[1:11].round(4)}")


if __name__ == '__main__':
    test_floor_polynomial()
    test_mixed_features()
    test_prime_indicator_polynomial()
    test_log_floor_identity()

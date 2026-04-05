"""
PROPOSAL 7: Neural Delta Oracle + Certification
================================================

IDEA: A completely different paradigm. Instead of computing delta(n) analytically,
LEARN it from data using a neural network, then CERTIFY the answer.

Key insight: even if we can't COMPUTE delta(n) in polylog time,
we might be able to PREDICT it, then VERIFY the prediction cheaply.

STEP 1 - PREDICTION:
Train a neural network f(n) that approximates delta(n) = p(n) - R^{-1}(n).
Input features: n, log(n), log(log(n)), fractional parts of n*gamma_k for
several zeta zeros gamma_k, etc.

The network needs to be accurate to within 1/2 (to round to the correct integer
part of the correction).

STEP 2 - CERTIFICATION:
Given x_candidate = round(R^{-1}(n) + f(n)):
1. Test if x_candidate is prime (AKS: O(log(x)^6), or probabilistic: O(log(x)^2))
2. Verify pi(x_candidate) = n (the hard part!)

For Step 2, we need a CHEAP VERIFICATION that pi(x) = n.
Option A: Verify pi(x_candidate) - pi(x_candidate - g) = 1 where g is the gap,
  and pi(x_candidate - g) is known from a previous computation.
Option B: Use a PRIMALITY CERTIFICATE that encodes pi(x).
Option C: Verify by checking no primes exist between x_candidate and the
  previous prime (only need primality testing of O(log(x)) numbers).

RADICAL OBSERVATION: If we know p(n-1), finding p(n) just requires
finding the next prime. The next prime after x is at most x + O(log(x)^2)
away (by Cramér's conjecture). Testing O(log(x)^2) numbers for primality
is O(log(x)^2 * log(x)^6) = O(log(x)^8) by AKS.

SO THE REAL PROBLEM IS BOOTSTRAPPING: if we could compute p(n0) for some
n0, then p(n0+1), p(n0+2), ... are each O(polylog) apart. But
accumulating from n0 to n takes O(n - n0) steps.

NEURAL SHORTCUT: If the neural net can predict delta(n) to within the
prime gap (which is O(log(x))), we get a CANDIDATE. Then:
1. Test primality of candidate: O(polylog)
2. If prime, count primes in [candidate-gap, candidate] to verify rank
3. This local counting is O(gap * polylog) = O(log(x) * polylog)

TOTAL: O(polylog) for prediction + O(polylog) for verification = O(polylog)
IF the neural net generalizes.

THE CATCH: Does delta(n) have learnable structure?
Test this empirically!
"""

import numpy as np
from sympy import prime, primepi, isprime, nextprime, prevprime, primerange
from mpmath import mp, mpf, li, log, exp
import math

mp.dps = 30

def compute_features(n, zeros_gamma=None):
    """
    Compute feature vector for predicting delta(n).

    Features motivated by the explicit formula:
    - n, log(n), log(log(n))
    - Fractional parts of gamma_k * log(n*log(n)) for zeta zeros gamma_k
    - These capture the oscillatory contributions of the zeros
    """
    if zeros_gamma is None:
        # First 20 zeta zero imaginary parts
        zeros_gamma = [
            14.134725, 21.022040, 25.010858, 30.424876, 32.935062,
            37.586178, 40.918719, 43.327073, 48.005151, 49.773832,
            52.970321, 56.446248, 59.347044, 60.831779, 65.112544,
            67.079811, 69.546402, 72.067158, 75.704691, 77.144840,
        ]

    ln_n = math.log(max(n, 2))
    ln_ln_n = math.log(max(ln_n, 1.01))

    # Approximate x ~ n * ln(n)
    x_approx = n * ln_n
    ln_x = math.log(max(x_approx, 2))

    features = [
        n,
        ln_n,
        ln_ln_n,
        n * ln_ln_n,
        math.sqrt(x_approx),
        1.0 / ln_n,
    ]

    # Oscillatory features from zeta zeros
    for gamma in zeros_gamma:
        phase = gamma * ln_x
        features.append(math.cos(phase))
        features.append(math.sin(phase))

    return np.array(features)

def build_training_data(n_range, zeros_gamma=None):
    """Build training data for the neural delta oracle."""
    from sympy import mobius

    X = []
    y = []

    for n in range(max(2, n_range[0]), n_range[1] + 1):
        p_n = int(prime(n))

        # Compute R^{-1}(n) approximately
        x = mpf(n) * log(mpf(n))
        for _ in range(20):
            R_val = mpf(0)
            for k in range(1, 30):
                mu_k = mobius(k)
                if mu_k != 0:
                    R_val += mpf(mu_k) / k * li(x ** (mpf(1)/k))
            R_prime = 1 / log(x)
            x = x + (mpf(n) - R_val) / R_prime
            if abs(mpf(n) - R_val) < mpf('1e-10'):
                break

        delta = p_n - float(x)
        features = compute_features(n, zeros_gamma)

        X.append(features)
        y.append(delta)

    return np.array(X), np.array(y)

def simple_regression_test(X_train, y_train, X_test, y_test):
    """
    Test a simple linear regression on the features.
    If even linear regression works, the structure is learnable.
    """
    from numpy.linalg import lstsq

    # Normalize features
    mean = X_train.mean(axis=0)
    std = X_train.std(axis=0)
    std[std == 0] = 1

    X_norm = (X_train - mean) / std
    X_test_norm = (X_test - mean) / std

    # Add bias
    X_aug = np.column_stack([X_norm, np.ones(len(X_norm))])
    X_test_aug = np.column_stack([X_test_norm, np.ones(len(X_test_norm))])

    # Solve least squares
    coeffs, residuals, rank, sv = lstsq(X_aug, y_train, rcond=None)

    # Predictions
    y_pred_train = X_aug @ coeffs
    y_pred_test = X_test_aug @ coeffs

    train_rmse = np.sqrt(np.mean((y_train - y_pred_train)**2))
    test_rmse = np.sqrt(np.mean((y_test - y_pred_test)**2))

    # How often would rounding give the correct delta?
    train_correct = np.mean(np.round(y_pred_train) == np.round(y_train))
    test_correct = np.mean(np.round(y_pred_test) == np.round(y_test))

    return {
        'train_rmse': train_rmse,
        'test_rmse': test_rmse,
        'train_round_accuracy': train_correct,
        'test_round_accuracy': test_correct,
        'coefficients': coeffs,
    }

def verification_cost_analysis(n_values):
    """
    Analyze the cost of verifying a candidate p(n).

    Given x_candidate:
    1. Is x_candidate prime? Cost: O(log(x)^6) via AKS
    2. Is it the RIGHT prime? Need to check no primes between p(n-1) and x_candidate

    The gap g = p(n) - p(n-1) determines verification cost.
    By Cramér's conjecture, g = O(log(p)^2).
    Checking g numbers for primality costs O(g * log(x)^2) (Miller-Rabin).
    """
    results = []
    for n in n_values:
        p_n = int(prime(n))
        p_prev = int(prime(n - 1)) if n > 1 else 1
        gap = p_n - p_prev
        log_p = math.log(p_n)

        # Verification cost (assuming candidate is correct)
        # Check gap numbers for primality
        verification_ops = gap * log_p**2  # Miller-Rabin cost
        cramer_expected_gap = log_p**2

        results.append({
            'n': n,
            'p_n': p_n,
            'gap': gap,
            'log_p': log_p,
            'cramer_expected': cramer_expected_gap,
            'verification_ops': verification_ops,
        })

    return results

def gap_prediction_test(n_range):
    """
    Separate test: can we predict the next prime gap from features?
    If yes, this gives another O(polylog) shortcut.
    """
    gaps = []
    features = []

    for n in range(max(3, n_range[0]), n_range[1] + 1):
        p_n = int(prime(n))
        p_prev = int(prime(n - 1))
        gap = p_n - p_prev

        # Features for gap prediction
        f = [
            math.log(p_prev),
            p_prev % 6,  # primes are ≡ 1 or 5 mod 6
            p_prev % 30,
            gap,  # previous gap
        ]

        gaps.append(gap)
        features.append(f)

    gaps = np.array(gaps)
    features = np.array(features)

    # Can we predict next gap from current features?
    # Shift: predict gaps[i+1] from features[i]
    X = features[:-1]
    y = gaps[1:]

    # Split
    split = int(0.7 * len(X))
    X_train, X_test = X[:split], X[split:]
    y_train, y_test = y[:split], y[split:]

    from numpy.linalg import lstsq
    X_aug = np.column_stack([X_train, np.ones(len(X_train))])
    X_test_aug = np.column_stack([X_test, np.ones(len(X_test))])

    coeffs, _, _, _ = lstsq(X_aug, y_train, rcond=None)

    y_pred = X_test_aug @ coeffs
    rmse = np.sqrt(np.mean((y_test - y_pred)**2))
    baseline_rmse = np.sqrt(np.mean((y_test - np.mean(y_train))**2))

    exact_match = np.mean(np.round(y_pred) == y_test)

    return {
        'rmse': rmse,
        'baseline_rmse': baseline_rmse,
        'improvement': 1 - rmse/baseline_rmse,
        'exact_match_rate': exact_match,
        'mean_gap': np.mean(gaps),
        'std_gap': np.std(gaps),
    }


if __name__ == '__main__':
    print("=" * 80)
    print("PROPOSAL 7: Neural Delta Oracle + Certification")
    print("=" * 80)

    print("\n--- Part A: Build training data (n=2..500) ---")
    X_train, y_train = build_training_data((2, 400))
    X_test, y_test = build_training_data((401, 500))
    print(f"  Training: {len(X_train)} samples, {X_train.shape[1]} features")
    print(f"  Test: {len(X_test)} samples")
    print(f"  delta stats: mean={np.mean(y_train):.2f}, std={np.std(y_train):.2f}")

    print("\n--- Part B: Linear regression on delta ---")
    result = simple_regression_test(X_train, y_train, X_test, y_test)
    print(f"  Train RMSE: {result['train_rmse']:.4f}")
    print(f"  Test RMSE: {result['test_rmse']:.4f}")
    print(f"  Train round accuracy: {result['train_round_accuracy']:.4f}")
    print(f"  Test round accuracy: {result['test_round_accuracy']:.4f}")

    print("\n--- Part C: Verification cost analysis ---")
    verify_results = verification_cost_analysis([100, 500, 1000, 5000, 10000])
    for r in verify_results:
        print(f"  n={r['n']:5d}: gap={r['gap']:3d}, log(p)={r['log_p']:.1f}, "
              f"Cramér_expected={r['cramer_expected']:.1f}, "
              f"verify_ops={r['verification_ops']:.0f}")

    print("\n--- Part D: Gap prediction test ---")
    gap_result = gap_prediction_test((3, 2000))
    print(f"  Gap RMSE: {gap_result['rmse']:.4f} (baseline: {gap_result['baseline_rmse']:.4f})")
    print(f"  Improvement: {gap_result['improvement']:.4f}")
    print(f"  Exact match rate: {gap_result['exact_match_rate']:.4f}")
    print(f"  Mean gap: {gap_result['mean_gap']:.2f}, std: {gap_result['std_gap']:.2f}")

    print("\n--- Part E: Feature importance (which zeta zeros matter most?) ---")
    # Look at coefficient magnitudes for the zero-based features
    coeffs = result['coefficients']
    feature_names = ['n', 'ln_n', 'ln_ln_n', 'n*ln_ln_n', 'sqrt_x', '1/ln_n']
    for i in range(20):
        feature_names.append(f'cos(g{i+1}*ln_x)')
        feature_names.append(f'sin(g{i+1}*ln_x)')
    feature_names.append('bias')

    # Sort by magnitude
    coeff_importance = [(abs(coeffs[i]), feature_names[i], coeffs[i])
                        for i in range(len(coeffs))]
    coeff_importance.sort(reverse=True)
    print("  Top 15 features by coefficient magnitude:")
    for mag, name, coeff in coeff_importance[:15]:
        print(f"    {name:20s}: coeff={coeff:+10.4f}")

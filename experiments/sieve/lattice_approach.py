"""
Session 6: Lattice-Based Approach to Prime Prediction

RADICAL IDEA: Use LLL lattice reduction to find integer relations
that express p(n) in terms of known quantities.

Approach 1: PSLQ/LLL Integer Relation Detection
  Given p(n), find integers a_i such that:
  a_0 * p(n) + a_1 * n*ln(n) + a_2 * n*ln(ln(n)) + a_3 * li^{-1}(n) + ... = 0
  If such a relation exists with SMALL a_i, we have a formula.

Approach 2: Simultaneous Approximation
  Find a single polynomial/expression that simultaneously approximates
  p(n) for many values of n. Use LLL to find the best coefficients.

Approach 3: Modular Lattice
  For each small prime q, p(n) mod q is determined by the position of
  p(n) among residue classes. Can LLL help reconstruct p(n) from these?

Approach 4: Diophantine Approximation
  Use LLL to find the best rational approximation to the "fractional part"
  of p(n) / some_function(n). If this has bounded partial quotients,
  we might predict it.

Approach 5: Feature Lattice
  Create a lattice where each basis vector corresponds to a number-theoretic
  function evaluated at n (li(n), R(n), n*ln(n), ...). Find short vectors
  that simultaneously approximate p(n) - R^{-1}(n) for many n.
"""

import numpy as np
from sympy import prime, primepi, log as slog, li as sli, floor as sfloor
from sympy import Rational, Matrix
import time

def compute_features(n):
    """Compute various number-theoretic features for index n."""
    if n < 2:
        return None
    ln_n = np.log(n)
    ln_ln_n = np.log(ln_n) if ln_n > 0 else 0

    features = {
        'n': n,
        'n_ln_n': n * ln_n,
        'n_ln_n_ln_ln_n': n * (ln_n + ln_ln_n),
        'cipolla_1': n * (ln_n + ln_ln_n - 1),
        'cipolla_2': n * (ln_n + ln_ln_n - 1 + (ln_ln_n - 2) / ln_n),
        'cipolla_3': n * (ln_n + ln_ln_n - 1 + (ln_ln_n - 2) / ln_n +
                         (ln_ln_n**2 - 6*ln_ln_n + 11) / (2 * ln_n**2)),
        'sqrt_n': np.sqrt(n),
        'n_over_ln_n': n / ln_n,
        'harmonic': np.log(n) + 0.5772156649,  # H_n ≈ ln(n) + γ
    }
    return features

def pslq_integer_relation(values, precision_bits=100):
    """
    Simple LLL-based integer relation detection.
    Given real numbers x_1, ..., x_m, find integers a_1, ..., a_m such that
    Σ a_i * x_i ≈ 0.
    """
    m = len(values)
    # Scale factor
    C = 10**15  # Large integer to scale real values

    # Build the LLL matrix
    # [1  0  0 ... 0  round(C*x_1)]
    # [0  1  0 ... 0  round(C*x_2)]
    # ...
    # [0  0  0 ... 1  round(C*x_m)]
    mat = np.zeros((m, m+1), dtype=object)  # Use Python ints to avoid overflow
    for i in range(m):
        mat[i, i] = int(1)
        mat[i, m] = int(round(C * float(values[i])))

    # Use sympy's LLL
    M = Matrix(mat.tolist())
    try:
        reduced = M.lll()
        # The shortest vector gives the integer relation
        # Check first few rows for short vectors
        results = []
        for i in range(min(3, reduced.rows)):
            row = [int(reduced[i, j]) for j in range(reduced.cols)]
            coeffs = row[:-1]
            residual = row[-1]
            dot = sum(c * v for c, v in zip(coeffs, values))
            results.append({
                'coeffs': coeffs,
                'residual_scaled': residual,
                'actual_residual': dot,
                'max_coeff': max(abs(c) for c in coeffs),
                'norm': np.sqrt(sum(c**2 for c in coeffs))
            })
        return results
    except Exception as e:
        return [{'error': str(e)}]

def experiment_1_integer_relations():
    """Find integer relations between p(n) and various functions of n."""
    print("="*70)
    print("EXPERIMENT 1: INTEGER RELATION DETECTION (PSLQ/LLL)")
    print("="*70)

    # For each n, try to find a relation:
    # a * p(n) = b * f1(n) + c * f2(n) + d * f3(n) + ...
    # Equivalently: a * p(n) + b * f1(n) + c * f2(n) + ... = 0

    for n in [100, 500, 1000]:
        p_n = int(prime(n))
        feats = compute_features(n)
        if feats is None:
            continue

        values = [
            float(p_n),           # p(n)
            feats['cipolla_2'],   # Cipolla approx
            feats['cipolla_3'],   # Higher order
            feats['n_ln_n'],      # n*ln(n)
            feats['sqrt_n'],      # sqrt(n)
            1.0                   # constant term
        ]

        labels = ['p(n)', 'cipolla_2', 'cipolla_3', 'n*ln(n)', 'sqrt(n)', '1']

        results = pslq_integer_relation(values)
        print(f"\nn = {n}, p({n}) = {p_n}")
        for r in results[:2]:
            if 'error' in r:
                print(f"  Error: {r['error']}")
                continue
            relation = ' + '.join(f'{c}*{l}' for c, l in zip(r['coeffs'], labels) if c != 0)
            print(f"  Relation: {relation} ≈ {r['actual_residual']:.6f}")
            print(f"  Max coeff: {r['max_coeff']}, norm: {r['norm']:.2f}")

def experiment_2_simultaneous_approximation():
    """
    Find a SINGLE set of coefficients that works for ALL n.

    We want: p(n) ≈ c_0 + c_1*f_1(n) + c_2*f_2(n) + ... for all n.

    Use least squares first, then LLL to find integer/rational coefficients.
    """
    print("\n" + "="*70)
    print("EXPERIMENT 2: SIMULTANEOUS APPROXIMATION")
    print("="*70)

    N_train = 5000
    N_test = 5000

    # Build feature matrix
    train_n = list(range(2, N_train + 2))
    test_n = list(range(N_train + 2, N_train + N_test + 2))

    def make_features(n_val):
        ln_n = np.log(n_val)
        ln_ln_n = np.log(max(ln_n, 1.01))
        return [
            n_val * ln_n,                        # n*ln(n)
            n_val * ln_ln_n,                      # n*ln(ln(n))
            n_val,                                # n
            n_val * (ln_ln_n - 2) / ln_n,         # first correction
            n_val * (ln_ln_n**2 - 6*ln_ln_n + 11) / (2*ln_n**2),  # second correction
            n_val * ln_ln_n**2 / ln_n**2,         # extra term
            np.sqrt(n_val) * ln_n,                # sqrt term
            1.0,                                  # constant
        ]

    feature_names = ['n*ln(n)', 'n*lnln(n)', 'n', 'n*(lnln-2)/ln',
                     'n*(lnln²-6lnln+11)/(2ln²)', 'n*lnln²/ln²',
                     'sqrt(n)*ln(n)', '1']

    # Build training data
    X_train = np.array([make_features(n) for n in train_n])
    y_train = np.array([int(prime(n)) for n in train_n], dtype=np.float64)

    # Least squares fit
    coeffs, residuals, rank, sv = np.linalg.lstsq(X_train, y_train, rcond=None)

    print("Least squares coefficients:")
    for name, c in zip(feature_names, coeffs):
        print(f"  {name}: {c:.10f}")

    # Training accuracy
    y_pred_train = X_train @ coeffs
    errors_train = y_pred_train - y_train
    print(f"\nTraining (n=2..{N_train+1}):")
    print(f"  Mean error: {np.mean(errors_train):.4f}")
    print(f"  Max |error|: {np.max(np.abs(errors_train)):.4f}")
    print(f"  RMS error: {np.sqrt(np.mean(errors_train**2)):.4f}")
    exact_train = np.sum(np.abs(errors_train) < 0.5)
    print(f"  Exact (round to correct): {exact_train}/{len(train_n)} ({100*exact_train/len(train_n):.1f}%)")

    # Test accuracy
    X_test = np.array([make_features(n) for n in test_n])
    y_test = np.array([int(prime(n)) for n in test_n], dtype=np.float64)
    y_pred_test = X_test @ coeffs
    errors_test = y_pred_test - y_test
    print(f"\nTest (n={N_train+2}..{N_train+N_test+1}):")
    print(f"  Mean error: {np.mean(errors_test):.4f}")
    print(f"  Max |error|: {np.max(np.abs(errors_test)):.4f}")
    print(f"  RMS error: {np.sqrt(np.mean(errors_test**2)):.4f}")
    exact_test = np.sum(np.abs(errors_test) < 0.5)
    print(f"  Exact (round to correct): {exact_test}/{len(test_n)} ({100*exact_test/len(test_n):.1f}%)")

    return coeffs

def experiment_3_correction_lattice():
    """
    Use LLL to find a formula for the CORRECTION δ(n) = p(n) - approx(n).

    If δ(n) can be expressed as a linear combination of known functions
    with small integer coefficients, LLL will find it.
    """
    print("\n" + "="*70)
    print("EXPERIMENT 3: CORRECTION TERM LATTICE REDUCTION")
    print("="*70)

    # Compute δ(n) for a range of n
    N = 200
    ns = list(range(10, N + 10))
    deltas = []
    for n in ns:
        p_n = int(prime(n))
        ln_n = np.log(n)
        ln_ln_n = np.log(ln_n)
        approx = n * (ln_n + ln_ln_n - 1 + (ln_ln_n - 2) / ln_n)
        deltas.append(p_n - approx)

    deltas = np.array(deltas)

    print(f"Correction δ(n) statistics for n={ns[0]}..{ns[-1]}:")
    print(f"  Mean: {np.mean(deltas):.4f}")
    print(f"  Std: {np.std(deltas):.4f}")
    print(f"  Range: [{np.min(deltas):.4f}, {np.max(deltas):.4f}]")

    # Try to express δ(n) as a function of n
    # Features: 1, ln(n), ln(ln(n)), n^{1/3}, sin(ln(n)*γ_1) where γ_1 is first zeta zero
    gamma_1 = 14.134725
    features_for_delta = []
    for n in ns:
        ln_n = np.log(n)
        ln_ln_n = np.log(ln_n)
        features_for_delta.append([
            1,
            ln_n,
            ln_ln_n,
            n**(1/3),
            np.sin(gamma_1 * ln_n),
            np.cos(gamma_1 * ln_n),
            np.sin(21.022 * ln_n),  # second zero
            np.cos(21.022 * ln_n),
        ])

    X = np.array(features_for_delta)
    coeffs, _, _, _ = np.linalg.lstsq(X, deltas, rcond=None)

    feature_names = ['1', 'ln(n)', 'ln(ln(n))', 'n^{1/3}',
                     'sin(γ₁·ln(n))', 'cos(γ₁·ln(n))',
                     'sin(γ₂·ln(n))', 'cos(γ₂·ln(n))']

    print("\nBest fit for δ(n):")
    for name, c in zip(feature_names, coeffs):
        if abs(c) > 0.001:
            print(f"  {c:+.6f} * {name}")

    residuals = deltas - X @ coeffs
    print(f"\nResidual after fitting:")
    print(f"  RMS: {np.sqrt(np.mean(residuals**2)):.4f}")
    print(f"  Max |residual|: {np.max(np.abs(residuals)):.4f}")
    exact = np.sum(np.abs(X @ coeffs - deltas) < 0.5)
    print(f"  Exact δ predictions: {exact}/{len(ns)} ({100*exact/len(ns):.1f}%)")

    # Now try on UNSEEN data
    test_ns = list(range(N + 10, N + 110))
    test_deltas = []
    for n in test_ns:
        p_n = int(prime(n))
        ln_n = np.log(n)
        ln_ln_n = np.log(ln_n)
        approx = n * (ln_n + ln_ln_n - 1 + (ln_ln_n - 2) / ln_n)
        test_deltas.append(p_n - approx)

    test_X = []
    for n in test_ns:
        ln_n = np.log(n)
        ln_ln_n = np.log(ln_n)
        test_X.append([1, ln_n, ln_ln_n, n**(1/3),
                       np.sin(gamma_1 * ln_n), np.cos(gamma_1 * ln_n),
                       np.sin(21.022 * ln_n), np.cos(21.022 * ln_n)])

    test_X = np.array(test_X)
    test_deltas = np.array(test_deltas)
    test_pred = test_X @ coeffs
    test_residuals = test_deltas - test_pred
    exact_test = np.sum(np.abs(test_residuals) < 0.5)
    print(f"\nTest set (n={test_ns[0]}..{test_ns[-1]}):")
    print(f"  RMS residual: {np.sqrt(np.mean(test_residuals**2)):.4f}")
    print(f"  Exact: {exact_test}/{len(test_ns)} ({100*exact_test/len(test_ns):.1f}%)")

def experiment_4_digit_extraction():
    """
    Can we compute individual BITS of p(n) independently?

    If p(n) = Σ b_k * 2^k, can we find b_k for each k separately?
    This would be a BBP-type formula for the prime sequence.
    """
    print("\n" + "="*70)
    print("EXPERIMENT 4: BIT-LEVEL ANALYSIS")
    print("="*70)

    # Analyze bit patterns of p(n) for n=1..1000
    primes = [int(prime(n)) for n in range(1, 1001)]

    # Look at each bit position
    max_bits = 15  # primes up to ~32000 need ~15 bits
    bit_probs = []
    for bit in range(max_bits):
        ones = sum(1 for p in primes if (p >> bit) & 1)
        prob = ones / len(primes)
        bit_probs.append(prob)
        if bit < 8:
            print(f"  Bit {bit}: P(1) = {prob:.4f}")

    # Are there correlations between consecutive bits?
    print("\nBit-bit correlations:")
    bit_arrays = np.array([[(p >> bit) & 1 for p in primes] for bit in range(max_bits)])
    for i in range(min(5, max_bits)):
        for j in range(i+1, min(5, max_bits)):
            corr = np.corrcoef(bit_arrays[i], bit_arrays[j])[0, 1]
            if abs(corr) > 0.05:
                print(f"  bits ({i},{j}): r = {corr:.4f}")

    # Look at p(n) mod 6 (all primes > 3 are 1 or 5 mod 6)
    print("\np(n) mod 6 distribution (n > 2):")
    mod6 = [p % 6 for p in primes[2:]]
    for r in range(6):
        count = mod6.count(r)
        if count > 0:
            print(f"  {r}: {count} ({100*count/len(mod6):.1f}%)")

    # Can we predict p(n) mod 6?
    # Dirichlet: asymptotically equal, but are there deviations?
    print("\nCumulative bias p(n) mod 6:")
    cum_1 = 0
    cum_5 = 0
    for i, p in enumerate(primes[2:], start=3):
        if p % 6 == 1:
            cum_1 += 1
        else:
            cum_5 += 1
        if (i % 200) == 0 or i == len(primes):
            total = cum_1 + cum_5
            print(f"  n={i}: mod6=1: {cum_1} ({100*cum_1/total:.1f}%), "
                  f"mod6=5: {cum_5} ({100*cum_5/total:.1f}%), bias={cum_1-cum_5}")

def experiment_5_recursive_structure():
    """
    Look for recursive/self-referential structure in the prime sequence.

    Idea: Can p(n) be expressed in terms of p(n-1), p(n-2), ...?
    Or more precisely, is the gap sequence g(n) = p(n+1) - p(n) predictable
    from previous gaps?
    """
    print("\n" + "="*70)
    print("EXPERIMENT 5: RECURSIVE / AUTOREGRESSIVE STRUCTURE")
    print("="*70)

    N = 10000
    primes = [int(prime(n)) for n in range(1, N + 2)]
    gaps = [primes[i+1] - primes[i] for i in range(len(primes)-1)]

    # Fit AR(k) model for various k
    gaps_arr = np.array(gaps, dtype=np.float64)

    for k in [1, 2, 3, 5, 10, 20]:
        # Build regression matrix
        X = np.column_stack([gaps_arr[k-1-i:N-1-i] for i in range(k)])
        y = gaps_arr[k:N]

        coeffs, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
        y_pred = X @ coeffs
        residuals = y - y_pred
        rms = np.sqrt(np.mean(residuals**2))
        r2 = 1 - np.var(residuals) / np.var(y)

        # Test: can we reconstruct primes from gaps?
        # If AR(k) were perfect, we'd have p(n+1) = p(n) + predicted_gap
        exact = np.sum(np.abs(np.round(y_pred) - y) < 0.5)

        print(f"  AR({k:2d}): RMS={rms:.3f}, R²={r2:.4f}, "
              f"exact gap predictions={exact}/{len(y)} ({100*exact/len(y):.1f}%)")

    # What about using the INDEX n as well?
    print("\n  Adding index features:")
    for k in [5, 10]:
        indices = np.arange(k, N, dtype=np.float64)
        ln_indices = np.log(indices)
        X = np.column_stack([
            gaps_arr[k-1-i:N-1-i] for i in range(k)
        ] + [ln_indices[:N-k], 1/ln_indices[:N-k]])
        y = gaps_arr[k:N]

        coeffs, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
        y_pred = X @ coeffs
        residuals = y - y_pred
        rms = np.sqrt(np.mean(residuals**2))
        r2 = 1 - np.var(residuals) / np.var(y)
        exact = np.sum(np.abs(np.round(y_pred) - y) < 0.5)

        print(f"  AR({k:2d})+idx: RMS={rms:.3f}, R²={r2:.4f}, "
              f"exact={exact}/{len(y)} ({100*exact/len(y):.1f}%)")

def main():
    print("Session 6: Lattice / LLL / Integer Relation Approaches")
    print(f"Started: {time.strftime('%H:%M:%S')}\n")

    t0 = time.time()

    experiment_1_integer_relations()
    experiment_2_simultaneous_approximation()
    experiment_3_correction_lattice()
    experiment_4_digit_extraction()
    experiment_5_recursive_structure()

    elapsed = time.time() - t0
    print(f"\n{'='*70}")
    print(f"Total time: {elapsed:.1f}s")

if __name__ == '__main__':
    main()

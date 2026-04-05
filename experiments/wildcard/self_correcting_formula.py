#!/usr/bin/env python3
"""
Self-Correcting Explicit Formula with Integer Constraints

IDEA: Use the explicit formula pi(x) = Li(x) - sum_rho Li(x^rho) + ...
with FEWER zeros than normally needed, then apply integer/monotonicity
constraints and primality testing to correct errors.

Hypothesis: Self-correction reduces the number of zeros needed for exact pi(x).

EXPERIMENT:
- Compute pi_approx(x) using N zeros for x in [2, X_MAX]
- Apply corrections: rounding, monotonicity, step-function, primality
- Measure accuracy vs number of zeros
- Key metric: minimum zeros for 100% accuracy up to x
"""

import sys
import os
import time
import numpy as np
from mpmath import mp, mpf, mpc, li, log, pi as mpi, sqrt, exp, loggamma, inf
from sympy import primepi, isprime, nextprime, prevprime
from collections import defaultdict

# High precision for Li(x^rho) computations
mp.dps = 50

# ---------------------------------------------------------------------------
# Load zeta zeros from data files
# ---------------------------------------------------------------------------
def load_zeros(max_zeros=500):
    """Load zeta zeros from data files, up to max_zeros."""
    zeros = []
    data_dir = "/apps/aplikacijos/prime-research/data"
    for fname in ["zeta_zeros_1000.txt", "zeta_zeros_500.txt",
                   "zeta_zeros_300.txt", "zeta_zeros_200.txt"]:
        fpath = os.path.join(data_dir, fname)
        if os.path.exists(fpath):
            with open(fpath) as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    try:
                        gamma = mpf(line)
                        zeros.append(gamma)
                    except:
                        pass
            if len(zeros) >= max_zeros:
                break
    # Deduplicate and sort
    zeros = sorted(set(zeros))
    return zeros[:max_zeros]


# ---------------------------------------------------------------------------
# Logarithmic integral Li(x) for complex argument via mpmath
# ---------------------------------------------------------------------------
def complex_li(z):
    """Compute Li(z) = li(z) - li(2) for complex z, where li is the
    logarithmic integral. Uses mpmath's li() which handles complex args."""
    try:
        return li(z, offset=True)  # offset=True gives Li(z) = li(z) - li(2)
    except:
        # Fallback: direct integration or series
        return li(z) - li(mpf(2))


# ---------------------------------------------------------------------------
# Explicit formula: pi_approx(x) using N zeros
# ---------------------------------------------------------------------------
def pi_explicit(x, zeros_imag, use_log_term=True):
    """
    Compute pi(x) via truncated explicit formula:
      pi(x) ~ Li(x) - sum_{rho} Li(x^rho) - log(2) + integral term

    Each zero gamma gives a pair of zeros rho = 1/2 + i*gamma and
    rho_bar = 1/2 - i*gamma. Their combined contribution is:
      Li(x^rho) + Li(x^{rho_bar}) = 2 * Re(Li(x^rho))

    Parameters:
        x: evaluation point (real, > 1)
        zeros_imag: list of imaginary parts of zeta zeros
        use_log_term: include the -log(2) + integral correction
    """
    x = mpf(x)
    if x < 2:
        return mpf(0)

    # Main term: Li(x)
    result = li(x, offset=True)

    # Zero contributions: -2 * sum Re(Li(x^rho))
    log_x = log(x)
    for gamma in zeros_imag:
        rho = mpc(mpf('0.5'), gamma)
        # x^rho = exp(rho * log(x))
        x_rho = exp(rho * log_x)
        li_val = complex_li(x_rho)
        result -= 2 * li_val.real

    # Correction terms
    if use_log_term:
        result -= log(mpf(2))
        # The integral from x to infinity of dt/(t*(t^2-1)*log(t)) is small
        # for large x; we approximate it as 0 for x >= 2

    return float(result)


# ---------------------------------------------------------------------------
# Self-correction algorithms
# ---------------------------------------------------------------------------
def round_to_integer(values):
    """Step 1: Round to nearest integer."""
    return [round(v) for v in values]


def enforce_monotonicity(values, xs):
    """Step 2: Enforce pi(x) is non-decreasing."""
    result = list(values)
    for i in range(1, len(result)):
        if result[i] < result[i-1]:
            result[i] = result[i-1]
    return result


def enforce_step_and_primality(values, xs):
    """
    Step 3: Enforce step-function constraint with primality testing.
    pi(x) can only increase at primes, and only by 1.
    For even x > 2, pi(x) = pi(x-1).
    """
    result = list(values)

    # First, identify where jumps occur
    for i in range(1, len(result)):
        x = xs[i]
        diff = result[i] - result[i-1]
        if diff > 0:
            # Jump detected - should only happen at primes
            if not isprime(x):
                # Not a prime, so pi shouldn't jump here
                # Decide: lower current or raise previous?
                # Use the raw approximation to decide
                result[i] = result[i-1]
            elif diff > 1:
                # Jump > 1 at a prime: set to exactly +1
                result[i] = result[i-1] + 1

    # Second pass: check if we missed any primes
    # (where pi should jump but didn't)
    for i in range(1, len(result)):
        x = xs[i]
        if isprime(x) and result[i] == result[i-1]:
            # We know pi should jump here by 1
            # But only if this doesn't conflict with later values
            # Increment this and all subsequent values
            for j in range(i, len(result)):
                result[j] += 1

    return result


def enforce_even_constraint(values, xs):
    """Step 4: For even x > 2, pi(x) = pi(x-1)."""
    result = list(values)
    for i in range(1, len(result)):
        x = xs[i]
        if x > 2 and x % 2 == 0:
            result[i] = result[i-1]
    return result


def full_self_correction(raw_values, xs):
    """Apply all correction steps in sequence."""
    v = round_to_integer(raw_values)
    v = enforce_monotonicity(v, xs)
    v = enforce_even_constraint(v, xs)
    # Note: we do NOT use enforce_step_and_primality by default since
    # it uses isprime which gives exact answers -- we want to measure
    # how much the other constraints help first
    return v


def full_correction_with_primality(raw_values, xs):
    """Apply all corrections including primality testing."""
    v = round_to_integer(raw_values)
    v = enforce_monotonicity(v, xs)
    v = enforce_step_and_primality(v, xs)
    return v


# ---------------------------------------------------------------------------
# Experiment runner
# ---------------------------------------------------------------------------
def compute_exact_pi(xs):
    """Compute exact pi(x) for all x in xs."""
    return [int(primepi(x)) for x in xs]


def accuracy(computed, exact):
    """Fraction of values that match exactly."""
    if len(computed) == 0:
        return 0.0
    matches = sum(1 for c, e in zip(computed, exact) if c == e)
    return matches / len(computed)


def max_error(computed, exact):
    """Maximum absolute error."""
    if len(computed) == 0:
        return 0
    return max(abs(c - e) for c, e in zip(computed, exact))


def mean_error(computed, exact):
    """Mean absolute error."""
    if len(computed) == 0:
        return 0.0
    return sum(abs(c - e) for c, e in zip(computed, exact)) / len(computed)


def run_experiment(x_max, zero_counts, verbose=True):
    """
    Run the self-correcting formula experiment.

    Args:
        x_max: compute pi(x) for x in [2, x_max]
        zero_counts: list of number of zeros to try
        verbose: print progress
    """
    xs = list(range(2, x_max + 1))

    if verbose:
        print(f"Computing exact pi(x) for x = 2 to {x_max}...")
    t0 = time.time()
    exact = compute_exact_pi(xs)
    t_exact = time.time() - t0
    if verbose:
        print(f"  Done in {t_exact:.2f}s. pi({x_max}) = {exact[-1]}")

    # Load all zeros we might need
    max_z = max(zero_counts)
    if verbose:
        print(f"Loading up to {max_z} zeta zeros...")
    all_zeros = load_zeros(max_z)
    if verbose:
        print(f"  Loaded {len(all_zeros)} zeros")

    results = {}

    for n_zeros in zero_counts:
        if n_zeros > len(all_zeros):
            if verbose:
                print(f"\nSkipping N={n_zeros} (only {len(all_zeros)} zeros available)")
            continue

        zeros = all_zeros[:n_zeros]
        if verbose:
            print(f"\n{'='*60}")
            print(f"Testing with N = {n_zeros} zeros")
            print(f"{'='*60}")

        # Compute raw approximations
        t0 = time.time()
        raw = []
        for i, x in enumerate(xs):
            val = pi_explicit(x, zeros)
            raw.append(val)
            if verbose and (i+1) % 500 == 0:
                print(f"  Computed {i+1}/{len(xs)} points...")
        t_compute = time.time() - t0

        # Method 1: Just round
        rounded = round_to_integer(raw)
        acc_round = accuracy(rounded, exact)
        err_round = max_error(rounded, exact)
        mean_err_round = mean_error(rounded, exact)

        # Method 2: Round + monotonicity + even constraint
        corrected = full_self_correction(raw, xs)
        acc_corr = accuracy(corrected, exact)
        err_corr = max_error(corrected, exact)
        mean_err_corr = mean_error(corrected, exact)

        # Method 3: Full correction with primality
        corrected_prime = full_correction_with_primality(raw, xs)
        acc_prime = accuracy(corrected_prime, exact)
        err_prime = max_error(corrected_prime, exact)
        mean_err_prime = mean_error(corrected_prime, exact)

        # Find where errors occur for rounding
        error_positions_round = [(xs[i], rounded[i], exact[i])
                                  for i in range(len(xs))
                                  if rounded[i] != exact[i]]
        error_positions_corr = [(xs[i], corrected[i], exact[i])
                                 for i in range(len(xs))
                                 if corrected[i] != exact[i]]
        error_positions_prime = [(xs[i], corrected_prime[i], exact[i])
                                  for i in range(len(xs))
                                  if corrected_prime[i] != exact[i]]

        results[n_zeros] = {
            'compute_time': t_compute,
            'round': {
                'accuracy': acc_round,
                'max_error': err_round,
                'mean_error': mean_err_round,
                'n_errors': len(error_positions_round),
                'first_errors': error_positions_round[:10],
            },
            'corrected': {
                'accuracy': acc_corr,
                'max_error': err_corr,
                'mean_error': mean_err_corr,
                'n_errors': len(error_positions_corr),
                'first_errors': error_positions_corr[:10],
            },
            'corrected_prime': {
                'accuracy': acc_prime,
                'max_error': err_prime,
                'mean_error': mean_err_prime,
                'n_errors': len(error_positions_prime),
                'first_errors': error_positions_prime[:10],
            },
        }

        if verbose:
            print(f"\n  Compute time: {t_compute:.2f}s")
            print(f"\n  Method 1 (round only):")
            print(f"    Accuracy: {acc_round*100:.2f}%")
            print(f"    Max error: {err_round}")
            print(f"    Mean error: {mean_err_round:.4f}")
            print(f"    Errors: {len(error_positions_round)}/{len(xs)}")
            if error_positions_round:
                print(f"    First error at x={error_positions_round[0][0]}: "
                      f"got {error_positions_round[0][1]}, expected {error_positions_round[0][2]}")

            print(f"\n  Method 2 (round + monotonicity + even):")
            print(f"    Accuracy: {acc_corr*100:.2f}%")
            print(f"    Max error: {err_corr}")
            print(f"    Mean error: {mean_err_corr:.4f}")
            print(f"    Errors: {len(error_positions_corr)}/{len(xs)}")

            print(f"\n  Method 3 (full correction + primality):")
            print(f"    Accuracy: {acc_prime*100:.2f}%")
            print(f"    Max error: {err_prime}")
            print(f"    Mean error: {mean_err_prime:.4f}")
            print(f"    Errors: {len(error_positions_prime)}/{len(xs)}")

    return results, exact, xs


def analyze_error_distribution(xs, raw, exact, n_zeros):
    """Analyze error distribution and correlations."""
    errors = [raw[i] - exact[i] for i in range(len(xs))]
    frac_errors = [raw[i] - round(raw[i]) for i in range(len(xs))]

    mean_e = np.mean(errors)
    std_e = np.std(errors)
    max_abs_e = max(abs(e) for e in errors)

    # Autocorrelation of errors
    errors_np = np.array(errors)
    if len(errors_np) > 10 and np.std(errors_np) > 0:
        autocorr_1 = np.corrcoef(errors_np[:-1], errors_np[1:])[0, 1]
        autocorr_5 = np.corrcoef(errors_np[:-5], errors_np[5:])[0, 1] if len(errors_np) > 10 else 0
    else:
        autocorr_1 = 0
        autocorr_5 = 0

    return {
        'mean_error': mean_e,
        'std_error': std_e,
        'max_abs_error': max_abs_e,
        'autocorr_lag1': autocorr_1,
        'autocorr_lag5': autocorr_5,
        'frac_within_0.5': sum(1 for e in errors if abs(e) < 0.5) / len(errors),
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    print("=" * 70)
    print("SELF-CORRECTING EXPLICIT FORMULA EXPERIMENT")
    print("=" * 70)
    print()

    # Parameters
    X_MAX = 1000  # Range [2, X_MAX]
    ZERO_COUNTS = [5, 10, 20, 30, 50, 100]

    # Run main experiment
    results, exact, xs = run_experiment(X_MAX, ZERO_COUNTS)

    # Additional analysis: error distribution for selected zero counts
    print("\n" + "=" * 70)
    print("ERROR DISTRIBUTION ANALYSIS")
    print("=" * 70)

    all_zeros = load_zeros(max(ZERO_COUNTS))
    for n_z in [10, 30, 50]:
        if n_z > len(all_zeros):
            continue
        zeros = all_zeros[:n_z]
        raw = [pi_explicit(x, zeros) for x in xs]
        dist = analyze_error_distribution(xs, raw, exact, n_z)
        print(f"\n  N={n_z} zeros:")
        print(f"    Mean error: {dist['mean_error']:.4f}")
        print(f"    Std error: {dist['std_error']:.4f}")
        print(f"    Max |error|: {dist['max_abs_error']:.4f}")
        print(f"    Autocorr(1): {dist['autocorr_lag1']:.4f}")
        print(f"    Autocorr(5): {dist['autocorr_lag5']:.4f}")
        print(f"    Fraction within 0.5: {dist['frac_within_0.5']*100:.1f}%")

    # Find minimum zeros for 100% accuracy (with each method)
    print("\n" + "=" * 70)
    print("MINIMUM ZEROS FOR 100% ACCURACY")
    print("=" * 70)

    for method_name in ['round', 'corrected', 'corrected_prime']:
        for n_z in sorted(results.keys()):
            if results[n_z][method_name]['accuracy'] == 1.0:
                print(f"  {method_name}: {n_z} zeros -> 100% accuracy up to x={X_MAX}")
                break
        else:
            best_n = max(results.keys())
            best_acc = results[best_n][method_name]['accuracy']
            print(f"  {method_name}: NOT achieved (best: {best_acc*100:.2f}% with {best_n} zeros)")

    # Extended test at higher range with best zero count
    print("\n" + "=" * 70)
    print("EXTENDED TEST: x up to 2000 with 100 zeros")
    print("=" * 70)
    ext_results, _, _ = run_experiment(2000, [100])

    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY TABLE")
    print("=" * 70)
    print(f"\n{'N zeros':>8} | {'Round Acc':>10} | {'Corr Acc':>10} | {'Prime Acc':>10} | "
          f"{'Round Err':>10} | {'Corr Err':>10} | {'Prime Err':>10} | {'Time':>8}")
    print("-" * 95)
    for n_z in sorted(results.keys()):
        r = results[n_z]
        print(f"{n_z:>8} | {r['round']['accuracy']*100:>9.2f}% | "
              f"{r['corrected']['accuracy']*100:>9.2f}% | "
              f"{r['corrected_prime']['accuracy']*100:>9.2f}% | "
              f"{r['round']['max_error']:>10} | "
              f"{r['corrected']['max_error']:>10} | "
              f"{r['corrected_prime']['max_error']:>10} | "
              f"{r['compute_time']:>7.1f}s")

    print("\nDone.")

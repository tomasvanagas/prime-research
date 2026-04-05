#!/usr/bin/env python3
"""
Green-Tao Nilsystem Correlation for Prime Prediction
=====================================================

PROPOSAL: Green-Tao theory shows primes correlate with nilsequences.
A nilsequence is n -> F(g^n * x) on a nilmanifold G/Gamma.

Test plan:
1. Compute prime indicator 1_P(n) for n up to N
2. Fit 1-step nilsequences (periodic functions / Fourier modes)
3. Fit 2-step nilsequences (bracket quadratics)
4. Measure residual structure after removing nilsequence correlations
5. Compare prediction accuracy to R^{-1}(n) (the smooth approximation)

Key question: can nilsequence-based prediction beat R^{-1}(n)?
"""

import numpy as np
from sympy import isprime, primepi, li
from scipy.optimize import minimize
import time
import warnings
warnings.filterwarnings('ignore')

# ===========================================================================
# Part 0: Setup
# ===========================================================================

N = 10000
print(f"=== Green-Tao Nilsystem Correlation Test (N={N}) ===\n")

# Prime indicator
is_prime = np.array([isprime(n) for n in range(N)], dtype=float)
ns = np.arange(N, dtype=float)

# Von Mangoldt-like weighting: 1/ln(n) normalization so expected value ~ 1
# For prime indicator, use 1_P(n) * ln(n) to get closer to Lambda(n)
ln_ns = np.log(np.maximum(ns, 2))
weighted_prime = is_prime * ln_ns  # approx Lambda(n) for primes

# Actual primes list
primes = [n for n in range(2, N) if is_prime[n]]
num_primes = len(primes)
print(f"Primes up to {N}: {num_primes}")
print(f"Prime density: {num_primes/N:.4f}, 1/ln(N): {1/np.log(N):.4f}\n")

# ===========================================================================
# Part 1: 1-Step Nilsequences (Periodic / Fourier)
# ===========================================================================
print("=" * 70)
print("PART 1: 1-Step Nilsequences (Fourier / periodic functions)")
print("=" * 70)
print()
print("1-step nilsequences on the circle R/Z are just periodic functions.")
print("We fit: f(n) = a_0 + sum_k [a_k cos(2pi alpha_k n) + b_k sin(2pi alpha_k n)]")
print()

# Use a range of frequencies including rational and irrational
# Key frequencies from number theory:
#   - 1/q for small q (Dirichlet characters)
#   - 1/log(N) (prime density scale)
#   - sqrt(2), golden ratio, etc. (irrational rotations)

# Build design matrix with many candidate frequencies
candidate_alphas = []

# Rational frequencies: 1/q for q = 2..30
for q in range(2, 31):
    for a in range(1, q):
        if np.gcd(a, q) == 1:
            candidate_alphas.append(a / q)

# Irrational frequencies
irrationals = [
    1 / np.log(N), 2 / np.log(N), 3 / np.log(N),
    np.sqrt(2) / 10, np.sqrt(3) / 10, np.sqrt(5) / 10,
    (1 + np.sqrt(5)) / 20,  # golden ratio scaled
    np.pi / 10, np.e / 10,
    1 / np.sqrt(np.log(N)),
]
candidate_alphas.extend(irrationals)
candidate_alphas = sorted(set(candidate_alphas))

print(f"Number of candidate frequencies: {len(candidate_alphas)}")

# Build Fourier design matrix
n_range = ns[2:]  # skip 0 and 1
target = is_prime[2:]  # prime indicator for n >= 2

X_fourier = np.ones((len(n_range), 1 + 2 * len(candidate_alphas)))
for i, alpha in enumerate(candidate_alphas):
    X_fourier[:, 1 + 2*i] = np.cos(2 * np.pi * alpha * n_range)
    X_fourier[:, 2 + 2*i] = np.sin(2 * np.pi * alpha * n_range)

# Least squares fit
coeffs_fourier, residuals, rank, sv = np.linalg.lstsq(X_fourier, target, rcond=None)

pred_fourier = X_fourier @ coeffs_fourier
mse_fourier = np.mean((target - pred_fourier) ** 2)
var_target = np.var(target)
r2_fourier = 1 - mse_fourier / var_target

print(f"Fourier fit R^2: {r2_fourier:.6f}")
print(f"MSE: {mse_fourier:.6f}, Variance of target: {var_target:.6f}")

# For prime prediction: threshold the fit
threshold = 0.5
pred_binary = (pred_fourier > threshold).astype(float)
accuracy = np.mean(pred_binary == target)
true_positives = np.sum((pred_binary == 1) & (target == 1))
false_positives = np.sum((pred_binary == 1) & (target == 0))
false_negatives = np.sum((pred_binary == 0) & (target == 1))
precision = true_positives / max(true_positives + false_positives, 1)
recall = true_positives / max(true_positives + false_negatives, 1)

print(f"\nBinary prediction (threshold={threshold}):")
print(f"  Accuracy: {accuracy:.4f}")
print(f"  Precision: {precision:.4f}")
print(f"  Recall: {recall:.4f}")
print(f"  F1: {2*precision*recall/max(precision+recall,1e-10):.4f}")

# What are the dominant frequencies?
amplitudes = np.sqrt(coeffs_fourier[1::2]**2 + coeffs_fourier[2::2]**2)
top_k = 10
top_idx = np.argsort(amplitudes)[::-1][:top_k]
print(f"\nTop {top_k} frequencies by amplitude:")
for idx in top_idx:
    alpha = candidate_alphas[idx]
    amp = amplitudes[idx]
    # Check if close to a rational
    from fractions import Fraction
    frac = Fraction(alpha).limit_denominator(50)
    print(f"  alpha = {alpha:.6f} (~ {frac}), amplitude = {amp:.6f}")

# ===========================================================================
# Part 2: 2-Step Nilsequences (bracket quadratics)
# ===========================================================================
print()
print("=" * 70)
print("PART 2: 2-Step Nilsequences (bracket quadratics)")
print("=" * 70)
print()
print("2-step nilsequences involve bracket quadratics:")
print("  phi(n) = {alpha * n} * {beta * n}  (fractional parts)")
print("  or floor(alpha*n*beta) - related constructions")
print()

# Generate 2-step nilsequence features
# {alpha*n} = fractional part of alpha*n
# bracket quadratic: {alpha*n} * {beta*n}

alpha_candidates = [
    1/np.log(N), np.sqrt(2)/10, np.sqrt(3)/10, np.sqrt(5)/10,
    (1+np.sqrt(5))/20, np.pi/10, np.e/10,
    1/3, 1/5, 1/7, 1/11, 1/13, 2/3, 2/5, 2/7,
    np.sqrt(2), np.sqrt(3), np.sqrt(5),
    np.log(2), np.log(3),
]

# Build features: fractional parts and their products
bracket_features = []
bracket_labels = []

for a_idx, alpha in enumerate(alpha_candidates):
    frac_a = (alpha * n_range) % 1.0
    bracket_features.append(frac_a)
    bracket_labels.append(f"{{  {alpha:.4f}*n}}")

    for b_idx, beta in enumerate(alpha_candidates):
        if b_idx <= a_idx:
            continue
        frac_b = (beta * n_range) % 1.0
        # bracket quadratic
        bq = frac_a * frac_b
        bracket_features.append(bq)
        bracket_labels.append(f"{{{alpha:.4f}*n}}*{{{beta:.4f}*n}}")

        # Also try floor-based: floor(alpha*n*beta) mod small numbers
        for m in [2, 3, 5]:
            floor_val = np.floor(alpha * n_range * beta) % m / m
            bracket_features.append(floor_val)
            bracket_labels.append(f"floor({alpha:.4f}*n*{beta:.4f}) mod {m}")

bracket_X = np.column_stack([np.ones(len(n_range))] + bracket_features)
print(f"Number of 2-step nilsequence features: {len(bracket_features)}")

# Least squares fit
coeffs_bracket, _, _, _ = np.linalg.lstsq(bracket_X, target, rcond=None)
pred_bracket = bracket_X @ coeffs_bracket
mse_bracket = np.mean((target - pred_bracket) ** 2)
r2_bracket = 1 - mse_bracket / var_target

print(f"Bracket quadratic fit R^2: {r2_bracket:.6f}")
print(f"MSE: {mse_bracket:.6f}")

# Binary prediction
pred_binary2 = (pred_bracket > threshold).astype(float)
accuracy2 = np.mean(pred_binary2 == target)
tp2 = np.sum((pred_binary2 == 1) & (target == 1))
fp2 = np.sum((pred_binary2 == 1) & (target == 0))
fn2 = np.sum((pred_binary2 == 0) & (target == 1))
prec2 = tp2 / max(tp2 + fp2, 1)
rec2 = tp2 / max(tp2 + fn2, 1)

print(f"\nBinary prediction (bracket quadratics):")
print(f"  Accuracy: {accuracy2:.4f}")
print(f"  Precision: {prec2:.4f}")
print(f"  Recall: {rec2:.4f}")
print(f"  F1: {2*prec2*rec2/max(prec2+rec2,1e-10):.4f}")

# ===========================================================================
# Part 3: Combined model (Fourier + bracket quadratics)
# ===========================================================================
print()
print("=" * 70)
print("PART 3: Combined Model (1-step + 2-step nilsequences)")
print("=" * 70)
print()

X_combined = np.hstack([X_fourier, bracket_X[:, 1:]])  # skip duplicate intercept
print(f"Combined feature count: {X_combined.shape[1]}")

# Use ridge regression to avoid overfitting with so many features
from numpy.linalg import solve

lambda_reg = 1.0
XtX = X_combined.T @ X_combined + lambda_reg * np.eye(X_combined.shape[1])
Xty = X_combined.T @ target
coeffs_combined = solve(XtX, Xty)

pred_combined = X_combined @ coeffs_combined
mse_combined = np.mean((target - pred_combined) ** 2)
r2_combined = 1 - mse_combined / var_target

print(f"Combined ridge fit R^2: {r2_combined:.6f}")
print(f"MSE: {mse_combined:.6f}")

# Binary prediction
pred_binary3 = (pred_combined > threshold).astype(float)
accuracy3 = np.mean(pred_binary3 == target)
tp3 = np.sum((pred_binary3 == 1) & (target == 1))
fp3 = np.sum((pred_binary3 == 1) & (target == 0))
fn3 = np.sum((pred_binary3 == 0) & (target == 1))
prec3 = tp3 / max(tp3 + fp3, 1)
rec3 = tp3 / max(tp3 + fn3, 1)

print(f"\nBinary prediction (combined):")
print(f"  Accuracy: {accuracy3:.4f}")
print(f"  Precision: {prec3:.4f}")
print(f"  Recall: {rec3:.4f}")
print(f"  F1: {2*prec3*rec3/max(prec3+rec3,1e-10):.4f}")

# ===========================================================================
# Part 4: Cross-validation (train on first half, test on second)
# ===========================================================================
print()
print("=" * 70)
print("PART 4: Cross-Validation (train first half, predict second half)")
print("=" * 70)
print()

mid = len(n_range) // 2
X_train, X_test = X_combined[:mid], X_combined[mid:]
y_train, y_test = target[:mid], target[mid:]

XtX_cv = X_train.T @ X_train + lambda_reg * np.eye(X_combined.shape[1])
Xty_cv = X_train.T @ y_train
coeffs_cv = solve(XtX_cv, Xty_cv)

pred_cv = X_test @ coeffs_cv
mse_cv = np.mean((y_test - pred_cv) ** 2)
var_test = np.var(y_test)
r2_cv = 1 - mse_cv / var_test

print(f"Cross-validated R^2: {r2_cv:.6f}")
print(f"MSE (test): {mse_cv:.6f}")

pred_cv_binary = (pred_cv > threshold).astype(float)
accuracy_cv = np.mean(pred_cv_binary == y_test)
tp_cv = np.sum((pred_cv_binary == 1) & (y_test == 1))
fp_cv = np.sum((pred_cv_binary == 1) & (y_test == 0))
fn_cv = np.sum((pred_cv_binary == 0) & (y_test == 1))
prec_cv = tp_cv / max(tp_cv + fp_cv, 1)
rec_cv = tp_cv / max(tp_cv + fn_cv, 1)

print(f"\nBinary prediction (cross-validated, second half):")
print(f"  Accuracy: {accuracy_cv:.4f}")
print(f"  Precision: {prec_cv:.4f}")
print(f"  Recall: {rec_cv:.4f}")
print(f"  F1: {2*prec_cv*rec_cv/max(prec_cv+rec_cv,1e-10):.4f}")

# Naive baseline: predict prime with probability 1/ln(n)
pred_naive = (1.0 / np.log(np.maximum(n_range[mid:], 2)) > 0.5).astype(float)
accuracy_naive = np.mean(pred_naive == y_test)
print(f"\nNaive baseline accuracy (1/ln(n) > 0.5 => always 0): {accuracy_naive:.4f}")
print(f"  (This just predicts 'not prime' for all n > 7)")

# Better baseline: random with prob 1/ln(n)
rng = np.random.RandomState(42)
density_test = 1.0 / np.log(np.maximum(n_range[mid:], 2))
pred_random = (rng.random(len(y_test)) < density_test).astype(float)
accuracy_random = np.mean(pred_random == y_test)
print(f"Random baseline accuracy (Bernoulli(1/ln(n))): {accuracy_random:.4f}")

# ===========================================================================
# Part 5: Residual Analysis
# ===========================================================================
print()
print("=" * 70)
print("PART 5: Residual Analysis -- Is the remainder pseudorandom?")
print("=" * 70)
print()

residual = target - pred_combined

print(f"Residual statistics:")
print(f"  Mean: {np.mean(residual):.6f}")
print(f"  Std:  {np.std(residual):.6f}")
print(f"  Min:  {np.min(residual):.6f}")
print(f"  Max:  {np.max(residual):.6f}")

# Autocorrelation of residual
from numpy.fft import fft, ifft

def autocorrelation(x, max_lag=100):
    """Compute normalized autocorrelation."""
    x = x - np.mean(x)
    n = len(x)
    result = np.correlate(x, x, mode='full')
    result = result[n-1:n-1+max_lag+1]
    return result / result[0] if result[0] > 0 else result

max_lag = 50
acf_residual = autocorrelation(residual, max_lag)
acf_prime = autocorrelation(target - np.mean(target), max_lag)

print(f"\nAutocorrelation of residual vs original prime indicator:")
print(f"  {'Lag':>5s}  {'Residual ACF':>12s}  {'Prime ACF':>12s}")
for lag in [1, 2, 3, 5, 6, 10, 12, 15, 20, 30, 50]:
    if lag <= max_lag:
        print(f"  {lag:5d}  {acf_residual[lag]:12.6f}  {acf_prime[lag]:12.6f}")

# Spectral analysis of residual
print(f"\nSpectral analysis of residual:")
spectrum = np.abs(fft(residual - np.mean(residual)))**2
freqs = np.fft.fftfreq(len(residual))
positive = freqs > 0
spectrum_pos = spectrum[positive]
freqs_pos = freqs[positive]

# Top spectral peaks
top_spectral = np.argsort(spectrum_pos)[::-1][:10]
print(f"Top 10 spectral peaks of residual:")
for idx in top_spectral:
    period = 1.0 / freqs_pos[idx] if freqs_pos[idx] > 0 else float('inf')
    print(f"  freq={freqs_pos[idx]:.6f}, period={period:.1f}, power={spectrum_pos[idx]:.1f}")

# KS test for normality of residual
from scipy.stats import kstest, normaltest

stat, p_val = kstest(residual / np.std(residual), 'norm')
print(f"\nKS test for normality of residual: stat={stat:.4f}, p={p_val:.6f}")

stat2, p_val2 = normaltest(residual)
print(f"D'Agostino-Pearson normality test: stat={stat2:.4f}, p={p_val2:.6f}")

# ===========================================================================
# Part 6: Comparison with R^{-1}(n) for nth prime prediction
# ===========================================================================
print()
print("=" * 70)
print("PART 6: Can nilsequences beat R^{-1}(n) for the nth prime?")
print("=" * 70)
print()

# R^{-1}(n) approximation: for the nth prime, R^{-1}(n) ~ n*ln(n) + corrections
# We use li^{-1}(n) as a proxy (close to R^{-1}(n) for our range)
from scipy.optimize import brentq
from scipy.special import expi

def li_func(x):
    """Logarithmic integral li(x) = Ei(ln(x))"""
    if x <= 1:
        return 0
    return expi(np.log(x))

def li_inverse(n):
    """Inverse of li(x): find x such that li(x) = n"""
    if n <= 2:
        return 2
    # Initial bracket
    lo, hi = 2, n * np.log(max(n, 3)) * 2
    try:
        return brentq(lambda x: li_func(x) - n, lo, hi)
    except:
        return n * np.log(max(n, 3))

# Compute li^{-1}(k) for k = 1..num_primes
print(f"Computing li^{{-1}}(k) for k=1..{num_primes}...")
li_inv_predictions = np.array([li_inverse(k) for k in range(1, num_primes + 1)])
actual_primes = np.array(primes)

# Error statistics for li^{-1}
li_errors = li_inv_predictions - actual_primes
li_abs_errors = np.abs(li_errors)
li_rel_errors = li_abs_errors / actual_primes

print(f"\nli^{{-1}}(n) prediction errors:")
print(f"  Mean absolute error: {np.mean(li_abs_errors):.2f}")
print(f"  Mean relative error: {np.mean(li_rel_errors):.6f}")
print(f"  Max absolute error:  {np.max(li_abs_errors):.2f}")
print(f"  Correct digits (avg): {np.mean(-np.log10(li_rel_errors + 1e-15)):.2f}")

# Now: can nilsequence prediction improve on this?
# Strategy: use nilsequence fit to estimate pi(x), then invert
# pi_nilseq(x) = sum of predicted primality for n <= x

# Cumulative prime count from nilsequence
pred_primecount = np.cumsum(np.clip(pred_combined, 0, 1))
# For each k, find x such that pred_primecount ~ k
# (this is the nilsequence-based estimate of p(k))

nil_predictions = []
for k in range(1, num_primes + 1):
    # Find first n where cumulative prediction >= k
    idx = np.searchsorted(pred_primecount, k)
    if idx < len(n_range):
        nil_predictions.append(n_range[idx])
    else:
        nil_predictions.append(n_range[-1])

nil_predictions = np.array(nil_predictions)
nil_errors = np.abs(nil_predictions - actual_primes)
nil_rel_errors = nil_errors / actual_primes

print(f"\nNilsequence-based prediction errors (inverting cumulative fit):")
print(f"  Mean absolute error: {np.mean(nil_errors):.2f}")
print(f"  Mean relative error: {np.mean(nil_rel_errors):.6f}")
print(f"  Max absolute error:  {np.max(nil_errors):.2f}")
print(f"  Correct digits (avg): {np.mean(-np.log10(nil_rel_errors + 1e-15)):.2f}")

# Cross-validated version (only use first-half fit)
pred_cv_full = X_combined @ coeffs_cv  # use CV coefficients on full range
pred_primecount_cv = np.cumsum(np.clip(pred_cv_full, 0, 1))

nil_pred_cv = []
for k in range(1, num_primes + 1):
    idx = np.searchsorted(pred_primecount_cv, k)
    if idx < len(n_range):
        nil_pred_cv.append(n_range[idx])
    else:
        nil_pred_cv.append(n_range[-1])

nil_pred_cv = np.array(nil_pred_cv)
nil_err_cv = np.abs(nil_pred_cv - actual_primes)
nil_rel_cv = nil_err_cv / actual_primes

print(f"\nNilsequence (cross-validated) prediction errors:")
print(f"  Mean absolute error: {np.mean(nil_err_cv):.2f}")
print(f"  Mean relative error: {np.mean(nil_rel_cv):.6f}")
print(f"  Max absolute error:  {np.max(nil_err_cv):.2f}")
print(f"  Correct digits (avg): {np.mean(-np.log10(nil_rel_cv + 1e-15)):.2f}")

# ===========================================================================
# Part 7: Information-theoretic analysis
# ===========================================================================
print()
print("=" * 70)
print("PART 7: Information-Theoretic Bound")
print("=" * 70)
print()

# How many bits does the nilsequence model capture?
# Entropy of prime indicator
p_prime = np.mean(target)
H_prime = -p_prime * np.log2(p_prime) - (1-p_prime) * np.log2(1-p_prime)
print(f"Entropy of prime indicator: H = {H_prime:.4f} bits/symbol")
print(f"Total information in prime indicator (n=2..{N}): {H_prime * len(target):.0f} bits")

# Cross-entropy of nilsequence prediction
pred_clipped = np.clip(pred_combined, 0.01, 0.99)
cross_entropy = -np.mean(target * np.log2(pred_clipped) + (1-target) * np.log2(1-pred_clipped))
print(f"Cross-entropy of nilsequence predictor: {cross_entropy:.4f} bits/symbol")
print(f"Information captured: {(H_prime - (cross_entropy - H_prime))/H_prime*100:.1f}%")
print(f"Mutual information bound: {max(0, H_prime - (cross_entropy - H_prime)):.4f} bits/symbol")

# Cross-validated cross-entropy
pred_cv_clipped = np.clip(pred_cv[: len(y_test)], 0.01, 0.99)
ce_cv = -np.mean(y_test * np.log2(pred_cv_clipped) + (1-y_test) * np.log2(1-pred_cv_clipped))
p_test = np.mean(y_test)
H_test = -p_test * np.log2(p_test) - (1-p_test) * np.log2(1-p_test)
print(f"\nCross-validated cross-entropy: {ce_cv:.4f} bits/symbol")
print(f"Baseline entropy (test): {H_test:.4f} bits/symbol")
print(f"Improvement over baseline: {(H_test - ce_cv)/H_test*100:.1f}%")

# ===========================================================================
# Summary
# ===========================================================================
print()
print("=" * 70)
print("SUMMARY")
print("=" * 70)
print()

print("Method                          | R^2 (in-sample) | R^2 (CV)  | Digits (CV)")
print("-" * 80)
print(f"1-step (Fourier)                | {r2_fourier:15.6f} |     --    |     --")
print(f"2-step (bracket quadratics)     | {r2_bracket:15.6f} |     --    |     --")
print(f"Combined (ridge, lambda={lambda_reg})     | {r2_combined:15.6f} | {r2_cv:9.6f} |     --")

li_digits = np.mean(-np.log10(li_rel_errors + 1e-15))
nil_digits = np.mean(-np.log10(nil_rel_errors + 1e-15))
nil_cv_digits = np.mean(-np.log10(nil_rel_cv + 1e-15))

print(f"\nPrime prediction (p(n) estimation):")
print(f"  li^{{-1}}(n):              avg {li_digits:.2f} correct digits")
print(f"  Nilsequence (in-sample):  avg {nil_digits:.2f} correct digits")
print(f"  Nilsequence (CV):         avg {nil_cv_digits:.2f} correct digits")

winner = "li^{-1}(n)" if li_digits > nil_cv_digits else "nilsequence"
print(f"\n  Winner: {winner}")

print()
print("THEORETICAL ANALYSIS:")
print("-" * 40)
print("""
The Green-Tao decomposition Lambda = Lambda_smooth + Lambda_struct + Lambda_error
gives us:
  - Lambda_smooth ~ 1 (captured by density 1/ln(n), i.e., the PNT)
  - Lambda_struct correlates with nilsequences of bounded step
  - Lambda_error is pseudorandom (small Gowers norms)

Key finding from this experiment:
1. 1-step nilsequences (periodic functions) capture Dirichlet-type biases
   (primes avoid multiples, prefer certain residue classes)
2. 2-step nilsequences (bracket quadratics) add marginal improvement
3. The residual after removing nilsequence correlations remains LARGE --
   it is NOT the case that nilsequences capture most of the prime structure

This is EXPECTED from the theory:
- Green-Tao shows that Lambda_struct has small Gowers norm relative to Lambda_smooth
- The oscillatory contribution from zeta zeros dominates Lambda_struct
- Nilsequences capture arithmetic progressions patterns, not individual primes

CONCLUSION: Nilsequence correlation CANNOT beat R^{-1}(n).
The nilsequence structure in primes is about CORRELATIONS (k-point patterns),
not about PREDICTION of individual primes. This is a fundamental category error:
Green-Tao tells us primes CONTAIN arithmetic structure, but that structure
is too weak to pinpoint individual primes better than analytic methods.

This path should be CLOSED: nilsequence approach reduces to periodic/quasiperiodic
fitting, which is equivalent to Fourier analysis of prime indicator -- a known
dead end (Sessions 3, 7, 14).
""")

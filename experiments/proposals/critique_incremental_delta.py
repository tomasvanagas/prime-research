#!/usr/bin/env python3
"""
Critique: Is delta(n) autocorrelation exploitable for incremental prime computation?

Claim: delta(n) = p(n) - R^{-1}(n) has autocorrelation r(1)~0.79 at lag 1,
so consecutive deltas share most information, enabling an incremental algorithm.

Tests:
1. Compute delta(n) for n=1..5000
2. AR(1) and AR(5) prediction of delta(n+1) from past deltas
3. Conditional entropy H(delta(n+1) | delta(n))
4. How prediction error scales with n
"""

import numpy as np
from sympy import prime, mobius
from scipy.special import expi  # for li(x) = Ei(ln(x))
import math
import warnings
warnings.filterwarnings('ignore')

# --- Riemann R function approximation ---
def li(x):
    """Logarithmic integral li(x) = Ei(ln(x))"""
    if x <= 1:
        return 0.0
    return expi(math.log(x))

def R_approx(x, terms=20):
    """Riemann R function: R(x) = sum_{k=1}^{terms} mu(k)/k * li(x^{1/k})"""
    if x <= 1:
        return 0.0
    result = 0.0
    for k in range(1, terms + 1):
        mu_k = int(mobius(k))
        if mu_k == 0:
            continue
        xk = x ** (1.0 / k)
        if xk > 1:
            result += mu_k / k * li(xk)
    return result

# --- Inverse R function via bisection ---
def R_inverse(n, terms=20):
    """Compute R^{-1}(n): the x such that R(x) = n, via bisection."""
    if n <= 0:
        return 0.0
    # Initial bracket: use n*ln(n) as rough estimate
    if n < 3:
        lo, hi = 2.0, 20.0
    else:
        ln_n = math.log(n)
        lo = max(2.0, n * (ln_n - 2))
        hi = n * (ln_n + 2)
    # Widen bracket if needed
    while R_approx(hi, terms) < n:
        hi *= 2
    while R_approx(lo, terms) > n:
        lo /= 2
    # Bisection
    for _ in range(100):
        mid = (lo + hi) / 2
        if R_approx(mid, terms) < n:
            lo = mid
        else:
            hi = mid
        if hi - lo < 1e-10:
            break
    return (lo + hi) / 2

# --- Main computation ---
N = 5000
print(f"Computing delta(n) = p(n) - R^{{-1}}(n) for n=1..{N}")
print("This may take a few minutes...")

primes_list = []
deltas = []
r_inv_values = []

for n in range(1, N + 1):
    p_n = prime(n)
    r_inv_n = R_inverse(n)
    delta_n = p_n - r_inv_n
    primes_list.append(p_n)
    r_inv_values.append(r_inv_n)
    deltas.append(delta_n)
    if n % 1000 == 0:
        print(f"  n={n}: p({n})={p_n}, R^{{-1}}({n})={r_inv_n:.2f}, delta={delta_n:.4f}")

deltas = np.array(deltas)

print(f"\n{'='*70}")
print("BASIC STATISTICS")
print(f"{'='*70}")
print(f"  mean(delta)   = {np.mean(deltas):.4f}")
print(f"  std(delta)    = {np.std(deltas):.4f}")
print(f"  min(delta)    = {np.min(deltas):.4f}")
print(f"  max(delta)    = {np.max(deltas):.4f}")
print(f"  median(delta) = {np.median(deltas):.4f}")

# --- Autocorrelation ---
print(f"\n{'='*70}")
print("AUTOCORRELATION STRUCTURE")
print(f"{'='*70}")

def autocorr(x, lag):
    """Sample autocorrelation at given lag."""
    n = len(x)
    xm = x - np.mean(x)
    c0 = np.sum(xm**2)
    if c0 == 0:
        return 0.0
    ck = np.sum(xm[:n-lag] * xm[lag:])
    return ck / c0

for lag in range(1, 11):
    r = autocorr(deltas, lag)
    print(f"  r({lag:2d}) = {r:.4f}")

r1 = autocorr(deltas, 1)

# --- AR(1) Predictor ---
print(f"\n{'='*70}")
print("AR(1) PREDICTOR: delta_hat(n+1) = mu + r(1)*(delta(n) - mu)")
print(f"{'='*70}")

mu = np.mean(deltas)
# AR(1): predict delta(n+1) = mu + r1*(delta(n) - mu)
ar1_pred = mu + r1 * (deltas[:-1] - mu)
ar1_errors = deltas[1:] - ar1_pred
ar1_rmse = np.sqrt(np.mean(ar1_errors**2))
ar1_mae = np.mean(np.abs(ar1_errors))

print(f"  RMSE          = {ar1_rmse:.4f}")
print(f"  MAE           = {ar1_mae:.4f}")
print(f"  max |error|   = {np.max(np.abs(ar1_errors)):.4f}")
print(f"  Fraction |error| < 0.5:  {np.mean(np.abs(ar1_errors) < 0.5):.4f}")
print(f"  Fraction |error| < 1.0:  {np.mean(np.abs(ar1_errors) < 1.0):.4f}")
print(f"  Fraction |error| < 2.0:  {np.mean(np.abs(ar1_errors) < 2.0):.4f}")

# --- AR(5) Predictor ---
print(f"\n{'='*70}")
print("AR(5) PREDICTOR (OLS fit)")
print(f"{'='*70}")

# Build design matrix for AR(5)
p_order = 5
X_ar = np.zeros((N - p_order, p_order))
for i in range(p_order):
    X_ar[:, i] = deltas[p_order - 1 - i : N - 1 - i]
y_ar = deltas[p_order:]

# OLS fit
XtX_inv = np.linalg.inv(X_ar.T @ X_ar)
beta = XtX_inv @ (X_ar.T @ y_ar)
print(f"  AR(5) coefficients: {beta}")

ar5_pred = X_ar @ beta
ar5_errors = y_ar - ar5_pred
ar5_rmse = np.sqrt(np.mean(ar5_errors**2))
ar5_mae = np.mean(np.abs(ar5_errors))

print(f"  RMSE          = {ar5_rmse:.4f}")
print(f"  MAE           = {ar5_mae:.4f}")
print(f"  max |error|   = {np.max(np.abs(ar5_errors)):.4f}")
print(f"  Fraction |error| < 0.5:  {np.mean(np.abs(ar5_errors) < 0.5):.4f}")
print(f"  Fraction |error| < 1.0:  {np.mean(np.abs(ar5_errors) < 1.0):.4f}")
print(f"  Fraction |error| < 2.0:  {np.mean(np.abs(ar5_errors) < 2.0):.4f}")

# --- Naive predictor (just guess mean) ---
print(f"\n{'='*70}")
print("BASELINE: Naive predictor delta_hat = mean(delta)")
print(f"{'='*70}")
naive_errors = deltas - mu
naive_rmse = np.sqrt(np.mean(naive_errors**2))
print(f"  RMSE          = {naive_rmse:.4f}")
print(f"  Fraction |error| < 0.5:  {np.mean(np.abs(naive_errors) < 0.5):.4f}")

# --- Conditional Entropy Estimation ---
print(f"\n{'='*70}")
print("CONDITIONAL ENTROPY H(delta(n+1) | delta(n))")
print(f"{'='*70}")

# Method 1: Gaussian assumption
# If (delta(n), delta(n+1)) is bivariate Gaussian with correlation r1,
# then H(Y|X) = 0.5 * log2(2*pi*e * sigma^2 * (1 - r1^2))
sigma = np.std(deltas)
var_cond_gaussian = sigma**2 * (1 - r1**2)
H_cond_gaussian = 0.5 * np.log2(2 * np.pi * np.e * var_cond_gaussian)
H_marginal_gaussian = 0.5 * np.log2(2 * np.pi * np.e * sigma**2)

print(f"\n  [Gaussian model]")
print(f"  sigma(delta)          = {sigma:.4f}")
print(f"  r(1)                  = {r1:.4f}")
print(f"  Conditional variance  = {var_cond_gaussian:.4f}")
print(f"  Conditional std dev   = {np.sqrt(var_cond_gaussian):.4f}")
print(f"  H(delta(n))           = {H_marginal_gaussian:.2f} bits")
print(f"  H(delta(n+1)|delta(n))= {H_cond_gaussian:.2f} bits")
print(f"  Mutual info I         = {H_marginal_gaussian - H_cond_gaussian:.2f} bits")
print(f"  Fraction of info shared= {(H_marginal_gaussian - H_cond_gaussian)/H_marginal_gaussian:.2f}")

# Method 2: Discretized empirical estimate
# Round deltas to nearest integer, compute empirical conditional entropy
print(f"\n  [Empirical discrete estimate (rounded to integers)]")
delta_int = np.round(deltas).astype(int)
from collections import Counter

# Joint counts
joint = Counter(zip(delta_int[:-1], delta_int[1:]))
marginal = Counter(delta_int[:-1])
n_pairs = len(delta_int) - 1

H_cond_emp = 0.0
for (x, y), count_xy in joint.items():
    p_xy = count_xy / n_pairs
    p_x = marginal[x] / n_pairs
    p_y_given_x = count_xy / marginal[x]
    if p_y_given_x > 0:
        H_cond_emp -= p_xy * np.log2(p_y_given_x)

H_marg_emp = 0.0
for x, count_x in marginal.items():
    p_x = count_x / n_pairs
    H_marg_emp -= p_x * np.log2(p_x)

print(f"  H(delta_int(n))            = {H_marg_emp:.2f} bits")
print(f"  H(delta_int(n+1)|delta_int(n)) = {H_cond_emp:.2f} bits")
print(f"  Mutual info I              = {H_marg_emp - H_cond_emp:.2f} bits")
print(f"  Unique marginal values     = {len(marginal)}")

# --- Scaling of prediction error with n ---
print(f"\n{'='*70}")
print("SCALING: How does AR(1) prediction error grow with n?")
print(f"{'='*70}")

# Split into windows and measure RMSE in each
window_size = 500
n_windows = (N - 1) // window_size
print(f"  Window size: {window_size}")
print(f"  {'Window':>8s} {'n range':>15s} {'RMSE':>10s} {'MAE':>10s} {'std(delta)':>12s} {'|err|<0.5':>10s}")

for w in range(n_windows):
    start = w * window_size
    end = (w + 1) * window_size
    window_errors = ar1_errors[start:end]
    window_deltas = deltas[start:end+1]
    rmse_w = np.sqrt(np.mean(window_errors**2))
    mae_w = np.mean(np.abs(window_errors))
    std_w = np.std(window_deltas)
    frac_w = np.mean(np.abs(window_errors) < 0.5)
    print(f"  {w+1:8d} {start+1:>6d}-{end:>6d} {rmse_w:10.4f} {mae_w:10.4f} {std_w:12.4f} {frac_w:10.4f}")

# --- Key diagnostic: is the autocorrelation from smooth trends? ---
print(f"\n{'='*70}")
print("DIAGNOSTIC: Is autocorrelation from a smooth trend or genuine memory?")
print(f"{'='*70}")

# Detrend: subtract a moving average, then recheck autocorrelation
from scipy.ndimage import uniform_filter1d
for window in [10, 50, 100, 500]:
    smoothed = uniform_filter1d(deltas, size=window, mode='nearest')
    detrended = deltas - smoothed
    r1_det = autocorr(detrended, 1)
    std_det = np.std(detrended)
    print(f"  After removing MA({window:3d}): r(1) = {r1_det:.4f}, std = {std_det:.4f}")

# --- Final Verdict ---
print(f"\n{'='*70}")
print("VERDICT")
print(f"{'='*70}")

can_round = np.mean(np.abs(ar1_errors) < 0.5)
print(f"""
  Autocorrelation r(1)         = {r1:.4f}
  AR(1) RMSE                   = {ar1_rmse:.4f}
  AR(5) RMSE                   = {ar5_rmse:.4f}

  For exact computation via rounding, we need |prediction error| < 0.5.
  AR(1) achieves this {can_round*100:.1f}% of the time.
  AR(5) achieves this {np.mean(np.abs(ar5_errors) < 0.5)*100:.1f}% of the time.

  Conditional entropy H(delta(n+1)|delta(n)) = {H_cond_gaussian:.2f} bits (Gaussian)
                                             ~= {H_cond_emp:.2f} bits (empirical discrete)

  CONCLUSION:
  - If conditional entropy >> 1 bit: each step requires genuinely new
    information that cannot be derived from previous deltas. Incremental
    approach FAILS.
  - If prediction error >> 0.5: cannot round to exact integer, so
    autocorrelation gives only approximate prediction, not exact.
    Incremental approach FAILS for exact computation.
  - The autocorrelation is {"EXPLOITABLE" if can_round > 0.95 and H_cond_gaussian < 1.0 else "NOT EXPLOITABLE"} for exact prime computation.
""")

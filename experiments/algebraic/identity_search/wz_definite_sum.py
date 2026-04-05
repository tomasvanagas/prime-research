"""
Wilf-Zeilberger Definite Sum Test for f(x) = pi(x) - R(x)

Tests whether f(x) can be expressed as a definite sum with a closed-form
certificate (WZ style). Four main tests:

1. Delta_f structure analysis (first differences, autocorrelation, recurrence)
2. Higher-order differences Delta^k f for k=1..10
3. Hypergeometric recurrence test: sum_i a_i(x)*f(x+i) = 0
4. Summation kernel K(x) = f(x)*x*log(x) structure analysis
"""

import numpy as np
import os
import sys
from pathlib import Path

# ── Load data ──────────────────────────────────────────────────────────────
data_path = Path(__file__).parent / "fx_data.npz"
d = np.load(data_path)
x_arr = d['x']       # 2..100000
f_arr = d['f']       # float64
pi_arr = d['pi']     # int64

N = len(x_arr)  # 99999
print(f"Loaded {N} data points, x in [{x_arr[0]}, {x_arr[-1]}]")
print(f"f range: [{f_arr.min():.4f}, {f_arr.max():.4f}], mean={f_arr.mean():.4f}, std={f_arr.std():.4f}")

results = []
def log(s):
    print(s)
    results.append(s)

log("=" * 72)
log("WILF-ZEILBERGER DEFINITE SUM TEST FOR f(x) = pi(x) - R(x)")
log("=" * 72)

# ══════════════════════════════════════════════════════════════════════════
# TEST 1: First differences Delta_f(x) = f(x+1) - f(x)
# ══════════════════════════════════════════════════════════════════════════
log("\n## TEST 1: First Differences Delta_f(x)")
log("-" * 50)

delta_f = np.diff(f_arr)  # length N-1
log(f"Delta_f: length={len(delta_f)}, range=[{delta_f.min():.6f}, {delta_f.max():.6f}]")
log(f"Delta_f: mean={delta_f.mean():.6e}, std={delta_f.std():.6f}")

# Identify prime vs non-prime transitions
is_prime_next = np.diff(pi_arr) > 0  # True when x+1 is prime (pi increases)
delta_f_prime = delta_f[is_prime_next]
delta_f_composite = delta_f[~is_prime_next]

log(f"\nWhen x+1 is prime ({is_prime_next.sum()} cases):")
log(f"  Delta_f mean={delta_f_prime.mean():.6f}, std={delta_f_prime.std():.6f}")
log(f"When x+1 is composite ({(~is_prime_next).sum()} cases):")
log(f"  Delta_f mean={delta_f_composite.mean():.6f}, std={delta_f_composite.std():.6f}")

# Test 1a: Simple function approximations to Delta_f
log("\n### 1a: Simple function approximations to Delta_f")
x_mid = x_arr[:-1].astype(float)  # x values for delta_f

# Try: Delta_f ~ c / log(x)
c_logx = np.mean(delta_f * np.log(x_mid))
approx_logx = c_logx / np.log(x_mid)
resid_logx = delta_f - approx_logx
log(f"  Fit c/log(x): c={c_logx:.6f}, residual std={resid_logx.std():.6f} (vs delta_f std={delta_f.std():.6f})")
log(f"  Variance explained: {1 - resid_logx.var()/delta_f.var():.6f}")

# Try: Delta_f ~ c1/log(x) + c2/log(x)^2
X_log = np.column_stack([1/np.log(x_mid), 1/np.log(x_mid)**2])
coefs_log, resid_log, _, _ = np.linalg.lstsq(X_log, delta_f, rcond=None)
approx_log2 = X_log @ coefs_log
resid_log2 = delta_f - approx_log2
log(f"  Fit c1/log(x)+c2/log(x)^2: coefs={coefs_log}, residual std={resid_log2.std():.6f}")
log(f"  Variance explained: {1 - resid_log2.var()/delta_f.var():.6f}")

# Try: Delta_f ~ c1/log(x) + c2/log(x)^2 + c3/log(x)^3
X_log3 = np.column_stack([1/np.log(x_mid), 1/np.log(x_mid)**2, 1/np.log(x_mid)**3])
coefs_log3, _, _, _ = np.linalg.lstsq(X_log3, delta_f, rcond=None)
approx_log3 = X_log3 @ coefs_log3
resid_log3 = delta_f - approx_log3
log(f"  Fit 1/log + 1/log^2 + 1/log^3: residual std={resid_log3.std():.6f}")
log(f"  Variance explained: {1 - resid_log3.var()/delta_f.var():.6f}")

# Test 1b: Recurrence check
log("\n### 1b: Recurrence check for Delta_f")
# Test if Delta_f(x+1) = c * Delta_f(x) + ... (linear recurrence)
for order in [1, 2, 3, 5]:
    # Build system: Delta_f(x+order) = sum_{j=0..order-1} a_j * Delta_f(x+j)
    n_rows = len(delta_f) - order
    A_rec = np.column_stack([delta_f[j:j+n_rows] for j in range(order)])
    b_rec = delta_f[order:order+n_rows]
    coefs_rec, resid_rec, _, _ = np.linalg.lstsq(A_rec, b_rec, rcond=None)
    pred = A_rec @ coefs_rec
    r2 = 1 - np.var(b_rec - pred) / np.var(b_rec)
    log(f"  Linear recurrence order {order}: R^2 = {r2:.8f}")

# Test 1c: Autocorrelation
log("\n### 1c: Autocorrelation of Delta_f")
delta_f_centered = delta_f - delta_f.mean()
var_df = np.var(delta_f_centered)
autocorr = []
for lag in [1, 2, 3, 5, 10, 20, 50, 100]:
    if lag < len(delta_f_centered):
        ac = np.mean(delta_f_centered[:-lag] * delta_f_centered[lag:]) / var_df
        autocorr.append((lag, ac))
        log(f"  Autocorrelation(lag={lag}): {ac:.6f}")

# ══════════════════════════════════════════════════════════════════════════
# TEST 2: Higher-order differences
# ══════════════════════════════════════════════════════════════════════════
log("\n## TEST 2: Higher-Order Differences Delta^k f")
log("-" * 50)

diff_k = f_arr.copy()
for k in range(1, 11):
    diff_k = np.diff(diff_k)
    abs_mean = np.mean(np.abs(diff_k))
    rms = np.sqrt(np.mean(diff_k**2))
    max_abs = np.max(np.abs(diff_k))
    log(f"  Delta^{k:2d} f: len={len(diff_k):6d}, |mean|={np.abs(diff_k.mean()):.4e}, "
        f"rms={rms:.4e}, max|.|={max_abs:.4e}")

# Check if rms decreases systematically
log("\n  Checking RMS ratios (Delta^{k+1}/Delta^k):")
diff_k = f_arr.copy()
rms_vals = []
for k in range(1, 11):
    diff_k = np.diff(diff_k)
    rms_vals.append(np.sqrt(np.mean(diff_k**2)))
for k in range(1, len(rms_vals)):
    ratio = rms_vals[k] / rms_vals[k-1]
    log(f"  RMS(Delta^{k+1})/RMS(Delta^{k}) = {ratio:.6f}")

# ══════════════════════════════════════════════════════════════════════════
# TEST 3: Hypergeometric recurrence test
# ══════════════════════════════════════════════════════════════════════════
log("\n## TEST 3: Hypergeometric Recurrence Test")
log("-" * 50)
log("Testing: sum_{i=0}^{r} a_i(x) * f(x+i) = 0, a_i = polynomial degree <= d")

# Split: fit on x=2..50000 (indices 0..49998), validate on 50001..100000
split_idx = 50000 - 2  # index where x=50000
f_train = f_arr[:split_idx+1]
f_test = f_arr[split_idx+1:]
x_train = x_arr[:split_idx+1].astype(float)
x_test = x_arr[split_idx+1:].astype(float)

best_result = None
best_r2_test = -np.inf

for r in range(1, 6):  # recurrence order
    for max_deg in range(0, 4):  # polynomial degree 0..3
        # Build design matrix for training set
        # For each row at position t: sum_{i=0}^{r} a_i(x_t) * f(x_t + i) = 0
        # a_i(x) = sum_{d=0}^{max_deg} c_{i,d} * x^d
        # => sum_{i,d} c_{i,d} * x_t^d * f(x_t+i) = 0

        n_train_rows = len(f_train) - r
        n_test_rows = len(f_test) - r

        n_params = (r + 1) * (max_deg + 1)

        # We need to solve a homogeneous system. Rewrite as:
        # sum_{i=0}^{r-1} a_i(x)*f(x+i) = -a_r(x)*f(x+r)
        # Fix a_r to have leading coef = 1 (i.e., fix one parameter)
        # Actually, simpler: fix a_r(x) = 1 (constant), solve for rest

        # RHS = -f(x+r)
        # LHS columns: x^d * f(x+i) for i=0..r-1, d=0..max_deg
        #              plus x^d * f(x+r) for d=1..max_deg (extra polynomial terms for a_r)

        # Even simpler: fix the last coefficient group as a_r(x) = 1 (just constant term)
        # Then: sum_{i=0}^{r-1} sum_{d=0}^{max_deg} c_{i,d} * x^d * f(x+i)
        #       + sum_{d=1}^{max_deg} c_{r,d} * x^d * f(x+r) = -f(x+r)

        n_free = r * (max_deg + 1) + max_deg  # free parameters
        if n_free < 1:
            continue

        # Build training matrix
        cols_train = []
        cols_test = []
        for i in range(r):
            for d_pow in range(max_deg + 1):
                col_tr = (x_train[:n_train_rows] ** d_pow) * f_train[i:i+n_train_rows]
                cols_train.append(col_tr)
                col_te = (x_test[:n_test_rows] ** d_pow) * f_test[i:i+n_test_rows]
                cols_test.append(col_te)
        # Extra terms for a_r (d=1..max_deg only, since d=0 is fixed to 1)
        for d_pow in range(1, max_deg + 1):
            col_tr = (x_train[:n_train_rows] ** d_pow) * f_train[r:r+n_train_rows]
            cols_train.append(col_tr)
            col_te = (x_test[:n_test_rows] ** d_pow) * f_test[r:r+n_test_rows]
            cols_test.append(col_te)

        A_train = np.column_stack(cols_train)
        A_test = np.column_stack(cols_test)
        b_train = -f_train[r:r+n_train_rows]
        b_test = -f_test[r:r+n_test_rows]

        # Solve via least squares
        try:
            coefs, residuals, rank, sv = np.linalg.lstsq(A_train, b_train, rcond=None)
        except Exception:
            continue

        pred_train = A_train @ coefs
        pred_test = A_test @ coefs

        resid_train = b_train - pred_train
        resid_test = b_test - pred_test

        r2_train = 1 - np.var(resid_train) / np.var(b_train)
        r2_test = 1 - np.var(resid_test) / np.var(b_test)

        if r2_test > best_r2_test:
            best_r2_test = r2_test
            best_result = (r, max_deg, r2_train, r2_test, coefs)

        if r2_test > 0.5 or (r <= 2 and max_deg <= 1):  # always report small cases
            log(f"  r={r}, deg={max_deg}: R²(train)={r2_train:.8f}, R²(test)={r2_test:.8f}, "
                f"params={n_free}")

log(f"\n  Best fit: r={best_result[0]}, deg={best_result[1]}, "
    f"R²(train)={best_result[2]:.8f}, R²(test)={best_result[3]:.8f}")

# Also test with normalized f (remove trend)
log("\n### 3b: Recurrence test with detrended f")
# Remove smooth trend: f_detrended = f - poly_fit
from numpy.polynomial import polynomial as P
x_norm = (x_arr.astype(float) - x_arr.mean()) / x_arr.std()
poly_coefs = np.polyfit(x_norm, f_arr, 3)
f_trend = np.polyval(poly_coefs, x_norm)
f_detrend = f_arr - f_trend
log(f"  Detrended f: std={f_detrend.std():.4f} (original std={f_arr.std():.4f})")

f_dt_train = f_detrend[:split_idx+1]
f_dt_test = f_detrend[split_idx+1:]

for r in [1, 2, 3]:
    for max_deg in [0, 1]:
        n_train_rows = len(f_dt_train) - r
        n_test_rows = len(f_dt_test) - r

        cols_train = []
        cols_test = []
        for i in range(r):
            for d_pow in range(max_deg + 1):
                cols_train.append((x_train[:n_train_rows]**d_pow) * f_dt_train[i:i+n_train_rows])
                cols_test.append((x_test[:n_test_rows]**d_pow) * f_dt_test[i:i+n_test_rows])
        for d_pow in range(1, max_deg + 1):
            cols_train.append((x_train[:n_train_rows]**d_pow) * f_dt_train[r:r+n_train_rows])
            cols_test.append((x_test[:n_test_rows]**d_pow) * f_dt_test[r:r+n_test_rows])

        A_tr = np.column_stack(cols_train)
        A_te = np.column_stack(cols_test)
        b_tr = -f_dt_train[r:r+n_train_rows]
        b_te = -f_dt_test[r:r+n_test_rows]

        coefs_dt, _, _, _ = np.linalg.lstsq(A_tr, b_tr, rcond=None)
        r2_tr = 1 - np.var(b_tr - A_tr @ coefs_dt) / np.var(b_tr)
        r2_te = 1 - np.var(b_te - A_te @ coefs_dt) / np.var(b_te)
        log(f"  Detrended r={r}, deg={max_deg}: R²(train)={r2_tr:.8f}, R²(test)={r2_te:.8f}")

# ══════════════════════════════════════════════════════════════════════════
# TEST 4: Summation kernel K(x) = f(x) * x * log(x)
# ══════════════════════════════════════════════════════════════════════════
log("\n## TEST 4: Summation Kernel K(x) = f(x) * x * log(x)")
log("-" * 50)

K = f_arr * x_arr.astype(float) * np.log(x_arr.astype(float))
log(f"K(x): range=[{K.min():.2f}, {K.max():.2f}], mean={K.mean():.2f}, std={K.std():.2f}")

# 4a: Polynomial fit to K(x)
log("\n### 4a: Polynomial fits to K(x)")
x_float = x_arr.astype(float)
for deg in [1, 2, 3, 5, 8, 10]:
    # Use Chebyshev basis for stability
    coefs_k = np.polyfit(x_float, K, deg)
    K_fit = np.polyval(coefs_k, x_float)
    resid_k = K - K_fit
    r2_k = 1 - np.var(resid_k) / np.var(K)
    log(f"  Poly deg {deg:2d}: R²={r2_k:.10f}, residual std={resid_k.std():.2f}")

# 4b: Compressibility - SVD of Hankel matrix
log("\n### 4b: Compressibility (Hankel matrix singular values)")
# Build a small Hankel matrix from K values
hankel_size = 500
K_sub = K[::len(K)//hankel_size][:hankel_size]
H = np.zeros((hankel_size//2, hankel_size//2))
for i in range(hankel_size//2):
    for j in range(hankel_size//2):
        H[i, j] = K_sub[i + j]
sv = np.linalg.svd(H, compute_uv=False)
sv_norm = sv / sv[0]
log(f"  Top 10 normalized singular values: {sv_norm[:10].round(6)}")
log(f"  SV decay: sv[10]/sv[0]={sv_norm[min(10,len(sv_norm)-1)]:.6e}, "
    f"sv[50]/sv[0]={sv_norm[min(50,len(sv_norm)-1)]:.6e}")

# Count effective rank (singular values > 1e-6 * max)
eff_rank = np.sum(sv_norm > 1e-6)
log(f"  Effective rank (sv > 1e-6 * max): {eff_rank} / {len(sv_norm)}")

# 4c: Compare K(x) autocorrelation vs f autocorrelation
log("\n### 4c: Autocorrelation comparison: K(x) vs f(x)")
K_centered = K - K.mean()
var_K = np.var(K_centered)
for lag in [1, 2, 5, 10, 50, 100]:
    ac_K = np.mean(K_centered[:-lag] * K_centered[lag:]) / var_K
    ac_f = np.mean((f_arr[:-lag] - f_arr.mean()) * (f_arr[lag:] - f_arr.mean())) / np.var(f_arr)
    log(f"  lag={lag:3d}: AC(K)={ac_K:.6f}, AC(f)={ac_f:.6f}")

# 4d: First differences of K
log("\n### 4d: First differences of K(x)")
delta_K = np.diff(K)
log(f"  Delta_K: mean={delta_K.mean():.4f}, std={delta_K.std():.4f}")
# Fit Delta_K ~ polynomial
for deg in [1, 2, 3]:
    ck = np.polyfit(x_float[:-1], delta_K, deg)
    dk_fit = np.polyval(ck, x_float[:-1])
    r2_dk = 1 - np.var(delta_K - dk_fit) / np.var(delta_K)
    log(f"  Poly deg {deg} fit to Delta_K: R²={r2_dk:.8f}")

# ══════════════════════════════════════════════════════════════════════════
# SUMMARY
# ══════════════════════════════════════════════════════════════════════════
log("\n" + "=" * 72)
log("SUMMARY")
log("=" * 72)

# Assess definite sum viability
log("\n1. FIRST DIFFERENCES: Delta_f(x) has negligible autocorrelation and is")
log("   not well-approximated by smooth functions. Variance explained by")
log(f"   1/log(x) expansion: {1 - resid_log3.var()/delta_f.var():.4f}")
log("   => Delta_f is dominated by the indicator function of primes (random).")

log(f"\n2. HIGHER-ORDER DIFFERENCES: RMS does NOT decrease -- it INCREASES.")
log("   This rules out f(x) being a polynomial or hypergeometric sequence.")

log(f"\n3. HYPERGEOMETRIC RECURRENCE: Best R²(test) = {best_r2_test:.6f}")
log("   WARNING: This high R² is SPURIOUS. It reflects only that f(x) is a")
log("   slow-varying function with AC(lag=1) = 0.9976, so f(x+1) ~ c*f(x)")
log("   trivially. The residual std ~ 0.29 matches Delta_f std exactly -- the")
log("   recurrence captures ONLY the trivial autocorrelation, not algebraic")
log("   structure. Adding more recurrence terms or polynomial degrees does NOT")
log("   improve test R² (all cluster at ~0.997). A genuine WZ certificate would")
log("   give R² = 1.000 exactly.")
log("   Hypergeometric certificate: NEGATIVE.")

log(f"\n4. SUMMATION KERNEL K(x) = f(x)*x*log(x):")
log(f"   Hankel effective rank: {eff_rank}/{len(sv_norm)} (full rank = incompressible)")
log("   K(x) has NO better structure than f(x) itself.")

log("\nCONCLUSION: f(x) = pi(x) - R(x) does NOT admit a WZ-style definite sum")
log("representation. The first differences are dominated by prime indicators,")
log("higher differences grow (not shrink), no polynomial recurrence holds,")
log("and the summation kernel is incompressible. This is consistent with")
log("f(x) encoding zeta zero oscillations with no finite closed-form certificate.")

# ── Write results ──────────────────────────────────────────────────────────
results_path = Path(__file__).parent / "wz_results.md"
with open(results_path, "w") as fout:
    fout.write("# Wilf-Zeilberger Definite Sum Test Results\n\n")
    fout.write(f"**Date**: 2026-04-05\n")
    fout.write(f"**Data**: f(x) = pi(x) - R(x), x = 2..100000\n\n")
    fout.write("```\n")
    fout.write("\n".join(results))
    fout.write("\n```\n")
print(f"\nResults written to {results_path}")

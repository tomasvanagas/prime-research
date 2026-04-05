#!/usr/bin/env python3
"""
Test if partial sums of the explicit formula satisfy any recurrence.

For zeta zeros gamma_k, define:
  S_K(x) = sum_{k=1}^{K} cos(gamma_k * log(x)) / sqrt(1/4 + gamma_k^2)

We test:
  1. Linear recurrences of order r=1..20
  2. Nonlinear (polynomial degree 2,3) recurrences
  3. Recurrences on differences d_K = S_K - S_{K-1}
  4. Convergence rate: geometric vs algebraic
  5. Convergence plot |S_K - S_1000| vs K
"""

import numpy as np
import os
import sys

DATA_PATH = "/apps/aplikacijos/prime-research/data/zeta_zeros_1000.txt"
OUT_DIR = "/apps/aplikacijos/prime-research/experiments/analytic/zeta_structure"

# ---- Load zeros ----
gammas = np.loadtxt(DATA_PATH)
N_ZEROS = len(gammas)
print(f"Loaded {N_ZEROS} zeta zeros, range [{gammas[0]:.4f}, {gammas[-1]:.4f}]")

X_VALUES = [100, 1000, 10000, 100000]

def compute_partial_sums(x, gammas):
    """Compute S_1, S_2, ..., S_N where S_K = sum_{k=1}^K cos(gamma_k*log(x))/sqrt(1/4+gamma_k^2)."""
    log_x = np.log(x)
    terms = np.cos(gammas * log_x) / np.sqrt(0.25 + gammas**2)
    return np.cumsum(terms)

# ---- 1. Compute partial sums for each x ----
all_sums = {}
for x in X_VALUES:
    all_sums[x] = compute_partial_sums(x, gammas)
    print(f"x={x:>6d}: S_1000 = {all_sums[x][-1]:.10f}, S_100 = {all_sums[x][99]:.10f}")

# ---- 2. Test linear recurrence of order r ----
def test_linear_recurrence(seq, max_order=20):
    """
    For each order r, fit S_{K+r} = a_1*S_{K+r-1} + ... + a_r*S_K via least squares.
    Return dict: order -> (coefficients, residual_norm, relative_residual).
    """
    n = len(seq)
    results = {}
    for r in range(1, max_order + 1):
        # Build matrix: each row is [S_{K+r-1}, S_{K+r-2}, ..., S_K]
        n_eqs = n - r
        if n_eqs < r + 10:
            break
        A = np.zeros((n_eqs, r))
        b = np.zeros(n_eqs)
        for i in range(n_eqs):
            for j in range(r):
                A[i, j] = seq[i + r - 1 - j]
            b[i] = seq[i + r]
        # Solve least squares
        coeffs, residuals, rank, sv = np.linalg.lstsq(A, b, rcond=None)
        predicted = A @ coeffs
        res_norm = np.linalg.norm(b - predicted)
        b_norm = np.linalg.norm(b)
        rel_res = res_norm / b_norm if b_norm > 1e-15 else float('inf')
        results[r] = (coeffs, res_norm, rel_res)
    return results

print("\n" + "="*80)
print("LINEAR RECURRENCE TEST on S_K")
print("="*80)

linear_results = {}
for x in X_VALUES:
    seq = all_sums[x]
    results = test_linear_recurrence(seq, max_order=20)
    linear_results[x] = results
    print(f"\nx = {x}:")
    print(f"  {'Order':>5s}  {'|residual|':>12s}  {'rel residual':>12s}")
    for r in sorted(results.keys()):
        coeffs, res_norm, rel_res = results[r]
        print(f"  {r:5d}  {res_norm:12.6e}  {rel_res:12.6e}")

# ---- 3. Test nonlinear recurrences ----
print("\n" + "="*80)
print("NONLINEAR RECURRENCE TEST: S_{K+1} = f(S_K, S_{K-1})")
print("="*80)

def test_nonlinear_recurrence(seq, degree=2):
    """
    Fit S_{K+1} = polynomial(S_K, S_{K-1}) of given degree.
    For degree 2: features are [1, S_K, S_{K-1}, S_K^2, S_K*S_{K-1}, S_{K-1}^2]
    For degree 3: add cubic terms.
    """
    n = len(seq)
    n_eqs = n - 2  # need S_{K-1}, S_K, S_{K+1}

    def build_features(s_k, s_km1, deg):
        feats = [1.0, s_k, s_km1]
        if deg >= 2:
            feats += [s_k**2, s_k * s_km1, s_km1**2]
        if deg >= 3:
            feats += [s_k**3, s_k**2 * s_km1, s_k * s_km1**2, s_km1**3]
        return feats

    A_rows = []
    b = []
    for i in range(1, n - 1):
        A_rows.append(build_features(seq[i], seq[i-1], degree))
        b.append(seq[i+1])

    A = np.array(A_rows)
    b = np.array(b)

    coeffs, residuals, rank, sv = np.linalg.lstsq(A, b, rcond=None)
    predicted = A @ coeffs
    res_norm = np.linalg.norm(b - predicted)
    b_norm = np.linalg.norm(b)
    rel_res = res_norm / b_norm if b_norm > 1e-15 else float('inf')
    return coeffs, res_norm, rel_res

nonlinear_results = {}
for x in X_VALUES:
    seq = all_sums[x]
    nonlinear_results[x] = {}
    print(f"\nx = {x}:")
    for deg in [2, 3]:
        coeffs, res_norm, rel_res = test_nonlinear_recurrence(seq, degree=deg)
        nonlinear_results[x][deg] = (coeffs, res_norm, rel_res)
        print(f"  Degree {deg}: |residual| = {res_norm:.6e}, rel = {rel_res:.6e}")

# ---- 4. Test recurrences on DIFFERENCES d_K = S_K - S_{K-1} ----
print("\n" + "="*80)
print("LINEAR RECURRENCE TEST on DIFFERENCES d_K = S_K - S_{K-1}")
print("="*80)

diff_linear_results = {}
for x in X_VALUES:
    seq = all_sums[x]
    diffs = np.diff(seq)  # d_K = S_K - S_{K-1} for K=2..N
    results = test_linear_recurrence(diffs, max_order=20)
    diff_linear_results[x] = results
    print(f"\nx = {x}:")
    print(f"  {'Order':>5s}  {'|residual|':>12s}  {'rel residual':>12s}")
    for r in sorted(results.keys()):
        coeffs, res_norm, rel_res = results[r]
        print(f"  {r:5d}  {res_norm:12.6e}  {rel_res:12.6e}")

# ---- 5. Convergence rate analysis ----
print("\n" + "="*80)
print("CONVERGENCE RATE ANALYSIS")
print("="*80)

convergence_results = {}
for x in X_VALUES:
    seq = all_sums[x]
    S_inf = seq[-1]  # S_1000 as proxy for S_inf
    errors = np.abs(seq[:-1] - S_inf)  # |S_K - S_1000| for K=1..999

    # Avoid zeros
    valid = errors > 1e-16
    K_vals = np.arange(1, N_ZEROS)[valid]
    err_vals = errors[valid]

    # Test geometric: log|error| ~ log(C) + K*log(r)
    # Fit on later portion (K > 100) to avoid transient
    mask_late = K_vals > 100
    if np.sum(mask_late) > 50:
        log_err = np.log(err_vals[mask_late])
        K_late = K_vals[mask_late]
        # Geometric fit: log|err| = a + b*K
        A_geo = np.column_stack([np.ones(len(K_late)), K_late.astype(float)])
        geo_coeffs = np.linalg.lstsq(A_geo, log_err, rcond=None)[0]
        geo_predicted = A_geo @ geo_coeffs
        geo_r2 = 1 - np.sum((log_err - geo_predicted)**2) / np.sum((log_err - np.mean(log_err))**2)
        geo_rate = np.exp(geo_coeffs[1])

        # Algebraic fit: log|err| = a + b*log(K)
        A_alg = np.column_stack([np.ones(len(K_late)), np.log(K_late.astype(float))])
        alg_coeffs = np.linalg.lstsq(A_alg, log_err, rcond=None)[0]
        alg_predicted = A_alg @ alg_coeffs
        alg_r2 = 1 - np.sum((log_err - alg_predicted)**2) / np.sum((log_err - np.mean(log_err))**2)
        alg_alpha = -alg_coeffs[1]
    else:
        geo_rate, geo_r2 = float('nan'), float('nan')
        alg_alpha, alg_r2 = float('nan'), float('nan')

    convergence_results[x] = {
        'geo_rate': geo_rate, 'geo_r2': geo_r2,
        'alg_alpha': alg_alpha, 'alg_r2': alg_r2,
    }
    print(f"\nx = {x}:")
    print(f"  Geometric: rate = {geo_rate:.6f}, R^2 = {geo_r2:.6f}")
    print(f"  Algebraic: alpha = {alg_alpha:.6f}, R^2 = {alg_r2:.6f}")
    better = "Algebraic" if alg_r2 > geo_r2 else "Geometric"
    print(f"  Better fit: {better}")

# ---- 6. Plot |S_K - S_1000| vs K ----
print("\n" + "="*80)
print("GENERATING CONVERGENCE PLOT")
print("="*80)

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    for idx, x in enumerate(X_VALUES):
        ax = axes[idx // 2][idx % 2]
        seq = all_sums[x]
        S_inf = seq[-1]
        errors = np.abs(seq[:-1] - S_inf)
        K_vals = np.arange(1, N_ZEROS)

        # Log-log plot
        valid = errors > 1e-16
        ax.semilogy(K_vals[valid], errors[valid], '.', markersize=1, alpha=0.5, label='|S_K - S_1000|')

        # Overlay algebraic fit
        cr = convergence_results[x]
        if not np.isnan(cr['alg_alpha']):
            K_fit = np.arange(100, 999)
            err_fit = np.exp(np.log(K_fit) * (-cr['alg_alpha'])) * np.exp(
                np.log(errors[valid][-1]) + cr['alg_alpha'] * np.log(K_vals[valid][-1])
            )
            # Simpler: just show the trend
            ax.set_title(f'x={x}, alpha={cr["alg_alpha"]:.2f}, geo_r={cr["geo_rate"]:.4f}')

        ax.set_xlabel('K (number of zeros)')
        ax.set_ylabel('|S_K - S_1000|')
        ax.grid(True, alpha=0.3)

    plt.suptitle('Partial Sum Convergence: |S_K - S_1000| vs K', fontsize=14)
    plt.tight_layout()
    plot_path = os.path.join(OUT_DIR, 'partial_sums_convergence.png')
    plt.savefig(plot_path, dpi=150)
    print(f"Plot saved to {plot_path}")
except ImportError:
    print("matplotlib not available, skipping plot")

# ---- 7. Additional: check if the sequence of individual terms has structure ----
print("\n" + "="*80)
print("INDIVIDUAL TERM ANALYSIS")
print("="*80)

for x in X_VALUES:
    log_x = np.log(x)
    terms = np.cos(gammas * log_x) / np.sqrt(0.25 + gammas**2)
    print(f"\nx = {x}:")
    print(f"  Mean term magnitude: {np.mean(np.abs(terms)):.6e}")
    print(f"  Std of terms: {np.std(terms):.6e}")
    print(f"  Max |term|: {np.max(np.abs(terms)):.6e}")
    print(f"  Autocorrelation lag-1: {np.corrcoef(terms[:-1], terms[1:])[0,1]:.6f}")
    # Check if terms decay
    # The envelope should decay like 1/gamma_k ~ 1/k (by Weyl law gamma_k ~ 2*pi*k/log(k))
    magnitudes = np.abs(terms)
    # Fit log(|term|) vs log(k)
    K_arr = np.arange(1, N_ZEROS + 1).astype(float)
    valid = magnitudes > 1e-20
    if np.sum(valid) > 100:
        A_dec = np.column_stack([np.ones(np.sum(valid)), np.log(K_arr[valid])])
        dec_coeffs = np.linalg.lstsq(A_dec, np.log(magnitudes[valid]), rcond=None)[0]
        print(f"  Envelope decay: |term_k| ~ k^{dec_coeffs[1]:.3f}")

# ---- 8. Summary ----
print("\n" + "="*80)
print("SUMMARY")
print("="*80)

print("\nLinear recurrence: checking if minimum relative residual < 0.01 for any order...")
for x in X_VALUES:
    min_rel = min(r[2] for r in linear_results[x].values())
    best_order = min(linear_results[x].keys(), key=lambda r: linear_results[x][r][2])
    print(f"  x={x:>6d}: best order={best_order}, rel_residual={min_rel:.6e}")
    if min_rel < 0.01:
        print(f"    *** POSSIBLE RECURRENCE DETECTED ***")
    else:
        print(f"    No recurrence (residual too large)")

print("\nNonlinear recurrence:")
for x in X_VALUES:
    for deg in [2, 3]:
        _, _, rel = nonlinear_results[x][deg]
        print(f"  x={x:>6d}, deg={deg}: rel_residual={rel:.6e}", end="")
        print("  *** POSSIBLE ***" if rel < 0.01 else "  No")

print("\nDifference recurrence:")
for x in X_VALUES:
    min_rel = min(r[2] for r in diff_linear_results[x].values())
    best_order = min(diff_linear_results[x].keys(), key=lambda r: diff_linear_results[x][r][2])
    print(f"  x={x:>6d}: best order={best_order}, rel_residual={min_rel:.6e}")

print("\nConvergence type:")
for x in X_VALUES:
    cr = convergence_results[x]
    better = "Algebraic" if cr['alg_r2'] > cr['geo_r2'] else "Geometric"
    print(f"  x={x:>6d}: {better} (alg R2={cr['alg_r2']:.4f}, geo R2={cr['geo_r2']:.4f})")

print("\nDone.")

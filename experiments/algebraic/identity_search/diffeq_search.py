#!/usr/bin/env python3
"""
Differential equation search for f(x) = pi(x) - R(x).

Tests whether the smoothed remainder f(x) satisfies any ODE or integral equation
in a regularized sense.

Sections:
1. Smoothed f(x) with Gaussian kernels
2. Linear ODE search via SVD
3. Multiplicative (Euler-type) ODE search
4. Non-linear ODE search
5. Volterra integral equation test
"""

import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.integrate import cumulative_trapezoid
import time, sys, os

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
DATA = np.load(os.path.join(os.path.dirname(__file__), "fx_data.npz"))
x_all = DATA["x"].astype(float)       # 2..100000
f_all = DATA["f"].astype(float)
pi_all = DATA["pi"].astype(float)

N = len(x_all)
dx = 1.0  # uniform spacing

results = []  # collect (section, text) for report

def report(section, text):
    results.append((section, text))
    print(f"[{section}] {text}")

# =========================================================================
# 1. Smoothed f(x) and derivatives
# =========================================================================
report("1", "Computing smoothed f(x) and derivatives...")

sigmas = [5, 10, 20, 50]

# High-order finite difference stencils (8th order central)
# Coefficients for 1st derivative (9-point stencil)
c1 = np.array([1/280, -4/105, 1/5, -4/5, 0, 4/5, -1/5, 4/105, -1/280])
# Coefficients for 2nd derivative (9-point stencil)
c2 = np.array([-1/560, 8/315, -1/5, 8/5, -205/72, 8/5, -1/5, 8/315, -1/560])
# Coefficients for 3rd derivative (7-point stencil)
c3 = np.array([-1/8, 1, -13/8, 0, 13/8, -1, 1/8])

def compute_derivatives(fs):
    """Compute f', f'', f''' using high-order stencils. Returns trimmed arrays."""
    n = len(fs)
    # Pad to avoid edge effects
    pad = 50
    # 1st derivative
    fp = np.convolve(fs, c1[::-1], mode='same') / dx
    # 2nd derivative
    fpp = np.convolve(fs, c2[::-1], mode='same') / dx**2
    # 3rd derivative
    fppp = np.convolve(fs, c3[::-1], mode='same') / dx**3
    # Trim edges where stencils are unreliable
    trim = 30
    return fp[trim:-trim], fpp[trim:-trim], fppp[trim:-trim]

smoothed = {}
for sigma in sigmas:
    fs = gaussian_filter1d(f_all, sigma=sigma)
    fp, fpp, fppp = compute_derivatives(fs)
    trim = 30
    xs = x_all[trim:-trim]
    fs_trim = fs[trim:-trim]
    smoothed[sigma] = {
        'x': xs, 'f': fs_trim, 'fp': fp, 'fpp': fpp, 'fppp': fppp
    }
    report("1", f"  sigma={sigma}: f range [{fs_trim.min():.3f}, {fs_trim.max():.3f}], "
           f"|f'| max={np.max(np.abs(fp)):.6f}, |f''| max={np.max(np.abs(fpp)):.6f}")

# =========================================================================
# 2. Linear ODE search via SVD
# =========================================================================
report("2", "Linear ODE search: sum P_j(x) * f^(j)(x) = 0")

def ode_search(xs, derivs, r_max, d_max, label):
    """
    derivs: list [f, f', f'', f'''] trimmed to same length as xs.
    r_max: max ODE order to test
    d_max: max polynomial degree in coefficients
    Returns list of (r, d, sigma_min, sigma_ratio, condition) tuples.
    """
    # Subsample for speed
    step = max(1, len(xs) // 2000)
    idx = np.arange(0, len(xs), step)
    xsub = xs[idx]
    dsub = [dd[idx] for dd in derivs]
    n_pts = len(xsub)

    # Normalize x to [0,1] for numerical stability
    x_norm = (xsub - xsub[0]) / (xsub[-1] - xsub[0])

    results_ode = []

    for r in range(1, r_max + 1):
        for d in range(0, d_max + 1):
            # Number of unknowns: (r+1) coefficients per polynomial, each degree d+1
            n_unknowns = (r + 1) * (d + 1)
            if n_unknowns > n_pts:
                continue

            # Build matrix A: each row is one x-value
            # Columns: for j=0..r, for k=0..d: x^k * f^(j)(x)
            A = np.zeros((n_pts, n_unknowns))
            col = 0
            for j in range(r + 1):
                for k in range(d + 1):
                    A[:, col] = (x_norm ** k) * dsub[j]
                    col += 1

            # SVD
            U, S, Vt = np.linalg.svd(A, full_matrices=False)
            sigma_min = S[-1]
            sigma_max = S[0]
            ratio = sigma_min / sigma_max if sigma_max > 0 else float('inf')

            # Also compute relative residual of best-fit null vector
            null_vec = Vt[-1]
            residual = np.linalg.norm(A @ null_vec) / np.linalg.norm(null_vec)
            rel_residual = residual / np.sqrt(n_pts)

            results_ode.append((r, d, sigma_min, ratio, rel_residual, n_unknowns))

    return results_ode

ode_results_all = {}
for sigma in sigmas:
    sd = smoothed[sigma]
    derivs = [sd['f'], sd['fp'], sd['fpp'], sd['fppp']]
    res = ode_search(sd['x'], derivs, r_max=3, d_max=3, label=f"sigma={sigma}")
    ode_results_all[sigma] = res

    report("2", f"\n  sigma={sigma}:")
    report("2", f"  {'r':>2} {'d':>2} {'sigma_min':>12} {'ratio':>12} {'rel_resid':>12} {'n_unk':>5}")
    for r, d, smin, ratio, rr, nunk in sorted(res, key=lambda t: t[4]):
        report("2", f"  {r:2d} {d:2d} {smin:12.4e} {ratio:12.4e} {rr:12.4e} {nunk:5d}")

# =========================================================================
# 3. Multiplicative (Euler-type) ODE
# =========================================================================
report("3", "Euler-type ODE tests")

euler_results = {}
for sigma in sigmas:
    sd = smoothed[sigma]
    xs, f, fp, fpp = sd['x'], sd['f'], sd['fp'], sd['fpp']
    n = len(xs)
    step = max(1, n // 3000)
    idx = np.arange(0, n, step)

    # Test 1: x*f'(x) + a*f(x) + b = 0
    # => [x*f', f, 1] @ [1, a, b] = 0
    A1 = np.column_stack([xs[idx] * fp[idx], f[idx], np.ones(len(idx))])
    U1, S1, Vt1 = np.linalg.svd(A1, full_matrices=False)
    null1 = Vt1[-1]
    res1 = np.linalg.norm(A1 @ null1) / (np.linalg.norm(null1) * np.sqrt(len(idx)))
    # Extract coefficients (normalize so coeff of x*f' = 1)
    if abs(null1[0]) > 1e-12:
        a_val = null1[1] / null1[0]
        b_val = null1[2] / null1[0]
    else:
        a_val = b_val = float('nan')

    # Test 2: x^2*f''(x) + a*x*f'(x) + b*f(x) + c = 0
    A2 = np.column_stack([xs[idx]**2 * fpp[idx], xs[idx] * fp[idx], f[idx], np.ones(len(idx))])
    U2, S2, Vt2 = np.linalg.svd(A2, full_matrices=False)
    null2 = Vt2[-1]
    res2 = np.linalg.norm(A2 @ null2) / (np.linalg.norm(null2) * np.sqrt(len(idx)))
    if abs(null2[0]) > 1e-12:
        a2 = null2[1] / null2[0]
        b2 = null2[2] / null2[0]
        c2 = null2[3] / null2[0]
    else:
        a2 = b2 = c2 = float('nan')

    euler_results[sigma] = {
        'order1': {'a': a_val, 'b': b_val, 'residual': res1,
                   'sigma_ratio': S1[-1]/S1[0]},
        'order2': {'a': a2, 'b': b2, 'c': c2, 'residual': res2,
                   'sigma_ratio': S2[-1]/S2[0]},
    }

    report("3", f"\n  sigma={sigma}:")
    report("3", f"    1st order Euler: x*f' + {a_val:.4f}*f + {b_val:.4f} = 0  "
           f"residual={res1:.4e}, sigma_ratio={S1[-1]/S1[0]:.4e}")
    report("3", f"    2nd order Euler: x^2*f'' + {a2:.4f}*x*f' + {b2:.4f}*f + {c2:.4f} = 0  "
           f"residual={res2:.4e}, sigma_ratio={S2[-1]/S2[0]:.4e}")

# =========================================================================
# 4. Non-linear ODE search: f'(x) = F(x, f(x))
# =========================================================================
report("4", "Non-linear ODE: f'(x) = F(x, f(x)), polynomial F")

nonlinear_results = {}
for sigma in sigmas:
    sd = smoothed[sigma]
    xs, f, fp = sd['x'], sd['f'], sd['fp']
    n = len(xs)
    step = max(1, n // 3000)
    idx = np.arange(0, n, step)

    x_sub = xs[idx]
    f_sub = f[idx]
    fp_sub = fp[idx]

    # Normalize
    x_norm = (x_sub - x_sub.min()) / (x_sub.max() - x_sub.min())
    f_norm = f_sub / (np.std(f_sub) + 1e-12)
    fp_norm = fp_sub

    nl_res = {}
    for deg_f in [1, 2, 3]:
        for deg_x in [0, 1, 2]:
            # Build design matrix: columns are x^i * f^j for i=0..deg_x, j=0..deg_f
            cols = []
            for i in range(deg_x + 1):
                for j in range(deg_f + 1):
                    cols.append((x_norm ** i) * (f_norm ** j))
            A = np.column_stack(cols)

            # Solve A @ c = fp_sub (least squares)
            c, residuals, rank, sv = np.linalg.lstsq(A, fp_sub, rcond=None)
            pred = A @ c
            rmse = np.sqrt(np.mean((pred - fp_sub)**2))
            rel_rmse = rmse / (np.std(fp_sub) + 1e-12)

            nl_res[(deg_f, deg_x)] = {'rmse': rmse, 'rel_rmse': rel_rmse}

    nonlinear_results[sigma] = nl_res

    report("4", f"\n  sigma={sigma}:")
    report("4", f"  {'deg_f':>5} {'deg_x':>5} {'RMSE':>12} {'rel_RMSE':>12}")
    for (df, dx), v in sorted(nl_res.items(), key=lambda t: t[1]['rel_rmse']):
        report("4", f"  {df:5d} {dx:5d} {v['rmse']:12.4e} {v['rel_rmse']:12.4e}")

# =========================================================================
# 5. Volterra integral equation test
# =========================================================================
report("5", "Volterra integral equation: f(x) = int_2^x K(x,t)*f(t)dt + g(x)")

volterra_results = {}
for sigma in sigmas:
    sd = smoothed[sigma]
    xs, f = sd['x'], sd['f']
    n = len(xs)

    # Subsample heavily for O(n^2) integrals
    step = max(1, n // 500)
    idx = np.arange(0, n, step)
    x_sub = xs[idx]
    f_sub = f[idx]
    m = len(x_sub)

    volt_res = {}

    # For each kernel K, compute I(x) = int_2^x K(x,t)*f(t)dt at subsampled points
    # Then fit: f(x) = alpha * I(x) + g(x), where g(x) is polynomial degree 0,1,2

    kernels = {
        '1/log(t)': lambda x, t: 1.0 / np.maximum(np.log(t), 0.7),
        '1/(x-t+1)': lambda x, t: 1.0 / (x - t + 1),
        '1/t': lambda x, t: 1.0 / t,
    }

    for kname, kfunc in kernels.items():
        # Compute integral at each point via trapezoidal rule
        I_vals = np.zeros(m)
        for i in range(1, m):
            t_pts = x_sub[:i+1]
            integrand = kfunc(x_sub[i], t_pts) * f_sub[:i+1]
            I_vals[i] = np.trapezoid(integrand, t_pts)

        # Fit: f = alpha*I + b0 + b1*(x-x0) + b2*(x-x0)^2
        x_c = (x_sub - x_sub[0]) / (x_sub[-1] - x_sub[0])
        for g_deg in [0, 1, 2]:
            cols = [I_vals]
            for k in range(g_deg + 1):
                cols.append(x_c ** k)
            A = np.column_stack(cols)

            # Skip first point (integral = 0 trivially)
            A_fit = A[1:]
            f_fit = f_sub[1:]

            c, _, _, _ = np.linalg.lstsq(A_fit, f_fit, rcond=None)
            pred = A_fit @ c
            rmse = np.sqrt(np.mean((pred - f_fit)**2))
            rel_rmse = rmse / (np.std(f_fit) + 1e-12)

            volt_res[(kname, g_deg)] = {
                'rmse': rmse, 'rel_rmse': rel_rmse, 'alpha': c[0]
            }

    volterra_results[sigma] = volt_res

    report("5", f"\n  sigma={sigma}:")
    report("5", f"  {'kernel':<15} {'g_deg':>5} {'RMSE':>12} {'rel_RMSE':>12} {'alpha':>12}")
    for (kn, gd), v in sorted(volt_res.items(), key=lambda t: t[1]['rel_rmse']):
        report("5", f"  {kn:<15} {gd:5d} {v['rmse']:12.4e} {v['rel_rmse']:12.4e} {v['alpha']:12.4e}")

# =========================================================================
# Summary: check if anything is close to zero
# =========================================================================
report("SUMMARY", "\n" + "="*70)
report("SUMMARY", "OVERALL ASSESSMENT")
report("SUMMARY", "="*70)

# Best linear ODE
best_linear = None
for sigma, res_list in ode_results_all.items():
    for r, d, smin, ratio, rr, nunk in res_list:
        if best_linear is None or rr < best_linear[2]:
            best_linear = (sigma, (r, d), rr, ratio)

report("SUMMARY", f"\nBest linear ODE: sigma={best_linear[0]}, order={best_linear[1][0]}, "
       f"deg={best_linear[1][1]}, rel_residual={best_linear[2]:.4e}, "
       f"sigma_ratio={best_linear[3]:.4e}")

# Best Euler
best_euler = None
for sigma, er in euler_results.items():
    for key, v in er.items():
        if best_euler is None or v['residual'] < best_euler[2]:
            best_euler = (sigma, key, v['residual'], v['sigma_ratio'])

report("SUMMARY", f"Best Euler ODE: sigma={best_euler[0]}, {best_euler[1]}, "
       f"residual={best_euler[2]:.4e}, sigma_ratio={best_euler[3]:.4e}")

# Best non-linear
best_nl = None
for sigma, nr in nonlinear_results.items():
    for (df, dx), v in nr.items():
        if best_nl is None or v['rel_rmse'] < best_nl[3]:
            best_nl = (sigma, df, dx, v['rel_rmse'])

report("SUMMARY", f"Best non-linear: sigma={best_nl[0]}, deg_f={best_nl[1]}, "
       f"deg_x={best_nl[2]}, rel_RMSE={best_nl[3]:.4e}")

# Best Volterra
best_volt = None
for sigma, vr in volterra_results.items():
    for (kn, gd), v in vr.items():
        if best_volt is None or v['rel_rmse'] < best_volt[3]:
            best_volt = (sigma, kn, gd, v['rel_rmse'])

report("SUMMARY", f"Best Volterra: sigma={best_volt[0]}, kernel={best_volt[1]}, "
       f"g_deg={best_volt[2]}, rel_RMSE={best_volt[3]:.4e}")

# Threshold for "interesting"
THRESHOLD = 0.01  # 1% relative residual
report("SUMMARY", f"\nThreshold for 'near-exact' DE: rel_residual < {THRESHOLD}")

any_found = False
if best_linear[2] < THRESHOLD:
    report("SUMMARY", "*** LINEAR ODE: CANDIDATE FOUND ***")
    any_found = True
if best_euler[2] < THRESHOLD:
    report("SUMMARY", "*** EULER ODE: CANDIDATE FOUND ***")
    any_found = True
if best_nl[3] < THRESHOLD:
    report("SUMMARY", "*** NON-LINEAR ODE: CANDIDATE FOUND ***")
    any_found = True
if best_volt[3] < THRESHOLD:
    report("SUMMARY", "*** VOLTERRA: CANDIDATE FOUND ***")
    any_found = True

if not any_found:
    report("SUMMARY", "NO differential equation found with residual < 1%.")
    report("SUMMARY", "f(x) = pi(x) - R(x) does NOT satisfy any simple ODE or integral equation.")
    report("SUMMARY", "This is consistent with f(x) encoding oscillatory zeta-zero contributions")
    report("SUMMARY", "that are not annihilated by any fixed-order linear/nonlinear DE.")

# =========================================================================
# Write results to markdown
# =========================================================================
md_path = os.path.join(os.path.dirname(__file__), "diffeq_results.md")
with open(md_path, "w") as fout:
    fout.write("# Differential Equation Search for f(x) = pi(x) - R(x)\n\n")
    fout.write("**Date:** 2026-04-05\n\n")
    fout.write("**Data:** x in [2, 100000], f(x) = pi(x) - R(x)\n\n")
    fout.write("**Method:** Gaussian smoothing + SVD-based ODE search + nonlinear fits + Volterra\n\n")

    current_section = None
    for sec, text in results:
        if sec != current_section:
            fout.write(f"\n## Section {sec}\n\n")
            current_section = sec
        fout.write(text + "\n")

    fout.write("\n---\n*Generated by diffeq_search.py*\n")

print(f"\nResults written to {md_path}")

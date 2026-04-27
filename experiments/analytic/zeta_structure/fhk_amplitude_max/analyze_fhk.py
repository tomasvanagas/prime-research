"""Post-process fhk_amplitude_max_results.json:

Produces additional diagnostics beyond the main script:
1. Per-T Anderson-Darling-like Gumbel goodness of fit.
2. FHK-vs-Selberg head-to-head model comparison via log-likelihood
   (treating each run as a sample from the respective model).
3. Bootstrap CIs on M_T mean and var.
4. Three-T linear-regression: <max> = a + b * loglog(T) + c * logloglog(T)
5. Histogram-style summary of M_T over all T pooled.
"""

import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def load(path):
    with open(path) as f:
        return json.load(f)


def gumbel_logpdf(x, loc, scale):
    z = (np.asarray(x) - loc) / scale
    return -z - np.exp(-z) - math.log(scale)


def normal_logpdf(x, loc, scale):
    z = (np.asarray(x) - loc) / scale
    return -0.5 * z * z - 0.5 * math.log(2 * math.pi) - math.log(scale)


def bootstrap_ci(values, stat_fn, B=2000, alpha=0.05, seed=42):
    rng = np.random.default_rng(seed)
    arr = np.asarray(values)
    n = len(arr)
    samples = np.empty(B)
    for b in range(B):
        idx = rng.integers(0, n, n)
        samples[b] = stat_fn(arr[idx])
    lo, hi = np.quantile(samples, [alpha / 2, 1 - alpha / 2])
    return float(lo), float(hi)


def main():
    json_path = os.path.join(HERE, "fhk_amplitude_max_results.json")
    if not os.path.exists(json_path):
        print(f"ERROR: {json_path} not found")
        sys.exit(1)

    out = []
    out.append("=" * 70)
    out.append("FHK |zeta| max statistics — post-hoc analysis")
    out.append("=" * 70)

    data = load(json_path)
    summary = data["summary"]
    raw = data["raw_per_T"]
    params = data["params"]
    out.append(f"K = {params['K']} windows per T, M = {params['M']} samples per window\n")

    # Per-T Gumbel vs Gauss fit
    pooled_M_T = []
    pooled_T = []
    for T_str in sorted(summary.keys(), key=lambda s: float(s)):
        T_base = float(T_str)
        log_T = math.log(T_base)
        loglog_T = math.log(log_T)
        logloglog_T = math.log(loglog_T)
        records = raw[T_str]
        maxes = np.array([r["max"] for r in records])
        M_T = maxes - loglog_T + 0.75 * logloglog_T
        pooled_M_T.extend(M_T.tolist())
        pooled_T.extend([T_base] * len(M_T))

        # Bootstrap CI on M_T mean and var
        mean_lo, mean_hi = bootstrap_ci(M_T, np.mean, seed=42)
        var_lo, var_hi = bootstrap_ci(M_T, lambda x: np.var(x, ddof=1), seed=43)

        # Goodness of fit: free Gumbel MLE vs free Gaussian MLE
        m = float(np.mean(M_T))
        s = float(np.std(M_T, ddof=1))
        # Method-of-moments Gumbel:
        gumbel_scale_mom = s * math.sqrt(6.0) / math.pi
        gumbel_loc_mom = m - 0.5772156649015329 * gumbel_scale_mom
        # Fixed-scale FHK Gumbel(loc, 1/2):
        gumbel_loc_fix = m - 0.5772156649015329 * 0.5
        # Gaussian:
        gauss_loc = m
        gauss_scale = s

        loglik_gumbel_free = float(np.sum(gumbel_logpdf(M_T, gumbel_loc_mom, gumbel_scale_mom)))
        loglik_gumbel_fix = float(np.sum(gumbel_logpdf(M_T, gumbel_loc_fix, 0.5)))
        loglik_gauss_free = float(np.sum(normal_logpdf(M_T, gauss_loc, gauss_scale)))

        out.append(f"\n--- T = {T_str} (K = {len(M_T)}) ---")
        out.append(f"  M_T mean = {m:+.4f}  95% CI [{mean_lo:+.4f}, {mean_hi:+.4f}]")
        out.append(f"  M_T var  = {s*s:.4f}  95% CI [{var_lo:.4f}, {var_hi:.4f}]")
        out.append(f"  FHK-Gumbel(1/2) predicted variance = {math.pi**2/24:.4f}")
        out.append(f"  Free Gumbel MLE: loc={gumbel_loc_mom:+.4f}, scale={gumbel_scale_mom:.4f}")
        out.append(f"  Log-likelihood (per-window means):")
        out.append(f"    free Gumbel(loc, scale)     = {loglik_gumbel_free / len(M_T):+.4f}")
        out.append(f"    FHK Gumbel(loc, 1/2)        = {loglik_gumbel_fix / len(M_T):+.4f}")
        out.append(f"    free Gaussian(mean, std)    = {loglik_gauss_free / len(M_T):+.4f}")
        # AIC / BIC: same parameter count for free Gumbel (2) vs free Gauss (2)
        # so log-likelihood difference is the head-to-head verdict.
        delta_LL_GumbelFree_vs_GaussFree = loglik_gumbel_free - loglik_gauss_free
        out.append(f"  ΔLL (free Gumbel − free Gaussian) = {delta_LL_GumbelFree_vs_GaussFree:+.3f}  "
                   f"({'GUMBEL' if delta_LL_GumbelFree_vs_GaussFree > 0 else 'GAUSS'} preferred)")
        # Likelihood ratio test (Vuong-style, 1-DoF z)
        # increments per sample
        inc_gumbel = gumbel_logpdf(M_T, gumbel_loc_mom, gumbel_scale_mom)
        inc_gauss = normal_logpdf(M_T, gauss_loc, gauss_scale)
        d = np.asarray(inc_gumbel) - np.asarray(inc_gauss)
        z_vuong = float(np.mean(d) / (np.std(d, ddof=1) / math.sqrt(len(d)))) if np.std(d) > 0 else float("nan")
        out.append(f"  Vuong z (Gumbel vs Gauss): {z_vuong:+.3f}  "
                   f"({'GUMBEL preferred at >2σ' if z_vuong > 2 else ('GAUSS preferred at >2σ' if z_vuong < -2 else 'inconclusive')})")
        # KS comparison
        sorted_M = np.sort(M_T)
        n = len(sorted_M)
        cdf_emp = (np.arange(1, n + 1) - 0.5) / n

        def gumbel_cdf(x, loc, scale):
            return np.exp(-np.exp(-(x - loc) / scale))

        def gauss_cdf(x, loc, scale):
            from math import erf
            return np.array([0.5 * (1 + erf((xi - loc) / (scale * math.sqrt(2)))) for xi in x])

        cdf_gum_free = gumbel_cdf(sorted_M, gumbel_loc_mom, gumbel_scale_mom)
        cdf_gum_fix = gumbel_cdf(sorted_M, gumbel_loc_fix, 0.5)
        cdf_gss = gauss_cdf(sorted_M, gauss_loc, gauss_scale)
        ks_gum_free = float(np.max(np.abs(cdf_emp - cdf_gum_free)))
        ks_gum_fix = float(np.max(np.abs(cdf_emp - cdf_gum_fix)))
        ks_gss = float(np.max(np.abs(cdf_emp - cdf_gss)))
        out.append(f"  KS to free Gumbel  = {ks_gum_free:.4f}")
        out.append(f"  KS to FHK Gumbel(1/2) = {ks_gum_fix:.4f}")
        out.append(f"  KS to free Gaussian = {ks_gss:.4f}")

    # Pooled regression
    out.append("\n--- Three-T regression: <max> = a + b * log log T + c * log log log T ---")
    T_strs = sorted(summary.keys(), key=lambda s: float(s))
    T_floats = [float(s) for s in T_strs]
    xs_loglog = np.array([math.log(math.log(t)) for t in T_floats])
    xs_logloglog = np.array([math.log(math.log(math.log(t))) for t in T_floats])
    ys_max = np.array([summary[t]["raw_max_mean"] for t in T_strs])
    ys_max_sem = np.array(
        [summary[t]["raw_max_std"] / math.sqrt(summary[t]["K_windows"]) for t in T_strs]
    )
    A = np.column_stack([np.ones(len(T_strs)), xs_loglog, xs_logloglog])
    coef, *_ = np.linalg.lstsq(A, ys_max, rcond=None)
    a_, b_, c_ = coef
    out.append(f"  Fitted: a = {a_:+.4f}, b = {b_:+.4f}, c = {c_:+.4f}")
    out.append(f"  FHK predicts:  b = +1.000, c = -0.750  (intercept = location of M_∞)")
    out.append(f"  Selberg-CLT:   b = +1.000, c =  0.000")
    pred_max_FHK = xs_loglog - 0.75 * xs_logloglog
    pred_max_Sel = xs_loglog
    out.append(f"\n  Per-T residuals from FHK (b=1, c=-3/4) form (free intercept):")
    intercept_FHK = float(np.mean(ys_max - pred_max_FHK))
    intercept_Sel = float(np.mean(ys_max - pred_max_Sel))
    res_FHK = ys_max - pred_max_FHK - intercept_FHK
    res_Sel = ys_max - pred_max_Sel - intercept_Sel
    for ts, rfhk, rsel, sem in zip(T_strs, res_FHK, res_Sel, ys_max_sem):
        out.append(f"    T={ts:>10s}: FHK residual = {rfhk:+.4f},  Selberg residual = {rsel:+.4f}  (sem ~ {sem:.4f})")
    rms_FHK = float(np.sqrt(np.mean(res_FHK ** 2)))
    rms_Sel = float(np.sqrt(np.mean(res_Sel ** 2)))
    out.append(f"  RMS residual: FHK form = {rms_FHK:.4f}, Selberg form = {rms_Sel:.4f}")
    out.append(f"  M_T common intercept (under FHK form): {intercept_FHK:+.4f}")

    # Pooled stats
    out.append("\n--- Pooled M_T across all T anchors ---")
    pooled = np.asarray(pooled_M_T)
    out.append(f"  N pooled = {len(pooled)}")
    out.append(f"  pooled mean = {pooled.mean():+.4f}  std = {pooled.std(ddof=1):.4f}")
    out.append(f"  pooled var  = {pooled.var(ddof=1):.4f}  vs FHK Gumbel(1/2) = {math.pi**2/24:.4f}")
    out.append(f"  pooled skewness = {((pooled - pooled.mean())**3).mean() / pooled.std()**3:.4f}  vs Gumbel = +1.139")
    out.append(f"  pooled excess kurtosis = {((pooled - pooled.mean())**4).mean() / pooled.std()**4 - 3:.4f}  vs Gumbel = +2.4")

    # Print to stdout and to file
    text = "\n".join(out)
    print(text)
    out_path = os.path.join(HERE, "fhk_amplitude_max_analysis.log")
    with open(out_path, "w") as f:
        f.write(text + "\n")
    print(f"\nWritten to {out_path}")


if __name__ == "__main__":
    main()

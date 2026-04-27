"""
ATTACK_VECTORS.md §C7 — Fyodorov-Hiary-Keating extreme-value statistics
of |zeta(1/2 + it)| over short windows.

Cross-domain technique: Gaussian multiplicative chaos / log-correlated
random fields. Channelled mathematician: Bourgain (extreme-value /
Fourier-analytic side of FHK).

PRE-STATED FALSIFICATION CRITERIA (also recorded in results.md):

F1 (E mode, B-grade): empirical M_T :=
    max_{t in [T, T+1]} log|zeta(1/2 + it)| - log log T + (3/4) log log log T
matches the FHK prediction within sample noise:
    - mean(M_T) consistent with Gumbel shift (no hidden T-dependence
      after subtracting the FHK normalisation)
    - var(M_T) consistent with pi^2/24 = 0.4112 (Gumbel scale 1/2)
    - Kolmogorov-Smirnov distance to fitted Gumbel < random-Gaussian
      surrogate distance.

F2 (A-grade): empirical M_T systematically deviates from FHK by > 5 sigma
in any of the three signatures (leading log log T constant, secondary
-3/4 log log log T correction, Gumbel tail) AND the deviation has a
structural arithmetic explanation.

F3 (I mode): empirical max systematically matches the *plain Selberg
CLT* extreme-value prediction (no -3/4 log log log T correction;
log|zeta| treated as pointwise N(0, (1/2) log log T) with effective
sample count log T per unit window) BETTER than FHK -- would mean the
freezing transition / log-correlation structure is not visible at
finite T.

CLOSURE LOGIC:
- F1 holds: B-grade close, refines E7.1 with first quantitative
  finite-T FHK confirmation; promotes CROSS_DOMAIN_TECHNIQUES "GMC /
  FHK extreme-value" PROPOSED -> USED (E).
- F2 holds: A-grade open, novel/ entry, new EDGE on amplitude side.
- F3 holds: B-grade close in the OPPOSITE direction -- structural
  distinction between log-correlated GMC (FHK) and naive-CLT regimes;
  closes "FHK is not finite-T-detectable, plain Selberg suffices."

Setup:
- Two T_base values: 1e4, 1e5.
- K = 80 unit-length windows per T_base, spaced by 10 to avoid
  overlap and inter-window correlations.
- M = 200 evenly-spaced samples per window (spacing 0.005, well below
  the mean zero spacing 2 pi / log(T/2 pi) ~ 0.85).
- Precision dps = 15 (sufficient for finding maxima of log|zeta|).

OUTPUT:
- raw per-window data: max, argmax_offset, mean, variance, second_max,
  log_abs_zeta percentiles.
- aggregated statistics per T_base.
- Gumbel / Gaussian / Selberg-CLT comparison summary.
"""

import json
import math
import os
import sys
import time

import mpmath
import numpy as np

mpmath.mp.dps = 15  # plenty for max-finding


HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS_JSON = os.path.join(HERE, "fhk_amplitude_max_results.json")
LOG_PATH = os.path.join(HERE, "fhk_amplitude_max.log")


def log_msg(s, log_file=None):
    print(s, flush=True)
    if log_file is not None:
        log_file.write(s + "\n")
        log_file.flush()


def log_abs_zeta(t):
    z = mpmath.zeta(mpmath.mpc("0.5", str(t)))
    av = abs(z)
    if av == 0:
        return -1e30
    return float(mpmath.log(av))


def sample_window(T_anchor, M=200):
    """Sample log|zeta(1/2 + it)| at M evenly-spaced points across
    [T_anchor, T_anchor + 1]."""
    dt = 1.0 / M
    samples = np.empty(M)
    for j in range(M):
        t = T_anchor + (j + 0.5) * dt
        samples[j] = log_abs_zeta(t)
    return samples


def per_window_stats(samples):
    M = len(samples)
    s_max = float(np.max(samples))
    s_argmax = int(np.argmax(samples))
    s_pos = (s_argmax + 0.5) / M  # offset within unit window in [0, 1]
    # second-max: largest sample whose argmax is at least 5 indices away
    # (avoid trivial near-neighbour pickup)
    mask = np.ones(M, dtype=bool)
    lo = max(0, s_argmax - 5)
    hi = min(M, s_argmax + 6)
    mask[lo:hi] = False
    if mask.sum() > 0:
        s_secmax = float(np.max(samples[mask]))
    else:
        s_secmax = float(np.partition(samples, -2)[-2])
    return {
        "max": s_max,
        "argmax_offset": s_pos,
        "second_max": s_secmax,
        "mean": float(np.mean(samples)),
        "var": float(np.var(samples)),
        "p99": float(np.quantile(samples, 0.99)),
        "p50": float(np.quantile(samples, 0.50)),
        "min": float(np.min(samples)),
    }


def run_anchor(T_base, K=80, M=200, window_spacing=10.0, log_file=None):
    """Run K unit windows starting at T_base + k * window_spacing."""
    results = []
    start = time.time()
    for k in range(K):
        T_anchor = T_base + window_spacing * k
        samples = sample_window(T_anchor, M)
        stats = per_window_stats(samples)
        stats["k"] = k
        stats["T_anchor"] = T_anchor
        results.append(stats)
        if (k + 1) % 5 == 0 or k == 0:
            elapsed = time.time() - start
            eta = elapsed / (k + 1) * (K - k - 1)
            log_msg(
                f"  T_base={T_base:.0e}: window {k+1}/{K}, "
                f"elapsed {elapsed:.0f}s, ETA {eta:.0f}s, "
                f"max={stats['max']:.3f}",
                log_file=log_file,
            )
    return results


def fhk_normalised(maxes, T_base):
    """Return M_T = max - log log T + (3/4) log log log T."""
    log_T = math.log(T_base)
    loglog_T = math.log(log_T)
    logloglog_T = math.log(loglog_T) if loglog_T > 0 else 0.0
    return [m - loglog_T + 0.75 * logloglog_T for m in maxes]


def selberg_normalised(maxes, T_base, K_eff_per_unit=None):
    """Plain Selberg-CLT prediction without the -3/4 log log log
    correction. log|zeta(1/2 + it)| ~ N(0, (1/2) log log T) pointwise.
    Treat each unit window as ~log T independent Gaussians (correlation
    length ~ 1/log T). Then max ~ sigma sqrt(2 log K_eff) where
    sigma^2 = (1/2) log log T."""
    log_T = math.log(T_base)
    loglog_T = math.log(log_T)
    sigma = math.sqrt(0.5 * loglog_T)
    if K_eff_per_unit is None:
        K_eff_per_unit = log_T
    pred_mean = sigma * math.sqrt(2.0 * math.log(K_eff_per_unit))
    return pred_mean


def gumbel_fit(values):
    """Method-of-moments fit of Gumbel(loc, scale=1/2) [FHK predicts
    scale=1/2, which gives variance pi^2/24]. We *also* report a free
    scale fit. Returns (loc_fixed_scale, free_scale, free_loc,
    KS_distance_to_fixed, KS_distance_to_free)."""
    arr = np.asarray(values, dtype=float)
    n = len(arr)

    # Free Gumbel(loc, scale) MOM fit:
    s_emp = float(np.std(arr, ddof=1))
    m_emp = float(np.mean(arr))
    free_scale = s_emp * math.sqrt(6.0) / math.pi
    euler_gamma = 0.5772156649015329
    free_loc = m_emp - euler_gamma * free_scale

    # Fixed scale=0.5 fit (FHK prediction):
    fixed_scale = 0.5
    loc_fixed = m_emp - euler_gamma * fixed_scale

    # KS distance to fitted distributions
    sorted_x = np.sort(arr)
    cdf_emp = (np.arange(1, n + 1) - 0.5) / n

    def gumbel_cdf(x, loc, scale):
        return np.exp(-np.exp(-(x - loc) / scale))

    cdf_fix = gumbel_cdf(sorted_x, loc_fixed, fixed_scale)
    cdf_free = gumbel_cdf(sorted_x, free_loc, free_scale)
    KS_fix = float(np.max(np.abs(cdf_emp - cdf_fix)))
    KS_free = float(np.max(np.abs(cdf_emp - cdf_free)))

    # Anderson-Darling-like statistic for tail emphasis (optional,
    # implemented as max in upper-half tail)
    upper = sorted_x[n // 2 :]
    cdf_emp_upper = cdf_emp[n // 2 :]
    cdf_fix_upper = gumbel_cdf(upper, loc_fixed, fixed_scale)
    KS_fix_upper = float(np.max(np.abs(cdf_emp_upper - cdf_fix_upper)))

    return {
        "n": n,
        "empirical_mean": m_emp,
        "empirical_std": s_emp,
        "empirical_var": s_emp * s_emp,
        "predicted_var_FHK_scale_half": (math.pi ** 2) / 24.0,
        "free_scale_MOM": free_scale,
        "free_loc_MOM": free_loc,
        "fixed_scale_FHK": fixed_scale,
        "loc_fixed_scale": loc_fixed,
        "KS_distance_fixed_scale_half": KS_fix,
        "KS_distance_free_scale": KS_free,
        "KS_distance_upper_half_fixed": KS_fix_upper,
    }


def gaussian_compare(values):
    arr = np.asarray(values, dtype=float)
    n = len(arr)
    m = float(np.mean(arr))
    s = float(np.std(arr, ddof=1))
    sorted_x = np.sort(arr)
    cdf_emp = (np.arange(1, n + 1) - 0.5) / n
    # Normal CDF
    cdf_norm = 0.5 * (1 + np.array([math.erf((x - m) / (s * math.sqrt(2))) for x in sorted_x]))
    KS_norm = float(np.max(np.abs(cdf_emp - cdf_norm)))
    return {
        "KS_distance_to_fitted_Gaussian": KS_norm,
        "skewness": float(((arr - m) ** 3).mean() / (s ** 3)),
        "excess_kurtosis": float(((arr - m) ** 4).mean() / (s ** 4) - 3),
    }


def aggregate(per_T_results):
    summary = {}
    for T_base, results in per_T_results.items():
        T_base_f = float(T_base)
        log_T = math.log(T_base_f)
        loglog_T = math.log(log_T)
        logloglog_T = math.log(loglog_T)
        maxes = [r["max"] for r in results]
        argmax_offsets = [r["argmax_offset"] for r in results]
        secmaxes = [r["second_max"] for r in results]
        means = [r["mean"] for r in results]
        variances = [r["var"] for r in results]

        # FHK-renormalised maxes
        M_T = [m - loglog_T + 0.75 * logloglog_T for m in maxes]
        # Selberg-renormalised (no log log log correction)
        Selberg_resid = [m - loglog_T for m in maxes]
        # Plain max minus log log T (which is the Selberg leading constant)
        # Difference of (FHK - Selberg) is exactly +0.75 log log log T.

        gumbel_stats = gumbel_fit(M_T)
        gaussian_stats_M = gaussian_compare(M_T)

        # Selberg-CLT predicted naive-Gaussian-extreme value (mean only)
        sigma = math.sqrt(0.5 * loglog_T)
        K_eff = log_T
        Selberg_max_mean = sigma * math.sqrt(2.0 * math.log(K_eff))
        Selberg_max_var = sigma * sigma  # rough, Mills-ratio order
        # Argmax uniformity test (KS vs uniform[0,1])
        sorted_off = np.sort(argmax_offsets)
        n = len(sorted_off)
        cdf_emp = (np.arange(1, n + 1) - 0.5) / n
        cdf_uniform = sorted_off
        KS_argmax_uniform = float(np.max(np.abs(cdf_emp - cdf_uniform)))

        summary[T_base] = {
            "T_base": T_base_f,
            "log_T": log_T,
            "loglog_T": loglog_T,
            "logloglog_T": logloglog_T,
            "FHK_subtracted_term": loglog_T - 0.75 * logloglog_T,
            "Selberg_subtracted_term": loglog_T,
            "K_windows": len(results),
            "raw_max_mean": float(np.mean(maxes)),
            "raw_max_std": float(np.std(maxes, ddof=1)),
            "raw_max_min": float(np.min(maxes)),
            "raw_max_max": float(np.max(maxes)),
            "M_T_mean": float(np.mean(M_T)),
            "M_T_std": float(np.std(M_T, ddof=1)),
            "M_T_var": float(np.var(M_T, ddof=1)),
            "M_T_min": float(np.min(M_T)),
            "M_T_max": float(np.max(M_T)),
            "Selberg_resid_mean": float(np.mean(Selberg_resid)),
            "Selberg_resid_std": float(np.std(Selberg_resid, ddof=1)),
            "Selberg_predicted_max_mean": Selberg_max_mean,
            "Selberg_excess_over_predicted": float(np.mean(maxes)) - Selberg_max_mean,
            "argmax_offset_mean": float(np.mean(argmax_offsets)),
            "argmax_offset_var": float(np.var(argmax_offsets, ddof=1)),
            "KS_argmax_to_uniform": KS_argmax_uniform,
            "second_max_mean": float(np.mean(secmaxes)),
            "second_max_to_max_gap_mean": float(
                np.mean([m - sm for m, sm in zip(maxes, secmaxes)])
            ),
            "log_abs_zeta_pointwise_mean": float(np.mean(means)),
            "log_abs_zeta_pointwise_var_mean": float(np.mean(variances)),
            "selberg_predicted_pointwise_var": 0.5 * loglog_T,
            "gumbel_fit": gumbel_stats,
            "gaussian_compare_M_T": gaussian_stats_M,
        }
    return summary


def main():
    log_file = open(LOG_PATH, "w")
    try:
        log_msg("=" * 60, log_file)
        log_msg("FHK |zeta(1/2 + it)| max-statistics — C7 attack vector", log_file)
        log_msg("Cross-domain: Gaussian multiplicative chaos (Saksman-Webb 2018)", log_file)
        log_msg("=" * 60, log_file)

        T_base_values = [1e4, 1e5, 1e6]
        K = int(os.environ.get("FHK_K", 100))
        M = int(os.environ.get("FHK_M", 200))
        log_msg(f"Parameters: K = {K} windows per T, M = {M} samples per window", log_file)
        log_msg(f"T_base values: {T_base_values}", log_file)
        log_msg("", log_file)

        per_T_raw = {}
        for T_base in T_base_values:
            log_msg(f"--- Running T_base = {T_base:.0e} ---", log_file)
            t0 = time.time()
            results = run_anchor(T_base, K=K, M=M, window_spacing=10.0, log_file=log_file)
            elapsed = time.time() - t0
            log_msg(f"  Done in {elapsed:.0f}s ({elapsed/K:.1f}s per window)", log_file)
            per_T_raw[str(T_base)] = results

        log_msg("", log_file)
        log_msg("Aggregating ...", log_file)
        summary = aggregate(per_T_raw)

        # Print summary table
        log_msg("", log_file)
        log_msg("=" * 60, log_file)
        log_msg("SUMMARY", log_file)
        log_msg("=" * 60, log_file)
        for T_str, stats in summary.items():
            log_msg(f"\nT_base = {T_str}", log_file)
            log_msg(f"  log log T               = {stats['loglog_T']:.4f}", log_file)
            log_msg(f"  log log log T           = {stats['logloglog_T']:.4f}", log_file)
            log_msg(f"  FHK subtraction         = {stats['FHK_subtracted_term']:.4f}", log_file)
            log_msg(f"  raw max mean            = {stats['raw_max_mean']:.4f} ± {stats['raw_max_std']:.4f}", log_file)
            log_msg(f"  M_T mean (= max - FHK)  = {stats['M_T_mean']:.4f} ± {stats['M_T_std']:.4f}", log_file)
            log_msg(f"  M_T var                 = {stats['M_T_var']:.4f}  (FHK Gumbel(1/2): {stats['gumbel_fit']['predicted_var_FHK_scale_half']:.4f})", log_file)
            log_msg(f"  Selberg pred max mean   = {stats['Selberg_predicted_max_mean']:.4f}", log_file)
            log_msg(f"  raw max - Selberg pred  = {stats['Selberg_excess_over_predicted']:.4f}", log_file)
            log_msg(f"  argmax mean             = {stats['argmax_offset_mean']:.4f}", log_file)
            log_msg(f"  KS to Uniform (argmax)  = {stats['KS_argmax_to_uniform']:.4f}", log_file)
            log_msg(f"  Gumbel free scale (MOM) = {stats['gumbel_fit']['free_scale_MOM']:.4f} (FHK: 0.5000)", log_file)
            log_msg(f"  KS to FHK Gumbel(loc, 1/2)         = {stats['gumbel_fit']['KS_distance_fixed_scale_half']:.4f}", log_file)
            log_msg(f"  KS to free Gumbel(loc, scale)      = {stats['gumbel_fit']['KS_distance_free_scale']:.4f}", log_file)
            log_msg(f"  KS to fitted Gaussian (M_T)        = {stats['gaussian_compare_M_T']['KS_distance_to_fitted_Gaussian']:.4f}", log_file)
            log_msg(f"  M_T skewness            = {stats['gaussian_compare_M_T']['skewness']:.4f}  (Gumbel: 1.139)", log_file)
            log_msg(f"  M_T excess kurtosis     = {stats['gaussian_compare_M_T']['excess_kurtosis']:.4f}  (Gumbel: 2.4)", log_file)
            log_msg(f"  pointwise log|zeta| var = {stats['log_abs_zeta_pointwise_var_mean']:.4f} (Selberg pred: {stats['selberg_predicted_pointwise_var']:.4f})", log_file)

        # Cross-T tests
        log_msg("", log_file)
        log_msg("Cross-T comparison (FHK predicts M_T mean ~ T-independent):", log_file)
        T_strs = sorted(summary.keys(), key=lambda s: float(s))
        T_floats = [float(s) for s in T_strs]
        if len(T_strs) >= 2:
            mt_means = [summary[t]["M_T_mean"] for t in T_strs]
            mt_stds = [summary[t]["M_T_std"] for t in T_strs]
            mt_vars = [summary[t]["M_T_var"] for t in T_strs]
            n_per = [summary[t]["K_windows"] for t in T_strs]
            sb_means = [summary[t]["Selberg_resid_mean"] for t in T_strs]
            sb_stds = [summary[t]["Selberg_resid_std"] for t in T_strs]

            for ts, m, s, n in zip(T_strs, mt_means, mt_stds, n_per):
                log_msg(
                    f"  T={ts}: M_T mean = {m:.4f}, sem = {s/math.sqrt(n):.4f}, "
                    f"var = {s*s:.4f}", log_file)

            log_msg("", log_file)
            log_msg("FHK predicts M_T -> Gumbel-distributed limit (T-independent).", log_file)
            log_msg(f"  Gumbel(loc, 1/2) variance = pi^2/24 = {math.pi**2/24:.4f}.", log_file)
            log_msg("", log_file)

            # Pairwise Z-test on M_T means
            log_msg("Pairwise Z(M_T mean diff) -- if FHK holds, all should be small:", log_file)
            for i in range(len(T_strs)):
                for j in range(i+1, len(T_strs)):
                    sd = math.sqrt(mt_stds[i]**2/n_per[i] + mt_stds[j]**2/n_per[j])
                    z = (mt_means[j] - mt_means[i]) / sd
                    log_msg(f"  Z(M_T(T={T_strs[j]}) - M_T(T={T_strs[i]})) = {z:.3f}", log_file)

            # Selberg-resid (= max - log log T) should drop by
            # +0.75 (log log log T2 - log log log T1) under FHK
            log_msg("", log_file)
            log_msg("Selberg-resid drops (under FHK should equal -0.75 d(log log log T)):", log_file)
            for i in range(len(T_strs)):
                for j in range(i+1, len(T_strs)):
                    T1 = T_floats[i]; T2 = T_floats[j]
                    obs = sb_means[j] - sb_means[i]
                    pred = -0.75 * (math.log(math.log(T2)) - math.log(math.log(T1)))
                    sd = math.sqrt(sb_stds[i]**2/n_per[i] + sb_stds[j]**2/n_per[j])
                    z_dev = (obs - pred) / sd
                    log_msg(
                        f"  T={T_strs[i]}->T={T_strs[j]}: observed drop = {obs:+.4f}, "
                        f"FHK pred = {pred:+.4f}, Z(obs - pred) = {z_dev:+.3f}",
                        log_file)

            # Linear regression test: max ~ a + b * log log T + c * log log log T
            # to identify empirically the leading and secondary coefficients.
            from numpy.polynomial import polynomial as P
            xs_loglog = np.array([math.log(math.log(t)) for t in T_floats])
            xs_logloglog = np.array([math.log(math.log(math.log(t))) for t in T_floats])
            ys_max = np.array([summary[t]["raw_max_mean"] for t in T_strs])
            # Fit max = a + b * loglog + c * logloglog if 3 points;
            # with only 3 T-values the system is exactly determined.
            log_msg("", log_file)
            log_msg("Empirical regression: <max> = a + b * log log T + c * log log log T", log_file)
            if len(T_strs) >= 3:
                A = np.column_stack([np.ones(len(T_strs)), xs_loglog, xs_logloglog])
                coef, *_ = np.linalg.lstsq(A, ys_max, rcond=None)
                a_, b_, c_ = coef
                log_msg(f"  fitted: a={a_:+.4f}, b={b_:+.4f} (FHK pred 1.0), c={c_:+.4f} (FHK pred -0.75)", log_file)
                log_msg("  Note: 3 points fits exactly; over-determination needs more T values.", log_file)
            else:
                # 2 points: solve b assuming c = -0.75 and intercept 0
                # Just report ratio
                log_msg("  (Fewer than 3 T values; skipping over-determined regression.)", log_file)

        # Write JSON
        out = {
            "params": {"K": K, "M": M, "T_base_values": T_base_values, "dps": 15},
            "raw_per_T": per_T_raw,
            "summary": summary,
        }
        with open(RESULTS_JSON, "w") as f:
            json.dump(out, f, indent=2, default=str)
        log_msg(f"\nWrote {RESULTS_JSON}", log_file)
    finally:
        log_file.close()


if __name__ == "__main__":
    main()

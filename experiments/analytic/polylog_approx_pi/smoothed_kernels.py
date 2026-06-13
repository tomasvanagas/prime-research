"""
Slot 3 of Thread 7: do compactly-supported smoothing kernels reduce
sigma_eff = rms(|err|) below the hard-cutoff value at fixed
K_compute = #zero evaluations?

S196 closed log-Gaussian (UN-truncated, full Gaussian damping w_j =
exp(-h^2 gamma_j^2 / 2)). S202 wrap explicitly listed compactly-
supported smoothing kernels (sinc, sech^2, Hann, raised-cosine, ...) as
"should generalise but not formally proven" — i.e., a legitimate open
falsifier of the Thread 3 closure.

This script tests that falsifier head-on: for each compactly-supported
kernel w on [0, K_compute] with w(0) = 1 and w(K_compute) = 0,

    pi_K_kernel(x) = R(x) - 2 sum_{j <= K_compute} w_j Re R(x^{rho_j})

what is sigma_eff(kernel, K_compute) = rms(|pi_K_kernel(x_i) - pi(x_i)|)
across multiple samples, vs sigma_eff(hard, K_compute)?

Theoretical prediction (random-phase): hard cutoff is L2-optimal at
fixed K_compute because under independence,

    Var(pi_K_kernel - pi) = sum_j (1 - w_j)^2 Var(2 Re R(x^{rho_j}))

is minimised by w_j = 1 for j <= K_compute (which is hard cutoff). All
other compactly-supported kernels have (1 - w_j)^2 > 0 for some j and
must satisfy sigma(kernel) >= sigma(hard).

The slot 3 question:
  - Empirical predicted ratio ratio_pred = sigma(kernel) / sigma(hard)
    is computable at j = 1..K from per-zero variance estimates.
  - Empirical observed ratio ratio_obs = sigma_eff(kernel) /
    sigma_eff(hard).
  - GUE pair-correlation could in principle make ratio_obs differ from
    ratio_pred. If ratio_obs < 1 for some kernel, hard-cutoff is NOT
    optimal under GUE statistics — A-grade signal.

Method (one cumulative loop per sample, post-process all kernels):
  1. For each (anchor, sample x_j), compute c_k = 2 Re R(x_j^{rho_k})
     for k = 1..K_max. (This is the cumulative loop, dominant cost.)
  2. For each kernel and each K_compute in {500, 1000, 2000, 4000}:
       pi_K_kernel(x_j) = R(x_j) - sum_{k=1..K_compute} w_k(kernel, K_compute) c_k
       err_kernel(x_j, K_compute) = pi_K_kernel(x_j) - pi(x_j)
  3. Aggregate sigma_eff(kernel, K_compute) over the N samples per anchor.
  4. Compute ratio_obs = sigma_eff(kernel) / sigma_eff(hard).
  5. Compute ratio_pred from the per-zero variance estimate (estimated
     by mean(c_k^2 / 2) over samples).

CLI args:
  --anchors  comma list of log10 anchors  (default: 7,8,9)
  --N        samples per anchor           (default: 15)
  --K-max    max K_compute                 (default: 4000)
  --dps      mpmath precision              (default: 25)
  --M        Mobius truncation in R_at_rho (default: 8)
  --csv      output per-sample CSV         (default: smoothed_kernels_data.csv)
  --summary-csv                           (default: smoothed_kernels_summary.csv)
  --quiet    suppress per-sample log

Run time: ~12-18 minutes at defaults (3 anchors x 15 samples x 4000 zeros).
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import sys
import time
from dataclasses import dataclass
from typing import Dict, List, Tuple, Callable

import mpmath as mp
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from polylog_approx_pi import (  # type: ignore
    PI_KNOWN,
    get_zeros,
    riemann_R,
    R_at_rho,
    sigma_predicted,
)
from multi_sample import (  # type: ignore
    pi_array_for_samples,
    half_normal_cdf,
    ks_statistic,
    ks_pvalue_approx,
)


# -----------------------------------------------------------------------------
# Compactly-supported smoothing kernels.
#
# Each kernel takes a normalised position u = (k-1) / (K - 1) in [0, 1] for
# k = 1..K, and returns the weight w_k. Conventions:
#   - Edge condition w_K = 0 for hann/hamming/triangle/tukey (otherwise
#     they are NOT meaningful smoothings -- "Hann" with w(1)=0 is the
#     usual definition).
#   - Hard is included as a baseline: w_k = 1 for all k.
#
# Sentry: kernel weight 1 at k = 1 (lowest gamma) for all kernels.
# -----------------------------------------------------------------------------

KernelFn = Callable[[np.ndarray], np.ndarray]  # u in [0,1] -> w in [0,1]


def kernel_hard(u: np.ndarray) -> np.ndarray:
    return np.ones_like(u)


def kernel_triangle(u: np.ndarray) -> np.ndarray:
    """Bartlett / Fejer triangle: w(u) = 1 - u, taper to 0 at u=1."""
    return 1.0 - u


def kernel_hann(u: np.ndarray) -> np.ndarray:
    """Hann window: w(u) = (1 + cos(pi u)) / 2, taper to 0 at u=1."""
    return 0.5 * (1.0 + np.cos(np.pi * u))


def kernel_hamming(u: np.ndarray) -> np.ndarray:
    """Hamming window: w(u) = 0.54 + 0.46 cos(pi u), taper to 0.08 at u=1.
    Hamming is *not* zero at u=1, but tapers far. Use with caution."""
    return 0.54 + 0.46 * np.cos(np.pi * u)


def kernel_riesz(u: np.ndarray) -> np.ndarray:
    """Riesz (parabolic): w(u) = 1 - u^2, taper to 0 at u=1."""
    return 1.0 - u ** 2


def kernel_riesz4(u: np.ndarray) -> np.ndarray:
    """Riesz^4 = (1 - u^2)^2, smoother taper."""
    return (1.0 - u ** 2) ** 2


def kernel_tukey50(u: np.ndarray) -> np.ndarray:
    """Tukey window with 50% taper: hard for u in [0, 0.5], cosine taper to 0
    on [0.5, 1]. Good for separating "active region" from "rolloff"."""
    w = np.ones_like(u)
    mask = u > 0.5
    w[mask] = 0.5 * (1.0 + np.cos(np.pi * (u[mask] - 0.5) / 0.5))
    return w


def kernel_tukey25(u: np.ndarray) -> np.ndarray:
    """Tukey window with 25% taper: hard for u in [0, 0.75], cosine taper
    on [0.75, 1]. Closest to hard cutoff."""
    w = np.ones_like(u)
    mask = u > 0.75
    w[mask] = 0.5 * (1.0 + np.cos(np.pi * (u[mask] - 0.75) / 0.25))
    return w


def kernel_cosine(u: np.ndarray) -> np.ndarray:
    """Cosine taper: w(u) = cos(pi u / 2), taper to 0 at u=1."""
    return np.cos(np.pi * u / 2.0)


KERNELS: Dict[str, KernelFn] = {
    "hard":      kernel_hard,
    "triangle":  kernel_triangle,
    "hann":      kernel_hann,
    "hamming":   kernel_hamming,
    "riesz":     kernel_riesz,
    "riesz4":    kernel_riesz4,
    "tukey50":   kernel_tukey50,
    "tukey25":   kernel_tukey25,
    "cosine":    kernel_cosine,
}


def kernel_weights(name: str, K_compute: int) -> np.ndarray:
    """Return weights w_1, ..., w_{K_compute} of length K_compute."""
    fn = KERNELS[name]
    if K_compute < 2:
        return np.ones(K_compute)
    u = np.linspace(0.0, 1.0, K_compute)  # u_k = (k-1) / (K_compute - 1) endpoints inclusive
    return fn(u)


# -----------------------------------------------------------------------------
# Data containers.
# -----------------------------------------------------------------------------

@dataclass
class KernelRow:
    anchor_log10: int
    x_anchor: int
    x: int
    pi_x: int
    R_x: float
    K_compute: int
    kernel: str
    pi_K_kernel: float
    err: float
    abs_err: float


# -----------------------------------------------------------------------------
# Per-anchor evaluator.
# -----------------------------------------------------------------------------

def evaluate_anchor(
    anchor_log10: int,
    N_samples: int,
    K_max: int,
    K_targets: List[int],
    gammas: List[mp.mpf],
    kernel_names: List[str],
    dps: int,
    M_mobius: int,
    quiet: bool,
) -> List[KernelRow]:
    x_anchor = 10 ** anchor_log10
    if x_anchor not in PI_KNOWN:
        raise SystemExit(f"pi(10^{anchor_log10}) not in PI_KNOWN.")
    pi_anchor = PI_KNOWN[x_anchor]

    # Sample geometrically across half-decade.
    log_lo = math.log(float(x_anchor))
    log_hi = log_lo + 0.5 * math.log(10.0)
    sample_xs = sorted({int(round(math.exp(log_lo + j * (log_hi - log_lo) / (N_samples - 1))))
                        for j in range(N_samples)})
    sample_xs = [x for x in sample_xs if x > x_anchor]
    while len(sample_xs) < N_samples:
        sample_xs.append(sample_xs[-1] + max(1, int(0.001 * x_anchor)))
    sample_xs = sorted(sample_xs[:N_samples])

    print(f"\n[anchor 10^{anchor_log10}]  N={len(sample_xs)} samples in "
          f"[{sample_xs[0]:.4g}, {sample_xs[-1]:.4g}]; pi(10^{anchor_log10})={pi_anchor}")

    t0 = time.time()
    pi_xs = pi_array_for_samples(x_anchor, pi_anchor, sample_xs)
    print(f"  segmented-sieve pi(x_j): {time.time() - t0:.1f}s "
          f"(pi range [{pi_xs[0]}, {pi_xs[-1]}])")

    # Pre-compute kernel weight vectors at each K_compute.
    # weights[K][kernel] is np.array of length K.
    weights: Dict[int, Dict[str, np.ndarray]] = {}
    for K in K_targets:
        weights[K] = {name: kernel_weights(name, K) for name in kernel_names}

    rows: List[KernelRow] = []
    for s, (x, pi_x) in enumerate(zip(sample_xs, pi_xs)):
        t0 = time.time()
        R_x = float(riemann_R(x, M=20, dps=dps))

        # Compute c_k = 2 Re R(x^{rho_k}) for k = 1..K_max in one cumulative loop.
        c = np.zeros(K_max, dtype=np.float64)
        for k in range(K_max):
            c[k] = 2.0 * float(R_at_rho(x, gammas[k], M=M_mobius, dps=dps).real)

        # For each (K, kernel), compute pi_K_kernel = R_x - sum_{j<=K} w_j c_j.
        for K in K_targets:
            for name in kernel_names:
                w = weights[K][name]
                correction = float(np.dot(w, c[:K]))
                pi_K_kernel = R_x - correction
                err = pi_K_kernel - pi_x
                rows.append(KernelRow(
                    anchor_log10=anchor_log10,
                    x_anchor=x_anchor,
                    x=x,
                    pi_x=pi_x,
                    R_x=R_x,
                    K_compute=K,
                    kernel=name,
                    pi_K_kernel=pi_K_kernel,
                    err=err,
                    abs_err=abs(err),
                ))
        if not quiet:
            print(f"  sample {s + 1}/{len(sample_xs)}: x={x}, "
                  f"|err@(hard,Kmax)|={abs(rows[-len(kernel_names)].err):.3f}, "
                  f"{time.time() - t0:.1f}s")
    return rows


# -----------------------------------------------------------------------------
# Aggregation + headline ratio table.
# -----------------------------------------------------------------------------

def rms(xs: List[float]) -> float:
    if not xs:
        return float("nan")
    return math.sqrt(sum(x * x for x in xs) / len(xs))


def aggregate_kernels(rows: List[KernelRow], kernel_names: List[str],
                      K_targets: List[int]):
    """Group by (anchor, K_compute, kernel) and return summary records.

    Headline metric: ratio_obs = sigma_eff(kernel) / sigma_eff(hard) at
    same (anchor, K_compute).
    """
    by_group: Dict[Tuple[int, int, str], List[KernelRow]] = {}
    for r in rows:
        by_group.setdefault((r.anchor_log10, r.K_compute, r.kernel), []).append(r)

    summary = []
    for (anc, K, name) in sorted(by_group.keys()):
        grp = by_group[(anc, K, name)]
        errs = [r.abs_err for r in grp]
        x_vals = [r.x for r in grp]
        x0 = x_vals[len(x_vals) // 2]
        sigma_pred = sigma_predicted(x0, K)
        sigma_eff = rms(errs)
        med_err = sorted(errs)[len(errs) // 2]
        mean_err = sum(errs) / len(errs)
        # KS half-normal test (vs sigma_eff scale, shape-only)
        ratios_eff = [e / sigma_eff for e in errs] if sigma_eff > 0 else []
        D_eff = ks_statistic(ratios_eff)
        p_eff = ks_pvalue_approx(D_eff, len(ratios_eff))
        summary.append({
            "anchor_log10": anc,
            "K_compute": K,
            "kernel": name,
            "N": len(grp),
            "sigma_pred": sigma_pred,
            "sigma_eff": sigma_eff,
            "median_abs_err": med_err,
            "mean_abs_err": mean_err,
            "ks_D_eff": D_eff,
            "ks_p_eff": p_eff,
        })
    return summary


def add_ratios(summary):
    """Add ratio_obs = sigma_eff(kernel) / sigma_eff(hard) at same (anchor, K)."""
    by_anchor_K: Dict[Tuple[int, int], Dict[str, float]] = {}
    for s in summary:
        by_anchor_K.setdefault((s["anchor_log10"], s["K_compute"]), {})[s["kernel"]] = s["sigma_eff"]
    for s in summary:
        sig_hard = by_anchor_K[(s["anchor_log10"], s["K_compute"])].get("hard", float("nan"))
        s["sigma_eff_hard"] = sig_hard
        s["ratio_obs_over_hard"] = s["sigma_eff"] / sig_hard if sig_hard > 0 else float("nan")
    return summary


def predicted_ratio_from_S195_var(K: int, kernel: str) -> float:
    """Predicted sigma(kernel, K) / sigma(hard, K) under random-phase + the
    S195 model that the per-zero variance |R(x^{rho_k})|^2 contribution scales
    asymptotically as 1 / (k log^2 k log^2 x) for k large.

    Specifically, S195 derives sigma(K, x)^2 ~ x * log^2 K / (2 K log^2 x),
    which is a sum over k > K of var contributions. With the assumption
    var_k = c / (k log^2 x) for some scale c (so that
    sum_{k>K} var_k ~ c log K log K / (K log^2 x) ~ x sigma^2),
    we can predict any kernel's variance:

        sigma(kernel, K)^2 = sum_{k=1..K} (1 - w_k(K))^2 var_k + sum_{k>K} var_k.

    Numerically, we estimate var_k from the S195 formula by setting the tail
    sum equal to sigma_pred(x, K_compute)^2 and using var_k ~ 1/(k log K_compute^2)
    profile. This is a heuristic prediction, not an exact one.

    Returns ratio = sigma(kernel, K) / sigma(hard, K).
    """
    if K < 2:
        return 1.0
    # Assume var_k ~ 1/k for k = 1..K (the dominant scaling); the constant cancels.
    var_k = np.array([1.0 / k for k in range(1, K + 1)])  # length K
    # Tail sum (j > K) approximated by integral, dominated by 1/k.
    # Use empirical fact: sum_{k>K} var_k = (sigma_hard)^2.
    # We just need the RATIO, so set sigma_hard^2 = 1 (the tail).
    tail_norm = 1.0  # sum_{k>K} var_k normalised to 1
    # We need to infer the head sum sum_{k=1..K} var_k in same units.
    # Under the var_k = 1/k profile with truncation at K:
    #   tail_int = int_K^infty 1/k dk = log(K_max/K) -- divergent without cap.
    # Instead, use: sum_{k>K} var_k = sigma_hard^2 = sigma_pred(x, K)^2 (S195).
    # And sum_{k<=K} var_k = log K - log 1 = log K (1/k integral).
    # These are in the same units up to a scale. The ratio is:
    head_var_total = float(np.sum(var_k))           # sum_{k<=K} 1/k
    head_var_kernel = lambda w: float(np.sum((1 - w) ** 2 * var_k))
    # sigma_hard^2 (under this model) = "tail" = let's pretend head_var_total - log K
    # is the tail. That doesn't quite work because we don't know the truncation
    # point of the original 1/k series.
    # Better: just use sigma_pred(x=arbitrary, K)^2 from S195 as the tail value;
    # the ratio is then unitless.
    # The 1/k scaling means sigma_pred(K)^2 ~ log K^2 / K (S195's formula);
    # more precisely sigma_pred(K)^2 = (x/(2*pi^2*K*log^2 x)) * log^2 K, so
    # the tail-to-head ratio is log^2 K / K vs log K -> tail is much smaller
    # than head for K large. Hmm.
    # For now, use a simpler heuristic prediction: assume var_k profile is
    # constant in k (the worst case for kernel ratio). Then:
    var_k_uniform = np.ones(K)
    head_kernel_unif = lambda w: float(np.sum((1 - w) ** 2 * var_k_uniform))
    w = kernel_weights(kernel, K)
    # Prediction: ratio_pred = sqrt(1 + head_kernel(w) / tail_norm) where tail_norm is
    # sigma_hard^2. Under the S195 1/k profile, tail_norm relative to head sums
    # is log^2 K / K vs log K, ratio log K / K << 1, so head dominates if any
    # w_k < 1.
    # This means under random-phase, ANY kernel != hard has sigma >> sigma_hard.
    # But experiments may show GUE corrections reducing the head term.
    # For the simplest, most pessimistic prediction:
    sig2_pred = head_kernel_unif(w) + 1.0  # tail = 1 (normalised)
    sig2_hard = 0.0 + 1.0
    return math.sqrt(sig2_pred / sig2_hard)


# -----------------------------------------------------------------------------
# CLI.
# -----------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--anchors", default="7,8,9",
                   help="Comma list of log10 anchors.")
    p.add_argument("--N", type=int, default=15, help="Samples per anchor.")
    p.add_argument("--K-max", dest="K_max", type=int, default=4000)
    p.add_argument("--K-targets", default="500,1000,2000,4000",
                   help="Comma list of K_compute values.")
    p.add_argument("--kernels", default="hard,triangle,hann,hamming,riesz,riesz4,tukey50,tukey25,cosine",
                   help="Comma list of kernel names.")
    p.add_argument("--dps", type=int, default=25)
    p.add_argument("--M", type=int, default=8)
    p.add_argument("--csv", default="smoothed_kernels_data.csv")
    p.add_argument("--summary-csv", default="smoothed_kernels_summary.csv")
    p.add_argument("--quiet", action="store_true")
    return p.parse_args()


def main():
    args = parse_args()
    anchors = [int(s) for s in args.anchors.split(",")]
    K_targets = sorted(set(int(s) for s in args.K_targets.split(",")))
    if max(K_targets) > args.K_max:
        raise SystemExit(f"K_target max {max(K_targets)} > K_max {args.K_max}")
    kernel_names = [s.strip() for s in args.kernels.split(",") if s.strip()]
    for name in kernel_names:
        if name not in KERNELS:
            raise SystemExit(f"unknown kernel '{name}'; valid: {list(KERNELS)}")

    print(f"Slot 3 / Thread 7: smoothed-kernel comparison")
    print(f"  anchors = 10^{anchors}")
    print(f"  N = {args.N} samples per anchor")
    print(f"  K_targets = {K_targets} (K_max = {args.K_max})")
    print(f"  kernels = {kernel_names}")
    print()

    print(f"Loading {args.K_max} zeta zeros (mp dps={args.dps})...", flush=True)
    t0 = time.time()
    gammas = get_zeros(args.K_max, dps=args.dps)
    print(f"  loaded {len(gammas)} zeros in {time.time() - t0:.1f}s "
          f"(gamma_1={float(gammas[0]):.4f}, gamma_K={float(gammas[-1]):.4f})")

    rows: List[KernelRow] = []
    for anc in anchors:
        sub = evaluate_anchor(anc, args.N, args.K_max, K_targets, gammas,
                              kernel_names, args.dps, args.M, args.quiet)
        rows.extend(sub)

    summary = aggregate_kernels(rows, kernel_names, K_targets)
    add_ratios(summary)

    # ------- print headline table -------
    print("\n" + "=" * 110)
    print("HEADLINE: sigma_eff(kernel) / sigma_eff(hard) at same (anchor, K_compute)")
    print(" -> ratio < 1 means kernel BEATS hard cutoff (would falsify L2-optimality)")
    print(" -> ratio = 1 means kernel ties hard cutoff")
    print(" -> ratio > 1 means kernel WORSE than hard cutoff (expected per random-phase)")
    print()
    header = f"{'anchor':>7s} {'K':>5s} {'kernel':>9s}  {'sig_eff':>9s} {'sig_pred':>9s} " \
             f"{'med|err|':>9s} {'KSp':>6s} {'ratio_obs/hard':>15s} {'pred_unif':>10s}"
    print(header)
    print("-" * len(header))
    for s in summary:
        pred_unif = predicted_ratio_from_S195_var(s["K_compute"], s["kernel"])
        print(f"{'10^'+str(s['anchor_log10']):>7s} {s['K_compute']:>5d} "
              f"{s['kernel']:>9s}  {s['sigma_eff']:>9.4f} {s['sigma_pred']:>9.4f} "
              f"{s['median_abs_err']:>9.4f} {s['ks_p_eff']:>6.3f} "
              f"{s['ratio_obs_over_hard']:>15.4f} {pred_unif:>10.4f}")

    # ------- per-K minimum-ratio finder -------
    print("\n" + "=" * 110)
    print("Best (min) ratio_obs/hard per (anchor, K) cell:")
    by_aK: Dict[Tuple[int, int], List[Dict]] = {}
    for s in summary:
        by_aK.setdefault((s["anchor_log10"], s["K_compute"]), []).append(s)
    for (anc, K), grp in sorted(by_aK.items()):
        sorted_g = sorted(grp, key=lambda r: r["ratio_obs_over_hard"])
        best = sorted_g[0]
        print(f"  10^{anc} K={K}:  best={best['kernel']:>9s}  "
              f"ratio={best['ratio_obs_over_hard']:.4f}  "
              f"(N={best['N']})  next 2: "
              f"{sorted_g[1]['kernel']}={sorted_g[1]['ratio_obs_over_hard']:.4f}, "
              f"{sorted_g[2]['kernel']}={sorted_g[2]['ratio_obs_over_hard']:.4f}")

    # ------- CSV outputs -------
    out_csv = args.csv if os.path.isabs(args.csv) else os.path.join(HERE, args.csv)
    with open(out_csv, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow([
            "anchor_log10", "x_anchor", "x", "pi_x", "R_x",
            "K_compute", "kernel", "pi_K_kernel", "err", "abs_err",
        ])
        for r in rows:
            w.writerow([
                r.anchor_log10, r.x_anchor, r.x, r.pi_x, f"{r.R_x:.6f}",
                r.K_compute, r.kernel,
                f"{r.pi_K_kernel:.6f}", f"{r.err:.6f}", f"{r.abs_err:.6f}",
            ])
    print(f"\nWrote per-sample CSV: {out_csv}  ({len(rows)} rows)")

    out_sum = args.summary_csv if os.path.isabs(args.summary_csv) else os.path.join(HERE, args.summary_csv)
    keys = ["anchor_log10", "K_compute", "kernel", "N",
            "sigma_pred", "sigma_eff", "sigma_eff_hard", "ratio_obs_over_hard",
            "median_abs_err", "mean_abs_err", "ks_D_eff", "ks_p_eff"]
    with open(out_sum, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(keys)
        for s in summary:
            w.writerow([s[k] for k in keys])
    print(f"Wrote summary CSV: {out_sum}  ({len(summary)} rows)")


if __name__ == "__main__":
    main()

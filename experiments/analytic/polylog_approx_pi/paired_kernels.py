"""
Slot 4 of Thread 7 (option a): do non-symmetric / position-correlated
("paired") weight kernels reduce sigma_eff = rms(|err|) below the
hard-cutoff value at fixed K_compute = #zero evaluations?

Slot 3 (S242) tested 8 symmetric compactly-supported smoothing kernels
(Hann, Hamming, Triangle, Riesz, Riesz4, Tukey25, Tukey50, Cosine) and
found *zero* of 96 (anchor, K, kernel) cells where kernel beats hard at
paired sign-test p < 0.05. The structural reason: under random-phase,
hard cutoff is L2-optimal at fixed K_compute, and GUE corrections
scale sigma_eff uniformly across symmetric kernels (factor ~0.74).

Slot 4 tests the only remaining kernel-family direction not yet covered:
**paired weights** w_{2j-1}, w_{2j} that explicitly use the position
parity of the zero index. The hypothesis: under GUE Wigner repulsion
(sin^2(pi t)/(pi t)^2 small-t structure on consecutive-zero spacing),
adjacent-zero pairs c_{2j-1} = 2 Re R(x^{rho_{2j-1}}) and c_{2j} = 2 Re
R(x^{rho_{2j}}) are anti-correlated. A paired kernel that exploits this
correlation might reduce sigma_eff below hard cutoff.

The 7 paired families tested here are:

  1. paired_hann       w_j = hann_{half-resolution}(ceil(j/2))
  2. paired_triangle   w_j = triangle_{half-resolution}(ceil(j/2))
  3. paired_riesz      w_j = riesz_{half-resolution}(ceil(j/2))
  4. antipair_03       w_{2j-1} = 1.3, w_{2j} = 0.7
  5. antipair_05       w_{2j-1} = 1.5, w_{2j} = 0.5
  6. half_int          hard for j < K, w_K = 0.5      (cutoff K - 0.5)
  7. boundary_pair     hard for j < K-1, w_{K-1}=0.75, w_K=0.25

A-grade signal for slot 4:
  * paired kernel != hard with sigma_eff_paired / sigma_eff_hard < 0.95
    at paired sign-test p < 0.05 in ANY (anchor, K_compute) cell.

If no such cell exists, slot 4 is a B-grade negative-shape closure of
the only remaining kernel-axis lever. Slot 5 would then be the
theoretical wrap to close Thread 7 as a partial-positive named-exponent
result with no unconditional algorithmic improvement available within
the linear partial-sum framework.

CLI args mirror smoothed_kernels.py.

Run time: ~12-18 minutes at defaults (3 anchors x 20 samples x 4000 zeros).
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
# Paired-kernel weight functions.
#
# All take K_compute and return np.array of length K_compute with entries
# w_1, ..., w_{K_compute}. Sentry: w_1 is always 1 (lowest gamma untouched).
# -----------------------------------------------------------------------------

KernelFn = Callable[[int], np.ndarray]  # K_compute -> w of length K_compute


def kernel_hard(K: int) -> np.ndarray:
    return np.ones(K)


def kernel_paired_hann(K: int) -> np.ndarray:
    """Hann at pair resolution: w_{2j-1} = w_{2j} = hann((j-1)/(K/2 - 1)),
    i.e. one Hann value per pair, applied to both members.
    For odd K, the last (singleton) zero is treated as a half-pair with
    its own weight at position (K+1)/2."""
    K_pairs = (K + 1) // 2
    if K_pairs < 2:
        return np.ones(K)
    u_pairs = np.linspace(0.0, 1.0, K_pairs)        # length K_pairs
    w_pairs = 0.5 * (1.0 + np.cos(np.pi * u_pairs))  # Hann on pairs
    # Repeat each pair-weight twice; clip to length K.
    w = np.repeat(w_pairs, 2)[:K]
    return w


def kernel_paired_triangle(K: int) -> np.ndarray:
    K_pairs = (K + 1) // 2
    if K_pairs < 2:
        return np.ones(K)
    u_pairs = np.linspace(0.0, 1.0, K_pairs)
    w_pairs = 1.0 - u_pairs                         # triangle on pairs
    w = np.repeat(w_pairs, 2)[:K]
    return w


def kernel_paired_riesz(K: int) -> np.ndarray:
    K_pairs = (K + 1) // 2
    if K_pairs < 2:
        return np.ones(K)
    u_pairs = np.linspace(0.0, 1.0, K_pairs)
    w_pairs = 1.0 - u_pairs ** 2                    # Riesz on pairs
    w = np.repeat(w_pairs, 2)[:K]
    return w


def kernel_antipair_03(K: int) -> np.ndarray:
    """w_{2j-1} = 1.3, w_{2j} = 0.7. Sum-per-pair preserved at 2.0."""
    w = np.empty(K)
    w[0::2] = 1.3
    w[1::2] = 0.7
    return w


def kernel_antipair_05(K: int) -> np.ndarray:
    """w_{2j-1} = 1.5, w_{2j} = 0.5. Stronger anti-pair perturbation."""
    w = np.empty(K)
    w[0::2] = 1.5
    w[1::2] = 0.5
    return w


def kernel_half_int(K: int) -> np.ndarray:
    """Hard for j < K, w_K = 0.5: half-integer cutoff at K - 0.5.
    Effectively integrates between K-1 and K zeros via boundary smoothing."""
    w = np.ones(K)
    if K >= 1:
        w[K - 1] = 0.5
    return w


def kernel_boundary_pair(K: int) -> np.ndarray:
    """Hard for j < K-1, w_{K-1} = 0.75, w_K = 0.25.
    Linear taper applied only to the boundary pair (K-1, K)."""
    w = np.ones(K)
    if K >= 2:
        w[K - 2] = 0.75
        w[K - 1] = 0.25
    return w


KERNELS: Dict[str, KernelFn] = {
    "hard":            kernel_hard,
    "paired_hann":     kernel_paired_hann,
    "paired_triangle": kernel_paired_triangle,
    "paired_riesz":    kernel_paired_riesz,
    "antipair_03":     kernel_antipair_03,
    "antipair_05":     kernel_antipair_05,
    "half_int":        kernel_half_int,
    "boundary_pair":   kernel_boundary_pair,
}


def kernel_weights(name: str, K_compute: int) -> np.ndarray:
    fn = KERNELS[name]
    return fn(K_compute)


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
# Per-anchor evaluator (one cumulative R_at_rho loop; post-process all kernels).
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

    log_lo = math.log(float(x_anchor))
    log_hi = log_lo + 0.5 * math.log(10.0)
    sample_xs = sorted({int(round(math.exp(log_lo + j * (log_hi - log_lo) / (N_samples - 1))))
                        for j in range(N_samples)})
    sample_xs = [x for x in sample_xs if x > x_anchor]
    while len(sample_xs) < N_samples:
        sample_xs.append(sample_xs[-1] + max(1, int(0.001 * x_anchor)))
    sample_xs = sorted(sample_xs[:N_samples])

    print(f"\n[anchor 10^{anchor_log10}]  N={len(sample_xs)} samples in "
          f"[{sample_xs[0]:.4g}, {sample_xs[-1]:.4g}]; "
          f"pi(10^{anchor_log10})={pi_anchor}", flush=True)

    t0 = time.time()
    pi_xs = pi_array_for_samples(x_anchor, pi_anchor, sample_xs)
    print(f"  segmented-sieve pi(x_j): {time.time() - t0:.1f}s "
          f"(pi range [{pi_xs[0]}, {pi_xs[-1]}])", flush=True)

    weights: Dict[int, Dict[str, np.ndarray]] = {}
    for K in K_targets:
        weights[K] = {name: kernel_weights(name, K) for name in kernel_names}

    rows: List[KernelRow] = []
    for s, (x, pi_x) in enumerate(zip(sample_xs, pi_xs)):
        t0 = time.time()
        R_x = float(riemann_R(x, M=20, dps=dps))

        c = np.zeros(K_max, dtype=np.float64)
        for k in range(K_max):
            c[k] = 2.0 * float(R_at_rho(x, gammas[k], M=M_mobius, dps=dps).real)

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
            hard_err = abs(rows[-len(kernel_names) - 0 + kernel_names.index("hard")].err)
            print(f"  sample {s + 1}/{len(sample_xs)}: x={x}, "
                  f"|err@(hard,Kmax)|={hard_err:.3f}, "
                  f"{time.time() - t0:.1f}s", flush=True)
    return rows


# -----------------------------------------------------------------------------
# Aggregation.
# -----------------------------------------------------------------------------

def rms(xs: List[float]) -> float:
    if not xs:
        return float("nan")
    return math.sqrt(sum(x * x for x in xs) / len(xs))


def aggregate_kernels(rows: List[KernelRow], kernel_names: List[str],
                      K_targets: List[int]):
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
    by_anchor_K: Dict[Tuple[int, int], Dict[str, float]] = {}
    for s in summary:
        by_anchor_K.setdefault((s["anchor_log10"], s["K_compute"]), {})[s["kernel"]] = s["sigma_eff"]
    for s in summary:
        sig_hard = by_anchor_K[(s["anchor_log10"], s["K_compute"])].get("hard", float("nan"))
        s["sigma_eff_hard"] = sig_hard
        s["ratio_obs_over_hard"] = s["sigma_eff"] / sig_hard if sig_hard > 0 else float("nan")
    return summary


# -----------------------------------------------------------------------------
# Paired-sign-test post-processor (built in to keep slot self-contained).
#
# For each (anchor, K, kernel != hard), and each sample i, ask whether
# |err_kernel(x_i)| < |err_hard(x_i)|. wins/N is the binomial estimate;
# one-sided p-value under Binom(N, 0.5) for "kernel beats hard".
# -----------------------------------------------------------------------------

def binom_one_sided_le(k: int, n: int, p: float = 0.5) -> float:
    """P(X <= k) under X ~ Binom(n, p). Used as p-value for k = wins, testing
    the null that kernel >= hard (so wins is small).  Equivalently
    P(kernel beats hard) p-value = P(Bin(n, 0.5) <= n - wins) by symmetry,
    BUT we want to report 'kernel beats hard' as: wins close to n is
    significant (small p)."""
    from math import comb
    return sum(comb(n, j) * (p ** j) * ((1 - p) ** (n - j)) for j in range(k + 1))


def paired_sign_test(rows: List[KernelRow], kernel_names: List[str],
                     K_targets: List[int]) -> List[Dict]:
    by_anc_K: Dict[Tuple[int, int], Dict[int, Dict[str, float]]] = {}
    for r in rows:
        d = by_anc_K.setdefault((r.anchor_log10, r.K_compute), {})
        s = d.setdefault(r.x, {})
        s[r.kernel] = r.abs_err

    out: List[Dict] = []
    for (anc, K), per_x in sorted(by_anc_K.items()):
        for name in kernel_names:
            if name == "hard":
                continue
            wins = 0
            ties = 0
            n = 0
            for x, errs in per_x.items():
                e_h = errs.get("hard")
                e_k = errs.get(name)
                if e_h is None or e_k is None:
                    continue
                n += 1
                if e_k < e_h:
                    wins += 1
                elif e_k == e_h:
                    ties += 1
            if n == 0:
                continue
            p_kernel_beats_hard = 1.0 - binom_one_sided_le(wins - 1, n) if wins > 0 else 1.0
            p_hard_beats_kernel = binom_one_sided_le(wins, n)
            out.append({
                "anchor_log10": anc,
                "K_compute": K,
                "kernel": name,
                "N": n,
                "wins_kernel": wins,
                "ties": ties,
                "p_kernel_beats_hard": p_kernel_beats_hard,
                "p_hard_beats_kernel": p_hard_beats_kernel,
            })
    return out


# -----------------------------------------------------------------------------
# CLI.
# -----------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--anchors", default="7,8,9")
    p.add_argument("--N", type=int, default=20)
    p.add_argument("--K-max", dest="K_max", type=int, default=4000)
    p.add_argument("--K-targets", default="500,1000,2000,4000")
    p.add_argument("--kernels",
                   default="hard,paired_hann,paired_triangle,paired_riesz,"
                           "antipair_03,antipair_05,half_int,boundary_pair")
    p.add_argument("--dps", type=int, default=25)
    p.add_argument("--M", type=int, default=8)
    p.add_argument("--csv", default="paired_kernels_data.csv")
    p.add_argument("--summary-csv", default="paired_kernels_summary.csv")
    p.add_argument("--paired-csv", default="paired_kernels_signtest.csv")
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
    if "hard" not in kernel_names:
        raise SystemExit("kernel set must include 'hard' baseline")

    print(f"Slot 4 / Thread 7: paired-kernel comparison")
    print(f"  anchors = 10^{anchors}")
    print(f"  N = {args.N} samples per anchor")
    print(f"  K_targets = {K_targets} (K_max = {args.K_max})")
    print(f"  kernels = {kernel_names}")
    print()

    print(f"Loading {args.K_max} zeta zeros (mp dps={args.dps})...", flush=True)
    t0 = time.time()
    gammas = get_zeros(args.K_max, dps=args.dps)
    print(f"  loaded {len(gammas)} zeros in {time.time() - t0:.1f}s "
          f"(gamma_1={float(gammas[0]):.4f}, gamma_K={float(gammas[-1]):.4f})",
          flush=True)

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
    print()
    header = (f"{'anchor':>7s} {'K':>5s} {'kernel':>16s}  "
              f"{'sig_eff':>9s} {'sig_pred':>9s} {'med|err|':>9s} "
              f"{'KSp':>6s} {'ratio_obs/hard':>15s}")
    print(header)
    print("-" * len(header))
    for s in summary:
        print(f"{'10^'+str(s['anchor_log10']):>7s} {s['K_compute']:>5d} "
              f"{s['kernel']:>16s}  {s['sigma_eff']:>9.4f} {s['sigma_pred']:>9.4f} "
              f"{s['median_abs_err']:>9.4f} {s['ks_p_eff']:>6.3f} "
              f"{s['ratio_obs_over_hard']:>15.4f}")

    print("\n" + "=" * 110)
    print("Best (min) ratio_obs/hard per (anchor, K) cell:")
    by_aK: Dict[Tuple[int, int], List[Dict]] = {}
    for s in summary:
        by_aK.setdefault((s["anchor_log10"], s["K_compute"]), []).append(s)
    for (anc, K), grp in sorted(by_aK.items()):
        non_hard = [r for r in grp if r["kernel"] != "hard"]
        sorted_g = sorted(non_hard, key=lambda r: r["ratio_obs_over_hard"])
        best = sorted_g[0] if sorted_g else None
        if best is not None:
            print(f"  10^{anc} K={K}:  best={best['kernel']:>16s}  "
                  f"ratio={best['ratio_obs_over_hard']:.4f}  "
                  f"(N={best['N']})")

    # ------- paired sign test -------
    paired = paired_sign_test(rows, kernel_names, K_targets)
    print("\n" + "=" * 110)
    print("PAIRED SIGN TEST (one-sided): wins_kernel/N + p-values")
    print(f"{'anchor':>7s} {'K':>5s} {'kernel':>16s}  "
          f"{'wins/N':>8s} {'p_k>h':>8s} {'p_h>k':>8s}")
    print("-" * 60)
    for r in paired:
        print(f"{'10^'+str(r['anchor_log10']):>7s} {r['K_compute']:>5d} "
              f"{r['kernel']:>16s}  {r['wins_kernel']}/{r['N']:<6d} "
              f"{r['p_kernel_beats_hard']:>8.4f} {r['p_hard_beats_kernel']:>8.4f}")

    # Headline: any kernel-beats-hard at p < 0.05?
    any_beat = [r for r in paired if r["p_kernel_beats_hard"] < 0.05]
    print("\n" + "=" * 110)
    if any_beat:
        print(f"** A-GRADE SIGNAL: {len(any_beat)} (anchor, K, kernel) cells "
              f"have kernel-beats-hard at p < 0.05:")
        for r in any_beat:
            print(f"  10^{r['anchor_log10']} K={r['K_compute']} {r['kernel']}: "
                  f"wins={r['wins_kernel']}/{r['N']}, p={r['p_kernel_beats_hard']:.4f}")
    else:
        print("** B-GRADE NEGATIVE-SHAPE: NO (anchor, K, kernel) cell has "
              "kernel-beats-hard at p < 0.05.")
        # Most "promising" cell (smallest p_kernel_beats_hard).
        best = min(paired, key=lambda r: r["p_kernel_beats_hard"])
        print(f"   Smallest p_kernel_beats_hard: {best['kernel']} at "
              f"10^{best['anchor_log10']} K={best['K_compute']}, "
              f"wins={best['wins_kernel']}/{best['N']}, "
              f"p={best['p_kernel_beats_hard']:.4f}")

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

    out_paired = args.paired_csv if os.path.isabs(args.paired_csv) else os.path.join(HERE, args.paired_csv)
    pkeys = ["anchor_log10", "K_compute", "kernel", "N",
             "wins_kernel", "ties", "p_kernel_beats_hard", "p_hard_beats_kernel"]
    with open(out_paired, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(pkeys)
        for r in paired:
            w.writerow([r[k] for k in pkeys])
    print(f"Wrote paired sign-test CSV: {out_paired}  ({len(paired)} rows)")


if __name__ == "__main__":
    main()

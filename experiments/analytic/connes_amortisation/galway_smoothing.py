"""
Slot 4/5 of commit thread "connes_amortisation" (S196).

Question: does Gaussian-window smoothing of the Riemann explicit formula
drop K* (the truncation index achieving |error| <= 0.5 with positive
hit-rate) from Theta(x) to polylog(x)?

Background (S195, slot 3): under the GUE random-phase heuristic the
unsmoothed truncation
   pi_K(x) = R(x) - 2 sum_{j<=K} Re R(x^rho_j)
has variance Var(pi(x) - pi_K(x)) ~ x log^2(K) / (2 pi^2 K log^2 x),
which forces K = Theta(x) for any positive in-distribution hit-rate.

Galway 2004 / Hiary 2011 use a *smoothed* explicit formula that gets
K = O(x^{1/2 + eps}). The S195 closure left open: does the smoothing
help in-distribution beyond the square-root barrier? Specifically,
could a polylog K work for distributional pi(x) +/- 1?

This script answers NO via the same heuristic. The argument:

Gaussian smoothing in log-coordinates with bandwidth h applies a
weight w_j(h) = exp(-h^2 gamma_j^2 / 2) to each zero. Define the
smoothed truncated approximation:

   pi_{K,h}(x) = R(x) - 2 sum_{j<=K} Re R(x^rho_j) * w_j(h).

Decompose the error:

   pi(x) - pi_{K,h}(x) = -2 sum_{j>K} Re R(x^rho_j)               (TAIL)
                       + 2 sum_{j<=K} Re R(x^rho_j) * (1 - w_j)    (BIAS)

Under iid uniform phases gamma_j log x mod 2pi:

   Var(TAIL)  ~ (2x/log^2 x) * sum_{j>K} 1/gamma_j^2
              ~ x log^2(K) / (2 pi^2 K log^2 x)

   Var(BIAS)  ~ (2x/log^2 x) * sum_{j<=K} (1 - w_j)^2 / gamma_j^2

For w_j = exp(-h^2 gamma_j^2 / 2):
   (1 - w_j) ~ h^2 gamma_j^2 / 2 if h gamma_j << 1
   (1 - w_j) ~ 1                  if h gamma_j >> 1

Asymptotically Var(BIAS) ~ x h log(1/h) / log^2 x for h small,
saturating at the unsmoothed Var(TAIL=0, K=infty)=0 when h=0.

Crucial observation: TAIL and BIAS sum over *disjoint* j ranges
(j > K vs j <= K), so under iid phases their variances add.

Total variance is bounded below by max(Var(TAIL), Var(BIAS)). For
ANY h, Var(TAIL)(K) ~ x log^2(K)/(K log^2 x) >= 1/4 unless K >= x/poly(log).
So smoothing CANNOT lower the K threshold below Theta(x/poly(log)).

This script verifies the heuristic prediction empirically:
  1. For 20 x values in a band, compute Re R(x^rho_j) per zero for
     j = 1..K_max.
  2. For each h in {0, 1e-5, ..., 1e-1}, sweep K = 1..K_max recording
     the filtered cumulative-sum error.
  3. Tabulate K*(h) := smallest K with hit-rate >= p.
  4. Verify K*(h) does NOT drop into polylog regime even for optimal
     h ~ log(x)/x.

Output:
  galway_smoothing_data.csv  -- per-(x, h, K) error grid
  K_star_table               -- K* vs h vs target hit-rate
  bias_table                 -- empirical bias variance vs predicted
"""

from __future__ import annotations

import csv
import math
import os as _os
import sys
import time
from pathlib import Path
from typing import Dict, List, Tuple

# Reuse infrastructure from connes_amortisation.py: zeros, R(x), R_at_rho,
# pi via sieve.
sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from connes_amortisation import (
    R_at_rho,
    build_prime_count_array,
    get_zeros,
    pi_from_array,
    riemann_R,
)

import numpy as np

# ----------------------------------------------------------------------------
# Per-zero contribution computation (the expensive step).
# ----------------------------------------------------------------------------


def per_zero_contributions(x: float, gammas: List[float]) -> np.ndarray:
    """For each gamma_j in gammas, compute c_j(x) := 2 Re R(x^rho_j).

    Returns a numpy array of length len(gammas).
    """
    K = len(gammas)
    out = np.empty(K, dtype=np.float64)
    for j in range(K):
        Rrho = R_at_rho(x, gammas[j])
        out[j] = 2.0 * Rrho.real
    return out


# ----------------------------------------------------------------------------
# Smoothed-truncated cumulative sums and error trajectories.
# ----------------------------------------------------------------------------


def smoothed_cumsum(c_x: np.ndarray, gammas_arr: np.ndarray, h: float) -> np.ndarray:
    """Compute cum_K = sum_{j<=K} c_j * w_j(h) for K = 1..len(c_x).

    w_j(h) = exp(-h^2 gamma_j^2 / 2).
    Returns array of length len(c_x).
    """
    if h <= 0.0:
        weights = np.ones_like(c_x)
    else:
        weights = np.exp(-(h * gammas_arr) ** 2 / 2.0)
    return np.cumsum(c_x * weights)


# ----------------------------------------------------------------------------
# Heuristic predictions.
# ----------------------------------------------------------------------------


def sigma_trunc_pred(x: float, K: int) -> float:
    """sigma^2_TAIL(K) ~ x log^2(K) / (2 pi^2 K log^2 x)."""
    if K <= 1 or x <= 1:
        return float("inf")
    return math.sqrt(x) * math.log(K) / (math.pi * math.sqrt(2 * K) * math.log(x))


def sigma_bias_pred(x: float, h: float, K: int, gammas_arr: np.ndarray) -> float:
    """sigma^2_BIAS(K, h) = (2x/log^2 x) sum_{j<=K} (1 - w_j(h))^2 / gamma_j^2.

    Computed directly from the explicit gamma_j values.
    """
    if x <= 1 or K < 1 or h <= 0.0:
        return 0.0
    gs = gammas_arr[:K]
    w = np.exp(-(h * gs) ** 2 / 2.0)
    var = 2.0 * x / math.log(x) ** 2 * np.sum((1.0 - w) ** 2 / gs ** 2)
    return math.sqrt(var)


def sigma_total_pred(x: float, K: int, h: float, gammas_arr: np.ndarray) -> float:
    """sigma_TAIL + sigma_BIAS in quadrature."""
    s_t = sigma_trunc_pred(x, K)
    s_b = sigma_bias_pred(x, h, K, gammas_arr)
    return math.sqrt(s_t ** 2 + s_b ** 2)


# ----------------------------------------------------------------------------
# Find K* (smallest K achieving threshold hit-rate)
# ----------------------------------------------------------------------------


def k_star_empirical(errs_per_x: np.ndarray, threshold: float = 0.5,
                     hit_target: float = 0.5) -> int:
    """errs_per_x: shape (n_x, K_max), |error| at each (x, K).

    Return smallest K such that fraction of x with |err| <= threshold >= hit_target,
    or -1 if no such K up to K_max.
    """
    n_x, K_max = errs_per_x.shape
    hit_per_K = (errs_per_x <= threshold).mean(axis=0)
    for K in range(1, K_max + 1):
        if hit_per_K[K - 1] >= hit_target:
            return K
    return -1


# ----------------------------------------------------------------------------
# Main experiment
# ----------------------------------------------------------------------------


def run_experiment(center: float = 1e5, n_samples: int = 20,
                   K_max: int = 400,
                   h_values: List[float] = None,
                   output_csv: str = None,
                   summary_csv: str = None):
    if h_values is None:
        # h = 0 is the unsmoothed baseline. Other h log-spaced including the
        # critical transition near h ~ 1/gamma_K_max ~ 4e-4 (for K=2000).
        h_values = [0.0, 1e-6, 1e-4, 3e-4, 5e-4, 1e-3, 2e-3, 5e-3, 1e-2,
                    3e-2, 1e-1]

    lo = center
    hi = center * math.sqrt(10.0)
    xs = list(np.geomspace(lo, hi, n_samples))
    sieve_limit = int(hi) + 100

    print(f"GALWAY-SMOOTHING: center={center}, n={n_samples}, "
          f"x in [{lo:.0f}, {hi:.0f}], K_max={K_max}, "
          f"|h|={len(h_values)}")

    print(f"Loading {K_max} zeros...")
    t0 = time.time()
    gammas = get_zeros(K_max, dps=15)
    gammas_arr = np.asarray(gammas, dtype=np.float64)
    print(f"  loaded in {time.time() - t0:.1f}s "
          f"(gamma_1={gammas[0]:.4f}, gamma_K={gammas[-1]:.4f})")

    print(f"Sieve up to {sieve_limit}...")
    t0 = time.time()
    cum = build_prime_count_array(sieve_limit)
    print(f"  done in {time.time() - t0:.1f}s")

    # Compute per-zero contributions for each x (the slow step).
    print(f"\nComputing per-zero contributions for {n_samples} x values...")
    contrib = np.empty((n_samples, K_max), dtype=np.float64)
    pi_xs = np.empty(n_samples, dtype=np.int64)
    R_xs = np.empty(n_samples, dtype=np.float64)
    for i, x in enumerate(xs):
        t1 = time.time()
        pi_xs[i] = int(cum[int(x)])
        R_xs[i] = riemann_R(x)
        contrib[i, :] = per_zero_contributions(x, gammas)
        print(f"  [{i+1}/{n_samples}] x={x:.0f}, pi(x)={pi_xs[i]}, "
              f"t={time.time() - t1:.1f}s")

    # Now build error grid: errs[i_h, i_x, K-1] = |err| at (h, x, K).
    print(f"\nBuilding error grid for {len(h_values)} h values, {n_samples} "
          f"x, K_max={K_max}...")
    errs = np.empty((len(h_values), n_samples, K_max), dtype=np.float64)
    for ih, h in enumerate(h_values):
        if h <= 0:
            weights = np.ones(K_max, dtype=np.float64)
        else:
            weights = np.exp(-(h * gammas_arr) ** 2 / 2.0)
        for i in range(n_samples):
            cumsum = np.cumsum(contrib[i, :] * weights)
            approx = R_xs[i] - cumsum  # shape (K_max,)
            errs[ih, i, :] = np.abs(approx - pi_xs[i])

    # Save the per-(h, x, K) data only for diagnostic K subset.
    if output_csv:
        with open(output_csv, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["h", "x", "K", "abs_err"])
            K_grid = sorted(set(
                [1, 5, 10, 25, 50, 100, 150, 200, 250, 300, K_max]
                + [int(round(math.log(xs[0]) ** p)) for p in (2, 3)]
            ))
            K_grid = [K for K in K_grid if 1 <= K <= K_max]
            for ih, h in enumerate(h_values):
                for i in range(n_samples):
                    for K in K_grid:
                        w.writerow([h, xs[i], K, errs[ih, i, K - 1]])
        print(f"Wrote {output_csv}")

    # K*(h, p) sweep: for each h, find smallest K with hit-rate >= p.
    print(f"\n=== K* vs h (smallest K with hit-rate >= p) ===")
    print(f"{'h':>10s} {'K*_p20':>10s} {'K*_p50':>10s} {'K*_p80':>10s} "
          f"{'predTailK*50':>14s}")
    print("-" * 60)
    summary_rows = []
    for ih, h in enumerate(h_values):
        K20 = k_star_empirical(errs[ih], hit_target=0.20)
        K50 = k_star_empirical(errs[ih], hit_target=0.50)
        K80 = k_star_empirical(errs[ih], hit_target=0.80)
        # Theoretical K* from sigma_tail alone (truncation):
        # P(|err|<=0.5) >= 0.5 with sigma_tail = 0.5/(sqrt(2)*erfinv(0.5)),
        # erfinv(0.5) ~ 0.4769.
        sig_target = 0.5 / (math.sqrt(2) * 0.476936)
        x_med = float(np.median(xs))
        # Solve sigma_tail(x_med, K) <= sig_target.
        K_pred = -1
        for K in range(2, K_max + 1):
            if sigma_trunc_pred(x_med, K) <= sig_target:
                K_pred = K
                break
        print(f"{h:>10.1e} {K20:>10d} {K50:>10d} {K80:>10d} "
              f"{K_pred:>14d}")
        summary_rows.append(dict(h=h, K_star_p20=K20, K_star_p50=K50,
                                 K_star_p80=K80, K_tail_pred=K_pred))

    # Bias variance comparison (predicted vs empirical).
    print(f"\n=== Bias variance: at K=K_max for each h, sigma_BIAS predicted vs empirical ===")
    print(f"{'h':>10s} {'predBIAS':>10s} {'empBIAS_rms':>12s} {'predTAIL':>10s} "
          f"{'empTOTAL_rms':>13s} {'predTOTAL':>11s}")
    print("-" * 75)
    bias_rows = []
    for ih, h in enumerate(h_values):
        x_med = float(np.median(xs))
        # Use very-large K (K=K_max) so TAIL is small; remaining error is BIAS.
        emp_rms_total = math.sqrt(np.mean(errs[ih, :, K_max - 1] ** 2))
        pred_bias = sigma_bias_pred(x_med, h, K_max, gammas_arr)
        pred_tail = sigma_trunc_pred(x_med, K_max)
        pred_total = math.sqrt(pred_bias ** 2 + pred_tail ** 2)
        # Empirical bias estimate: total^2 - tail_pred^2 (assuming tail
        # variance prediction is right).
        var_diff = max(0.0, emp_rms_total ** 2 - pred_tail ** 2)
        emp_bias_rms = math.sqrt(var_diff)
        print(f"{h:>10.1e} {pred_bias:>10.3f} {emp_bias_rms:>12.3f} "
              f"{pred_tail:>10.3f} {emp_rms_total:>13.3f} {pred_total:>11.3f}")
        bias_rows.append(dict(h=h, pred_bias=pred_bias,
                              emp_bias_rms=emp_bias_rms,
                              pred_tail=pred_tail,
                              emp_total_rms=emp_rms_total,
                              pred_total=pred_total))

    # Save summary CSV.
    if summary_csv:
        with open(summary_csv, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["h", "K_star_p20", "K_star_p50", "K_star_p80",
                        "K_tail_pred",
                        "pred_bias_sigma", "emp_bias_rms",
                        "pred_tail_sigma", "emp_total_rms", "pred_total_sigma"])
            for sr, br in zip(summary_rows, bias_rows):
                w.writerow([sr["h"], sr["K_star_p20"], sr["K_star_p50"],
                            sr["K_star_p80"], sr["K_tail_pred"],
                            br["pred_bias"], br["emp_bias_rms"],
                            br["pred_tail"], br["emp_total_rms"],
                            br["pred_total"]])
        print(f"\nWrote {summary_csv}")

    # Final scaling check: what is K* at p=50% as h varies? Does it
    # ever drop into polylog territory (K << sqrt(x))?
    print(f"\n=== Scaling diagnostic (x_med = {float(np.median(xs)):.3g}) ===")
    x_med = float(np.median(xs))
    sqrt_x = math.sqrt(x_med)
    log2x = math.log(x_med) ** 2
    log3x = math.log(x_med) ** 3
    print(f"sqrt(x) = {sqrt_x:.1f}   log^2(x) = {log2x:.1f}   "
          f"log^3(x) = {log3x:.1f}   K_max = {K_max}")
    print(f"Polylog candidates: K = log^2(x) = {int(log2x)}, "
          f"log^3(x) = {int(log3x)}")
    print(f"Square-root scale:  K = sqrt(x) = {int(sqrt_x)}, "
          f"sqrt(x) log(x) = {int(sqrt_x * math.log(x_med))}")
    print()
    for sr in summary_rows:
        h = sr["h"]
        K = sr["K_star_p50"]
        if K <= 0:
            label = "(no hit, K* > K_max)"
        elif K < log3x:
            label = "POLYLOG"
        elif K < sqrt_x:
            label = "sub-sqrt(x)"
        else:
            label = "sqrt(x) or above"
        print(f"  h={h:>9.1e}  K*_50 = {K:>4d}  ({label})")

    return summary_rows, bias_rows


def main():
    here = Path(__file__).parent
    run_experiment(
        center=1e5,
        n_samples=40,
        K_max=2000,
        output_csv=str(here / "galway_smoothing_data.csv"),
        summary_csv=str(here / "galway_smoothing_summary.csv"),
    )


if __name__ == "__main__":
    main()

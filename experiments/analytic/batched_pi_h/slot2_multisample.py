"""
slot2_multisample.py — Thread 8 (P2) commit slot 2

Multi-anchor empirical error characterisation of the Hardy-Littlewood
approximation HL_h(x) = S_h * li_2(x) for the prime-gap-counting
function pi_h(x) = #{p <= x : p, p+h both prime}.

Slot 1 (S418) established the structural dichotomy:
  - EXACT pi_h(x) per-h amortised cost = Theta(x/log x), not polylog.
  - APPROX HL_h(x) per-h cost = polylog, with empirical error <= 0.25 sqrt(x).

Slot 2 sharpens the approximate-regime characterisation. Single-script:

  For each anchor x_anchor in {10^6, 10^7, 10^8}:
    - Sample N=30 values x_j log-uniform in [x_anchor, x_anchor * 10^0.25].
    - For each h in {2, 6, 14, 22, 30, 42, 90, 198} (eight representative
      values spanning S_h in [1.32, 4.0]):
        - Compute exact pi_h(x_j) via single sieve + numpy AND-shift.
        - Compute HL_h(x_j) = S_h * li_2(x_j) via scipy.integrate.quad.
        - Record signed residual r = pi_h(x_j) - HL_h(x_j).
    - Aggregate per-(anchor, h):
        sigma_eff      = RMS(residual) over N samples
        sigma_Poisson  = sqrt(S_h * li_2(x_anchor))   (HL random-residual)
        sigma_eff/sqrt(x)
        Half-normal moment ratios median/sigma_eff, mean/sigma_eff
        KS p-value of |r|/sigma_eff vs half-Gaussian
        KS p-value of |r|/sigma_Poisson vs half-Gaussian

Output:
  slot2_data.csv      -- per-(anchor, x_j, h) raw residual
  slot2_summary.csv   -- per-(anchor, h) aggregate stats and KS p-values

The script is the sole computational artefact; results.md narrates.

Falsification statement: if the median KS p-value across all
(anchor, h) cells with N=30 is < 0.05 (i.e., uniformly rejects
half-Gaussian), then HL residuals are NOT half-Gaussian-shape and
the Thread-7 (P3) analogue does not transpose to the h-axis.
Conversely a stable sigma_eff/sigma_Poisson ratio across decades
with median KS p > 0.1 supports a Thread-7-shape A-grade target
for slot 4 (named-exponent corollary for pi_h).

Edges cited:
  - E1.5 (information density of pi)
  - S195 (variance formula machinery, transposed)
  - S224 (Correlation Dichotomy reference shape)
  - S418 slot-1 baseline (structural dichotomy this slot extends)
  - Thread 7 multi_sample.py (script structure template)
"""

from __future__ import annotations

import csv
import math
import os
import sys
import time
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------------
# Constants (matching slot 1)
# ---------------------------------------------------------------------------
C_2 = 0.6601618158468695739278121100145557784326233628

H_VALUES = [2, 6, 14, 22, 30, 42, 90, 198]
ANCHORS_LOG10 = [6, 7, 8]
N_SAMPLES = 30
WINDOW_LOG = 0.25 * math.log(10.0)  # quarter-decade

# ---------------------------------------------------------------------------
# Sieve + Hardy-Littlewood primitives
# ---------------------------------------------------------------------------


def sieve_eratosthenes(N: int) -> bytearray:
    is_p = bytearray(b"\x01") * (N + 1)
    is_p[0] = 0
    is_p[1] = 0
    bound = int(N ** 0.5) + 1
    for i in range(2, bound):
        if is_p[i]:
            start = i * i
            is_p[start : N + 1 : i] = b"\x00" * (((N - start) // i) + 1)
    return is_p


def li2(x: float) -> float:
    from scipy.integrate import quad
    v, _ = quad(lambda t: 1.0 / math.log(t) ** 2, 2.0, float(x), limit=400)
    return v


def singular_series(h: int) -> float:
    if h % 2 != 0:
        return 0.0
    n = h
    while n % 2 == 0:
        n //= 2
    factor = 1.0
    p = 3
    while p * p <= n:
        if n % p == 0:
            factor *= (p - 1) / (p - 2)
            while n % p == 0:
                n //= p
        p += 2
    if n > 1:
        factor *= (n - 1) / (n - 2)
    return 2.0 * C_2 * factor


# ---------------------------------------------------------------------------
# Sample design
# ---------------------------------------------------------------------------


def sample_xs(x_anchor: int, N: int) -> list[int]:
    """Log-uniform N samples in [x_anchor, x_anchor * exp(WINDOW_LOG)]."""
    log_lo = math.log(float(x_anchor))
    log_hi = log_lo + WINDOW_LOG
    raw = [int(round(math.exp(log_lo + j * (log_hi - log_lo) / (N - 1)))) for j in range(N)]
    # Dedup while preserving order; pad if dedup shrinks the list.
    seen = set()
    out = []
    for x in raw:
        if x not in seen:
            seen.add(x)
            out.append(x)
    while len(out) < N:
        # Fill with x_anchor + linear gap; rare for these decade scales.
        last = out[-1]
        out.append(last + max(1, x_anchor // 1000))
    return out[:N]


# ---------------------------------------------------------------------------
# Per-anchor evaluator
# ---------------------------------------------------------------------------


def evaluate_anchor(anchor_log10: int, N: int = N_SAMPLES) -> list[dict]:
    x_anchor = 10 ** anchor_log10
    xs = sample_xs(x_anchor, N)
    x_max = xs[-1]
    H_max = max(H_VALUES)

    print(f"\n[anchor x = 10^{anchor_log10}]  N={len(xs)}  x_max={x_max:,}  H_max={H_max}")

    t0 = time.perf_counter()
    is_p = sieve_eratosthenes(x_max + H_max)
    t_sieve = time.perf_counter() - t0
    print(f"  sieve [1, {x_max + H_max:,}] in {t_sieve:.2f}s")

    is_p_arr = np.frombuffer(is_p, dtype=np.uint8)

    t0 = time.perf_counter()
    li2_xs = [li2(x) for x in xs]
    t_li2 = time.perf_counter() - t0
    print(f"  li_2 evaluations: {t_li2:.2f}s")

    rows: list[dict] = []
    t_count = 0.0
    for h in H_VALUES:
        t1 = time.perf_counter()
        S_h = singular_series(h)
        # pair[i] = 1 iff i prime AND i+h prime, for i in [0, x_max].
        pair = is_p_arr[: x_max + 1] & is_p_arr[h : x_max + 1 + h]
        # cumcount[i] = #{j <= i : j prime AND j+h prime} = pi_h evaluated at i.
        cumcount = np.cumsum(pair, dtype=np.int64)
        for x_j, li2_xj in zip(xs, li2_xs):
            pi_h_xj = int(cumcount[x_j])
            HL = S_h * li2_xj
            r = pi_h_xj - HL
            sigma_pois = math.sqrt(S_h * li2_xj) if S_h > 0 else 0.0
            rows.append({
                "anchor_log10": anchor_log10,
                "x_anchor": x_anchor,
                "x": x_j,
                "h": h,
                "S_h": S_h,
                "pi_h": pi_h_xj,
                "HL": HL,
                "residual": r,
                "abs_residual": abs(r),
                "sigma_pred_poisson": sigma_pois,
                "sqrt_x": math.sqrt(x_j),
                "li2_x": li2_xj,
            })
        del pair, cumcount  # free memory before next h
        t_count += time.perf_counter() - t1
    print(f"  pair-count over {len(H_VALUES)} h-values: {t_count:.2f}s")
    return rows


# ---------------------------------------------------------------------------
# Statistics: KS test against half-Gaussian
# ---------------------------------------------------------------------------


def half_normal_cdf(z: float) -> float:
    return math.erf(z / math.sqrt(2)) if z >= 0 else 0.0


def ks_stat(ratios: list[float]) -> float:
    n = len(ratios)
    if n == 0:
        return float("nan")
    s = sorted(ratios)
    D = 0.0
    for i, z in enumerate(s):
        F_emp_lo = i / n
        F_emp_hi = (i + 1) / n
        F_th = half_normal_cdf(z)
        D = max(D, abs(F_emp_lo - F_th), abs(F_emp_hi - F_th))
    return D


def ks_pvalue(D: float, n: int) -> float:
    if n < 2:
        return float("nan")
    en = math.sqrt(n) * D + 0.12 * D + 0.11 * D / math.sqrt(n)
    s = 0.0
    for j in range(1, 50):
        term = 2 * (-1) ** (j - 1) * math.exp(-2 * j * j * en * en)
        s += term
        if abs(term) < 1e-12:
            break
    return min(1.0, max(0.0, s))


def aggregate(rows: list[dict]) -> list[dict]:
    """Per-(anchor, h) aggregate: KS across N=30 x-samples within the window."""
    groups: dict[tuple[int, int], list[dict]] = {}
    for r in rows:
        key = (r["anchor_log10"], r["h"])
        groups.setdefault(key, []).append(r)

    out = []
    for (anc, h), g in sorted(groups.items()):
        residuals = [row["residual"] for row in g]
        abs_res = [row["abs_residual"] for row in g]
        N = len(g)

        sigma_eff = math.sqrt(sum(r * r for r in residuals) / N)
        median_abs = sorted(abs_res)[N // 2]
        mean_abs = sum(abs_res) / N
        sig_pois = sorted([row["sigma_pred_poisson"] for row in g])[N // 2]
        sqrt_x_med = sorted([row["sqrt_x"] for row in g])[N // 2]
        # Within-window range diagnostic: if max|r|-min|r| << sigma_eff, the
        # residual is near-constant and KS shape-test is uninformative.
        rmin = min(residuals)
        rmax = max(residuals)
        amax = max(abs_res)
        amin = min(abs_res)

        ratios_eff = [a / sigma_eff for a in abs_res] if sigma_eff > 0 else []
        ratios_pois = [a / sig_pois for a in abs_res] if sig_pois > 0 else []
        D_eff = ks_stat(ratios_eff)
        p_eff = ks_pvalue(D_eff, len(ratios_eff))
        D_pois = ks_stat(ratios_pois)
        p_pois = ks_pvalue(D_pois, len(ratios_pois))

        out.append({
            "anchor_log10": anc,
            "h": h,
            "S_h": g[0]["S_h"],
            "N": N,
            "sigma_eff": sigma_eff,
            "sigma_pred_poisson": sig_pois,
            "sigma_eff_over_poisson": sigma_eff / sig_pois if sig_pois > 0 else float("nan"),
            "sigma_eff_over_sqrt_x": sigma_eff / sqrt_x_med if sqrt_x_med > 0 else float("nan"),
            "median_abs": median_abs,
            "mean_abs": mean_abs,
            "median_over_sigma_eff": median_abs / sigma_eff if sigma_eff > 0 else float("nan"),
            "mean_over_sigma_eff": mean_abs / sigma_eff if sigma_eff > 0 else float("nan"),
            "abs_range_over_sigma_eff": (amax - amin) / sigma_eff if sigma_eff > 0 else float("nan"),
            "signed_range": rmax - rmin,
            "ks_D_eff": D_eff,
            "ks_p_eff": p_eff,
            "ks_D_poisson": D_pois,
            "ks_p_poisson": p_pois,
        })
    return out


def aggregate_cross_h(rows: list[dict]) -> list[dict]:
    """Per-(anchor, x_j) aggregate: KS across the K=8 h-values at each window
    point. Residuals normalized by sigma_pred_poisson (h-dependent) so the
    pooled sample is dimensionless. N per cell = 8."""
    groups: dict[tuple[int, int], list[dict]] = {}
    for r in rows:
        key = (r["anchor_log10"], r["x"])
        groups.setdefault(key, []).append(r)

    out = []
    for (anc, x_j), g in sorted(groups.items()):
        N = len(g)
        # Normalized residuals: |r| / sigma_pred_poisson (h-dependent)
        norm_abs = [
            row["abs_residual"] / row["sigma_pred_poisson"]
            if row["sigma_pred_poisson"] > 0 else 0.0
            for row in g
        ]
        norm_signed = [
            row["residual"] / row["sigma_pred_poisson"]
            if row["sigma_pred_poisson"] > 0 else 0.0
            for row in g
        ]
        sigma_normalized = math.sqrt(sum(r * r for r in norm_signed) / N)
        # KS of |r|/sigma_pred_poisson directly (not re-normalized)
        D_pois = ks_stat(norm_abs)
        p_pois = ks_pvalue(D_pois, N)
        # KS of |r|/sigma_normalized (re-normalized within the cross-h ensemble)
        if sigma_normalized > 0:
            ratios_eff = [a / sigma_normalized for a in norm_abs]
        else:
            ratios_eff = []
        D_eff = ks_stat(ratios_eff)
        p_eff = ks_pvalue(D_eff, N)

        out.append({
            "anchor_log10": anc,
            "x": x_j,
            "N_h": N,
            "sigma_normalized": sigma_normalized,
            "mean_norm_abs": sum(norm_abs) / N,
            "median_norm_abs": sorted(norm_abs)[N // 2],
            "ks_D_h_pois": D_pois,
            "ks_p_h_pois": p_pois,
            "ks_D_h_eff": D_eff,
            "ks_p_h_eff": p_eff,
        })
    return out


def aggregate_pooled_per_anchor(rows: list[dict]) -> list[dict]:
    """Pool all (x_j, h) per anchor into one ensemble of N = N_x * N_h samples.
    Normalize each by sigma_pred_poisson, then KS against half-Gaussian."""
    groups: dict[int, list[dict]] = {}
    for r in rows:
        groups.setdefault(r["anchor_log10"], []).append(r)

    out = []
    for anc, g in sorted(groups.items()):
        N = len(g)
        norm_abs = [
            row["abs_residual"] / row["sigma_pred_poisson"]
            if row["sigma_pred_poisson"] > 0 else 0.0
            for row in g
        ]
        norm_signed = [
            row["residual"] / row["sigma_pred_poisson"]
            if row["sigma_pred_poisson"] > 0 else 0.0
            for row in g
        ]
        sigma_normalized = math.sqrt(sum(r * r for r in norm_signed) / N)
        D_pois = ks_stat(norm_abs)
        p_pois = ks_pvalue(D_pois, N)
        if sigma_normalized > 0:
            ratios_eff = [a / sigma_normalized for a in norm_abs]
        else:
            ratios_eff = []
        D_eff = ks_stat(ratios_eff)
        p_eff = ks_pvalue(D_eff, N)

        out.append({
            "anchor_log10": anc,
            "N_pooled": N,
            "sigma_normalized_eff": sigma_normalized,
            "mean_norm_abs": sum(norm_abs) / N,
            "median_norm_abs": sorted(norm_abs)[N // 2],
            "ks_D_pooled_pois": D_pois,
            "ks_p_pooled_pois": p_pois,
            "ks_D_pooled_eff": D_eff,
            "ks_p_pooled_eff": p_eff,
        })
    return out


# ---------------------------------------------------------------------------
# Main driver
# ---------------------------------------------------------------------------


def main():
    out_dir = Path(__file__).parent

    all_rows: list[dict] = []
    t_start = time.perf_counter()
    for anc in ANCHORS_LOG10:
        rows = evaluate_anchor(anc, N=N_SAMPLES)
        all_rows.extend(rows)
    t_total = time.perf_counter() - t_start
    print(f"\nTotal elapsed: {t_total:.1f}s, rows: {len(all_rows)}")

    # Per-row CSV
    raw_path = out_dir / "slot2_data.csv"
    with raw_path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow([
            "anchor_log10", "x_anchor", "x", "h", "S_h", "pi_h", "HL",
            "residual", "abs_residual", "sigma_pred_poisson", "sqrt_x", "li2_x",
        ])
        for r in all_rows:
            w.writerow([
                r["anchor_log10"], r["x_anchor"], r["x"], r["h"],
                f"{r['S_h']:.6f}", r["pi_h"], f"{r['HL']:.6f}",
                f"{r['residual']:.6f}", f"{r['abs_residual']:.6f}",
                f"{r['sigma_pred_poisson']:.6f}", f"{r['sqrt_x']:.6f}",
                f"{r['li2_x']:.6f}",
            ])
    print(f"Raw rows -> {raw_path}")

    # Aggregate summary
    summary = aggregate(all_rows)
    sum_path = out_dir / "slot2_summary.csv"
    with sum_path.open("w", newline="") as f:
        w = csv.writer(f)
        keys = list(summary[0].keys())
        w.writerow(keys)
        for r in summary:
            w.writerow([
                r[k] if isinstance(r[k], int) else f"{r[k]:.6f}"
                for k in keys
            ])
    print(f"Summary -> {sum_path}")

    # Print summary table
    print("\n" + "=" * 130)
    print(f"{'anc':>5s} {'h':>4s} {'S_h':>6s} {'sig_eff':>10s} {'sig_pois':>10s} "
          f"{'eff/pois':>9s} {'eff/sqrt_x':>11s} {'med/sig_eff':>12s} {'mean/sig_eff':>13s} "
          f"{'KSp_eff':>9s} {'KSp_pois':>9s}")
    print("-" * 130)
    for r in summary:
        print(f"{r['anchor_log10']:>5d} {r['h']:>4d} {r['S_h']:>6.3f} "
              f"{r['sigma_eff']:>10.3f} {r['sigma_pred_poisson']:>10.3f} "
              f"{r['sigma_eff_over_poisson']:>9.4f} "
              f"{r['sigma_eff_over_sqrt_x']:>11.4f} "
              f"{r['median_over_sigma_eff']:>12.4f} "
              f"{r['mean_over_sigma_eff']:>13.4f} "
              f"{r['ks_p_eff']:>9.4f} {r['ks_p_poisson']:>9.4f}")
    print("=" * 130)

    # Half-Gaussian theoretical reference: median/sigma=0.6745, mean/sigma=0.7979
    # Aggregate KS-p across all (anchor, h) cells.
    ks_p_eff_all = [r["ks_p_eff"] for r in summary]
    ks_p_pois_all = [r["ks_p_poisson"] for r in summary]
    median_ks_eff = sorted(ks_p_eff_all)[len(ks_p_eff_all) // 2]
    median_ks_pois = sorted(ks_p_pois_all)[len(ks_p_pois_all) // 2]
    print(f"\n  Theoretical half-Gaussian moments: median/sigma=0.6745, mean/sigma=0.7979")
    print(f"  median KS p-value across {len(ks_p_eff_all)} (anchor, h) cells: "
          f"eff={median_ks_eff:.4f}, poisson={median_ks_pois:.4f}")
    print(f"  cells with KS p_eff > 0.1:    {sum(1 for p in ks_p_eff_all if p > 0.1)}/{len(ks_p_eff_all)}")
    print(f"  cells with KS p_eff < 0.05:   {sum(1 for p in ks_p_eff_all if p < 0.05)}/{len(ks_p_eff_all)}")

    # Within-window range diagnostic: max|r|-min|r| over sigma_eff. If << 1
    # the residual is near-constant within the quarter-decade window.
    print("\n  Within-window range diagnostic: (max|r|-min|r|)/sigma_eff")
    for r in summary:
        print(f"    anc=10^{r['anchor_log10']}, h={r['h']:>3d}: "
              f"abs_range/sigma_eff = {r['abs_range_over_sigma_eff']:.3f}, "
              f"signed_range = {r['signed_range']:.1f}")

    # ========== Cross-h aggregate ==========
    cross_h = aggregate_cross_h(all_rows)
    crh_path = out_dir / "slot2_cross_h.csv"
    with crh_path.open("w", newline="") as f:
        w = csv.writer(f)
        keys = list(cross_h[0].keys())
        w.writerow(keys)
        for r in cross_h:
            w.writerow([
                r[k] if isinstance(r[k], int) else f"{r[k]:.6f}"
                for k in keys
            ])
    print(f"\nCross-h aggregate -> {crh_path}")

    # Aggregate across anchors: median KS p-value across all cross-h cells per anchor
    by_anc: dict[int, list[dict]] = {}
    for r in cross_h:
        by_anc.setdefault(r["anchor_log10"], []).append(r)
    print("\n  Cross-h KS summary (N_h=8 per cell, 30 cells per anchor):")
    print(f"  {'anc':>5s} {'N_cells':>9s} {'med KSp_eff':>12s} {'med KSp_pois':>12s} "
          f"{'cells p_eff>0.1':>17s} {'med sigma_norm':>15s}")
    for anc, cells in sorted(by_anc.items()):
        kspe = sorted([c["ks_p_h_eff"] for c in cells])
        ksp = sorted([c["ks_p_h_pois"] for c in cells])
        signorm = sorted([c["sigma_normalized"] for c in cells])
        med_kspe = kspe[len(kspe) // 2]
        med_ksp = ksp[len(ksp) // 2]
        med_sn = signorm[len(signorm) // 2]
        n_pass = sum(1 for p in kspe if p > 0.1)
        print(f"  {anc:>5d} {len(cells):>9d} {med_kspe:>12.4f} {med_ksp:>12.4f} "
              f"{n_pass:>17d} {med_sn:>15.4f}")

    # ========== Pooled per-anchor aggregate ==========
    pooled = aggregate_pooled_per_anchor(all_rows)
    pool_path = out_dir / "slot2_pooled.csv"
    with pool_path.open("w", newline="") as f:
        w = csv.writer(f)
        keys = list(pooled[0].keys())
        w.writerow(keys)
        for r in pooled:
            w.writerow([
                r[k] if isinstance(r[k], int) else f"{r[k]:.6f}"
                for k in keys
            ])
    print(f"\nPooled per-anchor -> {pool_path}")

    print("\n  Pooled per-anchor KS summary (N=240 per anchor):")
    print(f"  {'anc':>5s} {'N':>5s} {'sig_norm':>10s} {'mean_norm':>10s} {'med_norm':>10s} "
          f"{'KSp_eff':>9s} {'KSp_pois':>9s}")
    for r in pooled:
        print(f"  {r['anchor_log10']:>5d} {r['N_pooled']:>5d} "
              f"{r['sigma_normalized_eff']:>10.4f} {r['mean_norm_abs']:>10.4f} "
              f"{r['median_norm_abs']:>10.4f} "
              f"{r['ks_p_pooled_eff']:>9.4f} {r['ks_p_pooled_pois']:>9.4f}")


if __name__ == "__main__":
    main()

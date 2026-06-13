"""
slot2_cross_x.py — Thread 9 (P4) commit slot 2

Cross-x HL residual distribution analysis at fixed admissible tuple H,
with disjoint narrow polylog windows. The cross-x analogue of Thread 8
slot 2's cross-h ensemble (S419) on the new (x-axis at fixed H, varying
x) cell.

For each anchor x_anchor in {10^6, 10^7, 10^8}:
  - Sieve [x_anchor, x_anchor * 10^0.25 + h_max + w] once (one pass).
  - For each H in {(0,2), (0,2,6), (0,2,6,8)}:
      - Build pair_array[i] = AND_{h in H} is_prime[i+h].
      - Sample N = 200 disjoint x_i across [x_anchor, x_anchor * 10^0.25]
        with spacing >= w + h_max so windows are non-overlapping AND
        independent up to short-range correlation.
      - For each x_i, compute exact pi_H(x_i; w) via slice-sum on
        pair_array, plus the HL prediction
        HL_H(x_i; w) = C(H) * (li_k(x_i + w) - li_k(x_i)),
        approximating li_k via midpoint w / log^k(x_i + w/2).
      - Record signed residual r = pi_H - HL_H.
  - Aggregate per-(anchor, H, w_regime):
      sigma_eff      = RMS(residual) over N samples.
      sigma_pois     = sqrt(HL_H(x_anchor; w))     (HL random-residual).
      sigma_eff/sigma_pois suppression factor.
      Half-normal moment ratios median/sigma_eff, mean/sigma_eff.
      KS p-value of |r|/sigma_eff vs half-Gaussian.
      KS p-value of |r|/sigma_pois vs half-Gaussian.

Design choices that distinguish slot 2 from Thread 8 slot 2:
  - Disjoint windows (x_i >= x_{i-1} + w + h_max). Critical: this is the
    structural reason slot 2 cross-x works where Thread 8 cross-x failed.
    Thread 8 used a CUMULATIVE pi_h that drifts smoothly within quarter
    decade. Thread 9 uses WINDOWED pi_H whose disjoint-x samples are
    near-iid.
  - Two w-regimes per anchor:
      narrow:   w = log^2(x_anchor)  -- slot 1 baseline; HL_2 ~ 1.3 (discrete).
      wide:     w = 12 * log^2(x_anchor)  -- still polylog, HL_2 ~ 15 (continuous).
    The wide regime is a polylog scale chosen so HL_k=2 is approximately
    decade-stable (= S_2 * 12 = 15.85), enabling clean decade-stability
    tests for the suppression factor.
  - 3 H-tuples * 3 anchors * 2 w-regimes = 18 (anchor, H, w) cells, each
    with N = 200 disjoint x_i.

Output:
  slot2_data.csv       -- per-(anchor, H, w_regime, x) raw residual.
  slot2_summary.csv    -- per-(anchor, H, w_regime) aggregate.
  slot2_pooled.csv     -- pooled per anchor (for decade-stability check).
  slot2_run.log        -- run log.

Falsifiable hypotheses (slot 2):

  H4 (cross-x ensemble passes half-Gaussian KS): if the cross-x ensemble
      at fixed H is iid-like (disjoint windows), then the KS test of
      |r_i|/sigma_eff against half-Gaussian should pass at p > 0.1 for
      most cells. Falsified iff median KS_p_eff < 0.05 across all cells.

  H5 (sigma_eff/sigma_pois suppression factor exists, analogous to
      Thread 7 F_GUE = 0.755 and Thread 8 [0.36, 0.70]): empirical
      sigma_eff/sigma_pois should be < 1 (suppression below Poisson
      iid-Bernoulli prediction) AND > 0.3 (some residual variance).
      Falsified iff the ratio is OUT of [0.3, 1.05] in most cells, OR
      strongly decade-trending in a way that contradicts a fixed F_HL.

  H6 (decade-stability of suppression): the wide-regime ratio should be
      decade-stable (since HL_2 is decade-constant). Decade trend
      magnitude should be < 0.15 (analogous to Thread 7's flat F_GUE).

Edges cited:
  - E1.5 (information density of pi)
  - S195 / S224 (Thread 5 Correlation Dichotomy reference shape)
  - S418 / S419 / S421 (Thread 8 P2 EXACT/APPROX dichotomy and HL
    cross-h ensemble template — slot 2 uses cross-x analogue)
  - Hardy-Littlewood 1923 k-tuple conjecture (provides C(H))
  - S429 slot-1 baseline (provides HL_H(x; w) match at 6 cells; slot 2
    extends to N=200 cross-x distribution)
"""

from __future__ import annotations

import csv
import math
import os
import random
import time
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
ANCHORS_LOG10 = [6, 7, 8]
N_DISJOINT = 200
H_TUPLES = [
    ("k2_twin", (0, 2)),
    ("k3_admiss", (0, 2, 6)),
    ("k4_admiss", (0, 2, 6, 8)),
]
W_REGIMES = [
    ("narrow", 1),    # w = 1 * log^2(x)
    ("wide", 12),     # w = 12 * log^2(x), HL_2 ~ 15.85 (decade-stable)
]
WINDOW_LOG = 0.25 * math.log(10.0)  # quarter-decade


# ---------------------------------------------------------------------------
# Sieve + Hardy-Littlewood primitives
# ---------------------------------------------------------------------------

def sieve_eratosthenes(N: int) -> bytearray:
    """Sieve of Eratosthenes; returns bytearray of length N+1 with
    is_p[i] = 1 iff i prime."""
    is_p = bytearray(b"\x01") * (N + 1)
    is_p[0] = 0
    if N >= 1:
        is_p[1] = 0
    bound = int(N ** 0.5) + 1
    for i in range(2, bound):
        if is_p[i]:
            start = i * i
            is_p[start : N + 1 : i] = b"\x00" * (((N - start) // i) + 1)
    return is_p


def singular_series(H, P_max: int = 10000) -> float:
    """C(H) = product_{p prime} (1 - nu_H(p)/p) / (1 - 1/p)^k where
    nu_H(p) = #{distinct h mod p}. Truncates at P_max."""
    H = list(H)
    k = len(H)
    is_p = sieve_eratosthenes(P_max)
    val = 1.0
    for p in range(2, P_max + 1):
        if not is_p[p]:
            continue
        residues = set(h % p for h in H)
        nu = len(residues)
        val *= (1.0 - nu / p) / ((1.0 - 1.0 / p) ** k)
    return val


def hl_window(C_H: float, x: int, w: int, k: int) -> float:
    """Hardy-Littlewood prediction:
        pi_H(x; w) ~ C(H) * (li_k(x+w) - li_k(x))
                   ~ C(H) * w / log^k(x + w/2)
    Midpoint approximation; for narrow w this is accurate to O(w/x · log^{-k}).
    """
    x_mid = x + w / 2.0
    if x_mid < 2:
        return 0.0
    return C_H * w / (math.log(x_mid) ** k)


# ---------------------------------------------------------------------------
# Sample design: disjoint x_i across quarter-decade window
# ---------------------------------------------------------------------------

def sample_disjoint_xs(x_anchor: int, n: int, window_total: int, seed: int) -> list[int]:
    """Sample n disjoint window starts in [x_anchor, x_anchor * 10^0.25]
    with each window of width window_total = w + h_max non-overlapping.
    Uses uniform stratified sampling across the available range."""
    x_lo = x_anchor
    x_hi = int(x_anchor * math.exp(WINDOW_LOG))
    span = x_hi - x_lo
    # Maximum number of disjoint windows that fit
    max_disjoint = span // window_total
    if max_disjoint < n:
        # Should not happen for our anchor / w settings, but handle gracefully.
        n = int(max_disjoint)
    if n <= 0:
        return [x_lo]
    rng = random.Random(seed)
    # Stratified sampling: divide span into n bins, sample one window-start per bin
    bin_width = span / n
    xs = []
    for j in range(n):
        bin_lo = x_lo + int(j * bin_width)
        bin_hi = x_lo + int((j + 1) * bin_width) - window_total
        if bin_hi <= bin_lo:
            xs.append(bin_lo)
        else:
            xs.append(rng.randint(bin_lo, bin_hi))
    # Sort + verify disjoint
    xs.sort()
    cleaned = [xs[0]]
    for xi in xs[1:]:
        if xi >= cleaned[-1] + window_total:
            cleaned.append(xi)
    return cleaned


# ---------------------------------------------------------------------------
# Per-anchor evaluator
# ---------------------------------------------------------------------------

def evaluate_anchor(anchor_log10: int, log_fn) -> list[dict]:
    x_anchor = 10 ** anchor_log10
    x_max = int(x_anchor * math.exp(WINDOW_LOG))  # quarter-decade upper

    log_x_anc = math.log(x_anchor)
    # Determine w for each regime, and compute global h_max + w_max for sieve range
    w_by_regime = {label: max(2, int(c * log_x_anc ** 2)) for label, c in W_REGIMES}
    h_max_global = max(max(H) for _, H in H_TUPLES)
    w_max = max(w_by_regime.values())
    sieve_upper = x_max + w_max + h_max_global + 4

    log_fn(f"\n[anchor x = 10^{anchor_log10}]  x_anchor={x_anchor:,}  x_max={x_max:,}")
    log_fn(f"  w regimes: {w_by_regime}")
    log_fn(f"  sieve [1, {sieve_upper:,}]")

    t0 = time.perf_counter()
    is_p = sieve_eratosthenes(sieve_upper)
    is_p_arr = np.frombuffer(is_p, dtype=np.uint8)
    t_sieve = time.perf_counter() - t0
    log_fn(f"  sieved in {t_sieve:.2f}s")

    rows: list[dict] = []
    for k_label, H in H_TUPLES:
        k = len(H)
        C_H = singular_series(H)
        # Build pair_array[i] = AND_{h in H} is_p[i + h]
        # Length = sieve_upper - h_max + 1 (so we can index up to sieve_upper - h_max)
        L = sieve_upper - max(H)
        t1 = time.perf_counter()
        pair = is_p_arr[0:L].copy()
        for h in H[1:]:
            pair &= is_p_arr[h:h + L]
        t_pair = time.perf_counter() - t1
        log_fn(f"  H={H} (k={k}, C_H={C_H:.4f}): pair_array built in {t_pair:.2f}s, length={L}")

        for w_label, w_factor in W_REGIMES:
            w = w_by_regime[w_label]
            window_total = w + max(H) + 1
            seed = (anchor_log10 * 100 + sum(H) * 10 + w_factor) & 0xFFFFFFFF
            xs = sample_disjoint_xs(x_anchor, N_DISJOINT, window_total, seed)
            n_drawn = len(xs)
            # Compute pi_H(x_i; w) for each x_i: slice-sum of pair_array over [x_i, x_i + w]
            counts = np.empty(n_drawn, dtype=np.int64)
            for i, xi in enumerate(xs):
                # Window is [xi, xi+w] (length w+1), pi_H = #{j in [xi, xi+w] : pair[j]=1}
                if xi + w + 1 > L:
                    counts[i] = -1  # boundary issue; drop later
                else:
                    counts[i] = int(pair[xi:xi + w + 1].sum())
            valid_mask = counts >= 0
            xs_valid = [xi for xi, ok in zip(xs, valid_mask) if ok]
            counts_valid = counts[valid_mask]

            for xi, c in zip(xs_valid, counts_valid):
                HL = hl_window(C_H, xi, w, k)
                r = float(c) - HL
                rows.append({
                    "anchor_log10": anchor_log10,
                    "x_anchor": x_anchor,
                    "k_label": k_label,
                    "H": str(H),
                    "k": k,
                    "C_H": C_H,
                    "w_label": w_label,
                    "w": w,
                    "x": xi,
                    "pi_H": int(c),
                    "HL": HL,
                    "residual": r,
                    "abs_residual": abs(r),
                    "log_x": math.log(xi),
                    "sqrt_HL": math.sqrt(HL) if HL > 0 else 0.0,
                })
        # free pair memory
        del pair

    return rows


# ---------------------------------------------------------------------------
# KS test against half-Gaussian
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


# ---------------------------------------------------------------------------
# Aggregation
# ---------------------------------------------------------------------------

def aggregate_per_cell(rows: list[dict]) -> list[dict]:
    """Aggregate per-(anchor, k_label, w_label) cell. KS over N disjoint x_i."""
    groups: dict[tuple, list[dict]] = {}
    for r in rows:
        key = (r["anchor_log10"], r["k_label"], r["w_label"])
        groups.setdefault(key, []).append(r)
    out = []
    for (anc, k_lab, w_lab), g in sorted(groups.items()):
        N = len(g)
        residuals = [r["residual"] for r in g]
        abs_res = [r["abs_residual"] for r in g]
        sigma_eff = math.sqrt(sum(r * r for r in residuals) / N)
        # sigma_pois = sqrt(mean HL) — random-residual prediction
        mean_HL = sum(r["HL"] for r in g) / N
        sigma_pois = math.sqrt(max(mean_HL, 1e-12))
        median_abs = sorted(abs_res)[N // 2]
        mean_abs = sum(abs_res) / N
        # Counts distribution diagnostic
        pis = [r["pi_H"] for r in g]
        n_zero = sum(1 for p in pis if p == 0)
        max_count = max(pis)
        mean_count = sum(pis) / N
        # KS
        ratios_eff = [a / sigma_eff for a in abs_res] if sigma_eff > 0 else []
        ratios_pois = [a / sigma_pois for a in abs_res] if sigma_pois > 0 else []
        D_eff = ks_stat(ratios_eff)
        p_eff = ks_pvalue(D_eff, len(ratios_eff))
        D_pois = ks_stat(ratios_pois)
        p_pois = ks_pvalue(D_pois, len(ratios_pois))
        out.append({
            "anchor_log10": anc,
            "k_label": k_lab,
            "k": g[0]["k"],
            "C_H": g[0]["C_H"],
            "w_label": w_lab,
            "w": g[0]["w"],
            "N": N,
            "mean_HL": mean_HL,
            "mean_count": mean_count,
            "max_count": max_count,
            "n_zero": n_zero,
            "sigma_eff": sigma_eff,
            "sigma_pois": sigma_pois,
            "sigma_eff_over_pois": sigma_eff / sigma_pois if sigma_pois > 0 else float("nan"),
            "median_abs": median_abs,
            "mean_abs": mean_abs,
            "median_over_sigma_eff": median_abs / sigma_eff if sigma_eff > 0 else float("nan"),
            "mean_over_sigma_eff": mean_abs / sigma_eff if sigma_eff > 0 else float("nan"),
            "ks_D_eff": D_eff,
            "ks_p_eff": p_eff,
            "ks_D_pois": D_pois,
            "ks_p_pois": p_pois,
        })
    return out


def aggregate_pooled_per_anchor(rows: list[dict]) -> list[dict]:
    """Pool across all (k, w) at each anchor — for decade-stability check."""
    groups: dict[int, list[dict]] = {}
    for r in rows:
        groups.setdefault(r["anchor_log10"], []).append(r)
    out = []
    for anc, g in sorted(groups.items()):
        N = len(g)
        # Normalize residual by sqrt(HL) to make dimensionless
        norm_abs = [r["abs_residual"] / math.sqrt(max(r["HL"], 1e-12)) for r in g]
        norm_signed = [r["residual"] / math.sqrt(max(r["HL"], 1e-12)) for r in g]
        sigma_norm = math.sqrt(sum(r * r for r in norm_signed) / N)
        D_pois = ks_stat(norm_abs)
        p_pois = ks_pvalue(D_pois, N)
        if sigma_norm > 0:
            ratios_eff = [a / sigma_norm for a in norm_abs]
        else:
            ratios_eff = []
        D_eff = ks_stat(ratios_eff)
        p_eff = ks_pvalue(D_eff, N)
        out.append({
            "anchor_log10": anc,
            "N_pooled": N,
            "sigma_norm_eff": sigma_norm,
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
    log_path = out_dir / "slot2_run.log"
    log_lines: list[str] = []

    def log_fn(msg: str):
        print(msg)
        log_lines.append(msg)

    log_fn("=" * 78)
    log_fn("Thread 9 / P4 slot 2: cross-x HL residual distribution at fixed H")
    log_fn("=" * 78)
    log_fn(f"Anchors: {ANCHORS_LOG10}, N_disjoint per cell: {N_DISJOINT}")
    log_fn(f"H tuples: {[H for _, H in H_TUPLES]}")
    log_fn(f"w regimes: {W_REGIMES}")
    log_fn("")

    all_rows: list[dict] = []
    t_start = time.perf_counter()
    for anc in ANCHORS_LOG10:
        rs = evaluate_anchor(anc, log_fn)
        all_rows.extend(rs)
    t_total = time.perf_counter() - t_start
    log_fn(f"\nTotal elapsed: {t_total:.1f}s, raw rows: {len(all_rows)}")

    # Raw data CSV
    raw_path = out_dir / "slot2_data.csv"
    if all_rows:
        with raw_path.open("w", newline="") as f:
            keys = list(all_rows[0].keys())
            w = csv.DictWriter(f, fieldnames=keys)
            w.writeheader()
            for r in all_rows:
                w.writerow({k: (f"{r[k]:.6f}" if isinstance(r[k], float) else r[k]) for k in keys})
    log_fn(f"\nRaw -> {raw_path}  ({len(all_rows)} rows)")

    # Per-cell summary
    summary = aggregate_per_cell(all_rows)
    sum_path = out_dir / "slot2_summary.csv"
    if summary:
        with sum_path.open("w", newline="") as f:
            keys = list(summary[0].keys())
            w = csv.DictWriter(f, fieldnames=keys)
            w.writeheader()
            for r in summary:
                w.writerow({k: (f"{r[k]:.6f}" if isinstance(r[k], float) else r[k]) for k in keys})
    log_fn(f"Summary -> {sum_path}  ({len(summary)} cells)")

    # Pooled per anchor
    pooled = aggregate_pooled_per_anchor(all_rows)
    pool_path = out_dir / "slot2_pooled.csv"
    if pooled:
        with pool_path.open("w", newline="") as f:
            keys = list(pooled[0].keys())
            w = csv.DictWriter(f, fieldnames=keys)
            w.writeheader()
            for r in pooled:
                w.writerow({k: (f"{r[k]:.6f}" if isinstance(r[k], float) else r[k]) for k in keys})
    log_fn(f"Pooled -> {pool_path}  ({len(pooled)} anchor-rows)")

    # Print headline summary table
    log_fn("\n" + "=" * 130)
    log_fn(f"{'anc':>4s} {'k':>3s} {'w_lab':>7s} {'w':>6s} {'mean_HL':>9s} {'sig_eff':>9s} "
           f"{'sig_pois':>9s} {'eff/pois':>9s} {'med/sig_eff':>12s} {'mean/sig_eff':>13s} "
           f"{'KSp_eff':>9s} {'KSp_pois':>9s}")
    log_fn("-" * 130)
    for r in summary:
        log_fn(f"{r['anchor_log10']:>4d} {r['k']:>3d} {r['w_label']:>7s} {r['w']:>6d} "
               f"{r['mean_HL']:>9.3f} {r['sigma_eff']:>9.3f} {r['sigma_pois']:>9.3f} "
               f"{r['sigma_eff_over_pois']:>9.4f} {r['median_over_sigma_eff']:>12.4f} "
               f"{r['mean_over_sigma_eff']:>13.4f} {r['ks_p_eff']:>9.4f} {r['ks_p_pois']:>9.4f}")
    log_fn("=" * 130)

    log_fn(f"\nTheoretical half-Gaussian moments: median/sigma=0.6745, mean/sigma=0.7979.")

    # Decade stability of suppression factor for each (k, w) combination
    log_fn("\nDecade-stability of sigma_eff/sigma_pois  (per (k, w) across anchors):")
    by_kw: dict[tuple, list] = {}
    for r in summary:
        by_kw.setdefault((r["k_label"], r["w_label"]), []).append(r)
    for (k_lab, w_lab), cells in sorted(by_kw.items()):
        cells.sort(key=lambda c: c["anchor_log10"])
        ratios = [c["sigma_eff_over_pois"] for c in cells]
        anchors = [c["anchor_log10"] for c in cells]
        rng = max(ratios) - min(ratios)
        log_fn(f"  {k_lab:>10s} {w_lab:>7s}: anchors={anchors}, "
               f"sigma_eff/pois = {[f'{r:.3f}' for r in ratios]}, range={rng:.3f}")

    # KS-pass rate
    log_fn("\nKS-pass rate (vs half-Gaussian):")
    p_eff_all = [r["ks_p_eff"] for r in summary]
    p_pois_all = [r["ks_p_pois"] for r in summary]
    log_fn(f"  cells with KSp_eff > 0.10:  {sum(1 for p in p_eff_all if p > 0.1)} / {len(p_eff_all)}")
    log_fn(f"  cells with KSp_eff > 0.05:  {sum(1 for p in p_eff_all if p > 0.05)} / {len(p_eff_all)}")
    log_fn(f"  cells with KSp_pois > 0.10: {sum(1 for p in p_pois_all if p > 0.1)} / {len(p_pois_all)}")
    log_fn(f"  median KSp_eff:  {sorted(p_eff_all)[len(p_eff_all) // 2]:.4f}")
    log_fn(f"  median KSp_pois: {sorted(p_pois_all)[len(p_pois_all) // 2]:.4f}")

    # Write log file
    log_path.write_text("\n".join(log_lines) + "\n")
    print(f"\nLog -> {log_path}")


if __name__ == "__main__":
    main()

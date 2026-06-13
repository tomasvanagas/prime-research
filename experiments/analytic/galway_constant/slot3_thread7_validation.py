"""
P5 / Thread 10 (commit, slot 3) — Thread-7-shape validation via RMS scaling and
GUE-corrected K_emp prediction.

Slot 2 path (b) extrapolated K_emp(log10 x = 5.5) ≈ 51,778 under err ~ 1/√K
applied to the worst-case-of-N=30 statistic. That extrapolation matched
Thread-7-shape's prediction c_emp = 0.596 within 4%, but the err ~ 1/√K shape
itself was only weakly justified at slot 2.

Slot 3 strengthens the case via a different lens. Instead of fitting the noisy
worst-case-of-N |err|, it operates on the RMS over N samples, which is a
direct estimate of σ_eff(x, K) — exactly the quantity Thread 7 modelled. We
then compare σ_eff(x, K_max=20000) to Thread 7's σ_pred formula and predict
K_emp under the half-Gaussian-tail map K_emp ↔ σ_pred(x, K) · √(2 ln N) = ε.

Two contributions:

  1. **σ_eff/σ_pred ratio cross-check at K=20000.** For log10 x ∈ {5.0, 5.3,
     5.5} (slot 2's high-x extension), measure σ_eff = √(⟨err²⟩) over N=30
     and compare to σ_pred(x, 20000). The slot 2 extended_results.md
     reported the per-K ratio as drifting upward through 1.0; slot 3 fits
     a power-law σ_eff(K) = a/K^p over K=2000..20000 and reports both the
     fit and the K=K_max ratio.

  2. **GUE-corrected Thread-7-shape K_emp prediction.** Using factor 0.755
     ± 0.06 (Thread 7 typical-case mean) AND factor 1.0 (no correction —
     slot 1 reported worst-case GUE factor 0.796 ± 0.092), compute
     K_emp_pred at log10 x ∈ {4.0, 4.5, 5.0, 5.3, 5.5, 5.7, 6.0} via
     binary search on σ_pred(x, K) · √(2 ln 30) = ε = 1.

This is a pure analysis script — no new MPC compute required. Reads
slot1_summary.csv + slot2_finegrid_summary.csv + slot2_extended_summary.csv +
slot2_extended_traces.csv. Total runtime: < 1 s.

Outputs:
  - slot3_thread7_validation_summary.csv
  - slot3_thread7_validation_results.md (run separately, this script writes
    the data)
"""
from __future__ import annotations

import csv
import math
from collections import defaultdict
from typing import Dict, List, Tuple

import numpy as np


def sigma_pred(x: float, K: int) -> float:
    """Thread 7 σ_pred(x, K) = √x · log K / (π · √(2K) · log x)."""
    return math.sqrt(x) * math.log(K) / (math.pi * math.sqrt(2 * K) * math.log(x))


def find_K_emp(x: float, eps: float, N: int, factor: float) -> int:
    """Binary search for K such that factor · σ_pred(x, K) · √(2 ln N) ≈ ε.

    Models the worst-case-of-N statistic as a Gaussian with std factor · σ_pred,
    so the typical worst is factor · σ_pred · √(2 ln N).
    """
    f = math.sqrt(2 * math.log(N))
    target = eps / (f * factor)  # σ_pred we need
    lo, hi = 100, 10**8
    for _ in range(60):
        mid = (lo + hi) // 2
        if sigma_pred(x, mid) > target:
            lo = mid
        else:
            hi = mid
    return mid


def fit_rms_scaling(traces_csv: str, anchors: List[str], K_lo: int = 2000,
                    K_hi: int = 20000) -> Dict[str, Tuple[float, float, int]]:
    """Read traces, fit log RMS = log a - p log K over K_lo..K_hi.

    Returns dict {log10_x_str: (a, p, n_points)}.
    """
    buckets = defaultdict(list)
    with open(traces_csv) as f:
        rdr = csv.DictReader(f)
        for row in rdr:
            lx = row["log10_x_anchor"]
            K = int(row["K"])
            e = float(row["abs_err"])
            buckets[(lx, K)].append(e)

    results = {}
    for lx in anchors:
        Ks, rms = [], []
        for (lxx, K), es in sorted(buckets.items()):
            if lxx != lx or K < K_lo or K > K_hi:
                continue
            Ks.append(K)
            rms.append(math.sqrt(sum(e * e for e in es) / len(es)))
        if len(Ks) < 3:
            continue
        x_arr = np.log(np.array(Ks))
        y_arr = np.log(np.array(rms))
        coef = np.polyfit(x_arr, y_arr, 1)
        p = -float(coef[0])
        a = math.exp(float(coef[1]))
        results[lx] = (a, p, len(Ks))
    return results


def main() -> None:
    import os
    HERE = os.path.dirname(os.path.abspath(__file__))
    traces_csv = os.path.join(HERE, "slot2_extended_traces.csv")

    rms_fits = fit_rms_scaling(traces_csv, ["5.0", "5.3", "5.5"])

    rows: List[List[str]] = []

    print("=" * 78)
    print("Slot 3 — Part A: σ_eff vs σ_pred at K=20000 (Thread 7 cross-check)")
    print("=" * 78)
    print(f"{'log10 x':>8s} {'x':>10s} {'a':>8s} {'p':>6s} "
          f"{'σ_eff(20k)':>11s} {'σ_pred(20k)':>12s} {'ratio':>8s}")
    print("-" * 78)
    for lx in ["5.0", "5.3", "5.5"]:
        a, p, n = rms_fits[lx]
        x = 10 ** float(lx)
        sigma_eff_20k = a / 20000 ** p
        sigma_pred_20k = sigma_pred(x, 20000)
        ratio = sigma_eff_20k / sigma_pred_20k
        print(f"{lx:>8s} {x:>10.0f} {a:>8.3f} {p:>6.3f} "
              f"{sigma_eff_20k:>11.4f} {sigma_pred_20k:>12.4f} {ratio:>8.3f}")
        rows.append(["A", lx, f"{x:.0f}", f"{a:.4g}", f"{p:.4g}",
                     f"{sigma_eff_20k:.4g}", f"{sigma_pred_20k:.4g}",
                     f"{ratio:.4g}", "", ""])

    print()
    print("=" * 78)
    print("Slot 3 — Part B: GUE-corrected Thread-7 K_emp predictions")
    print("=" * 78)
    print(f"{'log10 x':>8s} {'x':>11s} "
          f"{'K_emp f=1':>10s} {'c_emp f=1':>10s} "
          f"{'K_emp f=.755':>13s} {'c_emp f=.755':>13s}")
    print("-" * 78)
    for lx in [4.0, 4.5, 5.0, 5.3, 5.5, 5.7, 6.0]:
        x = 10 ** lx
        K1 = find_K_emp(x, eps=1.0, N=30, factor=1.0)
        c1 = K1 / (math.sqrt(x) * math.log(x) ** 2)
        K2 = find_K_emp(x, eps=1.0, N=30, factor=0.755)
        c2 = K2 / (math.sqrt(x) * math.log(x) ** 2)
        print(f"{lx:>8.2f} {x:>11.0f} "
              f"{K1:>10d} {c1:>10.4f} "
              f"{K2:>13d} {c2:>13.4f}")
        rows.append(["B", f"{lx:.2f}", f"{x:.0f}",
                     "", "", "", "", "",
                     f"{K1}|{c1:.4f}", f"{K2}|{c2:.4f}"])

    print()
    print("=" * 78)
    print("Slot 3 — Part C: comparison vs slot 1+2 empirical c_emp")
    print("=" * 78)

    slot_summary = {}  # (log10_x_str, source) -> c_emp
    for path, source in [("slot1_summary.csv", "S1"),
                         ("slot2_finegrid_summary.csv", "S2f"),
                         ("slot2_extended_summary.csv", "S2x")]:
        full = os.path.join(HERE, path)
        with open(full) as f:
            rdr = csv.DictReader(f)
            for r in rdr:
                if r["eps"] != "1.0":
                    continue
                if r["c_emp"].strip() == "nan":
                    continue
                slot_summary[(r["log10_x_anchor"], source)] = float(r["c_emp"])

    slot1_anchors = ["4.00", "4.50", "5.00"]
    slot2f_anchors = [f"{a:.2f}" for a in [4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6,
                                            4.7, 4.8, 4.9]]
    slot2x_anchors = ["5.00", "5.30", "5.50"]
    print(f"{'log10 x':>8s} {'c_emp_S1':>10s} {'c_emp_S2f':>11s} "
          f"{'c_emp_S2x':>11s} {'c_emp_T7(.755)':>16s} {'c_emp_T7(1.0)':>14s}")
    print("-" * 78)
    seen = set()
    for lx in sorted(set(slot1_anchors) | set(slot2f_anchors) | set(slot2x_anchors)):
        x = 10 ** float(lx)
        c_emp_T7 = find_K_emp(x, 1.0, 30, 0.755) / (math.sqrt(x) * math.log(x) ** 2)
        c_emp_T7_no = find_K_emp(x, 1.0, 30, 1.0) / (math.sqrt(x) * math.log(x) ** 2)
        s1 = slot_summary.get((lx, "S1"), float("nan"))
        s2f = slot_summary.get((lx, "S2f"), float("nan"))
        s2x = slot_summary.get((lx, "S2x"), float("nan"))
        s1s = "  -      " if math.isnan(s1) else f"{s1:>10.4f}"
        s2fs = "   -       " if math.isnan(s2f) else f"{s2f:>11.4f}"
        s2xs = "   -       " if math.isnan(s2x) else f"{s2x:>11.4f}"
        print(f"{lx:>8s} {s1s} {s2fs} {s2xs} "
              f"{c_emp_T7:>16.4f} {c_emp_T7_no:>14.4f}")
        rows.append(["C", lx, f"{x:.0f}",
                     f"{s1:.4f}" if not math.isnan(s1) else "-",
                     f"{s2f:.4f}" if not math.isnan(s2f) else "-",
                     f"{s2x:.4f}" if not math.isnan(s2x) else "-",
                     "", "",
                     f"{c_emp_T7:.4f}",
                     f"{c_emp_T7_no:.4f}"])

    out_csv = os.path.join(HERE, "slot3_thread7_validation_summary.csv")
    with open(out_csv, "w") as f:
        w = csv.writer(f)
        w.writerow(["part", "log10_x", "x",
                    "rms_a_or_c_S1", "rms_p_or_c_S2f", "sigma_eff_20k_or_c_S2x",
                    "sigma_pred_20k", "ratio",
                    "K1_c1_or_c_T7_755", "K2_c2_or_c_T7_1"])
        w.writerows(rows)

    print(f"\nWrote {out_csv}.")


if __name__ == "__main__":
    main()

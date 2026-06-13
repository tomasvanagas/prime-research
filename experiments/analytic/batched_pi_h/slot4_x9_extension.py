"""
slot4_x9_extension.py — Thread 8 (P2) commit slot 4 (4b OPTIONAL)

Extends the slot3_q_truncation.py measurement to x = 10^9 anchor for a
third-decade scaling validation of the knee Q* = sqrt(x)/log(x) ~ 1525.

Reuses primitives from slot3_q_truncation.py. Runs ANCHOR=10^9 only
(slot 3 already covered 10^7 and 10^8). Outputs slot4_x9_*.csv with the
same structure as slot3 outputs.

Memory: sieve [1, 10^9 + h_max] is ~1 GB bytearray; pair AND-shift per h
is another ~1 GB transient numpy array. Total peak ~2.5 GB on a 26-h
ensemble. System has 122 GB total / 109 GB available, well within budget.

Time: at 10^8, slot 3 took 16s. At 10^9, expected ~3 minutes for sieve +
pair-count.
"""

from __future__ import annotations

import csv
import math
import time
from pathlib import Path

# Reuse all primitives from slot 3
from slot3_q_truncation import (
    H_VALUES,
    Q_VALUES_FINITE,
    Q_INF_LABEL,
    N_X_PER_ANCHOR,
    factor_h_odd_primes,
    s_q,
    s_q_cost_h,
    log_uniform_xs,
    evaluate_anchor,
    ks_stat,
    ks_pvalue,
)

ANCHOR_LOG10 = 9  # third-decade extension


def main() -> None:
    out_dir = Path(__file__).parent

    print(f"\n=== Slot 4 (4b) extension to x = 10^{ANCHOR_LOG10} ===\n")

    # Pre-compute factorisations and S_Q for each (h, Q).
    primes_by_h = {h: factor_h_odd_primes(h) for h in H_VALUES}
    max_p_by_h = {h: (max(primes_by_h[h]) if primes_by_h[h] else 0) for h in H_VALUES}
    s_q_table: dict[tuple[int, float], float] = {}
    for h in H_VALUES:
        s_q_table[(h, math.inf)] = s_q(h, math.inf, primes_by_h[h])
        for Q in Q_VALUES_FINITE:
            s_q_table[(h, float(Q))] = s_q(h, float(Q), primes_by_h[h])

    # Evaluate anchor 10^9.
    t_start = time.perf_counter()
    pi_h, li2_xs = evaluate_anchor(ANCHOR_LOG10)
    t_anchor = time.perf_counter() - t_start
    print(f"\nAnchor evaluation total: {t_anchor:.1f}s")

    # Build full data row list (anchor 10^9 only).
    all_data: list[dict] = []
    for x_j, li2_xj in sorted(li2_xs.items()):
        for h in H_VALUES:
            pi_h_xj = pi_h[(x_j, h)]
            S_inf_h = s_q_table[(h, math.inf)]
            row = {
                "anchor_log10": ANCHOR_LOG10,
                "x": x_j,
                "h": h,
                "max_p_h": max_p_by_h[h],
                "S_inf": S_inf_h,
                "pi_h": pi_h_xj,
                "li2_x": li2_xj,
                "sqrt_x": math.sqrt(x_j),
                "sigma_pred_pois": math.sqrt(S_inf_h * li2_xj) if S_inf_h > 0 else 0.0,
            }
            for Q in Q_VALUES_FINITE:
                S_Q_val = s_q_table[(h, float(Q))]
                HL_Q = S_Q_val * li2_xj
                r_Q = pi_h_xj - HL_Q
                row[f"S_{Q}"] = S_Q_val
                row[f"r_{Q}"] = r_Q
            HL_inf = S_inf_h * li2_xj
            row["r_inf"] = pi_h_xj - HL_inf
            all_data.append(row)

    raw_path = out_dir / "slot4_x9_data.csv"
    with raw_path.open("w", newline="") as f:
        w = csv.writer(f)
        keys = list(all_data[0].keys())
        w.writerow(keys)
        for r in all_data:
            w.writerow([
                r[k] if isinstance(r[k], int)
                else f"{r[k]:.6f}" if isinstance(r[k], float)
                else r[k]
                for k in keys
            ])
    print(f"\nRaw rows ({len(all_data)}) -> {raw_path}")

    # Cross-h aggregates per (anchor, x_j, Q)
    cross_h_groups: dict[tuple[int, int], list[dict]] = {}
    for r in all_data:
        key = (r["anchor_log10"], r["x"])
        cross_h_groups.setdefault(key, []).append(r)

    cross_h_rows = []
    for (anc, x_j), g in sorted(cross_h_groups.items()):
        N = len(g)
        for Q_label, Q_key in [(str(Q), f"r_{Q}") for Q in Q_VALUES_FINITE] + [("inf", "r_inf")]:
            residuals = [row[Q_key] for row in g]
            abs_residuals = [abs(r) for r in residuals]
            sigma_HL_Q = math.sqrt(sum(r * r for r in residuals) / N)
            mean_abs = sum(abs_residuals) / N
            median_abs = sorted(abs_residuals)[N // 2]

            if sigma_HL_Q > 0:
                ratios_eff = [a / sigma_HL_Q for a in abs_residuals]
                D_eff = ks_stat(ratios_eff)
                p_eff = ks_pvalue(D_eff, N)
            else:
                D_eff, p_eff = float("nan"), float("nan")

            cross_h_rows.append({
                "anchor_log10": anc,
                "x": x_j,
                "Q": Q_label,
                "N_h": N,
                "sigma_HL_Q": sigma_HL_Q,
                "mean_abs": mean_abs,
                "median_abs": median_abs,
                "ks_D_eff": D_eff,
                "ks_p_eff": p_eff,
            })

    cross_h_path = out_dir / "slot4_x9_cross_h.csv"
    with cross_h_path.open("w", newline="") as f:
        w = csv.writer(f)
        keys = list(cross_h_rows[0].keys())
        w.writerow(keys)
        for r in cross_h_rows:
            w.writerow([
                r[k] if isinstance(r[k], (int, str))
                else f"{r[k]:.6f}"
                for k in keys
            ])
    print(f"Cross-h aggregates ({len(cross_h_rows)}) -> {cross_h_path}")

    # Knee summary at anchor 10^9
    print("\n" + "=" * 130)
    print("x = 10^9 cross-h sigma_HL_Q vs Q  (knee predicted: max_p ~ sqrt(10^9)/log(10^9) ~ 1525)")
    print("=" * 130)
    print(f"{'anc':>5s} {'x':>13s} | "
          + " ".join([f"Q={q:<5d}" for q in Q_VALUES_FINITE])
          + f"  Q=inf       knee_Q  knee_max_p")
    knee_summary = []
    for (anc, x_j), g in sorted(cross_h_groups.items()):
        sigmas_by_Q = {}
        for Q in Q_VALUES_FINITE + [math.inf]:
            r_key = "r_inf" if Q == math.inf else f"r_{Q}"
            residuals = [row[r_key] for row in g]
            sigmas_by_Q[Q] = math.sqrt(sum(r * r for r in residuals) / len(residuals))
        sigma_inf = sigmas_by_Q[math.inf]
        knee = None
        for Q in Q_VALUES_FINITE:
            if sigmas_by_Q[Q] <= 1.05 * sigma_inf:
                knee = Q
                break
        knee_max_p = (
            max((max_p_by_h[h] for h in H_VALUES if max_p_by_h[h] <= knee), default=0)
            if knee else 0
        )
        knee_summary.append({
            "anchor_log10": anc,
            "x": x_j,
            "sigma_inf": sigma_inf,
            "knee_Q": knee,
            "knee_max_p": knee_max_p,
            **{f"sigma_Q_{Q}": sigmas_by_Q[Q] for Q in Q_VALUES_FINITE},
            "sigma_Q_inf": sigma_inf,
        })
        sigma_str = " ".join([f"{sigmas_by_Q[Q]:>7.0f}" for Q in Q_VALUES_FINITE])
        print(f"{anc:>5d} {x_j:>13d} | {sigma_str}  {sigma_inf:>7.0f}      {str(knee):>5s}    {knee_max_p:>3d}")

    knee_path = out_dir / "slot4_x9_knee.csv"
    with knee_path.open("w", newline="") as f:
        w = csv.writer(f)
        keys = list(knee_summary[0].keys())
        w.writerow(keys)
        for r in knee_summary:
            w.writerow([
                r[k] if isinstance(r[k], (int, str))
                else (f"{r[k]:.6f}" if r[k] is not None and isinstance(r[k], float) else r[k])
                for k in keys
            ])
    print(f"\nKnee summary -> {knee_path}")

    print(f"\nWall time: {time.perf_counter() - t_start:.1f}s\n")


if __name__ == "__main__":
    main()

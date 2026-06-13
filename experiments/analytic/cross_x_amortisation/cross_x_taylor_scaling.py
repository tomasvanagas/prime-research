"""
Cross-x amortisation, slot 2 supplementary: Taylor-P speedup as a function
of (K, M).

Companion to cross_x_batched_evaluator.py. The main script established
that Taylor-2 at (x=1e6, K=200, M=16) gives ~13x speedup with max error
2.3e-8 (well below the pi(x) +/- 0.5 threshold). This script tests the
asymptotic structure: does Taylor speedup scale with K (it should, by
constant multiplier), and how large can M be before Taylor-2 hits its
accuracy ceiling?

This data feeds the slot-2 falsifier conclusion: Taylor / Schoenhage
multipoint is a constant-factor speedup in K (no asymptotic alpha
reduction), with cluster-width-bounded M.

Usage:
  python3 cross_x_taylor_scaling.py
"""

from __future__ import annotations

import csv
import math
import os
import sys
import time
from typing import List, Tuple

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from cross_x_batched_evaluator import (  # noqa: E402
    batched_direct, batched_taylor, cluster_around, get_zeros,
)


def main():
    x = 1e6
    P = 2
    cases = [
        # (K, M values to test)
        (100, [4, 16, 64, 256]),
        (200, [4, 16, 64, 256]),
        (400, [4, 16, 64, 256]),
    ]
    K_max = max(K for K, _ in cases)
    print(f"Loading {K_max} zeros (x={x}, P={P})")
    gammas = get_zeros(K_max, dps=30)
    print(f"  loaded gamma_1={gammas[0]:.4f}, gamma_K_max={gammas[-1]:.4f}")
    rows = []
    print(f"\n{'K':>5}  {'M':>5}  {'T_dir':>8}  {'T_tay':>8}  "
          f"{'speedup':>8}  {'maxerr':>10}  {'per0err':>10}")
    for K, Ms in cases:
        for M in Ms:
            cl = cluster_around(x, M, mode="integer")
            d_a, T_d, _ = batched_direct(cl, K, gammas, include_R_x=False)
            t_a, T_ts, T_te, T_tt = batched_taylor(cl, K, gammas,
                                                   x0=x, P=P,
                                                   include_R_x=False)
            d = np.array(d_a); t = np.array(t_a)
            mx = float(np.max(np.abs(d - t)))
            per_zero = mx / K
            speedup = T_d / T_tt if T_tt > 0 else 0.0
            print(f"{K:>5}  {M:>5}  {T_d:>8.3f}  {T_tt:>8.3f}  "
                  f"{speedup:>7.2f}x  {mx:>10.3e}  {per_zero:>10.3e}")
            rows.append({
                "x": x, "K": K, "M": M, "P": P,
                "T_direct_total": T_d,
                "T_taylor_setup": T_ts,
                "T_taylor_eval": T_te,
                "T_taylor_total": T_tt,
                "speedup": speedup,
                "max_err": mx,
                "per_zero_err": per_zero,
                "passed_pi_threshold": int(mx <= 0.5),
            })

    out = os.path.join(HERE, "cross_x_taylor_scaling.csv")
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"\nWrote {out}")

    # Asymptotic Taylor speedup formula (slot-2 prediction):
    # T_dir = M * a * K, T_tay_setup = a' * K, T_tay_eval = M * b * K * P.
    # speedup = a / (a'/M + b * P).
    # As M -> infty, speedup -> a / (b * P) = constant, K-independent.
    #
    # Verify by checking: at fixed M, does speedup vary with K?
    print("\n=== K-scaling of Taylor speedup at fixed M ===")
    for M_target in [4, 16, 64, 256]:
        ratios = []
        for K, _ in cases:
            for r in rows:
                if r["M"] == M_target and r["K"] == K:
                    ratios.append((K, r["speedup"]))
                    break
        if len(ratios) >= 2:
            vals = [s for _, s in ratios]
            min_s, max_s = min(vals), max(vals)
            spread = (max_s - min_s) / min_s if min_s > 0 else 0
            print(f"  M={M_target}: speedups = "
                  f"{', '.join(f'K={K}:{s:.1f}x' for K, s in ratios)}  "
                  f"(spread={spread:.1%})")


if __name__ == "__main__":
    main()

"""§D6.c addendum — W-trick check.

Tests whether the residual HL-structure on mu_plus, mu_minus, sqfree
indicators (Q2 ≈ 1.0384) is also W-trick-killable, parallel to E2.13's
chi_P Q2 ≈ 1.003 at W=210.
"""

from __future__ import annotations

import json
import time

import numpy as np

from mu_weighted_chi_p_uk import panel_arrays, u2_pow4


def w_trick(f_full: np.ndarray, N_out: int, W: int, b: int) -> np.ndarray:
    out = np.zeros(N_out, dtype=np.float64)
    for n in range(N_out):
        idx = W * n + b
        if 0 <= idx < len(f_full):
            out[n] = f_full[idx]
    return out


def measure_short(name: str, f: np.ndarray) -> dict:
    N = len(f)
    mean_f = float(f.mean())
    l2sq = float((f * f).mean())
    u2_4 = u2_pow4(f)
    rec = {"N": N, "mean": mean_f, "L2_sq": l2sq, "u2_pow4": u2_4}
    if mean_f > 1e-12:
        rec["Q2"] = u2_4 / (mean_f ** 4)
    if l2sq > 1e-30:
        rec["Q2_norm"] = N * u2_4 / (2.0 * (l2sq ** 2))
    return rec


def main():
    N_short = 2048
    N_full = 220 * N_short + 100  # need at least W*N_short room for W up to 210
    print(f"Building panel at N_full = {N_full} ...")
    t0 = time.time()
    panel_full = panel_arrays(N_full)
    print(f"  done in {time.time()-t0:.2f}s")

    out = {"N_short": N_short, "N_full": N_full, "results": {}}

    # b = 1 is coprime to W for all W in {6, 30, 210}.
    for W in [6, 30, 210]:
        b = 1
        out["results"][f"W={W}"] = {}
        for name in ["chi_P", "sqfree", "mu_plus", "mu_minus",
                     "lam_plus", "lam_minus", "semi_primes"]:
            f_full = panel_full[name]
            # Need W*N_short + b <= N_full
            if W * N_short + b >= N_full:
                continue
            f_w = w_trick(f_full, N_short, W, b)
            rec = measure_short(name, f_w)
            out["results"][f"W={W}"][name] = rec
            q2_str = f"{rec.get('Q2', float('nan')):.4f}"
            print(f"  W={W:3d}  {name:14s}  mean={rec['mean']:.4f}  "
                  f"Q2={q2_str:>9s}  Q2_norm={rec['Q2_norm']:.4f}")
        print()

    # Bare baseline at N = N_short for comparison.
    out["bare_at_N_short"] = {}
    print(f"Bare values at N = {N_short} for comparison:")
    panel_short = panel_arrays(N_short)
    for name in ["chi_P", "sqfree", "mu_plus", "mu_minus",
                 "lam_plus", "lam_minus", "semi_primes"]:
        rec = measure_short(name, panel_short[name])
        out["bare_at_N_short"][name] = rec
        q2_str = f"{rec.get('Q2', float('nan')):.4f}"
        print(f"  bare {name:14s}  mean={rec['mean']:.4f}  Q2={q2_str:>9s}")

    with open("wtrick_check.json", "w") as f:
        json.dump(out, f, indent=2, default=str)
    print("\nWrote wtrick_check.json")


if __name__ == "__main__":
    main()

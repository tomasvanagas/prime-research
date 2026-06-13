"""
Extend Voronin shift scan to t in [1e5, 1e6] for the cleanest target
(exp(-s)), to test the empirical density model d(eps) ~ exp(-c/eps)
predicted by the [1e1, 1e5] scan in voronin_polylog.py.

If d(eps) ~ exp(-0.53 / eps) holds, we predict:
  ε = 0.06: d ~ exp(-8.83) ~ 1.5e-4    => 1 hit per ~7000 shifts
  ε = 0.05: d ~ exp(-10.60) ~ 2.5e-5   => 1 hit per ~40000 shifts
  ε = 0.04: d ~ exp(-13.25) ~ 1.7e-6   => 0 hits in 1000 shifts
With 1000 geometric shifts in [1e5, 1e6] (one decade, ~1000 trials),
we should marginally see eps=0.05-0.06 hits but not eps<0.04.

Polynomial alternative d(eps) ~ eps^k (k=2.05 from main scan):
  ε = 0.06: d ~ 0.06^2.05 ~ 0.0033    => ~3 hits per 1000 shifts
  ε = 0.05: d ~ 0.05^2.05 ~ 0.0023    => ~2 hits per 1000 shifts
  ε = 0.04: d ~ 0.04^2.05 ~ 0.0014    => ~1 hit per 1000 shifts
  ε = 0.03: d ~ 0.03^2.05 ~ 0.00074   => ~0-1 hits per 1000

If we see polynomial scaling, we'd expect ~6 hits at eps=0.04 in 1000 trials.
If we see exp(-c/eps) scaling, we'd expect 0 hits at eps=0.04.
"""

import json
import math
import time
from pathlib import Path

import mpmath
from mpmath import mp, mpc, mpf

mp.dps = 18

S0 = mpc("0.75", "0")
R = mpf("0.05")
N_BOUNDARY = 12
THETAS = [2 * math.pi * k / N_BOUNDARY for k in range(N_BOUNDARY)]
KPTS = [S0 + R * mpc(math.cos(th), math.sin(th)) for th in THETAS]


def f_exp_neg_s(s):
    return mpmath.exp(-s)


def main():
    n_shifts = 800
    t_lo, t_hi = 1.0e5, 1.0e6
    log_lo, log_hi = math.log10(t_lo), math.log10(t_hi)
    t_grid = [10 ** (log_lo + (log_hi - log_lo) * k / (n_shifts - 1))
              for k in range(n_shifts)]

    f_vals = [f_exp_neg_s(p) for p in KPTS]

    errs = []
    t0 = time.time()
    for i, t in enumerate(t_grid):
        sup_e = 0.0
        for p, fp in zip(KPTS, f_vals):
            z = mpmath.zeta(p + mpc(0, t))
            e = float(abs(z - fp))
            if e > sup_e:
                sup_e = e
        errs.append(sup_e)
        if (i + 1) % 25 == 0 or i == n_shifts - 1:
            elapsed = time.time() - t0
            avg = elapsed / (i + 1)
            eta = avg * (n_shifts - i - 1)
            print(f"  shift {i+1}/{n_shifts}  t={t:.0f}  "
                  f"elapsed={elapsed:.1f}s  ETA={eta:.0f}s",
                  flush=True)

    eps_grid = [0.5, 0.3, 0.2, 0.15, 0.10, 0.08, 0.06, 0.05, 0.04, 0.03]
    rows = []
    for eps in eps_grid:
        n_below = sum(1 for e in errs if e < eps)
        density = n_below / len(errs)
        t_first = None
        for t, e in zip(t_grid, errs):
            if e < eps:
                t_first = t
                break
        rows.append({"eps": eps, "n_below": n_below, "density": density,
                     "t_first": t_first})

    out = {
        "config": {
            "target": "exp(-s)",
            "n_shifts": n_shifts,
            "t_lo": t_lo,
            "t_hi": t_hi,
            "S0": str(S0),
            "R": str(R),
            "n_boundary": N_BOUNDARY,
        },
        "rows": rows,
        "t_grid": t_grid,
        "errs": errs,
        "min_err_overall": min(errs),
    }
    out_path = Path(__file__).parent / "voronin_extend_results.json"
    with open(out_path, "w") as fh:
        json.dump(out, fh, indent=2)

    print("\n========= extension scan summary (target exp(-s)) =========")
    print(f"min err overall: {min(errs):.4f}")
    print(f"{'eps':>8} {'n_below':>8} {'density':>10} {'t_first':>10}")
    for r in rows:
        tf = r["t_first"]
        tf_s = "None" if tf is None else f"{tf:.0f}"
        print(f"{r['eps']:>8.3f} {r['n_below']:>8} "
              f"{r['density']:>10.6f} {tf_s:>10}")

    # Compare to predictions.
    print("\n--- model comparison (n_hits in 800 shifts) ---")
    print(f"{'eps':>8} {'observed':>10} {'predicted_poly':>16} "
          f"{'predicted_exp(-c/eps)':>22}")
    for r in rows:
        pred_poly = (r["eps"] ** 2.05) * 800
        pred_exp = math.exp(-0.53 / r["eps"]) * 800
        print(f"{r['eps']:>8.3f} {r['n_below']:>10} "
              f"{pred_poly:>16.3f} {pred_exp:>22.6f}")


if __name__ == "__main__":
    main()

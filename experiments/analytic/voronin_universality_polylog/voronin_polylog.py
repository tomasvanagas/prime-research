"""
Voronin universality with effective shifts as a polylog approximator.

Attack vector ATTACK_VECTORS.md §B4.

Voronin (1975) proved: for any non-vanishing analytic f on a compact
K ⊂ {1/2 < Re(s) < 1} with connected complement, the set
   { t > 0 : sup_{s in K} |zeta(s + it) - f(s)| < eps }
has positive lower density.  Garunkstis 2003 gives the only known
effective bound on the smallest valid t:
   t* <= exp(exp(10 / eps^13))    -- astronomical.

QUESTION: is the empirical t*(f, eps) scaling polynomial or
super-polynomial in 1/eps for *natural* target f?  If polynomial,
evaluating zeta(s + i t*) at one shift gives a polylog reconstruction
of f, hence a polylog approximator of any quantity expressible via f.

EXPERIMENT:
  K = closed disc {|s - 0.75| <= 0.05}
  targets f:
    f1(s) = exp(s)             -- entire, non-vanishing
    f2(s) = exp(-s)            -- entire, non-vanishing
    f3(s) = 1 / (s - 0.5)      -- non-vanishing on K, simple pole outside
    f4(s) = 1 / (s - 0.4)      -- different pole location
    f5(s) = exp(s * s)         -- entire, slower variation
    f6(s) = exp(0.3 * s)       -- gentler exponential
  shift grid: 800 t-values, geometrically spaced in [10, 10^5]
  precision: mp.dps = 18 (sufficient at moderate t)

OUTPUT:
  - per-(f, eps) smallest t* where sup |zeta(s + it*) - f(s)| < eps
  - log10(t*) vs log10(1/eps) regression slope => polynomial-or-not
  - direct check of running min error curve
  - comparison vs Garunkstis worst-case bound
  - JSON results dump for analysis.

Cross-domain reference:
  Garunkstis 2003 "The effective universality theorem for the Riemann
                   zeta-function" (Acta Arith. 113).
  Steuding 2007  "Value-Distribution of L-Functions" Springer LNM 1877.
"""

import json
import math
import os
import time
from pathlib import Path

import mpmath
from mpmath import mp, mpc, mpf

mp.dps = 18  # 18 decimal digits is enough at t <= 10^5; speeds up zeta eval


# -------------------------- target functions --------------------------

def f_exp_s(s):
    return mpmath.exp(s)


def f_exp_neg_s(s):
    return mpmath.exp(-s)


def f_inv_s_05(s):
    return 1 / (s - mpf("0.5"))


def f_inv_s_04(s):
    return 1 / (s - mpf("0.4"))


def f_exp_ssq(s):
    return mpmath.exp(s * s)


def f_exp_03s(s):
    return mpmath.exp(mpf("0.3") * s)


TARGETS = {
    "exp(s)": f_exp_s,
    "exp(-s)": f_exp_neg_s,
    "1/(s-0.5)": f_inv_s_05,
    "1/(s-0.4)": f_inv_s_04,
    "exp(s^2)": f_exp_ssq,
    "exp(0.3 s)": f_exp_03s,
}


# -------------------------- compact set K --------------------------

S0 = mpc("0.75", "0")
R = mpf("0.05")
N_BOUNDARY = 12  # max modulus -> sup achieved on boundary

# Sample on the boundary of the disc.
THETAS = [2 * math.pi * k / N_BOUNDARY for k in range(N_BOUNDARY)]


def boundary_pts():
    return [S0 + R * mpc(math.cos(th), math.sin(th)) for th in THETAS]


KPTS = boundary_pts()


def f_values(target_fn):
    return [target_fn(p) for p in KPTS]


# -------------------------- core: sup error at shift t --------------------------

def zeta_on_boundary(t):
    """Return list of zeta(s + it) for s on the boundary of K."""
    return [mpmath.zeta(p + mpc(0, t)) for p in KPTS]


def sup_error_from_zeta(zvals, f_vals):
    """sup_{s in boundary K} |zeta(s + it) - f(s)| from precomputed zeta."""
    err = 0.0
    for z, fp in zip(zvals, f_vals):
        e = float(abs(z - fp))
        if e > err:
            err = e
    return err


# -------------------------- experiment --------------------------

def run(t_grid, targets):
    """For each shift t, record sup_error vs each target.

    Returns dict: target_name -> list of float errors (parallel to t_grid).
    """
    results = {name: [] for name in targets}
    f_vals = {name: f_values(fn) for name, fn in targets.items()}

    t0 = time.time()
    for i, t in enumerate(t_grid):
        zvals = zeta_on_boundary(t)
        for name, fn in targets.items():
            err = sup_error_from_zeta(zvals, f_vals[name])
            results[name].append(err)
        if (i + 1) % 25 == 0 or i == len(t_grid) - 1:
            elapsed = time.time() - t0
            avg = elapsed / (i + 1)
            remain = avg * (len(t_grid) - i - 1)
            print(
                f"  shift {i+1}/{len(t_grid)}  t={t:.1f}  "
                f"elapsed={elapsed:.1f}s  ETA={remain:.0f}s",
                flush=True,
            )
    return results


def smallest_t_below(t_grid, errs, eps):
    """Smallest t in grid such that error < eps; None if none."""
    for t, e in zip(t_grid, errs):
        if e < eps:
            return float(t), float(e)
    return None, None


def fit_loglog_slope(xs, ys):
    """OLS slope of ys vs xs in log space.  Returns (slope, intercept).

    slope > 0 means y grows with x; slope = polynomial-degree if both
    axes are log10.
    """
    n = len(xs)
    if n < 2:
        return None, None
    xm = sum(xs) / n
    ym = sum(ys) / n
    num = sum((xs[i] - xm) * (ys[i] - ym) for i in range(n))
    den = sum((xs[i] - xm) ** 2 for i in range(n))
    if den == 0:
        return None, None
    slope = num / den
    intercept = ym - slope * xm
    return slope, intercept


def main():
    out_dir = Path(__file__).parent
    out_dir.mkdir(exist_ok=True)

    # Geometrically spaced t-grid.  Coarser at first; can refine if signal.
    n_shifts = 800
    t_lo, t_hi = 10.0, 1.0e5
    log_lo, log_hi = math.log10(t_lo), math.log10(t_hi)
    t_grid = [10 ** (log_lo + (log_hi - log_lo) * k / (n_shifts - 1))
              for k in range(n_shifts)]

    print(f"Voronin universality scan: {n_shifts} shifts in "
          f"[{t_lo}, {t_hi}], {len(TARGETS)} targets, "
          f"{N_BOUNDARY} boundary points, mp.dps = {mp.dps}")
    print(f"Disc K: |s - {S0}| <= {R}")

    print("Running shift sweep...")
    results = run(t_grid, TARGETS)

    # For each target, scan eps grid; find smallest t with error < eps.
    eps_grid = [0.5, 0.4, 0.3, 0.25, 0.2, 0.15, 0.12, 0.10, 0.08, 0.06, 0.05]

    summary = {}
    for name, errs in results.items():
        running_min = []
        cur = float("inf")
        for e in errs:
            if e < cur:
                cur = e
            running_min.append(cur)

        rows = []
        log_inv_eps = []
        log_t = []
        for eps in eps_grid:
            t_star, err_at_t = smallest_t_below(t_grid, errs, eps)
            row = {
                "eps": eps,
                "t_star": t_star,
                "err_at_t_star": err_at_t,
            }
            rows.append(row)
            if t_star is not None:
                log_inv_eps.append(math.log10(1.0 / eps))
                log_t.append(math.log10(t_star))

        slope, intercept = fit_loglog_slope(log_inv_eps, log_t)
        # Garunkstis worst-case t bound for the smallest tested eps that
        # has data, capped to avoid overflow.
        garunkstis = []
        for eps in eps_grid:
            try:
                inner = 10.0 / (eps ** 13)
                if inner > 700:
                    garunkstis.append(float("inf"))
                else:
                    outer = math.exp(inner)
                    if outer > 700:
                        garunkstis.append(float("inf"))
                    else:
                        garunkstis.append(math.exp(outer))
            except OverflowError:
                garunkstis.append(float("inf"))

        # Also: for any t_star found, what is the *Voronin density*
        # within [10, t_star]? (positive-density check).
        finite_data = [
            (i, t_grid[i], errs[i])
            for i in range(len(t_grid))
            if errs[i] < eps_grid[-1]  # below smallest tested eps
        ]

        summary[name] = {
            "running_min_final": running_min[-1],
            "min_err_overall": min(errs),
            "min_t_overall": float(t_grid[errs.index(min(errs))]),
            "rows": rows,
            "loglog_slope": slope,
            "loglog_intercept": intercept,
            "garunkstis_bound_per_eps": garunkstis,
            "n_shifts_below_smallest_eps": len(finite_data),
        }

    out_path = out_dir / "voronin_polylog_results.json"
    with open(out_path, "w") as fh:
        json.dump(
            {
                "config": {
                    "n_shifts": n_shifts,
                    "t_lo": t_lo,
                    "t_hi": t_hi,
                    "S0": str(S0),
                    "R": str(R),
                    "n_boundary": N_BOUNDARY,
                    "mp_dps": mp.dps,
                    "eps_grid": eps_grid,
                },
                "summary": summary,
                "t_grid": t_grid,
                "errs_per_target": results,
            },
            fh,
            indent=2,
        )
    print(f"\nResults written to {out_path}")

    # Print a compact summary.
    print("\n" + "=" * 78)
    print("SUMMARY")
    print("=" * 78)
    for name, s in summary.items():
        print(f"\nTarget {name}:")
        print(f"  min error overall: {s['min_err_overall']:.4f} at t = {s['min_t_overall']:.1f}")
        print(f"  loglog slope log10(t*) vs log10(1/eps): "
              f"{s['loglog_slope'] if s['loglog_slope'] else 'NA'}")
        print(f"  N shifts with error < {eps_grid[-1]}: "
              f"{s['n_shifts_below_smallest_eps']}")
        for row in s["rows"]:
            t = row["t_star"]
            print(f"    eps = {row['eps']:.3f}: "
                  f"t* = {t if t is None else f'{t:.0f}'}, "
                  f"err = {row['err_at_t_star']}")


if __name__ == "__main__":
    main()

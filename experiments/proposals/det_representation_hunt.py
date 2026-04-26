"""
Proposal B3: Determinant-representation hunt for pi(x).

Many counting problems admit determinant representations:
    spanning trees     : matrix-tree theorem
    dimer tilings      : Kasteleyn (Pfaffian / determinant)
    non-intersecting paths : Lindström-Gessel-Viennot
    Frobenius traces of constructible sheaves : Grothendieck-Lefschetz

Hypothesis to test: does there exist a polynomial-size matrix M(x) with
*polylog-time computable* entries such that

    det M(x) = pi(x)   (or a simple algebraic variant)

We attack this empirically. Build a feature library

    F = { 1, log x, log log x, sqrt x, x/log x, R^{-1}(x), li(x), ... }

and at matrix size k = 2 enumerate all k^2 = 4 cell choices from F,
forming M(x) = [[f_a(x), f_b(x)], [f_c(x), f_d(x)]]. Compute det M(x_i)
on a target grid of x-values and score by L2 error vs pi(x_i).

The library is intentionally small (~12 features) so 12^4 = 20736
candidates fit in a few seconds.

Win condition: a (k=2 or k=3) candidate with residual << sqrt(N) over
N target values, AND that survives extrapolation to held-out points.

Usage:
    python3 det_representation_hunt.py
"""
from __future__ import annotations
import argparse
import itertools
import json
import time
import numpy as np
import sympy
from mpmath import mp, mpf, li


def Rinv_table(n_max: int) -> np.ndarray:
    """
    Return Rinv[n] for n in [0, n_max], using the convention
        Rinv(n) = mpmath.findroot(R(x) - n, ...)
    where R is Riemann's prime-counting function.

    For speed we instead use the asymptotic R^{-1}(n) ~ p_n approximation
    via Newton iteration on the logarithmic-integral.
    """
    mp.dps = 25
    out = np.zeros(n_max + 1)
    # use Riemann R via mpmath.zeta-based formula --- but that's slow;
    # instead, we use the inverse logarithmic integral, which gives ~50%
    # of digits and is the "smooth part" we'd subtract anyway.
    from mpmath import li, mpf, findroot
    for n in range(2, n_max + 1):
        # initial guess from PNT
        x0 = mpf(n) * mp.log(n)
        try:
            x = findroot(lambda y: li(y) - n, x0, tol=1e-3)
        except Exception:
            x = x0
        out[n] = float(x)
    out[0] = 0.0
    out[1] = 2.0
    return out


def pi_array(x_values: np.ndarray) -> np.ndarray:
    """exact pi(x) at given integer x's, via sympy."""
    return np.array([int(sympy.primepi(int(x))) for x in x_values], dtype=float)


def build_features(x_values: np.ndarray, rinv_lookup: dict[int, float] | None = None) -> dict[str, np.ndarray]:
    x = x_values.astype(float)
    feats = {
        "one":      np.ones_like(x),
        "x":        x,
        "logx":     np.log(x),
        "loglogx":  np.log(np.log(x)),
        "sqrtx":    np.sqrt(x),
        "x_over_logx": x / np.log(x),
        "log2x":    np.log(x) ** 2,
        "li":       np.array([float(li(xi)) for xi in x]),
        "x_inv":    1.0 / x,
        "logx_inv": 1.0 / np.log(x),
        "x_log_logx": x * np.log(np.log(x)),
        "x_minus_logx": x - np.log(x),
    }
    return feats


def hunt_k2(features: dict[str, np.ndarray], target: np.ndarray,
            top_k: int = 20) -> list[dict]:
    """
    For 2x2 matrices M = [[f_a, f_b], [f_c, f_d]],
    det M = f_a * f_d - f_b * f_c.
    """
    keys = list(features.keys())
    best: list[tuple[float, tuple[str, str, str, str], np.ndarray]] = []
    for choice in itertools.product(keys, repeat=4):
        a, b, c, d = choice
        det = features[a] * features[d] - features[b] * features[c]
        # L2 error normalised by length
        err = float(np.sqrt(np.mean((det - target) ** 2)))
        if not np.isfinite(err):
            continue
        best.append((err, choice, det))
        if len(best) > 4 * top_k:
            best.sort(key=lambda r: r[0])
            best = best[:top_k]
    best.sort(key=lambda r: r[0])
    out = []
    for err, choice, det in best[:top_k]:
        out.append({
            "rmse": err,
            "choice": list(choice),
            "det_first5": det[:5].tolist(),
            "target_first5": target[:5].tolist(),
        })
    return out


def hunt_k2_affine(features: dict[str, np.ndarray], target: np.ndarray,
                   top_k: int = 20) -> list[dict]:
    """
    Allow det M(x) + alpha * 1 + beta * (some_feature) = pi(x).
    Equivalently, fit per-candidate 2-parameter affine correction by
    least squares against (1, x).
    Reports RMSE *after* the affine correction.
    """
    keys = list(features.keys())
    x = features.get("x", None)
    if x is None:
        raise RuntimeError("missing 'x'")
    # design column-basis for the affine correction
    Phi = np.column_stack([np.ones_like(x), x, np.log(x)])
    target_vec = target.astype(float)
    best: list[tuple[float, tuple[str, str, str, str], np.ndarray, np.ndarray]] = []
    for choice in itertools.product(keys, repeat=4):
        a, b, c, d = choice
        det = features[a] * features[d] - features[b] * features[c]
        # least-squares: solve   det + Phi @ alpha = target
        # i.e.  Phi @ alpha = target - det
        try:
            alpha, *_ = np.linalg.lstsq(Phi, target_vec - det, rcond=None)
        except Exception:
            continue
        residual = target_vec - det - Phi @ alpha
        err = float(np.sqrt(np.mean(residual ** 2)))
        if not np.isfinite(err):
            continue
        best.append((err, choice, det, alpha))
        if len(best) > 4 * top_k:
            best.sort(key=lambda r: r[0])
            best = best[:top_k]
    best.sort(key=lambda r: r[0])
    out = []
    for err, choice, det, alpha in best[:top_k]:
        out.append({
            "rmse": err,
            "choice": list(choice),
            "alpha_const": float(alpha[0]),
            "alpha_x": float(alpha[1]),
            "alpha_logx": float(alpha[2]),
        })
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--x_max", type=int, default=10000)
    parser.add_argument("--n_grid", type=int, default=200)
    parser.add_argument("--top_k", type=int, default=15)
    args = parser.parse_args()

    print(f"[setup] training grid x in [10, {args.x_max}], "
          f"n_grid = {args.n_grid}")
    x_grid = np.linspace(10, args.x_max, args.n_grid).astype(int)
    x_grid = np.unique(x_grid)
    t0 = time.time()
    target = pi_array(x_grid)
    print(f"[setup] computed pi(x) on {len(x_grid)} points "
          f"in {time.time()-t0:.2f}s")

    print(f"[setup] building feature library at {len(x_grid)} points")
    features = build_features(x_grid)
    print(f"[setup] {len(features)} features: {list(features.keys())}")

    print(f"\n[hunt k=2 raw det]")
    t0 = time.time()
    results_raw = hunt_k2(features, target, top_k=args.top_k)
    print(f"[hunt] enumerated {len(features)**4} candidates "
          f"in {time.time()-t0:.2f}s")
    print(f"top {args.top_k} by RMSE:")
    for i, r in enumerate(results_raw):
        print(f"  {i:2d}. RMSE={r['rmse']:>12.3e}  M = "
              f"[[{r['choice'][0]:>14s},{r['choice'][1]:>14s}],"
              f" [{r['choice'][2]:>14s},{r['choice'][3]:>14s}]]")

    print(f"\n[hunt k=2 + affine (1, x, log x) correction]")
    t0 = time.time()
    results_aff = hunt_k2_affine(features, target, top_k=args.top_k)
    print(f"[hunt] enumerated {len(features)**4} candidates "
          f"in {time.time()-t0:.2f}s")
    print(f"top {args.top_k} by post-affine RMSE:")
    for i, r in enumerate(results_aff):
        print(f"  {i:2d}. RMSE={r['rmse']:>12.3e}  "
              f"alpha=({r['alpha_const']:+.3e},{r['alpha_x']:+.3e},"
              f"{r['alpha_logx']:+.3e})  M = "
              f"[[{r['choice'][0]:>13s},{r['choice'][1]:>13s}],"
              f" [{r['choice'][2]:>13s},{r['choice'][3]:>13s}]]")

    # holdout test: extrapolate the best 5 to x in [x_max+1000, x_max+5000]
    print(f"\n[holdout] extrapolating top 5 raw candidates "
          f"to x in [{args.x_max+1000}, {args.x_max+5000}]")
    holdout = np.linspace(args.x_max + 1000, args.x_max + 5000, 30).astype(int)
    holdout = np.unique(holdout)
    holdout_pi = pi_array(holdout)
    holdout_feats = build_features(holdout)
    extrap_table = []
    for r in results_raw[:5]:
        a, b, c, d = r["choice"]
        det = (holdout_feats[a] * holdout_feats[d]
               - holdout_feats[b] * holdout_feats[c])
        rmse = float(np.sqrt(np.mean((det - holdout_pi) ** 2)))
        extrap_table.append({"choice": r["choice"],
                             "training_rmse": r["rmse"],
                             "holdout_rmse": rmse,
                             "growth_ratio": rmse / max(r["rmse"], 1e-9)})
        print(f"  [[{a:>13s},{b:>13s}],[{c:>13s},{d:>13s}]]"
              f"  train={r['rmse']:.3e}  holdout={rmse:.3e}  "
              f"ratio={rmse/max(r['rmse'],1e-9):.2f}")

    out = {
        "args": vars(args),
        "feature_keys": list(features.keys()),
        "training_grid_size": int(len(x_grid)),
        "raw_top": results_raw,
        "affine_top": results_aff,
        "holdout_extrapolation": extrap_table,
    }
    out_path = (
        "/apps/aplikacijos/prime-research/experiments/proposals/"
        "det_representation_hunt_data.json"
    )
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2, default=float)
    print(f"\n[save] wrote {out_path}")


if __name__ == "__main__":
    main()

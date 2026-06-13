"""
Voronin density analysis: per-eps density of valid shifts vs eps.

Loads voronin_polylog_results.json (the 800-shift main scan), then for
each target f computes:
  - density d(eps) = #{t in scan : sup_err(t) < eps} / N_scan
  - density per t-decade (for measuring "log-Lebesgue density")
  - super-polynomial fit candidates:
      d(eps) ~ exp(-c / eps^a)            (Garunkstis-type, a ~ 13)
      d(eps) ~ exp(-c / eps)              (heuristic)
      t_first(eps) is the smallest t in scan with sup_err < eps
      log(t_first(eps)) vs (1/eps)        (exponential model)
      log(log(t_first(eps))) vs log(1/eps) (double-exponential model)

Output: voronin_density_results.md (markdown summary), voronin_density.json.
"""

import json
import math
from pathlib import Path

HERE = Path(__file__).parent


def fit_loglog(xs, ys):
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
    with open(HERE / "voronin_polylog_results.json") as fh:
        data = json.load(fh)

    t_grid = data["t_grid"]
    errs_per_target = data["errs_per_target"]
    cfg = data["config"]

    # Geometric t-grid bins by decade.
    decades = [(1, 2), (2, 3), (3, 4), (4, 5)]  # log10 ranges

    eps_fine = [0.5, 0.4, 0.3, 0.25, 0.2, 0.15, 0.12, 0.10, 0.08, 0.06,
                0.05, 0.04, 0.03, 0.025, 0.02]

    out = {"config": cfg, "per_target": {}}

    for name, errs in errs_per_target.items():
        per_eps = []
        for eps in eps_fine:
            n_below = sum(1 for e in errs if e < eps)
            density = n_below / len(errs)
            t_first = None
            for t, e in zip(t_grid, errs):
                if e < eps:
                    t_first = t
                    break

            # Per-decade hit count.
            per_dec = {}
            for lo, hi in decades:
                cnt = sum(
                    1
                    for t, e in zip(t_grid, errs)
                    if (10 ** lo <= t < 10 ** hi) and e < eps
                )
                per_dec[f"[10^{lo},10^{hi})"] = cnt

            per_eps.append(
                {
                    "eps": eps,
                    "n_below": n_below,
                    "density": density,
                    "t_first": t_first,
                    "per_decade": per_dec,
                }
            )

        # Polynomial-fit slope log10(t_first) vs log10(1/eps).
        xs_pp = []
        ys_pp = []
        for row in per_eps:
            if row["t_first"] is not None:
                xs_pp.append(math.log10(1.0 / row["eps"]))
                ys_pp.append(math.log10(row["t_first"]))
        slope_poly, _ = fit_loglog(xs_pp, ys_pp)

        # Exponential-fit: log(t_first) vs (1/eps)
        xs_exp = []
        ys_exp = []
        for row in per_eps:
            if row["t_first"] is not None:
                xs_exp.append(1.0 / row["eps"])
                ys_exp.append(math.log(row["t_first"]))
        slope_exp, intercept_exp = fit_loglog(xs_exp, ys_exp)

        # Double-exp: log(log(t_first)) vs log(1/eps)
        xs_de = []
        ys_de = []
        for row in per_eps:
            if row["t_first"] is not None and row["t_first"] > math.e:
                xs_de.append(math.log(1.0 / row["eps"]))
                ys_de.append(math.log(math.log(row["t_first"])))
        slope_de, intercept_de = fit_loglog(xs_de, ys_de)

        out["per_target"][name] = {
            "per_eps": per_eps,
            "min_err_overall": min(errs),
            "fit_polynomial_slope": slope_poly,
            "fit_exp_slope": slope_exp,
            "fit_exp_intercept": intercept_exp,
            "fit_doubleexp_slope": slope_de,
            "fit_doubleexp_intercept": intercept_de,
        }

    out_json = HERE / "voronin_density.json"
    with open(out_json, "w") as fh:
        json.dump(out, fh, indent=2)

    # Markdown summary.
    md = []
    md.append("# Voronin density analysis\n")
    md.append(f"- Scan: {cfg['n_shifts']} shifts in [{cfg['t_lo']}, "
              f"{cfg['t_hi']}], boundary points {cfg['n_boundary']}, "
              f"mp.dps {cfg['mp_dps']}\n")
    md.append(f"- Disc K: |s - {cfg['S0']}| <= {cfg['R']}\n\n")
    md.append("## Per-target scaling fits\n\n")
    md.append("Polynomial model: t_first(eps) ~ (1/eps)^slope_poly\n")
    md.append("Exponential model: log(t_first(eps)) ~ slope_exp / eps\n")
    md.append("Double-exp model: log(log(t_first(eps))) ~ slope_de * log(1/eps)\n\n")
    md.append("| target | min err | slope_poly | slope_exp | slope_de "
              "| n hits eps<0.10 | n hits eps<0.05 |\n")
    md.append("|---|---|---|---|---|---|---|\n")
    for name, r in out["per_target"].items():
        n10 = next(p["n_below"] for p in r["per_eps"] if p["eps"] == 0.10)
        n05 = next(p["n_below"] for p in r["per_eps"] if p["eps"] == 0.05)
        sp = r["fit_polynomial_slope"]
        se = r["fit_exp_slope"]
        sd = r["fit_doubleexp_slope"]
        md.append(
            f"| {name} | {r['min_err_overall']:.4f} | "
            f"{sp if sp is None else f'{sp:.2f}'} | "
            f"{se if se is None else f'{se:.3f}'} | "
            f"{sd if sd is None else f'{sd:.3f}'} | "
            f"{n10} | {n05} |\n"
        )

    md.append("\n## Per-target detail\n\n")
    for name, r in out["per_target"].items():
        md.append(f"### {name}\n\n")
        md.append("| eps | n_below | density | t_first | "
                  "[10^1,10^2) | [10^2,10^3) | [10^3,10^4) | [10^4,10^5) |\n")
        md.append("|---|---|---|---|---|---|---|---|\n")
        for p in r["per_eps"]:
            tf = p["t_first"]
            tf_s = "None" if tf is None else f"{tf:.0f}"
            d = p["per_decade"]
            md.append(
                f"| {p['eps']:.3f} | {p['n_below']} | "
                f"{p['density']:.5f} | {tf_s} | "
                f"{d['[10^1,10^2)']} | {d['[10^2,10^3)']} | "
                f"{d['[10^3,10^4)']} | {d['[10^4,10^5)']} |\n"
            )
        md.append("\n")

    out_md = HERE / "voronin_density_summary.md"
    with open(out_md, "w") as fh:
        fh.writelines(md)

    print(f"Wrote {out_json}")
    print(f"Wrote {out_md}")
    print()
    for name, r in out["per_target"].items():
        print(f"{name}: min_err={r['min_err_overall']:.4f}, "
              f"slope_poly={r['fit_polynomial_slope']}, "
              f"slope_exp={r['fit_exp_slope']}, "
              f"slope_de={r['fit_doubleexp_slope']}")


if __name__ == "__main__":
    main()

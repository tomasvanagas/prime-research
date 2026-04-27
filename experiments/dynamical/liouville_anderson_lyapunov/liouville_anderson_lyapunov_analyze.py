"""Diagnostics for G1 Liouville Anderson runs.

For each (N, seeds, energies) JSON result, compute:
  (a) sign-aware z-curve and integrated z^2 across energies;
  (b) seed-to-seed L^2 distance distribution (for the Rademacher
      ensemble) and rank of the lambda-curve in that distribution;
  (c) mean/std of gamma_lambda(E) ratio to Pastur-Figotin prediction
      inside the band;
  (d) two-point Chowla-style auto-correlation of lambda for h = 1..16
      (independent of the spectral test).

Stand-alone diagnostic: run after main experiment to land structural
analysis without rerunning the heavy curves.
"""
import argparse
import json
import os

import numpy as np


def analyse_result(path):
    with open(path) as f:
        d = json.load(f)
    energies = np.asarray(d["energies"])
    g_lam = np.asarray(d["gamma_lambda"])
    g_mean = np.asarray(d["gamma_rademacher_mean"])
    g_std = np.asarray(d["gamma_rademacher_std"])
    curves = np.asarray(d["gamma_rademacher_curves"])
    pf = np.asarray(d["pastur_figotin_pred"])

    # (a) Z-curve and chi^2 across energies (df = K).
    z = (g_lam - g_mean) / np.maximum(g_std, 1e-12)
    K = len(z)
    chi2 = float((z * z).sum())
    print(f"\n=== {os.path.basename(path)} ===")
    print(f"N={d['N']}, seeds={d['n_seeds']}, K={K}")
    print(f"max |z| = {np.abs(z).max():.3f} at E = {energies[np.abs(z).argmax()]:.3f}")
    print(f"sum z^2 = {chi2:.2f}  (df = {K}; mean of chi^2_K = {K})")
    p_chi2 = float(1.0 - chi2 / max(K, 1))  # one-sided rough
    print(f"  -> chi^2 / K = {chi2 / K:.3f} (Rademacher would expect ~ 1)")

    # (b) Rank of lambda L^2 distance in seed-to-seed distribution.
    seed_l2 = np.sqrt(((curves - g_mean[None, :]) ** 2).sum(axis=1))
    lam_l2 = float(np.sqrt(((g_lam - g_mean) ** 2).sum()))
    rank = int((seed_l2 < lam_l2).sum())
    print(f"L^2 distance from Rademacher mean:")
    print(f"  lambda      : {lam_l2:.5f}")
    print(f"  Rademacher seeds: median={np.median(seed_l2):.5f}, max={seed_l2.max():.5f}, min={seed_l2.min():.5f}")
    print(f"  rank of lambda among {len(seed_l2)} seeds: {rank} (so empirical p ~ {(len(seed_l2) - rank) / len(seed_l2):.3f})")

    # (c) Pastur-Figotin agreement inside the band [-2, 2].
    band_mask = ~np.isnan(pf)
    if band_mask.sum():
        ratio_lam = (g_lam[band_mask] / pf[band_mask])
        ratio_rad = (g_mean[band_mask] / pf[band_mask])
        print(f"Pastur-Figotin ratio (mean over {band_mask.sum()} in-band energies):")
        print(f"  gamma_lambda / gamma_PF       : {ratio_lam.mean():.4f}  (std {ratio_lam.std():.4f})")
        print(f"  gamma_Rademacher_mean / gamma_PF: {ratio_rad.mean():.4f}  (std {ratio_rad.std():.4f})")

    # (d) Top energies by |z|.
    order = np.argsort(-np.abs(z))[:5]
    print("Top 5 |z| energies:")
    for i in order:
        print(f"  E={energies[i]:+.4f}  z={z[i]:+.3f}  gamma_lambda={g_lam[i]:.5f}  rad_mean={g_mean[i]:.5f}±{g_std[i]:.5f}")

    return {
        "max_abs_z": float(np.abs(z).max()),
        "argmax_z_E": float(energies[np.abs(z).argmax()]),
        "chi2_over_K": chi2 / K,
        "lam_l2_rank": rank,
        "n_seeds": len(seed_l2),
    }


def chowla_auto(N, hmax=16):
    """Empirical Chowla-pair-correlation (1/N) Σ λ(n) λ(n+h) for h=0..hmax."""
    from liouville_anderson_lyapunov import liouville_pm1
    V = liouville_pm1(N)
    out = []
    for h in range(0, hmax + 1):
        if h == 0:
            c = float((V * V).mean())
        else:
            c = float((V[:N - h] * V[h:]).mean())
        out.append(c)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--paths", nargs="+", required=True)
    ap.add_argument("--chowla-N", type=int, default=10_000_000,
                    help="Chowla auto-correlation: largest N (cheap).")
    ap.add_argument("--chowla-hmax", type=int, default=16)
    args = ap.parse_args()

    summary = {}
    for p in args.paths:
        summary[os.path.basename(p)] = analyse_result(p)

    # Independent two-point Chowla test (cheap) on the bare λ.
    print("\n=== Chowla two-point auto-correlation of lambda ===")
    print(f"N = {args.chowla_N}")
    cor = chowla_auto(args.chowla_N, args.chowla_hmax)
    print("h    |  (1/N) Σ λ(n) λ(n+h)")
    for h, c in enumerate(cor):
        print(f" {h:3d} |  {c:+.6f}")
    print()
    # Standard error: for iid Rademacher, std of the empirical mean is
    # ~ 1/sqrt(N - h). Compute a z-score for each h >= 1.
    se = 1.0 / np.sqrt(args.chowla_N - args.chowla_hmax)
    print(f"Iid Rademacher SE ~ {se:.6f}")
    print("z = c / SE for h >= 1:")
    for h, c in enumerate(cor):
        if h == 0:
            continue
        print(f" h={h:3d}  z = {c / se:+.3f}")

    # Sum so we can report a chi^2-style aggregate.
    z_h = np.array([c / se for c in cor[1:]])
    print(f"\naggregate Σ z_h^2 over h=1..{args.chowla_hmax} = {(z_h * z_h).sum():.2f}  (Rademacher chi^2_{args.chowla_hmax} ~ {args.chowla_hmax})")

    out = {
        "summaries": summary,
        "chowla": {
            "N": args.chowla_N,
            "hmax": args.chowla_hmax,
            "se": float(se),
            "correlations": cor,
            "z_per_h": [float(c / se) for c in cor[1:]],
        },
    }
    out_path = os.path.join(os.path.dirname(__file__), "analysis_summary.json")
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nwrote {out_path}")


if __name__ == "__main__":
    main()

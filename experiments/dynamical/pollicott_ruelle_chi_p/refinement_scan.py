"""
Refinement-stability scan for the χ_P-weighted Gauss-map transfer
operator. Pollicott-Ruelle resonances are characterised by stability
under discretisation refinement; spurious eigenvalues drift.

Sweeps M_grid ∈ {30, 60, 90, 120, 160} and n_max ∈ {100, 200, 400, 800}
for each of {unweighted, χ_P, λ, Λ}. Reports top-5 eigenvalues per cell
and reports the L^∞ stability gap across cells.
"""

import argparse
import json
import time
from pathlib import Path

import numpy as np

from pollicott_ruelle_chi_p import (
    transfer_operator_matrix,
    top_eigenvalues,
    chi_P_weights,
    lambda_weights,
    Lambda_weights,
)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=str, default="refinement_scan.json")
    args = ap.parse_args()

    out_path = Path(__file__).parent / args.out

    M_grid_list = [30, 60, 90, 120, 160]
    n_max_list = [100, 200, 400, 800]
    K_top = 5

    weight_names = ["unweighted", "chi_P", "lambda", "Lambda"]
    weight_funcs = {
        "unweighted": lambda nm: (np.concatenate([[0], np.ones(nm)])),
        "chi_P": chi_P_weights,
        "lambda": lambda_weights,
        "Lambda": Lambda_weights,
    }

    record = {"M_grid_list": M_grid_list, "n_max_list": n_max_list,
              "K_top": K_top, "results": {}}

    for wname in weight_names:
        record["results"][wname] = {}
        print(f"\n=== {wname} ===")
        for n_max in n_max_list:
            w = weight_funcs[wname](n_max)
            for M in M_grid_list:
                t0 = time.time()
                L, _ = transfer_operator_matrix(w, M, n_max)
                e = top_eigenvalues(L, K_top)
                el = time.time() - t0
                key = f"M{M}_n{n_max}"
                record["results"][wname][key] = {
                    "M_grid": M, "n_max": n_max,
                    "real": [float(x.real) for x in e],
                    "imag": [float(x.imag) for x in e],
                    "abs": [float(abs(x)) for x in e],
                    "elapsed_s": el,
                }
                tops = ", ".join(f"{abs(x):.5f}" for x in e[:5])
                print(f"  {key:>11s}: |λ| = [{tops}]   ({el:.1f}s)")

    # Compute stability metrics: how much does |λ_k| change when we
    # double M_grid or n_max?
    print(f"\n=== Stability across (M_grid, n_max) ===")
    for wname in weight_names:
        results = record["results"][wname]
        # For each k ∈ {0..K_top-1}, gather the |λ_k| values across all
        # (M, n_max) cells. Compute std/mean ratio.
        for k in range(K_top):
            vals = np.array([results[key]["abs"][k]
                             for key in results.keys()])
            if vals.mean() > 1e-12:
                cv = vals.std() / vals.mean()  # coefficient of variation
            else:
                cv = float('nan')
            print(f"  {wname} |λ_{k}|: mean = {vals.mean():.6f}, "
                  f"std = {vals.std():.6f}, CV = {cv:.4f}")

    # The most stringent stability metric: change between
    # (M=160, n_max=800) and (M=120, n_max=400) cells. If a resonance
    # is real, this should be tiny (<1%); if spurious, large.
    print(f"\n=== Refinement convergence (M=120,n=400 → M=160,n=800) ===")
    for wname in weight_names:
        coarse = record["results"][wname]["M120_n400"]
        fine = record["results"][wname]["M160_n800"]
        for k in range(K_top):
            cv = coarse["abs"][k]
            fv = fine["abs"][k]
            rel = abs(fv - cv) / max(fv, 1e-12)
            print(f"  {wname} λ_{k}: coarse {cv:.7f} → fine {fv:.7f} "
                  f"(rel diff {rel:.4e})")

    with open(out_path, "w") as fp:
        json.dump(record, fp, indent=2)
    print(f"\n[saved] {out_path}")


if __name__ == "__main__":
    main()

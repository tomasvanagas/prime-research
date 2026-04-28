"""Scaling analysis for D5 CTQW: track ratio_max / equilib and top-eigenvector
overlap with v_p across x ∈ {32, 64, 96, 128, 192, 256, 384, 512}.

Falsification: if ratio_max stays bounded (no super-equilibrium concentration
that grows with x), then CTQW from |1> on D_x is no better than classical
equilibrium for prime detection — same conclusion as Szegedy (D4 closure).
"""
from __future__ import annotations

import argparse
import json
import math
import time
from pathlib import Path

import numpy as np
from scipy.linalg import eigh
from sympy import primerange

from ctqw_chi_p import (divisor_graph, coprime_graph, prime_indicator,
                        ctqw_amplitude_curve_eig)


def measure_one(x: int, t_grid: np.ndarray, graph_kind: str, hamiltonian: str):
    if graph_kind == "divisor":
        A = divisor_graph(x)
    else:
        A = coprime_graph(x)
    if hamiltonian == "A":
        H = A.copy()
    else:
        D = np.diag(A.sum(axis=1))
        H = D - A

    pi_x = sum(1 for _ in primerange(2, x + 1))
    e1 = np.zeros(x); e1[0] = 1.0
    v_p = prime_indicator(x)

    lam, U = eigh(H)
    phase = np.exp(-1j * np.outer(t_grid, lam))
    P = ctqw_amplitude_curve_eig(lam, U, e1, v_p, t_grid, phase)
    max_amp = float(np.max(P))
    mean_amp = float(np.mean(P))
    equilib = pi_x / x

    overlaps = U.T @ v_p          # |<u_k | v_p>|
    a1 = U.T @ e1                 # |<u_k | 1>|
    abs_overlap_v_p = np.abs(overlaps)
    abs_overlap_e1 = np.abs(a1)
    top_op = np.argsort(abs_overlap_v_p)[::-1][:5]
    return {
        "x": x, "graph": graph_kind, "ham": hamiltonian,
        "pi_x": pi_x, "equilib": equilib,
        "max_amp": max_amp, "mean_amp": mean_amp,
        "ratio_max": max_amp / equilib, "ratio_mean": mean_amp / equilib,
        "top_overlap_v_p": float(np.max(abs_overlap_v_p)),
        "top5_v_p_overlap": [
            {"k": int(i), "lambda": float(lam[i]),
             "overlap_v_p": float(abs_overlap_v_p[i]),
             "overlap_e1": float(abs_overlap_e1[i])}
            for i in top_op
        ],
        "lam_min": float(lam[0]), "lam_max": float(lam[-1]),
        "mean_degree": float(A.sum() / x),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--xs", type=int, nargs="+",
                    default=[32, 64, 96, 128, 192, 256, 384, 512, 768])
    ap.add_argument("--T-max", type=float, default=500.0)
    ap.add_argument("--T-pts", type=int, default=5001)
    ap.add_argument("--out", type=str, default="ctqw_scaling_results.json")
    args = ap.parse_args()

    t_grid = np.linspace(0.0, args.T_max, args.T_pts)
    rows = []
    for graph in ["divisor", "coprime"]:
        for ham in ["A"]:
            for x in args.xs:
                t0 = time.time()
                r = measure_one(x, t_grid, graph, ham)
                r["wall"] = time.time() - t0
                rows.append(r)
                print(f"{graph:>8s} H={ham} x={x:>4d}  pi={r['pi_x']:>4d}  "
                      f"ratio_max={r['ratio_max']:.4f}  ratio_mean={r['ratio_mean']:.4f}  "
                      f"top_overlap={r['top_overlap_v_p']:.4f}  ({r['wall']:.1f}s)")

    # Scaling fit: ratio_max ~ a + b * (1/log x) ?
    # Top-overlap ~ c * x^{-d} ?
    div_rows = [r for r in rows if r["graph"] == "divisor" and r["ham"] == "A"]
    xs = np.array([r["x"] for r in div_rows])
    ratios = np.array([r["ratio_max"] for r in div_rows])
    overlaps = np.array([r["top_overlap_v_p"] for r in div_rows])
    log_x = np.log(xs)

    # Fit ratio = a + b / log x
    A_fit = np.vstack([np.ones_like(xs), 1.0 / log_x]).T
    coef_r, *_ = np.linalg.lstsq(A_fit, ratios, rcond=None)
    # Fit overlap ~ c * x^{-d} -> log overlap = log c - d log x
    A_fit2 = np.vstack([np.ones_like(xs), log_x]).T
    coef_o, *_ = np.linalg.lstsq(A_fit2, np.log(overlaps), rcond=None)

    fit_summary = {
        "ratio_max_fit": {
            "form": "a + b / log(x)",
            "a": float(coef_r[0]),
            "b": float(coef_r[1]),
        },
        "top_overlap_fit": {
            "form": "c * x^{-d}",
            "log_c": float(coef_o[0]),
            "d": float(-coef_o[1]),
        },
    }
    print(f"\nFit (divisor, H=A):")
    print(f"  ratio_max = {coef_r[0]:.4f} + {coef_r[1]:+.4f} / log(x)  -> "
          f"x→∞ limit ≈ {coef_r[0]:.4f}")
    print(f"  top_overlap_v_p ≈ {math.exp(coef_o[0]):.4f} * x^{{ {coef_o[1]:+.4f} }}")

    Path(args.out).write_text(json.dumps({
        "args": vars(args), "rows": rows, "fit": fit_summary
    }, indent=2))
    print(f"\nWrote {args.out}")


if __name__ == "__main__":
    main()

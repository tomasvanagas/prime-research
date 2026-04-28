"""Supplementary scans for D5 CTQW: longer T window, Laplacian Hamiltonian
variant, and coprime graph. Reuses divisor_graph etc. from main module.
"""
from __future__ import annotations

import argparse
import json
import math
import time
from pathlib import Path

import numpy as np
from scipy.linalg import eigh

from ctqw_chi_p import (divisor_graph, coprime_graph, prime_indicator,
                        random_subset_indicator, cramer_odd_indicator,
                        ctqw_amplitude_curve_eig, random_regular_graph_adj)
from sympy import primerange


def degree_matrix(A: np.ndarray) -> np.ndarray:
    d = A.sum(axis=1)
    return np.diag(d)


def run_variant(x: int, t_grid: np.ndarray, n_seeds: int, rng: np.random.Generator,
                graph_kind: str, hamiltonian: str):
    if graph_kind == "divisor":
        A = divisor_graph(x)
    elif graph_kind == "coprime":
        A = coprime_graph(x)
    else:
        raise ValueError(graph_kind)
    if hamiltonian == "A":
        H = A.copy()
    elif hamiltonian == "L":  # Laplacian H = D - A
        H = degree_matrix(A) - A
    else:
        raise ValueError(hamiltonian)

    pi_x = sum(1 for _ in primerange(2, x + 1))
    e1 = np.zeros(x); e1[0] = 1.0
    v_p = prime_indicator(x)

    lam, U = eigh(H)
    phase_cache = np.exp(-1j * np.outer(t_grid, lam))

    P_prime = ctqw_amplitude_curve_eig(lam, U, e1, v_p, t_grid, phase_cache)
    max_prime = float(np.max(P_prime))
    tstar = float(t_grid[int(np.argmax(P_prime))])
    mean_prime = float(np.mean(P_prime))

    # C2 only — Cramér+odd-parity, the most stringent control
    rng_ = np.random.default_rng(int(rng.integers(1 << 30)))
    peaks = []
    means = []
    for s in range(n_seeds):
        v = cramer_odd_indicator(x, pi_x, rng_)
        P = ctqw_amplitude_curve_eig(lam, U, e1, v, t_grid, phase_cache)
        peaks.append(float(np.max(P)))
        means.append(float(np.mean(P)))
    peaks = np.array(peaks); means = np.array(means)
    z_peak = (max_prime - peaks.mean()) / max(peaks.std(ddof=1), 1e-300)
    z_mean = (mean_prime - means.mean()) / max(means.std(ddof=1), 1e-300)

    equilib = pi_x / x
    return {
        "x": x, "graph_kind": graph_kind, "hamiltonian": hamiltonian,
        "T_max": float(t_grid[-1]), "T_pts": len(t_grid),
        "pi_x": pi_x, "equilib": equilib,
        "max_amp": max_prime, "t_star": tstar,
        "ratio_max_over_equilib": max_prime / equilib,
        "mean_amp": mean_prime,
        "ratio_mean_over_equilib": mean_prime / equilib,
        "z_peak_C2": float(z_peak), "z_mean_C2": float(z_mean),
        "C2_peak_mean": float(peaks.mean()), "C2_peak_std": float(peaks.std(ddof=1)),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--xs", type=int, nargs="+", default=[64, 128, 256, 512])
    ap.add_argument("--T-max", type=float, default=1000.0)   # longer
    ap.add_argument("--T-pts", type=int, default=10001)
    ap.add_argument("--seeds", type=int, default=80)
    ap.add_argument("--out", type=str, default="ctqw_supp_results.json")
    ap.add_argument("--seed", type=int, default=20260427)
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    t_grid = np.linspace(0.0, args.T_max, args.T_pts)

    rows = []
    for graph in ["divisor", "coprime"]:
        for ham in ["A", "L"]:
            for x in args.xs:
                t0 = time.time()
                r = run_variant(x, t_grid, args.seeds, rng, graph, ham)
                r["wall"] = time.time() - t0
                rows.append(r)
                print(f"{graph:>8s}  H={ham}  x={x:>4d}  pi={r['pi_x']:>3d}  "
                      f"max={r['max_amp']:.4f}  ratio_max={r['ratio_max_over_equilib']:.3f}  "
                      f"t*={r['t_star']:>7.2f}  z_peak(C2)={r['z_peak_C2']:+.2f}  "
                      f"z_mean(C2)={r['z_mean_C2']:+.2f}  ({r['wall']:.1f}s)")
    Path(args.out).write_text(json.dumps({"args": vars(args), "rows": rows}, indent=2))
    print(f"\nWrote {args.out}")


if __name__ == "__main__":
    main()

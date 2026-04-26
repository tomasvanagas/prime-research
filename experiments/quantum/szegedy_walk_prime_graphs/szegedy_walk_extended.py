"""Extended Szegedy walk analysis: scaling laws and stationary distribution
prime bias.

Builds on szegedy_walk_prime_graphs.py with:
  - Wider x-sweep on each family (asymptotic confirmation).
  - Exact stationary prime mass; comparison to Mertens prediction pi^2/6.
  - Power-law fit on Cayley spectral gap to high precision.
  - Inspection of WHICH eigenvectors carry primality signal.

Output: scaling.json, stationary_bias.json, summary.log
"""
from __future__ import annotations

import math
import time
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh
from sympy import isprime, primerange

# Reuse builders from main module by importing in same dir
import sys
sys.path.insert(0, str(Path(__file__).parent))
from szegedy_walk_prime_graphs import (
    divisor_graph, coprime_graph, cayley_unit_graph,
    lazy_walk_transition, discriminant_matrix, spectral_gap_lazy,
)


def stationary_prime_mass(A: np.ndarray, vertex_list: list[int]) -> dict:
    """Stationary distribution of the lazy random walk = degree / volume.

    Returns mass on prime vertices, mass under uniform, ratio.
    """
    deg = A.sum(axis=1)
    vol = deg.sum()
    pi_st = deg / vol
    chi_P = np.array([1.0 if isprime(v) else 0.0 for v in vertex_list])
    n_p = int(chi_P.sum())
    n = len(vertex_list)
    stationary_p = float((pi_st * chi_P).sum())
    uniform_p = n_p / n if n > 0 else 0.0
    ratio = stationary_p / uniform_p if uniform_p > 0 else float("nan")
    return {
        "n_vertices": n,
        "n_primes": n_p,
        "uniform_prime_mass": uniform_p,
        "stationary_prime_mass": stationary_p,
        "ratio": ratio,
    }


def cayley_scaling_sweep(out: dict):
    """Run Cayley graph (Z/NZ)* with primes N (so |V| = N-1) over wide N.

    Goal: fit log(gap) vs log(N) to high precision.
    """
    sweep = [31, 61, 89, 127, 167, 211, 257, 307, 373, 449, 541, 673, 809, 1009]
    rows = []
    for N in sweep:
        A, vlist = cayley_unit_graph(N, [2, 3, 5])
        P, deg = lazy_walk_transition(A)
        D = discriminant_matrix(P)
        gap, eigvals, eigvecs = spectral_gap_lazy(D)
        # phase gap: 2 arcsin(sqrt(gap)) -- Szegedy spectral gap
        gap_clip = max(0.0, min(gap, 1.0))
        phase_gap = 2.0 * math.asin(math.sqrt(gap_clip))
        smix = 1.0 / phase_gap if phase_gap > 0 else float("inf")
        rows.append({
            "N": N,
            "n_vertices": A.shape[0],
            "spectral_gap": float(gap),
            "phase_gap": float(phase_gap),
            "szegedy_mixing": float(smix),
            "log_N": float(math.log(N)),
            "log_gap": float(math.log(max(gap, 1e-12))),
            "log_smix": float(math.log(max(smix, 1e-12))),
        })
    out["cayley_sweep"] = rows
    # fit log(gap) ~ a * log(N) + c
    arr_logN = np.array([r["log_N"] for r in rows])
    arr_loggap = np.array([r["log_gap"] for r in rows])
    arr_logsmix = np.array([r["log_smix"] for r in rows])
    a_gap, c_gap = np.polyfit(arr_logN, arr_loggap, 1)
    a_smix, c_smix = np.polyfit(arr_logN, arr_logsmix, 1)
    out["cayley_fit"] = {
        "gap_exponent": float(a_gap),
        "gap_const": float(c_gap),
        "smix_exponent": float(a_smix),
        "smix_const": float(c_smix),
        "interpretation": (
            f"gap(N) ~ N^{a_gap:.3f}; "
            f"szegedy_mixing(N) ~ N^{a_smix:.3f} -- {'POLYLOG' if abs(a_smix) < 0.05 else 'POLY(N)'}"
        ),
    }


def coprime_asymptotic_sweep(out: dict):
    """Coprime graph C_x at larger x; check if gap stays constant."""
    sweep = [50, 100, 200, 350, 500, 700, 1000]
    rows = []
    for x in sweep:
        A = coprime_graph(x)
        vlist = list(range(1, x + 1))
        P, deg = lazy_walk_transition(A)
        D = discriminant_matrix(P)
        gap, eigvals, eigvecs = spectral_gap_lazy(D)
        stat = stationary_prime_mass(A, vlist)
        # Asymptotic Mertens ratio: 6/pi^2 ~ 0.6079; reciprocal ~ 1.6449
        merten = math.pi ** 2 / 6.0
        rows.append({
            "x": x,
            "spectral_gap": float(gap),
            "phase_gap": 2 * math.asin(math.sqrt(max(0.0, min(gap, 1.0)))),
            "stationary_prime_ratio": stat["ratio"],
            "mertens_prediction": merten,
            "deviation": stat["ratio"] - merten,
            "stat_prime_mass": stat["stationary_prime_mass"],
            "uniform_prime_mass": stat["uniform_prime_mass"],
        })
    out["coprime_sweep"] = rows


def divisor_asymptotic_sweep(out: dict):
    """Divisor graph D_x at larger x; check if gap stays constant."""
    sweep = [50, 100, 200, 350, 500, 700, 1000]
    rows = []
    for x in sweep:
        A = divisor_graph(x)
        vlist = list(range(1, x + 1))
        P, deg = lazy_walk_transition(A)
        D = discriminant_matrix(P)
        gap, eigvals, eigvecs = spectral_gap_lazy(D)
        stat = stationary_prime_mass(A, vlist)
        rows.append({
            "x": x,
            "spectral_gap": float(gap),
            "phase_gap": 2 * math.asin(math.sqrt(max(0.0, min(gap, 1.0)))),
            "stationary_prime_ratio": stat["ratio"],
            "stat_prime_mass": stat["stationary_prime_mass"],
            "uniform_prime_mass": stat["uniform_prime_mass"],
        })
    out["divisor_sweep"] = rows


def secondary_eigenvector_primality(out: dict):
    """For each family at moderate size, find which eigenvector index has
    HIGHEST prime mass and which has LOWEST. Looking for non-trivial
    primality-localized eigenvectors beyond the trivial Perron one.
    """
    cases = [
        ("coprime", 250, None),
        ("divisor", 250, None),
        ("cayley_307", 307, [2, 3, 5]),
    ]
    rows = []
    for name, x, gens in cases:
        if name == "coprime":
            A = coprime_graph(x)
            vlist = list(range(1, x + 1))
        elif name == "divisor":
            A = divisor_graph(x)
            vlist = list(range(1, x + 1))
        else:
            A, vlist = cayley_unit_graph(x, gens)
        P, deg = lazy_walk_transition(A)
        D = discriminant_matrix(P)
        gap, eigvals, eigvecs = spectral_gap_lazy(D)
        chi_P = np.array([1.0 if isprime(v) else 0.0 for v in vlist])
        n = len(vlist)
        n_p = int(chi_P.sum())
        expected = n_p / n if n > 0 else 0.0
        # Mass on primes of |v_k|^2
        masses = []
        for k in range(min(50, eigvecs.shape[1])):
            v = eigvecs[:, k]
            prob = (v * v) / (v * v).sum()
            masses.append(float((prob * chi_P).sum()))
        masses_arr = np.array(masses)
        ratios = masses_arr / max(expected, 1e-12)
        argmax = int(np.argmax(ratios))
        argmin = int(np.argmin(ratios))
        rows.append({
            "case": name,
            "x": x,
            "expected_uniform_mass": expected,
            "eig0_mass": masses[0],
            "eig0_ratio": ratios[0],
            "max_prime_eig_index": argmax,
            "max_prime_eig_eigenvalue": float(eigvals[argmax]),
            "max_prime_eig_ratio": float(ratios[argmax]),
            "min_prime_eig_index": argmin,
            "min_prime_eig_eigenvalue": float(eigvals[argmin]),
            "min_prime_eig_ratio": float(ratios[argmin]),
            "n_eigenvectors_with_ratio_above_2": int(np.sum(ratios > 2.0)),
        })
    out["secondary_eigvecs"] = rows


def main():
    out_dir = Path(__file__).parent
    out: dict = {}
    log_lines: list[str] = []

    def log(msg: str):
        print(msg, flush=True)
        log_lines.append(msg)

    log("=" * 78)
    log("Extended Szegedy walk analysis")
    log("=" * 78)

    log("\n## Cayley sweep (N prime, gens={2,3,5})")
    cayley_scaling_sweep(out)
    for r in out["cayley_sweep"]:
        log(f"  N={r['N']:5d}  V={r['n_vertices']:5d}  gap={r['spectral_gap']:.5f}  "
            f"smix~{r['szegedy_mixing']:7.3f}")
    log(f"  fit: {out['cayley_fit']['interpretation']}")

    log("\n## Coprime asymptotic sweep")
    coprime_asymptotic_sweep(out)
    log(f"  Mertens prediction: pi^2/6 = {math.pi**2/6:.4f}")
    for r in out["coprime_sweep"]:
        log(f"  x={r['x']:5d}  gap={r['spectral_gap']:.5f}  "
            f"prime_ratio={r['stationary_prime_ratio']:.4f}  "
            f"dev={r['deviation']:+.4f}")

    log("\n## Divisor asymptotic sweep")
    divisor_asymptotic_sweep(out)
    for r in out["divisor_sweep"]:
        log(f"  x={r['x']:5d}  gap={r['spectral_gap']:.5f}  "
            f"prime_ratio={r['stationary_prime_ratio']:.4f}")

    log("\n## Secondary eigenvector primality search")
    secondary_eigenvector_primality(out)
    for r in out["secondary_eigvecs"]:
        log(f"  [{r['case']:>10s}] expected={r['expected_uniform_mass']:.4f}  "
            f"eig0_ratio={r['eig0_ratio']:.3f}  "
            f"max_eig#{r['max_prime_eig_index']} (ratio={r['max_prime_eig_ratio']:.3f})  "
            f"#eig with ratio>2: {r['n_eigenvectors_with_ratio_above_2']}")

    with open(out_dir / "scaling.json", "w") as f:
        json.dump(out, f, indent=2)
    with open(out_dir / "extended.log", "w") as f:
        f.write("\n".join(log_lines))


if __name__ == "__main__":
    main()

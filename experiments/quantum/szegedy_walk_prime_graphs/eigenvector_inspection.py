"""Inspect WHICH primes the divisor graph's high-ratio eigenvectors localize on.

Hypothesis: the >2x prime-mass eigenvectors are zero-modes of the star
sub-structure (vertex 1 as hub + primes p in (x/2, x] as leaves). If so,
they are a trivial graph-theoretic artifact of the divisor graph
construction (which already uses primality to define edges).
"""
from __future__ import annotations

import math
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh
from sympy import isprime

import sys
sys.path.insert(0, str(Path(__file__).parent))
from szegedy_walk_prime_graphs import (
    divisor_graph, lazy_walk_transition, discriminant_matrix, spectral_gap_lazy,
)


def main():
    out: dict = {}
    log_lines: list[str] = []

    def log(s):
        print(s, flush=True)
        log_lines.append(s)

    x = 250
    A = divisor_graph(x)
    vlist = list(range(1, x + 1))
    P, deg = lazy_walk_transition(A)
    D = discriminant_matrix(P)
    gap, eigvals, eigvecs = spectral_gap_lazy(D)

    chi_P = np.array([1.0 if isprime(v) else 0.0 for v in vlist])
    chi_P_large = np.array([1.0 if (isprime(v) and v > x / 2) else 0.0 for v in vlist])
    chi_P_small = np.array([1.0 if (isprime(v) and v <= x / 2) else 0.0 for v in vlist])

    expected = chi_P.sum() / len(vlist)
    n_large_primes = int(chi_P_large.sum())
    n_small_primes = int(chi_P_small.sum())
    log(f"x={x}, total primes={int(chi_P.sum())}, large primes p in (x/2, x]: {n_large_primes}, small: {n_small_primes}")

    rows = []
    for k in range(min(50, eigvecs.shape[1])):
        v = eigvecs[:, k]
        prob = (v * v) / (v * v).sum()
        m_all = float((prob * chi_P).sum())
        m_large = float((prob * chi_P_large).sum())
        m_small = float((prob * chi_P_small).sum())
        rows.append({
            "k": k,
            "eigenvalue": float(eigvals[k]),
            "ratio_all_primes": m_all / expected,
            "mass_large_primes": m_large,
            "mass_small_primes": m_small,
            "mass_total": float(prob.sum()),
        })

    log("\n  k | eig    | ratio | m_large | m_small | nonprime_mass")
    for r in rows[:50]:
        if r["ratio_all_primes"] > 1.5:
            nonp = 1.0 - r["mass_large_primes"] - r["mass_small_primes"]
            log(f"  {r['k']:2d} | {r['eigenvalue']:+.3f} | {r['ratio_all_primes']:5.2f} | "
                f"{r['mass_large_primes']:.3f}  | {r['mass_small_primes']:.3f}  | {nonp:.3f}")

    # Identify zero-eigenvalue cluster
    zero_count = int(np.sum(np.abs(eigvals) < 0.01))
    log(f"\n  #eigenvalues with |eig|<0.01: {zero_count}")

    # Statistics of nontrivial high-ratio eigenvectors
    high_ratio_indices = [r["k"] for r in rows if r["ratio_all_primes"] > 2.0]
    if high_ratio_indices:
        large_mass_mean = np.mean([rows[k]["mass_large_primes"] for k in high_ratio_indices])
        small_mass_mean = np.mean([rows[k]["mass_small_primes"] for k in high_ratio_indices])
        log(f"  high-ratio eigvecs: mean large-prime mass = {large_mass_mean:.3f}, "
            f"mean small-prime mass = {small_mass_mean:.3f}")
        log(f"  --> {'YES' if large_mass_mean > 0.5 else 'NO'}: "
            f"high-ratio eigvecs concentrate on LARGE PRIMES (graph leaves).")

    out["divisor_eigvecs_inspection"] = {
        "x": x,
        "n_large_primes": n_large_primes,
        "n_small_primes": n_small_primes,
        "expected_uniform_ratio_baseline": expected,
        "rows": rows,
        "high_ratio_indices": high_ratio_indices,
    }

    with open(Path(__file__).parent / "eigvec_inspection.json", "w") as f:
        json.dump(out, f, indent=2)
    with open(Path(__file__).parent / "eigvec_inspection.log", "w") as f:
        f.write("\n".join(log_lines))


if __name__ == "__main__":
    main()

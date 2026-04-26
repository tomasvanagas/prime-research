"""Check: do divisor graph "high-prime-ratio" eigenvectors localize on
DEGREE CLASS, with primes being incidental occupants of a degree class?

Hypothesis: eigenvectors with high prime-mass really localize on
{vertices of degree d}, where degree-2 vertices include both "primes
p with x/3 < p <= x/2" and "small prime squares p^2".
"""
from __future__ import annotations

import math
import json
from pathlib import Path

import numpy as np
from sympy import isprime, factorint
import sys
sys.path.insert(0, str(Path(__file__).parent))
from szegedy_walk_prime_graphs import (
    divisor_graph, lazy_walk_transition, discriminant_matrix, spectral_gap_lazy,
)


def main():
    log_lines: list[str] = []

    def log(s):
        print(s, flush=True)
        log_lines.append(s)

    x = 250
    A = divisor_graph(x)
    deg = A.sum(axis=1).astype(int)
    P, _ = lazy_walk_transition(A)
    D = discriminant_matrix(P)
    gap, eigvals, eigvecs = spectral_gap_lazy(D)
    vlist = list(range(1, x + 1))
    chi_P = np.array([1.0 if isprime(v) else 0.0 for v in vlist])

    # For each high-ratio eigenvector, dump the vertices it concentrates on.
    ks = [10, 14, 19, 25]
    for k in ks:
        v = eigvecs[:, k]
        prob = (v * v) / (v * v).sum()
        log(f"\n--- k={k}, eigenvalue={eigvals[k]:.4f} ---")
        # Top 12 vertices by |v|^2
        top = np.argsort(prob)[::-1][:15]
        for idx in top:
            n = idx + 1
            d = deg[idx]
            is_p = "PRIME" if chi_P[idx] else ""
            log(f"  vertex={n:3d}  prob={prob[idx]:.4f}  deg={d:3d}  "
                f"factorization={dict(factorint(n))}  {is_p}")

    # Now compute correlation of |v|^2 with the indicator chi_{deg=d} for
    # various d. If the eigenvector localizes on a degree class, the
    # correlation will be much higher than chi_P.
    log("\n--- Degree-class localization check ---")
    log("For k=14 (high prime ratio), what degree classes does it concentrate on?")
    v14 = eigvecs[:, 14]
    p14 = (v14 * v14) / (v14 * v14).sum()
    log(f"  Mass on each degree class:")
    for d_class in range(1, 12):
        chi_d = (deg == d_class).astype(float)
        n_d = int(chi_d.sum())
        if n_d == 0:
            continue
        m_d = float((p14 * chi_d).sum())
        ratio = m_d / (n_d / x)
        log(f"  deg={d_class:2d}, n_vertices_with_this_deg={n_d:3d}, "
            f"mass={m_d:.4f}, ratio_vs_uniform={ratio:.2f}")

    # Compute deg-class for primes
    log("\n--- Degree distribution of primes vs all vertices ---")
    for d_class in range(1, 8):
        n_d_total = int(np.sum(deg == d_class))
        n_d_primes = int(np.sum((deg == d_class) & (chi_P > 0)))
        log(f"  deg={d_class}: total={n_d_total}, primes={n_d_primes}, "
            f"frac_prime={n_d_primes/max(n_d_total,1):.2f}")

    with open(Path(__file__).parent / "degree_class.log", "w") as f:
        f.write("\n".join(log_lines))


if __name__ == "__main__":
    main()

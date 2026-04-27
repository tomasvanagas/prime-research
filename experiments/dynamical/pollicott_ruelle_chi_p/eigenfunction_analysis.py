"""
Eigenfunction analysis for the χ_P-weighted Gauss-map transfer operator.

Goals:
 1. Extract the leading eigenvector (the χ_P-invariant density analogue)
    and characterise its arithmetic content vs the unweighted Gauss
    measure 1/((1+x) log 2).
 2. Compare χ_P, λ, Λ leading eigenvectors against each other and against
    Cramér-model controls.
 3. Test the **trace identity**: tr(L_h^k) = Σ_n (cycle weights) ought to
    factor through Σ_p 1/p^{2k} (prime zeta) for h = χ_P.
"""

import argparse
import json
from pathlib import Path

import numpy as np
from sympy import primerange, primepi

from pollicott_ruelle_chi_p import (
    transfer_operator_matrix,
    eigenpairs_top,
    chi_P_weights,
    lambda_weights,
    Lambda_weights,
    chebyshev_lobatto,
    cramer_odd_weights,
)


def gauss_measure_density(x):
    """Unweighted Gauss-Kuzmin invariant density 1/((1+x) log 2)."""
    return 1.0 / ((1.0 + x) * np.log(2.0))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--M_grid", type=int, default=120)
    ap.add_argument("--n_max", type=int, default=400)
    ap.add_argument("--out", type=str, default="eigenfunction_analysis.json")
    args = ap.parse_args()
    out = Path(__file__).parent / args.out

    M = args.M_grid
    nm = args.n_max

    # === Build operators ===
    weights = {
        "unweighted": np.concatenate([[0.0], np.ones(nm)]),
        "chi_P": chi_P_weights(nm),
        "lambda": lambda_weights(nm),
        "Lambda": Lambda_weights(nm),
    }

    record = {"M_grid": M, "n_max": nm, "K_top": 5, "results": {}}

    x_nodes, _ = chebyshev_lobatto(M)
    record["x_nodes"] = x_nodes.tolist()

    for name, w in weights.items():
        L, _ = transfer_operator_matrix(w, M, nm)
        e, V = eigenpairs_top(L, K=5)
        # leading eigenvector v_0
        v0 = V[:, 0].real
        # normalise to unit L^2 with positive integral
        if v0.sum() < 0:
            v0 = -v0
        v0 = v0 / np.linalg.norm(v0)
        rec = {
            "eigenvalues_real": [float(z.real) for z in e],
            "eigenvalues_imag": [float(z.imag) for z in e],
            "eigenvalues_abs": [float(abs(z)) for z in e],
            "leading_eigvec_v0": v0.tolist(),
        }

        # Compare to unweighted Gauss density
        gm = gauss_measure_density(x_nodes)
        gm = gm / np.linalg.norm(gm)
        rec["v0_overlap_with_gauss_density"] = float(np.dot(v0, gm))

        # Trace identity check: tr(L_h^1) = Σ over fixed points
        # of the Gauss map of h(n)·|T'(y)|^{-1}. Fixed points of T(x)=1/x-n
        # are y_n = (-n + sqrt(n^2 + 4))/2 in [0,1]. There's one in each
        # cylinder [1/(n+1), 1/n].
        traces = []
        Lk = L.copy()
        for k in range(1, 7):
            tr = float(np.trace(Lk))
            traces.append(tr)
            if k < 6:
                Lk = Lk @ L
        rec["traces_L_to_k"] = traces

        # Predicted leading-order trace for h-weighted operator:
        # tr(L_h) ≈ Σ_n h(n) · y_n^2 / (1 - y_n^2 · |T'|) where T'(y_n) = -1/y_n^2.
        # Fixed-point trace formula: tr(L_h) = Σ_n h(n) / (n^2 + n - 1).
        # (This follows from y_n satisfying y_n(y_n + n) = 1, so 1/y_n = y_n + n,
        #  and at fixed point Jacobian = (x+n)^2 = 1/y_n^2; L_h(δ_{y_n}) = h(n)/y_n^2 · δ_{y_n};
        #  trace contribution at this fixed point is h(n) y_n^2 / (1 - y_n^2 / (1/y_n^2))
        #  = h(n) y_n^2 / (1 - y_n^4) for y_n ≠ ±1; but for periodic orbits the formula is
        #  tr(L_h^k) = Σ_{periodic-k orbits} h-product / |T'^k - 1| .)
        rec["pred_trace_fixed_pt_sum"] = float(sum(
            w[n] / (n**2 + n - 1) for n in range(1, nm+1) if w[n] != 0
        ))

        # χ_P-specific arithmetic predictions
        if name == "chi_P":
            # Σ_p 1/p^2 (prime zeta)
            P2 = sum(1.0 / p**2 for p in primerange(2, nm + 1))
            rec["prime_zeta_2"] = P2
            # Σ_p 1/(p^2 + p - 1)  (fixed-pt trace)
            P_fixed = sum(1.0 / (p**2 + p - 1) for p in primerange(2, nm + 1))
            rec["prime_fixed_pt_sum"] = P_fixed
            # Σ_p 1/(p(p+1))  (kernel sup over x in [0,1])
            P_inv_pp1 = sum(1.0 / (p * (p + 1)) for p in primerange(2, nm + 1))
            rec["prime_kernel_max_sum"] = P_inv_pp1

        record["results"][name] = rec
        e0 = e[0].real
        print(f"[{name}] λ_0 = {e0:.6f}, |λ_1|/|λ_0| = "
              f"{abs(e[1])/abs(e[0]):.4f}, "
              f"v0·(1/((1+x)log2)) overlap = {rec['v0_overlap_with_gauss_density']:.4f}, "
              f"tr(L^1) = {rec['traces_L_to_k'][0]:.6f}, "
              f"pred Σh(n)/(n²+n-1) = {rec['pred_trace_fixed_pt_sum']:.6f}")

    # χ_P-specific arithmetic structural test
    chi_rec = record["results"]["chi_P"]
    print(f"\n[χ_P arithmetic predictions]")
    print(f"  λ_0(χ_P)              = {chi_rec['eigenvalues_real'][0]:.7f}")
    print(f"  Σ_p 1/p²              = {chi_rec['prime_zeta_2']:.7f}")
    print(f"  Σ_p 1/(p²+p-1)        = {chi_rec['prime_fixed_pt_sum']:.7f}")
    print(f"  Σ_p 1/(p(p+1))        = {chi_rec['prime_kernel_max_sum']:.7f}")
    print(f"  tr(L_χP)              = {chi_rec['traces_L_to_k'][0]:.7f}")
    print(f"  pred Σ_p 1/(p²+p-1)   = {chi_rec['pred_trace_fixed_pt_sum']:.7f}")
    print(f"  λ_0 - Σ_p 1/(p(p+1))  = "
          f"{chi_rec['eigenvalues_real'][0] - chi_rec['prime_kernel_max_sum']:+.7f}")
    print(f"  λ_0 / Σ_p 1/p²        = "
          f"{chi_rec['eigenvalues_real'][0] / chi_rec['prime_zeta_2']:.7f}")

    # Now: does the χ_P-eigenvector resemble the unweighted Gauss density?
    # If yes → χ_P spectrum is just a damped version of Gauss-Kuzmin-Wirsing.
    # If no  → χ_P induces a different invariant density (arithmetic content).
    print(f"\n[Eigenvector overlap with unweighted Gauss density 1/((1+x)log2)]")
    for name in ["unweighted", "chi_P", "lambda", "Lambda"]:
        ov = record["results"][name]["v0_overlap_with_gauss_density"]
        print(f"  {name}: ⟨v_0, gauss⟩ = {ov:+.5f}  (1.0 = identical, 0 = orthogonal)")

    with open(out, "w") as fp:
        json.dump(record, fp, indent=2)
    print(f"\n[saved] {out}")


if __name__ == "__main__":
    main()

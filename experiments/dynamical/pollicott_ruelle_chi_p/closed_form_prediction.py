"""
CLOSED-FORM PREDICTION: λ_0^h ≈ Σ_n h(n) a_n where
  a_n = T_n / ||g||²,  g(x) = 1/((1+x) log 2),
  T_n = ∫₀¹ g(x) g(1/(x+n)) / (x+n)² dx
      = (1/(log 2)²) ∫₀¹ dx / [(1+x)(x+n)(x+n+1)]    (closed form)

The a_n's are derived from the Rayleigh quotient ⟨g, L_h g⟩ / ⟨g, g⟩
under the empirical observation that the leading eigenvector of L_h is
itself (to ~0.15% relative error in cosine similarity) the Gauss
density g.

Predictions:
  Σ_n a_n              = 1                    (consistency)
  Σ_n χ_P(n) a_n      ≈ |λ_0|^{χ_P}          (Cramér-like closure)
  Σ_n λ(n) a_n        ≈ λ_0^{λ}              (signed)
  Σ_n Λ(n) a_n        ≈ |λ_0|^{Λ}
"""

import argparse
import json
from pathlib import Path

import numpy as np
from sympy import primerange, isprime, factorint


def a_n_closed_form(n_max):
    """
    Exact a_n via partial fractions on (1+x)(x+n)(x+n+1).
    Returns array a[0..n_max], with a[0] = 0.
    """
    a = np.zeros(n_max + 1, dtype=np.float64)
    log2 = np.log(2.0)
    norm_g_sq = 0.5 / log2**2  # ||g||²

    # n=1 special case
    # T_1 = (1/log² 2) * (-ln 2 + 1/2 + ln(3/2))
    T1 = (1.0 / log2**2) * (-np.log(2.0) + 0.5 + np.log(1.5))
    a[1] = T1 / norm_g_sq

    # n >= 2 general formula
    for n in range(2, n_max + 1):
        # T_n = (1/log² 2) * [ ln 2 / (n(n-1))
        #                     - ln((n+1)/n) / (n-1)
        #                     + ln((n+2)/(n+1)) / n ]
        Tn = (1.0 / log2**2) * (
            np.log(2.0) / (n * (n - 1))
            - np.log((n + 1.0) / n) / (n - 1)
            + np.log((n + 2.0) / (n + 1)) / n
        )
        a[n] = Tn / norm_g_sq

    return a


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n_max", type=int, default=400)
    ap.add_argument("--out", type=str, default="closed_form_prediction.json")
    args = ap.parse_args()
    out_path = Path(__file__).parent / args.out

    n_max = args.n_max
    a = a_n_closed_form(n_max)

    print(f"[closed-form coefficients a_n via Gauss-Kuzmin Rayleigh quotient]")
    print(f"  n_max = {n_max}")
    print(f"  Σ a_n (n=1..{n_max}) = {a.sum():.10f}  (target: 1.000…)")
    # Check the tail: a_n ≈ ln 2 / (log² 2 · ||g||²) · 1/n³ for large n
    print(f"  a_1 = {a[1]:.6f}")
    print(f"  a_2 = {a[2]:.6f}")
    print(f"  a_3 = {a[3]:.6f}")
    print(f"  a_5 = {a[5]:.6f}")
    print(f"  a_10 = {a[10]:.6f}")
    print(f"  a_100 = {a[100]:.6e}")

    # === Predictions ===
    chi_P_pred = sum(a[p] for p in primerange(2, n_max + 1))
    print(f"\n[χ_P prediction]")
    print(f"  Σ_p prime ≤ {n_max} a_p = {chi_P_pred:.7f}")
    print(f"  vs measured |λ_0|^{{χ_P}} = 0.3596100  (M_grid=120, n_max=400)")
    print(f"  agreement: {chi_P_pred / 0.3596100:.5f}× ({(chi_P_pred - 0.3596100)*100:+.4f}% absolute)")

    # Liouville
    lam_arr = np.zeros(n_max + 1)
    for n in range(1, n_max + 1):
        Omega = sum(factorint(n).values())
        lam_arr[n] = (-1) ** Omega
    lam_pred = float(np.sum(lam_arr * a))
    print(f"\n[λ prediction]")
    print(f"  Σ_n λ(n) a_n = {lam_pred:.7f}")
    print(f"  vs measured λ_0^λ = +0.0902353")
    print(f"  agreement: {lam_pred / 0.0902353:.5f}×")

    # von Mangoldt
    Lam_arr = np.zeros(n_max + 1)
    for n in range(2, n_max + 1):
        f = factorint(n)
        if len(f) == 1:
            p = next(iter(f.keys()))
            Lam_arr[n] = float(np.log(p))
    Lam_pred = float(np.sum(Lam_arr * a))
    print(f"\n[Λ prediction]")
    print(f"  Σ_n Λ(n) a_n = {Lam_pred:.7f}")
    print(f"  vs measured |λ_0|^{{Λ}} = 0.4968265")
    print(f"  agreement: {Lam_pred / 0.4968265:.5f}×")

    # Cramér expectation: E[Σ_n Bern(1/log n) a_n] = Σ_n a_n / log n
    cra_expectation = sum(a[n] / np.log(n) for n in range(2, n_max + 1))
    print(f"\n[Cramér expectation]")
    print(f"  E[χ_P_Cramer |λ_0|] = Σ_n (1/log n) a_n = {cra_expectation:.7f}")
    print(f"  vs χ_P measured 0.3596100")
    print(f"  Cramér-χ_P agreement: {chi_P_pred / cra_expectation:.5f}×")

    # Save
    record = {
        "n_max": n_max,
        "a_n_sum": float(a.sum()),
        "a_n_first_30": [float(x) for x in a[:30]],
        "predictions": {
            "chi_P": float(chi_P_pred),
            "chi_P_measured_M120_N400": 0.3596100,
            "chi_P_relative_error": float((chi_P_pred - 0.3596100) / 0.3596100),
            "lambda": float(lam_pred),
            "lambda_measured_M120_N400": 0.0902353,
            "lambda_relative_error": float((lam_pred - 0.0902353) / 0.0902353),
            "Lambda": float(Lam_pred),
            "Lambda_measured_M120_N400": 0.4968265,
            "Lambda_relative_error": float((Lam_pred - 0.4968265) / 0.4968265),
            "cramer_expectation": float(cra_expectation),
        },
    }
    with open(out_path, "w") as fp:
        json.dump(record, fp, indent=2)
    print(f"\n[saved] {out_path}")


if __name__ == "__main__":
    main()

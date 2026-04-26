"""
Analyse the Gowers U^k experiment (gowers_uk_chi_p.py output) and
compare to the Hardy-Littlewood box singular series prediction.

Defines Q^k(f) := ||f||_{U^k}^{2^k} / E[f]^{2^k} (the "normalized" Gowers norm).

Hardy-Littlewood predicts:
  Q^k(chi_P) -> S_k (singular series of {0,1}^k cube) as N -> infinity.
  Q^k(chi_P_W,b) -> S^{(W)}_k = product_{p does not divide W} factor_p
                                  -> 1 as W -> infinity.
  Q^k(Bernoulli at p)            ≈ 1 + O(1/N).
"""

from __future__ import annotations
import argparse
import json
import math
import sys
from typing import Dict, List

from hardy_littlewood_box import singular_series_box, alpha_p


def Qk(rec: dict, k: int) -> float:
    """Q^k(f) = ||f||_{U^k}^{2^k} / mean^{2^k}."""
    p = rec["mean"]
    if k == 2:
        return rec["u2_pow4"] / (p ** 4)
    elif k == 3:
        if "u3_pow8" not in rec:
            return float("nan")
        return rec["u3_pow8"] / (p ** 8)
    else:
        raise ValueError(f"Unsupported k={k}")


def mean_std(xs: List[float]) -> tuple[float, float]:
    if not xs:
        return float("nan"), float("nan")
    m = sum(xs) / len(xs)
    if len(xs) < 2:
        return m, 0.0
    v = sum((x - m) ** 2 for x in xs) / (len(xs) - 1)
    return m, math.sqrt(v)


def w_trick_singular_series(W: int, k: int, P_max: int = 100) -> float:
    """
    S^{(W)}_k = product over primes p NOT dividing W of factor_p.
    """
    from sympy import primerange, isprime
    log_S = 0.0
    for p in primerange(2, P_max + 1):
        if W % p == 0:
            continue
        a = alpha_p(p, k)
        denom = (1 - 1.0 / p) ** (2 ** k)
        log_S += math.log(a / denom)
    return math.exp(log_S)


def report_run(N_data: dict, S2: float, S3: float, by_W: dict[int, dict]) -> str:
    out = []
    N = N_data["N"]
    out.append(f"\n=== N = {N} ===")
    out.append(f"chi_P density p = {N_data['p_chi']:.6f}, expected pi(N)/N = {N_data['p_chi']:.6f}")

    # chi_P
    chi = N_data["chi_P"]
    bern = N_data["bernoulli_chi_density"]
    Q2_chi = Qk(chi, 2)
    Q2_bern = [Qk(r, 2) for r in bern]
    Q2_bern_m, Q2_bern_s = mean_std(Q2_bern)

    out.append(f"\n  Q^2(chi_P)        = {Q2_chi:.4f}")
    out.append(f"  Q^2(Bernoulli)    = {Q2_bern_m:.4f} ± {Q2_bern_s:.4f} (n={len(Q2_bern)})")
    out.append(f"  HL prediction S_2 = {S2:.4f}")
    out.append(f"  ratio Q^2(chi)/S_2 = {Q2_chi / S2:.4f}    (1.0 = exact match)")

    if "u3_pow8" in chi:
        Q3_chi = Qk(chi, 3)
        Q3_bern = [Qk(r, 3) for r in bern]
        Q3_bern_m, Q3_bern_s = mean_std(Q3_bern)
        out.append(f"  Q^3(chi_P)        = {Q3_chi:.4f}")
        out.append(f"  Q^3(Bernoulli)    = {Q3_bern_m:.4f} ± {Q3_bern_s:.4f}")
        out.append(f"  HL prediction S_3 = {S3:.4f}")
        out.append(f"  ratio Q^3(chi)/S_3 = {Q3_chi / S3:.4f}")

    # Liouville
    if "liouville" in N_data:
        L = N_data["liouville"]
        bL = N_data["bernoulli_liouville_density"]
        Q2_L = Qk(L, 2)
        Q2_bL = [Qk(r, 2) for r in bL]
        Q2_bL_m, Q2_bL_s = mean_std(Q2_bL)
        out.append(f"\n  Liouville:")
        out.append(f"  Q^2(L)            = {Q2_L:.4f}  (predicted ~1; Liouville centered is Gowers-uniform)")
        out.append(f"  Q^2(Bernoulli@1/2)= {Q2_bL_m:.4f} ± {Q2_bL_s:.4f}")

    # W-tricks
    for W, info in by_W.items():
        key = f"w_trick_W{W}_b1"
        if key not in N_data:
            continue
        w_data = N_data[key]
        f = w_data["f"]
        bw = w_data["bernoulli_matched"]
        Q2_w = Qk(f, 2)
        Q2_bw = [Qk(r, 2) for r in bw]
        Q2_bw_m, Q2_bw_s = mean_std(Q2_bw)
        S2_W = info["S2"]
        out.append(f"\n  W-trick W={W}:")
        out.append(f"  Q^2(chi_W,b)      = {Q2_w:.4f}")
        out.append(f"  Q^2(Bernoulli)    = {Q2_bw_m:.4f} ± {Q2_bw_s:.4f}")
        out.append(f"  HL prediction S^(W)_2 = {S2_W:.4f}    (-> 1 as W -> infty)")
        out.append(f"  ratio Q^2(chi_W)/S^(W)_2 = {Q2_w / S2_W:.4f}")
        if "u3_pow8" in f:
            Q3_w = Qk(f, 3)
            S3_W = info["S3"]
            out.append(f"  Q^3(chi_W,b)      = {Q3_w:.4f}")
            out.append(f"  HL prediction S^(W)_3 = {S3_W:.4f}")
            out.append(f"  ratio Q^3(chi_W)/S^(W)_3 = {Q3_w / S3_W:.4f}")

    return "\n".join(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, help="JSON output from gowers_uk_chi_p.py")
    ap.add_argument("--P-max-HL", type=int, default=100)
    ap.add_argument("--out-md", default=None)
    args = ap.parse_args()

    with open(args.input) as f:
        data = json.load(f)

    # Compute HL singular series.
    print(f"Computing Hardy-Littlewood singular series (P_max = {args.P_max_HL})...")
    S2_data = singular_series_box(2, args.P_max_HL)
    S2 = S2_data["S_truncated"]
    S3_data = singular_series_box(3, min(args.P_max_HL, 47))  # alpha_p direct enum is slow above ~47
    S3 = S3_data["S_truncated"]
    print(f"  S_2 = {S2:.6f}")
    print(f"  S_3 = {S3:.6f}    (truncated at P_max <= 47 due to alpha_p cost)")
    print()

    by_W = {}
    for W in [6, 30, 210]:
        by_W[W] = {
            "S2": w_trick_singular_series(W, 2, args.P_max_HL),
            "S3": w_trick_singular_series(W, 3, min(args.P_max_HL, 47)),
        }
        print(f"  W={W}: S^(W)_2 = {by_W[W]['S2']:.6f}, S^(W)_3 = {by_W[W]['S3']:.6f}")

    print()

    chunks = []
    chunks.append("# Empirical Q^k(chi_P) vs Hardy-Littlewood prediction\n")
    chunks.append(f"S_2 = {S2:.6f}    (truncated, primes <= {args.P_max_HL})")
    chunks.append(f"S_3 = {S3:.6f}    (truncated, primes <= 47)")
    chunks.append("")
    chunks.append("W-tricked predictions:")
    for W, info in by_W.items():
        chunks.append(f"  W={W}: S^(W)_2 = {info['S2']:.4f}, S^(W)_3 = {info['S3']:.4f}")

    by_N = data["by_N"]
    for N_str, N_data in sorted(by_N.items(), key=lambda x: int(x[0])):
        chunk = report_run(N_data, S2, S3, by_W)
        chunks.append(chunk)
        print(chunk)

    if args.out_md:
        with open(args.out_md, "w") as f:
            f.write("\n".join(chunks))
        print(f"\nWrote {args.out_md}")


if __name__ == "__main__":
    main()

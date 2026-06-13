#!/usr/bin/env python3
"""Sequential major-arc stripping — does deficit close as Q grows?

For Q ∈ {1, 2, 3, 5, 7, 11, 13}, strip the major-arc subspaces V_q^prim
(spanned by characters mod q) from chi_P and re-measure log Mahler.

If Δ_∞ → 0 as Q → ∞: deficit is "sequentially major-arc-sourced" (good).
If Δ_∞ stays at -0.307: deficit is minor-arc / Polya-Carlson sourced.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
from parity_stripped_trinity import (
    chi_p_indicator,
    log_mahler,
    bernoulli_match_density,
)


def sqfree_le(Q: int) -> list[int]:
    out = []
    for q in range(1, Q + 1):
        # check squarefree
        n = q
        sqfree = True
        for p in range(2, int(n ** 0.5) + 1):
            if n % (p * p) == 0:
                sqfree = False
                break
        if sqfree:
            out.append(q)
    return out


def major_arc_strip(x: np.ndarray, q_list: list[int]) -> np.ndarray:
    """Project out V_q^prim subspaces for q in q_list (q=1 -> the constant
    vector; q>=2 -> primitive characters mod q via additive Fourier modes
    e^{2 pi i a n / q} for a coprime to q).

    Implementation: compute the additive-Fourier coefficient
        c_{a,q}(x) := (1/N) Σ_n x(n) e^{-2 pi i a n / q}
    for each (a,q) with gcd(a,q)=1, then subtract c_{a,q}·e^{2 pi i a n/q}
    (Re part for real x).
    """
    N = x.size
    n = np.arange(1, N + 1)
    out = x.astype(np.complex128).copy()
    for q in q_list:
        for a in range(0 if q == 1 else 1, q):
            if q > 1 and np.gcd(a, q) != 1:
                continue
            mode = np.exp(2j * np.pi * a * n / q)
            c = np.sum(out * np.conj(mode)) / N
            out = out - c * mode
    return out.real  # x was real, projection along real subspace


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--N_log2", type=int, default=14)
    p.add_argument("--Q_list", default="1,2,3,5,7,11,13,17,23,29")
    p.add_argument("--n_bern", type=int, default=20)
    p.add_argument("--seed", type=int, default=20260430)
    args = p.parse_args()

    rng = np.random.default_rng(args.seed)
    N = 2 ** args.N_log2
    chi = chi_p_indicator(N)
    pi_N = int(chi.sum())

    Q_values = [int(s) for s in args.Q_list.split(",")]
    log_m_chi_orig = log_mahler(chi)

    # Bernoulli baselines
    bern_log_m = []
    for _ in range(args.n_bern):
        b = bernoulli_match_density(N, pi_N, rng)
        bern_log_m.append(log_mahler(b))
    bern_mean = float(np.mean(bern_log_m))
    bern_std = float(np.std(bern_log_m, ddof=1))

    print(f"N=2^{args.N_log2}={N}, pi(N)={pi_N}")
    print(f"  log m(chi_P) = {log_m_chi_orig:.4f}")
    print(f"  BERN log m   = {bern_mean:.4f} ± {bern_std:.4f}")
    print(f"  Δ_orig       = {log_m_chi_orig - bern_mean:+.4f}")

    L2_chi = float(np.sqrt((chi * chi).sum()))
    rows = [
        {
            "Q": 0,
            "qs": [],
            "log_m_chi": log_m_chi_orig,
            "log_m_bern_mean": bern_mean,
            "log_m_bern_std": bern_std,
            "L2_chi": L2_chi,
            "L2_bern_mean": float(np.sqrt(pi_N)),
            "deficit": log_m_chi_orig - bern_mean,
            "shape_deficit": (log_m_chi_orig - np.log(L2_chi))
            - (bern_mean - np.log(np.sqrt(pi_N))),
        }
    ]

    # We strip cumulatively up to each q in Q_values.
    print(f"\nSequential cumulative stripping (chi AND BERN stripped):")
    header = "  Q   sqfree-q-list             log m(chi*)   log m(BERN*) σ_BERN*  Δ        shape-Δ  ‖chi*‖₂"
    print(header)
    # Cache stripped Bernoulli replicates (re-use for all Qs)
    rng2 = np.random.default_rng(args.seed + 1)
    bern_seqs = [
        bernoulli_match_density(N, pi_N, rng2) for _ in range(args.n_bern)
    ]
    for Q in Q_values:
        qs = sqfree_le(Q)
        x_strip = major_arc_strip(chi, qs)
        log_m_chi = log_mahler(x_strip)
        L2_chi_s = float(np.sqrt((x_strip * x_strip).sum()))
        bern_log_m_strip = []
        bern_L2_strip = []
        for b in bern_seqs:
            b_strip = major_arc_strip(b, qs)
            bern_log_m_strip.append(log_mahler(b_strip))
            bern_L2_strip.append(float(np.sqrt((b_strip * b_strip).sum())))
        bern_mean_s = float(np.mean(bern_log_m_strip))
        bern_std_s = float(np.std(bern_log_m_strip, ddof=1))
        bern_L2_s = float(np.mean(bern_L2_strip))
        deficit = log_m_chi - bern_mean_s
        shape_deficit = (log_m_chi - np.log(L2_chi_s)) - (
            bern_mean_s - np.log(bern_L2_s)
        )
        rows.append(
            {
                "Q": Q,
                "qs": qs,
                "log_m_chi": log_m_chi,
                "log_m_bern_mean": bern_mean_s,
                "log_m_bern_std": bern_std_s,
                "L2_chi": L2_chi_s,
                "L2_bern_mean": bern_L2_s,
                "deficit": deficit,
                "shape_deficit": shape_deficit,
            }
        )
        print(
            f"  Q={Q:>3} {str(qs)[:24]:<24} {log_m_chi:+8.4f}      "
            f"{bern_mean_s:+8.4f} ±{bern_std_s:.3f}  {deficit:+7.4f}  "
            f"{shape_deficit:+7.4f}  {L2_chi_s:7.2f}"
        )

    out = {
        "N": N,
        "pi_N": pi_N,
        "log_m_chi_orig": log_m_chi_orig,
        "bern_log_m_mean": bern_mean,
        "bern_log_m_std": bern_std,
        "rows": rows,
    }
    out_path = Path(__file__).parent / "sequential_strip_results.json"
    out_path.write_text(json.dumps(out, indent=2))
    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    main()

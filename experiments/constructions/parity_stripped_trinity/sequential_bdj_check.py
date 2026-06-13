#!/usr/bin/env python3
"""Sequential major-arc stripping for BDJ m_4 attribution.

Strip cumulative major arcs from chi_P standardised, build Toeplitz, measure m_4.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from parity_stripped_trinity import (
    chi_p_indicator,
    standardise,
    toeplitz_m4,
)
from sequential_strip_check import sqfree_le, major_arc_strip


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--N", type=int, default=1000)
    p.add_argument("--Q_list", default="0,1,2,3,5,7,11,13,17,23,29")
    args = p.parse_args()

    N = args.N
    chi = chi_p_indicator(N)
    pi_N = int(chi.sum())
    print(f"N={N}, pi(N)={pi_N}")
    print(f"Cumulative major-arc stripping for BDJ m_4:")
    print(f"  Q   q-list                                         m_4(chi*)   λ̃_max(chi*)   m_4 / (N/log²N)")

    log2_N = float(np.log(N))
    N_per_log2N = N / (log2_N ** 2)

    rows = []
    for Q in [int(s) for s in args.Q_list.split(",")]:
        if Q == 0:
            chi_strip = chi.copy()
            qs = []
        else:
            qs = sqfree_le(Q)
            chi_strip = major_arc_strip(chi, qs)
        try:
            eps = standardise(chi_strip)
        except ValueError:
            print(f"  Q={Q}: zero-variance after stripping, skipping")
            continue
        m4, lmax = toeplitz_m4(eps)
        ratio = m4 / N_per_log2N
        rows.append({"Q": Q, "qs": qs, "m4": m4, "lambda_max_scaled": lmax, "ratio_N_log2N": ratio})
        print(f"  Q={Q:>3} {str(qs)[:48]:<48} {m4:8.3f}    {lmax:8.3f}     {ratio:6.3f}")

    out_path = Path(__file__).parent / "sequential_bdj_results.json"
    out_path.write_text(json.dumps({"N": N, "pi_N": pi_N, "rows": rows}, indent=2))
    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    main()

"""
P5 / Thread 10 (commit, slot 2) — Galway worst-case K-constant tightening.

Slot 2 heavy-compute path: extend the zero database from K=8000 to K=20000
by combining `data/zeta_zeros_8000.txt` with the parallel chunks under
`data/zeros_chunks_slot2/chunk_*.txt`, then measure K_emp at
log10 x ∈ {5.0, 5.3, 5.5, 5.7} where the Galway-shape and Thread-7-shape
predictions diverge by factor ~3-5.

Decision rule:
  * If K_emp(log10 x = 5.5) ≤ ~19000:
    Galway-shape with c_emp ≈ const ≈ 0.20 confirmed across 1.5 decades.
  * If K=20000 still insufficient at log10 x = 5.5:
    Galway-shape c_emp ≈ const ruled out at log10 x = 5.5 → Thread-7-shape.

Imports helpers from slot1_worst_case_K.py and slot2_finegrid.py.
"""

from __future__ import annotations

import argparse
import csv
import glob
import math
import os
import random
import sys
import time
from typing import Dict, List, Tuple

import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
PR_ROOT = os.path.normpath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(PR_ROOT, "experiments", "analytic", "polylog_approx_pi"))
sys.path.insert(0, HERE)

from polylog_approx_pi import R_at_rho, riemann_R  # noqa: E402
from slot1_worst_case_K import (  # noqa: E402
    sieve_pi,
    cumulative_err_trace,
    sigma_predicted,
    K_typical,
    K_galway_grh_bound,
    find_K_emp,
)


def load_extended_zeros(K_target: int, dps: int = 18) -> List[mp.mpf]:
    """Combine zeta_zeros_8000.txt (lines 1..8000) with parallel chunks
    under data/zeros_chunks_slot2/. Each chunk file has lines '<idx> <value>',
    indices in [chunk_lo, chunk_hi]. Returns gamma_1..gamma_K_target.
    """
    mp.mp.dps = dps
    data_dir = os.path.normpath(os.path.join(HERE, "..", "..", "..", "data"))
    base_path = os.path.join(data_dir, "zeta_zeros_8000.txt")
    with open(base_path) as fh:
        base_zeros = [line.strip() for line in fh if line.strip()]
    assert len(base_zeros) >= 8000, f"base file has {len(base_zeros)} zeros"

    chunks_dir = os.path.join(data_dir, "zeros_chunks_slot2")
    chunk_files = sorted(glob.glob(os.path.join(chunks_dir, "chunk_*.txt")))
    extra: Dict[int, str] = {}
    for cf in chunk_files:
        with open(cf) as fh:
            for line in fh:
                parts = line.strip().split()
                if len(parts) == 2:
                    idx = int(parts[0])
                    extra[idx] = parts[1]
    print(f"loaded {len(base_zeros)} base zeros + {len(extra)} chunked zeros "
          f"(idx range {min(extra.keys()) if extra else '-'} .. "
          f"{max(extra.keys()) if extra else '-'})", flush=True)

    out: List[mp.mpf] = [mp.mpf(z) for z in base_zeros[:8000]]
    for k in range(8001, K_target + 1):
        if k not in extra:
            raise RuntimeError(f"missing zero idx {k} (K_target={K_target}, "
                               f"have up to {max(extra.keys()) if extra else 8000})")
        out.append(mp.mpf(extra[k]))
    assert len(out) == K_target, f"got {len(out)}, expected {K_target}"
    return out


def make_extended_K_milestones(K_max: int) -> List[int]:
    """Milestones tailored for K_max ~ 20000. 250 below 1000, 500 below 5000,
    1000 above 5000."""
    base = []
    for k in range(250, 1001, 250):
        base.append(k)
    for k in range(1500, 5001, 500):
        base.append(k)
    for k in range(6000, K_max + 1, 1000):
        base.append(k)
    if base[-1] != K_max:
        base.append(K_max)
    return sorted(set(base))


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--N", type=int, default=30)
    p.add_argument(
        "--anchors", default="5.0,5.3,5.5",
        help="comma list of log10 x anchors",
    )
    p.add_argument("--window", type=float, default=0.5)
    p.add_argument("--K_max", type=int, default=20000)
    p.add_argument("--eps_list", default="1,3")
    p.add_argument("--M", type=int, default=4)
    p.add_argument("--dps", type=int, default=18)
    p.add_argument("--csv_out", default=os.path.join(HERE, "slot2_extended_traces.csv"))
    p.add_argument("--summary_csv", default=os.path.join(HERE, "slot2_extended_summary.csv"))
    p.add_argument("--seed", type=int, default=0)
    return p.parse_args()


def main():
    args = parse_args()
    random.seed(args.seed)

    anchors = [float(s) for s in args.anchors.split(",")]
    eps_list = [float(s) for s in args.eps_list.split(",")]
    K_milestones = make_extended_K_milestones(args.K_max)
    print(f"K milestones (count={len(K_milestones)}): "
          f"{K_milestones[:5]}...{K_milestones[-5:]}", flush=True)

    max_x = int((10 ** max(anchors)) * (1 + args.window))
    sieve_n = max_x + 100
    print(f"Sieving pi(n) for n up to {sieve_n}...", flush=True)
    t0 = time.time()
    pi_table = sieve_pi(sieve_n)
    print(f"  sieve done in {time.time() - t0:.1f}s, pi({sieve_n}) = {pi_table[sieve_n]}",
          flush=True)

    print(f"Loading {args.K_max} zeta zeros (dps={args.dps})...", flush=True)
    t0 = time.time()
    gammas = load_extended_zeros(args.K_max, dps=args.dps)
    print(f"  loaded {len(gammas)} zeros in {time.time() - t0:.1f}s "
          f"(gamma_K = {float(gammas[-1]):.4f}).", flush=True)

    all_traces: Dict[float, List[Tuple[int, Dict[int, float]]]] = {}
    csv_rows: List[List] = []
    for log10_x in anchors:
        x_lo = int(10 ** log10_x)
        x_hi = int((10 ** log10_x) * (1 + args.window))
        x_samples = sorted(random.sample(range(x_lo, x_hi + 1), args.N))
        print(f"\n[anchor 10^{log10_x}, sampling {args.N} from "
              f"[{x_lo}, {x_hi}], K_max={args.K_max}]", flush=True)
        traces: List[Tuple[int, Dict[int, float]]] = []
        anchor_t0 = time.time()
        for i, x in enumerate(x_samples):
            pi_x = pi_table[x]
            errs = cumulative_err_trace(
                x, pi_x, gammas, args.K_max, args.M, args.dps, K_milestones
            )
            traces.append((x, errs))
            for K, e in errs.items():
                csv_rows.append([log10_x, x, pi_x, K, f"{e:.6f}"])
            if (i + 1) % 5 == 0 or i == 0:
                elapsed = time.time() - anchor_t0
                avg = elapsed / (i + 1)
                rem = avg * (args.N - i - 1)
                print(f"  [{i+1}/{args.N}] x={x} pi={pi_x} "
                      f"err@K=5000:{errs.get(5000, float('nan')):.3f} "
                      f"err@K_max:{errs[K_milestones[-1]]:.3f}  "
                      f"({elapsed:.0f}s, ETA {rem:.0f}s)", flush=True)
        all_traces[log10_x] = traces

    print("\n" + "=" * 100)
    print("K_emp = min K such that WORST-CASE |pi_K - pi(x)| <= eps over N samples")
    print(f"{'log10x':>8s} {'x':>8s} {'eps':>4s} "
          f"{'K_emp':>7s} {'K_typ':>6s} "
          f"{'K_pred_hG':>10s} {'K_gal(c=1)':>11s} "
          f"{'c_emp':>8s} {'c_emp_T7':>9s} {'sigma_K_emp':>11s}")
    print("-" * 100)

    summary_rows: List[List] = []
    for log10_x in anchors:
        x = 10 ** log10_x
        traces = [t[1] for t in all_traces[log10_x]]
        for eps in eps_list:
            K_emp, K_typ = find_K_emp(traces, eps, K_milestones)
            sigma_eq = eps / math.sqrt(2 * math.log(args.N))
            K_pred = K_typical(x, sigma_eq)
            K_galway = K_galway_grh_bound(x, c=1.0)
            sigma_at_K_emp = sigma_predicted(x, K_emp) if K_emp <= args.K_max else float("nan")
            c_emp = (K_emp / (math.sqrt(x) * math.log(x) ** 2) if K_emp <= args.K_max
                     else float("nan"))
            c_emp_T7_pred = 0.61 * K_pred / (math.sqrt(x) * math.log(x) ** 2)
            print(f"{log10_x:>8.2f} {x:>8.0f} {eps:>4.0f} "
                  f"{K_emp:>7d} {K_typ:>6d} "
                  f"{K_pred:>10d} {K_galway:>11d} "
                  f"{c_emp:>8.4f} {c_emp_T7_pred:>9.4f} {sigma_at_K_emp:>11.3f}")
            summary_rows.append([
                f"{log10_x:.2f}", x, eps, K_emp, K_typ, K_pred, K_galway,
                f"{c_emp:.6g}", f"{c_emp_T7_pred:.6g}", f"{sigma_at_K_emp:.6g}"
            ])

    print("-" * 100)

    with open(args.csv_out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["log10_x_anchor", "x", "pi_x", "K", "abs_err"])
        for r in csv_rows:
            w.writerow(r)
    print(f"\nWrote {len(csv_rows)} trace rows to {args.csv_out}")

    with open(args.summary_csv, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow([
            "log10_x_anchor", "x", "eps", "K_emp", "K_typ",
            "K_pred_half_gauss", "K_galway_c1", "c_emp", "c_emp_T7_pred",
            "sigma_at_K_emp"
        ])
        for r in summary_rows:
            w.writerow(r)
    print(f"Wrote {len(summary_rows)} summary rows to {args.summary_csv}")


if __name__ == "__main__":
    main()

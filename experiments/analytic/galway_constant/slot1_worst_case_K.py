"""
P5 / Thread 10 (commit, slot 1) — Galway worst-case K-constant tightening.

Target: OPEN_POSITIVE_TARGETS.md §P5 — what is the smallest K(x) such that the
truncated explicit-formula evaluator pi_K(x) attains pi(x) +/- eps WORST-CASE
over many random x in a window?

Galway 2004 (under GRH) gives K = O(x^{1/2 + epsilon}) worst case, with various
explicit constants in the prefactor for unconditional / GRH variants. The
Hardy-Littlewood / Riemann smooth approximation R(x) by itself gives error
~sqrt(x). Thread 7 (S240..S244) characterised the *typical* error as
sigma(K, x) ~ sqrt(x) log K / (pi sqrt(2K) log x) under Montgomery's random-
phase heuristic, and the empirical sigma_eff = 0.755 sigma_pred (GUE pair-
correlation factor).

This script measures the WORST-CASE K-frontier:

  K_emp(x, eps, N) := min K s.t. max_{i=1..N} |pi_K(x_i) - pi(x_i)| <= eps,

where x_i is sampled in a tight window around an anchor x. The outputs are:
  (a) the K_emp scaling vs x (does it look like x^{1/2}, x^{1/2}/log x, or x?)
  (b) the explicit prefactor c_emp(x) in K_emp ~ c_emp(x) * sqrt(x) * log^2 x.
  (c) an empirical comparison to the half-Gaussian-tail prediction
      K_emp(x, eps, N) ~ K_typical(x, eps / sqrt(2 log N)).

For Galway's conditional GRH bound K ~ c_Galway * x^{1/2 + epsilon}, we expect
c_emp <= c_Galway by a tight factor (c_Galway is loose in published proofs;
empirical is meant to test how loose).

Inputs (CLI):
  --N            samples per anchor (default 30)
  --anchors      comma list of anchor exponents log10 x (default 4.0, 4.5, 5.0)
  --window       multiplicative window radius (default 0.5 -> [x, 1.5x])
  --K_max        cap on K (default 8000)
  --eps_list     comma list of epsilon thresholds (default 1, 3, 10, 30, 100)
  --M            Mobius truncation in R_at_rho (default 4)
  --dps          mpmath dps (default 18)
  --csv_out      output csv path
  --seed         random seed (default 0)

Outputs:
  --csv_out: per-(anchor, x_i, K) trace CSV with |pi_K(x_i) - pi(x_i)|.
  Stdout: K_emp(anchor, eps) summary table.
"""

from __future__ import annotations

import argparse
import csv
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

from polylog_approx_pi import R_at_rho, get_zeros, riemann_R  # noqa: E402


# -----------------------------------------------------------------------------
# Sieve pi(n) for all n <= N.
# -----------------------------------------------------------------------------

def sieve_pi(N: int) -> List[int]:
    """Return cumulative pi(n) for n = 0..N via Eratosthenes."""
    is_p = bytearray([1]) * (N + 1)
    is_p[0] = 0
    if N >= 1:
        is_p[1] = 0
    for p in range(2, int(N ** 0.5) + 1):
        if is_p[p]:
            for j in range(p * p, N + 1, p):
                is_p[j] = 0
    pi = [0] * (N + 1)
    c = 0
    for n in range(N + 1):
        if is_p[n]:
            c += 1
        pi[n] = c
    return pi


# -----------------------------------------------------------------------------
# Cumulative pi_K(x) trace at K = 1..K_max for a single x.
# -----------------------------------------------------------------------------

def cumulative_err_trace(
    x: int,
    pi_x: int,
    gammas: List[mp.mpf],
    K_max: int,
    M: int,
    dps: int,
    K_milestones: List[int],
) -> Dict[int, float]:
    """Run pi_K(x) accumulating R_at_rho corrections; record |pi_K - pi(x)| at K_milestones.

    Returns {K -> |error|} for K in K_milestones.
    """
    R_x = riemann_R(x, M=20, dps=dps)
    correction = mp.mpf(0)
    out: Dict[int, float] = {}
    milestones_sorted = sorted(set(K_milestones))
    K_done = 0
    for K_target in milestones_sorted:
        K_target_capped = min(K_target, K_max)
        for k in range(K_done, K_target_capped):
            Rrho = R_at_rho(x, gammas[k], M=M, dps=dps)
            correction += 2 * Rrho.real
        K_done = K_target_capped
        pi_K = R_x - correction
        out[K_target] = abs(float(pi_K) - pi_x)
    return out


# -----------------------------------------------------------------------------
# K-milestones: log-uniform from K_min to K_max with extra fineness near low K.
# -----------------------------------------------------------------------------

def make_K_milestones(K_max: int) -> List[int]:
    base = [1, 2, 5, 10, 20, 50, 100, 200, 300, 500, 700,
            1000, 1500, 2000, 3000, 4000, 5000, 6000, 7000, 8000]
    return [k for k in base if k <= K_max]


# -----------------------------------------------------------------------------
# Theoretical predictions.
# -----------------------------------------------------------------------------

def sigma_predicted(x: float, K: int) -> float:
    """sigma(K, x) ~ sqrt(x) log K / (pi sqrt(2K) log x)  (Thread 7 random-phase)."""
    if K <= 1 or x <= 1:
        return float("inf")
    return float(
        mp.sqrt(mp.mpf(x)) * mp.log(K) / (mp.pi * mp.sqrt(2 * K) * mp.log(x))
    )


def K_typical(x: float, sigma_target: float) -> int:
    """Invert sigma(K, x) = sigma_target. Returns smallest K satisfying it."""
    if sigma_target <= 0:
        return 10 ** 9
    lo, hi = 2, 10 ** 9
    while lo < hi:
        mid = (lo + hi) // 2
        if sigma_predicted(x, mid) <= sigma_target:
            hi = mid
        else:
            lo = mid + 1
    return lo


def K_galway_grh_bound(x: float, c: float = 1.0) -> int:
    """Galway 2004-style explicit GRH worst-case bound K = c * x^{1/2} log^2 x.

    The published constant in Galway (and Buthe 2018 corollary B) is somewhat
    loose; c around 5..50 is typical. We report c=1 as the geometric anchor.
    """
    return max(2, int(round(c * (x ** 0.5) * (math.log(x) ** 2))))


# -----------------------------------------------------------------------------
# K_emp finder: smallest K in milestones with err <= eps WORST-CASE over samples.
# -----------------------------------------------------------------------------

def find_K_emp(
    err_traces: List[Dict[int, float]],
    eps: float,
    milestones: List[int],
) -> Tuple[int, int]:
    """Return (K_emp, K_emp_typical). K_emp = smallest milestone K with
    max_i err_traces[i][K] <= eps. K_emp_typical = smallest with median(err) <= eps.
    Returns (K_max + 1, K_max + 1) if not reached.

    `K_max + 1` value is a sentinel indicating exceeded budget.
    """
    K_emp = milestones[-1] + 1
    K_typ = milestones[-1] + 1
    for K in sorted(milestones):
        errs = sorted(t[K] for t in err_traces)
        worst = errs[-1]
        median = errs[len(errs) // 2]
        if K_emp > milestones[-1] and worst <= eps:
            K_emp = K
        if K_typ > milestones[-1] and median <= eps:
            K_typ = K
    return K_emp, K_typ


# -----------------------------------------------------------------------------
# Main.
# -----------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--N", type=int, default=30)
    p.add_argument("--anchors", default="4.0,4.5,5.0",
                   help="comma list of log10 x anchors")
    p.add_argument("--window", type=float, default=0.5)
    p.add_argument("--K_max", type=int, default=8000)
    p.add_argument("--eps_list", default="1,3,10,30,100,300")
    p.add_argument("--M", type=int, default=4)
    p.add_argument("--dps", type=int, default=18)
    p.add_argument("--csv_out", default=os.path.join(HERE, "slot1_traces.csv"))
    p.add_argument("--summary_csv", default=os.path.join(HERE, "slot1_summary.csv"))
    p.add_argument("--seed", type=int, default=0)
    return p.parse_args()


def main():
    args = parse_args()
    random.seed(args.seed)

    anchors = [float(s) for s in args.anchors.split(",")]
    eps_list = [float(s) for s in args.eps_list.split(",")]
    K_milestones = make_K_milestones(args.K_max)

    # Determine sieve range from largest anchor
    max_x = int((10 ** max(anchors)) * (1 + args.window))
    sieve_n = max_x + 100
    print(f"Sieving pi(n) for n up to {sieve_n}...", flush=True)
    t0 = time.time()
    pi_table = sieve_pi(sieve_n)
    print(f"  sieve done in {time.time() - t0:.1f}s, pi({sieve_n}) = {pi_table[sieve_n]}",
          flush=True)

    print(f"Loading {args.K_max} zeta zeros (dps={args.dps})...", flush=True)
    t0 = time.time()
    gammas = get_zeros(args.K_max, dps=args.dps)
    print(f"  loaded {len(gammas)} zeros in {time.time() - t0:.1f}s "
          f"(gamma_K = {float(gammas[-1]):.4f}).", flush=True)

    # ---- Sample x_i for each anchor and compute err traces ----
    all_traces: Dict[float, List[Tuple[int, Dict[int, float]]]] = {}
    csv_rows: List[List] = []
    for log10_x in anchors:
        x_lo = int(10 ** log10_x)
        x_hi = int((10 ** log10_x) * (1 + args.window))
        x_samples = sorted(random.sample(range(x_lo, x_hi + 1), args.N))
        print(f"\n[anchor 10^{log10_x}, sampling {args.N} from "
              f"[{x_lo}, {x_hi}], K_max = {args.K_max}]", flush=True)
        traces: List[Tuple[int, Dict[int, float]]] = []
        anchor_t0 = time.time()
        for i, x in enumerate(x_samples):
            pi_x = pi_table[x]
            t0 = time.time()
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
                      f"err@K=1000:{errs.get(1000, float('nan')):.2f} "
                      f"err@K_max:{errs[K_milestones[-1]]:.2f}  "
                      f"({elapsed:.0f}s, ETA {rem:.0f}s)", flush=True)
        all_traces[log10_x] = traces

    # ---- Find K_emp(anchor, eps) and K_typ(anchor, eps) ----
    print("\n" + "=" * 96)
    print("K_emp = min K such that WORST-CASE |pi_K(x_i) - pi(x_i)| <= eps over N samples")
    print(f"{'log10x':>8s} {'x_anchor':>10s} {'eps':>6s} "
          f"{'K_emp':>8s} {'K_typ':>8s} "
          f"{'K_pred(half-G)':>15s} {'K_galway(c=1)':>14s} "
          f"{'c_emp':>8s} {'sigma_K_emp':>11s}")
    print("-" * 96)

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
            print(f"{log10_x:>8.2f} {x:>10.0f} {eps:>6.0f} "
                  f"{K_emp:>8d} {K_typ:>8d} "
                  f"{K_pred:>15d} {K_galway:>14d} "
                  f"{c_emp:>8.4f} {sigma_at_K_emp:>11.3f}")
            summary_rows.append([
                f"{log10_x:.2f}", x, eps, K_emp, K_typ, K_pred, K_galway,
                f"{c_emp:.6g}", f"{sigma_at_K_emp:.6g}"
            ])

    print("-" * 96)
    print("legend:")
    print("  K_emp = smallest milestone K with worst-case err <= eps")
    print("  K_typ = smallest milestone K with MEDIAN err <= eps")
    print("  K_pred(half-G) = K such that sigma(K) = eps / sqrt(2 log N)  (theoretical worst")
    print("                   under Gaussian approx; sigma uses S195 random-phase formula)")
    print("  K_galway(c=1) = c * x^{1/2} log^2 x with c=1  (Galway-style explicit bound)")
    print("  c_emp = K_emp / (sqrt(x) * log^2 x)  -- empirical Galway-bound prefactor")

    # ---- CSV outputs ----
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
            "K_pred_half_gauss", "K_galway_c1", "c_emp", "sigma_at_K_emp"
        ])
        for r in summary_rows:
            w.writerow(r)
    print(f"Wrote {len(summary_rows)} summary rows to {args.summary_csv}")


if __name__ == "__main__":
    main()

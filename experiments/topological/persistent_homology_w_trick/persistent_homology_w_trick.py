#!/usr/bin/env python3
"""
D2.a — Persistent homology of W=210 W-tricked normalised prime gaps
====================================================================

Successor to D2 (S96, experiments/topological/persistent_homology_chi_p/).

Question (NOVELTY_CHALLENGES.md §D2.a): the W-trick at W=210 restores
Gowers U^2 uniformity of chi_P to within 0.1% (E2.13). Does it also
restore PH-Poisson behaviour of the gap sequence?

Protocol
--------

Same as D2 but the input sequence is filtered to a single residue class
b mod W = 210 (gcd(b, W) = 1):

  q_1 < q_2 < q_3 < ... = primes congruent to b mod 210

Cramer-normalisation on a single residue class:
  expected gap density = 1/(phi(W) log x), so
  E[q_{n+1} - q_n] ~ phi(W) log q_n.
  Take x_n := (q_{n+1} - q_n) / (phi(W) log q_n) — under Cramer +
  Dirichlet equidistribution this is asymptotically Exp(1).

Two falsifiers (pre-stated):

F_a (W-trick erases PH deviation):
    For all summary stats S in {T0, T1, T0_ratio_to_iid, T1_ratio_to_iid},
    z(W-tricked, B2) lies in [-3, +3] AND
    z(W-tricked, B1) lies in [-3, +3].
    -> E2.17 (S96 PH deviation) ~= E2.13 (Gowers HL singular series) in
       topological clothes. C-grade duplicate-plus refinement of E2.17.

F_b (W-trick partially erases the deviation):
    Either |z(W-tricked, B1)| or |z(W-tricked, B2)| reduced by >= 50%
    relative to S96 baseline (z = -10.31 / -8.70 / -4.20 / -11.99) but
    at least one |z| stays above 3.
    -> Partial reduction. PH detects a pure-singular-series component
       AND a residual structural component.  B-grade refinement of E2.17.

F_c (W-trick does NOT erase the deviation):
    All four S96 z-scores reproduced within 25% relative on the
    W-tricked sequence (e.g., T0 z(B2) stays in [-12, -6.5]; T1 z(B2)
    in [-15, -9]).
    -> PH detects an obstruction beyond singular-series cancellation.
    Distinct from E2.13. Genuine new content. B-grade negative-shape
    edge candidate.

Comparison anchors are the S96 main run at M=2000, d=3, x ~ 10^6:
    T0_PRIMES = 243.34, B2 mean = 277.43, std = 3.92 -> z = -8.70
    T1_PRIMES =  37.24, B2 mean =  56.09, std = 1.57 -> z = -11.99
    T0 z(B1) = -10.31, T1 z(B1) = -4.20

We pick THREE residue classes b in {1, 11, 13} (each coprime to 210)
and average the W-tricked summary across the three residues. This
controls for the (negligible) Dirichlet equidistribution bias and
gives a single comparison point.

Cross-domain refs:
- Carlsson 2009 BAMS 46.
- Edelsbrunner-Harer 2010.
- Bauer 2021 (ripser).
- Green-Tao 2008 (W-trick origin).
"""

import argparse
import json
import math
import time
from pathlib import Path

import numpy as np
from ripser import ripser


# ---------------------------------------------------------------------
# Prime sieve
# ---------------------------------------------------------------------

def segmented_sieve(N: int) -> np.ndarray:
    sieve = np.ones(N + 1, dtype=bool)
    sieve[:2] = False
    for i in range(2, int(N**0.5) + 1):
        if sieve[i]:
            sieve[i*i::i] = False
    return np.flatnonzero(sieve).astype(np.int64)


# ---------------------------------------------------------------------
# Persistence diagram summaries
# ---------------------------------------------------------------------

def diagram_summaries(dgm0, dgm1, tau_thr: float) -> dict:
    out = {}
    if len(dgm0):
        d0 = dgm0[np.isfinite(dgm0[:, 1])]
        if len(d0):
            life0 = d0[:, 1] - d0[:, 0]
            out["T0"] = float(life0.sum())
            out["L0"] = float(life0.max())
        else:
            out["T0"] = 0.0
            out["L0"] = 0.0
    else:
        out["T0"] = 0.0
        out["L0"] = 0.0
    if len(dgm1):
        life1 = dgm1[:, 1] - dgm1[:, 0]
        life1 = life1[np.isfinite(life1)]
        if len(life1):
            out["T1"] = float(life1.sum())
            out["L1"] = float(life1.max())
            out["N1"] = int((life1 > tau_thr).sum())
        else:
            out["T1"] = 0.0
            out["L1"] = 0.0
            out["N1"] = 0
    else:
        out["T1"] = 0.0
        out["L1"] = 0.0
        out["N1"] = 0
    return out


def takens_embed(x: np.ndarray, d: int, tau: int = 1) -> np.ndarray:
    n = len(x) - (d - 1) * tau
    if n <= 0:
        raise ValueError("not enough samples for embedding")
    y = np.zeros((n, d), dtype=float)
    for j in range(d):
        y[:, j] = x[j*tau : j*tau + n]
    return y


def ph01(cloud: np.ndarray, max_edge: float) -> dict:
    res = ripser(cloud, maxdim=1, thresh=max_edge)
    dgms = res["dgms"]
    return diagram_summaries(dgms[0], dgms[1], tau_thr=0.5)


# ---------------------------------------------------------------------
# Experiment
# ---------------------------------------------------------------------

def w_tricked_normalised_gaps(primes: np.ndarray, W: int, b: int,
                              start_x: int, M: int) -> tuple:
    """Return (M consecutive normalised gaps in residue b mod W,
    primes used).  The starting prime is the first p >= start_x with
    p % W == b.

    Normalisation: x_n = (q_{n+1} - q_n) / (phi(W) * log q_n).
    """
    phi_W = phi(W)
    # primes congruent to b mod W
    mask = (primes % W) == b
    sub = primes[mask]
    if len(sub) < M + 2:
        raise ValueError(f"only {len(sub)} primes ~ {b} mod {W}; need >= {M+2}")
    # find start
    start_idx = int(np.searchsorted(sub, start_x))
    if start_idx + M + 1 >= len(sub):
        raise ValueError(f"not enough primes mod {W} after start_x")
    q = sub[start_idx : start_idx + M + 1]   # M+1 primes -> M gaps
    g = np.diff(q).astype(float)
    log_q = np.log(q[:-1].astype(float))
    x = g / (phi_W * log_q)
    return x, q


def phi(n: int) -> int:
    """Euler totient."""
    result = n
    i = 2
    nn = n
    while i * i <= nn:
        if nn % i == 0:
            while nn % i == 0:
                nn //= i
            result -= result // i
        i += 1
    if nn > 1:
        result -= result // nn
    return result


def run(args):
    rng = np.random.default_rng(args.seed)
    t0 = time.time()

    # --- gather primes ---
    primes = segmented_sieve(args.N_max)
    n_primes = len(primes)
    print(f"[setup] sieve to {args.N_max:,}: {n_primes:,} primes "
          f"in {time.time()-t0:.1f}s")

    W = args.W
    phi_W = phi(W)
    print(f"[W-trick] W={W}, phi(W)={phi_W}, "
          f"residues used: {args.residues}")

    M = args.M
    d = args.d
    max_edge = args.max_edge
    K = args.K

    per_residue = []
    summaries_per_residue = []
    for b in args.residues:
        if math.gcd(b, W) != 1:
            print(f"  [skip] gcd({b},{W}) != 1; residue invalid")
            continue
        x, q = w_tricked_normalised_gaps(primes, W, b, args.start_x, M)
        print(f"[res b={b}] M={M}, q_start={int(q[0])}, q_end={int(q[-1])}, "
              f"mean(norm gap)={x.mean():.4f}, std={x.std():.4f}")

        cloud_p = takens_embed(x, d, tau=1)
        t1 = time.time()
        s_primes = ph01(cloud_p, max_edge=max_edge)
        print(f"  [PH] {time.time()-t1:.1f}s; summaries={s_primes}")

        # baselines (per residue, fresh draws to keep K-replicate stats valid)
        sums_B1 = []
        for k in range(K):
            xb = rng.exponential(scale=1.0, size=M).astype(float)
            cloud = takens_embed(xb, d, tau=1)
            sums_B1.append(ph01(cloud, max_edge=max_edge))
        sums_B2 = []
        for k in range(K):
            xb = rng.permutation(x)
            cloud = takens_embed(xb, d, tau=1)
            sums_B2.append(ph01(cloud, max_edge=max_edge))

        keys = ["T0", "L0", "T1", "L1", "N1"]
        def agg(rows):
            r = {}
            for k in keys:
                arr = np.array([rr[k] for rr in rows], dtype=float)
                r[k] = {"mean": float(arr.mean()),
                        "std":  float(arr.std(ddof=1)) if len(arr) > 1 else 0.0,
                        "min":  float(arr.min()),
                        "max":  float(arr.max())}
            return r
        B1 = agg(sums_B1)
        B2 = agg(sums_B2)
        z1 = {k: (s_primes[k] - B1[k]["mean"]) / max(B1[k]["std"], 1e-12)
              for k in keys}
        z2 = {k: (s_primes[k] - B2[k]["mean"]) / max(B2[k]["std"], 1e-12)
              for k in keys}
        rank1 = {k: float(sum(1 for r in sums_B1 if r[k] <= s_primes[k]))
                    / float(len(sums_B1)) for k in keys}
        rank2 = {k: float(sum(1 for r in sums_B2 if r[k] <= s_primes[k]))
                    / float(len(sums_B2)) for k in keys}

        per_residue.append({
            "b": b,
            "q_start": int(q[0]),
            "q_end": int(q[-1]),
            "x_mean": float(x.mean()),
            "x_std":  float(x.std()),
            "primes": s_primes,
            "B1": B1,
            "B2": B2,
            "z_B1": z1,
            "z_B2": z2,
            "rank_B1": rank1,
            "rank_B2": rank2,
        })
        summaries_per_residue.append({
            "b": b, "z_B1_T0": z1["T0"], "z_B2_T0": z2["T0"],
            "z_B1_T1": z1["T1"], "z_B2_T1": z2["T1"],
            "T0_PRIMES": s_primes["T0"], "T1_PRIMES": s_primes["T1"],
        })

    # ----- aggregate across residues -----
    keys = ["T0", "L0", "T1", "L1", "N1"]
    z_keys = ["z_B1", "z_B2"]
    pooled = {}
    for zk in z_keys:
        pooled[zk] = {}
        for k in keys:
            zs = [pr[zk][k] for pr in per_residue]
            pooled[zk][k] = {
                "mean": float(np.mean(zs)),
                "std":  float(np.std(zs, ddof=1)) if len(zs) > 1 else 0.0,
                "min":  float(np.min(zs)),
                "max":  float(np.max(zs)),
            }

    out = {
        "params": {
            "N_max": args.N_max,
            "W": W,
            "phi_W": phi_W,
            "residues": list(args.residues),
            "start_x": args.start_x,
            "M": M, "d": d, "tau": 1,
            "max_edge": max_edge,
            "K_baselines": K,
            "tau_thr_for_N1": 0.5,
            "seed": args.seed,
        },
        "per_residue": per_residue,
        "summaries_per_residue": summaries_per_residue,
        "pooled_z_across_residues": pooled,
        # S96 baseline (anchored, NOT recomputed here)
        "s96_baseline_main_run": {
            "M": 2000, "d": 3, "x": 1_000_000, "K": 20,
            "T0_PRIMES": 243.34, "T1_PRIMES": 37.24,
            "z_B1": {"T0": -10.31, "T1": -4.20},
            "z_B2": {"T0": -8.70,  "T1": -11.99},
            "ref": "experiments/topological/persistent_homology_chi_p/main_run_d3.json",
        },
        "elapsed_sec": time.time() - t0,
    }
    return out


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N-max", type=int, default=10_000_000,
                    help="upper bound for sieve")
    ap.add_argument("--W", type=int, default=210,
                    help="W-trick modulus (210 = 2*3*5*7)")
    ap.add_argument("--residues", type=int, nargs="+",
                    default=[1, 11, 13],
                    help="residue classes mod W (gcd(b,W)=1)")
    ap.add_argument("--start-x", type=int, default=1_000_000,
                    help="window start: nearest prime in residue b")
    ap.add_argument("--M", type=int, default=2000,
                    help="window length")
    ap.add_argument("--d", type=int, default=3,
                    help="Takens embedding dimension")
    ap.add_argument("--max-edge", type=float, default=4.0,
                    help="ripser max edge")
    ap.add_argument("--K", type=int, default=20,
                    help="baseline replicates per residue")
    ap.add_argument("--seed", type=int, default=20260427)
    ap.add_argument("--out",
                    default="w_trick_main.json")
    args = ap.parse_args()

    out = run(args)
    Path(args.out).write_text(json.dumps(out, indent=2))
    print("\n[done] wrote", args.out)
    print("\nPooled z across residues (W-tricked vs B1, IID Exp(1)):")
    for k in ["T0", "T1"]:
        m = out["pooled_z_across_residues"]["z_B1"][k]
        print(f"  {k}: mean = {m['mean']:+.2f}, range = "
              f"[{m['min']:+.2f}, {m['max']:+.2f}]")
    print("\nPooled z across residues (W-tricked vs B2, gap-permuted):")
    for k in ["T0", "T1"]:
        m = out["pooled_z_across_residues"]["z_B2"][k]
        print(f"  {k}: mean = {m['mean']:+.2f}, range = "
              f"[{m['min']:+.2f}, {m['max']:+.2f}]")
    print("\nS96 unconditioned baseline:")
    print(f"  T0 z(B1) = {out['s96_baseline_main_run']['z_B1']['T0']}, "
          f"T1 z(B1) = {out['s96_baseline_main_run']['z_B1']['T1']}")
    print(f"  T0 z(B2) = {out['s96_baseline_main_run']['z_B2']['T0']}, "
          f"T1 z(B2) = {out['s96_baseline_main_run']['z_B2']['T1']}")


if __name__ == "__main__":
    main()

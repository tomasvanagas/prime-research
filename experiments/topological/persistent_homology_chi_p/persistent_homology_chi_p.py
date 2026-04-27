#!/usr/bin/env python3
"""
D2 — Persistent homology of the prime gap sequence vs Poisson baseline
======================================================================

Cross-domain technique: persistent homology (Carlsson 2009 BAMS 46) on
Vietoris-Rips filtrations of Takens time-delay embeddings.

Question (ATTACK_VECTORS.md §D2): does the prime gap sequence carry
persistent topological features (in H_0 / H_1) that statistically
distinguish it from a matched-density Poisson process at the same
sample size?

Protocol
--------

1. Compute primes up to N_max via segmented sieve.
2. Build the gap sequence g_n = p_{n+1} - p_n.
3. Cramér-normalise: x_n = g_n / log(p_n) (asymptotically Exp(1) under
   the Cramér model).
4. Take a moving window of M consecutive normalised gaps starting near
   `start` and Takens-embed at delay tau=1 in dimension d:
       y_n = (x_n, x_{n+1}, ..., x_{n+d-1}) in R^d.
5. Compute persistent H_0, H_1 of the Vietoris-Rips filtration of the
   point cloud {y_n} via ripser.
6. Summary statistics on each persistence diagram:
       T0  = total persistence H_0 (sum of finite lifetimes)
       T1  = total persistence H_1
       L1  = max persistence H_1
       N1  = #{H_1 features with lifetime > tau_thr}  (tau_thr = 0.5)
       L0  = max finite lifetime H_0  (= largest connected-component
                                         merge distance)
7. Repeat the same computation on K Poisson baselines:
       (B1) IID Exp(1) with same M, d
       (B2) IID permutation of the empirical normalised-gap sample
            (preserves marginal, kills serial correlation)
8. Report (PRIMES, B1 mean +/- std, B2 mean +/- std) and bootstrap
   p-value: how many B1 / B2 replicates exceed the PRIMES summary.

Falsification (pre-registered)
------------------------------
F1: PRIMES summaries are within 2 sigma of BOTH B1 and B2 means for
    ALL of {T0, T1, L1, N1, L0} → primes are PH-indistinguishable from
    Poisson at this scale → 38th pseudorandomness measure (B-grade
    closure of D2).

F2: PRIMES is > 3 sigma above OR below B1 on at least one summary, but
    B2 (gap-permuted) lands in the same bin as PRIMES → the deviation
    is captured by the gap MARGINAL distribution, not serial structure
    → reduces D2 to the well-known gap-distribution non-Poisson tail
    (Goldston-Yıldırım, etc.) → C-grade duplicate refinement.

F3: PRIMES > 3 sigma above OR below BOTH B1 AND B2 → genuine serial-
    correlation persistent topological signature → B-grade negative-
    shape edge (PH detects a structure that gap-pair statistics miss).

Cross-domain ref:
- Carlsson 2009 "Topology and data" Bull. AMS 46(2), 255-308.
- Edelsbrunner-Harer 2010 "Computational Topology: An Introduction".
- ripser (Bauer 2021 J. Appl. Comput. Topol.).
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
    """Return all primes <= N as a sorted numpy array of int64."""
    sieve = np.ones(N + 1, dtype=bool)
    sieve[:2] = False
    for i in range(2, int(N**0.5) + 1):
        if sieve[i]:
            sieve[i*i::i] = False
    return np.flatnonzero(sieve).astype(np.int64)


# ---------------------------------------------------------------------
# Persistence diagram summaries
# ---------------------------------------------------------------------

def diagram_summaries(dgm0: np.ndarray, dgm1: np.ndarray,
                      tau_thr: float) -> dict:
    """
    dgm0, dgm1: ripser persistence diagrams (Nx2 arrays, possibly with
    inf in death column).  All finite features in H_0 measure when each
    cluster merged into the next; the unique infinite H_0 bar is the
    survivor of the global connected-component count.
    """
    out = {}
    # H_0
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
    # H_1
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


# ---------------------------------------------------------------------
# Takens embedding + ripser
# ---------------------------------------------------------------------

def takens_embed(x: np.ndarray, d: int, tau: int = 1) -> np.ndarray:
    """y_n = (x_n, x_{n+tau}, ..., x_{n+(d-1)*tau}) for n=0..N-(d-1)*tau-1."""
    n = len(x) - (d - 1) * tau
    if n <= 0:
        raise ValueError("not enough samples for embedding")
    y = np.zeros((n, d), dtype=float)
    for j in range(d):
        y[:, j] = x[j*tau : j*tau + n]
    return y


def ph01(cloud: np.ndarray, max_edge: float) -> dict:
    """Persistent H_0, H_1 summaries of cloud."""
    res = ripser(cloud, maxdim=1, thresh=max_edge)
    dgms = res["dgms"]
    return diagram_summaries(dgms[0], dgms[1], tau_thr=0.5)


# ---------------------------------------------------------------------
# Experiment
# ---------------------------------------------------------------------

def run(args):
    rng = np.random.default_rng(args.seed)
    t0 = time.time()

    # --- gather primes and gap window ---
    primes = segmented_sieve(args.N_max)
    n_primes = len(primes)
    print(f"[setup] sieve to {args.N_max:,}: {n_primes:,} primes "
          f"in {time.time()-t0:.1f}s")

    # gaps
    gaps = np.diff(primes).astype(float)            # length n_primes - 1
    log_p = np.log(primes[:-1].astype(float))
    norm = gaps / log_p                              # Cramér-normalised

    # window: M consecutive normalised gaps starting at the prime nearest
    # to args.start_x.
    start_idx = int(np.searchsorted(primes, args.start_x))
    M = args.M
    if start_idx + M >= len(norm):
        raise ValueError(f"window {M} starting at idx {start_idx} "
                         f"exceeds available gaps {len(norm)}")
    primes_window = norm[start_idx : start_idx + M]
    p_at_start = int(primes[start_idx])
    print(f"[window] M={M}, start prime ~ {p_at_start:,}, "
          f"mean(norm gap) = {primes_window.mean():.4f}, "
          f"std = {primes_window.std():.4f}")

    # --- PRIMES embedding + PH ---
    cloud_p = takens_embed(primes_window, args.d, tau=1)
    print(f"[primes] cloud {cloud_p.shape}; running ripser...")
    t1 = time.time()
    s_primes = ph01(cloud_p, max_edge=args.max_edge)
    print(f"[primes] PH in {time.time()-t1:.1f}s; summaries = {s_primes}")

    # --- B1: IID Exp(1) baselines ---
    print(f"[B1] running {args.K} IID Exp(1) baselines, M={M}, d={args.d}")
    sums_B1 = []
    for k in range(args.K):
        x = rng.exponential(scale=1.0, size=M).astype(float)
        cloud = takens_embed(x, args.d, tau=1)
        sums_B1.append(ph01(cloud, max_edge=args.max_edge))
        if (k + 1) % 5 == 0:
            print(f"   B1 {k+1}/{args.K}")

    # --- B2: gap-permuted baselines ---
    print(f"[B2] running {args.K} gap-permuted baselines")
    sums_B2 = []
    for k in range(args.K):
        x = rng.permutation(primes_window)
        cloud = takens_embed(x, args.d, tau=1)
        sums_B2.append(ph01(cloud, max_edge=args.max_edge))
        if (k + 1) % 5 == 0:
            print(f"   B2 {k+1}/{args.K}")

    # --- aggregate ---
    keys = ["T0", "L0", "T1", "L1", "N1"]
    def agg(rows):
        d = {}
        for k in keys:
            arr = np.array([r[k] for r in rows], dtype=float)
            d[k] = {"mean": float(arr.mean()),
                    "std":  float(arr.std(ddof=1)) if len(arr) > 1 else 0.0,
                    "min":  float(arr.min()),
                    "max":  float(arr.max())}
        return d

    B1 = agg(sums_B1)
    B2 = agg(sums_B2)

    # z-scores and one-sided p-values
    def z_p(primes_val, base):
        z = (primes_val - base["mean"]) / max(base["std"], 1e-12)
        # two-sided empirical p-value: how often |b - mean| >= |primes_val - mean|
        return z

    z1 = {k: z_p(s_primes[k], B1[k]) for k in keys}
    z2 = {k: z_p(s_primes[k], B2[k]) for k in keys}

    # rank tests: where in the K-replicate distribution does PRIMES sit?
    def rank(primes_val, rows, k):
        arr = np.array([r[k] for r in rows], dtype=float)
        return float((arr <= primes_val).sum()) / float(len(arr))

    r1 = {k: rank(s_primes[k], sums_B1, k) for k in keys}
    r2 = {k: rank(s_primes[k], sums_B2, k) for k in keys}

    out = {
        "params": {
            "N_max": args.N_max,
            "start_x": args.start_x,
            "p_at_start": p_at_start,
            "M": M,
            "d": args.d,
            "tau": 1,
            "max_edge": args.max_edge,
            "K_baselines": args.K,
            "tau_thr_for_N1": 0.5,
            "seed": args.seed,
        },
        "primes": s_primes,
        "B1_iid_exp1": B1,
        "B2_gap_permuted": B2,
        "z_score_vs_B1": z1,
        "z_score_vs_B2": z2,
        "rank_in_B1": r1,
        "rank_in_B2": r2,
        "raw_B1": sums_B1,
        "raw_B2": sums_B2,
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
    ap.add_argument("--start-x", type=int, default=1_000_000,
                    help="window start (nearest prime)")
    ap.add_argument("--M", type=int, default=2000,
                    help="window length (number of normalised gaps)")
    ap.add_argument("--d", type=int, default=3,
                    help="Takens embedding dimension")
    ap.add_argument("--max-edge", type=float, default=4.0,
                    help="ripser max edge length (filtration cap)")
    ap.add_argument("--K", type=int, default=20,
                    help="number of baseline replicates per arm")
    ap.add_argument("--seed", type=int, default=20260426)
    ap.add_argument("--out",
                    default="persistent_homology_chi_p_data.json")
    args = ap.parse_args()

    out = run(args)
    Path(args.out).write_text(json.dumps(out, indent=2))
    print("\n[done] wrote", args.out)
    print("PRIMES :", out["primes"])
    print("B1 mean:", {k: out["B1_iid_exp1"][k]["mean"] for k in out["primes"]})
    print("B2 mean:", {k: out["B2_gap_permuted"][k]["mean"] for k in out["primes"]})
    print("z(P;B1):", out["z_score_vs_B1"])
    print("z(P;B2):", out["z_score_vs_B2"])
    print("rank in B1:", out["rank_in_B1"])
    print("rank in B2:", out["rank_in_B2"])


if __name__ == "__main__":
    main()

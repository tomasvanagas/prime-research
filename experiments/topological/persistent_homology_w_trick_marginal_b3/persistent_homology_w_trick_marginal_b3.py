#!/usr/bin/env python3
"""
D2.a.1 — Persistent homology of W-tricked normalised prime gaps vs
continuous-marginal-matched baseline B3
=====================================================================

Successor to D2.a (S117, persistent_homology_w_trick.py).

Question (NOVELTY_CHALLENGES.md §D2.a.1, proposed S117): the W-trick
at W=210 erases the SERIAL component of the PH deficit on Cramer-
normalised prime gaps (z vs B2 collapses from |z|≈7-12 to |z|<=3).
Against B1 = IID Exp(1) the deficit persists (|z|≈7-12 with T1 sign-
flipping +). S117 attributed this to the W-tricked gap MARGINAL
itself differing from Exp(1) (gaps quantised to multiples of W=210 →
discrete spectrum).

But which sub-component of the marginal deviation drives the deficit?

  (i) the SHAPE of the marginal (different mean/variance/tail than
      Exp(1)); or
  (ii) the DISCRETENESS (support on quasi-grid of spacing ~ 0.318)?

Construction B3 = IID samples from a CONTINUOUS distribution matched
to the W-tricked empirical marginal (linearly-interpolated empirical
CDF, equivalently inverse-transform sampling). Both (i) and B3 share
the continuous envelope; (i) has discreteness, B3 does not.

Pre-stated falsifiers in results.md.

Cross-domain refs:
- Carlsson 2009 BAMS 46.
- Edelsbrunner-Harer 2010.
- Bauer 2021 (ripser).
- Devroye 1986 (inverse-transform sampling, §II.2.1).
"""

import argparse
import json
import math
import time
from pathlib import Path

import numpy as np
from ripser import ripser


def segmented_sieve(N: int) -> np.ndarray:
    sieve = np.ones(N + 1, dtype=bool)
    sieve[:2] = False
    for i in range(2, int(N**0.5) + 1):
        if sieve[i]:
            sieve[i*i::i] = False
    return np.flatnonzero(sieve).astype(np.int64)


def phi(n: int) -> int:
    """Euler totient."""
    result = n
    nn = n
    i = 2
    while i * i <= nn:
        if nn % i == 0:
            while nn % i == 0:
                nn //= i
            result -= result // i
        i += 1
    if nn > 1:
        result -= result // nn
    return result


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


def w_tricked_normalised_gaps(primes: np.ndarray, W: int, b: int,
                              start_x: int, M: int) -> tuple:
    phi_W = phi(W)
    mask = (primes % W) == b
    sub = primes[mask]
    if len(sub) < M + 2:
        raise ValueError(f"only {len(sub)} primes ~ {b} mod {W}; need >= {M+2}")
    start_idx = int(np.searchsorted(sub, start_x))
    if start_idx + M + 1 >= len(sub):
        raise ValueError(f"not enough primes mod {W} after start_x")
    q = sub[start_idx : start_idx + M + 1]
    g = np.diff(q).astype(float)
    log_q = np.log(q[:-1].astype(float))
    x = g / (phi_W * log_q)
    return x, q


def quantile_sample_from_ecdf(empirical_x: np.ndarray, n_samples: int,
                              rng: np.random.Generator) -> np.ndarray:
    """Inverse-transform sampling from the LINEARLY-INTERPOLATED empirical
    CDF of empirical_x.

    Algorithm: sort the input, draw u ~ Uniform(0, 1), map to position
    pos = u*(n-1) in [0, n-1], let i = floor(pos), frac = pos - i, return
    sorted_x[i] + frac*(sorted_x[i+1] - sorted_x[i]).

    Properties:
      - Distribution is CONTINUOUS (piecewise-uniform between sorted values).
      - Marginal CDF matches the empirical CDF at the sorted values.
      - Smooths the discrete spikes in the empirical density into a
        piecewise-constant density over [sorted_x[0], sorted_x[-1]].
      - Supports same range as empirical_x.

    This is the cleanest "continuous distribution matched to the
    empirical marginal" — no kernel choice, no bandwidth, no parametric
    family. Ref: Devroye 1986 *Non-Uniform Random Variate Generation*
    §II.2.1 (Inversion method).
    """
    sorted_x = np.sort(empirical_x)
    n = len(sorted_x)
    u = rng.uniform(0.0, 1.0, size=n_samples)
    pos = u * (n - 1)
    ilow = np.clip(np.floor(pos).astype(int), 0, n - 2)
    frac = pos - ilow
    return sorted_x[ilow] + frac * (sorted_x[ilow + 1] - sorted_x[ilow])


def aggregate_baseline(rows, keys):
    r = {}
    for k in keys:
        arr = np.array([rr[k] for rr in rows], dtype=float)
        r[k] = {"mean": float(arr.mean()),
                "std":  float(arr.std(ddof=1)) if len(arr) > 1 else 0.0,
                "min":  float(arr.min()),
                "max":  float(arr.max())}
    return r


def run(args):
    rng = np.random.default_rng(args.seed)
    t0 = time.time()

    primes = segmented_sieve(args.N_max)
    n_primes = len(primes)
    print(f"[setup] sieve to {args.N_max:,}: {n_primes:,} primes "
          f"in {time.time()-t0:.1f}s")

    W = args.W
    phi_W = phi(W)
    print(f"[W-trick] W={W}, phi(W)={phi_W}, residues used: {args.residues}")

    M = args.M
    d = args.d
    max_edge = args.max_edge
    K = args.K

    per_residue = []
    for b in args.residues:
        if math.gcd(b, W) != 1:
            print(f"  [skip] gcd({b},{W}) != 1; residue invalid")
            continue
        x, q = w_tricked_normalised_gaps(primes, W, b, args.start_x, M)
        print(f"[res b={b}] M={M}, q_start={int(q[0])}, q_end={int(q[-1])}, "
              f"mean={x.mean():.4f}, std={x.std():.4f}, "
              f"#unique={len(np.unique(np.round(x, 4)))}")

        cloud_p = takens_embed(x, d, tau=1)
        t1 = time.time()
        s_primes = ph01(cloud_p, max_edge=max_edge)
        print(f"  [PH primes] {time.time()-t1:.1f}s; T0={s_primes['T0']:.2f}, "
              f"T1={s_primes['T1']:.2f}")

        # B1: IID Exp(1)
        sums_B1 = []
        for k in range(K):
            xb = rng.exponential(scale=1.0, size=M).astype(float)
            cloud = takens_embed(xb, d, tau=1)
            sums_B1.append(ph01(cloud, max_edge=max_edge))

        # B2: gap-permutation (discrete marginal preserved exactly)
        sums_B2 = []
        for k in range(K):
            xb = rng.permutation(x)
            cloud = takens_embed(xb, d, tau=1)
            sums_B2.append(ph01(cloud, max_edge=max_edge))

        # B3 (NEW): IID inverse-CDF samples (continuous envelope)
        sums_B3 = []
        for k in range(K):
            xb = quantile_sample_from_ecdf(x, M, rng).astype(float)
            cloud = takens_embed(xb, d, tau=1)
            sums_B3.append(ph01(cloud, max_edge=max_edge))

        keys = ["T0", "L0", "T1", "L1", "N1"]
        B1 = aggregate_baseline(sums_B1, keys)
        B2 = aggregate_baseline(sums_B2, keys)
        B3 = aggregate_baseline(sums_B3, keys)

        def zof(B):
            return {k: (s_primes[k] - B[k]["mean"]) /
                       max(B[k]["std"], 1e-12)
                    for k in keys}

        def rankof(rows):
            return {k: float(sum(1 for r in rows if r[k] <= s_primes[k])) /
                       float(len(rows))
                    for k in keys}

        z1 = zof(B1)
        z2 = zof(B2)
        z3 = zof(B3)
        rank1 = rankof(sums_B1)
        rank2 = rankof(sums_B2)
        rank3 = rankof(sums_B3)

        per_residue.append({
            "b": b,
            "q_start": int(q[0]),
            "q_end": int(q[-1]),
            "x_mean": float(x.mean()),
            "x_std":  float(x.std()),
            "x_unique_count": int(len(np.unique(np.round(x, 6)))),
            "primes": s_primes,
            "B1": B1, "B2": B2, "B3": B3,
            "z_B1": z1, "z_B2": z2, "z_B3": z3,
            "rank_B1": rank1, "rank_B2": rank2, "rank_B3": rank3,
        })
        print(f"  [z] T0: B1={z1['T0']:+.2f} B2={z2['T0']:+.2f} "
              f"B3={z3['T0']:+.2f}")
        print(f"  [z] T1: B1={z1['T1']:+.2f} B2={z2['T1']:+.2f} "
              f"B3={z3['T1']:+.2f}")

    keys = ["T0", "L0", "T1", "L1", "N1"]
    z_keys = ["z_B1", "z_B2", "z_B3"]
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
        "pooled_z_across_residues": pooled,
        "s117_anchor_main_run": {
            "M": 1000, "d": 3, "x": 1_000_000, "K": 20,
            "residues": [1, 11, 13],
            "z_B1": {"T0": -9.07, "T1": +5.56},
            "z_B2": {"T0": -1.99, "T1": -0.67},
            "ref": "experiments/topological/persistent_homology_w_trick/w_trick_main.json",
        },
        "elapsed_sec": time.time() - t0,
    }
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N-max", type=int, default=10_000_000)
    ap.add_argument("--W", type=int, default=210)
    ap.add_argument("--residues", type=int, nargs="+", default=[1, 11, 13])
    ap.add_argument("--start-x", type=int, default=1_000_000)
    ap.add_argument("--M", type=int, default=1000)
    ap.add_argument("--d", type=int, default=3)
    ap.add_argument("--max-edge", type=float, default=4.0)
    ap.add_argument("--K", type=int, default=20)
    ap.add_argument("--seed", type=int, default=20260427)
    ap.add_argument("--out", default="b3_main.json")
    args = ap.parse_args()

    out = run(args)
    Path(args.out).write_text(json.dumps(out, indent=2))
    print("\n[done] wrote", args.out)
    print("\nPooled z across residues:")
    for k in ["T0", "T1"]:
        for zk in ["z_B1", "z_B2", "z_B3"]:
            m = out["pooled_z_across_residues"][zk][k]
            print(f"  {k} {zk}: mean={m['mean']:+.2f}, "
                  f"range=[{m['min']:+.2f}, {m['max']:+.2f}]")
    print("\nS117 anchor:")
    a = out["s117_anchor_main_run"]
    print(f"  T0 z(B1)={a['z_B1']['T0']}, z(B2)={a['z_B2']['T0']}")
    print(f"  T1 z(B1)={a['z_B1']['T1']}, z(B2)={a['z_B2']['T1']}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
D2.a.2 — Vary W and trace S^(W)_PH (PH analogue of E2.13's Gowers W-scan)
=========================================================================

Successor to D2.a (S117). NOVELTY_CHALLENGES.md §D2.a.2.

Question
--------
S117 showed: at W=210, the W-trick on Cramer-normalised prime gaps
collapses the serial-correlation component of E2.17 (z(B2) on T0/T1
goes from -7.45/-4.05 to <=|2.5|).

Vary W in {2, 6, 30, 210, 2310} (primorials). Predicted: a monotone
decay of |z(B2; W)| in W, parallel to the per-prime structure of the
Hardy-Littlewood singular series. This would give a "PH analogue of
S^(W)_2 / S^(W)_3" matching E2.13's Gowers W-scan.

Protocol
--------
For each W, pick the first min(3, phi(W)) residues b coprime to W
(b in {1}, {1,5}, {1,7,11}, {1,11,13}, {1,13,17} for W in
{2, 6, 30, 210, 2310} respectively). For each (W, b):

  1. Filter primes <= N_max to q congruent to b mod W.
  2. Take M consecutive primes starting at first q >= start_x.
  3. Cramer-normalise: x_n = (q_{n+1} - q_n) / (phi(W) * log q_n).
     E[x_n] -> 1 under Cramer + Dirichlet equidistribution; x_n is
     supported on a discrete grid k * W / (phi(W) * log q_n).
  4. Takens-embed at d = 3, tau = 1.
  5. ripser maxdim=1, thresh=4.0. Compute T0 (sum of finite H_0
     persistences) and T1 (sum of H_1 persistences).
  6. K=20 baselines:
       B1 = IID Exp(1) (continuous, no discreteness)
       B2 = gap-permutation (preserves W-tricked marginal, kills
            serial correlation)
  7. z(stat; B) = (PRIMES_W stat - mean(B stat)) / std(B stat).

Pool z across residues (simple mean) to obtain z(W) on T0 and T1.

Pre-stated falsifiers (in results.md before run).

Cross-domain refs
- Carlsson 2009 BAMS 46.
- Bauer 2021 (ripser).
- Green-Tao 2008 (W-trick origin).
- Gowers 2001 GAFA (Gowers norms).
- Hardy-Littlewood 1923 (singular series).
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


def phi(n: int) -> int:
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


def coprime_residues(W: int, k: int) -> list:
    """First k residues b in [1..W-1] with gcd(b, W) = 1."""
    out = []
    for b in range(1, W):
        if math.gcd(b, W) == 1:
            out.append(b)
            if len(out) == k:
                return out
    return out


# ---------------------------------------------------------------------
# PH summaries
# ---------------------------------------------------------------------

def diagram_summaries(dgm0, dgm1) -> dict:
    out = {}
    if len(dgm0):
        d0 = dgm0[np.isfinite(dgm0[:, 1])]
        if len(d0):
            life0 = d0[:, 1] - d0[:, 0]
            out["T0"] = float(life0.sum())
        else:
            out["T0"] = 0.0
    else:
        out["T0"] = 0.0
    if len(dgm1):
        life1 = dgm1[:, 1] - dgm1[:, 0]
        life1 = life1[np.isfinite(life1)]
        if len(life1):
            out["T1"] = float(life1.sum())
        else:
            out["T1"] = 0.0
    else:
        out["T1"] = 0.0
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
    return diagram_summaries(dgms[0], dgms[1])


# ---------------------------------------------------------------------
# Per-(W, b) experiment
# ---------------------------------------------------------------------

def w_tricked_normalised_gaps(primes: np.ndarray, W: int, b: int,
                              start_x: int, M: int) -> tuple:
    phi_W = phi(W)
    mask = (primes % W) == b
    sub = primes[mask]
    if len(sub) < M + 2:
        raise ValueError(
            f"only {len(sub)} primes ~ {b} mod {W}; need >= {M+2}"
        )
    start_idx = int(np.searchsorted(sub, start_x))
    if start_idx + M + 1 >= len(sub):
        raise ValueError(
            f"not enough primes mod {W} starting at {start_x}; "
            f"have {len(sub) - start_idx} after start, need {M + 1}"
        )
    q = sub[start_idx : start_idx + M + 1]
    g = np.diff(q).astype(float)
    log_q = np.log(q[:-1].astype(float))
    x = g / (phi_W * log_q)
    return x, q


def run_one(primes: np.ndarray, W: int, b: int, start_x: int, M: int,
            d: int, max_edge: float, K: int, rng) -> dict:
    x, q = w_tricked_normalised_gaps(primes, W, b, start_x, M)
    cloud_p = takens_embed(x, d, tau=1)
    s_primes = ph01(cloud_p, max_edge=max_edge)

    sums_B1 = []
    for _ in range(K):
        xb = rng.exponential(scale=1.0, size=M).astype(float)
        cloud = takens_embed(xb, d, tau=1)
        sums_B1.append(ph01(cloud, max_edge=max_edge))
    sums_B2 = []
    for _ in range(K):
        xb = rng.permutation(x)
        cloud = takens_embed(xb, d, tau=1)
        sums_B2.append(ph01(cloud, max_edge=max_edge))

    keys = ["T0", "T1"]
    def agg(rows):
        r = {}
        for kk in keys:
            arr = np.array([rr[kk] for rr in rows], dtype=float)
            r[kk] = {
                "mean": float(arr.mean()),
                "std":  float(arr.std(ddof=1)) if len(arr) > 1 else 0.0,
            }
        return r

    B1 = agg(sums_B1)
    B2 = agg(sums_B2)
    z1 = {kk: (s_primes[kk] - B1[kk]["mean"]) / max(B1[kk]["std"], 1e-12)
          for kk in keys}
    z2 = {kk: (s_primes[kk] - B2[kk]["mean"]) / max(B2[kk]["std"], 1e-12)
          for kk in keys}

    return {
        "W": W, "b": b,
        "phi_W": phi(W),
        "q_start": int(q[0]), "q_end": int(q[-1]),
        "x_mean": float(x.mean()), "x_std": float(x.std()),
        "primes": s_primes,
        "B1": B1, "B2": B2,
        "z_B1": z1, "z_B2": z2,
    }


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------

def run(args):
    rng = np.random.default_rng(args.seed)
    t0 = time.time()
    primes = segmented_sieve(args.N_max)
    print(f"[setup] sieve to {args.N_max:,}: {len(primes):,} primes "
          f"in {time.time()-t0:.1f}s")

    rows = []
    pooled = {}
    for W in args.W_list:
        ph_W = phi(W)
        residues = coprime_residues(W, min(args.n_residues, ph_W))
        print(f"\n[W={W}] phi={ph_W}, residues={residues}")
        per_W = []
        for b in residues:
            t1 = time.time()
            try:
                row = run_one(primes, W, b,
                              start_x=args.start_x,
                              M=args.M, d=args.d,
                              max_edge=args.max_edge, K=args.K,
                              rng=rng)
            except ValueError as exc:
                print(f"  [W={W}, b={b}] SKIP: {exc}")
                continue
            print(f"  [b={b}] q in [{row['q_start']}, {row['q_end']}], "
                  f"PRIMES T0={row['primes']['T0']:.2f}, "
                  f"T1={row['primes']['T1']:.2f}, "
                  f"z(B1)=({row['z_B1']['T0']:+.2f},{row['z_B1']['T1']:+.2f}), "
                  f"z(B2)=({row['z_B2']['T0']:+.2f},{row['z_B2']['T1']:+.2f}), "
                  f"{time.time()-t1:.1f}s")
            per_W.append(row)
            rows.append(row)

        # pool z across residues (simple mean), and primes mean
        if not per_W:
            continue
        keys = ["T0", "T1"]
        pool = {"W": W, "phi_W": ph_W, "n_residues": len(per_W)}
        for kk in keys:
            z1s = [r["z_B1"][kk] for r in per_W]
            z2s = [r["z_B2"][kk] for r in per_W]
            tps = [r["primes"][kk] for r in per_W]
            pool[f"z_B1_{kk}_mean"] = float(np.mean(z1s))
            pool[f"z_B1_{kk}_std"]  = (float(np.std(z1s, ddof=1))
                                       if len(z1s) > 1 else 0.0)
            pool[f"z_B2_{kk}_mean"] = float(np.mean(z2s))
            pool[f"z_B2_{kk}_std"]  = (float(np.std(z2s, ddof=1))
                                       if len(z2s) > 1 else 0.0)
            pool[f"PRIMES_{kk}_mean"] = float(np.mean(tps))
        pooled[W] = pool
        print(f"  [pool] z_B1=({pool['z_B1_T0_mean']:+.2f},"
              f"{pool['z_B1_T1_mean']:+.2f}), "
              f"z_B2=({pool['z_B2_T0_mean']:+.2f},"
              f"{pool['z_B2_T1_mean']:+.2f})")

    out = {
        "params": {
            "N_max": args.N_max,
            "W_list": list(args.W_list),
            "start_x": args.start_x,
            "M": args.M, "d": args.d, "tau": 1,
            "max_edge": args.max_edge,
            "K_baselines": args.K,
            "n_residues_max": args.n_residues,
            "seed": args.seed,
        },
        "per_residue_rows": rows,
        "pooled_by_W": pooled,
        "anchor_S96": {
            "M": 1000, "d": 3, "x": 1_000_000,
            "z_B1": {"T0": -5.96, "T1": -2.58},
            "z_B2": {"T0": -7.45, "T1": -4.05},
            "ref": "experiments/topological/persistent_homology_chi_p/",
        },
        "anchor_S117_W_210": {
            "M": 1000, "d": 3, "x": 1_000_000,
            "residues": [1, 11, 13],
            "z_B1": {"T0": -9.07, "T1": +5.56},
            "z_B2": {"T0": -1.99, "T1": -0.67},
            "ref": "experiments/topological/persistent_homology_w_trick/",
        },
        "elapsed_sec": time.time() - t0,
    }
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N-max", type=int, default=10_000_000)
    ap.add_argument("--W-list", type=int, nargs="+",
                    default=[2, 6, 30, 210, 2310])
    ap.add_argument("--start-x", type=int, default=1_000_000)
    ap.add_argument("--M", type=int, default=1000)
    ap.add_argument("--d", type=int, default=3)
    ap.add_argument("--max-edge", type=float, default=4.0)
    ap.add_argument("--K", type=int, default=20)
    ap.add_argument("--n-residues", type=int, default=3,
                    help="max residues per W (uses min(this, phi(W)))")
    ap.add_argument("--seed", type=int, default=20260427)
    ap.add_argument("--out", default="w_scan.json")
    args = ap.parse_args()

    out = run(args)
    Path(args.out).write_text(json.dumps(out, indent=2))
    print("\n[done] wrote", args.out)
    print("\nW-scan summary (pooled z across residues):")
    print(f"  {'W':>5} {'phi(W)':>6} {'n_b':>3}   "
          f"{'T0 z(B1)':>9} {'T0 z(B2)':>9}   "
          f"{'T1 z(B1)':>9} {'T1 z(B2)':>9}")
    for W in args.W_list:
        if W not in out["pooled_by_W"]:
            continue
        p = out["pooled_by_W"][W]
        print(f"  {W:>5} {p['phi_W']:>6} {p['n_residues']:>3}   "
              f"{p['z_B1_T0_mean']:>+9.2f} {p['z_B2_T0_mean']:>+9.2f}   "
              f"{p['z_B1_T1_mean']:>+9.2f} {p['z_B2_T1_mean']:>+9.2f}")
    print(f"  S96 unconditioned baseline (W=1):  "
          f"T0 z(B1)={out['anchor_S96']['z_B1']['T0']:+.2f}, "
          f"T0 z(B2)={out['anchor_S96']['z_B2']['T0']:+.2f}, "
          f"T1 z(B1)={out['anchor_S96']['z_B1']['T1']:+.2f}, "
          f"T1 z(B2)={out['anchor_S96']['z_B2']['T1']:+.2f}")


if __name__ == "__main__":
    main()

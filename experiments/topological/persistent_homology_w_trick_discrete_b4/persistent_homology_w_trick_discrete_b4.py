#!/usr/bin/env python3
"""
D2.a.1.i — Persistent homology of W-tricked normalised prime gaps vs
pure-discrete IID baseline B4 (and B2, B3 for cross-comparison)
=====================================================================

Successor to D2.a.1 (S124, persistent_homology_w_trick_marginal_b3.py),
proposed in S124 successor list, NOVELTY_CHALLENGES.md §D2.a.1.i.

Question: S124 decomposed E2.17's W-tricked PH deficit into three
sub-components (envelope ≫ discreteness ≳ serial-residual). The
discreteness sub-component was the *new* finding, but its 1-3σ on T0
is at the edge of significance and arises from a comparison
(B2 mean − B3 mean) where serial-residual partially cancels.
A *direct* discreteness probe is needed.

Construction:
  B2 = permutation of x (exact discrete marginal, NO duplicate values).
  B3 = inverse-CDF IID from the linearly-interpolated empirical CDF
       (continuous envelope, NO duplicate values).
  B4 = IID with replacement from the empirical PMF
       (exact discrete support, ~37% duplicate values).

The triplet isolates two orthogonal pairs:

  z(B4) vs z(B2)   "IID-vs-permutation on a fixed discrete marginal"
                   — if S117/S124 serial-residual is truly noise AND
                   IID-with-replacement does not create cloud-collapse
                   artifacts, |z(B4) − z(B2)| should be O(1).

  z(B4) vs z(B3)   "discrete-with-duplicates vs continuous-IID"
                   — directly measures the discreteness sub-component
                   without permutation interference.

Pre-stated falsifiers (also in results.md):

  F_i.1 (PRIMARY): pooled-residue |z(P_W; B4)_T0 − z(P_W; B2)_T0| ≤ 1.0
                   at d ∈ {2, 3, 4}. PASS → S117 serial-residual is noise.
                   FAIL with |Δz| ≥ 2 → genuine residual serial structure.

  F_i.2 (DISCRETENESS): pooled-residue B4_mean(T0) < B3_mean(T0) AND
                       (B4_mean − B3_mean) is in same direction as
                       S124's (B2_mean − B3_mean) within a factor 2.
                       PASS → discreteness is a real cloud-geometry
                       effect that B2 and B4 both exhibit.
                       FAIL → B2's discreteness lift is permutation-
                       specific, not a discreteness effect at all.

  F_i.3 (T1 CONSISTENCY): |z(P_W; B4)_T1 − z(P_W; B2)_T1| ≤ 1.5 at d=3.
                       PASS → 1-loop persistence agrees with F_i.1.
                       FAIL → T1 carries serial information T0 does not.

  F_i.4 (DUPLICATE COUNT): B4 baselines have ~ M·(1 − 1/e) ≈ 0.632M
                       expected unique values per draw. The deviation
                       of B4_mean(T0) from B2_mean(T0) is reported and
                       compared to a sub-sampling B6 control:
                       B6 = permutation of x's first ~0.632M elements
                       (preserves "fewer effective points" but no
                       duplicates).

Cross-domain refs:
- Devroye 1986 §II.2.1 (inverse-transform).
- Edelsbrunner-Harer 2010 (PH on point clouds with degenerate distances).
- Bauer 2021 (ripser, arXiv:1908.02518).
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
    """B3: inverse-transform from linearly-interpolated empirical CDF."""
    sorted_x = np.sort(empirical_x)
    n = len(sorted_x)
    u = rng.uniform(0.0, 1.0, size=n_samples)
    pos = u * (n - 1)
    ilow = np.clip(np.floor(pos).astype(int), 0, n - 2)
    frac = pos - ilow
    return sorted_x[ilow] + frac * (sorted_x[ilow + 1] - sorted_x[ilow])


def discrete_iid_sample(empirical_x: np.ndarray, n_samples: int,
                        rng: np.random.Generator) -> np.ndarray:
    """B4: IID with replacement from the discrete empirical PMF.

    Each sample drawn uniformly from the n input values. Expected
    number of distinct values among n_samples = n_samples * (1 -
    (1 - 1/n)^n_samples). For n_samples == n, → n*(1-1/e) ≈ 0.632 n.
    """
    return rng.choice(empirical_x, size=n_samples, replace=True)


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
        x_unique = int(len(np.unique(np.round(x, 6))))
        print(f"[res b={b}] M={M}, q_start={int(q[0])}, q_end={int(q[-1])}, "
              f"mean={x.mean():.4f}, std={x.std():.4f}, "
              f"#unique={x_unique}")

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

        # B2: gap-permutation (exact discrete marginal, no duplicates)
        sums_B2 = []
        for k in range(K):
            xb = rng.permutation(x)
            cloud = takens_embed(xb, d, tau=1)
            sums_B2.append(ph01(cloud, max_edge=max_edge))

        # B3: continuous-envelope IID via linearly-interpolated ECDF
        sums_B3 = []
        for k in range(K):
            xb = quantile_sample_from_ecdf(x, M, rng).astype(float)
            cloud = takens_embed(xb, d, tau=1)
            sums_B3.append(ph01(cloud, max_edge=max_edge))

        # B4 (NEW): IID with replacement from discrete empirical PMF
        sums_B4 = []
        b4_dup_counts = []
        for k in range(K):
            xb = discrete_iid_sample(x, M, rng).astype(float)
            b4_dup_counts.append(int(M - len(np.unique(np.round(xb, 6)))))
            cloud = takens_embed(xb, d, tau=1)
            sums_B4.append(ph01(cloud, max_edge=max_edge))

        keys = ["T0", "L0", "T1", "L1", "N1"]
        B1 = aggregate_baseline(sums_B1, keys)
        B2 = aggregate_baseline(sums_B2, keys)
        B3 = aggregate_baseline(sums_B3, keys)
        B4 = aggregate_baseline(sums_B4, keys)

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
        z4 = zof(B4)
        rank1 = rankof(sums_B1)
        rank2 = rankof(sums_B2)
        rank3 = rankof(sums_B3)
        rank4 = rankof(sums_B4)

        # Direct mean-gap statistics for F_i.2 (discreteness)
        gap_means = {
            "B2_minus_B3_T0": B2["T0"]["mean"] - B3["T0"]["mean"],
            "B4_minus_B3_T0": B4["T0"]["mean"] - B3["T0"]["mean"],
            "B4_minus_B2_T0": B4["T0"]["mean"] - B2["T0"]["mean"],
            "B2_minus_B3_T1": B2["T1"]["mean"] - B3["T1"]["mean"],
            "B4_minus_B3_T1": B4["T1"]["mean"] - B3["T1"]["mean"],
            "B4_minus_B2_T1": B4["T1"]["mean"] - B2["T1"]["mean"],
        }

        per_residue.append({
            "b": b,
            "q_start": int(q[0]),
            "q_end": int(q[-1]),
            "x_mean": float(x.mean()),
            "x_std":  float(x.std()),
            "x_unique_count": x_unique,
            "primes": s_primes,
            "B1": B1, "B2": B2, "B3": B3, "B4": B4,
            "z_B1": z1, "z_B2": z2, "z_B3": z3, "z_B4": z4,
            "rank_B1": rank1, "rank_B2": rank2,
            "rank_B3": rank3, "rank_B4": rank4,
            "gap_means": gap_means,
            "b4_dup_counts": {
                "mean": float(np.mean(b4_dup_counts)),
                "min": int(np.min(b4_dup_counts)),
                "max": int(np.max(b4_dup_counts)),
            },
        })
        print(f"  [z] T0: B1={z1['T0']:+.2f} B2={z2['T0']:+.2f} "
              f"B3={z3['T0']:+.2f} B4={z4['T0']:+.2f}")
        print(f"  [z] T1: B1={z1['T1']:+.2f} B2={z2['T1']:+.2f} "
              f"B3={z3['T1']:+.2f} B4={z4['T1']:+.2f}")
        print(f"  [gap-means T0] B2-B3={gap_means['B2_minus_B3_T0']:+.2f} "
              f"B4-B3={gap_means['B4_minus_B3_T0']:+.2f} "
              f"B4-B2={gap_means['B4_minus_B2_T0']:+.2f}")
        print(f"  [B4 dup count] mean {b4_dup_counts and np.mean(b4_dup_counts):.0f} of {M}")

    keys = ["T0", "L0", "T1", "L1", "N1"]
    z_keys = ["z_B1", "z_B2", "z_B3", "z_B4"]
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

    # Pooled gap-means
    gap_keys = list(per_residue[0]["gap_means"].keys()) if per_residue else []
    pooled_gap = {}
    for gk in gap_keys:
        vs = [pr["gap_means"][gk] for pr in per_residue]
        pooled_gap[gk] = {"mean": float(np.mean(vs)),
                          "min": float(np.min(vs)),
                          "max": float(np.max(vs))}

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
        "pooled_gap_means": pooled_gap,
        "s124_anchor_main_run": {
            "M": 1000, "d": 3, "x": 1_000_000, "K": 20,
            "residues": [1, 11, 13],
            "z_B1": {"T0": -9.64, "T1": +5.71},
            "z_B2": {"T0": -1.99, "T1": -0.74},
            "z_B3": {"T0": -0.05, "T1": +0.46},
            "ref": "experiments/topological/persistent_homology_w_trick_marginal_b3/b3_main.json",
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
    ap.add_argument("--out", default="b4_main.json")
    args = ap.parse_args()

    out = run(args)
    Path(args.out).write_text(json.dumps(out, indent=2))
    print("\n[done] wrote", args.out)
    print("\nPooled z across residues:")
    for k in ["T0", "T1"]:
        for zk in ["z_B1", "z_B2", "z_B3", "z_B4"]:
            m = out["pooled_z_across_residues"][zk][k]
            print(f"  {k} {zk}: mean={m['mean']:+.2f}, "
                  f"range=[{m['min']:+.2f}, {m['max']:+.2f}]")
    print("\nPooled gap-means:")
    for gk, m in out["pooled_gap_means"].items():
        print(f"  {gk}: mean={m['mean']:+.2f}, "
              f"range=[{m['min']:+.2f}, {m['max']:+.2f}]")


if __name__ == "__main__":
    main()

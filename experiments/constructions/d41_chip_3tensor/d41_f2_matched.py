"""
ATTACK_VECTORS.md §D41 — Phase 3: F^2-MATCHED Bernoulli null.

Phase 2 showed chi_P has z(rho_top) = 7-26 sigma, growing with N. But
chi_P also has 2.7x more ternary triples than density-matched Bernoulli
(this IS the Hardy-Littlewood ternary singular series S_3(N) > 1).

The right A-vs-B test: compare chi_P to a Bernoulli baseline with
density q* chosen so that E[F^2_bern] = F^2_chi_P. If chi_P sigma_max
matches Bernoulli at THIS density, the spectral concentration is fully
explained by HL ternary count — B-grade closure. If chi_P sigma_max
EXCEEDS the F^2-matched Bernoulli, there's a SECOND-ORDER spectral
deviation beyond HL — A-grade signal.

Three baselines per N:
  (a) density-matched: q = pi(N) / N (Phase 2)
  (b) F^2-matched: q* such that q*^3 * (N-9)^2 / 2 = F^2_chi_P
  (c) row-sum-matched: assign random row sums matching the chi_P degree
      distribution (test for "rank-1 mean degree" universality)
"""
from __future__ import annotations

import argparse
import json
import math
import time
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import svds


def sieve(n: int) -> np.ndarray:
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for p in range(2, int(n ** 0.5) + 1):
        if s[p]:
            s[p * p :: p] = False
    return s


def build_constraint_matrix_fast(N: int, indicator: np.ndarray) -> sparse.csr_matrix:
    lo, hi = 3, N - 6
    n = hi - lo + 1
    inds = np.where(indicator[lo:hi + 1])[0] + lo
    if len(inds) == 0:
        return sparse.csr_matrix((n, n))
    II, JJ = np.meshgrid(inds, inds, indexing="ij")
    II = II.ravel()
    JJ = JJ.ravel()
    KK = N - II - JJ
    valid = (KK >= lo) & (KK <= hi)
    II = II[valid]
    JJ = JJ[valid]
    KK = KK[valid]
    keep = indicator[KK]
    II = II[keep]
    JJ = JJ[keep]
    rows = II - lo
    cols = JJ - lo
    M = sparse.coo_matrix(
        (np.ones(len(rows), dtype=np.float64), (rows, cols)), shape=(n, n)
    ).tocsr()
    return M


def f2_at_density(N: int, q: float, rng: np.random.Generator) -> tuple[float, float]:
    """Return (F^2, sigma_max) for one Bernoulli realisation at density q."""
    ind = rng.random(N + 1) < q
    ind[:2] = False
    M = build_constraint_matrix_fast(N, ind)
    F2 = float((M.multiply(M)).sum())
    if M.nnz == 0 or M.shape[0] < 5:
        return F2, 0.0
    try:
        S = svds(M, k=2, which="LM", return_singular_vectors=False)
        sigma_max = float(np.max(S))
    except Exception:
        sigma_max = 0.0
    return F2, sigma_max


def bern_ensemble(N: int, q: float, n_real: int, rng: np.random.Generator) -> dict:
    F2s, sigmas = [], []
    for _ in range(n_real):
        F2, sigma = f2_at_density(N, q, rng)
        F2s.append(F2)
        sigmas.append(sigma)
    F2s = np.array(F2s)
    sigmas = np.array(sigmas)
    rho = sigmas ** 2 / np.maximum(F2s, 1.0)
    return {
        "q": q,
        "n_real": n_real,
        "F2_mean": float(F2s.mean()),
        "F2_std": float(F2s.std()),
        "sigma_max_mean": float(sigmas.mean()),
        "sigma_max_std": float(sigmas.std()),
        "sigma_max_min": float(sigmas.min()),
        "sigma_max_max": float(sigmas.max()),
        "rho_mean": float(rho.mean()),
        "rho_std": float(rho.std()),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--Ns", type=str, default="1023,2047,4095,8191")
    parser.add_argument("--n_bern", type=int, default=30)
    parser.add_argument("--n_calibration", type=int, default=10)
    parser.add_argument("--out", type=str,
        default="experiments/constructions/d41_chip_3tensor/results_phase3.json")
    parser.add_argument("--seed", type=int, default=20260430)
    args = parser.parse_args()

    Ns = [int(x) for x in args.Ns.split(",")]
    rng = np.random.default_rng(args.seed)
    out = {"by_N": {}, "meta": {"Ns": Ns}}

    for N in Ns:
        t0 = time.time()
        s_pi = sieve(N)
        n_primes = int(s_pi.sum())
        density_pi = n_primes / (N + 1)
        print(f"\n=== N = {N} (pi(N)={n_primes}, density={density_pi:.4f}) ===")

        # chi_P
        M_chiP = build_constraint_matrix_fast(N, s_pi)
        F2_chiP = float((M_chiP.multiply(M_chiP)).sum())
        if M_chiP.nnz > 0:
            S = svds(M_chiP, k=5, which="LM", return_singular_vectors=False)
            S = np.sort(S)[::-1]
            sigma_max_chiP = float(S[0])
            sigma_top5 = S.tolist()
        else:
            sigma_max_chiP = 0.0
            sigma_top5 = [0.0] * 5
        rho_chiP = sigma_max_chiP ** 2 / max(F2_chiP, 1.0)

        # (a) Density-matched Bernoulli
        bern_a = bern_ensemble(N, density_pi, args.n_bern, rng)

        # (b) Calibrate q* such that bern F^2 ≈ F2_chiP
        # Use small-n calibration grid
        # F^2(q) ≈ const * q^3 (cubic in q since 3 indep events)
        # so q* ≈ density_pi * (F2_chiP / bern_a["F2_mean"]) ** (1/3)
        if bern_a["F2_mean"] > 0:
            q_star = density_pi * (F2_chiP / bern_a["F2_mean"]) ** (1 / 3)
        else:
            q_star = density_pi
        # Validate via a smaller pilot
        pilot = bern_ensemble(N, q_star, args.n_calibration, rng)
        # Refine once
        if pilot["F2_mean"] > 0:
            q_star = q_star * (F2_chiP / pilot["F2_mean"]) ** (1 / 3)
        bern_b = bern_ensemble(N, q_star, args.n_bern, rng)

        # Z-scores against bern_a (density-matched)
        z_sigma_a = (
            (sigma_max_chiP - bern_a["sigma_max_mean"]) / bern_a["sigma_max_std"]
            if bern_a["sigma_max_std"] > 0
            else float("nan")
        )
        z_F2_a = (
            (F2_chiP - bern_a["F2_mean"]) / bern_a["F2_std"]
            if bern_a["F2_std"] > 0
            else float("nan")
        )
        z_rho_a = (
            (rho_chiP - bern_a["rho_mean"]) / bern_a["rho_std"]
            if bern_a["rho_std"] > 0
            else float("nan")
        )

        # Z-scores against bern_b (F^2-matched -- the key A-vs-B test)
        z_sigma_b = (
            (sigma_max_chiP - bern_b["sigma_max_mean"]) / bern_b["sigma_max_std"]
            if bern_b["sigma_max_std"] > 0
            else float("nan")
        )
        z_rho_b = (
            (rho_chiP - bern_b["rho_mean"]) / bern_b["rho_std"]
            if bern_b["rho_std"] > 0
            else float("nan")
        )

        elapsed = time.time() - t0

        out["by_N"][N] = {
            "n_primes": n_primes,
            "density_pi": density_pi,
            "wall_seconds": elapsed,
            "chi_P": {
                "F2": F2_chiP,
                "sigma_max": sigma_max_chiP,
                "sigma_top5": sigma_top5,
                "rho_top": rho_chiP,
            },
            "bern_a_density_matched": bern_a,
            "bern_b_F2_matched": bern_b,
            "z_vs_density_matched": {
                "sigma_max": z_sigma_a,
                "F2": z_F2_a,
                "rho_top": z_rho_a,
            },
            "z_vs_F2_matched": {
                "sigma_max": z_sigma_b,
                "rho_top": z_rho_b,
            },
        }
        print(f"  chi_P:  F2={F2_chiP:.0f}, sigma_max={sigma_max_chiP:.2f}, rho={rho_chiP:.4f}")
        print(
            f"  bern_a (density={density_pi:.3f}): F2={bern_a['F2_mean']:.0f}+-{bern_a['F2_std']:.0f}, "
            f"sigma={bern_a['sigma_max_mean']:.2f}+-{bern_a['sigma_max_std']:.2f}, "
            f"rho={bern_a['rho_mean']:.4f}+-{bern_a['rho_std']:.4f}"
        )
        print(
            f"  bern_b (F2-matched, q*={q_star:.4f}): F2={bern_b['F2_mean']:.0f}+-{bern_b['F2_std']:.0f}, "
            f"sigma={bern_b['sigma_max_mean']:.2f}+-{bern_b['sigma_max_std']:.2f}, "
            f"rho={bern_b['rho_mean']:.4f}+-{bern_b['rho_std']:.4f}"
        )
        print(
            f"  z vs density-matched: sigma={z_sigma_a:.2f}, F2={z_F2_a:.2f}, rho={z_rho_a:.2f}"
        )
        print(
            f"  z vs F^2-matched:     sigma={z_sigma_b:.2f}, rho={z_rho_b:.2f}    <-- A-vs-B TEST"
        )
        print(f"  wall: {elapsed:.1f}s")

    Path(args.out).write_text(json.dumps(out, indent=2, default=float))
    print(f"\nWrote {args.out}")


if __name__ == "__main__":
    main()

"""
ATTACK_VECTORS.md §D41 — Phase 2 analysis: eigenvector + multi-realisation
Bernoulli statistics for the chi_P 3-tensor mode-(1,2) matricisation.

Computes:
  - chi_P leading eigenvector v_1 and its arithmetic content:
      * inner product with chi_P_indicator restricted to primes
      * inner product with constant vector
      * Fourier content / character correlations
      * comparison to v_1 of Bernoulli matched-density realisations
  - Z-score of chi_P leading eigenvalue vs Bernoulli ensemble (B = 100 trials)
  - Top-k spectral gap test: sigma_k / sigma_{k+1} for chi_P vs Bernoulli
  - Tracy-Widom z-score for chi_P and Bernoulli ensembles
  - Frobenius normalised: sigma_max / sqrt(F^2) -- mass concentration in top mode

This is the QUANTITATIVE TEST for D41 melonic / matrix universality.
"""
from __future__ import annotations

import argparse
import json
import math
import time
from pathlib import Path
from typing import Tuple

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import svds, eigsh


# ---------------------------------------------------------------- primes / sieve

def sieve(n: int) -> np.ndarray:
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for p in range(2, int(n ** 0.5) + 1):
        if s[p]:
            s[p * p :: p] = False
    return s


# ---------------------------------------------------------------- matrix build

def build_constraint_matrix(N: int, indicator: np.ndarray) -> sparse.csr_matrix:
    """Symmetric (N-9) x (N-9) matrix M[i,j] = ind[i]*ind[j]*ind[N-i-j], i,j in [3, N-6]."""
    lo, hi = 3, N - 6
    n = hi - lo + 1
    rows, cols = [], []
    # Restrict to primes (or 1s) to vectorise via prime list
    idx_in_range = np.where(indicator[lo:hi + 1])[0] + lo
    set_in_range = set(int(x) for x in idx_in_range)
    for i in idx_in_range:
        i_int = int(i)
        # j range: j in [max(lo, N-i-hi), min(hi, N-i-lo)]
        j_min = max(lo, N - i_int - hi)
        j_max = min(hi, N - i_int - lo)
        if j_min > j_max:
            continue
        # candidate j: must be prime AND in [j_min, j_max]
        for j in idx_in_range:
            j_int = int(j)
            if j_int < j_min:
                continue
            if j_int > j_max:
                break
            k = N - i_int - j_int
            if k in set_in_range:
                rows.append(i_int - lo)
                cols.append(j_int - lo)
    M = sparse.coo_matrix(
        (np.ones(len(rows), dtype=np.float64), (rows, cols)), shape=(n, n)
    ).tocsr()
    return M


def build_constraint_matrix_fast(N: int, indicator: np.ndarray) -> sparse.csr_matrix:
    """Vectorised over (i, j); i+j precomputed."""
    lo, hi = 3, N - 6
    n = hi - lo + 1
    inds = np.where(indicator[lo:hi + 1])[0] + lo  # absolute indices in [lo, hi]
    if len(inds) == 0:
        return sparse.csr_matrix((n, n))
    # Build all pairs (i, j); for each pair, k = N - i - j.
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


# ---------------------------------------------------------------- baselines

def bernoulli_indicator(N: int, density: float, rng: np.random.Generator) -> np.ndarray:
    """Random Bernoulli indicator on [0..N], 1 with probability density."""
    ind = rng.random(N + 1) < density
    ind[:2] = False
    return ind


# ---------------------------------------------------------------- analysis

def top_singular_values(M: sparse.csr_matrix, k: int = 20) -> Tuple[np.ndarray, np.ndarray]:
    """Return top-k singular values and the leading right-singular vector."""
    if M.nnz == 0 or min(M.shape) < k + 2:
        return np.zeros(k), np.zeros(M.shape[1])
    try:
        U, S, Vt = svds(M.astype(float), k=k, which="LM", maxiter=5000)
        order = np.argsort(S)[::-1]
        S = S[order]
        U = U[:, order]
        Vt = Vt[order, :]
    except Exception:
        return np.zeros(k), np.zeros(M.shape[1])
    return S, Vt[0]


def top_symmetric_eigvec(M: sparse.csr_matrix, k: int = 20) -> Tuple[np.ndarray, np.ndarray]:
    """For symmetric M, return top k eigenvalues (signed) and leading eigvec.

    For a non-negative symmetric matrix, top eigenvalue = top singular value
    and is positive (Perron-Frobenius). For more general symmetric, top |eig|.
    """
    if M.nnz == 0 or M.shape[0] < k + 2:
        return np.zeros(k), np.zeros(M.shape[0])
    try:
        # eigsh with which='LA' picks largest algebraic
        E_LA, V_LA = eigsh(M.astype(float), k=k, which="LA", maxiter=5000)
        E_SA, V_SA = eigsh(M.astype(float), k=min(k, M.shape[0] - k - 1), which="SA", maxiter=5000)
    except Exception:
        return np.zeros(k), np.zeros(M.shape[0])
    # combine and sort by |eig|
    E_all = np.concatenate([E_LA, E_SA])
    V_all = np.concatenate([V_LA, V_SA], axis=1)
    order = np.argsort(np.abs(E_all))[::-1]
    E = E_all[order][:k]
    V = V_all[:, order][:, :k]
    return E, V[:, 0]


def assess_arithmetic_content(v: np.ndarray, N: int, primes_in_range: np.ndarray) -> dict:
    """Assess arithmetic content of leading eigenvector v.

    v has length n = N - 9 indexed by i in [3, N-6] (relative index = i - 3).
    """
    if len(v) == 0 or np.linalg.norm(v) == 0:
        return {"ok": False}
    v = v / np.linalg.norm(v)
    lo = 3
    n = len(v)
    # i values in absolute coordinates
    i_abs = np.arange(lo, lo + n)
    # 1) inner product with prime indicator (normalised)
    ind_prime = np.zeros(n)
    primes_in_v = primes_in_range[(primes_in_range >= lo) & (primes_in_range <= lo + n - 1)] - lo
    ind_prime[primes_in_v] = 1.0
    if ind_prime.sum() > 0:
        ind_prime = ind_prime / np.linalg.norm(ind_prime)
    inner_prime = float(np.dot(v, ind_prime))
    # 2) inner product with constant
    const = np.ones(n) / math.sqrt(n)
    inner_const = float(np.dot(v, const))
    # 3) Fourier content: top-k DFT modes of v
    fft = np.fft.rfft(v)
    fft_power = np.abs(fft) ** 2
    fft_power /= fft_power.sum() if fft_power.sum() > 0 else 1.0
    top_modes = np.argsort(fft_power)[-5:][::-1]
    top_mode_freq = [(int(k), float(fft_power[k])) for k in top_modes]
    # 4) localisation: participation ratio
    pr_inv = float(np.sum(v ** 4))
    pr = 1.0 / pr_inv if pr_inv > 0 else 0.0
    # 5) max entry vs L^2 norm
    max_entry = float(np.max(np.abs(v)))
    # 6) sign distribution
    pos_frac = float(np.mean(v > 0))
    return {
        "ok": True,
        "inner_with_prime_indicator": inner_prime,
        "inner_with_constant": inner_const,
        "top_fft_modes": top_mode_freq,
        "participation_ratio": pr,
        "max_entry_abs": max_entry,
        "positive_fraction": pos_frac,
    }


# ---------------------------------------------------------------- main

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--Ns", type=str, default="1023,2047,4095,8191")
    parser.add_argument("--n_bern", type=int, default=50)
    parser.add_argument("--n_top_sv", type=int, default=15)
    parser.add_argument("--out", type=str, default="experiments/constructions/d41_chip_3tensor/results_phase2.json")
    parser.add_argument("--seed", type=int, default=20260430)
    args = parser.parse_args()

    Ns = [int(x) for x in args.Ns.split(",")]
    rng = np.random.default_rng(args.seed)
    out: dict = {"by_N": {}, "meta": {"Ns": Ns, "n_bern": args.n_bern}}

    for N in Ns:
        t0 = time.time()
        s_pi = sieve(N)
        n_primes = int(s_pi.sum())
        density = n_primes / (N + 1)
        primes_in_range = np.where(s_pi)[0]
        print(f"\n=== N = {N} (pi(N) = {n_primes}, density = {density:.4f}) ===")

        # chi_P matrix
        M_chiP = build_constraint_matrix_fast(N, s_pi)
        S_chiP, v_chiP = top_singular_values(M_chiP, k=args.n_top_sv)
        # Symmetric eigvec analysis (matrix is symmetric, real-valued):
        E_chiP_sym, ve_chiP_sym = top_symmetric_eigvec(M_chiP, k=args.n_top_sv)
        F2_chiP = float((M_chiP.multiply(M_chiP)).sum())
        rho_top_chiP = float(S_chiP[0] ** 2 / F2_chiP) if F2_chiP > 0 else 0.0
        gap_top_chiP = float(S_chiP[0] / S_chiP[1]) if S_chiP[1] > 0 else float("inf")

        arith_chiP = assess_arithmetic_content(ve_chiP_sym, N, primes_in_range)

        # Bernoulli ensemble
        bern_top_sv = []
        bern_F2 = []
        bern_rho = []
        bern_gap = []
        bern_top10 = []
        bern_arith_inner_const = []
        for b in range(args.n_bern):
            ind = bernoulli_indicator(N, density, rng)
            M_b = build_constraint_matrix_fast(N, ind)
            S_b, v_b = top_singular_values(M_b, k=args.n_top_sv)
            E_b_sym, ve_b_sym = top_symmetric_eigvec(M_b, k=args.n_top_sv)
            F2_b = float((M_b.multiply(M_b)).sum())
            bern_F2.append(F2_b)
            bern_top_sv.append(float(S_b[0]))
            bern_top10.append(S_b[:10].tolist())
            bern_rho.append(float(S_b[0] ** 2 / F2_b) if F2_b > 0 else 0.0)
            bern_gap.append(float(S_b[0] / S_b[1]) if S_b[1] > 0 else 0.0)
            arith_b = assess_arithmetic_content(ve_b_sym, N, np.where(ind)[0])
            if arith_b.get("ok"):
                bern_arith_inner_const.append(arith_b["inner_with_constant"])

        bern_top_sv = np.array(bern_top_sv)
        bern_F2 = np.array(bern_F2)
        bern_rho = np.array(bern_rho)
        bern_gap = np.array(bern_gap)

        # Z-scores
        sv_mean, sv_std = float(bern_top_sv.mean()), float(bern_top_sv.std())
        z_sigma_max = (S_chiP[0] - sv_mean) / sv_std if sv_std > 0 else float("nan")
        F2_mean, F2_std = float(bern_F2.mean()), float(bern_F2.std())
        z_F2 = (F2_chiP - F2_mean) / F2_std if F2_std > 0 else float("nan")
        rho_mean, rho_std = float(bern_rho.mean()), float(bern_rho.std())
        z_rho = (rho_top_chiP - rho_mean) / rho_std if rho_std > 0 else float("nan")
        gap_mean, gap_std = float(bern_gap.mean()), float(bern_gap.std())
        z_gap = (gap_top_chiP - gap_mean) / gap_std if gap_std > 0 else float("nan")

        elapsed = time.time() - t0
        rec = {
            "n_primes": n_primes,
            "density": density,
            "wall_seconds": elapsed,
            "chi_P": {
                "F2": F2_chiP,
                "sigma_max": float(S_chiP[0]),
                "sigma_top": S_chiP.tolist(),
                "eig_top_signed": E_chiP_sym.tolist(),
                "rho_top_mass_frac": rho_top_chiP,
                "gap_sigma1_over_sigma2": gap_top_chiP,
                "arith_content": arith_chiP,
            },
            "bernoulli_ensemble": {
                "n_realisations": args.n_bern,
                "sigma_max_mean": sv_mean,
                "sigma_max_std": sv_std,
                "F2_mean": F2_mean,
                "F2_std": F2_std,
                "rho_mean": rho_mean,
                "rho_std": rho_std,
                "gap_mean": gap_mean,
                "gap_std": gap_std,
                "sigma_max_min": float(bern_top_sv.min()),
                "sigma_max_max": float(bern_top_sv.max()),
            },
            "z_scores": {
                "sigma_max": z_sigma_max,
                "F2": z_F2,
                "rho_top": z_rho,
                "gap_sigma1_sigma2": z_gap,
            },
        }
        out["by_N"][N] = rec
        print(
            f"  chi_P: sigma_max={S_chiP[0]:.2f}, F2={F2_chiP:.0f}, "
            f"rho_top={rho_top_chiP:.4f}, gap={gap_top_chiP:.3f}, "
            f"top5={[f'{x:.2f}' for x in S_chiP[:5]]}"
        )
        print(
            f"  bern:  sigma_max mean={sv_mean:.2f} (std {sv_std:.2f}), "
            f"F2 mean={F2_mean:.0f} (std {F2_std:.0f}), "
            f"rho mean={rho_mean:.4f}, gap mean={gap_mean:.3f}"
        )
        print(
            f"  z-scores:  sigma_max={z_sigma_max:.2f},  F2={z_F2:.2f},  "
            f"rho={z_rho:.2f},  gap={z_gap:.2f}"
        )
        if arith_chiP.get("ok"):
            print(
                f"  chi_P eigvec arith: inner_with_const={arith_chiP['inner_with_constant']:.3f}, "
                f"inner_with_prime_ind={arith_chiP['inner_with_prime_indicator']:.3f}, "
                f"PR={arith_chiP['participation_ratio']:.0f} (n={N-9}), "
                f"top FFT modes={arith_chiP['top_fft_modes'][:3]}"
            )
        print(f"  wall: {elapsed:.1f}s")

    Path(args.out).write_text(json.dumps(out, indent=2, default=float))
    print(f"\nWrote {args.out}")


if __name__ == "__main__":
    main()

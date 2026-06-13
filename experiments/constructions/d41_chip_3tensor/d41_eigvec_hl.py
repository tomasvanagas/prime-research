"""
ATTACK_VECTORS.md §D41 — Phase 4: leading eigenvector vs HL R_2(N-p).

If the chi_P leading eigenvector v satisfies v_p ~ sqrt(d_p) where
d_p = #{q prime : (p, q, N-p-q) is prime triple} = R_2(N-p) (Goldbach
2-rep at N-p), then the chi_P spectral structure is fully HL-predicted.
This is the standard Perron-Frobenius rank-1 approximation.

Empirical test: rank correlation of v_p^2 with d_p across primes p.
If correlation -> 1 at large N, HL fully explains the spectrum (B-grade).
If correlation < 1 with structured residual, A-grade.

Also: rank correlation of d_p (empirical) with HL S_2(N-p) prediction.
This separates "primes follow HL" (known) from "spectrum follows degree"
(the rank-1 hypothesis).
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy import sparse, stats
from scipy.sparse.linalg import eigsh


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
    II, JJ, KK = II[valid], JJ[valid], KK[valid]
    keep = indicator[KK]
    II, JJ = II[keep], JJ[keep]
    M = sparse.coo_matrix(
        (np.ones(len(II), dtype=np.float64), (II - lo, JJ - lo)), shape=(n, n)
    ).tocsr()
    return M


def hl_S2(M_target: int, p_max: int = 5000) -> float:
    """Hardy-Littlewood Goldbach singular series S_2(M) for even M."""
    if M_target <= 0 or M_target % 2:
        return 0.0
    s = sieve(p_max)
    primes = np.where(s)[0]
    C2 = 1.0
    for p in primes:
        if p < 3:
            continue
        C2 *= 1.0 - 1.0 / (p - 1) ** 2
    out = 2.0 * C2
    m = M_target
    for p in primes:
        if p * p > m:
            break
        if p > 2 and m % p == 0:
            out *= (p - 1) / (p - 2)
            while m % p == 0:
                m //= p
        elif m % p == 0:
            while m % p == 0:
                m //= p
    if m > 2:
        out *= (m - 1) / (m - 2)
    return out


def hl_r2_prediction(M_target: int) -> float:
    """HL prediction: r_2(M) ~ S_2(M) * M / log^2 M for even M >= 6."""
    if M_target < 6 or M_target % 2:
        return 0.0
    return hl_S2(M_target) * M_target / (math.log(M_target) ** 2)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--Ns", type=str, default="2047,8191,16383")
    parser.add_argument("--out", type=str,
        default="experiments/constructions/d41_chip_3tensor/results_phase4.json")
    args = parser.parse_args()

    Ns = [int(x) for x in args.Ns.split(",")]
    out = {"by_N": {}, "meta": {"Ns": Ns}}

    for N in Ns:
        s_pi = sieve(N)
        primes = np.where(s_pi)[0]
        primes_in_range = primes[(primes >= 3) & (primes <= N - 6)]
        M = build_constraint_matrix_fast(N, s_pi)
        # Compute leading eigenvector (symmetric matrix M)
        E, V = eigsh(M.astype(float), k=3, which="LA", maxiter=10000)
        order = np.argsort(np.abs(E))[::-1]
        E = E[order]
        V = V[:, order]
        v1 = V[:, 0]
        v1 = v1 / np.linalg.norm(v1)
        # Empirical degree: row sum of M restricted to primes
        row_sums = np.array(M.sum(axis=1)).ravel()  # length n = N - 9
        # HL R_2(N-p) prediction for each prime p
        hl_d = np.array([hl_r2_prediction(N - int(p)) for p in primes_in_range])
        # Empirical d_p = row sum at index (p - 3)
        emp_d = np.array([row_sums[int(p) - 3] for p in primes_in_range])
        # Eigvec mass at prime p
        v1_at_p = np.array([v1[int(p) - 3] for p in primes_in_range])
        v1_sq = v1_at_p ** 2
        # 1) correlation of v1^2 with d_p (rank-1 Perron-Frobenius hypothesis: v_p ~ sqrt(d_p))
        spear_v1sq_d = float(stats.spearmanr(v1_sq, emp_d).correlation)
        pear_v1sq_d = float(stats.pearsonr(v1_sq, emp_d)[0])
        # 2) correlation of empirical d_p with HL prediction
        spear_d_hl = float(stats.spearmanr(emp_d, hl_d).correlation)
        pear_d_hl = float(stats.pearsonr(emp_d, hl_d)[0])
        # 3) correlation of v1 with sqrt(emp_d)
        sqrt_d = np.sqrt(emp_d)
        if np.linalg.norm(sqrt_d) > 0:
            sqrt_d_norm = sqrt_d / np.linalg.norm(sqrt_d)
        else:
            sqrt_d_norm = sqrt_d
        # cosine similarity (sign-agnostic)
        v1_at_p_norm = v1_at_p / max(np.linalg.norm(v1_at_p), 1e-12)
        cos_v1_sqrtd = float(abs(np.dot(v1_at_p_norm, sqrt_d_norm)))
        # 4) for unsigned eigvec |v1|, cosine with HL-sqrt
        sqrt_hl = np.sqrt(np.maximum(hl_d, 0))
        if np.linalg.norm(sqrt_hl) > 0:
            sqrt_hl_norm = sqrt_hl / np.linalg.norm(sqrt_hl)
            cos_v1_sqrthl = float(abs(np.dot(v1_at_p_norm, sqrt_hl_norm)))
        else:
            cos_v1_sqrthl = 0.0

        out["by_N"][N] = {
            "n_primes": len(primes_in_range),
            "leading_eigvalue": float(E[0]),
            "second_eigvalue": float(E[1]) if len(E) > 1 else 0.0,
            "third_eigvalue": float(E[2]) if len(E) > 2 else 0.0,
            "abs_eig_ratio_2_to_1": float(abs(E[1]) / abs(E[0])) if abs(E[0]) > 0 else 0,
            # eigvec on primes
            "eigvec_mass_on_primes": float(np.sum(v1_sq)),
            "eigvec_localisation_PR": float(1.0 / np.sum(v1 ** 4)),
            # rank-1 PF test
            "spearman_v1sq_vs_emp_d": spear_v1sq_d,
            "pearson_v1sq_vs_emp_d": pear_v1sq_d,
            "cosine_v1_vs_sqrt_emp_d": cos_v1_sqrtd,
            # HL test
            "spearman_emp_d_vs_HL": spear_d_hl,
            "pearson_emp_d_vs_HL": pear_d_hl,
            "cosine_v1_vs_sqrt_HL": cos_v1_sqrthl,
            # summary stats
            "emp_d_mean": float(emp_d.mean()),
            "emp_d_std": float(emp_d.std()),
            "hl_d_mean": float(hl_d.mean()),
            "hl_d_std": float(hl_d.std()),
        }
        print(f"\n=== N = {N} ===")
        print(f"  leading eig: {E[0]:.2f},  second: {E[1]:.2f},  third: {E[2]:.2f}")
        print(f"  eig ratio |E2|/|E1|: {abs(E[1]) / abs(E[0]):.3f}")
        print(f"  eigvec mass on primes: {np.sum(v1_sq):.4f}")
        print(f"  PR localisation: {1.0 / np.sum(v1 ** 4):.1f}  (n_primes={len(primes_in_range)})")
        print(f"  Spearman corr(v1^2, emp_d) = {spear_v1sq_d:.4f}")
        print(f"  Spearman corr(emp_d, HL_d) = {spear_d_hl:.4f}  <-- HL fits primes")
        print(f"  Cosine(v1, sqrt(emp_d))    = {cos_v1_sqrtd:.4f}  <-- rank-1 PF")
        print(f"  Cosine(v1, sqrt(HL_d))     = {cos_v1_sqrthl:.4f}  <-- HL predicts spectrum")

    Path(args.out).write_text(json.dumps(out, indent=2, default=float))
    print(f"\nWrote {args.out}")


if __name__ == "__main__":
    main()

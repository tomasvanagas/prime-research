"""
ATTACK_VECTORS.md §D41 — Gurau-Witten melonic universality of large-N tensor
models on the chi_P 3-tensor.

Construction: T_{ijk}^N = 1[i+j+k = N] * chi_P(i) chi_P(j) chi_P(k)
indexed (i,j,k) in [3, N-6]^3.

Three matricizations / unfoldings are computed and analysed:

  (A) Mode-1 unfolding M_1: (N-9) x (N-9)^2,  M_1[i, (j,k)] = T[i,j,k].
      THEOREM (provable from the constraint i+j+k=N):
      M_1 M_1^T is diagonal with entries r_2(N-i) chi_P(i),
      where r_2(M) = #{(j,k) prime pair : j+k = M, both in valid range}.
      So mode-1 SVD is trivial; spectrum = sqrt(r_2(N-i)) over primes i.

  (B) Constraint-eliminated matricization M_2: (N-9) x (N-9),
      M_2[i, j] = chi_P(i) chi_P(j) chi_P(N-i-j) (with valid k=N-i-j range).
      Non-trivial; the natural rank-2 reduction. SVD bulk + edge tested
      against (i) Marchenko-Pastur, (ii) Tracy-Widom beta=2 edge,
      (iii) Bernoulli matched-density baseline, (iv) Liouville baseline.

  (C) HOSVD core / Tucker rank-1 leading approximation: tests whether
      a rank-1 outer-product approximation captures most of the tensor
      mass (would indicate matrix-rank-2 universality class) or the
      tensor is genuinely high-rank (a melonic-universality signature).

A-grade outcome: a deviation from MP+TW that is NOT explained by
Hardy-Littlewood ternary singular series and is consistent with
Gurau-Witten melonic universality at large N.

B-grade outcome: empirical match to MP+TW or to HL — first rank-3
universality measurement for chi_P, refining E2.13 by an explicit
rank-3 spectral identity.

Failure modes:
  (E)   spectrum is matrix-class (MP / TW) up to HL multiplicative
        prefactor.  A 41st pseudorandomness measure, this time
        rank-3, lands at the noise floor.
  (I)   the leading singular value matches HL R_3(N) prediction with
        no further structure (B-grade).
  (INC) tensor too sparse for clean large-N limit at accessible N.

Cross-domain refs (cited in d41_chip_3tensor_results.md):
  - Gurau 2011 arXiv:1011.2726
  - Witten 2016 arXiv:1610.09758
  - Klebanov-Tarnopolsky 2017 arXiv:1611.08915
  - Bonzom-Gurau-Riello-Rivasseau 2011 arXiv:1105.3122

Edges cited: E2.13 (chi_P U^k = HL singular series), E3.1 (Connes
operator), C7 / E7.1 (GUE up to order 6 — rank-2 zero correlations),
S74 (free cumulants — rank-2).
"""
from __future__ import annotations

import argparse
import json
import math
import time
from dataclasses import dataclass, asdict
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import svds


# ---------------------------------------------------------------- primes


def sieve(n: int) -> np.ndarray:
    """Boolean prime indicator chi_P[0..n], chi_P[i]=1 iff i is prime."""
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for p in range(2, int(n ** 0.5) + 1):
        if s[p]:
            s[p * p :: p] = False
    return s


def liouville(n: int) -> np.ndarray:
    """Liouville lambda(k) = (-1)^Omega(k) on [0..n]."""
    out = np.ones(n + 1, dtype=np.int8)
    out[0] = 0
    for k in range(2, n + 1):
        m, c = k, 0
        d = 2
        while d * d <= m:
            while m % d == 0:
                m //= d
                c += 1
            d += 1
        if m > 1:
            c += 1
        out[k] = -1 if (c & 1) else 1
    return out


# ---------------------------------------------------------------- HL constants


def hl_twin_constant(p_max: int = 10_000) -> float:
    """C_2 = prod_{p>=3} (1 - 1/(p-1)^2). Standard twin-prime constant."""
    s = sieve(p_max)
    primes = np.where(s)[0]
    c = 1.0
    for p in primes:
        if p < 3:
            continue
        c *= 1.0 - 1.0 / (p - 1) ** 2
    return c


def hl_ternary_singular_series(N: int, p_max: int = 10_000) -> float:
    """Vinogradov ternary singular series at N (odd):

        S_3(N) = prod_p [1 + 1/(p-1)^3] * prod_{p|N} [(p-1)^2 - 1] / [(p-1)^2 + (p-1)]
        more standard: S_3(N) = (1/2) prod_p [(1 - 1/(p-1)^2)] * prod_{p|N} [1 + 1/(p^2-3p+3)]
        (we use the Vaughan-style infinite product; small corrections immaterial here).

    Hardy-Littlewood-Vinogradov: R_3(N) ~ S_3(N) * N^2 / (2 * log^3 N) for odd N.

    For our purposes we just need the SCALING with N to compare to the
    leading singular value of M_2.
    """
    s = sieve(p_max)
    primes = np.where(s)[0]
    prod1 = 1.0
    for p in primes:
        if p < 2:
            continue
        prod1 *= 1.0 + 1.0 / (p - 1) ** 3
    prod2 = 1.0
    n_factors = []
    m = N
    for p in primes:
        if p * p > m:
            break
        if m % p == 0:
            n_factors.append(int(p))
            while m % p == 0:
                m //= p
    if m > 1:
        n_factors.append(int(m))
    for p in n_factors:
        prod2 *= 1.0 - 1.0 / ((p - 1) ** 2 - (p - 1) + 1)
    return 0.5 * prod1 * prod2


def hl_goldbach_singular_series(M: int, p_max: int = 10_000) -> float:
    """Hardy-Littlewood Goldbach 2-prime: r_2(M) ~ S_2(M) * M / log^2 M (even M).

    S_2(M) = 2 C_2 * prod_{p|M, p odd} (p-1)/(p-2).
    """
    if M % 2:
        return 0.0
    C2 = hl_twin_constant(p_max)
    out = 2.0 * C2
    s = sieve(p_max)
    primes = np.where(s)[0]
    m = M
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
    if m > 2 and m & 1:
        out *= (m - 1) / (m - 2)
    return out


# ---------------------------------------------------------------- tensor build


def build_constraint_matrix(N: int, indicator: np.ndarray) -> sparse.csr_matrix:
    """Return M_2 as defined above for a given indicator (boolean array length N+1)."""
    rows = []
    cols = []
    vals = []
    lo = 3
    hi = N - 6
    for i in range(lo, hi + 1):
        if not indicator[i]:
            continue
        # j range: j in [lo, hi]; k = N - i - j must be in [lo, hi].
        # so j in [max(lo, N-i-hi), min(hi, N-i-lo)]
        j_min = max(lo, N - i - hi)
        j_max = min(hi, N - i - lo)
        if j_min > j_max:
            continue
        for j in range(j_min, j_max + 1):
            if not indicator[j]:
                continue
            k = N - i - j
            if 0 <= k <= len(indicator) - 1 and indicator[k]:
                rows.append(i - lo)
                cols.append(j - lo)
                vals.append(1.0)
    n = hi - lo + 1
    M = sparse.coo_matrix((vals, (rows, cols)), shape=(n, n)).tocsr()
    return M


def build_bernoulli_baseline(
    N: int, density: float, rng: np.random.Generator
) -> sparse.csr_matrix:
    """Random Bernoulli indicator with given density on [0..N], then build M_2."""
    indicator = rng.random(N + 1) < density
    indicator[:2] = False
    return build_constraint_matrix(N, indicator)


def build_liouville_pos_baseline(N: int) -> sparse.csr_matrix:
    """Liouville-positive indicator: lambda(n)=+1 (squarefull-even-omega).

    Density ~ 1/2; gives a parity-matched control with no primality
    structure but the same 'density-1/2' scaling as random.
    """
    L = liouville(N)
    indicator = (L > 0).astype(bool)
    indicator[:2] = False
    return build_constraint_matrix(N, indicator)


# ---------------------------------------------------------------- diagonals


def mode1_diagonal_spectrum(
    N: int, indicator: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Compute the diagonal of M_1 M_1^T: r_2-like row counts.

    For row i, this is #{(j,k) : i+j+k=N, j,k in valid range, all three in indicator support}.
    Returns (i_values, diag_values).
    """
    lo = 3
    hi = N - 6
    n = hi - lo + 1
    i_vals = np.arange(lo, hi + 1, dtype=np.int64)
    diag = np.zeros(n, dtype=np.int64)
    # For each i, count primes j in [max(lo, N-i-hi), min(hi, N-i-lo)] with N-i-j prime.
    for idx, i in enumerate(i_vals):
        if not indicator[i]:
            continue
        j_min = max(lo, N - i - hi)
        j_max = min(hi, N - i - lo)
        if j_min > j_max:
            continue
        c = 0
        for j in range(j_min, j_max + 1):
            if indicator[j] and indicator[N - i - j]:
                c += 1
        diag[idx] = c
    return i_vals, diag


# ---------------------------------------------------------------- analysis


@dataclass
class SpectrumStats:
    label: str
    N: int
    aspect: float
    n_rows: int
    nnz: int
    density: float
    sigma_max: float
    sigma_top10: list
    fro_norm_sq: float
    rank_estimate: int
    bulk_match_mp: dict
    edge_match_tw: dict
    rank1_residual: dict


def marchenko_pastur_pdf(x: np.ndarray, gamma: float, sigma2: float) -> np.ndarray:
    """MP density for sample covariance with aspect gamma = p/n, variance sigma^2.

    Support [(1-sqrt(gamma))^2 sigma^2, (1+sqrt(gamma))^2 sigma^2].
    """
    a = (1 - math.sqrt(gamma)) ** 2 * sigma2
    b = (1 + math.sqrt(gamma)) ** 2 * sigma2
    out = np.zeros_like(x)
    mask = (x > a) & (x < b)
    out[mask] = (
        np.sqrt((b - x[mask]) * (x[mask] - a))
        / (2 * math.pi * gamma * x[mask] * sigma2)
    )
    return out


def tracy_widom_beta2_mean_var() -> tuple[float, float]:
    """TW_2 mean and variance (numerical)."""
    # Reference values: mean ~ -1.7711, variance ~ 0.8132
    return -1.7711, 0.8132


def fit_mp_to_eigenvalues(eigs: np.ndarray, gamma: float) -> dict:
    """Fit MP to bulk: estimate sigma^2 via mean(eigs) (since mean of MP = sigma^2).
    Then compute Wasserstein-1 (avg abs diff) between empirical CDF and MP CDF on a grid.
    """
    eigs = np.sort(eigs[eigs > 0])
    if len(eigs) < 50:
        return {"ok": False, "reason": "too few eigenvalues"}
    sigma2 = float(np.mean(eigs))
    if sigma2 <= 0:
        return {"ok": False, "reason": "non-positive sigma2"}
    a = (1 - math.sqrt(gamma)) ** 2 * sigma2
    b = (1 + math.sqrt(gamma)) ** 2 * sigma2
    grid = np.linspace(max(a, eigs.min()), min(b, eigs.max()), 1000)
    pdf_mp = marchenko_pastur_pdf(grid, gamma, sigma2)
    # Empirical histogram
    hist, edges = np.histogram(eigs, bins=50, density=True)
    centres = 0.5 * (edges[:-1] + edges[1:])
    pdf_mp_at = marchenko_pastur_pdf(centres, gamma, sigma2)
    # KS-style L1 distance
    cdf_mp = np.cumsum(pdf_mp) * (grid[1] - grid[0])
    cdf_mp /= cdf_mp[-1] if cdf_mp[-1] > 0 else 1.0
    eigs_in_supp = eigs[(eigs >= a) & (eigs <= b)]
    if len(eigs_in_supp) < 10:
        return {
            "ok": True,
            "sigma2": sigma2,
            "frac_in_mp_supp": len(eigs_in_supp) / len(eigs),
            "L1_dist": None,
            "ks_dist": None,
            "max_eig": float(eigs.max()),
            "mp_upper_edge": b,
            "comment": "too few eigenvalues in MP support",
        }
    eigs_in_supp = np.sort(eigs_in_supp)
    emp_cdf = np.arange(1, len(eigs_in_supp) + 1) / len(eigs_in_supp)
    mp_cdf_at_eigs = np.interp(eigs_in_supp, grid, cdf_mp)
    ks = float(np.max(np.abs(emp_cdf - mp_cdf_at_eigs)))
    L1 = float(np.mean(np.abs(emp_cdf - mp_cdf_at_eigs)))
    return {
        "ok": True,
        "sigma2": sigma2,
        "frac_in_mp_supp": len(eigs_in_supp) / len(eigs),
        "L1_dist": L1,
        "ks_dist": ks,
        "max_eig": float(eigs.max()),
        "mp_upper_edge": float(b),
        "edge_overshoot_ratio": float(eigs.max() / b),
    }


def estimate_tw_zscore(top_eig: float, gamma: float, sigma2: float, n: int) -> dict:
    """Apply MP+TW edge rescaling.

    For a sample covariance W = X X^T / m with X iid mean-0 variance-sigma^2,
    the largest eigenvalue lambda_max satisfies
       (lambda_max - mu_n) / sigma_n -> TW_2,
    where with N = n, m = n / gamma:
       mu_n = sigma^2 (1 + sqrt(gamma))^2,
       sigma_n = sigma^2 (1 + sqrt(gamma)) * (1/sqrt(gamma) + 1)^{1/3} * n^{-2/3}.
    """
    sg = math.sqrt(gamma)
    mu_n = sigma2 * (1 + sg) ** 2
    sigma_n = (
        sigma2
        * (1 + sg)
        * (1 / sg + 1) ** (1 / 3)
        * n ** (-2 / 3)
    )
    z = (top_eig - mu_n) / sigma_n
    tw_mean, tw_var = tracy_widom_beta2_mean_var()
    return {
        "mu_n": mu_n,
        "sigma_n": sigma_n,
        "z_score": z,
        "tw_mean": tw_mean,
        "tw_std": math.sqrt(tw_var),
        "delta_from_tw_mean_in_tw_std": (z - tw_mean) / math.sqrt(tw_var),
    }


def rank1_outer_residual(M: sparse.csr_matrix, k: int = 1) -> dict:
    """Compute leading-k SVD outer product approximation; report Frobenius residual ratio."""
    if min(M.shape) < k + 1 or M.nnz == 0:
        return {"ok": False, "reason": "matrix too small or empty"}
    fro = float(np.sqrt((M.multiply(M)).sum()))
    if fro == 0:
        return {"ok": False, "reason": "zero matrix"}
    try:
        U, S, Vt = svds(M.astype(float), k=k, which="LM")
    except Exception as e:
        return {"ok": False, "reason": str(e)}
    captured = float(np.sqrt(np.sum(S ** 2)))
    return {
        "ok": True,
        "fro": fro,
        "rank{}_captured".format(k): captured,
        "rank{}_residual_ratio".format(k): math.sqrt(
            max(0.0, fro ** 2 - sum(s ** 2 for s in S))
        ) / fro,
        "top_singular_values": S.tolist(),
    }


def analyse_matrix(label: str, N: int, M: sparse.csr_matrix, k_top: int = 30) -> SpectrumStats:
    n = M.shape[0]
    nnz = int(M.nnz)
    density = nnz / max(1, n * n)
    aspect = 1.0  # square
    fro2 = float((M.multiply(M)).sum())
    if nnz == 0 or fro2 == 0:
        return SpectrumStats(
            label=label,
            N=N,
            aspect=aspect,
            n_rows=n,
            nnz=nnz,
            density=density,
            sigma_max=0.0,
            sigma_top10=[0.0] * 10,
            fro_norm_sq=0.0,
            rank_estimate=0,
            bulk_match_mp={"ok": False, "reason": "empty"},
            edge_match_tw={"ok": False, "reason": "empty"},
            rank1_residual={"ok": False, "reason": "empty"},
        )
    # Top singular values
    k_top = min(k_top, n - 1)
    try:
        U, S, Vt = svds(M.astype(float), k=k_top, which="LM")
        S = np.sort(S)[::-1]
    except Exception:
        S = np.array([0.0])
    sigma_max = float(S[0])
    sigma_top10 = S[:10].tolist() if len(S) >= 10 else (S.tolist() + [0.0] * (10 - len(S)))

    # Eigenvalues of M M^T (sample covariance) are sigma_i^2.
    # Estimate bulk via the histogram of sigma_i^2 from a SHORTER svds + spectral truncation.
    # Strategy: sample bulk by computing eigenvalues of M M^T directly if n small enough;
    # otherwise estimate via Lanczos batch SVDs.
    if n <= 4500:
        Mdense = M.toarray()
        MMt = Mdense @ Mdense.T
        eigs = np.linalg.eigvalsh(MMt)
        eigs = eigs[eigs > 1e-10]
        rank_est = int(np.sum(eigs > 1e-10))
    else:
        # Approximate via Lanczos: compute many singular values via repeated svds
        # (bounded by n_top=200).
        try:
            U, S2, Vt = svds(M.astype(float), k=min(200, n - 2), which="LM")
            eigs = S2 ** 2
            rank_est = int(np.sum(eigs > 1e-8))
        except Exception:
            eigs = np.array([])
            rank_est = 0

    if len(eigs) > 50:
        gamma = aspect
        bulk = fit_mp_to_eigenvalues(eigs, gamma)
        if bulk.get("ok"):
            sigma2 = bulk["sigma2"]
            edge = estimate_tw_zscore(float(eigs.max()), gamma, sigma2, n)
        else:
            edge = {"ok": False, "reason": "no MP fit"}
    else:
        bulk = {"ok": False, "reason": "too few eigenvalues"}
        edge = {"ok": False, "reason": "too few eigenvalues"}

    rank1 = rank1_outer_residual(M, k=1)

    return SpectrumStats(
        label=label,
        N=N,
        aspect=aspect,
        n_rows=n,
        nnz=nnz,
        density=density,
        sigma_max=sigma_max,
        sigma_top10=sigma_top10,
        fro_norm_sq=fro2,
        rank_estimate=int(rank_est),
        bulk_match_mp=bulk,
        edge_match_tw=edge,
        rank1_residual=rank1,
    )


# ---------------------------------------------------------------- main


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--Ns", type=str, default="1024,2048,4096")
    parser.add_argument(
        "--out", type=str, default="experiments/constructions/d41_chip_3tensor"
    )
    parser.add_argument("--seed", type=int, default=20260430)
    args = parser.parse_args()

    Ns = [int(x) for x in args.Ns.split(",")]
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    results: dict = {"by_N": {}, "meta": {"Ns": Ns, "seed": args.seed}}

    for N in Ns:
        t0 = time.time()
        s_pi = sieve(N)
        n_primes = int(s_pi.sum())
        prime_density = n_primes / (N + 1)
        print(f"\n=== N = {N} (pi(N)={n_primes}, density={prime_density:.4f}) ===")

        # Build matrices
        M_chiP = build_constraint_matrix(N, s_pi)
        M_bern = build_bernoulli_baseline(N, prime_density, rng)
        M_liou = build_liouville_pos_baseline(N)

        # Mode-1 diagonal spectrum (provable)
        i_vals, diag_chiP = mode1_diagonal_spectrum(N, s_pi)

        # Analyse each
        stats_chiP = analyse_matrix("chi_P", N, M_chiP)
        stats_bern = analyse_matrix("bernoulli_density_matched", N, M_bern)
        stats_liou = analyse_matrix("liouville_positive", N, M_liou)

        # HL ternary leading-singular-value scaling
        hl_R3 = hl_ternary_singular_series(N) * (N ** 2) / (math.log(N) ** 3 if N > 2 else 1)
        sigma_max_pred_hl = math.sqrt(hl_R3 / max(1, stats_chiP.n_rows))
        # rough predicted leading sigma ~ sqrt(R_3(N) / n_rows) by mass-spread heuristic
        # if M is rank-1 with all 1/n entries equal density.

        elapsed = time.time() - t0
        results["by_N"][N] = {
            "n_primes": n_primes,
            "prime_density": prime_density,
            "wall_seconds": elapsed,
            "chi_P": asdict(stats_chiP),
            "bernoulli_baseline": asdict(stats_bern),
            "liouville_baseline": asdict(stats_liou),
            "mode1_diagonal_chiP": {
                "n_nonzero_rows": int(np.sum(diag_chiP > 0)),
                "max_diag": int(diag_chiP.max()),
                "mean_diag_over_primes": float(
                    np.mean(diag_chiP[diag_chiP > 0]) if (diag_chiP > 0).any() else 0
                ),
                "std_diag_over_primes": float(
                    np.std(diag_chiP[diag_chiP > 0]) if (diag_chiP > 0).any() else 0
                ),
                "histogram_bins_0_10_20_30_40_50": [
                    int(np.sum(diag_chiP == k)) for k in range(0, 51, 5)
                ],
            },
            "hl_ternary_predicted_R3": hl_R3,
            "sigma_max_pred_from_hl_heuristic": sigma_max_pred_hl,
        }
        # Print short summary
        print(
            f"  chi_P:    rank~{stats_chiP.rank_estimate}, "
            f"sigma_max={stats_chiP.sigma_max:.3f}, "
            f"density={stats_chiP.density:.4f}, "
            f"top SV (top 5)={[f'{x:.2f}' for x in stats_chiP.sigma_top10[:5]]}"
        )
        print(
            f"  bernoulli: rank~{stats_bern.rank_estimate}, "
            f"sigma_max={stats_bern.sigma_max:.3f}, "
            f"density={stats_bern.density:.4f}, "
            f"top SV (top 5)={[f'{x:.2f}' for x in stats_bern.sigma_top10[:5]]}"
        )
        print(
            f"  liouville: rank~{stats_liou.rank_estimate}, "
            f"sigma_max={stats_liou.sigma_max:.3f}, "
            f"density={stats_liou.density:.4f}"
        )
        if stats_chiP.bulk_match_mp.get("ok"):
            print(f"  chi_P MP: {stats_chiP.bulk_match_mp}")
        if stats_chiP.edge_match_tw.get("ok") if isinstance(stats_chiP.edge_match_tw, dict) else False:
            print(f"  chi_P TW edge: {stats_chiP.edge_match_tw}")
        if stats_chiP.rank1_residual.get("ok"):
            print(f"  chi_P rank-1 residual: {stats_chiP.rank1_residual}")
        print(
            f"  HL R_3(N) prediction: {hl_R3:.2f}, "
            f"empirical Frobenius^2 (chi_P 3-prime triple count): {stats_chiP.fro_norm_sq:.2f}"
        )
        print(f"  wall: {elapsed:.1f}s")

    out = out_dir / "results.json"
    out.write_text(json.dumps(results, indent=2, default=float))
    print(f"\nWrote {out}")


if __name__ == "__main__":
    main()

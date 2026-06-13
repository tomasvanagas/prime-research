"""
D3 attack: free-probability via Bryc-Dembo-Jiang Toeplitz/Hankel LSD on chi_P.

We embed the prime indicator chi_P[1..2N] as the data sequence underlying
real symmetric Toeplitz T_N (T_{ij} = chi_P(|i-j|+1) for 0<=i,j<N) and
real symmetric Hankel H_N (H_{ij} = chi_P(i+j+1)) matrices.

After centring (chi_P - p_N) and normalising entries to unit variance, the
classical Bryc-Dembo-Jiang 2006 (Ann. Probab. 34, 1-38) theorem says that for
i.i.d. mean-0 variance-1 entries, eigenvalues / sqrt(N) converge to universal
distributions gamma_T (Toeplitz) and gamma_H (Hankel) of unbounded support
with even moments
    m_{2k} = sum over "good" link/non-crossing pair partitions (Eulerian
             number combinatorics; see BDJ 2006 + Hammond-Miller 2005).
The first numerical values are
    Toeplitz:  m_2 = 1, m_4 = 8/3, m_6 ~ 11.67 (BDJ 2006)
    Hankel:    m_2 = 1, m_4 = 2, m_6 ~ 11
For Gaussian (semicircle is wrong here - that's Wigner) reference, we use
Bernoulli-matched-density NULL (i.i.d. Bern(p) - p, divided by sqrt(p(1-p))),
which is the exact 1st-2-moment-matched chi_P null.

Falsifier (pre-registered):
- Mode I (A-grade): chi_P empirical moment m_k differs from matched-Bernoulli
  null by >= 3 sigma at any k in {2, 4, 6, 8} for at least 2 distinct N >=
  1500, in the SAME direction for both T and H.
- Mode E (B-grade orthogonality): chi_P moments match null within +-2 sigma
  across all N tested. -> 15th orthogonal pseudorandomness category.
- Mode INC: variance dominates signal at the scales tested.

Cross-domain ingredient: Bryc-Dembo-Jiang 2006 universal Toeplitz/Hankel LSD,
real symmetric variant; Hammond-Miller 2005 Toeplitz "good" pair partition
moment expansion. Channelled mathematician: Bryc / Voiculescu (operator-
algebra free probability rigour).
"""
from __future__ import annotations

import argparse
import json
import math
import os
import time
from dataclasses import dataclass, asdict
from typing import Any

import numpy as np


# ---------------------------------------------------------------------------
# Sieve and chi_P sequence
# ---------------------------------------------------------------------------
def sieve_chi_p(M: int) -> np.ndarray:
    """Return chi_P[0..M-1] with chi_P[k] = 1 iff k+1 is prime, else 0.

    So chi_P sequence is indexed from 1 (i.e. chi_P at index 0 represents n=1).
    Convention: 1 is NOT prime.
    """
    sieve = np.ones(M + 1, dtype=bool)
    sieve[0] = sieve[1] = False
    for i in range(2, int(math.isqrt(M)) + 1):
        if sieve[i]:
            sieve[i * i :: i] = False
    return sieve[1 : M + 1].astype(np.float64)


# ---------------------------------------------------------------------------
# Matrix builders
# ---------------------------------------------------------------------------
def build_toeplitz(seq: np.ndarray, N: int) -> np.ndarray:
    """Symmetric Toeplitz, T_{ij} = seq[|i-j|], i,j in [0,N).

    Requires len(seq) >= N.
    """
    idx = np.abs(np.arange(N)[:, None] - np.arange(N)[None, :])
    return seq[idx]


def build_hankel(seq: np.ndarray, N: int) -> np.ndarray:
    """Symmetric Hankel, H_{ij} = seq[i+j], i,j in [0,N).

    Requires len(seq) >= 2N - 1.
    """
    idx = np.arange(N)[:, None] + np.arange(N)[None, :]
    return seq[idx]


# ---------------------------------------------------------------------------
# Spectral-moment computation
# ---------------------------------------------------------------------------
def normalised_eigenvalues(
    seq: np.ndarray,
    N: int,
    kind: str,
    p_mean: float,
    p_var: float,
) -> np.ndarray:
    """Centre `seq` by p_mean, scale to unit-variance, build matrix, return
    eigenvalues divided by sqrt(N) (BDJ scaling)."""
    if p_var <= 0:
        raise ValueError(f"non-positive entry variance p_var={p_var}")
    s_centered = (seq - p_mean) / math.sqrt(p_var)
    if kind == "T":
        M = build_toeplitz(s_centered, N)
    elif kind == "H":
        M = build_hankel(s_centered, N)
    else:
        raise ValueError(kind)
    w = np.linalg.eigvalsh(M)
    return w / math.sqrt(N)


def spectral_moments(eigs: np.ndarray, ks: tuple[int, ...]) -> dict[int, float]:
    """Empirical even moments m_k = mean_{lambda} lambda^k."""
    out: dict[int, float] = {}
    for k in ks:
        out[k] = float(np.mean(eigs ** k))
    return out


# ---------------------------------------------------------------------------
# Nulls: (B1) i.i.d. Bernoulli(p) matched-density, (B2) random permutation
# of chi_P (matches multiset of values exactly), same matrix-type
# ---------------------------------------------------------------------------
def bernoulli_null_eigs(
    rng: np.random.Generator,
    N: int,
    kind: str,
    p: float,
    n_seq: int,
) -> np.ndarray:
    """B1: sample i.i.d. Bern(p) - p sequence of length n_seq, build
    Toeplitz/Hankel of size N (after standardising), return eigenvalues/sqrt(N)."""
    raw = rng.binomial(1, p, size=n_seq).astype(np.float64)
    p_var = p * (1.0 - p)
    return normalised_eigenvalues(raw, N, kind, p_mean=p, p_var=p_var)


def permutation_null_eigs(
    rng: np.random.Generator,
    chi_seq: np.ndarray,
    N: int,
    kind: str,
    p: float,
) -> np.ndarray:
    """B2: random permutation of chi_seq; matches MULTISET of values exactly.
    Any deviation from chi_P here is from positional order (= prime arithmetic
    structure), not from value-distribution mismatch."""
    perm = rng.permutation(chi_seq)
    p_var = p * (1.0 - p)
    return normalised_eigenvalues(perm, N, kind, p_mean=p, p_var=p_var)


# ---------------------------------------------------------------------------
# Sweep
# ---------------------------------------------------------------------------
@dataclass
class NullStats:
    mean: dict[int, float]
    std: dict[int, float]
    z: dict[int, float]
    trials: int


@dataclass
class CellResult:
    N: int
    kind: str  # "T" or "H"
    p_N: float
    chi_moments: dict[int, float]
    bern: NullStats
    perm: NullStats
    chi_lambda_max: float
    chi_lambda_min: float
    chi_top5: list[float]
    chi_bot5: list[float]
    chi_bulk_moments: dict[int, float]  # excluding top-1 and bot-1
    seconds: float


def _null_stats(
    null_arrays: dict[int, np.ndarray],
    chi_m: dict[int, float],
    moments: tuple[int, ...],
    n_trials: int,
) -> NullStats:
    mean = {k: float(null_arrays[k].mean()) for k in moments}
    std = {k: float(null_arrays[k].std(ddof=1)) for k in moments}
    z = {
        k: (chi_m[k] - mean[k]) / std[k] if std[k] > 0 else 0.0
        for k in moments
    }
    return NullStats(mean=mean, std=std, z=z, trials=n_trials)


def run_cell(
    chi_seq_full: np.ndarray,
    N: int,
    kind: str,
    p_N: float,
    moments: tuple[int, ...],
    null_trials: int,
    rng: np.random.Generator,
) -> CellResult:
    t0 = time.time()
    n_seq_needed = N if kind == "T" else 2 * N - 1
    chi_seq_used = chi_seq_full[:n_seq_needed]
    chi_eigs = normalised_eigenvalues(
        chi_seq_used, N, kind, p_mean=p_N, p_var=p_N * (1.0 - p_N)
    )
    chi_m = spectral_moments(chi_eigs, moments)
    # Bulk = drop the (single) top and bot eigenvalue (parity major-arc spike
    # is the dominant outlier; for Hankel it's a +/- pair, so drop both ends).
    bulk = chi_eigs[1:-1]
    chi_bulk_m = spectral_moments(bulk, moments)

    bern_arrays = {k: np.empty(null_trials, dtype=np.float64) for k in moments}
    perm_arrays = {k: np.empty(null_trials, dtype=np.float64) for k in moments}
    for t in range(null_trials):
        e_b = bernoulli_null_eigs(rng, N, kind, p_N, n_seq_needed)
        m_b = spectral_moments(e_b, moments)
        e_p = permutation_null_eigs(rng, chi_seq_used, N, kind, p_N)
        m_p = spectral_moments(e_p, moments)
        for k in moments:
            bern_arrays[k][t] = m_b[k]
            perm_arrays[k][t] = m_p[k]

    bern = _null_stats(bern_arrays, chi_m, moments, null_trials)
    perm = _null_stats(perm_arrays, chi_m, moments, null_trials)
    return CellResult(
        N=N,
        kind=kind,
        p_N=p_N,
        chi_moments=chi_m,
        bern=bern,
        perm=perm,
        chi_lambda_max=float(chi_eigs.max()),
        chi_lambda_min=float(chi_eigs.min()),
        chi_top5=[float(x) for x in chi_eigs[-5:]],
        chi_bot5=[float(x) for x in chi_eigs[:5]],
        chi_bulk_moments=chi_bulk_m,
        seconds=time.time() - t0,
    )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--Ns",
        type=int,
        nargs="+",
        default=[500, 1000, 2000, 3500],
        help="matrix sizes",
    )
    ap.add_argument("--moments", type=int, nargs="+", default=[2, 4, 6, 8])
    ap.add_argument("--null-trials", type=int, default=50)
    ap.add_argument("--seed", type=int, default=20260428)
    ap.add_argument("--out-json", type=str, default=None)
    args = ap.parse_args()

    Ns = sorted(set(args.Ns))
    rng = np.random.default_rng(args.seed)
    moments = tuple(args.moments)

    # Pre-sieve enough chi_P
    M_needed = 2 * max(Ns) + 5
    chi_full = sieve_chi_p(M_needed)
    print(f"[sieve] generated chi_P on [1, {M_needed}], "
          f"prime count = {int(chi_full.sum())}")

    cells: list[CellResult] = []
    for N in Ns:
        n_seq = 2 * N - 1
        p_N = float(chi_full[:n_seq].mean())
        for kind in ("T", "H"):
            print(f"[cell] N={N} kind={kind} p_N={p_N:.6f} ... ", flush=True)
            r = run_cell(
                chi_full, N, kind, p_N, moments, args.null_trials, rng
            )
            cells.append(r)
            zs_b = ", ".join(f"m{k}:zB={r.bern.z[k]:+.2f}" for k in moments)
            zs_p = ", ".join(f"m{k}:zP={r.perm.z[k]:+.2f}" for k in moments)
            print(f"   chi:        {r.chi_moments}")
            print(f"   chi bulk:   {r.chi_bulk_moments}")
            print(f"   bern mean:  {r.bern.mean}")
            print(f"   perm mean:  {r.perm.mean}")
            print(f"   {zs_b}")
            print(f"   {zs_p}")
            top1_frac = (r.chi_lambda_max ** 4) / (N * r.chi_moments[4]) if 4 in r.chi_moments else 0.0
            sn_log = math.sqrt(N / math.log(N))
            print(f"   lambda_max = {r.chi_lambda_max:+.3f}, "
                  f"lambda_max / sqrt(N/log N) = {r.chi_lambda_max/sn_log:.3f}")
            print(f"   top-eig contrib to m_4: {top1_frac*100:.1f}%   "
                  f"({r.seconds:.1f}s)")

    # Summarize: max|z| across cells per moment k, separately for bern/perm
    def _max_z(attr: str) -> tuple[dict[int, float], dict[int, str]]:
        mz = {k: 0.0 for k in moments}
        am: dict[int, str] = {}
        for r in cells:
            stats = getattr(r, attr)
            for k in moments:
                if abs(stats.z[k]) > mz[k]:
                    mz[k] = abs(stats.z[k])
                    am[k] = f"N={r.N} kind={r.kind}"
        return mz, am

    max_z_bern, arg_bern = _max_z("bern")
    max_z_perm, arg_perm = _max_z("perm")

    print()
    print("=" * 60)
    print("Summary: max |z| per moment across all cells")
    for k in moments:
        print(f"  m_{k}: bern max|z| = {max_z_bern[k]:7.2f} at "
              f"{arg_bern[k]:<22} | perm max|z| = "
              f"{max_z_perm[k]:7.2f} at {arg_perm[k]}")

    n_tests = len(cells) * len(moments)
    bonf = 3 + math.sqrt(2 * math.log(n_tests))
    print(f"  total cells x moments = {n_tests}; Bonferroni-3sigma is "
          f"~ {bonf:.2f}-sigma raw threshold")

    out = {
        "config": vars(args),
        "cells": [asdict(c) for c in cells],
        "summary": {
            "bern_max_z_per_moment": max_z_bern,
            "bern_arg_max_z_per_moment": arg_bern,
            "perm_max_z_per_moment": max_z_perm,
            "perm_arg_max_z_per_moment": arg_perm,
            "n_cells_x_moments": n_tests,
            "bonferroni_3sigma_threshold": bonf,
        },
    }

    out_json = args.out_json or os.path.join(
        os.path.dirname(__file__),
        "bdj_toeplitz_hankel_chi_p_results.json",
    )
    # Convert int-keyed dicts to str for JSON
    def _stringify(o: Any) -> Any:
        if isinstance(o, dict):
            return {str(k): _stringify(v) for k, v in o.items()}
        if isinstance(o, list):
            return [_stringify(v) for v in o]
        return o

    with open(out_json, "w") as f:
        json.dump(_stringify(out), f, indent=2, default=str)
    print(f"\nWrote {out_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

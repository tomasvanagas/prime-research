"""
ATTACK_VECTORS §G2 — Gowers U^k norms of the Liouville function lambda(n).

Wild swing test of the Green-Tao Mobius / nilsequence orthogonality theorem
(arXiv:0807.1736) at *empirical* finite-N scale. Predictions:

  GT theorem (asymptotic):
    For every k >= 2 and A > 0,  ||lambda||_{U^k[Z/NZ]} = O(log^{-A} N).

  IID Rademacher (centred random +/-1) prediction:
    E ||rand||_{U^2[Z/NZ]}^4 ~ 2/N  (so ||rand||_{U^2} ~ N^{-1/4}).
    E ||rand||_{U^k[Z/NZ]}^{2^k} ~ k!/N^{k-1}  (asymptotically).

Three regimes, three grades:

  (E) lambda gives the SAME or SMALLER U^k norm than matched Rademacher random
      => Mobius-orthogonality empirically confirmed at scale, B-grade.
      Specifically: the ratio U^k(lambda) / U^k(rand) -> 0 or stays <= 1.

  (I) lambda gives a U^k norm LARGER than Rademacher random by a constant
      factor => same as random up to scale; tightens GT to constant rate.

  (NOVEL/A) lambda gives U^k with sustained ratio > 1 to random, OR exhibits
      sub-Rademacher decay rate, OR has a stable energy shift inconsistent
      with GT log^{-A} prediction.

Distinct from S87 (D6, chi_P Gowers): S87 used the indicator 1[lambda(n)=-1]
which has mean 1/2 and is a different object. G2 uses lambda(n) in {-1, 0, +1}
itself (lambda(n) = 0 for n with some prime square factor's contribution
nullified — actually lambda(n) = (-1)^Omega(n) is +/- only, never 0; lambda(0)
is undefined and we set to 0). The centring (mean-zero) is automatic for
lambda since mean -> 0 as N -> infty.

Implementation:
  - Liouville: linear sieve in O(N), works to N = 2^20 in seconds.
  - U^2 via FFT identity: ||f||_{U^2}^4 = (1/N^4) sum_k |fhat(k)|^4.
  - U^3 via recursion: ||f||_{U^3}^8 = (1/N) sum_h ||Delta_h f||_{U^2}^4
        with Delta_h f(x) = f(x) * f(x+h).  O(N^2 log N).
  - 50 Rademacher seeds for tight error bars on ratios.
  - chi_P (centred) reference for cross-comparison with S87.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from typing import Optional

import numpy as np
from sympy import primerange


# ----- Liouville sieve ------------------------------------------------------

def liouville_sieve(N: int) -> np.ndarray:
    """
    Linear sieve for lambda(n) = (-1)^Omega(n) for n in [0, N).
    Returns int8 array; lambda[0] = 0, lambda[1] = 1.
    """
    lam = np.zeros(N, dtype=np.int8)
    if N >= 2:
        lam[1] = 1
    smallest_prime = np.zeros(N, dtype=np.int32)
    primes = []
    for i in range(2, N):
        if smallest_prime[i] == 0:
            smallest_prime[i] = i
            primes.append(i)
            lam[i] = -1
        for p in primes:
            ip = i * p
            if ip >= N:
                break
            smallest_prime[ip] = p
            lam[ip] = -lam[i]
            if i % p == 0:
                break
    return lam


def prime_indicator(N: int) -> np.ndarray:
    chi = np.zeros(N, dtype=np.float64)
    for p in primerange(2, N):
        chi[p] = 1.0
    return chi


# ----- Gowers norms ---------------------------------------------------------

def u2_norm_pow4(f: np.ndarray) -> float:
    """||f||_{U^2}^4 = (1/N^4) sum_k |fhat(k)|^4 on Z/NZ."""
    N = len(f)
    F = np.fft.fft(f)
    return float(np.sum(np.abs(F) ** 4)) / (N ** 4)


def u3_norm_pow8(f: np.ndarray, max_h: Optional[int] = None) -> float:
    """
    ||f||_{U^3}^8 = (1/N) sum_h ||Delta_h f||_{U^2}^4
    where Delta_h f(x) = f(x) * f(x+h).

    O(N^2 log N) full; subsample h to control cost.
    """
    N = len(f)
    if max_h is None or max_h >= N:
        hs = list(range(N))
        denom = N
    else:
        step = max(1, N // max_h)
        hs = list(range(0, N, step))
        denom = len(hs)
    total = 0.0
    for h in hs:
        delta = f * np.roll(f, -h)
        total += u2_norm_pow4(delta)
    return total / denom


# ----- Controls -------------------------------------------------------------

def rademacher(N: int, rng: np.random.Generator) -> np.ndarray:
    """IID +/-1 Rademacher (mean 0, variance 1)."""
    return (2 * rng.integers(0, 2, size=N) - 1).astype(np.float64)


def matched_density_bernoulli(N: int, p: float, rng: np.random.Generator) -> np.ndarray:
    return rng.binomial(1, p, size=N).astype(np.float64)


# ----- Evaluator ------------------------------------------------------------

def evaluate(name: str, f: np.ndarray, do_u3: bool, max_h: Optional[int]) -> dict:
    N = len(f)
    mean = float(f.mean())
    f_centred = f - mean
    t0 = time.time()
    u2_4 = u2_norm_pow4(f_centred)
    t1 = time.time()
    u2 = u2_4 ** 0.25
    out = {
        "name": name,
        "N": N,
        "mean": mean,
        "u2_pow4_centred": u2_4,
        "u2_centred": u2,
        "u2_time_s": t1 - t0,
    }
    if do_u3:
        t2 = time.time()
        u3_8 = u3_norm_pow8(f_centred, max_h=max_h)
        t3 = time.time()
        u3 = u3_8 ** 0.125
        out["u3_pow8_centred"] = u3_8
        out["u3_centred"] = u3
        out["u3_time_s"] = t3 - t2
        out["u3_max_h"] = max_h
    return out


# ----- Per-N driver ---------------------------------------------------------

def run_N(N: int, do_u3: bool, n_seeds: int, max_h_u3: Optional[int],
          rng: np.random.Generator, include_chi: bool) -> dict:
    print(f"\n=== N = {N} (do_u3={do_u3}, max_h_u3={max_h_u3}, n_seeds={n_seeds}) ===")

    # Liouville  (lambda in {-1, 0, +1}; entries 0 only at n=0)
    lam = liouville_sieve(N).astype(np.float64)
    rec_lam = evaluate("lambda", lam, do_u3, max_h_u3)
    print(f"  lambda           u2={rec_lam['u2_centred']:.4e}  "
          f"u3={(rec_lam.get('u3_centred', float('nan'))):.4e}  "
          f"mean={rec_lam['mean']:.6f}")

    # mu(n) = lambda(n) * 1[n squarefree]  -- closely related to Liouville,
    # cleaner Mobius object, uses same sieve.
    # We do NOT include mu by default to keep cost bounded.

    # Rademacher controls
    rec_rad = []
    for s in range(n_seeds):
        r = rademacher(N, rng)
        rec_rad.append(evaluate(f"Rademacher_seed{s}", r, do_u3, max_h_u3))

    rad_u2_4 = np.array([r["u2_pow4_centred"] for r in rec_rad])
    print(f"  Rademacher x{n_seeds}  u2_pow4 mean={rad_u2_4.mean():.4e} "
          f"(theory 2/N = {2.0/N:.4e})  std={rad_u2_4.std(ddof=1):.4e}")

    rad_u3_8 = None
    if do_u3:
        rad_u3_8 = np.array([r["u3_pow8_centred"] for r in rec_rad])
        print(f"                  u3_pow8 mean={rad_u3_8.mean():.4e}  "
              f"std={rad_u3_8.std(ddof=1):.4e}")

    # Z-scores: lambda vs Rademacher pool
    z_u2_pow4 = (rec_lam["u2_pow4_centred"] - rad_u2_4.mean()) / max(rad_u2_4.std(ddof=1), 1e-300)
    z_u2 = (rec_lam["u2_centred"] - np.array([r["u2_centred"] for r in rec_rad]).mean()) / \
           max(np.array([r["u2_centred"] for r in rec_rad]).std(ddof=1), 1e-300)
    print(f"  ZSCORE  lambda vs Rademacher  U^2_pow4: {z_u2_pow4:+.2f}sigma  U^2: {z_u2:+.2f}sigma")
    z_u3_pow8 = None
    z_u3 = None
    if do_u3:
        z_u3_pow8 = (rec_lam["u3_pow8_centred"] - rad_u3_8.mean()) / max(rad_u3_8.std(ddof=1), 1e-300)
        z_u3 = (rec_lam["u3_centred"] - np.array([r["u3_centred"] for r in rec_rad]).mean()) / \
               max(np.array([r["u3_centred"] for r in rec_rad]).std(ddof=1), 1e-300)
        print(f"  ZSCORE  lambda vs Rademacher  U^3_pow8: {z_u3_pow8:+.2f}sigma  U^3: {z_u3:+.2f}sigma")

    out = {
        "N": N,
        "do_u3": do_u3,
        "n_seeds": n_seeds,
        "max_h_u3": max_h_u3,
        "lambda": rec_lam,
        "rademacher": rec_rad,
        "z_lambda_u2_pow4_vs_rad": z_u2_pow4,
        "z_lambda_u2_vs_rad": z_u2,
        "z_lambda_u3_pow8_vs_rad": z_u3_pow8,
        "z_lambda_u3_vs_rad": z_u3,
    }

    # chi_P (centred) at the same N for cross-comparison with S87 results.
    if include_chi:
        chi = prime_indicator(N)
        rec_chi = evaluate("chi_P", chi, do_u3, max_h_u3)
        # Matched-density Bernoulli control
        rec_bern = []
        for s in range(min(n_seeds, 10)):
            b = matched_density_bernoulli(N, rec_chi["mean"], rng)
            rec_bern.append(evaluate(f"BernoulliMatched_seed{s}", b, do_u3, max_h_u3))
        bern_u2_4 = np.array([r["u2_pow4_centred"] for r in rec_bern])
        z_chi_u2_pow4 = (rec_chi["u2_pow4_centred"] - bern_u2_4.mean()) / max(bern_u2_4.std(ddof=1), 1e-300)
        print(f"  chi_P            u2={rec_chi['u2_centred']:.4e}  "
              f"vs Bernoulli u2_pow4 mean={bern_u2_4.mean():.4e}  z={z_chi_u2_pow4:+.2f}sigma")
        out["chi_P"] = rec_chi
        out["bernoulli_chi"] = rec_bern
        out["z_chi_u2_pow4_vs_bern"] = z_chi_u2_pow4

    return out


def fit_decay_rate(N_list, u_pow_list):
    """
    Fit log(U) vs log(N) and log(U) vs log(log N).
    Return slopes plus the GT-side log^{-A} alternative fit metric.
    """
    N = np.array(N_list, dtype=np.float64)
    U = np.array(u_pow_list, dtype=np.float64)
    if np.any(U <= 0) or len(N) < 2:
        return {}
    logN = np.log(N)
    logU = np.log(U)
    A_loglog = np.vstack([logN, np.ones_like(logN)]).T
    slope_logN, intercept_logN = np.linalg.lstsq(A_loglog, logU, rcond=None)[0]

    log_logN = np.log(np.log(N))
    A_loglog_log = np.vstack([log_logN, np.ones_like(log_logN)]).T
    slope_loglogN, intercept_loglogN = np.linalg.lstsq(A_loglog_log, logU, rcond=None)[0]
    return {
        "logU_vs_logN_slope": float(slope_logN),
        "logU_vs_logN_intercept": float(intercept_logN),
        "logU_vs_loglogN_slope": float(slope_loglogN),
        "logU_vs_loglogN_intercept": float(intercept_loglogN),
        "N_list": N_list,
        "U_list": [float(u) for u in U],
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--N-list", type=int, nargs="+",
                    default=[4096, 16384, 65536, 262144])
    ap.add_argument("--N-list-u3", type=int, nargs="+",
                    default=[1024, 2048, 4096])
    ap.add_argument("--n-seeds", type=int, default=30)
    ap.add_argument("--seed", type=int, default=20260428)
    ap.add_argument("--max-h-u3", type=int, default=None,
                    help="Subsample h in U^3 (None = all h, slower).")
    ap.add_argument("--no-chi", action="store_true",
                    help="Skip chi_P comparison.")
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    out = {"config": vars(args), "by_N": {}}

    # Per-N runs.
    for N in args.N_list:
        do_u3 = N in args.N_list_u3
        # chi_P expensive at large N; skip beyond 65536.
        include_chi = (not args.no_chi) and (N <= 65536)
        res = run_N(
            N=N,
            do_u3=do_u3,
            n_seeds=args.n_seeds,
            max_h_u3=args.max_h_u3,
            rng=rng,
            include_chi=include_chi,
        )
        out["by_N"][N] = res
        with open(args.out, "w") as fh:
            json.dump(out, fh, indent=2, default=str)

    # Decay rate fits.
    summary = {}
    Ns = sorted(out["by_N"].keys())
    lam_u2_pow4 = [out["by_N"][N]["lambda"]["u2_pow4_centred"] for N in Ns]
    lam_u2 = [out["by_N"][N]["lambda"]["u2_centred"] for N in Ns]
    summary["lambda_u2_pow4_decay"] = fit_decay_rate(Ns, lam_u2_pow4)
    summary["lambda_u2_decay"] = fit_decay_rate(Ns, lam_u2)
    Ns_u3 = [N for N in Ns if "u3_centred" in out["by_N"][N]["lambda"]]
    if Ns_u3:
        lam_u3_pow8 = [out["by_N"][N]["lambda"]["u3_pow8_centred"] for N in Ns_u3]
        summary["lambda_u3_pow8_decay"] = fit_decay_rate(Ns_u3, lam_u3_pow8)
    out["summary"] = summary
    with open(args.out, "w") as fh:
        json.dump(out, fh, indent=2, default=str)

    print("\n--- summary ---")
    print(f"  lambda U^2_pow4 by N: {[f'{u:.3e}' for u in lam_u2_pow4]}")
    print(f"  Rademacher 2/N by N:  {[f'{2.0/N:.3e}' for N in Ns]}")
    if "logU_vs_logN_slope" in summary["lambda_u2_pow4_decay"]:
        print(f"  lambda U^2_pow4: log(U) = "
              f"{summary['lambda_u2_pow4_decay']['logU_vs_logN_slope']:.3f}*log(N) + "
              f"{summary['lambda_u2_pow4_decay']['logU_vs_logN_intercept']:.3f}")
        print(f"  Rademacher pred: log(U) = -1.000*log(N) + log(2)")
    print(f"\nWrote {args.out}")


if __name__ == "__main__":
    main()

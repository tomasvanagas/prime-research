"""
G2 extension. Adds:
  (a) Mobius mu(n) Gowers U^2 / U^3 norms (squarefree-restricted lambda).
  (b) Subsampled U^3 at larger N (16384, 65536) to test U^3 decay rate.
  (c) Inverse-theorem stress test: search for the largest |fhat(k)|^2
      among all k in [1, N/2] for lambda and mu, compare to predicted
      Rademacher max-Fourier-coefficient law.
"""

from __future__ import annotations

import argparse
import json
import time
from typing import Optional

import numpy as np

from liouville_gowers_uk import (
    liouville_sieve,
    u2_norm_pow4,
    u3_norm_pow8,
    rademacher,
)


def mobius_sieve(N: int) -> np.ndarray:
    """
    mu(n) = lambda(n) for squarefree n, else 0.
    Use linear-sieve smallest-prime to detect squarefree (n / spf is divisible by spf <=> not squarefree).
    """
    spf = np.zeros(N, dtype=np.int32)
    primes = []
    for i in range(2, N):
        if spf[i] == 0:
            spf[i] = i
            primes.append(i)
        for p in primes:
            ip = i * p
            if ip >= N:
                break
            spf[ip] = p
            if i % p == 0:
                break
    mu = np.zeros(N, dtype=np.int8)
    if N >= 2:
        mu[1] = 1
    for n in range(2, N):
        # squarefree iff n / spf[n] is not divisible by spf[n]
        m = n // spf[n]
        if m % spf[n] == 0:
            mu[n] = 0
            continue
        # mu(n) = (-1)^Omega(n) for squarefree n
        # Recurse: mu(n) = -mu(n/spf[n]) if (n/spf[n]) % spf[n] != 0
        mu[n] = -mu[m]
    return mu


def max_fourier_freq(f: np.ndarray, top_k: int = 5):
    """
    Return the top_k frequencies (excluding k=0) by |fhat(k)|, with values.
    Useful for inverse-theorem-style detection: any single Fourier mode
    with |fhat(k)|^2 = Omega(N) is a 1-step nilsequence correlation.
    """
    N = len(f)
    F = np.fft.fft(f - f.mean())
    mag = np.abs(F[1:N // 2 + 1])
    idxs = np.argsort(mag)[::-1][:top_k]
    return [(int(i + 1), float(mag[i]), float(mag[i] ** 2)) for i in idxs]


def evaluate_subsampled(name: str, f: np.ndarray, do_u3: bool,
                        max_h: Optional[int]) -> dict:
    N = len(f)
    mean = float(f.mean())
    f_c = f - mean
    t0 = time.time()
    u2_4 = u2_norm_pow4(f_c)
    t1 = time.time()
    out = {
        "name": name,
        "N": N,
        "mean": mean,
        "u2_pow4_centred": u2_4,
        "u2_centred": u2_4 ** 0.25,
        "u2_time_s": t1 - t0,
    }
    if do_u3:
        t2 = time.time()
        u3_8 = u3_norm_pow8(f_c, max_h=max_h)
        t3 = time.time()
        out["u3_pow8_centred"] = u3_8
        out["u3_centred"] = u3_8 ** 0.125
        out["u3_max_h"] = max_h
        out["u3_time_s"] = t3 - t2
    out["top5_fourier_modes"] = max_fourier_freq(f, top_k=5)
    out["fourier_max_amp_pow2"] = out["top5_fourier_modes"][0][2]
    return out


def run(N_list_u2_only, N_list_u3, n_seeds_u2, n_seeds_u3, max_h_u3, seed):
    rng = np.random.default_rng(seed)
    results = {}

    # Combine lists, run lambda and mu first; then Rademacher controls.
    all_N = sorted(set(N_list_u2_only) | set(N_list_u3))

    for N in all_N:
        do_u3 = N in N_list_u3
        # Cap max_h per N: full at small N, subsample at large N.
        if do_u3:
            cap = max_h_u3 if max_h_u3 else None
            if cap is not None and cap >= N:
                cap = None
        else:
            cap = None

        print(f"\n=== N = {N}  (do_u3={do_u3}, max_h={cap}) ===")
        lam = liouville_sieve(N).astype(np.float64)
        mu = mobius_sieve(N).astype(np.float64)
        rec_lam = evaluate_subsampled("lambda", lam, do_u3, cap)
        rec_mu = evaluate_subsampled("mu", mu, do_u3, cap)
        print(f"  lambda  u2={rec_lam['u2_centred']:.4e}  "
              f"u3={(rec_lam.get('u3_centred', 0)):.4e}  mean={rec_lam['mean']:.4e}  "
              f"max|fhat|^2={rec_lam['fourier_max_amp_pow2']:.2f}  N={N}")
        print(f"  mu      u2={rec_mu['u2_centred']:.4e}  "
              f"u3={(rec_mu.get('u3_centred', 0)):.4e}  mean={rec_mu['mean']:.4e}  "
              f"max|fhat|^2={rec_mu['fourier_max_amp_pow2']:.2f}")

        # Rademacher controls
        n_s = n_seeds_u3 if do_u3 else n_seeds_u2
        rec_rad = []
        for s in range(n_s):
            r = rademacher(N, rng)
            rec_rad.append(evaluate_subsampled(f"Rad_{s}", r, do_u3, cap))
        rad_u2_4 = np.array([r["u2_pow4_centred"] for r in rec_rad])
        rad_max = np.array([r["fourier_max_amp_pow2"] for r in rec_rad])
        print(f"  Rad x{n_s}  u2_pow4 mean={rad_u2_4.mean():.4e} std={rad_u2_4.std(ddof=1):.4e}  "
              f"max|fhat|^2 mean={rad_max.mean():.2f} std={rad_max.std(ddof=1):.2f}")
        # Z-scores
        z_lam = (rec_lam["u2_pow4_centred"] - rad_u2_4.mean()) / max(rad_u2_4.std(ddof=1), 1e-300)
        z_mu = (rec_mu["u2_pow4_centred"] - rad_u2_4.mean()) / max(rad_u2_4.std(ddof=1), 1e-300)
        z_lam_max = (rec_lam["fourier_max_amp_pow2"] - rad_max.mean()) / max(rad_max.std(ddof=1), 1e-300)
        z_mu_max = (rec_mu["fourier_max_amp_pow2"] - rad_max.mean()) / max(rad_max.std(ddof=1), 1e-300)
        print(f"  Z(lambda U^2_pow4 vs Rad) = {z_lam:+.2f}sigma   "
              f"Z(mu U^2_pow4 vs Rad) = {z_mu:+.2f}sigma")
        print(f"  Z(lambda max|fhat|^2 vs Rad) = {z_lam_max:+.2f}sigma   "
              f"Z(mu max|fhat|^2 vs Rad) = {z_mu_max:+.2f}sigma")
        if do_u3:
            rad_u3_8 = np.array([r["u3_pow8_centred"] for r in rec_rad])
            print(f"  Rad x{n_s}  u3_pow8 mean={rad_u3_8.mean():.4e} std={rad_u3_8.std(ddof=1):.4e}")
            z_lam_u3 = (rec_lam["u3_pow8_centred"] - rad_u3_8.mean()) / max(rad_u3_8.std(ddof=1), 1e-300)
            z_mu_u3 = (rec_mu["u3_pow8_centred"] - rad_u3_8.mean()) / max(rad_u3_8.std(ddof=1), 1e-300)
            print(f"  Z(lambda U^3_pow8 vs Rad) = {z_lam_u3:+.2f}sigma   "
                  f"Z(mu U^3_pow8 vs Rad) = {z_mu_u3:+.2f}sigma")
            results.setdefault("z_u3", {})[N] = {"lambda": z_lam_u3, "mu": z_mu_u3}
        results[N] = {
            "lambda": rec_lam,
            "mu": rec_mu,
            "rademacher": rec_rad,
            "z_lambda_u2_pow4": z_lam,
            "z_mu_u2_pow4": z_mu,
            "z_lambda_max_amp": z_lam_max,
            "z_mu_max_amp": z_mu_max,
        }
    return results


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--N-list-u2-only", type=int, nargs="+",
                    default=[16384, 65536, 262144, 1048576])
    ap.add_argument("--N-list-u3", type=int, nargs="+",
                    default=[1024, 2048, 4096, 8192])
    ap.add_argument("--n-seeds-u2", type=int, default=50)
    ap.add_argument("--n-seeds-u3", type=int, default=20)
    ap.add_argument("--max-h-u3", type=int, default=2048)
    ap.add_argument("--seed", type=int, default=20260428)
    args = ap.parse_args()

    out = run(
        N_list_u2_only=args.N_list_u2_only,
        N_list_u3=args.N_list_u3,
        n_seeds_u2=args.n_seeds_u2,
        n_seeds_u3=args.n_seeds_u3,
        max_h_u3=args.max_h_u3,
        seed=args.seed,
    )
    payload = {"config": vars(args), "by_N": out}
    with open(args.out, "w") as fh:
        json.dump(payload, fh, indent=2, default=str)
    print(f"\nWrote {args.out}")


if __name__ == "__main__":
    main()

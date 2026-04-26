"""
§D6 — Gowers U^k uniformity norms of the prime indicator chi_P.

Frontier attack from ATTACK_VECTORS.md §D6 (S85). Goal: compute the Gowers
U^2 and U^3 norms of the bare prime indicator chi_P : Z/NZ -> {0, 1} for
N up to 2^{16} (U^2) and 2^{12} (U^3), and compare against:

  (a) random Bernoulli at matched density p = pi(N)/N;
  (b) the Liouville indicator L(n) := 1[lambda(n) = -1];
  (c) random Bernoulli at density 1/2 (matched to Liouville);
  (d) W-tricked chi_P (Green-Tao W-trick) for several W.

The Green-Tao-Ziegler inverse theorem (arXiv:1009.3998) says: large
||f||_{U^{s+1}} <=> f correlates with an s-step nilsequence. The
project's 35+ pseudorandomness measures contain none from the Gowers
family; this is a fresh angle.

Implementation:
  - U^2 via FFT identity: ||f||_{U^2}^4 = (1/N^4) sum_k |fhat(k)|^4
  - U^3 via recursion: ||f||_{U^3}^8 = (1/N) sum_h ||Delta_h f||_{U^2}^4
        with Delta_h f(x) = f(x) f(x+h)
  - all on Z/NZ with cyclic shift

Falsification (pre-stated):
  H_random:    ||chi_P||_{U^k} matches random Bernoulli at same density,
               within sqrt(N) noise => 36th-37th pseudorandomness measure.
  H_HL:        ||chi_P||_{U^k} deviates from random by a CONSTANT
               multiplicative factor (Hardy-Littlewood singular series for
               the {0,1}^k-cube prime tuple) => no novel content beyond HL.
  H_novel:     ||chi_P||_{U^k} deviates from BOTH random and HL prediction,
               OR W-tricked chi_P fails Gowers-uniformity at some k <= 3.
               This is the A-grade scenario.

Cross-domain refs:
  Green-Tao,  arXiv:math/0606088
  Green-Tao,  arXiv:0807.1736   (Mobius and nilsequences)
  Green-Tao-Ziegler, arXiv:1009.3998 (U^{s+1}[N] inverse theorem)

Cite EDGES.md: this work either *adds* a 36th measure (extends the
pseudorandomness battery indexed by E1.10/E3.13/E7.1) or contradicts it.
"""

from __future__ import annotations
import argparse
import json
import os
import sys
import time
from typing import Tuple

import numpy as np
from sympy import primerange, factorint


# ----- helpers ---------------------------------------------------------------

def prime_indicator(N: int) -> np.ndarray:
    chi = np.zeros(N, dtype=np.float64)
    for p in primerange(2, N):
        chi[p] = 1.0
    return chi


def liouville_indicator(N: int) -> np.ndarray:
    """L(n) = 1 if Omega(n) is odd, else 0. L(0) = 0 by convention."""
    L = np.zeros(N, dtype=np.float64)
    for n in range(2, N):
        Omega = sum(factorint(n).values())
        if Omega % 2 == 1:
            L[n] = 1.0
    return L


def w_trick_chi(N: int, W: int, b: int) -> np.ndarray:
    """
    Returns f[n] = chi_P(W*n + b) for n in [0, N), with b coprime to W.
    The W-trick removes the local density bias coming from small primes p|W:
    on the residue class b mod W (gcd(b,W)=1), primes are NOT depleted by
    p|W, so f has 'cleaner' density phi(W)/W * (log scaling) inside [0, NW).
    """
    M = N * W + b + 10
    chi = np.zeros(M, dtype=np.float64)
    for p in primerange(2, M):
        chi[p] = 1.0
    out = np.zeros(N, dtype=np.float64)
    for n in range(N):
        idx = W * n + b
        if 0 <= idx < M:
            out[n] = chi[idx]
    return out


# ----- Gowers norms ---------------------------------------------------------

def u2_norm_pow4(f: np.ndarray) -> float:
    """||f||_{U^2}^4 = (1/N^4) sum_k |fhat(k)|^4   on Z/NZ."""
    N = len(f)
    F = np.fft.fft(f)
    return float(np.sum(np.abs(F) ** 4)) / (N ** 4)


def u3_norm_pow8(f: np.ndarray, max_h: int | None = None) -> float:
    """
    ||f||_{U^3}^8 = (1/N) sum_{h in Z/N} ||Delta_h f||_{U^2}^4
    where Delta_h f(x) = f(x) * f(x+h).
    O(N^2 log N) via FFT inside U^2.
    If max_h < N, we sample h uniformly from a subsample (faster, larger
    variance).
    """
    N = len(f)
    if max_h is None or max_h >= N:
        hs = range(N)
        denom = N
    else:
        # Use deterministic stride to subsample h cleanly.
        step = max(1, N // max_h)
        hs = range(0, N, step)
        denom = len(hs)
    total = 0.0
    for h in hs:
        delta_h_f = f * np.roll(f, -h)
        total += u2_norm_pow4(delta_h_f)
    return total / denom


# ----- Bernoulli reference -------------------------------------------------

def random_bernoulli(N: int, p: float, rng: np.random.Generator) -> np.ndarray:
    return rng.binomial(1, p, size=N).astype(np.float64)


# ----- main experiment -----------------------------------------------------

def evaluate_function(
    name: str,
    f: np.ndarray,
    do_u3: bool,
    max_h_u3: int | None,
) -> dict:
    N = len(f)
    p = float(f.mean())
    t0 = time.time()
    u2_4 = u2_norm_pow4(f)
    t1 = time.time()
    u2 = u2_4 ** 0.25
    print(f"  {name:34s}  N={N:6d}  p={p:.6f}  ||.||_{{U^2}}^4={u2_4:.6e}  ||.||_{{U^2}}={u2:.6e}  ({t1-t0:.2f}s)")

    out = {
        "name": name,
        "N": N,
        "mean": p,
        "u2_pow4": u2_4,
        "u2": u2,
        "u2_time_s": t1 - t0,
    }
    if do_u3:
        t2 = time.time()
        u3_8 = u3_norm_pow8(f, max_h=max_h_u3)
        t3 = time.time()
        u3 = u3_8 ** 0.125
        print(f"  {name:34s}  N={N:6d}            ||.||_{{U^3}}^8={u3_8:.6e}  ||.||_{{U^3}}={u3:.6e}  ({t3-t2:.2f}s)")
        out["u3_pow8"] = u3_8
        out["u3"] = u3
        out["u3_time_s"] = t3 - t2
        out["u3_max_h"] = max_h_u3
    return out


def derived_ratios(rec_chi: dict, rec_rand_list: list[dict]) -> dict:
    """Mean / std of chi_P / random ratios on U^2 and (if present) U^3."""
    ratios_u2 = [rec_chi["u2"] / r["u2"] for r in rec_rand_list]
    out = {
        "ratio_u2_chi_over_rand_mean": float(np.mean(ratios_u2)),
        "ratio_u2_chi_over_rand_std": float(np.std(ratios_u2, ddof=1)) if len(ratios_u2) > 1 else 0.0,
        "ratio_u2_pow4_chi_over_rand_mean": float(np.mean([rec_chi["u2_pow4"] / r["u2_pow4"] for r in rec_rand_list])),
    }
    if "u3" in rec_chi and "u3" in rec_rand_list[0]:
        ratios_u3 = [rec_chi["u3"] / r["u3"] for r in rec_rand_list]
        out["ratio_u3_chi_over_rand_mean"] = float(np.mean(ratios_u3))
        out["ratio_u3_chi_over_rand_std"] = float(np.std(ratios_u3, ddof=1)) if len(ratios_u3) > 1 else 0.0
        out["ratio_u3_pow8_chi_over_rand_mean"] = float(np.mean([rec_chi["u3_pow8"] / r["u3_pow8"] for r in rec_rand_list]))
    return out


def run(
    N: int,
    do_u3: bool,
    n_rand: int,
    seed: int,
    max_h_u3: int | None,
    do_liouville: bool,
    do_w_trick: bool,
) -> dict:
    print(f"\n=== N = {N}  (do_u3={do_u3}, n_rand={n_rand}, max_h_u3={max_h_u3}) ===")
    rng = np.random.default_rng(seed)

    chi = prime_indicator(N)
    p_chi = float(chi.mean())

    rec_chi = evaluate_function("chi_P", chi, do_u3, max_h_u3)

    # Random Bernoulli at matched density p_chi
    rec_rand = []
    for s in range(n_rand):
        b = random_bernoulli(N, p_chi, rng)
        rec_rand.append(evaluate_function(f"Bernoulli@p_chi seed={s}", b, do_u3, max_h_u3))

    res = {
        "N": N,
        "p_chi": p_chi,
        "do_u3": do_u3,
        "n_rand": n_rand,
        "max_h_u3": max_h_u3,
        "chi_P": rec_chi,
        "bernoulli_chi_density": rec_rand,
        "ratios_chi_vs_bernoulli": derived_ratios(rec_chi, rec_rand),
    }

    if do_liouville:
        print("\n  -- Liouville indicator (1 if Omega(n) odd) --")
        L = liouville_indicator(N)
        rec_L = evaluate_function("Liouville", L, do_u3, max_h_u3)
        # Random matched at p_L
        p_L = float(L.mean())
        rec_rand_L = []
        for s in range(min(n_rand, 5)):
            b = random_bernoulli(N, p_L, rng)
            rec_rand_L.append(evaluate_function(f"Bernoulli@p_L seed={s}", b, do_u3, max_h_u3))
        res["liouville"] = rec_L
        res["bernoulli_liouville_density"] = rec_rand_L
        res["ratios_liouville_vs_bernoulli"] = derived_ratios(rec_L, rec_rand_L)

    if do_w_trick:
        for W in [6, 30, 210]:
            # b = 1 is always coprime to W
            b = 1
            print(f"\n  -- W-trick chi_P with W={W}, b={b} --")
            f = w_trick_chi(N, W, b)
            rec_w = evaluate_function(f"chi_P_W{W}_b{b}", f, do_u3, max_h_u3)
            # Random at matched density
            p_w = float(f.mean())
            rec_rand_w = []
            for s in range(min(n_rand, 5)):
                b_arr = random_bernoulli(N, p_w, rng)
                rec_rand_w.append(evaluate_function(f"Bernoulli@p_w seed={s}", b_arr, do_u3, max_h_u3))
            key = f"w_trick_W{W}_b{b}"
            res[key] = {
                "f": rec_w,
                "bernoulli_matched": rec_rand_w,
                "ratios_w_vs_bernoulli": derived_ratios(rec_w, rec_rand_w),
            }
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True, help="Path to JSON results file.")
    ap.add_argument("--N-list", type=int, nargs="+", default=[1024, 4096, 16384, 65536])
    ap.add_argument("--N-list-u3", type=int, nargs="+", default=[1024, 2048, 4096])
    ap.add_argument("--n-rand", type=int, default=10)
    ap.add_argument("--seed", type=int, default=20260426)
    ap.add_argument("--no-u3", action="store_true")
    ap.add_argument("--max-h-u3", type=int, default=None,
                    help="Subsample h in U^3 (None = all h, slower).")
    ap.add_argument("--liouville-N-cap", type=int, default=8192,
                    help="Skip Liouville above this N (factorint cost).")
    ap.add_argument("--w-trick-N-cap", type=int, default=4096,
                    help="W-trick uses primes up to W*N + b; cap N.")
    ap.add_argument("--no-liouville", action="store_true")
    ap.add_argument("--no-w-trick", action="store_true")
    args = ap.parse_args()

    all_results = {
        "config": vars(args),
        "by_N": {},
    }

    # U^2 only: large N
    for N in args.N_list:
        do_u3 = (not args.no_u3) and (N in args.N_list_u3)
        res = run(
            N=N,
            do_u3=do_u3,
            n_rand=args.n_rand,
            seed=args.seed,
            max_h_u3=args.max_h_u3,
            do_liouville=(not args.no_liouville) and (N <= args.liouville_N_cap),
            do_w_trick=(not args.no_w_trick) and (N <= args.w_trick_N_cap),
        )
        all_results["by_N"][N] = res

        # Snapshot to disk after each N — never lose progress.
        with open(args.out, "w") as f:
            json.dump(all_results, f, indent=2, default=str)

    print(f"\nWrote {args.out}")


if __name__ == "__main__":
    main()

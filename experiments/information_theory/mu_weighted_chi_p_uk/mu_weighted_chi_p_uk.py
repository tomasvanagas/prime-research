"""
§D6.c — Mobius/Liouville-weighted prime-supported functions at U^2.

NOVELTY_CHALLENGES.md §D6.c (S93 follow-up). The literal target — the
"intermediate weighting" mu(n) * chi_P(n) — collapses trivially because
mu(p) = -1 on every prime p, so mu * chi_P = -chi_P pointwise on Z and
||mu * chi_P||_{U^2} = ||chi_P||_{U^2}. We document the collapse and
pivot to the natural broader question:

    Which prime / squarefree / Liouville / Mobius-related functions
    retain HL singular-series structure at U^2, and which collapse
    to Gowers-uniformity?

For each function f we report:
    mean(f), L2_sq(f) = E[f^2],
    ||f||_{U^2}^4 (raw),
    ||f - mean(f)||_{U^2}^4 (centered),
    Q2(f) := ||f||_{U^2}^4 / mean(f)^4   (only meaningful for indicators),
    Q2_norm(f) := N * ||f||_{U^2}^4 / (2 * L2_sq(f)^2)
        (ratio to Bernoulli-baseline; Q2_norm ~ 1 means Gowers-uniform).

Function panel:
    chi_P                = 1[n prime]                          (HL baseline)
    sqfree               = 1[n squarefree] = mu^2(n)           (HL structure?)
    mu_plus              = 1[mu(n) = +1] = squarefree, even omega
    mu_minus             = 1[mu(n) = -1] = squarefree, odd omega
    lam_plus             = 1[lambda(n) = +1] (Liouville positive, S87)
    lam_minus            = 1[lambda(n) = -1] (Liouville negative, S87)
    mu                   = mu(n) (signed)
    lam                  = lambda(n) (signed)
    mu_chi_P             = mu(n) * chi_P(n)   (= -chi_P, the trivial collapse)
    musq_chi_P           = mu^2(n) * chi_P(n) (= chi_P)
    semi_primes          = 1[Omega(n) = 2]   (HL structure on 2-almost-primes)
    chi_P_minus_mean     = chi_P - mean(chi_P) (centered chi_P)

Cite: E2.13 (Gowers norms of chi_P -> S_k, S87/S85). E2.14 (S88, Anderson),
E2.15 (S92, algebraic immunity), S87 (Liouville result).
Cross-domain ref: Sarnak's Mobius randomness conjecture; Green-Tao (2010)
arXiv:math/0606088 inverse-theorem framework.
"""

from __future__ import annotations

import argparse
import json
import math
import time
from pathlib import Path

import numpy as np


# ---------- arithmetic functions on [0, N) ---------------------------------

def smallest_prime_factor(N: int) -> np.ndarray:
    """spf[n] = smallest prime factor of n; spf[0] = spf[1] = 0."""
    spf = np.zeros(N, dtype=np.int64)
    for i in range(2, N):
        if spf[i] == 0:
            spf[i] = i
            for j in range(i * i, N, i):
                if spf[j] == 0:
                    spf[j] = i
    return spf


def factorise(spf: np.ndarray, n: int):
    """Yield (p, k) factor pairs of n using a precomputed spf array."""
    while n > 1:
        p = int(spf[n])
        k = 0
        while n % p == 0:
            n //= p
            k += 1
        yield p, k


def panel_arrays(N: int):
    """Build all panel functions on [0, N) using a single spf scan."""
    spf = smallest_prime_factor(N)

    chi_P = np.zeros(N, dtype=np.float64)
    mu = np.zeros(N, dtype=np.float64)
    lam = np.zeros(N, dtype=np.float64)
    sqfree = np.zeros(N, dtype=np.float64)
    Omega = np.zeros(N, dtype=np.int64)

    chi_P[2:] = (spf[2:] == np.arange(2, N))  # n prime iff spf[n] == n.
    # mu, lam, sqfree, Omega via factorisation.
    for n in range(2, N):
        # walk through prime power decomposition
        m = n
        big_omega = 0
        is_squarefree = 1
        sign_mu = 1
        while m > 1:
            p = int(spf[m])
            k = 0
            while m % p == 0:
                m //= p
                k += 1
            big_omega += k
            if k >= 2:
                is_squarefree = 0
            sign_mu *= -1  # one factor of -1 per *distinct* prime
        Omega[n] = big_omega
        sqfree[n] = float(is_squarefree)
        # mu(n) = (-1)^omega(n) on squarefree, 0 otherwise
        mu[n] = float(sign_mu * is_squarefree) if is_squarefree else 0.0
        # lambda(n) = (-1)^Omega(n)
        lam[n] = (-1.0) ** big_omega

    # mu(1) = 1, lam(1) = 1, sqfree(1) = 1
    if N > 1:
        mu[1] = 1.0
        lam[1] = 1.0
        sqfree[1] = 1.0

    # Indicators
    mu_plus = (mu == 1.0).astype(np.float64)
    mu_minus = (mu == -1.0).astype(np.float64)
    lam_plus = (lam == 1.0).astype(np.float64)
    lam_minus = (lam == -1.0).astype(np.float64)
    semi_primes = (Omega == 2).astype(np.float64)
    mu_chi_P = mu * chi_P
    musq_chi_P = (mu ** 2) * chi_P

    panel = {
        "chi_P":        chi_P,
        "sqfree":       sqfree,
        "mu_plus":      mu_plus,
        "mu_minus":     mu_minus,
        "lam_plus":     lam_plus,
        "lam_minus":    lam_minus,
        "mu":           mu,
        "lam":          lam,
        "mu_chi_P":     mu_chi_P,
        "musq_chi_P":   musq_chi_P,
        "semi_primes":  semi_primes,
    }
    return panel


# ---------- Gowers U^2 ------------------------------------------------------

def u2_pow4(f: np.ndarray) -> float:
    """||f||_{U^2}^4 = (1/N^4) sum_k |fhat(k)|^4 on Z/NZ."""
    N = len(f)
    F = np.fft.fft(f)
    return float(np.sum(np.abs(F) ** 4)) / (N ** 4)


# ---------- per-function record --------------------------------------------

def measure(name: str, f: np.ndarray, mu_chi_ref: np.ndarray) -> dict:
    N = len(f)
    mean_f = float(f.mean())
    l2sq = float((f * f).mean())
    u2_4_raw = u2_pow4(f)
    f_c = f - mean_f
    u2_4_centered = u2_pow4(f_c)

    rec = {
        "name": name,
        "N": N,
        "mean": mean_f,
        "L2_sq": l2sq,
        "u2_pow4_raw": u2_4_raw,
        "u2_pow4_centered": u2_4_centered,
    }

    # Q2 (only meaningful for non-negative indicator-like functions).
    if mean_f > 1e-12:
        rec["Q2"] = u2_4_raw / (mean_f ** 4)

    # Q2_norm = N * ||f||_{U^2}^4 / (2 * ||f||_2^4); ~1 if Gowers-uniform.
    if l2sq > 1e-30:
        rec["Q2_norm"] = N * u2_4_raw / (2.0 * (l2sq ** 2))
        rec["Q2_norm_centered"] = N * u2_4_centered / (2.0 * (l2sq ** 2))

    # Sanity: collapse check vs -chi_P.
    if name == "mu_chi_P":
        diff = float(np.abs(f + mu_chi_ref).max())
        rec["max_abs_diff_vs_neg_chi_P"] = diff

    if name == "musq_chi_P":
        diff = float(np.abs(f - mu_chi_ref).max())
        rec["max_abs_diff_vs_chi_P"] = diff

    return rec


def run_at_N(N: int) -> dict:
    print(f"\n{'='*72}\nN = {N}\n{'='*72}")
    t0 = time.time()
    panel = panel_arrays(N)
    t1 = time.time()
    print(f"  panel built in {t1-t0:.2f}s")

    chi_P = panel["chi_P"]
    out = {"N": N, "panel_build_s": t1 - t0, "records": {}}
    for name, f in panel.items():
        rec = measure(name, f, chi_P)
        out["records"][name] = rec
        mean_str = f"{rec['mean']:+.6f}"
        u2_str = f"{rec['u2_pow4_raw']:.4e}"
        cen_str = f"{rec['u2_pow4_centered']:.4e}"
        q2_str = f"{rec['Q2']:.4f}" if "Q2" in rec else "—"
        qn_str = f"{rec['Q2_norm']:.4f}" if "Q2_norm" in rec else "—"
        print(f"  {name:14s} mean={mean_str}  "
              f"||f||_U2^4={u2_str}  "
              f"centered={cen_str}  "
              f"Q2={q2_str:>9s}  Q2_norm={qn_str:>9s}")

    # Hardy-Littlewood reference values
    out["S2_chi_P"] = 2.300938
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--N-list", type=int, nargs="+",
                    default=[2048, 8192, 32768, 131072])
    args = ap.parse_args()

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)

    results = {"config": vars(args), "by_N": {}}
    for N in args.N_list:
        rec = run_at_N(N)
        results["by_N"][str(N)] = rec
        with open(args.out, "w") as f:
            json.dump(results, f, indent=2, default=str)

    print(f"\nWrote {args.out}")


if __name__ == "__main__":
    main()

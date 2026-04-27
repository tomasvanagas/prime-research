"""
§D6.b — Λ vs χ_P U^k comparison.

NOVELTY_CHALLENGES.md §D6.b. Compare Gowers U^2 / U^3 norms of:
  (a) bare χ_P (prime indicator)
  (b) bare Λ (log-weighted von Mangoldt)
  (c) ratios Q^k(Λ) / Q^k(χ_P)
  (d) centered variants χ_P/μ − 1 and Λ − 1
  (e) W-tricked variants for W ∈ {6, 30, 210}

Prediction (HL invariance):  Q^k(χ_P) ≈ Q^k(Λ) ≈ S_k  to leading order.

Pre-stated falsification — see lambda_vs_chi_p_uk_results.md.

Cite: E2.13 (Gowers norms of χ_P → S_k, S87).
Cross-domain ref: Gowers (2001); Green-Tao (2010, arXiv:math/0606088).
"""

from __future__ import annotations
import argparse
import json
import time
from pathlib import Path

import numpy as np
from sympy import primerange


# ---------- functions -------------------------------------------------------

def prime_indicator(N: int) -> np.ndarray:
    chi = np.zeros(N, dtype=np.float64)
    for p in primerange(2, N):
        chi[p] = 1.0
    return chi


def von_mangoldt(N: int) -> np.ndarray:
    """Λ(n) = log(p) if n = p^k, else 0; for n in [0, N)."""
    L = np.zeros(N, dtype=np.float64)
    for p in primerange(2, N):
        lp = float(np.log(p))
        pk = p
        while pk < N:
            L[pk] = lp
            pk *= p
    return L


def w_trick(f_full: np.ndarray, N_out: int, W: int, b: int) -> np.ndarray:
    """Returns g[n] = f_full[W*n + b] for n in [0, N_out)."""
    out = np.zeros(N_out, dtype=np.float64)
    for n in range(N_out):
        idx = W * n + b
        if 0 <= idx < len(f_full):
            out[n] = f_full[idx]
    return out


# ---------- Gowers norms ---------------------------------------------------

def u2_pow4(f: np.ndarray) -> float:
    """||f||_{U^2}^4 = (1/N^4) sum_k |fhat(k)|^4 on Z/NZ."""
    N = len(f)
    F = np.fft.fft(f)
    return float(np.sum(np.abs(F) ** 4)) / (N ** 4)


def u3_pow8(f: np.ndarray, max_h: int | None = None) -> float:
    """||f||_{U^3}^8 = (1/N) sum_h ||Δ_h f||_{U^2}^4, Δ_h f(x) = f(x)f(x+h)."""
    N = len(f)
    if max_h is None or max_h >= N:
        hs = range(N)
        denom = N
    else:
        step = max(1, N // max_h)
        hs = range(0, N, step)
        denom = len(hs)
    total = 0.0
    for h in hs:
        delta = f * np.roll(f, -h)
        total += u2_pow4(delta)
    return total / denom


# ---------- record ---------------------------------------------------------

def record(name: str, f: np.ndarray, do_u3: bool, max_h_u3: int | None,
           verbose: bool = True) -> dict:
    N = len(f)
    mu = float(f.mean())
    t0 = time.time()
    u2_4 = u2_pow4(f)
    t1 = time.time()
    rec = {
        "name": name,
        "N": N,
        "mean": mu,
        "u2_pow4": u2_4,
        "u2_time_s": t1 - t0,
    }
    if abs(mu) > 1e-12:
        rec["Q2"] = u2_4 / (mu ** 4)
    if do_u3:
        t2 = time.time()
        u3_8 = u3_pow8(f, max_h=max_h_u3)
        t3 = time.time()
        rec["u3_pow8"] = u3_8
        rec["u3_time_s"] = t3 - t2
        rec["u3_max_h"] = max_h_u3
        if abs(mu) > 1e-12:
            rec["Q3"] = u3_8 / (mu ** 8)
    if verbose:
        q2_str = f"{rec.get('Q2', float('nan')):.4f}" if 'Q2' in rec else "—"
        q3_str = f"{rec.get('Q3', float('nan')):.4f}" if 'Q3' in rec else "—"
        print(f"  {name:38s} N={N:6d} μ={mu:+.6f}  Q2={q2_str:>10s}  Q3={q3_str:>10s}  ({rec['u2_time_s']:.2f}s u2)")
    return rec


def centered_u2_pow4_from_Q(rec: dict) -> float | None:
    """Identity check: ||f/μ − 1||_{U^2}^4 = Q2(f) − 1."""
    if "Q2" in rec:
        return rec["Q2"] - 1.0
    return None


# ---------- runner ---------------------------------------------------------

def run_at_N(N: int, do_u3: bool, max_h_u3: int | None, do_w_trick: bool,
             do_centered: bool):
    print(f"\n{'='*70}\nN = {N}  (U^3: {do_u3}, W-trick: {do_w_trick})\n{'='*70}")

    chi = prime_indicator(N)
    lam = von_mangoldt(N)

    out = {"N": N, "do_u3": do_u3, "max_h_u3": max_h_u3, "records": {}}

    out["records"]["chi_P"] = record("chi_P", chi, do_u3, max_h_u3)
    out["records"]["Lambda"] = record("Lambda", lam, do_u3, max_h_u3)

    # Identity: ||(f/μ) − 1||_{U^k}^{2^k} = Q^k(f) − 1.  Check by direct compute.
    if do_centered:
        mu_chi = float(chi.mean())
        chi_c = chi / mu_chi - 1.0
        out["records"]["chi_P/mu - 1"] = record("chi_P/mu - 1", chi_c, do_u3, max_h_u3)

        mu_lam = float(lam.mean())
        lam_c = lam / mu_lam - 1.0
        out["records"]["Lambda/mu - 1"] = record("Lambda/mu - 1", lam_c, do_u3, max_h_u3)

        # also Λ - 1 (Λ has mean ≈ 1 by PNT, so this is approximately the same)
        lam_c1 = lam - 1.0
        out["records"]["Lambda - 1"] = record("Lambda - 1", lam_c1, do_u3, max_h_u3)

    # Comparison ratios at this N.
    qchi = out["records"]["chi_P"].get("Q2", float("nan"))
    qlam = out["records"]["Lambda"].get("Q2", float("nan"))
    out["Q2_ratio_Lambda_over_chi"] = qlam / qchi if qchi else float("nan")
    out["Q2_diff_pct"] = 100 * (qlam - qchi) / qchi if qchi else float("nan")
    print(f"  -- Q2(Λ)/Q2(χ_P) = {out['Q2_ratio_Lambda_over_chi']:.4f}, "
          f"diff = {out['Q2_diff_pct']:+.2f} %")
    if do_u3:
        qchi3 = out["records"]["chi_P"].get("Q3", float("nan"))
        qlam3 = out["records"]["Lambda"].get("Q3", float("nan"))
        out["Q3_ratio_Lambda_over_chi"] = qlam3 / qchi3 if qchi3 else float("nan")
        out["Q3_diff_pct"] = 100 * (qlam3 - qchi3) / qchi3 if qchi3 else float("nan")
        print(f"  -- Q3(Λ)/Q3(χ_P) = {out['Q3_ratio_Lambda_over_chi']:.4f}, "
              f"diff = {out['Q3_diff_pct']:+.2f} %")

    if do_w_trick:
        for W in [6, 30, 210]:
            b = 1
            full_len = W * N + b + 10
            chi_full = prime_indicator(full_len)
            lam_full = von_mangoldt(full_len)
            chi_w = w_trick(chi_full, N, W, b)
            lam_w = w_trick(lam_full, N, W, b)
            out["records"][f"chi_W{W}"] = record(f"chi_W{W}_b{b}", chi_w, do_u3, max_h_u3)
            out["records"][f"Lambda_W{W}"] = record(f"Lambda_W{W}_b{b}", lam_w, do_u3, max_h_u3)
            qcw = out["records"][f"chi_W{W}"].get("Q2", float("nan"))
            qlw = out["records"][f"Lambda_W{W}"].get("Q2", float("nan"))
            ratio = qlw / qcw if qcw else float("nan")
            print(f"  -- W={W:3d}: Q2(Λ_W)/Q2(χ_W) = {ratio:.4f}")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--N-list", type=int, nargs="+",
                    default=[1024, 4096, 16384, 65536])
    ap.add_argument("--N-list-u3", type=int, nargs="+",
                    default=[1024, 2048, 4096, 8192])
    ap.add_argument("--max-h-u3", type=int, default=None)
    ap.add_argument("--w-trick-N-cap", type=int, default=4096)
    ap.add_argument("--no-w-trick", action="store_true")
    ap.add_argument("--no-centered", action="store_true")
    args = ap.parse_args()

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)

    results = {"config": vars(args), "by_N": {}}
    for N in args.N_list:
        do_u3 = N in args.N_list_u3
        do_w_trick = (not args.no_w_trick) and (N <= args.w_trick_N_cap)
        rec = run_at_N(N, do_u3=do_u3, max_h_u3=args.max_h_u3,
                       do_w_trick=do_w_trick,
                       do_centered=not args.no_centered)
        results["by_N"][str(N)] = rec
        with open(args.out, "w") as f:
            json.dump(results, f, indent=2, default=str)

    print(f"\nWrote {args.out}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
F1.a refinement of E1.3 — cross-modulus generalisation of the bit-level
RH-scale anti-correlation valley discovered in S146.

S146 measured bit_J(p(n)) ≟ bit_J(round(Li^{-1}(n))) in base m=2 and
found a deep agreement-rate valley at J* = floor(log_2(p(N))/2). The
mechanism: under PNT-RH the predictor `Li^{-1}(n)` undershoots `p(n)`
by ~sqrt(p(n)) on the practical range; at bit-position J ~ J* the
rounding error has typical size 2^{J*}, flipping bit J* via a carry
with consistent sign.

This script measures the same predictor agreement in base m for
m ∈ {3, 5, 6, 30, 210}, asking whether the valley shifts to
J*(m) = floor(log_m(p(N))/2) — which would establish RH-shadow as
a generic predictor-mod-m-half-conductor phenomenon and refine E1.3.

Pre-stated falsification: see bit_J_pn_cross_modulus_results.md.
"""
from __future__ import annotations

import argparse
import json
import math
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent


def sieve_primes(limit: int) -> np.ndarray:
    """Sieve of Eratosthenes returning primes <= limit as np.uint64."""
    if limit < 2:
        return np.array([], dtype=np.uint64)
    sieve = np.ones(limit + 1, dtype=bool)
    sieve[:2] = False
    for i in range(2, int(limit ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i :: i] = False
    return np.nonzero(sieve)[0].astype(np.uint64)


def li_inverse_vectorised(n_arr: np.ndarray) -> np.ndarray:
    """Compute Li^{-1}(n) via Newton iteration on the asymptotic Li series.

    Identical numerics to S146's predictor; reused here for cross-modulus
    consistency.
    """
    n = n_arr.astype(np.float64)
    safe = np.maximum(n, 6.0)
    L0 = np.log(safe)
    LL = np.log(np.maximum(L0, 1.5))
    x = safe * (L0 + LL - 1.0)
    x = np.maximum(x, 2.0)
    for _ in range(40):
        x = np.maximum(x, 2.0)
        Lx = np.log(x)
        s = np.ones_like(x)
        term = np.ones_like(x)
        for k in range(1, 10):
            term *= k / Lx
            s += term
        Li = x / Lx * s
        Li_prime = 1.0 / Lx
        delta = (n - Li) / Li_prime
        x = x + delta
        x = np.maximum(x, 2.0)
    small = np.where(n_arr < 6)[0]
    base = np.array([2.0, 3.0, 5.0, 7.0, 11.0])
    for i in small:
        x[i] = base[int(n_arr[i]) - 1]
    return x


def digit_base_m(arr: np.ndarray, m: int, J: int) -> np.ndarray:
    """digit_J(x) := floor(x / m^J) mod m. (J=0 is the units digit.)"""
    if J == 0:
        return (arr % m).astype(np.int32)
    return ((arr // (m ** J)) % m).astype(np.int32)


def measure_for_m(
    primes: np.ndarray,
    li_round: np.ndarray,
    m: int,
    J_max: int,
) -> dict:
    """Per-J cross-modulus agreement & digit-distribution stats."""
    N = len(primes)
    last_p = int(primes[-1])
    log_m_p = math.log(last_p) / math.log(m)
    J_star_pred = int(log_m_p / 2)

    out: dict = {
        "m": int(m),
        "J_max": int(J_max),
        "N": int(N),
        "last_prime": int(last_p),
        "log_m_p": float(log_m_p),
        "J_star_predicted": int(J_star_pred),
        "baseline_uniform": 1.0 / m,
        "rows": [],
    }

    for J in range(J_max):
        d_p = digit_base_m(primes, m, J)
        d_l = digit_base_m(li_round, m, J)
        ag_li = float((d_p == d_l).mean())

        # Empirical digit distribution of p (informative for primorial m).
        bincount_p = np.bincount(d_p, minlength=m).astype(np.float64)
        probs_p = bincount_p / N
        max_p = float(probs_p.max())
        nonzero = probs_p[probs_p > 0]
        H_p = float(-np.sum(nonzero * np.log2(nonzero)))

        # Same for Li^{-1}.
        bincount_l = np.bincount(d_l, minlength=m).astype(np.float64)
        probs_l = bincount_l / N
        nonzero_l = probs_l[probs_l > 0]
        H_l = float(-np.sum(nonzero_l * np.log2(nonzero_l)))

        # Mean shift and shift histogram (digit_p - digit_l mod m).
        shift = (d_p - d_l) % m
        shift_hist = np.bincount(shift, minlength=m).astype(np.float64) / N

        out["rows"].append({
            "J": int(J),
            "ag_li": ag_li,
            "ag_const_majority": max_p,
            "H_p_bits": H_p,
            "H_l_bits": H_l,
            "shift_hist": shift_hist.tolist(),
        })

        log2m = math.log2(m)
        bar = "#" * max(0, int(ag_li * 60))
        baseline = 1.0 / m
        rel = ag_li / baseline
        marker = " <-- predicted J*" if J == J_star_pred else ""
        print(
            f"  m={m:3d} J={J:2d}  ag_li={ag_li:.4f}  "
            f"baseline=1/{m}={baseline:.4f}  rel={rel:.3f}  "
            f"H_p={H_p:.3f}/{log2m:.3f}b  {bar}{marker}",
            flush=True,
        )

    # Identify J* (argmin of ag_li) and characterise the dip.
    ag_arr = np.array([row["ag_li"] for row in out["rows"]])
    J_star_observed = int(np.argmin(ag_arr))
    out["J_star_observed"] = J_star_observed
    out["ag_li_at_J_star"] = float(ag_arr[J_star_observed])
    out["ag_li_relative"] = float(ag_arr[J_star_observed] / (1.0 / m))
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--limit", type=int, default=200_000_000)
    p.add_argument("--moduli", type=str, default="3,5,6,30,210")
    p.add_argument("--out", type=str,
                   default=str(HERE / "bit_J_pn_cross_modulus_results.json"))
    args = p.parse_args()

    moduli = [int(x) for x in args.moduli.split(",")]

    print(f"Sieving primes up to {args.limit:,}...", flush=True)
    t0 = time.time()
    primes = sieve_primes(args.limit)
    print(f"  Sieved {len(primes):,} primes in {time.time() - t0:.1f}s. "
          f"Last prime = {int(primes[-1]):,}", flush=True)

    n_arr = np.arange(1, len(primes) + 1, dtype=np.uint64)
    print("  Computing Li^{-1}(n)...", flush=True)
    t0 = time.time()
    li_pred = np.round(li_inverse_vectorised(n_arr)).astype(np.uint64)
    print(f"  Li^{{-1}} done in {time.time() - t0:.1f}s", flush=True)

    results: dict = {
        "limit": int(args.limit),
        "last_prime": int(primes[-1]),
        "N": int(len(primes)),
        "moduli": moduli,
        "modulus_results": [],
    }
    for m in moduli:
        J_max = int(math.log(int(primes[-1])) / math.log(m)) + 2
        print(f"\nMeasuring base-{m}, J = 0..{J_max - 1}  (predicted J* = "
              f"{int(math.log(int(primes[-1])) / math.log(m) / 2)}, "
              f"baseline 1/{m} = {1.0 / m:.4f})", flush=True)
        results["modulus_results"].append(
            measure_for_m(primes, li_pred, m, J_max)
        )

    with open(args.out, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nWrote {args.out}", flush=True)

    # Print headline summary table.
    print("\n=== HEADLINE: per-base dip summary ===", flush=True)
    print(f"{'m':>5}  {'log_m(p)':>10}  {'J*_pred':>8}  {'J*_obs':>8}  "
          f"{'ag(J*)':>8}  {'1/m':>8}  {'rel':>8}", flush=True)
    for r in results["modulus_results"]:
        print(f"{r['m']:>5}  {r['log_m_p']:>10.3f}  "
              f"{r['J_star_predicted']:>8d}  {r['J_star_observed']:>8d}  "
              f"{r['ag_li_at_J_star']:>8.4f}  {r['baseline_uniform']:>8.4f}  "
              f"{r['ag_li_relative']:>8.3f}", flush=True)


if __name__ == "__main__":
    main()

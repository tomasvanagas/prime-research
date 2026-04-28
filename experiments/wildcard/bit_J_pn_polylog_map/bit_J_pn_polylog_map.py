#!/usr/bin/env python3
"""
F1 refinement — bit_J(p(n)) computability map (LSB-side complement to E1.3).

E1.3 (S40, novel/carry_propagation_boundary.md) measures the agreement
rate between bit_J(p(n)) and bit_J(round(R^{-1}(n))) by bit position.
Top ~60% from MSB agree at >95%; bottom ~40% (including the LSB) drop
to ~50% (coin flip vs predictor). E1.3 calls this the "hard zone".

But agreement-with-R^{-1} conflates two distinct concepts:
    (a) R^{-1}(n) is a poor approximation at this bit position;
    (b) bit_J(p(n)) is intrinsically incompressible.

Counterexample: for n >= 2, bit_0(p(n)) = 1 (deterministic — primes are
odd). The bit is TRIVIALLY POLYLOG, but R^{-1}(n)'s LSB is essentially
uniform mod 2, so the agreement rate is ~50%. E1.3's predictor mis-
classifies bit_0 as "hard".

This script measures bit_J(p(n)) directly via three predictors:
    P_const(J): predict the empirically-majority bit at position J.
    P_PNT(J, n): predict bit_J( round(n * log n) ).
    P_li(J, n):  predict bit_J( round(Li^{-1}(n)) ).

Plus the intrinsic statistics (bias, peak discrepancy, lag-1 autocorr).

The map this produces:
    J = 0 (LSB):   trivial — P_const wins.
    J = 1..k:      Chebyshev/Dirichlet residue-class regime — bias is
                   O(sqrt(N)/log N) in a known direction.
    J = mid:       structural — neither P_const nor P_li dominates.
    J = top:       R^{-1}-predictable — P_li wins.

Pre-stated falsification:  see bit_J_pn_polylog_map_results.md.
"""

from __future__ import annotations
import json
import math
import time
from pathlib import Path

import numpy as np


HERE = Path(__file__).parent


def sieve_primes(limit: int) -> np.ndarray:
    """Sieve of Eratosthenes returning primes <= limit as np.uint32."""
    if limit < 2:
        return np.array([], dtype=np.uint32)
    sieve = np.ones(limit + 1, dtype=bool)
    sieve[:2] = False
    for i in range(2, int(limit ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i :: i] = False
    return np.nonzero(sieve)[0].astype(np.uint64)


def bit_J(arr: np.ndarray, J: int) -> np.ndarray:
    """Bit J of each element (LSB convention: J=0 is parity bit)."""
    return ((arr >> J) & 1).astype(np.int8)


def li_inverse(n: int) -> float:
    """Solve Li(x) = n for x via Newton's method on Li."""
    if n < 6:
        return [2.0, 3.0, 5.0, 7.0, 11.0][n - 1]
    # Asymptotic seed.
    x = float(n) * (math.log(n) + math.log(math.log(n)) - 1.0)
    for _ in range(60):
        # Compute Li(x) via series valid for large x.
        # Asymptotic expansion: Li(x) = x/log(x) * sum_{k>=0} k!/log^k(x)
        L = math.log(x)
        s = 1.0
        term = 1.0
        for k in range(1, 10):
            term *= k / L
            s += term
        Li = x / L * s
        Li_prime = 1.0 / L  # dLi/dx
        delta = (n - Li) / Li_prime
        if abs(delta) < 0.5:
            x += delta
            break
        x += delta
    return x


def li_inverse_vectorised(n_arr: np.ndarray) -> np.ndarray:
    """Compute Li^{-1}(n) for an array of n via Newton, vectorised."""
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
    # Hardcode tiny n.
    small = np.where(n_arr < 6)[0]
    base = np.array([2.0, 3.0, 5.0, 7.0, 11.0])
    for i in small:
        x[i] = base[int(n_arr[i]) - 1]
    return x


def measure_per_bit(primes: np.ndarray, J_max: int = 28) -> dict:
    """Compute statistics and predictor agreement for each J = 0..J_max-1."""
    N = len(primes)
    print(f"  N = {N} primes; max prime = {int(primes[-1])}", flush=True)

    # n indices (1-based) for predictor inputs.
    n_arr = np.arange(1, N + 1, dtype=np.uint64)
    pnt_pred = np.round(n_arr.astype(np.float64) * np.log(np.maximum(n_arr.astype(np.float64), 2.0))).astype(np.uint64)
    print("  Computing Li^{-1}(n) for all n...", flush=True)
    t0 = time.time()
    li_pred = np.round(li_inverse_vectorised(n_arr)).astype(np.uint64)
    print(f"  Li^{-1} done in {time.time()-t0:.1f}s", flush=True)

    out: dict = {"N": int(N), "J_max": int(J_max), "bits": []}
    for J in range(J_max):
        bj = bit_J(primes, J)
        bj_pnt = bit_J(pnt_pred, J)
        bj_li = bit_J(li_pred, J)

        bias = float(bj.mean()) - 0.5
        # Peak discrepancy of cumulative deviation from N/2.
        cum = np.cumsum(bj.astype(np.int64) - 0)
        # |sum_1^M b_J - M/2|
        M = np.arange(1, N + 1, dtype=np.float64)
        disc = cum.astype(np.float64) - M * 0.5
        peak_disc = float(np.max(np.abs(disc)))
        peak_disc_norm = peak_disc / math.sqrt(N)

        # Lag-1 autocorrelation.
        bm = bj.astype(np.float64) - 0.5
        denom = float(np.sum(bm * bm))
        if denom > 0:
            lag1 = float(np.sum(bm[:-1] * bm[1:]) / denom)
        else:
            lag1 = 0.0

        # Predictor agreement.
        agree_const_majority = max(float((bj == 0).mean()), float((bj == 1).mean()))
        agree_pnt = float((bj == bj_pnt).mean())
        agree_li = float((bj == bj_li).mean())

        # Chebyshev sign for J=1: bias toward p ≡ 3 (mod 4). For odd primes,
        # bit_1 of p in {1,3} mod 4 is 0 if p ≡ 1, 1 if p ≡ 3. Bias > 0 ⇒ favours ≡ 3.
        out["bits"].append({
            "J": J,
            "bias": bias,
            "peak_disc": peak_disc,
            "peak_disc_normalised_sqrtN": peak_disc_norm,
            "lag1_autocorr": lag1,
            "ag_const_majority": agree_const_majority,
            "ag_pnt": agree_pnt,
            "ag_li": agree_li,
        })

        bar = "#" * int(agree_li * 30)
        print(
            f"    J={J:2d}  bias={bias:+.4f}  peak/sqrtN={peak_disc_norm:6.2f}  "
            f"lag1={lag1:+.4f}  ag(const)={agree_const_majority:.3f}  "
            f"ag(PNT)={agree_pnt:.3f}  ag(Li^-1)={agree_li:.3f}  {bar}",
            flush=True,
        )

    return out


def main():
    import argparse

    p = argparse.ArgumentParser()
    p.add_argument("--limit", type=int, default=200_000_000,
                   help="sieve upper bound (primes up to this).")
    p.add_argument("--J-max", type=int, default=28)
    p.add_argument("--out", type=str, default=str(HERE / "bit_J_pn_results.json"))
    args = p.parse_args()

    print(f"Sieving primes up to {args.limit:,}...", flush=True)
    t0 = time.time()
    primes = sieve_primes(args.limit)
    print(f"  Sieved {len(primes):,} primes in {time.time()-t0:.1f}s. "
          f"Last prime = {int(primes[-1]):,}", flush=True)

    bits_used = math.ceil(math.log2(int(primes[-1]) + 1))
    print(f"  Last prime needs {bits_used} bits. Using J_max = {args.J_max}.", flush=True)

    print("\nMeasuring bit_J(p(n)) statistics and predictor agreement...", flush=True)
    res = measure_per_bit(primes, args.J_max)
    res["limit"] = int(args.limit)
    res["last_prime"] = int(primes[-1])
    res["bits_used"] = int(bits_used)

    with open(args.out, "w") as f:
        json.dump(res, f, indent=2)
    print(f"\nWrote {args.out}", flush=True)


if __name__ == "__main__":
    main()

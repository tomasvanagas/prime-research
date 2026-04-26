"""Proposal B — continued-fraction expansion of the prime constant.

  alpha = sum_{p prime} 2^{-p}    (binary digit n is 1 iff n prime)

If the regular CF [a_0; a_1, a_2, ...] of alpha has structural deviation
from Khintchine statistics (k-automaticity, low-frequency spectral
content, polynomial autocorrelation), then primality is polylog-
recoverable from CF data. If it is Khintchine-typical, this rules out
the CF route.
"""

from __future__ import annotations
import math
from collections import Counter
from pathlib import Path

import numpy as np
from mpmath import mp, mpf, floor


PRIME_LIMIT = 12000          # >= number of binary digits we want
TARGET_PARTIAL_QUOTIENTS = 1500   # how many a_n to extract


def sieve(n: int) -> list[int]:
    bs = bytearray(b"\x01") * (n + 1)
    bs[0] = bs[1] = 0
    for i in range(2, int(n**0.5) + 1):
        if bs[i]:
            bs[i * i :: i] = bytearray(len(bs[i * i :: i]))
    return [i for i in range(n + 1) if bs[i]]


def prime_constant(prec_bits: int) -> mpf:
    mp.prec = prec_bits + 64
    primes = sieve(prec_bits)
    half = mpf(1) / 2
    pw = half
    last_p = 1
    alpha = mpf(0)
    for p in primes:
        pw *= half ** (p - last_p)
        alpha += pw
        last_p = p
    return alpha


def continued_fraction(x: mpf, depth: int) -> list[int]:
    out = []
    for _ in range(depth):
        ai = int(floor(x))
        out.append(ai)
        frac = x - ai
        if frac == 0:
            break
        x = 1 / frac
        # if precision is exhausted (x explodes), bail.
        if mp.mag(x) > mp.prec - 16:
            break
    return out


def shannon_entropy(values: list[int]) -> float:
    c = Counter(values)
    n = sum(c.values())
    return -sum((v / n) * math.log2(v / n) for v in c.values())


def khintchine_geomean(a: list[int]) -> float:
    """log K_0 = (1/n) sum log a_i; should -> 0.987849... (log K_0 in nat) ≈
    Khintchine's constant K_0 = 2.685452..."""
    finite = [x for x in a if x > 0]
    return math.exp(sum(math.log(x) for x in finite) / len(finite))


def autocorrelation(a: list[int], lag: int) -> float:
    arr = np.array(a, dtype=float)
    arr = arr - arr.mean()
    n = len(arr) - lag
    num = np.sum(arr[: n] * arr[lag : lag + n])
    den = np.sum(arr * arr)
    return float(num / den) if den > 0 else 0.0


def power_spectrum(a: list[int]) -> np.ndarray:
    arr = np.log(np.maximum(np.array(a, dtype=float), 1.0))
    arr -= arr.mean()
    spec = np.abs(np.fft.rfft(arr)) ** 2
    return spec


def k_automaticity_proxy(a: list[int], k: int) -> float:
    """Heuristic test for k-automaticity: if a_n is k-automatic, then for
    each residue class r mod k^j the subsequence a_{kn+r} should match
    a global pattern. We compute Pearson correlation between the
    subsequence's empirical mean profile and the full sequence."""
    arr = np.array(a, dtype=float)
    # Compare the means of residue classes; for true k-automatic
    # sequences these will share a finite "shape" rather than match the
    # full distribution. We just record the inter-class variance.
    by_class = [arr[r::k] for r in range(k)]
    means = np.array([c.mean() for c in by_class])
    overall = arr.mean()
    inter_var = float(np.var(means))
    intra_var = float(np.mean([c.var() for c in by_class]))
    return inter_var / max(intra_var, 1e-12)


def main() -> None:
    print("# Proposal B — continued fraction of prime constant")
    bits = PRIME_LIMIT
    prec_bits = int(bits * 1.4)
    print(f"binary digits used: {bits}, mp.prec bits: {prec_bits}")
    alpha = prime_constant(prec_bits)
    print(f"alpha ~ {mp.nstr(alpha, 40)}")
    print()
    print(f"computing CF to depth {TARGET_PARTIAL_QUOTIENTS} ...")
    a = continued_fraction(alpha, TARGET_PARTIAL_QUOTIENTS)
    print(f"obtained {len(a)} partial quotients before precision exhaustion")
    print()
    print(f"first 30 partial quotients: {a[:30]}")
    print()
    # Khintchine geometric mean — for random reals this -> K_0 ≈ 2.6854.
    K_emp = khintchine_geomean(a[1:])      # skip a_0
    print(f"empirical Khintchine geometric mean: {K_emp:.4f}")
    print(f"  (theoretical for almost-every alpha:    2.6854...)")
    # Levy constant: (1/n) log q_n -> pi^2/(12 log 2) ≈ 1.18657.
    # We compute log q_n directly via the recurrence.
    q_prev, q_curr = 1, a[0] if a else 1
    log_q_growth = []
    for n in range(1, len(a)):
        q_new = a[n] * q_curr + q_prev
        q_prev, q_curr = q_curr, q_new
        if n % 50 == 0 and q_curr > 0:
            log_q_growth.append(math.log(q_curr) / n)
    if log_q_growth:
        print(f"empirical Levy constant (log q_n / n at n=last): "
              f"{log_q_growth[-1]:.4f}")
        print(f"  (theoretical:                                   1.1866)")
    print()
    # Distribution of small a_n values (Gauss-Kuzmin says
    # P(a_n = k) = log_2(1 + 1/(k(k+2))) ).
    counts = Counter(a)
    print("a_n distribution vs Gauss-Kuzmin (first 6 values):")
    print(f"{'k':>3} {'empirical':>10} {'GK pred':>10}")
    n_total = len(a)
    for k in range(1, 7):
        emp = counts.get(k, 0) / n_total
        gk = math.log2(1 + 1 / (k * (k + 2)))
        print(f"{k:>3} {emp:>10.4f} {gk:>10.4f}")
    print()
    # Autocorrelation at small lags.
    print("autocorrelation of a_n (lag 1..10):")
    for lag in range(1, 11):
        print(f"  lag {lag:>2}: {autocorrelation(a, lag):+.4f}")
    print()
    # Power spectrum slope (regress log|F|^2 vs log freq).
    spec = power_spectrum(a)
    spec = spec[1:]    # drop DC
    freqs = np.arange(1, len(spec) + 1)
    slope, _ = np.polyfit(np.log(freqs), np.log(spec + 1e-16), 1)
    print(f"power-spectrum slope (1/f^beta fit): beta = {-slope:.3f}")
    print(f"  (white noise: 0; pink: 1; brown: 2)")
    print()
    # k-automaticity proxy for k=2,3,5,7,10.
    print("k-automaticity proxy (high inter/intra variance ⇒ structure):")
    for k in (2, 3, 5, 7, 10):
        ratio = k_automaticity_proxy(a, k)
        print(f"  k={k:>2}: inter/intra = {ratio:.4f}")
    print()
    # Shannon entropy of the small-quotient distribution.
    H = shannon_entropy([x if x < 100 else 99 for x in a])
    print(f"Shannon entropy of a_n (bucketed): {H:.3f} bits")
    print()
    # Verdict.
    deviates = (
        abs(K_emp - 2.6854) > 0.5
        or any(abs(autocorrelation(a, lag)) > 0.1 for lag in range(1, 6))
        or -slope > 0.3
    )
    print(f"verdict: {'CANDIDATE STRUCTURE' if deviates else 'KHINTCHINE-TYPICAL'}")


if __name__ == "__main__":
    main()

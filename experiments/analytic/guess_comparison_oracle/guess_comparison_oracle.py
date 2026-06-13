#!/usr/bin/env python3
"""
guess_comparison_oracle.py — S491 follow-up: the "is my guess above or
below the real p(n)?" question, measured.

The comparison oracle [p(n) <= x?] is equivalent to [pi(x) >= n?], so
worst-case it is as hard as computing pi(x) (binary search recovers p(n)
from ~log x comparisons — E6.6 Aggarwal). But the *geography* of the
difficulty is the interesting part, and it is visible empirically:

  - against the BIASED canonical guess x = li^{-1}(n), the comparison
    bit is essentially constant (p(n) > guess for every n in any range
    reachable by computation — the li-vs-pi bias; Rubinstein-Sarnak 1994
    log-density of exceptions ~ 2.6e-7, first sign change near 1.4e316);
  - against the CENTERED guess x = R^{-1}(n) (inverse Riemann R), the
    comparison bit is an empirical fair coin. The better the guess, the
    less predictable the comparison — that IS the information barrier.

Also measured: |p(n) - guess| / sqrt(p(n)) stays O(polylog) for the
R-guess, i.e. the hard zone where comparisons cost a full pi(x)
computation has width ~sqrt(p(n)) — the familiar "last 50% of digits".

What would falsify the framing: the R-guess comparison bit showing a
stable predictable bias (would mean exploitable structure beyond R), or
li-guess sign flips at reachable scale (would contradict the known
Skewes-scale localisation).

Usage: python3 guess_comparison_oracle.py [--nmax 1000000] [--samples 40]
"""

import argparse

import numpy as np
from mpmath import mp, li, mpf
from mpmath import zeta as _zeta  # noqa: F401 (riemannr uses zeta internally)
import mpmath

mp.dps = 30


def sieve_primes(limit):
    s = np.ones(limit + 1, dtype=bool)
    s[:2] = False
    for p in range(2, int(limit ** 0.5) + 1):
        if s[p]:
            s[p * p::p] = False
    return np.flatnonzero(s)


def inverse(f, n, hi):
    """x with f(x) = n by bisection (f increasing)."""
    lo = mpf(2)
    hi = mpf(hi)
    for _ in range(120):
        mid = (lo + hi) / 2
        if f(mid) < n:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def main(nmax, samples):
    # primes up to a bound safely above p(nmax)
    import math
    bound = int(nmax * (math.log(nmax) + math.log(math.log(nmax))) * 1.1) + 100
    primes = sieve_primes(bound)
    assert len(primes) >= nmax

    ns = sorted(set(int(round(v)) for v in
                    np.logspace(3, np.log10(nmax), samples)))
    li_above = li_below = r_above = r_below = 0
    r_gap_norm = []
    print(f"{'n':>9} {'p(n)':>11} {'p-li_inv':>12} {'p-R_inv':>12} "
          f"{'|p-R_inv|/sqrt(p)':>18}")
    for n in ns:
        p = int(primes[n - 1])
        hi = 3 * bound
        g_li = inverse(lambda x: li(x), n, hi)
        g_r = inverse(lambda x: mpmath.riemannr(x), n, hi)
        d_li = p - float(g_li)
        d_r = p - float(g_r)
        li_above += d_li > 0
        li_below += d_li < 0
        r_above += d_r > 0
        r_below += d_r < 0
        r_gap_norm.append(abs(d_r) / p ** 0.5)
        print(f"{n:>9} {p:>11} {d_li:>12.1f} {d_r:>12.1f} "
              f"{abs(d_r)/p**0.5:>18.3f}")

    print(f"\nli-guess:  p(n) > li^-1(n) in {li_above}/{len(ns)} cases "
          f"(biased guess -> predictable comparison bit)")
    print(f"R-guess:   p(n) > R^-1(n)  in {r_above}/{len(ns)} cases, "
          f"below in {r_below}/{len(ns)} "
          f"(centered guess -> coin-flip comparison bit)")
    print(f"R-guess normalised gap |p - R^-1|/sqrt(p): "
          f"median {np.median(r_gap_norm):.3f}, max {max(r_gap_norm):.3f} "
          f"(hard-zone width ~ sqrt(p(n)) confirmed)")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmax", type=int, default=1000000)
    ap.add_argument("--samples", type=int, default=40)
    args = ap.parse_args()
    main(args.nmax, args.samples)

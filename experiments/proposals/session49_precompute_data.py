#!/usr/bin/env python3
"""
Session 49 — precompute primes p(n) and R^{-1}(n) for n=1..N and cache to npz.
This is a one-time helper that all three Proposal A/B/C scripts reuse.
"""
import math
import time
import os
import numpy as np
from sympy import primerange, mobius
from mpmath import mp, mpf, li, log as mlog, sqrt as msqrt

mp.dps = 18

CACHE = "/apps/aplikacijos/prime-research/experiments/proposals/session49_data.npz"
N = 8192


def first_N_primes(N):
    upper = max(20, int(N * (math.log(N) + math.log(math.log(N))) * 1.3))
    while True:
        ps = list(primerange(2, upper))
        if len(ps) >= N:
            return ps[:N]
        upper *= 2


def riemann_R(x):
    if x < 2:
        return mpf(0)
    s = mpf(0)
    k = 1
    while True:
        rt = mpf(x) ** (mpf(1) / k)
        if rt < mpf("1.05"):
            break
        mu = int(mobius(k))
        if mu != 0:
            s += mpf(mu) / k * li(rt)
        k += 1
        if k > 40:
            break
    return s


def riemann_R_prime_approx(x):
    """R'(x). Dominant term: 1/log(x). Accurate to ~1% for x > 100."""
    return 1.0 / float(mlog(mpf(x)))


def riemann_R_inverse(n, x_init=None):
    """Newton method seeded with x ~ n*log(n)*..."""
    n = mpf(n)
    if n < mpf(1):
        return mpf(2)
    if x_init is None:
        if n < 5:
            x = mpf(10)
        else:
            x = n * mlog(n) + n * mlog(mlog(n))
    else:
        x = mpf(x_init)
    for _ in range(20):
        fx = riemann_R(x) - n
        fpx = mpf(1) / mlog(x)
        step = fx / fpx
        x = x - step
        if abs(step) < mpf("1e-8"):
            break
    return float(x)


def main():
    if os.path.exists(CACHE):
        d = np.load(CACHE)
        if d["primes"].size >= N and d["rinv"].size >= N:
            print(f"Cache exists with {d['primes'].size} entries; skipping.")
            return

    print(f"Computing primes p(n) for n=1..{N}", flush=True)
    primes = first_N_primes(N)

    print("Computing R^{-1}(n) using Newton iteration with seed propagation", flush=True)
    rinv = np.zeros(N, dtype=np.float64)
    t0 = time.time()
    x_seed = None
    for i in range(N):
        rinv[i] = riemann_R_inverse(i + 1, x_init=x_seed)
        x_seed = rinv[i] + 1.0  # seed next iteration with the previous result
        if (i + 1) % 512 == 0:
            print(f"  n={i+1}, R^{{-1}}={rinv[i]:.2f}, t={time.time()-t0:.1f}s", flush=True)

    primes_arr = np.array(primes[:N], dtype=np.int64)
    delta = primes_arr.astype(np.float64) - rinv
    print(f"\ndelta range [{delta.min():.2f}, {delta.max():.2f}]  std={delta.std():.3f}")
    print(f"max |delta| = {np.abs(delta).max():.2f}")

    np.savez(CACHE, primes=primes_arr, rinv=rinv, delta=delta)
    print(f"Cached to {CACHE}")


if __name__ == "__main__":
    main()

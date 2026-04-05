#!/usr/bin/env python3
"""
Experiment 1: Compute f(x) = pi(x) - R(x) for x = 2..100000
Save the data for use by subsequent identity search experiments.

f(x) is the oscillatory residual of the prime counting function.
R(x) = sum_{k=1}^{inf} mu(k)/k * li(x^{1/k}) is Riemann's smooth approximation.
"""

import time
import json
import numpy as np
from mpmath import mpf, mp, li, log, sqrt, pi as MPI, power
from sympy import mobius as moebius

mp.dps = 30  # 30 digits sufficient for this range

XMAX = 100000

def sieve_pi(n):
    """Return array where pi_arr[x] = pi(x) for 0 <= x <= n."""
    is_prime = bytearray(b'\x01') * (n + 1)
    is_prime[0] = is_prime[1] = 0
    for i in range(2, int(n**0.5) + 1):
        if is_prime[i]:
            for j in range(i*i, n + 1, i):
                is_prime[j] = 0
    pi_arr = [0] * (n + 1)
    count = 0
    for i in range(n + 1):
        count += is_prime[i]
        pi_arr[i] = count
    return pi_arr, is_prime

def R_func(x, terms=80):
    """Riemann's R(x) = sum_{k=1}^{terms} mu(k)/k * li(x^{1/k})."""
    x = mpf(x)
    s = mpf(0)
    for k in range(1, terms + 1):
        mu_k = int(moebius(k))
        if mu_k == 0:
            continue
        xk = power(x, mpf(1) / k)
        if xk <= 1.0001:
            break
        s += mpf(mu_k) / k * li(xk)
    return s

print(f"Computing pi(x) sieve up to {XMAX}...")
t0 = time.time()
pi_arr, is_prime_arr = sieve_pi(XMAX)
print(f"  Sieve done in {time.time()-t0:.2f}s")
print(f"  pi({XMAX}) = {pi_arr[XMAX]}")

print(f"\nComputing R(x) for x=2..{XMAX}...")
t0 = time.time()

# Compute at sampled points for speed, plus store all for small range
# Full computation for x=2..100000 with mpmath is slow, so we:
# 1) Compute f(x) for ALL x=2..100000 (needed for some experiments)
# 2) Save to file

f_vals = {}  # x -> float(f(x))
R_vals = {}  # x -> float(R(x))

# For speed, compute R(x) in batches and report progress
batch_size = 5000
for batch_start in range(2, XMAX + 1, batch_size):
    batch_end = min(batch_start + batch_size, XMAX + 1)
    for x in range(batch_start, batch_end):
        Rx = R_func(x)
        fx = float(mpf(pi_arr[x]) - Rx)
        f_vals[x] = fx
        R_vals[x] = float(Rx)
    elapsed = time.time() - t0
    pct = (batch_end - 2) / (XMAX - 1) * 100
    print(f"  {pct:.0f}% done ({batch_end-1}/{XMAX}), elapsed {elapsed:.1f}s")

print(f"\nR(x) computation done in {time.time()-t0:.1f}s")

# Basic statistics
f_list = [f_vals[x] for x in range(2, XMAX + 1)]
print(f"\nf(x) statistics (x=2..{XMAX}):")
print(f"  min:  {min(f_list):.4f}")
print(f"  max:  {max(f_list):.4f}")
print(f"  mean: {np.mean(f_list):.6f}")
print(f"  std:  {np.std(f_list):.6f}")

# Save data
data = {
    'xmax': XMAX,
    'f_vals': {str(x): f_vals[x] for x in range(2, XMAX + 1)},
    'pi_vals': {str(x): pi_arr[x] for x in range(2, XMAX + 1)},
}

outfile = '/apps/aplikacijos/prime-research/experiments/algebraic/identity_search/fx_data.json'
print(f"\nSaving to {outfile}...")
with open(outfile, 'w') as fout:
    json.dump(data, fout)
print("Done.")

# Also save a compact numpy version
np.savez_compressed(
    '/apps/aplikacijos/prime-research/experiments/algebraic/identity_search/fx_data.npz',
    x=np.arange(2, XMAX + 1),
    f=np.array(f_list),
    pi=np.array([pi_arr[x] for x in range(2, XMAX + 1)]),
)
print("Saved numpy version too.")

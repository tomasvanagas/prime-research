#!/usr/bin/env python3
"""
Zero convergence analysis: How many zeta zeros does it take to
approximate delta(n) = p(n) - round(R^{-1}(n)) to within ±1?

This is the CORE experiment for Kt complexity: if K(n) zeros suffice,
then Kt(delta(n)|n) ~ K(n) * polylog. If K(n) ~ sqrt(p(n)), the barrier holds.

Uses the explicit formula: pi(x) = R(x) - sum_rho R(x^rho) + ...
where rho = 1/2 + i*gamma are nontrivial zeros of zeta.

For delta(n) specifically, we evaluate how the partial sum
S_K(x) = sum_{k=1}^{K} R(x^{rho_k}) + R(x^{bar{rho_k}})
converges to the true oscillatory correction.
"""

import numpy as np
import time
import sys
import os
from mpmath import mp, mpf, mpc, li, log, sqrt, pi as mpi, exp, fabs, re, im, power

mp.dps = 30

# Load zeta zeros
DATA_DIR = "/apps/aplikacijos/prime-research/data"
zeros_files = sorted([f for f in os.listdir(DATA_DIR) if 'zeros' in f.lower() or 'zero' in f.lower()])
print(f"Available data files: {zeros_files}")

# Try to load zeros
gammas = []
for fname in sorted(os.listdir(DATA_DIR)):
    fpath = os.path.join(DATA_DIR, fname)
    if os.path.isfile(fpath):
        try:
            with open(fpath, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line and not line.startswith('#'):
                        try:
                            val = float(line)
                            if val > 0:
                                gammas.append(val)
                        except ValueError:
                            pass
        except:
            pass

gammas = sorted(set(gammas))
print(f"Loaded {len(gammas)} unique zeta zeros")

if len(gammas) < 100:
    print("Not enough zeros in data files, computing first 200 zeros...")
    from mpmath import zetazero
    gammas = []
    for k in range(1, 201):
        z = zetazero(k)
        gammas.append(float(im(z)))
        if k % 50 == 0:
            print(f"  Computed {k} zeros...")
    print(f"Computed {len(gammas)} zeros")

gammas_mp = [mpf(str(g)) for g in gammas[:min(len(gammas), 1000)]]
NUM_ZEROS = len(gammas_mp)
print(f"Using {NUM_ZEROS} zeros for analysis")

# Load delta values
delta_all = np.load("/apps/aplikacijos/prime-research/experiments/information_theory/kt_delta_values.npy")
print(f"Loaded {len(delta_all)} delta values")

# Precompute primes for verification
from sympy import primerange
import math

N = len(delta_all)
upper = int(N * (math.log(N) + math.log(math.log(N)) + 2)) + 100
primes_list = list(primerange(2, upper + 1))[:N]
primes_arr = np.array(primes_list, dtype=np.int64)

# ---- Compute R(x) and R(x^rho) using mpmath ----
MU_MAX = 40
def compute_mobius(n):
    mu = [0] * (n + 1)
    mu[1] = 1
    is_prime_arr = [True] * (n + 1)
    ps = []
    for i in range(2, n + 1):
        if is_prime_arr[i]:
            ps.append(i)
            mu[i] = -1
        for p in ps:
            if i * p > n: break
            is_prime_arr[i * p] = False
            if i % p == 0:
                mu[i * p] = 0
                break
            else:
                mu[i * p] = -mu[i]
    return mu

MU = compute_mobius(MU_MAX)
MU_NZ = [(k, MU[k]) for k in range(1, MU_MAX + 1) if MU[k] != 0]

def R_mp(x):
    """Riemann R function with mpmath."""
    x = mpf(x)
    if x <= 1: return mpf(0)
    result = mpf(0)
    for k, mu_k in MU_NZ:
        xk = power(x, mpf(1)/k)
        if xk <= mpf('1.001'): break
        result += mpf(mu_k)/k * li(xk)
    return result

def R_mp_complex(x, rho):
    """R(x^rho) where rho = 1/2 + i*gamma, using mpmath."""
    x = mpf(x)
    if x <= 1: return mpc(0)
    result = mpc(0)
    for k, mu_k in MU_NZ:
        exponent = rho / k
        xk = power(x, exponent)
        if fabs(xk) < mpf('1.001'): break
        result += mpf(mu_k)/k * li(xk)
    return result

def oscillatory_sum_K(x, K):
    """
    Sum of R(x^rho) + R(x^{bar{rho}}) for first K zeros.
    This is the oscillatory correction to pi(x).
    Returns real value (imaginary parts cancel).
    """
    x = mpf(x)
    total = mpf(0)
    for k in range(min(K, NUM_ZEROS)):
        gamma = gammas_mp[k]
        rho = mpc(mpf('0.5'), gamma)
        rho_bar = mpc(mpf('0.5'), -gamma)
        contrib = R_mp_complex(x, rho) + R_mp_complex(x, rho_bar)
        total += re(contrib)
    return total

# ---- Main experiment: convergence of partial sums ----
print("\n" + "="*70)
print("EXPERIMENT 1: Zero convergence for specific n values")
print("="*70)

test_ns = [100, 500, 1000, 5000, 10000, 50000, 100000]
K_values = [1, 2, 5, 10, 20, 50, 100, 200, 500]
if NUM_ZEROS >= 1000:
    K_values.append(1000)

results = []
def log_result(msg):
    print(msg, flush=True)
    results.append(msg)

log_result(f"Zeros available: {NUM_ZEROS}")
log_result(f"Test n values: {test_ns}")
log_result(f"K values (number of zeros): {K_values}")
log_result("")

# For each test n, compute the explicit formula partial sum with K zeros
for n in test_ns:
    x = float(primes_arr[n-1])  # Use p(n) as the evaluation point for pi(x)
    # Actually, we should evaluate at R^{-1}(n) to get the correction
    # pi(R^{-1}(n)) = R(R^{-1}(n)) - sum_rho R(R^{-1}(n)^rho) + small terms
    # = n - sum_rho R(R^{-1}(n)^rho) + small terms
    # So delta(n) = p(n) - R^{-1}(n) ≈ correction involving zeros

    # But more directly: the explicit formula says
    # pi(x) = R(x) - sum_rho R(x^rho) - 1/ln(2) + integral term
    # So pi(p(n)) = n (since p(n) is the nth prime for counting purposes)
    # And R(p(n)) - sum_rho R(p(n)^rho) ≈ n
    # So sum_rho R(p(n)^rho) ≈ R(p(n)) - n

    true_delta = int(delta_all[n-1])
    R_x = float(R_mp(mpf(x)))
    oscillatory_true = R_x - n  # This should be close to the full zero sum

    log_result(f"\nn = {n}, p(n) = {x:.0f}, delta(n) = {true_delta}")
    log_result(f"  R(p(n)) = {R_x:.4f}, R(p(n)) - n = {oscillatory_true:.4f}")

    t0 = time.time()
    for K in K_values:
        if K > NUM_ZEROS:
            break
        osc_K = float(oscillatory_sum_K(mpf(x), K))
        error = oscillatory_true - osc_K
        elapsed = time.time() - t0
        log_result(f"  K={K:4d}: osc_sum={osc_K:+10.4f}, error={error:+10.4f}, |error|={abs(error):.4f}")
        if elapsed > 120:  # 2 min timeout per n
            log_result(f"  [TIMEOUT after {elapsed:.0f}s, stopping at K={K}]")
            break

    log_result(f"  Time for n={n}: {time.time()-t0:.1f}s")

# ---- Experiment 2: For each n, find minimum K such that |error| < 1 ----
print("\n" + "="*70)
print("EXPERIMENT 2: Minimum zeros needed for |error| < 1")
print("="*70)

log_result("\n" + "="*70)
log_result("EXPERIMENT 2: Minimum zeros needed for |error| < 1")
log_result("="*70)

# Test a range of n values
test_ns2 = list(range(100, 1001, 100)) + list(range(2000, 10001, 1000))
min_K_for_accuracy = []

for n in test_ns2:
    x = float(primes_arr[n-1])
    R_x = float(R_mp(mpf(x)))
    target = R_x - n

    cumsum = 0.0
    min_K = -1
    t0 = time.time()
    for K in range(1, min(NUM_ZEROS + 1, 501)):
        gamma = gammas_mp[K-1]
        rho = mpc(mpf('0.5'), gamma)
        rho_bar = mpc(mpf('0.5'), -gamma)
        contrib = float(re(R_mp_complex(mpf(x), rho) + R_mp_complex(mpf(x), rho_bar)))
        cumsum += contrib
        error = abs(target - cumsum)
        if error < 1.0 and min_K == -1:
            min_K = K
            break  # First time error < 1
        if time.time() - t0 > 30:  # 30s timeout per n
            break

    if min_K > 0:
        min_K_for_accuracy.append((n, x, min_K))
        log_result(f"  n={n:6d}, p(n)={x:10.0f}: K_min = {min_K:4d} zeros for |error| < 1")
    else:
        min_K_for_accuracy.append((n, x, -1))
        log_result(f"  n={n:6d}, p(n)={x:10.0f}: K_min > {min(NUM_ZEROS, 500)} (not reached)")

# ---- Experiment 3: Scaling of K_min ----
log_result("\n" + "="*70)
log_result("EXPERIMENT 3: Scaling analysis of K_min(n)")
log_result("="*70)

valid = [(n, x, K) for n, x, K in min_K_for_accuracy if K > 0]
if len(valid) > 5:
    ns_valid = np.array([v[0] for v in valid], dtype=float)
    xs_valid = np.array([v[1] for v in valid], dtype=float)
    ks_valid = np.array([v[2] for v in valid], dtype=float)

    # Fit K_min ~ n^alpha
    log_ns = np.log(ns_valid)
    log_ks = np.log(ks_valid)
    alpha, c = np.polyfit(log_ns, log_ks, 1)
    log_result(f"K_min ~ n^{alpha:.4f} (intercept={c:.4f})")

    # Fit K_min ~ x^alpha
    log_xs = np.log(xs_valid)
    alpha_x, c_x = np.polyfit(log_xs, log_ks, 1)
    log_result(f"K_min ~ x^{alpha_x:.4f} (intercept={c_x:.4f})")

    # Fit K_min ~ log(n)^alpha
    loglog_ns = np.log(np.log(ns_valid))
    alpha_ll, c_ll = np.polyfit(loglog_ns, log_ks, 1)
    log_result(f"K_min ~ log(n)^{alpha_ll:.4f} (intercept={c_ll:.4f})")

    # Fit K_min ~ sqrt(x)
    sqrt_xs = np.sqrt(xs_valid)
    corr_sqrt = np.corrcoef(sqrt_xs, ks_valid)[0, 1]
    log_result(f"Correlation(K_min, sqrt(x)): {corr_sqrt:.4f}")

    log_result(f"\nExpected for barrier: K_min ~ sqrt(x) ~ n^0.5 * polylog")
    log_result(f"Measured: K_min ~ n^{alpha:.4f}")
    log_result(f"If alpha ≈ 0.5: BARRIER CONFIRMED (need sqrt(x) zeros)")
    log_result(f"If alpha << 0.5: POSSIBLE SHORTCUT (fewer zeros suffice)")
else:
    log_result(f"Not enough valid data points ({len(valid)}) for scaling analysis")

# Save results
output_path = "/apps/aplikacijos/prime-research/experiments/information_theory/kt_complexity/zero_convergence_results.txt"
with open(output_path, 'w') as f:
    f.write('\n'.join(results))
print(f"\nResults saved to {output_path}")

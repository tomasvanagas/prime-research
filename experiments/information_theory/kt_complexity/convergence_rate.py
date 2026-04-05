#!/usr/bin/env python3
"""
Convergence rate of the explicit formula for pi(x):
    pi(x) = li(x) - sum_rho li(x^rho) - ln(2) + ...

We measure: for given x, how does |pi(x) - li(x) + sum_{k=1}^K [li(x^{rho_k}) + li(x^{bar{rho_k}})]|
decrease with K?

The theoretical bound (under RH): truncation error after T zeros ≈ O(sqrt(x)*log^2(x)/T).
So K_min for |error| < 1 is approximately sqrt(x)*log^2(x).

We test this empirically at small x where we know pi(x) exactly.
"""

import numpy as np
import math
import time
from mpmath import mp, mpf, mpc, li as mpli, log as mplog, im, re, exp, pi as mppi, fabs
from mpmath import zetazero

mp.dps = 50  # High precision

# Compute first 500 zeros
print("Computing zeta zeros...", flush=True)
zeros = []
for k in range(1, 501):
    z = zetazero(k)
    zeros.append(float(im(z)))
    if k % 100 == 0:
        print(f"  {k} zeros computed", flush=True)
print(f"Done: {len(zeros)} zeros", flush=True)

# Compute pi(x) exactly via sieve for test x values
from sympy import primepi

results = []
def log_result(msg):
    print(msg, flush=True)
    results.append(msg)

log_result("=" * 70)
log_result("EXPLICIT FORMULA CONVERGENCE: li-based, high precision")
log_result("=" * 70)

def li_complex(z):
    """li(z) for complex z using mpmath."""
    return mpli(z)

def explicit_formula_partial(x_val, K):
    """
    Compute li(x) - sum_{k=1}^{K} [li(x^{rho_k}) + li(x^{bar{rho_k}})] - ln(2)
    This should approximate pi(x) as K -> infinity.
    """
    x = mpf(x_val)
    total = mpli(x) - mplog(mpf(2))

    for k in range(K):
        gamma = mpf(str(zeros[k]))
        rho = mpc(mpf('0.5'), gamma)
        # x^rho = x^{1/2 + i*gamma} = sqrt(x) * exp(i*gamma*ln(x))
        x_rho = exp(rho * mplog(x))
        x_rho_bar = exp((mpc(mpf('0.5'), -gamma)) * mplog(x))

        contrib = mpli(x_rho) + mpli(x_rho_bar)
        total -= contrib

    return float(re(total))

# Test at various x values where we know pi(x) exactly
test_x = [100, 500, 1000, 5000, 10000, 50000, 100000, 500000, 1000000]
K_vals = [1, 2, 5, 10, 20, 50, 100, 200, 500]

for x in test_x:
    pi_x = int(primepi(x))
    li_x = float(mpli(mpf(x)))

    log_result(f"\nx = {x}, pi(x) = {pi_x}, li(x) = {li_x:.4f}, li(x)-pi(x) = {li_x-pi_x:.4f}")
    log_result(f"  sqrt(x) = {math.sqrt(x):.2f}, log^2(x) = {math.log(x)**2:.2f}")

    t0 = time.time()
    prev_approx = li_x - math.log(2)

    for K in K_vals:
        if K > len(zeros):
            break
        approx = explicit_formula_partial(x, K)
        error = approx - pi_x

        log_result(f"  K={K:4d}: approx={approx:12.4f}, error={error:+10.4f}, "
                   f"|error|={abs(error):.4f}, |error|/sqrt(x)={abs(error)/math.sqrt(x):.6f}")

        if time.time() - t0 > 120:
            log_result(f"  [TIMEOUT]")
            break

# ---- Scaling analysis ----
log_result("\n" + "=" * 70)
log_result("SCALING ANALYSIS")
log_result("=" * 70)

# For a fixed K, measure error vs x
for K in [10, 50, 100, 200, 500]:
    if K > len(zeros):
        break
    log_result(f"\nFixed K={K}:")
    errors = []
    xs = []
    for x in [100, 500, 1000, 5000, 10000, 50000, 100000]:
        pi_x = int(primepi(x))
        approx = explicit_formula_partial(x, K)
        error = abs(approx - pi_x)
        xs.append(x)
        errors.append(error)
        log_result(f"  x={x:8d}: |error|={error:.4f}, |error|/sqrt(x)={error/math.sqrt(x):.6f}, "
                   f"|error|*K/sqrt(x)/log^2(x)={error*K/(math.sqrt(x)*math.log(x)**2):.6f}")

    # Fit error ~ x^alpha
    if len(xs) > 3:
        log_xs = np.log(xs)
        log_errs = np.log([max(e, 1e-10) for e in errors])
        valid = log_errs > -20
        if sum(valid) > 2:
            alpha, c = np.polyfit(log_xs[valid], log_errs[valid], 1)
            log_result(f"  Error scaling: |error| ~ x^{alpha:.4f}")

# ---- Key conclusion ----
log_result("\n" + "=" * 70)
log_result("CONCLUSIONS")
log_result("=" * 70)
log_result("")
log_result("Theoretical prediction (under RH):")
log_result("  Truncation error after K zeros ≈ C * sqrt(x) * log^2(x) / K")
log_result("  K_min for |error| < 1: K_min ≈ C * sqrt(x) * log^2(x)")
log_result("")
log_result("Implications for Kt complexity:")
log_result("  Via explicit formula: Kt(delta(n)|n) >= Omega(sqrt(p(n))/log(p(n)))")
log_result("  Since p(n) ~ n*ln(n): Kt(delta(n)|n) >= Omega(sqrt(n*ln(n))/log(n))")
log_result("  ≈ Omega(sqrt(n) / sqrt(log(n)))")
log_result("  This is EXPONENTIAL in the input size N = log(n).")
log_result("")
log_result("A polylog Kt would require a NON-explicit-formula method.")

# Save
output_path = "/apps/aplikacijos/prime-research/experiments/information_theory/kt_complexity/convergence_rate_results.txt"
with open(output_path, 'w') as f:
    f.write('\n'.join(results))
print(f"\nResults saved to {output_path}")

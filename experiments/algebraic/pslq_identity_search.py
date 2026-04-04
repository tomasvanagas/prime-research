#!/usr/bin/env python3
"""
PSLQ / LLL Identity Search for f(x) = pi(x) - R(x)

Searches for algebraic, linear, recurrence, modular, functional,
and discrete-derivative identities involving the oscillatory residual.

Uses mpmath.pslq() for integer relation detection and basic LLL-style
polynomial relation searches.
"""

import sys
import math
import time
from collections import defaultdict

from mpmath import (
    mpf, mp, li, log, sqrt, sin, cos, pi as MPI, euler as EULER,
    pslq, fsum, power, fabs, nstr, matrix, inf
)
from sympy import mobius as moebius

# Precision
mp.dps = 50

# ─────────────────────────────────────────────────────────────────────
# 1. Compute pi(x) via sieve and R(x) via Gram series
# ─────────────────────────────────────────────────────────────────────

XMAX = 10000

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
    return pi_arr

def R_func(x, terms=100):
    """Riemann's R(x) = sum_{k=1}^{terms} mu(k)/k * li(x^{1/k})."""
    x = mpf(x)
    s = mpf(0)
    for k in range(1, terms + 1):
        mu_k = int(moebius(k))
        if mu_k == 0:
            continue
        xk = power(x, mpf(1) / k)
        if xk <= 1:
            break
        s += mpf(mu_k) / k * li(xk)
    return s

print("Computing pi(x) sieve up to", XMAX, "...")
t0 = time.time()
pi_arr = sieve_pi(XMAX)
print(f"  Sieve done in {time.time()-t0:.2f}s")

print("Computing R(x) for x=2..{} ...".format(XMAX))
t0 = time.time()
R_arr = [mpf(0)] * (XMAX + 1)
f_arr = [mpf(0)] * (XMAX + 1)  # f(x) = pi(x) - R(x)
for x in range(2, XMAX + 1):
    R_arr[x] = R_func(x)
    f_arr[x] = mpf(pi_arr[x]) - R_arr[x]
print(f"  R(x) done in {time.time()-t0:.2f}s")

# Quick stats
f_vals = [float(f_arr[x]) for x in range(2, XMAX + 1)]
print(f"  f(x) range: [{min(f_vals):.4f}, {max(f_vals):.4f}]")
print(f"  f(x) mean:  {sum(f_vals)/len(f_vals):.6f}")
print()

results = []  # (section, description, residual, detail)

def record(section, desc, residual, detail=""):
    results.append((section, desc, float(residual), detail))
    tag = "***INTERESTING***" if float(residual) < 1e-10 else ""
    print(f"  [{section}] {desc}: residual={float(residual):.6e} {tag}")

# ─────────────────────────────────────────────────────────────────────
# 2a. Linear relations via PSLQ at selected x values
# ─────────────────────────────────────────────────────────────────────
print("=" * 70)
print("SECTION A: Linear relations via PSLQ")
print("=" * 70)

gamma1_zero = mpf("14.134725141734693790457251983562470270784257115699")

test_x_values = [100, 200, 500, 1000, 2000, 5000, 10000]

for x in test_x_values:
    xm = mpf(x)
    fx = f_arr[x]

    # Build vector: [f(x), 1, log(x), sqrt(x), x^{1/3}, x^{1/4}, x^{1/5},
    #                li(x), li(sqrt(x)), li(x^{1/3}),
    #                sin(gamma1*log(x)), cos(gamma1*log(x))]
    vec = [
        fx,
        mpf(1),
        log(xm),
        sqrt(xm),
        power(xm, mpf(1)/3),
        power(xm, mpf(1)/4),
        power(xm, mpf(1)/5),
        li(xm),
        li(sqrt(xm)),
        li(power(xm, mpf(1)/3)),
        sin(gamma1_zero * log(xm)),
        cos(gamma1_zero * log(xm)),
    ]

    labels = ["f(x)", "1", "log(x)", "sqrt(x)", "x^{1/3}", "x^{1/4}", "x^{1/5}",
              "li(x)", "li(sqrt(x))", "li(x^{1/3})",
              "sin(g1*log(x))", "cos(g1*log(x))"]

    rel = pslq(vec, maxcoeff=10**6, maxsteps=2000)
    if rel is not None:
        # Check quality: coefficient of f(x) should be nonzero
        if rel[0] != 0:
            # Residual: the relation says sum(rel[i]*vec[i]) = 0
            residual = fabs(fsum(mpf(r)*v for r, v in zip(rel, vec)))
            terms = " + ".join(f"{r}*{l}" for r, l in zip(rel, labels) if r != 0)
            record("A-PSLQ", f"x={x}: {terms}", residual, str(rel))
        else:
            record("A-PSLQ", f"x={x}: relation found but coeff(f)=0 (spurious)", mpf(1), str(rel))
    else:
        record("A-PSLQ", f"x={x}: no relation found", mpf(1))

print()

# ─────────────────────────────────────────────────────────────────────
# 2b. Algebraic / polynomial relations via exhaustive search
# ─────────────────────────────────────────────────────────────────────
print("=" * 70)
print("SECTION B: Polynomial / algebraic relations")
print("=" * 70)

# Test: f(x)^d + linear combo of {1, sqrt(x), x^{1/3}, log(x)} = 0
# for degrees d=2,3,4 at multiple x, using PSLQ

for deg in [2, 3, 4]:
    print(f"\n  Degree {deg} polynomial in f(x):")
    for x in [500, 1000, 5000, 10000]:
        xm = mpf(x)
        fx = f_arr[x]
        # vector: [f^d, f^{d-1}, ..., f, 1, sqrt(x), x^{1/3}, log(x)]
        vec = [fx**k for k in range(deg, -1, -1)]
        vec += [sqrt(xm), power(xm, mpf(1)/3), log(xm)]

        labels = [f"f^{k}" for k in range(deg, -1, -1)]
        labels += ["sqrt(x)", "x^{1/3}", "log(x)"]

        rel = pslq(vec, maxcoeff=10**8, maxsteps=2000)
        if rel is not None and rel[0] != 0:
            residual = fabs(fsum(mpf(r)*v for r, v in zip(rel, vec)))
            terms = " + ".join(f"{r}*{l}" for r, l in zip(rel, labels) if r != 0)
            record("B-POLY", f"deg={deg}, x={x}: {terms}", residual, str(rel))
        else:
            tag = "(coeff=0)" if rel is not None else "(none)"
            record("B-POLY", f"deg={deg}, x={x}: no useful relation {tag}", mpf(1))

print()

# ─────────────────────────────────────────────────────────────────────
# 2c. Recurrence relations
# ─────────────────────────────────────────────────────────────────────
print("=" * 70)
print("SECTION C: Recurrence relations")
print("=" * 70)

# Test linear recurrence f(x) = a1*f(x-1) + a2*f(x-2) + ... + ak*f(x-k)
# over a range of x values (overdetermined system, least squares)

import numpy as np

f_float = np.array([float(f_arr[x]) for x in range(2, XMAX + 1)])

for order in [2, 3, 4, 5, 6, 8, 10]:
    # Build system: each row is [f(x-1), f(x-2), ..., f(x-order)] = f(x)
    n = len(f_float) - order
    A = np.zeros((n, order))
    b = np.zeros(n)
    for i in range(n):
        for j in range(order):
            A[i, j] = f_float[i + order - 1 - j]
        b[i] = f_float[i + order]

    coeffs, res, rank, sv = np.linalg.lstsq(A, b, rcond=None)
    predicted = A @ coeffs
    max_err = np.max(np.abs(predicted - b))
    rms_err = np.sqrt(np.mean((predicted - b)**2))
    rel_rms = rms_err / np.std(f_float)

    record("C-RECUR", f"order={order}: rms={rms_err:.6f}, max={max_err:.6f}, rel_rms={rel_rms:.6f}",
           mpf(rms_err), f"coeffs={coeffs.tolist()[:4]}...")

# Test: f(x) = a*f(x-1) + b*f(x-2) + c + d*log(x) + e/sqrt(x)
print("\n  Augmented recurrence (order 2 + log(x), 1/sqrt(x), 1):")
n = len(f_float) - 2
A = np.zeros((n, 5))
b = np.zeros(n)
for i in range(n):
    x = i + 4  # actual x value (starting from x=4 since f starts at x=2)
    A[i, 0] = f_float[i + 1]  # f(x-1)
    A[i, 1] = f_float[i]      # f(x-2)
    A[i, 2] = 1.0
    A[i, 3] = math.log(x)
    A[i, 4] = 1.0 / math.sqrt(x)
    b[i] = f_float[i + 2]

coeffs, res, rank, sv = np.linalg.lstsq(A, b, rcond=None)
predicted = A @ coeffs
rms_err = np.sqrt(np.mean((predicted - b)**2))
max_err = np.max(np.abs(predicted - b))
record("C-RECUR", f"aug_order2: rms={rms_err:.6f}, max={max_err:.6f}",
       mpf(rms_err), f"coeffs={coeffs.tolist()}")

print()

# ─────────────────────────────────────────────────────────────────────
# 2d. Modular identities
# ─────────────────────────────────────────────────────────────────────
print("=" * 70)
print("SECTION D: Modular identities")
print("=" * 70)

# f(x) is not integer in general, but pi(x) is integer and R(x) is real.
# We test: round(f(x)) mod m for patterns, and also pi(x) mod m vs x mod m.

# Since f(x) = pi(x) - R(x), and pi(x) is integer, let's look at
# floor(f(x)) mod m and fractional part patterns.

f_rounded = [round(float(f_arr[x])) for x in range(2, XMAX + 1)]
f_frac = [float(f_arr[x]) - round(float(f_arr[x])) for x in range(2, XMAX + 1)]

for m in [2, 3, 4, 5, 6, 7, 8, 12]:
    counts = defaultdict(int)
    for val in f_rounded:
        counts[val % m] += 1
    # Chi-squared test vs uniform
    expected = len(f_rounded) / m
    chi2 = sum((c - expected)**2 / expected for c in counts.values())
    dof = m - 1
    # p-value approximation: chi2/dof >> 1 means non-uniform
    record("D-MOD", f"round(f(x)) mod {m}: chi2/dof={chi2/dof:.4f}",
           mpf(1) if chi2/dof > 2 else mpf(chi2/dof),
           f"distribution={dict(sorted(counts.items()))}")

# Fractional part distribution
print("\n  Fractional part of f(x) statistics:")
f_frac_arr = np.array(f_frac)
print(f"    mean(frac) = {np.mean(f_frac_arr):.6f}")
print(f"    std(frac)  = {np.std(f_frac_arr):.6f}")
# Is fractional part clustered near specific values?
from collections import Counter
frac_bins = Counter(round(ff, 1) for ff in f_frac)
print(f"    frac part histogram (0.1 bins): {dict(sorted(frac_bins.items())[:10])}")

# Test: does f(x) mod 1 correlate with x mod small m?
for m in [2, 3, 5]:
    groups = defaultdict(list)
    for i, x in enumerate(range(2, XMAX + 1)):
        groups[x % m].append(f_frac[i])
    means = {k: np.mean(v) for k, v in groups.items()}
    spread = max(means.values()) - min(means.values())
    record("D-MOD", f"frac(f(x)) grouped by x mod {m}: spread={spread:.6f}",
           mpf(spread), f"group_means={means}")

print()

# ─────────────────────────────────────────────────────────────────────
# 2e. Functional equations: f(x*y) vs f(x), f(y)
# ─────────────────────────────────────────────────────────────────────
print("=" * 70)
print("SECTION E: Functional equations")
print("=" * 70)

# Test: f(xy) = a*f(x) + b*f(y) + c for fixed small y
for y in [2, 3, 5, 7]:
    xs = []
    fxy_vals = []
    fx_vals = []
    fy = float(f_arr[y])
    for x in range(2, XMAX // y + 1):
        xy = x * y
        if xy <= XMAX:
            xs.append(x)
            fxy_vals.append(float(f_arr[xy]))
            fx_vals.append(float(f_arr[x]))

    fxy_np = np.array(fxy_vals)
    # Fit: f(xy) = a*f(x) + b*f(y) + c = a*f(x) + (b*fy + c)
    # Simplified: f(xy) = a*f(x) + d  (since f(y) is constant)
    A = np.column_stack([np.array(fx_vals), np.ones(len(fx_vals))])
    coeffs, _, _, _ = np.linalg.lstsq(A, fxy_np, rcond=None)
    pred = A @ coeffs
    rms = np.sqrt(np.mean((pred - fxy_np)**2))
    record("E-FUNC", f"f({y}x) = a*f(x)+b: a={coeffs[0]:.6f}, b={coeffs[1]:.6f}, rms={rms:.6f}",
           mpf(rms))

# Test: f(x^2) vs f(x) (where x^2 <= XMAX)
xs2, fxs, fx2s = [], [], []
for x in range(2, int(XMAX**0.5) + 1):
    x2 = x * x
    if x2 <= XMAX:
        xs2.append(x)
        fxs.append(float(f_arr[x]))
        fx2s.append(float(f_arr[x2]))

fxs_np = np.array(fxs)
fx2s_np = np.array(fx2s)
A = np.column_stack([fxs_np, fxs_np**2, np.ones(len(fxs_np))])
coeffs, _, _, _ = np.linalg.lstsq(A, fx2s_np, rcond=None)
pred = A @ coeffs
rms = np.sqrt(np.mean((pred - fx2s_np)**2))
record("E-FUNC", f"f(x^2) = a*f(x)+b*f(x)^2+c: rms={rms:.6f}",
       mpf(rms), f"coeffs={coeffs.tolist()}")

print()

# ─────────────────────────────────────────────────────────────────────
# 2f. Discrete derivative patterns
# ─────────────────────────────────────────────────────────────────────
print("=" * 70)
print("SECTION F: Discrete derivative f(x+1)-f(x)")
print("=" * 70)

# Compute Df(x) = f(x+1) - f(x) and D2f = Df(x+1) - Df(x)
Df = np.array([f_float[i+1] - f_float[i] for i in range(len(f_float)-1)])
D2f = np.array([Df[i+1] - Df[i] for i in range(len(Df)-1)])

print(f"  Df stats: mean={np.mean(Df):.6f}, std={np.std(Df):.6f}")
print(f"  D2f stats: mean={np.mean(D2f):.6f}, std={np.std(D2f):.6f}")

# Df(x) is related to whether x+1 is prime minus the change in R(x).
# Since pi(x+1) - pi(x) = 1 if x+1 is prime, 0 otherwise,
# Df(x) = [x+1 is prime] - (R(x+1) - R(x))
# R(x+1) - R(x) ~ 1/log(x) for large x.
# So Df(x) ~ [x+1 is prime] - 1/log(x)

# Test: is Df(x) well-modeled by [x+1 is prime] - 1/log(x)?
xs_df = np.arange(2, XMAX)  # x values for Df
is_prime_xp1 = np.array([1.0 if pi_arr[x+1] > pi_arr[x] else 0.0 for x in range(2, XMAX)])
inv_logx = np.array([1.0 / math.log(x) for x in range(2, XMAX)])
model_Df = is_prime_xp1 - inv_logx
rms_Df_model = np.sqrt(np.mean((Df[:len(model_Df)] - model_Df)**2))
record("F-DISC", f"Df(x) vs [x+1 prime]-1/log(x): rms={rms_Df_model:.8f}",
       mpf(rms_Df_model))

# Better model: Df(x) = [x+1 prime] - (R(x+1)-R(x))
# This should be exactly 0 by definition! Let's verify.
DR = np.array([float(R_arr[x+1] - R_arr[x]) for x in range(2, min(XMAX, 2000))])
Df_exact = is_prime_xp1[:len(DR)] - DR
rms_exact = np.sqrt(np.mean(Df_exact**2))
max_exact = np.max(np.abs(Df_exact))
# This should match Df up to floating point
Df_check = Df[:len(DR)]
rms_check = np.sqrt(np.mean((Df_check - Df_exact)**2))
record("F-DISC", f"Df(x) vs [prime]-DR(x) consistency check (x<=2000): rms={rms_check:.2e}",
       mpf(rms_check))

# Autocorrelation of Df
acf = np.correlate(Df - np.mean(Df), Df - np.mean(Df), mode='full')
acf = acf[len(acf)//2:]  # positive lags
acf = acf / acf[0]  # normalize
print(f"\n  Autocorrelation of Df at lags 1-10: {[f'{acf[i]:.4f}' for i in range(1, 11)]}")

# PSLQ on Df at specific x
print("\n  PSLQ on Df(x) at selected x:")
for x in [100, 500, 1000, 5000]:
    xm = mpf(x)
    dfx = f_arr[x+1] - f_arr[x]
    is_p = mpf(1) if pi_arr[x+1] > pi_arr[x] else mpf(0)

    vec = [dfx, is_p, mpf(1)/log(xm), mpf(1)/log(xm)**2,
           mpf(1)/sqrt(xm), sin(gamma1_zero * log(xm)) / sqrt(xm)]
    labels_df = ["Df", "prime?", "1/log(x)", "1/log^2(x)",
                 "1/sqrt(x)", "sin(g1*logx)/sqrt(x)"]

    # Skip if any element is zero (PSLQ requirement)
    if any(v == 0 for v in vec):
        record("F-DISC", f"PSLQ Df x={x}: skipped (zero element)", mpf(1))
        continue

    rel = pslq(vec, maxcoeff=10**6, maxsteps=1000)
    if rel is not None and rel[0] != 0:
        residual = fabs(fsum(mpf(r)*v for r, v in zip(rel, vec)))
        terms = " + ".join(f"{r}*{l}" for r, l in zip(rel, labels_df) if r != 0)
        record("F-DISC", f"PSLQ Df x={x}: {terms}", residual)
    else:
        record("F-DISC", f"PSLQ Df x={x}: no relation", mpf(1))

print()

# ─────────────────────────────────────────────────────────────────────
# 2g. BONUS: Multi-point PSLQ -- test if f(x) at DIFFERENT x are related
# ─────────────────────────────────────────────────────────────────────
print("=" * 70)
print("SECTION G: Multi-point relations")
print("=" * 70)

# Are f(100), f(200), f(500), f(1000) rationally related?
multi_pts = [100, 200, 300, 500, 1000, 2000, 5000]
vec_multi = [f_arr[x] for x in multi_pts] + [mpf(1)]
labels_multi = [f"f({x})" for x in multi_pts] + ["1"]

rel = pslq(vec_multi, maxcoeff=10**8, maxsteps=3000)
if rel is not None:
    residual = fabs(fsum(mpf(r)*v for r, v in zip(rel, vec_multi)))
    terms = " + ".join(f"{r}*{l}" for r, l in zip(rel, labels_multi) if r != 0)
    record("G-MULTI", f"Multi-point: {terms}", residual)
else:
    record("G-MULTI", "No integer relation among f at selected points", mpf(1))

# Normalized: f(x)/sqrt(x) at multiple points
vec_norm = [f_arr[x] / sqrt(mpf(x)) for x in multi_pts] + [mpf(1)]
labels_norm = [f"f({x})/sqrt({x})" for x in multi_pts] + ["1"]

rel = pslq(vec_norm, maxcoeff=10**8, maxsteps=3000)
if rel is not None:
    residual = fabs(fsum(mpf(r)*v for r, v in zip(rel, vec_norm)))
    terms = " + ".join(f"{r}*{l}" for r, l in zip(rel, labels_norm) if r != 0)
    record("G-MULTI", f"Normalized multi-point: {terms}", residual)
else:
    record("G-MULTI", "No relation among f(x)/sqrt(x)", mpf(1))

print()

# ─────────────────────────────────────────────────────────────────────
# 2h. BONUS: Fourier analysis of f(x) -- look for dominant frequencies
# ─────────────────────────────────────────────────────────────────────
print("=" * 70)
print("SECTION H: Fourier analysis of f(x)")
print("=" * 70)

# f(x) in log-space: f(e^t) should show peaks at imaginary parts of zeros
# We sample f at equally spaced points and FFT

from numpy.fft import fft

# Sample f at x = 2..XMAX
f_signal = f_float.copy()
N = len(f_signal)
F = fft(f_signal)
freqs = np.fft.fftfreq(N, d=1.0)
magnitudes = np.abs(F[:N//2])
freqs_pos = freqs[:N//2]

# Top 10 frequencies
top_idx = np.argsort(magnitudes[1:])[-10:][::-1] + 1  # skip DC
print("  Top 10 Fourier peaks (in x-space):")
for idx in top_idx:
    period = 1.0 / freqs_pos[idx] if freqs_pos[idx] > 0 else float('inf')
    print(f"    freq={freqs_pos[idx]:.6f}, period={period:.2f}, magnitude={magnitudes[idx]:.4f}")

# Now in log-space: resample f(e^t) for t = log(2)..log(XMAX)
t_vals = np.linspace(math.log(2), math.log(XMAX), 4096)
f_log = np.interp(np.exp(t_vals), np.arange(2, XMAX + 1), f_float)
F_log = fft(f_log)
N_log = len(f_log)
dt = t_vals[1] - t_vals[0]
freqs_log = np.fft.fftfreq(N_log, d=dt)
mag_log = np.abs(F_log[:N_log//2])
freqs_log_pos = freqs_log[:N_log//2]

top_log_idx = np.argsort(mag_log[1:])[-10:][::-1] + 1
print("\n  Top 10 Fourier peaks in log-space (should show zeta zeros):")
print(f"  (First zeta zero gamma_1 = {float(gamma1_zero):.6f})")
for idx in top_log_idx:
    omega = freqs_log_pos[idx] * 2 * math.pi  # angular frequency
    print(f"    omega={omega:.6f} (freq={freqs_log_pos[idx]:.6f}), magnitude={mag_log[idx]:.4f}")

print()

# ─────────────────────────────────────────────────────────────────────
# SUMMARY
# ─────────────────────────────────────────────────────────────────────
print("=" * 70)
print("SUMMARY OF ALL RESULTS")
print("=" * 70)

results.sort(key=lambda r: r[2])

print(f"\n{'Rank':<5} {'Section':<12} {'Residual':<15} {'Description'}")
print("-" * 90)
for i, (sec, desc, res, detail) in enumerate(results):
    tag = " <<<INTERESTING>>>" if res < 1e-10 else ""
    print(f"{i+1:<5} {sec:<12} {res:<15.6e} {desc}{tag}")

# Check for any truly interesting results
interesting = [r for r in results if r[2] < 1e-10]
if interesting:
    print(f"\n*** FOUND {len(interesting)} POTENTIALLY INTERESTING RELATION(S)! ***")
    for sec, desc, res, detail in interesting:
        print(f"\n  {sec}: {desc}")
        print(f"  Residual: {res:.2e}")
        print(f"  Detail: {detail}")
else:
    print("\nNo relations with residual < 1e-10 found.")
    print("This is consistent with the barrier: f(x) encodes zeta zero phases")
    print("that are information-theoretically incompressible.")

# Save results to file
with open("experiments/algebraic/pslq_identity_results.txt", "w") as fout:
    fout.write("PSLQ / LLL Identity Search Results\n")
    fout.write(f"Date: 2026-04-04\n")
    fout.write(f"Range: x = 2..{XMAX}\n")
    fout.write(f"f(x) = pi(x) - R(x), R(x) = sum mu(k)/k * li(x^(1/k))\n\n")

    fout.write(f"f(x) range: [{min(f_vals):.4f}, {max(f_vals):.4f}]\n")
    fout.write(f"f(x) mean:  {sum(f_vals)/len(f_vals):.6f}\n\n")

    fout.write("ALL RESULTS (sorted by residual):\n")
    fout.write(f"{'Rank':<5} {'Section':<12} {'Residual':<15} {'Description'}\n")
    fout.write("-" * 90 + "\n")
    for i, (sec, desc, res, detail) in enumerate(results):
        tag = " <<<INTERESTING>>>" if res < 1e-10 else ""
        fout.write(f"{i+1:<5} {sec:<12} {res:<15.6e} {desc}{tag}\n")

    fout.write("\n\nDETAILED RESULTS:\n\n")
    for sec, desc, res, detail in results:
        fout.write(f"[{sec}] {desc}\n")
        fout.write(f"  Residual: {res:.6e}\n")
        if detail:
            fout.write(f"  Detail: {detail}\n")
        fout.write("\n")

    if interesting:
        fout.write(f"\n*** {len(interesting)} INTERESTING RELATION(S) ***\n")
        for sec, desc, res, detail in interesting:
            fout.write(f"\n  {sec}: {desc}\n  Residual: {res:.2e}\n  Detail: {detail}\n")
    else:
        fout.write("\nCONCLUSION: No algebraic/linear/recurrence identity found for f(x).\n")
        fout.write("f(x) = pi(x) - R(x) appears to have no simple closed-form expression.\n")
        fout.write("This supports the information-theoretic barrier: the oscillatory\n")
        fout.write("residual encodes ~10^48 zeta zero phases that are incompressible.\n")

print("\nResults saved to experiments/algebraic/pslq_identity_results.txt")
print("Done.")

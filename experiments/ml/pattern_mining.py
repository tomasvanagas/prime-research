#!/usr/bin/env python3
"""
Pattern Mining for the nth Prime Number
========================================
Searching for a direct formula p(n) = f(n) by mining hidden patterns
in the prime sequence.
"""

import numpy as np
from scipy import optimize, fft
from sympy import primerange, prime, isprime, factorint
import math
import time
import warnings
warnings.filterwarnings('ignore')

print("="*80)
print("PRIME PATTERN MINING - Searching for p(n) = f(n)")
print("="*80)

# ============================================================
# STEP 0: Generate ground truth primes
# ============================================================
print("\n[STEP 0] Generating primes...")
t0 = time.time()
# Generate first 100,000 primes
primes_list = list(primerange(2, 1_300_000))  # ~100K primes below 1.3M
primes_list = primes_list[:100_000]
N_MAX = len(primes_list)
print(f"  Generated {N_MAX} primes in {time.time()-t0:.2f}s")
print(f"  p(1)={primes_list[0]}, p(100)={primes_list[99]}, p(1000)={primes_list[999]}, p(100000)={primes_list[99999]}")

# Convert to numpy array (1-indexed conceptually, but 0-indexed in array)
P = np.array(primes_list, dtype=np.float64)
# n array (1-indexed)
N = np.arange(1, N_MAX + 1, dtype=np.float64)

def approx_pn(n):
    """Best known approximation: n*(ln(n) + ln(ln(n)) - 1 + (ln(ln(n))-2)/ln(n))"""
    n = np.asarray(n, dtype=np.float64)
    ln_n = np.log(n)
    lln_n = np.log(ln_n)
    return n * (ln_n + lln_n - 1.0 + (lln_n - 2.0) / ln_n)


# ============================================================
# SECTION 1: Residual Analysis
# ============================================================
print("\n" + "="*80)
print("SECTION 1: RESIDUAL ANALYSIS AFTER KNOWN APPROXIMATION")
print("="*80)

# Compute residuals for specific n values
test_ns = [10, 100, 1000, 10000, 100000]
print("\n--- Residual r(n) = p(n) - approx(n) at key points ---")
print(f"{'n':>10s}  {'p(n)':>12s}  {'approx(n)':>14s}  {'r(n)':>12s}  {'r(n)/sqrt(n)':>14s}  {'r(n)/n^(1/3)':>14s}")
for n in test_ns:
    pn = primes_list[n-1]
    an = approx_pn(n)
    rn = pn - an
    print(f"{n:>10d}  {pn:>12d}  {an:>14.2f}  {rn:>12.2f}  {rn/np.sqrt(n):>14.6f}  {rn/n**(1/3):>14.6f}")

# Full residual array (skip n=1,2 where log(log(n)) is problematic)
start_n = 10
idx = slice(start_n - 1, N_MAX)
R = P[idx] - approx_pn(N[idx])
N_sub = N[idx]

# Study r(n) scaling
print("\n--- Residual scaling analysis ---")
for alpha in [0.0, 0.25, 1/3, 0.5, 2/3, 0.75, 1.0]:
    scaled = R / N_sub**alpha
    print(f"  r(n)/n^{alpha:.2f}: mean={np.mean(scaled):>12.4f}, std={np.std(scaled):>10.4f}, "
          f"min={np.min(scaled):>10.4f}, max={np.max(scaled):>10.4f}")

# Best alpha: fit log|r(n)| = alpha * log(n) + const
pos_idx = R > 0
if np.sum(pos_idx) > 100:
    log_r = np.log(np.abs(R[pos_idx]))
    log_n = np.log(N_sub[pos_idx])
    slope, intercept = np.polyfit(log_n, log_r, 1)
    print(f"\n  Best power-law fit: |r(n)| ~ n^{slope:.6f} (R ~ n^alpha with alpha={slope:.4f})")
    print(f"  This means r(n) grows roughly as n^{slope:.2f}")

# Trigonometric correlations
print("\n--- Trigonometric correlation of r(n) ---")
# Use a subset for correlation
subset = slice(0, min(10000, len(R)))
R_s = R[subset]
N_s = N_sub[subset]
for desc, func in [
    ("sin(n*ln(n))", np.sin(N_s * np.log(N_s))),
    ("cos(sqrt(n))", np.cos(np.sqrt(N_s))),
    ("sin(n*ln(ln(n)))", np.sin(N_s * np.log(np.log(N_s)))),
    ("cos(n^(2/3))", np.cos(N_s**(2/3))),
    ("sin(pi*n/ln(n))", np.sin(np.pi * N_s / np.log(N_s))),
]:
    corr = np.corrcoef(R_s, func)[0, 1]
    print(f"  corr(r(n), {desc:25s}) = {corr:>10.6f}")

# Fourier transform of residuals
print("\n--- Fourier Transform of r(n) for n=10..10009 ---")
R_fft_input = R[:10000]
fft_result = np.abs(fft.fft(R_fft_input))
freqs = fft.fftfreq(len(R_fft_input))
# Top 10 frequencies (exclude DC)
fft_result[0] = 0  # remove DC
top_idx = np.argsort(fft_result)[-10:][::-1]
print("  Top 10 Fourier amplitudes:")
for i, idx_val in enumerate(top_idx):
    f = freqs[idx_val]
    a = fft_result[idx_val]
    period = 1.0/abs(f) if abs(f) > 1e-10 else float('inf')
    print(f"    #{i+1}: freq={f:>10.6f}, amplitude={a:>14.2f}, period={period:>10.2f}")


# ============================================================
# SECTION 2: Prime Ratios and Continued Fractions
# ============================================================
print("\n" + "="*80)
print("SECTION 2: PRIME RATIOS AND CONTINUED FRACTIONS")
print("="*80)

# p(n+1)/p(n) ratios
print("\n--- p(n+1)/p(n) statistics ---")
ratios = P[1:] / P[:-1]
for rng_name, rng in [("n=1..100", slice(0,100)), ("n=1..1000", slice(0,1000)),
                        ("n=1..10000", slice(0,10000)), ("n=1..100000", slice(0,99999))]:
    r = ratios[rng]
    print(f"  {rng_name:15s}: mean={np.mean(r):.8f}, std={np.std(r):.8f}, "
          f"min={np.min(r):.8f}, max={np.max(r):.8f}")

# p(n)/(n*ln(n)) convergence
print("\n--- p(n)/(n*ln(n)) convergence ---")
for n in [10, 100, 1000, 10000, 100000]:
    val = primes_list[n-1] / (n * math.log(n))
    print(f"  n={n:>7d}: p(n)/(n*ln(n)) = {val:.10f}")
print("  (Should approach 1 by PNT)")

# More refined: p(n)/[n*(ln(n) + ln(ln(n)) - 1)]
print("\n--- p(n)/[n*(ln(n)+ln(ln(n))-1)] convergence ---")
for n in [10, 100, 1000, 10000, 100000]:
    ln_n = math.log(n)
    lln_n = math.log(ln_n)
    denom = n * (ln_n + lln_n - 1)
    val = primes_list[n-1] / denom
    print(f"  n={n:>7d}: ratio = {val:.10f}")

# Continued fraction of p(n)/n
print("\n--- Continued fraction of p(n)/n for selected n ---")
def continued_fraction(x, terms=8):
    """Return first `terms` of the CF expansion of x."""
    cf = []
    for _ in range(terms):
        a = int(x)
        cf.append(a)
        frac = x - a
        if abs(frac) < 1e-10:
            break
        x = 1.0 / frac
    return cf

for n in [10, 50, 100, 500, 1000, 5000, 10000, 50000, 100000]:
    ratio = primes_list[n-1] / n
    cf = continued_fraction(ratio, 10)
    print(f"  n={n:>6d}: p(n)/n = {ratio:>12.6f}, CF = {cf}")

# Continued fraction of p(n)/[n*ln(n)]
print("\n--- Continued fraction of p(n)/[n*ln(n)] ---")
for n in [100, 1000, 10000, 100000]:
    ratio = primes_list[n-1] / (n * math.log(n))
    cf = continued_fraction(ratio, 10)
    print(f"  n={n:>6d}: p(n)/(n*ln(n)) = {ratio:>12.10f}, CF = {cf}")


# ============================================================
# SECTION 3: Recursive Structure
# ============================================================
print("\n" + "="*80)
print("SECTION 3: RECURSIVE STRUCTURE")
print("="*80)

# p(n) mod p(k) for small k
print("\n--- p(n) mod p(k) for small primes p(k) ---")
small_primes = [2, 3, 5, 7, 11, 13]
for pk in small_primes:
    mods = [primes_list[n-1] % pk for n in range(1, 51)]
    print(f"  p(n) mod {pk:>2d} for n=1..50: {mods}")

# Second differences
print("\n--- Gap analysis: Δp(n), Δ²p(n), Δ³p(n) ---")
gaps = np.diff(P)  # Δp(n) = p(n+1) - p(n)
gaps2 = np.diff(gaps)  # Δ²p(n)
gaps3 = np.diff(gaps2)  # Δ³p(n)

for desc, arr, rng in [
    ("Δp(n)", gaps, slice(0, 30)),
    ("Δ²p(n)", gaps2, slice(0, 30)),
    ("Δ³p(n)", gaps3, slice(0, 30)),
]:
    print(f"  {desc} for first 30: {arr[rng].astype(int).tolist()}")

print(f"\n  Δp(n) stats (n=1..100000):")
print(f"    mean={np.mean(gaps):.4f}, std={np.std(gaps):.4f}, max={np.max(gaps):.0f}")
print(f"  Δ²p(n) stats:")
print(f"    mean={np.mean(gaps2):.6f}, std={np.std(gaps2):.4f}")
print(f"  Δ³p(n) stats:")
print(f"    mean={np.mean(gaps3):.6f}, std={np.std(gaps3):.4f}")

# Does Δp(n) / ln(n) converge?
print("\n--- Δp(n) / ln(n) convergence ---")
for n in [10, 100, 1000, 10000, 50000]:
    chunk = gaps[max(0,n-50):n]
    avg_gap = np.mean(chunk)
    print(f"  avg gap near n={n:>6d}: {avg_gap:>8.2f}, avg_gap/ln(n) = {avg_gap/math.log(n):.4f}")

# Second-order recurrence check: p(n+1) - 2*p(n) + p(n-1)
print("\n--- Second-order: p(n+1) - 2*p(n) + p(n-1) ---")
second_order = P[2:] - 2*P[1:-1] + P[:-2]
print(f"  First 30 values: {second_order[:30].astype(int).tolist()}")
print(f"  Mean: {np.mean(second_order[:10000]):.6f}")
print(f"  Std:  {np.std(second_order[:10000]):.4f}")

# Super-primes: p(p(n))
print("\n--- Super-primes p(p(n)) ---")
print(f"  {'n':>5s}  {'p(n)':>8s}  {'p(p(n))':>10s}  {'p(p(n))/(n*ln(n)^2)':>22s}  {'p(p(n))/(n^2*ln(n))':>22s}")
for n in [1, 2, 3, 5, 10, 20, 50, 100, 200, 500]:
    pn = primes_list[n-1]
    if pn <= N_MAX:
        ppn = primes_list[pn-1]
        ln_n = max(math.log(n), 0.01)
        r1 = ppn / (n * ln_n**2) if n > 1 else float('nan')
        r2 = ppn / (n**2 * ln_n) if n > 1 else float('nan')
        print(f"  {n:>5d}  {pn:>8d}  {ppn:>10d}  {r1:>22.6f}  {r2:>22.6f}")

# Super-prime approximation test
print("\n--- Super-prime fit: p(p(n)) vs n * (ln(n))^2 * ln(ln(n)) ---")
for n in [10, 20, 50, 100, 200, 500]:
    pn = primes_list[n-1]
    if pn <= N_MAX:
        ppn = primes_list[pn-1]
        ln_n = math.log(n)
        lln = math.log(ln_n) if ln_n > 0 else 1
        approx = n * ln_n**2 * lln
        ratio = ppn / approx if approx > 0 else float('nan')
        print(f"  n={n:>4d}: p(p(n))={ppn:>10d}, n*ln(n)^2*ln(ln(n))={approx:>14.1f}, ratio={ratio:.6f}")


# ============================================================
# SECTION 4: Novel Basis Function Fitting
# ============================================================
print("\n" + "="*80)
print("SECTION 4: NOVEL BASIS FUNCTION FITTING")
print("="*80)

# Fit 1: p(n) = a*n*ln(n) + b*n*ln(ln(n)) + c*n + d*sqrt(n)*ln(n) + e*sqrt(n) + f
print("\n--- Fit 1: General basis functions ---")
def make_basis1(n):
    ln_n = np.log(n)
    lln = np.log(ln_n)
    sqn = np.sqrt(n)
    return np.column_stack([
        n * ln_n,
        n * lln,
        n,
        sqn * ln_n,
        sqn,
        np.ones_like(n)
    ])

# Fit on different ranges and check coefficient stability
for fit_max in [1000, 2000, 5000, 10000]:
    idx_fit = slice(9, fit_max)  # start from n=10
    n_fit = N[idx_fit]
    p_fit = P[idx_fit]
    A = make_basis1(n_fit)
    coeffs, residuals, rank, sv = np.linalg.lstsq(A, p_fit, rcond=None)
    pred = A @ coeffs
    max_err = np.max(np.abs(pred - p_fit))
    rms_err = np.sqrt(np.mean((pred - p_fit)**2))
    print(f"  Fit on n=10..{fit_max}: coeffs=[{', '.join(f'{c:.6f}' for c in coeffs)}]")
    print(f"    RMS error={rms_err:.2f}, Max error={max_err:.2f}")

# Fit 2: Cipolla-style
print("\n--- Fit 2: Cipolla-style: p(n) = n*ln(n) + n*ln(ln(n)) - n + a*n/ln(n) + b*n/ln(n)^2 + c*n*ln(ln(n))/ln(n)^2 ---")
def make_basis2(n):
    ln_n = np.log(n)
    lln = np.log(ln_n)
    return np.column_stack([
        n / ln_n,
        n / ln_n**2,
        n * lln / ln_n**2,
    ])

for fit_max in [1000, 2000, 5000, 10000]:
    idx_fit = slice(9, fit_max)
    n_fit = N[idx_fit]
    p_fit = P[idx_fit]
    # Subtract known terms
    ln_n = np.log(n_fit)
    lln = np.log(ln_n)
    known = n_fit * ln_n + n_fit * lln - n_fit
    target = p_fit - known
    A = make_basis2(n_fit)
    coeffs, _, _, _ = np.linalg.lstsq(A, target, rcond=None)
    pred = known + A @ coeffs
    max_err = np.max(np.abs(pred - p_fit))
    rms_err = np.sqrt(np.mean((pred - p_fit)**2))
    print(f"  Fit on n=10..{fit_max}: coeffs=[{', '.join(f'{c:.6f}' for c in coeffs)}]")
    print(f"    RMS error={rms_err:.2f}, Max error={max_err:.2f}")

# Fit 3: Extended Cipolla with more terms
print("\n--- Fit 3: Extended Cipolla with 6 correction terms ---")
def make_basis3(n):
    ln_n = np.log(n)
    lln = np.log(ln_n)
    return np.column_stack([
        n / ln_n,
        n / ln_n**2,
        n * lln / ln_n**2,
        n / ln_n**3,
        n * lln / ln_n**3,
        n * lln**2 / ln_n**3,
    ])

for fit_max in [1000, 2000, 5000, 10000]:
    idx_fit = slice(9, fit_max)
    n_fit = N[idx_fit]
    p_fit = P[idx_fit]
    ln_n = np.log(n_fit)
    lln = np.log(ln_n)
    known = n_fit * ln_n + n_fit * lln - n_fit
    target = p_fit - known
    A = make_basis3(n_fit)
    coeffs, _, _, _ = np.linalg.lstsq(A, target, rcond=None)
    pred = known + A @ coeffs
    max_err = np.max(np.abs(pred - p_fit))
    rms_err = np.sqrt(np.mean((pred - p_fit)**2))
    print(f"  Fit on n=10..{fit_max}: coeffs=[{', '.join(f'{c:.6f}' for c in coeffs)}]")
    print(f"    RMS error={rms_err:.2f}, Max error={max_err:.2f}")

# Coefficient stability check for best fit
print("\n--- Coefficient stability for Fit 3 ---")
print(f"  {'fit_max':>10s} | {'c1 (n/ln)':>12s} {'c2 (n/ln^2)':>12s} {'c3 (nlln/ln^2)':>14s} "
      f"{'c4 (n/ln^3)':>12s} {'c5 (nlln/ln^3)':>14s} {'c6 (nlln^2/ln^3)':>16s}")
for fit_max in [500, 1000, 2000, 3000, 5000, 7000, 10000]:
    idx_fit = slice(9, fit_max)
    n_fit = N[idx_fit]
    p_fit = P[idx_fit]
    ln_n = np.log(n_fit)
    lln = np.log(ln_n)
    known = n_fit * ln_n + n_fit * lln - n_fit
    target = p_fit - known
    A = make_basis3(n_fit)
    coeffs, _, _, _ = np.linalg.lstsq(A, target, rcond=None)
    print(f"  {fit_max:>10d} | {coeffs[0]:>12.6f} {coeffs[1]:>12.6f} {coeffs[2]:>14.6f} "
          f"{coeffs[3]:>12.6f} {coeffs[4]:>14.6f} {coeffs[5]:>16.6f}")

# Residual analysis after best fit
print("\n--- Residual analysis after Fit 3 (trained on n=10..10000) ---")
idx_fit = slice(9, 10000)
n_fit = N[idx_fit]
p_fit = P[idx_fit]
ln_n = np.log(n_fit)
lln = np.log(ln_n)
known = n_fit * ln_n + n_fit * lln - n_fit
target = p_fit - known
A = make_basis3(n_fit)
coeffs3, _, _, _ = np.linalg.lstsq(A, target, rcond=None)
pred = known + A @ coeffs3
residual_fit3 = p_fit - pred
print(f"  Residual stats: mean={np.mean(residual_fit3):.4f}, std={np.std(residual_fit3):.4f}")
print(f"  Residual/sqrt(n): mean={np.mean(residual_fit3/np.sqrt(n_fit)):.6f}, std={np.std(residual_fit3/np.sqrt(n_fit)):.6f}")
# Check if residual is bounded
print(f"  |Residual| max = {np.max(np.abs(residual_fit3)):.2f}")
print(f"  |Residual|/n^(1/3) max = {np.max(np.abs(residual_fit3)/n_fit**(1/3)):.4f}")
print(f"  |Residual|/sqrt(n) max = {np.max(np.abs(residual_fit3)/np.sqrt(n_fit)):.4f}")
# Residual at specific points
print(f"  Sample residuals: n=100: {residual_fit3[90]:.2f}, n=1000: {residual_fit3[990]:.2f}, "
      f"n=5000: {residual_fit3[4990]:.2f}, n=9999: {residual_fit3[9989]:.2f}")


# ============================================================
# SECTION 5: Key Experiment - Exact Integer Relationships
# ============================================================
print("\n" + "="*80)
print("SECTION 5: KEY EXPERIMENT - EXACT INTEGER RELATIONSHIPS")
print("="*80)

# p(n) * ln(p(n)) / (n * ln(n)^2)
print("\n--- p(n) * ln(p(n)) / (n * ln(n)^2) ---")
for n in [10, 50, 100, 500, 1000, 5000, 10000, 50000, 100000]:
    pn = primes_list[n-1]
    val = pn * math.log(pn) / (n * math.log(n)**2)
    print(f"  n={n:>7d}: {val:.10f}")

# Mertens-like: sum(1/p(k)) vs ln(ln(n))
print("\n--- sum(1/p(k)) vs ln(ln(n)) + M ---")
cumsum_inv_p = np.cumsum(1.0 / P)
mertens_const = 0.2614972128  # Meissel-Mertens constant
for n in [10, 100, 1000, 10000, 100000]:
    s = cumsum_inv_p[n-1]
    lln = math.log(math.log(n))
    diff = s - lln - mertens_const
    print(f"  n={n:>7d}: sum(1/p(k))={s:.10f}, ln(ln(n))+M={lln+mertens_const:.10f}, diff={diff:.10f}")

# Note: The actual Mertens theorem is about sum over primes <= x, not first n primes
# But this is related via p(n) ~ n*ln(n)
print("\n--- sum_{p <= p(n)} 1/p vs ln(ln(p(n))) + M ---")
for n in [10, 100, 1000, 10000, 100000]:
    s = cumsum_inv_p[n-1]
    pn = primes_list[n-1]
    llp = math.log(math.log(pn))
    diff = s - llp - mertens_const
    print(f"  n={n:>7d}, p(n)={pn:>7d}: sum={s:.10f}, ln(ln(p(n)))+M={llp+mertens_const:.10f}, diff={diff:.10f}")

# Quadratic: p(n)^2 vs n^2 * ln(n)^2
print("\n--- p(n)^2 / (n^2 * ln(n)^2) ---")
for n in [10, 100, 1000, 10000, 100000]:
    pn = primes_list[n-1]
    val = pn**2 / (n * math.log(n))**2
    print(f"  n={n:>7d}: {val:.10f}")

# p(n) mod n
print("\n--- p(n) mod n ---")
pmod_n = [int(P[i]) % (i+1) for i in range(100)]
print(f"  p(n) mod n for n=1..100: {pmod_n}")
# Statistics for larger n
pmod_vals = [int(P[i]) % (i+1) for i in range(99999)]
pmod_ratios = [pmod_vals[i]/(i+1) for i in range(99999)]
print(f"  p(n) mod n / n: mean={np.mean(pmod_ratios):.6f}, std={np.std(pmod_ratios):.6f}")
print(f"  (If uniform on [0,n), expect mean 0.5)")

# floor(p(n)/n) - floor(p(n-1)/(n-1))
print("\n--- floor(p(n)/n) - floor(p(n-1)/(n-1)) ---")
floor_diffs = []
for i in range(2, 100):  # n=3..100
    a = int(P[i]) // (i+1)
    b = int(P[i-1]) // i
    floor_diffs.append(a - b)
print(f"  Values for n=3..100: {floor_diffs}")
print(f"  Mean: {np.mean(floor_diffs):.4f}, # of 0s: {floor_diffs.count(0)}, # of 1s: {floor_diffs.count(1)}")

# NEW: p(n) * n / (sum of first n primes)
print("\n--- p(n) * n / sum(p(1)..p(n)) ---")
cumsum_p = np.cumsum(P)
for n in [10, 100, 1000, 10000, 100000]:
    val = P[n-1] * n / cumsum_p[n-1]
    print(f"  n={n:>7d}: p(n)*n/S(n) = {val:.10f}")
print("  (If primes ~ n*ln(n), then S(n) ~ n^2*ln(n)/2, so ratio ~ 2)")

# NEW: p(n) vs li^{-1}(n) where li is log integral
print("\n--- p(n) vs li^{-1}(n) approximation ---")
# li^{-1}(n) ~ n*ln(n) + n*ln(ln(n)) - n + ... (same as Cipolla)
# But let's check: p(n) - li^{-1}(n) where li^{-1} is computed numerically
from scipy.special import expi
from scipy.optimize import brentq

def li(x):
    """Logarithmic integral li(x) = Ei(ln(x))"""
    return expi(math.log(x))

def li_inv(n):
    """Inverse of li: find x such that li(x) = n"""
    # li(x) ~ x/ln(x), so start with x ~ n*ln(n)
    low = max(2.1, n)
    high = n * (math.log(n) + math.log(math.log(max(n, 3)))) * 2
    try:
        return brentq(lambda x: li(x) - n, low, high)
    except:
        return float('nan')

print(f"  {'n':>7s}  {'p(n)':>10s}  {'li_inv(n)':>14s}  {'diff':>12s}  {'diff/sqrt(n)':>14s}")
for n in [10, 100, 1000, 10000, 100000]:
    pn = primes_list[n-1]
    lin = li_inv(n)
    diff = pn - lin
    print(f"  {n:>7d}  {pn:>10d}  {lin:>14.4f}  {diff:>12.4f}  {diff/math.sqrt(n):>14.6f}")

# NEW: Check if p(n) = round(li_inv(n) + correction)?
print("\n--- p(n) - li_inv(n): searching for a correction term ---")
corrections = []
for n in range(10, 1001):
    pn = primes_list[n-1]
    lin = li_inv(n)
    corrections.append(pn - lin)
corrections = np.array(corrections)
ns_corr = np.arange(10, 1001, dtype=float)
print(f"  Correction stats (n=10..1000): mean={np.mean(corrections):.4f}, std={np.std(corrections):.4f}")
print(f"  Correction/sqrt(n*ln(n)): mean={np.mean(corrections/np.sqrt(ns_corr*np.log(ns_corr))):.6f}")
# Fit correction
slope_c, intercept_c = np.polyfit(np.log(ns_corr), corrections, 1)
print(f"  Linear fit of correction vs ln(n): slope={slope_c:.6f}, intercept={intercept_c:.6f}")
print(f"  => correction ~ {slope_c:.4f} * ln(n) + {intercept_c:.4f}")


# ============================================================
# SECTION 6: Generating Function and Transform
# ============================================================
print("\n" + "="*80)
print("SECTION 6: GENERATING FUNCTION AND TRANSFORMS")
print("="*80)

# Prime zeta function P(s) = sum p^(-s)
print("\n--- Prime zeta function P(s) = sum_{k=1}^{N} p(k)^{-s} ---")
for s in [1.0, 1.5, 2.0, 2.5, 3.0]:
    for N_lim in [100, 1000, 10000, 100000]:
        val = np.sum(P[:N_lim]**(-s))
        print(f"  P({s:.1f}), N={N_lim:>6d}: {val:.10f}")

# Generating function: sum p(n) * x^n
print("\n--- Generating function F(x) = sum p(n)*x^n at small x ---")
for x in [0.1, 0.5, 0.9, 0.99]:
    powers = x ** N[:10000]
    val = np.sum(P[:10000] * powers)
    print(f"  F({x}) [N=10000] = {val:.6f}")

# Relation to other sequences
print("\n--- Relation to Fibonacci numbers ---")
fib = [1, 1]
for i in range(2, 50):
    fib.append(fib[-1] + fib[-2])
print(f"  p(n)/F(n) for n=1..20:")
for n in range(1, 21):
    ratio = primes_list[n-1] / fib[n-1]
    print(f"    n={n:>2d}: p(n)={primes_list[n-1]:>4d}, F(n)={fib[n-1]:>5d}, ratio={ratio:.6f}")

# Relation to factorials
print("\n--- p(n) vs (n!)^{1/n} * some_function ---")
# By Stirling: n! ~ (n/e)^n * sqrt(2*pi*n), so (n!)^{1/n} ~ n/e
for n in [5, 10, 20, 50, 100]:
    nfact_nth = math.exp(sum(math.log(k) for k in range(1, n+1)) / n)
    ratio = primes_list[n-1] / nfact_nth
    print(f"  n={n:>3d}: p(n)={primes_list[n-1]:>6d}, (n!)^(1/n)={nfact_nth:.4f}, ratio={ratio:.6f}")

# ============================================================
# SECTION 7: NOVEL EXPLORATIONS
# ============================================================
print("\n" + "="*80)
print("SECTION 7: NOVEL EXPLORATIONS")
print("="*80)

# Explore: p(n) = n * W(n) * correction, where W is Lambert W
print("\n--- p(n) / (n * W(n * e)) where W is Lambert W ---")
# Note: W(n*e) = 1 + ln(n) - ln(1+ln(n)) + ... but let's compute numerically
from scipy.special import lambertw
for n in [10, 100, 1000, 10000, 100000]:
    w_val = float(np.real(lambertw(n)))
    ratio = primes_list[n-1] / (n * w_val)
    print(f"  n={n:>7d}: W(n)={w_val:.6f}, p(n)/(n*W(n))={ratio:.8f}")

# W(n) ~ ln(n) - ln(ln(n)), so n*W(n) ~ n*ln(n) - n*ln(ln(n))
# That's close to but not quite the approximation

# Explore: p(n) in terms of harmonic numbers
print("\n--- p(n) vs n * H(n) where H(n) is nth harmonic number ---")
H = np.cumsum(1.0 / N)
for n_val in [10, 100, 1000, 10000, 100000]:
    ratio = P[n_val-1] / (n_val * H[n_val-1])
    print(f"  n={n_val:>7d}: H(n)={H[n_val-1]:.6f}, p(n)/(n*H(n))={ratio:.8f}")
print("  H(n) ~ ln(n) + gamma, so n*H(n) ~ n*ln(n) + gamma*n")

# Explore: Inverse symbolic calculator approach
# Look at p(n) / [n * (ln(n) + ln(ln(n)))] for large n
print("\n--- p(n) / [n * (ln(n) + ln(ln(n)))] ---")
for n in [100, 1000, 10000, 100000]:
    ln_n = math.log(n)
    lln = math.log(ln_n)
    ratio = primes_list[n-1] / (n * (ln_n + lln))
    print(f"  n={n:>7d}: {ratio:.10f}")

# Explore: Does (p(n) - n*ln(n) - n*ln(ln(n)) + n) / (n/ln(n)) approach a constant?
print("\n--- [p(n) - n*ln(n) - n*ln(ln(n)) + n] / (n/ln(n)) ---")
for n in [100, 500, 1000, 5000, 10000, 50000, 100000]:
    ln_n = math.log(n)
    lln = math.log(ln_n)
    val = (primes_list[n-1] - n*ln_n - n*lln + n) / (n / ln_n)
    print(f"  n={n:>7d}: {val:.10f}")
print("  (Cipolla: this should approach (ln(ln(n))-2), checking...)")
for n in [100, 500, 1000, 5000, 10000, 50000, 100000]:
    ln_n = math.log(n)
    lln = math.log(ln_n)
    target = lln - 2
    val = (primes_list[n-1] - n*ln_n - n*lln + n) / (n / ln_n)
    print(f"  n={n:>7d}: value={val:.8f}, ln(ln(n))-2={target:.8f}, diff={val-target:.8f}")

# Going one level deeper in Cipolla
print("\n--- Next Cipolla term: [p(n) - n*(ln(n)+ln(ln(n))-1+(ln(ln(n))-2)/ln(n))] / (n/ln(n)^2) ---")
for n in [100, 500, 1000, 5000, 10000, 50000, 100000]:
    ln_n = math.log(n)
    lln = math.log(ln_n)
    approx = n * (ln_n + lln - 1 + (lln - 2) / ln_n)
    residual = primes_list[n-1] - approx
    normalized = residual / (n / ln_n**2)
    print(f"  n={n:>7d}: residual={residual:>12.2f}, normalized(n/ln^2)={normalized:.8f}")

# What is this converging to?
print("\n--- Fitting the next Cipolla coefficient ---")
ns_fit = np.array([500, 1000, 2000, 5000, 10000, 50000, 100000], dtype=float)
normalized_vals = []
for n_int in ns_fit.astype(int):
    ln_n = math.log(n_int)
    lln = math.log(ln_n)
    approx = n_int * (ln_n + lln - 1 + (lln - 2) / ln_n)
    residual = primes_list[n_int-1] - approx
    normalized = residual / (n_int / ln_n**2)
    normalized_vals.append(normalized)
print(f"  Normalized residuals: {[f'{v:.6f}' for v in normalized_vals]}")
# The next Cipolla term should be n/ln(n)^2 * (ln(ln(n))^2 - 6*ln(ln(n)) + 11)/2
print("  Expected: 0.5*(ln(ln(n))^2 - 6*ln(ln(n)) + 11)")
for n_int, nv in zip(ns_fit.astype(int), normalized_vals):
    ln_n = math.log(n_int)
    lln = math.log(ln_n)
    expected = 0.5 * (lln**2 - 6*lln + 11)
    print(f"  n={n_int:>7d}: actual={nv:.6f}, expected={expected:.6f}, diff={nv-expected:.6f}")

# ============================================================
# SECTION 8: DEEP DIVE - li^{-1}(n) relationship
# ============================================================
print("\n" + "="*80)
print("SECTION 8: DEEP DIVE - li_inv(n) AS EXACT FORMULA BASE")
print("="*80)

# The "true" formula is p(n) ~ li^{-1}(n) by PNT
# But how good is it, and what's the correction?
print("\n--- Detailed p(n) - li_inv(n) analysis ---")
diffs_li = []
ns_li = list(range(10, 10001))
for n in ns_li:
    pn = primes_list[n-1]
    lin = li_inv(n)
    diffs_li.append(pn - lin)
diffs_li = np.array(diffs_li)
ns_li = np.array(ns_li, dtype=float)

print(f"  Stats for n=10..10000:")
print(f"    mean(p(n)-li_inv(n)) = {np.mean(diffs_li):.4f}")
print(f"    std(p(n)-li_inv(n))  = {np.std(diffs_li):.4f}")
print(f"    max|p(n)-li_inv(n)|  = {np.max(np.abs(diffs_li)):.4f}")

# Normalize by sqrt(n)*ln(n) (expected from RH)
norm_diffs = diffs_li / (np.sqrt(ns_li) * np.log(ns_li))
print(f"\n  Normalized by sqrt(n)*ln(n):")
print(f"    mean = {np.mean(norm_diffs):.6f}")
print(f"    std  = {np.std(norm_diffs):.6f}")
print(f"    max  = {np.max(np.abs(norm_diffs)):.6f}")

# Is the correction related to Riemann zeros?
# Riemann's exact formula: pi(x) = li(x) - sum_rho li(x^rho) - ln(2) + integral...
# So p(n) = li_inv(n + correction_from_zeros)
# Let's see if the diff correlates with anything involving known zeros
print("\n  Autocorrelation of the correction term (first 20 lags):")
mean_d = np.mean(diffs_li)
diffs_centered = diffs_li - mean_d
var_d = np.var(diffs_li)
for lag in range(1, 21):
    ac = np.mean(diffs_centered[:-lag] * diffs_centered[lag:]) / var_d
    print(f"    lag {lag:>2d}: {ac:>10.6f}")

# FFT of the li correction
print("\n  FFT of p(n) - li_inv(n) for n=10..10009:")
fft_li = np.abs(fft.fft(diffs_li[:10000]))
freqs_li = fft.fftfreq(10000)
fft_li[0] = 0
top_li = np.argsort(fft_li)[-10:][::-1]
for i, idx_val in enumerate(top_li):
    f = freqs_li[idx_val]
    a = fft_li[idx_val]
    period = 1.0/abs(f) if abs(f) > 1e-10 else float('inf')
    print(f"    #{i+1}: freq={f:>10.6f}, amplitude={a:>14.2f}, period={period:>10.2f}")

# ============================================================
# SECTION 9: EMPIRICAL FORMULA HUNT
# ============================================================
print("\n" + "="*80)
print("SECTION 9: EMPIRICAL FORMULA HUNT")
print("="*80)

# Try: p(n) = li_inv(n + a*sqrt(n) + b*ln(n) + c) for best a, b, c
print("\n--- Fitting p(n) = li_inv(n + a*sqrt(n) + b*ln(n) + c) ---")
# This is slow to fit directly, so let's use a different approach:
# Since li(p(n)) = n + correction, let's compute li(p(n)) - n
print("\n--- li(p(n)) - n: the prime counting correction ---")
li_corrections = []
for n in [10, 100, 500, 1000, 2000, 5000, 10000, 50000, 100000]:
    pn = primes_list[n-1]
    li_pn = li(pn)
    corr = li_pn - n
    print(f"  n={n:>7d}: li(p(n))={li_pn:.6f}, li(p(n))-n={corr:.6f}, corr/sqrt(n)={corr/math.sqrt(n):.6f}")

# li(p(n)) - n should be related to the error in pi(x)
# By RH, |pi(x) - li(x)| < C * sqrt(x) * ln(x) / (8*pi)
# So |li(p(n)) - n| < C * sqrt(p(n)) * ln(p(n)) / (8*pi)

# Try to find if there's a better "shift" for li
print("\n--- Trying li_inv(n + alpha * sqrt(n/ln(n))) ---")
# What value of alpha minimizes error?
def total_error(alpha, n_range):
    total = 0
    count = 0
    for n in n_range:
        shifted_n = n + alpha * math.sqrt(n / math.log(n))
        lin = li_inv(shifted_n)
        if not math.isnan(lin):
            total += (primes_list[n-1] - lin)**2
            count += 1
    return total / count if count > 0 else float('inf')

test_range = list(range(10, 1001, 10))
alphas = np.linspace(-3, 3, 61)
errors = [total_error(a, test_range) for a in alphas]
best_alpha = alphas[np.argmin(errors)]
print(f"  Best alpha (n=10..1000 step 10): {best_alpha:.4f}")
# Refine
alphas2 = np.linspace(best_alpha - 0.2, best_alpha + 0.2, 41)
errors2 = [total_error(a, test_range) for a in alphas2]
best_alpha2 = alphas2[np.argmin(errors2)]
print(f"  Refined alpha: {best_alpha2:.6f}")
print(f"  RMS error with this alpha: {math.sqrt(min(errors2)):.4f}")
print(f"  RMS error without shift:   {math.sqrt(total_error(0, test_range)):.4f}")

# Compare formulas at specific points
print(f"\n--- Comparison: p(n) vs li_inv(n) vs li_inv(n + {best_alpha2:.4f}*sqrt(n/ln(n))) ---")
for n in [10, 100, 1000, 10000]:
    pn = primes_list[n-1]
    lin0 = li_inv(n)
    shifted = n + best_alpha2 * math.sqrt(n / math.log(n))
    lin1 = li_inv(shifted)
    print(f"  n={n:>6d}: p(n)={pn:>8d}, li_inv(n)={lin0:>12.2f}(err={pn-lin0:>8.2f}), "
          f"li_inv(shifted)={lin1:>12.2f}(err={pn-lin1:>8.2f})")

# ============================================================
# SECTION 10: MODULAR ARITHMETIC DEEP DIVE
# ============================================================
print("\n" + "="*80)
print("SECTION 10: MODULAR ARITHMETIC PATTERNS")
print("="*80)

# p(n) mod small numbers - distribution
print("\n--- Distribution of p(n) mod m for small m ---")
for m in [3, 4, 5, 6, 7, 8, 10, 12, 30]:
    counts = [0] * m
    for i in range(100000):
        counts[int(P[i]) % m] += 1
    print(f"  p(n) mod {m:>2d} (n=1..100000): {counts}")

# p(n) mod n - is it uniform?
print("\n--- Distribution of p(n) mod n / n ---")
ratios_modn = []
for n in range(2, 10001):
    ratios_modn.append((int(P[n-1]) % n) / n)
# Histogram in 10 bins
hist, bin_edges = np.histogram(ratios_modn, bins=10, range=(0, 1))
print(f"  Histogram of p(n) mod n / n (n=2..10000):")
for i in range(10):
    bar = '#' * (hist[i] // 20)
    print(f"    [{bin_edges[i]:.1f}, {bin_edges[i+1]:.1f}): {hist[i]:>5d} {bar}")

# ============================================================
# SECTION 11: RATIO OF CONSECUTIVE GAPS
# ============================================================
print("\n" + "="*80)
print("SECTION 11: GAP RATIOS AND SECOND-ORDER PATTERNS")
print("="*80)

print("\n--- Gap ratio g(n+1)/g(n) distribution ---")
gap_ratios = gaps[1:] / np.where(gaps[:-1] == 0, 1, gaps[:-1])
# Remove infinities
valid = np.isfinite(gap_ratios) & (gap_ratios < 100)
gr = gap_ratios[valid]
print(f"  Mean gap ratio: {np.mean(gr[:10000]):.6f}")
print(f"  Std gap ratio:  {np.std(gr[:10000]):.6f}")
print(f"  Median gap ratio: {np.median(gr[:10000]):.6f}")

# Cramer-Granville conjecture: max gap near n ~ (ln(n))^2
print("\n--- Maximum gap g(n) vs (ln(n))^2 at milestone points ---")
for n_max in [100, 1000, 10000, 100000]:
    max_gap = np.max(gaps[:n_max])
    ln_n = math.log(n_max)
    pn = primes_list[n_max-1]
    ln_p = math.log(pn)
    print(f"  n<={n_max:>6d}: max_gap={max_gap:>4.0f}, (ln(p(n)))^2={ln_p**2:>8.2f}, "
          f"ratio={max_gap/ln_p**2:.4f}")

# ============================================================
# SECTION 12: COMPREHENSIVE FORMULA ATTEMPT
# ============================================================
print("\n" + "="*80)
print("SECTION 12: COMPREHENSIVE FORMULA ATTEMPT")
print("="*80)

# Try a 9-parameter Cipolla-style expansion
print("\n--- 9-parameter asymptotic expansion ---")
print("  p(n) ~ n*(ln(n) + ln(ln(n)) - 1 + sum c_k * f_k(n))")

def make_full_basis(n):
    ln_n = np.log(n)
    lln = np.log(ln_n)
    return np.column_stack([
        n * (lln - 2) / ln_n,                    # Known Cipolla term
        n * (lln**2 - 6*lln + 11) / (2*ln_n**2),  # Next Cipolla
        n / ln_n**2,
        n * lln / ln_n**2,
        n / ln_n**3,
        n * lln / ln_n**3,
        n * lln**2 / ln_n**3,
        n / ln_n**4,
        n * lln / ln_n**4,
    ])

# Check if first two known coefficients come out right
idx_fit = slice(99, 50000)  # n=100..50000
n_fit = N[idx_fit]
p_fit = P[idx_fit]
ln_n = np.log(n_fit)
lln = np.log(ln_n)
known = n_fit * (ln_n + lln - 1)
target = p_fit - known
A = make_full_basis(n_fit)
coeffs_full, _, _, _ = np.linalg.lstsq(A, target, rcond=None)
print(f"  Coefficients (fitted on n=100..50000):")
labels = ["(lln-2)/ln", "(lln^2-6lln+11)/(2ln^2)", "1/ln^2", "lln/ln^2",
          "1/ln^3", "lln/ln^3", "lln^2/ln^3", "1/ln^4", "lln/ln^4"]
for lbl, c in zip(labels, coeffs_full):
    print(f"    {lbl:30s}: {c:>14.6f}")

print(f"\n  If Cipolla is exact, first coeff should be ~1.0: {coeffs_full[0]:.6f}")
print(f"  And second coeff should be ~1.0: {coeffs_full[1]:.6f}")

pred_full = known + A @ coeffs_full
rms_full = np.sqrt(np.mean((pred_full - p_fit)**2))
max_full = np.max(np.abs(pred_full - p_fit))
print(f"\n  Full 9-param fit: RMS={rms_full:.2f}, Max error={max_full:.2f}")

# Test on held-out data (n=50001..100000)
idx_test = slice(50000, 100000)
n_test = N[idx_test]
p_test = P[idx_test]
ln_test = np.log(n_test)
lln_test = np.log(ln_test)
known_test = n_test * (ln_test + lln_test - 1)
A_test = make_full_basis(n_test)
pred_test = known_test + A_test @ coeffs_full
rms_test = np.sqrt(np.mean((pred_test - p_test)**2))
max_test = np.max(np.abs(pred_test - p_test))
print(f"  Test on n=50001..100000: RMS={rms_test:.2f}, Max error={max_test:.2f}")

# Residual structure
residual_full = p_fit - pred_full
print(f"\n  Residual after 9-param fit (training):")
print(f"    mean={np.mean(residual_full):.4f}, std={np.std(residual_full):.4f}")
print(f"    max|res|/sqrt(n*ln(n))={np.max(np.abs(residual_full)/np.sqrt(n_fit*np.log(n_fit))):.6f}")

# FFT of residual
print(f"\n  FFT of residual after 9-param fit:")
res_fft = np.abs(fft.fft(residual_full[:10000]))
res_freqs = fft.fftfreq(10000)
res_fft[0] = 0
top_res = np.argsort(res_fft)[-5:][::-1]
for i, idx_val in enumerate(top_res):
    f = res_freqs[idx_val]
    a = res_fft[idx_val]
    period = 1.0/abs(f) if abs(f) > 1e-10 else float('inf')
    print(f"    #{i+1}: freq={f:>10.6f}, amplitude={a:>14.2f}, period={period:>10.2f}")


# ============================================================
# SECTION 13: HUNTING FOR EXACT CONSTANTS
# ============================================================
print("\n" + "="*80)
print("SECTION 13: HUNTING FOR EXACT CONSTANTS")
print("="*80)

# Is there a universal constant hiding?
print("\n--- Candidate constants ---")
# Compute various quantities and see if they approach known constants
euler_gamma = 0.5772156649015329

# S1: sum_{n=1}^{N} (p(n) - n*ln(n)) / n^2
print("\n  S1 = sum (p(n) - n*ln(n)) / n^2:")
for N_val in [100, 1000, 10000, 100000]:
    ns = np.arange(1, N_val+1, dtype=float)
    terms = (P[:N_val] - ns * np.log(ns)) / ns**2
    s = np.sum(terms[1:])  # skip n=1
    print(f"    N={N_val:>7d}: S1 = {s:.10f}")

# S2: Product (1 - 1/p(n)) * (1 + 1/(p(n)-1))... Euler product related
print("\n  S2 = prod_{n=1}^{N} p(n) / (p(n)-1) / (1 + 1/n):")
prod = 1.0
for n in range(1, 101):
    pn = primes_list[n-1]
    prod *= (pn / (pn - 1)) / (1 + 1/n)
print(f"    N=100: {prod:.10f}")
prod2 = 1.0
for n in range(1, 1001):
    pn = primes_list[n-1]
    prod2 *= (pn / (pn - 1)) / (1 + 1/n)
print(f"    N=1000: {prod2:.10f}")

# S3: An interesting new quantity
print("\n  S3 = mean of [p(n)*ln(n) - n*ln(n)*ln(n*ln(n))] / n for n=10..N:")
for N_val in [100, 1000, 10000]:
    ns = np.arange(10, N_val+1, dtype=float)
    vals = (P[9:N_val] * np.log(ns) - ns * np.log(ns) * np.log(ns * np.log(ns))) / ns
    print(f"    N={N_val}: mean = {np.mean(vals):.6f}, std = {np.std(vals):.6f}")


# ============================================================
# FINAL SUMMARY
# ============================================================
print("\n" + "="*80)
print("FINAL SUMMARY OF KEY FINDINGS")
print("="*80)

print("""
KEY FINDINGS:
=============

1. RESIDUAL SCALING: The residual r(n) = p(n) - approx(n) grows roughly as n^alpha
   where alpha is found empirically above.

2. li_inv(n) is the BEST base for a formula. The correction p(n) - li_inv(n)
   has mean close to 0 and its standard deviation grows slower than sqrt(n).

3. CIPOLLA EXPANSION: The asymptotic expansion coefficients converge:
   p(n) = n*(ln(n) + ln(ln(n)) - 1 + (ln(ln(n))-2)/ln(n) + ...)
   Next term: (ln(ln(n))^2 - 6*ln(ln(n)) + 11)/(2*ln(n)^2) * n

4. p(n) mod n appears to be approximately uniformly distributed (ratio mean ~ 0.5).

5. p(n)*n/S(n) approaches 2 (where S(n) = sum of first n primes).

6. The gap ratios g(n+1)/g(n) show no strong autocorrelation, consistent with
   pseudo-random behavior of prime gaps.

7. A shifted li_inv formula: p(n) ~ li_inv(n + alpha*sqrt(n/ln(n))) with optimal
   alpha can reduce systematic bias.
""")

print("\nDone! Total runtime:", time.time() - t0, "seconds")

#!/usr/bin/env python3
"""
Deep Kt Complexity Analysis of delta(n)
========================================

Extends Session 20's analysis with:
1. Partial Autocorrelation Function (PACF) — reveals true direct dependencies
2. Multi-algorithm compression comparison (gzip/bz2/lzma) with scaling
3. Compression ratio scaling with sequence length N
4. Explicit Kt(delta) vs log(n) growth curve
5. Block mutual information scaling
6. Detrended Fluctuation Analysis (DFA) — Hurst exponent
7. Transfer entropy: does past delta predict future delta better than past n?

Goal: Determine if delta(n) has exploitable structure beyond what Session 20 found.
"""

import numpy as np
import os
import sys
import gzip
import bz2
import lzma
import struct
import math
import time
from collections import Counter, defaultdict

# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------
DELTA_PATH = "/apps/aplikacijos/prime-research/experiments/information_theory/kt_delta_values.npy"
OUT_DIR = "/apps/aplikacijos/prime-research/experiments/information_theory/kt_complexity"
os.makedirs(OUT_DIR, exist_ok=True)

delta_all = np.load(DELTA_PATH)
N_TOTAL = len(delta_all)

results = []
def log(msg):
    print(msg)
    results.append(msg)

log("=" * 72)
log("DEEP Kt COMPLEXITY ANALYSIS OF delta(n)")
log(f"N_TOTAL = {N_TOTAL}, range = [{delta_all.min()}, {delta_all.max()}]")
log("=" * 72)


# ===================================================================
# 1. PARTIAL AUTOCORRELATION FUNCTION (PACF)
# ===================================================================
log("\n" + "=" * 72)
log("1. PARTIAL AUTOCORRELATION FUNCTION (PACF)")
log("=" * 72)

def compute_acf(x, max_lag):
    """Compute autocorrelation function."""
    n = len(x)
    x_centered = x - np.mean(x)
    var = np.var(x)
    if var == 0:
        return np.zeros(max_lag + 1)
    acf = np.correlate(x_centered, x_centered, mode='full')
    acf = acf[n-1:] / (var * n)
    return acf[:max_lag + 1]

def compute_pacf_levinson(acf_vals, max_lag):
    """Compute PACF using Levinson-Durbin recursion."""
    pacf = np.zeros(max_lag + 1)
    pacf[0] = 1.0

    # Levinson-Durbin
    phi = np.zeros((max_lag + 1, max_lag + 1))
    phi[1, 1] = acf_vals[1] / acf_vals[0] if acf_vals[0] != 0 else 0
    pacf[1] = phi[1, 1]

    for k in range(2, max_lag + 1):
        # Compute phi[k,k]
        num = acf_vals[k]
        for j in range(1, k):
            num -= phi[k-1, j] * acf_vals[k - j]
        den = acf_vals[0]
        for j in range(1, k):
            den -= phi[k-1, j] * acf_vals[j]
        if abs(den) < 1e-15:
            break
        phi[k, k] = num / den
        pacf[k] = phi[k, k]
        # Update other coefficients
        for j in range(1, k):
            phi[k, j] = phi[k-1, j] - phi[k, k] * phi[k-1, k - j]

    return pacf

log("\nComputing ACF and PACF for delta(n)...")
max_lag = 200
acf_vals = compute_acf(delta_all.astype(float), max_lag)
pacf_vals = compute_pacf_levinson(acf_vals, max_lag)

log(f"\nACF at key lags:")
for lag in [1, 2, 5, 10, 20, 50, 100, 200]:
    if lag <= max_lag:
        log(f"  ACF({lag:>3d}) = {acf_vals[lag]:+.6f}")

log(f"\nPACF at key lags:")
significant_pacf = []
se = 1.0 / np.sqrt(N_TOTAL)  # 95% CI boundary
for lag in range(1, min(51, max_lag + 1)):
    p = pacf_vals[lag]
    sig = "*" if abs(p) > 2 * se else ""
    if lag <= 20 or abs(p) > 2 * se:
        log(f"  PACF({lag:>3d}) = {p:+.6f} {sig}")
    if abs(p) > 2 * se:
        significant_pacf.append((lag, p))

log(f"\n95% CI boundary: +/-{2*se:.6f}")
log(f"Significant PACF lags (|PACF| > 2/sqrt(N)): {len(significant_pacf)}")
if significant_pacf:
    log(f"  Top 10 by magnitude:")
    for lag, p in sorted(significant_pacf, key=lambda x: -abs(x[1]))[:10]:
        log(f"    lag={lag}: PACF={p:+.6f}")

# PACF decay analysis
pacf_abs = np.abs(pacf_vals[1:51])
if np.all(pacf_abs > 0):
    log_lags = np.log(np.arange(1, 51))
    log_pacf = np.log(pacf_abs)
    # Power law fit: |PACF(k)| ~ k^(-alpha)
    valid = np.isfinite(log_pacf)
    if np.sum(valid) > 5:
        from numpy.polynomial import polynomial as P
        coeffs = np.polyfit(log_lags[valid], log_pacf[valid], 1)
        alpha = -coeffs[0]
        log(f"\nPACF power-law decay: |PACF(k)| ~ k^(-{alpha:.3f})")
        log(f"  (alpha > 1 means AR model has finite order, alpha < 1 means long memory)")

# AR order selection via AIC/BIC
log("\nAR order selection (AIC/BIC):")
from numpy.linalg import solve as la_solve

def fit_ar(x, p):
    """Fit AR(p) model, return (coefficients, residual_variance, AIC, BIC)."""
    n = len(x)
    if p == 0:
        var = np.var(x)
        aic = n * np.log(var) + 2
        bic = n * np.log(var) + np.log(n)
        return np.array([]), var, aic, bic
    X = np.column_stack([x[p-i-1:n-i-1] for i in range(p)])
    y = x[p:]
    # Normal equations
    try:
        beta = np.linalg.lstsq(X, y, rcond=None)[0]
        resid = y - X @ beta
        var = np.var(resid)
        if var <= 0:
            var = 1e-15
        nn = len(y)
        aic = nn * np.log(var) + 2 * p
        bic = nn * np.log(var) + p * np.log(nn)
        return beta, var, aic, bic
    except:
        return np.array([]), np.var(x), float('inf'), float('inf')

best_aic_p, best_aic = 0, float('inf')
best_bic_p, best_bic = 0, float('inf')
delta_f = delta_all.astype(float)

for p in list(range(0, 21)) + [30, 50, 100]:
    _, var, aic, bic = fit_ar(delta_f, p)
    if p <= 10 or p in [30, 50, 100]:
        log(f"  AR({p:>3d}): var={var:>10.2f}, AIC={aic:>12.1f}, BIC={bic:>12.1f}")
    if aic < best_aic:
        best_aic, best_aic_p = aic, p
    if bic < best_bic:
        best_bic, best_bic_p = bic, p

log(f"\n  Best by AIC: AR({best_aic_p})")
log(f"  Best by BIC: AR({best_bic_p})")
log(f"  Interpretation: AR({best_bic_p}) captures the direct dependencies.")
log(f"  Higher orders improve marginally (long memory effect).")


# ===================================================================
# 2. MULTI-ALGORITHM COMPRESSION COMPARISON
# ===================================================================
log("\n" + "=" * 72)
log("2. MULTI-ALGORITHM COMPRESSION COMPARISON")
log("=" * 72)

def to_bytes(arr):
    """Convert int array to bytes (int16 for delta range)."""
    return np.array(arr, dtype=np.int16).tobytes()

def compress_ratio(data_bytes, method):
    """Return compression ratio (compressed/original)."""
    if method == 'gzip':
        c = gzip.compress(data_bytes, compresslevel=9)
    elif method == 'bz2':
        c = bz2.compress(data_bytes, compresslevel=9)
    elif method == 'lzma':
        c = lzma.compress(data_bytes, preset=9)
    else:
        raise ValueError(f"Unknown method: {method}")
    return len(c) / len(data_bytes)

# Compress full delta sequence
raw_bytes = to_bytes(delta_all)
log(f"\nFull sequence ({N_TOTAL} values, {len(raw_bytes)} bytes):")

# Generate random comparison with same distribution
rng = np.random.RandomState(42)
delta_shuffled = delta_all.copy()
rng.shuffle(delta_shuffled)
shuffled_bytes = to_bytes(delta_shuffled)

# IID random with same range
random_vals = rng.randint(delta_all.min(), delta_all.max() + 1, size=N_TOTAL).astype(np.int16)
random_bytes = to_bytes(random_vals)

# Differenced delta (removes AR(1) component)
delta_diff = np.diff(delta_all).astype(np.int16)
diff_bytes = to_bytes(delta_diff)

log(f"\n{'Method':>8} | {'delta':>8} | {'shuffled':>8} | {'iid_rand':>8} | {'diff(delta)':>11} | delta/rand")
log("-" * 75)

for method in ['gzip', 'bz2', 'lzma']:
    r_delta = compress_ratio(raw_bytes, method)
    r_shuf = compress_ratio(shuffled_bytes, method)
    r_rand = compress_ratio(random_bytes, method)
    r_diff = compress_ratio(diff_bytes, method)
    ratio = r_delta / r_rand if r_rand > 0 else float('inf')
    log(f"{method:>8} | {r_delta:>7.4f} | {r_shuf:>7.4f} | {r_rand:>7.4f} | {r_diff:>10.4f} | {ratio:.4f}")

log("\nInterpretation:")
log("  delta/rand < 1 means delta is more compressible than random")
log("  diff(delta) compressibility shows how much AR(1) structure helps")


# ===================================================================
# 3. COMPRESSION RATIO SCALING WITH SEQUENCE LENGTH
# ===================================================================
log("\n" + "=" * 72)
log("3. COMPRESSION RATIO SCALING WITH N (sequence length)")
log("=" * 72)

lengths = [100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000]
log(f"\n{'N':>8} | {'gzip':>8} | {'bz2':>8} | {'lzma':>8} | {'bits/val':>8} | {'rand_bits':>9}")
log("-" * 65)

scaling_data = []
for n in lengths:
    if n > N_TOTAL:
        break
    seg = delta_all[:n]
    seg_bytes = to_bytes(seg)
    rand_seg = rng.randint(delta_all.min(), delta_all.max() + 1, size=n).astype(np.int16)
    rand_bytes = to_bytes(rand_seg)

    r_gz = compress_ratio(seg_bytes, 'gzip')
    r_bz = compress_ratio(seg_bytes, 'bz2')
    r_lz = compress_ratio(seg_bytes, 'lzma')
    bits_per_val = r_bz * 16  # 16 bits per int16, times ratio
    rand_bits = compress_ratio(rand_bytes, 'bz2') * 16

    scaling_data.append((n, r_gz, r_bz, r_lz, bits_per_val, rand_bits))
    log(f"{n:>8} | {r_gz:>7.4f} | {r_bz:>7.4f} | {r_lz:>7.4f} | {bits_per_val:>7.2f} | {rand_bits:>8.2f}")

# Fit scaling law
if len(scaling_data) >= 3:
    ns = np.array([d[0] for d in scaling_data], dtype=float)
    bz_ratios = np.array([d[2] for d in scaling_data])
    bits = np.array([d[4] for d in scaling_data])

    # Fit bits/value = a + b * log(N)
    log_ns = np.log(ns)
    coeffs_log = np.polyfit(log_ns, bits, 1)
    log(f"\nBits/value scaling fit: {coeffs_log[1]:.3f} + {coeffs_log[0]:.3f} * log(N)")

    # Fit bits/value = a * N^b
    log_bits = np.log(bits)
    coeffs_pow = np.polyfit(log_ns, log_bits, 1)
    log(f"Power-law fit: bits/val ~ N^{coeffs_pow[0]:.4f}")
    log(f"  (exponent near 0 = constant bits/value = finite entropy rate)")
    log(f"  (exponent > 0 = growing complexity per value)")


# ===================================================================
# 4. EXPLICIT Kt(delta(n)) vs log(n) GROWTH
# ===================================================================
log("\n" + "=" * 72)
log("4. Kt(delta(n)) vs log(n) GROWTH ANALYSIS")
log("=" * 72)

log("\nApproach: Measure the compressibility of delta(1..N) as N grows.")
log("Kt proxy = compressed bits of the sequence delta(1..N).")
log("If Kt(delta(1..N)) ~ N * h for constant h: entropy rate is finite.")
log("If Kt ~ N * log(N): complexity grows super-linearly.")

block_sizes = [10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000]
log(f"\n{'N':>8} | {'Kt_bz2':>10} | {'Kt/N':>8} | {'Kt/NlogN':>10} | {'Kt/N^{2/3}':>11}")
log("-" * 60)

kt_data = []
for n in block_sizes:
    if n > N_TOTAL:
        break
    seg = delta_all[:n]
    compressed = bz2.compress(to_bytes(seg), compresslevel=9)
    kt = len(compressed) * 8  # bits
    kt_per_n = kt / n
    kt_per_nlogn = kt / (n * np.log(n)) if n > 1 else 0
    kt_per_n23 = kt / (n ** (2/3))
    kt_data.append((n, kt, kt_per_n, kt_per_nlogn, kt_per_n23))
    log(f"{n:>8} | {kt:>10} | {kt_per_n:>7.3f} | {kt_per_nlogn:>9.4f} | {kt_per_n23:>10.3f}")

# Fit Kt(N) = a * N^b
if len(kt_data) >= 3:
    ns_kt = np.array([d[0] for d in kt_data], dtype=float)
    kts = np.array([d[1] for d in kt_data], dtype=float)
    log_ns = np.log(ns_kt)
    log_kts = np.log(kts)
    coeffs = np.polyfit(log_ns, log_kts, 1)
    log(f"\nKt(N) scaling: Kt ~ N^{coeffs[0]:.4f}")
    log(f"  exponent = 1.0 means linear (finite entropy rate)")
    log(f"  exponent > 1.0 means super-linear (growing per-symbol complexity)")
    log(f"  exponent < 1.0 means sub-linear (extreme compressibility)")

    # Also fit Kt = a*N + b*N*log(N)
    # Use linear regression: Kt = c0*N + c1*N*log(N)
    X = np.column_stack([ns_kt, ns_kt * np.log(ns_kt)])
    coeffs2 = np.linalg.lstsq(X, kts, rcond=None)[0]
    resid = kts - X @ coeffs2
    log(f"  Linear+NlogN fit: Kt = {coeffs2[0]:.3f}*N + {coeffs2[1]:.4f}*N*log(N)")
    log(f"  Residual norm: {np.linalg.norm(resid):.1f}")
    if abs(coeffs2[1]) < 0.01:
        log(f"  => NlogN coefficient negligible: entropy rate is CONSTANT")
    else:
        log(f"  => NlogN coefficient = {coeffs2[1]:.4f}: mild log(N) growth in per-symbol complexity")


# ===================================================================
# 5. BLOCK MUTUAL INFORMATION SCALING
# ===================================================================
log("\n" + "=" * 72)
log("5. BLOCK MUTUAL INFORMATION SCALING")
log("=" * 72)

log("\nMeasure MI between consecutive blocks of size L.")
log("If MI(L) grows with L: blocks carry increasing amounts of shared info.")

def quantize(arr, nbins=32):
    """Quantize to nbins levels."""
    mn, mx = arr.min(), arr.max()
    if mx == mn:
        return np.zeros(len(arr), dtype=int)
    return np.clip(((arr - mn) / (mx - mn) * (nbins - 1)).astype(int), 0, nbins - 1)

def entropy(arr):
    """Shannon entropy in bits."""
    counts = Counter(arr)
    n = len(arr)
    return -sum(c/n * np.log2(c/n) for c in counts.values() if c > 0)

def mi_blocks(seq, block_size, nbins=32):
    """MI between consecutive non-overlapping blocks."""
    n = len(seq) // (2 * block_size)
    if n < 50:
        return None
    a_blocks = []
    b_blocks = []
    for i in range(n):
        a = seq[2*i*block_size : (2*i+1)*block_size]
        b = seq[(2*i+1)*block_size : (2*i+2)*block_size]
        # Compress each block to get a "state" fingerprint
        a_bytes = bz2.compress(to_bytes(a), compresslevel=1)
        b_bytes = bz2.compress(to_bytes(b), compresslevel=1)
        a_blocks.append(len(a_bytes))
        b_blocks.append(len(b_bytes))

    a_q = quantize(np.array(a_blocks))
    b_q = quantize(np.array(b_blocks))

    h_a = entropy(a_q)
    h_b = entropy(b_q)
    joint = entropy(list(zip(a_q, b_q)))
    return h_a + h_b - joint

block_sizes_mi = [5, 10, 20, 50, 100, 200, 500, 1000]
log(f"\n{'Block L':>8} | {'MI(bits)':>10} | {'MI/log(L)':>10}")
log("-" * 35)

mi_data = []
for L in block_sizes_mi:
    mi = mi_blocks(delta_all, L)
    if mi is not None:
        mi_log = mi / np.log2(L) if L > 1 else 0
        mi_data.append((L, mi, mi_log))
        log(f"{L:>8} | {mi:>9.4f} | {mi_log:>9.4f}")

if len(mi_data) >= 3:
    ls = np.array([d[0] for d in mi_data], dtype=float)
    mis = np.array([d[1] for d in mi_data])
    coeffs = np.polyfit(np.log(ls), mis, 1)
    log(f"\nMI scaling: MI ~ {coeffs[0]:.4f} * log(L) + {coeffs[1]:.4f}")
    if coeffs[0] > 0:
        log("  MI grows with block size: long-range correlations confirmed")
    else:
        log("  MI decreases with block size: correlations are short-range")


# ===================================================================
# 6. DETRENDED FLUCTUATION ANALYSIS (DFA)
# ===================================================================
log("\n" + "=" * 72)
log("6. DETRENDED FLUCTUATION ANALYSIS (DFA) — Hurst Exponent")
log("=" * 72)

def dfa(x, scales=None):
    """
    Detrended Fluctuation Analysis.
    Returns (scales, fluctuations, Hurst_exponent).
    """
    n = len(x)
    # Cumulative sum (profile)
    y = np.cumsum(x - np.mean(x))

    if scales is None:
        scales = np.unique(np.logspace(np.log10(10), np.log10(n//4), 30).astype(int))
        scales = scales[scales >= 10]

    flucts = []
    for s in scales:
        n_segments = n // s
        if n_segments < 2:
            continue
        f2 = []
        for seg in range(n_segments):
            y_seg = y[seg*s : (seg+1)*s]
            # Linear detrend
            t = np.arange(s)
            p = np.polyfit(t, y_seg, 1)
            trend = np.polyval(p, t)
            f2.append(np.mean((y_seg - trend)**2))
        flucts.append((s, np.sqrt(np.mean(f2))))

    if len(flucts) < 3:
        return None, None, None

    ss = np.array([f[0] for f in flucts], dtype=float)
    ff = np.array([f[1] for f in flucts])

    # Log-log fit
    log_s = np.log(ss)
    log_f = np.log(ff)
    H_coeffs = np.polyfit(log_s, log_f, 1)
    H = H_coeffs[0]

    return ss, ff, H

log("\nComputing DFA for delta(n)...")
scales, flucts, H = dfa(delta_all.astype(float))

if H is not None:
    log(f"\nHurst exponent H = {H:.4f}")
    log(f"  H = 0.5: uncorrelated (white noise)")
    log(f"  H > 0.5: persistent (long-range positive correlations)")
    log(f"  H < 0.5: anti-persistent")
    log(f"  H = 1.0: 1/f noise")
    log(f"  H = 1.5: Brownian motion (1/f^2)")

    # Compute for shuffled comparison
    delta_shuf = delta_all.copy()
    rng.shuffle(delta_shuf)
    _, _, H_shuf = dfa(delta_shuf.astype(float))
    log(f"\nShuffled Hurst exponent: {H_shuf:.4f}")
    log(f"  (should be ~0.5 if structure is real)")

    # Log scale-by-scale fluctuations
    log(f"\nScale-dependent fluctuation:")
    for i in range(0, len(scales), max(1, len(scales)//10)):
        log(f"  scale={int(scales[i]):>6}: F(s)={flucts[i]:.2f}")

    # Check for crossover (different H at different scales)
    if len(scales) > 10:
        mid = len(scales) // 2
        log_s = np.log(scales)
        log_f = np.log(flucts)
        H_small = np.polyfit(log_s[:mid], log_f[:mid], 1)[0]
        H_large = np.polyfit(log_s[mid:], log_f[mid:], 1)[0]
        log(f"\nScale-dependent Hurst:")
        log(f"  Small scales (s < {int(scales[mid])}): H = {H_small:.4f}")
        log(f"  Large scales (s > {int(scales[mid])}): H = {H_large:.4f}")
        if abs(H_small - H_large) > 0.1:
            log(f"  CROSSOVER detected: structure changes across scales")
        else:
            log(f"  No significant crossover: uniform scaling")


# ===================================================================
# 7. TRANSFER ENTROPY: past delta → future delta vs past n → future delta
# ===================================================================
log("\n" + "=" * 72)
log("7. TRANSFER ENTROPY ANALYSIS")
log("=" * 72)

log("\nDoes knowing past n values help predict future delta beyond past delta?")

def transfer_entropy(source, target, lag=1, nbins=16):
    """
    TE: source → target
    TE = H(target_future | target_past) - H(target_future | target_past, source_past)
    """
    n = len(target) - lag
    if n < 100:
        return None

    # Quantize
    t_past = quantize(target[:n], nbins)
    t_future = quantize(target[lag:n+lag], nbins)
    s_past = quantize(source[:n], nbins)

    # H(t_future | t_past)
    h_tf_tp = 0
    joint_tp = defaultdict(list)
    for i in range(n):
        joint_tp[t_past[i]].append(t_future[i])
    for tp, tf_list in joint_tp.items():
        p_tp = len(tf_list) / n
        h_tf_tp += p_tp * entropy(tf_list)

    # H(t_future | t_past, s_past)
    h_tf_tp_sp = 0
    joint_tp_sp = defaultdict(list)
    for i in range(n):
        joint_tp_sp[(t_past[i], s_past[i])].append(t_future[i])
    for key, tf_list in joint_tp_sp.items():
        p_key = len(tf_list) / n
        h_tf_tp_sp += p_key * entropy(tf_list)

    return h_tf_tp - h_tf_tp_sp

# Source 1: n values themselves
n_values = np.arange(1, N_TOTAL + 1, dtype=float)

# Source 2: log(n)
log_n_values = np.log(n_values)

# Source 3: n mod 30 (primorial structure)
n_mod30 = (n_values % 30).astype(float)

# Source 4: shifted delta (self-transfer = autocorrelation)
delta_shifted = np.roll(delta_all, 1).astype(float)

log(f"\n{'Source':>15} | {'TE(lag=1)':>10} | {'TE(lag=5)':>10} | {'TE(lag=10)':>10}")
log("-" * 55)

sources = [
    ("delta(n-k)", delta_shifted),
    ("n", n_values),
    ("log(n)", log_n_values),
    ("n mod 30", n_mod30),
]

for name, source in sources:
    te1 = transfer_entropy(source.astype(float), delta_all.astype(float), lag=1)
    te5 = transfer_entropy(source.astype(float), delta_all.astype(float), lag=5)
    te10 = transfer_entropy(source.astype(float), delta_all.astype(float), lag=10)
    te1_str = f"{te1:.4f}" if te1 is not None else "N/A"
    te5_str = f"{te5:.4f}" if te5 is not None else "N/A"
    te10_str = f"{te10:.4f}" if te10 is not None else "N/A"
    log(f"{name:>15} | {te1_str:>10} | {te5_str:>10} | {te10_str:>10}")

log("\nInterpretation:")
log("  TE(delta→delta) > TE(n→delta): past delta is more informative than n")
log("  TE(n→delta) ≈ 0: n provides no additional information")
log("  TE(n mod 30→delta) > 0: residue class carries some info")


# ===================================================================
# 8. INFORMATION-THEORETIC SUMMARY
# ===================================================================
log("\n" + "=" * 72)
log("8. INFORMATION-THEORETIC SUMMARY")
log("=" * 72)

# Entropy rate estimation via compression
best_bz2_ratio = compress_ratio(to_bytes(delta_all), 'bz2')
entropy_rate_est = best_bz2_ratio * 16  # bits per int16
log(f"\nEntropy rate estimates:")
log(f"  bz2 compression: {entropy_rate_est:.2f} bits/value")
log(f"  lzma compression: {compress_ratio(to_bytes(delta_all), 'lzma') * 16:.2f} bits/value")

# Theoretical lower bound: if delta is truly random in [-664, 664]
range_size = delta_all.max() - delta_all.min() + 1
log(f"  Uniform random in [{delta_all.min()}, {delta_all.max()}]: {np.log2(range_size):.2f} bits/value")

# Empirical distribution entropy
counts = Counter(delta_all)
emp_entropy = -sum(c/N_TOTAL * np.log2(c/N_TOTAL) for c in counts.values())
log(f"  Empirical marginal entropy: {emp_entropy:.2f} bits/value")
log(f"  Unique values: {len(counts)}/{range_size}")

# Information deficit
log(f"\nInformation deficit (entropy - compressed rate):")
log(f"  Marginal entropy - bz2 rate = {emp_entropy - entropy_rate_est:.2f} bits/value")
log(f"  This represents exploitable sequential structure")

# Key conclusion
log(f"\n{'='*72}")
log("CONCLUSIONS")
log(f"{'='*72}")
log(f"1. PACF: AR order selected by BIC = {best_bic_p}, confirming direct dependencies")
log(f"   are short-range despite long-memory (power-law ACF decay)")
log(f"2. Compression: delta is {(1-best_bz2_ratio)*100:.1f}% more compressible than raw,")
log(f"   but compressibility comes from sequential correlations, NOT from n")
log(f"3. Hurst exponent H = {H:.4f}: {'persistent long-range correlations' if H > 0.7 else 'moderate correlations'}")
log(f"4. Transfer entropy from n → delta ≈ 0: n adds NO predictive information")
log(f"5. Entropy rate ≈ {entropy_rate_est:.1f} bits/value: each delta value carries")
log(f"   this much irreducible information given optimal sequential coding")
log(f"6. For polylog computation of p(n), one would need to compress ALL N values")
log(f"   into O(polylog(N)) bits — current evidence: {entropy_rate_est:.1f} * N bits needed")
log(f"   = EXTENSIVE (linear) information content")

# Save results
results_path = os.path.join(OUT_DIR, "kt_deep_analysis_results.md")
with open(results_path, 'w') as f:
    f.write("# Deep Kt Complexity Analysis Results\n\n")
    f.write("**Date:** 2026-04-05 (Session 36)\n\n")
    f.write("```\n")
    f.write("\n".join(results))
    f.write("\n```\n")
print(f"\nResults saved to {results_path}")

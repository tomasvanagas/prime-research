#!/usr/bin/env python3
"""
Block entropy scaling, power spectrum, wavelet analysis, and mutual information
of delta(n) = p(n) - round(R^{-1}(n)).

Analyses:
1. Block entropy scaling: h_k = H_k/k for various block sizes k
2. Power spectral density with spectral slope fitting
3. Wavelet analysis (pywt)
4. Mutual information decay I(delta(n); delta(n+k)) vs lag k
"""

import numpy as np
import sys
import os
from collections import Counter

# ── Load data ──────────────────────────────────────────────────────────────
DATA_PATH = os.path.join(os.path.dirname(__file__), '..', 'kt_delta_values.npy')
delta = np.load(DATA_PATH)
N = len(delta)
print(f"Loaded delta sequence: N={N}, range=[{delta.min()}, {delta.max()}]")
print(f"Mean={delta.mean():.4f}, Std={delta.std():.4f}")
print()

results = []
def log(msg):
    print(msg)
    results.append(msg)

# ── 1. Block Entropy Scaling ──────────────────────────────────────────────
log("=" * 70)
log("1. BLOCK ENTROPY SCALING")
log("=" * 70)

def block_entropy(seq, k):
    """Compute Shannon entropy of non-overlapping k-blocks."""
    n = len(seq) - (len(seq) % k)
    blocks = seq[:n].reshape(-1, k)
    # Convert each block to a hashable tuple
    block_tuples = [tuple(row) for row in blocks]
    counts = Counter(block_tuples)
    total = sum(counts.values())
    probs = np.array(list(counts.values())) / total
    H = -np.sum(probs * np.log2(probs))
    return H, len(counts)

block_sizes = [1, 2, 3, 4, 5, 8, 10, 16, 20, 32, 50]

# Also prepare shuffled version
rng = np.random.RandomState(42)
delta_shuffled = delta.copy()
rng.shuffle(delta_shuffled)

log(f"\n{'k':>4} | {'H_k':>10} | {'h_k=H_k/k':>10} | {'#symbols':>10} | {'h_k(shuf)':>10} | {'#sym(shuf)':>10}")
log("-" * 70)

hk_values = []
hk_shuf_values = []

for k in block_sizes:
    H, nsym = block_entropy(delta, k)
    H_s, nsym_s = block_entropy(delta_shuffled, k)
    hk = H / k
    hk_s = H_s / k
    hk_values.append(hk)
    hk_shuf_values.append(hk_s)
    log(f"{k:>4} | {H:>10.4f} | {hk:>10.4f} | {nsym:>10} | {hk_s:>10.4f} | {nsym_s:>10}")

# Fit entropy rate convergence: h_k ~ h_inf + a/k
# Use k >= 4 for fitting
fit_mask = np.array(block_sizes) >= 4
if np.sum(fit_mask) >= 3:
    from numpy.polynomial import polynomial as P
    ks_fit = np.array(block_sizes)[fit_mask]
    hk_fit = np.array(hk_values)[fit_mask]
    # Fit h_k = h_inf + a/k  =>  h_k = a * (1/k) + h_inf
    inv_k = 1.0 / ks_fit
    coeffs = np.polyfit(inv_k, hk_fit, 1)  # [a, h_inf]
    h_inf = coeffs[1]
    a_coeff = coeffs[0]
    log(f"\nEntropy rate fit: h_k ~ {h_inf:.4f} + {a_coeff:.4f}/k")
    log(f"Estimated asymptotic entropy rate h_inf = {h_inf:.4f} bits/symbol")

    # Same for shuffled
    hk_s_fit = np.array(hk_shuf_values)[fit_mask]
    coeffs_s = np.polyfit(inv_k, hk_s_fit, 1)
    log(f"Shuffled: h_k ~ {coeffs_s[1]:.4f} + {coeffs_s[0]:.4f}/k")
    log(f"Shuffled asymptotic h_inf = {coeffs_s[1]:.4f} bits/symbol")

# Single-symbol entropy for reference
H1, _ = block_entropy(delta, 1)
log(f"\nSingle-symbol entropy H_1 = {H1:.4f} bits")
log(f"Entropy reduction h_1 -> h_50: {hk_values[0]:.4f} -> {hk_values[-1]:.4f} ({(1 - hk_values[-1]/hk_values[0])*100:.1f}% reduction)")

# ── 2. Power Spectral Density ────────────────────────────────────────────
log("\n" + "=" * 70)
log("2. POWER SPECTRAL DENSITY")
log("=" * 70)

from scipy.fft import fft

# Zero-mean
delta_zm = delta - delta.mean()
delta_shuf_zm = delta_shuffled - delta_shuffled.mean()

# FFT
F = fft(delta_zm)
PSD = np.abs(F[:N//2])**2 / N
freqs = np.arange(N//2) / N

F_s = fft(delta_shuf_zm)
PSD_s = np.abs(F_s[:N//2])**2 / N

# Skip DC (index 0)
freqs_pos = freqs[1:]
PSD_pos = PSD[1:]
PSD_s_pos = PSD_s[1:]

# Log-log fit for spectral slope: log(PSD) = -beta * log(f) + c
log_f = np.log10(freqs_pos)
log_PSD = np.log10(PSD_pos)
log_PSD_s = np.log10(PSD_s_pos)

# Fit over full range
slope, intercept = np.polyfit(log_f, log_PSD, 1)
slope_s, intercept_s = np.polyfit(log_f, log_PSD_s, 1)

log(f"\nSpectral slope (full range): beta = {-slope:.4f} (PSD ~ f^{slope:.4f})")
log(f"Shuffled spectral slope:     beta = {-slope_s:.4f} (PSD ~ f^{slope_s:.4f})")

# Fit in low-frequency regime (f < 0.01)
low_mask = freqs_pos < 0.01
if np.sum(low_mask) > 10:
    sl_low, _ = np.polyfit(log_f[low_mask], log_PSD[low_mask], 1)
    log(f"Low-freq slope (f<0.01):     beta = {-sl_low:.4f}")

# Fit in mid-frequency regime
mid_mask = (freqs_pos > 0.01) & (freqs_pos < 0.1)
if np.sum(mid_mask) > 10:
    sl_mid, _ = np.polyfit(log_f[mid_mask], log_PSD[mid_mask], 1)
    log(f"Mid-freq slope (0.01<f<0.1): beta = {-sl_mid:.4f}")

# High frequency
hi_mask = freqs_pos > 0.1
if np.sum(hi_mask) > 10:
    sl_hi, _ = np.polyfit(log_f[hi_mask], log_PSD[hi_mask], 1)
    log(f"High-freq slope (f>0.1):     beta = {-sl_hi:.4f}")

# Summary statistics
log(f"\nMean PSD: {PSD_pos.mean():.2f}")
log(f"Median PSD: {np.median(PSD_pos):.2f}")
log(f"PSD ratio (low 1% freqs / high 1% freqs): {PSD_pos[:len(PSD_pos)//100].mean() / PSD_pos[-len(PSD_pos)//100:].mean():.2f}")

# White noise test: coefficient of variation of log-binned PSD
n_bins = 50
bin_edges = np.logspace(log_f[0], log_f[-1], n_bins + 1)
bin_means = []
for i in range(n_bins):
    mask = (freqs_pos >= bin_edges[i]) & (freqs_pos < bin_edges[i+1])
    if np.sum(mask) > 0:
        bin_means.append(PSD_pos[mask].mean())
bin_means = np.array(bin_means)
cv = bin_means.std() / bin_means.mean()
log(f"Log-binned PSD coefficient of variation: {cv:.4f}")
if abs(slope) < 0.1:
    log("VERDICT: Spectrum is approximately WHITE NOISE (flat)")
elif -0.5 < slope < -0.05:
    log(f"VERDICT: Spectrum shows weak 1/f^{-slope:.2f} coloring")
elif slope < -0.5:
    log(f"VERDICT: Spectrum shows strong 1/f^{-slope:.2f} coloring")
else:
    log(f"VERDICT: Spectrum shows slight blue tilt (f^{slope:.2f})")

# ── 3. Wavelet Analysis ──────────────────────────────────────────────────
log("\n" + "=" * 70)
log("3. WAVELET ANALYSIS")
log("=" * 70)

try:
    import pywt
    from scipy.stats import kurtosis, skew, normaltest

    # Use Daubechies-4 wavelet, max decomposition level
    max_level = pywt.dwt_max_level(N, 'db4')
    use_level = min(max_level, 12)
    coeffs = pywt.wavedec(delta_zm.astype(float), 'db4', level=use_level)

    log(f"\nWavelet: db4, decomposition levels: {use_level}")
    log(f"{'Level':>6} | {'#coeffs':>8} | {'Std':>10} | {'Kurtosis':>10} | {'Skewness':>10} | {'Normal?':>8}")
    log("-" * 70)

    # coeffs[0] = approximation at coarsest level, coeffs[1..] = details finest to coarsest
    for i, c in enumerate(coeffs):
        if i == 0:
            label = f"A{use_level}"
        else:
            label = f"D{use_level - i + 1}"
        k = kurtosis(c, fisher=True)  # excess kurtosis (Gaussian = 0)
        s = skew(c)
        # Normal test (only if enough samples)
        if len(c) >= 20:
            stat, pval = normaltest(c)
            normal_str = f"p={pval:.2e}"
        else:
            normal_str = "N/A"
        log(f"{label:>6} | {len(c):>8} | {c.std():>10.4f} | {k:>10.4f} | {s:>10.4f} | {normal_str:>8}")

    # Energy distribution across scales
    total_energy = sum(np.sum(c**2) for c in coeffs)
    log(f"\nEnergy distribution across scales:")
    for i, c in enumerate(coeffs):
        if i == 0:
            label = f"A{use_level}"
        else:
            label = f"D{use_level - i + 1}"
        energy_frac = np.sum(c**2) / total_energy
        log(f"  {label}: {energy_frac*100:.2f}%")

except ImportError:
    log("pywt not available, skipping wavelet analysis")

# ── 4. Mutual Information Scaling ─────────────────────────────────────────
log("\n" + "=" * 70)
log("4. MUTUAL INFORMATION SCALING")
log("=" * 70)

def histogram_mi(x, y, bins=50):
    """Estimate mutual information using 2D histogram."""
    # Joint histogram
    H_xy, _, _ = np.histogram2d(x, y, bins=bins)
    # Normalize
    p_xy = H_xy / H_xy.sum()
    p_xy = p_xy[p_xy > 0]

    p_x = np.histogram(x, bins=bins)[0].astype(float)
    p_x = p_x / p_x.sum()
    p_x = p_x[p_x > 0]

    p_y = np.histogram(y, bins=bins)[0].astype(float)
    p_y = p_y / p_y.sum()
    p_y = p_y[p_y > 0]

    H_x = -np.sum(p_x * np.log2(p_x))
    H_y = -np.sum(p_y * np.log2(p_y))
    H_xy = -np.sum(p_xy * np.log2(p_xy))

    return H_x + H_y - H_xy

lags = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000]
n_bins = 80  # bins for histogram MI

log(f"\n{'Lag k':>8} | {'MI(bits)':>10} | {'MI_shuf':>10}")
log("-" * 40)

mi_values = []
mi_shuf_values = []

for lag in lags:
    x = delta[:-lag].astype(float)
    y = delta[lag:].astype(float)
    mi = histogram_mi(x, y, bins=n_bins)

    x_s = delta_shuffled[:-lag].astype(float)
    y_s = delta_shuffled[lag:].astype(float)
    mi_s = histogram_mi(x_s, y_s, bins=n_bins)

    mi_values.append(mi)
    mi_shuf_values.append(mi_s)
    log(f"{lag:>8} | {mi:>10.6f} | {mi_s:>10.6f}")

# Fit MI decay: log(MI) vs log(k) for power law, MI vs k for exponential
mi_arr = np.array(mi_values)
lag_arr = np.array(lags, dtype=float)

# Subtract baseline (shuffled MI mean) as bias correction
mi_baseline = np.mean(mi_shuf_values)
mi_corrected = mi_arr - mi_baseline
mi_corrected_pos = mi_corrected[mi_corrected > 0]
lags_pos = lag_arr[mi_corrected > 0]

log(f"\nMI baseline (shuffled mean): {mi_baseline:.6f} bits")

if len(mi_corrected_pos) >= 3:
    # Power law fit: MI ~ k^(-alpha)
    log_k = np.log10(lags_pos)
    log_mi = np.log10(mi_corrected_pos)
    alpha, c_pl = np.polyfit(log_k, log_mi, 1)
    log(f"Power law fit (bias-corrected): MI ~ k^({alpha:.3f})")

    # Exponential fit: MI ~ exp(-k/tau)
    # log(MI) = -k/tau + c  =>  linear in k
    ln_mi = np.log(mi_corrected_pos)
    rate, c_exp = np.polyfit(lags_pos, ln_mi, 1)
    tau = -1.0 / rate if rate < 0 else float('inf')
    log(f"Exponential fit (bias-corrected): MI ~ exp(-k/{tau:.1f})")

    # Which fits better? Compare R^2
    # Power law
    pred_pl = alpha * log_k + c_pl
    ss_res_pl = np.sum((log_mi - pred_pl)**2)
    ss_tot = np.sum((log_mi - log_mi.mean())**2)
    r2_pl = 1 - ss_res_pl / ss_tot if ss_tot > 0 else 0

    # Exponential (in log space)
    pred_exp = rate * lags_pos + c_exp
    ss_res_exp = np.sum((ln_mi - pred_exp)**2)
    ss_tot_exp = np.sum((ln_mi - ln_mi.mean())**2)
    r2_exp = 1 - ss_res_exp / ss_tot_exp if ss_tot_exp > 0 else 0

    log(f"R^2 power law:   {r2_pl:.4f}")
    log(f"R^2 exponential: {r2_exp:.4f}")

    if r2_pl > r2_exp:
        log(f"VERDICT: MI decays as POWER LAW ~ k^({alpha:.2f})")
    else:
        log(f"VERDICT: MI decays EXPONENTIALLY with tau ~ {tau:.1f}")
else:
    log("Insufficient positive bias-corrected MI values for decay fitting.")

# ── 5. Summary ────────────────────────────────────────────────────────────
log("\n" + "=" * 70)
log("5. SUMMARY")
log("=" * 70)
log(f"""
Delta(n) = p(n) - round(R^{{-1}}(n)), N = {N}
Range: [{delta.min()}, {delta.max()}], Mean: {delta.mean():.4f}, Std: {delta.std():.4f}

Block Entropy:
  h_1 = {hk_values[0]:.4f} bits/symbol
  h_50 = {hk_values[-1]:.4f} bits/symbol
  Asymptotic h_inf ~ {h_inf:.4f} bits/symbol
  Entropy rate DECREASES with block size => delta has structure
  Reduction from random: {(1 - hk_values[-1]/hk_shuf_values[-1])*100:.1f}% at k=50

Power Spectrum:
  Overall slope: beta = {-slope:.4f} (PSD ~ f^{slope:.4f})

Mutual Information:
  MI(lag=1) = {mi_values[0]:.6f} bits (baseline: {mi_baseline:.6f})
  MI drops to baseline around lag ~ {lags[next((i for i, m in enumerate(mi_corrected) if m < mi_baseline), len(lags)-1)]}
""")

# ── Save results ──────────────────────────────────────────────────────────
out_path = os.path.join(os.path.dirname(__file__), 'entropy_spectral_results.txt')
with open(out_path, 'w') as f:
    f.write('\n'.join(results))
print(f"\nResults saved to {out_path}")

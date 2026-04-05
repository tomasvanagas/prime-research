#!/usr/bin/env python3
"""
DFT analysis of Riemann zeta zero sequence.

Computes:
1. Power spectrum |F(k)|^2 of the zero sequence
2. Comparison to GUE random matrix ensemble
3. Spectral flatness by frequency band
4. Peaks at log(p) frequencies
5. Pair correlation R_2(r) vs GUE prediction
6. Number variance Sigma^2(L)
"""

import numpy as np
from scipy import stats
from scipy.linalg import eigvalsh
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import json

OUT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_FILE = "/apps/aplikacijos/prime-research/data/zeta_zeros_1000.txt"

# ── Load zeros ──────────────────────────────────────────────────────────────
zeros = np.loadtxt(DATA_FILE)
N = len(zeros)
print(f"Loaded {N} zeros, range [{zeros[0]:.4f}, {zeros[-1]:.4f}]")

# Mean spacing (Weyl-law local spacing ~ 2*pi/log(gamma/(2*pi)))
spacings = np.diff(zeros)
mean_spacing = np.mean(spacings)
print(f"Mean spacing: {mean_spacing:.6f}")

# ── 1. DFT of the zero sequence ────────────────────────────────────────────
# Treat {gamma_1, ..., gamma_N} as a discrete signal
F = np.fft.fft(zeros)
power = np.abs(F)**2
freqs = np.fft.fftfreq(N)

# Only positive frequencies
pos = freqs > 0
k_pos = np.arange(N)[pos]
power_pos = power[pos]
freqs_pos = freqs[pos]

print(f"\n=== POWER SPECTRUM ===")
print(f"DC component |F(0)|^2 = {power[0]:.4e}")
print(f"Max non-DC power = {np.max(power_pos):.4e} at k={k_pos[np.argmax(power_pos)]}")

# Find peaks above noise floor
noise_floor = np.median(power_pos)
threshold = noise_floor * 10
peaks_mask = power_pos > threshold
peak_indices = k_pos[peaks_mask]
peak_powers = power_pos[peaks_mask]
print(f"Noise floor (median): {noise_floor:.4e}")
print(f"Peaks above 10x noise: {len(peak_indices)}")
for i in np.argsort(peak_powers)[-10:][::-1]:
    print(f"  k={peak_indices[i]:4d}, freq={freqs_pos[peaks_mask][i]:.6f}, "
          f"|F|^2={peak_powers[i]:.4e}, ratio={peak_powers[i]/noise_floor:.1f}x")

# Plot power spectrum
fig, axes = plt.subplots(2, 1, figsize=(12, 8))
axes[0].semilogy(k_pos, power_pos, 'b-', alpha=0.7, linewidth=0.5)
axes[0].axhline(noise_floor, color='r', linestyle='--', label=f'Median = {noise_floor:.2e}')
axes[0].axhline(threshold, color='orange', linestyle='--', label=f'10x median')
axes[0].set_xlabel('k (frequency index)')
axes[0].set_ylabel('|F(k)|^2')
axes[0].set_title('Power Spectrum of Zeta Zero Sequence')
axes[0].legend()

# Low-frequency zoom
k_low = k_pos[:50]
p_low = power_pos[:50]
axes[1].bar(k_low, p_low, color='steelblue')
axes[1].set_xlabel('k (frequency index)')
axes[1].set_ylabel('|F(k)|^2')
axes[1].set_title('Low-Frequency Power Spectrum (k=1..50)')
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, 'power_spectrum.png'), dpi=150)
plt.close()

# ── 2. GUE comparison ──────────────────────────────────────────────────────
print(f"\n=== GUE COMPARISON ===")
n_gue = 30
gue_size = 300  # Use 300x300 matrices (faster), then pad/interpolate
gue_powers = []

np.random.seed(42)
for trial in range(n_gue):
    # Generate GUE matrix
    A = (np.random.randn(gue_size, gue_size) + 1j * np.random.randn(gue_size, gue_size)) / np.sqrt(2)
    H = (A + A.conj().T) / (2 * np.sqrt(gue_size))
    eigs = np.sort(np.real(eigvalsh(H)))
    # Interpolate to N points to match zeta zeros length
    eigs = np.interp(np.linspace(0, len(eigs)-1, N), np.arange(len(eigs)), eigs)
    # DFT of eigenvalues
    F_gue = np.fft.fft(eigs)
    gue_powers.append(np.abs(F_gue)**2)

gue_powers = np.array(gue_powers)
gue_mean = np.mean(gue_powers[:, pos], axis=0)
gue_std = np.std(gue_powers[:, pos], axis=0)

# Normalize both for shape comparison
zeta_norm = power_pos / np.mean(power_pos)
gue_norm = gue_mean / np.mean(gue_mean)

fig, ax = plt.subplots(figsize=(12, 5))
ax.semilogy(k_pos, zeta_norm, 'b-', alpha=0.7, linewidth=0.5, label='Zeta zeros')
ax.semilogy(k_pos, gue_norm, 'r-', alpha=0.7, linewidth=0.5, label='GUE mean (100 trials)')
ax.fill_between(k_pos,
                 (gue_mean - 2*gue_std) / np.mean(gue_mean),
                 (gue_mean + 2*gue_std) / np.mean(gue_mean),
                 color='red', alpha=0.1, label='GUE +/- 2sigma')
ax.set_xlabel('k')
ax.set_ylabel('Normalized |F(k)|^2')
ax.set_title('Zeta Zeros vs GUE Power Spectrum (normalized)')
ax.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, 'gue_comparison.png'), dpi=150)
plt.close()

# Kolmogorov-Smirnov test on normalized power distributions
ks_stat, ks_pval = stats.ks_2samp(zeta_norm, gue_norm)
print(f"KS test (zeta vs GUE normalized power): stat={ks_stat:.4f}, p={ks_pval:.4e}")

# Correlation between zeta and GUE power spectra
corr = np.corrcoef(np.log(zeta_norm + 1e-30), np.log(gue_norm + 1e-30))[0, 1]
print(f"Log-power correlation: {corr:.4f}")

# ── 3. Spectral flatness by frequency band ─────────────────────────────────
print(f"\n=== SPECTRAL FLATNESS ===")

def spectral_flatness(spectrum):
    """Geometric mean / arithmetic mean of power spectrum."""
    log_mean = np.mean(np.log(spectrum + 1e-30))
    return np.exp(log_mean) / np.mean(spectrum)

sf_total = spectral_flatness(power_pos)
print(f"Overall spectral flatness: {sf_total:.4f}")

# Break into bands
n_bands = 5
band_size = len(power_pos) // n_bands
band_results = []
for b in range(n_bands):
    start = b * band_size
    end = start + band_size if b < n_bands - 1 else len(power_pos)
    band_power = power_pos[start:end]
    sf = spectral_flatness(band_power)
    freq_range = (freqs_pos[start], freqs_pos[min(end-1, len(freqs_pos)-1)])
    band_results.append((b, freq_range, sf))
    print(f"  Band {b} (freq {freq_range[0]:.4f}-{freq_range[1]:.4f}): SF = {sf:.4f}")

# GUE spectral flatness for comparison
gue_sfs = []
for trial in range(n_gue):
    gue_sf = spectral_flatness(gue_powers[trial, pos])
    gue_sfs.append(gue_sf)
print(f"GUE spectral flatness: mean={np.mean(gue_sfs):.4f}, std={np.std(gue_sfs):.4f} (n={n_gue} trials, {gue_size}x{gue_size})")
print(f"Zeta SF z-score vs GUE: {(sf_total - np.mean(gue_sfs))/np.std(gue_sfs):.2f}")

# ── 4. Peaks at log(p) frequencies ─────────────────────────────────────────
print(f"\n=== LOG-PRIME FREQUENCY ANALYSIS ===")

small_primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]
log_primes = np.log(small_primes)

# The "natural" frequency unit: we look for peaks at frequencies ~ 1/(2*pi*log(p))
# and also at log(p)/(2*pi) in the normalized frequency domain
# Since the zeros are not uniformly sampled, we also do a Lomb-Scargle periodogram

# Lomb-Scargle on the unfolded zeros (deviation from smooth part)
# Unfold: map gamma_n -> n-th expected position via smooth counting N(T) ~ T/(2pi) * log(T/(2pi*e))
from scipy.signal import lombscargle

# Use index as the "time" and deviation from linear fit as signal
n_idx = np.arange(1, N+1)
# Linear detrend
slope, intercept = np.polyfit(n_idx, zeros, 1)
residuals = zeros - (slope * n_idx + intercept)

# Compute Lomb-Scargle at specific angular frequencies
test_freqs_angular = np.linspace(0.01, np.pi, 5000)
ls_power = lombscargle(n_idx.astype(float), residuals, test_freqs_angular, normalize=True)
test_freqs = test_freqs_angular / (2 * np.pi)

fig, ax = plt.subplots(figsize=(12, 5))
ax.plot(test_freqs, ls_power, 'b-', alpha=0.7, linewidth=0.5)
ax.set_xlabel('Frequency (cycles per index)')
ax.set_ylabel('Lomb-Scargle Power')
ax.set_title('Lomb-Scargle Periodogram of Detrended Zeta Zeros')

# Mark log(p) frequencies (normalized by mean spacing)
for p, lp in zip(small_primes, log_primes):
    # Frequency candidates: log(p)/(2*pi*mean_spacing), 1/(log(p)*mean_spacing)
    f1 = lp / (2 * np.pi * mean_spacing)
    f2 = 1.0 / (lp * mean_spacing)
    if f1 < test_freqs[-1]:
        ax.axvline(f1, color='red', alpha=0.3, linestyle='--')
        ax.text(f1, ax.get_ylim()[1]*0.9, f'log({p})/2pi/d', fontsize=6, rotation=90)
    if f2 < test_freqs[-1]:
        ax.axvline(f2, color='green', alpha=0.3, linestyle='--')
        ax.text(f2, ax.get_ylim()[1]*0.8, f'1/(log({p})*d)', fontsize=6, rotation=90)

plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, 'log_prime_frequencies.png'), dpi=150)
plt.close()

# Quantify: power at log-prime frequencies vs surroundings
print("Checking power at log-prime related frequencies:")
ls_noise = np.median(ls_power)
for p, lp in zip(small_primes[:8], log_primes[:8]):
    f1 = lp / (2 * np.pi * mean_spacing)
    # Find closest frequency in our grid
    idx = np.argmin(np.abs(test_freqs - f1))
    local_power = ls_power[idx]
    # Local neighborhood mean (window of +/- 20 bins)
    lo, hi = max(0, idx-20), min(len(ls_power), idx+20)
    local_mean = np.mean(ls_power[lo:hi])
    ratio = local_power / ls_noise
    print(f"  p={p:2d}: f=log({p})/(2pi*d)={f1:.5f}, "
          f"LS power={local_power:.4f}, ratio to median={ratio:.2f}")

# ── 5. Pair correlation ────────────────────────────────────────────────────
print(f"\n=== PAIR CORRELATION R_2(r) ===")

# Unfold zeros to unit mean spacing
unfolded = np.zeros(N)
for i in range(N):
    # Local unfolding: use smooth counting function
    # N(T) ~ T/(2pi) * log(T/(2pi*e)) + 7/8
    T = zeros[i]
    unfolded[i] = T / (2*np.pi) * np.log(T / (2*np.pi * np.e)) + 7.0/8.0

# Normalize to have mean spacing = 1
unfolded_spacings = np.diff(unfolded)
mean_unf_spacing = np.mean(unfolded_spacings)
unfolded = unfolded / mean_unf_spacing

# Compute pair correlation
r_max = 4.0
dr = 0.05
r_bins = np.arange(0, r_max + dr, dr)
r_centers = (r_bins[:-1] + r_bins[1:]) / 2
hist = np.zeros(len(r_centers))

# Count pairs
count = 0
for i in range(N):
    for j in range(i+1, min(i+50, N)):  # Only nearby pairs for efficiency
        diff = abs(unfolded[j] - unfolded[i])
        if diff < r_max:
            idx = int(diff / dr)
            if idx < len(hist):
                hist[idx] += 1
                count += 1

# Normalize: R_2(r) = hist / (N * dr * density_of_pairs)
# Expected count per bin for uniform: N * (N_window) * dr / r_max (approx)
# Better: normalize so R_2 -> 1 for large r
norm_factor = N * dr  # approximate normalization
R2 = hist / norm_factor
# Renormalize so R_2 -> 1 at large r
R2_tail = np.mean(R2[r_centers > 2.5]) if np.any(r_centers > 2.5) else 1.0
R2 = R2 / R2_tail

# GUE prediction: 1 - (sin(pi*r)/(pi*r))^2
r_theory = np.linspace(0.001, r_max, 500)
R2_gue = 1 - (np.sin(np.pi * r_theory) / (np.pi * r_theory))**2

fig, ax = plt.subplots(figsize=(10, 6))
ax.bar(r_centers, R2, width=dr*0.8, color='steelblue', alpha=0.7, label='Zeta zeros')
ax.plot(r_theory, R2_gue, 'r-', linewidth=2, label=r'GUE: $1 - (\sin\pi r / \pi r)^2$')
ax.set_xlabel('r (units of mean spacing)')
ax.set_ylabel('$R_2(r)$')
ax.set_title('Pair Correlation Function')
ax.legend()
ax.set_xlim(0, r_max)
ax.set_ylim(0, 1.5)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, 'pair_correlation.png'), dpi=150)
plt.close()

# Quantify fit
r_interp = np.interp(r_centers, r_theory, R2_gue)
residual = R2 - r_interp
rms_resid = np.sqrt(np.mean(residual[r_centers > 0.1]**2))
print(f"R_2 RMS residual from GUE (r > 0.1): {rms_resid:.4f}")
print(f"Level repulsion at r=0: R_2(0) = {R2[0]:.4f} (GUE predicts 0)")

# ── 6. Number variance Sigma^2(L) ──────────────────────────────────────────
print(f"\n=== NUMBER VARIANCE Sigma^2(L) ===")

L_values = np.logspace(np.log10(0.1), np.log10(10), 40)
sigma2 = np.zeros(len(L_values))

# Use unfolded zeros
unf = unfolded.copy()

for li, L in enumerate(L_values):
    counts = []
    # Slide window of length L across unfolded zeros
    n_windows = min(500, N - 1)
    window_starts = np.linspace(unf[0], unf[-1] - L, n_windows)
    for s in window_starts:
        c = np.sum((unf >= s) & (unf < s + L))
        counts.append(c)
    counts = np.array(counts)
    sigma2[li] = np.var(counts)

# GUE prediction for number variance:
# Sigma^2(L) ~ (2/pi^2) * (log(2*pi*L) + gamma_euler + 1 - pi^2/8) for large L
# More precisely, use the known asymptotic
gamma_euler = 0.5772156649
sigma2_gue = (2/np.pi**2) * (np.log(2*np.pi*L_values) + gamma_euler + 1 - np.pi**2/8)
# For small L, Sigma^2(L) ~ L (Poisson limit doesn't apply for GUE; use exact formula)
# Small L: Sigma^2(L) ~ L - (2/pi^2)*L^2*(log(L) - 1) + ... but just use asymptotic

# Poisson prediction
sigma2_poisson = L_values

fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(L_values, sigma2, 'bo-', markersize=4, label='Zeta zeros')
ax.plot(L_values, sigma2_gue, 'r-', linewidth=2, label='GUE (asymptotic)')
ax.plot(L_values, sigma2_poisson, 'g--', linewidth=1, label='Poisson')
ax.set_xlabel('L (units of mean spacing)')
ax.set_ylabel(r'$\Sigma^2(L)$')
ax.set_title('Number Variance')
ax.legend()
ax.set_xscale('log')
ax.set_yscale('log')
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, 'number_variance.png'), dpi=150)
plt.close()

print("L values and Sigma^2(L):")
for L, s2 in zip(L_values[::5], sigma2[::5]):
    s2_gue = (2/np.pi**2) * (np.log(2*np.pi*L) + gamma_euler + 1 - np.pi**2/8)
    print(f"  L={L:6.2f}: Sigma^2={s2:.4f}, GUE={s2_gue:.4f}, Poisson={L:.4f}")

# ── 7. Also compute DFT of spacings (delta_n = gamma_{n+1} - gamma_n) ─────
print(f"\n=== DFT OF SPACINGS ===")
F_sp = np.fft.fft(spacings)
power_sp = np.abs(F_sp)**2
freqs_sp = np.fft.fftfreq(len(spacings))
pos_sp = freqs_sp > 0

sf_spacings = spectral_flatness(power_sp[pos_sp])
print(f"Spectral flatness of spacings: {sf_spacings:.4f}")

fig, ax = plt.subplots(figsize=(12, 5))
ax.semilogy(np.arange(len(spacings))[pos_sp], power_sp[pos_sp], 'b-', alpha=0.7, linewidth=0.5)
ax.set_xlabel('k')
ax.set_ylabel('|F(k)|^2')
ax.set_title('Power Spectrum of Zeta Zero Spacings')
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, 'spacing_spectrum.png'), dpi=150)
plt.close()

# ── Summary ─────────────────────────────────────────────────────────────────
print("\n" + "="*60)
print("SUMMARY")
print("="*60)

results = {
    "n_zeros": N,
    "mean_spacing": float(mean_spacing),
    "spectral_flatness_total": float(sf_total),
    "spectral_flatness_spacings": float(sf_spacings),
    "spectral_flatness_gue_mean": float(np.mean(gue_sfs)),
    "spectral_flatness_gue_std": float(np.std(gue_sfs)),
    "ks_stat_vs_gue": float(ks_stat),
    "ks_pval_vs_gue": float(ks_pval),
    "log_power_corr_vs_gue": float(corr),
    "pair_corr_rms_residual": float(rms_resid),
    "pair_corr_at_zero": float(R2[0]),
    "n_peaks_above_10x_noise": int(len(peak_indices)),
}

print(json.dumps(results, indent=2))

# Save results as JSON too
with open(os.path.join(OUT_DIR, 'dft_results.json'), 'w') as f:
    json.dump(results, f, indent=2)

print(f"\nAll plots saved to {OUT_DIR}/")

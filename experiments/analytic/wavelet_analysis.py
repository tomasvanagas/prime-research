"""
SESSION 6: Deep Multi-Scale Wavelet Analysis of δ(n) = p(n) - approx(n)
=========================================================================
Investigates whether the correction term between p(n) and various
approximations has exploitable structure at any scale.

Analyses:
  1. Three approximations: Li^{-1}(n), R^{-1}(n), Cipolla-style
  2. Wavelet decomposition (Daubechies, Haar, Symlet)
  3. Power spectrum / spectral density -- is it f^{-α}?
  4. Chebyshev polynomial expansion of the correction
  5. Periodicity search against known constants (ln2, π, zeta zeros)
  6. Hurst exponent -- random walk or long-range correlated?
  7. Detrended fluctuation analysis (DFA)
"""

import numpy as np
import time
import sys
import math
from collections import defaultdict

# ─── Imports ───
from scipy import signal, optimize, stats
from scipy.interpolate import interp1d
from scipy.special import expi  # Ei(x) for Li(x)
import pywt

# ═══════════════════════════════════════════════════════════════
# PART 0: Generate primes and compute δ(n) for three approximations
# FAST versions using scipy only (no mpmath per-call)
# ═══════════════════════════════════════════════════════════════

def sieve_primes(limit):
    """Sieve of Eratosthenes up to limit."""
    is_prime = bytearray([1]) * (limit + 1)
    is_prime[0] = is_prime[1] = 0
    for i in range(2, int(limit**0.5) + 1):
        if is_prime[i]:
            is_prime[i*i::i] = bytearray(len(is_prime[i*i::i]))
    return [i for i in range(2, limit + 1) if is_prime[i]]

def li_func(x):
    """Logarithmic integral li(x) = Ei(ln(x)) for x > 1."""
    if x <= 1:
        return -1e30
    return expi(math.log(x))

LI_2 = li_func(2.0)  # li(2) ≈ 1.0451

def Li_func(x):
    """Li(x) = li(x) - li(2) = integral from 2 to x of 1/ln(t) dt."""
    return li_func(x) - LI_2

def inverse_li_fast(n):
    """Fast Li^{-1}(n): solve Li(x) = n via Newton's method with float arithmetic."""
    if n <= 0:
        return 2.0
    # Initial guess: n * ln(n)
    x = max(2.1, n * math.log(max(n, 2)))
    for _ in range(50):
        val = Li_func(x)
        deriv = 1.0 / math.log(x)
        dx = (n - val) / deriv
        x = x + dx
        if abs(dx) < 1e-10:
            break
    return x

def R_func(x, n_terms=100):
    """
    Riemann R function: R(x) = sum_{k=1}^{inf} mu(k)/k * li(x^{1/k})
    Computed using the Gram series for speed.
    R(x) = 1 + sum_{k=1}^{inf} (ln x)^k / (k * k! * zeta(k+1))
    """
    if x <= 1:
        return 0.0
    lnx = math.log(x)
    # Use the Gram series representation
    # R(x) = 1 + sum_{n=1}^{N} (ln x)^n / (n * n! * zeta(n+1))
    # Precomputed zeta values for small arguments aren't easily available,
    # so use the Mobius function approach with limited terms
    # mu(1)=1, mu(2)=-1, mu(3)=-1, mu(4)=0, mu(5)=-1, mu(6)=1, ...
    mu = [0, 1, -1, -1, 0, -1, 1, -1, 0, 0, 1,
          -1, 0, -1, 1, 1, 0, -1, 0, -1, 0,
          1, 1, -1, 0, 0, 1, 0, 0, -1, -1,
          -1, 0, 1, 1, 1, 0, -1, 1, 1, 0,
          -1, -1, -1, 0, 0, -1, -1, 0, 0, 1]
    result = 0.0
    for k in range(1, len(mu)):
        if mu[k] == 0:
            continue
        xk = x ** (1.0/k)
        if xk <= 1.0001:
            continue
        result += mu[k] / k * li_func(xk)
    return result

def inverse_R_fast(n):
    """Fast R^{-1}(n): solve R(x) = n via Newton's method."""
    if n <= 0:
        return 2.0
    x = max(2.1, n * math.log(max(n, 2)))
    for _ in range(80):
        val = R_func(x)
        deriv = 1.0 / math.log(x) if x > 1 else 1.0
        dx = (n - val) / deriv
        # Dampen large steps
        if abs(dx) > x * 0.5:
            dx = 0.5 * x * (1 if dx > 0 else -1)
        x = x + dx
        if x < 2:
            x = 2.1
        if abs(dx) < 1e-10:
            break
    return x

def cipolla_approx(n):
    """Cipolla's asymptotic expansion for p(n)."""
    if n < 6:
        return [2, 3, 5, 7, 11][n-1]
    ln_n = math.log(n)
    ln_ln_n = math.log(ln_n)
    return n * (ln_n + ln_ln_n - 1.0 + (ln_ln_n - 2.0) / ln_n)


N_MAX = 100000
print(f"Generating primes up to cover p({N_MAX})...")
t0 = time.time()

primes_list = sieve_primes(1500000)
p = primes_list[:N_MAX]
assert len(p) == N_MAX, f"Only got {len(p)} primes"
print(f"  Sieve done in {time.time()-t0:.2f}s. p({N_MAX}) = {p[-1]}")

# ─── Compute approximations (fast) ───
print("Computing all three approximations for n=1..100000...")
t0 = time.time()

delta_li = np.zeros(N_MAX)
delta_R = np.zeros(N_MAX)
delta_cip = np.zeros(N_MAX)

for i in range(N_MAX):
    n = i + 1
    pn = p[i]
    delta_cip[i] = pn - cipolla_approx(n)
    if n <= 1:
        delta_li[i] = pn - 2.0
        delta_R[i] = pn - 2.0
    else:
        delta_li[i] = pn - inverse_li_fast(n)
        delta_R[i] = pn - inverse_R_fast(n)

    if (i+1) % 20000 == 0:
        elapsed = time.time() - t0
        print(f"  n={i+1}/{N_MAX} done ({elapsed:.1f}s)")

print(f"All approximations computed in {time.time()-t0:.1f}s")

# Save raw deltas
np.savez('/apps/aplikacijos/prime-research/session6_experiments/deltas.npz',
         delta_li=delta_li, delta_R=delta_R, delta_cip=delta_cip,
         primes=np.array(p))

results = []
results.append("# Session 6: Wavelet Analysis of δ(n) = p(n) - approx(n)")
results.append("")
results.append(f"N_MAX = {N_MAX}, p({N_MAX}) = {p[-1]}")
results.append("")

# ═══════════════════════════════════════════════════════════════
# PART 1: Basic Statistics
# ═══════════════════════════════════════════════════════════════

results.append("## Part 1: Basic Statistics of δ(n)")
results.append("")
for name, delta in [("Li^{-1}", delta_li), ("R^{-1}", delta_R), ("Cipolla", delta_cip)]:
    results.append(f"### {name}")
    results.append(f"  Mean:    {np.mean(delta):.4f}")
    results.append(f"  Std:     {np.std(delta):.4f}")
    results.append(f"  Min:     {np.min(delta):.4f}")
    results.append(f"  Max:     {np.max(delta):.4f}")
    results.append(f"  Median:  {np.median(delta):.4f}")
    # Skewness and kurtosis
    results.append(f"  Skew:    {stats.skew(delta):.4f}")
    results.append(f"  Kurt:    {stats.kurtosis(delta):.4f}")
    results.append("")

# ═══════════════════════════════════════════════════════════════
# PART 2: Power Spectrum Analysis
# ═══════════════════════════════════════════════════════════════

results.append("## Part 2: Power Spectrum Analysis")
results.append("")

for name, delta in [("Li^{-1}", delta_li), ("R^{-1}", delta_R), ("Cipolla", delta_cip)]:
    # Detrend first
    detrended = signal.detrend(delta)

    # Compute power spectral density using Welch's method
    freqs, psd = signal.welch(detrended, fs=1.0, nperseg=min(8192, len(delta)//4),
                               noverlap=None, scaling='density')

    # Fit power law: PSD ~ f^{-α} (in log-log space)
    # Skip DC and very low frequencies
    mask = freqs > 0.001
    log_f = np.log10(freqs[mask])
    log_psd = np.log10(psd[mask])

    # Linear fit in log-log space
    slope, intercept, r_value, p_value, std_err = stats.linregress(log_f, log_psd)
    alpha = -slope

    results.append(f"### {name} Power Spectrum")
    results.append(f"  Spectral slope α = {alpha:.4f} (PSD ~ f^{{-{alpha:.4f}}})")
    results.append(f"  R² of fit: {r_value**2:.6f}")
    results.append(f"  Interpretation:")
    if abs(alpha) < 0.3:
        results.append(f"    α ≈ 0 → WHITE NOISE (no exploitable frequency structure)")
    elif 0.3 <= alpha < 0.8:
        results.append(f"    α ≈ 0.5 → 1/f^0.5 (weakly correlated, fractional noise)")
    elif 0.8 <= alpha < 1.3:
        results.append(f"    α ≈ 1 → 1/f NOISE (pink noise, long-range correlations)")
    elif 1.3 <= alpha < 1.8:
        results.append(f"    α ≈ 1.5 → Between pink and brown noise")
    elif 1.8 <= alpha < 2.3:
        results.append(f"    α ≈ 2 → BROWNIAN / RED NOISE (random walk!)")
    else:
        results.append(f"    α = {alpha:.2f} → unusual spectral slope")
    results.append("")

    # Check for spectral peaks (periodic components)
    peak_indices, peak_props = signal.find_peaks(psd, height=np.median(psd)*5,
                                                  prominence=np.median(psd)*3)
    if len(peak_indices) > 0:
        results.append(f"  Spectral peaks found at frequencies:")
        for pi_idx in peak_indices[:10]:
            f_peak = freqs[pi_idx]
            period = 1.0/f_peak if f_peak > 0 else float('inf')
            results.append(f"    f={f_peak:.6f} (period={period:.2f}), power={psd[pi_idx]:.2f}")
        results.append("")
    else:
        results.append(f"  No significant spectral peaks found.")
        results.append("")

# ═══════════════════════════════════════════════════════════════
# PART 3: Wavelet Decomposition
# ═══════════════════════════════════════════════════════════════

results.append("## Part 3: Wavelet Decomposition")
results.append("")

# Use the R^{-1} correction as our primary target (best approximation)
delta_main = delta_R.copy()

for wavelet_name in ['db4', 'haar', 'sym8']:
    results.append(f"### Wavelet: {wavelet_name}")

    # Multi-level decomposition
    max_level = pywt.dwt_max_level(len(delta_main), wavelet_name)
    levels_to_use = min(max_level, 16)

    coeffs = pywt.wavedec(delta_main, wavelet_name, level=levels_to_use)

    # Energy at each level
    total_energy = np.sum(delta_main**2)
    results.append(f"  Total signal energy: {total_energy:.2f}")
    results.append(f"  Decomposition levels: {levels_to_use}")
    results.append(f"  Level | Coeffs | Energy | % of Total | Scale (approx n)")
    results.append(f"  ------|--------|--------|------------|------------------")

    for level_idx, c in enumerate(coeffs):
        energy = np.sum(c**2)
        pct = 100.0 * energy / total_energy if total_energy > 0 else 0
        # Scale: level 0 = approximation (largest scale), then details from coarse to fine
        if level_idx == 0:
            scale_desc = f"~{N_MAX} (approximation)"
        else:
            scale = N_MAX / (2**level_idx)
            scale_desc = f"~{scale:.0f}"
        results.append(f"  {level_idx:5d} | {len(c):6d} | {energy:12.2f} | {pct:8.4f}%  | {scale_desc}")
    results.append("")

    # Reconstruct from only the coarsest levels and measure residual
    for keep_levels in [1, 2, 3, 5, 8]:
        if keep_levels > levels_to_use:
            continue
        # Zero out detail coefficients beyond keep_levels
        coeffs_truncated = [coeffs[0]] + [np.zeros_like(c) if i >= keep_levels
                                            else c for i, c in enumerate(coeffs[1:], 0)]
        # Actually: coeffs[0] = approx, coeffs[1] = coarsest detail, coeffs[-1] = finest
        # Keep approx + first (keep_levels-1) detail levels
        coeffs_truncated = list(coeffs)
        for i in range(keep_levels, len(coeffs_truncated)):
            coeffs_truncated[i] = np.zeros_like(coeffs_truncated[i])

        reconstructed = pywt.waverec(coeffs_truncated, wavelet_name)[:N_MAX]
        residual = delta_main - reconstructed
        rmse = np.sqrt(np.mean(residual**2))
        r2 = 1.0 - np.sum(residual**2) / np.sum((delta_main - np.mean(delta_main))**2)
        results.append(f"  Keep {keep_levels} levels: RMSE={rmse:.4f}, R²={r2:.6f}")
    results.append("")

# ═══════════════════════════════════════════════════════════════
# PART 4: Wavelet Coefficient Statistics per Scale
# ═══════════════════════════════════════════════════════════════

results.append("## Part 4: Wavelet Coefficient Distribution Analysis")
results.append("")
results.append("Testing if wavelet coefficients at each scale are Gaussian (random)")
results.append("vs. structured (heavy-tailed, skewed, etc.)")
results.append("")

coeffs = pywt.wavedec(delta_R, 'db4', level=14)
for level_idx in range(1, min(len(coeffs), 15)):
    c = coeffs[level_idx]
    if len(c) < 20:
        continue
    # Shapiro-Wilk test for normality (on subsample if too large)
    sample = c[:5000] if len(c) > 5000 else c
    if len(sample) >= 8:
        sw_stat, sw_p = stats.shapiro(sample)
    else:
        sw_stat, sw_p = 0, 0

    # Anderson-Darling
    ad_result = stats.anderson(sample, dist='norm')

    # Excess kurtosis (0 = Gaussian)
    kurt = stats.kurtosis(c)
    skew = stats.skew(c)

    scale = N_MAX / (2**level_idx)
    results.append(f"  Level {level_idx} (scale ~{scale:.0f}): "
                   f"kurt={kurt:.3f}, skew={skew:.3f}, "
                   f"Shapiro p={sw_p:.2e}, "
                   f"{'GAUSSIAN' if sw_p > 0.05 else 'NON-GAUSSIAN'}")

results.append("")

# ═══════════════════════════════════════════════════════════════
# PART 5: Hurst Exponent / DFA
# ═══════════════════════════════════════════════════════════════

results.append("## Part 5: Hurst Exponent and Long-Range Dependence")
results.append("")

def hurst_rs(ts, min_window=10, max_window=None):
    """Estimate Hurst exponent via R/S analysis."""
    n = len(ts)
    if max_window is None:
        max_window = n // 4

    window_sizes = []
    rs_values = []

    ws = min_window
    while ws <= max_window:
        window_sizes.append(ws)
        rs_list = []
        for start in range(0, n - ws + 1, ws):
            segment = ts[start:start+ws]
            mean_seg = np.mean(segment)
            cumdev = np.cumsum(segment - mean_seg)
            R = np.max(cumdev) - np.min(cumdev)
            S = np.std(segment, ddof=1)
            if S > 0:
                rs_list.append(R / S)
        if rs_list:
            rs_values.append(np.mean(rs_list))
        ws = int(ws * 1.5)

    log_ws = np.log(window_sizes[:len(rs_values)])
    log_rs = np.log(rs_values)
    slope, intercept, r, p, se = stats.linregress(log_ws, log_rs)
    return slope, r**2

def dfa(ts, min_window=10, max_window=None, order=1):
    """Detrended Fluctuation Analysis."""
    n = len(ts)
    if max_window is None:
        max_window = n // 4

    # Cumulative sum (profile)
    profile = np.cumsum(ts - np.mean(ts))

    window_sizes = []
    fluctuations = []

    ws = min_window
    while ws <= max_window:
        window_sizes.append(ws)
        n_windows = n // ws
        if n_windows == 0:
            break
        fluct_sq = []
        for i in range(n_windows):
            segment = profile[i*ws:(i+1)*ws]
            x = np.arange(ws)
            coeffs_poly = np.polyfit(x, segment, order)
            trend = np.polyval(coeffs_poly, x)
            fluct_sq.append(np.mean((segment - trend)**2))
        fluctuations.append(np.sqrt(np.mean(fluct_sq)))
        ws = int(ws * 1.3)

    log_ws = np.log(window_sizes[:len(fluctuations)])
    log_f = np.log(fluctuations)
    slope, intercept, r, p, se = stats.linregress(log_ws, log_f)
    return slope, r**2

for name, delta in [("Li^{-1}", delta_li), ("R^{-1}", delta_R), ("Cipolla", delta_cip)]:
    H_rs, r2_rs = hurst_rs(delta)
    H_dfa, r2_dfa = dfa(delta)

    results.append(f"### {name}")
    results.append(f"  Hurst (R/S):  H = {H_rs:.4f} (R² = {r2_rs:.4f})")
    results.append(f"  Hurst (DFA):  H = {H_dfa:.4f} (R² = {r2_dfa:.4f})")

    if H_dfa < 0.4:
        interp = "ANTI-PERSISTENT (mean-reverting)"
    elif H_dfa < 0.6:
        interp = "UNCORRELATED (random walk / white noise)"
    elif H_dfa < 0.8:
        interp = "WEAKLY PERSISTENT (some long-range correlation)"
    else:
        interp = "STRONGLY PERSISTENT (highly correlated, trending)"
    results.append(f"  Interpretation: {interp}")
    results.append("")

# ═══════════════════════════════════════════════════════════════
# PART 6: Chebyshev Polynomial Expansion
# ═══════════════════════════════════════════════════════════════

results.append("## Part 6: Chebyshev Polynomial Expansion")
results.append("")

# Map n to [-1, 1] for Chebyshev -- subsample for speed
CHEB_SUBSAMPLE = 10  # Use every 10th point
idx_cheb = np.arange(0, N_MAX, CHEB_SUBSAMPLE)
n_cheb = len(idx_cheb)
x_cheb = 2.0 * idx_cheb / N_MAX - 1.0
x_cheb_full = 2.0 * np.arange(N_MAX) / N_MAX - 1.0

for name, delta in [("R^{-1}", delta_R)]:
    results.append(f"### {name} (fitted on {n_cheb} subsampled points, evaluated on all {N_MAX})")
    delta_sub = delta[idx_cheb]

    for degree in [5, 10, 20, 50, 100, 200]:
        coeffs_cheb = np.polynomial.chebyshev.chebfit(x_cheb, delta_sub, degree)
        fitted_full = np.polynomial.chebyshev.chebval(x_cheb_full, coeffs_cheb)
        residual = delta - fitted_full
        rmse = np.sqrt(np.mean(residual**2))
        r2 = 1.0 - np.sum(residual**2) / np.sum((delta - np.mean(delta))**2)
        results.append(f"  Degree {degree:4d}: RMSE={rmse:.4f}, R²={r2:.6f}")

    # How does RMSE scale with degree?
    degrees = list(range(5, 201, 5))
    rmses = []
    for deg in degrees:
        c = np.polynomial.chebyshev.chebfit(x_cheb, delta_sub, deg)
        f = np.polynomial.chebyshev.chebval(x_cheb_full, c)
        rmses.append(np.sqrt(np.mean((delta - f)**2)))

    # Fit power law: RMSE ~ degree^{-beta}
    log_d = np.log(degrees)
    log_r = np.log(rmses)
    slope, intercept, r, p_val, se = stats.linregress(log_d, log_r)
    results.append(f"  RMSE decay rate: RMSE ~ degree^{{{slope:.3f}}} (R²={r**2:.4f})")
    results.append(f"  (For smooth functions expect fast decay; for noise expect slow)")
    results.append("")

# ═══════════════════════════════════════════════════════════════
# PART 7: Periodicity Search Against Known Constants
# ═══════════════════════════════════════════════════════════════

results.append("## Part 7: Periodicity Search Against Known Constants")
results.append("")

# Compute autocorrelation
def autocorrelation(x, max_lag):
    """Compute normalized autocorrelation using FFT (fast)."""
    x_centered = x - np.mean(x)
    n = len(x)
    var = np.var(x)
    if var == 0:
        return np.zeros(max_lag)
    # FFT-based autocorrelation
    fft_x = np.fft.rfft(x_centered, n=2*n)
    acf_full = np.fft.irfft(np.abs(fft_x)**2)[:max_lag]
    acf_full = acf_full / (var * n)
    return acf_full

# Known constants to test as periods
import math
known_constants = {
    'ln(2)': math.log(2),           # 0.6931...
    'π': math.pi,                    # 3.1416...
    'π²/6': math.pi**2/6,           # ζ(2) = 1.6449...
    'e': math.e,                     # 2.7183...
    'ln(10)': math.log(10),          # 2.3026...
    'sqrt(2)': math.sqrt(2),         # 1.4142...
    'γ (Euler)': 0.5772156649,
    '1/ln(2)': 1/math.log(2),       # 1.4427...
    '2π': 2*math.pi,                 # 6.2832...
    'ζ zero 1 (14.13...)': 14.134725,
    'ζ zero 2 (21.02...)': 21.022040,
    'ζ zero 3 (25.01...)': 25.010858,
    '2π/ln(2)': 2*math.pi/math.log(2),  # 9.0647...
    '2π/ln(3)': 2*math.pi/math.log(3),  # 5.7187...
}

# For each constant, check if the FFT has power at frequency 1/constant
# or if autocorrelation shows peaks at integer multiples

delta_detrended = signal.detrend(delta_R)
fft_vals = np.fft.rfft(delta_detrended)
fft_power = np.abs(fft_vals)**2
fft_freqs = np.fft.rfftfreq(len(delta_detrended), d=1.0)

# Background power level (median)
bg_power = np.median(fft_power[1:])

results.append("Testing if FFT power at f=1/constant is significantly above background:")
results.append("")

for const_name, const_val in sorted(known_constants.items()):
    # Test frequencies related to this constant
    test_freqs = [1.0/const_val, const_val/N_MAX, 1.0/(const_val * math.log(2))]

    for tf in [1.0/const_val]:
        if tf >= 0.5:  # Nyquist
            continue
        # Find nearest FFT bin
        idx = np.argmin(np.abs(fft_freqs - tf))
        if idx == 0:
            continue
        power_at_freq = fft_power[idx]
        # Also check a neighborhood
        neighborhood = slice(max(1, idx-2), min(len(fft_power), idx+3))
        max_nearby = np.max(fft_power[neighborhood])
        ratio = max_nearby / bg_power
        sig = "***" if ratio > 10 else "**" if ratio > 5 else "*" if ratio > 3 else ""
        results.append(f"  {const_name:25s}: f={tf:.6f}, power/bg={ratio:.2f} {sig}")

results.append("")

# Autocorrelation analysis
results.append("Autocorrelation analysis (R^{-1} correction):")
results.append("")
max_lag = 1000
acf = autocorrelation(delta_R, max_lag)
# Find peaks in autocorrelation
acf_peaks, _ = signal.find_peaks(acf[1:], height=0.02, distance=3)
acf_peaks += 1  # Offset for skipping lag=0

if len(acf_peaks) > 0:
    results.append(f"  Autocorrelation peaks (lag, acf value):")
    for peak_lag in acf_peaks[:20]:
        results.append(f"    lag={peak_lag}, acf={acf[peak_lag]:.4f}")
else:
    results.append(f"  No significant autocorrelation peaks found.")
results.append("")

# ═══════════════════════════════════════════════════════════════
# PART 8: Multiscale Entropy Analysis
# ═══════════════════════════════════════════════════════════════

results.append("## Part 8: Multiscale Analysis - Variance Scaling")
results.append("")
results.append("If δ(n) is white noise: Var(mean of k samples) ~ 1/k")
results.append("If δ(n) has structure: different scaling")
results.append("")

for name, delta in [("R^{-1}", delta_R)]:
    block_sizes = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000]
    variances = []

    for bs in block_sizes:
        n_blocks = len(delta) // bs
        if n_blocks < 5:
            continue
        block_means = [np.mean(delta[i*bs:(i+1)*bs]) for i in range(n_blocks)]
        variances.append((bs, np.var(block_means)))

    results.append(f"### {name} Variance Scaling")
    results.append(f"  Block Size | Var(block mean) | Ratio to white noise")
    results.append(f"  -----------|-----------------|--------------------")

    base_var = variances[0][1] if variances else 1
    for bs, v in variances:
        expected_white = base_var / bs  # White noise prediction
        ratio = v / expected_white if expected_white > 0 else 0
        results.append(f"  {bs:10d} | {v:15.6f} | {ratio:.4f}")

    # Fit scaling: Var ~ k^{-γ}
    if len(variances) > 3:
        log_bs = np.log([v[0] for v in variances])
        log_var = np.log([v[1] for v in variances])
        slope, intercept, r, p_val, se = stats.linregress(log_bs, log_var)
        results.append(f"  Variance scaling: Var ~ k^{{{slope:.4f}}} (R²={r**2:.4f})")
        results.append(f"  White noise would give slope = -1.0")
        results.append(f"  Random walk would give slope = -1.0 with different prefactor")
        results.append(f"  Persistent process: slope > -1.0")
    results.append("")

# ═══════════════════════════════════════════════════════════════
# PART 9: Wavelet Scalogram - Energy Distribution
# ═══════════════════════════════════════════════════════════════

results.append("## Part 9: Continuous Wavelet Transform Scalogram Summary")
results.append("")

# Use CWT with Ricker (Mexican hat) wavelet at various scales
scales = np.logspace(0.5, 4, 50)  # scales from ~3 to 10000
# Subsample for speed
subsample = delta_R[::10]  # every 10th point
cwt_matrix, cwt_freqs = pywt.cwt(subsample, scales[:30], 'mexh')

# Energy at each scale
results.append(f"  Scale | Energy Fraction | Interpretation")
results.append(f"  ------|----------------|---------------")
total_cwt_energy = np.sum(cwt_matrix**2)
for i, s in enumerate(scales[:30]):
    energy = np.sum(cwt_matrix[i]**2)
    pct = 100.0 * energy / total_cwt_energy
    real_scale = s * 10  # Account for subsampling
    results.append(f"  {real_scale:8.1f} | {pct:12.4f}%   | {'<-- HIGH' if pct > 8 else ''}")

results.append("")

# ═══════════════════════════════════════════════════════════════
# PART 10: Information-Theoretic Measures
# ═══════════════════════════════════════════════════════════════

results.append("## Part 10: Predictability Test (Conditional Entropy)")
results.append("")
results.append("Can δ(n) be predicted from δ(n-1), δ(n-2), ...?")
results.append("")

# Bin the corrections and compute conditional entropy
def binned_entropy(data, n_bins=50):
    """Compute entropy of binned data."""
    counts, _ = np.histogram(data, bins=n_bins)
    probs = counts / len(data)
    probs = probs[probs > 0]
    return -np.sum(probs * np.log2(probs))

for name, delta in [("R^{-1}", delta_R)]:
    H_marginal = binned_entropy(delta)

    # Conditional: bin (δ(n), δ(n-1)) jointly
    n_bins = 50
    for lag in [1, 2, 5, 10]:
        joint = np.column_stack([delta[lag:], delta[:-lag]])
        H_joint, _, _ = np.histogram2d(joint[:,0], joint[:,1], bins=n_bins)
        H_joint = H_joint / H_joint.sum()
        H_joint = H_joint[H_joint > 0]
        H_xy = -np.sum(H_joint * np.log2(H_joint))

        # H(Y|X) = H(X,Y) - H(X)
        H_cond = H_xy - binned_entropy(delta[:-lag], n_bins)
        mi = H_marginal - H_cond  # Mutual information

        results.append(f"  {name} lag-{lag}: H(δ)={H_marginal:.3f} bits, "
                       f"H(δ|δ_{{-{lag}}})={H_cond:.3f} bits, MI={mi:.4f} bits")
    results.append("")

# ═══════════════════════════════════════════════════════════════
# PART 11: Subtracting Smooth Component - Is Residual White?
# ═══════════════════════════════════════════════════════════════

results.append("## Part 11: After Removing Smooth Component")
results.append("")

# Try subtracting various smooth fits and check if residual is white noise
for name, delta in [("R^{-1}", delta_R)]:
    ns = np.arange(1, N_MAX+1, dtype=float)

    # Fit: δ(n) ≈ a * sqrt(n) * ln(n) + b * sqrt(n) + c * ln(n) + d
    # This is motivated by the known error term of R^{-1}
    X = np.column_stack([
        np.sqrt(ns) * np.log(ns),
        np.sqrt(ns),
        np.log(ns),
        np.ones(N_MAX),
        np.log(ns)**2,
        ns**0.25,
    ])

    # Least squares fit
    coeffs_fit, residuals, rank, sv = np.linalg.lstsq(X, delta, rcond=None)
    smooth_fit = X @ coeffs_fit
    residual = delta - smooth_fit

    results.append(f"### {name} with smooth component removed")
    results.append(f"  Fit: a*sqrt(n)*ln(n) + b*sqrt(n) + c*ln(n) + d + e*ln²(n) + f*n^0.25")
    results.append(f"  Coefficients: {[f'{c:.6f}' for c in coeffs_fit]}")
    results.append(f"  Residual RMSE: {np.sqrt(np.mean(residual**2)):.4f}")
    results.append(f"  Residual as fraction of δ std: {np.std(residual)/np.std(delta):.4f}")
    results.append("")

    # Is the residual white noise?
    # 1. Ljung-Box test (via autocorrelation)
    acf_resid = autocorrelation(residual, 50)
    n_r = len(residual)
    # Q statistic
    Q = n_r * (n_r + 2) * np.sum(acf_resid[1:50]**2 / (n_r - np.arange(1, 50)))
    # Under H0 (white noise), Q ~ chi²(48)
    p_ljung = 1 - stats.chi2.cdf(Q, 48)

    results.append(f"  Ljung-Box Q={Q:.2f}, p={p_ljung:.2e} ({'WHITE NOISE' if p_ljung > 0.05 else 'NOT white noise'})")

    # 2. Runs test
    median_r = np.median(residual)
    runs = 1 + np.sum(np.diff((residual > median_r).astype(int)) != 0)
    n_above = np.sum(residual > median_r)
    n_below = np.sum(residual <= median_r)
    expected_runs = 1 + 2*n_above*n_below/(n_above+n_below)
    var_runs = 2*n_above*n_below*(2*n_above*n_below - n_above - n_below)/((n_above+n_below)**2 * (n_above+n_below-1))
    z_runs = (runs - expected_runs) / np.sqrt(var_runs)
    p_runs = 2 * (1 - stats.norm.cdf(abs(z_runs)))

    results.append(f"  Runs test: z={z_runs:.4f}, p={p_runs:.4e} ({'RANDOM' if p_runs > 0.05 else 'STRUCTURED'})")
    results.append("")

    # Power spectrum of residual
    freqs_r, psd_r = signal.welch(residual, fs=1.0, nperseg=min(8192, len(residual)//4))
    mask_r = freqs_r > 0.001
    slope_r, _, r_r, _, _ = stats.linregress(np.log10(freqs_r[mask_r]), np.log10(psd_r[mask_r]))
    results.append(f"  Residual spectral slope: α = {-slope_r:.4f}")
    results.append("")

# ═══════════════════════════════════════════════════════════════
# PART 12: Cross-Scale Correlations (Wavelet Modulus Maxima)
# ═══════════════════════════════════════════════════════════════

results.append("## Part 12: Wavelet Cross-Scale Structure")
results.append("")

coeffs_full = pywt.wavedec(delta_R, 'db4', level=14)
# Check if adjacent scales have correlated coefficient magnitudes
results.append("Correlation between |coefficients| at adjacent wavelet scales:")
for i in range(1, min(len(coeffs_full)-1, 13)):
    c1 = np.abs(coeffs_full[i])
    c2 = np.abs(coeffs_full[i+1])
    # Resample to same length
    min_len = min(len(c1), len(c2))
    c1_rs = signal.resample(c1, min_len)
    c2_rs = signal.resample(c2, min_len)
    corr = np.corrcoef(c1_rs, c2_rs)[0,1]
    results.append(f"  Level {i} vs {i+1}: r = {corr:.4f}")

results.append("")

# ═══════════════════════════════════════════════════════════════
# FINAL VERDICT
# ═══════════════════════════════════════════════════════════════

results.append("## CONCLUSIONS")
results.append("")
results.append("### Summary of Findings")
results.append("")
results.append("(Auto-generated based on numerical results above)")
results.append("")

# Write results
output_path = '/apps/aplikacijos/prime-research/session6_experiments/wavelet_results.md'
with open(output_path, 'w') as f:
    f.write('\n'.join(results))
print(f"\nResults written to {output_path}")
print("\n".join(results[-40:]))

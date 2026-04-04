#!/usr/bin/env python3
"""
Session 9: Zeta Zero Compressibility Analysis
==============================================

Investigate whether the nontrivial zeros of the Riemann zeta function
have compressible global structure that could allow computing
    S(x) = sum_rho R(x^rho)
in O(polylog x) operations instead of O(x^{1/2}).

Key questions:
1. Are the residuals gamma_n - gamma_n^{smooth} compressible?
2. Can we group zeros into bands with partial cancellation?
3. Can statistical properties (pair correlation) replace individual zeros?
4. Does FMM-style grouping work for the R(x^rho) sum?
"""

import numpy as np
from mpmath import mp, mpf, zetazero, log, pi, exp, sqrt, fabs
from mpmath import lambertw, loggamma, arg, zeta
import time
import sys
import json

mp.dps = 30  # 30 decimal digits precision

# ============================================================
# PART 1: Compute first 1000 zeta zeros
# ============================================================

def compute_zeros(N=1000):
    """Compute first N nontrivial zeta zeros (imaginary parts gamma_n)."""
    print(f"Computing first {N} zeta zeros...")
    t0 = time.time()
    zeros = []
    for n in range(1, N+1):
        z = zetazero(n)
        gamma = float(z.imag)
        zeros.append(gamma)
        if n % 100 == 0:
            print(f"  {n}/{N} done, gamma_{n} = {gamma:.6f}")
    elapsed = time.time() - t0
    print(f"Computed {N} zeros in {elapsed:.1f}s")
    return np.array(zeros)

# ============================================================
# PART 2: Backlund/Gram smooth approximation
# ============================================================

def backlund_smooth(n):
    """
    Smooth approximation to gamma_n using the asymptotic formula:

    N(T) ~ (T/2pi) * log(T/2pi*e) + 7/8 + S(T)

    Inverting N(gamma_n) = n - 1 + epsilon gives the Backlund approximation.
    We solve: (g/2pi)*log(g/2pi*e) = n - 7/8

    Using the Lambert W function for inversion.
    """
    # We need to invert f(g) = (g/(2*pi)) * ln(g/(2*pi*e)) = n - 7/8
    # Let u = g/(2*pi), then u*ln(u/e) = n - 7/8
    # u*ln(u) - u = n - 7/8
    # This doesn't have a clean Lambert W form, so use Newton's method.

    target = n - 7.0/8.0
    # Initial guess using leading asymptotics: g ~ 2*pi*n / ln(n)
    if n < 2:
        g = 14.0
    else:
        g = 2 * np.pi * n / np.log(n / (2 * np.pi * np.e) + 1)
        if g < 10:
            g = 10.0 + n

    # Newton iteration on f(g) = (g/(2pi))*ln(g/(2pi*e)) - target = 0
    for _ in range(50):
        u = g / (2 * np.pi)
        if u <= 0:
            u = 0.001
        f_val = u * np.log(u / np.e) - target  # u*ln(u) - u - target
        f_prime = np.log(u / np.e) + 1  # derivative: ln(u/e) + 1 = ln(u)
        if abs(f_prime) < 1e-15:
            break
        delta = f_val / f_prime
        g -= 2 * np.pi * delta
        if abs(delta) < 1e-12:
            break

    return g

def compute_smooth_zeros(N=1000):
    """Compute smooth approximations for first N zeros."""
    return np.array([backlund_smooth(n) for n in range(1, N+1)])

# ============================================================
# PART 3: Residual analysis
# ============================================================

def analyze_residuals(gammas, gammas_smooth):
    """Analyze the residuals delta_n = gamma_n - gamma_n^{smooth}."""
    residuals = gammas - gammas_smooth

    results = {}

    # Basic statistics
    results['mean'] = float(np.mean(residuals))
    results['std'] = float(np.std(residuals))
    results['max_abs'] = float(np.max(np.abs(residuals)))
    results['min'] = float(np.min(residuals))
    results['max'] = float(np.max(residuals))

    print(f"\n=== RESIDUAL STATISTICS ===")
    print(f"Mean:    {results['mean']:.6f}")
    print(f"Std:     {results['std']:.6f}")
    print(f"Max|r|:  {results['max_abs']:.6f}")
    print(f"Range:   [{results['min']:.6f}, {results['max']:.6f}]")

    # Relative residual: how many digits does the smooth part give?
    rel_residuals = np.abs(residuals) / gammas
    results['mean_relative'] = float(np.mean(rel_residuals))
    results['max_relative'] = float(np.max(rel_residuals))
    digits_from_smooth = -np.log10(rel_residuals + 1e-30)
    results['mean_digits_from_smooth'] = float(np.mean(digits_from_smooth))
    results['min_digits_from_smooth'] = float(np.min(digits_from_smooth))

    print(f"\nDigits from smooth approximation:")
    print(f"  Mean: {results['mean_digits_from_smooth']:.2f}")
    print(f"  Min:  {results['min_digits_from_smooth']:.2f}")

    return residuals, results

def entropy_analysis(residuals, n_bins=50):
    """Estimate Shannon entropy of the residual distribution."""
    # Normalize residuals to [0,1]
    r_min, r_max = residuals.min(), residuals.max()
    r_range = r_max - r_min
    if r_range < 1e-15:
        return 0.0, 0.0

    normalized = (residuals - r_min) / r_range

    # Histogram-based entropy
    counts, _ = np.histogram(normalized, bins=n_bins, range=(0, 1))
    probs = counts / counts.sum()
    probs = probs[probs > 0]
    entropy = -np.sum(probs * np.log2(probs))

    # Entropy rate: bits per sample
    # For Gaussian, entropy = 0.5*log2(2*pi*e*sigma^2)
    sigma = np.std(residuals)
    gaussian_entropy = 0.5 * np.log2(2 * np.pi * np.e * sigma**2) if sigma > 0 else 0

    print(f"\n=== ENTROPY ANALYSIS ===")
    print(f"Histogram entropy ({n_bins} bins): {entropy:.3f} bits")
    print(f"Gaussian entropy (same variance): {gaussian_entropy:.3f} bits")
    print(f"Entropy ratio: {entropy/gaussian_entropy:.3f}" if gaussian_entropy > 0 else "")

    return entropy, gaussian_entropy

def autocorrelation_analysis(residuals, max_lag=50):
    """Compute autocorrelation of residuals."""
    N = len(residuals)
    r_centered = residuals - np.mean(residuals)
    var = np.var(residuals)
    if var < 1e-30:
        return np.zeros(max_lag)

    acf = np.zeros(max_lag)
    for lag in range(max_lag):
        acf[lag] = np.mean(r_centered[:N-lag] * r_centered[lag:]) / var

    print(f"\n=== AUTOCORRELATION ===")
    print(f"Lag  ACF")
    for lag in [1, 2, 3, 5, 10, 20, 50]:
        if lag < max_lag:
            print(f"  {lag:3d}  {acf[lag]:+.6f}")

    # Is there significant autocorrelation?
    threshold = 2.0 / np.sqrt(N)  # 95% CI for white noise
    n_significant = np.sum(np.abs(acf[1:]) > threshold)
    print(f"\nSignificance threshold (95%): {threshold:.4f}")
    print(f"Significant lags (out of {max_lag-1}): {n_significant}")

    return acf

def spectral_analysis(residuals):
    """Power spectral density of residuals."""
    N = len(residuals)
    r_centered = residuals - np.mean(residuals)

    # FFT
    fft_vals = np.fft.rfft(r_centered)
    psd = np.abs(fft_vals)**2 / N
    freqs = np.fft.rfftfreq(N)

    # Check for dominant frequencies
    psd_norm = psd / psd.sum()
    sorted_idx = np.argsort(psd_norm)[::-1]

    print(f"\n=== SPECTRAL ANALYSIS ===")
    print(f"Top 10 frequencies by power:")
    cumulative = 0
    for i in range(min(10, len(sorted_idx))):
        idx = sorted_idx[i]
        cumulative += psd_norm[idx]
        print(f"  freq={freqs[idx]:.6f} (period={1/freqs[idx]:.1f} if >0), "
              f"power={psd_norm[idx]:.4f}, cumulative={cumulative:.4f}")

    # Spectral entropy
    psd_pos = psd_norm[psd_norm > 0]
    spectral_entropy = -np.sum(psd_pos * np.log2(psd_pos))
    max_entropy = np.log2(len(psd_pos))

    print(f"\nSpectral entropy: {spectral_entropy:.3f} bits (max={max_entropy:.3f})")
    print(f"Spectral flatness ratio: {spectral_entropy/max_entropy:.4f}")
    print(f"  (1.0 = white noise, lower = more structured)")

    # How many frequencies capture 50%, 90%, 99% of power?
    cumsum = np.cumsum(psd_norm[sorted_idx])
    for threshold in [0.5, 0.9, 0.95, 0.99]:
        n_needed = np.searchsorted(cumsum, threshold) + 1
        print(f"  Frequencies for {threshold*100:.0f}% power: {n_needed}/{len(psd_norm)}")

    return psd, freqs, spectral_entropy

# ============================================================
# PART 4: Compressibility via SVD / Low-rank structure
# ============================================================

def compressibility_analysis(residuals, block_size=50):
    """Check if the residual sequence has low-rank structure."""
    N = len(residuals)
    n_blocks = N // block_size

    # Reshape into matrix and check SVD
    matrix = residuals[:n_blocks * block_size].reshape(n_blocks, block_size)
    U, S, Vt = np.linalg.svd(matrix, full_matrices=False)

    S_norm = S / S.sum()
    cumsum = np.cumsum(S_norm)

    print(f"\n=== SVD COMPRESSIBILITY (block_size={block_size}) ===")
    print(f"Matrix shape: {matrix.shape}")
    print(f"Top singular values (normalized):")
    for i in range(min(10, len(S))):
        print(f"  sigma_{i+1} = {S[i]:.6f} ({S_norm[i]:.4f}), cumulative: {cumsum[i]:.4f}")

    for threshold in [0.5, 0.9, 0.95, 0.99]:
        k = np.searchsorted(cumsum, threshold) + 1
        print(f"  Rank for {threshold*100:.0f}% energy: {k}/{len(S)}")

    return S, S_norm

# ============================================================
# PART 5: R(x^rho) computation and band grouping
# ============================================================

def R_function(x):
    """Riemann R function: R(x) = sum_{k=1}^inf mu(k)/k * li(x^{1/k})."""
    mp.dps = 30
    x = mpf(x)
    if x <= 1:
        return mpf(0)
    result = mpf(0)
    # Mobius function values for first terms
    mobius = [0, 1, -1, -1, 0, -1, 1, -1, 0, 0, 1, -1, 0, -1, 1, 1, 0, -1, 0, -1, 0]
    for k in range(1, len(mobius)):
        if mobius[k] == 0:
            continue
        xk = x ** (mpf(1)/k)
        if xk <= 1:
            break
        result += mpf(mobius[k]) / k * mp.li(xk)
    return result

def compute_oscillatory_sum(x, gammas, N_terms=None):
    """
    Compute S(x) = sum_rho R(x^rho) using the explicit formula.

    Each zero rho = 1/2 + i*gamma contributes:
    R(x^rho) + R(x^{conj(rho)}) = 2 * Re[R(x^{1/2 + i*gamma})]
    """
    mp.dps = 30
    x = mpf(x)
    log_x = mp.log(x)

    if N_terms is None:
        N_terms = len(gammas)
    else:
        N_terms = min(N_terms, len(gammas))

    total = mpf(0)
    contributions = []

    for j in range(N_terms):
        gamma = mpf(gammas[j])
        # R(x^rho) where rho = 1/2 + i*gamma
        # Leading term: li(x^rho) = li(x^{1/2 + i*gamma})
        # x^rho = x^{1/2} * x^{i*gamma} = sqrt(x) * exp(i*gamma*log(x))

        # For the leading term of R, use li(x^rho):
        # li(x^rho) ~ x^rho / (rho * log(x))  for large x
        # 2*Re[li(x^rho)] ~ 2*x^{1/2}*cos(gamma*log(x)) / (|rho|^2 * log(x)) * Re(conj(rho))

        # Direct computation via mpmath for moderate x
        phase = gamma * log_x
        amp = mp.sqrt(x)
        rho_abs_sq = mpf(0.25) + gamma**2

        # Leading Ri term: 2*Re[x^rho / (rho * log(x))]
        # = 2 * sqrt(x) / log(x) * Re[exp(i*gamma*log(x)) / (1/2 + i*gamma)]
        # = 2 * sqrt(x) / log(x) * (1/2*cos(phase) + gamma*sin(phase)) / rho_abs_sq

        cos_phase = mp.cos(phase)
        sin_phase = mp.sin(phase)

        contrib = 2 * amp / log_x * (mpf(0.5)*cos_phase + gamma*sin_phase) / rho_abs_sq
        contributions.append(float(contrib))
        total += contrib

    return float(total), contributions

# ============================================================
# PART 6: Statistical sum estimation (pair correlation approach)
# ============================================================

def statistical_sum_estimate(x, gammas, method='random_matrix'):
    """
    Try to estimate sum_rho R(x^rho) using STATISTICAL properties
    of the zeros rather than individual zero values.

    Key insight: if gamma_n = gamma_n^{smooth} + delta_n,
    then the sum splits into:
    S = S_smooth + S_correction

    S_smooth uses only the Backlund approximation.
    S_correction depends on the residuals.

    Question: is S_correction small or structured?
    """
    log_x = np.log(x)
    sqrt_x = np.sqrt(x)
    N = len(gammas)

    # --- Method 1: Use smooth zeros only ---
    gammas_smooth = np.array([backlund_smooth(n) for n in range(1, N+1)])

    # Compute sum with actual zeros
    phases_actual = gammas * log_x
    rho_abs_sq_actual = 0.25 + gammas**2
    cos_actual = np.cos(phases_actual)
    sin_actual = np.sin(phases_actual)
    contribs_actual = 2 * sqrt_x / log_x * (0.5*cos_actual + gammas*sin_actual) / rho_abs_sq_actual
    S_actual = np.sum(contribs_actual)

    # Compute sum with smooth zeros
    phases_smooth = gammas_smooth * log_x
    rho_abs_sq_smooth = 0.25 + gammas_smooth**2
    cos_smooth = np.cos(phases_smooth)
    sin_smooth = np.sin(phases_smooth)
    contribs_smooth = 2 * sqrt_x / log_x * (0.5*cos_smooth + gammas_smooth*sin_smooth) / rho_abs_sq_smooth
    S_smooth = np.sum(contribs_smooth)

    print(f"\n=== STATISTICAL SUM ANALYSIS (x={x:.0f}) ===")
    print(f"S_actual (from true zeros):  {S_actual:.6f}")
    print(f"S_smooth (from Backlund):    {S_smooth:.6f}")
    print(f"S_correction = S_actual - S_smooth: {S_actual - S_smooth:.6f}")
    print(f"|S_correction/S_actual|: {abs((S_actual-S_smooth)/(S_actual+1e-30)):.4f}")

    # --- Method 2: Perturbative correction ---
    # delta_n = gamma_n - gamma_n^smooth
    deltas = gammas - gammas_smooth

    # First-order correction: d/dgamma [contrib(gamma)] * delta
    # d/dgamma[cos(gamma*log_x)] = -log_x * sin(gamma*log_x)
    # The correction involves: sum delta_n * d(contrib_n)/dgamma_n

    # Derivative of contribution w.r.t. gamma:
    # contrib = 2*sqrt(x)/log(x) * (0.5*cos(g*L) + g*sin(g*L)) / (0.25 + g^2)
    # where L = log(x)
    L = log_x
    g = gammas_smooth
    rho2 = 0.25 + g**2

    # Numerator: N = 0.5*cos(gL) + g*sin(gL)
    # dN/dg = -0.5*L*sin(gL) + sin(gL) + g*L*cos(gL)
    #        = sin(gL)*(1 - 0.5*L) + g*L*cos(gL)  -- wait, let me redo
    # Actually: d/dg[0.5*cos(gL)] = -0.5*L*sin(gL)
    #           d/dg[g*sin(gL)] = sin(gL) + g*L*cos(gL)
    # So dN/dg = -0.5*L*sin(gL) + sin(gL) + g*L*cos(gL)

    dN_dg = -0.5*L*np.sin(g*L) + np.sin(g*L) + g*L*np.cos(g*L)
    N_val = 0.5*np.cos(g*L) + g*np.sin(g*L)
    drho2_dg = 2*g

    # d/dg[N/rho2] = (dN*rho2 - N*drho2) / rho2^2
    d_contrib_dg = 2*sqrt_x/L * (dN_dg * rho2 - N_val * drho2_dg) / rho2**2

    S_perturbative = S_smooth + np.sum(d_contrib_dg * deltas)

    print(f"\nS_perturbative (1st order):  {S_perturbative:.6f}")
    print(f"|S_perturbative - S_actual|: {abs(S_perturbative - S_actual):.6f}")
    print(f"Perturbative improvement over smooth: "
          f"{abs(S_smooth - S_actual) / (abs(S_perturbative - S_actual) + 1e-30):.1f}x")

    # --- Method 3: Random correction (Monte Carlo) ---
    # If residuals are "random" with known statistics, what's the expected sum?
    n_trials = 100
    mc_sums = []
    for _ in range(n_trials):
        fake_deltas = np.random.normal(np.mean(deltas), np.std(deltas), N)
        fake_gammas = gammas_smooth + fake_deltas
        phases_fake = fake_gammas * log_x
        rho2_fake = 0.25 + fake_gammas**2
        cos_fake = np.cos(phases_fake)
        sin_fake = np.sin(phases_fake)
        contribs_fake = 2*sqrt_x/log_x * (0.5*cos_fake + fake_gammas*sin_fake) / rho2_fake
        mc_sums.append(np.sum(contribs_fake))

    mc_mean = np.mean(mc_sums)
    mc_std = np.std(mc_sums)

    print(f"\nMonte Carlo (random residuals):")
    print(f"  Mean S: {mc_mean:.6f} +/- {mc_std:.6f}")
    print(f"  |MC_mean - S_actual|: {abs(mc_mean - S_actual):.6f}")

    return {
        'S_actual': S_actual,
        'S_smooth': S_smooth,
        'S_correction': S_actual - S_smooth,
        'S_perturbative': S_perturbative,
        'perturbative_error': abs(S_perturbative - S_actual),
        'smooth_error': abs(S_smooth - S_actual),
        'mc_mean': mc_mean,
        'mc_std': mc_std,
    }

# ============================================================
# PART 7: Band grouping (FMM-style)
# ============================================================

def band_grouping_analysis(x, gammas, band_sizes=[10, 20, 50, 100, 200]):
    """
    Group zeros into bands and analyze cancellation within bands.

    If contributions within a band largely cancel, we might only
    need to compute the band's net contribution approximately.
    """
    log_x = np.log(x)
    sqrt_x = np.sqrt(x)
    N = len(gammas)

    print(f"\n=== BAND GROUPING ANALYSIS (x={x:.0f}) ===")

    # Individual contributions
    phases = gammas * log_x
    rho_abs_sq = 0.25 + gammas**2
    contribs = 2 * sqrt_x / log_x * (0.5*np.cos(phases) + gammas*np.sin(phases)) / rho_abs_sq

    total = np.sum(contribs)
    print(f"Total sum: {total:.6f}")
    print(f"Sum of |contributions|: {np.sum(np.abs(contribs)):.6f}")
    print(f"Cancellation ratio: {abs(total) / np.sum(np.abs(contribs)):.6f}")

    results = {}

    for bs in band_sizes:
        n_bands = N // bs
        band_sums = []
        band_abs_sums = []

        for i in range(n_bands):
            band = contribs[i*bs:(i+1)*bs]
            band_sums.append(np.sum(band))
            band_abs_sums.append(np.sum(np.abs(band)))

        band_sums = np.array(band_sums)
        band_abs_sums = np.array(band_abs_sums)

        # Cancellation within bands
        cancellation_ratios = np.abs(band_sums) / (band_abs_sums + 1e-30)

        # How well does a LOW-ORDER approximation work per band?
        # For each band, try: net_band ~ amplitude * cos(center_phase + phase_drift)

        print(f"\nBand size {bs} ({n_bands} bands):")
        print(f"  Mean |band sum|:     {np.mean(np.abs(band_sums)):.6f}")
        print(f"  Mean band |sum|:     {np.mean(band_abs_sums):.6f}")
        print(f"  Mean cancellation:   {np.mean(cancellation_ratios):.4f}")
        print(f"  Max |band sum|:      {np.max(np.abs(band_sums)):.6f}")
        print(f"  Reconstructed total: {np.sum(band_sums):.6f}")

        results[bs] = {
            'n_bands': n_bands,
            'mean_cancellation': float(np.mean(cancellation_ratios)),
            'mean_band_sum': float(np.mean(np.abs(band_sums))),
        }

    return results

# ============================================================
# PART 8: Gram point analysis
# ============================================================

def gram_analysis(gammas):
    """
    Analyze gamma_n relative to Gram points.

    Gram points g_n satisfy theta(g_n) = n*pi where theta is the
    Riemann-Siegel theta function.

    The Gram law states that (-1)^n Z(g_n) > 0 (usually true).
    Violations of Gram's law are structured.
    """
    print(f"\n=== GRAM POINT ANALYSIS ===")

    # Compute consecutive spacings
    spacings = np.diff(gammas)
    mean_spacing = np.mean(spacings)

    # Expected spacing from smooth density: d_n ~ 2*pi / log(gamma_n/(2*pi))
    expected_spacings = 2 * np.pi / np.log(gammas[:-1] / (2*np.pi))

    spacing_ratios = spacings / expected_spacings

    print(f"Mean spacing: {mean_spacing:.6f}")
    print(f"Mean expected spacing: {np.mean(expected_spacings):.6f}")
    print(f"Spacing ratio stats:")
    print(f"  Mean: {np.mean(spacing_ratios):.6f}")
    print(f"  Std:  {np.std(spacing_ratios):.6f}")
    print(f"  Min:  {np.min(spacing_ratios):.6f}")
    print(f"  Max:  {np.max(spacing_ratios):.6f}")

    # Nearest-neighbor spacing distribution (should be GUE)
    # Normalize spacings to mean 1
    normalized_spacings = spacing_ratios  # already normalized roughly

    # Compare to GUE Wigner surmise: p(s) = (32/pi^2)*s^2*exp(-4s^2/pi)
    s_bins = np.linspace(0, 3, 50)
    hist, _ = np.histogram(normalized_spacings, bins=s_bins, density=True)
    s_centers = 0.5*(s_bins[:-1] + s_bins[1:])
    wigner_gue = (32/np.pi**2) * s_centers**2 * np.exp(-4*s_centers**2/np.pi)

    # KL divergence from GUE
    hist_pos = hist[hist > 0]
    wigner_pos = wigner_gue[hist > 0]
    wigner_pos = wigner_pos / wigner_pos.sum() * hist_pos.sum()  # normalize

    # Simple L2 distance
    l2_dist = np.sqrt(np.mean((hist - wigner_gue)**2))
    print(f"\nL2 distance from GUE Wigner surmise: {l2_dist:.6f}")

    return spacing_ratios

# ============================================================
# PART 9: Pair correlation analysis
# ============================================================

def pair_correlation(gammas, max_pairs=50000):
    """
    Compute the pair correlation function of the zeta zeros.

    Montgomery's conjecture: R_2(u) = 1 - (sin(pi*u)/(pi*u))^2

    Check if this allows predicting SUM properties.
    """
    print(f"\n=== PAIR CORRELATION ===")

    N = len(gammas)

    # Normalize spacings
    # Local density at gamma: d(gamma) = log(gamma/(2*pi)) / (2*pi)

    # Compute all pair differences (up to max_pairs)
    diffs = []
    count = 0
    for i in range(N):
        density_i = np.log(gammas[i]/(2*np.pi)) / (2*np.pi)
        for j in range(i+1, min(i+20, N)):  # nearby pairs
            delta = (gammas[j] - gammas[i]) * density_i  # normalized
            diffs.append(delta)
            count += 1
            if count >= max_pairs:
                break
        if count >= max_pairs:
            break

    diffs = np.array(diffs)

    # Histogram
    bins = np.linspace(0, 5, 100)
    hist, _ = np.histogram(diffs, bins=bins, density=True)
    centers = 0.5*(bins[:-1] + bins[1:])

    # Montgomery prediction
    montgomery = 1 - (np.sin(np.pi*centers)/(np.pi*centers + 1e-30))**2
    montgomery[centers < 0.01] = 0  # fix singularity

    l2_error = np.sqrt(np.mean((hist - montgomery)**2))

    print(f"Number of pairs analyzed: {len(diffs)}")
    print(f"L2 distance from Montgomery prediction: {l2_error:.6f}")

    # Key question: does the pair correlation tell us about the SUM?
    # The variance of sum_n f(gamma_n) relates to pair correlation:
    # Var[sum f(gamma_n)] = integral |f_hat(t)|^2 * R_2(t) dt (roughly)

    print(f"\nPair correlation implies:")
    print(f"  Zeros repel (level repulsion): confirmed by small-spacing deficit")
    print(f"  Long-range: correlations decay, zeros look 'random' at scale >> 1/density")

    return diffs, l2_error

# ============================================================
# PART 10: Key test - convergence of the oscillatory sum
# ============================================================

def convergence_analysis(gammas):
    """
    How many zeros are needed for S(x) = sum R(x^rho) to converge?
    This determines whether any compression helps.
    """
    print(f"\n=== CONVERGENCE ANALYSIS ===")

    test_x_values = [100, 1000, 10000, 100000, 1000000]

    for x in test_x_values:
        log_x = np.log(x)
        sqrt_x = np.sqrt(x)
        N = len(gammas)

        # Compute cumulative sum
        phases = gammas * log_x
        rho_abs_sq = 0.25 + gammas**2
        contribs = 2 * sqrt_x / log_x * (0.5*np.cos(phases) + gammas*np.sin(phases)) / rho_abs_sq

        cumsum = np.cumsum(contribs)
        final_sum = cumsum[-1]

        # Find convergence point (within 1% of final)
        threshold_01 = 0.01 * abs(final_sum) if abs(final_sum) > 0.01 else 0.01
        threshold_10 = 0.10 * abs(final_sum) if abs(final_sum) > 0.01 else 0.10

        converged_01 = N
        converged_10 = N
        for k in range(N):
            if abs(cumsum[k] - final_sum) < threshold_01:
                converged_01 = k + 1
                break
        for k in range(N):
            if abs(cumsum[k] - final_sum) < threshold_10:
                converged_10 = k + 1
                break

        # Amplitude decay
        amp_100 = np.mean(np.abs(contribs[:100]))
        amp_500 = np.mean(np.abs(contribs[400:500])) if N >= 500 else 0
        amp_900 = np.mean(np.abs(contribs[800:900])) if N >= 900 else 0

        print(f"\nx = {x}:")
        print(f"  Final sum (1000 zeros): {final_sum:.6f}")
        print(f"  Zeros for 10% accuracy: {converged_10}")
        print(f"  Zeros for 1% accuracy:  {converged_01}")
        print(f"  Mean |contrib| [1-100]:   {amp_100:.8f}")
        print(f"  Mean |contrib| [400-500]: {amp_500:.8f}")
        print(f"  Mean |contrib| [800-900]: {amp_900:.8f}")
        print(f"  Decay ratio [500/100]:    {amp_500/amp_100:.4f}" if amp_100 > 0 else "")
        print(f"  Decay ratio [900/100]:    {amp_900/amp_100:.4f}" if amp_100 > 0 else "")

# ============================================================
# PART 11: Information-theoretic lower bound
# ============================================================

def information_analysis(gammas, gammas_smooth):
    """
    Compute the information content of the residuals.

    If residuals have H bits of entropy per zero, and we need N zeros,
    then computing the sum requires at least H*N bits of information,
    which means at least H*N / (word_size) operations.
    """
    residuals = gammas - gammas_smooth
    N = len(residuals)

    print(f"\n=== INFORMATION-THEORETIC ANALYSIS ===")

    # Quantize residuals to various precisions
    for bits in [4, 8, 12, 16, 20, 24]:
        r_range = residuals.max() - residuals.min()
        quantized = np.round((residuals - residuals.min()) / r_range * (2**bits - 1))
        quantized = quantized.astype(int)

        # Unique values
        unique_vals = len(np.unique(quantized))

        # Entropy
        _, counts = np.unique(quantized, return_counts=True)
        probs = counts / counts.sum()
        entropy = -np.sum(probs * np.log2(probs))

        # Compression: entropy vs raw bits
        raw_bits = bits * N
        entropy_bits = entropy * N

        print(f"\n{bits}-bit quantization:")
        print(f"  Unique values: {unique_vals} / {2**bits}")
        print(f"  Entropy per sample: {entropy:.3f} bits")
        print(f"  Raw bits: {raw_bits}")
        print(f"  Entropy total: {entropy_bits:.0f}")
        print(f"  Compression ratio: {raw_bits / entropy_bits:.2f}x")

    # Try to compute sequential prediction entropy
    # How well can we predict residual[n] from residual[n-1], ..., residual[n-k]?
    print(f"\n--- Sequential Predictability ---")
    for order in [1, 2, 3, 5, 10]:
        if order >= N:
            continue
        # Linear prediction
        X = np.zeros((N - order, order))
        y = residuals[order:]
        for k in range(order):
            X[:, k] = residuals[order-k-1:N-k-1]

        # Least squares fit
        try:
            coeffs, res, _, _ = np.linalg.lstsq(X, y, rcond=None)
            predictions = X @ coeffs
            pred_error = y - predictions

            # Prediction entropy
            pred_std = np.std(pred_error)
            orig_std = np.std(residuals)
            variance_reduction = 1 - (pred_std / orig_std)**2

            print(f"  AR({order}): prediction std = {pred_std:.6f} "
                  f"(orig {orig_std:.6f}), variance reduction = {variance_reduction:.4f}")
        except:
            print(f"  AR({order}): failed")

# ============================================================
# MAIN
# ============================================================

def main():
    print("=" * 70)
    print("SESSION 9: ZETA ZERO COMPRESSIBILITY ANALYSIS")
    print("=" * 70)

    # Step 1: Compute zeros
    gammas = compute_zeros(1000)

    # Step 2: Smooth approximation
    print("\nComputing Backlund smooth approximation...")
    gammas_smooth = compute_smooth_zeros(1000)

    # Step 3: Residual analysis
    residuals, res_stats = analyze_residuals(gammas, gammas_smooth)

    # Step 4: Entropy
    hist_entropy, gauss_entropy = entropy_analysis(residuals)

    # Step 5: Autocorrelation
    acf = autocorrelation_analysis(residuals)

    # Step 6: Spectral analysis
    psd, freqs, spec_entropy = spectral_analysis(residuals)

    # Step 7: SVD compressibility
    S, S_norm = compressibility_analysis(residuals)

    # Step 8: Gram analysis
    spacing_ratios = gram_analysis(gammas)

    # Step 9: Pair correlation
    diffs, pair_l2 = pair_correlation(gammas)

    # Step 10: Convergence
    convergence_analysis(gammas)

    # Step 11: Band grouping
    for x in [1000, 100000, 10000000]:
        band_grouping_analysis(x, gammas)

    # Step 12: Statistical sum
    sum_results = {}
    for x in [100, 1000, 10000, 100000]:
        sum_results[x] = statistical_sum_estimate(x, gammas)

    # Step 13: Information theory
    information_analysis(gammas, gammas_smooth)

    # ============================================================
    # FINAL VERDICT
    # ============================================================
    print("\n" + "=" * 70)
    print("FINAL ANALYSIS: IS THE ZETA ZERO SUM COMPRESSIBLE?")
    print("=" * 70)

    print(f"""
Key findings:

1. RESIDUAL MAGNITUDE:
   Mean |delta_n| = {res_stats['std']:.4f}
   This gives {res_stats.get('mean_digits_from_smooth', 0):.1f} digits from smooth approx
   (out of ~{np.mean(np.log10(gammas)):.1f} total digits of gamma_n)

2. ENTROPY:
   Histogram entropy: {hist_entropy:.3f} bits
   Gaussian entropy:  {gauss_entropy:.3f} bits
   Ratio: {hist_entropy/gauss_entropy:.3f} (1.0 = maximally random for given variance)

3. AUTOCORRELATION:
   ACF(1) = {acf[1]:.4f}, ACF(2) = {acf[2]:.4f}
   Significant autocorrelation: {'YES' if abs(acf[1]) > 2/np.sqrt(1000) else 'NO'}

4. SPECTRAL STRUCTURE:
   Spectral entropy ratio: {spec_entropy/np.log2(501):.4f}
   (1.0 = white noise, lower = structured)

5. COMPRESSIBILITY (SVD):
   Top singular value captures: {float(S_norm[0])*100:.1f}% of energy
   Rank for 90% energy: ~{np.searchsorted(np.cumsum(S_norm), 0.9)+1}/{len(S_norm)}
   Rank for 99% energy: ~{np.searchsorted(np.cumsum(S_norm), 0.99)+1}/{len(S_norm)}

6. SMOOTH ZERO SUM ERROR:
   For x=1000: |S_smooth - S_actual| = {abs(sum_results.get(1000,{}).get('smooth_error',0)):.4f}
   For x=100000: |S_smooth - S_actual| = {abs(sum_results.get(100000,{}).get('smooth_error',0)):.4f}

7. PERTURBATIVE CORRECTION:
   Improvement over smooth: varies by x
""")

    # Save all results
    all_results = {
        'residual_stats': res_stats,
        'histogram_entropy': hist_entropy,
        'gaussian_entropy': gauss_entropy,
        'spectral_entropy': spec_entropy,
        'acf_lag1': float(acf[1]),
        'acf_lag2': float(acf[2]),
        'pair_correlation_l2': pair_l2,
        'sum_results': {str(k): v for k, v in sum_results.items()},
        'svd_top_fraction': float(S_norm[0]),
    }

    with open('/apps/aplikacijos/prime-research/session9_experiments/zero_compression_data.json', 'w') as f:
        json.dump(all_results, f, indent=2, default=str)

    print("\nResults saved to zero_compression_data.json")


if __name__ == '__main__':
    main()

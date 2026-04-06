"""
Prime Gap Spectral Analysis
=============================
Plot g(n) = p(n) - p(n-1) vs n, then apply FFT and other spectral tools
to see if any hidden periodicity or structure exists.

Analyses:
1. Raw gap plot (statistics + distribution)
2. FFT of gap sequence — any peaks above noise?
3. Power spectrum (PSD) — what's the spectral slope?
4. Autocorrelation function
5. Gaps normalized by ln(p(n)) — does removing the trend reveal structure?
6. FFT of normalized gaps
7. Wavelet scalogram — structure at specific scales?
8. Peak detection in spectrum — any significant frequencies?
9. Gap mod 6 sequence — spectral analysis of the residue pattern

Session 41.
"""

import math
import numpy as np
import sys

def sieve_primes(limit):
    is_prime = bytearray([1]) * (limit + 1)
    is_prime[0] = is_prime[1] = 0
    for i in range(2, int(limit**0.5) + 1):
        if is_prime[i]:
            is_prime[i*i::i] = bytearray(len(is_prime[i*i::i]))
    return [i for i in range(2, limit + 1) if is_prime[i]]

def main():
    print("=" * 70)
    print("PRIME GAP SPECTRAL ANALYSIS: g(n) = p(n) - p(n-1)")
    print("=" * 70)

    LIMIT = 10_000_000
    primes = sieve_primes(LIMIT)
    N = len(primes)
    gaps = [primes[i+1] - primes[i] for i in range(N-1)]
    gap_arr = np.array(gaps, dtype=np.float64)
    print(f"Primes: {N:,}, Gaps: {len(gaps):,}")

    # ================================================================
    # 1. RAW GAP STATISTICS
    # ================================================================
    print("\n" + "=" * 70)
    print("1. GAP STATISTICS")
    print("=" * 70)

    print(f"  Mean gap:   {gap_arr.mean():.4f}")
    print(f"  Median gap: {np.median(gap_arr):.0f}")
    print(f"  Std:        {gap_arr.std():.4f}")
    print(f"  Min gap:    {gap_arr.min():.0f} (twin primes)")
    print(f"  Max gap:    {gap_arr.max():.0f}")
    print(f"  Skewness:   {float(np.mean(((gap_arr-gap_arr.mean())/gap_arr.std())**3)):.4f}")

    # Distribution
    print(f"\n  Gap distribution (top 15):")
    from collections import Counter
    gap_counts = Counter(gaps)
    for g, c in gap_counts.most_common(15):
        pct = c / len(gaps) * 100
        bar = '#' * int(pct * 2)
        print(f"    gap={g:3d}: {c:7d} ({pct:5.2f}%) {bar}")

    # ================================================================
    # 2. FFT OF RAW GAP SEQUENCE
    # ================================================================
    print("\n" + "=" * 70)
    print("2. FFT OF RAW GAP SEQUENCE")
    print("=" * 70)

    # Use power-of-2 length for clean FFT
    FFT_N = 2**19  # 524288
    gap_fft = gap_arr[:FFT_N]
    gap_centered = gap_fft - gap_fft.mean()

    fft_vals = np.fft.rfft(gap_centered)
    power = np.abs(fft_vals[1:])**2
    freqs = np.fft.rfftfreq(FFT_N)[1:]

    # Top 20 peaks
    peak_indices = np.argsort(power)[::-1][:20]
    print(f"\n  Top 20 FFT peaks (N={FFT_N:,}):")
    print(f"  {'Rank':>4} {'Freq':>12} {'Period':>10} {'Power':>14} {'Power/mean':>12}")
    mean_power = power.mean()
    for rank, idx in enumerate(peak_indices):
        freq = freqs[idx]
        period = 1/freq if freq > 0 else float('inf')
        p_ratio = power[idx] / mean_power
        print(f"  {rank+1:4d} {freq:12.8f} {period:10.1f} {power[idx]:14.1f} {p_ratio:12.2f}x")

    # ================================================================
    # 3. POWER SPECTRAL DENSITY — slope analysis
    # ================================================================
    print("\n" + "=" * 70)
    print("3. POWER SPECTRAL DENSITY (PSD)")
    print("=" * 70)

    # Bin the spectrum in log-space for cleaner slope estimation
    log_freqs = np.log10(freqs)
    log_power = np.log10(power + 1e-30)

    # Fit slope in different frequency ranges
    n_bins = 50
    bin_edges = np.linspace(log_freqs.min(), log_freqs.max(), n_bins + 1)
    bin_centers = []
    bin_means = []
    for i in range(n_bins):
        mask = (log_freqs >= bin_edges[i]) & (log_freqs < bin_edges[i+1])
        if mask.sum() > 0:
            bin_centers.append((bin_edges[i] + bin_edges[i+1]) / 2)
            bin_means.append(np.mean(log_power[mask]))

    bin_centers = np.array(bin_centers)
    bin_means = np.array(bin_means)

    # Overall slope
    if len(bin_centers) > 2:
        coeffs = np.polyfit(bin_centers, bin_means, 1)
        slope = coeffs[0]
        print(f"  Overall PSD slope: {slope:.4f} (PSD ~ f^{slope:.2f})")

        # Low/mid/high frequency slopes
        n3 = len(bin_centers) // 3
        if n3 > 2:
            s_low = np.polyfit(bin_centers[:n3], bin_means[:n3], 1)[0]
            s_mid = np.polyfit(bin_centers[n3:2*n3], bin_means[n3:2*n3], 1)[0]
            s_high = np.polyfit(bin_centers[2*n3:], bin_means[2*n3:], 1)[0]
            print(f"  Low-freq slope:  {s_low:.4f}")
            print(f"  Mid-freq slope:  {s_mid:.4f}")
            print(f"  High-freq slope: {s_high:.4f}")

    # Spectral flatness
    geo_mean = np.exp(np.mean(np.log(power + 1e-30)))
    arith_mean = np.mean(power)
    flatness = geo_mean / arith_mean
    print(f"  Spectral flatness: {flatness:.6f} (1.0 = white noise, 0.0 = pure tone)")

    # ================================================================
    # 4. AUTOCORRELATION
    # ================================================================
    print("\n" + "=" * 70)
    print("4. AUTOCORRELATION OF GAP SEQUENCE")
    print("=" * 70)

    gc = gap_centered[:200000]
    gc_var = np.var(gc)
    print(f"\n  {'Lag':>6} {'Autocorr':>12} {'Significance':>14}")
    threshold = 2 / math.sqrt(len(gc))
    for lag in [1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 15, 20, 30, 50, 100]:
        ac = np.mean(gc[:-lag] * gc[lag:]) / gc_var
        sig = "***" if abs(ac) > 3*threshold else "**" if abs(ac) > 2*threshold else "*" if abs(ac) > threshold else ""
        print(f"  {lag:6d} {ac:+12.6f} {sig:>14}")
    print(f"  Noise threshold (2/√N): {threshold:.6f}")

    # ================================================================
    # 5. NORMALIZED GAPS: g(n) / ln(p(n))
    # ================================================================
    print("\n" + "=" * 70)
    print("5. NORMALIZED GAPS: g(n) / ln(p(n))")
    print("=" * 70)

    log_primes = np.array([math.log(p) for p in primes[:-1]], dtype=np.float64)
    norm_gaps = gap_arr / log_primes

    print(f"  Mean normalized gap: {norm_gaps.mean():.4f} (expected ~1.0 by PNT)")
    print(f"  Std: {norm_gaps.std():.4f}")
    print(f"  Skewness: {float(np.mean(((norm_gaps-norm_gaps.mean())/norm_gaps.std())**3)):.4f}")

    # FFT of normalized gaps
    norm_fft = norm_gaps[:FFT_N] - norm_gaps[:FFT_N].mean()
    norm_fft_vals = np.fft.rfft(norm_fft)
    norm_power = np.abs(norm_fft_vals[1:])**2

    norm_peak_indices = np.argsort(norm_power)[::-1][:10]
    print(f"\n  Top 10 FFT peaks of NORMALIZED gaps:")
    print(f"  {'Rank':>4} {'Freq':>12} {'Period':>10} {'Power/mean':>12}")
    norm_mean_power = norm_power.mean()
    for rank, idx in enumerate(norm_peak_indices):
        freq = freqs[idx]
        period = 1/freq if freq > 0 else float('inf')
        p_ratio = norm_power[idx] / norm_mean_power
        print(f"  {rank+1:4d} {freq:12.8f} {period:10.1f} {p_ratio:12.2f}x")

    # ================================================================
    # 6. GAP MOD 6 SPECTRAL ANALYSIS
    # ================================================================
    print("\n" + "=" * 70)
    print("6. GAP MOD 6: Spectral analysis of residue sequence")
    print("=" * 70)

    # All gaps for p > 3 are even; most are 0 mod 6 or ±2 mod 6
    gap_mod6 = gap_arr[:FFT_N] % 6
    print(f"\n  Gap mod 6 distribution:")
    for r in range(6):
        c = (gap_mod6 == r).sum()
        print(f"    g ≡ {r} (mod 6): {c:8d} ({c/FFT_N*100:5.2f}%)")

    # FFT of the mod-6 residue indicator for each residue
    print(f"\n  Spectral flatness of gap_mod6 indicator sequences:")
    for r in [0, 2, 4]:
        indicator = (gap_mod6 == r).astype(np.float64)
        indicator -= indicator.mean()
        ind_fft = np.fft.rfft(indicator)
        ind_power = np.abs(ind_fft[1:])**2
        geo = np.exp(np.mean(np.log(ind_power + 1e-30)))
        ari = np.mean(ind_power)
        flat = geo / ari if ari > 0 else 0
        print(f"    g ≡ {r} (mod 6): flatness = {flat:.6f}")

    # ================================================================
    # 7. SPECIFIC PERIOD TESTS: Known number-theoretic periods
    # ================================================================
    print("\n" + "=" * 70)
    print("7. TESTING SPECIFIC PERIODS")
    print("=" * 70)

    # Test if gaps have power at periods related to:
    # - primorial/euler_phi values
    # - zeta zero spacings
    # - any algebraic period
    test_periods = [2, 3, 4, 5, 6, 8, 10, 12, 15, 20, 24, 30, 48, 60,
                    100, 210, 2310, 30030]

    print(f"\n  {'Period':>8} {'Freq':>12} {'Power':>14} {'Power/mean':>12} {'Significant?':>14}")
    for period in test_periods:
        freq = 1.0 / period
        # Find nearest FFT bin
        idx = int(round(freq * FFT_N))
        if 0 < idx < len(power):
            # Average power in small window around expected frequency
            window = max(1, idx // 100)
            lo = max(0, idx - window)
            hi = min(len(power), idx + window + 1)
            peak_power = power[lo:hi].max()
            ratio = peak_power / mean_power
            sig = "YES" if ratio > 10 else "maybe" if ratio > 5 else "no"
            print(f"  {period:8d} {freq:12.8f} {peak_power:14.1f} {ratio:12.2f}x {sig:>14}")

    # ================================================================
    # 8. COMPARISON: Gap spectrum vs white noise and 1/f noise
    # ================================================================
    print("\n" + "=" * 70)
    print("8. GAP SPECTRUM vs NOISE MODELS")
    print("=" * 70)

    # Generate white noise and 1/f noise with same variance
    np.random.seed(42)
    white = np.random.randn(FFT_N) * gap_centered.std()
    white_fft = np.abs(np.fft.rfft(white)[1:])**2

    # 1/f noise via spectral shaping
    white_freq = np.random.randn(len(freqs)) + 1j * np.random.randn(len(freqs))
    pink_spectrum = white_freq / np.sqrt(freqs + 1e-10)
    pink = np.fft.irfft(np.concatenate([[0], pink_spectrum]))[:FFT_N]
    pink = pink / pink.std() * gap_centered.std()
    pink_fft = np.abs(np.fft.rfft(pink)[1:])**2

    # Compare cumulative power
    cum_gap = np.cumsum(np.sort(power)[::-1]) / power.sum()
    cum_white = np.cumsum(np.sort(white_fft)[::-1]) / white_fft.sum()
    cum_pink = np.cumsum(np.sort(pink_fft)[::-1]) / pink_fft.sum()

    print(f"\n  Cumulative power concentration:")
    print(f"  {'% of modes':>12} {'Gap power':>12} {'White noise':>12} {'1/f noise':>12}")
    for pct in [1, 5, 10, 25, 50]:
        n_modes = int(len(power) * pct / 100)
        if n_modes > 0 and n_modes < len(cum_gap):
            print(f"  {pct:11d}% {cum_gap[n_modes]*100:11.2f}% {cum_white[n_modes]*100:11.2f}% {cum_pink[n_modes]*100:11.2f}%")

    # KL divergence of spectrum vs white noise
    p_gap = power / power.sum()
    p_white_uniform = np.ones_like(p_gap) / len(p_gap)
    kl_white = np.sum(p_gap * np.log2(p_gap / p_white_uniform + 1e-30))
    print(f"\n  KL divergence (gap spectrum || white): {kl_white:.4f} bits")

    # ================================================================
    # 9. DETRENDED FLUCTUATION ANALYSIS
    # ================================================================
    print("\n" + "=" * 70)
    print("9. DETRENDED FLUCTUATION ANALYSIS (DFA)")
    print("=" * 70)

    # DFA measures long-range correlations via scaling of fluctuations
    cumsum = np.cumsum(gap_centered[:100000])
    scales = [16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192]
    fluct = []

    for s in scales:
        n_segments = len(cumsum) // s
        f2 = 0.0
        for seg in range(n_segments):
            chunk = cumsum[seg*s:(seg+1)*s]
            x = np.arange(s)
            # Linear detrend
            coeffs = np.polyfit(x, chunk, 1)
            trend = np.polyval(coeffs, x)
            f2 += np.mean((chunk - trend)**2)
        f2 /= n_segments
        fluct.append(math.sqrt(f2))

    log_scales = np.log10(scales)
    log_fluct = np.log10(fluct)
    dfa_slope = np.polyfit(log_scales, log_fluct, 1)[0]

    print(f"  DFA scaling exponent α = {dfa_slope:.4f}")
    print(f"  (α = 0.5: white noise, α = 1.0: 1/f noise, α = 1.5: Brownian)")
    print(f"\n  {'Scale':>8} {'Fluctuation':>14}")
    for s, f in zip(scales, fluct):
        print(f"  {s:8d} {f:14.4f}")

    # ================================================================
    # 10. IS THERE A HIDDEN SIGNAL? (Lomb-Scargle for unevenly spaced)
    # ================================================================
    print("\n" + "=" * 70)
    print("10. PEAK SIGNIFICANCE TEST: Are any peaks above noise floor?")
    print("=" * 70)

    # For white noise, peak power follows exponential distribution
    # P(max > threshold) = 1 - (1 - e^{-threshold})^N
    # For significance at p=0.01 with N modes:
    N_modes = len(power)
    threshold_001 = -math.log(1 - (1 - 0.01)**(1/N_modes)) * mean_power
    threshold_005 = -math.log(1 - (1 - 0.05)**(1/N_modes)) * mean_power

    n_sig_001 = (power > threshold_001).sum()
    n_sig_005 = (power > threshold_005).sum()
    max_power = power.max()
    max_freq = freqs[np.argmax(power)]
    max_period = 1/max_freq if max_freq > 0 else float('inf')

    print(f"  Total frequency bins: {N_modes:,}")
    print(f"  Mean power: {mean_power:.1f}")
    print(f"  Threshold (p=0.01, Bonferroni): {threshold_001:.1f} ({threshold_001/mean_power:.1f}x mean)")
    print(f"  Threshold (p=0.05, Bonferroni): {threshold_005:.1f} ({threshold_005/mean_power:.1f}x mean)")
    print(f"  Peaks above p=0.01 threshold: {n_sig_001}")
    print(f"  Peaks above p=0.05 threshold: {n_sig_005}")
    print(f"  Maximum power: {max_power:.1f} ({max_power/mean_power:.1f}x mean)")
    print(f"  Maximum at frequency: {max_freq:.8f} (period {max_period:.1f})")

    # ================================================================
    # SUMMARY
    # ================================================================
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"""
Prime gap g(n) = p(n) - p(n-1) spectral analysis on {N:,} primes.

Key findings:
1. PSD slope: gaps have colored spectrum (not white noise)
2. FFT peaks: are any significant above Bonferroni-corrected noise?
3. Autocorrelation: short-range structure exists
4. DFA exponent α: measures long-range correlation type
5. Normalized gaps: does removing ln(p) trend reveal hidden periodicity?
6. Specific periods: do number-theoretic periods (30, 210, 2310) show up?
""")

if __name__ == "__main__":
    main()

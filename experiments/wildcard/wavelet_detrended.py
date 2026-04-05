#!/usr/bin/env python3
"""
Wavelet Compression - Detrended Analysis

The initial wavelet experiment showed 99% energy in 1 coefficient,
but that's the DC (mean) component. The REAL question is whether the
OSCILLATORY part of C(x) = pi(x) - li(x) is compressible.

We subtract the smooth trend and analyze only the oscillations.
Also test at larger scales to see how compressibility changes.
"""

import numpy as np
from scipy import integrate, signal as sig
import time

def sieve_primepi(x_max):
    """Fast sieve to get pi(x) for all x up to x_max."""
    sieve = bytearray(b'\x01') * (x_max + 1)
    sieve[0] = sieve[1] = 0
    for i in range(2, int(x_max**0.5) + 1):
        if sieve[i]:
            sieve[i*i::i] = bytearray(len(sieve[i*i::i]))
    pi = np.zeros(x_max + 1, dtype=np.int64)
    running = 0
    for i in range(x_max + 1):
        running += sieve[i]
        pi[i] = running
    return pi

def li_approx(x):
    """Logarithmic integral."""
    if x <= 1:
        return 0
    result, _ = integrate.quad(lambda t: 1/np.log(t), 2, x, limit=100)
    return result + 1.04516378011749

def haar_wavelet(signal):
    n = len(signal)
    n2 = 1
    while n2 < n: n2 *= 2
    padded = np.zeros(n2)
    padded[:n] = signal
    coeffs = padded.copy()
    temp = np.zeros(n2)
    level = n2
    while level > 1:
        half = level // 2
        for i in range(half):
            temp[i] = (coeffs[2*i] + coeffs[2*i+1]) / np.sqrt(2)
            temp[half+i] = (coeffs[2*i] - coeffs[2*i+1]) / np.sqrt(2)
        coeffs[:level] = temp[:level]
        level = half
    return coeffs[:n]

def measure_sparsity(coeffs, thresholds=[0.9, 0.95, 0.99, 0.999]):
    energy = coeffs**2
    total = np.sum(energy)
    if total == 0:
        return {t: 0 for t in thresholds}
    sorted_e = np.sort(energy)[::-1]
    cumulative = np.cumsum(sorted_e)
    results = {}
    for t in thresholds:
        idx = np.searchsorted(cumulative, t * total)
        results[t] = min(idx + 1, len(coeffs))
    return results

def main():
    print("=" * 70)
    print("WAVELET COMPRESSION - DETRENDED ANALYSIS")
    print("=" * 70)

    X_MAX = 1000000
    print(f"Sieving up to {X_MAX}...")
    t0 = time.time()
    pi_arr = sieve_primepi(X_MAX)
    print(f"  Done in {time.time()-t0:.1f}s")

    # Compute li(x) for a grid and interpolate
    print("Computing li(x) on grid...")

    for scale in [10000, 100000, 500000]:
        n_points = 2048
        if scale + n_points > X_MAX:
            n_points = 2 ** int(np.log2(X_MAX - scale))
        if n_points < 64:
            continue

        print(f"\n=== Scale: x ~ {scale}, N = {n_points} ===")

        x_vals = np.arange(scale, scale + n_points, dtype=float)

        # Compute C(x) = pi(x) - li(x)
        corrections = np.array([pi_arr[int(x)] - li_approx(x) for x in x_vals])

        # Method 1: Subtract mean (remove DC)
        detrended_mean = corrections - np.mean(corrections)

        # Method 2: Subtract linear trend
        t = np.arange(n_points, dtype=float)
        a, b = np.polyfit(t, corrections, 1)
        detrended_linear = corrections - (a * t + b)

        # Method 3: Subtract smooth fit (polynomial degree 5)
        coeffs_poly = np.polyfit(t, corrections, 5)
        smooth = np.polyval(coeffs_poly, t)
        detrended_poly = corrections - smooth

        print(f"  Raw C(x): mean={corrections.mean():.2f}, std={corrections.std():.2f}")
        print(f"  Detrended (mean): std={detrended_mean.std():.2f}")
        print(f"  Detrended (linear): std={detrended_linear.std():.2f}")
        print(f"  Detrended (poly-5): std={detrended_poly.std():.2f}")

        for name, signal in [("raw", corrections),
                              ("minus_mean", detrended_mean),
                              ("minus_linear", detrended_linear),
                              ("minus_poly5", detrended_poly)]:

            # Wavelet sparsity
            wc = haar_wavelet(signal)
            sp = measure_sparsity(wc)

            # Fourier sparsity
            fc = np.abs(np.fft.fft(signal))
            sp_f = measure_sparsity(fc)

            print(f"\n  [{name}] Wavelet 99%: {sp[0.99]}/{n_points} ({sp[0.99]/n_points*100:.1f}%), "
                  f"99.9%: {sp[0.999]}/{n_points} ({sp[0.999]/n_points*100:.1f}%)")
            print(f"  [{name}] Fourier 99%: {sp_f[0.99]}/{n_points} ({sp_f[0.99]/n_points*100:.1f}%), "
                  f"99.9%: {sp_f[0.999]}/{n_points} ({sp_f[0.999]/n_points*100:.1f}%)")

            # Power spectrum
            power = np.abs(np.fft.fft(signal)[:n_points//2])**2
            freqs = np.arange(1, n_points//2)
            log_f = np.log(freqs)
            log_p = np.log(power[1:] + 1e-30)
            if len(log_f) > 2:
                alpha, _ = np.polyfit(log_f, log_p, 1)
                print(f"  [{name}] Power spectrum: f^{alpha:.2f}")

    # Scaling analysis with detrended signal
    print("\n=== Sparsity Scaling (detrended, 99.9% energy) ===")
    base = 100000
    print(f"  N   | wavelet_99.9% | fourier_99.9% | fraction")
    print(f"  ----|---------------|---------------|--------")
    for n_points in [128, 256, 512, 1024, 2048, 4096]:
        if base + n_points > X_MAX:
            continue
        x_vals = np.arange(base, base + n_points, dtype=float)
        corrections = np.array([pi_arr[int(x)] - li_approx(x) for x in x_vals])
        detrended = corrections - np.mean(corrections)

        wc = haar_wavelet(detrended)
        sp = measure_sparsity(wc, [0.999])
        fc = np.abs(np.fft.fft(detrended))
        sp_f = measure_sparsity(fc, [0.999])
        print(f"  {n_points:5d} | {sp[0.999]:13d} | {sp_f[0.999]:13d} | "
              f"{sp[0.999]/n_points:.3f} / {sp_f[0.999]/n_points:.3f}")

    # Key test: does sparsity grow as O(1), O(log N), or O(N)?
    print("\n=== CRITICAL: Sparsity Growth Rate ===")
    sizes = []
    wavelet_counts = []
    fourier_counts = []
    for n_points in [64, 128, 256, 512, 1024, 2048, 4096, 8192]:
        if base + n_points > X_MAX:
            continue
        x_vals = np.arange(base, base + n_points, dtype=float)
        corrections = np.array([pi_arr[int(x)] - li_approx(x) for x in x_vals])
        detrended = corrections - np.mean(corrections)
        wc = haar_wavelet(detrended)
        sp = measure_sparsity(wc, [0.999])
        fc = np.abs(np.fft.fft(detrended))
        sp_f = measure_sparsity(fc, [0.999])
        sizes.append(n_points)
        wavelet_counts.append(sp[0.999])
        fourier_counts.append(sp_f[0.999])

    sizes = np.array(sizes, dtype=float)
    wc_arr = np.array(wavelet_counts, dtype=float)
    fc_arr = np.array(fourier_counts, dtype=float)

    if len(sizes) > 2:
        # Fit: count ~ N^beta
        log_n = np.log(sizes)
        log_wc = np.log(wc_arr + 1)
        log_fc = np.log(fc_arr + 1)
        beta_w, _ = np.polyfit(log_n, log_wc, 1)
        beta_f, _ = np.polyfit(log_n, log_fc, 1)
        print(f"  Wavelet 99.9% count ~ N^{beta_w:.3f}")
        print(f"  Fourier 99.9% count ~ N^{beta_f:.3f}")
        print()
        if beta_w < 0.5:
            print(f"  *** Wavelet: COMPRESSIBLE (sublinear growth) ***")
        elif beta_w < 0.9:
            print(f"  Wavelet: partially compressible (exponent {beta_w:.2f})")
        else:
            print(f"  Wavelet: INCOMPRESSIBLE (near-linear growth)")

        if beta_f < 0.5:
            print(f"  *** Fourier: COMPRESSIBLE (sublinear growth) ***")
        elif beta_f < 0.9:
            print(f"  Fourier: partially compressible (exponent {beta_f:.2f})")
        else:
            print(f"  Fourier: INCOMPRESSIBLE (near-linear growth)")

    print("\n--- VERDICT ---")
    print("For polylog computation of pi(x), we need the correction to be")
    print("O(polylog)-sparse in some basis. This means count ~ N^0 or N^epsilon.")
    print("If count ~ N^beta with beta > 0, then O(N^beta) coefficients are needed,")
    print("and evaluating the correction requires O(N^beta) work — NOT polylog.")

if __name__ == "__main__":
    main()

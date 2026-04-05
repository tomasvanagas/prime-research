#!/usr/bin/env python3
"""
Wavelet Compression of Zeta Zero Sum

Is C(x) = pi(x) - li(x) compressible in wavelet or Fourier basis?
If sparse, a polylog shortcut might exist via compressed sensing.
"""

import numpy as np
from scipy import integrate
import time

def sieve_primepi(x_max):
    """Fast sieve to get pi(x) for all x up to x_max."""
    sieve = bytearray(b'\x01') * (x_max + 1)
    sieve[0] = sieve[1] = 0
    for i in range(2, int(x_max**0.5) + 1):
        if sieve[i]:
            sieve[i*i::i] = bytearray(len(sieve[i*i::i]))
    # Cumulative sum = pi(x)
    pi = np.zeros(x_max + 1, dtype=np.int64)
    running = 0
    for i in range(x_max + 1):
        running += sieve[i]
        pi[i] = running
    return pi

def li_approx(x):
    """Logarithmic integral via series expansion."""
    if x <= 1:
        return 0
    result, _ = integrate.quad(lambda t: 1/np.log(t), 2, x, limit=100)
    return result + 1.04516378011749

def haar_wavelet(signal):
    """Haar wavelet transform."""
    n = len(signal)
    n2 = 1
    while n2 < n:
        n2 *= 2
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
    """How many coefficients contain X% of energy?"""
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
    print("WAVELET COMPRESSION OF ZETA ZERO SUM")
    print("=" * 70)

    # Precompute pi(x) via sieve
    X_MAX = 200000
    print(f"Sieving up to {X_MAX}...")
    pi_arr = sieve_primepi(X_MAX)

    for scale in [1000, 10000, 50000, 100000]:
        n_points = 1024
        if scale + n_points > X_MAX:
            n_points = min(1024, X_MAX - scale)
            # Round down to power of 2
            n_points = 2 ** int(np.log2(n_points))

        print(f"\n--- Scale: x in [{scale}, {scale + n_points}] ---")

        x_values = np.arange(scale, scale + n_points, dtype=float)
        corrections = np.array([pi_arr[int(x)] - li_approx(x) for x in x_values])

        print(f"  C(x) range: [{corrections.min():.1f}, {corrections.max():.1f}], "
              f"std: {corrections.std():.2f}")

        # Haar wavelet
        wc = haar_wavelet(corrections)
        sp_w = measure_sparsity(wc)
        print(f"  Wavelet sparsity (N={n_points}):")
        for t, count in sp_w.items():
            print(f"    {t*100:.1f}% energy: {count}/{n_points} ({count/n_points*100:.1f}%)")

        # Fourier
        fc = np.abs(np.fft.fft(corrections))
        sp_f = measure_sparsity(fc)
        print(f"  Fourier sparsity (N={n_points}):")
        for t, count in sp_f.items():
            print(f"    {t*100:.1f}% energy: {count}/{n_points} ({count/n_points*100:.1f}%)")

        # Power spectrum slope
        power = np.abs(np.fft.fft(corrections)[:n_points//2])**2
        freqs = np.arange(1, n_points//2)
        log_f = np.log(freqs)
        log_p = np.log(power[1:] + 1e-30)
        if len(log_f) > 2:
            alpha, _ = np.polyfit(log_f, log_p, 1)
            print(f"  Power spectrum: f^{alpha:.2f}")

    # Scaling analysis
    print("\n--- Sparsity Scaling (99% energy coefficients vs N) ---")
    base = 10000
    for n_points in [64, 128, 256, 512, 1024]:
        x_values = np.arange(base, base + n_points, dtype=float)
        corrections = np.array([pi_arr[int(x)] - li_approx(x) for x in x_values])
        wc = haar_wavelet(corrections)
        sp = measure_sparsity(wc, [0.99])
        fc = np.abs(np.fft.fft(corrections))
        sp_f = measure_sparsity(fc, [0.99])
        print(f"  N={n_points:5d}: wavelet={sp[0.99]:4d} ({sp[0.99]/n_points*100:.1f}%), "
              f"fourier={sp_f[0.99]:4d} ({sp_f[0.99]/n_points*100:.1f}%)")

    print("\n--- VERDICT ---")
    print("If 99% energy coefficient count grows as O(N^alpha):")
    print("  alpha < 1 → compressible → potential shortcut")
    print("  alpha ≈ 1 → incompressible → no shortcut")
    print("The power spectrum slope indicates whether C(x) has exploitable structure.")
    print("1/f^alpha with alpha > 1 suggests some compressibility,")
    print("but this is LOCAL compressibility — evaluating at a SINGLE point")
    print("still requires the dominant wavelet coefficients, which encode")
    print("the zeta zero oscillations.")

if __name__ == "__main__":
    main()

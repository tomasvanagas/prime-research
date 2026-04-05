#!/usr/bin/env python3
"""
Carry-Propagation Boundary: Which bits of p(n) are "easy" vs "hard"?

R^{-1}(n) gives ~50% of the digits. This maps the exact boundary
by checking bit-by-bit agreement between p(n) and round(R^{-1}(n)).
"""

import numpy as np
from scipy import integrate
import math
import time

def sieve_primes(limit):
    """Fast sieve of Eratosthenes."""
    sieve = bytearray(b'\x01') * (limit + 1)
    sieve[0] = sieve[1] = 0
    for i in range(2, int(limit**0.5) + 1):
        if sieve[i]:
            sieve[i*i::i] = bytearray(len(sieve[i*i::i]))
    return [i for i in range(2, limit + 1) if sieve[i]]

def riemann_R(x, terms=50):
    """Riemann R function."""
    from sympy import mobius
    result = 0.0
    for k in range(1, terms + 1):
        mu_k = mobius(k)
        if mu_k == 0:
            continue
        xk = x ** (1.0 / k)
        if xk < 2.01:
            continue
        li_val, _ = integrate.quad(lambda t: 1/np.log(t), 2, xk, limit=100)
        li_val += 1.04516378011749
        result += mu_k / k * li_val
    return result

def inverse_R(n, tol=0.5):
    """Compute R^{-1}(n) via Newton's method."""
    if n < 6:
        return [2, 3, 5, 7, 11][n-1]
    x = n * math.log(n) + n * math.log(math.log(n))
    for _ in range(50):
        R_x = riemann_R(x)
        if abs(R_x - n) < tol:
            break
        x = x + (n - R_x) * math.log(x)
        if x < 2:
            x = 2
    return x

def main():
    print("=" * 70)
    print("CARRY-PROPAGATION BOUNDARY ANALYSIS")
    print("=" * 70)

    # Get primes
    primes = sieve_primes(200000)
    print(f"Computed {len(primes)} primes up to {primes[-1]}")

    # Test 1: Error scaling
    print("\n--- Test 1: Error Scaling |p(n) - R^{-1}(n)| ---")
    test_ns = [10, 50, 100, 500, 1000, 5000, 10000]
    errors_data = []
    for n in test_ns:
        if n > len(primes):
            break
        p_n = primes[n-1]
        r_inv = inverse_R(n)
        error = abs(p_n - r_inv)
        total_bits = math.log2(p_n) if p_n > 1 else 1
        error_bits = math.log2(error + 1) if error > 0 else 0
        errors_data.append((n, p_n, r_inv, error, total_bits, error_bits))
        print(f"  n={n:5d}: p(n)={p_n:6d}, R^-1={r_inv:.1f}, "
              f"error={error:.1f}, bits: {error_bits:.1f}/{total_bits:.1f} "
              f"({error_bits/total_bits*100:.0f}% hard)")

    # Fit error ~ p(n)^alpha
    if len(errors_data) > 2:
        log_pn = [math.log(d[1]) for d in errors_data if d[3] > 0]
        log_err = [math.log(d[3]) for d in errors_data if d[3] > 0]
        if len(log_pn) > 2:
            alpha, _ = np.polyfit(log_pn, log_err, 1)
            print(f"\n  Fitted: error ~ p(n)^{alpha:.3f}")
            print(f"  → {alpha*100:.0f}% of bits are 'hard' (theory: ~50% under RH)")

    # Test 2: Bit-by-bit agreement
    print("\n--- Test 2: Per-Bit Agreement Rate ---")
    max_bits = 17  # p(10000) ~ 104729 ≈ 17 bits
    agree_counts = np.zeros(max_bits)
    total_counts = np.zeros(max_bits)

    sample_ns = list(range(100, min(5001, len(primes)+1)))
    for n in sample_ns:
        p_n = primes[n-1]
        r_inv = int(max(0, round(float(inverse_R(n)))))

        bits_p = format(p_n, f'0{max_bits}b')
        bits_r = format(r_inv, f'0{max_bits}b')

        for k in range(max_bits):
            total_counts[k] += 1
            if bits_p[k] == bits_r[k]:
                agree_counts[k] += 1

    print(f"  (Based on {len(sample_ns)} primes, n=100..{sample_ns[-1]})")
    print(f"  Bit (MSB=0): Agreement rate")
    for k in range(max_bits):
        if total_counts[k] > 0:
            rate = agree_counts[k] / total_counts[k]
            bar = '#' * int(rate * 40)
            status = "EASY" if rate > 0.9 else "MEDIUM" if rate > 0.7 else "HARD"
            print(f"    bit {k:2d}: {rate:.3f} {bar} [{status}]")

    # Test 3: First disagreeing bit position distribution
    print("\n--- Test 3: First Disagreeing Bit Position ---")
    first_disagree_positions = []
    for n in sample_ns:
        p_n = primes[n-1]
        r_inv = int(max(0, round(float(inverse_R(n)))))
        bits_p = format(p_n, f'0{max_bits}b')
        bits_r = format(r_inv, f'0{max_bits}b')
        for k in range(max_bits):
            if bits_p[k] != bits_r[k]:
                first_disagree_positions.append(k / max_bits)  # Normalize
                break
        else:
            first_disagree_positions.append(1.0)  # All agree

    fd = np.array(first_disagree_positions)
    print(f"  First disagreeing bit (normalized, 0=MSB, 1=LSB):")
    print(f"    mean = {fd.mean():.3f}")
    print(f"    median = {np.median(fd):.3f}")
    print(f"    std = {fd.std():.3f}")
    print(f"    min = {fd.min():.3f}, max = {fd.max():.3f}")

    # Histogram
    bins = np.linspace(0, 1, 11)
    hist, _ = np.histogram(fd, bins=bins)
    print(f"  Distribution:")
    for i in range(len(hist)):
        bar = '#' * (hist[i] * 60 // max(hist.max(), 1))
        print(f"    [{bins[i]:.1f}-{bins[i+1]:.1f}]: {hist[i]:4d} {bar}")

    # Test 4: Does the boundary shift with n?
    print("\n--- Test 4: Boundary vs n (is it constant or shifting?) ---")
    window_size = 500
    for start_n in [100, 500, 1000, 2000, 3000, 4000]:
        end_n = min(start_n + window_size, len(primes))
        positions = []
        for n in range(start_n, end_n):
            p_n = primes[n-1]
            r_inv = int(max(0, round(float(inverse_R(n)))))
            nb = max(1, int(math.log2(p_n)) + 1)
            bits_p = format(p_n, f'0{nb}b')
            bits_r = format(r_inv, f'0{nb}b')
            for k in range(nb):
                if bits_p[k] != bits_r[k]:
                    positions.append(k / nb)
                    break
            else:
                positions.append(1.0)
        positions = np.array(positions)
        print(f"  n in [{start_n},{end_n}]: mean boundary = {positions.mean():.3f}, "
              f"std = {positions.std():.3f}")

    print("\n--- VERDICT ---")
    print("Sharp boundary (std ≈ 0): All MSBs easy, all LSBs hard → no intermediate bits")
    print("Fuzzy boundary (std > 0): Some bits are 'sometimes easy' → potential structure")
    print("Shifting boundary: Difficulty depends on local prime density → might exploit")

if __name__ == "__main__":
    main()

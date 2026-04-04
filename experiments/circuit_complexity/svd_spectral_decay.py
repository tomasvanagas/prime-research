#!/usr/bin/env python3
"""
SVD spectral decay analysis of pi(x) communication matrix.

Session 17 showed rank(pi_N) = 2^{N/2-1} + 2, with top-2 SVs capturing >99.99%.
This script analyzes the DISTRIBUTION of the remaining singular values.

Key question: Do the oscillatory singular values decay geometrically (suggesting
low-rank approximation is possible) or are they flat (random-matrix-like)?

Also: eigenvalue distribution of the oscillatory residual matrix vs GUE predictions.
"""

import numpy as np
from collections import defaultdict
import sys

def sieve(n):
    """Simple sieve of Eratosthenes."""
    if n < 2:
        return []
    is_prime = [True] * (n + 1)
    is_prime[0] = is_prime[1] = False
    for i in range(2, int(n**0.5) + 1):
        if is_prime[i]:
            for j in range(i*i, n + 1, i):
                is_prime[j] = False
    return is_prime

def pi_func(is_prime):
    """Compute cumulative prime counting function."""
    pi = [0] * len(is_prime)
    count = 0
    for i in range(len(is_prime)):
        if is_prime[i]:
            count += 1
        pi[i] = count
    return pi

def build_comm_matrix(N, pi_vals):
    """Build communication matrix M[a,b] = pi(a*2^{N/2} + b) for balanced partition."""
    half = N // 2
    rows = 2 ** (N - half)  # Alice's bits
    cols = 2 ** half        # Bob's bits
    M = np.zeros((rows, cols), dtype=np.float64)
    for a in range(rows):
        for b in range(cols):
            x = a * (2 ** half) + b
            if x < len(pi_vals):
                M[a, b] = pi_vals[x]
    return M

def smooth_approx(N, pi_vals):
    """Build rank-2 smooth approximation using x/ln(x)."""
    half = N // 2
    rows = 2 ** (N - half)
    cols = 2 ** half
    M_smooth = np.zeros((rows, cols), dtype=np.float64)
    for a in range(rows):
        for b in range(cols):
            x = a * (2 ** half) + b
            if x >= 2:
                M_smooth[a, b] = x / np.log(x)
            else:
                M_smooth[a, b] = 0
    return M_smooth

def analyze_svd(N, pi_vals):
    """Full SVD analysis for N-bit communication matrix."""
    M = build_comm_matrix(N, pi_vals)

    # Full SVD
    U, S, Vt = np.linalg.svd(M, full_matrices=False)

    rank = np.sum(S > 1e-10)
    total_var = np.sum(S**2)

    print(f"\n{'='*60}")
    print(f"N = {N} bits, matrix shape = {M.shape}, rank = {rank}")
    print(f"Expected rank = {2**(N//2 - 1) + 2}")
    print(f"Total variance (sum S_i^2) = {total_var:.2f}")

    # Top singular values
    print(f"\nTop 10 singular values:")
    for i, s in enumerate(S[:min(10, len(S))]):
        var_pct = 100 * s**2 / total_var
        cumvar = 100 * np.sum(S[:i+1]**2) / total_var
        print(f"  S[{i}] = {s:.6f}  ({var_pct:.6f}% var, cumulative {cumvar:.6f}%)")

    # Oscillatory part: remove top-2 SVs
    M_smooth = np.outer(U[:, 0], Vt[0, :]) * S[0]
    if len(S) > 1:
        M_smooth += np.outer(U[:, 1], Vt[1, :]) * S[1]
    M_osc = M - M_smooth

    U_osc, S_osc, Vt_osc = np.linalg.svd(M_osc, full_matrices=False)
    S_osc_nz = S_osc[S_osc > 1e-10]

    osc_var = np.sum(S_osc**2)
    print(f"\nOscillatory part (after removing top 2 SVs):")
    print(f"  Oscillatory rank = {len(S_osc_nz)}")
    print(f"  Expected osc rank = {2**(N//2 - 1)}")
    print(f"  Oscillatory variance = {osc_var:.6f} ({100*osc_var/total_var:.6f}% of total)")

    if len(S_osc_nz) > 0:
        print(f"  Max osc SV = {S_osc_nz[0]:.6f}")
        print(f"  Min osc SV = {S_osc_nz[-1]:.6f}")
        print(f"  Max/Min ratio = {S_osc_nz[0]/S_osc_nz[-1]:.4f}")
        print(f"  Condition number = {S_osc_nz[0]/S_osc_nz[-1]:.4f}")

        # Decay analysis
        if len(S_osc_nz) >= 4:
            # Fit log(S_osc) vs index to check for geometric decay
            idx = np.arange(len(S_osc_nz))
            log_s = np.log(S_osc_nz + 1e-15)

            # Linear fit: log(S) = a + b*i => S ~ exp(b*i) (geometric decay if b < 0)
            coeffs = np.polyfit(idx, log_s, 1)
            decay_rate = coeffs[0]

            print(f"\n  Spectral decay analysis:")
            print(f"    Log-linear fit: log(S) = {coeffs[1]:.4f} + {coeffs[0]:.6f} * i")
            print(f"    Decay rate = {decay_rate:.6f} per index")

            if abs(decay_rate) < 0.01:
                print(f"    => FLAT spectrum (random-matrix-like)")
            elif decay_rate < -0.1:
                print(f"    => GEOMETRIC decay (low-rank approximation possible!)")
            else:
                print(f"    => MILD decay")

            # Fit log(S) vs log(i) for power-law decay
            if len(S_osc_nz) >= 4:
                idx_nz = idx[1:]  # skip i=0
                log_idx = np.log(idx_nz)
                log_s_nz = log_s[1:]
                coeffs_pl = np.polyfit(log_idx, log_s_nz, 1)
                print(f"    Power-law fit: log(S) = {coeffs_pl[1]:.4f} + {coeffs_pl[0]:.4f} * log(i)")
                print(f"    Power-law exponent = {coeffs_pl[0]:.4f}")

        # Distribution statistics
        print(f"\n  Distribution of oscillatory SVs:")
        print(f"    Mean = {np.mean(S_osc_nz):.6f}")
        print(f"    Std = {np.std(S_osc_nz):.6f}")
        print(f"    CV (std/mean) = {np.std(S_osc_nz)/np.mean(S_osc_nz):.4f}")
        print(f"    Skewness = {float(np.mean(((S_osc_nz - np.mean(S_osc_nz))/np.std(S_osc_nz))**3)):.4f}" if np.std(S_osc_nz) > 0 else "    Skewness = N/A")

        # Compare to Marchenko-Pastur (random matrix prediction)
        n_rows, n_cols = M_osc.shape
        gamma = n_rows / n_cols if n_cols > 0 else 1
        sigma2 = osc_var / (n_rows * n_cols) if n_rows * n_cols > 0 else 0
        mp_lower = sigma2 * (1 - np.sqrt(gamma))**2
        mp_upper = sigma2 * (1 + np.sqrt(gamma))**2
        print(f"\n  Marchenko-Pastur prediction (if random):")
        print(f"    sigma^2 = {sigma2:.6f}")
        print(f"    MP range: [{np.sqrt(mp_lower):.6f}, {np.sqrt(mp_upper):.6f}]")
        print(f"    Actual range: [{S_osc_nz[-1]:.6f}, {S_osc_nz[0]:.6f}]")

        # Cumulative variance of oscillatory part
        print(f"\n  Cumulative variance of oscillatory SVs:")
        cum_var = np.cumsum(S_osc_nz**2) / osc_var * 100
        for pct in [50, 75, 90, 95, 99]:
            idx_pct = np.searchsorted(cum_var, pct)
            if idx_pct < len(cum_var):
                print(f"    {pct}% captured by top {idx_pct+1}/{len(S_osc_nz)} osc SVs")

    # All SVs for detailed output
    return {
        'N': N,
        'rank': rank,
        'S': S,
        'S_osc': S_osc_nz,
        'total_var': total_var,
        'osc_var': osc_var
    }


def cross_N_analysis(results):
    """Compare spectral properties across different N values."""
    print(f"\n{'='*60}")
    print("CROSS-N COMPARISON")
    print(f"{'='*60}")

    print(f"\n{'N':>4} {'rank':>6} {'osc_rank':>8} {'max_osc':>10} {'min_osc':>10} {'ratio':>8} {'CV':>8} {'decay_rate':>12}")
    print("-" * 80)

    for r in results:
        N = r['N']
        S_osc = r['S_osc']
        if len(S_osc) > 0:
            ratio = S_osc[0] / S_osc[-1]
            cv = np.std(S_osc) / np.mean(S_osc) if np.mean(S_osc) > 0 else 0

            # Decay rate
            if len(S_osc) >= 2:
                idx = np.arange(len(S_osc))
                coeffs = np.polyfit(idx, np.log(S_osc + 1e-15), 1)
                dr = coeffs[0]
            else:
                dr = 0

            print(f"{N:>4} {r['rank']:>6} {len(S_osc):>8} {S_osc[0]:>10.4f} {S_osc[-1]:>10.4f} {ratio:>8.2f} {cv:>8.4f} {dr:>12.6f}")

    # Check if decay rate changes with N
    print(f"\nKey question: Does the max/min ratio grow with N?")
    print("  If ratio → const: oscillatory SVs are uniformly spread (random-like)")
    print("  If ratio → ∞: SVs concentrate, possible low-rank approx of oscillatory part")

    # Check scaling of max oscillatory SV
    Ns = [r['N'] for r in results if len(r['S_osc']) > 0]
    max_svs = [r['S_osc'][0] for r in results if len(r['S_osc']) > 0]
    if len(Ns) >= 3:
        log_N = np.log([2**n for n in Ns])
        log_sv = np.log(max_svs)
        coeffs = np.polyfit(log_N, log_sv, 1)
        print(f"\n  Max osc SV scales as x^{coeffs[0]:.4f} (x = 2^N)")
        print(f"  Expected: ~x^0 (constant) for random, ~x^{0.5} for structured")


def low_rank_approximation_error(N, pi_vals):
    """Test how well low-rank approximations work for the oscillatory part."""
    M = build_comm_matrix(N, pi_vals)
    U, S, Vt = np.linalg.svd(M, full_matrices=False)

    print(f"\nLow-rank approximation error for N={N}:")
    print(f"{'rank_k':>8} {'max_err':>10} {'mean_err':>10} {'exact_pct':>10}")

    total_entries = M.shape[0] * M.shape[1]

    for k in [2, 3, 4, 5, 6, 8, 10, 15, 20]:
        if k > len(S):
            break
        M_approx = sum(S[i] * np.outer(U[:, i], Vt[i, :]) for i in range(k))
        M_rounded = np.round(M_approx)

        errors = np.abs(M - M_rounded)
        max_err = np.max(errors)
        mean_err = np.mean(errors)
        exact = np.sum(errors < 0.5)

        print(f"{k:>8} {max_err:>10.4f} {mean_err:>10.6f} {100*exact/total_entries:>9.2f}%")


def main():
    # Sieve up to 2^20 + some buffer
    max_x = 2**20 + 100
    print(f"Sieving primes up to {max_x}...")
    is_prime = sieve(max_x)
    pi_vals = pi_func(is_prime)
    print(f"pi({max_x}) = {pi_vals[max_x-1]}")

    results = []

    for N in [4, 6, 8, 10, 12, 14, 16, 18, 20]:
        if 2**N > max_x:
            break
        r = analyze_svd(N, pi_vals)
        results.append(r)

    cross_N_analysis(results)

    # Low-rank approximation tests
    for N in [10, 12, 14, 16, 18, 20]:
        if 2**N > max_x:
            break
        low_rank_approximation_error(N, pi_vals)

    # Final: dump all oscillatory SVs for the largest feasible N
    largest = results[-1]
    N = largest['N']
    S_osc = largest['S_osc']
    print(f"\n{'='*60}")
    print(f"ALL oscillatory SVs for N={N} ({len(S_osc)} values):")
    for i, s in enumerate(S_osc):
        if i < 50 or i >= len(S_osc) - 5:
            print(f"  S_osc[{i}] = {s:.8f}")
        elif i == 50:
            print(f"  ... ({len(S_osc) - 55} more) ...")

if __name__ == '__main__':
    main()

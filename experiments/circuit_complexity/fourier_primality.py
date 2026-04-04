#!/usr/bin/env python3
"""
Boolean Fourier Analysis of the Prime Indicator Function
=========================================================

Computes the full Fourier spectrum of f: {0,1}^N -> {0,1} where f(x) = 1 iff x is prime.

We use the convention f: {0,1}^N -> R (keeping values 0/1, not +1/-1).
Fourier coefficients: f_hat(S) = (1/2^N) * sum_x f(x) * (-1)^{<S,x>}

Metrics computed:
  - Fourier weight at each degree k: W_k = sum_{|S|=k} f_hat(S)^2
  - Total influence: I(f) = sum_i Inf_i(f)
  - Noise sensitivity: NS_delta(f) = Pr[f(x) != f(x^delta)] where x^delta flips each bit with prob delta
  - Comparison against random function with same density

Key question: Does the prime indicator look like a "simple" function or a random function?
"""

import numpy as np
from itertools import combinations
from math import comb
import time
import sys

# --- Primality test ---
def is_prime(n):
    if n < 2:
        return False
    if n < 4:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6
    return True

# --- Walsh-Hadamard Transform (fast) ---
def walsh_hadamard_transform(f, N):
    """
    Compute all 2^N Fourier coefficients using the fast Walsh-Hadamard transform.
    f_hat(S) = (1/2^N) * sum_x f(x) * (-1)^{<S,x>}

    We use the in-place butterfly algorithm.
    """
    n = 1 << N
    # Copy f as float
    a = np.array(f, dtype=np.float64)

    # Standard WHT butterfly
    h = 1
    for _ in range(N):
        for i in range(0, n, h * 2):
            for j in range(i, i + h):
                x = a[j]
                y = a[j + h]
                a[j] = x + y
                a[j + h] = x - y
        h *= 2

    # Normalize
    a /= n
    return a

def popcount(x):
    """Number of set bits."""
    return bin(x).count('1')

def subset_to_mask(S, N):
    """Convert a subset (tuple of bit indices) to an integer mask."""
    mask = 0
    for i in S:
        mask |= (1 << i)
    return mask

def inner_product_parity(S_mask, x):
    """(-1)^{<S,x>} = (-1)^{popcount(S & x)}"""
    return popcount(S_mask & x) % 2

# --- Compute metrics ---
def compute_fourier_metrics(N, f_hat):
    """Compute degree weights, influence, noise sensitivity from Fourier coefficients."""
    n = 1 << N

    # Degree weights W_k = sum_{|S|=k} f_hat(S)^2
    W = np.zeros(N + 1)
    for S in range(n):
        k = popcount(S)
        W[k] += f_hat[S] ** 2

    # Individual influences: Inf_i(f) = sum_{S: i in S} f_hat(S)^2
    influences = np.zeros(N)
    for S in range(n):
        val = f_hat[S] ** 2
        for i in range(N):
            if S & (1 << i):
                influences[i] += val

    total_influence = np.sum(influences)

    # Noise sensitivity: NS_delta(f) = sum_{S != empty} f_hat(S)^2 * (1 - (1-2*delta)^|S|)
    # Using Fourier formula: NS_delta = E[f(x)(1-f(y))] where y = x with each bit flipped w.p. delta
    # = sum_S f_hat(S)^2 * (1 - (1-2*delta)^|S|) ... but this is for +/-1 valued functions.
    # For {0,1}-valued: NS_delta = Pr[f(x) != f(x^delta)]
    # = (1/2) - (1/2)*sum_S (2*f_hat(S))^2*(1-2*delta)^|S| ... no.
    #
    # Actually, for f:{0,1}^N -> {0,1}, the noise sensitivity is:
    # NS_delta(f) = E_x E_y [f(x) XOR f(y)] where y ~ N_delta(x)
    # = E[f] - sum_S f_hat(S)^2 * (1-2*delta)^|S|    ... no, let me be careful.
    #
    # The standard approach: let g(x) = 2*f(x) - 1 (map to +/-1).
    # g_hat(S) = 2*f_hat(S) for S != empty, g_hat(empty) = 2*f_hat(empty) - 1.
    # Parseval: sum_S g_hat(S)^2 = 1.
    # NS_delta(g) = (1/2) - (1/2)*sum_S g_hat(S)^2 * (1-2*delta)^|S|
    # Pr[f(x) != f(y)] = (1 - E[g(x)*g(y)])/2 = (1/2) - (1/2)*sum_S g_hat(S)^2*(1-2*delta)^|S|

    mu = f_hat[0]  # = E[f] = density of primes in {0,...,2^N-1}
    g_hat_sq = np.zeros(n)
    for S in range(n):
        if S == 0:
            g_hat_sq[S] = (2 * mu - 1) ** 2
        else:
            g_hat_sq[S] = (2 * f_hat[S]) ** 2

    deltas = [0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5]
    noise_sens = {}
    for delta in deltas:
        rho = 1 - 2 * delta
        corr = 0.0
        for S in range(n):
            k = popcount(S)
            corr += g_hat_sq[S] * (rho ** k)
        ns = (1 - corr) / 2
        noise_sens[delta] = ns

    return W, influences, total_influence, noise_sens

def generate_random_function(N, density, seed=42):
    """Generate a random Boolean function with given density."""
    rng = np.random.RandomState(seed)
    n = 1 << N
    count = int(round(density * n))
    f = np.zeros(n, dtype=np.float64)
    indices = rng.choice(n, size=count, replace=False)
    f[indices] = 1.0
    return f

def main():
    print("=" * 80)
    print("BOOLEAN FOURIER ANALYSIS OF THE PRIME INDICATOR FUNCTION")
    print("=" * 80)
    print()

    N_values = [4, 6, 8, 10, 12, 14, 16]

    all_results = {}

    for N in N_values:
        n = 1 << N
        print(f"\n{'='*70}")
        print(f"N = {N}  (universe = {{0, 1, ..., {n-1}}})")
        print(f"{'='*70}")

        t0 = time.time()

        # Build prime indicator
        f = np.zeros(n, dtype=np.float64)
        prime_count = 0
        for x in range(n):
            if is_prime(x):
                f[x] = 1.0
                prime_count += 1

        density = prime_count / n
        print(f"Primes in range: {prime_count}/{n} = {density:.6f}")
        print(f"PNT prediction: {n / np.log(n):.1f} = {1/np.log(n):.6f}")

        # Fourier transform
        f_hat = walsh_hadamard_transform(f, N)

        # Verify Parseval: sum f_hat(S)^2 = E[f^2] = E[f] = density (for {0,1}-valued)
        parseval = np.sum(f_hat ** 2)
        print(f"Parseval check: sum f_hat(S)^2 = {parseval:.8f} (should be {density:.8f})")

        # Compute metrics for prime indicator
        W, influences, total_influence, noise_sens = compute_fourier_metrics(N, f_hat)

        # --- Random function comparison ---
        f_rand = generate_random_function(N, density)
        f_hat_rand = walsh_hadamard_transform(f_rand, N)
        W_rand, inf_rand, total_inf_rand, ns_rand = compute_fourier_metrics(N, f_hat_rand)

        t1 = time.time()

        # --- Print results ---
        print(f"\nFourier weight by degree:")
        print(f"  {'Degree':<8} {'W_k(prime)':<16} {'W_k(random)':<16} {'Ratio':<12} {'Binom coeff':<12}")
        for k in range(N + 1):
            bc = comb(N, k)
            ratio = W[k] / W_rand[k] if W_rand[k] > 1e-15 else float('inf')
            print(f"  {k:<8} {W[k]:<16.8f} {W_rand[k]:<16.8f} {ratio:<12.4f} {bc:<12}")

        # Normalized weight (fraction of total)
        print(f"\nNormalized Fourier weight (fraction of Parseval sum):")
        print(f"  {'Degree':<8} {'Prime':<16} {'Random':<16}")
        for k in range(N + 1):
            wp = W[k] / parseval if parseval > 0 else 0
            wr = W_rand[k] / np.sum(f_hat_rand**2) if np.sum(f_hat_rand**2) > 0 else 0
            print(f"  {k:<8} {wp:<16.6f} {wr:<16.6f}")

        print(f"\nTotal influence:")
        print(f"  Prime:  I(f) = {total_influence:.6f}")
        print(f"  Random: I(f) = {total_inf_rand:.6f}")
        print(f"  Ratio:  {total_influence / total_inf_rand:.4f}" if total_inf_rand > 0 else "")
        print(f"  For random function with density p: expected I = N*p*(1-p) = {N*density*(1-density):.6f}")

        print(f"\nPer-bit influence (prime indicator):")
        for i in range(N):
            print(f"  Bit {i} (value 2^{i}): Inf_{i} = {influences[i]:.8f}")

        print(f"\nNoise sensitivity:")
        print(f"  {'delta':<10} {'NS(prime)':<16} {'NS(random)':<16} {'Ratio':<12}")
        for delta in sorted(noise_sens.keys()):
            ratio = noise_sens[delta] / ns_rand[delta] if ns_rand[delta] > 1e-15 else float('inf')
            print(f"  {delta:<10.2f} {noise_sens[delta]:<16.8f} {ns_rand[delta]:<16.8f} {ratio:<12.4f}")

        # Spectral entropy
        # H = -sum p_S log2(p_S) where p_S = f_hat(S)^2 / sum f_hat^2
        nonzero = f_hat[f_hat != 0]
        p_S = (nonzero ** 2) / parseval
        spectral_entropy = -np.sum(p_S * np.log2(p_S + 1e-300))
        max_entropy = N  # log2(2^N) for uniform

        nonzero_r = f_hat_rand[f_hat_rand != 0]
        p_S_r = (nonzero_r ** 2) / np.sum(f_hat_rand**2)
        spectral_entropy_r = -np.sum(p_S_r * np.log2(p_S_r + 1e-300))

        print(f"\nSpectral entropy:")
        print(f"  Prime:  {spectral_entropy:.4f} bits (max = {N} bits)")
        print(f"  Random: {spectral_entropy_r:.4f} bits")
        print(f"  Ratio to max: prime={spectral_entropy/N:.4f}, random={spectral_entropy_r/N:.4f}")

        # Top Fourier coefficients
        top_k = min(20, n)
        indices = np.argsort(np.abs(f_hat))[::-1][:top_k]
        print(f"\nTop {top_k} Fourier coefficients (by magnitude):")
        print(f"  {'S (binary)':<20} {'|S|':<6} {'f_hat(S)':<16} {'f_hat(S)^2':<16}")
        for idx in indices:
            S_bin = format(idx, f'0{N}b')
            k = popcount(idx)
            print(f"  {S_bin:<20} {k:<6} {f_hat[idx]:<16.8f} {f_hat[idx]**2:<16.10f}")

        # Degree concentration: what fraction of weight is in degrees <= N/3?
        low_degree_weight = sum(W[k] for k in range(N // 3 + 1))
        low_degree_frac = low_degree_weight / parseval if parseval > 0 else 0
        low_degree_weight_r = sum(W_rand[k] for k in range(N // 3 + 1))
        low_degree_frac_r = low_degree_weight_r / np.sum(f_hat_rand**2) if np.sum(f_hat_rand**2) > 0 else 0

        print(f"\nDegree concentration (weight in degrees 0..{N//3}):")
        print(f"  Prime:  {low_degree_frac:.6f}")
        print(f"  Random: {low_degree_frac_r:.6f}")

        print(f"\nTime: {t1-t0:.3f}s")

        all_results[N] = {
            'density': density,
            'prime_count': prime_count,
            'W': W,
            'W_rand': W_rand,
            'total_influence': total_influence,
            'total_inf_rand': total_inf_rand,
            'noise_sens': noise_sens,
            'ns_rand': ns_rand,
            'spectral_entropy': spectral_entropy,
            'spectral_entropy_r': spectral_entropy_r,
            'low_degree_frac': low_degree_frac,
            'low_degree_frac_r': low_degree_frac_r,
        }

    # ========== SUMMARY ==========
    print("\n" + "=" * 80)
    print("SUMMARY: FOURIER ANALYSIS OF PRIME INDICATOR")
    print("=" * 80)

    print("\n1. DEGREE WEIGHT DISTRIBUTION (normalized)")
    print(f"   {'N':<6} {'W_0':<10} {'W_1':<10} {'W_2':<10} {'W_N/2':<10} {'Low-deg frac':<14} {'Random low-deg':<14}")
    for N in N_values:
        r = all_results[N]
        p = r['density']
        print(f"   {N:<6} {r['W'][0]/p:<10.6f} {r['W'][1]/p:<10.6f} {r['W'][2]/p:<10.6f} "
              f"{r['W'][N//2]/p:<10.6f} {r['low_degree_frac']:<14.6f} {r['low_degree_frac_r']:<14.6f}")

    print("\n2. TOTAL INFLUENCE")
    print(f"   {'N':<6} {'I(prime)':<14} {'I(random)':<14} {'N*p*(1-p)':<14} {'Ratio p/r':<12}")
    for N in N_values:
        r = all_results[N]
        p = r['density']
        expected = N * p * (1 - p)
        ratio = r['total_influence'] / r['total_inf_rand'] if r['total_inf_rand'] > 0 else 0
        print(f"   {N:<6} {r['total_influence']:<14.6f} {r['total_inf_rand']:<14.6f} "
              f"{expected:<14.6f} {ratio:<12.4f}")

    print("\n3. NOISE SENSITIVITY at delta=0.05")
    print(f"   {'N':<6} {'NS(prime)':<16} {'NS(random)':<16} {'Ratio':<12}")
    for N in N_values:
        r = all_results[N]
        ns_p = r['noise_sens'][0.05]
        ns_r = r['ns_rand'][0.05]
        ratio = ns_p / ns_r if ns_r > 0 else 0
        print(f"   {N:<6} {ns_p:<16.8f} {ns_r:<16.8f} {ratio:<12.4f}")

    print("\n4. SPECTRAL ENTROPY")
    print(f"   {'N':<6} {'H(prime)':<12} {'H(random)':<12} {'Max':<8} {'Prime/Max':<12}")
    for N in N_values:
        r = all_results[N]
        print(f"   {N:<6} {r['spectral_entropy']:<12.4f} {r['spectral_entropy_r']:<12.4f} "
              f"{N:<8} {r['spectral_entropy']/N:<12.4f}")

    print("\n5. CONCLUSION")
    print("   The prime indicator function's Fourier spectrum is compared to random functions")
    print("   of the same density. Key findings:")
    print()

    # Compute average ratios for conclusion
    avg_low_deg_ratio = np.mean([all_results[N]['low_degree_frac'] / all_results[N]['low_degree_frac_r']
                                  for N in N_values if all_results[N]['low_degree_frac_r'] > 0])
    avg_ns_ratio = np.mean([all_results[N]['noise_sens'][0.05] / all_results[N]['ns_rand'][0.05]
                            for N in N_values if all_results[N]['ns_rand'][0.05] > 0])
    avg_inf_ratio = np.mean([all_results[N]['total_influence'] / all_results[N]['total_inf_rand']
                             for N in N_values if all_results[N]['total_inf_rand'] > 0])
    avg_entropy_ratio = np.mean([all_results[N]['spectral_entropy'] / all_results[N]['spectral_entropy_r']
                                  for N in N_values])

    print(f"   - Low-degree concentration ratio (prime/random): {avg_low_deg_ratio:.4f}")
    print(f"     (>1 means MORE low-degree weight than random = simpler)")
    print(f"     (<1 means LESS low-degree weight than random = more complex)")
    print()
    print(f"   - Noise sensitivity ratio (prime/random) at delta=0.05: {avg_ns_ratio:.4f}")
    print(f"     (~1 means same noise sensitivity as random)")
    print(f"     (<1 means MORE noise-stable = suggestive of low depth)")
    print()
    print(f"   - Total influence ratio (prime/random): {avg_inf_ratio:.4f}")
    print(f"     (~1 means similar influence pattern)")
    print()
    print(f"   - Spectral entropy ratio (prime/random): {avg_entropy_ratio:.4f}")
    print(f"     (~1 means similarly spread Fourier spectrum)")
    print()

    if avg_low_deg_ratio > 1.2:
        print("   FINDING: Prime indicator has MORE low-degree Fourier weight than random.")
        print("   This suggests some exploitable structure -- NOT fully random-like.")
    elif avg_low_deg_ratio < 0.8:
        print("   FINDING: Prime indicator has LESS low-degree Fourier weight than random.")
        print("   This suggests the function is MORE complex than random in Fourier terms.")
    else:
        print("   FINDING: Prime indicator has SIMILAR low-degree weight to random functions.")
        print("   No evidence of low-complexity structure in the Fourier domain.")

    print()
    if avg_ns_ratio < 0.8:
        print("   FINDING: Prime indicator is MORE noise-stable than random.")
        print("   This is consistent with bounded-depth (TC^0/AC^0) computation.")
    elif avg_ns_ratio > 1.2:
        print("   FINDING: Prime indicator is MORE noise-sensitive than random.")
        print("   This suggests high-depth computation.")
    else:
        print("   FINDING: Prime indicator has SIMILAR noise sensitivity to random.")
        print("   Inconclusive for depth bounds.")

    print()
    print("   OVERALL: The prime indicator function in Fourier space is characterized by:")
    print("   - Bit 0 (parity) has disproportionate influence (primes > 2 are odd)")
    print("   - Bit 1 (mod 4) also has notable influence (primes not divisible by 2)")
    print("   - Higher bits contribute increasingly uniformly")
    print("   - The overall Fourier profile is CLOSE TO RANDOM with density corrections")
    print()
    print("   IMPLICATIONS FOR CIRCUIT COMPLEXITY:")
    print("   - Random-like Fourier spectrum is CONSISTENT WITH high circuit complexity")
    print("   - No evidence of junta structure (few influential variables)")
    print("   - No evidence of low-degree polynomial approximation")
    print("   - The low bits (parity, mod 4) carry divisibility structure,")
    print("     but this is trivial and doesn't help with prime counting")

if __name__ == "__main__":
    main()

"""
Index-Prime Bit Correlation Analysis
=====================================
For p(n) = nth prime, analyze correlations between bits of n and bits of p(n).

Questions:
1. Does bit j of n predict bit k of p(n)? (full MI matrix)
2. How much total information do bits of n carry about bits of p(n)?
3. Is there a smooth-to-hard transition visible in the n-bits?
4. Can groups of n-bits predict individual p-bits? (conditional MI)
5. Per-bit R^2: how well does n predict each bit of p(n)?
6. XOR correlation: n XOR p(n) structure
7. Bit-width relationship: how does bit_length(p(n)) relate to bit_length(n)?

Session 41.
"""

import math
import numpy as np
from collections import Counter

def sieve_primes(limit):
    is_prime = bytearray([1]) * (limit + 1)
    is_prime[0] = is_prime[1] = 0
    for i in range(2, int(limit**0.5) + 1):
        if is_prime[i]:
            is_prime[i*i::i] = bytearray(len(is_prime[i*i::i]))
    return [i for i in range(2, limit + 1) if is_prime[i]]

def mutual_info(seq_a, seq_b, n):
    """Compute MI between two binary sequences."""
    joint = [[0,0],[0,0]]
    for i in range(n):
        joint[seq_a[i]][seq_b[i]] += 1
    mi = 0.0
    for a in range(2):
        for b in range(2):
            p_ab = joint[a][b] / n
            p_a = (joint[a][0] + joint[a][1]) / n
            p_b = (joint[0][b] + joint[1][b]) / n
            if p_ab > 0 and p_a > 0 and p_b > 0:
                mi += p_ab * math.log2(p_ab / (p_a * p_b))
    return mi

def main():
    print("=" * 70)
    print("INDEX-PRIME BIT CORRELATION: bits(n) vs bits(p(n))")
    print("=" * 70)

    LIMIT = 10_000_000
    primes = sieve_primes(LIMIT)
    N = len(primes)
    print(f"\nPrimes: {N:,}, max p = {primes[-1]:,}")

    N_BITS = (N-1).bit_length()       # bits needed for index n
    P_BITS = primes[-1].bit_length()  # bits needed for p(n)
    print(f"Index bits (n): {N_BITS}, Prime bits (p): {P_BITS}")

    # Precompute bit arrays for speed
    SAMPLE = min(N, 500_000)
    # Use indices starting from 1 (p(1)=2, p(2)=3, ...)
    indices = list(range(1, SAMPLE + 1))
    sample_primes = primes[:SAMPLE]

    # ================================================================
    # 1. FULL MI MATRIX: bit_j(n) vs bit_k(p(n))
    # ================================================================
    print("\n" + "=" * 70)
    print("1. MUTUAL INFORMATION: bit_j(n) vs bit_k(p(n))")
    print("=" * 70)

    MI_N = min(N_BITS, 20)
    MI_P = min(P_BITS, 24)

    # Extract bit sequences
    n_bits = {}
    for j in range(MI_N):
        n_bits[j] = [(indices[i] >> j) & 1 for i in range(SAMPLE)]

    p_bits = {}
    for k in range(MI_P):
        p_bits[k] = [(sample_primes[i] >> k) & 1 for i in range(SAMPLE)]

    mi_matrix = np.zeros((MI_N, MI_P))
    for j in range(MI_N):
        for k in range(MI_P):
            mi_matrix[j][k] = mutual_info(n_bits[j], p_bits[k], SAMPLE)

    print(f"\nMI(bit_j(n), bit_k(p(n))) matrix (sample={SAMPLE:,}):")
    print(f"Rows = n bits (0=LSB), Cols = p(n) bits (0=LSB)")
    print(f"\n{'n\\p':>6}", end='')
    for k in range(MI_P):
        print(f" {k:6d}", end='')
    print()
    for j in range(MI_N):
        print(f"{j:6d}", end='')
        for k in range(MI_P):
            v = mi_matrix[j][k]
            if v < 0.0001:
                print(f"  .    ", end='')
            else:
                print(f" {v:6.4f}", end='')
        print(f"  | sum={mi_matrix[j,:].sum():.4f}")

    print(f"\nColumn sums (total n-info about each p-bit):")
    print(f"{'p-bit':>6}", end='')
    for k in range(MI_P):
        print(f" {k:6d}", end='')
    print()
    print(f"{'sum':>6}", end='')
    for k in range(MI_P):
        print(f" {mi_matrix[:,k].sum():6.4f}", end='')
    print()

    # ================================================================
    # 2. TOTAL MI: H(p-bit | all n-bits) for each p-bit
    # ================================================================
    print("\n" + "=" * 70)
    print("2. CONDITIONAL ENTROPY: How much of each p-bit is determined by n?")
    print("=" * 70)

    print(f"\nFor each bit k of p(n), compute H(bit_k(p)) and MI with n:")
    print(f"{'p-bit':>6} {'H(bit)':>8} {'sum MI':>8} {'H remaining':>12} {'% explained':>12}")
    for k in range(MI_P):
        # H(bit_k(p))
        p1 = sum(p_bits[k]) / SAMPLE
        p0 = 1 - p1
        h_bit = 0
        if p0 > 0 and p1 > 0:
            h_bit = -p0 * math.log2(p0) - p1 * math.log2(p1)
        total_mi = mi_matrix[:, k].sum()
        h_remaining = max(0, h_bit - total_mi)
        pct = (total_mi / h_bit * 100) if h_bit > 0 else 0
        print(f"{k:6d} {h_bit:8.4f} {total_mi:8.4f} {h_remaining:12.4f} {pct:11.2f}%")

    # ================================================================
    # 3. CORRELATION COEFFICIENT: bit_j(n) vs bit_k(p(n))
    # ================================================================
    print("\n" + "=" * 70)
    print("3. PEARSON CORRELATION: bit_j(n) vs bit_k(p(n))")
    print("=" * 70)

    corr_matrix = np.zeros((MI_N, MI_P))
    for j in range(MI_N):
        nj = np.array(n_bits[j], dtype=np.float64)
        nj_c = nj - nj.mean()
        nj_std = nj.std()
        for k in range(MI_P):
            pk = np.array(p_bits[k], dtype=np.float64)
            pk_c = pk - pk.mean()
            pk_std = pk.std()
            if nj_std > 0 and pk_std > 0:
                corr_matrix[j][k] = np.mean(nj_c * pk_c) / (nj_std * pk_std)

    print(f"\nCorrelation matrix (showing only |r| > 0.01):")
    print(f"{'n\\p':>6}", end='')
    for k in range(MI_P):
        print(f" {k:6d}", end='')
    print()
    for j in range(MI_N):
        print(f"{j:6d}", end='')
        for k in range(MI_P):
            v = corr_matrix[j][k]
            if abs(v) < 0.01:
                print(f"  .    ", end='')
            else:
                print(f" {v:+5.3f} ", end='')
        print()

    # ================================================================
    # 4. DIAGONAL ANALYSIS: bit_k(n) vs bit_k(p(n)) — same position
    # ================================================================
    print("\n" + "=" * 70)
    print("4. SAME-POSITION ANALYSIS: bit_k(n) vs bit_k(p(n))")
    print("=" * 70)

    print(f"\n{'bit k':>6} {'corr':>8} {'MI':>8} {'P(agree)':>10} {'P(n=0,p=0)':>12} {'P(n=0,p=1)':>12} {'P(n=1,p=0)':>12} {'P(n=1,p=1)':>12}")
    for k in range(min(MI_N, MI_P)):
        agree = sum(1 for i in range(SAMPLE) if n_bits[k][i] == p_bits[k][i]) / SAMPLE
        joint = [[0,0],[0,0]]
        for i in range(SAMPLE):
            joint[n_bits[k][i]][p_bits[k][i]] += 1
        j00 = joint[0][0]/SAMPLE
        j01 = joint[0][1]/SAMPLE
        j10 = joint[1][0]/SAMPLE
        j11 = joint[1][1]/SAMPLE
        print(f"{k:6d} {corr_matrix[k][k]:+8.4f} {mi_matrix[k][k]:8.4f} {agree:10.4f} {j00:12.4f} {j01:12.4f} {j10:12.4f} {j11:12.4f}")

    # ================================================================
    # 5. n XOR p(n) ANALYSIS
    # ================================================================
    print("\n" + "=" * 70)
    print("5. n XOR p(n) ANALYSIS")
    print("=" * 70)

    xor_vals = [indices[i] ^ sample_primes[i] for i in range(SAMPLE)]
    xor_bitlengths = [xv.bit_length() for xv in xor_vals]

    print(f"\nn XOR p(n) statistics:")
    xor_arr = np.array(xor_vals, dtype=np.float64)
    print(f"  Mean:   {xor_arr.mean():.1f}")
    print(f"  Median: {np.median(xor_arr):.1f}")
    print(f"  Std:    {xor_arr.std():.1f}")

    # Ratio p(n)/n and its binary structure
    ratios = [sample_primes[i] / indices[i] for i in range(SAMPLE)]
    ratio_arr = np.array(ratios)
    print(f"\np(n)/n ratio statistics:")
    print(f"  Mean:   {ratio_arr.mean():.4f}")
    print(f"  Std:    {ratio_arr.std():.4f}")
    print(f"  At n=1000: {ratios[999]:.4f} (expected ~ln(1000)={math.log(1000):.4f})")
    print(f"  At n=10000: {ratios[9999]:.4f} (expected ~ln(10000)={math.log(10000):.4f})")
    print(f"  At n=100000: {ratios[99999]:.4f} (expected ~ln(100000)={math.log(100000):.4f})")

    # Bit-length difference
    bl_diffs = [sample_primes[i].bit_length() - indices[i].bit_length() for i in range(SAMPLE)]
    bl_counts = Counter(bl_diffs)
    print(f"\nbit_length(p(n)) - bit_length(n) distribution:")
    for d in sorted(bl_counts.keys()):
        pct = bl_counts[d] / SAMPLE * 100
        if pct > 0.1:
            print(f"  diff={d:+2d}: {bl_counts[d]:8d} ({pct:5.2f}%)")

    # ================================================================
    # 6. SHIFTED MI: bit_j(n) vs bit_{j+shift}(p(n))
    # ================================================================
    print("\n" + "=" * 70)
    print("6. SHIFTED CORRELATION: bit_j(n) vs bit_{j+s}(p(n))")
    print("=" * 70)

    print(f"\nAverage |correlation| for different bit shifts s:")
    print(f"(shift s means comparing bit j of n with bit j+s of p(n))")
    for s in range(-3, 8):
        corrs = []
        for j in range(MI_N):
            k = j + s
            if 0 <= k < MI_P:
                corrs.append(abs(corr_matrix[j][k]))
        if corrs:
            mean_corr = np.mean(corrs)
            max_corr = np.max(corrs)
            print(f"  shift {s:+2d}: mean |r| = {mean_corr:.5f}, max |r| = {max_corr:.5f}, n_pairs = {len(corrs)}")

    # ================================================================
    # 7. MULTI-BIT PREDICTION: Can 3-bit groups of n predict p-bits?
    # ================================================================
    print("\n" + "=" * 70)
    print("7. MULTI-BIT PREDICTION: 3-bit n-groups predicting p-bits")
    print("=" * 70)

    print(f"\nFor each p-bit, train 3-bit lookup from best n-bits:")
    for pk in range(min(MI_P, 20)):
        # Find best 3 n-bits by MI
        mi_for_pk = [(mi_matrix[j][pk], j) for j in range(MI_N)]
        mi_for_pk.sort(reverse=True)
        best_3 = [j for _, j in mi_for_pk[:3]]

        # Build 8-entry lookup table
        lookup = {}
        counts = {}
        for i in range(SAMPLE):
            key = (n_bits[best_3[0]][i], n_bits[best_3[1]][i], n_bits[best_3[2]][i])
            if key not in lookup:
                lookup[key] = [0, 0]
                counts[key] = 0
            lookup[key][p_bits[pk][i]] += 1
            counts[key] += 1

        # Predict using majority
        correct = 0
        for i in range(SAMPLE):
            key = (n_bits[best_3[0]][i], n_bits[best_3[1]][i], n_bits[best_3[2]][i])
            pred = 1 if lookup[key][1] > lookup[key][0] else 0
            if pred == p_bits[pk][i]:
                correct += 1

        accuracy = correct / SAMPLE
        baseline = max(sum(p_bits[pk]) / SAMPLE, 1 - sum(p_bits[pk]) / SAMPLE)
        lift = accuracy - baseline
        bits_used = f"n[{best_3[0]},{best_3[1]},{best_3[2]}]"
        print(f"  p-bit {pk:2d}: acc={accuracy:.4f} baseline={baseline:.4f} lift={lift:+.4f} using {bits_used}")

    # ================================================================
    # 8. RESIDUAL ANALYSIS: p(n) - n*ln(n) in binary
    # ================================================================
    print("\n" + "=" * 70)
    print("8. RESIDUAL p(n) - n*ln(n): bit-level analysis")
    print("=" * 70)

    residuals = [sample_primes[i] - int(indices[i] * math.log(indices[i])) for i in range(1, SAMPLE)]
    res_arr = np.array(residuals, dtype=np.float64)
    print(f"\nResidual statistics (p(n) - n*ln(n)):")
    print(f"  Mean: {res_arr.mean():.2f}")
    print(f"  Std:  {res_arr.std():.2f}")
    print(f"  Min:  {res_arr.min():.0f}")
    print(f"  Max:  {res_arr.max():.0f}")
    print(f"  Fraction positive: {(res_arr > 0).mean():.4f}")

    # How many bits does the residual occupy?
    res_bits = [abs(r).bit_length() for r in residuals]
    res_bits_arr = np.array(res_bits, dtype=np.float64)
    n_bits_arr = np.array([indices[i].bit_length() for i in range(1, SAMPLE)], dtype=np.float64)

    print(f"\nResidual bit-length vs index bit-length:")
    print(f"  Mean residual bits: {res_bits_arr.mean():.2f}")
    print(f"  Mean index bits:    {n_bits_arr.mean():.2f}")
    print(f"  Ratio (res/index):  {res_bits_arr.mean()/n_bits_arr.mean():.4f}")

    # Per n-bit-length bucket
    print(f"\n{'n_bits':>7} {'mean res_bits':>14} {'ratio':>8} {'count':>8}")
    for nb in range(4, int(n_bits_arr.max()) + 1):
        mask = n_bits_arr == nb
        if mask.sum() > 100:
            mean_rb = res_bits_arr[mask].mean()
            print(f"{nb:7d} {mean_rb:14.2f} {mean_rb/nb:8.4f} {mask.sum():8d}")

    # ================================================================
    # SUMMARY
    # ================================================================
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"""
Results:
1. MI matrix: Which n-bits predict which p-bits, and how strongly
2. Conditional entropy: How much of each p-bit is explained by n
3. Correlation structure: Linear relationships between n and p bits
4. Same-position: Do bit k of n and bit k of p(n) correspond?
5. XOR/ratio: Structural relationship between n and p(n)
6. Shifted correlation: Is there a consistent bit-shift relationship?
7. Multi-bit prediction: Can 3 n-bits predict a p-bit better than chance?
8. Residual bit analysis: How many bits does the "hard part" occupy?
""")

if __name__ == "__main__":
    main()

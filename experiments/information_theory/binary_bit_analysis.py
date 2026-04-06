"""
Binary Bit Analysis of Consecutive Primes
==========================================
Examines unstudied binary perspectives:
1. XOR structure: p(n) XOR p(n+1) — where do bits flip?
2. Hamming distance patterns between consecutive primes
3. Per-bit-position statistics — P(bit_k = 1) across all primes
4. Run-length structure in binary representations
5. Bit transition matrices — Markov model per bit position
6. Spectral analysis per bit position — FFT of bit_k(p(n)) as n varies
7. Bit-position mutual information — which bits predict which?

Session 41: Novel binary perspective analysis.
"""

import math
import sys
import numpy as np
from collections import Counter, defaultdict

def sieve_primes(limit):
    """Simple sieve of Eratosthenes."""
    is_prime = bytearray([1]) * (limit + 1)
    is_prime[0] = is_prime[1] = 0
    for i in range(2, int(limit**0.5) + 1):
        if is_prime[i]:
            is_prime[i*i::i] = bytearray(len(is_prime[i*i::i]))
    return [i for i in range(2, limit + 1) if is_prime[i]]

def to_bits(n, width):
    """Convert integer to binary array (MSB first)."""
    return [(n >> (width - 1 - i)) & 1 for i in range(width)]

def main():
    print("=" * 70)
    print("BINARY BIT ANALYSIS OF CONSECUTIVE PRIMES")
    print("=" * 70)

    # Generate primes
    LIMIT = 10_000_000
    primes = sieve_primes(LIMIT)
    N = len(primes)
    print(f"\nPrimes up to {LIMIT:,}: {N:,}")
    print(f"Largest prime: {primes[-1]:,} ({primes[-1].bit_length()} bits)")

    BIT_WIDTH = primes[-1].bit_length()  # max bits needed
    print(f"Bit width for analysis: {BIT_WIDTH}")

    # ================================================================
    # 1. XOR STRUCTURE: p(n) XOR p(n+1)
    # ================================================================
    print("\n" + "=" * 70)
    print("1. XOR STRUCTURE: p(n) XOR p(n+1)")
    print("=" * 70)

    xor_vals = [primes[i] ^ primes[i+1] for i in range(N-1)]

    # Which bits flip most often?
    bit_flip_count = [0] * BIT_WIDTH
    for xv in xor_vals:
        for b in range(BIT_WIDTH):
            if (xv >> b) & 1:
                bit_flip_count[b] += 1

    print(f"\nPer-bit flip frequency (bit 0 = LSB):")
    print(f"{'Bit':>4} {'Flips':>10} {'Rate':>8} {'Expected(random)':>16}")
    for b in range(min(BIT_WIDTH, 24)):
        rate = bit_flip_count[b] / (N - 1)
        print(f"{b:4d} {bit_flip_count[b]:10d} {rate:8.4f} {'~0.50':>16}")

    # XOR bit length distribution (how many bits change?)
    xor_bitlengths = [xv.bit_length() for xv in xor_vals]
    print(f"\nXOR bit-length (highest flipped bit) distribution:")
    bl_counts = Counter(xor_bitlengths)
    for bl in sorted(bl_counts.keys()):
        pct = bl_counts[bl] / (N-1) * 100
        if pct > 0.1:
            print(f"  {bl:2d} bits: {bl_counts[bl]:8d} ({pct:5.2f}%)")

    # ================================================================
    # 2. HAMMING DISTANCE PATTERNS
    # ================================================================
    print("\n" + "=" * 70)
    print("2. HAMMING DISTANCE BETWEEN CONSECUTIVE PRIMES")
    print("=" * 70)

    hamming_dists = [bin(xv).count('1') for xv in xor_vals]
    hamming_arr = np.array(hamming_dists, dtype=np.float64)

    print(f"\nHamming distance statistics:")
    print(f"  Mean:   {hamming_arr.mean():.4f}")
    print(f"  Median: {np.median(hamming_arr):.1f}")
    print(f"  Std:    {hamming_arr.std():.4f}")
    print(f"  Min:    {hamming_arr.min():.0f}")
    print(f"  Max:    {hamming_arr.max():.0f}")

    print(f"\nHamming distance distribution:")
    hd_counts = Counter(hamming_dists)
    for hd in sorted(hd_counts.keys()):
        pct = hd_counts[hd] / (N-1) * 100
        bar = '#' * int(pct * 2)
        print(f"  {hd:2d}: {hd_counts[hd]:8d} ({pct:5.2f}%) {bar}")

    # Autocorrelation of Hamming distances
    hd_centered = hamming_arr - hamming_arr.mean()
    hd_var = np.var(hd_centered)
    print(f"\nHamming distance autocorrelation:")
    for lag in [1, 2, 3, 5, 10, 20, 50, 100]:
        if lag < len(hd_centered):
            ac = np.mean(hd_centered[:-lag] * hd_centered[lag:]) / hd_var
            print(f"  lag {lag:3d}: {ac:+.6f}")

    # ================================================================
    # 3. PER-BIT-POSITION STATISTICS
    # ================================================================
    print("\n" + "=" * 70)
    print("3. PER-BIT-POSITION STATISTICS: P(bit_k = 1) ACROSS PRIMES")
    print("=" * 70)

    bit_ones = [0] * BIT_WIDTH
    for p in primes:
        for b in range(BIT_WIDTH):
            if (p >> b) & 1:
                bit_ones[b] += 1

    print(f"\n{'Bit':>4} {'P(1)':>8} {'Deviation from 0.5':>20} {'Significance':>14}")
    for b in range(min(BIT_WIDTH, 24)):
        p1 = bit_ones[b] / N
        dev = p1 - 0.5
        # z-score for binomial
        z = dev / (0.5 / math.sqrt(N))
        sig = "***" if abs(z) > 3.29 else "**" if abs(z) > 2.58 else "*" if abs(z) > 1.96 else ""
        print(f"{b:4d} {p1:8.5f} {dev:+20.5f} {z:+8.2f} {sig}")

    # ================================================================
    # 4. RUN-LENGTH STRUCTURE
    # ================================================================
    print("\n" + "=" * 70)
    print("4. RUN-LENGTH STRUCTURE IN BINARY REPRESENTATIONS")
    print("=" * 70)

    # Analyze run lengths of 0s and 1s in binary representations
    run_lengths_0 = []
    run_lengths_1 = []
    sample_size = min(N, 500_000)  # sample for speed

    for p in primes[:sample_size]:
        bits = bin(p)[2:]  # binary string without '0b'
        current = bits[0]
        run = 1
        for c in bits[1:]:
            if c == current:
                run += 1
            else:
                if current == '0':
                    run_lengths_0.append(run)
                else:
                    run_lengths_1.append(run)
                current = c
                run = 1
        if current == '0':
            run_lengths_0.append(run)
        else:
            run_lengths_1.append(run)

    rl0 = np.array(run_lengths_0, dtype=np.float64)
    rl1 = np.array(run_lengths_1, dtype=np.float64)

    print(f"\nRun-length statistics (sample of {sample_size:,} primes):")
    print(f"  0-runs: mean={rl0.mean():.3f}, std={rl0.std():.3f}, max={rl0.max():.0f}")
    print(f"  1-runs: mean={rl1.mean():.3f}, std={rl1.std():.3f}, max={rl1.max():.0f}")

    # For truly random binary strings, run length follows geometric(0.5)
    # mean = 2, P(run >= k) = 2^{-(k-1)}
    print(f"\n  Expected for random: mean=2.000, geometric distribution")

    print(f"\n  0-run distribution vs geometric:")
    rl0_counts = Counter(run_lengths_0)
    total_rl0 = len(run_lengths_0)
    for k in range(1, 12):
        observed = rl0_counts.get(k, 0) / total_rl0
        expected = 0.5 ** k  # geometric
        ratio = observed / expected if expected > 0 else 0
        print(f"    len={k:2d}: obs={observed:.5f} exp={expected:.5f} ratio={ratio:.4f}")

    print(f"\n  1-run distribution vs geometric:")
    rl1_counts = Counter(run_lengths_1)
    total_rl1 = len(run_lengths_1)
    for k in range(1, 12):
        observed = rl1_counts.get(k, 0) / total_rl1
        expected = 0.5 ** k
        ratio = observed / expected if expected > 0 else 0
        print(f"    len={k:2d}: obs={observed:.5f} exp={expected:.5f} ratio={ratio:.4f}")

    # ================================================================
    # 5. BIT TRANSITION MATRICES (Markov model per bit position)
    # ================================================================
    print("\n" + "=" * 70)
    print("5. BIT TRANSITION MATRICES: P(bit_k(p(n+1)) | bit_k(p(n)))")
    print("=" * 70)

    print(f"\n{'Bit':>4} {'P(0→0)':>8} {'P(0→1)':>8} {'P(1→0)':>8} {'P(1→1)':>8} {'Persistence':>12} {'Chi2':>10}")
    for b in range(min(BIT_WIDTH, 24)):
        trans = [[0,0],[0,0]]
        for i in range(N-1):
            b_curr = (primes[i] >> b) & 1
            b_next = (primes[i+1] >> b) & 1
            trans[b_curr][b_next] += 1

        total = N - 1
        t00 = trans[0][0] / max(trans[0][0] + trans[0][1], 1)
        t01 = trans[0][1] / max(trans[0][0] + trans[0][1], 1)
        t10 = trans[1][0] / max(trans[1][0] + trans[1][1], 1)
        t11 = trans[1][1] / max(trans[1][0] + trans[1][1], 1)
        persistence = (trans[0][0] + trans[1][1]) / total
        # Chi-squared for independence
        n00, n01, n10, n11 = trans[0][0], trans[0][1], trans[1][0], trans[1][1]
        row0 = n00 + n01
        row1 = n10 + n11
        col0 = n00 + n10
        col1 = n01 + n11
        chi2 = 0
        for r, c, obs in [(row0,col0,n00),(row0,col1,n01),(row1,col0,n10),(row1,col1,n11)]:
            exp = r * c / total
            if exp > 0:
                chi2 += (obs - exp)**2 / exp
        print(f"{b:4d} {t00:8.4f} {t01:8.4f} {t10:8.4f} {t11:8.4f} {persistence:12.4f} {chi2:10.1f}")

    # ================================================================
    # 6. SPECTRAL ANALYSIS PER BIT POSITION
    # ================================================================
    print("\n" + "=" * 70)
    print("6. SPECTRAL ANALYSIS: FFT OF bit_k(p(n)) SEQUENCES")
    print("=" * 70)

    # Use first 2^17 = 131072 primes for clean FFT
    FFT_N = min(2**17, N)
    print(f"\nUsing first {FFT_N:,} primes for FFT analysis")

    print(f"\n{'Bit':>4} {'DC (mean)':>10} {'Peak freq':>10} {'Peak power':>11} {'Spectral flatness':>18} {'Entropy':>8}")
    for b in range(min(BIT_WIDTH, 20)):
        bit_seq = np.array([(primes[i] >> b) & 1 for i in range(FFT_N)], dtype=np.float64)
        bit_seq -= bit_seq.mean()  # remove DC

        fft_vals = np.fft.rfft(bit_seq)
        power = np.abs(fft_vals[1:])**2  # exclude DC
        power_norm = power / power.sum() if power.sum() > 0 else power

        # Spectral flatness: geometric mean / arithmetic mean
        log_power = np.log(power_norm + 1e-30)
        geo_mean = np.exp(np.mean(log_power))
        arith_mean = np.mean(power_norm)
        flatness = geo_mean / arith_mean if arith_mean > 0 else 0

        # Spectral entropy
        p_spec = power_norm[power_norm > 0]
        entropy = -np.sum(p_spec * np.log2(p_spec))
        max_entropy = np.log2(len(power_norm))

        peak_idx = np.argmax(power) + 1
        peak_freq = peak_idx / FFT_N

        print(f"{b:4d} {bit_seq.mean() + ((primes[0] >> b) & 1):10.5f} {peak_freq:10.6f} {power[peak_idx-1]:11.1f} {flatness:18.6f} {entropy/max_entropy:8.4f}")

    # ================================================================
    # 7. INTER-BIT MUTUAL INFORMATION
    # ================================================================
    print("\n" + "=" * 70)
    print("7. INTER-BIT MUTUAL INFORMATION: I(bit_j; bit_k) WITHIN p(n)")
    print("=" * 70)

    MI_BITS = min(16, BIT_WIDTH)
    sample = min(N, 500_000)

    # Build joint distributions
    mi_matrix = np.zeros((MI_BITS, MI_BITS))
    for j in range(MI_BITS):
        for k in range(j+1, MI_BITS):
            joint = [[0,0],[0,0]]
            for i in range(sample):
                bj = (primes[i] >> j) & 1
                bk = (primes[i] >> k) & 1
                joint[bj][bk] += 1

            # MI calculation
            mi = 0.0
            for a in range(2):
                for b_val in range(2):
                    p_ab = joint[a][b_val] / sample
                    p_a = (joint[a][0] + joint[a][1]) / sample
                    p_b = (joint[0][b_val] + joint[1][b_val]) / sample
                    if p_ab > 0 and p_a > 0 and p_b > 0:
                        mi += p_ab * math.log2(p_ab / (p_a * p_b))
            mi_matrix[j][k] = mi
            mi_matrix[k][j] = mi

    print(f"\nMutual information matrix (bits), sample={sample:,}:")
    print(f"{'':>4}", end='')
    for k in range(MI_BITS):
        print(f" {k:6d}", end='')
    print()
    for j in range(MI_BITS):
        print(f"{j:4d}", end='')
        for k in range(MI_BITS):
            if j == k:
                print(f"   --- ", end='')
            else:
                print(f" {mi_matrix[j][k]:6.4f}", end='')
        print()

    # ================================================================
    # 8. XOR PREDICTABILITY: Can bit_k(p(n) XOR p(n+1)) be predicted?
    # ================================================================
    print("\n" + "=" * 70)
    print("8. XOR PREDICTABILITY: AUTOCORRELATION OF XOR BITS")
    print("=" * 70)

    print(f"\nAutocorrelation of bit_k(p(n) XOR p(n+1)) at lags 1-5:")
    print(f"{'Bit':>4}", end='')
    for lag in [1, 2, 3, 5, 10]:
        print(f" {'lag'+str(lag):>8}", end='')
    print(f" {'verdict':>10}")

    AC_N = min(N-1, 200_000)
    for b in range(min(BIT_WIDTH, 20)):
        xor_bits = np.array([(xor_vals[i] >> b) & 1 for i in range(AC_N)], dtype=np.float64)
        xor_centered = xor_bits - xor_bits.mean()
        xor_var = np.var(xor_centered)

        acs = []
        for lag in [1, 2, 3, 5, 10]:
            if xor_var > 0:
                ac = np.mean(xor_centered[:-lag] * xor_centered[lag:]) / xor_var
            else:
                ac = 0.0
            acs.append(ac)

        # Is any autocorrelation significant? (threshold ~ 2/sqrt(N))
        threshold = 2.0 / math.sqrt(AC_N)
        significant = any(abs(ac) > threshold for ac in acs)
        verdict = "STRUCT" if significant else "RANDOM"

        print(f"{b:4d}", end='')
        for ac in acs:
            print(f" {ac:+8.5f}", end='')
        print(f" {verdict:>10}")

    print(f"\n  Significance threshold (2/sqrt(N)): {2/math.sqrt(AC_N):.6f}")

    # ================================================================
    # 9. CONSECUTIVE PRIME BIT DIFFERENCE PATTERNS
    # ================================================================
    print("\n" + "=" * 70)
    print("9. GAP-CONDITIONED BIT ANALYSIS")
    print("=" * 70)

    # For primes with gap=2 (twin primes), gap=4, gap=6, etc.
    # How do bits change?
    gap_bit_flips = defaultdict(lambda: [0] * BIT_WIDTH)
    gap_counts = Counter()

    for i in range(N-1):
        gap = primes[i+1] - primes[i]
        if gap <= 30:
            gap_counts[gap] += 1
            xv = primes[i] ^ primes[i+1]
            for b in range(BIT_WIDTH):
                if (xv >> b) & 1:
                    gap_bit_flips[gap][b] += 1

    print(f"\nBit flip rates conditioned on gap size:")
    print(f"{'Gap':>4} {'Count':>8} {'Mean HD':>8} {'Bits flipped (0-7)':>40}")
    for gap in sorted(gap_counts.keys()):
        if gap_counts[gap] < 100:
            continue
        count = gap_counts[gap]
        flips = gap_bit_flips[gap]
        mean_hd = sum(flips) / count
        bit_rates = [flips[b] / count for b in range(8)]
        rates_str = " ".join(f"{r:.3f}" for r in bit_rates)
        print(f"{gap:4d} {count:8d} {mean_hd:8.3f} {rates_str:>40}")

    # ================================================================
    # SUMMARY
    # ================================================================
    print("\n" + "=" * 70)
    print("SUMMARY OF FINDINGS")
    print("=" * 70)

    print("""
Key questions answered:
1. XOR structure: Do bit flips between consecutive primes show patterns?
2. Hamming distance: Is the number of changed bits predictable?
3. Per-bit statistics: Are certain bit positions biased?
4. Run lengths: Do primes have unusual binary run-length structure?
5. Transitions: Are bit transitions Markovian or structured?
6. Spectral: Do individual bit sequences have frequency structure?
7. Inter-bit MI: Do bits within a prime carry information about each other?
8. XOR predictability: Can we predict which bits will flip next?
9. Gap-conditioned: Does knowing the gap help predict bit changes?
""")


if __name__ == "__main__":
    main()

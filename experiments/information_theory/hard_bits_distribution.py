"""
Hard Bits Distribution Analysis
================================
The bottom ~50% of bits of p(n) are the "hard" part that R^{-1}(n) gets wrong.
Question: Is the distribution of these hard bits uniform, or are some patterns
more likely? If non-uniform, we could prioritize bruteforce search accordingly.

Analyses:
1. Distribution of delta(n) = p(n) - round(R^{-1}(n)) — is it symmetric? peaked?
2. Per-bit bias in the hard half
3. Joint distribution of hard bit pairs/triples
4. Delta mod small numbers — any residue biases?
5. Conditional distribution: given the EASY bits, what's the distribution of HARD bits?
6. Ordering optimization: if we search by most-likely delta first, how much do we save?
7. Delta sign bias and magnitude distribution
8. Hard-bit patterns conditioned on n mod small numbers

Session 41.
"""

import math
import numpy as np
from collections import Counter, defaultdict

def sieve_primes(limit):
    is_prime = bytearray([1]) * (limit + 1)
    is_prime[0] = is_prime[1] = 0
    for i in range(2, int(limit**0.5) + 1):
        if is_prime[i]:
            is_prime[i*i::i] = bytearray(len(is_prime[i*i::i]))
    return [i for i in range(2, limit + 1) if is_prime[i]]

def R_inverse_approx(n):
    """Approximate R^{-1}(n) using n*ln(n) + corrections."""
    if n < 2:
        return 2
    ln_n = math.log(n)
    ln_ln_n = math.log(ln_n) if ln_n > 1 else 0.1
    # Cipolla-like approximation for R^{-1}
    x = n * (ln_n + ln_ln_n - 1 + (ln_ln_n - 2) / ln_n
             - (ln_ln_n**2 - 6*ln_ln_n + 11) / (2 * ln_n**2))
    return x

def main():
    print("=" * 70)
    print("HARD BITS DISTRIBUTION ANALYSIS")
    print("=" * 70)

    LIMIT = 10_000_000
    primes = sieve_primes(LIMIT)
    N = len(primes)
    print(f"Primes: {N:,}, max = {primes[-1]:,}")

    # Compute delta(n) = p(n) - round(R^{-1}(n))
    SAMPLE = min(N, 500_000)
    # Skip first 100 primes (small-number effects)
    START = 100

    deltas = []
    r_inv_vals = []
    for i in range(START, START + SAMPLE):
        n = i + 1  # 1-indexed
        r_inv = R_inverse_approx(n)
        delta = primes[i] - round(r_inv)
        deltas.append(delta)
        r_inv_vals.append(r_inv)

    delta_arr = np.array(deltas, dtype=np.float64)
    print(f"Sample: {SAMPLE:,} primes (indices {START+1} to {START+SAMPLE})")

    # ================================================================
    # 1. DELTA DISTRIBUTION: Is it symmetric? Gaussian? Heavy-tailed?
    # ================================================================
    print("\n" + "=" * 70)
    print("1. DISTRIBUTION OF delta(n) = p(n) - round(R^{-1}(n))")
    print("=" * 70)

    print(f"\n  Mean:     {delta_arr.mean():+.2f}")
    print(f"  Median:   {np.median(delta_arr):+.1f}")
    print(f"  Std:      {delta_arr.std():.2f}")
    print(f"  Skewness: {float(np.mean(((delta_arr - delta_arr.mean())/delta_arr.std())**3)):.4f}")
    print(f"  Kurtosis: {float(np.mean(((delta_arr - delta_arr.mean())/delta_arr.std())**4)) - 3:.4f}")
    print(f"  Min:      {delta_arr.min():.0f}")
    print(f"  Max:      {delta_arr.max():.0f}")
    print(f"  % positive: {(delta_arr > 0).mean()*100:.2f}%")
    print(f"  % negative: {(delta_arr < 0).mean()*100:.2f}%")
    print(f"  % zero:     {(delta_arr == 0).mean()*100:.2f}%")

    # Percentile distribution
    print(f"\n  Percentile distribution of |delta|:")
    for pct in [25, 50, 75, 90, 95, 99, 99.9]:
        val = np.percentile(np.abs(delta_arr), pct)
        print(f"    {pct:5.1f}%ile: |delta| <= {val:.0f}")

    # ================================================================
    # 2. DELTA MOD SMALL NUMBERS: Any residue biases?
    # ================================================================
    print("\n" + "=" * 70)
    print("2. DELTA MOD SMALL NUMBERS: Residue biases")
    print("=" * 70)

    for mod in [2, 3, 4, 5, 6, 8, 10, 12, 16, 30]:
        residues = [(d % mod) for d in deltas]
        counts = Counter(residues)
        expected = SAMPLE / mod
        chi2 = sum((counts.get(r, 0) - expected)**2 / expected for r in range(mod))
        max_dev = max(abs(counts.get(r, 0)/SAMPLE - 1/mod) for r in range(mod))

        # Show distribution for small mods
        if mod <= 6:
            dist = " ".join(f"{r}:{counts.get(r,0)/SAMPLE:.4f}" for r in range(mod))
            print(f"  mod {mod:2d}: chi2={chi2:10.1f} max_dev={max_dev:.5f}  [{dist}]")
        else:
            print(f"  mod {mod:2d}: chi2={chi2:10.1f} max_dev={max_dev:.5f}")

    # ================================================================
    # 3. PER-BIT BIAS IN THE HARD HALF
    # ================================================================
    print("\n" + "=" * 70)
    print("3. PER-BIT BIAS IN THE HARD HALF OF p(n)")
    print("=" * 70)

    # Determine bit width and the hard/easy boundary
    max_p = max(primes[START:START+SAMPLE])
    BW = max_p.bit_length()
    HARD_BOUNDARY = BW // 2  # bottom half = hard bits

    print(f"  Bit width: {BW}, hard bits: 0-{HARD_BOUNDARY-1}, easy bits: {HARD_BOUNDARY}-{BW-1}")

    print(f"\n  Per-bit P(1) in the HARD half of p(n):")
    print(f"  {'Bit':>4} {'P(1)':>8} {'Dev from 0.5':>14} {'z-score':>10}")
    for b in range(HARD_BOUNDARY):
        ones = sum(1 for i in range(SAMPLE) if (primes[START+i] >> b) & 1)
        p1 = ones / SAMPLE
        dev = p1 - 0.5
        z = dev / (0.5 / math.sqrt(SAMPLE))
        sig = "***" if abs(z) > 3.29 else ""
        print(f"  {b:4d} {p1:8.5f} {dev:+14.5f} {z:+10.2f} {sig}")

    # ================================================================
    # 4. JOINT DISTRIBUTION OF HARD BIT PAIRS
    # ================================================================
    print("\n" + "=" * 70)
    print("4. JOINT DISTRIBUTION: HARD BIT PAIRS")
    print("=" * 70)

    print(f"  Mutual information between hard bit pairs of p(n):")
    print(f"  {'(j,k)':>8} {'MI(bits)':>10} {'P(00)':>8} {'P(01)':>8} {'P(10)':>8} {'P(11)':>8}")
    for j in range(min(HARD_BOUNDARY, 8)):
        for k in range(j+1, min(HARD_BOUNDARY, 8)):
            joint = [[0,0],[0,0]]
            for i in range(SAMPLE):
                bj = (primes[START+i] >> j) & 1
                bk = (primes[START+i] >> k) & 1
                joint[bj][bk] += 1
            mi = 0.0
            for a in range(2):
                for b in range(2):
                    p_ab = joint[a][b] / SAMPLE
                    p_a = (joint[a][0] + joint[a][1]) / SAMPLE
                    p_b = (joint[0][b] + joint[1][b]) / SAMPLE
                    if p_ab > 0 and p_a > 0 and p_b > 0:
                        mi += p_ab * math.log2(p_ab / (p_a * p_b))
            if mi > 0.0001:
                p00 = joint[0][0]/SAMPLE
                p01 = joint[0][1]/SAMPLE
                p10 = joint[1][0]/SAMPLE
                p11 = joint[1][1]/SAMPLE
                print(f"  ({j},{k}){'':<3} {mi:10.6f} {p00:8.4f} {p01:8.4f} {p10:8.4f} {p11:8.4f}")

    # Count how many pairs have MI > threshold
    sig_pairs = 0
    total_pairs = 0
    for j in range(HARD_BOUNDARY):
        for k in range(j+1, HARD_BOUNDARY):
            total_pairs += 1
            joint = [[0,0],[0,0]]
            for idx in range(min(SAMPLE, 50000)):  # subsample for speed
                bj = (primes[START+idx] >> j) & 1
                bk = (primes[START+idx] >> k) & 1
                joint[bj][bk] += 1
            n_sub = min(SAMPLE, 50000)
            mi = 0.0
            for a in range(2):
                for b in range(2):
                    p_ab = joint[a][b] / n_sub
                    p_a = (joint[a][0] + joint[a][1]) / n_sub
                    p_b = (joint[0][b] + joint[1][b]) / n_sub
                    if p_ab > 0 and p_a > 0 and p_b > 0:
                        mi += p_ab * math.log2(p_ab / (p_a * p_b))
            if mi > 0.001:
                sig_pairs += 1
    print(f"\n  Pairs with MI > 0.001: {sig_pairs}/{total_pairs}")

    # ================================================================
    # 5. DELTA DISTRIBUTION CONDITIONED ON n mod 6, n mod 30
    # ================================================================
    print("\n" + "=" * 70)
    print("5. DELTA DISTRIBUTION CONDITIONED ON n mod 6 AND n mod 30")
    print("=" * 70)

    for mod in [6, 30]:
        print(f"\n  Conditioned on n mod {mod}:")
        print(f"  {'residue':>8} {'count':>8} {'mean δ':>10} {'std δ':>10} {'median δ':>10} {'% positive':>12}")
        for r in range(mod):
            subset = [deltas[i] for i in range(SAMPLE) if (START + 1 + i) % mod == r]
            if len(subset) > 100:
                sa = np.array(subset)
                print(f"  {r:8d} {len(subset):8d} {sa.mean():+10.2f} {sa.std():10.2f} {np.median(sa):+10.1f} {(sa > 0).mean()*100:11.2f}%")

    # ================================================================
    # 6. SEARCH ORDERING OPTIMIZATION
    # ================================================================
    print("\n" + "=" * 70)
    print("6. SEARCH ORDERING: How much can we save by searching likely deltas first?")
    print("=" * 70)

    abs_deltas = np.abs(delta_arr)
    sorted_abs = np.sort(abs_deltas)

    # If we search symmetrically outward from R^{-1}(n): 0, ±1, ±2, ...
    # How many candidates until we hit the right one?

    # Simulate: for each prime, how many candidates in {R^{-1} + k} do we check
    # if we search in order of most-likely delta?

    # First: histogram of delta values
    delta_ints = [int(d) for d in deltas]
    delta_counts = Counter(delta_ints)
    total = len(delta_ints)

    # Most common deltas
    print(f"\n  Top 20 most common delta values:")
    for val, cnt in delta_counts.most_common(20):
        print(f"    delta = {val:+6d}: {cnt:6d} ({cnt/total*100:.3f}%)")

    # Search strategy: order candidates by frequency of delta
    # For each observed delta, what rank would it have in frequency-sorted order?
    sorted_deltas = [v for v, c in delta_counts.most_common()]
    delta_rank = {v: i+1 for i, v in enumerate(sorted_deltas)}

    ranks = [delta_rank[d] for d in delta_ints]
    rank_arr = np.array(ranks, dtype=np.float64)

    print(f"\n  If we search deltas in frequency order:")
    print(f"    Mean rank (candidates checked): {rank_arr.mean():.1f}")
    print(f"    Median rank: {np.median(rank_arr):.0f}")
    print(f"    90th percentile rank: {np.percentile(rank_arr, 90):.0f}")
    print(f"    99th percentile rank: {np.percentile(rank_arr, 99):.0f}")
    print(f"    Total distinct deltas: {len(delta_counts)}")

    # Compare with naive search (outward from 0)
    naive_ranks = [abs(d) + 1 for d in delta_ints]  # search 0, ±1, ±2, ...
    naive_arr = np.array(naive_ranks, dtype=np.float64)
    print(f"\n  If we search outward from 0 (naive ±1, ±2, ...):")
    print(f"    Mean candidates checked: {naive_arr.mean():.1f}")
    print(f"    Median: {np.median(naive_arr):.0f}")
    print(f"    90th percentile: {np.percentile(naive_arr, 90):.0f}")

    # Even smarter: search outward but skip non-prime candidates
    # How many of the candidates near R^{-1} are prime?
    # Just measure: what fraction of integers near p(n) are prime?
    prime_density_near = []
    for i in range(0, min(SAMPLE, 100000)):
        p = primes[START + i]
        window = 100
        count_primes = sum(1 for j in range(max(2, p - window), p + window + 1)
                          if j < LIMIT and j > 1 and all(j % d != 0 for d in range(2, min(int(j**0.5)+1, 1000))))
        # Rough density using PNT
        prime_density_near.append(1 / math.log(p))

    avg_density = np.mean(prime_density_near[:10000])
    print(f"\n  Local prime density near p(n): ~{avg_density:.5f} (1/ln(p))")
    print(f"  Expected candidates to check (search only primes): {1/avg_density:.1f}")

    # ================================================================
    # 7. MAGNITUDE DISTRIBUTION OF |delta|
    # ================================================================
    print("\n" + "=" * 70)
    print("7. MAGNITUDE DISTRIBUTION: |delta(n)| scaling")
    print("=" * 70)

    # Bin by n-size and check if |delta| scales predictably
    bin_edges = [100, 500, 1000, 5000, 10000, 50000, 100000, 500000]
    print(f"\n  {'n range':>20} {'mean |δ|':>10} {'std |δ|':>10} {'mean |δ|/√p':>14} {'mean |δ|/p^0.4':>16}")
    for bi in range(len(bin_edges) - 1):
        lo, hi = bin_edges[bi], bin_edges[bi+1]
        subset_d = []
        subset_p = []
        for i in range(SAMPLE):
            n = START + 1 + i
            if lo <= n < hi:
                subset_d.append(abs(deltas[i]))
                subset_p.append(primes[START + i])
        if len(subset_d) > 10:
            sd = np.array(subset_d)
            sp = np.array(subset_p, dtype=np.float64)
            mean_d = sd.mean()
            std_d = sd.std()
            ratio_sqrt = (sd / np.sqrt(sp)).mean()
            ratio_04 = (sd / sp**0.4).mean()
            print(f"  {lo:>8d}-{hi:<8d} {mean_d:10.1f} {std_d:10.1f} {ratio_sqrt:14.6f} {ratio_04:16.6f}")

    # ================================================================
    # 8. CONDITIONAL ON n mod 30: DOES IT HELP?
    # ================================================================
    print("\n" + "=" * 70)
    print("8. CONDITIONAL SEARCH: Does n mod 30 narrow the delta distribution?")
    print("=" * 70)

    for mod in [6, 30]:
        all_std = delta_arr.std()
        cond_stds = []
        for r in range(mod):
            subset = [deltas[i] for i in range(SAMPLE) if (START + 1 + i) % mod == r]
            if len(subset) > 100:
                cond_stds.append(np.std(subset))
        mean_cond_std = np.mean(cond_stds)
        reduction = (1 - mean_cond_std / all_std) * 100
        print(f"  mod {mod:2d}: unconditional std={all_std:.1f}, mean conditional std={mean_cond_std:.1f}, reduction={reduction:.2f}%")

    # ================================================================
    # 9. EVEN/ODD PATTERN IN HARD BITS
    # ================================================================
    print("\n" + "=" * 70)
    print("9. IS delta(n) EVEN OR ODD MORE OFTEN?")
    print("=" * 70)

    even_count = sum(1 for d in deltas if d % 2 == 0)
    odd_count = SAMPLE - even_count
    print(f"  Even deltas: {even_count} ({even_count/SAMPLE*100:.2f}%)")
    print(f"  Odd deltas:  {odd_count} ({odd_count/SAMPLE*100:.2f}%)")

    # Mod 6 of delta (since primes > 3 are 1 or 5 mod 6)
    print(f"\n  delta mod 6 distribution:")
    d6 = Counter(d % 6 for d in deltas)
    for r in range(6):
        print(f"    delta ≡ {r} (mod 6): {d6.get(r,0)/SAMPLE*100:.2f}%")

    # ================================================================
    # SUMMARY
    # ================================================================
    print("\n" + "=" * 70)
    print("SUMMARY: Is there a statistical edge in the hard bits?")
    print("=" * 70)


if __name__ == "__main__":
    main()

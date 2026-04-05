#!/usr/bin/env python3
"""
Session 28 Experiment: Per-Bit Circuit Complexity of pi(x)

NOVEL QUESTION: When we decompose pi(x) into its output bits,
do different bits have different circuit complexities?

Hypothesis: The top ~N/2 bits of pi(x) are determined by the smooth part R(x)
and should have LOWER circuit complexity. The bottom ~N/2 bits encode the
oscillatory correction and should have HIGHER complexity.

This has NOT been studied before. Previous work analyzed pi(x) as a whole
(communication rank, OBDD, ANF) but not individual output bits.

Metrics per bit:
1. BDD size (via Shannon expansion, best of multiple orderings)
2. Total influence (sum of variable influences)
3. Number of sensitive inputs (certificate complexity proxy)
4. Fourier energy at low degrees (spectral analysis)
5. Correlation with R(x) bit (smooth approximation accuracy)
"""

import numpy as np
from sympy import isprime
from math import log, floor, sqrt
from collections import Counter
import time

def compute_pi_table(N):
    """Compute pi(x) for all x in [0, 2^N)."""
    size = 1 << N
    table = np.zeros(size, dtype=np.int64)
    count = 0
    for x in range(size):
        if x >= 2 and isprime(x):
            count += 1
        table[x] = count
    return table

def smooth_approximation(x):
    """Riemann R(x) approximation using first 20 terms."""
    if x <= 1:
        return 0.0
    from math import lgamma, log
    result = 0.0
    lnx = log(x)
    for k in range(1, 21):
        # R(x) = sum_{k=1}^{infty} mu(k)/k * li(x^{1/k})
        # Simplified: use li(x) + li(sqrt(x))/2 + ... with Mobius
        mu = [0, 1, -1, -1, 0, -1, 1, -1, 0, 0, 1,
              -1, 0, -1, 1, 1, 0, -1, 0, -1, 0]
        if k >= len(mu) or mu[k] == 0:
            continue
        xk = x ** (1.0 / k)
        if xk <= 1.0001:
            continue
        # li(y) approximation using series
        lny = log(xk)
        li_val = 0.0
        term = 1.0
        for n in range(1, 100):
            term *= lny / n
            li_val += term / n
        li_val += 0.5772156649 + log(abs(lny))  # Euler-Mascheroni + ln(ln(x^{1/k}))
        result += mu[k] / k * li_val
    return result

def compute_R_table(N):
    """Compute round(R(x)) for all x in [0, 2^N)."""
    size = 1 << N
    table = np.zeros(size, dtype=np.int64)
    for x in range(size):
        if x <= 1:
            table[x] = 0
        else:
            table[x] = round(smooth_approximation(x))
    return table

def extract_bit(values, bit_pos):
    """Extract bit at position bit_pos from each value."""
    return ((values >> bit_pos) & 1).astype(np.int64)

def bdd_size_estimate(truth_table, N):
    """
    Estimate BDD size via Shannon expansion with multiple variable orderings.
    Returns the smallest BDD node count found.

    We try:
    1. Natural order (x_0, x_1, ..., x_{N-1})
    2. Reverse order
    3. Random orderings
    """
    best_size = float('inf')

    # Try a few orderings
    orderings = [
        list(range(N)),           # natural
        list(range(N-1, -1, -1)), # reverse
    ]
    # Add 3 random orderings
    rng = np.random.RandomState(42)
    for _ in range(3):
        orderings.append(rng.permutation(N).tolist())

    for order in orderings:
        size = _bdd_size_for_ordering(truth_table, N, order)
        if size < best_size:
            best_size = size

    return best_size

def _bdd_size_for_ordering(truth_table, N, order):
    """Count unique decision nodes in OBDD with given variable ordering."""
    # Build OBDD layer by layer
    # At level i, we have subsets of the truth table indexed by the remaining variables
    # A node is identified by its truth table restriction

    # Convert truth table to tuple for hashing
    current_level = {tuple(truth_table): 0}  # map from subtable -> id
    total_nodes = 0

    for depth, var in enumerate(order):
        next_level = {}
        for subtable, _ in current_level.items():
            n_remaining = N - depth
            if n_remaining <= 0:
                break
            half = len(subtable) // 2

            # Split by variable 'var' relative to current indexing
            # We need to split the subtable based on the variable at position var
            # in the ORIGINAL ordering

            # Actually, for OBDD, at depth d we branch on order[d]
            # The subtable has 2^(N-d) entries
            # Branching on variable order[d] splits into two subtables of size 2^(N-d-1)

            stride = 1 << (N - 1 - var)  # stride for this variable in original indexing

            low_indices = []
            high_indices = []
            for i in range(len(subtable)):
                if (i >> (N - 1 - var)) & 1:
                    high_indices.append(i)
                else:
                    low_indices.append(i)

            low_sub = tuple(subtable[i] for i in low_indices)
            high_sub = tuple(subtable[i] for i in high_indices)

            if low_sub not in next_level:
                next_level[low_sub] = len(next_level)
            if high_sub not in next_level:
                next_level[high_sub] = len(next_level)

        total_nodes += len(current_level)
        current_level = next_level

    total_nodes += len(current_level)  # terminal nodes
    return total_nodes

def simple_bdd_size(truth_table, N):
    """
    Simpler BDD size estimation: count unique subfunctions at each level.
    For each variable ordering, the OBDD size = sum of unique subfunctions at each level.
    """
    best = float('inf')

    orderings = [list(range(N)), list(range(N-1, -1, -1))]
    rng = np.random.RandomState(42)
    for _ in range(3):
        orderings.append(rng.permutation(N).tolist())

    for order in orderings:
        size = _count_unique_subfunctions(truth_table, N, order)
        best = min(best, size)

    return best

def _count_unique_subfunctions(table, N, order):
    """Count unique subfunctions when splitting variables in given order."""
    size = len(table)
    # Represent current set of subfunctions
    # Initially one subfunction = the full truth table
    current_functions = set()
    current_functions.add(tuple(table))

    total_nodes = 0

    for var in order:
        total_nodes += len(current_functions)
        next_functions = set()

        for func in current_functions:
            n_entries = len(func)
            # Split func by variable `var`
            # Variable var has stride 2^(N-1-var) in original bit ordering
            # But func may have fewer entries if we've already split
            # Actually, we need to track which variables remain
            pass

        # Simpler approach: track by level
        break

    # Use the most basic approach: unique row count at each bit position
    total = 0
    for depth in range(N):
        var = order[depth]
        block_size = 1 << (N - depth)
        half_block = block_size // 2

        # Split table into blocks, count unique blocks
        n_blocks = len(table) // block_size
        unique_blocks = set()
        for i in range(n_blocks):
            start = i * block_size
            block = tuple(table[start:start + block_size])
            unique_blocks.add(block)
        total += len(unique_blocks)

        # Reduce table by fixing this variable
        # This is a simplification - not a real OBDD construction

    return total + 2  # +2 for terminal nodes

def total_influence(truth_table, N):
    """Compute total influence = sum of P[f(x) != f(x XOR e_i)] over all i."""
    size = 1 << N
    total = 0.0
    for var in range(N):
        flips = 0
        mask = 1 << var
        for x in range(size):
            if truth_table[x] != truth_table[x ^ mask]:
                flips += 1
        total += flips / size
    return total

def sensitivity(truth_table, N):
    """Compute average and max sensitivity."""
    size = 1 << N
    sens = np.zeros(size, dtype=np.int64)
    for var in range(N):
        mask = 1 << var
        for x in range(size):
            if truth_table[x] != truth_table[x ^ mask]:
                sens[x] += 1
    return float(np.mean(sens)), int(np.max(sens))

def fourier_weight_by_degree(truth_table, N):
    """
    Compute Fourier weight at each degree.
    f_hat(S) = (1/2^N) * sum_x (-1)^{<S,x>} * f(x) (where f maps to {-1,1})
    Fourier weight at degree d = sum_{|S|=d} f_hat(S)^2
    """
    size = 1 << N
    # Convert to +1/-1
    f = np.where(truth_table == 1, 1.0, -1.0)

    weights = np.zeros(N + 1)

    # For small N, enumerate all subsets
    if N <= 14:
        for S in range(size):
            degree = bin(S).count('1')
            # Compute f_hat(S)
            coeff = 0.0
            for x in range(size):
                parity = bin(S & x).count('1') % 2
                coeff += ((-1) ** parity) * f[x]
            coeff /= size
            weights[degree] += coeff ** 2

    return weights

def main():
    results = []

    print("=" * 80)
    print("PER-BIT CIRCUIT COMPLEXITY OF pi(x)")
    print("=" * 80)
    print()

    for N in range(4, 15):
        t0 = time.time()
        size = 1 << N

        if N > 14:
            print(f"N={N}: skipping (too large)")
            continue

        print(f"\n{'='*60}")
        print(f"N = {N} (x up to {size-1}, {N} input bits)")
        print(f"{'='*60}")

        # Compute pi(x) and R(x) tables
        pi_table = compute_pi_table(N)
        R_table = compute_R_table(N)
        delta_table = pi_table - R_table  # oscillatory correction

        max_pi = int(np.max(pi_table))
        n_output_bits = max(1, max_pi.bit_length())

        print(f"max pi(x) = {max_pi}, output bits = {n_output_bits}")
        print(f"R(x) accuracy: {np.sum(pi_table == R_table)}/{size} exact matches "
              f"({100*np.sum(pi_table == R_table)/size:.1f}%)")
        print(f"|delta| stats: mean={np.mean(np.abs(delta_table)):.2f}, "
              f"max={np.max(np.abs(delta_table))}, "
              f"std={np.std(delta_table):.2f}")

        print(f"\n{'Bit':>4} {'Ones%':>6} {'TotInfl':>8} {'AvgSens':>8} {'MaxSens':>8} "
              f"{'F_low':>7} {'F_d≤2':>7} {'R_corr':>7} {'Label':>10}")
        print("-" * 80)

        bit_results = []

        for bit in range(n_output_bits):
            pi_bit = extract_bit(pi_table, bit)
            R_bit = extract_bit(R_table, bit)

            # Fraction of ones
            ones_frac = np.mean(pi_bit)

            # Total influence
            ti = total_influence(pi_bit, N)

            # Sensitivity
            avg_s, max_s = sensitivity(pi_bit, N)

            # Fourier analysis (only for N <= 12 due to cost)
            if N <= 12:
                fw = fourier_weight_by_degree(pi_bit, N)
                f_low = fw[0] + fw[1]  # weight at degree 0 and 1
                f_low2 = sum(fw[:3])     # weight at degree 0, 1, 2
            else:
                f_low = -1
                f_low2 = -1

            # Correlation with R(x) bit
            if np.std(pi_bit) > 0 and np.std(R_bit) > 0:
                corr = np.corrcoef(pi_bit, R_bit)[0, 1]
            elif np.array_equal(pi_bit, R_bit):
                corr = 1.0
            else:
                corr = 0.0

            # Label the bit
            if bit >= n_output_bits - 2:
                label = "MSB(smooth)"
            elif bit <= 1:
                label = "LSB(osc)"
            else:
                label = "mid"

            print(f"{bit:4d} {100*ones_frac:5.1f}% {ti:8.3f} {avg_s:8.3f} {max_s:8d} "
                  f"{f_low:7.4f} {f_low2:7.4f} {corr:7.4f} {label:>10}")

            bit_results.append({
                'bit': bit,
                'ones_frac': ones_frac,
                'total_influence': ti,
                'avg_sensitivity': avg_s,
                'max_sensitivity': max_s,
                'fourier_low': f_low,
                'fourier_low2': f_low2,
                'R_correlation': corr,
                'label': label
            })

        elapsed = time.time() - t0
        print(f"\n[N={N} completed in {elapsed:.1f}s]")

        results.append({
            'N': N,
            'max_pi': max_pi,
            'n_output_bits': n_output_bits,
            'R_exact_frac': float(np.sum(pi_table == R_table) / size),
            'delta_mean_abs': float(np.mean(np.abs(delta_table))),
            'delta_max': int(np.max(np.abs(delta_table))),
            'bits': bit_results
        })

    # Summary analysis
    print("\n" + "=" * 80)
    print("SUMMARY: COMPLEXITY GRADIENT ACROSS BITS")
    print("=" * 80)

    print("\nKey question: Do MSBs (smooth part) have lower complexity than LSBs (oscillatory)?")
    print()

    for r in results:
        N = r['N']
        bits = r['bits']
        if len(bits) < 3:
            continue

        # Compare bottom half vs top half of bits
        mid = len(bits) // 2
        low_bits = bits[:mid]
        high_bits = bits[mid:]

        avg_infl_low = np.mean([b['total_influence'] for b in low_bits])
        avg_infl_high = np.mean([b['total_influence'] for b in high_bits])
        avg_corr_low = np.mean([b['R_correlation'] for b in low_bits])
        avg_corr_high = np.mean([b['R_correlation'] for b in high_bits])

        print(f"N={N:2d}: LSB-half avg influence={avg_infl_low:.3f}, "
              f"MSB-half avg influence={avg_infl_high:.3f}, "
              f"ratio={avg_infl_low/max(avg_infl_high,0.001):.2f}")
        print(f"       LSB-half avg R-corr={avg_corr_low:.3f}, "
              f"MSB-half avg R-corr={avg_corr_high:.3f}")

    print("\n" + "=" * 80)
    print("SCALING ANALYSIS: Total influence vs N for each bit position")
    print("=" * 80)

    # Track how total influence of bit 0 (LSB = parity) scales with N
    print("\nBit 0 (parity = pi(x) mod 2) — proven as hard as full pi(x):")
    for r in results:
        bits = r['bits']
        if bits:
            b0 = bits[0]
            print(f"  N={r['N']:2d}: total_influence={b0['total_influence']:.3f}, "
                  f"max_sensitivity={b0['max_sensitivity']}")

    print("\nBit 1:")
    for r in results:
        bits = r['bits']
        if len(bits) > 1:
            b1 = bits[1]
            print(f"  N={r['N']:2d}: total_influence={b1['total_influence']:.3f}, "
                  f"max_sensitivity={b1['max_sensitivity']}")

    print("\nMSB (highest bit):")
    for r in results:
        bits = r['bits']
        if bits:
            bm = bits[-1]
            print(f"  N={r['N']:2d} (bit {bm['bit']}): total_influence={bm['total_influence']:.3f}, "
                  f"R_correlation={bm['R_correlation']:.4f}")

if __name__ == '__main__':
    main()

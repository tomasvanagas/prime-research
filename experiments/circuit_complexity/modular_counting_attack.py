#!/usr/bin/env python3
"""
Modular Counting Attack: Wheel Decomposition Circuit Complexity
================================================================

Novel question: Does decomposing pi(x) by residue class mod M (primorial)
give SMALLER circuits for each per-class counting function?

Key distinction from prior CRT work (Sessions 3,7,9,13,20-22,24):
  - We are NOT computing p(n) mod q
  - We ARE asking: is pi_r(x) = #{n <= x : n ≡ r mod M, n prime} a simpler
    function of (x div M) than pi(x) is of x?

Experiments:
1. Mixed-radix conditional entropy of pi(x) increments
2. Per-class circuit complexity (BDD/gate count proxy)
3. Cross-class mutual information
4. Divide-and-conquer scaling analysis

Author: Claude + human
Date: 2026-04-05
"""

import math
import time
import sys
from collections import defaultdict, Counter
from functools import lru_cache

# ─── Utility functions ───

def sieve(limit):
    """Sieve of Eratosthenes, returns boolean array."""
    if limit < 2:
        return [False] * (limit + 1)
    is_prime = [True] * (limit + 1)
    is_prime[0] = is_prime[1] = False
    for i in range(2, int(limit**0.5) + 1):
        if is_prime[i]:
            for j in range(i*i, limit + 1, i):
                is_prime[j] = False
    return is_prime

def pi_array(limit):
    """Compute cumulative prime-counting function."""
    isp = sieve(limit)
    pi = [0] * (limit + 1)
    for i in range(1, limit + 1):
        pi[i] = pi[i-1] + (1 if isp[i] else 0)
    return pi, isp

def euler_phi(n):
    """Euler's totient function."""
    result = n
    p = 2
    temp = n
    while p * p <= temp:
        if temp % p == 0:
            while temp % p == 0:
                temp //= p
            result -= result // p
        p += 1
    if temp > 1:
        result -= result // temp
    return result

def coprime_residues(M):
    """Return sorted list of residues coprime to M."""
    return sorted([r for r in range(1, M + 1) if math.gcd(r, M) == 1])

def entropy(counts):
    """Shannon entropy from a Counter/dict of counts."""
    total = sum(counts.values())
    if total == 0:
        return 0.0
    h = 0.0
    for c in counts.values():
        if c > 0:
            p = c / total
            h -= p * math.log2(p)
    return h

def mutual_information(joint_counts, marginal1_counts, marginal2_counts):
    """Mutual information I(X;Y) = H(X) + H(Y) - H(X,Y)."""
    hx = entropy(marginal1_counts)
    hy = entropy(marginal2_counts)
    hxy = entropy(joint_counts)
    return max(0.0, hx + hy - hxy)

# ─── Primorials ───
PRIMORIALS = {
    6: [2, 3],
    30: [2, 3, 5],
    210: [2, 3, 5, 7],
    2310: [2, 3, 5, 7, 11],
    30030: [2, 3, 5, 7, 11, 13],
}

# ═══════════════════════════════════════════════════════════════════════════════
# EXPERIMENT 1: Mixed-Radix Conditional Entropy
# ═══════════════════════════════════════════════════════════════════════════════

def experiment1_mixed_radix_entropy():
    """
    For each M and each bit-size N:
      - For each coprime residue r mod M, consider x ≡ r (mod M)
      - Compute H(pi(x+M) - pi(x) | x ≡ r mod M, x div M = q)
      - This is the entropy of the pi-increment over one full wheel period,
        conditioned on which period q we're in
      - Compare to unconditional H(pi(x+M) - pi(x))
    """
    print("=" * 72)
    print("EXPERIMENT 1: Mixed-Radix Conditional Entropy of pi(x)")
    print("=" * 72)
    print()

    results = {}

    for N in [10, 12, 14, 16]:
        limit = (1 << N) - 1
        pi_arr, isp = pi_array(limit + 5000)  # extra room for increments

        for M in [6, 30, 210, 2310]:
            if M > limit // 4:
                continue

            coprime_res = coprime_residues(M)
            phi_M = len(coprime_res)

            # Unconditional: distribution of pi(x+M) - pi(x) for x in [M, limit-M]
            uncond_counts = Counter()
            for x in range(M, min(limit - M, len(pi_arr) - M - 1)):
                delta = pi_arr[x + M] - pi_arr[x]
                uncond_counts[delta] += 1
            H_uncond = entropy(uncond_counts)

            # Per-class conditional entropy
            # For each residue r, group by q = x div M, measure entropy of increment
            per_class_entropies = {}
            per_class_counts = {}

            for r in coprime_res:
                # x = q*M + r for q = 0, 1, 2, ...
                # Increment = pi(x + M) - pi(x) = primes in (x, x+M]
                increments_by_q = {}
                overall_counts = Counter()

                max_q = (limit - M) // M
                for q in range(1, max_q + 1):
                    x = q * M + r
                    if x + M >= len(pi_arr):
                        break
                    delta = pi_arr[x + M] - pi_arr[x]
                    increments_by_q[q] = delta
                    overall_counts[delta] += 1

                # Conditional entropy H(delta | q) using binning
                # Group q into blocks of size B and measure entropy within each block
                B = max(16, max_q // 64)
                total_weight = 0
                weighted_H = 0.0
                for block_start in range(1, max_q + 1, B):
                    block_counts = Counter()
                    for q in range(block_start, min(block_start + B, max_q + 1)):
                        if q in increments_by_q:
                            block_counts[increments_by_q[q]] += 1
                    if sum(block_counts.values()) > 0:
                        w = sum(block_counts.values())
                        weighted_H += w * entropy(block_counts)
                        total_weight += w

                H_cond = weighted_H / total_weight if total_weight > 0 else 0
                H_marginal = entropy(overall_counts)
                per_class_entropies[r] = (H_marginal, H_cond)
                per_class_counts[r] = sum(overall_counts.values())

            # Average conditional entropy across classes
            total_samples = sum(per_class_counts.values())
            avg_H_marginal = sum(per_class_counts[r] * per_class_entropies[r][0]
                                  for r in coprime_res) / total_samples if total_samples > 0 else 0
            avg_H_cond = sum(per_class_counts[r] * per_class_entropies[r][1]
                              for r in coprime_res) / total_samples if total_samples > 0 else 0

            results[(N, M)] = {
                'H_uncond': H_uncond,
                'avg_H_marginal': avg_H_marginal,
                'avg_H_cond': avg_H_cond,
                'phi_M': phi_M,
                'num_classes': len(coprime_res),
                'reduction': (H_uncond - avg_H_cond) / H_uncond * 100 if H_uncond > 0 else 0,
            }

            print(f"N={N:2d} bits, M={M:5d} (φ(M)={phi_M:4d}):")
            print(f"  H_uncond(Δπ)        = {H_uncond:.4f} bits")
            print(f"  avg H_marginal(Δπ|r)= {avg_H_marginal:.4f} bits")
            print(f"  avg H_cond(Δπ|r,q)  = {avg_H_cond:.4f} bits")
            print(f"  Entropy reduction   = {results[(N,M)]['reduction']:.1f}%")
            print()

    return results


# ═══════════════════════════════════════════════════════════════════════════════
# EXPERIMENT 2: Per-Class Circuit Complexity (BDD Size Proxy)
# ═══════════════════════════════════════════════════════════════════════════════

def truth_table_complexity(truth_table):
    """
    Proxy for circuit complexity: count the number of "changes" in the
    truth table (transitions 0->1 or 1->0) and measure compressibility.

    Returns dict with multiple complexity measures:
    - transitions: number of bit-changes in truth table
    - run_count: number of runs of same value
    - entropy_rate: entropy of the difference sequence
    - autocorrelation_1: first-order autocorrelation
    """
    if not truth_table:
        return {'transitions': 0, 'run_count': 0, 'entropy_rate': 0, 'autocorr_1': 0}

    n = len(truth_table)
    transitions = sum(1 for i in range(1, n) if truth_table[i] != truth_table[i-1])

    # Runs
    run_count = 1
    for i in range(1, n):
        if truth_table[i] != truth_table[i-1]:
            run_count += 1

    # Entropy rate of difference sequence
    diffs = [truth_table[i] - truth_table[i-1] for i in range(1, n)]
    diff_counts = Counter(diffs)
    ent_rate = entropy(diff_counts)

    # Autocorrelation
    mean = sum(truth_table) / n
    var = sum((t - mean)**2 for t in truth_table) / n
    if var > 0 and n > 1:
        autocorr = sum((truth_table[i] - mean) * (truth_table[i+1] - mean)
                        for i in range(n-1)) / ((n-1) * var)
    else:
        autocorr = 0

    return {
        'transitions': transitions,
        'run_count': run_count,
        'entropy_rate': ent_rate,
        'autocorr_1': autocorr,
        'normalized_transitions': transitions / max(n - 1, 1),
    }


def experiment2_per_class_circuits():
    """
    For each M, build the per-class function pi_r(x) as a function of q = x div M.
    Measure complexity of each class's truth table (pi_r mod 2 as a function of q bits).

    KEY QUESTION: Are per-class functions simpler?
    """
    print("=" * 72)
    print("EXPERIMENT 2: Per-Class Circuit Complexity")
    print("=" * 72)
    print()

    results = {}

    for N in [10, 12, 14, 16]:
        limit = (1 << N) - 1
        pi_arr, isp = pi_array(limit)

        # Full pi(x) mod 2 complexity
        full_tt = [pi_arr[x] % 2 for x in range(2, limit + 1)]
        full_complexity = truth_table_complexity(full_tt)

        for M in [6, 30, 210]:
            if M > limit // 8:
                continue

            coprime_res = coprime_residues(M)

            class_complexities = {}
            for r in coprime_res:
                # pi_r(q) = #{n <= q*M+r : n ≡ r mod M, n prime}
                # Actually, build cumulative count for this class
                tt = []
                count = 0
                max_q = limit // M
                for q in range(0, max_q + 1):
                    x = q * M + r
                    if x <= limit and x >= 2 and isp[x]:
                        count += 1
                    tt.append(count % 2)  # pi_r mod 2

                class_complexities[r] = truth_table_complexity(tt)

            # Average complexity across classes
            avg_transitions = sum(c['normalized_transitions'] for c in class_complexities.values()) / len(class_complexities)
            avg_entropy_rate = sum(c['entropy_rate'] for c in class_complexities.values()) / len(class_complexities)
            avg_autocorr = sum(c['autocorr_1'] for c in class_complexities.values()) / len(class_complexities)

            results[(N, M)] = {
                'full_norm_trans': full_complexity['normalized_transitions'],
                'avg_class_norm_trans': avg_transitions,
                'full_entropy_rate': full_complexity['entropy_rate'],
                'avg_class_entropy_rate': avg_entropy_rate,
                'full_autocorr': full_complexity['autocorr_1'],
                'avg_class_autocorr': avg_autocorr,
                'trans_ratio': avg_transitions / full_complexity['normalized_transitions'] if full_complexity['normalized_transitions'] > 0 else 0,
                'num_classes': len(coprime_res),
            }

            print(f"N={N:2d} bits, M={M:4d} ({len(coprime_res)} classes):")
            print(f"  Full pi(x) mod 2:  transitions={full_complexity['normalized_transitions']:.4f}, "
                  f"H_rate={full_complexity['entropy_rate']:.4f}, autocorr={full_complexity['autocorr_1']:.4f}")
            print(f"  Avg class pi_r mod2: transitions={avg_transitions:.4f}, "
                  f"H_rate={avg_entropy_rate:.4f}, autocorr={avg_autocorr:.4f}")
            print(f"  Transition ratio (class/full): {results[(N,M)]['trans_ratio']:.4f}")

            # Show best and worst classes
            best_r = min(class_complexities, key=lambda r: class_complexities[r]['normalized_transitions'])
            worst_r = max(class_complexities, key=lambda r: class_complexities[r]['normalized_transitions'])
            print(f"  Best class r={best_r}: transitions={class_complexities[best_r]['normalized_transitions']:.4f}")
            print(f"  Worst class r={worst_r}: transitions={class_complexities[worst_r]['normalized_transitions']:.4f}")
            print()

    return results


# ═══════════════════════════════════════════════════════════════════════════════
# EXPERIMENT 3: Cross-Class Mutual Information
# ═══════════════════════════════════════════════════════════════════════════════

def experiment3_cross_class_mi():
    """
    For pairs of coprime residue classes r1, r2 mod M:
    Measure I(pi_{r1}(x) ; pi_{r2}(x)) where x ranges over multiples of M.

    If classes are nearly independent, pi(x) decomposes into independent subproblems.
    """
    print("=" * 72)
    print("EXPERIMENT 3: Cross-Class Mutual Information")
    print("=" * 72)
    print()

    results = {}

    for N in [12, 14, 16]:
        limit = (1 << N) - 1
        pi_arr, isp = pi_array(limit)

        for M in [6, 30, 210]:
            if M > limit // 8:
                continue

            coprime_res = coprime_residues(M)
            max_q = limit // M

            # Build per-class increments: for each q, delta_r = is_prime(q*M + r)?
            class_increments = {}
            for r in coprime_res:
                inc = []
                for q in range(1, max_q + 1):
                    x = q * M + r
                    if x <= limit:
                        inc.append(1 if isp[x] else 0)
                    else:
                        inc.append(0)
                class_increments[r] = inc

            # Pairwise mutual information
            mi_values = []
            num_pairs = 0
            max_pairs = min(len(coprime_res), 20)  # Limit pairs for speed
            selected = coprime_res[:max_pairs]

            for i in range(len(selected)):
                for j in range(i + 1, len(selected)):
                    r1, r2 = selected[i], selected[j]
                    inc1 = class_increments[r1]
                    inc2 = class_increments[r2]

                    # Joint and marginal distributions
                    min_len = min(len(inc1), len(inc2))
                    joint = Counter()
                    m1 = Counter()
                    m2 = Counter()
                    for k in range(min_len):
                        joint[(inc1[k], inc2[k])] += 1
                        m1[inc1[k]] += 1
                        m2[inc2[k]] += 1

                    mi = mutual_information(joint, m1, m2)
                    mi_values.append(mi)
                    num_pairs += 1

            avg_mi = sum(mi_values) / len(mi_values) if mi_values else 0
            max_mi = max(mi_values) if mi_values else 0

            # For context: marginal entropy of a single class
            sample_r = coprime_res[0]
            sample_counts = Counter(class_increments[sample_r])
            H_marginal = entropy(sample_counts)

            results[(N, M)] = {
                'avg_mi': avg_mi,
                'max_mi': max_mi,
                'H_marginal': H_marginal,
                'mi_over_H': avg_mi / H_marginal if H_marginal > 0 else 0,
                'num_pairs': num_pairs,
            }

            print(f"N={N:2d} bits, M={M:4d} ({len(coprime_res)} classes, {num_pairs} pairs tested):")
            print(f"  Marginal H(is_prime in class): {H_marginal:.6f} bits")
            print(f"  Avg mutual info I(r1;r2):      {avg_mi:.6f} bits")
            print(f"  Max mutual info:               {max_mi:.6f} bits")
            print(f"  I/H ratio:                     {results[(N,M)]['mi_over_H']:.6f}")
            independence = "NEAR-INDEPENDENT" if avg_mi < 0.01 else "DEPENDENT"
            print(f"  Assessment: {independence}")
            print()

    return results


# ═══════════════════════════════════════════════════════════════════════════════
# EXPERIMENT 4: Divide-and-Conquer Scaling
# ═══════════════════════════════════════════════════════════════════════════════

def experiment4_scaling():
    """
    Does per-class circuit size grow SLOWER than full circuit size?

    Measure:
    - Full pi(x) mod 2 truth table complexity for N bits
    - Per-class truth table complexity for log2(N/M) input bits
    - Total cost = φ(M) * per-class cost
    - Is this less than full cost?
    """
    print("=" * 72)
    print("EXPERIMENT 4: Divide-and-Conquer Scaling Analysis")
    print("=" * 72)
    print()

    results = {}

    for N in [10, 12, 14, 16]:
        limit = (1 << N) - 1
        pi_arr, isp = pi_array(limit)

        # Full complexity
        full_tt = [pi_arr[x] % 2 for x in range(2, limit + 1)]
        full_comp = truth_table_complexity(full_tt)

        print(f"--- N = {N} bits (x up to {limit}) ---")
        print(f"Full pi(x) mod 2: {full_comp['transitions']} transitions, "
              f"norm={full_comp['normalized_transitions']:.4f}")

        for M in [6, 30, 210, 2310]:
            if M > limit // 8:
                continue

            coprime_res = coprime_residues(M)
            phi_M = len(coprime_res)
            max_q = limit // M
            q_bits = math.ceil(math.log2(max_q + 1)) if max_q > 0 else 1

            total_class_transitions = 0
            total_class_entries = 0

            for r in coprime_res:
                tt = []
                count = 0
                for q in range(0, max_q + 1):
                    x = q * M + r
                    if x <= limit and x >= 2 and isp[x]:
                        count += 1
                    tt.append(count % 2)

                comp = truth_table_complexity(tt)
                total_class_transitions += comp['transitions']
                total_class_entries += len(tt)

            # Compare: total transitions across all classes vs full
            ratio = total_class_transitions / full_comp['transitions'] if full_comp['transitions'] > 0 else 0

            results[(N, M)] = {
                'full_transitions': full_comp['transitions'],
                'total_class_transitions': total_class_transitions,
                'ratio': ratio,
                'phi_M': phi_M,
                'q_bits': q_bits,
                'input_bits_saved': N - q_bits,
            }

            print(f"  M={M:5d}: φ(M)={phi_M:4d}, q_bits={q_bits:2d} (saved {N - q_bits} input bits)")
            print(f"    Total class transitions: {total_class_transitions} vs full: {full_comp['transitions']}")
            print(f"    Ratio (class_total / full): {ratio:.4f}")
            verdict = "WIN" if ratio < 1.0 else "LOSS"
            print(f"    Verdict: {verdict}")

        print()

    return results


# ═══════════════════════════════════════════════════════════════════════════════
# EXPERIMENT 5: Per-Class Entropy as Function of M (Scaling Law)
# ═══════════════════════════════════════════════════════════════════════════════

def experiment5_entropy_scaling():
    """
    How does the per-class conditional entropy scale with M?
    Key theoretical prediction:
    - Prime density in class r mod M is ~ 1/(φ(M) * ln(x/M))
    - Entropy contribution per class should decrease as ~1/φ(M)
    - But the NUMBER of classes is φ(M)
    - So total entropy = φ(M) * H_per_class -- does this stay constant, grow, or shrink?
    """
    print("=" * 72)
    print("EXPERIMENT 5: Per-Class Entropy Scaling with M")
    print("=" * 72)
    print()

    N = 16
    limit = (1 << N) - 1
    pi_arr, isp = pi_array(limit)

    print(f"Fixed N={N} bits (x up to {limit})")
    print()

    for M in [2, 6, 30, 210, 2310]:
        if M > limit // 8:
            continue

        coprime_res = coprime_residues(M)
        phi_M = len(coprime_res)
        max_q = limit // M

        # Per-class: entropy of prime indicator within each class
        class_entropies = []
        class_densities = []

        for r in coprime_res:
            prime_count = 0
            total = 0
            for q in range(1, max_q + 1):
                x = q * M + r
                if x <= limit:
                    total += 1
                    if isp[x]:
                        prime_count += 1

            density = prime_count / total if total > 0 else 0
            class_densities.append(density)

            if density > 0 and density < 1:
                h = -density * math.log2(density) - (1 - density) * math.log2(1 - density)
            else:
                h = 0
            class_entropies.append(h)

        avg_density = sum(class_densities) / len(class_densities) if class_densities else 0
        avg_H = sum(class_entropies) / len(class_entropies) if class_entropies else 0
        total_H = phi_M * avg_H

        # For comparison: overall prime density around limit
        overall_density = pi_arr[limit] / limit

        print(f"M={M:5d} (φ(M)={phi_M:4d}):")
        print(f"  Overall prime density:       {overall_density:.6f}")
        print(f"  Avg per-class prime density:  {avg_density:.6f} (expected ~{1/math.log(limit/M):.6f})")
        print(f"  Avg per-class entropy:        {avg_H:.6f} bits")
        print(f"  Total entropy (φ(M) × avg):   {total_H:.4f} bits")
        print(f"  Density × φ(M):               {avg_density * phi_M:.6f} (should ≈ overall × M/1)")
        print()


# ═══════════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    print("MODULAR COUNTING ATTACK: Wheel Decomposition Circuit Complexity")
    print("=" * 72)
    print()

    t0 = time.time()

    r1 = experiment1_mixed_radix_entropy()
    print()

    r2 = experiment2_per_class_circuits()
    print()

    r3 = experiment3_cross_class_mi()
    print()

    r4 = experiment4_scaling()
    print()

    r5 = experiment5_entropy_scaling()
    print()

    elapsed = time.time() - t0
    print(f"Total runtime: {elapsed:.1f}s")

    # ─── Summary ───
    print()
    print("=" * 72)
    print("SUMMARY")
    print("=" * 72)
    print()

    print("Exp 1 - Entropy reduction from mixed-radix conditioning:")
    for (N, M), v in sorted(r1.items()):
        print(f"  N={N:2d}, M={M:5d}: reduction={v['reduction']:.1f}%")
    print()

    print("Exp 2 - Circuit complexity ratio (per-class / full):")
    for (N, M), v in sorted(r2.items()):
        print(f"  N={N:2d}, M={M:4d}: trans_ratio={v['trans_ratio']:.4f}")
    print()

    print("Exp 3 - Cross-class mutual information:")
    for (N, M), v in sorted(r3.items()):
        print(f"  N={N:2d}, M={M:4d}: I/H={v['mi_over_H']:.6f} ({'INDEP' if v['mi_over_H'] < 0.01 else 'DEP'})")
    print()

    print("Exp 4 - Divide-and-conquer scaling:")
    for (N, M), v in sorted(r4.items()):
        print(f"  N={N:2d}, M={M:5d}: ratio={v['ratio']:.4f} ({'WIN' if v['ratio'] < 1.0 else 'LOSS'})")
    print()

    print("KEY FINDINGS:")
    print("(see experiments/circuit_complexity/modular_counting_attack_results.md)")

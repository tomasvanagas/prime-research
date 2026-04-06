"""
Higher-Dimensional Prime Structure Analysis
=============================================
The Ulam spiral (2D) revealed mod-4 structure. Question: do higher-dimensional
arrangements reveal MORE structure?

Approaches:
1. CRT embedding: map n → (n mod 2, n mod 3, n mod 5, n mod 7, ...) = k-dimensional lattice
   This is the "natural" higher-dimensional structure of integers.
2. Multi-modular MI: does the JOINT distribution of (n mod p1, n mod p2, ...) predict
   primality better than each modulus alone?
3. k-tuple structure: do prime constellations (twin, triplet, etc.) carry info?
4. Higher-dimensional spirals: 3D, 4D, 5D spiral embeddings
5. Fractal / self-similarity: do primes at different scales share structure?
6. The key test: after exhausting ALL modular info, what's left?

Session 41.
"""

import math
import numpy as np
from collections import Counter, defaultdict
from itertools import product

def sieve_primes(limit):
    is_prime = bytearray([1]) * (limit + 1)
    is_prime[0] = is_prime[1] = 0
    for i in range(2, int(limit**0.5) + 1):
        if is_prime[i]:
            is_prime[i*i::i] = bytearray(len(is_prime[i*i::i]))
    return is_prime

def mutual_info_categorical(labels_a, labels_b, n):
    """MI between two categorical sequences."""
    joint = defaultdict(int)
    count_a = defaultdict(int)
    count_b = defaultdict(int)
    for i in range(n):
        joint[(labels_a[i], labels_b[i])] += 1
        count_a[labels_a[i]] += 1
        count_b[labels_b[i]] += 1
    mi = 0.0
    for (a, b), c_ab in joint.items():
        p_ab = c_ab / n
        p_a = count_a[a] / n
        p_b = count_b[b] / n
        if p_ab > 0:
            mi += p_ab * math.log2(p_ab / (p_a * p_b))
    return mi

def entropy(labels, n):
    """Shannon entropy of a categorical sequence."""
    counts = Counter(labels)
    h = 0.0
    for c in counts.values():
        p = c / n
        if p > 0:
            h -= p * math.log2(p)
    return h

def main():
    print("=" * 70)
    print("HIGHER-DIMENSIONAL PRIME STRUCTURE ANALYSIS")
    print("=" * 70)

    LIMIT = 2_000_000
    is_prime = sieve_primes(LIMIT)
    N = LIMIT - 1  # numbers 2..LIMIT

    primes_list = [i for i in range(2, LIMIT + 1) if is_prime[i]]
    N_PRIMES = len(primes_list)
    prime_density = N_PRIMES / N
    h_prime = -prime_density * math.log2(prime_density) - (1-prime_density) * math.log2(1-prime_density)
    print(f"Primes: {N_PRIMES:,} / {N:,} = {prime_density:.6f}")
    print(f"H(is_prime) = {h_prime:.6f} bits")

    # ================================================================
    # 1. CRT EMBEDDING: Cumulative information from modular dimensions
    # ================================================================
    print("\n" + "=" * 70)
    print("1. CRT EMBEDDING: Adding dimensions one prime at a time")
    print("=" * 70)

    small_primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43]
    SAMPLE = min(N, 500_000)
    numbers = list(range(2, SAMPLE + 2))
    is_p_seq = [1 if is_prime[n] else 0 for n in numbers]

    print(f"\nCumulative MI as we add modular dimensions:")
    print(f"{'Dims':>5} {'Moduli':>30} {'Product':>12} {'MI (bits)':>10} {'% of H':>8} {'Marginal':>10}")

    prev_mi = 0
    for k in range(1, len(small_primes) + 1):
        mods = small_primes[:k]
        modulus = 1
        for m in mods:
            modulus *= m

        # Create joint label: tuple of residues
        labels = [tuple(n % m for m in mods) for n in numbers]
        mi = mutual_info_categorical(labels, is_p_seq, SAMPLE)
        marginal = mi - prev_mi
        pct = mi / h_prime * 100

        mods_str = ",".join(str(m) for m in mods)
        print(f"{k:5d} {mods_str:>30} {modulus:>12,} {mi:10.6f} {pct:7.2f}% {marginal:10.6f}")
        prev_mi = mi

    # ================================================================
    # 2. INTERACTION INFORMATION: Do moduli have synergies?
    # ================================================================
    print("\n" + "=" * 70)
    print("2. INTERACTION INFORMATION: Do pairs of moduli have synergies?")
    print("=" * 70)

    print(f"\nII(mod_p1; mod_p2; is_prime) = MI(joint) - MI(p1) - MI(p2)")
    print(f"Positive = synergy (together know more), Negative = redundancy")
    print(f"\n{'(p1,p2)':>10} {'MI(p1)':>8} {'MI(p2)':>8} {'MI(joint)':>10} {'II':>10} {'Type':>10}")

    for i in range(min(8, len(small_primes))):
        for j in range(i+1, min(8, len(small_primes))):
            p1, p2 = small_primes[i], small_primes[j]

            lab1 = [n % p1 for n in numbers]
            lab2 = [n % p2 for n in numbers]
            lab_joint = [(n % p1, n % p2) for n in numbers]

            mi1 = mutual_info_categorical(lab1, is_p_seq, SAMPLE)
            mi2 = mutual_info_categorical(lab2, is_p_seq, SAMPLE)
            mi_joint = mutual_info_categorical(lab_joint, is_p_seq, SAMPLE)

            ii = mi_joint - mi1 - mi2
            typ = "SYNERGY" if ii > 0.0001 else "REDUNDANT" if ii < -0.0001 else "ADDITIVE"
            print(f"({p1:2d},{p2:2d}){'':<3} {mi1:8.5f} {mi2:8.5f} {mi_joint:10.5f} {ii:+10.5f} {typ:>10}")

    # ================================================================
    # 3. TRIPLE INTERACTIONS: Do 3-way combinations reveal more?
    # ================================================================
    print("\n" + "=" * 70)
    print("3. TRIPLE INTERACTIONS: 3-modulus synergies")
    print("=" * 70)

    triples = [(2,3,5), (2,3,7), (2,5,7), (3,5,7), (2,3,11), (5,7,11), (2,3,5,7)]
    print(f"{'Moduli':>15} {'MI(joint)':>10} {'Sum MI(each)':>14} {'Synergy':>10}")
    for mods in triples:
        lab_joint = [tuple(n % m for m in mods) for n in numbers]
        mi_joint = mutual_info_categorical(lab_joint, is_p_seq, SAMPLE)

        mi_sum = sum(
            mutual_info_categorical([n % m for n in numbers], is_p_seq, SAMPLE)
            for m in mods
        )
        synergy = mi_joint - mi_sum
        mods_str = ",".join(str(m) for m in mods)
        print(f"{mods_str:>15} {mi_joint:10.5f} {mi_sum:14.5f} {synergy:+10.5f}")

    # ================================================================
    # 4. THE LIMIT: How much total modular information exists?
    # ================================================================
    print("\n" + "=" * 70)
    print("4. THE LIMIT: Total modular information vs total entropy")
    print("=" * 70)

    # Use progressively larger primorial moduli
    primorials = []
    prod = 1
    for p in small_primes:
        prod *= p
        primorials.append((prod, p))

    print(f"\n{'Primorial':>15} {'Up to p=':>8} {'MI':>10} {'% of H':>8} {'Remaining H':>12}")
    for modulus, max_p in primorials:
        if modulus > SAMPLE:
            # Can't measure MI reliably when modulus > sample size
            # Estimate using Euler product
            remaining_fraction = 1.0
            for p in small_primes:
                if p <= max_p:
                    remaining_fraction *= (1 - 1/p)
            # MI from sieving = H(prime) - H(prime | coprime to primorial)
            # After sieving, density among survivors = prime_density / remaining_fraction
            # if survivor is coprime to primorial
            est_mi = h_prime  # placeholder
            print(f"{modulus:>15,} {max_p:>8} {'(too large to measure)':>30}")
            continue

        labels = [n % modulus for n in numbers]
        mi = mutual_info_categorical(labels, is_p_seq, SAMPLE)
        pct = mi / h_prime * 100
        remaining = h_prime - mi
        print(f"{modulus:>15,} {max_p:>8} {mi:10.6f} {pct:7.2f}% {remaining:12.6f}")

    # Theoretical limit from Euler product
    print(f"\n  Theoretical analysis:")
    print(f"  Sieving by primes up to y removes all composites with smallest factor ≤ y.")
    print(f"  After sieving by all primes up to √x, only primes remain (+ 1).")
    print(f"  The Legendre sieve inclusion-exclusion IS the modular structure.")

    # ================================================================
    # 5. RESIDUAL AFTER FULL MODULAR SIEVE: What structure remains?
    # ================================================================
    print("\n" + "=" * 70)
    print("5. RESIDUAL: Among mod-30030 survivors, any structure?")
    print("=" * 70)

    # Numbers coprime to 30030 = 2*3*5*7*11*13
    def coprime_to_30030(n):
        return all(n % p != 0 for p in [2, 3, 5, 7, 11, 13])

    survivors = [(n, is_prime[n]) for n in numbers if coprime_to_30030(n)]
    n_surv = len(survivors)
    surv_primes = sum(1 for _, ip in survivors if ip)
    surv_density = surv_primes / n_surv

    print(f"  Survivors coprime to 30030: {n_surv:,} ({n_surv/SAMPLE*100:.1f}%)")
    print(f"  Primes among survivors: {surv_primes:,} ({surv_density*100:.2f}%)")
    print(f"  (Expected: {1/math.log(SAMPLE/2)*100:.2f}%)")

    # Among survivors, check higher moduli for structure
    surv_numbers = [n for n, _ in survivors]
    surv_is_prime = [ip for _, ip in survivors]

    for mod in [17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61]:
        lab = [n % mod for n in surv_numbers]
        mi = mutual_info_categorical(lab, surv_is_prime, n_surv)
        h_surv = entropy(surv_is_prime, n_surv)
        pct = mi / h_surv * 100 if h_surv > 0 else 0

        # Expected MI from sieving out mod-p composites
        # Among survivors, those divisible by p are composite
        # Fraction divisible by p: ~1/p (since coprime to 30030, and p > 13)
        expected_mi_approx = (1/mod) * math.log2(mod)  # rough upper bound

        print(f"  mod {mod:2d}: MI = {mi:.6f} bits ({pct:.3f}% of survivor H={h_surv:.4f}), expected ~{1/mod:.4f} fraction sieved")

    # ================================================================
    # 6. DOES COMBINING DIMENSIONS EVER BEAT THEIR SUM?
    # ================================================================
    print("\n" + "=" * 70)
    print("6. HIGHER-DIMENSIONAL LATTICE: Genuine non-linear structure?")
    print("=" * 70)

    # The CRT maps n → (n mod 2, n mod 3, ..., n mod p_k) into a k-dimensional lattice.
    # Question: is there structure in the GEOMETRY of this lattice (distances, clusters)
    # that carries information beyond the individual residues?

    # Test: among numbers with SAME residue class mod 30, does their position
    # in the mod-7 × mod-11 × mod-13 lattice predict primality?
    print(f"\n  Among numbers ≡ 1 (mod 30), checking 3D lattice (mod 7 × 11 × 13):")

    subset_1mod30 = [(n, is_prime[n]) for n in numbers if n % 30 == 1]
    n_sub = len(subset_1mod30)

    # 3D label
    labels_3d = [(n % 7, n % 11, n % 13) for n, _ in subset_1mod30]
    is_p_sub = [ip for _, ip in subset_1mod30]

    mi_3d = mutual_info_categorical(labels_3d, is_p_sub, n_sub)

    # Compare with sum of individual MIs
    mi_7 = mutual_info_categorical([n % 7 for n, _ in subset_1mod30], is_p_sub, n_sub)
    mi_11 = mutual_info_categorical([n % 11 for n, _ in subset_1mod30], is_p_sub, n_sub)
    mi_13 = mutual_info_categorical([n % 13 for n, _ in subset_1mod30], is_p_sub, n_sub)
    mi_sum = mi_7 + mi_11 + mi_13

    print(f"  MI(7×11×13 jointly) = {mi_3d:.6f}")
    print(f"  MI(7) + MI(11) + MI(13) = {mi_sum:.6f}")
    print(f"  Synergy = {mi_3d - mi_sum:+.6f}")
    print(f"  (Positive = geometric structure helps, zero = purely additive)")

    # 5D test
    labels_5d = [(n % 7, n % 11, n % 13, n % 17, n % 19) for n, _ in subset_1mod30]
    mi_5d = mutual_info_categorical(labels_5d, is_p_sub, n_sub)
    mi_sum5 = mi_7 + mi_11 + mi_13
    mi_17 = mutual_info_categorical([n % 17 for n, _ in subset_1mod30], is_p_sub, n_sub)
    mi_19 = mutual_info_categorical([n % 19 for n, _ in subset_1mod30], is_p_sub, n_sub)
    mi_sum5 += mi_17 + mi_19

    print(f"\n  MI(7×11×13×17×19 jointly) = {mi_5d:.6f}")
    print(f"  Sum of individual MIs = {mi_sum5:.6f}")
    print(f"  Synergy = {mi_5d - mi_sum5:+.6f}")

    # ================================================================
    # 7. EUCLIDEAN DISTANCE IN CRT SPACE: Do primes cluster?
    # ================================================================
    print("\n" + "=" * 70)
    print("7. DO PRIMES CLUSTER IN CRT SPACE?")
    print("=" * 70)

    # Map primes and composites to CRT coordinates and measure clustering
    mods_for_dist = [7, 11, 13, 17, 19]
    dim = len(mods_for_dist)

    # Sample primes and composites among mod-30 coprime
    coprime30 = [n for n in range(31, min(SAMPLE, 200000)) if n % 30 in [1,7,11,13,17,19,23,29]]
    prime_coords = []
    comp_coords = []
    for n in coprime30:
        coord = tuple(n % m for m in mods_for_dist)
        if is_prime[n]:
            prime_coords.append(coord)
        else:
            comp_coords.append(coord)

    # Count occupancy of each lattice cell
    prime_cells = Counter(prime_coords)
    comp_cells = Counter(comp_coords)
    all_cells = set(prime_cells.keys()) | set(comp_cells.keys())

    # For each cell: prime fraction
    cell_fracs = []
    cell_sizes = []
    for cell in all_cells:
        np_ = prime_cells.get(cell, 0)
        nc_ = comp_cells.get(cell, 0)
        total = np_ + nc_
        if total > 20:  # enough data
            cell_fracs.append(np_ / total)
            cell_sizes.append(total)

    frac_arr = np.array(cell_fracs)
    overall_frac = len(prime_coords) / (len(prime_coords) + len(comp_coords))

    print(f"  CRT lattice: {dim}D with moduli {mods_for_dist}")
    print(f"  Total cells occupied: {len(all_cells)}")
    print(f"  Cells with >20 entries: {len(cell_fracs)}")
    print(f"  Overall prime fraction: {overall_frac:.4f}")
    print(f"  Per-cell prime fraction: mean={frac_arr.mean():.4f}, std={frac_arr.std():.4f}")
    print(f"  Min cell fraction: {frac_arr.min():.4f}")
    print(f"  Max cell fraction: {frac_arr.max():.4f}")
    print(f"  Coefficient of variation: {frac_arr.std()/frac_arr.mean():.4f}")

    # How many cells have fraction = 0 (all composite)?
    zero_cells = sum(1 for f in cell_fracs if f == 0)
    print(f"  Cells with 0% primes: {zero_cells}/{len(cell_fracs)}")
    print(f"  (These are residue classes divisible by one of {mods_for_dist})")

    # Among cells where no modulus divides, is there variation?
    nondiv_fracs = []
    for cell in all_cells:
        # Check if any coordinate is 0 (meaning n ≡ 0 mod that prime)
        if 0 not in cell:
            np_ = prime_cells.get(cell, 0)
            nc_ = comp_cells.get(cell, 0)
            total = np_ + nc_
            if total > 20:
                nondiv_fracs.append(np_ / total)

    if nondiv_fracs:
        nd_arr = np.array(nondiv_fracs)
        expected_std = math.sqrt(overall_frac * (1-overall_frac) / np.mean(cell_sizes))
        print(f"\n  Among cells not divisible by any of {mods_for_dist}:")
        print(f"    Cells: {len(nondiv_fracs)}")
        print(f"    Mean fraction: {nd_arr.mean():.4f}")
        print(f"    Std fraction:  {nd_arr.std():.4f}")
        print(f"    Expected std (binomial noise): {expected_std:.4f}")
        print(f"    Ratio (actual/expected std): {nd_arr.std()/expected_std:.3f}")
        print(f"    (Ratio ≈ 1.0 means NO structure beyond noise)")

    # ================================================================
    # 8. THE ULTIMATE TEST: Dimensionality vs information scaling
    # ================================================================
    print("\n" + "=" * 70)
    print("8. SCALING: How does information grow with dimensions?")
    print("=" * 70)

    print(f"\n{'k (dims)':>8} {'Moduli used':>30} {'MI (bits)':>10} {'MI/dim':>8} {'Cumulative %H':>14}")

    cumulative_primes_list = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]
    for k in range(1, len(cumulative_primes_list) + 1):
        mods = cumulative_primes_list[:k]
        modulus = 1
        for m in mods:
            modulus *= m

        if modulus > SAMPLE * 2:
            break

        labels = [tuple(n % m for m in mods) for n in numbers]
        mi = mutual_info_categorical(labels, is_p_seq, SAMPLE)
        mi_per_dim = mi / k
        pct = mi / h_prime * 100

        mods_str = ",".join(str(m) for m in mods)
        print(f"{k:8d} {mods_str:>30} {mi:10.6f} {mi_per_dim:8.5f} {pct:13.2f}%")

    # Extrapolate: what's the theoretical limit?
    print(f"\n  Theoretical limit:")
    print(f"  After sieving by ALL primes up to √x, only primes + {1} remain.")
    print(f"  So in the limit of infinite dimensions: MI → H(is_prime) = {h_prime:.6f} bits")
    print(f"  But this IS the sieve of Eratosthenes — checking all primes up to √x.")
    print(f"  Each new dimension (prime p) adds ~log2(p/(p-1)) bits ≈ 1/(p·ln2) bits.")

    # ================================================================
    # SUMMARY
    # ================================================================
    print("\n" + "=" * 70)
    print("SUMMARY: Does higher dimensionality help?")
    print("=" * 70)

if __name__ == "__main__":
    main()

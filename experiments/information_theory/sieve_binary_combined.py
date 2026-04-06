"""
Sieve + Binary Combined Analysis
==================================
Can we combine:
  - R^{-1}(n) giving top ~50% of bits (smooth part)
  - Wheel mod 30030 constraining bottom ~15 bits to 5760 patterns
  - Binary structure of primes

to squeeze MORE determined bits than either alone?

Key idea: If R^{-1} fixes the top half and the wheel constrains the bottom ~15 bits,
maybe only a narrow "gap" of truly unknown bits remains in the MIDDLE.

Analyses:
1. How many bottom bits does mod 30030 actually constrain?
2. Combined: R^{-1} top bits + wheel bottom bits → how many unknown bits remain?
3. Do the wheel residues have non-uniform binary patterns?
4. Cross-information: does knowing the wheel class help predict middle bits?
5. Practical search space: R^{-1} bracket × wheel filter → how many candidates?
6. Can we determine the wheel class from R^{-1}(n)?
7. Bit-level entropy after ALL known constraints applied

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
    if n < 2: return 2
    ln_n = math.log(n)
    ln_ln_n = math.log(ln_n) if ln_n > 1 else 0.1
    return n * (ln_n + ln_ln_n - 1 + (ln_ln_n - 2) / ln_n
                - (ln_ln_n**2 - 6*ln_ln_n + 11) / (2 * ln_n**2))

def main():
    print("=" * 70)
    print("SIEVE + BINARY COMBINED ANALYSIS")
    print("=" * 70)

    LIMIT = 10_000_000
    primes = sieve_primes(LIMIT)
    N = len(primes)
    print(f"Primes: {N:,}, max = {primes[-1]:,}")

    BW = primes[-1].bit_length()
    HALF = BW // 2
    print(f"Bit width: {BW}, half = {HALF}")

    # Wheel mod 30030
    WHEEL = 30030
    coprime_residues = [r for r in range(WHEEL) if all(r % p != 0 for p in [2,3,5,7,11,13])]
    N_RESIDUES = len(coprime_residues)
    print(f"Wheel mod {WHEEL}: {N_RESIDUES} coprime residues out of {WHEEL}")

    # ================================================================
    # 1. HOW MANY BOTTOM BITS DOES MOD 30030 CONSTRAIN?
    # ================================================================
    print("\n" + "=" * 70)
    print("1. BITS CONSTRAINED BY WHEEL MOD 30030")
    print("=" * 70)

    # 30030 ≈ 2^14.87, so the wheel constrains the bottom ~15 bits
    wheel_bits = math.log2(WHEEL)
    residue_bits = math.log2(N_RESIDUES)
    constrained_bits = wheel_bits - residue_bits

    print(f"  log2(30030) = {wheel_bits:.2f} bits")
    print(f"  log2(5760 residues) = {residue_bits:.2f} bits")
    print(f"  Bits constrained by wheel: {constrained_bits:.2f}")
    print(f"  (Out of ~15 bottom bits, wheel eliminates ~{constrained_bits:.1f})")

    # Per-bit entropy of the bottom 15 bits AMONG coprime residues
    print(f"\n  Per-bit entropy of coprime residues (0=fully constrained, 1=random):")
    for b in range(16):
        ones = sum(1 for r in coprime_residues if (r >> b) & 1)
        p1 = ones / N_RESIDUES
        if p1 > 0 and p1 < 1:
            h = -p1 * math.log2(p1) - (1-p1) * math.log2(1-p1)
        elif p1 == 0 or p1 == 1:
            h = 0.0
        else:
            h = 1.0
        print(f"    bit {b:2d}: P(1) = {p1:.4f}, H = {h:.4f} {'← CONSTRAINED' if h < 0.9 else ''}")

    # ================================================================
    # 2. COMBINED: R^{-1} + WHEEL → HOW MANY UNKNOWN BITS?
    # ================================================================
    print("\n" + "=" * 70)
    print("2. COMBINED CONSTRAINTS: R^{-1}(n) + wheel mod 30030")
    print("=" * 70)

    SAMPLE = min(N, 200_000)
    START = 1000  # skip small primes

    # For each prime, compute:
    # - R^{-1} estimate → which top bits are correct?
    # - Wheel residue → which bottom-bit patterns are possible?
    # - How many candidates remain?

    total_r_inv_correct_bits = []
    total_combined_candidates = []
    total_wheel_candidates = []
    total_naive_window = []

    for i in range(START, START + SAMPLE):
        p = primes[i]
        n = i + 1
        r_est = R_inverse_approx(n)
        r_int = round(r_est)

        # How many top bits does R^{-1} get right?
        xor = p ^ r_int
        if xor == 0:
            correct_top = BW
        else:
            correct_top = BW - xor.bit_length()
        total_r_inv_correct_bits.append(correct_top)

        # R^{-1} gives a window: [r_int - error, r_int + error]
        error = abs(p - r_int)

        # Within that window, how many are coprime to 30030?
        # Window size = 2 * error
        window = 2 * error + 1
        # Expected coprime fraction = phi(30030)/30030 = 5760/30030
        wheel_fraction = N_RESIDUES / WHEEL
        expected_coprime = window * wheel_fraction
        total_wheel_candidates.append(expected_coprime)
        total_naive_window.append(window)

        # Among coprime candidates, how many are actually prime?
        # Prime density near p: 1/ln(p)
        prime_density = 1 / math.log(p)
        expected_primes = expected_coprime * prime_density / wheel_fraction
        total_combined_candidates.append(expected_primes)

    r_bits_arr = np.array(total_r_inv_correct_bits, dtype=np.float64)
    wheel_cand_arr = np.array(total_wheel_candidates, dtype=np.float64)
    naive_arr = np.array(total_naive_window, dtype=np.float64)
    combined_arr = np.array(total_combined_candidates, dtype=np.float64)

    print(f"\n  R^{{-1}} correct top bits: mean={r_bits_arr.mean():.1f}/{BW}, median={np.median(r_bits_arr):.0f}")
    print(f"  R^{{-1}} search window: mean={naive_arr.mean():.0f} integers")
    print(f"  After wheel filter: mean={wheel_cand_arr.mean():.0f} candidates ({N_RESIDUES/WHEEL*100:.1f}% of window)")
    print(f"  Expected primes in window: mean={combined_arr.mean():.1f}")
    print(f"\n  Bit accounting for {BW}-bit primes:")
    print(f"    Top bits from R^{{-1}}:       {r_bits_arr.mean():.1f} bits determined")
    print(f"    Bottom bits from wheel:      {constrained_bits:.1f} bits constrained")
    print(f"    Overlap (bottom bits R^{{-1}} also gets): ~0 (R^{{-1}} fails on bottom)")
    print(f"    TOTAL determined:            ~{r_bits_arr.mean() + constrained_bits:.1f} of {BW}")
    print(f"    REMAINING unknown:           ~{BW - r_bits_arr.mean() - constrained_bits:.1f} bits")

    # ================================================================
    # 3. DO WHEEL RESIDUES HAVE EXPLOITABLE BINARY PATTERNS?
    # ================================================================
    print("\n" + "=" * 70)
    print("3. BINARY PATTERNS OF WHEEL RESIDUES")
    print("=" * 70)

    # Among the 5760 coprime residues mod 30030, what binary patterns appear?
    # Group by bottom 4 bits
    bottom4 = Counter(r & 0xF for r in coprime_residues)
    print(f"\n  Bottom 4 bits of coprime residues mod 30030:")
    print(f"  {'Pattern':>8} {'Count':>8} {'Fraction':>10}")
    for pattern in sorted(bottom4.keys()):
        c = bottom4[pattern]
        print(f"  {pattern:4d} ({bin(pattern):>6s}) {c:8d} {c/N_RESIDUES:10.4f}")

    # Bottom 3 bits (since all primes > 2 are odd, bit 0 = 1)
    bottom3 = Counter((r >> 1) & 0x7 for r in coprime_residues)
    print(f"\n  Bits 1-3 of coprime residues:")
    for pattern in sorted(bottom3.keys()):
        c = bottom3[pattern]
        print(f"  bits[3:1]={bin(pattern):>5s}: {c:6d} ({c/N_RESIDUES*100:.1f}%)")

    # ================================================================
    # 4. CAN WE DETERMINE WHEEL CLASS FROM R^{-1}?
    # ================================================================
    print("\n" + "=" * 70)
    print("4. CAN R^{-1}(n) PREDICT THE WHEEL RESIDUE CLASS?")
    print("=" * 70)

    # R^{-1}(n) mod 30030 vs actual p(n) mod 30030
    r_inv_residues = []
    actual_residues = []
    match_count = 0

    for i in range(START, START + SAMPLE):
        p = primes[i]
        n = i + 1
        r_est = round(R_inverse_approx(n))
        r_res = r_est % WHEEL
        p_res = p % WHEEL
        r_inv_residues.append(r_res)
        actual_residues.append(p_res)
        if r_res == p_res:
            match_count += 1

    match_rate = match_count / SAMPLE
    print(f"  R^{{-1}}(n) mod 30030 == p(n) mod 30030: {match_count}/{SAMPLE} ({match_rate*100:.2f}%)")
    print(f"  Random baseline: {1/N_RESIDUES*100:.4f}%")

    # How close? Distance between R^{-1} residue and actual residue
    residue_diffs = [(actual_residues[i] - r_inv_residues[i]) % WHEEL for i in range(SAMPLE)]
    # Map to [-WHEEL/2, WHEEL/2]
    residue_diffs = [d if d <= WHEEL//2 else d - WHEEL for d in residue_diffs]
    rd_arr = np.array(residue_diffs, dtype=np.float64)
    print(f"  Residue difference (p - R^{{-1}}) mod 30030:")
    print(f"    Mean: {rd_arr.mean():+.1f}")
    print(f"    Std:  {rd_arr.std():.1f}")
    print(f"    (For reference, |delta| mean ≈ {naive_arr.mean()/2:.0f})")

    # ================================================================
    # 5. PRACTICAL SEARCH: Combined constraints
    # ================================================================
    print("\n" + "=" * 70)
    print("5. PRACTICAL COMBINED SEARCH")
    print("=" * 70)

    # Strategy: Use R^{-1} for window, wheel to filter, primality test to confirm
    # How many primality tests needed?

    print(f"\n  Search strategy comparison:")
    print(f"  {'Method':>40} {'Candidates':>12} {'Reduction':>10}")
    print(f"  {'R^{-1} window (all integers)':>40} {naive_arr.mean():12.0f} {'1x':>10}")
    print(f"  {'R^{-1} + wheel 30030':>40} {wheel_cand_arr.mean():12.0f} {naive_arr.mean()/wheel_cand_arr.mean():9.1f}x")
    print(f"  {'R^{-1} + primality test only':>40} {combined_arr.mean():12.1f} {naive_arr.mean()/combined_arr.mean():9.1f}x")

    # What if we use smaller wheels too?
    for wheel_mod, phi in [(6, 2), (30, 8), (210, 48), (2310, 480), (30030, 5760)]:
        frac = phi / wheel_mod
        candidates = naive_arr.mean() * frac
        print(f"  {'R^{-1} + wheel ' + str(wheel_mod):>40} {candidates:12.0f} {naive_arr.mean()/candidates:9.1f}x")

    # ================================================================
    # 6. THE MIDDLE BITS: After R^{-1} top and wheel bottom, what's left?
    # ================================================================
    print("\n" + "=" * 70)
    print("6. THE MIDDLE BITS: Entropy of undetermined bits")
    print("=" * 70)

    # For each prime, determine:
    # - Which bits R^{-1} gets right (from MSB down)
    # - Which bits wheel constrains (from LSB up)
    # - The "gap" of unknown middle bits

    middle_entropies = []
    for i in range(START, min(START + 50000, N)):
        p = primes[i]
        n = i + 1
        r_est = round(R_inverse_approx(n))

        # Top bits correct from R^{-1}
        xor = p ^ r_est
        if xor == 0:
            top_correct = BW
        else:
            top_correct = BW - xor.bit_length()

        # Bottom bits constrained by wheel: ~2.4 bits
        bottom_constrained = constrained_bits

        # Middle unknown
        middle = BW - top_correct - bottom_constrained
        middle_entropies.append(max(0, middle))

    mid_arr = np.array(middle_entropies)
    print(f"  Total bits: {BW}")
    print(f"  Top bits (R^{{-1}}): mean {r_bits_arr[:50000].mean():.1f}")
    print(f"  Bottom bits (wheel): {constrained_bits:.1f}")
    print(f"  Middle unknown bits: mean {mid_arr.mean():.1f}")
    print(f"  Middle unknown bits: median {np.median(mid_arr):.0f}")
    print(f"  Candidates = 2^(middle): mean {2**mid_arr.mean():.0f}")

    # ================================================================
    # 7. PER-BIT ANALYSIS: Which bits are determined by EITHER source?
    # ================================================================
    print("\n" + "=" * 70)
    print("7. PER-BIT: Which bits are determined by R^{-1} vs wheel vs neither?")
    print("=" * 70)

    # For each bit position, measure:
    # (a) P(bit correct from R^{-1})
    # (b) Entropy reduction from knowing wheel class
    # (c) Combined

    SUB = min(50000, SAMPLE)
    print(f"\n  {'Bit':>4} {'R^-1 correct':>13} {'Wheel H red':>12} {'Combined P(correct)':>20} {'Status':>10}")

    for b in range(BW):
        # (a) R^{-1} correct for this bit
        r_correct = 0
        for i in range(START, START + SUB):
            p = primes[i]
            r_est = round(R_inverse_approx(i + 1))
            if (p >> b) & 1 == (r_est >> b) & 1:
                r_correct += 1
        r_acc = r_correct / SUB

        # (b) Wheel entropy reduction for this bit
        # Among 5760 coprime residues, what's P(bit_b = 1)?
        wheel_ones = sum(1 for r in coprime_residues if (r >> b) & 1)
        wheel_p1 = wheel_ones / N_RESIDUES
        if wheel_p1 > 0 and wheel_p1 < 1:
            wheel_h = -wheel_p1 * math.log2(wheel_p1) - (1-wheel_p1) * math.log2(1-wheel_p1)
        else:
            wheel_h = 0.0
        wheel_reduction = 1.0 - wheel_h

        # Status
        if r_acc > 0.95:
            status = "R^-1"
        elif wheel_reduction > 0.1:
            status = "WHEEL"
        elif r_acc > 0.6:
            status = "partial"
        else:
            status = "UNKNOWN"

        combined = max(r_acc, 0.5 + wheel_reduction/2)

        print(f"  {b:4d} {r_acc:13.4f} {wheel_reduction:12.4f} {combined:20.4f} {status:>10}")

    # ================================================================
    # 8. SEARCH SPACE AFTER ALL CONSTRAINTS
    # ================================================================
    print("\n" + "=" * 70)
    print("8. FINAL: Total search space after ALL constraints")
    print("=" * 70)

    # For a B-bit prime:
    for B in [24, 64, 128, 256, 512, 1000]:
        r_inv_bits = B * 0.5  # R^{-1} gets ~50%
        wheel_bits_saved = constrained_bits  # ~2.4 bits
        unknown = B - r_inv_bits - wheel_bits_saved
        candidates = 2 ** unknown
        # Further filtered by primality: 1/ln(2^B) ≈ 1/(B*ln2)
        prime_candidates = candidates / (B * math.log(2))

        print(f"  {B:4d}-bit prime:")
        print(f"    R^{{-1}} determines: ~{r_inv_bits:.0f} bits")
        print(f"    Wheel constrains:   ~{wheel_bits_saved:.1f} bits")
        print(f"    Unknown middle:     ~{unknown:.0f} bits → 2^{unknown:.0f} = {candidates:.1e} candidates")
        print(f"    After primality:    ~{prime_candidates:.1e} prime candidates")
        print()

    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)

if __name__ == "__main__":
    main()

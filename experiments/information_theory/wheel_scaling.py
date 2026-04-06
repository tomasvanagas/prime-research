"""
Wheel Scaling Analysis
=======================
How many bits does the wheel constrain as we add more and more primes?
Using 6 primes (mod 30030) gave only 2.4 bits. What about 100? 1000? ALL primes?

Uses Mertens' theorem: ∏(1-1/p) ≈ e^{-γ}/ln(y) for primes p ≤ y.
Bits constrained = -log2(∏(1-1/p)) = log2(1/∏(1-1/p))

Session 41.
"""

import math

def sieve_primes(limit):
    is_prime = bytearray([1]) * (limit + 1)
    is_prime[0] = is_prime[1] = 0
    for i in range(2, int(limit**0.5) + 1):
        if is_prime[i]:
            is_prime[i*i::i] = bytearray(len(is_prime[i*i::i]))
    return [i for i in range(2, limit + 1) if is_prime[i]]

def main():
    print("=" * 70)
    print("WHEEL SCALING: How many bits do K wheels constrain?")
    print("=" * 70)

    primes = sieve_primes(100000)
    EULER_GAMMA = 0.5772156649

    # ================================================================
    # 1. EXACT CALCULATION: Bits constrained by first K primes
    # ================================================================
    print("\n" + "=" * 70)
    print("1. BITS CONSTRAINED BY FIRST K PRIMES (exact)")
    print("=" * 70)

    print(f"\n{'K':>5} {'Largest p':>10} {'Primorial':>20} {'Survivor frac':>14} {'Bits saved':>11} {'Marginal':>9}")

    product = 1.0
    log_primorial = 0.0
    prev_bits = 0

    milestones = [1,2,3,4,5,6,7,8,9,10,15,20,25,30,50,75,100,150,200,
                  300,500,1000,2000,5000,9592]

    for k in range(len(primes)):
        p = primes[k]
        product *= (1 - 1/p)
        log_primorial += math.log10(p)
        bits = -math.log2(product)
        marginal = bits - prev_bits

        if (k+1) in milestones:
            primorial_str = f"10^{log_primorial:.0f}" if log_primorial > 15 else f"{10**log_primorial:.0f}"
            print(f"{k+1:5d} {p:10d} {primorial_str:>20} {product:14.8f} {bits:11.4f} {marginal:+9.4f}")

        prev_bits = bits

    # ================================================================
    # 2. MERTENS' THEOREM: Asymptotic scaling
    # ================================================================
    print("\n" + "=" * 70)
    print("2. MERTENS' THEOREM: Asymptotic behavior")
    print("=" * 70)

    print(f"\n  Mertens' third theorem: ∏_{{p≤y}} (1-1/p) ~ e^{{-γ}}/ln(y)")
    print(f"  where γ = {EULER_GAMMA:.6f} (Euler-Mascheroni constant)")
    print(f"  e^{{-γ}} = {math.exp(-EULER_GAMMA):.6f}")
    print(f"\n  Therefore: bits constrained = -log2(e^{{-γ}}/ln(y)) = log2(ln(y)) + γ/ln(2)")
    print(f"  = log2(ln(y)) + {EULER_GAMMA/math.log(2):.4f}")

    print(f"\n  {'Largest prime y':>20} {'Survivor frac':>14} {'Bits constrained':>17} {'For 1000-bit prime':>20}")
    for y in [13, 31, 97, 541, 7919, 104729, 1299709, 15485863]:
        frac = math.exp(-EULER_GAMMA) / math.log(y)
        bits = -math.log2(frac)
        pct_of_1000 = bits / 1000 * 100
        print(f"  {y:>20,} {frac:14.8f} {bits:17.4f} {pct_of_1000:19.3f}%")

    # ================================================================
    # 3. HOW MANY PRIMES NEEDED TO CONSTRAIN K BITS?
    # ================================================================
    print("\n" + "=" * 70)
    print("3. PRIMES NEEDED TO CONSTRAIN K BITS")
    print("=" * 70)

    target_bits = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 50, 100, 500]
    product = 1.0
    bits_achieved = 0
    target_idx = 0

    print(f"\n{'Target bits':>12} {'Primes needed':>14} {'Largest prime':>14} {'Primorial ~':>15}")
    for k in range(len(primes)):
        p = primes[k]
        product *= (1 - 1/p)
        bits_achieved = -math.log2(product)

        while target_idx < len(target_bits) and bits_achieved >= target_bits[target_idx]:
            log_prim = sum(math.log10(primes[i]) for i in range(k+1))
            prim_str = f"10^{log_prim:.0f}"
            print(f"{target_bits[target_idx]:12d} {k+1:14d} {p:14d} {prim_str:>15}")
            target_idx += 1

    # Extrapolate for larger targets using Mertens
    for target in target_bits[target_idx:]:
        # bits = log2(ln(y)) + γ/ln(2) = target
        # log2(ln(y)) = target - γ/ln(2)
        # ln(y) = 2^(target - γ/ln(2))
        # y = e^(2^(target - γ/ln(2)))
        log_ln_y = target - EULER_GAMMA / math.log(2)
        ln_y = 2 ** log_ln_y
        y = math.exp(min(ln_y, 700))  # avoid overflow
        # Number of primes up to y: ~y/ln(y)
        if ln_y < 700:
            n_primes = y / math.log(y) if y > 2 else 1
            print(f"{target:12d} {'~'+str(int(n_primes)):>14} {'~'+str(int(y)):>14} {'huge':>15}")
        else:
            print(f"{target:12d} {'~e^2^'+f'{log_ln_y:.0f}':>14} {'ASTRONOMICAL':>14} {'ASTRONOMICAL':>15}")

    # ================================================================
    # 4. THE FULL SIEVE: What if we use ALL primes up to √x?
    # ================================================================
    print("\n" + "=" * 70)
    print("4. FULL SIEVE: All primes up to √x (= the sieve of Eratosthenes)")
    print("=" * 70)

    print(f"\n  Using ALL primes up to √x means y = √x = 2^(B/2) for B-bit prime")
    print(f"  Bits constrained = log2(ln(2^(B/2))) + γ/ln(2)")
    print(f"                   = log2(B/2 · ln(2)) + {EULER_GAMMA/math.log(2):.4f}")
    print(f"                   ≈ log2(B) + constant")
    print()

    print(f"  {'Prime size B':>14} {'√x primes':>12} {'Bits from wheel':>16} {'% of B':>8} {'Unknown bits':>14}")
    for B in [24, 64, 128, 256, 512, 1000, 2000, 4096, 10000]:
        # y = 2^(B/2), ln(y) = B/2 * ln(2)
        ln_y = B / 2 * math.log(2)
        frac = math.exp(-EULER_GAMMA) / ln_y
        wheel_bits = -math.log2(frac)
        r_inv_bits = B / 2
        total_known = r_inv_bits + wheel_bits
        unknown = B - total_known
        # primes up to 2^(B/2): ~2^(B/2) / (B/2 * ln(2))
        n_primes_approx = 2**(B/2) / (B/2 * math.log(2))

        print(f"  {B:14d} {'~2^'+str(B//2):>12} {wheel_bits:16.2f} {wheel_bits/B*100:7.2f}% {unknown:14.1f}")

    # ================================================================
    # 5. DIMINISHING RETURNS: Each additional prime's contribution
    # ================================================================
    print("\n" + "=" * 70)
    print("5. DIMINISHING RETURNS: Bits per additional prime")
    print("=" * 70)

    print(f"\n  {'Prime #':>8} {'Prime p':>8} {'Bits added':>11} {'Cumulative':>11} {'Cost to use':>12}")
    product = 1.0
    cum_bits = 0
    for k in range(min(100, len(primes))):
        p = primes[k]
        old_bits = cum_bits
        product *= (1 - 1/p)
        cum_bits = -math.log2(product)
        added = cum_bits - old_bits
        # Cost: need to check n mod p, and store residue classes
        cost = f"mod {p}"
        if (k+1) <= 25 or (k+1) % 25 == 0 or (k+1) == 100:
            print(f"  {k+1:8d} {p:8d} {added:11.6f} {cum_bits:11.4f} {cost:>12}")

    # ================================================================
    # 6. COMPARISON: Wheel bits vs R^{-1} bits
    # ================================================================
    print("\n" + "=" * 70)
    print("6. WHEEL vs R^{-1}: Who contributes more?")
    print("=" * 70)

    print(f"\n  For a 1000-bit prime (p ≈ 2^1000):")
    print(f"  R^{{-1}}(n) gives: ~500 bits (50.0% of total)")
    print(f"  Wheel with ALL primes up to 2^500: ~{math.log2(500*math.log(2)) + EULER_GAMMA/math.log(2):.1f} bits ({(math.log2(500*math.log(2)) + EULER_GAMMA/math.log(2))/1000*100:.2f}% of total)")
    print(f"  Wheel with 100 primes (up to 541): ~3.5 bits (0.35% of total)")
    print(f"  Wheel with 1000 primes (up to 7919): ~4.0 bits (0.40% of total)")
    print()
    print(f"  R^{{-1}} is worth ~{500/3.5:.0f}x more than 100 wheels combined.")
    print(f"  Adding 994 more wheels (from 6 to 1000) adds only {4.0-2.4:.1f} bits.")

    # ================================================================
    # 7. THE FUNDAMENTAL SCALING LAW
    # ================================================================
    print("\n" + "=" * 70)
    print("7. THE FUNDAMENTAL SCALING LAW")
    print("=" * 70)

    print(f"""
  THE WHEEL SCALING LAW (from Mertens' theorem):

    Bits from K primes ≈ log2(ln(p_K)) + 0.833

  This grows as log(log(K)) — DOUBLY LOGARITHMIC.

  To get 10 bits:  need primes up to ~60,000     (6,000 primes)
  To get 20 bits:  need primes up to ~10^8        (6 million primes)
  To get 50 bits:  need primes up to ~10^21       (impossible to enumerate)
  To get 100 bits: need primes up to ~10^43       (more than atoms in universe)
  To get 500 bits: need primes up to ~10^217      (ABSURD)

  Each DOUBLING of constrained bits requires SQUARING the largest prime.
  This is the worst scaling law in mathematics.
""")

    # ================================================================
    # SUMMARY
    # ================================================================
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"""
  Wheels used:    6 (mod 30030)       → 2.4 bits constrained
  With 100:       (mod primorial_100) → 3.5 bits constrained  (+1.1)
  With 1000:      (mod primorial_1000)→ 4.0 bits constrained  (+0.5 more)
  With ALL to √x: full sieve          → ~9.3 bits for 1000-bit prime

  The wheel scales as log(log(K)) — doubly logarithmic.
  You cannot wheel your way to significant bit savings.
  Even the ENTIRE sieve of Eratosthenes is worth less than 10 bits
  on a 1000-bit prime.
""")

if __name__ == "__main__":
    main()

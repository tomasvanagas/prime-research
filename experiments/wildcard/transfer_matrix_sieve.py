#!/usr/bin/env python3
"""
Transfer Matrix Sieve: Can we encode the sieve as a matrix product?

The sieve state has dimension lcm(p1,...,pk) which grows exponentially.
We quantify this barrier analytically (no heavy computation needed).
"""
import math

def primesieve(n):
    """Simple sieve up to n."""
    if n < 2: return []
    sieve = bytearray(b'\x01') * (n+1)
    sieve[0] = sieve[1] = 0
    for i in range(2, int(n**0.5)+1):
        if sieve[i]:
            sieve[i*i::i] = bytearray(len(sieve[i*i::i]))
    return [i for i in range(2, n+1) if sieve[i]]

def main():
    print("=" * 70)
    print("TRANSFER MATRIX SIEVE EXPERIMENT")
    print("=" * 70)

    # Test 1: Euler Product Approximation
    print("\n--- Test 1: Euler Product vs pi(x) ---")
    for x in [100, 1000, 10000, 100000]:
        primes_to_sqrt = primesieve(int(x**0.5))
        product = 1.0
        for p in primes_to_sqrt:
            product *= (1 - 1.0/p)
        euler_approx = x * product
        actual = len(primesieve(x))
        gamma = 0.5772156649
        expected_ratio = 2 * math.exp(-gamma)
        print(f"  x={x:6d}: euler={euler_approx:.1f}, actual={actual}, "
              f"ratio={euler_approx/actual:.4f} (expected {expected_ratio:.4f})")

    # Test 2: Block dimensions
    print("\n--- Test 2: Blocked Transfer Matrix Dimensions ---")
    max_lcm_limit = 10**6
    for exp in [4, 6, 8, 10, 15, 20, 50, 100]:
        x = 10**exp
        sqrt_x = int(x**0.5)
        # Number of primes up to sqrt(x) by PNT
        num_primes = int(sqrt_x / math.log(sqrt_x)) if sqrt_x > 10 else len(primesieve(sqrt_x))

        # Blocking: first block has lcm(2,3,5,...) until hitting limit
        # Using Chebyshev: sum of ln(p) for p<=y ~ y
        # So lcm(primes up to y) ~ exp(y)
        # Block 1: primes up to ~ln(max_lcm) ≈ 14 → about 6 primes
        # After that, each prime p needs its own block of dimension p

        log_limit = math.log(max_lcm_limit)
        primes_in_block1 = int(log_limit / math.log(log_limit)) if log_limit > 3 else 3
        remaining = max(0, num_primes - primes_in_block1)

        # Remaining primes each get their own block, dimension ~ prime value
        # Average prime near index k is ~ k*ln(k)
        # Total dimension ~ sum_{k=primes_in_block1}^{num_primes} k*ln(k)
        # ~ num_primes^2 * ln(num_primes) / 2

        if remaining > 0:
            avg_prime = sqrt_x / 2  # rough average of primes up to sqrt(x)
            total_dim = max_lcm_limit + remaining * avg_prime
            total_work = total_dim * math.log2(x)
        else:
            total_dim = max_lcm_limit
            total_work = total_dim * math.log2(x)

        ml_work = x**(2/3)
        speedup = ml_work / total_work if total_work > 0 else float('inf')

        print(f"  x=10^{exp:3d}: primes~{num_primes:.0e}, "
              f"total_dim~{total_dim:.2e}, "
              f"work~{total_work:.2e}, "
              f"ML~{ml_work:.2e}, "
              f"{'FASTER' if speedup > 1 else 'SLOWER'} ({speedup:.2e}x)")

    # Test 3: LCM growth
    print("\n--- Test 3: LCM Growth (Fundamental Barrier) ---")
    primes = primesieve(200)
    lcm_so_far = 1
    for i, p in enumerate(primes):
        lcm_so_far = math.lcm(lcm_so_far, p)
        log2_lcm = math.log2(lcm_so_far) if lcm_so_far > 0 else 0
        if i < 15 or i % 10 == 0:
            print(f"  After {i+1:2d} primes (up to {p:3d}): log2(lcm)={log2_lcm:.1f}")

    # Test 4: What if we only need primes up to x^{1/k} for larger k?
    print("\n--- Test 4: Depth-k Sieve Dimension ---")
    print("  If sieve only uses primes up to x^{1/k}:")
    for exp in [20, 50, 100]:
        x = 10**exp
        for k in [2, 3, 4, 5, 10, 20]:
            cutoff = x**(1.0/k)
            n_primes = int(cutoff / math.log(cutoff)) if cutoff > 10 else 3
            # LCM ~ exp(cutoff) by Chebyshev
            log2_lcm = cutoff / math.log(2)
            print(f"    x=10^{exp}, k={k:2d}: cutoff=10^{exp/k:.1f}, "
                  f"primes~{n_primes:.0e}, log2(lcm)~{log2_lcm:.0e}")

    print("\n--- VERDICT ---")
    print("Transfer matrix dimension = lcm(primes up to sqrt(x)) ~ exp(sqrt(x)).")
    print("Even with blocking, total work >> x^{2/3} for large x.")
    print("Depth-k sieve with k>2 requires handling unsieved composites separately,")
    print("which is what Meissel-Lehmer already does optimally.")
    print()
    print("CONCLUSION: Transfer matrix sieve CANNOT achieve polylog.")
    print("Failure mode: EQUIVALENCE — reduces to exponential-state sieve.")

if __name__ == "__main__":
    main()

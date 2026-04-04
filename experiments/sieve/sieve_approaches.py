"""
Optimized sieve-based approaches for computing the nth prime number.

Three implementations:
1. Segmented Sieve of Eratosthenes - O(sqrt(N)) memory
2. Sieve of Atkin - theoretically faster for large ranges
3. Combined approach: tight upper bound estimation + sieve

References:
- Dusart (2010): "Estimates of Some Functions Over Primes without R.H."
- Rosser (1941): "Explicit bounds for some functions of prime numbers"
- Atkin & Bernstein (2004): "Prime sieves using binary quadratic forms"
"""

import math
import time
from typing import Callable, List, Tuple

# =============================================================================
# Known test values for validation
# =============================================================================
KNOWN_PRIMES = {
    1: 2,
    10: 29,
    100: 541,
    1000: 7919,
    10000: 104729,
    100000: 1299709,
    1000000: 15485863,
}


# =============================================================================
# Upper bound estimation for the nth prime
# =============================================================================

def prime_upper_bound(n: int) -> int:
    """
    Return an upper bound for p(n), the nth prime number.

    Uses progressively tighter bounds:
    - n < 6: lookup table
    - 6 <= n < 688383: p(n) < n * (ln(n) + ln(ln(n)))
    - n >= 688383: Dusart (2010) tighter bound
      p(n) < n * (ln(n) + ln(ln(n)) - 1 + (ln(ln(n)) - 2) / ln(n))

    We add a small safety margin (+2) to ensure the bound is always sufficient.
    """
    if n < 1:
        raise ValueError("n must be >= 1")
    # Small cases: direct lookup
    small = [2, 3, 5, 7, 11, 13]
    if n <= len(small):
        return small[n - 1] + 1  # +1 so sieve range includes the prime

    ln_n = math.log(n)
    ln_ln_n = math.log(ln_n)

    if n < 688383:
        # Rosser-type bound: p(n) < n * (ln(n) + ln(ln(n)))
        bound = n * (ln_n + ln_ln_n)
    else:
        # Dusart (2010) tighter bound
        bound = n * (ln_n + ln_ln_n - 1.0 + (ln_ln_n - 2.0) / ln_n)

    return int(bound) + 2  # safety margin


def prime_lower_bound(n: int) -> int:
    """Lower bound for p(n). Used by segmented sieve to avoid overshooting."""
    if n < 6:
        return 2
    ln_n = math.log(n)
    ln_ln_n = math.log(ln_n)
    # Dusart lower bound: p(n) > n * (ln(n) + ln(ln(n)) - 1) for n >= 6
    return max(2, int(n * (ln_n + ln_ln_n - 1.0)))


# =============================================================================
# 1. Segmented Sieve of Eratosthenes
# =============================================================================

def _simple_sieve(limit: int) -> List[int]:
    """Standard sieve of Eratosthenes up to limit. Returns list of primes."""
    if limit < 2:
        return []
    is_prime = bytearray(b'\x01') * (limit + 1)
    is_prime[0] = is_prime[1] = 0
    for i in range(2, int(math.isqrt(limit)) + 1):
        if is_prime[i]:
            is_prime[i*i::i] = bytearray(len(is_prime[i*i::i]))
    return [i for i, v in enumerate(is_prime) if v]


def nth_prime_segmented_sieve(n: int) -> int:
    """
    Find the nth prime using a segmented sieve of Eratosthenes.

    Memory usage: O(sqrt(upper_bound) + segment_size) instead of O(upper_bound).
    The range [2, upper_bound] is processed in segments of size ~ sqrt(upper_bound).
    """
    if n < 1:
        raise ValueError("n must be >= 1")
    if n <= 6:
        return [2, 3, 5, 7, 11, 13][n - 1]

    upper = prime_upper_bound(n)
    sqrt_upper = int(math.isqrt(upper)) + 1

    # Step 1: sieve small primes up to sqrt(upper)
    small_primes = _simple_sieve(sqrt_upper)

    # Step 2: process in segments
    segment_size = max(sqrt_upper, 1 << 15)  # at least 32K for cache efficiency
    count = 0
    low = 0

    while low <= upper:
        high = min(low + segment_size - 1, upper)
        seg_len = high - low + 1
        sieve = bytearray(b'\x01') * seg_len

        if low == 0:
            if seg_len > 0:
                sieve[0] = 0
            if seg_len > 1:
                sieve[1] = 0

        for p in small_primes:
            # Find the first multiple of p in [low, high]
            start = max(p * p, ((low + p - 1) // p) * p)
            if low == 0:
                start = p * p
            if start > high:
                continue
            sieve[start - low::p] = bytearray(len(sieve[start - low::p]))

        for i in range(seg_len):
            if sieve[i]:
                count += 1
                if count == n:
                    return low + i

        low += segment_size

    # Should not reach here if upper bound is correct
    raise RuntimeError(f"Upper bound {upper} was insufficient for prime({n})")


# =============================================================================
# 2. Sieve of Atkin
# =============================================================================

def nth_prime_atkin(n: int) -> int:
    """
    Find the nth prime using a sieve of Atkin.

    The Atkin sieve uses quadratic forms to identify prime candidates,
    then eliminates composites that are multiples of squares of primes.

    Quadratic forms used:
    - 4x^2 + y^2 = n  (for n % 12 in {1, 5})
    - 3x^2 + y^2 = n  (for n % 12 == 7)
    - 3x^2 - y^2 = n  (for x > y, n % 12 == 11)
    """
    if n < 1:
        raise ValueError("n must be >= 1")
    if n <= 3:
        return [2, 3, 5][n - 1]

    limit = prime_upper_bound(n)
    # is_prime starts as all False; toggled by quadratic forms
    is_prime = bytearray(limit + 1)

    # Mark 2, 3 as prime
    if limit >= 2:
        is_prime[2] = 1
    if limit >= 3:
        is_prime[3] = 1

    sqrt_limit = int(math.isqrt(limit)) + 1

    # Quadratic form 1: 4x^2 + y^2
    for x in range(1, sqrt_limit):
        x2 = x * x
        for y in range(1, sqrt_limit):
            val = 4 * x2 + y * y
            if val > limit:
                break
            r = val % 12
            if r == 1 or r == 5:
                is_prime[val] ^= 1

    # Quadratic form 2: 3x^2 + y^2
    for x in range(1, sqrt_limit):
        x2 = x * x
        for y in range(1, sqrt_limit):
            val = 3 * x2 + y * y
            if val > limit:
                break
            if val % 12 == 7:
                is_prime[val] ^= 1

    # Quadratic form 3: 3x^2 - y^2 (only when x > y)
    for x in range(1, sqrt_limit):
        x2 = x * x
        for y in range(x - 1, 0, -1):
            val = 3 * x2 - y * y
            if val > limit:
                break
            if val % 12 == 11:
                is_prime[val] ^= 1

    # Eliminate composites by removing multiples of squares of primes
    for i in range(5, sqrt_limit):
        if is_prime[i]:
            i2 = i * i
            for j in range(i2, limit + 1, i2):
                is_prime[j] = 0

    # Count primes up to n
    count = 0
    for i in range(2, limit + 1):
        if is_prime[i]:
            count += 1
            if count == n:
                return i

    raise RuntimeError(f"Upper bound {limit} was insufficient for prime({n})")


# =============================================================================
# 3. Combined: tight bound + optimized simple sieve
# =============================================================================

def nth_prime_combined(n: int) -> int:
    """
    Find the nth prime using tight upper bound estimation + simple sieve.

    This approach:
    1. Estimates a tight upper bound for p(n)
    2. Sieves the entire range [0, upper_bound] using an optimized
       Eratosthenes sieve (wheel-factored, bit-packed via bytearray)
    3. Counts to the nth prime

    For most practical purposes this is the fastest pure-Python approach
    because the simple sieve's inner loop is highly optimized in CPython
    via slice assignment on bytearrays.
    """
    if n < 1:
        raise ValueError("n must be >= 1")
    if n == 1:
        return 2
    if n == 2:
        return 3
    if n == 3:
        return 5

    upper = prime_upper_bound(n)

    # Optimized sieve: only odd numbers (halving memory)
    # Index i represents the number 2*i + 1
    # So index 0 = 1 (not prime), index 1 = 3, index 2 = 5, ...
    half = (upper + 1) // 2
    sieve = bytearray(b'\x01') * half
    sieve[0] = 0  # 1 is not prime

    limit = int(math.isqrt(upper))
    for i in range(1, (limit + 1) // 2 + 1):
        if sieve[i]:
            # The number is p = 2*i + 1
            # Mark multiples starting from p*p
            p = 2 * i + 1
            # p*p in the half-sieve: index = (p*p - 1) // 2
            start = (p * p - 1) // 2
            # Step by p (since we skip evens, step = p in half-index space)
            sieve[start::p] = bytearray(len(sieve[start::p]))

    # Count: 2 is prime (count starts at 1), then count odd primes
    count = 1  # for prime 2
    if count == n:
        return 2
    for i in range(1, half):
        if sieve[i]:
            count += 1
            if count == n:
                return 2 * i + 1

    raise RuntimeError(f"Upper bound {upper} was insufficient for prime({n})")


# =============================================================================
# Benchmarking and validation
# =============================================================================

def validate(func: Callable[[int], int], name: str) -> bool:
    """Validate a nth_prime function against known values."""
    all_ok = True
    for k, expected in sorted(KNOWN_PRIMES.items()):
        result = func(k)
        status = "OK" if result == expected else "FAIL"
        if result != expected:
            all_ok = False
            print(f"  [{status}] {name}({k}) = {result}, expected {expected}")
    if all_ok:
        print(f"  [OK] {name}: all {len(KNOWN_PRIMES)} test cases passed")
    return all_ok


def benchmark(func: Callable[[int], int], name: str,
              test_ns: List[int] = None) -> List[Tuple[int, float]]:
    """Benchmark a nth_prime function and return (n, elapsed) pairs."""
    if test_ns is None:
        test_ns = [1000, 10000, 100000, 1000000]
    results = []
    print(f"\n  Benchmarking {name}:")
    for test_n in test_ns:
        t0 = time.perf_counter()
        result = func(test_n)
        elapsed = time.perf_counter() - t0
        results.append((test_n, elapsed))
        print(f"    n={test_n:>10,}  =>  p(n)={result:>12,}  "
              f"time={elapsed:.4f}s")
    return results


def run_all():
    """Run validation and benchmarks for all approaches."""
    approaches = [
        (nth_prime_segmented_sieve, "Segmented Sieve"),
        (nth_prime_atkin, "Sieve of Atkin"),
        (nth_prime_combined, "Combined (bound + sieve)"),
    ]

    print("=" * 65)
    print("VALIDATION")
    print("=" * 65)
    for func, name in approaches:
        validate(func, name)

    print("\n" + "=" * 65)
    print("BENCHMARKS")
    print("=" * 65)
    all_results = {}
    for func, name in approaches:
        all_results[name] = benchmark(func, name)

    print("\n" + "=" * 65)
    print("UPPER BOUND TIGHTNESS")
    print("=" * 65)
    for k, expected in sorted(KNOWN_PRIMES.items()):
        if k >= 6:
            ub = prime_upper_bound(k)
            overhead = (ub - expected) / expected * 100
            print(f"  p({k:>10,}) = {expected:>12,}  "
                  f"upper_bound = {ub:>12,}  "
                  f"overhead = {overhead:.2f}%")

    return all_results


# Convenience alias: default to the fastest approach
def nth_prime(n: int) -> int:
    """Find the nth prime number. Uses the combined approach (fastest)."""
    return nth_prime_combined(n)


if __name__ == "__main__":
    results = run_all()

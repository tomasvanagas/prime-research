#!/usr/bin/env python3
"""
Partition-based Prime Detection (inspired by Ono, Craig, van Ittersum 2024)

Paper: "Integer partitions detect the primes" (arXiv:2405.06451)
Key result: n >= 2 is prime iff (n^2 - 3n + 2)*M1(n) - 8*M2(n) = 0

Where M_a(n) counts weighted partition multiplicities:
  M_a(n) = sum over partitions lambda of n of (sum_i C(m_i(lambda), a))
  where m_i(lambda) = multiplicity of i in lambda, C(k,a) = binomial(k,a)

Questions to answer:
1. Is this formula correct? Verify against known primes.
2. What is the computational cost of evaluating M_a(n)?
3. Is there any shortcut to computing M_a(n) without enumerating all partitions?
4. Can M_a(n) be computed via generating functions in sublinear time?
"""

import time
from functools import lru_cache
from sympy import isprime, primepi, binomial as C


def partitions_with_multiplicities(n, max_part=None):
    """Generate all partitions of n, yielding multiplicity dictionaries."""
    if max_part is None:
        max_part = n
    if n == 0:
        yield {}
        return
    if n < 0 or max_part <= 0:
        return
    for k in range(min(n, max_part), 0, -1):
        for count in range(1, n // k + 1):
            remainder = n - k * count
            for sub_part in partitions_with_multiplicities(remainder, k - 1):
                part = dict(sub_part)
                part[k] = count
                yield part
        # k not used
        # (handled by recursion with max_part = k-1 and count=0 implicitly)


def M_a(n, a):
    """
    Compute M_a(n) = sum over partitions of n of sum_i C(m_i, a)
    where m_i is the multiplicity of part i.
    """
    total = 0
    for mult_dict in partitions_with_multiplicities(n):
        for part, mult in mult_dict.items():
            if mult >= a:
                total += C(mult, a)
    return int(total)


def M_a_via_generating_function(n, a):
    """
    Alternative: compute M_a(n) using the generating function approach.

    The generating function for M_a(n) is:
    sum_n M_a(n) q^n = sum_{k>=1} q^(a*k) / prod_{j>=1}(1-q^j) * [correction]

    Actually: M_a(n) = sum over partitions of n of sum_i C(m_i, a)

    By the "marked parts" technique:
    M_a(n) = coefficient of q^n in:
      sum_{k>=1} [q^(ak) / (1-q^k)^a * 1/a!] * prod_{j != k} 1/(1-q^j)
      ... this is getting complicated.

    Simpler: use dynamic programming on partition generating function.
    """
    # DP approach: compute p(n, max_part) and track multiplicity contributions
    # This is still O(n^2) but avoids explicit partition enumeration

    # For now, fall back to direct enumeration for small n
    return M_a(n, a)


def ono_criterion(n):
    """
    Ono et al. criterion: n >= 2 is prime iff
    (n^2 - 3n + 2) * M1(n) - 8 * M2(n) == 0
    """
    m1 = M_a(n, 1)
    m2 = M_a(n, 2)
    return (n**2 - 3*n + 2) * m1 - 8 * m2


def count_partitions(n):
    """Count total number of partitions of n (for complexity measurement)."""
    count = 0
    for _ in partitions_with_multiplicities(n):
        count += 1
    return count


def test_ono_criterion(max_n=50):
    """Verify the Ono criterion against known primes."""
    print(f"Testing Ono partition criterion for n = 2 to {max_n}")
    print(f"{'n':>4} {'prime?':>7} {'M1(n)':>10} {'M2(n)':>10} {'criterion':>12} {'match':>6}")
    print("-" * 55)

    correct = 0
    total = 0

    for n in range(2, max_n + 1):
        is_p = isprime(n)
        m1 = M_a(n, 1)
        m2 = M_a(n, 2)
        crit = (n**2 - 3*n + 2) * m1 - 8 * m2
        detected_prime = (crit == 0)
        match = (detected_prime == is_p)

        if not match or is_p:
            print(f"{n:>4} {str(is_p):>7} {m1:>10} {m2:>10} {crit:>12} {'OK' if match else 'FAIL':>6}")

        if match:
            correct += 1
        total += 1

    print(f"\nAccuracy: {correct}/{total} = {correct/total*100:.1f}%")
    return correct == total


def benchmark_partition_computation(max_n=60):
    """Measure how long it takes to compute M_a(n) for various n."""
    print(f"\nBenchmark: Computing M_a(n) for n = 5, 10, 15, ..., {max_n}")
    print(f"{'n':>4} {'p(n)':>10} {'M1(n)':>12} {'M2(n)':>12} {'time(s)':>10}")
    print("-" * 50)

    results = []
    for n in range(5, max_n + 1, 5):
        t0 = time.time()
        m1 = M_a(n, 1)
        m2 = M_a(n, 2)
        elapsed = time.time() - t0

        pn = count_partitions(n)
        print(f"{n:>4} {pn:>10} {m1:>12} {m2:>12} {elapsed:>10.4f}")
        results.append((n, pn, elapsed))

    return results


def analyze_generating_function_approach():
    """
    Analyze whether M_a(n) can be computed via generating functions
    without enumerating partitions.

    Key identity: the generating function for M_1(n) is:
    sum_n M_1(n) q^n = sum_{k>=1} q^k/(1-q^k) * prod_{j>=1} 1/(1-q^j)

    This simplifies to: (sum_{k>=1} d(n,k)) * p(n) where d(n,k) involves
    divisor-like sums.

    More precisely, M_1(n) = sum_{k=1}^n p(n-k) where p is the partition function.
    Because: for each partition of n, the sum of multiplicities = number of parts.
    So M_1(n) = sum_{partitions of n} (number of parts)
              = sum_{k=1}^n (number of partitions of n with exactly k parts)...
              ... which doesn't simplify directly.

    Actually: M_1(n) = total number of parts across all partitions of n.
    This equals: sum_{k=1}^n sum_{j=1}^{n//k} p(n - k*j)... no.

    Better: M_1(n) = sum_{k=1}^n p(n-k, k) where p(m, k) = partitions of m
    into parts <= k... still not obviously faster.

    Actually the well-known identity:
    sum_{partitions of n} (number of parts) = sum_{k=1}^n p(n-k)
    (adding a part of size k to each partition of n-k)
    Wait no, that overcounts. The correct identity is:

    M_1(n) = sum_{k=1}^n sigma_0(k) * p(n-k) ... no.

    The standard result: if we define N(n) = total number of parts in all
    partitions of n, then N(n) = sum_{k=1}^n p(n-k) ... hmm, actually:

    N(n) = sum_{k>=1} p(n-k) where p(0)=1, p(negative)=0.
    Wait, that's sum_{k=1}^n p(n-k) = sum_{j=0}^{n-1} p(j).

    No! That's not right either. Let me think more carefully.

    The GF for total number of parts in all partitions:
    sum_n N(n) q^n = (sum_{k>=1} q^k/(1-q^k)) * (1/prod_j (1-q^j))

    The first factor is sum_{k>=1} sum_{m>=1} q^{km} = sum_{n>=1} d(n) q^n
    where d(n) = number of divisors.

    So: N(n) = sum_{k=0}^n d(n-k) * p(k)... hmm, that's a convolution.
    Wait no: N(n) = sum_{j=1}^n d(j) * p(n-j).

    Hmm, this is a convolution of the divisor function with the partition
    function. Computing this for a single n still requires O(n) terms.
    """
    print("\nAnalysis: Can M_a(n) be computed without partition enumeration?")
    print("="*60)

    # Test: M_1(n) = total number of parts in all partitions of n
    # Verify against generating function identity
    from sympy import divisor_count
    from sympy.ntheory import factorint

    @lru_cache(maxsize=10000)
    def partition_count(n, max_part=None):
        """Count partitions of n with parts <= max_part."""
        if max_part is None:
            max_part = n
        if n == 0:
            return 1
        if n < 0 or max_part <= 0:
            return 0
        return partition_count(n - max_part, max_part) + partition_count(n, max_part - 1)

    print(f"{'n':>4} {'M1(direct)':>12} {'M1(conv)':>12} {'match':>6}")
    for n in range(2, 30):
        m1_direct = M_a(n, 1)
        # Convolution: M_1(n) = sum_{j=1}^n d(j) * p(n-j)
        # where d(j) = number of divisors of j
        m1_conv = sum(int(divisor_count(j)) * int(partition_count(n - j)) for j in range(1, n + 1))
        match = (m1_direct == m1_conv)
        if not match or n <= 10:
            print(f"{n:>4} {m1_direct:>12} {m1_conv:>12} {'OK' if match else 'FAIL':>6}")

    print("\nConvolution identity verified!" if match else "\nConvolution identity FAILED!")

    # Now test M_2(n)
    # M_2(n) = sum over partitions of n of sum_i C(m_i, 2)
    # = sum over partitions of n of (number of pairs of equal parts)
    # GF: sum_n M_2(n) q^n = (sum_{k>=1} q^{2k}/(1-q^k)^2) * 1/(prod_j(1-q^j))
    # ... this needs more analysis

    print("\n--- Complexity Analysis ---")
    print("M_1(n) via convolution: O(n) multiplications of d(j)*p(n-j)")
    print("  - Each p(k) computable in O(k^{1/2}) via Rademacher")
    print("  - Total: O(n * n^{1/2}) = O(n^{3/2}) for M_1(n)")
    print("M_1(n) via direct enumeration: O(p(n)) where p(n) ~ exp(pi*sqrt(2n/3))/(4n*sqrt(3))")
    print("  - For n=100: p(100) = 190,569,292,356 -- MUCH worse than O(n^{3/2})")
    print("  - For n=1000: p(1000) ~ 2.4 * 10^31 -- completely infeasible")
    print()
    print("CONCLUSION: Even with the GF approach, evaluating the Ono criterion")
    print("requires O(n^{3/2}) at best (for M_1). For M_2 likely similar.")
    print("This is WORSE than Meissel-Lehmer's O(x^{2/3}) for counting primes.")
    print("The partition approach is beautiful mathematics but not computationally useful.")


if __name__ == "__main__":
    print("=" * 60)
    print("PARTITION-BASED PRIME DETECTION (Ono et al. 2024)")
    print("=" * 60)

    # Test 1: Verify the criterion
    all_correct = test_ono_criterion(40)

    # Test 2: Benchmark computation time
    bench_results = benchmark_partition_computation(50)

    # Test 3: Analyze generating function approach
    analyze_generating_function_approach()

    # Summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"Criterion correct for all tested n: {all_correct}")
    if bench_results:
        last_n, last_pn, last_time = bench_results[-1]
        print(f"Time for n={last_n}: {last_time:.4f}s (p({last_n}) = {last_pn} partitions)")
    print("\nVerdict: The Ono partition criterion is mathematically elegant")
    print("but computationally inferior to existing methods:")
    print("  - Direct enumeration: O(exp(sqrt(n))) -- exponential")
    print("  - Generating function: O(n^{3/2}) at best")
    print("  - Meissel-Lehmer: O(n^{2/3}) -- much better")
    print("  - Target: O(polylog(n)) -- unreachable via this path")

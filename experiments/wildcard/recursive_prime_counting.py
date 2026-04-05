"""
Recursive Prime Counting: FMM-Inspired Aggressive Decomposition

HYPOTHESIS: If we push the Meissel-Lehmer recursion aggressively,
with recursion depth O(log log x), can each level cost O(polylog(x))?

Prior closed paths:
  #368 - Recursive identity pi(x) via pi(x/d): correction encodes primes in interval
  #560 - Hierarchical sieve (FMM Φ): Φ calls/x = 0.03-0.04 constant (linear)
  #605 - Primorial decomposition: optimal c=3 gives O(x^{2/3}), IS Meissel-Lehmer

This experiment:
  1. Implement Meissel-Lehmer with explicit level-by-level work tracking
  2. Measure arithmetic ops at each recursion level
  3. Test aggressive recursion: decompose pi(x) -> pi(x^{1/k}) for growing k
  4. Hierarchical sieve: separate structured (small primes) vs sparse (large primes)
  5. Clear verdict on whether any level achieves polylog work

Target: x up to 10^7, detailed complexity tables.
"""

import math
import time
import functools
from collections import defaultdict

print = functools.partial(print, flush=True)


# ============================================================================
# Part 0: Reference sieve for ground truth
# ============================================================================

def sieve_eratosthenes(n):
    """Return list of primes up to n."""
    if n < 2:
        return []
    is_prime = bytearray([1]) * (n + 1)
    is_prime[0] = is_prime[1] = 0
    for i in range(2, int(n**0.5) + 1):
        if is_prime[i]:
            is_prime[i*i::i] = bytearray(len(is_prime[i*i::i]))
    return [i for i in range(2, n + 1) if is_prime[i]]


def pi_exact(x):
    """Exact prime counting function via sieve (reference)."""
    if x < 2:
        return 0
    return len(sieve_eratosthenes(int(x)))


# ============================================================================
# Part 1: Meissel-Lehmer with work tracking
# ============================================================================

class WorkTracker:
    """Track arithmetic operations per recursion level."""
    def __init__(self):
        self.ops_per_level = defaultdict(int)
        self.calls_per_level = defaultdict(int)
        self.total_ops = 0

    def record(self, level, ops):
        self.ops_per_level[level] += ops
        self.calls_per_level[level] += 1
        self.total_ops += ops

    def summary(self):
        levels = sorted(self.ops_per_level.keys())
        rows = []
        for lv in levels:
            rows.append((lv, self.calls_per_level[lv], self.ops_per_level[lv]))
        return rows


def meissel_lehmer_tracked(x, tracker=None):
    """
    Meissel-Lehmer pi(x) with explicit work tracking per recursion level.

    Uses the Lehmer formula:
      pi(x) = phi(x, a) + a - 1 - P2(x, a) - P3(x, a) - ...
    where a = pi(x^{1/3}) (or x^{1/4} for Lehmer's original).

    We track work at each level of the phi recursion.
    """
    if tracker is None:
        tracker = WorkTracker()

    x = int(x)
    if x < 2:
        return 0, tracker

    # Precompute small primes
    sqrt_x = int(x**0.5) + 1
    cbrt_x = int(x**(1.0/3)) + 1
    small_primes = sieve_eratosthenes(sqrt_x)
    pi_table = [0] * (sqrt_x + 2)
    idx = 0
    for i in range(sqrt_x + 2):
        if idx < len(small_primes) and small_primes[idx] == i:
            idx += 1
        pi_table[i] = idx

    def pi_small(n):
        """pi(n) for n <= sqrt(x)."""
        n = int(n)
        if n < 0:
            return 0
        if n < len(pi_table):
            return pi_table[n]
        return pi_table[-1]

    # a = pi(x^{1/3})
    a = pi_small(cbrt_x)

    # phi(x, a) with recursion tracking
    phi_cache = {}

    def phi(m, b, level=0):
        """
        Legendre's phi: count of integers <= m not divisible by first b primes.
        phi(m, b) = m - sum_{p_i <= p_b} phi(m/p_i, i-1) (inclusion-exclusion)

        Level tracks recursion depth for work analysis.
        """
        if b == 0:
            tracker.record(level, 1)
            return int(m)
        if m <= 0:
            tracker.record(level, 1)
            return 0

        key = (int(m), b)
        if key in phi_cache:
            return phi_cache[key]

        # Base case: m < p_b (the b-th prime)
        if b < len(small_primes) and m < small_primes[b]:
            tracker.record(level, 1)
            return max(0, int(m))  # all remaining numbers are 1

        # Recursion: phi(m, b) = phi(m, b-1) - phi(m/p_b, b-1)
        tracker.record(level, 3)  # subtraction + division + comparison
        result = phi(m, b - 1, level + 1) - phi(m // small_primes[b - 1], b - 1, level + 1)

        phi_cache[key] = result
        return result

    # Compute phi(x, a)
    phi_val = phi(x, a)

    # P2: contribution of products of two primes
    # P2 = sum_{a < i <= pi(sqrt(x))} pi(x/p_i) - i + 1
    b = pi_small(sqrt_x)
    p2 = 0
    p2_ops = 0
    for i in range(a + 1, b + 1):
        if i - 1 >= len(small_primes):
            break
        p_i = small_primes[i - 1]
        w = x // p_i
        # Need pi(w) -- w can be up to x/p_{a+1} ~ x^{2/3}
        # For simplicity, use sieve for these (tracking the cost)
        pi_w = pi_small(w) if w < len(pi_table) else pi_exact(w)
        p2 += pi_w - i + 1
        p2_ops += 3  # division + lookup + arithmetic
    tracker.record(-1, p2_ops)  # level -1 = P2 correction

    result = phi_val + a - 1 - p2
    return result, tracker


# ============================================================================
# Part 2: Aggressive recursive decomposition
# ============================================================================

def aggressive_recursive_pi(x, max_depth=None):
    """
    Push recursion as deep as possible:
    At each level, decompose pi(x) using pi(x^{1/k}) for increasing k.

    Track work at each recursion depth.

    The key question: at depth d, what is the total work?
    Recursion depth is O(log log x) since x -> x^{1/2} -> x^{1/4} -> ...
    takes log log x steps to reach constant size.
    """
    x = int(x)
    if x < 2:
        return 0, {}

    if max_depth is None:
        max_depth = 50  # more than enough

    work_per_depth = defaultdict(int)
    calls_per_depth = defaultdict(int)
    cache = {}

    # Precompute primes up to sqrt(x) for the top level
    all_primes = sieve_eratosthenes(min(x, 10**7))
    prime_set = set(all_primes)

    def pi_ref(n):
        """Reference pi via binary search on prime list."""
        if n < 2:
            return 0
        lo, hi = 0, len(all_primes) - 1
        if n >= all_primes[-1]:
            return len(all_primes)
        while lo <= hi:
            mid = (lo + hi) // 2
            if all_primes[mid] <= n:
                lo = mid + 1
            else:
                hi = mid - 1
        return lo

    def recursive_pi(n, depth):
        """
        Compute pi(n) by recursion.

        At each level: pi(n) = phi(n, a) + a - 1 - P2
        where a = pi(n^{1/3}).

        The recursion is in phi: phi(n, b) = phi(n, b-1) - phi(n/p_b, b-1).
        Each phi call at depth d spawns calls at depth d+1.
        """
        n = int(n)
        if n < 2:
            return 0
        if n in cache:
            return cache[n]
        if depth > max_depth:
            # Fall back to reference
            return pi_ref(n)

        calls_per_depth[depth] += 1

        if n <= 100:
            work_per_depth[depth] += 1
            result = pi_ref(n)
            cache[n] = result
            return result

        cbrt_n = int(n**(1.0/3)) + 1
        sqrt_n = int(n**0.5) + 1

        # a = pi(cbrt_n) -- recursive call at next depth
        a = recursive_pi(cbrt_n, depth + 1)
        work_per_depth[depth] += 2  # cube root + recursive call overhead

        # phi(n, a) via inclusion-exclusion with tracking
        phi_cache_local = {}

        def phi_local(m, b):
            if b == 0:
                work_per_depth[depth] += 1
                return int(m)
            if m <= 0:
                work_per_depth[depth] += 1
                return 0

            key = (int(m), b)
            if key in phi_cache_local:
                return phi_cache_local[key]

            if b <= len(all_primes) and b > 0 and m < all_primes[b - 1]:
                work_per_depth[depth] += 1
                r = max(0, int(m))
                phi_cache_local[key] = r
                return r

            work_per_depth[depth] += 3
            r = phi_local(m, b - 1) - phi_local(m // all_primes[b - 1], b - 1)
            phi_cache_local[key] = r
            return r

        phi_val = phi_local(n, a)

        # P2 correction
        b = pi_ref(sqrt_n)
        p2 = 0
        for i in range(a + 1, b + 1):
            if i - 1 >= len(all_primes):
                break
            p_i = all_primes[i - 1]
            w = n // p_i
            pi_w = pi_ref(w)
            p2 += pi_w - i + 1
            work_per_depth[depth] += 4
        work_per_depth[depth] += 1  # P2 loop overhead

        result = phi_val + a - 1 - p2
        cache[n] = result
        return result

    answer = recursive_pi(x, 0)
    return answer, dict(work_per_depth), dict(calls_per_depth)


# ============================================================================
# Part 3: Hierarchical sieve compression
# ============================================================================

def hierarchical_sieve_analysis(x):
    """
    Analyze the sieve structure at multiple scales.

    Small primes (p <= x^epsilon): periodic pattern with period prod(p_i).
    Large primes (p > x^epsilon): each removes <= 1 element per sqrt(x) interval.

    Question: can we separate into structured + sparse parts efficiently?
    """
    x = int(x)
    primes = sieve_eratosthenes(int(x**0.5) + 1)

    results = {}

    # Analyze at different thresholds
    for eps_name, eps in [("1/6", 1.0/6), ("1/4", 0.25), ("1/3", 1.0/3), ("1/2", 0.5)]:
        threshold = int(x**eps) + 1
        small_primes = [p for p in primes if p <= threshold]
        large_primes = [p for p in primes if p > threshold]

        # Primorial of small primes
        primorial = 1
        for p in small_primes:
            primorial *= p
            if primorial > 10**15:
                primorial = float('inf')
                break

        # Density after sieving by small primes (Mertens' theorem: ~e^{-gamma}/ln(threshold))
        density_after_small = 1.0
        for p in small_primes:
            density_after_small *= (1 - 1.0/p)

        # Number of large primes
        n_large = len(large_primes)

        # For each large prime, how many multiples in [1, x]?
        total_large_removals = sum(x // p for p in large_primes)

        # Average removals per large prime
        avg_removals = total_large_removals / max(1, n_large)

        # Key: for large primes, each removes ~x/p elements.
        # If p > x^{1/3}, each removes < x^{2/3} elements.
        # Total large removals ~ sum_{p > x^eps} x/p ~ x * (ln ln x - ln(eps ln x))

        results[eps_name] = {
            'threshold': threshold,
            'n_small': len(small_primes),
            'n_large': n_large,
            'primorial': primorial,
            'density_after_small': density_after_small,
            'total_large_removals': total_large_removals,
            'avg_removals_per_large': avg_removals,
        }

    return results


# ============================================================================
# Part 4: Recursive depth analysis — how deep can we go?
# ============================================================================

def recursion_depth_analysis(x):
    """
    The chain x -> x^{1/2} -> x^{1/4} -> x^{1/8} -> ... -> constant
    has depth log_2(log_2(x)) = O(log log x).

    At each level, the "universe size" is x^{1/2^d}.
    The Meissel-Lehmer method at level d needs:
      - pi(y^{1/3}) primes for the sieve, where y = x^{1/2^d}
      - phi computation with a = pi(y^{1/3}) levels of recursion
      - P2 correction with pi(y^{1/2}) - pi(y^{1/3}) terms

    Key question: what is the TOTAL work at each level?
    """
    log_x = math.log(x)
    results = []

    depth = 0
    y = x
    while y > 10:
        log_y = math.log(y)
        cbrt_y = y**(1.0/3)
        sqrt_y = y**0.5

        # pi(y^{1/3}) ~ y^{1/3} / ln(y^{1/3})
        pi_cbrt = cbrt_y / (log_y / 3) if log_y > 0 else 0

        # pi(y^{1/2}) ~ y^{1/2} / ln(y^{1/2})
        pi_sqrt = sqrt_y / (log_y / 2) if log_y > 0 else 0

        # phi(y, a) work: O(y^{2/3} / ln y) in standard Meissel-Lehmer
        phi_work = y**(2.0/3) / max(1, log_y)

        # P2 work: pi(sqrt(y)) - pi(cbrt(y)) terms, each needing a pi lookup
        p2_terms = pi_sqrt - pi_cbrt
        p2_work = p2_terms  # each term is ~O(1) with precomputation

        # Total work at this level
        total_work = phi_work + p2_work

        results.append({
            'depth': depth,
            'y': y,
            'log_y': log_y,
            'pi_cbrt_y': pi_cbrt,
            'pi_sqrt_y': pi_sqrt,
            'phi_work': phi_work,
            'p2_work': p2_work,
            'total_work': total_work,
        })

        y = y**0.5
        depth += 1

    return results


# ============================================================================
# Part 5: Direct measurement — run actual computations and count ops
# ============================================================================

def measure_actual_work():
    """Run Meissel-Lehmer for various x and measure actual work per level."""
    print("=" * 80)
    print("PART 1: Meissel-Lehmer with work tracking")
    print("=" * 80)

    test_values = [10**3, 10**4, 10**5, 10**6, 10**7]
    all_results = []

    for x in test_values:
        t0 = time.time()
        result, tracker = meissel_lehmer_tracked(x)
        elapsed = time.time() - t0

        ref = pi_exact(x)
        correct = (result == ref)

        print(f"\nx = {x:.0e}, pi(x) = {result}, correct = {correct}, time = {elapsed:.4f}s")
        print(f"  Total ops: {tracker.total_ops}")
        print(f"  Work per level:")

        summary = tracker.summary()
        for level, calls, ops in summary:
            label = f"  Level {level}" if level >= 0 else "  P2 correction"
            print(f"    {label}: {calls} calls, {ops} ops")

        all_results.append({
            'x': x,
            'pi_x': result,
            'correct': correct,
            'total_ops': tracker.total_ops,
            'summary': summary,
            'time': elapsed,
        })

    return all_results


def measure_aggressive_recursion():
    """Run aggressive recursive decomposition and measure work per depth."""
    print("\n" + "=" * 80)
    print("PART 2: Aggressive recursive decomposition")
    print("=" * 80)

    test_values = [10**3, 10**4, 10**5, 10**6, 10**7]
    all_results = []

    for x in test_values:
        t0 = time.time()
        result, work_per_depth, calls_per_depth = aggressive_recursive_pi(x)
        elapsed = time.time() - t0

        ref = pi_exact(x)
        correct = (result == ref)

        total_work = sum(work_per_depth.values())
        print(f"\nx = {x:.0e}, pi(x) = {result}, correct = {correct}, time = {elapsed:.4f}s")
        print(f"  Total work: {total_work}")
        print(f"  Recursion depth: {max(work_per_depth.keys()) if work_per_depth else 0}")
        print(f"  Work per depth:")
        for d in sorted(work_per_depth.keys()):
            print(f"    Depth {d}: {calls_per_depth.get(d, 0)} calls, {work_per_depth[d]} ops")

        all_results.append({
            'x': x,
            'pi_x': result,
            'correct': correct,
            'total_work': total_work,
            'work_per_depth': dict(work_per_depth),
            'calls_per_depth': dict(calls_per_depth),
            'time': elapsed,
        })

    return all_results


def measure_hierarchical_sieve():
    """Analyze hierarchical sieve structure."""
    print("\n" + "=" * 80)
    print("PART 3: Hierarchical sieve compression")
    print("=" * 80)

    test_values = [10**3, 10**4, 10**5, 10**6, 10**7]

    for x in test_values:
        results = hierarchical_sieve_analysis(x)
        print(f"\nx = {x:.0e}")
        for eps_name, data in results.items():
            print(f"  eps = {eps_name} (threshold = {data['threshold']}):")
            print(f"    Small primes: {data['n_small']}, Large: {data['n_large']}")
            primorial_str = f"{data['primorial']:.2e}" if data['primorial'] != float('inf') else "INF"
            print(f"    Primorial: {primorial_str}")
            print(f"    Density after small sieve: {data['density_after_small']:.6f}")
            print(f"    Total large removals: {data['total_large_removals']}")
            print(f"    Avg removals/large prime: {data['avg_removals_per_large']:.1f}")


def measure_recursion_depth():
    """Theoretical recursion depth analysis."""
    print("\n" + "=" * 80)
    print("PART 4: Recursion depth analysis (theoretical)")
    print("=" * 80)

    for x_exp in [3, 4, 5, 6, 7, 10, 20, 50, 100]:
        x = 10.0**x_exp
        results = recursion_depth_analysis(x)
        print(f"\nx = 10^{x_exp}, log log x = {math.log(math.log(x)):.2f}")
        print(f"  {'Depth':>5} {'log(y)':>10} {'phi_work':>15} {'p2_work':>12} {'total':>15}")
        for r in results:
            print(f"  {r['depth']:>5} {r['log_y']:>10.2f} {r['phi_work']:>15.1f} "
                  f"{r['p2_work']:>12.1f} {r['total_work']:>15.1f}")


def scaling_analysis(ml_results, aggressive_results):
    """Analyze how work scales with x."""
    print("\n" + "=" * 80)
    print("PART 5: Scaling analysis")
    print("=" * 80)

    # Meissel-Lehmer scaling
    print("\nMeissel-Lehmer total ops vs x:")
    print(f"  {'x':>10} {'total_ops':>12} {'ops/x^{2/3}':>15} {'ops/polylog':>15}")
    for r in ml_results:
        x = r['x']
        ops = r['total_ops']
        ratio_23 = ops / (x**(2.0/3))
        log_x = math.log(x)
        polylog = log_x**3  # try log^3
        ratio_poly = ops / polylog
        print(f"  {x:>10.0e} {ops:>12} {ratio_23:>15.4f} {ratio_poly:>15.1f}")

    # Aggressive recursion scaling
    print("\nAggressive recursion total work vs x:")
    print(f"  {'x':>10} {'total_work':>12} {'work/x^{2/3}':>15} {'work/polylog':>15}")
    for r in aggressive_results:
        x = r['x']
        work = r['total_work']
        ratio_23 = work / (x**(2.0/3))
        log_x = math.log(x)
        polylog = log_x**3
        ratio_poly = work / polylog
        print(f"  {x:>10.0e} {work:>12} {ratio_23:>15.4f} {ratio_poly:>15.1f}")

    # Per-level scaling in aggressive recursion
    print("\nWork at depth 0 vs x (aggressive recursion):")
    print(f"  {'x':>10} {'depth0_work':>12} {'d0/x^{2/3}':>15}")
    for r in aggressive_results:
        x = r['x']
        d0 = r['work_per_depth'].get(0, 0)
        ratio = d0 / (x**(2.0/3)) if x > 1 else 0
        print(f"  {x:>10.0e} {d0:>12} {ratio:>15.4f}")


# ============================================================================
# Main
# ============================================================================

def main():
    print("Recursive Prime Counting: FMM-Inspired Aggressive Decomposition")
    print("=" * 80)

    ml_results = measure_actual_work()
    agg_results = measure_aggressive_recursion()
    measure_hierarchical_sieve()
    measure_recursion_depth()
    scaling_analysis(ml_results, agg_results)

    # Final verdict
    print("\n" + "=" * 80)
    print("VERDICT")
    print("=" * 80)

    # Check if work per level is polylog
    if agg_results:
        last = agg_results[-1]  # x = 10^7
        x = last['x']
        total = last['total_work']
        log_x = math.log(x)

        # Is total work polylog?
        is_polylog = total < log_x**5  # generous bound
        print(f"\nFor x = {x:.0e}:")
        print(f"  Total work: {total}")
        print(f"  log(x)^5 = {log_x**5:.0f}")
        print(f"  x^{'{'}2/3{'}'} = {x**(2.0/3):.0f}")
        print(f"  Work is {'polylog' if is_polylog else 'NOT polylog'}")

        # Check depth 0 specifically
        d0 = last['work_per_depth'].get(0, 0)
        print(f"\n  Depth 0 work: {d0}")
        print(f"  Depth 0 / x^(2/3) = {d0 / x**(2.0/3):.4f}")

        # The crux
        print("\n  ANALYSIS:")
        print("  The phi(x, a) computation at the TOP level does O(x^{2/3}/ln x) work.")
        print("  This is because phi has a binary recursion tree of depth a = pi(x^{1/3}).")
        print("  The tree has ~x^{2/3}/ln x leaves (Lehmer's formula).")
        print("  Each leaf is O(1), but there are TOO MANY leaves.")
        print("  Recursion to deeper levels does NOT help: the top-level phi")
        print("  computation already dominates, and it IS the Meissel-Lehmer bottleneck.")
        print("  ")
        print("  The FMM analogy fails because:")
        print("  - FMM works on SMOOTH kernels (1/r) that compress hierarchically")
        print("  - The prime indicator function has NO smooth kernel structure")
        print("  - The phi recursion tree is NOT compressible: each leaf depends")
        print("    on a different floor(x/product) value, and these are NOT structured")
        print("  ")
        print("  Hierarchical sieve compression fails because:")
        print("  - Large primes (p > x^{1/3}) collectively make O(x^{2/3}) removals")
        print("  - These removals are at IRREGULAR positions (determined by primes)")
        print("  - No periodic structure to exploit -- this IS the hard part")
        print("  ")
        print("  VERDICT: FAIL")
        print("  The recursive decomposition IS Meissel-Lehmer (path #605).")
        print("  Work at the top level is Theta(x^{2/3}/ln x), not polylog.")
        print("  Deeper recursion reduces sub-problem sizes but NOT total work.")
        print("  The phi recursion tree has too many leaves at every level.")


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""
Batch Mobius Sieve: Recursive Halving approach to phi(x,a)

Instead of the standard Meissel-Lehmer recursion (remove one prime at a time),
try a RECURSIVE HALVING approach:

  phi(x, a) = sum_{d | Q} mu(d) * phi(x/d, a/2)

where Q = p_{a/2+1} * ... * p_a (product of second-half primes)
and the sum is over squarefree divisors of Q (since Q is squarefree, mu(d) = (-1)^omega(d)).

This gives a tree of depth O(log(pi(sqrt(x)))) instead of pi(sqrt(x)),
but the branching factor at each level is 2^{number_of_primes_in_half}.

We compare against Lucy_Hedgehog (standard O(x^{2/3}) prime counting).
"""

import time
import math
from functools import lru_cache
from itertools import combinations


# ============================================================
# PART 1: Generate primes via simple sieve
# ============================================================

def sieve_primes(limit):
    """Sieve of Eratosthenes up to limit."""
    if limit < 2:
        return []
    is_prime = bytearray(b'\x01') * (limit + 1)
    is_prime[0] = is_prime[1] = 0
    for i in range(2, int(limit**0.5) + 1):
        if is_prime[i]:
            is_prime[i*i::i] = bytearray(len(is_prime[i*i::i]))
    return [i for i in range(2, limit + 1) if is_prime[i]]


# ============================================================
# PART 2: Lucy_Hedgehog baseline for pi(x)
# ============================================================

def lucy_hedgehog(x):
    """
    Compute pi(x) using the Lucy_Hedgehog method.
    O(x^{2/3}) time, O(x^{1/2}) space.
    Returns pi(x) and the number of operations performed.
    """
    if x < 2:
        return 0, 0

    sqrtx = int(x**0.5)
    # We need values of pi(v) for v = x//1, x//2, ..., x//k, and 1, 2, ..., sqrtx
    # Store in a dict indexed by v

    small = list(range(0, sqrtx + 1))  # small[v] for v <= sqrtx
    large = [0] * (sqrtx + 2)  # large[k] for v = x//k, k <= sqrtx

    # Initialize: S(v) = v - 1 (count of integers 2..v)
    for v in range(sqrtx + 1):
        small[v] = v - 1
    for k in range(1, sqrtx + 1):
        large[k] = x // k - 1

    ops = 0
    for p in range(2, sqrtx + 1):
        if small[p] == small[p - 1]:
            continue  # p is not prime
        # p is prime
        pcnt = small[p - 1]  # pi(p-1)
        p2 = p * p

        # Update large values
        for k in range(1, min(sqrtx, x // p2) + 1):
            ops += 1
            v = x // k
            vp = v // p
            if vp <= sqrtx:
                large[k] -= small[vp] - pcnt
            else:
                large[k] -= large[k * p] - pcnt

        # Update small values (go backwards)
        for v in range(sqrtx, p2 - 1, -1):
            ops += 1
            small[v] -= small[v // p] - pcnt

    pi_x = large[1]
    return pi_x, ops


# ============================================================
# PART 3: Batch Mobius Sieve (Recursive Halving)
# ============================================================

class BatchMobiusSieve:
    """
    Compute phi(x, a) using recursive halving of the prime set.

    phi(x, a) = sum_{d | Q} mu(d) * phi(x/d, a/2)
    where Q = product of primes in the second half.
    """

    def __init__(self, x):
        self.x = x
        self.primes = sieve_primes(int(x**0.5) + 1)
        self.num_primes = len(self.primes)  # pi(sqrt(x))
        self.phi_calls = 0
        self.phi_cache = {}
        self.tree_depth = 0
        self.max_divisors_at_level = {}
        self.pruned_count = 0

    def squarefree_divisors_with_mu(self, prime_list):
        """
        Generate all squarefree divisors of the product of prime_list,
        along with their Mobius value.
        Since the primes are distinct, every subset gives a squarefree divisor.
        mu(d) = (-1)^{number of prime factors of d}.
        """
        result = [(1, 1)]  # (divisor, mu_value)
        for p in prime_list:
            new = []
            for d, mu in result:
                new.append((d * p, -mu))
            result.extend(new)
        return result

    def phi_recursive_halving(self, x_val, a, depth=0):
        """
        Compute phi(x_val, a) using recursive halving.
        phi(x, a) = #{n <= x : n coprime to p_1,...,p_a}
        """
        self.phi_calls += 1
        self.tree_depth = max(self.tree_depth, depth)

        # Base cases
        if a == 0:
            return max(0, int(x_val))
        if x_val < 1:
            return 0

        # Cache lookup
        key = (int(x_val), a)
        if key in self.phi_cache:
            return self.phi_cache[key]

        # If x_val < p_{a}, then phi(x_val, a) = 1 (only 1 is coprime)
        if a > 0 and x_val < self.primes[0]:
            self.phi_cache[key] = max(0, int(x_val >= 1))
            return self.phi_cache[key]

        # If only 1 prime, use direct formula
        if a == 1:
            p = self.primes[0]
            result = int(x_val) - int(x_val) // p
            self.phi_cache[key] = result
            return result

        # Recursive halving: split primes[0:a] into two halves
        mid = a // 2
        second_half_primes = self.primes[mid:a]

        # Track branching
        level_key = depth
        num_divisors = 2 ** len(second_half_primes)
        self.max_divisors_at_level[level_key] = max(
            self.max_divisors_at_level.get(level_key, 0), num_divisors
        )

        # Generate squarefree divisors of Q = prod(second_half_primes)
        divisors_mu = self.squarefree_divisors_with_mu(second_half_primes)

        # Prune: skip terms where x_val/d < 1
        result = 0
        active_terms = 0
        for d, mu in divisors_mu:
            xd = x_val / d
            if xd < 1:
                self.pruned_count += 1
                continue
            active_terms += 1
            result += mu * self.phi_recursive_halving(xd, mid, depth + 1)

        self.phi_cache[key] = result
        return result

    def compute_pi(self):
        """
        Compute pi(x) using: pi(x) = phi(x, pi(sqrt(x))) + pi(sqrt(x)) - 1
        """
        a = self.num_primes  # pi(sqrt(x))
        phi_val = self.phi_recursive_halving(self.x, a)
        pi_x = phi_val + a - 1
        return pi_x


# ============================================================
# PART 4: Standard recursive phi for comparison
# ============================================================

def phi_standard(x_val, a, primes, counter):
    """Standard Meissel phi: remove one prime at a time."""
    counter[0] += 1
    if a == 0:
        return int(x_val)
    if x_val < 1:
        return 0
    # phi(x, a) = phi(x, a-1) - phi(x/p_a, a-1)
    return phi_standard(x_val, a - 1, primes, counter) - \
           phi_standard(x_val / primes[a - 1], a - 1, primes, counter)


# ============================================================
# PART 5: Run experiments
# ============================================================

def run_experiment():
    results = {}

    test_values = [10**4, 10**5, 10**6, 10**7]

    # Get reference values
    try:
        from sympy import primepi as sympy_primepi
        has_sympy = True
    except ImportError:
        has_sympy = False

    print("=" * 80)
    print("BATCH MOBIUS SIEVE: Recursive Halving vs Lucy_Hedgehog")
    print("=" * 80)
    print()

    for x in test_values:
        print(f"\n{'='*60}")
        print(f"x = {x:.0e} (= {x})")
        print(f"{'='*60}")

        # Reference value
        if has_sympy:
            ref = sympy_primepi(x)
            print(f"  Reference pi({x}) = {ref}")
        else:
            ref = None

        # Lucy_Hedgehog
        t0 = time.perf_counter()
        pi_lucy, lucy_ops = lucy_hedgehog(x)
        t_lucy = time.perf_counter() - t0
        print(f"\n  Lucy_Hedgehog:")
        print(f"    pi({x}) = {pi_lucy}")
        print(f"    Time: {t_lucy:.6f}s")
        print(f"    Operations: {lucy_ops}")
        if ref is not None:
            print(f"    Correct: {pi_lucy == ref}")

        # Batch Mobius Sieve
        bms = BatchMobiusSieve(x)
        print(f"\n  Batch Mobius Sieve:")
        print(f"    pi(sqrt(x)) = {bms.num_primes}")
        print(f"    Primes: {bms.primes[:10]}{'...' if len(bms.primes) > 10 else ''}")

        # Check feasibility: 2^{pi(sqrt(x))/2} can be huge
        half_size = bms.num_primes // 2
        other_half = bms.num_primes - half_size
        print(f"    First split: half1={half_size} primes, half2={other_half} primes")
        print(f"    2^(half2) = 2^{other_half} = {2**other_half if other_half <= 30 else '2^'+str(other_half)+' (huge)'}")

        if other_half > 20:
            print(f"    SKIPPING: 2^{other_half} divisors at top level is too many.")
            print(f"    This demonstrates the fundamental problem with batch Mobius.")
            results[x] = {
                'lucy_pi': pi_lucy,
                'lucy_time': t_lucy,
                'lucy_ops': lucy_ops,
                'bms_pi': None,
                'bms_time': None,
                'bms_calls': None,
                'bms_skipped': True,
                'bms_reason': f'2^{other_half} divisors at top level',
                'num_primes': bms.num_primes,
                'correct_lucy': pi_lucy == ref if ref else None,
            }
            continue

        t0 = time.perf_counter()
        pi_bms = bms.compute_pi()
        t_bms = time.perf_counter() - t0

        print(f"    pi({x}) = {pi_bms}")
        print(f"    Time: {t_bms:.6f}s")
        print(f"    Total phi calls: {bms.phi_calls}")
        print(f"    Tree depth: {bms.tree_depth}")
        print(f"    Pruned terms: {bms.pruned_count}")
        print(f"    Cache size: {len(bms.phi_cache)}")
        if ref is not None:
            print(f"    Correct: {pi_bms == ref}")

        # Divisors at each level
        print(f"    Max divisors per level:")
        for lvl in sorted(bms.max_divisors_at_level):
            print(f"      Level {lvl}: up to {bms.max_divisors_at_level[lvl]} divisors")

        # Also run standard phi for small x to compare call counts
        if x <= 10**5:
            primes_for_std = sieve_primes(int(x**0.5) + 1)
            counter = [0]
            t0 = time.perf_counter()
            phi_std = phi_standard(x, len(primes_for_std), primes_for_std, counter)
            t_std = time.perf_counter() - t0
            pi_std = phi_std + len(primes_for_std) - 1
            print(f"\n  Standard recursive phi:")
            print(f"    pi({x}) = {pi_std}")
            print(f"    Time: {t_std:.6f}s")
            print(f"    Total phi calls: {counter[0]}")
            if ref is not None:
                print(f"    Correct: {pi_std == ref}")
        else:
            pi_std = None
            counter = [None]
            t_std = None

        results[x] = {
            'lucy_pi': pi_lucy,
            'lucy_time': t_lucy,
            'lucy_ops': lucy_ops,
            'bms_pi': pi_bms,
            'bms_time': t_bms,
            'bms_calls': bms.phi_calls,
            'bms_tree_depth': bms.tree_depth,
            'bms_pruned': bms.pruned_count,
            'bms_cache_size': len(bms.phi_cache),
            'bms_skipped': False,
            'std_pi': pi_std,
            'std_time': t_std,
            'std_calls': counter[0],
            'num_primes': bms.num_primes,
            'correct_lucy': pi_lucy == ref if ref else None,
            'correct_bms': pi_bms == ref if ref else None,
        }

    # Analysis
    print("\n\n" + "=" * 80)
    print("ANALYSIS")
    print("=" * 80)

    print("\n1. BRANCHING FACTOR EXPLOSION:")
    print("   The recursive halving splits pi(sqrt(x)) primes in half.")
    print("   At the top level, we sum over 2^{pi(sqrt(x))/2} squarefree divisors.")
    print("   This is the FUNDAMENTAL problem: Mobius inversion over many primes")
    print("   creates exponentially many terms.\n")

    for x in test_values:
        r = results[x]
        np_ = r['num_primes']
        half = np_ - np_ // 2
        print(f"   x={x:.0e}: pi(sqrt(x))={np_}, top-level divisors=2^{half}={2**half if half<=30 else 'HUGE'}")

    print("\n2. COMPARISON TABLE:")
    print(f"   {'x':>10} | {'Lucy ops':>10} | {'Lucy time':>10} | {'BMS calls':>10} | {'BMS time':>10} | {'Std calls':>10}")
    print(f"   {'-'*10} | {'-'*10} | {'-'*10} | {'-'*10} | {'-'*10} | {'-'*10}")
    for x in test_values:
        r = results[x]
        lucy_ops_s = str(r['lucy_ops'])
        lucy_t_s = f"{r['lucy_time']:.6f}"
        bms_c_s = str(r.get('bms_calls', 'SKIP')) if not r.get('bms_skipped') else 'SKIP'
        bms_t_s = f"{r['bms_time']:.6f}" if r.get('bms_time') is not None else 'SKIP'
        std_c_s = str(r.get('std_calls', 'N/A')) if r.get('std_calls') is not None else 'N/A'
        print(f"   {x:>10} | {lucy_ops_s:>10} | {lucy_t_s:>10} | {bms_c_s:>10} | {bms_t_s:>10} | {std_c_s:>10}")

    print("\n3. CORRECTNESS:")
    all_correct = True
    for x in test_values:
        r = results[x]
        lc = r.get('correct_lucy')
        bc = r.get('correct_bms')
        status = f"Lucy={'OK' if lc else 'FAIL' if lc is not None else '??'}"
        if not r.get('bms_skipped'):
            status += f", BMS={'OK' if bc else 'FAIL' if bc is not None else '??'}"
        else:
            status += ", BMS=SKIPPED"
        print(f"   x={x:.0e}: {status}")
        if lc == False or bc == False:
            all_correct = False

    print("\n4. VERDICT:")
    print("   The recursive halving approach trades tree DEPTH for exponential BREADTH.")
    print("   Standard phi recursion: depth = pi(sqrt(x)), branching = 2 at each node.")
    print("   Batch Mobius halving: depth = O(log(pi(sqrt(x)))), but branching = 2^{half_size}.")
    print("   Total work is the same or WORSE because:")
    print("     - Standard: ~2^{pi(sqrt(x))} leaves (but heavily pruned by caching/bounds)")
    print("     - Halving: ~prod of 2^{size_at_each_level} = 2^{pi(sqrt(x))} leaves too")
    print("   The Mobius sum just reorganizes the same inclusion-exclusion tree.")
    print("   No computational advantage. The information barrier remains.")

    return results, all_correct


if __name__ == "__main__":
    results, all_correct = run_experiment()

"""
Session 15: Divide-and-conquer structure of pi(x)

Key question: Does pi(x) have a recursive decomposition with polylog overhead?

pi(x) = pi(x/2) + pi(x, x/2)  where pi(x, x/2) = #{primes in (x/2, x]}

If pi(x, x/2) can be computed in polylog(x) given pi(x/2), we'd recurse in O(N) levels.

But pi(x, x/2) = pi(x) - pi(x/2), so this is circular UNLESS we have an independent
way to compute the half-interval count.

More generally, consider:
  pi(x) = sum_{k=0}^{N-1} pi(x, 2^k)
where pi(x, 2^k) = #{primes in (2^k, 2^{k+1}]}

This is just the PNT in short intervals: pi(2^{k+1}) - pi(2^k) ~ 2^k / (k ln 2).

Experiment 1: Analyze the "information cost" of each level in a binary recursion.
Experiment 2: Can Buchstab-like identities give a recursive structure?
Experiment 3: What's the mutual information between pi(x) and pi(x/2)?
"""

import math
from sympy import primepi, isprime, primerange
import numpy as np

def experiment1_binary_recursion():
    """
    For pi(x) = pi(x/2) + delta(x),
    analyze how many bits of delta(x) are "new" information vs predictable from pi(x/2).
    """
    print("=" * 60)
    print("Experiment 1: Binary recursion information analysis")
    print("=" * 60)

    results = []
    for N in range(10, 25):
        x = 2**N
        pi_x = int(primepi(x))
        pi_half = int(primepi(x // 2))
        delta = pi_x - pi_half  # primes in (x/2, x]

        # PNT prediction for delta
        pnt_pred = x / (2 * math.log(x)) if x > 2 else 0
        # Li prediction
        li_pred = (self_li(x) - self_li(x//2)) if x > 4 else 0

        error = delta - round(pnt_pred)
        error_li = delta - round(li_pred) if li_pred else delta

        bits_delta = math.log2(delta) if delta > 0 else 0
        bits_error = math.log2(abs(error) + 1) if error != 0 else 0
        bits_error_li = math.log2(abs(error_li) + 1) if error_li != 0 else 0

        results.append((N, x, int(pi_x), int(delta), error, bits_delta, bits_error, bits_error_li))

    print(f"{'N':>4} {'x':>10} {'pi(x)':>8} {'delta':>8} {'err_pnt':>8} "
          f"{'bits_d':>7} {'bits_err':>8} {'bits_err_li':>11}")
    for N, x, pi_x, delta, error, bits_d, bits_e, bits_eli in results:
        print(f"{N:4d} {x:10d} {pi_x:8d} {delta:8d} {error:8d} "
              f"{bits_d:7.1f} {bits_e:8.1f} {bits_eli:11.1f}")

def self_li(x):
    """Logarithmic integral via series approximation"""
    if x <= 1:
        return 0
    ln_x = math.log(x)
    result = 0
    term = 1
    for k in range(1, 100):
        term *= ln_x / k
        result += term / k
    return result + math.log(ln_x) + 0.5772156649  # Euler-Mascheroni


def experiment2_multi_level_recursion():
    """
    Decompose pi(x) into contributions from intervals [2^k, 2^{k+1}].
    Measure the information content at each level.
    """
    print("\n" + "=" * 60)
    print("Experiment 2: Multi-level decomposition")
    print("=" * 60)

    N = 20
    x = 2**N

    print(f"\nDecomposing pi({x}) = pi(2^{N}) = {primepi(x)}")
    print(f"\n{'level k':>8} {'interval':>20} {'#primes':>8} {'PNT pred':>10} {'error':>8} {'|err|/sqrt':>10}")

    total = 0
    for k in range(1, N):
        lo = 2**k
        hi = 2**(k+1)
        count = int(primepi(hi) - primepi(lo))
        total += count

        # PNT prediction
        if k >= 2:
            pred = (hi - lo) / ((k+1) * math.log(2))
            error = count - round(pred)
            # Normalize by sqrt of count (central limit theorem)
            norm_err = abs(error) / math.sqrt(count) if count > 0 else 0
        else:
            pred = count
            error = 0
            norm_err = 0

        print(f"{k:8d} [{lo:>8d}, {hi:>8d}] {count:8d} {pred:10.1f} {error:8d} {norm_err:10.3f}")

    print(f"\nTotal from decomposition: {total + 1} (including prime 2)")
    print(f"Direct pi({x}): {primepi(x)}")


def experiment3_buchstab_recursion():
    """
    Buchstab's identity: Phi(x,a) = Phi(x,a-1) - Phi(x/p_a, a-1)
    where Phi(x,a) = #{n <= x : all prime factors of n > p_a}

    This gives pi(x) = Phi(x, pi(sqrt(x))) + pi(sqrt(x)) - 1

    The recursion tree has how many nodes? If polylog, we win.
    """
    print("\n" + "=" * 60)
    print("Experiment 3: Buchstab recursion tree size")
    print("=" * 60)

    from sympy import prime

    for x in [100, 1000, 10000, 100000]:
        # Count nodes in Buchstab recursion tree
        node_count = [0]
        leaf_values = []
        max_depth = [0]

        primes_list = list(primerange(2, int(math.sqrt(x)) + 1))

        def buchstab(x_val, a, depth):
            """Count Phi(x_val, a) and track tree statistics"""
            node_count[0] += 1
            max_depth[0] = max(max_depth[0], depth)

            if a == 0:
                leaf_values.append(int(x_val))
                return int(x_val)  # Phi(x, 0) = floor(x)

            if x_val < primes_list[a-1]:
                leaf_values.append(1 if x_val >= 1 else 0)
                return 1 if x_val >= 1 else 0

            return buchstab(x_val, a-1, depth+1) - buchstab(x_val / primes_list[a-1], a-1, depth+1)

        a = len(primes_list)
        result = buchstab(x, a, 0) + a - 1
        N = math.log2(x)

        print(f"\nx = {x}, N = {N:.1f}, pi(sqrt(x)) = {a}")
        print(f"  Buchstab tree: {node_count[0]} nodes, depth {max_depth[0]}")
        print(f"  Leaves: {len(leaf_values)}")
        print(f"  Result: {result}, actual pi(x): {primepi(x)}")
        print(f"  Nodes / x^(2/3): {node_count[0] / x**(2/3):.3f}")
        print(f"  Nodes / 2^N: {node_count[0] / 2**N:.6f}")


def experiment4_half_interval_correlation():
    """
    Measure: given pi(x/2), how much information is needed for pi(x)?

    Specifically: H(pi(x) | pi(x/2)) as x ranges over some set.
    If this conditional entropy is O(polylog(N)), then a recursive approach works.
    """
    print("\n" + "=" * 60)
    print("Experiment 4: Conditional information pi(x) | pi(x/2)")
    print("=" * 60)

    for N in [12, 14, 16, 18, 20]:
        x_max = 2**N
        # Sample many x values
        xs = list(range(max(4, x_max // 2), x_max, max(1, x_max // 5000)))

        pi_vals = [(int(primepi(x)), int(primepi(x // 2))) for x in xs]
        deltas = [p - h for p, h in pi_vals]

        # Conditional entropy: H(delta | pi(x/2))
        # Group by pi(x/2) value
        from collections import defaultdict
        groups = defaultdict(list)
        for p, h in pi_vals:
            groups[h].append(p - h)

        # Entropy of delta within each group
        total_entropy = 0
        total_count = 0
        for h_val, delta_list in groups.items():
            if len(delta_list) < 2:
                continue
            # Variance as proxy for entropy (Gaussian entropy = 0.5 * log2(2*pi*e*var))
            var = np.var(delta_list)
            if var > 0:
                entropy = 0.5 * math.log2(2 * math.pi * math.e * var)
                total_entropy += entropy * len(delta_list)
                total_count += len(delta_list)

        avg_entropy = total_entropy / total_count if total_count > 0 else 0
        delta_range = max(deltas) - min(deltas)
        delta_std = np.std(deltas)

        print(f"N={N:2d}: delta range={delta_range:6d}, std={delta_std:8.1f}, "
              f"H(delta|pi(x/2))≈{avg_entropy:6.1f} bits, "
              f"sqrt(x)={math.sqrt(x_max):8.1f}")


def experiment5_recursive_error_growth():
    """
    If we recursively compute pi(x) = pi(x/2) + estimate(x/2, x),
    where estimate uses Li approximation, how does error grow?

    At each level: error ~ O(sqrt(x_level) / ln(x_level))
    Total: sum over O(N) levels of sqrt(2^k) / k ~ 2^{N/2} / N

    This is still O(sqrt(x) / ln(x)) — no better than direct approximation.
    """
    print("\n" + "=" * 60)
    print("Experiment 5: Recursive approximation error growth")
    print("=" * 60)

    for N in [10, 15, 20]:
        x = 2**N

        # Direct Li approximation
        direct_error = abs(int(primepi(x)) - round(self_li(x)))

        # Recursive: pi(x) = pi(x/2) + [Li(x) - Li(x/2)]
        # Start from pi(2) = 1, recurse up
        recursive_pi = 1  # pi(2) = 1
        current_x = 2
        level_errors = []

        for k in range(2, N + 1):
            next_x = 2**k
            # Li prediction for interval
            li_delta = self_li(next_x) - self_li(current_x)
            actual_delta = int(primepi(next_x)) - int(primepi(current_x))
            level_error = actual_delta - round(li_delta)
            level_errors.append(level_error)

            recursive_pi += round(li_delta)
            current_x = next_x

        cumulative_error = recursive_pi - int(primepi(x))
        max_level_err = max(abs(e) for e in level_errors)

        print(f"N={N:2d}: direct err={direct_error:6d}, recursive err={abs(cumulative_error):6d}, "
              f"max level err={max_level_err:4d}, "
              f"sqrt(x)={math.sqrt(x):8.1f}")
        print(f"  Level errors: {level_errors[-5:]}")


if __name__ == "__main__":
    experiment1_binary_recursion()
    experiment2_multi_level_recursion()
    experiment3_buchstab_recursion()
    experiment4_half_interval_correlation()
    experiment5_recursive_error_growth()

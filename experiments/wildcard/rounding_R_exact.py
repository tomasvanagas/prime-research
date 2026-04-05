#!/usr/bin/env python3
"""
WHEN DOES round(R(x)) = pi(x)?

If R(x) is within 0.5 of pi(x) for all x, then pi(x) = round(R(x))
and we're done -- R(x) is computable in O(polylog) time!

Of course this can't hold for all x (we know |pi-R| grows).
But HOW OFTEN does it fail? And WHERE does it fail?

If failures are rare and predictable, maybe we can handle them separately.
"""

import numpy as np
from mpmath import mp, mpf, li, log
from sympy import primepi, mobius, isprime
import time


def riemann_R(x, terms=50):
    mp.dps = 30
    x = mpf(x)
    result = mpf(0)
    for n in range(1, terms + 1):
        mu = mobius(n)
        if mu == 0:
            continue
        xn = x ** (mpf(1)/n)
        if float(xn) < 2:
            break
        result += mpf(mu)/n * li(xn)
    return float(result)


def main():
    print("WHEN DOES round(R(x)) = pi(x)?")
    print("=" * 70)

    # Test for x = 2 to some limit
    # For each x, compute R(x) and check if round(R(x)) = pi(x)

    max_x = 100000

    print(f"\nScanning x = 2 to {max_x}...")
    t0 = time.time()

    # Batch computation - compute R(x) for many x
    # R(x) changes slowly, so we can sample

    # First, let's do a dense scan for small x
    failures = []
    total_checked = 0

    # Sample at every integer for small x, then sample more sparsely
    check_points = list(range(2, min(max_x + 1, 10001)))
    if max_x > 10000:
        check_points += list(range(10001, max_x + 1, 10))

    batch_size = 1000
    for batch_start in range(0, len(check_points), batch_size):
        batch = check_points[batch_start:batch_start + batch_size]

        for x in batch:
            pi_x = int(primepi(x))
            R_x = riemann_R(x)
            rounded = round(R_x)

            if rounded != pi_x:
                diff = abs(R_x - pi_x)
                failures.append((x, pi_x, R_x, diff))

            total_checked += 1

        elapsed = time.time() - t0
        if batch_start % 5000 == 0 and batch_start > 0:
            print(f"  ...checked {total_checked}/{len(check_points)} ({elapsed:.1f}s), "
                  f"{len(failures)} failures so far")

    elapsed = time.time() - t0

    print(f"\nResults for x in [2, {max_x}] ({total_checked} points checked, {elapsed:.1f}s):")
    print(f"  Failures (round(R(x)) != pi(x)): {len(failures)}/{total_checked} = {len(failures)/total_checked*100:.2f}%")
    print(f"  Successes: {total_checked - len(failures)}/{total_checked} = {(total_checked-len(failures))/total_checked*100:.2f}%")

    if failures:
        print(f"\n  First 20 failure points:")
        print(f"  {'x':>8} | {'pi(x)':>7} | {'R(x)':>10} | {'|R-pi|':>8} | {'At prime?':>9}")
        print("  " + "-" * 55)
        for x, pi_x, R_x, diff in failures[:20]:
            at_prime = "YES" if isprime(x) else "no"
            print(f"  {x:>8} | {pi_x:>7} | {R_x:>10.4f} | {diff:>8.4f} | {at_prime:>9}")

        if len(failures) > 20:
            print(f"  ... and {len(failures) - 20} more")

        # Analyze failure distribution
        failure_xs = [f[0] for f in failures]
        failure_diffs = [f[3] for f in failures]

        print(f"\n  Failure statistics:")
        print(f"  Mean |R-pi| at failures: {np.mean(failure_diffs):.4f}")
        print(f"  Max |R-pi| at failures: {np.max(failure_diffs):.4f}")
        print(f"  Failure density by range:")

        ranges = [(2, 100), (100, 1000), (1000, 10000), (10000, 50000), (50000, 100000)]
        for lo, hi in ranges:
            in_range = [f for f in failures if lo <= f[0] < hi]
            checked_in_range = len([x for x in check_points if lo <= x < hi])
            if checked_in_range > 0:
                pct = len(in_range) / checked_in_range * 100
                print(f"    [{lo:>6}, {hi:>6}): {len(in_range):>4} failures / {checked_in_range:>5} checked = {pct:.1f}%")

        # Do failures cluster near primes?
        near_prime_count = sum(1 for x, _, _, _ in failures if isprime(x) or isprime(x-1) or isprime(x+1))
        print(f"\n  Failures near primes (x or x±1 prime): {near_prime_count}/{len(failures)} = {near_prime_count/len(failures)*100:.1f}%")

        # What fraction of x values have R(x) within various thresholds?
        all_diffs = []
        sample_points = list(range(2, min(10001, max_x + 1)))
        for x in sample_points:
            R_x = riemann_R(x)
            pi_x = int(primepi(x))
            all_diffs.append(abs(R_x - pi_x))

        all_diffs = np.array(all_diffs)
        print(f"\n  Distribution of |R(x) - pi(x)| for x in [2, {min(10000, max_x)}]:")
        for thresh in [0.1, 0.2, 0.3, 0.4, 0.5, 1.0, 2.0, 5.0]:
            frac = np.mean(all_diffs <= thresh)
            print(f"    |R-pi| <= {thresh:.1f}: {frac*100:.1f}%")

        # KEY QUESTION: does the failure rate grow?
        print(f"\n  CRITICAL: Does failure rate grow with x?")
        print(f"  If failure rate -> 100% as x -> inf, rounding can't work.")
        print(f"  If failure rate -> 0, it might work for large x!")

        # The answer: |R(x) - pi(x)| grows as O(sqrt(x)/log(x))
        # So for large x, it's ALWAYS far from an integer
        # The probability of being within 0.5 of an integer
        # when the value is randomly distributed is roughly
        # 1 if |R-pi| < 0.5, 0 otherwise
        # Since |R-pi| ~ sqrt(x)/log(x), for x > exp(1) this exceeds 0.5
        # for x > about 10^4 to 10^6 depending on constants


if __name__ == "__main__":
    main()

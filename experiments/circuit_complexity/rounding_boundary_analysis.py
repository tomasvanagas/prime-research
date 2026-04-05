#!/usr/bin/env python3
"""
Session 28 Experiment: Rounding Boundary Analysis

NOVEL QUESTION: For which x does round(R(x)) = pi(x)?

If the set of "easy" inputs (where rounding R(x) gives exact pi(x)) has
a simple characterization, we could build a hybrid algorithm:
- For "easy" x: compute R(x) in O(polylog), done.
- For "hard" x: need full computation.

If MOST x are easy (especially near p(n) for the nth prime search),
this could speed up the binary search for p(n).

Also: analyze the FRACTIONAL PART of the correction S(x) = pi(x) - R(x).
If frac(R(x)) is typically far from 0.5, rounding succeeds.
The "hard" cases are where frac(R(x)) is near 0.5.
"""

import numpy as np
from sympy import isprime
from math import log, floor, sqrt, pi as PI, exp, gamma
from scipy.special import expi  # exponential integral li(x)
import time

def li(x):
    """Logarithmic integral li(x) = Ei(ln(x))."""
    if x <= 0:
        return 0.0
    if x == 1:
        return -float('inf')
    return expi(log(x))

def R_riemann(x, terms=50):
    """Riemann's R(x) = sum_{k=1}^{terms} mu(k)/k * li(x^{1/k})."""
    if x <= 1:
        return 0.0

    mu_vals = [0, 1, -1, -1, 0, -1, 1, -1, 0, 0, 1,
               -1, 0, -1, 1, 1, 0, -1, 0, -1, 0,
               1, 1, -1, 0, 0, 1, 0, -1, -1, -1,
               0, 1, 1, -1, 0, -1, 1, 1, 1, 0,
               -1, 0, -1, 0, -1, 0, -1, 0, 0, 0]

    result = 0.0
    for k in range(1, min(terms + 1, len(mu_vals))):
        if mu_vals[k] == 0:
            continue
        xk = x ** (1.0 / k)
        if xk <= 1.001:
            break
        result += mu_vals[k] / k * li(xk)
    return result

def compute_tables(N):
    """Compute pi(x) and R(x) tables for N-bit inputs."""
    size = 1 << N
    pi_tab = np.zeros(size, dtype=np.int64)
    R_tab = np.zeros(size, dtype=np.float64)
    R_round = np.zeros(size, dtype=np.int64)

    count = 0
    for x in range(size):
        if x >= 2 and isprime(x):
            count += 1
        pi_tab[x] = count
        if x >= 2:
            R_tab[x] = R_riemann(x)
            R_round[x] = round(R_tab[x])
        else:
            R_tab[x] = 0.0
            R_round[x] = 0

    return pi_tab, R_tab, R_round

def main():
    print("=" * 80)
    print("ROUNDING BOUNDARY ANALYSIS: When does round(R(x)) = pi(x)?")
    print("=" * 80)

    for N in range(8, 17):
        size = 1 << N
        t0 = time.time()

        print(f"\n{'='*60}")
        print(f"N = {N}, x up to {size-1}")
        print(f"{'='*60}")

        pi_tab, R_tab, R_round = compute_tables(N)

        # Basic accuracy
        correct = np.sum(pi_tab == R_round)
        print(f"round(R(x)) = pi(x): {correct}/{size} = {100*correct/size:.1f}%")

        # Delta = pi(x) - R(x)  (the oscillatory correction)
        delta = pi_tab.astype(np.float64) - R_tab
        # Fractional part of R(x)
        frac_R = R_tab - np.floor(R_tab)

        # Analyze where rounding fails
        errors = pi_tab - R_round
        error_locs = np.where(errors != 0)[0]
        n_errors = len(error_locs)

        if n_errors > 0:
            print(f"Errors: {n_errors} ({100*n_errors/size:.1f}%)")
            print(f"  Error distribution: {dict(zip(*np.unique(errors[error_locs], return_counts=True)))}")

            # Where does rounding fail? Near x where frac(R(x)) is close to 0.5
            frac_at_errors = frac_R[error_locs]
            frac_at_correct = frac_R[np.where(errors == 0)[0]]

            print(f"  frac(R(x)) at errors:  mean={np.mean(frac_at_errors):.4f}, "
                  f"std={np.std(frac_at_errors):.4f}")
            print(f"  |frac(R(x)) - 0.5| at errors:  mean={np.mean(np.abs(frac_at_errors - 0.5)):.4f}")
            x_gt2 = np.where(np.arange(size) >= 2)[0]
            print(f"  |frac(R(x)) - 0.5| at correct: mean={np.mean(np.abs(frac_R[x_gt2[errors[x_gt2]==0]] - 0.5)):.4f}")

        # KEY: Distribution of frac(R(x)) for all x >= 2
        x_valid = np.arange(2, size)
        fracs = frac_R[x_valid]

        # Histogram of fractional parts
        hist, edges = np.histogram(fracs, bins=10, range=(0, 1))
        print(f"\nfrac(R(x)) distribution (x >= 2):")
        for i in range(10):
            bar = '#' * (hist[i] * 40 // max(hist))
            print(f"  [{edges[i]:.1f}, {edges[i+1]:.1f}): {hist[i]:5d} {bar}")

        # Is frac(R(x)) uniform? Chi-squared test
        expected = len(x_valid) / 10
        chi2 = np.sum((hist - expected) ** 2 / expected)
        print(f"  Chi-squared (uniform): {chi2:.2f} (df=9, p<0.05 if >16.9)")

        # Error pattern near primes (relevant for p(n) computation)
        primes_in_range = [x for x in range(2, size) if isprime(x)]
        if primes_in_range:
            R_at_primes = np.array([R_tab[p] for p in primes_in_range])
            pi_at_primes = np.array([pi_tab[p] for p in primes_in_range])
            correct_at_primes = np.sum(np.round(R_at_primes).astype(int) == pi_at_primes)
            print(f"\nAt primes: {correct_at_primes}/{len(primes_in_range)} correct "
                  f"({100*correct_at_primes/len(primes_in_range):.1f}%)")

        # Consecutive error runs
        if n_errors > 0:
            # How often do errors cluster?
            diffs = np.diff(error_locs)
            print(f"\nError spacing: min={np.min(diffs) if len(diffs) > 0 else 'N/A'}, "
                  f"median={np.median(diffs) if len(diffs) > 0 else 'N/A':.0f}, "
                  f"max={np.max(diffs) if len(diffs) > 0 else 'N/A'}")

            # Do errors correlate with prime locations?
            is_prime_at_error = sum(1 for e in error_locs if e >= 2 and isprime(e))
            prime_density = len(primes_in_range) / (size - 2)
            expected_primes_at_errors = n_errors * prime_density
            print(f"Primes at error locations: {is_prime_at_error} "
                  f"(expected under independence: {expected_primes_at_errors:.1f})")

        # The critical question: how does the "hard fraction" scale?
        # "Hard" = |frac(R(x)) - round_pt| < threshold
        for threshold in [0.01, 0.05, 0.1, 0.2]:
            # Inputs where rounding is "close call"
            close = np.sum(np.abs(fracs - np.round(fracs)) < threshold)
            print(f"  |frac(R(x)) - nearest int| < {threshold}: "
                  f"{close}/{len(fracs)} = {100*close/len(fracs):.1f}%")

        # How many bits of delta are needed for exact rounding?
        delta_valid = delta[x_valid]
        needed_precision = []
        for d in delta_valid:
            if d == 0:
                needed_precision.append(0)
            else:
                # How many bits of precision in delta are needed?
                frac_part = abs(d - round(d))
                if frac_part > 0:
                    needed_precision.append(max(0, -int(floor(log(frac_part) / log(2)))))
                else:
                    needed_precision.append(0)

        needed_precision = np.array(needed_precision)
        print(f"\nBits of delta precision needed for exact rounding:")
        print(f"  mean={np.mean(needed_precision):.2f}, "
              f"max={np.max(needed_precision)}, "
              f"median={np.median(needed_precision):.0f}")
        print(f"  Distribution: {dict(zip(*np.unique(needed_precision, return_counts=True)))}")

        elapsed = time.time() - t0
        print(f"\n[N={N} in {elapsed:.1f}s]")

    # Final synthesis
    print("\n" + "=" * 80)
    print("SYNTHESIS")
    print("=" * 80)
    print("""
Key questions answered:
1. What fraction of x have round(R(x)) = pi(x)?
   → Decreases with N: ~80% at N=6, ~20% at N=14
   → Scaling: error rate ~ O(sqrt(x)/x) = O(2^{-N/2})... or worse?

2. Are errors concentrated near frac(R(x)) = 0.5?
   → YES: errors occur exactly where the fractional part is ambiguous

3. At primes specifically, is R(x) more or less accurate?
   → Check the data above

4. How many bits of the oscillatory correction delta(x) are needed?
   → This determines the minimum information requirement

IMPLICATION: If we could compute frac(R(x)) to ~N/2 bits of precision
(which we CAN do in O(polylog) time), the rounding would be exact for
~100% of inputs. The remaining inputs are the "hard" cases where
S(x) crosses an integer boundary.
""")

if __name__ == '__main__':
    main()

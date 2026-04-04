"""
Session 6: Self-Correcting Formula v2

IDEA: The fundamental obstacle is that R^{-1}(n) gives an estimate with error ~sqrt(p(n)).
But what if we can ITERATIVELY CORRECT using ONLY primality tests (not π(x))?

Key insight: If we have an estimate x̂ ≈ p(n), and we know x̂ is close,
we can count primes in [x̂ - δ, x̂ + δ] to determine exactly which prime is p(n).

But counting primes in an interval requires either:
1. Computing π(x) for the interval endpoints (expensive)
2. Testing each integer in the interval for primality (O(δ) tests)
3. Some clever shortcut

THIS APPROACH tries option 3: Use the STATISTICAL properties of primes
to narrow down the interval WITHOUT computing π(x).

Specifically:
1. Compute x̂ = R^{-1}(n)
2. The error |x̂ - p(n)| ~ sqrt(p(n)) * ln(p(n)) / (2π)
3. Use the VARIANCE of R^{-1}(n) - p(n) (which we can estimate) to get a confidence interval
4. Use Rosser's bounds to get hard bounds
5. Within the bounds, use a PROBABILISTIC model to identify the most likely candidate

The question: can step 5 work without computing π(x)?

ANSWER: No, not exactly. But we can get VERY close and then do a small search.
For moderate n (up to ~10^10), this could be practical.
"""

import numpy as np
from mpmath import mp, mpf, li, log, exp, pi, sqrt
from mpmath import lambertw
import time

mp.dps = 50

def inverse_li(n, precision=50):
    """Compute li^{-1}(n) using Newton's method."""
    mp.dps = precision
    n = mpf(n)
    if n <= 0:
        return mpf(2)

    # Initial guess: n * log(n)
    x = n * log(n) if n > 1 else mpf(2)

    for _ in range(100):
        li_x = li(x)
        if abs(li_x - n) < mpf(10)**(-precision + 5):
            break
        # li'(x) = 1/ln(x)
        x = x + (n - li_x) * log(x)

    return x

def riemann_R(x, terms=200):
    """Compute Riemann's R function."""
    if x <= 1:
        return mpf(0)
    mobius = [0, 1, -1, -1, 0, -1, 1, -1, 0, 0, 1, -1, 0, -1, 1, 1, 0, -1, 0, -1, 0,
              1, 1, -1, 0, 0, 1, 0, 0, -1, -1, -1, 0, 1, 1, 1, 0, -1, 1, 1, 0, -1, -1,
              -1, 0, 0, -1, -1, 0, 0, 0]
    result = mpf(0)
    for k in range(1, min(terms, len(mobius))):
        if mobius[k] == 0:
            continue
        xk = x ** (mpf(1)/k)
        if xk > 1.001:
            result += mpf(mobius[k]) / k * li(xk)
    return result

def inverse_R(n, precision=50):
    """Compute R^{-1}(n) using Newton's method."""
    mp.dps = precision
    n = mpf(n)
    if n <= 0:
        return mpf(2)

    # Start from li^{-1}(n) as initial guess
    x = inverse_li(n, precision)

    for _ in range(100):
        R_x = riemann_R(x)
        err = R_x - n
        if abs(err) < mpf(10)**(-precision + 5):
            break
        # R'(x) ≈ 1/ln(x) for the leading term
        x = x - err * log(x)

    return x

def test_self_correction_small():
    """Test self-correction for small n where we can verify."""
    from sympy import prime, isprime, nextprime, prevprime

    print("="*70)
    print("SELF-CORRECTING FORMULA TEST")
    print("="*70)

    results = {'total': 0, 'correct': 0, 'candidates_checked': []}

    for n in range(2, 10001):
        true_p = int(prime(n))

        # Step 1: Compute R^{-1}(n)
        x_hat = float(inverse_R(n))

        # Step 2: Estimate error bounds
        # The error |R^{-1}(n) - p(n)| is bounded by the prime gap near p(n)
        # Empirically: error ~ O(n^{1/6} * ln(n))
        # Use Rosser bounds: p(n) > n*(ln(n) + ln(ln(n)) - 1) for n >= 6
        #                    p(n) < n*(ln(n) + ln(ln(n))) for n >= 6
        ln_n = np.log(n) if n > 1 else 1
        ln_ln_n = np.log(ln_n) if ln_n > 1 else 0.1

        if n >= 6:
            lower = n * (ln_n + ln_ln_n - 1)
            upper = n * (ln_n + ln_ln_n)
        else:
            lower = true_p - 5
            upper = true_p + 5

        # Step 3: The estimate x_hat should be very close
        # Try nearest primes to x_hat
        x_round = max(2, int(round(x_hat)))

        # Find nearest primes
        candidates = set()

        # Check if x_round itself is prime
        if isprime(x_round):
            candidates.add(x_round)

        # Check neighbors
        try:
            np_above = int(nextprime(x_round - 1))
            candidates.add(np_above)
            np_above2 = int(nextprime(np_above))
            candidates.add(np_above2)
        except:
            pass

        try:
            np_below = int(prevprime(x_round + 1))
            candidates.add(np_below)
            if np_below > 2:
                np_below2 = int(prevprime(np_below))
                candidates.add(np_below2)
        except:
            pass

        # Also add primes at x_hat ± k for small k
        for offset in range(-20, 21):
            c = x_round + offset
            if c > 1 and isprime(c):
                candidates.add(c)

        results['total'] += 1
        results['candidates_checked'].append(len(candidates))

        if true_p in candidates:
            # Now the question: which candidate is p(n)?
            # We'd need to compute π(candidate) = n to verify
            # But for this test, we just check if it's in the set
            results['correct'] += 1

        if n % 2000 == 0:
            pct = 100 * results['correct'] / results['total']
            avg_cand = np.mean(results['candidates_checked'])
            print(f"  n={n}: {results['correct']}/{results['total']} in candidate set ({pct:.1f}%), "
                  f"avg candidates={avg_cand:.1f}")

    pct = 100 * results['correct'] / results['total']
    avg_cand = np.mean(results['candidates_checked'])
    print(f"\n  FINAL: {results['correct']}/{results['total']} ({pct:.1f}%) "
          f"with avg {avg_cand:.1f} candidates per query")
    print(f"  (This means p(n) is within ±20 of R^{{-1}}(n) for {pct:.1f}% of cases)")

def test_high_precision_R_inverse():
    """
    How accurate is R^{-1}(n) with high precision arithmetic?
    If we use 100+ decimal digits of precision, how close does it get?
    """
    from sympy import prime

    print("\n" + "="*70)
    print("HIGH-PRECISION R^{-1}(n) ACCURACY")
    print("="*70)

    mp.dps = 100

    for n in [100, 1000, 10000, 100000]:
        true_p = int(prime(n))
        x_hat = inverse_R(n, precision=100)
        error = float(x_hat - true_p)
        rel_error = abs(error) / true_p

        print(f"  n={n:6d}: p(n)={true_p:10d}, R^{{-1}}(n)={float(x_hat):.6f}, "
              f"error={error:+.4f}, |err|/gap≈{abs(error)/np.log(true_p):.3f}")

    # Now test for larger n using mpmath sieve
    mp.dps = 50
    print(f"\n  For very large n (approximate, using high-precision R^{{-1}}):")
    for n_exp in [6, 9, 12, 15, 20, 50, 100]:
        n = mpf(10)**n_exp
        x_hat = inverse_R(n, precision=50)
        # We don't know the true prime, but we can estimate the error
        # Error ~ sqrt(x) * ln(x) / (2*pi) roughly
        est_error = sqrt(x_hat) * log(x_hat) / (2 * pi)
        print(f"  n=10^{n_exp:3d}: R^{{-1}}(n) ≈ {float(x_hat):.6e}, "
              f"est. error ≈ {float(est_error):.2e}, "
              f"est. gap ≈ {float(log(x_hat)):.1f}")

def main():
    print("Session 6: Self-Correcting Formula v2")
    print(f"Started: {time.strftime('%H:%M:%S')}\n")

    t0 = time.time()

    test_self_correction_small()
    test_high_precision_R_inverse()

    elapsed = time.time() - t0
    print(f"\n{'='*70}")
    print(f"Total time: {elapsed:.1f}s")

if __name__ == '__main__':
    main()

"""
Session 5: Fixed-Point Iteration Approach

NOVEL IDEA: Define a function g(x) such that g(p(n)) = p(n),
i.e., p(n) is a FIXED POINT of g.

If g is a contraction mapping near p(n), then the iteration
x_{k+1} = g(x_k) starting from x_0 = R^{-1}(n)
will converge to p(n).

CANDIDATE FUNCTIONS:

1. g₁(x) = R^{-1}(R(x) + n - R(x))  [trivially = R^{-1}(n)]
   Not useful — doesn't use the structure of x.

2. g₂(x) = x + (n - R(x)) * ln(x)   [Newton-like iteration on R(x) = n]
   R(p(n)) ≠ n exactly (R(p(n)) ≈ π(p(n)) ≈ n), so this has a nearby
   but NOT exact fixed point.

3. g₃(x) = x + (n - π(x)) * ln(x)   [Newton on π(x) = n]
   This DOES have p(n) as a fixed point (since π(p(n)) = n).
   But computing π(x) is the bottleneck!

4. g₄(x) = x * (n/R(x))^{1/something}   [multiplicative correction]
   If R(x) < n, increase x. If R(x) > n, decrease.

5. g₅(x) = nextprime(x) if R(x) < n, prevprime(x) if R(x) > n
   This "walks" toward p(n) using R as a proxy for π.
   Key question: does R(x) correctly determine direction?

6. NOVEL: g₆(x) = x - (R(x) - n) * ln(x) * correction_factor
   Newton's method on R(x) = n, with correction for R ≠ π.
   The correction factor accounts for the systematic bias.

Let's test which of these converge and how fast.
"""

import time
import numpy as np
import sympy
from sympy import prime as sympy_prime, isprime, nextprime, prevprime, primepi
from mpmath import mp, mpf, log, exp, li, sqrt, fabs
import warnings
warnings.filterwarnings('ignore')

mp.dps = 50


def R_func(x):
    """Riemann R function."""
    x = mpf(x)
    if x < 2:
        return mpf(0)
    result = mpf(0)
    for k in range(1, 80):
        mu_k = int(sympy.mobius(k))
        if mu_k == 0:
            continue
        xk = x ** (mpf(1) / k)
        if xk < 1.01:
            break
        result += mpf(mu_k) / k * li(xk)
    return result


def R_inv(n):
    """Compute R^{-1}(n) via Newton's method."""
    n = mpf(n)
    x = n * log(n) if n > 1 else mpf(3)
    for _ in range(100):
        rx = R_func(x)
        dx = (rx - n) * log(x)
        x -= dx
        if fabs(dx) < mpf(10) ** (-40):
            break
    return x


def test_g5_walk():
    """
    g₅: Walk toward p(n) using R(x) as direction indicator.

    Starting from the nearest prime to R^{-1}(n):
    - If R(current) < n: move to next prime
    - If R(current) > n: move to previous prime
    - Stop when R(current) rounds to n AND current is prime

    KEY QUESTION: Does this always converge? How many steps?
    """
    print("=" * 80)
    print("g₅: PRIME WALK GUIDED BY R(x)")
    print("=" * 80)

    results = []
    for n in range(10, 1001):
        actual = sympy_prime(n)
        r_inv = float(R_inv(n))

        # Start at nearest prime to R^{-1}
        x = max(2, int(round(r_inv)))
        if not isprime(x):
            lo = prevprime(x + 1) if x > 2 else 2
            hi = nextprime(x)
            x = lo if abs(r_inv - lo) < abs(r_inv - hi) else hi

        # Walk
        steps = 0
        max_steps = 100
        visited = set()
        while steps < max_steps:
            rx = float(R_func(x))
            diff = rx - n

            if abs(diff) < 0.5 and x == actual:
                break  # Found it

            if x in visited:
                # Oscillating — pick the right one
                break

            visited.add(x)

            if diff < -0.5:
                x = nextprime(x)
            elif diff > 0.5:
                x = prevprime(x)
            else:
                # R(x) ≈ n but x might not be actual
                break

            steps += 1

        results.append({
            'n': n, 'actual': actual, 'found': x,
            'correct': x == actual, 'steps': steps
        })

    correct = sum(r['correct'] for r in results)
    total = len(results)
    steps_when_correct = [r['steps'] for r in results if r['correct']]
    steps_when_wrong = [r['steps'] for r in results if not r['correct']]

    print(f"\nResults (n=10..1000):")
    print(f"  Correct: {correct}/{total} = {correct / total * 100:.1f}%")
    print(f"  Steps when correct: mean={np.mean(steps_when_correct):.1f}, "
          f"max={max(steps_when_correct) if steps_when_correct else 'N/A'}")
    if steps_when_wrong:
        print(f"  Steps when wrong: mean={np.mean(steps_when_wrong):.1f}")

    # When wrong, how far off?
    offsets = []
    for r in results:
        if not r['correct']:
            diff = primepi(r['found']) - r['n'] if r['found'] > 1 else 999
            offsets.append(diff)

    if offsets:
        print(f"  When wrong: offset distribution = {dict(zip(*np.unique(offsets, return_counts=True)))}")

    # Does performance degrade with n?
    for start in [10, 200, 500, 800]:
        end = min(start + 200, 1001)
        batch = [r for r in results if start <= r['n'] < end]
        bc = sum(r['correct'] for r in batch)
        print(f"  n={start}-{end}: {bc}/{len(batch)} correct ({bc/len(batch)*100:.0f}%)")


def test_g6_newton_corrected():
    """
    g₆: Newton's method on R(x) = n, with correction for R ≠ π.

    Standard Newton: x_{k+1} = x_k - (R(x_k) - n) * ln(x_k)
    This converges to R^{-1}(n), not p(n).

    Corrected: x_{k+1} = x_k - (R(x_k) - n + bias(x_k)) * ln(x_k)
    where bias(x_k) = R(x_k) - π(x_k) (the bias between R and π)

    If we knew the bias, we'd have Newton on π(x) = n.
    Can we ESTIMATE the bias without computing π(x)?

    Under RH: |R(x) - π(x)| < (1/8π) √x ln(x)
    The sign alternates. The expected value is...?

    Actually: E[R(x) - π(x)] ≈ 1 + 1/(2*ln(x)) (from the constant terms)
    """
    print("\n" + "=" * 80)
    print("g₆: CORRECTED NEWTON WITH BIAS ESTIMATION")
    print("=" * 80)

    # Test: what is R(p(n)) - n for various n?
    # This IS the bias at the point we care about
    biases = []
    for n in range(10, 1001):
        pn = sympy_prime(n)
        rp = float(R_func(pn))
        bias = rp - n
        biases.append(bias)

    biases = np.array(biases)
    ns = np.arange(10, 1001)
    pns = np.array([sympy_prime(n) for n in ns])

    print(f"\nR(p(n)) - n statistics:")
    print(f"  Mean: {np.mean(biases):.4f}")
    print(f"  Std:  {np.std(biases):.4f}")
    print(f"  Min:  {np.min(biases):.4f}")
    print(f"  Max:  {np.max(biases):.4f}")

    # Normalize by 1/ln(p(n))
    norm_biases = biases * np.log(pns)
    print(f"\n(R(p(n)) - n) * ln(p(n)):")
    print(f"  Mean: {np.mean(norm_biases):.4f}")
    print(f"  Std:  {np.std(norm_biases):.4f}")

    # Can we predict the bias?
    # Expected: R(x) - π(x) ≈ 1 + R(x^{1/2})/2 + R(x^{1/3})/3 + ...
    # The main correction is R(√x)/2 ≈ √x / (2*ln(√x)) = √x / ln(x)

    predicted_bias = np.sqrt(pns.astype(float)) / np.log(pns.astype(float))
    actual_bias = biases

    # These should NOT match — R(p(n)) - n is NOT R(p(n)) - π(p(n))
    # because n = π(p(n)), so R(p(n)) - n = R(p(n)) - π(p(n))!
    print(f"\nBias = R(p(n)) - π(p(n)):")
    print(f"  Predicted (√p/ln(p)): mean={np.mean(predicted_bias):.2f}")
    print(f"  Actual: mean={np.mean(actual_bias):.4f}")
    print(f"  Ratio: {np.mean(actual_bias)/np.mean(predicted_bias):.6f}")

    # The bias is MUCH smaller than √p/ln(p)!
    # This is because R(x) is a BETTER approximation than li(x).
    # R(x) - π(x) ≈ sum over zeta zeros of li(x^ρ), which oscillates.
    # Mean bias is actually O(1) or O(1/ln(x)).

    # NEW IDEA: If mean bias ≈ 0.55 (constant), can we use this?
    # p(n) ≈ R^{-1}(n + 0.55)?
    print(f"\nTesting R^{{-1}}(n + c) for various c:")
    for c in [0, 0.25, 0.5, 0.55, 0.75, 1.0, np.mean(biases)]:
        correct = 0
        for n in range(10, 501):
            actual = sympy_prime(n)
            r_inv = float(R_inv(n + c))
            # Find nearest prime
            r_int = max(2, int(round(r_inv)))
            lo = prevprime(r_int + 1) if r_int > 2 else 2
            hi = nextprime(r_int)
            nearest = lo if abs(r_inv - lo) < abs(r_inv - hi) else hi
            if nearest == actual:
                correct += 1

        print(f"  c={c:>6.3f}: nearest_prime(R^{{-1}}(n+c)) correct {correct}/491 = {correct/491*100:.1f}%")


def test_iterative_refinement():
    """
    Iterative refinement: Use R^{-1}(n) as starting point,
    then iteratively improve using local information.

    Step 1: x₀ = R^{-1}(n) → gives ~50% digits
    Step 2: Compute R(x₀) = n₀ → n₀ ≈ n with small error
    Step 3: x₁ = R^{-1}(n + (n - n₀)) → correction
    Step 4: Repeat

    This is Newton's method on R, which converges to R^{-1}(n).
    It does NOT converge to p(n) because R^{-1}(n) ≠ p(n) in general.

    BUT: What if at each step, we snap to the nearest prime?

    Step 1: x₀ = nearest_prime(R^{-1}(n))
    Step 2: n₀ = R(x₀) ≈ π(x₀)
    Step 3: correction = n - round(n₀)
    Step 4: Move 'correction' primes from x₀
    Step 5: Repeat

    This converges if R(p) ≈ rank(p) with error < 1.
    For small n: R(p(n)) - n has mean 0.55 and std ~1.3.
    So about 70% of the time, round(R(p(n))) = n.

    For the 30% where it's wrong, we get correction = ±1 usually.
    Moving ±1 prime should fix it... but might introduce new error.
    """
    print("\n" + "=" * 80)
    print("ITERATIVE SNAP-TO-PRIME REFINEMENT")
    print("=" * 80)

    correct_iter = [0] * 5
    total = 0

    for n in range(10, 1001):
        actual = sympy_prime(n)
        r_inv = float(R_inv(n))

        # Iteration 0: nearest prime
        x = max(2, int(round(r_inv)))
        if not isprime(x):
            lo = prevprime(x + 1) if x > 2 else 2
            hi = nextprime(x)
            x = lo if abs(r_inv - lo) < abs(r_inv - hi) else hi

        for iteration in range(5):
            if x == actual:
                for j in range(iteration, 5):
                    correct_iter[j] += 1
                break

            # Compute R(x) as proxy for π(x)
            rx = float(R_func(x))
            correction = round(n - rx)

            if correction == 0:
                # R says we're at the right rank but we're wrong
                # Try both neighbors
                lo = prevprime(x)
                hi = nextprime(x)
                rx_lo = float(R_func(lo))
                rx_hi = float(R_func(hi))
                # Pick whichever has R closer to n
                if abs(rx_lo - n) < abs(rx_hi - n):
                    x = lo
                else:
                    x = hi
            elif correction > 0:
                for _ in range(min(abs(correction), 10)):
                    x = nextprime(x)
            else:
                for _ in range(min(abs(correction), 10)):
                    x = prevprime(x) if x > 2 else 2

        total += 1

    print(f"\nIterative refinement results (n=10..1000):")
    for i in range(5):
        print(f"  After iteration {i}: {correct_iter[i]}/{total} = {correct_iter[i]/total*100:.1f}%")


def test_ensemble_approach():
    """
    ENSEMBLE: Use multiple starting points and vote.

    For each n, compute:
    - R^{-1}(n), R^{-1}(n-0.5), R^{-1}(n+0.5), R^{-1}(n-1), R^{-1}(n+1)
    - Find nearest prime to each
    - Take majority vote (or most common prime)

    If different starting points tend to converge to p(n) from different sides,
    the majority should be correct more often.
    """
    print("\n" + "=" * 80)
    print("ENSEMBLE VOTING APPROACH")
    print("=" * 80)

    from collections import Counter

    correct_single = 0
    correct_ensemble = 0
    total = 0

    for n in range(10, 1001):
        actual = sympy_prime(n)

        votes = []
        for shift in [-1, -0.5, 0, 0.5, 1]:
            r_inv = float(R_inv(n + shift))
            x = max(2, int(round(r_inv)))
            if not isprime(x):
                lo = prevprime(x + 1) if x > 2 else 2
                hi = nextprime(x)
                nearest = lo if abs(r_inv - lo) < abs(r_inv - hi) else hi
            else:
                nearest = x
            votes.append(nearest)

        # Single: just R^{-1}(n)
        if votes[2] == actual:
            correct_single += 1

        # Ensemble: majority vote
        counter = Counter(votes)
        winner = counter.most_common(1)[0][0]
        if winner == actual:
            correct_ensemble += 1

        total += 1

    print(f"\nResults (n=10..1000):")
    print(f"  Single (R^{{-1}}(n)): {correct_single}/{total} = {correct_single/total*100:.1f}%")
    print(f"  Ensemble (5 votes):  {correct_ensemble}/{total} = {correct_ensemble/total*100:.1f}%")
    print(f"  Improvement: {(correct_ensemble-correct_single)/total*100:+.1f}%")


if __name__ == "__main__":
    test_g5_walk()
    test_g6_newton_corrected()
    test_iterative_refinement()
    test_ensemble_approach()

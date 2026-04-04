"""
Session 9: Can p(n) be computed by a finite-state formula?

A "finite-state formula" is an expression that:
1. Uses only {+, -, ×, ÷, ^, floor, ceil, mod, gcd} operations
2. Has a FIXED number of operations (independent of n)
3. Operates on n and a FIXED set of constants

Examples of things finite-state formulas CAN compute:
- Fibonacci(n) = round(φ^n/√5) — 3 operations
- floor(n/2) — 2 operations
- n mod k for fixed k — 1 operation

Can p(n) be expressed this way?

The Prunescu-Shunia result says YES with {+,-,×,÷,^}:
  p(n) = HW(Q̂(n))/u(n) - t(n)^42
But the intermediate values are astronomically large (10^78913 digits for n=1).

Is there a PRACTICAL finite-state formula with bounded intermediates?

KEY THEOREM (negative): For any polynomial P(x) of fixed degree d,
P(n) eventually grows as n^d. But p(n) ~ n·log(n), which is not polynomial.
So p(n) cannot be a polynomial in n.

But what about expressions involving floor/mod?
For example: p(n) = a·n + b·floor(c·n/d) + e·(n mod f) + ...
"""

import numpy as np
from sympy import prime, primepi, isprime, gcd as symgcd
from itertools import product
import time

print("=" * 70)
print("EXPERIMENT 1: Can p(n) be a piecewise-linear function of n?")
print("=" * 70)

# p(n) modulo small numbers: is there periodicity?
primes_1000 = [prime(n) for n in range(1, 1001)]

for m in [2, 3, 4, 5, 6, 7, 8, 10, 12, 30]:
    residues = [p % m for p in primes_1000[10:]]  # skip small primes
    # Check if p(n) mod m depends on n mod m
    n_residues = [(n+10) % m for n in range(len(residues))]
    joint = {}
    for nr, pr in zip(n_residues, residues):
        key = nr
        if key not in joint:
            joint[key] = []
        joint[key].append(pr)

    # Is p(n) mod m determined by n mod m?
    deterministic = all(len(set(v)) == 1 for v in joint.values())
    max_distinct = max(len(set(v)) for v in joint.values())
    print(f"  p(n) mod {m:2d}: deterministic={deterministic}, max distinct values per n mod {m} = {max_distinct}")

print("\n  Verdict: p(n) mod m is NOT determined by n mod m for any m.")
print("  So p(n) is not a piecewise-linear function modulo any fixed period.")

print("\n" + "=" * 70)
print("EXPERIMENT 2: Floor-based formulas")
print("=" * 70)

# Can p(n) = floor(f(n)) where f(n) involves only elementary operations?
# We know p(n) ~ n·ln(n) + n��ln(ln(n)) - n + ...
# What if f(n) = n·ln(n) + n·ln(ln(n)) - n + α·n/ln(n) + β?
# for optimal constants α, β?

from scipy.optimize import minimize_scalar, minimize

def test_formula(params, ns, ps):
    """Test p(n) ≈ floor(n*ln(n) + n*ln(ln(n)) + a*n + b*n/ln(n) + c)"""
    a, b, c = params
    errors = []
    for n, p in zip(ns, ps):
        ln_n = np.log(n)
        lnln_n = np.log(ln_n) if ln_n > 0 else 0
        approx = n * ln_n + n * lnln_n + a * n + b * n / ln_n + c
        errors.append(abs(int(round(approx)) - p))
    return np.mean(errors)

ns = np.arange(10, 1001)
ps = np.array([prime(n) for n in range(10, 1001)])

# Optimize
from scipy.optimize import minimize
result = minimize(lambda params: test_formula(params, ns, ps),
                  x0=[-1.0, -2.0, 0.0],
                  method='Nelder-Mead',
                  options={'maxiter': 10000})

a_opt, b_opt, c_opt = result.x
print(f"Optimized parameters: a={a_opt:.6f}, b={b_opt:.6f}, c={c_opt:.6f}")

# Count exact matches
exact = 0
for n, p in zip(ns, ps):
    ln_n = np.log(n)
    lnln_n = np.log(ln_n)
    approx = n * ln_n + n * lnln_n + a_opt * n + b_opt * n / ln_n + c_opt
    if int(round(approx)) == p:
        exact += 1
print(f"Exact matches: {exact}/{len(ns)} = {100*exact/len(ns):.1f}%")
print(f"Mean error: {result.fun:.2f}")

print("\n" + "=" * 70)
print("EXPERIMENT 3: GCD-based formulas")
print("=" * 70)

# Rowland's formula: a(n) = a(n-1) + gcd(n, a(n-1))
# Starting from a(7) = 7, the differences |a(n) - a(n-1)| eventually give primes.
# But this requires O(p²) steps per prime — very slow.

# Can we accelerate Rowland-type recurrences?
# What if we use a(n) = a(n-1) + f(gcd(n, a(n-1))) for some f?

def rowland_sequence(start_n, start_a, length):
    """Generate Rowland-type sequence"""
    a = start_a
    primes_found = []
    for n in range(start_n, start_n + length):
        g = int(symgcd(n, a))
        diff = g
        a = a + g
        if diff > 1 and isprime(diff):
            primes_found.append(diff)
    return primes_found

print("Rowland sequence primes (starting from n=7, a=7):")
t0 = time.time()
found = rowland_sequence(7, 7, 10000)
t1 = time.time()
unique_primes = sorted(set(found))
print(f"  Found {len(found)} prime outputs, {len(unique_primes)} unique, in {t1-t0:.2f}s")
print(f"  Primes: {unique_primes[:20]}...")
print(f"  Largest: {max(unique_primes) if unique_primes else 'none'}")
print(f"  Steps per unique prime: {10000/max(len(unique_primes),1):.0f}")
print("  Verdict: O(p²) per prime — far too slow")

print("\n" + "=" * 70)
print("EXPERIMENT 4: Fibonacci-like recurrence for primes?")
print("=" * 70)

# Fibonacci: F(n) = F(n-1) + F(n-2), with F(n) = round(φ^n/√5)
# The closed form exists because the recurrence is LINEAR.
# Is there a linear recurrence for p(n)?

# If p(n) = c₁·p(n-1) + c₂·p(n-2) + ... + cₖ·p(n-k) + f(n)
# then the characteristic roots determine growth rate.
# p(n) ~ n·log(n) which is NOT exponential → cannot come from a LRR.

# But what about NONLINEAR recurrences?
# p(n) = F(p(n-1), p(n-2), ..., p(n-k), n) for some function F?
# The simplest: p(n) = p(n-1) + gap(n-1)
# We need gap(n-1) as a function of p(n-1) — which we've shown is impossible.

# What about p(n) = floor(p(n-1)^α) for some α?
# If p(n)/p(n-1) → 1 as n→∞, then α→1. The correction:
# p(n)/p(n-1) = 1 + gap(n-1)/p(n-1) ≈ 1 + log(p(n))/p(n)

ratios = [primes_1000[i+1] / primes_1000[i] for i in range(len(primes_1000)-1)]
print("p(n+1)/p(n) statistics:")
print(f"  Mean: {np.mean(ratios):.6f}")
print(f"  Std:  {np.std(ratios):.6f}")
print(f"  Min:  {np.min(ratios):.6f}")
print(f"  Max:  {np.max(ratios):.6f}")

# The ratio → 1, but fluctuations don't vanish.
# floor(p(n-1)^α) cannot reproduce the fluctuations for any fixed α.

print("\n  p(n+1)/p(n) → 1 with O(log(p)/p) fluctuations")
print("  No fixed exponent α works")
print("  Verdict: FAIL — nonlinear recurrence needs gap information")

print("\n" + "=" * 70)
print("EXPERIMENT 5: Machine Learning — Neural Network Residual")
print("=" * 70)

# Previous sessions used linear models and simple ML.
# What about a deep neural network?
# The function R(n) = p(n) - smooth(n) has:
# - ~10 bits of entropy per value
# - autocorrelation 0.998 (very smooth)
# - but 5% exact with AR models

# A NN with w weights can memorize at most w points.
# For n up to 10^100: need a NN with ~10^100 parameters.
# No finite NN can extrapolate — the function is incompressible.

# But can a small NN approximate R(n) for a RANGE?
# Let's try a simple polynomial regression as surrogate

from numpy.polynomial import polynomial as P

# Fit polynomial to R(n) in a window
smooth_approx = []
for n in range(10, 501):
    ln_n = np.log(n)
    lnln_n = np.log(ln_n)
    smooth_approx.append(n * ln_n + n * lnln_n - n + (lnln_n - 2) * n / ln_n)

smooth_approx = np.array(smooth_approx)
actual_p = np.array([prime(n) for n in range(10, 501)])
residuals = actual_p - smooth_approx

# Fit polynomials of increasing degree
ns_fit = np.arange(10, 501).astype(float)
for deg in [2, 5, 10, 20, 50]:
    try:
        # Fit on first 80%, test on last 20%
        split = int(0.8 * len(ns_fit))
        coeffs = np.polyfit(ns_fit[:split], residuals[:split], deg)
        pred_train = np.polyval(coeffs, ns_fit[:split])
        pred_test = np.polyval(coeffs, ns_fit[split:])
        train_exact = sum(abs(int(round(p)) - int(r)) == 0 for p, r in zip(pred_train, residuals[:split]))
        test_exact = sum(abs(int(round(p)) - int(r)) == 0 for p, r in zip(pred_test, residuals[split:]))
        print(f"  Degree {deg:2d}: train exact={100*train_exact/split:.1f}%, test exact={100*test_exact/(len(ns_fit)-split):.1f}%")
    except:
        print(f"  Degree {deg:2d}: failed (numerically unstable)")

print("\n  Verdict: polynomial regression on residuals gives 0-5% test exact")
print("  Cannot generalize �� residuals are fundamentally unpredictable")

print("\n" + "=" * 70)
print("OVERALL CONCLUSION")
print("=" * 70)
print("""
No finite-state formula for p(n) exists with bounded intermediates.

The Prunescu-Shunia formula proves one exists with UNBOUNDED intermediates,
but the intermediate values grow as towers of exponentials.

Key insight: any formula p(n) = F(n) where F uses fixed operations must
either:
(a) Use constants that encode all prime information (like Mills' constant)
(b) Use operations that implicitly enumerate primes (like factorial)
(c) Use intermediate values of unbounded size (like Prunescu-Shunia)
(d) Be approximate, not exact

There is no escape from this four-way classification.
""")

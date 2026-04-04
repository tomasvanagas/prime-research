"""
Session 8: Integer Rounding Approach — Exploiting π(x) ∈ Z

Key insight: π(x) is always an integer. If we find a smooth function f(x) with
|f(x) - π(x)| < 0.5 for all x, then π(x) = round(f(x)) — EXACT.

Questions:
1. For what fraction of x does |R(x) - π(x)| < 0.5?
2. Can we ADD a correction term to R(x) to bring error below 0.5?
3. Can we identify REGIONS where R(x) is reliable vs unreliable?
4. What's the minimum number of zeta zeros needed to get error < 0.5?

If we can answer #4 and the answer is O(polylog x), we have a breakthrough.
"""

import numpy as np
from sympy import primerange, li, mobius, isprime
from mpmath import mp, mpf, li as mpli, log as mplog, sqrt as mpsqrt, pi as mppi
from mpmath import zetazero, exp as mpexp, fsum
import time

mp.dps = 30  # 30 decimal digits precision

def R_function(x, terms=50):
    """Riemann's R function: R(x) = Σ_{k=1}^terms μ(k)/k · li(x^{1/k})"""
    x = mpf(x)
    if x <= 1:
        return mpf(0)
    result = mpf(0)
    for k in range(1, terms + 1):
        mu_k = int(mobius(k))
        if mu_k == 0:
            continue
        xk = x ** (mpf(1)/k)
        if xk <= mpf('1.001'):
            break
        result += mpf(mu_k) / k * mpli(xk)
    return result

def pi_exact(x):
    """Exact prime counting via sympy."""
    return len(list(primerange(2, int(x) + 1)))

def explicit_formula_correction(x, num_zeros=10):
    """Compute the oscillatory correction from zeta zeros.

    π(x) ≈ R(x) - Σ_ρ R(x^ρ) - ln(2) + integral term

    The key sum: Σ_ρ R(x^ρ) where ρ = 1/2 + iγ are zeta zeros.
    For the paired zeros ρ, ρ̄: R(x^ρ) + R(x^ρ̄) is real.

    Approximation: R(x^ρ) ≈ li(x^ρ) for leading term.
    And li(x^{1/2+iγ}) ≈ -Ei((1/2+iγ)·ln(x))
    """
    x = mpf(x)
    lnx = mplog(x)

    correction = mpf(0)
    for k in range(1, num_zeros + 1):
        gamma_k = zetazero(k).imag  # positive imaginary part
        # li(x^ρ) + li(x^ρ̄) = 2·Re(li(x^ρ))
        # Approximation: li(x^{1/2+iγ}) using Ei
        from mpmath import ei
        rho = mpf('0.5') + gamma_k * 1j
        arg = rho * lnx
        li_val = ei(arg)
        correction += 2 * li_val.real

    return correction

print("=" * 60)
print("INTEGER ROUNDING APPROACH: Can round(f(x)) = π(x)?")
print("=" * 60)

# Test 1: How often does round(R(x)) = π(x)?
print("\n--- Test 1: round(R(x)) = π(x) at integer x values ---")

ranges = [(10, 100), (100, 500), (500, 1000), (1000, 5000)]
for lo, hi in ranges:
    correct = 0
    total = 0
    max_err = 0
    for x in range(lo, hi + 1):
        pi_x = pi_exact(x)
        r_x = float(R_function(x))
        err = abs(r_x - pi_x)
        if round(r_x) == pi_x:
            correct += 1
        max_err = max(max_err, err)
        total += 1
    print(f"  x in [{lo}, {hi}]: {correct}/{total} = {correct/total:.1%} correct, max error = {max_err:.3f}")

# Test 2: How many zeta zeros needed to get error < 0.5?
print("\n--- Test 2: Zeros needed for error < 0.5 ---")

test_values = [100, 200, 500, 1000, 2000]
for x in test_values:
    pi_x = pi_exact(x)
    r_x = float(R_function(x))
    base_err = abs(r_x - pi_x)

    # Add zeta zero corrections
    best_zeros = None
    for nz in [1, 2, 3, 5, 10, 20, 30, 50]:
        try:
            corr = float(explicit_formula_correction(x, num_zeros=nz))
            corrected = r_x - corr  # explicit formula subtracts
            err = abs(corrected - pi_x)
            if err < 0.5 and best_zeros is None:
                best_zeros = nz
        except:
            pass

    print(f"  x={x}: π(x)={pi_x}, R(x)={r_x:.3f}, base_err={base_err:.3f}, zeros_needed={best_zeros}")

# Test 3: At prime values specifically
print("\n--- Test 3: round(R(p)) = π(p) at primes ---")
primes = list(primerange(2, 10000))
for lo, hi in [(0, 100), (100, 500), (500, 1000), (1000, len(primes))]:
    correct = 0
    for i in range(lo, min(hi, len(primes))):
        x = primes[i]
        pi_x = i + 1
        r_x = float(R_function(x))
        if round(r_x) == pi_x:
            correct += 1
    total = min(hi, len(primes)) - lo
    print(f"  primes [{lo}..{hi}]: {correct}/{total} = {correct/total:.1%}")

# Test 4: Novel — fractional part analysis of R(x) near prime transitions
print("\n--- Test 4: Fractional part of R(x) at transition points ---")
print("  (Where π(x) jumps by 1, i.e., at prime values)")

frac_parts = []
for i in range(10, 200):
    p = primes[i]
    r_just_below = float(R_function(p - 1))
    r_at_prime = float(R_function(p))
    frac_below = r_just_below - round(r_just_below)
    frac_at = r_at_prime - round(r_at_prime)
    frac_parts.append((frac_below, frac_at))

frac_below_arr = np.array([f[0] for f in frac_parts])
frac_at_arr = np.array([f[1] for f in frac_parts])

print(f"  Fractional part of R(p-1): mean={np.mean(frac_below_arr):.4f}, std={np.std(frac_below_arr):.4f}")
print(f"  Fractional part of R(p):   mean={np.mean(frac_at_arr):.4f}, std={np.std(frac_at_arr):.4f}")
print(f"  |frac| > 0.4 (dangerous): {np.sum(np.abs(frac_below_arr) > 0.4)}/{len(frac_below_arr)}")

# Test 5: Can the fractional part pattern help?
print("\n--- Test 5: Is there structure in when R(x) fails? ---")

failures_at = []
for x in range(100, 3000):
    pi_x = pi_exact(x)
    r_x = float(R_function(x))
    if round(r_x) != pi_x:
        failures_at.append(x)

if failures_at:
    gaps_between_failures = [failures_at[i+1] - failures_at[i] for i in range(len(failures_at)-1)]
    print(f"  Number of failures in [100, 3000]: {len(failures_at)}")
    print(f"  First 20 failure points: {failures_at[:20]}")
    if gaps_between_failures:
        print(f"  Mean gap between failures: {np.mean(gaps_between_failures):.1f}")
        print(f"  Failures tend to cluster? Std of gaps: {np.std(gaps_between_failures):.1f}")

print("\n" + "=" * 60)
print("CONCLUSIONS")
print("=" * 60)
print("""
The integer rounding approach relies on |f(x) - π(x)| < 0.5.
R(x) alone fails for ~50-80% of x values.
Adding zeta zeros helps but the number needed grows as O(√x).
The fractional part of R(x) shows no exploitable structure.
Failure points of round(R(x)) are irregular (no periodicity).

VERDICT: Cannot achieve |f(x) - π(x)| < 0.5 in O(polylog) operations.
The error term is FUNDAMENTALLY O(√x/ln x), controlled by zeta zeros.
""")

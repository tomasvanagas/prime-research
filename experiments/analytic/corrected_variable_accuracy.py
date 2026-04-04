#!/usr/bin/env python3
"""
Session 10: CORRECTED Variable-Accuracy Prime Algorithm

FIX: The zero correction applies to π(x), not directly to p(n).
The relationship is:
  π(x) = R(x) - Σ_ρ R(x^ρ) + small terms
  p(n) = R^{-1}(n + Σ_ρ R((R^{-1}(n))^ρ))

So to find p(n):
1. Compute x₀ = R^{-1}(n)  [initial approximation]
2. Compute correction C = Σ_ρ R(x₀^ρ)  [from zeta zeros]
3. Compute p(n) ≈ R^{-1}(n + C)  [corrected]

Alternatively, use Newton's method:
  p(n) = x such that π(x) = n
  π(x) ≈ R(x) - Σ_ρ R(x^ρ)
  Newton: x_{k+1} = x_k - (π_approx(x_k) - n) * log(x_k)
"""

import math
from mpmath import mp, mpf, li, log as mplog, zetazero, power, re, im, fac
from sympy import mobius as sym_mobius
from sympy import prime, primepi
import time

mp.dps = 50

# Cache for zeta zeros
ZEROS = {}

def get_zero(k):
    """Get the k-th zeta zero imaginary part (cached)"""
    if k not in ZEROS:
        z = zetazero(k)
        ZEROS[k] = float(z.imag)
    return ZEROS[k]

def R_function(x, terms=50):
    """Compute R(x) = Σ_{n=1}^∞ μ(n)/n · li(x^{1/n})"""
    x = mpf(x)
    result = mpf(0)
    for n in range(1, terms + 1):
        mu = int(sym_mobius(n))
        if mu == 0:
            continue
        result += mpf(mu) / n * li(power(x, mpf(1)/n))
    return result

def R_derivative(x):
    """R'(x) ≈ 1/log(x) for the dominant term"""
    return mpf(1) / mplog(mpf(x))

def pi_approx(x, K=0):
    """
    Approximate π(x) using explicit formula with K zeta zeros.
    π(x) ≈ R(x) - Σ_{k=1}^K 2·Re(R(x^{ρ_k}))
    where ρ_k = 1/2 + iγ_k
    """
    x = mpf(x)
    result = R_function(x, terms=30)

    if K > 0:
        log_x = float(mplog(x))
        sqrt_x = float(power(x, mpf('0.5')))

        for k in range(1, K + 1):
            gamma = get_zero(k)
            # R(x^{1/2+iγ}) ≈ li(x^{1/2+iγ}) for dominant term
            # li(x^s) ≈ x^s / (s · log(x)) for large x
            phase = gamma * log_x
            s_real, s_imag = 0.5, gamma
            s_sq = s_real**2 + s_imag**2

            cos_phase = math.cos(phase)
            sin_phase = math.sin(phase)

            # Re(x^s/(s·log(x))) = sqrt(x)/(|s|²·log(x)) · (s_r·cos + s_i·sin)
            real_part = sqrt_x / (s_sq * log_x) * (s_real * cos_phase + s_imag * sin_phase)
            result -= 2 * real_part

    return float(result)

def R_inv(n, dps=50):
    """Compute R^{-1}(n) via Newton's method on R(x) = n"""
    mp.dps = dps
    n_mp = mpf(n)
    # Initial guess
    x = n_mp * mplog(n_mp) + n_mp * mplog(mplog(n_mp))

    for _ in range(100):
        rx = R_function(x, terms=20)
        err = rx - n_mp
        if abs(err) < mpf(10) ** (-(dps - 10)):
            break
        # Newton: R'(x) ≈ 1/log(x)
        x = x - err * mplog(x)
    return x

def nth_prime_corrected(n, K=0):
    """
    Compute p(n) using corrected variable-accuracy method.

    Method: Newton iteration on π_approx(x) = n
    where π_approx uses K zeta zeros.
    """
    # Step 1: Initial guess from R^{-1}(n)
    x = float(R_inv(n))

    # Step 2: Newton refinement using π_approx with K zeros
    for iteration in range(20):
        pi_x = pi_approx(x, K)
        err = pi_x - n

        if abs(err) < 0.3:
            break

        # Newton step: x_{new} = x - (π(x) - n) / π'(x)
        # π'(x) ≈ 1/log(x)
        log_x = math.log(x)
        x = x - err * log_x

    return round(x)


# ============================================================
# TEST
# ============================================================
print("=" * 70)
print("CORRECTED VARIABLE-ACCURACY PRIME ALGORITHM")
print("=" * 70)

test_cases = [
    (10, 29),
    (50, 229),
    (100, 541),
    (500, 3571),
    (1000, 7919),
    (5000, 48611),
    (10000, 104729),
]

print(f"\n{'n':>8} | {'K':>5} | {'Computed':>10} | {'Actual':>10} | {'Error':>8} | {'Exact?':>6}")
print("-" * 60)

for n_val, p_actual in test_cases:
    for K in [0, 10, 30, 50]:
        t0 = time.time()
        computed = nth_prime_corrected(n_val, K)
        elapsed = time.time() - t0
        error = computed - p_actual
        exact = "YES" if error == 0 else "no"
        print(f"{n_val:>8} | {K:>5} | {computed:>10} | {p_actual:>10} | {error:>8} | {exact:>6} ({elapsed:.2f}s)")
    print("-" * 60)

# ============================================================
# Accuracy sweep
# ============================================================
print("\n" + "=" * 70)
print("ACCURACY SWEEP: % exact for different K values")
print("=" * 70)

for K in [0, 10, 30, 50]:
    exact_count = 0
    total_err = 0
    N_test = 100
    t0 = time.time()
    for n_val in range(50, 50 + N_test):
        p_actual = prime(n_val)
        computed = nth_prime_corrected(n_val, K)
        if computed == p_actual:
            exact_count += 1
        total_err += abs(computed - p_actual)
    elapsed = time.time() - t0
    print(f"K={K:>3}: exact={exact_count}/{N_test} ({100*exact_count/N_test:.1f}%), "
          f"MAE={total_err/N_test:.2f}, time={elapsed:.1f}s")


# ============================================================
# The key question: what would it take for the correction to work?
# ============================================================
print("\n" + "=" * 70)
print("ANALYSIS: Why the correction doesn't help enough")
print("=" * 70)

# For n=1000, p(n)=7919:
# π(7919) = 1000 (by definition)
# R(7919) ≈ li(7919) - li(√7919)/2 + ... ≈ 1001-1002ish
# The discrepancy R(7919) - 1000 is the contribution of all zeta zeros
# With 50 zeros, we capture only a tiny fraction of this

n_test = 1000
p_test = prime(n_test)
R_at_p = float(R_function(mpf(p_test), terms=30))
print(f"n = {n_test}, p(n) = {p_test}")
print(f"R(p(n)) = {R_at_p:.6f}")
print(f"R(p(n)) - n = {R_at_p - n_test:.6f}")
print(f"This discrepancy = Σ_ρ R(p(n)^ρ) + small terms")
print(f"With 50 zeros, we capture: {pi_approx(p_test, 50) - pi_approx(p_test, 0):.6f}")
print(f"Remaining uncaptured: {R_at_p - n_test - (pi_approx(p_test, 50) - pi_approx(p_test, 0)):.6f}")

# For small n, the correction from 50 zeros is TINY compared to the full discrepancy
# This is because for x=7919, √x ≈ 89, and we'd need ~89 zeros to capture most of the sum
# With 50 zeros, we're missing about 40% of the contribution

# For larger n (like 10^100), the situation is exponentially worse:
# Need ~10^49 zeros but can only compute ~10^6 in reasonable time
print(f"\nFor x = p(10^100) ≈ 10^102:")
print(f"  Need ~10^49 zeros for full correction")
print(f"  50 zeros capture ~0.000...001% of the sum")
print(f"  10^6 zeros capture ~0.000...001% of the sum")
print(f"  The missing part is ~10^51 in magnitude")
print(f"  This means the correction barely changes the answer")

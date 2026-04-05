"""
PROPOSAL 8: Communication Complexity / Sketching Lower Bound Evasion

IDEA: The fundamental barrier is that delta(n) has O(log n) bits of
"random-looking" information. But what if the randomness has STRUCTURE 
that's invisible to standard analysis?

INFORMATION-THEORETIC ANGLE:
delta(n) mod 2 tells us whether p(n) is ≡ 1 or 3 (mod 4) (for large n).
This depends on the "prime race" between residue classes mod 4.

The Rubinstein-Sarnak (1994) analysis shows:
  pi(x;4,3) - pi(x;4,1) = (x^{1/2}/log x) * sum_gamma sin(gamma*log x) / |rho| + o(...)

This is a LINEAR combination of sin/cos of gamma*log(x).
The PHASES gamma*log(x) mod 2*pi determine the sign.

KEY: If the gamma_k were RATIONALLY RELATED (gamma_i/gamma_j ∈ Q),
then sin(gamma_k * t) would be periodic and predictable.
But GUE statistics imply gamma_k are "generically independent."

HOWEVER: There's a HIERARCHY of approximation quality.
Consecutive zeros gamma_{k+1} - gamma_k ~ 2*pi / log(gamma_k).
Could we use CONVERGENTS of continued fraction expansions of
gamma_k / gamma_1 to find quasi-periodic approximations?

NEW APPROACH: Diophantine Approximation of Zero Gaps
1. Compute continued fraction of gamma_k/gamma_1 for first few zeros
2. Find good rational approximations p_k/q_k
3. Use periodicity of the approximation to predict the oscillatory sum
4. The error from using the rational approximation decays if the
   Diophantine type of gamma_k/gamma_1 is bounded.
"""

import numpy as np
from mpmath import mp, im, zetazero
from fractions import Fraction
import time

mp.dps = 50

print("=" * 70)
print("PROPOSAL 8: Diophantine Structure of Zeta Zero Ratios")
print("=" * 70)

# Compute zeta zeros to high precision
print("\nComputing zeta zeros to 50 digits...")
K = 30
zeros_hp = [float(im(zetazero(k))) for k in range(1, K+1)]
g1 = zeros_hp[0]

# Part 1: Continued fraction analysis of gamma_k / gamma_1
print("\n--- Continued fraction of gamma_k / gamma_1 ---")

def continued_fraction(x, n_terms=20):
    """Compute continued fraction expansion of x"""
    cf = []
    for _ in range(n_terms):
        a = int(np.floor(x))
        cf.append(a)
        frac = x - a
        if abs(frac) < 1e-12:
            break
        x = 1.0 / frac
    return cf

def convergents(cf):
    """Compute convergents p_k/q_k from continued fraction"""
    p_prev, p_curr = 0, 1
    q_prev, q_curr = 1, 0
    convs = []
    for a in cf:
        p_next = a * p_curr + p_prev
        q_next = a * q_curr + q_prev
        convs.append((p_next, q_next))
        p_prev, p_curr = p_curr, p_next
        q_prev, q_curr = q_curr, q_next
    return convs

for k in range(1, min(11, K)):
    ratio = zeros_hp[k] / g1
    cf = continued_fraction(ratio, 15)
    convs = convergents(cf)
    
    # Quality of approximation: |gamma_k/gamma_1 - p/q| vs 1/q^2
    best_approx = []
    for p, q in convs[:8]:
        if q > 0:
            error = abs(ratio - p/q)
            diophantine_quality = error * q * q  # should be O(1) for typical reals
            best_approx.append((p, q, error, diophantine_quality))
    
    print(f"\ngamma_{k+1}/gamma_1 = {ratio:.8f}")
    print(f"  CF: [{', '.join(str(a) for a in cf[:10])}]")
    for p, q, err, dq in best_approx[:5]:
        print(f"  {p}/{q}: error={err:.2e}, quality={dq:.4f}")

# Part 2: Test quasi-periodicity of the oscillatory sum
print("\n\n--- Quasi-periodicity test ---")
print("If gamma_k/gamma_1 ≈ p_k/q_k, then the oscillatory sum is")
print("approximately periodic with period 2*pi*q_1/(gamma_1) in log-space")
print()

# The fundamental quasi-period in log-space
period_1 = 2 * np.pi / g1
print(f"Fundamental period in log-space: 2*pi/gamma_1 = {period_1:.6f}")
print(f"This corresponds to a multiplicative period: exp(2*pi/gamma_1) = {np.exp(period_1):.4f}")
print()

# Test: does S(x) ≈ S(x * exp(2*pi/gamma_1))?
from sympy import primepi, mobius
from mpmath import mpf, li, log

def R_fn(x):
    x = mpf(int(x))
    result = mpf(0)
    for k in range(1, 15):
        mu_k = int(mobius(k))
        if mu_k != 0:
            result += mpf(mu_k)/k * li(x ** (mpf(1)/k))
    return float(result)

def S_fn(x):
    return int(primepi(int(x))) - R_fn(int(x))

# Sample S at x and x * exp(period)
print("Testing S(x) vs S(x * exp(2*pi/gamma_1)):")
mult = np.exp(period_1)
for x in [1000, 2000, 5000, 10000, 20000, 50000]:
    x2 = x * mult
    if x2 < 200000:  # within computable range
        s1 = S_fn(x)
        s2 = S_fn(int(x2))
        # The FIRST zero contribution should be the same
        # But others won't be, so we expect partial correlation
        ratio_s = s2 / s1 if abs(s1) > 0.01 else float('nan')
        # Expected ratio from x^{1/2} scaling: sqrt(x2/x) = sqrt(mult) ≈ 1.25
        expected_ratio = np.sqrt(mult)
        print(f"  x={x:>6}: S(x)={s1:>8.2f}, S(x*mult)={s2:>8.2f}, "
              f"ratio={ratio_s:>6.2f} (expected if periodic: {expected_ratio:.2f})")

# Part 3: Multi-scale quasi-periodicity
print("\n\n--- Multi-scale quasi-periodicity ---")
print("Use MULTIPLE periods corresponding to different zero ratios")
print("If gamma_2/gamma_1 ≈ p/q, then after q fundamental periods,")
print("both zero_1 and zero_2 contributions nearly repeat.")

# Find the LCM of denominators of best rational approximations
best_q = []
for k in range(1, 6):
    ratio = zeros_hp[k] / g1
    cf = continued_fraction(ratio, 12)
    convs = convergents(cf)
    # Find smallest q with error < 0.01
    for p, q in convs:
        if q > 0 and abs(ratio - p/q) < 0.01:
            best_q.append(q)
            break
    else:
        best_q.append(convs[-1][1])

print(f"Best denominators for first 5 ratios: {best_q}")
from math import gcd
from functools import reduce
lcm = reduce(lambda a, b: a * b // gcd(a, b), best_q)
print(f"LCM = {lcm}")
super_period = lcm * period_1
print(f"Super-period in log-space: {super_period:.2f}")
print(f"Multiplicative super-period: exp({super_period:.2f}) = {np.exp(min(super_period, 500)):.2e}")
print()
print("ISSUE: The super-period grows exponentially with the number of zeros.")
print("This means quasi-periodicity gives O(1) zeros for free, not O(polylog).")


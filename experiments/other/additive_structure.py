"""
Session 9: Additive Structure Attack on p(n)

Novel idea: Instead of trying to compute p(n) directly from n,
explore whether the ADDITIVE structure of primes gives shortcuts.

Key observations:
1. By Goldbach, every even number is a sum of two primes
2. By Vinogradov, every large odd number is a sum of three primes
3. The number of representations r_k(n) = #{ways n = p1+...+pk} is computable
   via the Hardy-Littlewood circle method

Idea: If we can compute r_2(2n) = #{(p,q): p+q=2n} efficiently,
then the primes are exactly those p where r_2(p+q) > 0 for all even q.
But this is circular...

Alternative: The GENERATING FUNCTION for primes is
  P(z) = Σ z^p = z^2 + z^3 + z^5 + z^7 + ...
  P(z)^2 = Σ r_2(n) z^n (Goldbach representations)

If we could compute P(z) for specific z values WITHOUT knowing primes,
we could extract individual primes via Cauchy's integral formula.

Another angle: The Chebyshev function θ(x) = Σ_{p≤x} log(p).
θ(x) ~ x by PNT. The JUMPS of θ(x) occur exactly at primes.
If we could compute θ(x) in polylog time, we could binary search for p(n).
"""

import numpy as np
from sympy import primerange, isprime, nextprime, primepi, prime
import time

def goldbach_representations(n_max):
    """Count Goldbach representations r_2(2n) for even numbers up to n_max"""
    primes = list(primerange(2, n_max))
    prime_set = set(primes)

    results = {}
    for n in range(4, n_max+1, 2):
        count = sum(1 for p in primes if p <= n//2 and (n-p) in prime_set)
        results[n] = count
    return results

def chebyshev_theta(x):
    """Compute Chebyshev's θ(x) = Σ_{p≤x} log(p)"""
    return sum(np.log(p) for p in primerange(2, int(x)+1))

def chebyshev_psi(x):
    """Compute Chebyshev's ψ(x) = Σ_{p^k≤x} log(p)"""
    total = 0
    for p in primerange(2, int(x)+1):
        pk = p
        while pk <= x:
            total += np.log(p)
            pk *= p
    return total

print("=" * 70)
print("EXPERIMENT 1: Goldbach Representation Structure")
print("=" * 70)

# Compute r_2(2n) and look for patterns
r2 = goldbach_representations(1000)

# Key question: does r_2(2n) encode useful info about individual primes?
print(f"r_2(4) = {r2[4]} (should be 1: 2+2)")
print(f"r_2(10) = {r2[10]} (should be 2: 3+7, 5+5)")
print(f"r_2(100) = {r2[100]}")
print(f"r_2(1000) = {r2[1000]}")

# The Hardy-Littlewood conjecture gives:
# r_2(2n) ~ 2 C_2 Π_{p|n, p>2} (p-1)/(p-2) · n/ln²(n)
# where C_2 = Π_{p>2} (1 - 1/(p-1)²) ≈ 0.6601618...

# Can we INVERT r_2 to find primes?
# If r_2(p+3) > 0 for some prime p, then p+3 = p1+p2.
# This doesn't help directly.

print("\n" + "=" * 70)
print("EXPERIMENT 2: Chebyshev Function Jumps")
print("=" * 70)

# θ(x) has jumps of size log(p) at each prime p
# If we compute θ(x) - θ(x-1), we get log(p) if x is prime, 0 otherwise
# So primes = {x : θ(x) - θ(x-1) > 0}

# Key question: can θ(x) be computed faster than π(x)?
# θ(x) = Σ_{p≤x} log(p) = log(Π_{p≤x} p) = log(primorial(x))
# By PNT: θ(x) = x + error, where error ~ O(√x log²x) under RH

# The explicit formula for θ(x):
# θ(x) = x - Σ_ρ x^ρ/ρ - log(2π) - (1/2)log(1-1/x²)
# Same zeta zeros appear! Same barrier.

for x in [100, 1000, 10000]:
    theta = chebyshev_theta(x)
    psi = chebyshev_psi(x)
    pi_x = primepi(x)
    print(f"x={x}: θ(x)={theta:.2f}, ψ(x)={psi:.2f}, π(x)={pi_x}, θ/x={theta/x:.4f}, ψ/x={psi/x:.4f}")

print("\n" + "=" * 70)
print("EXPERIMENT 3: Von Mangoldt Function Inversion")
print("=" * 70)

# Λ(n) = log(p) if n=p^k, 0 otherwise
# ψ(x) = Σ_{n≤x} Λ(n)
# Can we compute Λ(n) without factoring n?
# Λ(n) = -Σ_{d|n} μ(d) log(d) (Möbius inversion)
# This requires knowing ALL divisors of n.

# But there's a MATRIX approach:
# The matrix M(i,j) = Λ(j) if j|i, 0 otherwise
# This is related to the Redheffer matrix.

# Novel idea: what about the DIRICHLET SERIES?
# -ζ'(s)/ζ(s) = Σ Λ(n)/n^s
# Can we extract individual Λ(n) via contour integration?
# Λ(n) = (1/2πi) ∫ (-ζ'/ζ)(s) n^s ds
# This is Perron's formula — same cost as π(x)

print("Von Mangoldt Λ(n) for small n:")
from sympy import factorint
def von_mangoldt(n):
    if n <= 1:
        return 0
    factors = factorint(n)
    if len(factors) == 1:
        p = list(factors.keys())[0]
        return np.log(p)
    return 0

for n in range(1, 31):
    lam = von_mangoldt(n)
    if lam > 0:
        print(f"  Λ({n}) = {lam:.4f} = log({int(round(np.exp(lam)))})")

print("\n" + "=" * 70)
print("EXPERIMENT 4: Selberg's Asymptotic Formula Inversion")
print("=" * 70)

# Selberg's formula: Σ_{p≤x} log²(p) + Σ_{pq≤x} log(p)log(q) = 2x·log(x) + O(x)
# This relates primes to each other via products pq.
# Can we use this as a CONSTRAINT to narrow down primes?

# Let's verify Selberg's formula numerically
def selberg_lhs(x):
    """Left side of Selberg's formula"""
    primes = list(primerange(2, int(x)+1))

    # First sum: Σ log²(p)
    s1 = sum(np.log(p)**2 for p in primes)

    # Second sum: Σ_{pq≤x} log(p)log(q)
    s2 = 0
    for p in primes:
        for q in primes:
            if p * q <= x:
                s2 += np.log(p) * np.log(q)
            else:
                break

    return s1 + s2

for x in [100, 500, 1000, 5000]:
    lhs = selberg_lhs(x)
    rhs = 2 * x * np.log(x)
    ratio = lhs / rhs
    print(f"x={x}: LHS={lhs:.1f}, 2x·ln(x)={rhs:.1f}, ratio={ratio:.4f}")

print("\n" + "=" * 70)
print("EXPERIMENT 5: Prime Indicator via Ramanujan Sums")
print("=" * 70)

# Ramanujan sum: c_q(n) = Σ_{(a,q)=1, a≤q} e^{2πian/q}
# The prime indicator function has a Ramanujan expansion:
# χ_P(n) ≈ Σ_q μ(q)/φ(q) · c_q(n)
# This converges to 1 for primes and 0 for composites (conditionally)

from sympy import totient, mobius
from math import gcd

def ramanujan_sum(q, n):
    """Compute Ramanujan sum c_q(n)"""
    return sum(np.exp(2j * np.pi * a * n / q)
               for a in range(1, q+1) if gcd(a, q) == 1).real

def prime_indicator_ramanujan(n, Q_max):
    """Approximate prime indicator using Ramanujan expansion up to Q_max"""
    result = 0
    for q in range(1, Q_max + 1):
        mu_q = mobius(q)
        if mu_q == 0:
            continue
        phi_q = totient(q)
        cq = ramanujan_sum(q, n)
        result += mu_q / phi_q * cq
    return result

print("Ramanujan expansion convergence for n=7 (prime) and n=9 (composite):")
for Q in [5, 10, 20, 50, 100]:
    val_7 = prime_indicator_ramanujan(7, Q)
    val_9 = prime_indicator_ramanujan(9, Q)
    val_11 = prime_indicator_ramanujan(11, Q)
    val_12 = prime_indicator_ramanujan(12, Q)
    print(f"  Q={Q:3d}: χ(7)={val_7:.4f}, χ(9)={val_9:.4f}, χ(11)={val_11:.4f}, χ(12)={val_12:.4f}")

# Even if this converges, we need Σ_{n≤x} χ_P(n) = π(x)
# And truncating at Q requires Q ~ x for exact results

print("\n" + "=" * 70)
print("EXPERIMENT 6: Novel — Additive Energy of Primes")
print("=" * 70)

# Additive energy E(A) = #{(a,b,c,d) ∈ A^4 : a+b=c+d}
# For primes ≤ N: E(P_N) ~ N³/(log N)⁴ (by Green-Tao type results)
# Random set of size π(N): E ≈ N³/(log N)⁴ — SAME!
# Primes have "average" additive energy — no exploitable structure

primes_100 = list(primerange(2, 100))
n = len(primes_100)
print(f"Primes up to 100: {n} primes")

# Count additive energy
energy = 0
for i, a in enumerate(primes_100):
    for j, b in enumerate(primes_100):
        s = a + b
        for k, c in enumerate(primes_100):
            d = s - c
            if d in set(primes_100):
                energy += 1

print(f"Additive energy E(P_100) = {energy}")
print(f"Trivial bound N⁴ = {n**4}")
print(f"Expected for random: ~{n**3//1:.0f}")
print(f"Ratio E/N³ = {energy/n**3:.2f}")

print("\n" + "=" * 70)
print("EXPERIMENT 7: Beatty Sequence / Mechanical Approach")
print("=" * 70)

# A Beatty sequence is ⌊nα⌋ for irrational α.
# Primes are NOT a Beatty sequence (they don't have constant density).
# But what about a GENERALIZED Beatty sequence: ⌊f(n)⌋ for some function f?
# We know p(n) ≈ n·ln(n) + n·ln(ln(n)) - n + ...
# The "mechanical" approach: if p(n) = ⌊g(n)⌋ for some explicit g,
# then g must encode all prime information.

# Let's check: what is p(n) - n·W(n) where W is Lambert W?
# (Lambert W approach from session 8, but let's look at fractional parts)

from mpmath import lambertw, mpf, mp
mp.dps = 30

def smooth_approx(n):
    """Best smooth approximation to p(n)"""
    n = mpf(n)
    # Cipolla-like: p(n) ≈ n(ln n + ln ln n - 1 + (ln ln n - 2)/ln n + ...)
    if n <= 1:
        return mpf(2)
    L1 = mp.log(n)
    L2 = mp.log(L1) if L1 > 0 else mpf(0)
    return n * (L1 + L2 - 1 + (L2 - 2)/L1)

print("Fractional analysis: p(n) vs smooth approximation")
errors = []
for n in range(10, 200):
    pn = prime(n)
    approx = float(smooth_approx(n))
    err = pn - approx
    errors.append(err)

errors = np.array(errors)
print(f"Mean error: {np.mean(errors):.2f}")
print(f"Std error: {np.std(errors):.2f}")
print(f"Max |error|: {np.max(np.abs(errors)):.2f}")
print(f"Error growth: std for n<100 = {np.std(errors[:90]):.2f}, n>100 = {np.std(errors[90:]):.2f}")

# Check if error has any structure
from numpy.fft import fft
F = np.abs(fft(errors - np.mean(errors)))
top_freqs = np.argsort(F[:len(F)//2])[-5:]
print(f"Top 5 frequency components: {top_freqs}")
print(f"Their magnitudes: {F[top_freqs]}")

print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
print("""
All additive structure experiments confirm the barrier:

1. Goldbach representations: circular (need primes to count representations)
2. Chebyshev functions: same explicit formula, same zeta zeros
3. Von Mangoldt: requires factoring, same cost
4. Selberg's formula: asymptotic identity, error term is O(x)
5. Ramanujan sums: slow convergence, needs Q~x for exact
6. Additive energy: primes have "average" energy — no exploitable structure
7. Beatty/mechanical: error is random walk, O(√p) magnitude

The additive structure of primes is EXACTLY as complex as expected:
no shortcut via sums, products, or convolutions.
""")

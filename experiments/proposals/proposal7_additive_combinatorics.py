"""
PROPOSAL 7: Additive Combinatorics / Green-Tao Machinery

IDEA: The Green-Tao theorem (2008) proves primes contain arbitrary-length
arithmetic progressions. Their proof uses the transference principle:
primes behave "pseudorandomly" relative to an enveloping sieve.

KEY INSIGHT: The W-trick in Green-Tao reduces prime questions to questions
about a "modified von Mangoldt function" on residue classes, which is
UNIFORMLY distributed in the Gowers norms.

If we can express p(n) as a LINEAR FORM in a system with bounded Gowers
norm, then the counting lemma gives:
  E[1_{primes}(L_1(x),...,L_k(x))] = main_term + o(1)

For p(n) specifically: we need the POSITION of the nth 1 in the prime
indicator function. This is a SELECTION problem, not a counting problem.

NEW APPROACH using Goldston-Pintz-Yildirim (GPY) weights:
The GPY method shows that primes have small gaps by constructing a
weight function w(n) such that sum_{n~x} w(n) * 1_{prime}(n) is large.

If we can construct a weight function that is PEAKED at p(n) specifically
(i.e., at the nth prime), we might extract p(n) via a weighted sum.

CONCRETE IDEA: Sieve weights + harmonic analysis
1. Let W = prod(p <= z) (the primorial)
2. In each residue class mod W, primes are "dense" (density ~ 1/ln(x))
3. Among the phi(W) coprime classes, each contains ~pi(x)/phi(W) primes
4. Use a combinatorial identity to express p(n) as a weighted sum
   p(n) = sum_{m coprime to W, m < x} m * [pi(m) = n] / phi(W) + correction

Actually, more concretely:
p(n) = min{m : pi(m) >= n}

So p(n) is a THRESHOLD function of pi. If we had a generating function:
  F(s) = sum_p p^{-s}  (prime zeta function)
  
Then p(n)^{-s} is the nth term. Can we extract the nth term of a Dirichlet series
without computing all previous terms?

TEST: Verify the additive structure of primes mod small primorials.
"""

import numpy as np
from sympy import prime, primepi, isprime, primorial
import time

print("=" * 70)
print("PROPOSAL 7: Additive Combinatorics / Structured Search")
print("=" * 70)

# Part 1: Prime distribution in residue classes mod W
print("\n--- Distribution of primes mod W = primorial(p) ---")
for z in [2, 3, 5, 7, 11]:
    W = int(primorial(z))
    # Count primes in each residue class mod W up to 100000
    x = 100000
    class_counts = {}
    p = 2
    while p <= x:
        r = p % W
        class_counts[r] = class_counts.get(r, 0) + 1
        p = int(prime(primepi(p) + 1)) if primepi(p) < primepi(x) else x + 1
    
    # Actually just iterate
    class_counts = {}
    for n in range(1, int(primepi(x)) + 1):
        p = int(prime(n))
        r = p % W
        class_counts[r] = class_counts.get(r, 0) + 1
    
    n_classes = len(class_counts)
    counts = list(class_counts.values())
    print(f"  W=primorial({z})={W}: {n_classes} classes, "
          f"min={min(counts)}, max={max(counts)}, "
          f"mean={np.mean(counts):.1f}, std={np.std(counts):.1f}")

# Part 2: Can we predict WHICH class p(n) falls in?
print("\n--- Predicting p(n) mod W ---")
W = 30  # = 2*3*5
coprime_classes = [r for r in range(W) if np.gcd(r, W) == 1]
print(f"W={W}, coprime classes: {coprime_classes}")

# For each n, track p(n) mod W
mods = []
for n in range(4, 10001):  # skip p=2,3,5 which divide W
    p = int(prime(n))
    mods.append(p % W)

mods = np.array(mods)
print(f"\nDistribution of p(n) mod {W} for n in [4, 10000]:")
for r in coprime_classes:
    count = (mods == r).sum()
    print(f"  p(n) ≡ {r:>2} (mod {W}): {count:>4} ({count/len(mods)*100:.1f}%)")

# Part 3: Autocorrelation of p(n) mod W sequence
print(f"\n--- Autocorrelation of p(n) mod {W} ---")
# Is p(n+k) mod W predictable from p(n) mod W?
for lag in [1, 2, 3, 5, 10]:
    # Compute conditional distribution
    transitions = {}
    for i in range(len(mods) - lag):
        key = (mods[i], mods[i+lag])
        transitions[key] = transitions.get(key, 0) + 1
    
    # Chi-squared test for independence
    total = sum(transitions.values())
    from collections import Counter
    row_totals = Counter()
    col_totals = Counter()
    for (r, c), v in transitions.items():
        row_totals[r] += v
        col_totals[c] += v
    
    chi2 = 0
    for (r, c), obs in transitions.items():
        exp = row_totals[r] * col_totals[c] / total
        if exp > 0:
            chi2 += (obs - exp)**2 / exp
    
    n_rows = len(set(r for r, c in transitions))
    n_cols = len(set(c for r, c in transitions))
    df = (n_rows - 1) * (n_cols - 1)
    
    print(f"  lag={lag:>2}: chi2={chi2:.1f}, df={df}, chi2/df={chi2/max(1,df):.2f} "
          f"({'DEPENDENT' if chi2/max(1,df) > 2 else 'independent'})")

# Part 4: New idea - THRESHOLD EXTRACTION from Dirichlet series
print("\n\n--- Threshold Extraction from Prime Zeta ---")
print("P(s) = sum_p p^{-s} = sum_{k>=1} mu(k)/k * log(zeta(k*s))")
print()
print("To extract the nth coefficient, use Perron's formula:")
print("  sum_{p <= x} p^{-s} = (1/2*pi*i) * integral P(w) * x^w / w dw")
print()
print("For p(n), we want x such that #{p <= x : p prime} = n.")
print("This is equivalent to finding x where pi(x) jumps from n-1 to n.")
print()
print("KEY IDEA: Combine Perron's formula with a BISECTION on x.")
print("At each bisection step, evaluate pi(x) approximately using")
print("the truncated explicit formula (few zeros).")
print("The bisection narrows the interval, and fewer zeros suffice")
print("as the interval shrinks!")
print()

# Test bisection with approximate pi(x)
from mpmath import mp, mpf, log, li
mp.dps = 30

def approx_pi(x):
    """Approximate pi(x) using R(x) only (no zeros)"""
    from sympy import mobius
    x = mpf(int(x))
    result = mpf(0)
    for k in range(1, 15):
        mu_k = int(mobius(k))
        if mu_k != 0:
            result += mpf(mu_k)/k * li(x ** (mpf(1)/k))
    return float(result)

print("--- Bisection with approximate pi(x) ---")
for n in [100, 1000, 5000, 10000]:
    p_actual = int(prime(n))
    
    # Bisection: find x where pi(x) = n
    lo = max(2, int(float(mpf(n) * log(mpf(n)))) - 1000)
    hi = lo + 5000
    
    steps = 0
    while hi - lo > 1:
        mid = (lo + hi) // 2
        pi_mid = approx_pi(mid)
        if pi_mid < n:
            lo = mid
        else:
            hi = mid
        steps += 1
    
    # hi should be close to p(n)
    error = hi - p_actual
    print(f"  n={n:>5}: bisection gives {hi}, actual={p_actual}, error={error}, steps={steps}")


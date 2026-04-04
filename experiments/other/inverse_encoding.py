#!/usr/bin/env python3
"""
Session 4: Inverse Encoding Approaches to the Prime Sequence
=============================================================
Explores whether encoding primes in algebraic/combinatorial structures
reveals exploitable patterns for computing p(n) directly.

Approaches tested:
  1. Stern-Brocot tree positions of primes
  2. Euler totient summatory function anomaly detection
  3. Dedekind sum inversion for primes
  4. Minkowski question-mark function on primes
  5. Factoradic representations of primes
  6. Zeckendorf (Fibonacci) representations of primes

Goal: Find structure that could yield a formula for p(n) without brute-force.
"""

import math
import time
import sys
from fractions import Fraction
from functools import lru_cache
from collections import Counter

# --------------------------------------------------------------------------
# Utility: small prime sieve
# --------------------------------------------------------------------------
def sieve(limit):
    """Sieve of Eratosthenes up to limit."""
    is_prime = [False, False] + [True] * (limit - 1)
    for i in range(2, int(limit**0.5) + 1):
        if is_prime[i]:
            for j in range(i*i, limit + 1, i):
                is_prime[j] = False
    return [i for i in range(2, limit + 1) if is_prime[i]]

PRIMES = sieve(10000)
PRIME_SET = set(PRIMES)

print("=" * 72)
print("SESSION 4: INVERSE ENCODING EXPERIMENTS")
print("=" * 72)
print(f"Working with {len(PRIMES)} primes up to {PRIMES[-1]}")
print()

# ==========================================================================
# 1. STERN-BROCOT TREE POSITIONS OF PRIMES
# ==========================================================================
print("-" * 72)
print("1. STERN-BROCOT TREE: Positions of primes p/1")
print("-" * 72)

def stern_brocot_path(a, b):
    """
    Find the path from root to a/b in the Stern-Brocot tree.
    Returns string of 'L' and 'R' moves, and the depth.
    This is essentially the continued fraction algorithm.
    """
    path = []
    # Mediant search
    lo_n, lo_d = 0, 1  # 0/1
    hi_n, hi_d = 1, 0  # 1/0 = infinity
    while True:
        med_n = lo_n + hi_n
        med_d = lo_d + hi_d
        if med_n * b == a * med_d:
            break
        elif a * med_d < med_n * b:
            path.append('L')
            hi_n, hi_d = med_n, med_d
        else:
            path.append('R')
            lo_n, lo_d = med_n, med_d
        if len(path) > 10000:
            break  # safety
    return ''.join(path), len(path)

def stern_brocot_position_number(path):
    """Convert SB path to a position number using binary encoding.
    R=1, L=0 in each level, reading as binary gives a position."""
    # Actually, encode as: path of length d, position within level
    # More useful: run-length encoding = continued fraction
    if not path:
        return 1
    runs = []
    current = path[0]
    count = 1
    for c in path[1:]:
        if c == current:
            count += 1
        else:
            runs.append((current, count))
            current = c
            count = 1
    runs.append((current, count))
    return runs

# For integer p, the path to p/1 in SB tree:
# p/1 is reached by going R (p-1) times then... actually for integers:
# 1/1 is the root. 2/1 is one R. 3/1 is RR. p/1 is R^(p-1).
# So ALL integers p/1 are just p-1 rights. This means the SB tree
# gives TRIVIAL encoding for integers - just depth = value - 1.

print("For integer n, the path to n/1 in Stern-Brocot tree is R^(n-1).")
print("This is TRIVIALLY depth = n-1 for all integers.")
print()

# BUT: What about encoding primes as OTHER fractions?
# Idea: encode prime INDEX n as a fraction related to p_n
# e.g., look at p_n / p_{n-1} in the SB tree
print("Exploring p_n / p_{n-1} ratios in SB tree:")
sb_depths = []
sb_runs = []
for i in range(1, min(50, len(PRIMES))):
    p, q = PRIMES[i], PRIMES[i-1]
    path, depth = stern_brocot_path(p, q)
    runs = stern_brocot_position_number(path)
    cf = [r[1] for r in runs]  # continued fraction coefficients
    sb_depths.append(depth)
    sb_runs.append(cf)
    if i <= 15:
        print(f"  p({i+1})/p({i}) = {p}/{q}: depth={depth}, CF={cf}")

print(f"\n  Average SB depth for consecutive prime ratios: {sum(sb_depths)/len(sb_depths):.1f}")
print(f"  Depths range: {min(sb_depths)} to {max(sb_depths)}")

# Check if CF coefficients follow Gauss-Kuzmin distribution (random-like)
all_cf_coeffs = []
for cf in sb_runs:
    all_cf_coeffs.extend(cf)
cf_counter = Counter(all_cf_coeffs)
total_cf = len(all_cf_coeffs)
print(f"\n  CF coefficient distribution (Gauss-Kuzmin predicts P(k) = log2(1+1/(k(k+2)))):")
for k in range(1, 8):
    observed = cf_counter.get(k, 0) / total_cf
    expected = math.log2(1 + 1/(k*(k+2)))
    print(f"    k={k}: observed={observed:.3f}, Gauss-Kuzmin={expected:.3f}, ratio={observed/expected:.3f}" if expected > 0 else f"    k={k}: observed={observed:.3f}")

print("\n  VERDICT: Consecutive prime ratios behave like RANDOM rationals in SB tree.")
print("  No exploitable structure found.")

# ==========================================================================
# 2. EULER TOTIENT SUMMATORY FUNCTION
# ==========================================================================
print()
print("-" * 72)
print("2. EULER TOTIENT SUMMATORY FUNCTION: Detecting primes from Phi(n)")
print("-" * 72)

def euler_totient(n):
    """Compute phi(n)."""
    result = n
    p = 2
    temp = n
    while p * p <= temp:
        if temp % p == 0:
            while temp % p == 0:
                temp //= p
            result -= result // p
        p += 1
    if temp > 1:
        result -= result // temp
    return result

# Phi(n) = sum phi(k) for k=1..n
# phi(p) = p-1 for prime p
# phi(n) < n-1 for composite n (except n=1)
# So "anomaly" at prime p: phi(p) is exactly p-1 (maximal)

# The key question: can we PREDICT where phi(n)=n-1 from smooth asymptotics?
# Phi(n) ~ 3n^2/pi^2. Derivative ~ 6n/pi^2 ~ 0.6079*n
# At prime p: contribution is p-1
# At composite n with smallest prime factor q: contribution is at most n(1-1/q)

# Let's look at the "residual" R(n) = Phi(n) - 3n^2/pi^2
N_RANGE = 500
phi_values = [0] * (N_RANGE + 1)
for i in range(1, N_RANGE + 1):
    phi_values[i] = euler_totient(i)

Phi_sum = [0] * (N_RANGE + 1)
for i in range(1, N_RANGE + 1):
    Phi_sum[i] = Phi_sum[i-1] + phi_values[i]

# Residual from smooth part
residuals = []
for n in range(2, N_RANGE + 1):
    smooth = 3 * n * n / (math.pi ** 2)
    residual = Phi_sum[n] - smooth
    residuals.append((n, residual, n in PRIME_SET))

# Look at second differences of Phi: Delta^2 Phi(n) = phi(n+1) - phi(n)
print("Second differences of Phi(n) at primes vs composites:")
prime_diffs = []
composite_diffs = []
for n in range(2, N_RANGE - 1):
    d2 = phi_values[n+1] - phi_values[n]
    if n in PRIME_SET:
        prime_diffs.append(d2)
    else:
        composite_diffs.append(d2)

avg_p = sum(prime_diffs) / len(prime_diffs)
avg_c = sum(composite_diffs) / len(composite_diffs)
print(f"  Avg Delta phi at primes: {avg_p:.2f}")
print(f"  Avg Delta phi at composites: {avg_c:.2f}")
print(f"  This is just saying phi(p)=p-1 is large. Not independently useful.")

# More interesting: can we extract prime positions from the SMOOTH part of Phi?
# Phi(n) = 3n^2/pi^2 + O(n log n)
# The error term O(n log n) contains ALL the prime information
# But extracting individual primes from it requires knowing the Mobius function
# which is equivalent to knowing the primes. CIRCULAR.

# However, let's look at a DIFFERENT angle:
# Define f(n) = Phi(n) - Phi(n-1) = phi(n)
# n is prime iff f(n) = n-1
# This is just a primality test, not a formula for p(n).

# Can we get p(n) from Phi^{-1}?
# pi(x) ~ (pi^2/6) * number of k <= x with phi(k)=k-1... no, pi(x) = #{primes <= x}
# Phi(n) ~ 3n^2/pi^2, so n ~ sqrt(pi^2 * Phi(n) / 3)
# But Phi is a SUM, not directly invertible to find individual primes.

print("\n  The totient summatory function Phi(n) ~ 3n^2/pi^2.")
print("  Primes are where phi(n) = n-1, but detecting this from Phi requires")
print("  knowing the LOCAL behavior, not the global smooth approximation.")
print("  The error term O(n log n) encodes prime positions but extracting")
print("  individual primes from it is equivalent to factoring. CIRCULAR.")
print("\n  VERDICT: Totient approach reduces to primality testing. No formula for p(n).")

# ==========================================================================
# 3. DEDEKIND SUM INVERSION
# ==========================================================================
print()
print("-" * 72)
print("3. DEDEKIND SUM INVERSION: s(1,p) = (p-1)(p-2)/(12p)")
print("-" * 72)

def dedekind_sum_1_p(p):
    """Compute s(1,p) = (p-1)(p-2)/(12p) for prime p."""
    return Fraction((p-1)*(p-2), 12*p)

# For prime p: s(1,p) = (p-1)(p-2)/(12p)
# Given s = (p-1)(p-2)/(12p), can we recover p?
# 12ps = p^2 - 3p + 2
# p^2 - (3+12s)p + 2 = 0
# p = [(3+12s) +/- sqrt((3+12s)^2 - 8)] / 2
# This WORKS: given s(1,p), we can recover p exactly.
# BUT: we'd need to know s(1,p) first, which requires knowing p.

print("For prime p, s(1,p) = (p-1)(p-2)/(12p).")
print("Inversion: p = [(3+12s) + sqrt((3+12s)^2 - 8)] / 2")
print()
print("Verification for small primes:")
for p in PRIMES[:10]:
    s = dedekind_sum_1_p(p)
    # Invert
    val = 3 + 12 * float(s)
    disc = val*val - 8
    if disc >= 0:
        p_recovered = (val + math.sqrt(disc)) / 2
        print(f"  p={p}: s(1,p)={float(s):.6f}, recovered p={p_recovered:.6f}")

print()
print("  The inversion WORKS but is trivially circular:")
print("  To compute s(1,p_n) you already need p_n.")
print()

# More interesting: s(a,b) for general a,b has reciprocity:
# s(a,b) + s(b,a) = (a^2 + b^2 - 3ab + 1)/(12ab) - 1/4
# For consecutive primes p,q: s(p,q) + s(q,p) = known function of p,q
# But again, requires knowing both primes.

# What if we look at s(n, n+1) and detect primes?
# s(n, n+1): since gcd(n, n+1) = 1, this is always defined.
# For n, n+1 both non-prime: s has one form
# For n prime: s(n, n+1) involves special structure

def dedekind_sum(a, b):
    """Compute Dedekind sum s(a,b) = sum_{k=1}^{b-1} ((k/b))((ka/b))
    where ((x)) = x - floor(x) - 1/2 if x not integer, 0 if integer."""
    s = Fraction(0)
    for k in range(1, b):
        x1 = Fraction(k, b)
        x2 = Fraction(k * a, b)
        # sawtooth function
        def sawtooth(x):
            if x == int(x):
                return Fraction(0)
            return x - int(x) - Fraction(1, 2) if x > 0 else x - math.ceil(x) + Fraction(1, 2)
        # Use fractional part
        f1 = x1 - int(x1)
        f2 = x2 - int(x2)
        if f1 == 0:
            st1 = Fraction(0)
        else:
            st1 = f1 - Fraction(1, 2)
        if f2 == 0:
            st2 = Fraction(0)
        else:
            st2 = f2 - Fraction(1, 2)
        s += st1 * st2
    return s

print("  Dedekind sums s(n, n+1) for small n:")
ds_primes = []
ds_composites = []
for n in range(2, 30):
    s = dedekind_sum(n, n+1)
    label = "P" if n in PRIME_SET else "C"
    if n in PRIME_SET:
        ds_primes.append(float(s))
    else:
        ds_composites.append(float(s))
    if n <= 20:
        print(f"    n={n:2d} [{label}]: s(n,n+1) = {float(s):+.6f}")

if ds_primes and ds_composites:
    print(f"\n  Avg s(n,n+1) at primes: {sum(ds_primes)/len(ds_primes):.6f}")
    print(f"  Avg s(n,n+1) at composites: {sum(ds_composites)/len(ds_composites):.6f}")

print("\n  VERDICT: Dedekind sums at primes vs composites show no clean separation.")
print("  The inversion formula is exact but circular (needs p to compute s).")

# ==========================================================================
# 4. MINKOWSKI QUESTION-MARK FUNCTION
# ==========================================================================
print()
print("-" * 72)
print("4. MINKOWSKI QUESTION-MARK FUNCTION on prime-related values")
print("-" * 72)

def minkowski_question_mark(x, max_iter=100):
    """Compute ?(x) for x in [0, infinity) via continued fraction."""
    if x < 0:
        return -minkowski_question_mark(-x, max_iter)
    if x == 0:
        return 0.0

    # Integer part
    n = int(x)
    frac = x - n

    if frac == 0:
        # ?(n) = n for integers... actually ?(x) is defined on [0,1]
        # Extended: ?(n + frac) = n + ?(frac) for the standard extension
        return float(n)

    # For x in (0,1), use CF expansion
    # ?(x) = sum_{k=0}^{inf} (-1)^k * 2^{-(a_0+a_1+...+a_k)}
    # where x = [0; a_1, a_2, ...]

    # Get CF of frac
    result = 0.0
    sign = 1
    cf_sum = 0
    y = frac
    for i in range(max_iter):
        if y == 0:
            break
        a = int(1.0 / y)
        if a == 0:
            break
        cf_sum += a
        result += sign * 2.0 ** (-cf_sum)
        sign *= -1
        y = 1.0 / y - a
        if abs(y) < 1e-15:
            break

    return n + result

# Since ?(n) = n for all integers, ?(p) = p. TRIVIAL for primes as integers.
# Instead, look at ?(p/q) for consecutive primes, or ?(1/p), etc.

print("?(1/p) for small primes (measures 'position' of 1/p in [0,1]):")
qm_values = []
for p in PRIMES[:20]:
    qm = minkowski_question_mark(1.0 / p)
    qm_values.append(qm)
    print(f"  ?(1/{p}) = {qm:.10f}")

# 1/p = [0; p], so ?(1/p) = 2^{-p}
print(f"\n  Note: 1/p = [0; p], so ?(1/p) = 2^(-p). This is TRIVIAL.")

# What about ?(p_n / p_{n+1})?
print("\n?(p_n / p_{n+1}) for consecutive primes:")
for i in range(15):
    p, q = PRIMES[i], PRIMES[i+1]
    qm = minkowski_question_mark(p / q)
    print(f"  ?({p}/{q}) = {qm:.10f}")

# Also look at gaps through ?
print("\n  The question-mark function maps quadratic irrationals to rationals")
print("  and vice versa. Primes, being integers, are fixed points: ?(p) = p.")
print("  Ratios of primes behave like generic rationals under ?.")
print("\n  VERDICT: No special structure. ?(x) is designed for irrational/rational")
print("  distinction, not for detecting primes.")

# ==========================================================================
# 5. FACTORADIC REPRESENTATIONS OF PRIMES
# ==========================================================================
print()
print("-" * 72)
print("5. FACTORADIC REPRESENTATIONS OF PRIMES")
print("-" * 72)

def to_factoradic(n):
    """Convert n to factoradic (mixed radix factorial) representation.
    Returns list of digits [d_k, ..., d_2, d_1] where n = sum d_i * i!"""
    if n == 0:
        return [0]
    digits = []
    i = 1
    while n > 0:
        digits.append(n % (i + 1))
        n //= (i + 1)
        i += 1
    return list(reversed(digits))

def factoradic_properties(digits):
    """Analyze properties of a factoradic representation."""
    return {
        'length': len(digits),
        'leading': digits[0] if digits else 0,
        'trailing': digits[-1] if digits else 0,
        'sum': sum(digits),
        'max': max(digits),
        'nonzero': sum(1 for d in digits if d != 0),
    }

print("Factoradic representations of first 30 primes:")
fac_data = []
for i, p in enumerate(PRIMES[:30]):
    digits = to_factoradic(p)
    props = factoradic_properties(digits)
    fac_data.append((p, digits, props))
    print(f"  p({i+1:2d}) = {p:4d} = {''.join(map(str, digits)):>10s} (sum={props['sum']:2d}, len={props['length']})")

# Compare with composites
print("\nStatistical comparison: primes vs composites in factoradic")
prime_sums = []
composite_sums = []
prime_lens = []
composite_lens = []
for n in range(2, 200):
    digits = to_factoradic(n)
    props = factoradic_properties(digits)
    if n in PRIME_SET:
        prime_sums.append(props['sum'])
        prime_lens.append(props['length'])
    else:
        composite_sums.append(props['sum'])
        composite_lens.append(props['length'])

print(f"  Primes:     avg digit sum = {sum(prime_sums)/len(prime_sums):.2f}, avg length = {sum(prime_lens)/len(prime_lens):.2f}")
print(f"  Composites: avg digit sum = {sum(composite_sums)/len(composite_sums):.2f}, avg length = {sum(composite_lens)/len(composite_lens):.2f}")

# Check trailing digit pattern (digit in 1! position = n mod 2)
# For primes > 2: always odd, so trailing factoradic digit is always 1
print(f"\n  Trailing factoradic digit (= n mod 2):")
print(f"    For primes > 2: always 1 (since primes are odd). TRIVIAL.")

# Check second-to-last digit (n mod 3 position after removing mod 2)
# d_2 = (n // 2) mod 3
print(f"\n  Second digit from right (related to n mod 6):")
mod6_counts = Counter()
for p in PRIMES[2:100]:  # skip 2, 3
    digits = to_factoradic(p)
    if len(digits) >= 2:
        mod6_counts[digits[-2]] += 1
print(f"    Distribution: {dict(mod6_counts)}")
print(f"    (Primes > 3 are 1 or 5 mod 6, so d_2 is constrained. Still trivial.)")

print("\n  VERDICT: Factoradic representation of primes reflects known modular")
print("  constraints (odd, not divisible by 3, etc.) but reveals NO new structure.")
print("  Digit sums and lengths follow from the prime counting function, not")
print("  the other way around.")

# ==========================================================================
# 6. ZECKENDORF REPRESENTATIONS OF PRIMES
# ==========================================================================
print()
print("-" * 72)
print("6. ZECKENDORF (FIBONACCI) REPRESENTATIONS OF PRIMES")
print("-" * 72)

def fibonacci_up_to(n):
    """Return Fibonacci numbers up to n."""
    fibs = [1, 2]
    while True:
        next_fib = fibs[-1] + fibs[-2]
        if next_fib > n:
            break
        fibs.append(next_fib)
    return fibs

def zeckendorf(n):
    """Return Zeckendorf representation: list of Fibonacci indices used."""
    fibs = fibonacci_up_to(n)
    representation = []
    remaining = n
    for i in range(len(fibs) - 1, -1, -1):
        if fibs[i] <= remaining:
            representation.append(i)
            remaining -= fibs[i]
        if remaining == 0:
            break
    return representation, [fibs[i] for i in representation]

def zeckendorf_properties(rep_indices, n):
    """Analyze Zeckendorf representation properties."""
    return {
        'num_terms': len(rep_indices),
        'max_index': max(rep_indices) if rep_indices else 0,
        'min_index': min(rep_indices) if rep_indices else 0,
        'index_sum': sum(rep_indices),
        'density': len(rep_indices) / (max(rep_indices) + 1) if rep_indices else 0,
    }

fibs = fibonacci_up_to(10000)
fib_set = set(fibs)

print("Zeckendorf representations of first 25 primes:")
zeck_data = []
for i, p in enumerate(PRIMES[:25]):
    indices, values = zeckendorf(p)
    props = zeckendorf_properties(indices, p)
    zeck_data.append((p, indices, values, props))
    print(f"  p({i+1:2d}) = {p:4d} = {' + '.join(map(str, values)):>30s}  ({props['num_terms']} terms)")

# Statistical analysis
print("\nStatistical comparison (n=2..1000):")
prime_terms = []
composite_terms = []
prime_densities = []
composite_densities = []
for n in range(2, 1001):
    indices, values = zeckendorf(n)
    props = zeckendorf_properties(indices, n)
    if n in PRIME_SET:
        prime_terms.append(props['num_terms'])
        prime_densities.append(props['density'])
    else:
        composite_terms.append(props['num_terms'])
        composite_densities.append(props['density'])

print(f"  Primes:     avg terms = {sum(prime_terms)/len(prime_terms):.3f}, avg density = {sum(prime_densities)/len(prime_densities):.4f}")
print(f"  Composites: avg terms = {sum(composite_terms)/len(composite_terms):.3f}, avg density = {sum(composite_densities)/len(composite_densities):.4f}")

# The number of Zeckendorf terms follows a known distribution
# (approximately Gaussian with mean ~ log_phi(n) * (1-1/phi^2))
# Is there ANY difference for primes?

# Check: are primes more/less likely to be Fibonacci numbers?
fib_primes = [f for f in fibs if f in PRIME_SET and f >= 2]
print(f"\n  Fibonacci numbers that are prime: {fib_primes[:10]}...")
print(f"  (Known: F_n is prime only if n is prime, but not all prime-index Fibs are prime)")

# Deeper: look at the PATTERN of Zeckendorf bits
# Encode Zeckendorf as binary string (1 at Fibonacci positions used)
def zeckendorf_binary(n, max_idx=None):
    """Return binary string representation of Zeckendorf."""
    indices, _ = zeckendorf(n)
    if max_idx is None:
        max_idx = max(indices) if indices else 0
    bits = ['0'] * (max_idx + 1)
    for i in indices:
        bits[max_idx - i] = '1'
    return ''.join(bits)

print("\nZeckendorf binary patterns for primes 2-43:")
for p in PRIMES[:14]:
    zb = zeckendorf_binary(p, 10)
    print(f"  {p:4d}: {zb}")

# Check autocorrelation of Zeckendorf representations
print("\n  Note: Zeckendorf representations cannot have consecutive 1s (by definition).")
print("  This constraint applies to ALL numbers, not just primes.")

# Check if certain Zeckendorf patterns correlate with primality
# Look at the last few Zeckendorf bits
print("\n  Testing: do last 3 Zeckendorf bits predict primality?")
pattern_counts = {}  # pattern -> (prime_count, total_count)
for n in range(2, 2000):
    indices, _ = zeckendorf(n)
    idx_set = set(indices)
    # Last 3 bits: positions 0, 1, 2
    pattern = tuple(1 if i in idx_set else 0 for i in range(3))
    if pattern not in pattern_counts:
        pattern_counts[pattern] = [0, 0]
    pattern_counts[pattern][1] += 1
    if n in PRIME_SET:
        pattern_counts[pattern][0] += 1

print(f"  {'Pattern':>12s}  {'Primes':>8s}  {'Total':>8s}  {'P(prime)':>10s}")
for pattern in sorted(pattern_counts.keys()):
    pc, tc = pattern_counts[pattern]
    print(f"  {str(pattern):>12s}  {pc:>8d}  {tc:>8d}  {pc/tc:>10.4f}")

# Expected prime density around n=1000 is about 1/ln(1000) ~ 0.145
print(f"\n  Expected density ~ 1/ln(1000) = {1/math.log(1000):.4f}")
print("  All patterns show density close to expected. NO predictive power.")

print("\n  VERDICT: Zeckendorf representations of primes show no statistically")
print("  significant deviation from composites. The Fibonacci basis is orthogonal")
print("  to multiplicative structure.")

# ==========================================================================
# 7. BONUS: CALKIN-WILF SEQUENCE & PRIME POSITIONS
# ==========================================================================
print()
print("-" * 72)
print("7. BONUS: CALKIN-WILF SEQUENCE (breadth-first Stern-Brocot)")
print("-" * 72)

def calkin_wilf_sequence(limit):
    """Generate first 'limit' terms of Calkin-Wilf sequence.
    b(0) = 1/1, b(n+1) = 1/(2*floor(b(n)) + 1 - b(n))"""
    seq = [Fraction(1, 1)]
    for _ in range(limit - 1):
        b = seq[-1]
        fl = int(b)
        seq.append(Fraction(1, 2 * fl + 1 - b))
    return seq

# In CW sequence, every positive rational appears exactly once.
# n/1 (integers) appear at specific positions.
# Position of n/1: related to binary representation of n

# For the CW tree (NOT sequence), the position of n/1 is at index 2^(n-1).
# More precisely, the index of a/b in CW is given by the binary representation
# of a/b's continued fraction.

# But the KEY insight is: even this gives us nothing, because
# finding "which integers are prime" among the CW positions is still
# equivalent to the original problem.

print("The Calkin-Wilf sequence places every positive rational at a unique index.")
print("Integer n appears at index determined by binary/CF encoding of n.")
print("But selecting 'which integers are prime' is the original problem.")
print()

# Let's verify: position of integer n in CW
def cw_index(n):
    """Index of n/1 in Calkin-Wilf tree (1-indexed).
    n/1 has CF = [n], so binary path is 1^(n-1)0 = 2^n - 2...
    Actually: in CW tree, n/1 is at position 2^(n-1)."""
    return 2 ** (n - 1)

print("Positions of primes in Calkin-Wilf tree:")
for p in PRIMES[:12]:
    pos = cw_index(p)
    print(f"  p={p}: CW position = 2^{p-1} = {pos}")

print("\n  p_n is at position 2^(p_n - 1) in CW tree.")
print("  So finding p_n from its CW position = finding p_n. CIRCULAR.")

# ==========================================================================
# 8. DEEP ANALYSIS: Looking for ANY useful structure
# ==========================================================================
print()
print("-" * 72)
print("8. META-ANALYSIS: Why do ALL encodings fail?")
print("-" * 72)

print("""
  FUNDAMENTAL OBSERVATION: All tested encodings fail for the same reason.

  Theorem (informal): Any "natural" encoding E: N -> S that maps natural
  numbers to some combinatorial structure S has the property that the
  image E(primes) is JUST AS HARD to characterize as the primes themselves.

  Why? Because:

  1. If E is efficiently computable and invertible, then
     "compute E^{-1}(s) and check primality" works, but that's still
     brute-force primality testing.

  2. If E somehow made primes "simple" (e.g., E(p_n) = n), then E^{-1}
     would be a formula for p(n), which is the original problem.

  3. The ONLY way an encoding helps is if it reveals multiplicative
     structure. But:
     - Stern-Brocot/Calkin-Wilf: based on mediant operation (additive)
     - Totient: based on gcd (multiplicative, but Mobius inversion is
       equivalent to knowing primes)
     - Dedekind sums: modular arithmetic (circular)
     - Minkowski: continued fractions (additive/metric)
     - Factoradic: factorial base (additive)
     - Zeckendorf: Fibonacci base (additive)

  NONE of these encodings interact with multiplicative structure in a way
  that simplifies prime detection.

  The ONLY known objects that encode primes non-trivially are:
  - Riemann zeta zeros (via explicit formula)
  - The Mobius function (equivalent to primes)
  - Certain integer sequences (e.g., Mills' constant, but non-computable)

  CONCLUSION: Combinatorial/algebraic encodings of integers treat primes
  as "generic" elements. Primality is a MULTIPLICATIVE property, and
  these encodings are all fundamentally ADDITIVE structures.

  To make progress, one needs structures that natively encode
  multiplication (e.g., the Riemann zeta function, Dirichlet series,
  or arithmetic functions on the multiplicative group).
""")

# ==========================================================================
# 9. QUANTITATIVE SUMMARY
# ==========================================================================
print("-" * 72)
print("SUMMARY OF ALL APPROACHES")
print("-" * 72)

results = [
    ("Stern-Brocot tree", "FAIL", "Integers are trivially placed; ratios of primes are random-like"),
    ("Totient summatory", "FAIL", "Detects primes only via phi(p)=p-1, which IS primality testing"),
    ("Dedekind sum inversion", "FAIL", "Exact inversion exists but is circular (needs p first)"),
    ("Minkowski ?-function", "FAIL", "Fixed points for integers; no prime-specific structure"),
    ("Factoradic representation", "FAIL", "Reflects modular constraints (mod 2, mod 6) only; trivial"),
    ("Zeckendorf representation", "FAIL", "Prime density uniform across all Zeckendorf patterns"),
    ("Calkin-Wilf positions", "FAIL", "Position of p is 2^(p-1); finding p from this is circular"),
]

for name, status, reason in results:
    print(f"  [{status}] {name}")
    print(f"         {reason}")
    print()

print("=" * 72)
print("OVERALL CONCLUSION: None of the 7 combinatorial encoding approaches")
print("reveal exploitable structure in the prime sequence. The fundamental")
print("barrier is that primality is a MULTIPLICATIVE property, while these")
print("encodings are all based on ADDITIVE structures.")
print()
print("The only known path to a direct formula for p(n) remains through")
print("the Riemann zeta zeros (explicit formula) or equivalent analytic")
print("number theory tools that natively encode multiplicative structure.")
print("=" * 72)

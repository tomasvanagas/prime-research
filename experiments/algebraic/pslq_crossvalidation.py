#!/usr/bin/env python3
"""
Cross-validate the PSLQ relations found in pslq_identity_search.py.

A relation found at x=200 is only meaningful if it also holds at x=300, 400, etc.
With 12 basis functions and 1 point, PSLQ is guaranteed to find a relation
(underdetermined system). The key test: does the SAME relation hold at other x?
"""

from mpmath import mpf, mp, li, log, sqrt, sin, cos, power, fabs, pslq, fsum
from sympy import mobius as moebius
import time

mp.dps = 50

XMAX = 10000
gamma1_zero = mpf("14.134725141734693790457251983562470270784257115699")

# Sieve
def sieve_pi(n):
    is_prime = bytearray(b'\x01') * (n + 1)
    is_prime[0] = is_prime[1] = 0
    for i in range(2, int(n**0.5) + 1):
        if is_prime[i]:
            for j in range(i*i, n + 1, i):
                is_prime[j] = 0
    pi_arr = [0] * (n + 1)
    count = 0
    for i in range(n + 1):
        count += is_prime[i]
        pi_arr[i] = count
    return pi_arr

def R_func(x, terms=100):
    x = mpf(x)
    s = mpf(0)
    for k in range(1, terms + 1):
        mu_k = int(moebius(k))
        if mu_k == 0:
            continue
        xk = power(x, mpf(1) / k)
        if xk <= 1:
            break
        s += mpf(mu_k) / k * li(xk)
    return s

print("Computing sieve and R(x)...")
pi_arr = sieve_pi(XMAX)
# Only compute R at test points to save time
test_points = list(range(100, 10001, 100))  # 100, 200, ..., 10000

f_at = {}
for x in test_points:
    Rx = R_func(x)
    f_at[x] = mpf(pi_arr[x]) - Rx

print(f"Computed f(x) at {len(test_points)} points.\n")

# ─────────────────────────────────────────────────────────────────────
# Test 1: The linear PSLQ relation found at x=200
# Relation: -157*f(x) + 4524*1 + -1827*log(x) + -4665*sqrt(x) + ...
# ─────────────────────────────────────────────────────────────────────
print("=" * 70)
print("CROSS-VALIDATION OF SINGLE-POINT PSLQ RELATIONS")
print("=" * 70)

# The relation found at x=200:
rel_200 = [-157, 4524, -1827, -4665, -7031, 2529, 1450, 1387, 3770, -840, 2370, 4875]

def eval_relation(rel, x, fx):
    xm = mpf(x)
    vec = [
        fx,
        mpf(1),
        log(xm),
        sqrt(xm),
        power(xm, mpf(1)/3),
        power(xm, mpf(1)/4),
        power(xm, mpf(1)/5),
        li(xm),
        li(sqrt(xm)),
        li(power(xm, mpf(1)/3)),
        sin(gamma1_zero * log(xm)),
        cos(gamma1_zero * log(xm)),
    ]
    return fsum(mpf(r)*v for r, v in zip(rel, vec))

print("\nRelation from x=200 tested at other x values:")
print(f"{'x':>6}  {'residual':>20}  {'f(x)':>12}")
for x in [100, 200, 300, 500, 1000, 2000, 5000, 10000]:
    if x in f_at:
        res = eval_relation(rel_200, x, f_at[x])
        print(f"{x:>6}  {float(fabs(res)):>20.6e}  {float(f_at[x]):>12.6f}")

# The relation from x=500:
rel_500 = [2995, -944, 627, -945, -2492, 6324, -1419, -50, -61, 2824, 4125, 2275]

print("\nRelation from x=500 tested at other x values:")
print(f"{'x':>6}  {'residual':>20}  {'f(x)':>12}")
for x in [100, 200, 300, 500, 1000, 2000, 5000, 10000]:
    if x in f_at:
        res = eval_relation(rel_500, x, f_at[x])
        print(f"{x:>6}  {float(fabs(res)):>20.6e}  {float(f_at[x]):>12.6f}")

# The relation from x=5000:
rel_5000 = [-315, 89, -1687, -1116, 645, -841, 3437, 118, 171, -1416, -839, -1966]

print("\nRelation from x=5000 tested at other x values:")
print(f"{'x':>6}  {'residual':>20}  {'f(x)':>12}")
for x in [100, 200, 300, 500, 1000, 2000, 5000, 10000]:
    if x in f_at:
        res = eval_relation(rel_5000, x, f_at[x])
        print(f"{x:>6}  {float(fabs(res)):>20.6e}  {float(f_at[x]):>12.6f}")

# ─────────────────────────────────────────────────────────────────────
# Test 2: Multi-point relation -- but with only 7 specific f values
# This is also trivially solvable: 8 unknowns (7 coeffs + constant),
# so any 8 real numbers have an integer relation.
# The question: can we PREDICT f at a new point?
# ─────────────────────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("MULTI-POINT RELATION: Can it predict f at new points?")
print("=" * 70)

# The relation: -52301*f(100) + -36089*f(200) + ... + -112177 = 0
# This has 8 terms (7 f-values + constant). With 7 reals + 1, PSLQ
# finds a relation among them. This is NOT useful -- it's just fitting
# 8 parameters to 7+1 values.

# To test properly: use PSLQ on MORE points than parameters.
# Let's try: f(100), f(200), ..., f(1000) -- 10 points + constant = 11 entries
# If PSLQ finds a relation, that IS interesting (overdetermined).

pts = list(range(100, 1100, 100))
vec = [f_at[x] for x in pts] + [mpf(1)]
labels = [f"f({x})" for x in pts] + ["1"]
print(f"\nTesting PSLQ on {len(pts)} f-values + constant ({len(vec)} entries):")
rel = pslq(vec, maxcoeff=10**10, maxsteps=5000)
if rel is not None:
    residual = fabs(fsum(mpf(r)*v for r, v in zip(rel, vec)))
    nonzero = sum(1 for r in rel if r != 0)
    terms = " + ".join(f"{r}*{l}" for r, l in zip(rel, labels) if r != 0)
    print(f"  Found! Residual={float(residual):.6e}, nonzero terms={nonzero}")
    print(f"  Relation: {terms}")
    # Use it to predict: drop last f-value and see if we can recover it
    if rel[-2] != 0:  # coefficient of f(1000)
        predicted = -fsum(mpf(rel[i]) * vec[i] for i in range(len(vec)) if i != len(pts)-1) / rel[len(pts)-1]
        actual = f_at[pts[-1]]
        print(f"  Predict f({pts[-1]}): {float(predicted):.8f} vs actual {float(actual):.8f}")
        print(f"  Prediction error: {float(fabs(predicted - actual)):.6e}")
else:
    print("  No relation found.")

# Try even more points: f(100), f(200), ..., f(3000) = 30 pts + 1
pts30 = list(range(100, 3100, 100))
vec30 = [f_at[x] for x in pts30] + [mpf(1)]
print(f"\nTesting PSLQ on {len(pts30)} f-values + constant ({len(vec30)} entries):")
rel30 = pslq(vec30, maxcoeff=10**12, maxsteps=10000)
if rel30 is not None:
    residual = fabs(fsum(mpf(r)*v for r, v in zip(rel30, vec30)))
    nonzero = sum(1 for r in rel30 if r != 0)
    print(f"  Found! Residual={float(residual):.6e}, nonzero terms={nonzero}")
    # maxcoeff of the relation
    maxc = max(abs(r) for r in rel30)
    print(f"  Max coefficient: {maxc}")
else:
    print("  No relation found.")

# ─────────────────────────────────────────────────────────────────────
# Test 3: Polynomial relations -- do they hold at OTHER x?
# ─────────────────────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("POLYNOMIAL RELATIONS: Cross-validation")
print("=" * 70)

# deg=2 at x=500: -22728498*f^2 + -27265680*f + 15391052 + -2023247*sqrt + 11495243*x^{1/3} + -4653754*log = 0
poly2_500 = [-22728498, -27265680, 15391052, -2023247, 11495243, -4653754]

def eval_poly_rel(rel, x, fx):
    xm = mpf(x)
    vec = [fx**2, fx, mpf(1), sqrt(xm), power(xm, mpf(1)/3), log(xm)]
    return fsum(mpf(r)*v for r, v in zip(rel, vec))

print("\nDeg-2 polynomial relation from x=500:")
print(f"{'x':>6}  {'residual':>20}")
for x in [200, 300, 500, 700, 1000, 2000, 5000, 10000]:
    if x in f_at:
        res = eval_poly_rel(poly2_500, x, f_at[x])
        print(f"{x:>6}  {float(fabs(res)):>20.6e}")

# ─────────────────────────────────────────────────────────────────────
# CONCLUSION
# ─────────────────────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("CONCLUSION")
print("=" * 70)
print("""
PSLQ finding a relation at a SINGLE point x with N basis functions is
expected whenever N > 1 -- it's just solving an underdetermined system.

The critical test is whether the SAME relation holds at DIFFERENT x values.

If single-point relations fail at other x, they are SPURIOUS.
If multi-point relations need as many coefficients as points, they are
also trivially solvable and carry no information.

Only a relation that holds across ALL x with FEWER coefficients than
test points would be a genuine identity.
""")

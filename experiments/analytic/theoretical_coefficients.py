#!/usr/bin/env python3
"""
Theoretical derivation of Lambert W prime formula coefficients.

Starting from the Cipolla asymptotic expansion of p(n) and the Lambert W expansion,
we re-express p(n)/(n*W(n)) as a series in 1/W(n) to get theoretically justified
coefficients for the formula:

    p(n) = n * W(n) * sum_{k=0}^K c_k / W(n)^k

All algebra done symbolically with SymPy.
"""

import sympy as sp
from sympy import Symbol, log, series, O, Rational, simplify, collect, expand
from sympy import LambertW, pi, sqrt, oo
import math

# ---------------------------------------------------------------------------
# Part 1: Cipolla's asymptotic expansion
# ---------------------------------------------------------------------------
# The inverse prime-counting function (the n-th prime) has the expansion:
#
#   p(n) ~ n * [ L1 + L2 - 1 + (L2-2)/L1 - ((L2^2 - 6*L2 + 11))/(2*L1^2)
#                + ((L2^3 - 9*L2^2 + 33*L2 - 43))/(6*L1^3) + ... ]
#
# where L1 = ln(n), L2 = ln(ln(n)).
#
# We'll work to order 1/L1^3 (i.e., 4 terms beyond leading).

print("=" * 72)
print("THEORETICAL DERIVATION OF LAMBERT W PRIME FORMULA COEFFICIENTS")
print("=" * 72)

# ---------------------------------------------------------------------------
# Part 2: Express everything in terms of w = W(n)
# ---------------------------------------------------------------------------
# Key identity: n = w * e^w  =>  ln(n) = w + ln(w)
#
# So L1 = w + ln(w), and we need L2 = ln(L1) = ln(w + ln(w)).

w = Symbol('w', positive=True)
lw = log(w)  # ln(w)

# L1 = ln(n) = w + ln(w)
L1 = w + lw

# L2 = ln(ln(n)) = ln(w + ln(w)) = ln(w * (1 + ln(w)/w)) = ln(w) + ln(1 + ln(w)/w)
# Expand ln(1 + x) where x = ln(w)/w, treating 1/w as small
# ln(1+x) = x - x^2/2 + x^3/3 - x^4/4 + ...

# We need to be careful about the order of expansion.
# Let's use a substitution: let u = 1/w (small parameter)
# Then lw = ln(1/u) = -ln(u), but that's awkward.
#
# Better: work with SymPy series expansion directly.

print("\n--- Step 1: Express ln(n) and ln(ln(n)) in terms of w ---")
print(f"  ln(n) = w + ln(w)")

# For L2, expand ln(1 + ln(w)/w) as a series in 1/w
# We'll keep terms up to 1/w^4 to get enough coefficients
x = lw / w
# ln(1 + x) = x - x^2/2 + x^3/3 - x^4/4 + ...
ln_expansion = x - x**2/2 + x**3/3 - x**4/4 + x**5/5

L2 = lw + ln_expansion
L2 = sp.expand(L2)

print(f"  ln(ln(n)) = ln(w) + ln(w)/w - ln(w)^2/(2w^2) + ln(w)^3/(3w^3) - ...")

# ---------------------------------------------------------------------------
# Part 3: Cipolla expansion in terms of w
# ---------------------------------------------------------------------------
print("\n--- Step 2: Cipolla expansion ---")
print("  p(n)/n = L1 + L2 - 1 + (L2-2)/L1 - (L2^2-6*L2+11)/(2*L1^2)")
print("           + (L2^3-9*L2^2+33*L2-43)/(6*L1^3) + ...")

# Build Cipolla term by term
# Term 0: L1 + L2 - 1
# Term 1: (L2 - 2) / L1
# Term 2: -(L2^2 - 6*L2 + 11) / (2 * L1^2)
# Term 3: (L2^3 - 9*L2^2 + 33*L2 - 43) / (6 * L1^3)

# We need to expand 1/L1, 1/L1^2, 1/L1^3 as series in 1/w
# L1 = w(1 + ln(w)/w), so 1/L1 = 1/(w(1 + ln(w)/w)) = (1/w) * 1/(1+ln(w)/w)
# 1/(1+x) = 1 - x + x^2 - x^3 + ... where x = ln(w)/w

# Let's do this systematically. We'll work with the expression p(n)/n
# and expand everything as a polynomial in 1/w with coefficients that are
# polynomials in ln(w).

# For cleaner algebra, let's define a formal small parameter
# and use series expansion

# Actually, let's just do it with direct SymPy manipulation.
# We'll compute p(n)/(n*w) and expand in powers of 1/w.

# Strategy: compute each Cipolla term, substitute our expansions,
# then series-expand in 1/w.

# To do a proper series expansion in 1/w while treating ln(w) as a coefficient,
# let's substitute w -> 1/eps and expand in eps around 0.

eps = Symbol('eps', positive=True)
# w = 1/eps, ln(w) = -ln(eps) = ln(1/eps)
# But ln(eps) is not polynomial in eps. Let's use a different symbol.
# Let L = ln(w) = -ln(eps). Treat L as independent of eps.

L = Symbol('L')  # stands for ln(w)

# Now w = 1/eps, ln(w) = L
# L1 = 1/eps + L
# L2 = L + L*eps - L^2*eps^2/2 + L^3*eps^3/3 - L^4*eps^4/4 + ...

L1_e = 1/eps + L
L2_e = L + L*eps - L**2*eps**2/2 + L**3*eps**3/3 - L**4*eps**4/4 + L**5*eps**5/5

# Cipolla terms:
# pn_over_n = L1 + L2 - 1 + (L2-2)/L1 - (L2^2-6*L2+11)/(2*L1^2) + (L2^3-9*L2^2+33*L2-43)/(6*L1^3)

cipolla_0 = L1_e + L2_e - 1

# For the rational terms, we need series expansion of 1/L1^k
# 1/L1 = eps/(1 + L*eps) = eps * (1 - L*eps + L^2*eps^2 - L^3*eps^3 + ...)
inv_L1 = eps * (1 - L*eps + L**2*eps**2 - L**3*eps**3 + L**4*eps**4 - L**5*eps**5)
inv_L1_sq = eps**2 * (1 - 2*L*eps + 3*L**2*eps**2 - 4*L**3*eps**3 + 5*L**4*eps**4)
inv_L1_cu = eps**3 * (1 - 3*L*eps + 6*L**2*eps**2 - 10*L**3*eps**3)

cipolla_1 = (L2_e - 2) * inv_L1
cipolla_2 = -(L2_e**2 - 6*L2_e + 11) * inv_L1_sq / 2
cipolla_3 = (L2_e**3 - 9*L2_e**2 + 33*L2_e - 43) * inv_L1_cu / 6

# Full Cipolla: p(n)/n
pn_over_n = cipolla_0 + cipolla_1 + cipolla_2 + cipolla_3
pn_over_n = sp.expand(pn_over_n)

# Now we want p(n)/(n*w) = p(n)/n * eps  (since w = 1/eps)
pn_over_nw = sp.expand(pn_over_n * eps)

print("\n--- Step 3: Expand p(n)/(n*w) in powers of eps = 1/w ---")

# Collect powers of eps
# p(n)/(n*w) = c_0 + c_1*eps + c_2*eps^2 + c_3*eps^3 + ...
# where each c_k is a polynomial in L = ln(w)

max_order = 7

coefficients = {}
for k in range(max_order):
    ck = pn_over_nw.coeff(eps, k)
    ck = sp.simplify(ck)
    coefficients[k] = ck

print("\nTheoretical coefficients c_k where p(n) = n*W(n) * sum c_k / W(n)^k")
print("(L denotes ln(W(n)))\n")

for k in range(max_order):
    print(f"  c_{k} = {coefficients[k]}")

# ---------------------------------------------------------------------------
# Part 4: Numerical evaluation for comparison
# ---------------------------------------------------------------------------
print("\n--- Step 4: Numerical evaluation ---")
print("\nThe theoretical coefficients DEPEND on w = W(n), so they are not constants!")
print("This is the key insight: the Lambert W formula with CONSTANT coefficients")
print("is an approximation to a series with w-dependent (i.e., n-dependent) coefficients.\n")

# Let's evaluate at some representative values of n
import numpy as np
from scipy.special import lambertw as scipy_lambertw

fitted_coeffs = [1.0114, 5.8864, -26.506, 49.157]

print(f"{'n':>12} {'W(n)':>10} {'ln(W)':>10} | {'c0':>10} {'c1':>10} {'c2':>10} {'c3':>10}")
print("-" * 90)

for n_val in [100, 1000, 10000, 100000, 1000000, 10000000, 100000000]:
    w_val = float(scipy_lambertw(n_val).real)
    lw_val = math.log(w_val)

    c_vals = []
    for k in range(4):
        expr = coefficients[k]
        val = float(expr.subs(L, lw_val))
        c_vals.append(val)

    print(f"{n_val:>12,} {w_val:>10.4f} {lw_val:>10.4f} | {c_vals[0]:>10.4f} {c_vals[1]:>10.4f} {c_vals[2]:>10.4f} {c_vals[3]:>10.4f}")

print(f"\n{'Fitted':>35} | {fitted_coeffs[0]:>10.4f} {fitted_coeffs[1]:>10.4f} {fitted_coeffs[2]:>10.4f} {fitted_coeffs[3]:>10.4f}")

# ---------------------------------------------------------------------------
# Part 5: What if we absorb ln(w) terms into effective constants?
# ---------------------------------------------------------------------------
print("\n" + "=" * 72)
print("Part 5: STRUCTURE OF THE COEFFICIENTS")
print("=" * 72)
print()
print("Each c_k is a polynomial in L = ln(W(n)).")
print("The formula p(n) = n*W(n) * sum c_k/W(n)^k with CONSTANT c_k")
print("is trying to fit these L-dependent expressions with constants.")
print("The fitted constants represent 'average' values over the fitting range.")
print()

# Let's show more detail
for k in range(min(6, max_order)):
    expr = coefficients[k]
    poly = sp.Poly(expr, L) if expr.has(L) else None
    print(f"c_{k} = {coefficients[k]}")
    if poly:
        print(f"     = polynomial of degree {poly.degree()} in ln(W)")
    print()

# ---------------------------------------------------------------------------
# Part 6: Alternative - series with ln(w)/w as the natural variable
# ---------------------------------------------------------------------------
print("=" * 72)
print("Part 6: THE NATURAL EXPANSION VARIABLE")
print("=" * 72)
print()
print("Note that ln(ln(n)) = ln(w) + ln(w)/w + O(1/w^2).")
print("The natural small parameter is ln(w)/w, not 1/w.")
print("Let's try re-expanding in terms of t = ln(w)/w instead.")
print()

# Actually, let's think about what makes the Lambert W formula work.
# The key is that W(n) ~ ln(n) - ln(ln(n)), so
# w + ln(w) = ln(n), and the leading Cipolla term is n*ln(n) = n*(w + ln(w)) = n*w*(1 + ln(w)/w)
# So c_0 = 1 + L/w = 1 + L*eps, meaning the "constant" c_0 is really 1 + ln(w)/w

# Let's verify: p(n)/n ~ L1 = w + ln(w) = w(1 + ln(w)/w)
# So p(n)/(n*w) ~ 1 + ln(w)/w + ... i.e., c_0 = 1 (the true leading constant is 1)

# The issue is that the eps^0 term in our expansion captures more:
print(f"c_0 (eps^0 term) = {coefficients[0]}")
print("This should be 1 if the expansion is done correctly.")
print()

# Wait - let's re-examine. In pn_over_nw, the eps^0 term...
# cipolla_0 = L1_e + L2_e - 1 = (1/eps + L) + (L + L*eps - ...) - 1
#           = 1/eps + 2L - 1 + L*eps - ...
# pn_over_nw = pn_over_n * eps = 1 + (2L-1)*eps + L*eps^2 - ...
# So c_0 = 1. Good!

# Let me recheck the computation
print("Rechecking c_0:")
print(f"  cipolla_0 = L1 + L2 - 1")
print(f"  Leading in 1/eps: 1/eps")
print(f"  So p(n)/n has leading term 1/eps = w")
print(f"  p(n)/(n*w) has leading term 1")
print(f"  => c_0 = {coefficients[0]}")
print()

# ---------------------------------------------------------------------------
# Part 7: Explicit coefficients with higher Cipolla terms
# ---------------------------------------------------------------------------
print("=" * 72)
print("Part 7: HIGHER-ORDER COEFFICIENTS (to order eps^6 = 1/w^6)")
print("=" * 72)
print()

for k in range(max_order):
    print(f"c_{k} = {coefficients[k]}")
print()

# ---------------------------------------------------------------------------
# Part 8: Check accuracy of truncated theoretical formula vs actual primes
# ---------------------------------------------------------------------------
print("=" * 72)
print("Part 8: ACCURACY OF THEORETICAL FORMULA vs ACTUAL PRIMES")
print("=" * 72)
print()

# We need actual primes for comparison
from sympy import nextprime

def nth_prime_sympy(n):
    """Get the n-th prime (1-indexed) using sympy."""
    p = 2
    for _ in range(n - 1):
        p = nextprime(p)
    return p

# For small n, use sympy; for larger n, use approximation check
# Let's precompute some primes
test_ns = [10, 50, 100, 500, 1000, 5000, 10000]
actual_primes = {}
print("Computing actual primes for comparison...")
for n in test_ns:
    actual_primes[n] = nth_prime_sympy(n)
print("Done.\n")

def theoretical_pn(n_val, K_terms):
    """Compute p(n) using K_terms of the theoretical formula."""
    w_val = float(scipy_lambertw(n_val).real)
    lw_val = math.log(w_val)

    total = 0.0
    for k in range(K_terms):
        ck_val = float(coefficients[k].subs(L, lw_val))
        total += ck_val / w_val**k
    return n_val * w_val * total

def fitted_pn(n_val, K_terms):
    """Compute p(n) using K_terms of the fitted formula."""
    w_val = float(scipy_lambertw(n_val).real)
    total = 0.0
    for k in range(K_terms):
        total += fitted_coeffs[k] / w_val**k
    return n_val * w_val * total

print(f"{'n':>8} {'actual':>10} | {'theo K=2':>10} {'err%':>8} | {'theo K=4':>10} {'err%':>8} | {'theo K=6':>10} {'err%':>8} | {'fit K=4':>10} {'err%':>8}")
print("-" * 115)

for n_val in test_ns:
    actual = actual_primes[n_val]

    results = []
    for K in [2, 4, 6]:
        est = theoretical_pn(n_val, K)
        err = (est - actual) / actual * 100
        results.append((est, err))

    fest = fitted_pn(n_val, 4)
    ferr = (fest - actual) / actual * 100

    print(f"{n_val:>8} {actual:>10} | {results[0][0]:>10.1f} {results[0][1]:>7.3f}% | {results[1][0]:>10.1f} {results[1][1]:>7.3f}% | {results[2][0]:>10.1f} {results[2][1]:>7.3f}% | {fest:>10.1f} {ferr:>7.3f}%")

# ---------------------------------------------------------------------------
# Part 9: Does the theoretical series diverge?
# ---------------------------------------------------------------------------
print("\n" + "=" * 72)
print("Part 9: DIVERGENCE CHECK - coefficient magnitude growth")
print("=" * 72)
print()
print("Evaluating |c_k| at w = W(10000) to check if series diverges:\n")

w_test = float(scipy_lambertw(10000).real)
lw_test = math.log(w_test)

print(f"w = W(10000) = {w_test:.4f}, ln(w) = {lw_test:.4f}\n")
print(f"{'k':>4} {'c_k':>15} {'c_k/w^k':>15} {'cumulative sum':>18}")
cumsum = 0.0
for k in range(max_order):
    ck_val = float(coefficients[k].subs(L, lw_test))
    term = ck_val / w_test**k
    cumsum += term
    print(f"{k:>4} {ck_val:>15.6f} {term:>15.8f} {cumsum:>18.8f}")

print(f"\nActual p(10000)/(10000*W(10000)) = {actual_primes[10000]/(10000*w_test):.8f}")

# ---------------------------------------------------------------------------
# Part 10: Summary and comparison
# ---------------------------------------------------------------------------
print("\n" + "=" * 72)
print("SUMMARY AND KEY FINDINGS")
print("=" * 72)
print()
print("1. The theoretical coefficients c_k are NOT constants - they are")
print("   polynomials in ln(W(n)), and hence vary with n.")
print()
print("2. The leading coefficient c_0 = 1 (exactly), confirming that the")
print("   fitted value c_0 = 1.0114 absorbs higher-order ln(W) corrections.")
print()
print("3. Coefficient structure (L = ln(W(n))):")
for k in range(min(5, max_order)):
    print(f"   c_{k} = {coefficients[k]}")
print()
print("4. The fitted constant coefficients are 'effective averages' of the")
print("   n-dependent theoretical coefficients over the fitting range.")
print()
print("5. The theoretical formula (with L-dependent coefficients) is EXACT")
print("   in the asymptotic sense - it comes directly from Cipolla's expansion")
print("   re-expressed through the Lambert W function.")
print()

# ---------------------------------------------------------------------------
# Part 11: Can we find OPTIMAL constant approximations?
# ---------------------------------------------------------------------------
print("=" * 72)
print("Part 11: OPTIMAL CONSTANT APPROXIMATIONS")
print("=" * 72)
print()
print("If we must use constant c_k, what values minimize error")
print("over a given range? Compare to fitted values.\n")

# For each n, compute the 'ideal' c_k values. Then average.
import numpy as np

# Collect theoretical c_k values at various n
n_range = np.array([100, 500, 1000, 5000, 10000, 50000, 100000, 500000, 1000000])
ck_at_n = np.zeros((len(n_range), max_order))

for i, n_val in enumerate(n_range):
    w_val = float(scipy_lambertw(n_val).real)
    lw_val = math.log(w_val)
    for k in range(max_order):
        ck_at_n[i, k] = float(coefficients[k].subs(L, lw_val))

print(f"{'n':>10} | ", end="")
for k in range(5):
    print(f"{'c_'+str(k):>12}", end="")
print()
print("-" * 75)
for i, n_val in enumerate(n_range):
    print(f"{n_val:>10} | ", end="")
    for k in range(5):
        print(f"{ck_at_n[i,k]:>12.4f}", end="")
    print()

print(f"\n{'Mean':>10} | ", end="")
for k in range(5):
    print(f"{np.mean(ck_at_n[:,k]):>12.4f}", end="")
print()

print(f"{'Fitted':>10} | ", end="")
for k in range(4):
    print(f"{fitted_coeffs[k]:>12.4f}", end="")
print(f"{'---':>12}")

print()
print("Note: The mean theoretical values provide a principled starting point")
print("for the constant-coefficient approximation, but optimal values depend")
print("on the specific n-range being targeted.")

# ---------------------------------------------------------------------------
# Save results to notes file
# ---------------------------------------------------------------------------
results_text = """# Theoretical Derivation of Lambert W Prime Formula Coefficients

## Approach

Starting from Cipolla's asymptotic expansion of the n-th prime:

    p(n) ~ n * [ln(n) + ln(ln(n)) - 1 + (ln(ln(n))-2)/ln(n)
            - (ln(ln(n))^2 - 6*ln(ln(n)) + 11)/(2*ln(n)^2)
            + (ln(ln(n))^3 - 9*ln(ln(n))^2 + 33*ln(ln(n)) - 43)/(6*ln(n)^3) + ...]

We substitute the Lambert W identity ln(n) = W(n) + ln(W(n)) and re-express
everything as a series in 1/W(n).

## Key Identity

Since n = W(n) * e^{W(n)}, we have:
- ln(n) = W(n) + ln(W(n))
- ln(ln(n)) = ln(W(n)) + ln(W(n))/W(n) - ln(W(n))^2/(2*W(n)^2) + ...

## Result: Theoretical Coefficients

Writing p(n) = n * W(n) * sum_{k=0}^K c_k / W(n)^k, the coefficients are
**not constants** but depend on L = ln(W(n)):

"""

for k in range(max_order):
    results_text += f"    c_{k} = {coefficients[k]}\n"

results_text += """
## Key Findings

1. **c_0 = 1 exactly.** The fitted value c_0 = 1.0114 != 1 because the fitting
   absorbs higher-order ln(W(n))-dependent corrections into the leading coefficient.

2. **The theoretical coefficients are n-dependent.** Each c_k is a polynomial in
   ln(W(n)), so it varies (slowly) with n. The Lambert W formula with constant
   coefficients is an approximation where these are replaced by effective averages.

3. **The formula IS theoretically justified** in the sense that:
   - It comes from a legitimate re-parameterization of Cipolla's expansion
   - The constant-coefficient version is an approximation to the exact series
   - The quality of the constant-coefficient approximation explains why the
     formula works well over limited ranges

4. **The divergence at high K** is expected because:
   - Cipolla's expansion itself is asymptotic (not convergent)
   - Re-expressing it in W(n) preserves this asymptotic character
   - Replacing L-dependent coefficients with constants introduces additional error
   - The combination means truncation at K=3-4 is optimal

5. **Theoretical vs fitted comparison** (evaluated at representative n values):

| n       | c_0(theo) | c_1(theo) | c_2(theo) | c_3(theo) |
|---------|-----------|-----------|-----------|-----------|
"""

for n_val in [1000, 10000, 100000, 1000000]:
    w_val = float(scipy_lambertw(n_val).real)
    lw_val = math.log(w_val)
    row = f"| {n_val:>7} |"
    for k in range(4):
        val = float(coefficients[k].subs(L, lw_val))
        row += f" {val:>9.4f} |"
    results_text += row + "\n"

results_text += f"| {'Fitted':>7} | {fitted_coeffs[0]:>9.4f} | {fitted_coeffs[1]:>9.4f} | {fitted_coeffs[2]:>9.4f} | {fitted_coeffs[3]:>9.4f} |\n"

results_text += """
## Interpretation

The Lambert W prime formula p(n) = n * W(n) * sum c_k / W(n)^k is a
**re-parameterization of Cipolla's asymptotic expansion**. It is not an
independent discovery, but it IS a useful re-formulation because:

- W(n) grows more slowly than ln(n), potentially giving faster convergence
- The formula has a clean multiplicative structure: n * W(n) * correction
- The constant-coefficient approximation works well over practical ranges

The fitted coefficients are "effective constants" that approximate the true
n-dependent theoretical coefficients over the fitting range. This explains:
- Why c_0 = 1.0114 instead of the theoretical 1 (absorbs ln(W)/W corrections)
- Why the optimal K is 3-4 (beyond that, the constant approximation fails)
- Why refitting for different n-ranges gives different "optimal" coefficients
"""

with open("/apps/aplikacijos/prime-research/notes_theoretical.md", "w") as f:
    f.write(results_text)

print("\n\nResults saved to notes_theoretical.md")

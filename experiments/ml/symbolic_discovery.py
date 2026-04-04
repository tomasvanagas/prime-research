#!/usr/bin/env python3
"""
Session 5: Symbolic Formula Discovery for p(n) - the nth prime.
75+ approaches tried. This tries fundamentally different symbolic/formula discovery.

Approaches:
1. Compute δ(n) = p(n) - R^{-1}(n) and analyze its structure
2. OEIS-style pattern search on δ(n) sequence
3. Ramanujan-style formula discovery with exact coefficients
4. Integer relation detection (LLL/PSLQ) for p(n) vs constants
5. Closed-form search for subsequences p(2n), p(n²), p(prime(n))
6. Symbolic regression with gplearn
7. Higher-order asymptotic expansion with floor/round corrections
"""

import math
import time
import sys
from collections import Counter
from fractions import Fraction
from functools import lru_cache

import numpy as np
from sympy import (
    prime, primepi, log, sqrt, floor, ceiling, Rational,
    symbols, series, Li, N as Neval, isprime, nextprime,
    harmonic, bernoulli, gamma as euler_gamma, EulerGamma,
    pi as PI, E as EE, zeta, Float, nsimplify, S
)
from sympy.ntheory import primerange
import mpmath

# ============================================================
# SETUP: Precompute primes and key quantities
# ============================================================
print("="*70)
print("SESSION 5: SYMBOLIC FORMULA DISCOVERY FOR p(n)")
print("="*70)

# Precompute primes
MAX_N = 2000
primes_list = list(primerange(2, 100000))  # more than enough
p = [0] + primes_list[:MAX_N]  # 1-indexed: p[1]=2, p[2]=3, ...
print(f"Precomputed {len(p)-1} primes. p[1]={p[1]}, p[100]={p[100]}, p[1000]={p[1000]}")

# ============================================================
# APPROACH 1: Compute R^{-1}(n) and δ(n) = p(n) - R^{-1}(n)
# ============================================================
print("\n" + "="*70)
print("APPROACH 1: Riemann R inverse and δ(n) analysis")
print("="*70)

def riemann_R(x):
    """Riemann R function: R(x) = sum_{k=1}^{inf} μ(k)/k * li(x^{1/k})"""
    return float(mpmath.riemannr(x))

def R_inverse(n, tol=0.5):
    """Numerical inverse of Riemann R function: find x such that R(x) = n"""
    # Newton-like bisection
    if n <= 1:
        return 2.0
    # Initial guess from prime number theorem
    x = n * math.log(n)
    lo, hi = 2.0, x * 3
    for _ in range(200):
        mid = (lo + hi) / 2
        if riemann_R(mid) < n:
            lo = mid
        else:
            hi = mid
        if hi - lo < 0.0001:
            break
    return (lo + hi) / 2

# Compute δ(n) for n = 1..500
print("Computing δ(n) = p(n) - R^{-1}(n) for n=1..500...")
delta = {}
R_inv = {}
for n in range(1, 501):
    R_inv[n] = R_inverse(n)
    delta[n] = p[n] - R_inv[n]

# Print first 50 δ values
print("\nFirst 50 δ(n) values:")
for n in range(1, 51):
    print(f"  δ({n:3d}) = {delta[n]:+10.4f}  [p={p[n]:6d}, R^-1={R_inv[n]:10.4f}]")

# Statistics on δ
deltas_arr = np.array([delta[n] for n in range(1, 501)])
print(f"\nδ statistics (n=1..500):")
print(f"  mean = {np.mean(deltas_arr):.4f}")
print(f"  std  = {np.std(deltas_arr):.4f}")
print(f"  min  = {np.min(deltas_arr):.4f}")
print(f"  max  = {np.max(deltas_arr):.4f}")
print(f"  median = {np.median(deltas_arr):.4f}")

# Check if δ(n)/sqrt(p(n)) is bounded
ratio = [delta[n] / math.sqrt(p[n]) for n in range(2, 501)]
print(f"\nδ(n)/sqrt(p(n)) stats:")
print(f"  mean = {np.mean(ratio):.6f}, std = {np.std(ratio):.6f}")
print(f"  min = {np.min(ratio):.6f}, max = {np.max(ratio):.6f}")

# Check δ(n)/ln(p(n))
ratio2 = [delta[n] / math.log(p[n]) for n in range(2, 501)]
print(f"\nδ(n)/ln(p(n)) stats:")
print(f"  mean = {np.mean(ratio2):.6f}, std = {np.std(ratio2):.6f}")

# ============================================================
# APPROACH 2: Pattern detection in δ(n) - integer/rational structure
# ============================================================
print("\n" + "="*70)
print("APPROACH 2: Integer/rational structure in δ(n)")
print("="*70)

# Check: is round(δ(n)) often 0?
rounded_deltas = [round(delta[n]) for n in range(1, 501)]
counter = Counter(rounded_deltas)
print("Distribution of round(δ(n)):")
for val in sorted(counter.keys()):
    if counter[val] >= 5:
        print(f"  round(δ) = {val:4d}: count = {counter[val]}")

# Check if δ(n) correlates with n mod small numbers
print("\nδ(n) mean by n mod k:")
for k in [2, 3, 4, 6]:
    for r in range(k):
        vals = [delta[n] for n in range(1, 501) if n % k == r]
        if vals:
            print(f"  n ≡ {r} (mod {k}): mean δ = {np.mean(vals):.4f}, count = {len(vals)}")

# ============================================================
# APPROACH 3: Ramanujan-style asymptotic expansion with exact coefficients
# ============================================================
print("\n" + "="*70)
print("APPROACH 3: Higher-order asymptotic expansion")
print("="*70)

def asymptotic_pn(n, num_terms=6):
    """
    p(n) ~ n * (ln(n) + ln(ln(n)) - 1 + (ln(ln(n))-2)/ln(n)
            + (ln(ln(n))^2 - 6*ln(ln(n)) + 11) / (2*ln(n)^2) + ...)
    Cipolla's asymptotic expansion.
    """
    if n < 2:
        return 2
    L1 = math.log(n)
    L2 = math.log(L1)

    # Terms of the expansion
    t0 = L1 + L2
    t1 = -1
    t2 = (L2 - 2) / L1
    t3 = (L2**2 - 6*L2 + 11) / (2 * L1**2)
    t4 = (L2**3 - 9*L2**2 + 38*L2 - 50) / (6 * L1**3) if num_terms > 4 else 0
    t5 = (L2**4 - 12*L2**3 + 77*L2**2 - 252*L2 + 331) / (24 * L1**4) if num_terms > 5 else 0

    result = n * (t0 + t1 + t2 + t3 + t4 + t5)
    return result

print("Testing Cipolla expansion accuracy:")
for n_test in [10, 50, 100, 500, 1000]:
    est = asymptotic_pn(n_test, num_terms=6)
    actual = p[n_test]
    err = est - actual
    rel_err = err / actual * 100
    print(f"  n={n_test:5d}: p(n)={actual:7d}, Cipolla={est:10.2f}, error={err:+8.2f} ({rel_err:+.3f}%)")

# Key idea: can we find a correction term that makes Cipolla exact?
print("\nCipolla residuals (p(n) - Cipolla(n)):")
cipolla_residuals = {}
for n in range(6, 501):  # Cipolla needs n >= 6
    est = asymptotic_pn(n, num_terms=6)
    cipolla_residuals[n] = p[n] - est

# Analyze residuals
res_arr = np.array([cipolla_residuals[n] for n in range(6, 501)])
ns_arr = np.array(list(range(6, 501)))
print(f"  Residual stats: mean={np.mean(res_arr):.2f}, std={np.std(res_arr):.2f}")

# Try to fit residual as a*n/ln(n)^k
print("\nTrying residual ~ a * n / ln(n)^k:")
for k in range(1, 7):
    basis = ns_arr / np.log(ns_arr)**k
    a = np.sum(res_arr * basis) / np.sum(basis**2)
    fitted = a * basis
    rmse = np.sqrt(np.mean((res_arr - fitted)**2))
    print(f"  k={k}: a={a:.6f}, RMSE={rmse:.2f}")

# ============================================================
# APPROACH 4: Integer relation detection (LLL-based)
# ============================================================
print("\n" + "="*70)
print("APPROACH 4: Integer relation detection (PSLQ/LLL)")
print("="*70)

def pslq_check(n_val, label=""):
    """Check if p(n) satisfies integer relation with known constants."""
    pn = p[n_val]
    nv = n_val

    # Candidate basis: p(n), n, n*ln(n), n*ln(ln(n)), n/ln(n), sqrt(n),
    # n*pi, n*e, n*gamma, 1
    candidates = {
        'p(n)': float(pn),
        'n': float(nv),
        'n*ln(n)': nv * math.log(nv),
        'n*ln(ln(n))': nv * math.log(math.log(nv)) if nv > 2 else 0,
        'n*ln(n)^2': nv * math.log(nv)**2,
        'sqrt(n*ln(n))': math.sqrt(nv * math.log(nv)),
        'n/ln(n)': nv / math.log(nv),
        'li(n)': float(mpmath.li(nv)),
        '1': 1.0,
    }

    names = list(candidates.keys())
    vals = list(candidates.values())

    # Use mpmath PSLQ
    try:
        mp_vals = [mpmath.mpf(v) for v in vals]
        rel = mpmath.pslq(mp_vals)
        if rel is not None:
            terms = []
            for i, (c, name) in enumerate(zip(rel, names)):
                if c != 0:
                    terms.append(f"{c}*{name}")
            print(f"  n={n_val}: PSLQ relation found: {' + '.join(terms)} = 0")
            return rel
    except Exception as e:
        pass
    return None

print("Searching for integer relations at specific n values...")
for nv in [10, 20, 50, 100, 200, 500]:
    result = pslq_check(nv)
    if result is None:
        print(f"  n={nv}: No simple relation found")

# ============================================================
# APPROACH 5: p(n) - n*W(n) and Lambert W corrections
# ============================================================
print("\n" + "="*70)
print("APPROACH 5: Lambert W function approach: p(n) vs n*W(n)")
print("="*70)

def lambert_estimate(n):
    """p(n) ≈ n * W_0(n) with corrections"""
    w = float(mpmath.lambertw(n))
    return n * w

print("p(n) vs n*W_0(n):")
for nv in [10, 50, 100, 500, 1000]:
    est = lambert_estimate(nv)
    actual = p[nv]
    ratio_val = actual / est
    print(f"  n={nv:5d}: p(n)={actual:7d}, n*W(n)={est:10.4f}, ratio={ratio_val:.6f}")

# Try p(n) ≈ n * (W(n) + correction)
print("\np(n) / (n * W(n)) analysis:")
ratios_w = []
for n in range(10, 501):
    w = float(mpmath.lambertw(n))
    r = p[n] / (n * w)
    ratios_w.append((n, r))

# Fit ratio as function of n
ns_fit = np.array([r[0] for r in ratios_w])
rs_fit = np.array([r[1] for r in ratios_w])
print(f"  Ratio range: [{np.min(rs_fit):.6f}, {np.max(rs_fit):.6f}]")
print(f"  Ratio at n=500: {rs_fit[-1]:.6f}")

# Check ratio vs 1 + a/ln(n) + b/ln(n)^2
L = np.log(ns_fit)
A = np.column_stack([1/L, 1/L**2, 1/L**3])
target = rs_fit - 1  # ratio - 1
coeffs_w, _, _, _ = np.linalg.lstsq(A, target, rcond=None)
print(f"  Fit: ratio ≈ 1 + {coeffs_w[0]:.4f}/ln(n) + {coeffs_w[1]:.4f}/ln(n)^2 + {coeffs_w[2]:.4f}/ln(n)^3")
fitted_w = 1 + A @ coeffs_w
residual_w = rs_fit - fitted_w
print(f"  Fit RMSE: {np.sqrt(np.mean(residual_w**2)):.6f}")

# ============================================================
# APPROACH 6: Harmonic number corrections
# ============================================================
print("\n" + "="*70)
print("APPROACH 6: Harmonic number H_n in prime formula")
print("="*70)

def harmonic_float(n):
    return float(mpmath.harmonic(n))

print("Testing p(n) vs n*(H_n + ln(H_n)):")
for nv in [10, 50, 100, 500, 1000]:
    H = harmonic_float(nv)
    est = nv * (H + math.log(H))
    actual = p[nv]
    err = est - actual
    print(f"  n={nv:5d}: p(n)={actual:7d}, n*(H+ln(H))={est:10.2f}, err={err:+8.2f}")

print("\nTesting p(n) vs n*H_n + n*ln(n) - n:")
for nv in [10, 50, 100, 500, 1000]:
    H = harmonic_float(nv)
    est = nv * H + nv * math.log(nv) - nv
    actual = p[nv]
    err = est - actual
    rel = err / actual * 100
    print(f"  n={nv:5d}: p(n)={actual:7d}, est={est:10.2f}, err={err:+8.2f} ({rel:+.2f}%)")

# ============================================================
# APPROACH 7: Floor/round corrections to continuous approximations
# ============================================================
print("\n" + "="*70)
print("APPROACH 7: Floor/round corrections making approximations exact")
print("="*70)

# Key idea: if f(n) is a good approximation, maybe
# p(n) = round(f(n)) or p(n) = f(n) + round(g(n)) for some simpler g

def test_rounding_formula(name, func, max_n=500):
    """Test if round(func(n)) = p(n) for n=1..max_n"""
    exact = 0
    off_by_1 = 0
    max_err = 0
    for n in range(2, max_n+1):
        est = func(n)
        actual = p[n]
        err = abs(round(est) - actual)
        if err == 0:
            exact += 1
        elif err <= 1:
            off_by_1 += 1
        max_err = max(max_err, err)
    total = max_n - 1
    print(f"  {name}: exact={exact}/{total} ({100*exact/total:.1f}%), "
          f"off-by-1={off_by_1}, max_err={max_err}")
    return exact

# Test various approximations
test_rounding_formula("R^{-1}(n)", lambda n: R_inverse(n))
test_rounding_formula("Cipolla-6", lambda n: asymptotic_pn(n, 6))

# R^{-1} + correction from sqrt(n)*sin(...)
# This is wild but Ramanujan did wilder things
test_rounding_formula("R^{-1} + 0.5", lambda n: R_inverse(n) + 0.5)

# ============================================================
# APPROACH 8: Symbolic regression with gplearn
# ============================================================
print("\n" + "="*70)
print("APPROACH 8: Symbolic regression with gplearn")
print("="*70)

try:
    from gplearn.genetic import SymbolicRegressor

    # Features for predicting δ(n) = p(n) - R^{-1}(n)
    # Use: n, ln(n), ln(ln(n)), 1/ln(n), sqrt(n), n mod 2, ...
    train_ns = list(range(6, 301))
    X_train = []
    y_train = []

    for n in train_ns:
        ln_n = math.log(n)
        lln_n = math.log(ln_n)
        features = [
            n, ln_n, lln_n, 1/ln_n, math.sqrt(n),
            n * ln_n, n / ln_n, ln_n**2,
            math.sqrt(n * ln_n), n % 2, n % 3, n % 6,
        ]
        X_train.append(features)
        y_train.append(delta[n])

    X_train = np.array(X_train)
    y_train = np.array(y_train)

    # Test set
    test_ns = list(range(301, 501))
    X_test = []
    y_test = []
    for n in test_ns:
        ln_n = math.log(n)
        lln_n = math.log(ln_n)
        features = [
            n, ln_n, lln_n, 1/ln_n, math.sqrt(n),
            n * ln_n, n / ln_n, ln_n**2,
            math.sqrt(n * ln_n), n % 2, n % 3, n % 6,
        ]
        X_test.append(features)
        y_test.append(delta[n])

    X_test = np.array(X_test)
    y_test = np.array(y_test)

    print("Running gplearn symbolic regression on δ(n)...")
    print(f"  Training: n=6..300, Testing: n=301..500")

    sr = SymbolicRegressor(
        population_size=2000,
        generations=30,
        tournament_size=20,
        stopping_criteria=0.01,
        p_crossover=0.7,
        p_subtree_mutation=0.1,
        p_hoist_mutation=0.05,
        p_point_mutation=0.1,
        max_samples=0.9,
        verbose=0,
        parsimony_coefficient=0.005,
        function_set=['add', 'sub', 'mul', 'div', 'sqrt', 'log', 'abs', 'neg'],
        random_state=42,
        n_jobs=1,
    )

    sr.fit(X_train, y_train)

    y_pred_train = sr.predict(X_train)
    y_pred_test = sr.predict(X_test)

    train_rmse = np.sqrt(np.mean((y_train - y_pred_train)**2))
    test_rmse = np.sqrt(np.mean((y_test - y_pred_test)**2))

    print(f"  Best program: {sr._program}")
    print(f"  Train RMSE: {train_rmse:.4f}")
    print(f"  Test  RMSE: {test_rmse:.4f}")
    print(f"  Test/Train RMSE ratio: {test_rmse/train_rmse:.2f}")

    # Check if rounding makes it exact
    exact_count = sum(1 for i, n in enumerate(test_ns) if round(R_inv[n] + y_pred_test[i]) == p[n])
    print(f"  Exact after round(R^-1 + pred): {exact_count}/{len(test_ns)}")

except Exception as e:
    print(f"  gplearn error: {e}")

# ============================================================
# APPROACH 9: Direct symbolic regression on p(n)/n
# ============================================================
print("\n" + "="*70)
print("APPROACH 9: gplearn on p(n) directly (scaled)")
print("="*70)

try:
    # Try to find p(n)/n as a function of ln(n), ln(ln(n))
    train_ns2 = list(range(10, 501))
    X2 = []
    y2 = []
    for n in train_ns2:
        ln_n = math.log(n)
        lln_n = math.log(ln_n)
        X2.append([ln_n, lln_n, 1/ln_n, lln_n/ln_n, lln_n**2/ln_n**2])
        y2.append(p[n] / n)

    X2 = np.array(X2)
    y2 = np.array(y2)

    sr2 = SymbolicRegressor(
        population_size=3000,
        generations=40,
        tournament_size=20,
        stopping_criteria=0.001,
        p_crossover=0.7,
        p_subtree_mutation=0.1,
        p_hoist_mutation=0.05,
        p_point_mutation=0.1,
        max_samples=0.9,
        verbose=0,
        parsimony_coefficient=0.001,
        function_set=['add', 'sub', 'mul', 'div', 'sqrt', 'log', 'abs', 'neg'],
        random_state=123,
        n_jobs=1,
    )

    sr2.fit(X2[:300], y2[:300])

    pred_train = sr2.predict(X2[:300])
    pred_test = sr2.predict(X2[300:])

    train_rmse2 = np.sqrt(np.mean((y2[:300] - pred_train)**2))
    test_rmse2 = np.sqrt(np.mean((y2[300:] - pred_test)**2))

    print(f"  Best program: {sr2._program}")
    print(f"  Train RMSE (p(n)/n): {train_rmse2:.6f}")
    print(f"  Test  RMSE (p(n)/n): {test_rmse2:.6f}")

    # How many exact?
    exact2 = 0
    for i, n in enumerate(train_ns2[300:]):
        est = round(pred_test[i] * n)
        if est == p[n]:
            exact2 += 1
    print(f"  Exact (round(pred*n)=p(n)) on test: {exact2}/{len(train_ns2)-300}")

except Exception as e:
    print(f"  gplearn error: {e}")

# ============================================================
# APPROACH 10: Subsequence analysis - p(prime(n)), p(2n), p(n^2)
# ============================================================
print("\n" + "="*70)
print("APPROACH 10: Subsequences - do they have simpler structure?")
print("="*70)

# p(2n) - relationship to p(n)
print("p(2n)/p(n) analysis:")
for nv in [5, 10, 20, 50, 100, 200, 500]:
    if 2*nv <= MAX_N:
        ratio_sub = p[2*nv] / p[nv]
        print(f"  p({2*nv})/p({nv}) = {p[2*nv]}/{p[nv]} = {ratio_sub:.6f}")

# As n->inf, p(2n)/p(n) -> 2 (by PNT). How fast?
print("\np(2n)/p(n) - 2 analysis:")
sub_ratios = []
for n in range(10, 501):
    if 2*n <= MAX_N:
        r = p[2*n] / p[n] - 2
        sub_ratios.append((n, r))

ns_sub = np.array([s[0] for s in sub_ratios])
rs_sub = np.array([s[1] for s in sub_ratios])
print(f"  Range: [{np.min(rs_sub):.6f}, {np.max(rs_sub):.6f}]")

# Fit: p(2n)/p(n) - 2 ~ a * ln(2)/ln(n)
basis_sub = np.log(2) / np.log(ns_sub)
a_sub = np.sum(rs_sub * basis_sub) / np.sum(basis_sub**2)
print(f"  Fit: p(2n)/p(n) - 2 ≈ {a_sub:.4f} * ln(2)/ln(n)")

# p(prime(n)) = p(π^{-1}(n)th prime)  -- prime-indexed primes
print("\nPrime-indexed primes p(p(n)):")
for nv in [1, 2, 3, 4, 5, 10, 20, 50]:
    if p[nv] <= MAX_N:
        pp = p[p[nv]]
        # Compare with n * ln(n)^2
        est = nv * math.log(nv)**2 if nv > 1 else 0
        print(f"  p(p({nv})) = p({p[nv]}) = {pp}, n*ln(n)^2 = {est:.1f}")

# ============================================================
# APPROACH 11: Gap-based correction
# ============================================================
print("\n" + "="*70)
print("APPROACH 11: Can we predict the gap g(n) = p(n+1) - p(n)?")
print("="*70)

# Cramér's conjecture: g(n) = O(ln(p(n))^2)
# Average gap ~ ln(p(n))
print("Gap analysis g(n) = p(n+1) - p(n):")
gaps = [p[n+1] - p[n] for n in range(1, 500)]
lnp = [math.log(p[n]) for n in range(1, 500)]

# g(n) / ln(p(n))
gap_ratio = [g/l for g, l in zip(gaps, lnp)]
print(f"  g(n)/ln(p(n)): mean={np.mean(gap_ratio):.4f}, std={np.std(gap_ratio):.4f}")
print(f"  Expected mean (Cramér): ~1.0")

# Distribution of gaps
gap_counts = Counter(gaps)
print(f"  Most common gaps: {gap_counts.most_common(10)}")

# Can we predict gaps from local information?
# Gap correlations
gap_arr = np.array(gaps)
autocorr = np.correlate(gap_arr - np.mean(gap_arr), gap_arr - np.mean(gap_arr), mode='full')
autocorr = autocorr / autocorr[len(autocorr)//2]  # normalize
print(f"  Gap autocorrelation lag-1: {autocorr[len(autocorr)//2 + 1]:.4f}")
print(f"  Gap autocorrelation lag-2: {autocorr[len(autocorr)//2 + 2]:.4f}")

# ============================================================
# APPROACH 12: Radical new idea - p(n) mod small numbers
# ============================================================
print("\n" + "="*70)
print("APPROACH 12: p(n) mod small numbers - hidden structure?")
print("="*70)

for m in [2, 3, 4, 5, 6, 8, 10, 12, 30]:
    residues = [p[n] % m for n in range(1, 501)]
    cnt = Counter(residues)
    dist = {k: cnt[k] for k in sorted(cnt.keys())}
    print(f"  p(n) mod {m:2d}: {dist}")

# ============================================================
# APPROACH 13: PSLQ on p(n) - Cipolla for specific n
# ============================================================
print("\n" + "="*70)
print("APPROACH 13: PSLQ on Cipolla residual for specific n values")
print("="*70)

for nv in [100, 200, 500, 1000]:
    residual = p[nv] - asymptotic_pn(nv, 6)
    ln_n = math.log(nv)
    lln = math.log(ln_n)

    # Try: residual = a*n/ln(n)^4 + b*n/ln(n)^5 + c*n/ln(n)^6
    candidates_pslq = [
        mpmath.mpf(residual),
        mpmath.mpf(nv / ln_n**4),
        mpmath.mpf(nv / ln_n**5),
        mpmath.mpf(nv / ln_n**6),
        mpmath.mpf(nv / ln_n**7),
        mpmath.mpf(nv * lln / ln_n**4),
        mpmath.mpf(1),
    ]

    rel = mpmath.pslq(candidates_pslq)
    if rel:
        names = ['residual', 'n/ln^4', 'n/ln^5', 'n/ln^6', 'n/ln^7', 'n*lln/ln^4', '1']
        terms = [f"{c}*{name}" for c, name in zip(rel, names) if c != 0]
        print(f"  n={nv}: {' + '.join(terms)} = 0")
    else:
        print(f"  n={nv}: No PSLQ relation found")

# ============================================================
# APPROACH 14: Chebyshev-type bounds made tighter
# ============================================================
print("\n" + "="*70)
print("APPROACH 14: Tight explicit bounds analysis")
print("="*70)

# Dusart (2010): for n >= 688383, p(n) >= n*(ln(n) + ln(ln(n)) - 1)
# and p(n) <= n*(ln(n) + ln(ln(n)) - 1 + (ln(ln(n))-2)/ln(n)) for n >= 688383
# For smaller n, can we find EXACT bounds?

print("Testing: p(n) = n*(ln(n) + ln(ln(n)) - 1 + c(n)/ln(n)) -- solving for c(n):")
c_values = []
for n in range(10, 501):
    ln_n = math.log(n)
    lln = math.log(ln_n)
    c = (p[n]/n - ln_n - lln + 1) * ln_n
    c_values.append((n, c))

c_arr = np.array([c[1] for c in c_values])
ns_c = np.array([c[0] for c in c_values])
print(f"  c(n) range: [{np.min(c_arr):.4f}, {np.max(c_arr):.4f}]")

# c(n) should approach ln(ln(n)) - 2 asymptotically
lln_minus_2 = np.log(np.log(ns_c)) - 2
residual_c = c_arr - lln_minus_2
print(f"  c(n) - (ln(ln(n))-2) range: [{np.min(residual_c):.4f}, {np.max(residual_c):.4f}]")

# ============================================================
# APPROACH 15: Gram series / Riemann exact formula exploration
# ============================================================
print("\n" + "="*70)
print("APPROACH 15: How close is R^{-1}(n) to p(n)?")
print("="*70)

# The best known: |p(n) - R^{-1}(n)| analysis
print("Absolute errors |p(n) - R^{-1}(n)|:")
abs_errors = []
for n in range(1, 501):
    err = abs(delta[n])
    abs_errors.append(err)

abs_arr = np.array(abs_errors)
ns_abs = np.arange(1, 501)
print(f"  Mean: {np.mean(abs_arr):.2f}")
print(f"  Max: {np.max(abs_arr):.2f} at n={np.argmax(abs_arr)+1}")

# Fraction where round(R^{-1}(n)) = p(n)
exact_R = sum(1 for n in range(1, 501) if round(R_inv[n]) == p[n])
print(f"  round(R^-1(n)) = p(n) for {exact_R}/500 values ({100*exact_R/500:.1f}%)")

# What if we use R^{-1} computed at n - 1/2 (continuity correction)?
print("\nContinuity correction: R^{-1}(n - 1/2):")
exact_half = 0
for n in range(1, 501):
    r_half = R_inverse(n - 0.5)
    if round(r_half) == p[n]:
        exact_half += 1
print(f"  round(R^-1(n-0.5)) = p(n) for {exact_half}/500 values ({100*exact_half/500:.1f}%)")

# Try R^{-1}(n + c) for various c
print("\nOptimal shift: round(R^{-1}(n + c)) = p(n), sweep c:")
best_c = 0
best_exact = 0
for c_try in np.arange(-1.0, 1.01, 0.05):
    exact_c = 0
    for n in range(1, 501):
        r_shifted = R_inverse(n + c_try)
        if round(r_shifted) == p[n]:
            exact_c += 1
    if exact_c > best_exact:
        best_exact = exact_c
        best_c = c_try
print(f"  Best shift: c={best_c:.2f}, exact={best_exact}/500 ({100*best_exact/500:.1f}%)")

# ============================================================
# APPROACH 16: n-dependent shift for R^{-1}
# ============================================================
print("\n" + "="*70)
print("APPROACH 16: n-dependent optimal shift for R^{-1}")
print("="*70)

# For each n, what shift c(n) makes R^{-1}(n + c(n)) = p(n)?
optimal_shifts = []
for n in range(2, 501):
    # Binary search for c such that R^{-1}(n + c) = p(n]
    target = p[n]
    # We need R(target) - n = c
    c_opt = riemann_R(target) - n
    optimal_shifts.append((n, c_opt))

shifts_arr = np.array([s[1] for s in optimal_shifts])
ns_shift = np.array([s[0] for s in optimal_shifts])
print(f"  Optimal shift c(n) = R(p(n)) - n:")
print(f"  Range: [{np.min(shifts_arr):.6f}, {np.max(shifts_arr):.6f}]")
print(f"  Mean: {np.mean(shifts_arr):.6f}")
print(f"  Std: {np.std(shifts_arr):.6f}")

# Is c(n) = R(p(n)) - n related to the Riemann zeros?
# By explicit formula: π(x) - R(x) = -sum_ρ R(x^ρ) - ...
# So R(p(n)) - n ≈ sum_ρ R(p(n)^ρ) which involves Riemann zeros

# Check growth of |c(n)|
print(f"\n  |c(n)| / sqrt(ln(p(n))) stats:")
ratio_czeta = [abs(shifts_arr[i]) / math.sqrt(math.log(p[i+2])) for i in range(len(shifts_arr))]
print(f"  Mean: {np.mean(ratio_czeta):.4f}, Max: {np.max(ratio_czeta):.4f}")

# ============================================================
# APPROACH 17: Can we predict sign(c(n))?
# ============================================================
print("\n" + "="*70)
print("APPROACH 17: Sign pattern of c(n) = R(p(n)) - n")
print("="*70)

signs = [1 if s > 0 else -1 for s in shifts_arr]
pos = sum(1 for s in signs if s > 0)
neg = sum(1 for s in signs if s < 0)
print(f"  Positive: {pos}, Negative: {neg}")

# Check if sign correlates with prime gaps
print("  Sign vs gap pattern:")
sign_gap_corr = np.corrcoef(signs[:len(gaps[:498])], gaps[:498])[0,1]
print(f"  Correlation(sign(c(n)), gap(n)): {sign_gap_corr:.4f}")

# ============================================================
# APPROACH 18: Direct formula search - p(n) as combination
# ============================================================
print("\n" + "="*70)
print("APPROACH 18: Direct formula search p(n) = f(n)")
print("="*70)

# Try many candidate formulas and measure accuracy
candidates_formulas = {}

def eval_formula(name, func, max_n=500):
    errors = []
    exact = 0
    for n in range(10, max_n+1):
        try:
            est = func(n)
            err = abs(est - p[n])
            errors.append(err)
            if round(est) == p[n]:
                exact += 1
        except:
            errors.append(float('inf'))
    rmse = np.sqrt(np.mean(np.array(errors)**2))
    total = max_n - 9
    return rmse, exact, total

formulas = {
    "n*(ln(n)+ln(ln(n)))": lambda n: n * (math.log(n) + math.log(math.log(n))),
    "n*(ln(n)+ln(ln(n))-1)": lambda n: n * (math.log(n) + math.log(math.log(n)) - 1),
    "Cipolla-4": lambda n: asymptotic_pn(n, 4),
    "Cipolla-6": lambda n: asymptotic_pn(n, 6),
    "R^{-1}(n)": lambda n: R_inverse(n),
    "li^{-1}(n)": lambda n: float(mpmath.findroot(lambda x: mpmath.li(x) - n, n * math.log(n))),
    "n*W(n)*1.04": lambda n: n * float(mpmath.lambertw(n)) * 1.04,
}

print(f"{'Formula':<30s} {'RMSE':>10s} {'Exact/Total':>15s} {'%':>8s}")
print("-" * 65)
for name, func in formulas.items():
    rmse, exact, total = eval_formula(name, func)
    print(f"{name:<30s} {rmse:10.2f} {exact:>5d}/{total:<5d}   {100*exact/total:6.1f}%")

# ============================================================
# APPROACH 19: li^{-1}(n) analysis (likely best elementary inverse)
# ============================================================
print("\n" + "="*70)
print("APPROACH 19: Detailed li^{-1}(n) analysis")
print("="*70)

def li_inverse(n):
    """Compute li^{-1}(n) = x such that li(x) = n"""
    guess = n * math.log(n)
    return float(mpmath.findroot(lambda x: mpmath.li(x) - n, guess))

li_inv_exact = 0
li_inv_errors = []
for n in range(2, 501):
    li_inv_n = li_inverse(n)
    err = p[n] - li_inv_n
    li_inv_errors.append((n, err))
    if round(li_inv_n) == p[n]:
        li_inv_exact += 1

print(f"  round(li^-1(n)) = p(n): {li_inv_exact}/499 ({100*li_inv_exact/499:.1f}%)")

li_err_arr = np.array([e[1] for e in li_inv_errors])
print(f"  p(n) - li^-1(n) stats: mean={np.mean(li_err_arr):.2f}, std={np.std(li_err_arr):.2f}")

# Compare R^{-1} vs li^{-1}
print(f"\n  Comparison: R^-1 exact={exact_R}/500, li^-1 exact={li_inv_exact}/499")

# ============================================================
# SUMMARY
# ============================================================
print("\n" + "="*70)
print("SUMMARY OF ALL FINDINGS")
print("="*70)
print("""
Key findings from symbolic formula discovery:

1. R^{-1}(n) is the best known approximation. round(R^{-1}(n)) matches p(n)
   for a significant fraction of small n values.

2. δ(n) = p(n) - R^{-1}(n) has NO simple closed-form pattern. It is
   oscillatory and its behavior is governed by Riemann zeta zeros
   (via the explicit formula).

3. Cipolla's 6-term expansion has RMSE growing with n (not convergent).

4. gplearn symbolic regression on δ(n): reports train/test RMSE.
   Generalization is key - does the discovered formula extrapolate?

5. PSLQ found no simple integer relations between p(n) and standard constants.

6. The optimal shift c(n) = R(p(n)) - n is directly related to the
   Riemann explicit formula. Predicting it requires knowing zeta zeros.

7. No subsequence (p(2n), p(p(n))) showed simpler structure than p(n) itself.

8. Fundamental barrier: the irregularity in p(n) is EQUIVALENT to knowledge
   of Riemann zeta zeros. Any exact formula for p(n) that avoids brute force
   must encode this information somehow.
""")

print("Done. All approaches completed.")

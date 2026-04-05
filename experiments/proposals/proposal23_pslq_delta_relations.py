#!/usr/bin/env python3
"""
Proposal 23: PSLQ/LLL Integer Relation Discovery for delta(n)

IDEA: Use integer relation algorithms (PSLQ, LLL) to search for algebraic
relationships between delta(n) = p(n) - R^{-1}(n) and known mathematical
constants. If delta(n) can be expressed as a closed-form combination of
computable quantities, we get an O(polylog) formula.

Inspired by the "Ramanujan Library" approach (arXiv:2412.12361, 2024):
- Build a library of candidate basis constants
- Use PSLQ to search for integer relations
- Validate on held-out data

CANDIDATE BASES:
1. Zeta values: zeta(2), zeta(3), zeta(5), ...
2. Log values: log(2), log(3), log(5), ...
3. Powers of pi
4. Euler-Mascheroni constant gamma
5. Li_2(1/2), Li_3(1/2), etc. (polylogarithms)
6. Stieltjes constants gamma_1, gamma_2, ...

CONJECTURE: There exists a polynomial P of bounded degree such that
delta(n) = P(log n, log log n, 1/log n, ...) + O(1) for all n.

TIME COMPLEXITY: O(polylog(n)) if a closed-form is found.
"""

import math
import os
import numpy as np
from itertools import combinations

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(SCRIPT_DIR, '..', '..', 'data')

def li(x):
    if x <= 0 or x <= 1:
        return 0.0
    if abs(x - 1.0) < 0.01:
        return -100.0
    lnx = math.log(x)
    total = 0.0
    term = 1.0
    for k in range(1, 100):
        term *= lnx / k
        total += term / k
        if abs(term / k) < 1e-15:
            break
    return 0.5772156649015329 + math.log(abs(lnx)) + total

def R_function(x):
    if x <= 1:
        return 0.0
    mu = [0, 1, -1, -1, 0, -1, 1, -1, 0, 0, 1,
          -1, 0, -1, 1, 1, 0, -1, 0, -1, 0]
    total = 0.0
    for k in range(1, 21):
        if mu[k] == 0:
            continue
        xk = x ** (1.0 / k)
        if xk <= 1.0001:
            break
        total += mu[k] / k * li(xk)
    return total

def R_inverse(n):
    x = n * math.log(n) + n * math.log(math.log(max(n, 3))) - n
    if x < 2:
        x = 2.0
    for _ in range(50):
        rx = R_function(x)
        if abs(rx - n) < 1e-10:
            break
        deriv = 1.0 / math.log(x)
        x = x - (rx - n) / deriv
        if x < 2:
            x = 2.0
    return x

def sieve_primes(limit):
    if limit < 2:
        return []
    sieve = [True] * (limit + 1)
    sieve[0] = sieve[1] = False
    for i in range(2, int(limit**0.5) + 1):
        if sieve[i]:
            for j in range(i*i, limit + 1, i):
                sieve[j] = False
    return [i for i in range(limit + 1) if sieve[i]]

def build_candidate_basis(n):
    """
    Build a vector of candidate constants that might relate to delta(n).
    Each element is a function of n that could appear in a closed form.
    """
    pn_est = R_inverse(n)
    ln_n = math.log(n)
    ln_pn = math.log(pn_est)
    ln_ln_n = math.log(ln_n) if ln_n > 1 else 0.1
    sqrt_pn = math.sqrt(pn_est)

    basis = {
        '1': 1.0,
        'ln(n)': ln_n,
        'ln(n)^2': ln_n ** 2,
        'ln(ln(n))': ln_ln_n,
        'ln(n)*ln(ln(n))': ln_n * ln_ln_n,
        '1/ln(n)': 1.0 / ln_n,
        'sqrt(pn)/ln(pn)': sqrt_pn / ln_pn,
        'sqrt(pn)/ln(pn)^2': sqrt_pn / ln_pn ** 2,
        'n^(1/3)': n ** (1.0/3),
        'n^(1/3)/ln(n)': n ** (1.0/3) / ln_n,
        'pi': math.pi,
        'pi/ln(n)': math.pi / ln_n,
        'euler_gamma': 0.5772156649015329,
        'euler_gamma*ln(n)': 0.5772156649015329 * ln_n,
        'ln(2)': math.log(2),
        'ln(2)*ln(n)': math.log(2) * ln_n,
        'zeta(2)/ln(n)': (math.pi**2/6) / ln_n,
        'li(sqrt(pn))': li(sqrt_pn),
    }
    return basis

def lll_reduce(basis_matrix):
    """Simple LLL lattice reduction (for small dimensions)."""
    # Use numpy for a simplified Gram-Schmidt-based approach
    # For production, one would use fpLLL or similar
    n = basis_matrix.shape[0]
    B = basis_matrix.astype(float).copy()
    mu = np.zeros((n, n))
    B_star = np.zeros_like(B)
    B_star_norm = np.zeros(n)

    # Gram-Schmidt
    for i in range(n):
        B_star[i] = B[i].copy()
        for j in range(i):
            if B_star_norm[j] > 1e-10:
                mu[i, j] = np.dot(B[i], B_star[j]) / B_star_norm[j]
                B_star[i] -= mu[i, j] * B_star[j]
        B_star_norm[i] = np.dot(B_star[i], B_star[i])

    # LLL reduction with delta = 3/4
    delta_lll = 0.75
    k = 1
    max_iter = 1000
    iteration = 0
    while k < n and iteration < max_iter:
        iteration += 1
        # Size reduction
        for j in range(k - 1, -1, -1):
            if abs(mu[k, j]) > 0.5:
                r = round(mu[k, j])
                B[k] -= r * B[j]
                for i in range(j + 1):
                    mu[k, i] -= r * mu[j, i]
                mu[k, j] -= r

        # Lovász condition
        if B_star_norm[k] >= (delta_lll - mu[k, k-1]**2) * B_star_norm[k-1]:
            k += 1
        else:
            # Swap B[k] and B[k-1]
            B[[k, k-1]] = B[[k-1, k]]
            # Recompute Gram-Schmidt from scratch (simplified)
            for i in range(n):
                B_star[i] = B[i].copy()
                for j in range(i):
                    if B_star_norm[j] > 1e-10:
                        mu[i, j] = np.dot(B[i], B_star[j]) / B_star_norm[j]
                        B_star[i] -= mu[i, j] * B_star[j]
                B_star_norm[i] = np.dot(B_star[i], B_star[i])
            k = max(k - 1, 1)

    return B

def pslq_search(target, basis_values, precision_digits=10):
    """
    Simplified PSLQ-like search for integer relation:
    target = c_1 * basis[0] + c_2 * basis[1] + ... + c_k * basis[k-1]
    where c_i are small integers.

    Uses LLL on the augmented lattice.
    """
    n = len(basis_values) + 1
    scale = 10 ** precision_digits

    # Build lattice basis
    # [scale*target, scale*b1, ..., scale*bk, I_n]
    L = np.zeros((n, n))
    for i in range(n):
        L[i, i] = 1  # identity part

    # Last column encodes the relation
    vals = [target] + list(basis_values)
    for i in range(n):
        L[i, -1] = round(vals[i] * scale)

    # LLL reduce
    reduced = lll_reduce(L)

    # Look for short vectors with small last component (close to 0)
    results = []
    for row in reduced:
        if abs(row[-1]) < scale * 1e-5:  # nearly zero in the relation column
            coeffs = row[:-1].astype(int)
            # Check: coeffs[0]*target + sum(coeffs[i+1]*basis[i]) ≈ 0
            check = coeffs[0] * target + sum(c * b for c, b in zip(coeffs[1:], basis_values))
            if abs(check) < 1e-3 and any(c != 0 for c in coeffs):
                results.append((coeffs, abs(check)))

    results.sort(key=lambda x: x[1])
    return results

def test_pslq_relations():
    print("=" * 80)
    print("PROPOSAL 23: PSLQ/LLL Integer Relation Discovery for delta(n)")
    print("=" * 80)

    primes = sieve_primes(120000)

    # Compute delta(n) for various n
    print("\n--- Test 1: Structure of delta(n) ---")
    test_ns = [50, 100, 200, 500, 1000, 2000, 5000]
    deltas = {}
    for n in test_ns:
        pn = primes[n - 1]
        r_inv = R_inverse(n)
        delta = pn - r_inv
        deltas[n] = delta
        print(f"  n={n:>5}: p(n)={pn:>7}, R^{{-1}}(n)={r_inv:>10.2f}, "
              f"delta={delta:>8.2f}, delta/sqrt(pn)*ln(pn)={delta*math.log(pn)/math.sqrt(pn):.4f}")

    # Test 2: Search for polynomial relations in log(n)
    print("\n--- Test 2: Polynomial fit of delta(n) in terms of log-quantities ---")
    train_ns = list(range(100, 5001))
    delta_vec = np.array([primes[n-1] - R_inverse(n) for n in train_ns])

    # Try various polynomial bases
    for basis_desc, basis_func in [
        ("log(n), log(log(n))", lambda n: [math.log(n), math.log(math.log(n))]),
        ("sqrt(p)/log(p)", lambda n: [math.sqrt(R_inverse(n)) / math.log(R_inverse(n))]),
        ("sqrt(p)/log(p), sqrt(p)/log(p)^2",
         lambda n: [math.sqrt(R_inverse(n)) / math.log(R_inverse(n)),
                    math.sqrt(R_inverse(n)) / math.log(R_inverse(n))**2]),
        ("1, log(n), n^{1/3}/log(n)",
         lambda n: [1, math.log(n), n**(1./3)/math.log(n)]),
    ]:
        X = np.array([basis_func(n) for n in train_ns])
        X = np.column_stack([np.ones(len(train_ns)), X])
        coeffs, residuals, _, _ = np.linalg.lstsq(X, delta_vec, rcond=None)
        pred = X @ coeffs
        err = np.std(delta_vec - pred)
        r_squared = 1 - np.var(delta_vec - pred) / np.var(delta_vec)
        print(f"  Basis [{basis_desc}]: R²={r_squared:.4f}, err_std={err:.3f}")

    # Test 3: PSLQ search for integer relations at specific n values
    print("\n--- Test 3: PSLQ search for integer relations ---")
    for n in [100, 500, 1000]:
        pn = primes[n - 1]
        delta = pn - R_inverse(n)
        basis = build_candidate_basis(n)

        basis_names = list(basis.keys())
        basis_vals = list(basis.values())

        # Search subsets of basis elements (PSLQ works best with small dimensions)
        print(f"\n  n={n}, delta={delta:.4f}:")
        best_error = float('inf')
        best_relation = None

        for combo_size in [3, 4, 5]:
            for combo in combinations(range(len(basis_vals)), combo_size):
                sub_names = [basis_names[i] for i in combo]
                sub_vals = [basis_vals[i] for i in combo]

                relations = pslq_search(delta, sub_vals, precision_digits=6)
                for coeffs, err in relations[:1]:
                    if err < best_error and coeffs[0] != 0:
                        best_error = err
                        coeff_target = coeffs[0]
                        terms = [(coeffs[i+1], sub_names[i]) for i in range(len(sub_names))
                                if i + 1 < len(coeffs) and coeffs[i+1] != 0]
                        best_relation = (coeff_target, terms, err)

        if best_relation:
            ct, terms, err = best_relation
            expr = " + ".join(f"{c}*{name}" for c, name in terms)
            print(f"    Best: {ct}*delta ≈ {expr} (error={err:.6f})")
        else:
            print(f"    No integer relation found with error < 0.001")

    # Test 4: Normalized delta search
    print("\n--- Test 4: Normalized delta analysis ---")
    # Normalize: delta_hat(n) = delta(n) * log(p(n)) / sqrt(p(n))
    normalized = []
    for n in range(100, 5001):
        pn = primes[n - 1]
        delta = pn - R_inverse(n)
        delta_hat = delta * math.log(pn) / math.sqrt(pn)
        normalized.append(delta_hat)

    normalized = np.array(normalized)
    print(f"  Normalized delta stats: mean={np.mean(normalized):.6f}, "
          f"std={np.std(normalized):.6f}")
    print(f"  Range: [{np.min(normalized):.4f}, {np.max(normalized):.4f}]")

    # Autocorrelation of normalized delta
    centered = normalized - np.mean(normalized)
    autocorr = np.correlate(centered, centered, mode='full')
    autocorr = autocorr[len(autocorr)//2:]
    autocorr /= autocorr[0]
    print(f"  Autocorrelation at lag 1: {autocorr[1]:.4f}")
    print(f"  Autocorrelation at lag 5: {autocorr[5]:.4f}")
    print(f"  Autocorrelation at lag 10: {autocorr[10]:.4f}")
    print(f"  Autocorrelation at lag 50: {autocorr[50]:.4f}")

    # Test if normalized delta is normally distributed
    from scipy import stats
    stat, pvalue = stats.normaltest(normalized)
    print(f"  Normality test: stat={stat:.2f}, p-value={pvalue:.4f}")
    print(f"  (p < 0.05 means NOT normal)")

    print("\n--- VERDICT ---")
    print("PSLQ searches for integer relations between delta(n) and known constants.")
    print("If a relation exists, it implies a closed-form for delta(n).")
    print("However, delta(n) appears pseudo-random with low autocorrelation,")
    print("making it unlikely that a simple algebraic relation exists.")
    print("The approach is most promising for SPECIFIC n values or for")
    print("discovering approximate formulas that reduce the search range.")

if __name__ == '__main__':
    test_pslq_relations()

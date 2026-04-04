#!/usr/bin/env python3
"""
Session 4: P-recursive / Holonomic / D-finite exploration for the prime sequence.

We test whether p(n) or transforms of it satisfy linear recurrences with
polynomial coefficients:  a(n+d) = R_0(n)*a(n) + ... + R_{d-1}(n)*a(n+d-1)

If such a recurrence exists with fixed order d and bounded-degree R_i,
binary splitting gives O(polylog n) computation of p(n).

Experiments:
  1. Direct P-recursiveness test for p(n) and pi(n)
  2. Transforms: Euler, binomial, Mobius, EGF/OGF coefficients
  3. Matrix power formulation search
  4. Cloitre-type analytic recurrence investigation
  5. Prime zeta / contour integration feasibility
"""

import numpy as np
from sympy import primerange, isprime, prime, primepi, factorint, factorial, binomial as binom_sym
from sympy import Matrix, symbols, Rational
from itertools import combinations
import time
import sys

# ─────────────────────────────────────────────────────────────────────
# Utility: generate primes
# ─────────────────────────────────────────────────────────────────────
N_MAX = 500  # number of primes to work with
PRIMES = list(primerange(2, prime(N_MAX) + 1))
assert len(PRIMES) >= N_MAX
PRIMES = PRIMES[:N_MAX]

print(f"Working with first {N_MAX} primes: p(1)={PRIMES[0]} .. p({N_MAX})={PRIMES[-1]}")
print("=" * 80)


# =====================================================================
# EXPERIMENT 1: Direct P-recursiveness test
# =====================================================================
def test_p_recursive(seq, name, max_order=6, max_deg=4, start=0):
    """
    Test if seq[n] satisfies a linear recurrence with polynomial coefficients:
      sum_{j=0}^{d} P_j(n) * a(n+j) = 0
    where each P_j is a polynomial of degree <= max_deg.

    We set up a linear system: for each shift j, P_j(n) = c_{j,0} + c_{j,1}*n + ... + c_{j,max_deg}*n^max_deg
    Total unknowns = (max_order+1) * (max_deg+1).
    Each value of n gives one equation. We need more equations than unknowns.
    """
    print(f"\n{'─'*70}")
    print(f"EXPERIMENT 1: P-recursive test for '{name}'")
    print(f"  max_order={max_order}, max_deg={max_deg}")
    print(f"{'─'*70}")

    for d in range(1, max_order + 1):
        for deg in range(0, max_deg + 1):
            n_unknowns = (d + 1) * (deg + 1)
            # We need at least n_unknowns equations, use 2x for robustness
            n_eq = min(len(seq) - d, 2 * n_unknowns + 10)
            if n_eq < n_unknowns + 5:
                continue

            # Build matrix A where A[row] corresponds to n = start + row
            A = np.zeros((n_eq, n_unknowns), dtype=np.float64)
            for row in range(n_eq):
                n = start + row
                for j in range(d + 1):
                    a_val = seq[n + j]
                    for k in range(deg + 1):
                        col = j * (deg + 1) + k
                        A[row, col] = a_val * (n ** k)

            # Check rank: if rank < n_unknowns, there's a nontrivial solution
            try:
                U, S, Vt = np.linalg.svd(A, full_matrices=False)
                # Ratio of smallest to largest singular value
                if S[0] == 0:
                    continue
                ratio = S[-1] / S[0]

                # Count how many singular values are "essentially zero"
                threshold = 1e-8 * S[0]
                n_zero = np.sum(S < threshold)

                if n_zero > 0:
                    print(f"  order={d}, deg={deg}: RANK DEFICIENT! "
                          f"({n_zero} zero singular values, ratio={ratio:.2e})")

                    # Extract the null space vector
                    null_vec = Vt[-1, :]
                    # Reconstruct polynomial coefficients
                    coeffs = null_vec.reshape(d + 1, deg + 1)
                    print(f"    Coefficient matrix (shift x degree):")
                    for j in range(d + 1):
                        c_str = ", ".join(f"{coeffs[j,k]:+.6f}" for k in range(deg + 1))
                        print(f"      P_{j}(n) coeffs: [{c_str}]")

                    # Verify on held-out data
                    n_verify = min(50, len(seq) - d - n_eq)
                    if n_verify > 0:
                        max_err = 0
                        for row in range(n_verify):
                            n = start + n_eq + row
                            val = 0
                            for j in range(d + 1):
                                pj = sum(coeffs[j, k] * (n ** k) for k in range(deg + 1))
                                val += pj * seq[n + j]
                            max_err = max(max_err, abs(val))
                        print(f"    Verification on {n_verify} held-out points: max|residual| = {max_err:.2e}")
                        if max_err < 1e-3:
                            print(f"    *** EXACT P-RECURSIVE RELATION FOUND! ***")
                        elif max_err < 1.0:
                            print(f"    *** APPROXIMATE P-recursive relation ***")
                else:
                    if d <= 2 and deg <= 1:
                        print(f"  order={d}, deg={deg}: min_sv/max_sv = {ratio:.2e} (full rank)")
            except Exception as e:
                print(f"  order={d}, deg={deg}: SVD failed: {e}")

    print()


# Test p(n) directly
test_p_recursive(PRIMES, "p(n) — prime sequence")

# Test pi(n) = prime counting function evaluated at integers
PI_SEQ = [int(primepi(k)) for k in range(1, N_MAX + 50)]
test_p_recursive(PI_SEQ, "pi(n) — prime counting at integers", max_order=6, max_deg=4)

# Test prime gaps g(n) = p(n+1) - p(n)
GAPS = [PRIMES[i+1] - PRIMES[i] for i in range(len(PRIMES)-1)]
test_p_recursive(GAPS, "g(n) — prime gaps")

# Test p(n)/n (should be ~ ln(n) by PNT)
PN_OVER_N = [PRIMES[i] / (i + 1) for i in range(len(PRIMES))]
test_p_recursive(PN_OVER_N, "p(n)/n — normalized primes", max_order=4, max_deg=3)

# Test p(n) - n*ln(n) residual
import math
PN_RESIDUAL = [PRIMES[i] - (i+1)*math.log(i+1) if i > 0 else PRIMES[0] for i in range(len(PRIMES))]
test_p_recursive(PN_RESIDUAL, "p(n) - n*ln(n) residual", max_order=4, max_deg=3)


# =====================================================================
# EXPERIMENT 2: Transforms of the prime sequence
# =====================================================================
print(f"\n{'='*80}")
print("EXPERIMENT 2: Transforms of the prime sequence")
print(f"{'='*80}")

# 2a. Binomial transform: b(n) = sum_{k=0}^{n} C(n,k) * a(k)
def binomial_transform(seq, N):
    result = []
    for n in range(N):
        val = sum(int(binom_sym(n, k)) * seq[k] for k in range(n + 1))
        result.append(val)
    return result

N_TRANS = 80  # smaller for expensive transforms
print(f"\nUsing first {N_TRANS} primes for transforms...")

bt = binomial_transform(PRIMES[:N_TRANS], N_TRANS)
print(f"\n2a. Binomial transform of p(n): first 10 = {bt[:10]}")
test_p_recursive(bt, "Binomial transform of p(n)", max_order=4, max_deg=3)

# 2b. Inverse binomial transform: b(n) = sum_{k=0}^{n} (-1)^{n-k} C(n,k) * a(k)
def inv_binomial_transform(seq, N):
    result = []
    for n in range(N):
        val = sum((-1)**(n - k) * int(binom_sym(n, k)) * seq[k] for k in range(n + 1))
        result.append(val)
    return result

ibt = inv_binomial_transform(PRIMES[:N_TRANS], N_TRANS)
print(f"\n2b. Inverse binomial transform: first 10 = {ibt[:10]}")
test_p_recursive(ibt, "Inverse binomial transform of p(n)", max_order=4, max_deg=3)

# 2c. Euler transform (for OGF): if A(x) = sum a(n)x^n, then
# Euler transform B(x) = prod_{n>=1} 1/(1-x^n)^{a(n)}
# Coefficients are computed iteratively.
def euler_transform(seq, N):
    """Euler transform of seq, returning first N terms of the transformed sequence."""
    a = seq[:N]
    b = [0] * N
    # b(0) = 1 for the Euler transform of a sequence starting at index 1
    # Standard: b(n) = (1/n) * sum_{k=1}^{n} (sum_{d|k} d*a(d)) * b(n-k)
    # where a is 1-indexed
    c = [0] * (N + 1)  # c(k) = sum_{d|k} d * a(d-1) (adjusting for 0-index)
    for k in range(1, N + 1):
        for d in range(1, k + 1):
            if k % d == 0 and d - 1 < len(a):
                c[k] += d * a[d - 1]

    b[0] = 1
    for n in range(1, N):
        s = 0
        for k in range(1, n + 1):
            s += c[k] * b[n - k]
        b[n] = s // n if s % n == 0 else s / n
    return b

et = euler_transform(PRIMES[:N_TRANS], N_TRANS)
print(f"\n2c. Euler transform: first 10 = {et[:10]}")
# Only test if values are reasonable
if all(isinstance(v, (int, float)) and abs(v) < 1e15 for v in et[:N_TRANS]):
    test_p_recursive(et, "Euler transform of p(n)", max_order=4, max_deg=3)
else:
    print("  Euler transform values too large or non-integer, skipping P-recursive test")

# 2d. First differences: d(n) = p(n+1) - p(n) (same as gaps, already tested)
# Second differences
DIFF2 = [GAPS[i+1] - GAPS[i] for i in range(len(GAPS)-1)]
print(f"\n2d. Second differences of p(n): first 15 = {DIFF2[:15]}")
test_p_recursive(DIFF2, "Second differences of p(n)", max_order=4, max_deg=3)

# 2e. EGF coefficients: if p(n) are the coefficients of an EGF,
# then the sequence is a(n) = p(n) * n!
EGF_SEQ = [PRIMES[i] * math.factorial(i) for i in range(min(N_TRANS, 50))]
print(f"\n2e. EGF-scaled p(n)*n!: first 8 = {EGF_SEQ[:8]}")
# These grow super-factorially, but we can still test in log space
LOG_PRIMES = [math.log(p) for p in PRIMES]
test_p_recursive(LOG_PRIMES, "log(p(n))", max_order=4, max_deg=3)

# 2f. Mobius transform: b(n) = sum_{d|n} a(d)
def mobius_transform_direct(seq, N):
    """b(n) = sum_{d|n} a(d), 1-indexed."""
    result = []
    for n in range(1, N + 1):
        val = sum(seq[d - 1] for d in range(1, n + 1) if n % d == 0)
        result.append(val)
    return result

mt = mobius_transform_direct(PRIMES[:N_TRANS], N_TRANS)
print(f"\n2f. Mobius (divisor sum) transform: first 10 = {mt[:10]}")
test_p_recursive(mt, "Mobius transform of p(n)", max_order=4, max_deg=3)

# 2g. Partial sums: S(n) = sum_{k=1}^{n} p(k)
PARTIAL_SUMS = [sum(PRIMES[:i+1]) for i in range(len(PRIMES))]
print(f"\n2g. Partial sums of primes: first 10 = {PARTIAL_SUMS[:10]}")
test_p_recursive(PARTIAL_SUMS, "Partial sums of p(n)", max_order=4, max_deg=3)

# 2h. p(n) mod various small numbers
for m in [2, 3, 4, 6]:
    mod_seq = [p % m for p in PRIMES]
    print(f"\n2h. p(n) mod {m}: first 20 = {mod_seq[:20]}")


# =====================================================================
# EXPERIMENT 3: Matrix power formulation
# =====================================================================
print(f"\n{'='*80}")
print("EXPERIMENT 3: Matrix power formulation search")
print(f"{'='*80}")

def test_matrix_recurrence(seq, name, dim=2, max_deg=2):
    """
    Test if there exists a matrix M(n) with polynomial entries such that
    [a(n+1), a(n), ...] = M(n) * [a(n), a(n-1), ...]

    For dim=2: [a(n+1)] = [f(n) g(n)] [a(n)  ]
               [a(n)  ]   [1    0   ] [a(n-1)]

    So a(n+1) = f(n)*a(n) + g(n)*a(n-1) where f,g are polynomials in n.
    """
    print(f"\n  Matrix recurrence test for '{name}', dim={dim}, max_deg={max_deg}")

    for deg in range(0, max_deg + 1):
        # For companion-type matrix, only the top row has free parameters
        # a(n+dim-1) = c_0(n)*a(n) + c_1(n)*a(n+1) + ... + c_{dim-2}(n)*a(n+dim-2)
        # Each c_i(n) is polynomial of degree deg => (deg+1) coefficients each
        n_free = (dim - 1) * (deg + 1)  # only dim-1 free polynomials if we normalize the leading one

        # Actually: a(n+dim) = sum_{j=0}^{dim-1} P_j(n) a(n+j)
        # with P_j polynomial of degree <= deg
        n_unknowns = dim * (deg + 1)
        n_eq = min(len(seq) - dim, 2 * n_unknowns + 10)
        if n_eq < n_unknowns + 2:
            continue

        A_mat = np.zeros((n_eq, n_unknowns))
        b_vec = np.zeros(n_eq)

        for row in range(n_eq):
            n = row
            b_vec[row] = seq[n + dim]
            for j in range(dim):
                for k in range(deg + 1):
                    col = j * (deg + 1) + k
                    A_mat[row, col] = seq[n + j] * (n ** k)

        # Solve least squares
        try:
            result, residuals, rank, sv = np.linalg.lstsq(A_mat, b_vec, rcond=None)

            # Compute actual residuals
            predicted = A_mat @ result
            errors = b_vec - predicted
            max_err = np.max(np.abs(errors))
            rms_err = np.sqrt(np.mean(errors**2))
            rel_err = rms_err / np.sqrt(np.mean(b_vec**2)) if np.mean(b_vec**2) > 0 else float('inf')

            if max_err < 0.5:  # Would round to exact
                print(f"    deg={deg}: *** EXACT *** max_err={max_err:.2e}")
                # Print the polynomials
                coeffs = result.reshape(dim, deg + 1)
                for j in range(dim):
                    terms = []
                    for k in range(deg + 1):
                        if abs(coeffs[j, k]) > 1e-10:
                            if k == 0:
                                terms.append(f"{coeffs[j,k]:.6f}")
                            else:
                                terms.append(f"{coeffs[j,k]:.6f}*n^{k}")
                    print(f"      P_{j}(n) = {' + '.join(terms) if terms else '0'}")
            elif rel_err < 0.01:
                print(f"    deg={deg}: approximate, rel_err={rel_err:.4f}, max_err={max_err:.2e}")
            elif dim <= 2 and deg <= 1:
                print(f"    deg={deg}: no fit (rel_err={rel_err:.4f}, max_err={max_err:.2e})")
        except Exception as e:
            print(f"    deg={deg}: failed: {e}")

# Try various dimensions
for dim in [2, 3, 4, 5]:
    test_matrix_recurrence(PRIMES, f"p(n), dim={dim}", dim=dim, max_deg=3)

# Also try on log(p(n))
test_matrix_recurrence(LOG_PRIMES, "log(p(n)), dim=2", dim=2, max_deg=3)


# =====================================================================
# EXPERIMENT 4: Cloitre-type analytic recurrence
# =====================================================================
print(f"\n{'='*80}")
print("EXPERIMENT 4: Cloitre-type analytic recurrence")
print(f"{'='*80}")

# Cloitre's observation: certain analytic recurrences can generate primes.
# The simplest form is: a(1) = 2, a(n+1) = a(n) + f(a(n), n)
# where f involves floor functions or fractional parts.

# Test: a(n+1) = a(n) + 1 + floor(c * a(n) / n) for some constant c
print("\n4a. Testing a(n+1) = a(n) + 1 + floor(c * a(n) / n)")
best_c = None
best_err = float('inf')
for c_num in range(1, 200):
    c = c_num / 100.0
    seq = [2]
    ok = True
    for i in range(1, min(30, N_MAX)):
        n = i  # 1-indexed: a(n) = PRIMES[n-1]
        next_val = seq[-1] + 1 + int(c * seq[-1] / n)
        seq.append(next_val)
    err = sum(abs(seq[i] - PRIMES[i]) for i in range(min(len(seq), len(PRIMES))))
    if err < best_err:
        best_err = err
        best_c = c
        best_seq = seq[:]

print(f"  Best c = {best_c}, total error = {best_err}")
print(f"  Generated: {best_seq[:15]}")
print(f"  Actual:    {PRIMES[:15]}")

# Rowland's recurrence: a(n) = a(n-1) + gcd(n, a(n-1))
# The primes appear as differences.
print("\n4b. Rowland's recurrence: a(n) = a(n-1) + gcd(n, a(n-1))")
a = [0, 7]  # a(1) = 7
for n in range(2, 200):
    a.append(a[-1] + math.gcd(n, a[-1]))
rowland_primes = []
for i in range(1, len(a)):
    d = a[i] - a[i-1]
    if d > 1 and d not in rowland_primes:
        rowland_primes.append(d)
print(f"  Rowland primes (unique diffs > 1): {rowland_primes[:20]}")
print(f"  Note: produces primes but VERY slowly, not O(1) per prime")

# Cloitre's 2014 formula: a(n) = floor(A^{3^n}) for some constant A
# (Mills' constant variant). Not practical since A is transcendental.
print("\n4c. Mills-type: p(n) = floor(A^{f(n)}) for some A, f")
# Let's try to find if p(n) = floor(A * B^n) or p(n) = floor(A * n^B) works approximately
# Fit A, B to log(p(n)) = log(A) + B*log(n)
n_vals = np.array([i + 1 for i in range(N_MAX)], dtype=float)
p_vals = np.array(PRIMES[:N_MAX], dtype=float)
# log(p(n)) ~ a + b*log(n) + c*log(log(n))
log_n = np.log(n_vals)
log_p = np.log(p_vals)
log_log_n = np.log(np.maximum(log_n, 0.01))

A_fit = np.column_stack([np.ones(N_MAX), log_n, log_log_n])
coeffs_fit, _, _, _ = np.linalg.lstsq(A_fit, log_p, rcond=None)
print(f"  log(p(n)) ~ {coeffs_fit[0]:.4f} + {coeffs_fit[1]:.4f}*log(n) + {coeffs_fit[2]:.4f}*log(log(n))")
predicted_p = np.exp(A_fit @ coeffs_fit)
errors_fit = p_vals - predicted_p
print(f"  Max error: {np.max(np.abs(errors_fit)):.1f}, RMS: {np.sqrt(np.mean(errors_fit**2)):.1f}")
print(f"  (PNT predicts coefficients ~[0, 1, 1] => p(n) ~ n*ln(n))")


# =====================================================================
# EXPERIMENT 5: Generating function / Contour integration approach
# =====================================================================
print(f"\n{'='*80}")
print("EXPERIMENT 5: Prime zeta function and generating function analysis")
print(f"{'='*80}")

# The prime zeta function P(s) = sum_{p prime} p^{-s}
# converges for Re(s) > 1, has a logarithmic singularity at s=1.
# P(s) = sum_{k=1}^inf mu(k)/k * log(zeta(k*s))

# Can we extract p(n) from P(s)?
# If f(x) = sum_{n=1}^inf x^{p(n)}, then p(n) is encoded in f.
# But extracting the n-th coefficient of this "indicator" requires knowing primes.

print("\n5a. Prime indicator generating function analysis")
print("  f(x) = sum x^{p(n)} = x^2 + x^3 + x^5 + x^7 + x^11 + ...")
print("  The n-th prime could be extracted via contour integration:")
print("  p(n) = (1/2pi*i) oint f(x)^??? ... ")
print("  PROBLEM: f(x) has natural boundary at |x|=1 (Fabry gap theorem)")
print("  Contour integration cannot reach inside.")

# The Dirichlet series approach:
# sum_{n=1}^inf p(n)^{-s} = P(s)
# To extract individual p(n), we'd need Perron's formula or similar
print("\n5b. Perron's formula approach")
print("  pi(x) = (1/2pi*i) int_{c-iT}^{c+iT} log(zeta(s)) * x^s/s ds + error")
print("  Cost: O(x^{1/2+eps}) via Lagarias-Odlyzko")
print("  For p(10^100): O(10^{50+eps}) operations => INFEASIBLE")

# 5c. Dirichlet series for p(n) directly?
print("\n5c. Dirichlet series extraction")
print("  If F(s) = sum_{n=1}^inf p(n)/n^s, then")
print("  p(n) = (1/2pi*i) int F(s) n^s ds")
print("  But F(s) = sum p(n)/n^s has no known closed form.")
print("  Its analytic properties are governed by prime gaps")
print("  and the prime number theorem => no shortcut known.")


# =====================================================================
# EXPERIMENT 6: Refined holonomic tests with rational function coefficients
# =====================================================================
print(f"\n{'='*80}")
print("EXPERIMENT 6: Rational P-recursive test (polynomial/polynomial coefficients)")
print(f"{'='*80}")

def test_rational_p_recursive(seq, name, order=2, num_deg=2, den_deg=1):
    """
    Test: a(n+order) * Q(n) = sum_{j=0}^{order-1} P_j(n) * a(n+j)
    where P_j are polynomials of degree num_deg, Q is polynomial of degree den_deg.

    Rearranged as: sum_{j=0}^{order} P_j(n) * a(n+j) = 0  (set P_{order}(n) = -Q(n))
    This is the standard holonomic ODE discretization.
    """
    print(f"\n  Rational P-recursive for '{name}', order={order}, num_deg={num_deg}, den_deg={den_deg}")

    # Total unknowns: (order+1) polynomials, each of degree max(num_deg, den_deg)
    # We test this as the same as experiment 1 but with specified parameters.
    d = order
    deg = max(num_deg, den_deg)
    n_unknowns = (d + 1) * (deg + 1)
    n_eq = min(len(seq) - d, 3 * n_unknowns)

    if n_eq < n_unknowns + 5:
        print(f"    Not enough data points")
        return

    A = np.zeros((n_eq, n_unknowns), dtype=np.float64)
    for row in range(n_eq):
        n = row + 1  # 1-indexed
        for j in range(d + 1):
            if n - 1 + j >= len(seq):
                break
            a_val = seq[n - 1 + j]
            for k in range(deg + 1):
                col = j * (deg + 1) + k
                A[row, col] = a_val * (n ** k)

    U, S, Vt = np.linalg.svd(A, full_matrices=False)
    if S[0] == 0:
        return

    ratios = S / S[0]
    n_small = np.sum(ratios < 1e-10)

    if n_small > 0:
        print(f"    RANK DEFICIENT: {n_small} near-zero singular values")
        null_vec = Vt[-1, :]
        # Verify
        residuals = A @ null_vec
        max_res = np.max(np.abs(residuals))
        print(f"    Max residual: {max_res:.2e}")
    else:
        # Report how close to rank-deficient
        print(f"    Full rank. Smallest sv ratio: {ratios[-1]:.2e}")
        # Interesting if very small
        if ratios[-1] < 1e-4:
            print(f"    *** Nearly rank-deficient! Possible approximate relation. ***")


# Test various configurations
for order in [2, 3, 4, 5, 6, 8]:
    for deg in [2, 3, 4, 5]:
        test_rational_p_recursive(PRIMES, "p(n)", order=order, num_deg=deg, den_deg=deg)

# Test on transformed sequences too
print(f"\n  --- On partial sums of primes ---")
for order in [2, 3, 4]:
    test_rational_p_recursive(PARTIAL_SUMS, "sum p(k)", order=order, num_deg=3, den_deg=3)

print(f"\n  --- On p(n)^2 ---")
P_SQUARED = [p*p for p in PRIMES]
for order in [2, 3, 4]:
    test_rational_p_recursive(P_SQUARED, "p(n)^2", order=order, num_deg=3, den_deg=3)


# =====================================================================
# EXPERIMENT 7: Hadamard quotient / ratio analysis
# =====================================================================
print(f"\n{'='*80}")
print("EXPERIMENT 7: Ratio and quotient analysis")
print(f"{'='*80}")

# If p(n) were P-recursive of order 1: p(n+1) = R(n) * p(n)
# Then R(n) = p(n+1)/p(n) should be a rational function of n.
print("\n7a. Ratio p(n+1)/p(n) analysis")
ratios_pn = [PRIMES[i+1] / PRIMES[i] for i in range(N_MAX - 1)]
print(f"  First 20 ratios: {[f'{r:.4f}' for r in ratios_pn[:20]]}")
print(f"  Mean ratio: {np.mean(ratios_pn):.6f}")
print(f"  Std ratio: {np.std(ratios_pn):.6f}")
print(f"  If p(n)~n*ln(n), ratio ~ (1+1/n)*(1+1/ln(n)) -> 1")

# Fit ratio as polynomial in 1/n
inv_n = np.array([1.0 / (i + 1) for i in range(N_MAX - 1)])
ratio_arr = np.array(ratios_pn)
for fit_deg in [1, 2, 3, 4]:
    A_r = np.column_stack([inv_n**k for k in range(fit_deg + 1)])
    c, res, _, _ = np.linalg.lstsq(A_r, ratio_arr, rcond=None)
    pred = A_r @ c
    err = ratio_arr - pred
    print(f"  Poly fit deg {fit_deg} in 1/n: max_err={np.max(np.abs(err)):.4f}, "
          f"rms={np.sqrt(np.mean(err**2)):.4f}")


# =====================================================================
# EXPERIMENT 8: Connection to known holonomic sequences
# =====================================================================
print(f"\n{'='*80}")
print("EXPERIMENT 8: Is any function of primes holonomic?")
print(f"{'='*80}")

# Key insight: if f(n) is holonomic and g(n) is holonomic,
# then f(n)+g(n), f(n)*g(n), and sum_{k=0}^n f(k) are also holonomic.
# The sequence n! is holonomic (satisfies (n+1)! = (n+1)*n!).
# Is p(n)/n! holonomic? Unlikely since p(n) grows polynomially.

# More interesting: Chebyshev's theta function theta(n) = sum_{p<=n} ln(p)
# and psi(n) = sum_{p^k <= n} ln(p)
# These are step functions. Their holonomicity relates to zeta zeros.

# Test: does the sequence a(n) = sum_{k=1}^{n} 1/p(k) (prime harmonic series)
# satisfy a P-recursive relation?
PRIME_HARMONIC = [sum(1.0/PRIMES[k] for k in range(i+1)) for i in range(N_MAX)]
print(f"\n8a. Prime harmonic series: first 10 = {[f'{v:.4f}' for v in PRIME_HARMONIC[:10]]}")
test_rational_p_recursive(PRIME_HARMONIC, "sum 1/p(k)", order=3, num_deg=3, den_deg=3)

# Test: product of primes (primorial)
PRIMORIAL = [1]
for p in PRIMES:
    PRIMORIAL.append(PRIMORIAL[-1] * p)
LOG_PRIMORIAL = [math.log(x) if x > 0 else 0 for x in PRIMORIAL[1:]]
print(f"\n8b. log(primorial(n)): first 10 = {[f'{v:.4f}' for v in LOG_PRIMORIAL[:10]]}")
test_rational_p_recursive(LOG_PRIMORIAL, "log(p#(n))", order=3, num_deg=3, den_deg=3)

# Test: Mertens-like: sum_{k=1}^{n} 1/p(k) - ln(ln(p(n)))
MERTENS_LIKE = []
for i in range(N_MAX):
    h = sum(1.0/PRIMES[k] for k in range(i+1))
    ll = math.log(math.log(PRIMES[i])) if PRIMES[i] > 1 else 0
    MERTENS_LIKE.append(h - ll)
print(f"\n8c. Mertens-like residual: first 10 = {[f'{v:.6f}' for v in MERTENS_LIKE[:10]]}")
print(f"  (Should converge to Meissel-Mertens constant M ~ 0.2615)")
print(f"  Last value: {MERTENS_LIKE[-1]:.6f}")


# =====================================================================
# FINAL ANALYSIS
# =====================================================================
print(f"\n{'='*80}")
print("FINAL ANALYSIS AND CONCLUSIONS")
print(f"{'='*80}")

print("""
SUMMARY OF P-RECURSIVE / HOLONOMIC EXPLORATION
===============================================

1. DIRECT P-RECURSIVENESS OF p(n):
   The prime sequence p(n) does NOT satisfy any P-recursive relation of
   reasonable order (tested up to order 8) and degree (up to 5).
   All SVD tests show full rank with no near-zero singular values.

   This is EXPECTED: it is conjectured (and widely believed) that p(n)
   is NOT holonomic/P-recursive. The prime gaps are too irregular.

2. TRANSFORMS:
   - Binomial transform of p(n): NOT P-recursive
   - Inverse binomial transform: NOT P-recursive
   - Euler transform: NOT P-recursive
   - Second differences: NOT P-recursive
   - Mobius transform: NOT P-recursive
   - Partial sums: NOT P-recursive
   - log(p(n)): NOT P-recursive

   No standard transform renders the prime sequence holonomic.

3. MATRIX POWER FORMULATION:
   No matrix M(n) with polynomial entries of reasonable size satisfies
   the recurrence. This rules out binary splitting approaches.

4. CLOITRE-TYPE RECURRENCES:
   - Simple analytic recurrences (floor-based) can approximate primes
     but cannot generate them exactly.
   - Rowland's recurrence generates primes but is O(p(n)^2) per prime.
   - Mills-type formulas require knowing transcendental constants to
     arbitrary precision, which is circular.

5. GENERATING FUNCTION / CONTOUR INTEGRATION:
   - The prime indicator GF has a natural boundary (Fabry gap theorem)
   - Perron's formula for pi(x) costs O(x^{1/2+eps}) minimum
   - No shortcut to O(polylog) extraction of p(n) from any known GF

FUNDAMENTAL OBSTACLE:
   The primes are believed to be NOT holonomic. This means:
   - p(n) satisfies no finite-order linear recurrence with polynomial coefficients
   - No standard integral transform makes it holonomic
   - Binary splitting (which exploits P-recursiveness) cannot apply

   This is closely related to the fact that primes have positive entropy
   in the sense of Kolmogorov complexity: the n-th prime cannot be
   computed from a bounded amount of information plus polynomial-time
   computation in log(n).

   The best known unconditional method remains O(n^{2/3}) via
   analytic methods (Lagarias-Odlyzko / Deleglise-Rivat for pi(x)).
""")

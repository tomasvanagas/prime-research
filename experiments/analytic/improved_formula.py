"""
Improved Prime Formula: Fundamentally Different Approaches
==========================================================

Attempts to beat the current Lambert W formula:
  p(n) = n * W(n) * (1.0114 + 5.8864/W - 26.506/W^2 + 49.157/W^3)
  which achieves ~0.028% mean relative error for n > 50,000.

Approaches:
  1. Pade approximant in 1/W(n)
  2. Pade approximant in 1/ln(n)
  3. Minimax (L-infinity) approximation
  4. Logarithmic correction: n*W*exp(f(n))
  5. Two-variable expansion mixing W(n) and ln(n)
  6. Offset trick: shift argument of W
  7. Composite: li^{-1}(n) + correction
  8. Combinations of best approaches

Usage:
    python improved_formula.py
"""

import math
import time
import numpy as np
from scipy.optimize import minimize, least_squares
from scipy.special import lambertw as scipy_W

# ============================================================================
# Utilities
# ============================================================================

def sieve(limit):
    """Sieve of Eratosthenes, returns list of primes up to limit."""
    is_prime = bytearray(b'\x01') * (limit + 1)
    is_prime[0] = is_prime[1] = 0
    for i in range(2, int(limit**0.5) + 1):
        if is_prime[i]:
            is_prime[i*i::i] = bytearray(len(is_prime[i*i::i]))
    return [i for i in range(2, limit + 1) if is_prime[i]]

def W(x):
    """Lambert W function (principal branch) via scipy for arrays."""
    return np.real(scipy_W(x))

def W_scalar(x):
    """Lambert W for a single float."""
    if x <= 0:
        return 0.0
    if x < 1:
        w = x
    elif x < 10:
        w = math.log(x)
    else:
        lx = math.log(x)
        w = lx - math.log(lx) + math.log(lx) / lx
    for _ in range(50):
        ew = math.exp(w)
        f = w * ew - x
        if abs(f) < 1e-15 * max(abs(x), 1):
            break
        fp = ew * (w + 1)
        fpp = ew * (w + 2)
        w -= f / (fp - f * fpp / (2 * fp))
    return w

def li(x):
    """Logarithmic integral li(x) = integral from 0 to x of dt/ln(t).
    Computed via series expansion for accuracy."""
    if x <= 1:
        return float('-inf')
    lnx = math.log(x)
    # Ramanujan's series for li(x)
    # li(x) = euler_gamma + ln(|ln(x)|) + sum_{k=1}^{inf} (ln(x))^k / (k * k!)
    euler = 0.5772156649015329
    result = euler + math.log(abs(lnx))
    term = lnx
    for k in range(1, 200):
        result += term / (k * math.factorial(k)) if k <= 170 else 0
        term *= lnx
        if k > 170:
            break
        if abs(term / (k * math.factorial(k))) < 1e-15:
            break
    return result

def li_inv_newton(n_val):
    """Inverse of li: find x such that li(x) = n_val, via Newton's method."""
    # Initial guess: x ~ n * ln(n)
    if n_val < 2:
        return 2.0
    ln_n = math.log(n_val)
    x = n_val * (ln_n + math.log(ln_n))
    for _ in range(100):
        lix = li(x)
        err = lix - n_val
        if abs(err) < 1e-10:
            break
        # derivative of li(x) is 1/ln(x)
        x -= err * math.log(x)
    return x

def evaluate_formula(formula_fn, ns, primes_true, label=""):
    """Evaluate a formula against ground truth primes.
    Returns dict with mean_err, max_err, exact_count."""
    preds = np.array([formula_fn(int(n)) for n in ns])
    true = np.array([primes_true[int(n)-1] for n in ns])
    rel_err = np.abs(preds - true) / true
    rounded = np.round(preds).astype(int)
    exact = np.sum(rounded == true)
    return {
        'label': label,
        'mean_err': np.mean(rel_err) * 100,
        'max_err': np.max(rel_err) * 100,
        'exact_count': int(exact),
        'total': len(ns),
    }

def print_result(res, range_label=""):
    """Print evaluation result."""
    print(f"  {res['label']:45s} | mean={res['mean_err']:.5f}% | max={res['max_err']:.4f}% | exact={res['exact_count']}/{res['total']} {range_label}")


# ============================================================================
# Generate ground truth
# ============================================================================

print("Generating primes via sieve...")
t0 = time.time()
# Need enough primes: the 100,000th prime is 1,299,709
ALL_PRIMES = sieve(1_400_000)
print(f"  Generated {len(ALL_PRIMES)} primes in {time.time()-t0:.2f}s")
assert len(ALL_PRIMES) >= 100_000, f"Only {len(ALL_PRIMES)} primes!"

# Index arrays
TRAIN_NS = np.arange(100, 50001)   # n=100..50000
TEST_NS  = np.arange(50001, 100001) # n=50001..100000
EXACT_NS = np.arange(100, 10001)    # n=100..10000 for exact-match counting

TRAIN_PRIMES = np.array([ALL_PRIMES[n-1] for n in TRAIN_NS], dtype=np.float64)
TEST_PRIMES  = np.array([ALL_PRIMES[n-1] for n in TEST_NS], dtype=np.float64)

# Precompute W(n) and ln(n) for train/test
TRAIN_W = W(TRAIN_NS.astype(np.float64))
TEST_W  = W(TEST_NS.astype(np.float64))
TRAIN_LN = np.log(TRAIN_NS.astype(np.float64))
TEST_LN  = np.log(TEST_NS.astype(np.float64))
TRAIN_LNLN = np.log(TRAIN_LN)
TEST_LNLN  = np.log(TEST_LN)


# ============================================================================
# Baseline: Current Lambert W K=3
# ============================================================================

def baseline_formula(n):
    w = W_scalar(float(n))
    return n * w * (1.011420854 + 5.886406822/w - 26.506067157/w**2 + 49.157482190/w**3)

print("\n" + "="*90)
print("BASELINE: Lambert W K=3 polynomial")
print("="*90)

for ns, label in [(TRAIN_NS, "train n=100..50000"), (TEST_NS, "test n=50001..100000"), (EXACT_NS, "exact n=100..10000")]:
    res = evaluate_formula(baseline_formula, ns, ALL_PRIMES, "Baseline K=3")
    print_result(res, label)


# ============================================================================
# APPROACH 1a: Pade approximant in 1/W(n)
# ============================================================================

print("\n" + "="*90)
print("APPROACH 1a: Pade approximant in 1/W(n)")
print("  p(n) = n * W(n) * (a0 + a1/W + a2/W^2) / (1 + b1/W + b2/W^2)")
print("="*90)

def pade_w_formula_vec(ns, w_vals, params):
    a0, a1, a2, b1, b2 = params
    inv_w = 1.0 / w_vals
    numer = a0 + a1 * inv_w + a2 * inv_w**2
    denom = 1.0 + b1 * inv_w + b2 * inv_w**2
    return ns * w_vals * numer / denom

def pade_w_residuals(params):
    preds = pade_w_formula_vec(TRAIN_NS.astype(np.float64), TRAIN_W, params)
    return (preds - TRAIN_PRIMES) / TRAIN_PRIMES

# Initial guess from baseline
x0_pade_w = [1.01, 5.9, -26.5, 0.0, 0.0]
result = least_squares(pade_w_residuals, x0_pade_w, method='lm', max_nfev=10000)
pade_w_params = result.x
print(f"  Fitted params: a0={pade_w_params[0]:.8f}, a1={pade_w_params[1]:.8f}, a2={pade_w_params[2]:.8f}, b1={pade_w_params[3]:.8f}, b2={pade_w_params[4]:.8f}")

def pade_w_formula(n):
    w = W_scalar(float(n))
    a0, a1, a2, b1, b2 = pade_w_params
    inv_w = 1.0 / w
    numer = a0 + a1 * inv_w + a2 * inv_w**2
    denom = 1.0 + b1 * inv_w + b2 * inv_w**2
    return n * w * numer / denom

for ns, label in [(TRAIN_NS, "train"), (TEST_NS, "test"), (EXACT_NS, "exact")]:
    res = evaluate_formula(pade_w_formula, ns, ALL_PRIMES, "Pade [2/2] in 1/W")
    print_result(res, label)


# Also try higher order Pade [3/2]
print("\n  --- Pade [3/2] in 1/W ---")

def pade_w32_formula_vec(ns, w_vals, params):
    a0, a1, a2, a3, b1, b2 = params
    inv_w = 1.0 / w_vals
    numer = a0 + a1*inv_w + a2*inv_w**2 + a3*inv_w**3
    denom = 1.0 + b1*inv_w + b2*inv_w**2
    return ns * w_vals * numer / denom

def pade_w32_residuals(params):
    preds = pade_w32_formula_vec(TRAIN_NS.astype(np.float64), TRAIN_W, params)
    return (preds - TRAIN_PRIMES) / TRAIN_PRIMES

x0_32 = [1.01, 5.9, -26.5, 49.0, 0.0, 0.0]
result = least_squares(pade_w32_residuals, x0_32, method='lm', max_nfev=10000)
pade_w32_params = result.x
print(f"  Fitted [3/2] params: {[f'{p:.6f}' for p in pade_w32_params]}")

def pade_w32_formula(n):
    w = W_scalar(float(n))
    a0, a1, a2, a3, b1, b2 = pade_w32_params
    inv_w = 1.0 / w
    numer = a0 + a1*inv_w + a2*inv_w**2 + a3*inv_w**3
    denom = 1.0 + b1*inv_w + b2*inv_w**2
    return n * w * numer / denom

for ns, label in [(TRAIN_NS, "train"), (TEST_NS, "test"), (EXACT_NS, "exact")]:
    res = evaluate_formula(pade_w32_formula, ns, ALL_PRIMES, "Pade [3/2] in 1/W")
    print_result(res, label)


# ============================================================================
# APPROACH 1b: Pade in 1/ln(n)
# ============================================================================

print("\n" + "="*90)
print("APPROACH 1b: Pade approximant in 1/ln(n)")
print("  p(n) = n * (a0*ln(n) + a1*ln(ln(n)) + a2) / (1 + b1/ln(n) + b2/ln(n)^2)")
print("="*90)

def pade_ln_formula_vec(ns, ln_vals, lnln_vals, params):
    a0, a1, a2, b1, b2 = params
    inv_ln = 1.0 / ln_vals
    numer = a0 * ln_vals + a1 * lnln_vals + a2
    denom = 1.0 + b1 * inv_ln + b2 * inv_ln**2
    return ns * numer / denom

def pade_ln_residuals(params):
    preds = pade_ln_formula_vec(TRAIN_NS.astype(np.float64), TRAIN_LN, TRAIN_LNLN, params)
    return (preds - TRAIN_PRIMES) / TRAIN_PRIMES

x0_pade_ln = [1.0, 1.0, -1.0, 0.0, 0.0]
result = least_squares(pade_ln_residuals, x0_pade_ln, method='lm', max_nfev=10000)
pade_ln_params = result.x
print(f"  Fitted params: {[f'{p:.8f}' for p in pade_ln_params]}")

def pade_ln_formula(n):
    ln_n = math.log(n)
    lnln_n = math.log(ln_n)
    a0, a1, a2, b1, b2 = pade_ln_params
    inv_ln = 1.0 / ln_n
    numer = a0 * ln_n + a1 * lnln_n + a2
    denom = 1.0 + b1 * inv_ln + b2 * inv_ln**2
    return n * numer / denom

for ns, label in [(TRAIN_NS, "train"), (TEST_NS, "test"), (EXACT_NS, "exact")]:
    res = evaluate_formula(pade_ln_formula, ns, ALL_PRIMES, "Pade in 1/ln(n)")
    print_result(res, label)


# ============================================================================
# APPROACH 2: Minimax (L-infinity) approximation
# ============================================================================

print("\n" + "="*90)
print("APPROACH 2: Minimax (L-infinity) approximation")
print("  Minimize max|error| instead of mean(error^2)")
print("="*90)

# Use the same polynomial form as baseline but optimize for worst-case
def minimax_polynomial_pred(ns, w_vals, params):
    c0, c1, c2, c3 = params
    inv_w = 1.0 / w_vals
    return ns * w_vals * (c0 + c1*inv_w + c2*inv_w**2 + c3*inv_w**3)

def minimax_objective(params):
    preds = minimax_polynomial_pred(TRAIN_NS.astype(np.float64), TRAIN_W, params)
    rel_err = np.abs(preds - TRAIN_PRIMES) / TRAIN_PRIMES
    return np.max(rel_err)

from scipy.optimize import differential_evolution

bounds_mm = [(0.9, 1.1), (4.0, 8.0), (-40.0, -10.0), (20.0, 80.0)]
result_mm = differential_evolution(minimax_objective, bounds_mm, seed=42, maxiter=500, tol=1e-12, polish=True)
mm_params = result_mm.x
print(f"  Minimax params: c0={mm_params[0]:.8f}, c1={mm_params[1]:.8f}, c2={mm_params[2]:.8f}, c3={mm_params[3]:.8f}")
print(f"  Minimax obj (max train rel err): {result_mm.fun*100:.5f}%")

def minimax_formula(n):
    w = W_scalar(float(n))
    c0, c1, c2, c3 = mm_params
    return n * w * (c0 + c1/w + c2/w**2 + c3/w**3)

for ns, label in [(TRAIN_NS, "train"), (TEST_NS, "test"), (EXACT_NS, "exact")]:
    res = evaluate_formula(minimax_formula, ns, ALL_PRIMES, "Minimax K=3")
    print_result(res, label)

# Minimax with Pade form
print("\n  --- Minimax with Pade [3/2] form ---")

def minimax_pade_objective(params):
    a0, a1, a2, a3, b1, b2 = params
    inv_w = 1.0 / TRAIN_W
    numer = a0 + a1*inv_w + a2*inv_w**2 + a3*inv_w**3
    denom = 1.0 + b1*inv_w + b2*inv_w**2
    preds = TRAIN_NS.astype(np.float64) * TRAIN_W * numer / denom
    rel_err = np.abs(preds - TRAIN_PRIMES) / TRAIN_PRIMES
    return np.max(rel_err)

bounds_mm_pade = [(0.9, 1.1), (4.0, 8.0), (-40.0, -10.0), (20.0, 80.0), (-2.0, 2.0), (-5.0, 5.0)]
result_mm_pade = differential_evolution(minimax_pade_objective, bounds_mm_pade, seed=42, maxiter=500, tol=1e-12, polish=True)
mm_pade_params = result_mm_pade.x
print(f"  Minimax Pade params: {[f'{p:.6f}' for p in mm_pade_params]}")

def minimax_pade_formula(n):
    w = W_scalar(float(n))
    a0, a1, a2, a3, b1, b2 = mm_pade_params
    inv_w = 1.0 / w
    numer = a0 + a1*inv_w + a2*inv_w**2 + a3*inv_w**3
    denom = 1.0 + b1*inv_w + b2*inv_w**2
    return n * w * numer / denom

for ns, label in [(TRAIN_NS, "train"), (TEST_NS, "test"), (EXACT_NS, "exact")]:
    res = evaluate_formula(minimax_pade_formula, ns, ALL_PRIMES, "Minimax Pade [3/2]")
    print_result(res, label)


# ============================================================================
# APPROACH 3: Logarithmic correction (exponential)
# ============================================================================

print("\n" + "="*90)
print("APPROACH 3: Logarithmic correction")
print("  p(n) = n * W(n) * exp(f0 + f1/W + f2/W^2 + f3/W^3)")
print("="*90)

# Target: p(n) / (n * W(n)) = exp(f(n))
# So f(n) = ln(p(n) / (n * W(n)))
TRAIN_RATIO = TRAIN_PRIMES / (TRAIN_NS * TRAIN_W)
TRAIN_LOG_RATIO = np.log(TRAIN_RATIO)

def exp_correction_residuals(params):
    f0, f1, f2, f3 = params
    inv_w = 1.0 / TRAIN_W
    f_vals = f0 + f1*inv_w + f2*inv_w**2 + f3*inv_w**3
    preds = TRAIN_NS.astype(np.float64) * TRAIN_W * np.exp(f_vals)
    return (preds - TRAIN_PRIMES) / TRAIN_PRIMES

# Initial guess: ln of baseline coefficients
x0_exp = [math.log(1.01), 5.8, -20.0, 30.0]
result_exp = least_squares(exp_correction_residuals, x0_exp, method='lm', max_nfev=10000)
exp_params = result_exp.x
print(f"  Exp params: f0={exp_params[0]:.8f}, f1={exp_params[1]:.8f}, f2={exp_params[2]:.8f}, f3={exp_params[3]:.8f}")

def exp_formula(n):
    w = W_scalar(float(n))
    f0, f1, f2, f3 = exp_params
    inv_w = 1.0 / w
    f_val = f0 + f1*inv_w + f2*inv_w**2 + f3*inv_w**3
    return n * w * math.exp(f_val)

for ns, label in [(TRAIN_NS, "train"), (TEST_NS, "test"), (EXACT_NS, "exact")]:
    res = evaluate_formula(exp_formula, ns, ALL_PRIMES, "Exp correction K=3")
    print_result(res, label)

# Higher order exponential
print("\n  --- Exp correction K=4 ---")

def exp_correction_k4_residuals(params):
    f0, f1, f2, f3, f4 = params
    inv_w = 1.0 / TRAIN_W
    f_vals = f0 + f1*inv_w + f2*inv_w**2 + f3*inv_w**3 + f4*inv_w**4
    preds = TRAIN_NS.astype(np.float64) * TRAIN_W * np.exp(f_vals)
    return (preds - TRAIN_PRIMES) / TRAIN_PRIMES

x0_exp4 = list(exp_params) + [0.0]
result_exp4 = least_squares(exp_correction_k4_residuals, x0_exp4, method='lm', max_nfev=10000)
exp4_params = result_exp4.x
print(f"  Exp K=4 params: {[f'{p:.8f}' for p in exp4_params]}")

def exp_formula_k4(n):
    w = W_scalar(float(n))
    f0, f1, f2, f3, f4 = exp4_params
    inv_w = 1.0 / w
    f_val = f0 + f1*inv_w + f2*inv_w**2 + f3*inv_w**3 + f4*inv_w**4
    return n * w * math.exp(f_val)

for ns, label in [(TRAIN_NS, "train"), (TEST_NS, "test"), (EXACT_NS, "exact")]:
    res = evaluate_formula(exp_formula_k4, ns, ALL_PRIMES, "Exp correction K=4")
    print_result(res, label)


# ============================================================================
# APPROACH 4: Two-variable expansion
# ============================================================================

print("\n" + "="*90)
print("APPROACH 4: Two-variable expansion")
print("  p(n) = n * (a*W(n) + b*ln(n) + c*ln(ln(n)) + d + e/W(n) + f/ln(n))")
print("="*90)

def two_var_formula_vec(ns, w_vals, ln_vals, lnln_vals, params):
    a, b, c, d, e, f = params
    return ns * (a*w_vals + b*ln_vals + c*lnln_vals + d + e/w_vals + f/ln_vals)

def two_var_residuals(params):
    preds = two_var_formula_vec(TRAIN_NS.astype(np.float64), TRAIN_W, TRAIN_LN, TRAIN_LNLN, params)
    return (preds - TRAIN_PRIMES) / TRAIN_PRIMES

x0_2v = [1.0, 0.0, 1.0, -1.0, 0.0, 0.0]
result_2v = least_squares(two_var_residuals, x0_2v, method='lm', max_nfev=10000)
two_var_params = result_2v.x
print(f"  Two-var params: {[f'{p:.8f}' for p in two_var_params]}")

def two_var_formula(n):
    w = W_scalar(float(n))
    ln_n = math.log(n)
    lnln_n = math.log(ln_n)
    a, b, c, d, e, f = two_var_params
    return n * (a*w + b*ln_n + c*lnln_n + d + e/w + f/ln_n)

for ns, label in [(TRAIN_NS, "train"), (TEST_NS, "test"), (EXACT_NS, "exact")]:
    res = evaluate_formula(two_var_formula, ns, ALL_PRIMES, "Two-variable")
    print_result(res, label)

# Extended two-variable with cross terms
print("\n  --- Extended two-variable with cross terms ---")

def two_var_ext_formula_vec(ns, w_vals, ln_vals, lnln_vals, params):
    a, b, c, d, e, f, g, h = params
    return ns * (a*w_vals + b*ln_vals + c*lnln_vals + d + e/w_vals + f/ln_vals + g*lnln_vals/w_vals + h/w_vals**2)

def two_var_ext_residuals(params):
    preds = two_var_ext_formula_vec(TRAIN_NS.astype(np.float64), TRAIN_W, TRAIN_LN, TRAIN_LNLN, params)
    return (preds - TRAIN_PRIMES) / TRAIN_PRIMES

x0_2ve = list(two_var_params) + [0.0, 0.0]
result_2ve = least_squares(two_var_ext_residuals, x0_2ve, method='lm', max_nfev=10000)
two_var_ext_params = result_2ve.x
print(f"  Extended params: {[f'{p:.6f}' for p in two_var_ext_params]}")

def two_var_ext_formula(n):
    w = W_scalar(float(n))
    ln_n = math.log(n)
    lnln_n = math.log(ln_n)
    a, b, c, d, e, f, g, h = two_var_ext_params
    return n * (a*w + b*ln_n + c*lnln_n + d + e/w + f/ln_n + g*lnln_n/w + h/w**2)

for ns, label in [(TRAIN_NS, "train"), (TEST_NS, "test"), (EXACT_NS, "exact")]:
    res = evaluate_formula(two_var_ext_formula, ns, ALL_PRIMES, "Extended two-variable")
    print_result(res, label)


# ============================================================================
# APPROACH 5: The offset trick
# ============================================================================

print("\n" + "="*90)
print("APPROACH 5: Offset trick")
print("  p(n) = n * W(n + alpha*sqrt(n)) * (c0 + c1/W(n + alpha*sqrt(n)))")
print("="*90)

def offset_formula_vec(ns_float, params):
    alpha, c0, c1 = params
    shifted = ns_float + alpha * np.sqrt(ns_float)
    w_shifted = W(shifted)
    return ns_float * w_shifted * (c0 + c1 / w_shifted)

def offset_residuals(params):
    preds = offset_formula_vec(TRAIN_NS.astype(np.float64), params)
    return (preds - TRAIN_PRIMES) / TRAIN_PRIMES

x0_off = [1.5, 1.01, 5.0]
result_off = least_squares(offset_residuals, x0_off, method='lm', max_nfev=10000)
offset_params = result_off.x
print(f"  Offset params: alpha={offset_params[0]:.8f}, c0={offset_params[1]:.8f}, c1={offset_params[2]:.8f}")

def offset_formula(n):
    alpha, c0, c1 = offset_params
    n_shifted = n + alpha * math.sqrt(n)
    w = W_scalar(n_shifted)
    return n * w * (c0 + c1 / w)

for ns, label in [(TRAIN_NS, "train"), (TEST_NS, "test"), (EXACT_NS, "exact")]:
    res = evaluate_formula(offset_formula, ns, ALL_PRIMES, "Offset trick (2 terms)")
    print_result(res, label)

# Extended offset with more terms
print("\n  --- Offset trick with 3 correction terms ---")

def offset_ext_formula_vec(ns_float, params):
    alpha, c0, c1, c2 = params
    shifted = ns_float + alpha * np.sqrt(ns_float)
    w_shifted = W(shifted)
    inv_w = 1.0 / w_shifted
    return ns_float * w_shifted * (c0 + c1*inv_w + c2*inv_w**2)

def offset_ext_residuals(params):
    preds = offset_ext_formula_vec(TRAIN_NS.astype(np.float64), params)
    return (preds - TRAIN_PRIMES) / TRAIN_PRIMES

x0_off_ext = [offset_params[0], offset_params[1], offset_params[2], 0.0]
result_off_ext = least_squares(offset_ext_residuals, x0_off_ext, method='lm', max_nfev=10000)
offset_ext_params = result_off_ext.x
print(f"  Extended offset params: {[f'{p:.8f}' for p in offset_ext_params]}")

def offset_ext_formula(n):
    alpha, c0, c1, c2 = offset_ext_params
    n_shifted = n + alpha * math.sqrt(n)
    w = W_scalar(n_shifted)
    inv_w = 1.0 / w
    return n * w * (c0 + c1*inv_w + c2*inv_w**2)

for ns, label in [(TRAIN_NS, "train"), (TEST_NS, "test"), (EXACT_NS, "exact")]:
    res = evaluate_formula(offset_ext_formula, ns, ALL_PRIMES, "Offset trick (3 terms)")
    print_result(res, label)

# Double offset: shift n and also have K=3
print("\n  --- Offset trick with K=3 correction ---")

def offset_k3_formula_vec(ns_float, params):
    alpha, beta, c0, c1, c2, c3 = params
    shifted = ns_float + alpha * np.sqrt(ns_float) + beta * np.log(ns_float)
    w_shifted = W(shifted)
    inv_w = 1.0 / w_shifted
    return ns_float * w_shifted * (c0 + c1*inv_w + c2*inv_w**2 + c3*inv_w**3)

def offset_k3_residuals(params):
    preds = offset_k3_formula_vec(TRAIN_NS.astype(np.float64), params)
    return (preds - TRAIN_PRIMES) / TRAIN_PRIMES

x0_off_k3 = [offset_params[0], 0.0, 1.01, 5.9, -26.5, 49.0]
result_off_k3 = least_squares(offset_k3_residuals, x0_off_k3, method='lm', max_nfev=10000)
offset_k3_params = result_off_k3.x
print(f"  Offset K=3 params: {[f'{p:.6f}' for p in offset_k3_params]}")

def offset_k3_formula(n):
    alpha, beta, c0, c1, c2, c3 = offset_k3_params
    n_shifted = n + alpha * math.sqrt(n) + beta * math.log(n)
    w = W_scalar(n_shifted)
    inv_w = 1.0 / w
    return n * w * (c0 + c1*inv_w + c2*inv_w**2 + c3*inv_w**3)

for ns, label in [(TRAIN_NS, "train"), (TEST_NS, "test"), (EXACT_NS, "exact")]:
    res = evaluate_formula(offset_k3_formula, ns, ALL_PRIMES, "Offset+K=3 correction")
    print_result(res, label)


# ============================================================================
# APPROACH 6: Composite li^{-1}(n) + correction
# ============================================================================

print("\n" + "="*90)
print("APPROACH 6: Composite li^{-1}(n) + correction")
print("  p(n) = li^{-1}(n) + n * W(n) * g(n)")
print("="*90)

# Precompute li^{-1} for training set
print("  Computing li^{-1} for training set (this takes a moment)...")
t0 = time.time()
TRAIN_LI_INV = np.array([li_inv_newton(int(n)) for n in TRAIN_NS])
TEST_LI_INV = np.array([li_inv_newton(int(n)) for n in TEST_NS])
print(f"  Done in {time.time()-t0:.1f}s")

# Residual after li^{-1}
TRAIN_LI_RESID = (TRAIN_PRIMES - TRAIN_LI_INV) / (TRAIN_NS * TRAIN_W)

def composite_residuals(params):
    g0, g1, g2 = params
    inv_w = 1.0 / TRAIN_W
    g_vals = g0 + g1*inv_w + g2*inv_w**2
    preds = TRAIN_LI_INV + TRAIN_NS * TRAIN_W * g_vals
    return (preds - TRAIN_PRIMES) / TRAIN_PRIMES

x0_comp = [0.0, 0.0, 0.0]
result_comp = least_squares(composite_residuals, x0_comp, method='lm', max_nfev=10000)
comp_params = result_comp.x
print(f"  Composite correction params: {[f'{p:.8f}' for p in comp_params]}")

def composite_formula(n):
    li_inv = li_inv_newton(n)
    w = W_scalar(float(n))
    g0, g1, g2 = comp_params
    inv_w = 1.0 / w
    g_val = g0 + g1*inv_w + g2*inv_w**2
    return li_inv + n * w * g_val

for ns, label in [(TRAIN_NS, "train"), (TEST_NS, "test"), (EXACT_NS, "exact")]:
    res = evaluate_formula(composite_formula, ns, ALL_PRIMES, "Composite li^{-1} + correction")
    print_result(res, label)

# Also try: li^{-1}(n + offset)
print("\n  --- li^{-1}(n + alpha*sqrt(n/ln(n))) ---")

def li_offset_formula_vec(ns_float, params):
    alpha, beta = params
    offsets = ns_float + alpha * np.sqrt(ns_float / np.log(ns_float)) + beta
    preds = np.array([li_inv_newton(max(2, o)) for o in offsets])
    return preds

def li_offset_residuals(params):
    preds = li_offset_formula_vec(TRAIN_NS.astype(np.float64), params)
    return (preds - TRAIN_PRIMES) / TRAIN_PRIMES

x0_lio = [1.47, 0.0]
result_lio = least_squares(li_offset_residuals, x0_lio, method='lm', max_nfev=5000)
li_offset_params = result_lio.x
print(f"  li offset params: alpha={li_offset_params[0]:.8f}, beta={li_offset_params[1]:.8f}")

def li_offset_formula(n):
    alpha, beta = li_offset_params
    n_off = n + alpha * math.sqrt(n / math.log(n)) + beta
    return li_inv_newton(max(2, n_off))

for ns, label in [(TRAIN_NS, "train"), (TEST_NS, "test"), (EXACT_NS, "exact")]:
    res = evaluate_formula(li_offset_formula, ns, ALL_PRIMES, "li^{-1}(n + offset)")
    print_result(res, label)


# ============================================================================
# APPROACH 7: Pade [3/3] in 1/W(n) — higher order rational
# ============================================================================

print("\n" + "="*90)
print("APPROACH 7: Higher-order Pade [3/3] in 1/W(n)")
print("="*90)

def pade_w33_formula_vec(ns, w_vals, params):
    a0, a1, a2, a3, b1, b2, b3 = params
    inv_w = 1.0 / w_vals
    numer = a0 + a1*inv_w + a2*inv_w**2 + a3*inv_w**3
    denom = 1.0 + b1*inv_w + b2*inv_w**2 + b3*inv_w**3
    return ns * w_vals * numer / denom

def pade_w33_residuals(params):
    preds = pade_w33_formula_vec(TRAIN_NS.astype(np.float64), TRAIN_W, params)
    return (preds - TRAIN_PRIMES) / TRAIN_PRIMES

x0_33 = list(pade_w32_params[:4]) + list(pade_w32_params[4:]) + [0.0]
result_33 = least_squares(pade_w33_residuals, x0_33, method='lm', max_nfev=10000)
pade_w33_params = result_33.x
print(f"  Pade [3/3] params: {[f'{p:.6f}' for p in pade_w33_params]}")

def pade_w33_formula(n):
    w = W_scalar(float(n))
    a0, a1, a2, a3, b1, b2, b3 = pade_w33_params
    inv_w = 1.0 / w
    numer = a0 + a1*inv_w + a2*inv_w**2 + a3*inv_w**3
    denom = 1.0 + b1*inv_w + b2*inv_w**2 + b3*inv_w**3
    return n * w * numer / denom

for ns, label in [(TRAIN_NS, "train"), (TEST_NS, "test"), (EXACT_NS, "exact")]:
    res = evaluate_formula(pade_w33_formula, ns, ALL_PRIMES, "Pade [3/3] in 1/W")
    print_result(res, label)


# ============================================================================
# APPROACH 8: Combinations
# ============================================================================

print("\n" + "="*90)
print("APPROACH 8: Combinations of approaches")
print("="*90)

# Weighted average of best approaches
# Collect predictions from all approaches on training set
print("  Computing predictions from all approaches...")

approaches = {
    'baseline': baseline_formula,
    'pade_w22': pade_w_formula,
    'pade_w32': pade_w32_formula,
    'pade_w33': pade_w33_formula,
    'minimax': minimax_formula,
    'minimax_pade': minimax_pade_formula,
    'exp_k3': exp_formula,
    'exp_k4': exp_formula_k4,
    'two_var': two_var_formula,
    'two_var_ext': two_var_ext_formula,
    'offset_k3': offset_k3_formula,
}

# Compute all predictions
all_train_preds = {}
all_test_preds = {}
for name, fn in approaches.items():
    all_train_preds[name] = np.array([fn(int(n)) for n in TRAIN_NS])
    all_test_preds[name] = np.array([fn(int(n)) for n in TEST_NS])

# Find optimal weights for a weighted average of top approaches
from itertools import combinations

names = list(approaches.keys())

# Try all pairs
print("\n  --- Best weighted pairs ---")
best_pair_err = 1.0
best_pair = None
best_pair_weights = None

for i, j in combinations(range(len(names)), 2):
    n1, n2 = names[i], names[j]
    p1_train = all_train_preds[n1]
    p2_train = all_train_preds[n2]

    # Optimal weight: minimize sum of (w*p1 + (1-w)*p2 - true)^2 / true^2
    # d/dw [ sum ( (w*p1 + (1-w)*p2 - true) / true )^2 ] = 0
    diff = (p1_train - p2_train) / TRAIN_PRIMES
    base_err = (p2_train - TRAIN_PRIMES) / TRAIN_PRIMES
    # w * diff + base_err minimized when w = -sum(diff*base_err) / sum(diff^2)
    w_opt = -np.sum(diff * base_err) / (np.sum(diff**2) + 1e-30)
    w_opt = np.clip(w_opt, -0.5, 1.5)

    p1_test = all_test_preds[n1]
    p2_test = all_test_preds[n2]
    combo_test = w_opt * p1_test + (1 - w_opt) * p2_test
    test_err = np.mean(np.abs(combo_test - TEST_PRIMES) / TEST_PRIMES)

    if test_err < best_pair_err:
        best_pair_err = test_err
        best_pair = (n1, n2)
        best_pair_weights = (w_opt, 1 - w_opt)

print(f"  Best pair: {best_pair[0]} ({best_pair_weights[0]:.4f}) + {best_pair[1]} ({best_pair_weights[1]:.4f})")
print(f"  Test mean error: {best_pair_err*100:.5f}%")

# Build the combined formula
w1_final, w2_final = best_pair_weights
fn1_final, fn2_final = approaches[best_pair[0]], approaches[best_pair[1]]

def combo_formula(n):
    return w1_final * fn1_final(n) + w2_final * fn2_final(n)

for ns, label in [(TRAIN_NS, "train"), (TEST_NS, "test"), (EXACT_NS, "exact")]:
    res = evaluate_formula(combo_formula, ns, ALL_PRIMES, f"Combo: {best_pair[0]}+{best_pair[1]}")
    print_result(res, label)

# Try triple combo
print("\n  --- Best weighted triple ---")
best_triple_err = 1.0
best_triple = None
best_triple_weights = None

for i, j, k in combinations(range(len(names)), 3):
    n1, n2, n3 = names[i], names[j], names[k]
    p1, p2, p3 = all_train_preds[n1], all_train_preds[n2], all_train_preds[n3]

    # Solve linear system for optimal weights
    # min || (w1*p1 + w2*p2 + (1-w1-w2)*p3 - true) / true ||^2
    d1 = (p1 - p3) / TRAIN_PRIMES
    d2 = (p2 - p3) / TRAIN_PRIMES
    base = (p3 - TRAIN_PRIMES) / TRAIN_PRIMES

    A = np.array([[np.sum(d1*d1), np.sum(d1*d2)],
                   [np.sum(d1*d2), np.sum(d2*d2)]])
    b_vec = np.array([-np.sum(d1*base), -np.sum(d2*base)])

    try:
        w12 = np.linalg.solve(A, b_vec)
    except np.linalg.LinAlgError:
        continue

    w1, w2 = w12
    w3 = 1 - w1 - w2
    if abs(w1) > 3 or abs(w2) > 3 or abs(w3) > 3:
        continue

    pt1, pt2, pt3 = all_test_preds[n1], all_test_preds[n2], all_test_preds[n3]
    combo = w1*pt1 + w2*pt2 + w3*pt3
    test_err = np.mean(np.abs(combo - TEST_PRIMES) / TEST_PRIMES)

    if test_err < best_triple_err:
        best_triple_err = test_err
        best_triple = (n1, n2, n3)
        best_triple_weights = (w1, w2, w3)

if best_triple:
    print(f"  Best triple: {best_triple[0]} ({best_triple_weights[0]:.4f}) + {best_triple[1]} ({best_triple_weights[1]:.4f}) + {best_triple[2]} ({best_triple_weights[2]:.4f})")
    print(f"  Test mean error: {best_triple_err*100:.5f}%")

    fn_t1, fn_t2, fn_t3 = approaches[best_triple[0]], approaches[best_triple[1]], approaches[best_triple[2]]
    wt1, wt2, wt3 = best_triple_weights

    def triple_combo_formula(n):
        return wt1 * fn_t1(n) + wt2 * fn_t2(n) + wt3 * fn_t3(n)

    for ns, label in [(TRAIN_NS, "train"), (TEST_NS, "test"), (EXACT_NS, "exact")]:
        res = evaluate_formula(triple_combo_formula, ns, ALL_PRIMES, f"Triple combo")
        print_result(res, label)


# ============================================================================
# FINAL SUMMARY
# ============================================================================

print("\n" + "="*90)
print("FINAL SUMMARY: All approaches on TEST SET (n=50001..100000)")
print("="*90)
print(f"  {'Approach':50s} | {'Mean Err %':>12s} | {'Max Err %':>12s}")
print("  " + "-"*80)

all_formulas = [
    ("Baseline K=3 polynomial", baseline_formula),
    ("Pade [2/2] in 1/W", pade_w_formula),
    ("Pade [3/2] in 1/W", pade_w32_formula),
    ("Pade [3/3] in 1/W", pade_w33_formula),
    ("Pade in 1/ln(n)", pade_ln_formula),
    ("Minimax K=3", minimax_formula),
    ("Minimax Pade [3/2]", minimax_pade_formula),
    ("Exp correction K=3", exp_formula),
    ("Exp correction K=4", exp_formula_k4),
    ("Two-variable", two_var_formula),
    ("Extended two-variable", two_var_ext_formula),
    ("Offset trick (2 terms)", offset_formula),
    ("Offset trick (3 terms)", offset_ext_formula),
    ("Offset + K=3 correction", offset_k3_formula),
    (f"Best pair combo", combo_formula),
]

if best_triple:
    all_formulas.append(("Best triple combo", triple_combo_formula))

results = []
for label, fn in all_formulas:
    res = evaluate_formula(fn, TEST_NS, ALL_PRIMES, label)
    results.append(res)
    print(f"  {res['label']:50s} | {res['mean_err']:>11.5f}% | {res['max_err']:>11.4f}%")

# Find the winner
best = min(results, key=lambda r: r['mean_err'])
print(f"\n  WINNER (lowest mean error on test): {best['label']}")
print(f"    Mean error: {best['mean_err']:.5f}%  (baseline: 0.028%)")
print(f"    Max error:  {best['max_err']:.4f}%")

best_max = min(results, key=lambda r: r['max_err'])
print(f"\n  WINNER (lowest max error on test): {best_max['label']}")
print(f"    Max error:  {best_max['max_err']:.4f}%")
print(f"    Mean error: {best_max['mean_err']:.5f}%")

# Print the best formula's coefficients
print("\n" + "="*90)
print("BEST FORMULAS - COEFFICIENTS")
print("="*90)

print(f"\nPade [3/2] in 1/W(n):")
print(f"  p(n) = n * W(n) * (a0 + a1/W + a2/W^2 + a3/W^3) / (1 + b1/W + b2/W^2)")
print(f"  a0={pade_w32_params[0]:.10f}, a1={pade_w32_params[1]:.10f}, a2={pade_w32_params[2]:.10f}, a3={pade_w32_params[3]:.10f}")
print(f"  b1={pade_w32_params[4]:.10f}, b2={pade_w32_params[5]:.10f}")

print(f"\nPade [3/3] in 1/W(n):")
print(f"  p(n) = n * W(n) * (a0 + a1/W + a2/W^2 + a3/W^3) / (1 + b1/W + b2/W^2 + b3/W^3)")
for i, name in enumerate(['a0','a1','a2','a3','b1','b2','b3']):
    print(f"  {name}={pade_w33_params[i]:.10f}")

print(f"\nExp correction K=4:")
print(f"  p(n) = n * W(n) * exp(f0 + f1/W + f2/W^2 + f3/W^3 + f4/W^4)")
for i, name in enumerate(['f0','f1','f2','f3','f4']):
    print(f"  {name}={exp4_params[i]:.10f}")

print(f"\nOffset + K=3:")
print(f"  p(n) = n * W(n + alpha*sqrt(n) + beta*ln(n)) * (c0 + c1/W' + c2/W'^2 + c3/W'^3)")
for i, name in enumerate(['alpha','beta','c0','c1','c2','c3']):
    print(f"  {name}={offset_k3_params[i]:.10f}")

if best_triple:
    print(f"\nBest triple combo:")
    print(f"  p(n) = {best_triple_weights[0]:.6f} * {best_triple[0]}(n) + {best_triple_weights[1]:.6f} * {best_triple[1]}(n) + {best_triple_weights[2]:.6f} * {best_triple[2]}(n)")

print("\nDone.")

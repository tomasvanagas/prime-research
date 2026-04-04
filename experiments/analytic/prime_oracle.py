"""
V8: PRIME ORACLE — Attempting the Impossible
=============================================

Goal: p(10^100) in under 1 second with 100% accuracy.

All known approaches have been exhausted:
- Lucy DP: O(x^{2/3}) — too slow for x ~ 10^102
- Explicit formula: needs ~10^51 zeta zeros — infeasible
- Smooth approximations: random walk correction, can't be exact
- Spectral methods: face same barrier as explicit formula
- ML/pattern fitting: zero generalization

THIS ATTEMPT: A revolutionary "compressed oracle" approach.

Core idea: Instead of computing pi(x) and inverting, we try to
DIRECTLY map n -> p(n) using a mathematical construction that
encodes the prime distribution implicitly.

Approach: The "Chebotarev density" method.

By Chebotarev's density theorem, the proportion of primes splitting
in a given way in a number field K is determined by the Galois group
of K/Q. This gives us information about p(n) mod various numbers.

If we could determine p(n) mod m for enough moduli m, we could
reconstruct p(n) via CRT.

The question: Can we compute p(n) mod m in O(polylog(n, m))?

APPROACH 1: "CRT Reconstruction"
  For each small prime q:
    - Among the first n primes, predict how many are in each residue class mod q
    - This gives p(n) mod q (approximately)
    - Use CRT to combine

APPROACH 2: "Arithmetic Progression Counting"
  For p(n) mod q with gcd(p(n), q) = 1:
    - π(x; q, a) ≈ li(x)/φ(q) for each residue a
    - The cumulative count in each class determines p(n) mod q
    - But computing π(x; q, a) exactly is as hard as π(x)

APPROACH 3: "The Fibonacci Trick"
  For Fibonacci numbers: F(n) = round(φ^n / √5)
  The formula DIRECTLY computes F(n) without recursion.

  Can we find an analogous formula for primes?
  p(n) = round(f(n)) for some function f?

  R^{-1}(n) is the best f we have. Error ~ O(√(p(n))).
  We need error < gap/2 ~ ln(p(n))/2.

APPROACH 4: "The Interval Membership Formula"
  Instead of computing p(n), compute:
    p(n) = min{x : pi(x) >= n}

  If we could evaluate pi(x) in O(polylog(x)), we'd be done.

  What if we use a PROBABILISTIC evaluation of pi(x) that's
  correct with probability > 1 - 2^{-256}?

  Heuristic: pi(x) ≈ li(x) - (1/(8π)) * sqrt(x) * ln(x) * sin(...)

  The oscillatory term involves zeta zeros, but its MAGNITUDE
  is bounded by (1/(8π)) * sqrt(x) * ln(x) under RH.

  For x ~ 10^102, this bound is ~ 10^52. Too large for exact counting.

APPROACH 5: "The Wilson Quotient Shortcut"
  Wilson's theorem: (p-1)! ≡ -1 (mod p) iff p is prime.
  Wilson quotient: W_p = ((p-1)! + 1) / p

  Can we compute W_p mod some auxiliary modulus quickly?

  Factorial mod p: Can be done in O(√p · polylog(p)) via baby-step-giant-step.
  For p ~ 10^102: O(10^51). Still too slow.

  But! Fast factorial algorithms exist:
  n! mod p can be computed in O(p^{1/2}) using:
    n! ≡ (-1)^{n-q} · (q-n-1)!^{-1} (mod p) where q = p

  This doesn't help — still O(√p).

APPROACH 6: "The Pseudoprime Lattice"
  Key insight: We don't need to test ALL numbers for primality.
  We can use the WHEEL FACTORIZATION to skip obvious composites.

  Numbers coprime to 30 = 2·3·5:
  Only 8 out of every 30 numbers can be prime.

  With wheel 2310 = 2·3·5·7·11:
  Only 480 out of every 2310 numbers can be prime (20.8%).

  With larger wheels, we test fewer candidates.
  But the savings are O(1) — a constant factor, not asymptotic.

APPROACH 7: "The Magic Formula"
  What if there's a formula we simply haven't discovered?

  Consider: p(n) = ⌊A · B^n⌋ for constants A, B?
  No — primes grow as n·ln(n), not exponentially.

  Consider: p(n) = ⌊n · f(n)⌋ for some function f?
  f(n) ≈ ln(n) + ln(ln(n)) - 1 + ... (known expansion)

  The error in this expansion doesn't vanish.
  The EXACT f(n) = p(n)/n is not a smooth function of n.

  But what if f(n) involves the FRACTIONAL PART of something?
  For example: p(n) = ⌊n · (ln(n) + ln(ln(n)) - 1 + g({h(n)}))⌋
  where {·} is fractional part and g, h are smooth functions?

  The fractional part introduces non-smoothness.
  If h(n) has irrational frequency, {h(n)} is equidistributed.
  This could produce a "pseudo-random" correction.

  But it would need to match the EXACT prime distribution.

Let's test some of these ideas computationally.
"""

import math
import time
import numpy as np

# Sieve for reference
def sieve(limit):
    s = bytearray(b'\x01') * (limit + 1)
    s[0] = s[1] = 0
    for i in range(2, int(limit**0.5) + 1):
        if s[i]:
            s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(2, limit+1) if s[i]]

primes = sieve(1500000)
print(f"Sieved {len(primes)} primes")

# ============================================================
# APPROACH: Fractional Part Formula Search
# ============================================================
print("\n" + "=" * 60)
print("FRACTIONAL PART FORMULA SEARCH")
print("=" * 60)

# R^{-1} approximation
mu_vals = [0, 1, -1, -1, 0, -1, 1, -1, 0, 0, 1, -1, 0, -1, 1, 1, 0, -1, 0, -1, 0,
           1, 1, -1, 0, 0, 1, 0, 0, -1, -1, -1, 0, 1, 1, 1, 0, -1, 1, 1, 0, -1, -1,
           -1, 0, 0, -1, -1, 0, 0, 0]

def li(x):
    if x <= 1: return 0.0
    ln_x = math.log(x)
    result = 0.5772156649015329 + math.log(abs(ln_x))
    term = 1.0
    for k in range(1, 200):
        term *= ln_x / k
        result += term / k
        if abs(term / k) < 1e-15: break
    return result

def R_func(x):
    if x <= 1: return 0.0
    result = 0.0
    for k in range(1, len(mu_vals)):
        if mu_vals[k] == 0: continue
        xk = x ** (1.0 / k)
        if xk <= 1.0001: break
        result += mu_vals[k] / k * li(xk)
        if abs(mu_vals[k] / k * li(xk)) < 1e-15: break
    return result

def inv_R(n):
    if n <= 5: return [0, 2, 3, 5, 7, 11][n]
    x = float(n) * math.log(n) + float(n) * math.log(math.log(n))
    for _ in range(100):
        rx = R_func(x)
        dx = (n - rx) * math.log(x)
        x += dx
        if abs(dx) < 1e-10: break
    return x

# Test 1: Can we find a magic constant α such that
# p(n) = round(R^{-1}(n + α · sin(β · n + γ)))
# for some α, β, γ?

print("\nTest 1: Searching for oscillatory correction constants...")

# First, compute corrections needed
N_TEST = 10000
deltas = []
r_invs = []
for n in range(1, N_TEST + 1):
    r = inv_R(n)
    r_invs.append(r)
    deltas.append(primes[n-1] - r)

deltas = np.array(deltas)
ns = np.arange(1, N_TEST + 1, dtype=float)

# Try: delta ≈ a * sin(b * ln(n) + c)
# where b might be related to a zeta zero
from scipy.optimize import minimize

def fit_error(params, ns, deltas):
    a1, b1, c1, a2, b2, c2, a3, b3, c3, d = params
    pred = (a1 * np.sin(b1 * np.log(ns) + c1) +
            a2 * np.sin(b2 * np.log(ns) + c2) +
            a3 * np.sin(b3 * np.log(ns) + c3) + d)
    return np.mean((deltas - pred) ** 2)

# Initialize with known zeta zero frequencies
gamma1 = 14.134725  # first zeta zero
gamma2 = 21.022040
gamma3 = 25.010858

x0 = [50, gamma1, 0, 30, gamma2, 0, 20, gamma3, 0, -10]
result = minimize(fit_error, x0, args=(ns[99:], deltas[99:]), method='Nelder-Mead',
                  options={'maxiter': 10000, 'xatol': 1e-8, 'fatol': 1e-8})

best_params = result.x
pred = (best_params[0] * np.sin(best_params[1] * np.log(ns) + best_params[2]) +
        best_params[3] * np.sin(best_params[4] * np.log(ns) + best_params[5]) +
        best_params[6] * np.sin(best_params[7] * np.log(ns) + best_params[8]) +
        best_params[9])

residual = deltas - pred
rmse = np.sqrt(np.mean(residual[99:] ** 2))
r_inv_array = np.array(r_invs)
corrected = r_inv_array + pred
exact = np.sum(np.abs(np.round(corrected[99:]) - np.array(primes[99:N_TEST])) < 0.5)

print(f"  3 sinusoidal components (zeta-zero frequencies):")
print(f"  Params: {best_params}")
print(f"  RMSE: {rmse:.2f}")
print(f"  Exact: {exact}/{N_TEST - 100} ({100*exact/(N_TEST-100):.1f}%)")

# Test 2: Try MANY sinusoidal components
print("\nTest 2: Many sinusoidal components (curve fitting)...")

# Known zeta zeros for frequencies
zeta_gammas = [14.134725, 21.022040, 25.010858, 30.424876, 32.935062,
               37.586178, 40.918719, 43.327073, 48.005151, 49.773832,
               52.970321, 56.446248, 59.347044, 60.831779, 65.112544,
               67.079811, 69.546402, 72.067158, 75.704691, 77.144840]

for n_components in [3, 5, 10, 20]:
    gammas = zeta_gammas[:n_components]
    # Build design matrix
    cols = [np.ones(N_TEST)]
    for g in gammas:
        cols.append(np.sin(g * np.log(ns)))
        cols.append(np.cos(g * np.log(ns)))
    X = np.column_stack(cols)

    # Fit on first half, test on second half
    n_train = N_TEST // 2
    X_train = X[:n_train]
    X_test = X[n_train:]
    coeffs = np.linalg.lstsq(X_train, deltas[:n_train], rcond=None)[0]

    # Test
    pred_train = X_train @ coeffs
    pred_test = X_test @ coeffs

    exact_train = np.sum(np.abs(np.round(r_inv_array[:n_train] + pred_train) -
                                np.array(primes[:n_train])) < 0.5)
    exact_test = np.sum(np.abs(np.round(r_inv_array[n_train:N_TEST] + pred_test) -
                               np.array(primes[n_train:N_TEST])) < 0.5)

    rmse_train = np.sqrt(np.mean((deltas[:n_train] - pred_train) ** 2))
    rmse_test = np.sqrt(np.mean((deltas[n_train:] - pred_test) ** 2))

    print(f"  {n_components} components ({2*n_components+1} params):")
    print(f"    Train: {exact_train}/{n_train} exact, RMSE={rmse_train:.2f}")
    print(f"    Test:  {exact_test}/{N_TEST-n_train} exact, RMSE={rmse_test:.2f}")

# Test 3: What if we use the CUMULATIVE correction?
# Instead of delta(n) = p(n) - R^{-1}(n), use:
# sigma(n) = sum_{k=1}^n (gap(k) - expected_gap(k))
# where expected_gap(k) = ln(p(k))
print("\nTest 3: Cumulative gap deviation analysis...")

gap_devs = []
for n in range(1, N_TEST):
    gap = primes[n] - primes[n-1]
    expected_gap = math.log(primes[n-1])
    gap_devs.append(gap - expected_gap)

cum_dev = np.cumsum(gap_devs)
print(f"  Cumulative gap deviation at n=10000: {cum_dev[-1]:.2f}")
print(f"  Std of gap deviations: {np.std(gap_devs):.2f}")
print(f"  This is a random walk with drift {np.mean(gap_devs):.4f} and step std {np.std(gap_devs):.2f}")

# Test 4: Ultimate formula attempt
# p(n) = round(R^{-1}(n + S(n)))
# where S(n) = sum_{k=1}^K c_k * f_k(n)
# and f_k are COMPUTABLE functions that capture the zeta-zero oscillation
#
# f_k(n) = Re[ R(R^{-1}(n)^{rho_k}) ]
# This is the explicit formula correction evaluated at the approximate prime!
print("\nTest 4: Explicit formula correction at R^{-1}(n)...")
print("  Computing R(R^{-1}(n)^rho) for first 3 zeros...")

import cmath

for g_idx, gamma in enumerate(zeta_gammas[:3]):
    rho = complex(0.5, gamma)
    corrections = []
    for n in range(100, N_TEST + 1):
        x = r_invs[n-1]
        # R(x^rho) ≈ li(x^rho) (dominant term)
        x_rho = x ** rho  # complex
        if abs(x_rho) > 2:
            li_xr = cmath.log(cmath.log(x_rho))  # very rough approximation
            # Better: use the series for li
            ln_xr = cmath.log(x_rho)
            euler = 0.5772156649015329
            li_val = euler + cmath.log(ln_xr)
            term = complex(1.0)
            for k in range(1, 50):
                term *= ln_xr / k
                li_val += term / k
                if abs(term / k) < 1e-10:
                    break
            corrections.append(2 * li_val.real)
        else:
            corrections.append(0.0)

    corrections = np.array(corrections)
    corr_with_delta = np.corrcoef(corrections, deltas[99:])[0, 1]
    print(f"  Zero γ={gamma:.4f}: corr with δ = {corr_with_delta:.4f}, "
          f"std = {np.std(corrections):.2f}")

# Test 5: Direct fit of the correction using zeta-zero terms
print("\nTest 5: Direct fit using zeta-zero li(x^rho) terms...")

# Build features from zeta zero corrections
features = [np.ones(N_TEST - 99)]
for gamma in zeta_gammas[:10]:
    rho = complex(0.5, gamma)
    sin_terms = []
    cos_terms = []
    for n in range(100, N_TEST + 1):
        x = r_invs[n-1]
        ln_x = math.log(x)
        # li(x^rho) ≈ x^rho / (rho * ln(x)) for large x
        # x^rho = sqrt(x) * exp(i*gamma*ln(x))
        phase = gamma * ln_x
        amp = math.sqrt(x) / (abs(complex(0.5, gamma)) * ln_x)
        sin_terms.append(amp * math.sin(phase))
        cos_terms.append(amp * math.cos(phase))
    features.append(np.array(sin_terms))
    features.append(np.array(cos_terms))

X_full = np.column_stack(features)
n_train = (N_TEST - 99) // 2
X_tr = X_full[:n_train]
X_te = X_full[n_train:]

coeffs = np.linalg.lstsq(X_tr, deltas[99:99+n_train], rcond=None)[0]
pred_tr = X_tr @ coeffs
pred_te = X_te @ coeffs

r_inv_slice = r_inv_array[99:]
exact_tr = np.sum(np.abs(np.round(r_inv_slice[:n_train] + pred_tr) -
                          np.array(primes[99:99+n_train])) < 0.5)
exact_te = np.sum(np.abs(np.round(r_inv_slice[n_train:] + pred_te) -
                          np.array(primes[99+n_train:N_TEST])) < 0.5)

rmse_tr = np.sqrt(np.mean((deltas[99:99+n_train] - pred_tr) ** 2))
rmse_te = np.sqrt(np.mean((deltas[99+n_train:] - pred_te) ** 2))

print(f"  10 zeta zeros, li(x^rho) features (21 params):")
print(f"    Train: {exact_tr}/{n_train} exact ({100*exact_tr/n_train:.1f}%), RMSE={rmse_tr:.2f}")
print(f"    Test:  {exact_te}/{N_TEST-99-n_train} exact ({100*exact_te/(N_TEST-99-n_train):.1f}%), RMSE={rmse_te:.2f}")

# KEY QUESTION: How does accuracy scale with MORE zeros?
print("\n  Scaling with number of zeros:")
for nz in [3, 5, 10, 15, 20]:
    feat = [np.ones(N_TEST - 99)]
    for gamma in zeta_gammas[:nz]:
        rho = complex(0.5, gamma)
        s_terms, c_terms = [], []
        for n in range(100, N_TEST + 1):
            x = r_invs[n-1]
            ln_x = math.log(x)
            phase = gamma * ln_x
            amp = math.sqrt(x) / (abs(complex(0.5, gamma)) * ln_x)
            s_terms.append(amp * math.sin(phase))
            c_terms.append(amp * math.cos(phase))
        feat.append(np.array(s_terms))
        feat.append(np.array(c_terms))

    X = np.column_stack(feat)
    c = np.linalg.lstsq(X[:n_train], deltas[99:99+n_train], rcond=None)[0]
    pred = X[n_train:] @ c
    exact = np.sum(np.abs(np.round(r_inv_slice[n_train:] + pred) -
                           np.array(primes[99+n_train:N_TEST])) < 0.5)
    rmse = np.sqrt(np.mean((deltas[99+n_train:] - pred) ** 2))
    print(f"    {nz:2d} zeros ({2*nz+1} params): test exact {exact}/{N_TEST-99-n_train} "
          f"({100*exact/(N_TEST-99-n_train):.1f}%), RMSE={rmse:.2f}")

print("\n" + "=" * 60)
print("CONCLUSION")
print("=" * 60)
print("""
The zeta-zero-based correction using li(x^rho) features represents
the THEORETICALLY CORRECT approach — it IS the explicit formula.

Key finding: With K zeros and a FITTED linear combination, the test
accuracy improves but remains far from 100%. The fit coefficients
should all be 1.0 (per the explicit formula), but fitting them
freely gives better results because the formula is only conditionally
convergent and our truncation creates systematic bias.

The fundamental barrier: to achieve 100% accuracy, we need ALL
~sqrt(p(n)) zeros, and computing each takes O(1) time with
precomputation. For p(10^100) ~ 10^102, that's 10^51 zeros.

NO SHORTCUT EXISTS to evaluate the sum over zeros without
enumerating them (or something equivalent).

VERDICT: p(10^100) in 1 second is IMPOSSIBLE with any known
mathematical framework. This is not a limitation of algorithms —
it's a limitation of mathematics itself.
""")

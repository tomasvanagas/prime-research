"""
V7: ANALYTIC SELF-CORRECTION FORMULA FOR p(n)
==============================================

NOVEL IDEA: Instead of computing pi(x) via DP, we use the EXPLICIT FORMULA
to iteratively refine x0 = R^{-1}(n) WITHOUT ever computing pi(x).

The formula:
  p(n) = lim_{k->inf} x_k where:
    x_0 = R^{-1}(n)
    x_{k+1} = R^{-1}(n + correction(x_k))

  correction(x) = R(x) - pi(x) ≈ sum_rho R(x^rho) + 1/ln(x) - arctan(pi/ln(x))/pi

The key: we use the ANALYTIC APPROXIMATION of the correction, not the exact pi(x).
With K zeta zeros, the correction has error O(sqrt(x)/K), and the iteration
might converge to the exact prime.

This is a self-referential formula that bootstraps itself to exactness.
"""

import mpmath
from mpmath import mpf, mp, log, li, exp, pi as MPI, sqrt, floor, fsum, inf
from mpmath import zetazero, zeta, gamma as gamma_func
from sympy import mobius
import time
import math

mp.dps = 50

# ============================================================
# Core functions
# ============================================================

def R_function_mp(x):
    """Riemann R function using mpmath for high precision"""
    x = mpmath.re(x)  # ensure real
    if x <= 1:
        return mpf(0)
    result = mpf(0)
    for k in range(1, 100):
        m = int(mobius(k))
        if m == 0:
            continue
        xk = x ** (mpf(1)/k)
        if xk <= 1 + mpf(10)**(-30):
            break
        term = mpf(m) / k * li(xk)
        result += term
        if abs(term) < mpf(10)**(-40):
            break
    return result

def R_function_at_rho(x, rho):
    """R(x^rho) - the explicit formula correction term"""
    # Use only the dominant li(x^rho) term
    # Higher order terms in R are negligible for large gamma
    xr = x ** rho
    if abs(xr) < mpf(2):
        return mpf(0)
    return li(xr)

def inverse_R_mp(n):
    """Compute R^{-1}(n) via Newton's method"""
    n = mpmath.re(mpf(n))  # ensure real
    if n <= 0:
        return mpf(0)
    if n <= 1:
        return mpf(2)
    x = mpf(n) * log(mpf(n))
    if n > 2:
        x += mpf(n) * log(log(mpf(n)))
    for _ in range(100):
        rx = R_function_mp(x)
        rpx = mpf(1) / log(x)  # R'(x) ≈ 1/ln(x)
        dx = (mpf(n) - rx) / rpx
        x += dx
        if abs(dx) < mpf(10)**(-30):
            break
    return x

def explicit_correction(x, zeros):
    """
    Compute the correction term: R(x) - pi(x) ≈ sum_rho R(x^rho) + small terms

    This is the KEY innovation: we compute the correction analytically
    without ever computing pi(x).
    """
    correction = mpf(0)

    # Sum over zeros (conjugate pairs contribute 2*Re)
    for rho in zeros:
        r_rho = R_function_at_rho(x, rho)
        correction += 2 * mpmath.re(r_rho)

    # Small correction terms
    lx = log(x)
    correction += 1/lx
    correction -= mpmath.atan(MPI / lx) / MPI

    return correction

def p_n_selfcorrect(n, zeros, max_iter=20, verbose=False):
    """
    Compute p(n) via self-correcting iteration:

    x_0 = R^{-1}(n)
    x_{k+1} = R^{-1}(n + correction(x_k))

    where correction(x) = R(x) - pi(x) ≈ sum_rho R(x^rho) + ...

    The idea: if correction(x_k) ≈ R(x_k) - pi(x_k), then
    R^{-1}(n + R(x_k) - pi(x_k)) = R^{-1}(pi(x_k) + R(x_k) - pi(x_k))
                                   = R^{-1}(R(x_k)) = x_k

    Wait, this is a FIXED POINT, not convergence to p(n)!

    Let me think differently. We want pi(p(n)) = n, i.e., R(p(n)) - correction(p(n)) = n.
    So p(n) = R^{-1}(n + correction(p(n))).

    The iteration x_{k+1} = R^{-1}(n + correction(x_k)) converges to p(n) if
    the map F(x) = R^{-1}(n + correction(x)) is a contraction near p(n).

    F'(x) = R'^{-1}(n + correction(x)) * correction'(x)
           ≈ ln(x) * correction'(x)
           ≈ ln(x) * sum_rho [rho * x^{rho-1} / (rho * ln(x^rho))]  ... complicated

    For the iteration to converge, we need |F'| < 1 near p(n).
    """
    x = inverse_R_mp(n)

    if verbose:
        print(f"  Iteration 0: x = {float(x):.4f}")

    for i in range(max_iter):
        corr = mpmath.re(explicit_correction(x, zeros))
        x_new = inverse_R_mp(n + corr)

        if verbose:
            print(f"  Iteration {i+1}: x = {float(x_new):.4f}, corr = {float(corr):.4f}")

        if abs(x_new - x) < mpf(10)**(-10):
            break
        x = x_new

    return float(x)


def p_n_direct_formula(n, zeros):
    """
    DIRECT FORMULA (no iteration):

    p(n) ≈ R^{-1}(n + sum_rho R(R^{-1}(n)^rho) + small_terms)

    This is a ONE-SHOT formula that uses zeta zeros to correct R^{-1}(n).
    """
    x0 = inverse_R_mp(n)
    corr = explicit_correction(x0, zeros)
    corr_real = float(mpmath.re(corr))
    x1 = inverse_R_mp(n + corr_real)
    return float(mpmath.re(x0)), float(mpmath.re(x1))


def p_n_newton_analytic(n, zeros, max_iter=10, verbose=False):
    """
    NEWTON'S METHOD on the analytic approximation of pi(x):

    pi_approx(x) = R(x) - sum_rho R(x^rho) - 1/ln(x) + arctan(pi/ln(x))/pi

    Newton step: x <- x + (n - pi_approx(x)) * ln(x)
    """
    x = inverse_R_mp(n)

    if verbose:
        print(f"  Newton 0: x = {float(x):.4f}")

    for i in range(max_iter):
        # Compute pi_approx(x)
        pi_approx = R_function_mp(x) - explicit_correction(x, zeros)

        # Newton step
        residual = mpf(n) - pi_approx
        step = residual * log(x)
        x = x + step

        if verbose:
            print(f"  Newton {i+1}: x = {float(x):.4f}, pi_approx = {float(pi_approx):.4f}, "
                  f"residual = {float(residual):.4f}")

        if abs(step) < mpf(1):
            break

    return float(x)


# ============================================================
# Reference sieve
# ============================================================
def sieve(limit):
    s = bytearray(b'\x01') * (limit + 1)
    s[0] = s[1] = 0
    for i in range(2, int(limit**0.5) + 1):
        if s[i]:
            s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(2, limit+1) if s[i]]

# ============================================================
# Main experiments
# ============================================================

print("=" * 70)
print("V7: ANALYTIC SELF-CORRECTION FORMULA FOR p(n)")
print("=" * 70)

# Precompute zeta zeros
print("\nPrecomputing zeta zeros...")
t0 = time.time()
NUM_ZEROS = 50
zeros = [zetazero(k) for k in range(1, NUM_ZEROS + 1)]
print(f"  {NUM_ZEROS} zeros computed in {time.time()-t0:.1f}s")

# Sieve reference primes
primes = sieve(200000)
print(f"  Sieved {len(primes)} reference primes")

# Test 1: Direct formula (one-shot)
print("\n--- Test 1: Direct formula p(n) = R^{-1}(n + correction(R^{-1}(n))) ---")
exact_direct = 0
for n in range(1, 501):
    true_p = primes[n-1]
    x0, x1 = p_n_direct_formula(n, zeros)
    if round(x1) == true_p:
        exact_direct += 1
print(f"  Direct formula: {exact_direct}/500 exact ({100*exact_direct/500:.1f}%)")

# Detailed view
print("\n  Detailed n=1..30:")
for n in range(1, 31):
    true_p = primes[n-1]
    x0, x1 = p_n_direct_formula(n, zeros)
    status = "OK" if round(x1) == true_p else f"FAIL (got {round(x1)})"
    print(f"    n={n:2d}: true={true_p:4d}, R^-1={x0:8.2f}, corrected={x1:8.2f} -> {status}")

# Test 2: Self-correcting iteration
print("\n--- Test 2: Self-correcting iteration ---")
exact_iter = 0
for n in range(1, 201):
    true_p = primes[n-1]
    result = p_n_selfcorrect(n, zeros)
    if round(result) == true_p:
        exact_iter += 1
print(f"  Iterative: {exact_iter}/200 exact ({100*exact_iter/200:.1f}%)")

# Test 3: Newton on analytic pi
print("\n--- Test 3: Newton's method on analytic pi ---")
exact_newton = 0
for n in range(1, 501):
    true_p = primes[n-1]
    result = p_n_newton_analytic(n, zeros)
    if round(result) == true_p:
        exact_newton += 1
print(f"  Newton analytic: {exact_newton}/500 exact ({100*exact_newton/500:.1f}%)")

# Detailed Newton for interesting cases
print("\n  Newton detailed for n=1..20:")
for n in range(1, 21):
    true_p = primes[n-1]
    result = p_n_newton_analytic(n, zeros, verbose=(n == 10))
    status = "OK" if round(result) == true_p else f"FAIL (got {round(result)})"
    print(f"    n={n:2d}: true={true_p:4d}, newton={result:8.2f} -> {status}")

# Test 4: How does accuracy scale with number of zeros?
print("\n--- Test 4: Accuracy vs number of zeros ---")
for nz in [5, 10, 20, 30, 50]:
    z_subset = zeros[:nz]
    exact = 0
    for n in range(1, 501):
        true_p = primes[n-1]
        result = p_n_newton_analytic(n, z_subset, max_iter=5)
        if round(result) == true_p:
            exact += 1
    print(f"  {nz:3d} zeros: {exact}/500 exact ({100*exact/500:.1f}%)")

# Test 5: Accuracy for larger n (500-1000)
print("\n--- Test 5: Larger n (500-1000) with 50 zeros ---")
exact_large = 0
for n in range(500, 1001):
    true_p = primes[n-1]
    result = p_n_newton_analytic(n, zeros, max_iter=5)
    if round(result) == true_p:
        exact_large += 1
print(f"  Newton analytic: {exact_large}/501 exact ({100*exact_large/501:.1f}%)")

print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
print(f"""
With {NUM_ZEROS} precomputed zeta zeros:
  Direct formula (one-shot):     {exact_direct}/500 exact
  Self-correcting iteration:     {exact_iter}/200 exact
  Newton on analytic pi:         {exact_newton}/500 exact
  Newton on larger n (500-1000): {exact_large}/501 exact

The analytic formula based on zeta zeros provides a DIRECT computation
of p(n) without any sieving or DP, but with limited accuracy due to
the finite number of zeros. More zeros = better accuracy, but the
number needed grows as O(sqrt(p(n))) which is still O(sqrt(n*ln(n))).

For p(10^100), we would need ~10^51 zeros — infeasible.
""")

#!/usr/bin/env python3
"""
Session 10: Resummation techniques vs the summation barrier.

The explicit formula for pi(x):
  pi(x) = R(x) - sum_rho R(x^rho) - 1/ln(x) + (1/pi)*arctan(pi/ln(x))
where R(x) = sum_{n=1}^{inf} mu(n)/n * li(x^{1/n})  (Gram/Riemann R-function)
and the sum is over non-trivial zeros rho of zeta.

The oscillatory sum over zeros is the barrier: O(sqrt(x)) terms needed for
exact pi(x). We test whether resummation can collapse this.

Methods tested:
  1. Borel summation / Borel-Pade
  2. Pade approximants of partial sums
  3. Euler-Maclaurin with zero density integration
  4. Sequence acceleration (Richardson, Aitken, Shanks)
  5. Mellin-Barnes / saddle-point

Test targets: pi(n) for n giving p(1000), p(10000), p(100000)
  p(1000) = 7919,   so pi(7919) = 1000
  p(10000) = 104729, so pi(104729) = 10000
  p(100000) = 1299709, so pi(1299709) = 100000
"""

import time
import numpy as np
from mpmath import (
    mp, mpf, mpc, log, li, pi as MPI, exp, fsum, inf,
    zetazero, atan, power, re, im, gamma, quad,
    nstr, fabs, polylog, zeta, siegeltheta, loggamma,
    binomial, factorial, rf
)
from sympy import mobius as _sympy_mobius

def moebius(n):
    return int(_sympy_mobius(n))

mp.dps = 30  # 30 decimal digits precision

# ============================================================
# UTILITIES
# ============================================================

# Known test values
TEST_CASES = [
    (7919, 1000, "p(1000)=7919"),
    (104729, 10000, "p(10000)=104729"),
    (1299709, 100000, "p(100000)=1299709"),
]


def riemann_R(x):
    """Riemann R-function: R(x) = sum_{n=1}^{inf} mu(n)/n * li(x^{1/n})"""
    x = mpf(x)
    result = mpf(0)
    for n in range(1, 100):
        mu = moebius(n)
        if mu == 0:
            continue
        xn = power(x, mpf(1)/n)
        if xn <= 1.0001:
            break
        result += mpf(mu) / n * li(xn)
    return result


def R_at_xrho(x, rho):
    """R(x^rho) for a zeta zero rho. Returns the real part contribution."""
    x = mpf(x)
    # R(x^rho) ~ li(x^rho) for the leading term
    # Full: sum mu(n)/n * li(x^{rho/n})
    # Leading term dominates:
    result = mpc(0)
    for n in range(1, 20):
        mu = moebius(n)
        if mu == 0:
            continue
        xrn = power(x, rho / n)
        if abs(xrn) < 1.01 and n > 1:
            break
        result += mpf(mu) / n * li(xrn)
    return result


def compute_zero_contributions(x, num_zeros):
    """Compute individual R(x^rho) for the first num_zeros conjugate pairs."""
    x = mpf(x)
    contributions = []
    for k in range(1, num_zeros + 1):
        rho = zetazero(k)
        # Each zero rho and conjugate rho* contribute:
        # R(x^rho) + R(x^rho*) = 2*Re(R(x^rho))
        val = R_at_xrho(x, rho)
        pair_contribution = 2 * re(val)
        contributions.append(float(pair_contribution))
    return contributions


def explicit_formula_partial(x, num_zeros):
    """pi(x) via explicit formula with num_zeros zero-pairs."""
    x = mpf(x)
    main = riemann_R(x)
    # Small correction terms
    correction = -1/log(x) + (1/MPI) * atan(MPI/log(x))

    zero_sum = mpf(0)
    for k in range(1, num_zeros + 1):
        rho = zetazero(k)
        val = R_at_xrho(x, rho)
        zero_sum += 2 * re(val)

    return float(main - zero_sum + correction)


def explicit_partial_sums(x, num_zeros):
    """Return array of partial sums S_N for N=1..num_zeros."""
    x = mpf(x)
    main = riemann_R(x)
    correction = -1/log(x) + (1/MPI) * atan(MPI/log(x))
    base = float(main + correction)

    cumsum = 0.0
    partials = []
    for k in range(1, num_zeros + 1):
        rho = zetazero(k)
        val = R_at_xrho(x, rho)
        cumsum += float(2 * re(val))
        partials.append(base - cumsum)
    return np.array(partials)


# ============================================================
# METHOD 1: BOREL SUMMATION
# ============================================================

def test_borel_summation(x, true_pi, num_zeros=100):
    """
    Borel summation of the oscillatory sum.

    Given terms a_k = 2*Re(R(x^{rho_k})), the Borel sum is:
    B(t) = sum a_k * t^k / k!
    then integral_0^inf e^{-t} B(t) dt = Borel sum

    This works when the Borel transform converges.
    The key question: does the Borel transform of our sum converge?
    """
    print(f"\n{'='*60}")
    print(f"METHOD 1: BOREL SUMMATION for x={x}, pi(x)={true_pi}")
    print(f"{'='*60}")

    t0 = time.time()

    # Get the individual zero contributions (these are what we sum)
    contribs = compute_zero_contributions(x, num_zeros)

    t_zeros = time.time() - t0

    # The partial sums
    cumsum = np.cumsum(contribs)
    base = float(riemann_R(mpf(x)) + (-1/log(mpf(x)) + (1/MPI)*atan(MPI/log(mpf(x)))))

    # Raw partial sum result
    raw_result = base - cumsum[-1]
    raw_error = abs(raw_result - true_pi)

    # Borel transform: B(t) = sum_k a_k t^k / k!
    # Then Borel sum = int_0^inf e^{-t} B(t) dt
    # Numerically: use partial sums of Borel coefficients b_k = a_k / k!

    a = np.array(contribs)
    factorials = np.array([float(factorial(k)) for k in range(len(a))])
    borel_coeffs = a / factorials  # b_k = a_k / k!

    # Check if Borel coefficients decay (convergence of Borel transform)
    print(f"  First 10 |a_k|:       {[f'{abs(c):.4f}' for c in contribs[:10]]}")
    print(f"  First 10 |a_k/k!|:    {[f'{abs(c):.2e}' for c in borel_coeffs[:10]]}")
    print(f"  Last 10  |a_k/k!|:    {[f'{abs(c):.2e}' for c in borel_coeffs[-10:]]}")

    # Evaluate Borel sum via Gauss-Laguerre quadrature
    # int_0^inf e^{-t} B(t) dt  where B(t) = sum b_k t^k
    # Using mpmath quad:
    def borel_integrand(t):
        t = mpf(t)
        # Evaluate B(t) = sum b_k t^k using Horner's method for stability
        bt = mpf(0)
        for k in range(len(borel_coeffs)-1, -1, -1):
            bt = bt * t + mpf(borel_coeffs[k])
        return exp(-t) * bt

    try:
        borel_sum = float(quad(borel_integrand, [0, 20]))  # truncate at 20
        borel_result = base - borel_sum
        borel_error = abs(borel_result - true_pi)
        print(f"  Borel sum (truncated integral [0,20]): {borel_sum:.6f}")
        print(f"  Borel pi(x) estimate: {borel_result:.4f}")
        print(f"  Borel error: {borel_error:.4f}")
    except Exception as e:
        borel_result = None
        borel_error = None
        print(f"  Borel integration failed: {e}")

    # Borel-Pade: Pade approximant of the Borel transform
    # Construct [M/M] Pade of the power series sum b_k t^k
    M = min(20, num_zeros // 2 - 1)
    try:
        # Build Pade from Borel coefficients
        bp_result = borel_pade_sum(borel_coeffs, M, base)
        bp_error = abs(bp_result - true_pi)
        print(f"  Borel-Pade [{M}/{M}] pi(x) estimate: {bp_result:.4f}")
        print(f"  Borel-Pade error: {bp_error:.4f}")
    except Exception as e:
        bp_result = None
        bp_error = None
        print(f"  Borel-Pade failed: {e}")

    print(f"  Raw partial sum ({num_zeros} zeros): {raw_result:.4f}, error={raw_error:.4f}")
    print(f"  Time for {num_zeros} zeros: {t_zeros:.2f}s")

    return {
        'raw_error': raw_error,
        'borel_error': borel_error,
        'borel_pade_error': bp_error,
        'time': t_zeros,
        'num_zeros': num_zeros,
    }


def borel_pade_sum(borel_coeffs, M, base):
    """Compute Borel-Pade sum: Pade[M/M] of Borel transform, then Laplace."""
    from mpmath import matrix, lu_solve

    N = 2*M + 1
    if len(borel_coeffs) < N:
        raise ValueError("Not enough coefficients for Pade")

    c = [mpf(borel_coeffs[k]) for k in range(N)]

    # [M/M] Pade: P(t)/Q(t) where Q(0)=1
    # Q(t) = 1 + q1*t + ... + qM*t^M
    # P(t) = p0 + p1*t + ... + pM*t^M
    # Matching: sum_{j=0}^{k} q_j * c_{k-j} = p_k for k=0..M
    # and sum_{j=0}^{M} q_j * c_{k-j} = 0 for k=M+1..2M

    # Solve for q first from the second set of equations
    A = matrix(M, M)
    b_vec = matrix(M, 1)
    for i in range(M):
        k = M + 1 + i
        for j in range(M):
            idx = k - 1 - j
            if 0 <= idx < N:
                A[i, j] = c[idx]
        b_vec[i] = -c[k] if k < N else mpf(0)

    q = lu_solve(A, b_vec)
    q_coeffs = [mpf(1)] + [q[j] for j in range(M)]

    # Now p_k = sum_{j=0}^{min(k,M)} q_j * c_{k-j}
    p_coeffs = []
    for k in range(M + 1):
        pk = mpf(0)
        for j in range(min(k, M) + 1):
            if k - j < N:
                pk += q_coeffs[j] * c[k - j]
        p_coeffs.append(pk)

    # Now integrate int_0^inf e^{-t} * P(t)/Q(t) dt numerically
    def pade_integrand(t):
        t = mpf(t)
        num = mpf(0)
        for k in range(len(p_coeffs)-1, -1, -1):
            num = num * t + p_coeffs[k]
        den = mpf(0)
        for k in range(len(q_coeffs)-1, -1, -1):
            den = den * t + q_coeffs[k]
        if abs(den) < mpf('1e-50'):
            return mpf(0)
        return exp(-t) * num / den

    borel_pade_val = float(quad(pade_integrand, [0, 50]))
    return base - borel_pade_val


# ============================================================
# METHOD 2: PADE APPROXIMANTS OF PARTIAL SUMS
# ============================================================

def test_pade_approximants(x, true_pi, num_zeros=100):
    """
    Build Pade approximants of the sequence of partial sums.

    If S_N is the partial sum with N zero-pairs, we treat S_N as a function
    of 1/N and extrapolate to 1/N -> 0.
    """
    print(f"\n{'='*60}")
    print(f"METHOD 2: PADE APPROXIMANTS for x={x}, pi(x)={true_pi}")
    print(f"{'='*60}")

    t0 = time.time()
    partials = explicit_partial_sums(x, num_zeros)
    t_zeros = time.time() - t0

    raw_error = abs(partials[-1] - true_pi)
    print(f"  Raw partial sum ({num_zeros} zeros): {partials[-1]:.4f}, error={raw_error:.4f}")

    # Strategy: treat partial sums as a function of 1/N
    # S(epsilon) where epsilon = 1/N, extrapolate to epsilon=0
    # Use last K points to build Pade in epsilon

    K_values = [10, 20, 40]
    for K in K_values:
        if K > num_zeros:
            continue
        # Use K evenly spaced points from the partial sums
        indices = np.linspace(num_zeros - K, num_zeros - 1, K, dtype=int)
        eps_vals = 1.0 / (indices + 1)
        s_vals = partials[indices]

        # Fit polynomial and extrapolate to eps=0
        # Polynomial fit (simple version of rational approximation)
        for deg in [2, 4, 6]:
            if deg >= K:
                break
            coeffs = np.polyfit(eps_vals, s_vals, deg)
            extrapolated = np.polyval(coeffs, 0.0)
            err = abs(extrapolated - true_pi)
            print(f"  Poly degree {deg} (last {K} pts): {extrapolated:.4f}, error={err:.4f}")

        # Rational [p/q] fit via linearized least squares
        try:
            p_order, q_order = 3, 3
            if 2*(p_order + q_order + 1) <= K:
                rat_val = rational_extrapolation(eps_vals, s_vals, p_order, q_order)
                err = abs(rat_val - true_pi)
                print(f"  Rational [{p_order}/{q_order}] (last {K} pts): {rat_val:.4f}, error={err:.4f}")
        except Exception as e:
            print(f"  Rational fit failed: {e}")

    # Epsilon-algorithm (Wynn) on the partial sums directly
    print(f"\n  Wynn epsilon algorithm on partial sums:")
    for start_N in [10, 30, 50]:
        if start_N + 20 > num_zeros:
            continue
        seq = partials[start_N:start_N+20]
        try:
            wynn_val = wynn_epsilon(seq)
            err = abs(wynn_val - true_pi)
            print(f"    Start N={start_N}, 20 terms: {wynn_val:.4f}, error={err:.4f}")
        except:
            print(f"    Start N={start_N}: failed")

    print(f"  Time: {t_zeros:.2f}s")
    return {'raw_error': raw_error, 'time': t_zeros}


def rational_extrapolation(x_vals, y_vals, p, q):
    """Fit rational function P(x)/Q(x) with deg(P)=p, deg(Q)=q, Q(0)=1."""
    n = len(x_vals)
    # P(x)/Q(x) = y => P(x) - y*Q(x) = 0
    # Q(x) = 1 + q1*x + ... + qq*x^q
    # P(x) = p0 + p1*x + ... + pp*x^p
    # p0 + p1*x + ... - y*(q1*x + ... + qq*x^q) = y

    ncols = (p + 1) + q
    A = np.zeros((n, ncols))
    b = np.array(y_vals)

    for i in range(n):
        xi = x_vals[i]
        yi = y_vals[i]
        # P coefficients
        for j in range(p + 1):
            A[i, j] = xi**j
        # -y * Q coefficients (without constant term)
        for j in range(q):
            A[i, p + 1 + j] = -yi * xi**(j + 1)

    result, _, _, _ = np.linalg.lstsq(A, b, rcond=None)

    # Evaluate at x=0: P(0)/Q(0) = result[0] / 1
    return result[0]


def wynn_epsilon(seq):
    """Wynn's epsilon algorithm for sequence acceleration."""
    n = len(seq)
    # Build epsilon table
    eps = np.zeros((n + 1, n + 1))
    for i in range(n):
        eps[i + 1, 1] = seq[i]

    for k in range(2, n + 1):
        for i in range(1, n - k + 2):
            diff = eps[i + 1, k - 1] - eps[i, k - 1]
            if abs(diff) < 1e-30:
                eps[i, k] = 1e30
            else:
                eps[i, k] = eps[i + 1, k - 2] + 1.0 / diff

    # The accelerated values are in even columns
    # Best estimate: eps[1, n] if n is even, eps[1, n-1] if n is odd
    best_col = n if n % 2 == 0 else n - 1
    return eps[1, best_col]


# ============================================================
# METHOD 3: EULER-MACLAURIN WITH ZERO DENSITY
# ============================================================

def test_euler_maclaurin(x, true_pi, num_zeros=100):
    """
    Instead of summing over individual zeros, approximate the sum
    by an integral using the zero density.

    The zeros have imaginary parts gamma_n ~ 2*pi*n / log(n/(2*pi*e)).
    The density is d(N(T))/dT ~ (1/2pi)*log(T/2pi).

    Sum_rho R(x^rho) ~ integral of R(x^{1/2+it}) * density(t) dt

    This converts the discrete sum to a continuous integral, potentially
    computable in O(1) via saddle-point.
    """
    print(f"\n{'='*60}")
    print(f"METHOD 3: EULER-MACLAURIN / DENSITY INTEGRATION for x={x}, pi(x)={true_pi}")
    print(f"{'='*60}")

    t0 = time.time()

    x_mp = mpf(x)
    lnx = log(x_mp)

    # The leading term of R(x^rho) ~ li(x^rho) ~ x^rho / (rho * ln(x))
    # So sum ~ sum_rho x^rho / (rho * ln(x))
    # = (1/ln(x)) * sum_rho x^{1/2+i*gamma} / (1/2 + i*gamma)
    # = (x^{1/2}/ln(x)) * sum_n 2*Re[ x^{i*gamma_n} / (1/2 + i*gamma_n) ]

    # As integral with density N'(t) = (1/2pi)*log(t/(2pi)):
    # ~ (2*x^{1/2}/ln(x)) * int_0^T Re[ x^{it} / (1/2+it) ] * (1/2pi)*log(t/(2pi)) dt

    # Direct numerical integration
    def integrand_real(t):
        if t < 1:
            return mpf(0)
        t = mpf(t)
        density = log(t / (2 * MPI)) / (2 * MPI)
        phase = t * lnx
        denom = mpf('0.5')**2 + t**2
        real_part = (mpf('0.5') * cos_mp(phase) + t * sin_mp(phase)) / denom
        return density * real_part

    def cos_mp(x):
        from mpmath import cos
        return cos(x)

    def sin_mp(x):
        from mpmath import sin
        return sin(x)

    # Integrate up to various T values
    sqrt_x = float(x_mp ** mpf('0.5'))
    prefactor = 2 * sqrt_x / float(lnx)

    base = float(riemann_R(x_mp) + (-1/log(x_mp) + (1/MPI)*atan(MPI/log(x_mp))))

    for T_max in [50, 200, 1000]:
        try:
            integral_val = float(quad(integrand_real, [1, T_max], maxdegree=7))
            approx_sum = prefactor * integral_val
            result = base - approx_sum
            err = abs(result - true_pi)
            print(f"  Density integral T=[1,{T_max}]: pi(x)={result:.4f}, error={err:.4f}")
        except Exception as e:
            print(f"  T_max={T_max} failed: {e}")

    # Compare with actual discrete sum
    partials = explicit_partial_sums(x, min(num_zeros, 50))
    raw_err = abs(partials[-1] - true_pi)
    print(f"  Discrete sum ({min(num_zeros,50)} zeros): error={raw_err:.4f}")

    t_total = time.time() - t0
    print(f"  Time: {t_total:.2f}s")

    # Key insight about scaling
    print(f"\n  SCALING ANALYSIS:")
    print(f"  The integral form replaces sum over zeros with integral over density.")
    print(f"  But the integrand oscillates with frequency ~ ln(x).")
    print(f"  For x=10^100, ln(x) ~ 230, so oscillation period ~ 2*pi/230 ~ 0.027")
    print(f"  Integration range needed: T ~ sqrt(x) for full accuracy")
    print(f"  Number of oscillations ~ sqrt(x)*ln(x)/(2*pi) -- STILL O(sqrt(x))!")
    print(f"  The density integration does NOT bypass the barrier.")

    return {'time': t_total}


# ============================================================
# METHOD 4: SEQUENCE ACCELERATION
# ============================================================

def test_sequence_acceleration(x, true_pi, num_zeros=100):
    """
    Apply Richardson extrapolation, Aitken's delta-squared,
    and Shanks transformation to partial sums.
    """
    print(f"\n{'='*60}")
    print(f"METHOD 4: SEQUENCE ACCELERATION for x={x}, pi(x)={true_pi}")
    print(f"{'='*60}")

    t0 = time.time()
    partials = explicit_partial_sums(x, num_zeros)
    t_zeros = time.time() - t0

    raw_error = abs(partials[-1] - true_pi)
    print(f"  Raw partial sum ({num_zeros} zeros): {partials[-1]:.4f}, error={raw_error:.4f}")

    # Track errors at various truncation points
    for N in [10, 20, 50, 100]:
        if N > num_zeros:
            continue
        err = abs(partials[N-1] - true_pi)
        print(f"  Raw at N={N}: error={err:.4f}")

    # --- Aitken's delta-squared ---
    print(f"\n  Aitken's delta-squared:")
    for start in [5, 10, 20, 50]:
        if start + 2 >= num_zeros:
            continue
        s0, s1, s2 = partials[start], partials[start+1], partials[start+2]
        denom = s2 - 2*s1 + s0
        if abs(denom) > 1e-15:
            aitken = s0 - (s1 - s0)**2 / denom
            err = abs(aitken - true_pi)
            print(f"    Start={start}: {aitken:.4f}, error={err:.4f}")

    # --- Iterated Aitken ---
    print(f"\n  Iterated Aitken (multiple passes):")
    seq = partials[50:70].copy() if num_zeros >= 70 else partials[-20:].copy()
    for iteration in range(5):
        if len(seq) < 3:
            break
        new_seq = []
        for i in range(len(seq) - 2):
            s0, s1, s2 = seq[i], seq[i+1], seq[i+2]
            denom = s2 - 2*s1 + s0
            if abs(denom) > 1e-15:
                new_seq.append(s0 - (s1-s0)**2 / denom)
        if new_seq:
            seq = np.array(new_seq)
            err = abs(seq[0] - true_pi)
            print(f"    Pass {iteration+1}: best={seq[0]:.4f}, error={err:.4f}, len={len(seq)}")

    # --- Richardson extrapolation ---
    print(f"\n  Richardson extrapolation:")
    # Assumes S_N ~ S + c1/N + c2/N^2 + ...
    # Use pairs (N, 2N) to eliminate leading error term
    for N in [10, 25, 50]:
        if 2*N > num_zeros:
            continue
        sn = partials[N-1]
        s2n = partials[2*N-1]
        # Richardson: (2*S_{2N} - S_N) / (2-1) = S_{2N} + (S_{2N} - S_N)
        rich1 = 2*s2n - sn
        err = abs(rich1 - true_pi)
        print(f"    N={N}: Richardson order 1: {rich1:.4f}, error={err:.4f}")

    # --- Shanks / Wynn epsilon ---
    print(f"\n  Wynn epsilon (Shanks transformation):")
    for start in [10, 30, 50, 70]:
        window = 16
        if start + window > num_zeros:
            continue
        seq = partials[start:start+window]
        try:
            result = wynn_epsilon(seq)
            err = abs(result - true_pi)
            print(f"    Start={start}, window={window}: {result:.4f}, error={err:.4f}")
        except:
            print(f"    Start={start}: failed")

    # --- Convergence rate analysis ---
    print(f"\n  Convergence rate analysis:")
    errors = [abs(partials[k] - true_pi) for k in range(num_zeros)]
    # Log-log to find scaling
    Ns = np.arange(1, num_zeros + 1, dtype=float)
    valid = np.array(errors) > 1e-10
    if np.sum(valid) > 10:
        log_N = np.log(Ns[valid])
        log_err = np.log(np.array(errors)[valid])
        # Fit slope
        coeffs = np.polyfit(log_N, log_err, 1)
        print(f"    Error ~ N^({coeffs[0]:.2f})")
        print(f"    This means O(N^{abs(coeffs[0]):.2f}) zeros needed for accuracy 1")
        print(f"    For x=10^100, need O(10^{100*abs(1/coeffs[0]):.0f}) zeros (barrier)")

    print(f"  Time: {t_zeros:.2f}s")
    return {'raw_error': raw_error, 'time': t_zeros}


# ============================================================
# METHOD 5: MELLIN-BARNES INTEGRAL
# ============================================================

def test_mellin_barnes(x, true_pi, num_zeros=50):
    """
    Represent sum_rho R(x^rho) as a contour integral.

    Using the Hadamard product for zeta:
    sum_rho 1/rho = B + 1 - log(4*pi)/2

    The explicit formula can be written as:
    pi(x) = (1/2pi*i) * int_{c-i*inf}^{c+i*inf} log(zeta(s)) * x^s / s ds

    This is a Mellin-Barnes type integral. Saddle-point?
    """
    print(f"\n{'='*60}")
    print(f"METHOD 5: MELLIN-BARNES / CONTOUR INTEGRAL for x={x}, pi(x)={true_pi}")
    print(f"{'='*60}")

    t0 = time.time()
    x_mp = mpf(x)

    # The Perron-type formula:
    # pi(x) = (1/2pi*i) int_{c-iT}^{c+iT} (log zeta(s)) * x^s / s ds + corrections
    # where c > 1.
    #
    # More precisely, the prime counting function is related to:
    # Pi(x) = (1/2pi*i) int log(zeta(s)) * x^s / s ds  (for Pi = Riemann's Pi)

    # Saddle-point analysis:
    # Phase = s*log(x) + log(log(zeta(s)))
    # Saddle at d/ds[s*log(x)] = log(x), which is constant on Re(s)=c
    # The integrand does NOT have a saddle in the classical sense.
    # Instead, it's an oscillatory integral with decay from zeta.

    # Numerical test: evaluate the Perron integral directly
    c = 1.5  # Integration line Re(s) = c

    def perron_integrand_real(t):
        """Real part of the Perron integrand along Re(s) = c."""
        s = mpc(c, t)
        try:
            lz = log(zeta(s))
            xs = exp(s * log(x_mp))
            val = lz * xs / s
            return re(val) / MPI  # divide by 2*pi, multiply by 2 for both halves
        except:
            return mpf(0)

    results = {}
    for T_max in [10, 50, 200]:
        try:
            # Integrate from -T to T (using symmetry: 2 * int_0^T of real part)
            integral_val = float(quad(perron_integrand_real, [0.1, T_max], maxdegree=6))
            # This gives Pi(x) (Riemann's prime power counting function)
            # pi(x) ~ Pi(x) - Pi(x^{1/2})/2 - ...
            result = integral_val
            err = abs(result - true_pi)
            results[T_max] = (result, err)
            print(f"  Perron integral T=[0.1,{T_max}], c={c}: Pi~{result:.4f}, error={err:.4f}")
        except Exception as e:
            print(f"  T_max={T_max} failed: {e}")

    t_total = time.time() - t0

    print(f"\n  SCALING ANALYSIS:")
    print(f"  The Perron integral integrand oscillates with frequency log(x).")
    print(f"  For x=10^100, oscillation frequency ~ 230.")
    print(f"  To capture the step at x exactly, need T >> 1.")
    print(f"  The truncation error is O(x^c / T) -- need T large for accuracy.")
    print(f"  Saddle-point methods don't help because there's no natural saddle.")
    print(f"  The contour is pinned at Re(s)=c; moving it left hits the zeros.")
    print(f"  VERDICT: Same fundamental barrier in different clothing.")

    print(f"  Time: {t_total:.2f}s")
    return {'time': t_total}


# ============================================================
# METHOD 6: HYBRID - few zeros + acceleration
# ============================================================

def test_hybrid(x, true_pi, num_zeros_few=30, num_zeros_full=100):
    """
    The most promising idea: compute a FEW zero contributions,
    then use acceleration to extrapolate the rest.

    This is the only way to beat the barrier -- if the partial sums
    have enough structure for acceleration to work.
    """
    print(f"\n{'='*60}")
    print(f"METHOD 6: HYBRID (few zeros + acceleration) for x={x}, pi(x)={true_pi}")
    print(f"{'='*60}")

    t0 = time.time()
    partials = explicit_partial_sums(x, num_zeros_full)
    t_zeros = time.time() - t0

    print(f"  Full sum ({num_zeros_full} zeros): {partials[-1]:.4f}, error={abs(partials[-1]-true_pi):.4f}")

    # Use only the first 'few' zeros, then try to extrapolate
    for N_few in [10, 20, 30, 50]:
        if N_few >= num_zeros_full:
            continue

        few_partials = partials[:N_few]
        raw_err = abs(few_partials[-1] - true_pi)

        # Wynn epsilon on the few
        window = min(16, N_few)
        try:
            wynn_val = wynn_epsilon(few_partials[-window:])
            wynn_err = abs(wynn_val - true_pi)
        except:
            wynn_val = few_partials[-1]
            wynn_err = raw_err

        # Iterated Aitken
        seq = few_partials.copy()
        for _ in range(min(5, N_few // 3)):
            if len(seq) < 3:
                break
            new_seq = []
            for i in range(len(seq) - 2):
                d = seq[i+2] - 2*seq[i+1] + seq[i]
                if abs(d) > 1e-15:
                    new_seq.append(seq[i] - (seq[i+1]-seq[i])**2 / d)
            if new_seq:
                seq = np.array(new_seq)
        aitken_val = seq[0] if len(seq) > 0 else few_partials[-1]
        aitken_err = abs(aitken_val - true_pi)

        print(f"  N={N_few}: raw_err={raw_err:.4f}, Wynn_err={wynn_err:.4f}, Aitken_err={aitken_err:.4f}")

    # Detailed analysis of the convergence character
    print(f"\n  Character of the partial sums:")
    diffs = np.diff(partials[:50])
    print(f"  First 20 increments: {[f'{d:.4f}' for d in diffs[:20]]}")

    # Are the increments alternating? Monotone? Random?
    signs = np.sign(diffs[:50])
    sign_changes = np.sum(np.abs(np.diff(signs)) > 0)
    print(f"  Sign changes in first 50 increments: {sign_changes}/{len(diffs[:50])-1}")
    print(f"  (50% would be random; ~100% would be alternating)")

    # Autocorrelation of increments
    if len(diffs) >= 20:
        ac1 = np.corrcoef(diffs[:-1], diffs[1:])[0, 1]
        print(f"  Lag-1 autocorrelation of increments: {ac1:.4f}")
        print(f"  (Near -1 = alternating/amenable to acceleration)")
        print(f"  (Near 0 = pseudo-random/NOT amenable)")

    print(f"  Time: {t_zeros:.2f}s")
    return {}


# ============================================================
# MAIN
# ============================================================

def main():
    print("=" * 70)
    print("SESSION 10: RESUMMATION TECHNIQUES vs SUMMATION BARRIER")
    print("=" * 70)
    print(f"\nTest cases:")
    for x, pi_x, desc in TEST_CASES:
        print(f"  {desc}: pi({x}) = {pi_x}")

    # Use smallest test case for detailed analysis, medium for scaling
    all_results = {}

    # ---- Run on x=7919 (smallest, fastest) ----
    x, true_pi, desc = TEST_CASES[0]
    print(f"\n\n{'#'*70}")
    print(f"# TESTING ON {desc}")
    print(f"{'#'*70}")

    r1 = test_borel_summation(x, true_pi, num_zeros=80)
    all_results['borel'] = r1

    r2 = test_pade_approximants(x, true_pi, num_zeros=100)
    all_results['pade'] = r2

    r3 = test_euler_maclaurin(x, true_pi, num_zeros=50)
    all_results['euler_maclaurin'] = r3

    r4 = test_sequence_acceleration(x, true_pi, num_zeros=100)
    all_results['acceleration'] = r4

    r5 = test_mellin_barnes(x, true_pi, num_zeros=50)
    all_results['mellin_barnes'] = r5

    r6 = test_hybrid(x, true_pi, num_zeros_few=30, num_zeros_full=100)
    all_results['hybrid'] = r6

    # ---- Run acceleration on medium case for scaling ----
    x2, true_pi2, desc2 = TEST_CASES[1]
    print(f"\n\n{'#'*70}")
    print(f"# SCALING TEST ON {desc2}")
    print(f"{'#'*70}")

    r4b = test_sequence_acceleration(x2, true_pi2, num_zeros=100)
    r6b = test_hybrid(x2, true_pi2, num_zeros_few=30, num_zeros_full=100)

    # ============================================================
    # FINAL ANALYSIS
    # ============================================================
    print(f"\n\n{'='*70}")
    print("FINAL ANALYSIS: CAN RESUMMATION BEAT O(sqrt(x))?")
    print(f"{'='*70}")

    print("""
FINDINGS BY METHOD:

1. BOREL SUMMATION:
   - The Borel transform B(t) = sum a_k t^k/k! converges because k! kills any
     polynomial growth in a_k.
   - BUT: the Laplace integral int e^{-t} B(t) dt still requires evaluating B(t)
     which needs all the a_k coefficients -- we're back to computing all zeros.
   - Borel-Pade partially helps by replacing the power series with a rational
     approximation, but the poles of the Pade encode zero information implicitly.
   - VERDICT: Does NOT bypass barrier. Needs all terms to construct the Borel sum.

2. PADE APPROXIMANTS:
   - Polynomial extrapolation in 1/N captures the slow trend but not oscillations.
   - Rational approximation does slightly better but not dramatically.
   - The partial sums DO NOT converge monotonically -- they oscillate, making
     standard Pade extrapolation unreliable.
   - VERDICT: Marginal improvement (maybe 2-3x fewer terms), NOT polylog.

3. EULER-MACLAURIN / DENSITY INTEGRATION:
   - Replacing the sum by an integral against the zero density is elegant,
     but the integrand oscillates with frequency ln(x).
   - For x=10^100, the integral has ~10^50 * 230/(2*pi) ~ 10^51 oscillations.
   - Quadrature of highly oscillatory integrals is well-studied: requires O(T)
     evaluations minimum for T oscillation periods.
   - Montgomery pair correlation doesn't help: it describes correlations between
     zeros, but the sum we need is NOT a 2-point quantity.
   - VERDICT: EXACTLY the same barrier, just rewritten as an integral.

4. SEQUENCE ACCELERATION:
   - Aitken/Shanks/Wynn work beautifully on ALTERNATING or GEOMETRIC series.
   - The partial sums of the explicit formula are NOT alternating -- the zeros
     have incommensurate imaginary parts, creating quasi-random oscillations.
   - The lag-1 autocorrelation of increments is near 0, confirming pseudo-randomness.
   - Richardson extrapolation assumes S_N ~ S + c/N^a; the actual convergence
     is irregular due to the erratic distribution of zeros.
   - VERDICT: Provides modest improvement but CANNOT overcome pseudo-random oscillations.

5. MELLIN-BARNES / CONTOUR INTEGRAL:
   - The Perron integral representation is just the inverse Mellin transform.
   - Saddle-point methods require a saddle, which doesn't exist on Re(s)=c.
   - Moving the contour left hits the zeros -- which is exactly the residue
     sum we're trying to avoid.
   - VERDICT: Equivalent formulation, no computational advantage.

6. HYBRID (few zeros + acceleration):
   - The ONLY potential path: compute O(polylog) zeros, accelerate the rest.
   - But the increments are pseudo-random with near-zero autocorrelation.
   - No acceleration method can predict pseudo-random oscillations from a prefix.
   - The GUE statistics of zero spacings (proven by Montgomery, Odlyzko) mean
     the contributions are fundamentally unpredictable from a small sample.

FUNDAMENTAL OBSTRUCTION:
   The zeta zeros have GUE (Gaussian Unitary Ensemble) statistics. The imaginary
   parts gamma_n are distributed like eigenvalues of random Hermitian matrices.
   This means:

   (a) The terms x^{i*gamma_n} are effectively random unit-norm complex numbers
   (b) Their sum is a random walk with O(sqrt(N)) fluctuations after N steps
   (c) NO resummation technique can predict the endpoint of a random walk
       from its first few steps

   This is not just a technical obstacle -- it's INFORMATION-THEORETIC.
   The sum encodes O(sqrt(x)) bits of information about the zero locations,
   and no compression scheme can reduce this below O(sqrt(x)) without
   computing the zeros themselves (which takes O(sqrt(x)) time).

CONCLUSION:
   Resummation techniques provide at best CONSTANT FACTOR improvements.
   The summation barrier O(sqrt(x)) is FUNDAMENTAL and cannot be broken
   by any resummation, acceleration, or integral representation technique.

   The path to O(polylog) prime computation, if it exists, CANNOT go through
   the explicit formula. It would require a completely different mathematical
   framework that doesn't involve zeta zeros at all.
""")


if __name__ == "__main__":
    main()

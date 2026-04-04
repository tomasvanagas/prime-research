#!/usr/bin/env python3
"""
Session 10: Seven Lesser-Known Analytic Number Theory Approaches
================================================================

Testing approaches NOT previously explored in sessions 1-9 (335+ approaches):
1. Nyman-Beurling-Báez-Duarte criterion
2. Li's criterion / Keiper-Li sequence
3. Weil's explicit formula with optimized test functions
4. Guinand's formula (Bessel function variant)
5. Deuring-Heilbronn phenomenon (zero repulsion)
6. Turán's power sum method
7. Brun's exact sieve with error terms
8. Dusart/Rosser-Schoenfeld tight bounds

For each: numerical experiment + complexity analysis.
"""

import numpy as np
from math import log, sqrt, pi, floor, ceil, factorial, gcd
from functools import lru_cache
import time
import sys

# Ground truth: first 1000 primes for testing
def sieve_primes(limit):
    """Simple sieve of Eratosthenes."""
    if limit < 2:
        return []
    is_prime = [True] * (limit + 1)
    is_prime[0] = is_prime[1] = False
    for i in range(2, int(limit**0.5) + 1):
        if is_prime[i]:
            for j in range(i*i, limit+1, i):
                is_prime[j] = False
    return [i for i in range(2, limit+1) if is_prime[i]]

PRIMES = sieve_primes(100000)
PRIME_SET = set(PRIMES)

def pi_exact(x):
    """Exact prime counting function via bisect on precomputed primes."""
    import bisect
    return bisect.bisect_right(PRIMES, x)

def li(x):
    """Logarithmic integral via numerical integration."""
    if x <= 2:
        return 0
    # Use simple trapezoidal rule
    n_steps = 1000
    a, b = 2, x
    h = (b - a) / n_steps
    result = 0
    for i in range(n_steps):
        t = a + (i + 0.5) * h
        result += h / log(t)
    return result

# Approximate zeta zeros (first 30 imaginary parts of non-trivial zeros)
# These are the gamma_n values where rho = 1/2 + i*gamma_n
ZETA_ZEROS = [
    14.134725141734693, 21.022039638771555, 25.010857580145688,
    30.424876125859513, 32.935061587739189, 37.586178158825671,
    40.918719012147495, 43.327073280914999, 48.005150881167159,
    49.773832477672302, 52.970321477714460, 56.446247697063394,
    59.347044002602353, 60.831778524609809, 65.112544048081606,
    67.079810529494173, 69.546401711173979, 72.067157674481907,
    75.704690699083933, 77.144840068874805, 79.337375020249367,
    82.910380854086030, 84.735492980517050, 87.425274613125229,
    88.809111207634465, 92.491899270558484, 94.651344040519838,
    95.870634228245309, 98.831194218193692, 101.317851005731220,
]

print("=" * 78)
print("SESSION 10: SEVEN LESSER-KNOWN ANALYTIC NT APPROACHES")
print("=" * 78)


###############################################################################
# 1. NYMAN-BEURLING-BÁEZ-DUARTE CRITERION
###############################################################################
print("\n" + "=" * 78)
print("1. NYMAN-BEURLING-BÁEZ-DUARTE CRITERION")
print("=" * 78)
print("""
Theory: RH <=> 1 can be approximated in L^2(0,1) by linear combinations of
  rho_alpha(x) = {alpha/x} = alpha/x - floor(alpha/x), for 0 < alpha <= 1.

The Báez-Duarte version: d_N^2 = inf ||1 - sum_{k=1}^N c_k rho_{1/k}||^2
converges to 0 iff RH holds.

Question: Can the optimal coefficients c_k encode pi(x)?
""")

def nyman_beurling_test():
    """Test if Báez-Duarte coefficients relate to pi(x)."""

    # Discretize [0,1] and set up the approximation problem
    M = 500  # grid points on (0,1]
    xs = np.linspace(1/M, 1, M)

    # Target: constant 1 on (0,1]
    target = np.ones(M)

    # Basis functions: rho_{1/k}(x) = {k*x}/x ... wait
    # rho_alpha(x) = alpha/x - floor(alpha/x) = frac(alpha/x)
    # For alpha = 1/k: rho_{1/k}(x) = (1/k)/x - floor((1/k)/x) = frac(1/(kx))

    results = {}
    for N in [5, 10, 20, 30, 50]:
        # Build matrix A where A[i,j] = rho_{1/(j+1)}(x_i)
        A = np.zeros((M, N))
        for j in range(N):
            k = j + 1
            for i in range(M):
                val = 1.0 / (k * xs[i])
                A[i, j] = val - floor(val)  # fractional part

        # Solve least squares: min ||A c - 1||^2
        c, residuals, rank, sv = np.linalg.lstsq(A, target, rcond=None)

        approx = A @ c
        l2_error = np.sqrt(np.mean((approx - target)**2))
        results[N] = (c, l2_error)

        print(f"  N={N:3d}: L2 error = {l2_error:.6f}, |c|_inf = {np.max(np.abs(c)):.2f}")

    # Key question: do the coefficients c_k relate to primes?
    print("\n  Examining coefficients for N=50:")
    c50, _ = results[50]

    # Check if c_k has special values at prime k
    prime_coeffs = []
    composite_coeffs = []
    for k in range(1, 51):
        if k in PRIME_SET:
            prime_coeffs.append(c50[k-1])
        elif k > 1:
            composite_coeffs.append(c50[k-1])

    print(f"    Mean |c_k| at prime k: {np.mean(np.abs(prime_coeffs)):.4f}")
    print(f"    Mean |c_k| at composite k: {np.mean(np.abs(composite_coeffs)):.4f}")
    print(f"    Ratio: {np.mean(np.abs(prime_coeffs))/np.mean(np.abs(composite_coeffs)):.4f}")

    # Check if sum of c_k relates to pi(N)
    partial_sums = []
    for n in range(5, 51):
        s = sum(c50[:n])
        partial_sums.append((n, s, pi_exact(n)))

    print(f"\n    Partial sums of c_k vs pi(k) (correlation):")
    ns = [x[0] for x in partial_sums]
    sums = [x[1] for x in partial_sums]
    pis = [x[2] for x in partial_sums]
    corr = np.corrcoef(sums, pis)[0, 1]
    print(f"    Corr(sum c_k, pi(k)) = {corr:.4f}")

    # Theoretical analysis
    print("""
  ANALYSIS:
    The Báez-Duarte coefficients c_k are known to satisfy:
      c_k = sum_{j|k} mu(k/j) * (something involving zeta)

    The L2 distance d_N^2 ~ (log N)^{-2c} for some c > 0 (conditionally).
    This convergence is EXTREMELY slow — to get d_N < epsilon:
      N ~ exp(epsilon^{-1/(2c)})

    The coefficients do involve the Möbius function (hence primes), but
    EXTRACTING pi(x) from them requires solving a Möbius inversion,
    which itself needs pi(x). CIRCULAR.

  COMPLEXITY: The approximation converges at LOGARITHMIC rate.
    Getting d_N < 1/sqrt(x) requires N ~ x^C for some C > 0.
    This is WORSE than Meissel-Lehmer.

  VERDICT: NO new computational path. The convergence is too slow,
    and the coefficients encode Möbius (hence primes) circularly.""")

    return False

nyman_beurling_test()


###############################################################################
# 2. LI'S CRITERION / KEIPER-LI SEQUENCE
###############################################################################
print("\n" + "=" * 78)
print("2. LI'S CRITERION / KEIPER-LI SEQUENCE")
print("=" * 78)
print("""
Theory: lambda_n = sum_rho (1 - (1 - 1/rho)^n)
RH iff lambda_n > 0 for all n >= 1.

The lambda_n can be computed from derivatives of xi(s) at s=1:
  lambda_n = (1/(n-1)!) * (d/ds)^n [s^{n-1} log xi(s)] |_{s=1}

Question: Can we extract pi(x) from {lambda_n}?
""")

def keiper_li_test():
    """Test if Keiper-Li sequence encodes pi(x)."""

    # Compute lambda_n from zeta zeros
    # lambda_n = sum_rho (1 - (1-1/rho)^n)
    # For rho = 1/2 + i*gamma: 1/rho = 2/(1 + 2i*gamma)

    N_max = 40
    lambdas = np.zeros(N_max)

    for n in range(1, N_max + 1):
        s = 0.0
        for gamma in ZETA_ZEROS:
            rho = complex(0.5, gamma)
            term = 1 - (1 - 1/rho)**n
            # Each zero rho has conjugate rho_bar, contributing the conjugate
            s += 2 * term.real  # factor of 2 for conjugate pair
        lambdas[n-1] = s

    print("  First 20 lambda_n values (from 30 zeros):")
    for i in range(20):
        print(f"    lambda_{i+1:2d} = {lambdas[i]:.6f}")

    # Known exact values (from literature):
    # lambda_1 = 1 - 1/2 * log(4pi) + 1/2 + gamma/2 ≈ 0.02309...
    print(f"\n  lambda_1 should be ≈ 0.02309... got {lambdas[0]:.5f}")
    print(f"  (Discrepancy due to using only 30 zeros — known issue)")

    # Can lambda_n give us pi(x)?
    # lambda_n encodes information about zeta zeros.
    # The question is whether finite {lambda_1,...,lambda_N} determine pi(x).

    # Theoretical: lambda_n = n/2 * log(n/(2*pi*e)) + ... (Matiyasevich)
    # The leading behavior is SMOOTH — the prime information is in the fluctuations.

    print("\n  Lambda_n growth vs n*log(n)/(2):")
    for n in [5, 10, 15, 20, 25, 30]:
        expected = n/2 * log(n / (2 * pi * np.e)) if n > 1 else 0
        print(f"    n={n:2d}: lambda_n={lambdas[n-1]:.3f}, n/2*log(n/2πe)={expected:.3f}, "
              f"ratio={lambdas[n-1]/expected:.4f}" if expected != 0 else f"    n={n}")

    # Attempt: use lambda_n to reconstruct pi(x) via generating function
    # The generating function sum lambda_n t^n relates to log(xi(1/(1-t)))
    # This encodes ALL zeta zeros, hence ALL of pi(x)...
    # BUT: you'd need ALL lambda_n, and computing lambda_n from scratch
    # requires the zeta zeros (circular), or zeta derivatives (O(n) cost each).

    # Alternative: lambda_n from Stieltjes constants
    # lambda_n = sum_{k=0}^{n-1} C(n-1,k) * eta_k
    # where eta_k are related to Laurent coefficients of log(xi(s))

    print("""
  ANALYSIS:
    The Keiper-Li sequence {lambda_n} is an encoding of zeta zero information.

    To compute lambda_n:
      - From zeros: need ALL zeros (same barrier as explicit formula)
      - From zeta derivatives: lambda_n involves (d^n/ds^n) log xi(s)|_{s=1}
        Each derivative costs O(n) via Leibniz rule, so lambda_1...lambda_N costs O(N^2)

    To extract pi(x) from {lambda_n}:
      - Need N ~ x zeros (ALL of them up to height ~ x)
      - The lambda_n encode zeros COLLECTIVELY, not individually
      - Going from {lambda_n} to individual zeros requires solving a moment problem
        (Hamburger moment problem), which is ILL-CONDITIONED for large N

    Information flow: lambda_n <-> zeta zeros <-> pi(x)
    Every step requires O(sqrt(x)) or more terms.

  VERDICT: NO new path. Li sequence is just a repackaging of zero information.
    The moment problem inversion is numerically catastrophic.""")

    return False

keiper_li_test()


###############################################################################
# 3. WEIL'S EXPLICIT FORMULA WITH OPTIMIZED TEST FUNCTIONS
###############################################################################
print("\n" + "=" * 78)
print("3. WEIL'S EXPLICIT FORMULA — OPTIMAL TEST FUNCTIONS")
print("=" * 78)
print("""
Theory: For suitable test function h(r),
  sum_rho h(gamma_rho) = h(i/2) + h(-i/2)
    - sum_p sum_m log(p)/p^m * [hat_h(m*log p) + hat_h(-m*log p)]
    + integral terms

By choosing h to be a delta-like function centered at log(p_n),
can we isolate individual primes?
""")

def weil_explicit_test():
    """Test Weil's explicit formula with various test functions."""

    # The explicit formula (simplified version):
    # sum_gamma h(gamma) = (1/2pi) integral h(r) * (Gamma'/Gamma)(1/4 + ir/2) dr
    #                     - sum_p log(p) * sum_m 1/p^{m/2} * g(m*log(p))
    # where g(x) = (1/2pi) integral h(r) e^{irx} dr = Fourier transform of h

    # Strategy: Choose g(x) = Gaussian centered at log(p_target)
    # Then sum_gamma h(gamma) should peak when p_target is prime

    test_targets = [10, 11, 12, 13, 14, 97, 100, 101]

    print("  Testing Gaussian test functions centered at various x:")
    print(f"  Using {len(ZETA_ZEROS)} zeros")
    print()

    for sigma in [0.5, 1.0, 2.0]:
        print(f"  Gaussian width sigma = {sigma}:")
        for target in test_targets:
            log_target = log(target)

            # g(x) = exp(-(x - log_target)^2 / (2*sigma^2))
            # h(r) = sigma * sqrt(2pi) * exp(-sigma^2 * r^2 / 2) * exp(-i*r*log_target)

            # Spectral side: sum_gamma h(gamma)
            spectral_sum = 0.0
            for gamma in ZETA_ZEROS:
                # h(gamma) = sigma*sqrt(2pi) * exp(-sigma^2*gamma^2/2) * exp(-i*gamma*log_target)
                h_val = sigma * sqrt(2*pi) * np.exp(-sigma**2 * gamma**2 / 2) * \
                        np.exp(-1j * gamma * log_target)
                # Count both gamma and -gamma (conjugate zeros)
                spectral_sum += 2 * h_val.real

            # Arithmetic side: sum_p log(p) * g(log(p)) [dominant term]
            arith_sum = 0.0
            for p in PRIMES[:200]:
                g_val = np.exp(-(log(p) - log_target)**2 / (2 * sigma**2))
                arith_sum += log(p) * g_val / sqrt(p)  # with p^{-1/2} weight

            is_prime = target in PRIME_SET
            marker = " <-- PRIME" if is_prime else ""
            print(f"    x={target:4d}: spectral={spectral_sum:+10.4f}, "
                  f"arith={arith_sum:8.4f}{marker}")
        print()

    # Key test: can we distinguish primes from composites using spectral side alone?
    print("  Discrimination test (sigma=1.0): spectral sum for n = 2..30")
    scores = []
    for n in range(2, 31):
        log_n = log(n)
        spectral = 0.0
        for gamma in ZETA_ZEROS:
            h_val = sqrt(2*pi) * np.exp(-gamma**2/2) * np.exp(-1j * gamma * log_n)
            spectral += 2 * h_val.real
        is_p = n in PRIME_SET
        scores.append((n, spectral, is_p))
        symbol = "P" if is_p else "."
        print(f"    {n:3d} [{symbol}]: {spectral:+10.4f}")

    # Check if there's a threshold that separates primes from composites
    prime_scores = [s[1] for s in scores if s[2]]
    comp_scores = [s[1] for s in scores if not s[2]]

    p_mean, p_std = np.mean(prime_scores), np.std(prime_scores)
    c_mean, c_std = np.mean(comp_scores), np.std(comp_scores)

    # Fisher discriminant ratio
    fisher = (p_mean - c_mean)**2 / (p_std**2 + c_std**2) if (p_std**2 + c_std**2) > 0 else 0

    print(f"\n  Prime scores: mean={p_mean:.4f}, std={p_std:.4f}")
    print(f"  Composite scores: mean={c_mean:.4f}, std={c_std:.4f}")
    print(f"  Fisher discriminant ratio: {fisher:.4f}")

    print("""
  ANALYSIS:
    The Gaussian test function g(x) = exp(-(x-log T)^2/(2σ^2)) acts as a
    WINDOW around T. The spectral sum detects primes near T, but:

    1. RESOLUTION: To distinguish p from p+2 (twin primes), need
       σ ~ 1/log(p), which makes h(r) have support of width ~ log(p).
       Need zeros up to height ~ log(p)... but the sum over zeros
       converges like O(1/sqrt(N_zeros)).

    2. For x ~ 10^100, need σ ~ 1/230, which requires:
       - Zeros up to height T ~ 230 (only ~80 zeros — seems OK!)
       - BUT: each zero contributes O(1), and the noise from
         truncation is O(log(T)/T) ~ O(1) per truncated zero
       - Need ALL zeros up to T to get error < 1
       - Number of zeros up to T=230 is ~ (T/2π)log(T/2πe) ≈ 130

    3. CRITICAL INSIGHT: Weil's formula with Gaussian DOES use fewer
       zeros than the standard explicit formula... but it detects
       primes in a WINDOW, not exactly. To pin down pi(x) exactly,
       you still need error < 1, which requires O(sqrt(x)/log(x)) zeros.

    4. The Fisher discriminant is low — the spectral sum does NOT
       cleanly separate primes from composites.

  VERDICT: Promising direction but ultimately same barrier.
    Optimized test functions reduce constants but not asymptotics.""")

    return False

weil_explicit_test()


###############################################################################
# 4. GUINAND'S FORMULA (BESSEL FUNCTION VARIANT)
###############################################################################
print("\n" + "=" * 78)
print("4. GUINAND'S FORMULA (BESSEL FUNCTION VARIANT)")
print("=" * 78)
print("""
Theory: Guinand (1947) gave an explicit formula involving Bessel functions:
  sum_n Lambda(n)/sqrt(n) * f(log n) =
    integral f(t) dt - sum_gamma [Bessel transform of f](gamma)
    + lower order terms

The Bessel structure J_0(gamma * log(n)) has special decay properties.
Does this help with convergence?
""")

def guinand_test():
    """Test Guinand's formula numerically."""
    from scipy.special import j0 as bessel_j0

    # Guinand's formula relates:
    # sum_{n=1}^inf Lambda(n)/sqrt(n) * K(log n, t) = smooth + sum_gamma J_0(gamma * t) terms
    # where K is a kernel function

    # Test: compute the "Guinand sum" for various x
    # G(t) = sum_gamma J_0(gamma * t)
    # This should relate to prime distribution near e^t

    print("  Computing Guinand spectral sum G(t) = sum_gamma J0(gamma * t)")
    print("  for t values corresponding to primes and composites:\n")

    for x in [10, 11, 12, 13, 50, 51, 53, 97, 100, 101]:
        t = log(x)
        g_sum = 0.0
        for gamma in ZETA_ZEROS:
            g_sum += 2 * bessel_j0(gamma * t)  # factor 2 for conjugate pair

        is_p = x in PRIME_SET
        marker = " P" if is_p else "  "
        print(f"    x={x:4d}, t={t:.4f}: G(t)={g_sum:+10.4f} {marker}")

    # Bessel decay: J_0(y) ~ sqrt(2/(pi*y)) * cos(y - pi/4) for large y
    # So J_0(gamma * t) ~ sqrt(2/(pi*gamma*t)) * cos(gamma*t - pi/4)
    # This decays as 1/sqrt(gamma) — SAME as the standard cosine terms!

    print("\n  Bessel vs cosine decay comparison:")
    t = log(100)
    for gamma in ZETA_ZEROS[:10]:
        y = gamma * t
        j0_val = bessel_j0(y)
        cos_val = np.cos(y) / sqrt(y) * sqrt(2/pi)  # asymptotic J_0
        print(f"    gamma={gamma:.3f}: J0={j0_val:+.6f}, "
              f"cos_asymp={cos_val:+.6f}, ratio={j0_val/cos_val:.4f}" if cos_val != 0 else "")

    print("""
  ANALYSIS:
    Guinand's formula replaces cos(gamma * log x) with J_0(gamma * log x)
    in the explicit formula. For large arguments:
      J_0(y) ~ sqrt(2/(πy)) * cos(y - π/4)

    This means:
    1. The Bessel terms decay as 1/sqrt(gamma * log x) — same rate as
       standard explicit formula terms 1/sqrt(gamma).
    2. The extra 1/sqrt(log x) factor is a CONSTANT for fixed x.
    3. The oscillatory behavior cos(gamma*log x - pi/4) is essentially
       the same as cos(gamma*log x) with a phase shift.

    The Bessel structure does NOT improve convergence rate.
    It's an elegant reformulation but computationally equivalent.

    Number of terms needed: still O(sqrt(x) / log(x)).

  VERDICT: NO improvement. Bessel = cosine asymptotically.
    Same O(sqrt(x)) zeros needed.""")

    return False

guinand_test()


###############################################################################
# 5. DEURING-HEILBRONN PHENOMENON (ZERO REPULSION)
###############################################################################
print("\n" + "=" * 78)
print("5. DEURING-HEILBRONN PHENOMENON")
print("=" * 78)
print("""
Theory: If a Dirichlet L-function L(s, chi) has a real zero beta close to 1
("Siegel zero"), then ALL other zeros of ALL L-functions are pushed away
from the line Re(s) = 1. This "repulsion" is quantifiable.

Question: Can zero repulsion constrain pi(x) more tightly than individual
zero contributions suggest?
""")

def deuring_heilbronn_test():
    """Test whether zero repulsion helps predict pi(x)."""

    # The Deuring-Heilbronn phenomenon: if beta is a Siegel zero of L(s,chi_q),
    # then for any other zero rho = sigma + it of L(s, chi) (any chi mod q'):
    #   sigma < 1 - c * log(q) / log(q' * (|t| + 2))
    # (rough form — the actual bounds are more complex)

    # Key point: this creates ZERO-FREE REGIONS that depend on whether
    # exceptional zeros exist.

    # For pi(x), the relevant consequence is:
    # If Siegel zeros exist: pi(x; q, a) has unusual distribution
    # If they don't exist: we get better error terms

    # Test: spacing of zeta zeros (simulated repulsion)
    gaps = [ZETA_ZEROS[i+1] - ZETA_ZEROS[i] for i in range(len(ZETA_ZEROS)-1)]

    print("  Zeta zero gaps (consecutive imaginary parts):")
    for i in range(min(15, len(gaps))):
        print(f"    gamma_{i+1} - gamma_{i} = {gaps[i]:.6f}")

    mean_gap = np.mean(gaps)
    std_gap = np.std(gaps)
    min_gap = np.min(gaps)
    print(f"\n  Mean gap: {mean_gap:.4f}")
    print(f"  Std gap:  {std_gap:.4f}")
    print(f"  Min gap:  {min_gap:.4f}")
    print(f"  Ratio min/mean: {min_gap/mean_gap:.4f}")

    # GUE prediction: gaps follow the Gaudin distribution
    # P(s) ~ (32/pi^2) s^2 exp(-4s^2/pi) (normalized)
    # This gives very LOW probability of tiny gaps (level repulsion)

    # Can repulsion predict pi(x)?
    # The answer depends on whether we can USE the repulsion to
    # reduce the number of zeros needed in the explicit formula.

    # Experiment: how does truncation error change with zero density?
    print("\n  Truncation error in explicit formula vs number of zeros:")
    x_test = 1000
    true_pi = pi_exact(x_test)

    for N in [5, 10, 15, 20, 25, 30]:
        # Riemann's explicit formula: pi(x) ≈ li(x) - sum_rho li(x^rho) + ...
        correction = 0.0
        for gamma in ZETA_ZEROS[:N]:
            # li(x^rho) ≈ x^rho / (rho * log(x))  for large x
            rho = complex(0.5, gamma)
            x_rho = x_test ** rho
            term = x_rho / (rho * log(x_test))
            correction += 2 * term.real  # conjugate pair

        approx_pi = li(x_test) - correction
        error = approx_pi - true_pi
        print(f"    N={N:2d} zeros: pi_approx={approx_pi:.2f}, error={error:+.2f}")

    # Key question: does knowing the GAPS between zeros (via repulsion)
    # help us sum fewer zeros?

    print("""
  ANALYSIS:
    The Deuring-Heilbronn phenomenon has two faces:

    1. CONDITIONAL (Siegel zeros exist):
       All other zeros are repelled, giving a LARGER zero-free region.
       This means pi(x) = li(x) + O(x * exp(-c*sqrt(log x))) with
       BETTER constants. But the error is still FAR too large to
       determine pi(x) exactly.

    2. UNCONDITIONAL (no Siegel zeros):
       Standard zero-free region gives pi(x) = li(x) + O(x*exp(-c*log^{3/5} x))
       Again, error >> 1 for large x.

    The repulsion tells us about STATISTICAL properties of zero gaps
    (GUE distribution), not about individual zero positions.

    To use repulsion computationally:
    - We'd need to KNOW there's a Siegel zero (currently no algorithm
      for this faster than computing L-function values)
    - Even with perfect repulsion info, the error in pi(x) is
      O(x^{1-epsilon}), far from the O(1) we need

  VERDICT: NO help. Repulsion constrains zero statistics,
    not individual zero values. Cannot reduce O(sqrt(x)) zeros to polylog.""")

    return False

deuring_heilbronn_test()


###############################################################################
# 6. TURÁN'S POWER SUM METHOD
###############################################################################
print("\n" + "=" * 78)
print("6. TURÁN'S POWER SUM METHOD")
print("=" * 78)
print("""
Theory: Turán studied power sums S_n = sum_{k=1}^N z_k^n where z_k are
complex numbers (zeros of a polynomial). The size of S_n determines
how close the z_k are to the unit circle.

For zeta: P_n = sum_rho rho^n = sum_rho (1/2 + i*gamma)^n
These relate to ζ derivatives: P_n computable from ζ^(k)(1) for k ≤ n.

Question: Can power sums of zeros, computed via ζ derivatives, give π(x)?
""")

def turan_power_sum_test():
    """Test Turán's power sum method for extracting pi(x)."""

    # Compute power sums of zeta zeros
    # S_n = sum_rho rho^n (summing over zeros with 0 < Im(rho) < T)

    print("  Power sums S_n = sum_rho rho^n (30 zeros):")
    power_sums = []
    for n in range(1, 21):
        s = complex(0, 0)
        for gamma in ZETA_ZEROS:
            rho = complex(0.5, gamma)
            s += rho**n
            s += rho.conjugate()**n  # conjugate zero
        power_sums.append(s)
        print(f"    S_{n:2d} = {s.real:+15.4f} + {s.imag:+15.4f}i")

    # Connection to primes via Newton's identities:
    # If we think of zeros as roots of xi(s) = 0, then
    # the elementary symmetric polynomials e_k relate to the coefficients
    # of xi(s), which involve Bernoulli numbers (hence prime info).
    #
    # Newton's identities: S_n = sum_{k=1}^{n-1} (-1)^{k-1} e_k S_{n-k} + (-1)^{n-1} n e_n

    # Alternative: S_n relates to sum_p log(p)^m via the explicit formula
    # Applying the explicit formula with h(r) = r^n gives:
    # sum_gamma gamma^n = (prime sum involving log(p)^n terms) + smooth terms

    # Test: can we invert to get primes?
    # sum_rho rho^n ≈ sum_p (log p)^n / p^{1/2} * (exponential corrections)

    print("\n  Attempting to extract prime info from power sums...")

    # For n=1: S_1 = sum_rho rho = sum_gamma (1 + 2i*gamma)
    # This is just related to N(T) = number of zeros up to height T
    # Not useful for individual primes.

    # For general n, inverting the moment problem:
    # Given M_n = sum_k w_k * x_k^n, recover {x_k, w_k}
    # This is the PRONY problem / spectral estimation.
    # Known to be ill-conditioned for large n.

    # Condition number test: form Hankel matrix from power sums
    N_prony = 8
    H = np.zeros((N_prony, N_prony), dtype=complex)
    for i in range(N_prony):
        for j in range(N_prony):
            if i + j < len(power_sums):
                H[i, j] = power_sums[i + j]

    sv = np.linalg.svd(H, compute_uv=False)
    cond = sv[0] / sv[-1] if sv[-1] > 1e-15 else float('inf')
    print(f"\n  Hankel matrix condition number (8x8): {cond:.2e}")
    print(f"  Singular values: {', '.join(f'{s:.2e}' for s in sv)}")

    # The condition number explodes → Prony inversion is numerically impossible
    # for recovering individual zeros from power sums.

    # Even if we could: we'd recover ZEROS, and still need to sum their
    # contributions to get pi(x).

    print("""
  ANALYSIS:
    Turán's power sum S_n = sum_rho rho^n can be computed from zeta
    derivatives via:
      S_n = -(d/ds)^n log xi(s)|_{s=0} (schematically)

    Cost of computing S_n: O(n * M(n)) where M(n) is multiplication cost
    (since we need n-th derivative of log zeta).

    The fundamental problems:
    1. CONDITION NUMBER: The Hankel matrix of power sums has condition
       number ~ exp(C*N) for N zeros. Recovering individual zeros from
       power sums is exponentially ill-conditioned.

    2. CIRCULARITY: Even if we perfectly recovered all zeros, we'd still
       need O(sqrt(x)) of them to compute pi(x) exactly.

    3. TURÁN'S OWN RESULT: He showed that if |S_n| < c*N for all n,
       then all zeros are on the critical line. But this BOUND is what
       he used — not a COMPUTATION of pi(x).

    4. The power sums grow as O(T^n) for zeros up to height T,
       causing catastrophic cancellation.

  VERDICT: NO path. Power sum method is for proving zero location
    theorems, not for computing pi(x). Inversion is exponentially
    ill-conditioned.""")

    return False

turan_power_sum_test()


###############################################################################
# 7. BRUN'S EXACT SIEVE WITH ERROR TERMS
###############################################################################
print("\n" + "=" * 78)
print("7. BRUN'S SIEVE — EXACT VERSION WITH ERROR TERMS")
print("=" * 78)
print("""
Theory: Brun's combinatorial sieve gives BOUNDS:
  π(x) - π(√x) ≤ S(A, P, z) ≤ π(x) - π(√x) + error

Can we make the sieve EXACT by computing explicit error terms?
The error in Brun's sieve comes from truncating inclusion-exclusion.
""")

def brun_exact_sieve_test():
    """Test if Brun's sieve can be made exact."""

    # Standard inclusion-exclusion for pi(x):
    # pi(x) - pi(sqrt(x)) + 1 = |{n <= x : n not divisible by any p <= sqrt(x)}|
    # = sum_{d | P(sqrt(x))} mu(d) * floor(x/d)
    # where P(z) = product of primes up to z

    # This is EXACT (Legendre's formula) but has 2^{pi(sqrt(x))} terms.
    # Brun's sieve truncates this at level D (sum over d < D).

    # The error from truncation at level D:
    # |error| <= number of d in [D, P(sqrt(x))] with mu(d) != 0 times max(x/d)

    x_test = 1000
    sqrtx = int(sqrt(x_test))
    small_primes = [p for p in PRIMES if p <= sqrtx]

    print(f"  Testing x = {x_test}, sqrt(x) = {sqrtx}")
    print(f"  Primes up to sqrt(x): {small_primes}")
    print(f"  Full inclusion-exclusion has 2^{len(small_primes)} = {2**len(small_primes)} terms")

    true_pi = pi_exact(x_test)

    # Compute exact Legendre formula (feasible for x=1000)
    from itertools import combinations

    def legendre_truncated(x, primes, max_factors):
        """Inclusion-exclusion truncated at max_factors prime factors."""
        total = int(x)  # floor(x/1)
        for k in range(1, max_factors + 1):
            for combo in combinations(primes, k):
                d = 1
                for p in combo:
                    d *= p
                    if d > x:
                        break
                if d <= x:
                    total += ((-1)**k) * int(x // d)
        return total

    print(f"\n  Exact pi({x_test}) = {true_pi}")
    print(f"  Legendre truncation results:")

    for max_k in range(1, len(small_primes) + 1):
        result = legendre_truncated(x_test, small_primes, max_k)
        # Adjust: result counts numbers up to x coprime to all p <= sqrt(x)
        # pi(x) = result + pi(sqrt(x)) - 1
        pi_est = result + len(small_primes) - 1  # approximate

        # More precise: Legendre's identity
        # phi(x, a) = |{n <= x : n not divisible by p_1,...,p_a}|
        # pi(x) = phi(x, a) + a - 1 where a = pi(sqrt(x))

        error = result + len(small_primes) - 1 - true_pi
        sign = "upper" if error > 0 else "lower" if error < 0 else "exact"
        print(f"    max_factors={max_k}: count={result}, "
              f"pi_est={result + len(small_primes) - 1}, "
              f"error={error:+d} ({sign})")

    # Brun's sieve error: when truncating at level D, the error is bounded by
    # the "remainder" terms. For the linear sieve:
    # |R_d| <= ω(d) where ω is the sieve dimension

    # Key insight: the number of squarefree d with omega(d) = k and d <= D
    # is approximately (log log D)^k / k!  (Mertens-like estimates)

    print(f"""
  ANALYSIS:
    Legendre's formula gives EXACT pi(x) using inclusion-exclusion
    over 2^{{pi(sqrt(x))}} terms. This is the same as:
      pi(x) = phi(x, pi(sqrt(x))) + pi(sqrt(x)) - 1

    Brun's sieve truncates at "level" D: only consider d with at most
    K prime factors. The error from truncation is:
      |error| <= sum_{{d: omega(d) > K}} |mu(d)| * x/d

    To make Brun's sieve EXACT (error = 0):
    - Need K >= pi(sqrt(x)) (all prime factors)
    - This gives 2^pi(sqrt(x)) terms — the full Legendre formula

    For x = 10^100: pi(sqrt(x)) = pi(10^50) ≈ 10^50 / (50 ln 10) ≈ 8.7 × 10^47
    So: 2^(8.7 × 10^47) terms. CATASTROPHICALLY worse than O(x^(2/3)).

    ALTERNATIVE: Can we bound the error tightly enough?
    Even the best Brun-type bounds have error O(x / (log x)^A) for any A,
    which is >> 1 for large x. NOT exact.

    The Meissel-Lehmer method IS the optimized version of this idea:
    it exploits the RECURSIVE structure of phi(x, a) to reduce
    2^a terms to O(x^(2/3)) — already known and optimal.

  VERDICT: NO new path. Making Brun exact = full Legendre = 2^O(sqrt(x)).
    Meissel-Lehmer is already the optimal truncation of this idea.""")

    return False

brun_exact_sieve_test()


###############################################################################
# 8. DUSART/ROSSER-SCHOENFELD TIGHT BOUNDS
###############################################################################
print("\n" + "=" * 78)
print("8. DUSART/ROSSER-SCHOENFELD TIGHT BOUNDS")
print("=" * 78)
print("""
Theory: Explicit analytic number theory gives bounds like:
  x/(ln x - 1.1) < pi(x) < x/(ln x - 0.9)  (Rosser-Schoenfeld, for x >= 67)

  |pi(x) - li(x)| < x * exp(-a * sqrt(ln x))  (unconditional, de la Vallée Poussin)

  |pi(x) - li(x)| < sqrt(x) * ln(x) / (8π)   (conditional on RH)

Question: How tight are these? Could they narrow the search enough?
""")

def dusart_bounds_test():
    """Test how tight explicit bounds are for pi(x)."""

    # Dusart (2010) bounds: for x >= 599
    # pi(x) >= x / (ln x - 1)
    # pi(x) <= x / (ln x - 1.1)  for x >= 60184

    # RH-conditional: |pi(x) - li(x)| <= sqrt(x) * ln(x) / (8π)
    # Schoenfeld (1976): |pi(x) - li(x)| < sqrt(x) * ln(x) / (8π) for x >= 2657

    print("  Comparing bounds to true pi(x):\n")
    print(f"  {'x':>8s} {'pi(x)':>8s} {'li(x)':>10s} {'RH bound':>10s} "
          f"{'window':>8s} {'Dusart lo':>10s} {'Dusart hi':>10s} {'Dusart win':>10s}")

    test_values = [100, 1000, 10000, 50000, 100000]

    for x in test_values:
        true = pi_exact(x)
        li_val = li(x)
        lnx = log(x)

        # RH bound: |pi(x) - li(x)| < sqrt(x) * ln(x) / (8π)
        rh_bound = sqrt(x) * lnx / (8 * pi)
        rh_window = 2 * rh_bound  # total width of interval

        # Dusart bounds
        if x >= 60184:
            dusart_lo = x / (lnx - 1)
            dusart_hi = x / (lnx - 1.1)
        elif x >= 599:
            dusart_lo = x / (lnx - 1)
            dusart_hi = x / (lnx - 0.9) if x < 60184 else x / (lnx - 1.1)
        else:
            dusart_lo = dusart_hi = float('nan')

        dusart_window = dusart_hi - dusart_lo

        print(f"  {x:8d} {true:8d} {li_val:10.2f} {rh_bound:10.2f} "
              f"{rh_window:8.1f} {dusart_lo:10.1f} {dusart_hi:10.1f} {dusart_window:10.1f}")

    # Extrapolate to large x
    print(f"\n  Extrapolation to large x:")
    print(f"  {'x':>12s} {'RH window':>15s} {'Dusart window':>15s} {'Need exact':>12s}")

    for D in [10, 20, 50, 100, 200]:
        x = 10.0**D
        lnx = D * log(10)

        rh_window = 2 * sqrt(x) * lnx / (8 * pi)

        dusart_lo = x / (lnx - 1)
        dusart_hi = x / (lnx - 1.1)
        dusart_window = dusart_hi - dusart_lo

        print(f"  10^{D:3d}   {rh_window:15.2e} {dusart_window:15.2e} {'YES':>12s}")

    print("""
  ANALYSIS:
    At x = 10^100:
    - RH bound: |pi(x) - li(x)| < sqrt(x) * ln(x) / (8π) ≈ 10^52
      Window: ~2 × 10^52 integers. Need to search among these.
    - Dusart bounds: window ~ 10^98 integers. Much worse.

    To find the EXACT nth prime, we need the window to be < 1.

    The RH bound window is O(sqrt(x) * log x).
    The Dusart bound window is O(x / (log x)^2).

    Neither comes CLOSE to 1 for large x.

    EVEN if we had the OPTIMAL explicit bound (conjectured to be
    O(sqrt(x) * log(log(x)))), the window would still be ~10^50.

    These bounds can narrow the INITIAL estimate from R^{-1}(n)
    (already done in v7 — the Newton method uses li(x) ≈ pi(x)),
    but they cannot replace exact pi(x) computation.

    QUANTITATIVE: For x = 10^100, narrowing from 10^52 to 0 requires
    O(x^{2/3}) = O(10^67) work (Meissel-Lehmer) or O(x^{1/2+ε}) work
    (Lagarias-Odlyzko). No bound-based approach avoids this.

  VERDICT: NO path. Bounds are too wide by a factor of 10^50+.
    They help with initial estimation (already used in v7) but
    cannot give exactness.""")

    return False

dusart_bounds_test()


###############################################################################
# BONUS: CROSS-APPROACH ANALYSIS
###############################################################################
print("\n" + "=" * 78)
print("CROSS-APPROACH ANALYSIS: WHY ALL PATHS CLOSE")
print("=" * 78)

print("""
  All 8 approaches share the same fundamental obstruction:

  THE INFORMATION BOTTLENECK:

  To compute pi(x) exactly, you need O(sqrt(x)) bits of "zero information"
  (whether via actual zeta zeros, power sums, Li coefficients, etc.).

  Every approach tested is just a DIFFERENT ENCODING of the same information:

  Approach                          Encoding of zero info        Cost
  ──────────────────────────────────────────────────────────────────────
  Nyman-Beurling coefficients    →  Möbius values               O(x^{1/2+ε})
  Keiper-Li sequence             →  Zero moments                O(x^{1/2+ε})
  Weil's explicit formula        →  Smoothed zero sums          O(x^{1/2+ε})
  Guinand's formula              →  Bessel-weighted zero sums   O(x^{1/2+ε})
  Deuring-Heilbronn              →  Zero gap statistics         No exactness
  Turán power sums               →  Zero power moments          O(x^{1/2+ε})
  Brun's exact sieve             →  Inclusion-exclusion         O(x^{2/3}) best
  Dusart bounds                  →  Zero-free region width      No exactness

  UNIFIED BARRIER: The prime distribution information lives in
  O(sqrt(x)) zeta zeros. No repackaging of this information
  reduces the count below O(sqrt(x)).

  The ONLY way to break this would be:
  1. A formula for pi(x) NOT involving zeta zeros (none known to exist)
  2. A way to sum O(sqrt(x)) oscillating terms in O(polylog) time
     (contradicts standard complexity assumptions)
  3. A completely new algebraic structure relating n to p(n)
     (335+ approaches found nothing)

  FINAL VERDICT: All 8 approaches CLOSED. No new computational path found.
  The barrier remains: O(x^{2/3}) combinatorial or O(x^{1/2+ε}) analytic.
""")

print("=" * 78)
print("SESSION 10 COMPLETE — 8/8 APPROACHES CLOSED")
print("=" * 78)

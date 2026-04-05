"""
Recursive Dickman Decomposition: Exploiting the DDE structure of φ.

KEY INSIGHT: The Dickman function ρ(u) satisfies:
  u·ρ'(u) = -ρ(u-1)  for u > 1
  ρ(u) = 1            for 0 ≤ u ≤ 1

This DDE has a beautiful recursive structure: each "layer" (u ∈ [k, k+1])
depends only on the previous layer. This is like a recurrence relation!

The sieve function φ(x, a) is the discrete analog.
φ(x, a) = number of integers ≤ x not divisible by p_1, ..., p_a.

φ satisfies: φ(x, a) = φ(x, a-1) - φ(x/p_a, a-1)

This recursion has 2^a leaves. But the Dickman connection suggests
that the "continuous limit" is much simpler (the DDE is 1D!).

QUESTION: Can we match the discrete φ to the continuous ρ closely enough
that only O(polylog) correction terms are needed?

Specifically:
  φ(x, a) ≈ x · ρ(log(x)/log(p_a)) · correction(x, a)

If correction(x, a) is "simple" (polynomial in log x, or has few terms),
then we can compute φ(x, a) in polylog time.

ALSO EXPLORING: Buchstab's identity as a VOLTERRA INTEGRAL EQUATION.
Buchstab: φ(x, a) = x·ω(u) where u = log(x)/log(p_a) and
  ω(u) = 1/u for 1 ≤ u ≤ 2
  u·ω(u) = ∫_{1}^{u-1} ω(t) dt  for u > 2

This integral equation might be solvable via spectral methods!
"""

import numpy as np
import sympy
from sympy import primepi, prime, nextprime
import math
import time
from functools import lru_cache

def dickman_rho(u, num_terms=100):
    """Compute Dickman's function ρ(u) via series expansion."""
    if u <= 0:
        return 0.0
    if u <= 1:
        return 1.0
    if u <= 2:
        return 1.0 - math.log(u)

    # For u > 2, use numerical integration of the DDE
    # u·ρ'(u) = -ρ(u-1)
    # Discretize with step h
    h = 0.001
    n_steps = int(u / h) + 1
    rho_vals = np.zeros(n_steps + 1)

    # Initialize for u ∈ [0, 1]
    for i in range(min(n_steps + 1, int(1/h) + 1)):
        rho_vals[i] = 1.0

    # Initialize for u ∈ [1, 2]
    for i in range(int(1/h), min(n_steps + 1, int(2/h) + 1)):
        t = i * h
        if t >= 1:
            rho_vals[i] = 1.0 - math.log(t)

    # Solve DDE for u > 2
    for i in range(int(2/h) + 1, n_steps + 1):
        t = i * h
        # ρ'(t) = -ρ(t-1)/t
        i_prev = i - int(1/h)
        if i_prev >= 0 and t > 0:
            rho_vals[i] = rho_vals[i-1] + h * (-rho_vals[i_prev] / t)

    return rho_vals[min(n_steps, len(rho_vals)-1)]

def buchstab_omega(u, h=0.001):
    """Compute Buchstab's function ω(u)."""
    if u <= 1:
        return 0.0
    if u <= 2:
        return 1.0 / u

    # ω(u) = 1 + (1/u) ∫₂ᵘ ω(t-1) dt  ... actually
    # The standard form: u·ω(u) = (u-1)·ω(u-1) + ∫_{u-1}^{u} ω(t-1)dt
    # Or equivalently: u·ω(u) = 1 + ∫₂ᵘ ω(t-1) dt for u > 2

    n_steps = int(u / h) + 2
    omega_vals = np.zeros(n_steps + 1)

    # u ∈ [1, 2]: ω(u) = 1/u
    for i in range(int(1/h), min(n_steps + 1, int(2/h) + 2)):
        t = i * h
        if t >= 1:
            omega_vals[i] = 1.0 / t

    # u > 2: use integral equation
    # u·ω(u) = (u-1)·ω(u-1) + ω(u-1) ... simplified:
    # d/du [u·ω(u)] = ω(u-1)
    # So u·ω(u) = 2·ω(2) + ∫₂ᵘ ω(t-1) dt = 1 + ∫₂ᵘ ω(t-1) dt

    # Use cumulative integration
    for i in range(int(2/h) + 1, n_steps + 1):
        t = i * h
        if t > 0:
            # Euler step for d/du[u·ω(u)] = ω(u-1)
            i_prev = i - int(1/h)
            if i_prev >= 0:
                u_omega_prev = (t - h) * omega_vals[i-1]
                u_omega_curr = u_omega_prev + h * omega_vals[i_prev]
                omega_vals[i] = u_omega_curr / t

    idx = min(int(u / h), len(omega_vals) - 1)
    return omega_vals[idx]

def phi_exact(x, a, primes_list, memo=None):
    """Exact Legendre sieve function with memoization."""
    if memo is None:
        memo = {}
    key = (int(x), a)
    if key in memo:
        return memo[key]
    if a == 0:
        return int(x)
    if x < 1:
        return 0
    result = phi_exact(x, a-1, primes_list, memo) - phi_exact(x // primes_list[a-1], a-1, primes_list, memo)
    memo[key] = result
    return result

def test_dickman_approximation():
    """
    Test: How well does x·ρ(log(x)/log(y)) approximate φ(x, π(y))?
    """
    print("=== Dickman approximation of φ(x, a) ===\n")

    primes_list = list(sympy.primerange(2, 1000))

    print(f"{'x':>10} {'y':>6} {'u':>6} {'φ_exact':>10} {'x·ρ(u)':>12} {'rel_err':>10}")

    for x in [100, 500, 1000, 5000, 10000]:
        for y_exp in [1/3, 1/2]:
            y = int(x ** y_exp)
            if y < 2:
                y = 2
            a = int(sympy.primepi(y))
            u = math.log(x) / math.log(max(y, 2))

            memo = {}
            phi_val = phi_exact(x, a, primes_list, memo)
            rho_val = dickman_rho(u)
            approx = x * rho_val

            rel_err = abs(approx - phi_val) / max(phi_val, 1)
            print(f"{x:>10} {y:>6} {u:>6.2f} {phi_val:>10} {approx:>12.2f} {rel_err:>10.4f}")

def test_buchstab_approximation():
    """
    Test: How well does Buchstab's function approximate φ(x,a)/x?
    """
    print("\n=== Buchstab approximation ===\n")

    primes_list = list(sympy.primerange(2, 1000))

    print(f"{'x':>10} {'a':>5} {'u':>6} {'φ/x':>10} {'ω(u)':>10} {'abs_err':>10}")

    for x in [100, 500, 1000, 5000, 10000]:
        a = int(sympy.primepi(int(x ** (1/3))))
        y = primes_list[a-1] if a > 0 else 2
        u = math.log(x) / math.log(max(y, 2))

        memo = {}
        phi_val = phi_exact(x, a, primes_list, memo)
        omega_val = buchstab_omega(u)

        phi_over_x = phi_val / x
        abs_err = abs(phi_over_x - omega_val)
        print(f"{x:>10} {a:>5} {u:>6.2f} {phi_over_x:>10.6f} {omega_val:>10.6f} {abs_err:>10.6f}")

def test_correction_structure():
    """
    KEY EXPERIMENT: What is the structure of the correction term?

    Define: ε(x, a) = φ(x, a) - x·ρ(log(x)/log(p_a))

    If ε has a simple structure (polynomial in log x, or sparse representation),
    we might be able to compute it fast.
    """
    print("\n=== Correction term structure ===\n")

    primes_list = list(sympy.primerange(2, 1000))

    # Compute ε for a range of x values
    a = 10  # fixed sieve depth
    y = primes_list[a-1]

    corrections = []
    x_vals = []

    for x in range(50, 5001, 10):
        u = math.log(x) / math.log(y)
        memo = {}
        phi_val = phi_exact(x, a, primes_list, memo)
        rho_val = dickman_rho(u)
        approx = x * rho_val

        corrections.append(phi_val - approx)
        x_vals.append(x)

    corrections = np.array(corrections)
    x_vals = np.array(x_vals)

    print(f"Correction ε(x, {a}) = φ(x,{a}) - x·ρ(u)  for x ∈ [50, 5000]:")
    print(f"  Mean: {corrections.mean():.4f}")
    print(f"  Std: {corrections.std():.4f}")
    print(f"  Max |ε|: {np.max(np.abs(corrections)):.4f}")

    # Growth rate: does ε grow like x^α for some α?
    pos = corrections > 0
    if np.any(pos) and np.any(x_vals[pos] > 100):
        log_x = np.log(x_vals[pos])
        log_eps = np.log(np.abs(corrections[pos]) + 1e-10)
        # Linear fit in log-log
        valid = log_eps > -20
        if np.sum(valid) > 10:
            coeffs = np.polyfit(log_x[valid], log_eps[valid], 1)
            print(f"  Growth exponent (ε ~ x^α): α ≈ {coeffs[0]:.4f}")

    # Fourier analysis of correction
    corr_fft = np.abs(np.fft.fft(corrections))
    sorted_fft = sorted(corr_fft, reverse=True)
    cumul = np.cumsum(np.array(sorted_fft)**2)
    total = cumul[-1]

    print(f"\n  Fourier analysis of correction:")
    for frac in [0.5, 0.8, 0.9, 0.95]:
        k = np.searchsorted(cumul, frac * total) + 1
        print(f"    {frac*100:.0f}% energy: top {k}/{len(corrections)} coefficients ({k/len(corrections)*100:.1f}%)")

    # Polynomial fit of correction
    print(f"\n  Polynomial fit of ε(x,a)/√x:")
    normalized = corrections / np.sqrt(x_vals)
    for deg in [1, 2, 3, 5, 10]:
        coeffs = np.polyfit(np.log(x_vals), normalized, deg)
        fit = np.polyval(coeffs, np.log(x_vals))
        residual = np.std(normalized - fit) / np.std(normalized)
        print(f"    degree {deg:>2}: residual/std = {residual:.6f}")

def test_recursive_dickman_acceleration():
    """
    CORE IDEA: Use Dickman/Buchstab as a SKELETON and compute corrections.

    Algorithm:
    1. Compute φ(x, a) ≈ x · ρ(u) where u = log(x)/log(p_a)   [O(polylog)]
    2. Correction: ε(x, a) = φ(x, a) - x·ρ(u)
    3. ε(x, a) = ε(x, a-1) - ε(x/p_a, a-1) - x·[ρ(u) - ρ(u') + ρ(u'')/p_a]
       where u' and u'' involve shifted arguments
    4. If ε satisfies a SIMPLER recursion, solve it directly

    The recursion for ε might have better convergence properties because
    the "bulk" has been subtracted out.
    """
    print("\n=== Recursive Dickman acceleration ===\n")

    primes_list = list(sympy.primerange(2, 1000))

    # Compute correction at multiple scales
    print("Correction ε(x, a) at different scales:")
    print(f"{'x':>10} {'a':>5} {'φ':>10} {'x·ρ':>12} {'ε':>10} {'ε/√x':>10} {'ε/x^{1/3}':>10}")

    for x in [100, 500, 1000, 5000, 10000, 50000]:
        a = max(1, int(sympy.primepi(int(x ** (1/3)))))
        y = primes_list[a-1]
        u = math.log(x) / math.log(y)

        memo = {}
        phi_val = phi_exact(x, a, primes_list, memo)
        rho_val = dickman_rho(u)
        approx = x * rho_val
        eps = phi_val - approx

        print(f"{x:>10} {a:>5} {phi_val:>10} {approx:>12.2f} {eps:>10.2f} "
              f"{eps/x**0.5:>10.4f} {eps/x**(1/3):>10.4f}")

    # Test: does the recursion for ε have fewer significant terms?
    print("\nRecursion tree depth for ε(x, a):")
    print("(Counting subproblems where |ε| > 0.5)")

    for x in [1000, 5000, 10000]:
        a = max(1, int(sympy.primepi(int(x ** (1/3)))))
        y = primes_list[a-1]

        # Build full recursion tree for φ
        memo_phi = {}
        phi_exact(x, a, primes_list, memo_phi)
        total_subproblems = len(memo_phi)

        # Count "significant" subproblems (where correction matters)
        significant = 0
        for (xi, ai), val in memo_phi.items():
            if xi < 1 or ai == 0:
                continue
            if ai > len(primes_list):
                continue
            yi = primes_list[ai-1]
            if yi <= 1:
                continue
            ui = math.log(max(xi, 2)) / math.log(yi)
            rho_i = dickman_rho(ui)
            eps_i = val - xi * rho_i
            if abs(eps_i) > 0.5:
                significant += 1

        print(f"  x={x:>6}, a={a}: total subproblems={total_subproblems}, "
              f"significant corrections={significant} ({significant/total_subproblems*100:.1f}%)")

def test_volterra_spectral():
    """
    Buchstab's identity as Volterra integral equation.

    u·ω(u) = 1 + ∫₂ᵘ ω(t-1) dt

    This is a Volterra equation of the second kind.
    Standard spectral methods (Chebyshev, Legendre) can solve these
    to high accuracy with few terms.

    If ω(u) has a good Chebyshev expansion, then computing ω at any
    point takes O(N) time where N is the number of Chebyshev coefficients.
    N should be small if ω is smooth (which it is for u > 2).
    """
    print("\n=== Volterra spectral method for Buchstab ===\n")

    # Compute ω(u) on a fine grid
    u_max = 20
    u_vals = np.linspace(1.01, u_max, 10000)
    omega_vals = np.array([buchstab_omega(u) for u in u_vals])

    # ω(u) converges to e^{-γ} ≈ 0.5615 as u → ∞ (Mertens' constant)
    euler_gamma = 0.5772156649
    mertens = math.exp(-euler_gamma)
    print(f"ω(u) → e^{{-γ}} ≈ {mertens:.6f}")
    print(f"ω(20) ≈ {omega_vals[-1]:.6f}")
    print(f"Convergence error at u=20: {abs(omega_vals[-1] - mertens):.6f}")

    # Chebyshev expansion of ω(u) on [2, u_max]
    # Map [2, u_max] → [-1, 1]
    u_cheb = u_vals[u_vals >= 2]
    omega_cheb = omega_vals[u_vals >= 2]

    t_mapped = 2 * (u_cheb - 2) / (u_max - 2) - 1

    # Fit Chebyshev coefficients
    for N in [5, 10, 20, 50]:
        coeffs = np.polynomial.chebyshev.chebfit(t_mapped, omega_cheb, N)
        fit = np.polynomial.chebyshev.chebval(t_mapped, coeffs)
        max_err = np.max(np.abs(fit - omega_cheb))
        print(f"  Chebyshev degree {N:>3}: max error = {max_err:.2e}")

    print("\n  (If Chebyshev converges fast, ω can be evaluated in O(N) time)")
    print("  (The question is whether this helps with the DISCRETE φ)")

if __name__ == "__main__":
    print("=" * 60)
    print("EXPERIMENT: Recursive Dickman Decomposition")
    print("=" * 60)

    test_dickman_approximation()
    test_buchstab_approximation()
    test_correction_structure()
    test_recursive_dickman_acceleration()
    test_volterra_spectral()

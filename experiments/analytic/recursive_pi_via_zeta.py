"""
Session 6: Novel Recursive Approach via Zeta Functional Equation

RADICAL IDEA: The functional equation ζ(s) = 2^s π^{s-1} sin(πs/2) Γ(1-s) ζ(1-s)
relates ζ at s and 1-s. Can we exploit this self-referential structure?

More specifically:
  log ζ(s) = Σ_p Σ_k p^{-ks}/k for Re(s) > 1

This is the Euler product. The LOG of zeta is the prime zeta function (up to k terms).

APPROACH: Telescoping product / recursive formula
  π(x) = li(x) - Σ_ρ li(x^ρ) + correction

The zeros ρ satisfy N(T) ~ T/(2π) * log(T/(2πe)). Can we use
ZERO-FREE REGIONS to bound the contribution of distant zeros
and get a formula that's exact for specific n?

KEY NEW IDEA: What if we don't need ALL zeros, but can compute
their COLLECTIVE effect via a contour integral that avoids
individual zero computation?

The Perron formula: π(x) = (1/2πi) ∫ -ζ'(s)/ζ(s) * x^s/s ds

This integral picks up residues at poles and zeros. If we evaluate
it numerically using a contour that AVOIDS all zeros, we get a
formula that doesn't require knowing individual zeros!

The challenge: ζ'/ζ has poles at every zero, so any contour that
avoids them must go around them, which effectively counts them.

But what about a SHIFTED contour? If we integrate along Re(s) = σ > 1,
we get a formula involving ζ values on the convergent side.
"""

import numpy as np
from mpmath import mp, mpf, mpc, pi, log, exp, zeta, zetazero, gamma
from mpmath import quad, im, re, fabs, loggamma, sqrt, li, cos, sin
import time

mp.dps = 30

def perron_integral_shifted(x, sigma=2.0, T=50):
    """
    Compute π(x) via Perron's formula with a shifted contour at Re(s) = σ > 1.

    ψ(x) = (1/2πi) ∫_{σ-iT}^{σ+iT} -ζ'(s)/ζ(s) * x^s/s ds

    where ψ(x) = Σ_{p^k ≤ x} ln(p) is the Chebyshev function.
    Then π(x) can be recovered from ψ via Möbius inversion.
    """
    x = mpf(x)

    def integrand(t):
        s = mpc(sigma, t)
        try:
            # -ζ'(s)/ζ(s) * x^s / s
            z = zeta(s)
            if abs(z) < 1e-50:
                return mpc(0, 0)
            # Numerical derivative: ζ'(s) ≈ (ζ(s+h) - ζ(s-h))/(2h)
            h = mpc(0, mpf('1e-8'))
            zeta_deriv = (zeta(s + h) - zeta(s - h)) / (2 * h)
            return (-zeta_deriv / z) * exp(s * log(x)) / s
        except Exception:
            return mpc(0, 0)

    # Integrate from -T to T
    result = quad(lambda t: integrand(t), [mpf(-T), mpf(T)])
    psi_x = re(result) / (2 * pi)

    return psi_x

def test_perron():
    """Test the Perron integral approach."""
    print("="*70)
    print("PERRON INTEGRAL WITH SHIFTED CONTOUR")
    print("="*70)

    from sympy import primepi

    for x in [20, 50, 100]:
        true_pi = int(primepi(x))

        for T in [20, 50, 100]:
            for sigma in [1.5, 2.0, 3.0]:
                psi_est = perron_integral_shifted(x, sigma=sigma, T=T)
                # ψ(x) ≈ x - Σ_ρ x^ρ/ρ - ...
                # π(x) ≈ ψ(x)/ln(x) + ψ(x^{1/2})/(2*ln(x)) + ...
                # Rough: π(x) ≈ ψ(x)/ln(x)
                pi_est = float(psi_est) / float(log(mpf(x)))

                print(f"  x={x}, σ={sigma}, T={T}: "
                      f"ψ̂={float(psi_est):.2f}, π̂≈{pi_est:.2f}, true π={true_pi}")

def experiment_prime_generating_function():
    """
    Another approach: the prime generating function
    F(z) = Σ_n p(n) * z^n

    If we know F(z) near z=0, we can extract p(n) via Taylor coefficients.

    But what IS F(z)? Can we compute it without knowing the primes?

    Using the relation: F(z) = Σ_k (Σ_n [n-th prime] z^n)
    This doesn't help directly.

    Alternative: prime INDICATOR generating function
    χ(z) = Σ_p z^p (sum over primes)

    We know χ(z) = Σ_p z^p, and from the Euler product:
    log(1/(1-z^p)) appears in Π_p 1/(1-z^p) = "partition function for prime parts"

    Can we extract χ from the partition generating function?
    """
    print("\n" + "="*70)
    print("PRIME GENERATING FUNCTION ANALYSIS")
    print("="*70)

    from sympy import prime

    # Compute the prime generating function F(z) for small z
    N = 100
    primes = [int(prime(n)) for n in range(1, N + 1)]

    # F(z) = Σ p(n) z^n
    # At z = 1/2: F(1/2) = Σ p(n) / 2^n
    # This converges since p(n) ~ n*ln(n) < 2^n for large n

    z_values = [0.1, 0.3, 0.5, 0.7, 0.9]
    for z in z_values:
        F_z = sum(p * z**n for n, p in enumerate(primes, 1))
        print(f"  F({z}) = {F_z:.6f}")

    # The radius of convergence of F(z) is 1 (since p(n) grows polynomially)
    # Near z=1: F(z) ~ Σ n*ln(n) * z^n ≈ -z*d/dz[Σ ln(n)*z^n]
    # This diverges as z→1, giving F(z) ~ 1/(1-z)^2 * ln(1/(1-z))

    # Key question: does F(z) have a nice analytic continuation beyond |z|=1?
    # Answer: NO, because p(n) is "random-like". The generating function
    # has the unit circle as a natural boundary (lacunary-like behavior).

    print("\n  Growth of F(z) near z=1:")
    for z in [0.9, 0.95, 0.99, 0.999]:
        F_z = sum(p * z**n for n, p in enumerate(primes, 1))
        expected = 1/(1-z)**2 * np.log(1/(1-z))
        print(f"  F({z:.3f}) = {F_z:.1f}, "
              f"approx 1/(1-z)²·ln(1/(1-z)) = {expected:.1f}, "
              f"ratio = {F_z/expected:.3f}")

def experiment_moebius_transform_shortcut():
    """
    The Möbius inversion formula: if g(n) = Σ_{d|n} f(d), then f(n) = Σ_{d|n} μ(n/d) g(d).

    Apply this to: Λ(n) = Σ_{d|n} μ(d) * ln(n/d) (von Mangoldt)

    We have: ψ(x) = Σ_{n≤x} Λ(n) = x - Σ_ρ x^ρ/ρ - ln(2π)
    And: θ(x) = Σ_{p≤x} ln(p) = ψ(x) - ψ(x^{1/2}) - ψ(x^{1/3}) - ...
    And: ��(x) = Σ_k μ(k)/k * J(x^{1/k}) where J(x) = li(x) + Σ li(x^ρ) + ...

    New idea: What if we compute ψ(x) mod small integers?

    ψ(x) mod m = Σ_{n≤x} Λ(n) mod m

    If Λ(n) has structure mod m, this might be faster.
    """
    print("\n" + "="*70)
    print("MÖBIUS SHORTCUT: ψ(x) MOD m")
    print("="*70)

    from sympy import factorint, mobius, primepi

    # Compute ψ(x) and π(x) for comparison
    def von_mangoldt(n):
        """Λ(n) = ln(p) if n = p^k, else 0."""
        if n <= 1:
            return 0
        factors = factorint(n)
        if len(factors) == 1:
            p = list(factors.keys())[0]
            return float(np.log(p))
        return 0

    def chebyshev_psi(x):
        """ψ(x) = Σ_{n≤x} Λ(n)."""
        return sum(von_mangoldt(n) for n in range(2, int(x) + 1))

    # Test ψ(x) modular structure
    for x in [100, 500]:
        psi = chebyshev_psi(x)
        true_pi = int(primepi(x))

        print(f"\n  x = {x}: ψ(x) = {psi:.4f}, π(x) = {true_pi}")
        print(f"  ψ(x)/x = {psi/x:.6f} (should → 1)")
        print(f"  ψ(x) - x = {psi - x:.4f} (should be small)")

        # Can we compute ψ(x) mod m for small m?
        # This doesn't make sense for real-valued ψ, but we can look at
        # floor(ψ(x)) mod m
        psi_int = int(round(psi))
        for m in [2, 3, 5, 7, 10]:
            print(f"    floor(ψ) mod {m} = {psi_int % m}")

def experiment_optimal_approximation():
    """
    Instead of finding an exact formula, find the BEST possible
    approximation using a bounded number of operations.

    Question: Given K operations (additions, multiplications, special
    function evaluations), what is the MINIMUM possible error in
    approximating p(n)?

    Approach: Use optimization to find the best formula tree
    with K nodes.
    """
    print("\n" + "="*70)
    print("OPTIMAL K-OPERATION APPROXIMATION")
    print("="*70)

    from sympy import prime
    from scipy.optimize import minimize

    N = 5000
    ns = np.arange(2, N + 2, dtype=np.float64)
    true_primes = np.array([int(prime(n)) for n in range(2, N + 2)], dtype=np.float64)

    # Parametric families of increasing complexity

    # Family 1: n * (ln(n) + ln(ln(n)) + a + b/ln(n) + c*ln(ln(n))/ln(n) + d/ln(n)^2)
    def family1(params, ns):
        a, b, c, d = params
        ln_n = np.log(ns)
        ln_ln_n = np.log(ln_n)
        return ns * (ln_n + ln_ln_n + a + b/ln_n + c*ln_ln_n/ln_n + d/ln_n**2)

    def loss1(params):
        pred = family1(params, ns)
        return np.mean((pred - true_primes)**2)

    res1 = minimize(loss1, [-1, -2, 1, 6], method='Nelder-Mead',
                    options={'maxiter': 10000, 'xatol': 1e-10})
    pred1 = family1(res1.x, ns)
    err1 = pred1 - true_primes
    exact1 = np.sum(np.abs(err1) < 0.5)

    print(f"\nFamily 1 (4 params): a={res1.x[0]:.6f}, b={res1.x[1]:.6f}, "
          f"c={res1.x[2]:.6f}, d={res1.x[3]:.6f}")
    print(f"  RMS error: {np.sqrt(np.mean(err1**2)):.2f}")
    print(f"  Max |error|: {np.max(np.abs(err1)):.2f}")
    print(f"  Exact: {exact1}/{N} ({100*exact1/N:.1f}%)")

    # Family 2: Same + extra terms
    def family2(params, ns):
        a, b, c, d, e, f = params
        ln_n = np.log(ns)
        ln_ln_n = np.log(ln_n)
        return ns * (ln_n + ln_ln_n + a + b/ln_n + c*ln_ln_n/ln_n +
                     d/ln_n**2 + e*ln_ln_n**2/ln_n**2 + f*ln_ln_n/ln_n**2)

    def loss2(params):
        pred = family2(params, ns)
        return np.mean((pred - true_primes)**2)

    res2 = minimize(loss2, [-1, -2, 1, 6, 0, 0], method='Nelder-Mead',
                    options={'maxiter': 20000, 'xatol': 1e-10})
    pred2 = family2(res2.x, ns)
    err2 = pred2 - true_primes
    exact2 = np.sum(np.abs(err2) < 0.5)

    print(f"\nFamily 2 (6 params):")
    print(f"  Params: {[f'{x:.6f}' for x in res2.x]}")
    print(f"  RMS error: {np.sqrt(np.mean(err2**2)):.2f}")
    print(f"  Max |error|: {np.max(np.abs(err2)):.2f}")
    print(f"  Exact: {exact2}/{N} ({100*exact2/N:.1f}%)")

    # Family 3: Add oscillatory terms
    def family3(params, ns):
        a, b, c, d, e, f, amp1, freq1, amp2, freq2 = params
        ln_n = np.log(ns)
        ln_ln_n = np.log(ln_n)
        base = ns * (ln_n + ln_ln_n + a + b/ln_n + c*ln_ln_n/ln_n +
                     d/ln_n**2 + e*ln_ln_n**2/ln_n**2 + f*ln_ln_n/ln_n**2)
        osc = amp1 * np.sqrt(ns) * np.sin(freq1 * ln_n) + \
              amp2 * np.sqrt(ns) * np.cos(freq2 * ln_n)
        return base + osc

    def loss3(params):
        pred = family3(params, ns)
        return np.mean((pred - true_primes)**2)

    # Use first zeta zeros as initial frequencies
    res3 = minimize(loss3,
                    list(res2.x) + [0.1, 14.13, 0.1, 21.02],
                    method='Nelder-Mead',
                    options={'maxiter': 50000, 'xatol': 1e-10})
    pred3 = family3(res3.x, ns)
    err3 = pred3 - true_primes
    exact3 = np.sum(np.abs(err3) < 0.5)

    print(f"\nFamily 3 (10 params, with oscillatory):")
    print(f"  RMS error: {np.sqrt(np.mean(err3**2)):.2f}")
    print(f"  Max |error|: {np.max(np.abs(err3)):.2f}")
    print(f"  Exact: {exact3}/{N} ({100*exact3/N:.1f}%)")
    print(f"  Best freq1={res3.x[7]:.4f} (cf γ₁=14.13)")
    print(f"  Best freq2={res3.x[9]:.4f} (cf γ₂=21.02)")

    # Test generalization
    test_ns = np.arange(N + 2, N + 1002, dtype=np.float64)
    test_primes = np.array([int(prime(n)) for n in range(N + 2, N + 1002)], dtype=np.float64)
    test_pred3 = family3(res3.x, test_ns)
    test_err3 = test_pred3 - test_primes
    test_exact3 = np.sum(np.abs(test_err3) < 0.5)

    print(f"\n  Generalization test (n={N+2}..{N+1001}):")
    print(f"    RMS error: {np.sqrt(np.mean(test_err3**2)):.2f}")
    print(f"    Exact: {test_exact3}/1000 ({100*test_exact3/1000:.1f}%)")

def main():
    print("Session 6: Recursive π via Zeta & Novel Approaches")
    print(f"Started: {time.strftime('%H:%M:%S')}\n")

    t0 = time.time()

    test_perron()
    experiment_prime_generating_function()
    experiment_moebius_transform_shortcut()
    experiment_optimal_approximation()

    elapsed = time.time() - t0
    print(f"\n{'='*70}")
    print(f"Total time: {elapsed:.1f}s")

if __name__ == '__main__':
    main()

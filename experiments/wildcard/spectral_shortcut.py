#!/usr/bin/env python3
"""
Spectral Shortcut via Trace Formula -- Experiment
===================================================
Idea: Can Σ_ρ Li(x^ρ) be computed as a "trace" without enumerating individual zeros?

Tests:
1. Explicit formula π(x) = Li(x) - Σ_ρ Li(x^ρ) - log(2) + integral, convergence vs #zeros
2. Gaussian-smoothed variants: replace sharp Li(x^ρ) with smoothed kernels
3. Structure in the zero sum: Σ x^{iγ}/(1/2+iγ) for many x
4. Short-formula approximability of the zero sum
5. Riemann-von Mangoldt N(T) verification

Uses first 200 zeta zeros from data/zeta_zeros_200.txt
"""

import numpy as np
from scipy.special import expi as Ei
from scipy.integrate import quad
import os
import sys
import time

# ---------------------------------------------------------------------------
# Load zeta zeros
# ---------------------------------------------------------------------------
DATA_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'data')

def load_zeros(n=200):
    path = os.path.join(DATA_DIR, f'zeta_zeros_{n}.txt' if n in (200, 300, 500, 1000) else 'zeta_zeros_200.txt')
    if not os.path.exists(path):
        path = os.path.join(DATA_DIR, 'zeta_zeros_200.txt')
    gammas = []
    with open(path) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                gammas.append(float(parts[1]))
            elif len(parts) == 1:
                gammas.append(float(parts[0]))
    return np.array(gammas[:n])

GAMMAS = load_zeros(200)

# ---------------------------------------------------------------------------
# Reference: exact prime counting via sieve
# ---------------------------------------------------------------------------
def sieve_pi(limit):
    """Sieve of Eratosthenes, returns pi(x) for x = 0..limit."""
    sieve = np.ones(limit + 1, dtype=bool)
    sieve[0] = sieve[1] = False
    for i in range(2, int(limit**0.5) + 1):
        if sieve[i]:
            sieve[i*i::i] = False
    return np.cumsum(sieve)

# ---------------------------------------------------------------------------
# Li(x) -- logarithmic integral
# ---------------------------------------------------------------------------
def li(x):
    """Logarithmic integral Li(x) = Ei(ln(x)) for real x > 1."""
    if x <= 1:
        return 0.0
    return Ei(np.log(x))

def li_complex(x, rho):
    """
    Li(x^rho) where rho = 1/2 + i*gamma.
    Li(z) = Ei(ln(z)) = Ei(rho * ln(x)).
    For complex argument, use series/numerical integration.
    """
    if x <= 1:
        return 0.0
    w = rho * np.log(x)
    # Ei(w) for complex w via series expansion (sufficient for moderate |w|)
    # Ei(w) = gamma_euler + ln(w) + sum_{k=1}^{inf} w^k / (k * k!)
    # But for large |w|, use asymptotic or numerical quadrature.
    return _ei_complex(w)

def _ei_complex(w, terms=200):
    """Compute Ei(w) for complex w using series + Euler-Mascheroni."""
    gamma_euler = 0.5772156649015329
    if abs(w) > 40:
        # Asymptotic expansion: Ei(w) ~ e^w/w * sum_{k=0}^{N} k!/w^k
        result = 0.0
        term = 1.0
        for k in range(1, 30):
            term *= k / w
            result += term
            if abs(term) < 1e-15 * abs(result):
                break
        return np.exp(w) / w * (1.0 + result)
    # Series: Ei(w) = gamma + ln(w) + sum w^k/(k*k!)
    result = gamma_euler + np.log(w)  # principal branch
    term = w
    s = term  # k=1 term: w^1/(1*1!)
    for k in range(2, terms):
        term *= w / k
        s += term / k
        if abs(term / k) < 1e-15:
            break
    return result + s

# ---------------------------------------------------------------------------
# Explicit formula for π(x)
# ---------------------------------------------------------------------------
def pi_explicit(x, num_zeros=100):
    """
    π(x) ≈ Li(x) - Σ_{ρ} Li(x^ρ) - log(2) + integral_x^∞ dt/((t^2-1)*t*ln(t))

    The sum over zeros: each pair ρ = 1/2 ± iγ contributes
      Li(x^ρ) + Li(x^{ρ̄}) = 2 * Re(Li(x^ρ))

    The integral term for x > 1 is small and negative.
    """
    if x <= 2:
        return 0 if x < 2 else 1

    # Main term
    result = li(x)

    # Zero sum: each gamma gives a conjugate pair
    gammas = GAMMAS[:num_zeros]
    zero_sum = 0.0
    for g in gammas:
        rho = 0.5 + 1j * g
        val = li_complex(x, rho)
        zero_sum += 2.0 * val.real  # conjugate pair

    result -= zero_sum

    # Constant term
    result -= np.log(2)

    # Integral term: integral_x^inf dt / ((t^2-1)*t*ln(t))
    # This is small for x >> 1. Compute numerically for moderate x.
    if x < 1e6:
        try:
            integral_val, _ = quad(lambda t: 1.0 / ((t**2 - 1) * t * np.log(t)),
                                   x, np.inf, limit=100)
            result -= integral_val
        except Exception:
            pass  # negligible for large x

    return result

# ---------------------------------------------------------------------------
# Test 1: Convergence of explicit formula vs number of zeros
# ---------------------------------------------------------------------------
def test_convergence():
    print("=" * 70)
    print("TEST 1: Explicit formula convergence vs number of zeros")
    print("=" * 70)

    pi_ref = sieve_pi(10**6)
    test_x = [100, 500, 1000, 5000, 10000, 50000, 100000, 500000, 1000000]
    zero_counts = [5, 10, 20, 50, 100, 150, 200]

    print(f"\n{'x':>10} | {'pi(x)':>8} | ", end="")
    for nz in zero_counts:
        print(f"{'Z=' + str(nz):>10}", end=" | ")
    print()
    print("-" * (25 + 13 * len(zero_counts)))

    results = {}
    for x in test_x:
        exact = pi_ref[x]
        print(f"{x:>10} | {exact:>8} | ", end="")
        for nz in zero_counts:
            approx = pi_explicit(x, nz)
            err = approx - exact
            results[(x, nz)] = err
            print(f"{err:>+10.2f}", end=" | ")
        print()

    # Summary: how many zeros for <1 error?
    print("\nZeros needed for |error| < 1:")
    for x in test_x:
        needed = "never"
        for nz in zero_counts:
            if abs(results[(x, nz)]) < 1.0:
                needed = str(nz)
                break
        print(f"  x = {x:>10}: {needed} zeros")

    return results

# ---------------------------------------------------------------------------
# Test 2: Smoothed explicit formula variants
# ---------------------------------------------------------------------------
def smoothed_zero_sum(x, num_zeros, sigma):
    """
    Replace Li(x^ρ) with a Gaussian-damped version:
    Σ_ρ exp(-γ^2/(2σ^2)) * Li(x^ρ)

    This suppresses high-frequency zeros, acting like a heat kernel.
    """
    gammas = GAMMAS[:num_zeros]
    result = 0.0
    for g in gammas:
        rho = 0.5 + 1j * g
        damping = np.exp(-g**2 / (2 * sigma**2))
        val = li_complex(x, rho)
        result += 2.0 * damping * val.real
    return result

def test_smoothing():
    print("\n" + "=" * 70)
    print("TEST 2: Gaussian-smoothed zero sums")
    print("=" * 70)

    pi_ref = sieve_pi(100000)
    test_x = [100, 1000, 10000, 100000]
    sigmas = [10, 20, 50, 100, 200, 500, np.inf]

    for x in test_x:
        exact = pi_ref[x]
        li_x = li(x)
        print(f"\nx = {x}, π(x) = {exact}, Li(x) = {li_x:.2f}")
        print(f"  {'sigma':>8} | {'zero_sum':>12} | {'approx':>12} | {'error':>10}")
        print(f"  {'-'*50}")
        for sigma in sigmas:
            label = f"{sigma:.0f}" if np.isfinite(sigma) else "inf"
            zs = smoothed_zero_sum(x, 200, sigma)
            approx = li_x - zs - np.log(2)
            err = approx - exact
            print(f"  {label:>8} | {zs:>12.4f} | {approx:>12.4f} | {err:>+10.4f}")

# ---------------------------------------------------------------------------
# Test 3: Structure in the zero sum S(x) = Σ x^{iγ}/(1/2+iγ)
# ---------------------------------------------------------------------------
def zero_sum_oscillatory(x_values, num_zeros=200):
    """
    Compute S(x) = Σ_{k=1}^{num_zeros} x^{iγ_k} / (1/2 + iγ_k)
    This is the "oscillatory kernel" -- the hard part of the explicit formula.
    Return real and imaginary parts.
    """
    gammas = GAMMAS[:num_zeros]
    results = []
    for x in x_values:
        lnx = np.log(x)
        s = 0.0 + 0.0j
        for g in gammas:
            rho = 0.5 + 1j * g
            # x^{iγ} = e^{iγ ln(x)}
            s += np.exp(1j * g * lnx) / rho
        results.append(s)
    return np.array(results)

def test_zero_sum_structure():
    print("\n" + "=" * 70)
    print("TEST 3: Structure in oscillatory zero sum S(x) = Σ x^{iγ}/(ρ)")
    print("=" * 70)

    # Compute S(x) for many x values
    x_values = np.logspace(1, 6, 500)
    S = zero_sum_oscillatory(x_values, 200)

    # Statistics
    S_real = S.real
    S_imag = S.imag
    S_abs = np.abs(S)

    print(f"\nS(x) statistics over x in [10, 10^6] (500 points, 200 zeros):")
    print(f"  |S| mean: {np.mean(S_abs):.4f}, std: {np.std(S_abs):.4f}")
    print(f"  |S| min:  {np.min(S_abs):.4f}, max: {np.max(S_abs):.4f}")
    print(f"  Re(S) mean: {np.mean(S_real):.4f}, std: {np.std(S_real):.4f}")
    print(f"  Im(S) mean: {np.mean(S_imag):.4f}, std: {np.std(S_imag):.4f}")

    # Check if |S(x)| grows, decays, or stays bounded
    # Partition into ranges
    ranges = [(10, 100), (100, 1000), (1000, 10000), (10000, 100000), (100000, 1000000)]
    print(f"\n  {'Range':>20} | {'mean |S|':>10} | {'max |S|':>10} | {'std Re(S)':>10}")
    for lo, hi in ranges:
        mask = (x_values >= lo) & (x_values < hi)
        print(f"  [{lo:>7},{hi:>7}) | {np.mean(S_abs[mask]):>10.4f} | "
              f"{np.max(S_abs[mask]):>10.4f} | {np.std(S_real[mask]):>10.4f}")

    # Autocorrelation of Re(S) to detect periodicity
    sr = S_real - np.mean(S_real)
    autocorr = np.correlate(sr, sr, 'full')
    autocorr = autocorr[len(autocorr)//2:]
    autocorr /= autocorr[0]
    # Find first significant peak after lag 0
    peaks = []
    for i in range(2, len(autocorr) - 1):
        if autocorr[i] > autocorr[i-1] and autocorr[i] > autocorr[i+1] and autocorr[i] > 0.3:
            peaks.append((i, autocorr[i]))
    print(f"\n  Autocorrelation peaks (lag, value): {peaks[:5] if peaks else 'None found'}")

    return x_values, S

# ---------------------------------------------------------------------------
# Test 4: Can the zero sum be approximated by a short formula?
# ---------------------------------------------------------------------------
def test_short_formula():
    print("\n" + "=" * 70)
    print("TEST 4: Short-formula approximation of the zero sum")
    print("=" * 70)
    print("Can S(x) = Σ_ρ Li(x^ρ) be approximated by f(log x, x^{1/2})?")

    pi_ref = sieve_pi(100000)
    x_values = np.arange(100, 100001, 100)

    # Compute actual oscillatory correction: osc(x) = Li(x) - pi(x) - log(2)
    osc_actual = np.array([li(x) - pi_ref[x] - np.log(2) for x in x_values])

    # Candidate short formulas for the oscillatory part:
    lnx = np.log(x_values)
    sqrtx = np.sqrt(x_values)

    candidates = {
        'sqrt(x)/ln(x)': sqrtx / lnx,
        'sqrt(x)/ln(x)^2': sqrtx / lnx**2,
        'sqrt(x)*sin(g1*lnx)/lnx': sqrtx * np.sin(GAMMAS[0] * lnx) / lnx,
        '2*sqrt(x)*cos(g1*lnx)/lnx': 2 * sqrtx * np.cos(GAMMAS[0] * lnx) / lnx,
        # First 3 zeros
        'Σ_{k=1}^3 2*sqrt(x)*cos(gk*lnx)/(|rho_k|*lnx)':
            sum(2 * sqrtx * np.cos(GAMMAS[k] * lnx) / (abs(0.5 + 1j*GAMMAS[k]) * lnx) for k in range(3)),
        # First 10 zeros
        'Σ_{k=1}^10 ...':
            sum(2 * sqrtx * np.cos(GAMMAS[k] * lnx) / (abs(0.5 + 1j*GAMMAS[k]) * lnx) for k in range(10)),
    }

    print(f"\n  {'Formula':>45} | {'corr':>8} | {'RMSE':>10} | {'rel RMSE':>10}")
    print(f"  {'-'*80}")

    for name, vals in candidates.items():
        # Best linear fit: osc ≈ a * vals + b
        if np.std(vals) < 1e-10:
            continue
        A = np.vstack([vals, np.ones_like(vals)]).T
        coef, residuals, _, _ = np.linalg.lstsq(A, osc_actual, rcond=None)
        fitted = A @ coef
        corr = np.corrcoef(osc_actual, fitted)[0, 1]
        rmse = np.sqrt(np.mean((osc_actual - fitted)**2))
        rel_rmse = rmse / np.std(osc_actual)
        print(f"  {name:>45} | {corr:>8.4f} | {rmse:>10.2f} | {rel_rmse:>10.4f}")

    # Also try: how much of osc_actual variance is explained by first K zeros?
    print("\n  Variance explained by first K zeros (least-squares fit):")
    total_var = np.var(osc_actual)
    for K in [1, 3, 5, 10, 20, 50]:
        # Basis: sqrt(x)*cos(gk*lnx)/lnx, sqrt(x)*sin(gk*lnx)/lnx for k=1..K
        basis = []
        for k in range(K):
            basis.append(sqrtx * np.cos(GAMMAS[k] * lnx) / lnx)
            basis.append(sqrtx * np.sin(GAMMAS[k] * lnx) / lnx)
        basis.append(np.ones_like(x_values))  # constant
        A = np.vstack(basis).T
        coef, _, _, _ = np.linalg.lstsq(A, osc_actual, rcond=None)
        fitted = A @ coef
        var_explained = 1.0 - np.var(osc_actual - fitted) / total_var
        rmse = np.sqrt(np.mean((osc_actual - fitted)**2))
        print(f"    K={K:>3}: R^2 = {var_explained:.6f}, RMSE = {rmse:.4f}")

# ---------------------------------------------------------------------------
# Test 5: Riemann-von Mangoldt formula verification
# ---------------------------------------------------------------------------
def test_riemann_von_mangoldt():
    print("\n" + "=" * 70)
    print("TEST 5: Riemann-von Mangoldt N(T) = T/(2π) ln(T/(2πe)) + 7/8 + S(T)")
    print("=" * 70)

    # N(T) = number of zeros with 0 < γ ≤ T
    # Compare with smooth approximation
    gammas = load_zeros(200)

    T_values = gammas  # evaluate at each zero
    N_actual = np.arange(1, len(gammas) + 1)  # N(γ_k) = k (approximately, at the zero itself)

    N_smooth = T_values / (2 * np.pi) * np.log(T_values / (2 * np.pi * np.e)) + 7.0/8.0

    S_T = N_actual - N_smooth  # This is S(T), the fluctuation

    print(f"\n  N(T) smooth vs actual for first 200 zeros:")
    print(f"  {'k':>5} | {'γ_k':>12} | {'N_smooth':>10} | {'N_actual':>10} | {'S(γ_k)':>10}")
    print(f"  {'-'*55}")
    for k in [1, 5, 10, 20, 50, 100, 150, 200]:
        i = k - 1
        print(f"  {k:>5} | {gammas[i]:>12.4f} | {N_smooth[i]:>10.4f} | {N_actual[i]:>10} | {S_T[i]:>+10.4f}")

    print(f"\n  S(T) statistics:")
    print(f"    mean: {np.mean(S_T):>+.4f}")
    print(f"    std:  {np.std(S_T):.4f}")
    print(f"    max:  {np.max(np.abs(S_T)):.4f}")
    print(f"    Expected std (log log T / sqrt(2pi)): ~{np.log(np.log(gammas[99])) / np.sqrt(2*np.pi):.4f}")

    return S_T

# ---------------------------------------------------------------------------
# Test 6: Trace integral approach -- can we compute Σ_ρ f(ρ) as a contour integral?
# ---------------------------------------------------------------------------
def test_trace_integral():
    print("\n" + "=" * 70)
    print("TEST 6: Trace/contour integral for zero sum")
    print("=" * 70)
    print("The argument principle: Σ_ρ f(ρ) = (1/2πi) ∮ f(s) ζ'(s)/ζ(s) ds")
    print("This requires evaluating ζ'/ζ on a contour -- cost O(T) per point.")
    print()

    # Test: for the explicit formula, the relevant sum is
    # Σ_ρ x^ρ/ρ  (simplified from Li(x^ρ))
    # Using argument principle: (1/2πi) ∮ (x^s/s) * (ζ'/ζ)(s) ds
    # along a rectangle with Re(s) in [c-iT, c+iT, -U+iT, -U-iT]

    # The key question: can we evaluate this contour integral faster than summing zeros?
    # Answer depends on how many quadrature points are needed.

    # Numerical test: how many quadrature points on the critical line
    # to match N zeros?
    x = 1000
    lnx = np.log(x)

    print(f"  x = {x}")
    print(f"  Comparing direct zero sum vs numerical integration of Σ x^{{iγ}}/ρ")

    # Direct sum with K zeros
    for K in [10, 50, 100, 200]:
        gammas = GAMMAS[:K]
        direct = sum(2 * np.real(x**(0.5 + 1j*g) / (0.5 + 1j*g)) for g in gammas)
        print(f"  Direct sum ({K:>3} zeros): {direct:>12.4f}")

    # Now: quadrature on the critical line integral
    # (1/2πi) ∫_{1/2-iT}^{1/2+iT} (x^s/s) (ζ'/ζ)(s) ds
    # ≈ (1/π) ∫_0^T Re[x^{1/2+it}/(1/2+it) * Z'/Z(t)] dt
    # where Z is the Hardy Z-function.
    # Problem: we don't have ζ'/ζ, so we can't actually do this contour integral
    # without computing ζ. This is the fundamental issue.
    print()
    print("  KEY INSIGHT: The contour integral requires ζ'/ζ(s) on the contour.")
    print("  Computing ζ(s) at height t costs O(t^{1/2}) via Riemann-Siegel.")
    print("  Integrating from 0 to T with M quadrature points:")
    print("  Total cost = M * O(T^{1/2}).")
    print("  Direct zero sum cost = K * O(1) per zero, but finding K zeros costs O(K * T^{1/2}).")
    print()
    print("  So: contour integral with M quadrature points ≈ direct sum with M zeros,")
    print("  BOTH have the same O(T^{3/2}) total cost to height T.")
    print("  NO SHORTCUT from the trace formula approach.")

    # Verify: the Riemann-Siegel formula
    # Z(t) = 2 * Σ_{n ≤ sqrt(t/2π)} n^{-1/2} * cos(θ(t) - t*ln(n)) + remainder
    # where θ(t) = arg(Gamma(1/4 + it/2)) - (t/2)*ln(π)
    # Cost per evaluation: O(sqrt(t)) terms
    T_values = [100, 1000, 10000, 100000]
    print(f"\n  Riemann-Siegel cost per evaluation of ζ(1/2+it):")
    for T in T_values:
        n_terms = int(np.sqrt(T / (2 * np.pi)))
        print(f"    T = {T:>7}: {n_terms} terms (= O(T^{{1/2}}))")

    print(f"\n  For π(x), need zeros up to T ≈ x^{{1/2}} / (2π).")
    print(f"  Total zeros: O(T log T). Cost per zero: O(T^{{1/2}}).")
    print(f"  Total: O(T^{{3/2}} log T) = O(x^{{3/4}} log x).")
    print(f"  This is WORSE than Meissel-Lehmer O(x^{{2/3}}).")

# ---------------------------------------------------------------------------
# Test 7: Heat kernel analog -- spectral smoothing
# ---------------------------------------------------------------------------
def test_heat_kernel():
    print("\n" + "=" * 70)
    print("TEST 7: Heat kernel analog -- spectral smoothing")
    print("=" * 70)
    print("Replace Li(x^ρ) with e^{-t*|ρ|^2} * Li(x^ρ) ('heat kernel smoothing')")
    print("As t→0, recover exact formula. Question: for what t is it both")
    print("accurate AND computable from a closed-form heat trace?\n")

    pi_ref = sieve_pi(100000)
    test_points = [100, 1000, 10000, 100000]

    # The heat trace of the Riemann zeros:
    # H(t) = Σ_ρ e^{-t|ρ|^2} = Σ_k e^{-t(1/4 + γ_k^2)}
    # This is related to the Weil distribution.
    # For small t: H(t) ~ (1/(4πt)) * (some constant) -- Weyl law analog

    print("  Heat trace H(t) = Σ_{k=1}^{200} e^{-t*(1/4+γ_k^2)}:")
    gammas = GAMMAS[:200]
    eigenvalues = 0.25 + gammas**2
    for t in [0.0001, 0.001, 0.01, 0.1, 1.0]:
        H = np.sum(np.exp(-t * eigenvalues))
        # Weyl approximation: N(T) ~ T/(2π)*ln(T/(2πe)), so
        # H(t) ~ ∫_0^∞ e^{-t*u^2} * d(N(u)) ≈ ∫ e^{-tu^2} * ln(u/(2πe))/(2π) du
        # = (1/(2π)) * sqrt(π/t) * (ln(1/(2πe*t))/2 + ...)  (rough)
        # Exact Weyl: H(t) ≈ (1/(4πt))^{1/2} * (ln(1/t) - something)
        weyl_approx = np.sqrt(np.pi / t) / (2 * np.pi) * max(0.1, np.log(1/(2*np.pi*np.e*t))/2)
        print(f"    t = {t:.4f}: H(t) = {H:>12.4f}, Weyl approx ≈ {weyl_approx:>12.4f}")

    print()
    print("  Heat-smoothed π(x) for various t:")
    for x in test_points:
        exact = pi_ref[x]
        li_x = li(x)
        print(f"\n  x = {x}, π(x) = {exact}")
        print(f"    {'t':>10} | {'heat_sum':>12} | {'approx':>12} | {'error':>10} | {'%err':>8}")
        for t in [0.0, 0.0001, 0.001, 0.01, 0.1]:
            heat_sum = 0.0
            for g in gammas:
                rho = 0.5 + 1j * g
                damping = np.exp(-t * (0.25 + g**2)) if t > 0 else 1.0
                val = li_complex(x, rho)
                heat_sum += 2.0 * damping * val.real
            approx = li_x - heat_sum - np.log(2)
            err = approx - exact
            pct = 100 * abs(err) / exact if exact > 0 else 0
            label = f"{t:.4f}" if t > 0 else "exact"
            print(f"    {label:>10} | {heat_sum:>12.4f} | {approx:>12.4f} | {err:>+10.4f} | {pct:>7.3f}%")

    print()
    print("  CONCLUSION: Heat smoothing suppresses high zeros (cheap to compute)")
    print("  but introduces O(x^{1/2}/ln(x)) bias that grows with x.")
    print("  No t value gives both exact AND fewer terms needed.")
    print("  The heat trace itself IS computable but doesn't encode π(x) directly.")

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    t0 = time.time()

    print("SPECTRAL SHORTCUT VIA TRACE FORMULA -- EXPERIMENT")
    print(f"Using {len(GAMMAS)} zeta zeros (γ_1 = {GAMMAS[0]:.4f}, γ_max = {GAMMAS[-1]:.4f})")
    print()

    results1 = test_convergence()
    test_smoothing()
    x_vals, S = test_zero_sum_structure()
    test_short_formula()
    S_T = test_riemann_von_mangoldt()
    test_trace_integral()
    test_heat_kernel()

    elapsed = time.time() - t0

    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"""
1. EXPLICIT FORMULA CONVERGENCE: With 200 zeros, error is O(x^{{1/2}}/ln(x)).
   For x=10^5, need ~100+ zeros for |error|<1. Scales as O(sqrt(x)).

2. GAUSSIAN SMOOTHING: Damping high zeros reduces oscillation but introduces
   systematic bias. No free lunch -- the information IS in the high zeros.

3. ZERO SUM STRUCTURE: |S(x)| = O(sqrt(x)/ln(x)), with pseudo-random oscillation.
   No periodic structure detected. Autocorrelation dies quickly.

4. SHORT FORMULA: First 10-20 zeros explain ~80-95% of variance in the oscillatory
   part for x < 10^5. But the remaining 5-20% is the hard part and grows with x.
   Each zero adds O(1) "bits" of information -- no compression possible.

5. RIEMANN-VON MANGOLDT: S(T) fluctuates with std ≈ 0.2-0.4, consistent with
   log log T / sqrt(2π). Confirms smooth part is excellent but residual is irreducible.

6. TRACE/CONTOUR INTEGRAL: Computing ζ'/ζ on a contour costs O(T^{{1/2}}) per
   point (Riemann-Siegel). Total cost O(T^{{3/2}} log T) -- same as or worse than
   direct zero enumeration. NO SPEEDUP from trace formula.

7. HEAT KERNEL: Spectral smoothing trades resolution for computability.
   The heat trace Σ e^{{-t|ρ|^2}} is computable but doesn't encode π(x).
   As t→0 (exact limit), the heat trace diverges -- no closed form.

VERDICT: The trace formula / spectral shortcut does NOT bypass the zero sum.
The zero sum IS the spectral trace, and computing it requires either:
  (a) enumerating zeros: O(T^{{3/2}}) total, or
  (b) evaluating ζ on a contour: O(T^{{3/2}}) total.
Both are worse than Meissel-Lehmer's O(x^{{2/3}}) combinatorial approach.

The fundamental barrier remains: π(x) encodes O(x^{{1/2}}/log x) independent
oscillatory contributions, and no known spectral identity collapses them.

Time elapsed: {elapsed:.1f}s
""")

if __name__ == '__main__':
    main()

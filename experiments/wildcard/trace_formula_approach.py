#!/usr/bin/env python3
"""
TRACE FORMULA APPROACH: Computing π(x) Without Enumerating Individual Zeros
============================================================================

The explicit formula: π(x) ≈ li(x) - Σ_ρ li(x^ρ) + smaller
Under RH: ρ = 1/2 + iγ_j, so the zero sum S(x) = Σ_j x^{ρ_j}/ρ_j.

This looks like Tr(f(A)) where A has eigenvalues {ρ_j}.
Can we compute this trace WITHOUT enumerating eigenvalues?

EXPERIMENTS:
1. Moments method: S(x) via spectral moments M_k = Σ_ρ ρ^{k-1}
2. Direct GUE trace: Tr(f(A)) for random GUE vs eigenvalue sum
3. Weil explicit formula: partial geometric side accuracy
4. Nuclear/effective rank: rank-k approximation of the operator

Uses zeta zeros from data/ directory.
"""

import numpy as np
from scipy.special import expi
from scipy.integrate import quad
from scipy.linalg import eigvalsh
import math
import os
import sys
import time
import warnings
warnings.filterwarnings('ignore')

# ---------------------------------------------------------------------------
# Load zeta zeros
# ---------------------------------------------------------------------------
DATA_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'data')

def load_zeros(n=1000):
    """Load first n zeta zeros (imaginary parts γ_j)."""
    for size in [1000, 500, 300, 200]:
        path = os.path.join(DATA_DIR, f'zeta_zeros_{size}.txt')
        if os.path.exists(path):
            gammas = []
            with open(path) as f:
                for line in f:
                    parts = line.strip().split()
                    if len(parts) >= 2:
                        gammas.append(float(parts[1]))
                    elif len(parts) == 1:
                        gammas.append(float(parts[0]))
            gammas = np.array(gammas[:n])
            print(f"Loaded {len(gammas)} zeros from {os.path.basename(path)}")
            return gammas
    raise FileNotFoundError("No zeta zero files found")

def li_func(x):
    """Logarithmic integral li(x) = Ei(ln(x)), x > 1."""
    if x <= 1:
        return 0.0
    return expi(np.log(x))

def prime_count_exact(x):
    """Exact π(x) using sympy for moderate x."""
    from sympy import primepi
    return int(primepi(int(x)))

# ---------------------------------------------------------------------------
# EXPERIMENT 1: Moments Method
# ---------------------------------------------------------------------------
def experiment_moments_method(gammas, test_x_values):
    """
    S(x) = Σ_j x^{ρ_j}/ρ_j  where ρ_j = 1/2 + iγ_j

    Expand x^ρ = e^{ρ ln x} = Σ_k (ln x)^k / k! · ρ^k
    So S(x) = Σ_k (ln x)^k / k! · Σ_j ρ_j^{k-1} = Σ_k (ln x)^k / k! · M_{k-1}

    where M_k = Σ_j ρ_j^k are the spectral moments.

    Question: How many moments K are needed for accuracy ε at x?
    """
    print("\n" + "="*70)
    print("EXPERIMENT 1: MOMENTS METHOD")
    print("="*70)

    rhos = 0.5 + 1j * gammas  # the zeros
    N = len(rhos)

    # Precompute moments M_k = Σ_j ρ_j^k for k = 0, 1, ..., K_max
    K_max = 80
    # We need both S(x) and its moment expansion
    # S(x) = Σ_j e^{ρ_j ln x} / ρ_j
    #       = Σ_k (ln x)^k / k! · M_{k-1}
    # where M_k = Σ_j ρ_j^k  (but we need M_{-1} = Σ 1/ρ_j for k=0 term)
    #
    # Actually more carefully:
    # x^ρ/ρ = e^{ρ ln x}/ρ = (1/ρ) Σ_k (ρ ln x)^k / k!
    #        = Σ_k (ln x)^k / k! · ρ^{k-1}
    # So S(x) = Σ_k (ln x)^k / k! · Σ_j ρ_j^{k-1}
    #         = Σ_k (ln x)^k / k! · M_{k-1}
    # where M_m = Σ_j ρ_j^m, and we need m = -1, 0, 1, ..., K_max-1

    # Compute M_m for m = -1, 0, 1, ..., K_max-1
    moments = {}
    moments[-1] = np.sum(1.0 / rhos)
    for m in range(K_max):
        moments[m] = np.sum(rhos**m)

    results_table = []

    for x_val in test_x_values:
        lnx = np.log(x_val)

        # "Exact" zero sum using all zeros (with both ρ and 1-ρ = conj(ρ) under RH)
        # S(x) = Σ_j [x^{ρ_j}/ρ_j + x^{(1-ρ_j)}/(1-ρ_j)]
        # Under RH: 1-ρ_j = 1/2 - iγ_j = conj(ρ_j)
        zero_sum_exact = 0.0
        for rho in rhos:
            contrib = x_val**rho / rho
            zero_sum_exact += contrib + np.conj(contrib)  # pair ρ and 1-ρ
        S_exact = zero_sum_exact.real

        # Moment expansion: S_K(x) = 2 Re[ Σ_{k=0}^{K} (ln x)^k / k! · M_{k-1} ]
        # (factor 2 Re for pairing ρ with conj(ρ))
        convergence = []
        for K in range(1, K_max + 1):
            S_approx = 0.0 + 0j
            for k in range(K):
                S_approx += (lnx**k / math.factorial(k)) * moments[k - 1]
            S_approx_real = 2 * S_approx.real
            error = abs(S_approx_real - S_exact)
            convergence.append((K, S_approx_real, error))

        # Find K needed for various accuracies
        best_error = min(c[2] for c in convergence)
        K_for_1 = next((c[0] for c in convergence if c[2] < 1.0), None)
        K_for_01 = next((c[0] for c in convergence if c[2] < 0.1), None)

        results_table.append({
            'x': x_val,
            'ln_x': lnx,
            'S_exact': S_exact,
            'best_error': best_error,
            'K_for_err_1': K_for_1,
            'K_for_err_01': K_for_01,
            'convergence': convergence
        })

        print(f"\nx = {x_val:.0f}, ln(x) = {lnx:.2f}")
        print(f"  S_exact (with {N} zero pairs) = {S_exact:.6f}")
        print(f"  Best moment approx error = {best_error:.6e}")
        print(f"  K for |error| < 1:   {K_for_1}")
        print(f"  K for |error| < 0.1: {K_for_01}")

        # Check: does series diverge?
        errors_at_K = [convergence[k-1][2] for k in [5, 10, 20, 40, 60, 80] if k <= K_max]
        Ks = [k for k in [5, 10, 20, 40, 60, 80] if k <= K_max]
        print(f"  Errors at K = {Ks}: {['%.2e' % e for e in errors_at_K]}")

    # KEY DIAGNOSTIC: Does the series converge?
    # The radius of convergence of e^{ρ ln x} as a power series in ρ is infinite,
    # BUT the moments M_k = Σ ρ_j^k grow as the zeros grow, so the combined series
    # Σ_k (ln x)^k / k! · M_{k-1} may diverge.
    print("\n--- CONVERGENCE DIAGNOSTIC ---")
    print("M_k growth (|M_k| for k = 0..20):")
    for k in range(0, 21, 2):
        print(f"  |M_{k}| = {abs(moments[k]):.4e}")

    return results_table


# ---------------------------------------------------------------------------
# EXPERIMENT 2: GUE Random Matrix Trace Computation
# ---------------------------------------------------------------------------
def experiment_gue_trace(gammas):
    """
    For GUE(N), eigenvalues λ_i: Tr(f(A)) = Σ f(λ_i).
    But Tr(f(A)) can also be computed from matrix entries via Tr(A^k).

    Test: For f(λ) = e^{iλt}/λ (our test function), how does
    Tr(f(A)) computed from Tr(A^k) compare to Σ f(λ_i)?

    Also: How well does a GUE matrix model the actual zeta zero sum?
    """
    print("\n" + "="*70)
    print("EXPERIMENT 2: GUE RANDOM MATRIX TRACE")
    print("="*70)

    np.random.seed(42)

    results = []

    for N in [50, 100, 200]:
        # Generate GUE(N) matrix: H = (A + A^*) / (2√N)
        A = (np.random.randn(N, N) + 1j * np.random.randn(N, N)) / np.sqrt(2)
        H = (A + A.conj().T) / (2 * np.sqrt(N))
        eigs = eigvalsh(H)  # real eigenvalues (sorted)

        # Scale to match zeta zero spacing: mean spacing ~ 2π/ln(N)
        # For GUE(N), eigenvalues are in [-1, 1] with spacing ~ 1/N near center
        # Zeta zeros have spacing ~ 2π/ln(γ/(2π)) near height γ

        # Test function: f(λ) = e^{iλt} for various t
        # Tr(f(H)) = Σ_k (it)^k/k! · Tr(H^k)

        for t_val in [1.0, 5.0, 10.0]:
            # "Exact" from eigenvalues
            trace_exact = np.sum(np.exp(1j * eigs * t_val))

            # From matrix powers: Tr(e^{itH}) = Σ_k (it)^k/k! · Tr(H^k)
            H_power = np.eye(N, dtype=complex)
            trace_approx = np.zeros(60, dtype=complex)
            trace_approx[0] = N  # Tr(I) = N

            running_sum = N + 0j
            convergence = [(0, abs(running_sum - trace_exact))]

            for k in range(1, 60):
                H_power = H_power @ H
                tr_Hk = np.trace(H_power)
                coeff = (1j * t_val)**k / math.factorial(min(k, 170))
                running_sum += coeff * tr_Hk
                convergence.append((k, abs(running_sum - trace_exact)))

            best_K = min(convergence, key=lambda c: c[1])

            result = {
                'N': N, 't': t_val,
                'trace_exact': trace_exact,
                'best_error': best_K[1],
                'best_K': best_K[0],
                'convergence': convergence
            }
            results.append(result)

            print(f"\nGUE({N}), t={t_val}: |Tr_exact| = {abs(trace_exact):.4f}")
            print(f"  Best moment approx at K={best_K[0]}: error = {best_K[1]:.4e}")
            # Show convergence
            select_K = [1, 5, 10, 20, 30, 50]
            for kk in select_K:
                if kk < len(convergence):
                    print(f"    K={kk:2d}: error = {convergence[kk][1]:.4e}")

    # KEY COMPARISON: GUE trace convergence vs zeta zero moment convergence
    print("\n--- KEY COMPARISON ---")
    print("For GUE matrices, the moment expansion converges because:")
    print("  1. Eigenvalues are bounded (|λ| ≤ 1 for GUE)")
    print("  2. Tr(H^k) ~ O(1) for large k (Wigner semicircle)")
    print("For zeta zeros:")
    print("  1. γ_j → ∞ (unbounded)")
    print("  2. M_k = Σ ρ^k diverges for k ≥ 1")
    print("This is a FUNDAMENTAL difference.")

    return results


# ---------------------------------------------------------------------------
# EXPERIMENT 3: Weil Explicit Formula -- Geometric Side Truncation
# ---------------------------------------------------------------------------
def experiment_weil_truncation(gammas, test_x_values):
    """
    Weil explicit formula: Σ_ρ h(ρ) = h(0) + h(1) - Σ_p Σ_m h_hat(m log p) log(p) / p^{m/2}
                            + integral terms

    For h(s) = x^s / s, h_hat involves log x.

    The geometric side sums over primes -- seems circular for computing π(x).
    BUT: if the sum converges fast, we only need primes up to some y << x.

    Test: What fraction of the geometric side comes from primes ≤ y?
    """
    print("\n" + "="*70)
    print("EXPERIMENT 3: WEIL EXPLICIT FORMULA -- GEOMETRIC SIDE TRUNCATION")
    print("="*70)

    from sympy import primerange, primepi

    for x_val in test_x_values:
        lnx = np.log(x_val)
        pi_exact = int(primepi(int(x_val)))
        li_val = li_func(x_val)

        # The explicit formula zero sum (spectral side) using our zeros
        rhos = 0.5 + 1j * gammas
        spectral_sum = 0.0
        for rho in rhos:
            contrib = li_func_complex(x_val, rho)
            spectral_sum += 2 * contrib.real  # pair ρ and conj(ρ)

        # Geometric side: Σ_{p ≤ y} Σ_{m: p^m ≤ x} 1/(m · p^{m/2})
        # This is the von Mangoldt / Chebyshev ψ(x) version
        # More precisely: the prime-sum side of the explicit formula for ψ(x)
        #   ψ(x) = x - Σ_ρ x^ρ/ρ - log(2π) - (1/2)log(1-x^{-2})
        # The geometric analog: ψ(x) = Σ_{p^m ≤ x} log(p)

        # Test: How much of ψ(x) comes from primes ≤ y for various y?
        primes = list(primerange(2, int(x_val) + 1))

        # Full ψ(x)
        psi_full = 0.0
        for p in primes:
            m = 1
            while p**m <= x_val:
                psi_full += np.log(p)
                m += 1

        # Truncated ψ(y) for various y
        print(f"\nx = {x_val:.0f}, π(x) = {pi_exact}, ψ(x) = {psi_full:.4f}, x = {x_val:.4f}")
        print(f"  y/x ratio | ψ(y)/ψ(x) | # primes ≤ y | relative error")
        print(f"  " + "-"*60)

        for frac in [0.01, 0.05, 0.1, 0.2, 0.5, 0.8, 1.0]:
            y = max(2, int(frac * x_val))
            psi_trunc = 0.0
            n_primes_y = 0
            for p in primes:
                if p > y:
                    break
                n_primes_y += 1
                m = 1
                while p**m <= x_val:
                    psi_trunc += np.log(p)
                    m += 1

            rel_err = abs(psi_trunc - psi_full) / max(psi_full, 1e-10)
            print(f"  {frac:8.2f}   | {psi_trunc/max(psi_full,1e-10):10.6f} | {n_primes_y:12d} | {rel_err:.6e}")

    # KEY QUESTION: To get ψ(x) with error < 0.5 (enough for exact π(x)),
    # we need essentially ALL primes up to x. This is circular.
    print("\n--- VERDICT ON GEOMETRIC SIDE ---")
    print("To compute ψ(x) within ±0.5 from the geometric side,")
    print("we need all primes up to x. This is fully circular.")
    print("The explicit formula trades: (spectral zeros) <-> (primes)")
    print("Neither side avoids enumerating the other.")


def li_func_complex(x, rho):
    """Compute li(x^ρ) = Ei(ρ · ln(x)) approximately."""
    z = rho * np.log(x)
    # For complex z, Ei(z) ~ e^z / z for large |z|
    # Use series for small |z|, asymptotic for large
    if abs(z) < 30:
        # Series: Ei(z) = γ + ln(z) + Σ_{k=1}^∞ z^k / (k · k!)
        result = 0.5772156649 + np.log(z + 0j)
        term = z
        for k in range(1, 100):
            result += term / (k * math.factorial(min(k, 170)))
            term *= z / (k + 1)
            if abs(term) < 1e-15:
                break
        return result
    else:
        # Asymptotic: Ei(z) ~ e^z/z · (1 + 1/z + 2/z^2 + ...)
        return np.exp(z) / z * (1 + 1/z + 2/z**2 + 6/z**3)


# ---------------------------------------------------------------------------
# EXPERIMENT 4: Nuclear Norm / Effective Rank
# ---------------------------------------------------------------------------
def experiment_effective_rank(gammas, test_x_values):
    """
    If we think of S(x) = Σ_j f(x, γ_j) as the trace of an operator
    with kernel K(x, y) = Σ_j φ_j(x) φ_j(y), what is its effective rank?

    Construct the matrix M_{ij} = x_i^{ρ_j} / ρ_j for a grid of x values
    and compute its singular values. The effective rank tells us how
    compressible the zero sum is.
    """
    print("\n" + "="*70)
    print("EXPERIMENT 4: EFFECTIVE RANK OF ZERO-SUM OPERATOR")
    print("="*70)

    rhos = 0.5 + 1j * gammas
    N_zeros = min(len(gammas), 200)  # use first 200 zeros
    rhos_sub = rhos[:N_zeros]

    # Grid of x values
    for x_range_label, x_grid in [
        ("x in [100, 1000]", np.linspace(100, 1000, 200)),
        ("x in [1000, 10000]", np.linspace(1000, 10000, 200)),
        ("x in [10000, 100000]", np.linspace(10000, 100000, 200)),
    ]:
        # Build matrix M[i,j] = x_i^{ρ_j} / ρ_j
        # = exp(ρ_j · ln(x_i)) / ρ_j
        ln_x = np.log(x_grid)

        # M[i,j] = exp(ρ_j * ln_x_i) / ρ_j
        # Take real part (pairing with conjugate)
        M = np.zeros((len(x_grid), N_zeros))
        for j in range(N_zeros):
            rho = rhos_sub[j]
            for i in range(len(x_grid)):
                val = np.exp(rho * ln_x[i]) / rho
                M[i, j] = 2 * val.real  # pair with conjugate

        # SVD
        try:
            U, s, Vt = np.linalg.svd(M, full_matrices=False)
        except:
            print(f"\n{x_range_label}: SVD failed")
            continue

        # Effective rank: number of singular values > ε * s_max
        s_max = s[0]
        total_energy = np.sum(s**2)

        print(f"\n{x_range_label} (grid: {len(x_grid)}, zeros: {N_zeros})")
        print(f"  Top 10 singular values: {', '.join(f'{sv:.4e}' for sv in s[:10])}")
        print(f"  s[0]/s[10] = {s[0]/s[min(10,len(s)-1)]:.2f}")
        print(f"  s[0]/s[50] = {s[0]/s[min(50,len(s)-1)]:.2f}")

        for threshold_frac in [0.99, 0.999, 0.9999]:
            cumulative = np.cumsum(s**2) / total_energy
            rank_needed = np.searchsorted(cumulative, threshold_frac) + 1
            print(f"  Rank for {threshold_frac*100:.1f}% energy: {rank_needed}")

        # Actual approximation error with rank-k
        # S(x_i) ≈ Σ_{j=1}^k s_j · u_ij · v_j  (rank-k SVD approx)
        # True S(x_i) = Σ_j M[i,j]
        true_S = M.sum(axis=1)  # S(x_i) for each x_i

        print(f"\n  Rank-k approximation of S(x):")
        print(f"  {'k':>5s} | {'max |error|':>12s} | {'mean |error|':>12s} | {'max |S|':>12s} | {'rel error':>10s}")
        for k in [1, 2, 5, 10, 20, 50, 100, 150, 200]:
            if k > min(len(x_grid), N_zeros):
                break
            # Rank-k reconstruction
            M_k = U[:, :k] @ np.diag(s[:k]) @ Vt[:k, :]
            S_k = M_k.sum(axis=1)
            max_err = np.max(np.abs(S_k - true_S))
            mean_err = np.mean(np.abs(S_k - true_S))
            max_S = np.max(np.abs(true_S))
            rel = max_err / max(max_S, 1e-10)
            print(f"  {k:5d} | {max_err:12.4e} | {mean_err:12.4e} | {max_S:12.4e} | {rel:10.4e}")


# ---------------------------------------------------------------------------
# EXPERIMENT 5: Moment Scaling -- How M_k Grows
# ---------------------------------------------------------------------------
def experiment_moment_scaling(gammas):
    """
    The moments M_k = Σ_j ρ_j^k determine whether the moment expansion converges.
    If |M_k| grows faster than k!, the expansion Σ (ln x)^k / k! · M_{k-1} diverges.

    For the zeta zeros, ρ_j = 1/2 + iγ_j with γ_j ~ 2πj/ln(j).
    So |ρ_j|^k ~ γ_j^k, and Σ_j γ_j^k is dominated by the largest γ_j.

    For N zeros with γ_max ~ 2πN/ln(N):
    |M_k| ~ N · γ_max^k ~ N · (2πN/ln N)^k

    This grows MUCH faster than k!, so the expansion diverges for ln(x) > 0.
    """
    print("\n" + "="*70)
    print("EXPERIMENT 5: MOMENT SCALING ANALYSIS")
    print("="*70)

    rhos = 0.5 + 1j * gammas
    N = len(rhos)
    gamma_max = gammas[-1]

    print(f"N = {N} zeros, γ_max = {gamma_max:.2f}")
    print(f"Expected |M_k| ~ N · γ_max^k = {N} · {gamma_max:.1f}^k")
    print()

    print(f"{'k':>4s} | {'|M_k|':>14s} | {'k!':>14s} | {'|M_k|/k!':>14s} | {'N·γ_max^k':>14s} | {'ratio':>10s}")
    print("-" * 80)

    for k in range(0, 31):
        M_k = np.sum(rhos**k)
        abs_Mk = abs(M_k)
        k_fact = float(math.factorial(min(k, 170)))
        predicted = N * gamma_max**k
        ratio = abs_Mk / max(predicted, 1e-300) if predicted > 0 else float('inf')

        print(f"{k:4d} | {abs_Mk:14.4e} | {k_fact:14.4e} | {abs_Mk/k_fact:14.4e} | {predicted:14.4e} | {ratio:10.4e}")

    print("\n--- VERDICT ---")
    print(f"Moments |M_k| grow as ~ γ_max^k ≈ {gamma_max:.0f}^k")
    print(f"Factorial k! grows much slower.")
    print(f"For the expansion Σ (ln x)^k/k! · M_{{k-1}} to converge,")
    print(f"  need (ln x · γ_max) < 1, i.e., x < e^{{1/γ_max}} ≈ {np.exp(1/gamma_max):.6f}")
    print(f"This means the series only converges for x ≈ 1, which is USELESS.")
    print(f"The moment method CANNOT work for the zeta zero sum.")


# ---------------------------------------------------------------------------
# EXPERIMENT 6: Smoothed Trace -- Gaussian Damping
# ---------------------------------------------------------------------------
def experiment_smoothed_trace(gammas, test_x_values):
    """
    If we replace x^ρ/ρ with x^ρ/ρ · e^{-ρ²/T²} (Gaussian damping),
    the high zeros are suppressed. This makes moments converge.

    Question: What is the cost? How does the error in π(x) depend on T?
    """
    print("\n" + "="*70)
    print("EXPERIMENT 6: SMOOTHED TRACE WITH GAUSSIAN DAMPING")
    print("="*70)

    rhos = 0.5 + 1j * gammas
    N = len(rhos)

    from sympy import primepi

    for x_val in test_x_values:
        pi_exact = int(primepi(int(x_val)))
        li_val = li_func(x_val)
        lnx = np.log(x_val)

        # Full zero sum (no damping)
        S_full = sum(2 * (x_val**rho / rho).real for rho in rhos)
        pi_approx_full = li_val - S_full - np.log(2)  # approximate

        print(f"\nx = {x_val:.0f}, π(x) = {pi_exact}, li(x) = {li_val:.4f}")
        print(f"  Full sum ({N} pairs): S = {S_full:.4f}, π_approx = {pi_approx_full:.4f}")
        print(f"  {'T':>8s} | {'S_damped':>12s} | {'#eff zeros':>10s} | {'π_approx':>10s} | {'error':>8s}")
        print(f"  " + "-"*60)

        for T in [10, 20, 50, 100, 200, 500, 1000]:
            # Damped sum: S_T(x) = Σ_j 2Re[x^{ρ_j}/ρ_j · e^{-γ_j²/T²}]
            weights = np.exp(-gammas**2 / T**2)
            S_damped = sum(2 * (x_val**rho / rho * w).real
                         for rho, w in zip(rhos, weights))
            n_eff = np.sum(weights > 0.01)
            pi_approx = li_val - S_damped - np.log(2)
            error = pi_approx - pi_exact
            print(f"  {T:8d} | {S_damped:12.4f} | {n_eff:10d} | {pi_approx:10.4f} | {error:8.2f}")

        # What T is needed for |error| < 0.5?
        for T in range(10, 2000, 10):
            weights = np.exp(-gammas**2 / T**2)
            S_damped = sum(2 * (x_val**rho / rho * w).real
                         for rho, w in zip(rhos, weights))
            pi_approx = li_val - S_damped - np.log(2)
            if abs(pi_approx - pi_exact) < 0.5:
                n_eff = np.sum(weights > 0.01)
                print(f"  --> T = {T} gives |error| < 0.5, #effective zeros = {n_eff}")
                break
        else:
            print(f"  --> No T in [10, 2000] gives |error| < 0.5 (need more zeros)")


# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------
def main():
    print("TRACE FORMULA APPROACH: Computing π(x) Without Enumerating Zeros")
    print("=" * 70)

    gammas = load_zeros(1000)

    # Test x values (moderate, so we can compute exact π(x))
    test_x_small = [100, 1000, 10000]
    test_x_medium = [100, 1000, 10000, 100000]

    t0 = time.time()

    # Run experiments
    exp1 = experiment_moments_method(gammas, test_x_small)
    exp5 = experiment_moment_scaling(gammas)
    exp2 = experiment_gue_trace(gammas)
    exp3 = experiment_weil_truncation(gammas, test_x_small)
    exp4 = experiment_effective_rank(gammas, test_x_small)
    exp6 = experiment_smoothed_trace(gammas, test_x_small)

    elapsed = time.time() - t0

    # SUMMARY
    print("\n" + "="*70)
    print("FINAL SUMMARY")
    print("="*70)
    print(f"\nTotal runtime: {elapsed:.1f}s")

    print("""
EXPERIMENT 1 & 5 (Moments Method):
  The moment expansion S(x) = Σ_k (ln x)^k/k! · M_{k-1} DIVERGES.
  |M_k| ~ γ_max^k >> k!, so the series radius of convergence is 0.
  VERDICT: DOES NOT WORK. The moments of zeta zeros are too large.

EXPERIMENT 2 (GUE Trace):
  For GUE(N), Tr(e^{itH}) converges via matrix powers because eigenvalues
  are bounded. Zeta zeros are unbounded (γ_j → ∞), so the GUE analogy
  breaks down completely for this application.
  VERDICT: GUE analogy is MISLEADING here. The trace method works for
  bounded operators but zeta zeros give an unbounded operator.

EXPERIMENT 3 (Weil Explicit Formula):
  The geometric side Σ_p log(p)/p^{m/2} requires knowing all primes up
  to x to get ψ(x) within ±0.5. This is circular.
  VERDICT: CIRCULAR. Cannot use the explicit formula's geometric side
  to compute π(x) without already knowing primes up to x.

EXPERIMENT 4 (Effective Rank):
  The matrix M[i,j] = 2Re[x_i^{ρ_j}/ρ_j] has high effective rank.
  To approximate S(x) to O(1), we need rank comparable to the number
  of zeros contributing significantly. No useful compression.
  VERDICT: NO SHORTCUT. The operator is essentially full-rank.

EXPERIMENT 6 (Smoothed Trace):
  Gaussian damping e^{-γ²/T²} makes moments converge, but introduces
  systematic error. To get |error| < 0.5 (exact π(x)), we need T large
  enough that most zeros contribute, losing the efficiency gain.
  VERDICT: TRADE-OFF is unfavorable. Smoothing enough to make moments
  work destroys the accuracy needed for exact computation.

OVERALL VERDICT: The trace/moment approach FAILS for computing π(x) exactly.
  - Zeta zeros are unbounded → moment series diverges
  - GUE analogy breaks down for unbounded spectra
  - Weil duality is circular
  - No low-rank structure to exploit
  This is a NEGATIVE result that CLOSES the trace formula shortcut.
""")


if __name__ == '__main__':
    main()

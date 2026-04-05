"""
Sublinear Correction: Pushing below O(x^{2/3})

The experiments show δ(n) is partially compressible (5x better than random).
Question: Can we exploit this to compute π(x) in o(x^{2/3})?

Key observation from spectral compression experiment:
- DCT captures 99% of δ(n) energy in 10.4% of coefficients
- The correction ε in the Dickman decomposition is 95% energy in 2.2%

This suggests a HYBRID algorithm:
1. Compute R(x) in O(polylog) — the smooth part
2. Compute the "spectrally dominant" correction using O(x^α) work for α < 2/3
3. Round to get π(x)

The question is whether the dominant correction terms can be identified
without computing all of them.

IDEA: The Lagarias-Odlyzko analytic method gives O(x^{1/2+ε}).
It uses the Perron formula with T = x^{1/2+ε}.
The Galway variant (2004) improved constants.
Can we combine DCT compression with the analytic method?

ALSO: The Deléglise-Rivat O(x^{2/3}/ln²x) uses combinatorial sieving.
The Platt-Trudgian verification of RH up to 3×10^12 zeros enables
conditional algorithms. Under RH, the analytic method works.

What's the BEST achievable complexity using both analytic and combinatorial
methods, with the compressibility of δ as a structural assumption?
"""

import numpy as np
import sympy
from sympy import primepi, prime
from mpmath import mp, li, zetazero, log, exp
import math
import time

def compute_correction_spectrum(N=10000):
    """
    Compute the DCT spectrum of δ(n) = π(n) - R(n) and analyze
    which coefficients are needed for exact rounding.
    """
    mp.dps = 20

    print("=== Correction spectrum analysis ===\n")

    # Compute δ(n)
    deltas = []
    for n in range(2, N + 1):
        pi_n = int(sympy.primepi(n))
        R_n = float(li(n) - 0.5 * li(n**0.5))
        deltas.append(pi_n - R_n)

    deltas = np.array(deltas)

    from scipy.fft import dct, idct

    # DCT of δ
    D = dct(deltas, type=2, norm='ortho')
    magnitudes = np.abs(D)

    # Find minimum number of coefficients for exact rounding
    # We need: |δ(n) - δ_approx(n)| < 0.5 for all n
    sorted_indices = np.argsort(magnitudes)[::-1]

    for target_error in [0.5, 0.3, 0.1]:
        for k in range(1, len(D) + 1):
            # Keep top k coefficients
            D_approx = np.zeros_like(D)
            D_approx[sorted_indices[:k]] = D[sorted_indices[:k]]
            delta_approx = idct(D_approx, type=2, norm='ortho')
            max_err = np.max(np.abs(deltas - delta_approx))
            if max_err < target_error:
                print(f"  N={N}: {k}/{len(D)} DCT coefficients ({k/len(D)*100:.1f}%) "
                      f"for max_error < {target_error}")
                break
        else:
            print(f"  N={N}: ALL coefficients needed for max_error < {target_error}")

    # Analyze which coefficients are needed
    # Are they low-frequency or spread across the spectrum?
    k_for_rounding = None
    for k in range(1, len(D) + 1):
        D_approx = np.zeros_like(D)
        D_approx[sorted_indices[:k]] = D[sorted_indices[:k]]
        delta_approx = idct(D_approx, type=2, norm='ortho')
        max_err = np.max(np.abs(deltas - delta_approx))
        if max_err < 0.5:
            k_for_rounding = k
            break

    if k_for_rounding:
        needed_indices = sorted(sorted_indices[:k_for_rounding])
        print(f"\n  Indices of needed coefficients (for rounding):")
        print(f"  Min index: {min(needed_indices)}")
        print(f"  Max index: {max(needed_indices)}")
        print(f"  Median index: {np.median(needed_indices):.0f}")
        print(f"  Mean index: {np.mean(needed_indices):.0f}")
        print(f"  Are they concentrated? Top 10%: {sum(1 for i in needed_indices if i < len(D)*0.1)}")
        print(f"  Top 50%: {sum(1 for i in needed_indices if i < len(D)*0.5)}")

    return deltas, D

def test_hybrid_algorithm(N=5000):
    """
    Test a HYBRID approach:
    1. R(x) approximation (polylog)
    2. Partial zero sum with K zeros
    3. DCT correction of the residual

    Can we get exact π(x) with K << √x zeros + O(polylog) DCT work?
    """
    mp.dps = 20

    print("\n=== Hybrid algorithm test ===\n")

    # Precompute zeta zeros
    max_K = 50
    gammas = [float(zetazero(k).imag) for k in range(1, max_K + 1)]

    # Compute δ(n) for training
    deltas_train = []
    ns = np.arange(2, N + 1, dtype=float)
    for n_val in ns:
        n = int(n_val)
        pi_n = int(sympy.primepi(n))
        R_n = float(li(n) - 0.5 * li(n**0.5))
        deltas_train.append(pi_n - R_n)
    deltas_train = np.array(deltas_train)

    # Compute zero sum contribution for various K
    for K in [5, 10, 20, 30, 50]:
        zero_correction = np.zeros(len(ns))
        for k in range(K):
            gamma_k = gammas[k]
            for i, x in enumerate(ns):
                log_x = np.log(x)
                sqrt_x = np.sqrt(x)
                phase = gamma_k * log_x
                rho_abs = np.sqrt(0.25 + gamma_k**2)
                zero_correction[i] += 2 * sqrt_x * np.cos(phase) / (rho_abs * log_x)

        # Residual after K zeros
        residual = deltas_train + zero_correction  # δ = -zero_sum + small, so residual = δ + zero_sum

        max_residual = np.max(np.abs(residual))
        mean_residual = np.mean(np.abs(residual))

        # Can the residual be captured by DCT compression?
        from scipy.fft import dct, idct
        R = dct(residual, type=2, norm='ortho')
        magnitudes = np.abs(R)
        sorted_idx = np.argsort(magnitudes)[::-1]

        # Find k_dct for max error < 0.5
        k_dct = None
        for kd in range(1, len(R) + 1):
            R_approx = np.zeros_like(R)
            R_approx[sorted_idx[:kd]] = R[sorted_idx[:kd]]
            res_approx = idct(R_approx, type=2, norm='ortho')
            if np.max(np.abs(residual - res_approx)) < 0.5:
                k_dct = kd
                break

        print(f"  K={K:>3} zeros: max_residual={max_residual:.2f}, "
              f"mean_residual={mean_residual:.2f}, "
              f"DCT coeffs for rounding: {k_dct if k_dct else '>ALL'}")

    print("\n  Key: if K zeros + k_dct DCT coefficients << N, we have a speedup")
    print("  But computing DCT coefficients of the residual requires knowing the residual,")
    print("  which requires knowing π(x) — circular unless the DCT coefficients are PREDICTABLE")

def test_dct_coefficient_predictability(max_N=5000, step=500):
    """
    KEY EXPERIMENT: Are the DCT coefficients of δ(n) predictable from
    shorter sequences?

    If we can compute the DCT of δ for n ≤ N₁ < N, can we predict the
    DCT coefficients for N₁ < n ≤ N?

    This would enable a truly sublinear algorithm:
    1. Compute δ(n) for n ≤ N₁ (via Meissel-Lehmer, cost O(N₁^{2/3}))
    2. Predict DCT for n > N₁
    3. Round to get π(x)
    """
    mp.dps = 20
    from scipy.fft import dct

    print("\n=== DCT coefficient predictability ===\n")

    # Compute δ(n) for full range
    deltas = []
    for n in range(2, max_N + 1):
        pi_n = int(sympy.primepi(n))
        R_n = float(li(n) - 0.5 * li(n**0.5))
        deltas.append(pi_n - R_n)
    deltas = np.array(deltas)

    # Compute DCT at different truncation points
    # Check: do the first K DCT coefficients stabilize?
    print(f"  DCT coefficient stability (first 10 coefficients):")
    print(f"  {'N':>6}  {'c[0]':>8}  {'c[1]':>8}  {'c[2]':>8}  {'c[3]':>8}  {'c[4]':>8}")

    for N in range(step, max_N + 1, step):
        D = dct(deltas[:N-1], type=2, norm='ortho')
        # Normalize by sqrt(N) to compare across different N
        D_norm = D / np.sqrt(N)
        print(f"  {N:>6}  {D_norm[0]:>8.4f}  {D_norm[1]:>8.4f}  {D_norm[2]:>8.4f}  "
              f"{D_norm[3]:>8.4f}  {D_norm[4]:>8.4f}")

    # Check correlation between DCT on [2, N/2] and [2, N]
    print(f"\n  Correlation of DCT coefficients at different scales:")
    N1 = max_N // 2
    N2 = max_N
    D1 = dct(deltas[:N1-1], type=2, norm='ortho')
    D2 = dct(deltas[:N2-1], type=2, norm='ortho')

    # Subsample D2 to match D1 length
    D2_sub = D2[:len(D1)]

    correlation = np.corrcoef(D1, D2_sub)[0, 1]
    print(f"  Correlation between DCT(δ, N={N1}) and DCT(δ, N={N2}): {correlation:.6f}")

    # Individual coefficient correlations
    for k in [0, 1, 2, 5, 10, 50]:
        if k < len(D1):
            vals_at_k = []
            for N in range(step, max_N + 1, step):
                D = dct(deltas[:N-1], type=2, norm='ortho')
                vals_at_k.append(D[k] / np.sqrt(N))
            vals = np.array(vals_at_k)
            # Is this converging?
            if len(vals) > 2:
                trend = np.polyfit(range(len(vals)), vals, 1)
                print(f"  c[{k}]/√N: mean={vals.mean():.4f}, std={vals.std():.4f}, "
                      f"trend={trend[0]:.6f}/step")

def test_zero_sum_as_dct():
    """
    FINAL INSIGHT: The zero sum IS a type of transform.

    Σ_ρ R(x^ρ) = Σ_k [a_k sin(γ_k log x) + b_k cos(γ_k log x)] / ...

    This is a Fourier-like expansion in log-frequency space.
    The "frequencies" are the zeta zeros γ_k.
    The "coefficients" a_k, b_k are determined by the zeros.

    If we change variables to u = log(x), this becomes:
    Σ_k c_k e^{iγ_k u} · envelope(u)

    This is a FINITE RATE OF INNOVATION (FRI) signal!
    FRI theory says: if a signal has K innovations per unit interval,
    it can be reconstructed from 2K+1 samples.

    For the zero sum with K zeros: 2K+1 samples of π(x) suffice!
    But we need π(x) values to get the samples — circular again.

    UNLESS: we can get approximate samples from the Meissel-Lehmer method
    at a few carefully chosen points, then interpolate.
    """
    print("\n=== Finite Rate of Innovation analysis ===\n")

    mp.dps = 20

    # The zero sum has K "innovations" (one per zero)
    # In the range [log(a), log(b)], the effective number of oscillations
    # from zero k is γ_k * (log(b) - log(a)) / (2π)

    gammas = [float(zetazero(k).imag) for k in range(1, 51)]

    log_a = np.log(100)
    log_b = np.log(10000)
    interval = log_b - log_a

    print(f"  Interval: [{100}, {10000}], log-length = {interval:.2f}")
    print(f"\n  Oscillations per zero in this interval:")

    total_oscillations = 0
    for k, gamma in enumerate(gammas[:20]):
        osc = gamma * interval / (2 * np.pi)
        total_oscillations += osc
        print(f"    γ_{k+1} = {gamma:.2f}: {osc:.1f} oscillations")

    print(f"\n  Total oscillations from 20 zeros: {total_oscillations:.0f}")
    print(f"  By FRI theory: need ≥ {2*int(total_oscillations)+1} samples")
    print(f"  Interval length: {10000-100} integers")
    print(f"  Sampling density: {(2*int(total_oscillations)+1)/(10000-100)*100:.1f}%")

    # For x ~ 10^100:
    # γ_k ~ k * 2π/ln(k) for large k (by zero density)
    # Total oscillations from K zeros ≈ Σ_{k=1}^{K} γ_k * Δ / (2π)
    # ≈ K² * Δ / ln(K) for Δ = log-interval width
    # Need K ~ √x zeros, Δ ~ log(x)
    # Total oscillations ~ x * log(x) / ln(√x) ~ x
    # So sampling density ~ 100% — no savings via FRI

    print(f"\n  Scaling to x ~ 10^100:")
    print(f"  Need K ~ √x ~ 10^50 zeros")
    print(f"  Total oscillations ~ 10^100")
    print(f"  → Need ~10^100 samples (= all integers)")
    print(f"  FRI provides NO savings at scale")

if __name__ == "__main__":
    print("=" * 60)
    print("EXPERIMENT: Sublinear Correction via Spectral Methods")
    print("=" * 60)

    deltas, D = compute_correction_spectrum(5000)
    test_hybrid_algorithm(3000)
    test_dct_coefficient_predictability(5000)
    test_zero_sum_as_dct()

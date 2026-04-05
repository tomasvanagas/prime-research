#!/usr/bin/env python3
"""
Experiment: Compressibility of the Meissel-Lehmer sieve function.

The Lucy_Hedgehog / Meissel-Lehmer algorithm computes:
  S(v, p) = |{n ≤ v : n has no prime factor ≤ p}|

for v in the set V = {⌊x/k⌋ : k=1,...,√x} ∪ {1,...,√x} (about 2√x values).

Each prime p updates: S(v,p) = S(v,p-1) - [S(v/p,p-1) - S(p-1,p-1)]

THE KEY QUESTION: Is S(v, p) as a function of v "compressible"?
- If S(·, p) is piecewise polynomial of degree d, updates preserve this
  and the per-prime cost drops from O(√x) to O(d * breakpoints).
- If S(·, p) has a Fourier-sparse representation, FFT could help.
- If S(·, p) is low-rank in some 2D sense, SVD could help.

We test these on small-to-medium x values (up to 10^8).
"""

import numpy as np
from sympy import primerange, primepi
import time


def lucy_hedgehog(x):
    """
    Standard Lucy_Hedgehog algorithm for π(x).
    Returns the full S(v, final_p) array for analysis.
    """
    sqrtx = int(x**0.5)

    # V = set of values ⌊x/k⌋ for k=1,...,x, deduplicated
    # Stored as: small[v] for v ≤ √x, large[k] for v = ⌊x/k⌋ where k ≤ √x

    small = [0] * (sqrtx + 2)  # S(v, p) for v = 0, 1, ..., sqrtx
    large = [0] * (sqrtx + 2)  # S(⌊x/k⌋, p) for k = 1, 2, ..., sqrtx

    # Initialize: S(v, 1) = v - 1 (all integers from 2 to v)
    for i in range(sqrtx + 2):
        small[i] = max(i - 1, 0)
    for k in range(1, sqrtx + 1):
        large[k] = x // k - 1

    def get_S(v):
        if v <= sqrtx:
            return small[v]
        else:
            return large[x // v]

    def set_S(v, val):
        if v <= sqrtx:
            small[v] = val
        else:
            large[x // v] = val

    # Sieve: for each prime p ≤ √x
    # S(v, p) = S(v, p-1) - [S(⌊v/p⌋, p-1) - S(p-1, p-1)]
    history = []  # Store snapshots for analysis

    # Snapshot initial state
    history.append(('init', list(small), list(large)))

    for p in range(2, sqrtx + 1):
        if small[p] == small[p-1]:
            continue  # p is not prime (already sieved)

        # p is prime
        p_count = small[p - 1]  # π(p-1) = S(p-1, p-1)

        # Update large values
        for k in range(1, min(sqrtx, x // (p * p)) + 1):
            v = x // k
            vp = v // p
            large[k] -= get_S(vp) - p_count

        # Update small values
        for v in range(sqrtx, p*p - 1, -1):
            vp = v // p
            small[v] -= get_S(vp) - p_count

        history.append((p, list(small), list(large)))

    return small[sqrtx] if sqrtx <= sqrtx else large[1], history, small, large


def analyze_sieve_shape(x):
    """Analyze how S(v, p) changes shape as p increases."""
    print(f"\n{'='*70}")
    print(f"SIEVE FUNCTION SHAPE ANALYSIS: x = {x}")
    print(f"{'='*70}")

    sqrtx = int(x**0.5)
    result, history, final_small, final_large = lucy_hedgehog(x)

    print(f"  √x = {sqrtx}, π(x) = {result}")
    print(f"  Number of sieve steps: {len(history) - 1}")

    # Analyze each snapshot
    for step_name, small, large in history[:1] + history[1::max(1, len(history)//5)] + history[-1:]:
        # Construct the full S(v) function for analysis
        # v values and S values
        vs = list(range(1, sqrtx + 1))
        Svs = [small[v] for v in vs]

        # Also the large values
        for k in range(sqrtx, 0, -1):
            v = x // k
            if v > sqrtx:
                vs.append(v)
                Svs.append(large[k])

        vs = np.array(vs, dtype=float)
        Svs = np.array(Svs, dtype=float)

        # Sort by v
        order = np.argsort(vs)
        vs = vs[order]
        Svs = Svs[order]

        # Fit polynomial
        # S(v) should be approximately v * euler_product for large v
        # Normalize: S(v) / v → density
        mask = vs > 10
        density = np.zeros_like(Svs, dtype=float)
        density[mask] = Svs[mask] / vs[mask]

        # Check: is S(v) well-approximated by a polynomial in v?
        if len(vs[mask]) > 5:
            for deg in [1, 2, 3, 5]:
                try:
                    coeffs = np.polyfit(vs[mask], Svs[mask], deg)
                    pred = np.polyval(coeffs, vs[mask])
                    max_err = np.max(np.abs(Svs[mask] - pred))
                    rmse = np.sqrt(np.mean((Svs[mask] - pred)**2))
                except:
                    max_err = rmse = float('inf')

        # Also check: is S(v) well-described by C * v / log(v)?
        log_vs = np.log(vs[mask])
        pred_logform = vs[mask] / log_vs
        scale = np.sum(Svs[mask] * pred_logform) / np.sum(pred_logform**2)
        residual_log = np.max(np.abs(Svs[mask] - scale * pred_logform))

        if step_name == 'init':
            print(f"\n  After init (S(v,1) = v-1):")
        else:
            print(f"\n  After sieving by p={step_name}:")
        print(f"    Range of S: [{Svs[0]:.0f}, {Svs[-1]:.0f}]")
        print(f"    Density S(v)/v at v={vs[-1]:.0f}: {density[-1]:.4f}")
        if len(vs[mask]) > 5:
            print(f"    Poly deg 1 max error: {max_err:.2f}")
            print(f"    v/log(v) fit residual: {residual_log:.2f}, scale={scale:.4f}")


def test_piecewise_polynomial():
    """
    KEY TEST: Is S(v, p) piecewise polynomial?

    If S(v, p) is a polynomial of degree d on intervals between
    breakpoints, and sieving by the next prime PRESERVES this
    (with a bounded number of new breakpoints), then the total
    cost is O(breakpoints × primes × d).
    """
    print(f"\n{'='*70}")
    print("PIECEWISE POLYNOMIAL TEST")
    print(f"{'='*70}")

    for x in [1000, 10000, 100000]:
        sqrtx = int(x**0.5)

        # Compute S(v, p) for all v = 1, ..., x and each prime p ≤ √x
        # (Brute force for small x)

        if x > 10000:
            # Use Lucy Hedgehog for larger x
            _, history, _, _ = lucy_hedgehog(x)

            for step_name, small, large in [history[0], history[len(history)//2], history[-1]]:
                # Just analyze the "small" part: S(v, p) for v = 1, ..., √x
                vs = np.arange(1, sqrtx + 1)
                Svs = np.array([small[v] for v in vs], dtype=float)

                # Look for breakpoints: where does the second difference change sign?
                if len(Svs) > 4:
                    diffs2 = np.diff(Svs, n=2)
                    sign_changes = np.sum(diffs2[:-1] * diffs2[1:] < 0)

                    # Higher-order differences
                    diffs3 = np.diff(Svs, n=3) if len(Svs) > 5 else np.array([0])
                    diffs4 = np.diff(Svs, n=4) if len(Svs) > 6 else np.array([0])

                    p_label = f"p={step_name}" if step_name != 'init' else "init"
                    print(f"\n  x={x}, {p_label}, v = 1..{sqrtx}")
                    print(f"    2nd diff sign changes: {sign_changes} / {len(diffs2)}")
                    print(f"    3rd diff |max|: {np.max(np.abs(diffs3)):.2f}")
                    print(f"    4th diff |max|: {np.max(np.abs(diffs4)):.2f}")

        else:
            # Full computation for small x
            is_composite = [False] * (x + 1)
            primes_list = []

            for p in range(2, sqrtx + 1):
                if is_composite[p]:
                    continue
                primes_list.append(p)
                for m in range(p*p, x+1, p):
                    is_composite[m] = True

            # Compute S(v, p) for all v and each prime stage
            S = np.zeros((x + 1,), dtype=int)
            for v in range(2, x + 1):
                S[v] = v - 1  # S(v, 1) = v - 1

            for p_idx, p in enumerate(primes_list):
                p_count = S[p - 1]
                new_S = S.copy()
                for v in range(x, p*p - 1, -1):
                    new_S[v] = S[v] - (S[v // p] - p_count)
                S = new_S

                # Analyze S(v) after this prime
                vs = np.arange(1, x + 1)
                Svs = S[1:].astype(float)

                # Differences
                d1 = np.diff(Svs)
                d2 = np.diff(d1)

                # d1 should be 0 or 1 eventually (since S is a step function counting integers)
                # But S(v, p) is NOT a step function - it's the count of smooth numbers.

                # How many distinct values does d1 take?
                unique_d1 = len(np.unique(d1))
                unique_d2 = len(np.unique(d2))

                if p_idx < 3 or p_idx == len(primes_list) - 1:
                    print(f"\n  x={x}, after p={p}:")
                    print(f"    S range: [{S[1]}, {S[x]}]")
                    print(f"    unique 1st diffs: {unique_d1}")
                    print(f"    unique 2nd diffs: {unique_d2}")
                    print(f"    first 20 diffs: {d1[:20].astype(int).tolist()}")


def fourier_analysis_of_sieve():
    """
    Test: Does S(v, p) have a sparse Fourier representation?
    If the DFT of S(v, p) has few large coefficients, FFT-based
    methods could speed up the sieve.
    """
    print(f"\n{'='*70}")
    print("FOURIER ANALYSIS OF SIEVE FUNCTION")
    print(f"{'='*70}")

    for x in [10000, 100000]:
        sqrtx = int(x**0.5)
        _, history, _, _ = lucy_hedgehog(x)

        for step_name, small, large in [history[-1]]:
            vs = np.arange(1, sqrtx + 1)
            Svs = np.array([small[v] for v in vs], dtype=float)

            # Detrend: subtract linear fit
            slope, intercept = np.polyfit(vs, Svs, 1)
            residual = Svs - (slope * vs + intercept)

            # DFT of residual
            fft_res = np.fft.rfft(residual)
            magnitudes = np.abs(fft_res)
            total_energy = np.sum(magnitudes**2)

            # How many Fourier modes capture 90%, 99% of energy?
            sorted_mag = np.sort(magnitudes**2)[::-1]
            cumulative = np.cumsum(sorted_mag) / total_energy
            k90 = np.searchsorted(cumulative, 0.90) + 1
            k99 = np.searchsorted(cumulative, 0.99) + 1

            # Sparsity ratio
            n_modes = len(magnitudes)
            threshold = 0.01 * magnitudes.max()
            n_significant = np.sum(magnitudes > threshold)

            print(f"\n  x={x}, final sieve (p up to {sqrtx}):")
            print(f"    Signal length: {len(Svs)}")
            print(f"    DFT modes: {n_modes}")
            print(f"    Modes for 90% energy: {k90} ({k90/n_modes*100:.1f}%)")
            print(f"    Modes for 99% energy: {k99} ({k99/n_modes*100:.1f}%)")
            print(f"    Significant modes (>1% of max): {n_significant} ({n_significant/n_modes*100:.1f}%)")
            print(f"    Top 5 mode magnitudes: {sorted_mag[:5]}")
            print(f"    → Fourier {'SPARSE' if n_significant < n_modes * 0.1 else 'NOT sparse'}")


def rank_of_sieve_updates():
    """
    Test: Does the UPDATE ΔS(v) = S(v,p) - S(v,p-1) have low rank?
    If so, each sieve step could be compressed.
    """
    print(f"\n{'='*70}")
    print("RANK OF SIEVE UPDATES")
    print(f"{'='*70}")

    x = 10000
    sqrtx = int(x**0.5)
    _, history, _, _ = lucy_hedgehog(x)

    # Collect all update vectors
    updates = []
    prev_small = None
    for step_name, small, large in history:
        if prev_small is not None and step_name != 'init':
            delta = np.array(small) - np.array(prev_small)
            if np.any(delta != 0):
                updates.append(delta)
        prev_small = small

    if len(updates) > 1:
        # Stack into matrix
        U = np.array(updates, dtype=float)
        print(f"\n  x={x}: {U.shape[0]} update vectors, each of length {U.shape[1]}")

        # SVD
        sv = np.linalg.svd(U, compute_uv=False)
        sv_norm = sv / sv[0] if sv[0] > 0 else sv
        energy = np.cumsum(sv**2) / np.sum(sv**2) if np.sum(sv**2) > 0 else np.zeros_like(sv)
        k90 = np.searchsorted(energy, 0.90) + 1
        k99 = np.searchsorted(energy, 0.99) + 1
        rank = np.sum(sv > sv[0] * 1e-10)

        print(f"    Rank: {rank} / {min(U.shape)}")
        print(f"    Singular values (top 10): {sv[:10]}")
        print(f"    k for 90% energy: {k90}")
        print(f"    k for 99% energy: {k99}")
        print(f"    → Updates are {'LOW-RANK' if k99 < rank * 0.5 else 'full-rank'}")


def hyperbolic_lattice_test():
    """
    WILD IDEA: The set of v-values in Lucy_Hedgehog is {⌊x/k⌋ : k=1,...,√x}.
    This is a "hyperbolic lattice" (related to the divisor problem).

    On this lattice, the sieve update S(v,p) = S(v,p-1) - S(v/p,p-1) + c
    is a SHIFT operation (v → v/p).

    If we use a DIFFERENT parameterization of the lattice
    (e.g., log-space), the shift becomes translation, and
    translation is diagonalized by Fourier transform!

    Can we solve the sieve in Fourier-log space?
    """
    print(f"\n{'='*70}")
    print("HYPERBOLIC LATTICE / LOG-FOURIER TEST")
    print(f"{'='*70}")

    for x in [10000, 100000]:
        sqrtx = int(x**0.5)
        _, history, _, _ = lucy_hedgehog(x)

        # Get final S(v) on the hyperbolic lattice
        final_step = history[-1]
        small = final_step[1]
        large = final_step[2]

        # Construct (v, S(v)) on the full lattice
        lattice_v = []
        lattice_S = []

        for v in range(1, sqrtx + 1):
            lattice_v.append(v)
            lattice_S.append(small[v])

        for k in range(sqrtx, 0, -1):
            v = x // k
            if v > sqrtx:
                lattice_v.append(v)
                lattice_S.append(large[k])

        lattice_v = np.array(lattice_v, dtype=float)
        lattice_S = np.array(lattice_S, dtype=float)

        # In log space
        log_v = np.log(lattice_v + 1)

        # Fit S(v) in log-space
        # S(v) ≈ C * v / log(v) for the final sieve (prime counting)
        mask = lattice_v > 10
        log_S = np.log(lattice_S[mask] + 1)

        # Check if S(v) vs v relationship in log-log is linear
        # log S = α log v + β → S = e^β v^α
        alpha, beta = np.polyfit(np.log(lattice_v[mask]), log_S, 1)
        pred = np.exp(beta) * lattice_v[mask]**alpha
        max_rel_err = np.max(np.abs(lattice_S[mask] - pred) / (lattice_S[mask] + 1))

        print(f"\n  x={x}:")
        print(f"    Lattice size: {len(lattice_v)} points")
        print(f"    Log-log fit: S(v) ≈ {np.exp(beta):.4f} * v^{alpha:.4f}")
        print(f"    Max relative error: {max_rel_err:.4f}")

        # The key: in the sieve update S(v,p) = S(v,p-1) - S(v/p,p-1) + c,
        # the operation v → v/p in log-space is log(v) → log(v) - log(p).
        # This is a TRANSLATION in log-space!

        # If S(v, p-1) were a trigonometric polynomial in log(v),
        # then S(v/p, p-1) would be the same polynomial shifted by log(p).
        # The update would just modify Fourier coefficients.

        # Test: is the residual (S(v) - power law) sparse in Fourier-log space?
        residual = lattice_S[mask] - pred
        if len(residual) > 10:
            # Resample to uniform log-spacing
            n_uniform = 256
            log_v_uniform = np.linspace(np.log(lattice_v[mask][0]),
                                        np.log(lattice_v[mask][-1]),
                                        n_uniform)
            res_interp = np.interp(log_v_uniform, np.log(lattice_v[mask]), residual)

            fft_res = np.fft.rfft(res_interp)
            mag = np.abs(fft_res)
            total_energy = np.sum(mag**2)
            sorted_mag = np.sort(mag**2)[::-1]
            cum_energy = np.cumsum(sorted_mag) / total_energy if total_energy > 0 else np.zeros_like(sorted_mag)
            k90 = np.searchsorted(cum_energy, 0.90) + 1
            k99 = np.searchsorted(cum_energy, 0.99) + 1

            print(f"    Fourier of residual in log-space:")
            print(f"      Modes for 90% energy: {k90}/{n_uniform//2 + 1}")
            print(f"      Modes for 99% energy: {k99}/{n_uniform//2 + 1}")
            print(f"      → {'SPARSE' if k99 < 20 else 'NOT sparse'} in log-Fourier space")


if __name__ == "__main__":
    print("SIEVE FUNCTION COMPRESSION EXPERIMENTS")
    print("Can the Meissel-Lehmer sieve function be compressed?")

    analyze_sieve_shape(10000)
    analyze_sieve_shape(100000)

    test_piecewise_polynomial()

    fourier_analysis_of_sieve()

    rank_of_sieve_updates()

    hyperbolic_lattice_test()

    print(f"\n{'='*70}")
    print("SYNTHESIS")
    print(f"{'='*70}")
    print("""
    If ANY of these compression methods work, the sieve can be
    accelerated beyond O(x^{2/3}):

    1. Piecewise polynomial: If breakpoints grow polylog → O(polylog * primes)
    2. Fourier sparse: If k99 = O(polylog) → FFT-based sieve in O(polylog²)
    3. Low-rank updates: If rank = O(polylog) → compress update matrices
    4. Log-Fourier: If translation structure enables FFT → O(polylog per prime)

    The "translation in log-space" idea is the most promising:
    it exploits the MULTIPLICATIVE structure of the sieve
    (v → v/p is multiplication by 1/p, which is translation in log-space).
    """)

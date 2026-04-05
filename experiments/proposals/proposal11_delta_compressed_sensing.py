"""
PROPOSAL 3: Compressed Sensing on Delta(n) = p(n) - R^{-1}(n)
==============================================================

IDEA: delta(n) = p(n) - R^{-1}(n) has only O(log n) bits of information.
But computing it takes O(x^{2/3}) time. The question is whether delta(n)
has STRUCTURE that allows faster recovery.

KEY OBSERVATION: delta(n) is determined by the oscillatory part of the
explicit formula. In the frequency domain (over the zeta zeros), it's
a superposition of all zeros' contributions.

But delta(n) AS A FUNCTION OF n might have structure:
- It's piecewise constant (changes only at prime gaps)
- Its distribution is approximately Gaussian (by CLT over zeros)
- It might be compressible in a wavelet or other basis

COMPRESSED SENSING APPROACH:
1. Compute delta(n) exactly for O(polylog(N)) strategically chosen n values
2. Fit a model that predicts delta(N) from these samples
3. The model exploits the fact that delta is a function of the zero-sum

ALTERNATIVE: INTERPOLATION ON delta(n)
- delta(n) varies slowly most of the time (RH => |delta(n)| << x^{1/2+eps}/log(x))
- Abrupt changes happen at "prime gap events"
- If we can identify the O(polylog) gap events in [N-K, N+K], we can interpolate

ANOTHER ANGLE: RECURSIVE HALVING
- Compute p(n) by binary search on pi(x) = n
- Binary search takes O(log(x)) steps
- Each step: is pi(x_mid) >= n?
- If we can answer this COMPARISON query in O(polylog), done.
- Approximate pi(x_mid) ~ R(x_mid) gives ~50% of bits
- Need: a way to determine sign(pi(x_mid) - n) from residual info

COMPLEXITY: O(polylog) if delta has polylog effective degrees of freedom
ASSUMPTION: delta(n) is compressible in some basis, or the comparison
oracle pi(x) >= n can be answered from local information.
"""

import numpy as np
from sympy import prime, primepi, isprime, nextprime, prevprime
from mpmath import mp, mpf, li, log, exp, pi as mpi
import math

mp.dps = 30

def R_inv_approx(n):
    """Approximate R^{-1}(n) using Newton's method on R(x) = n."""
    from mpmath import lambertw
    # Initial guess: n * ln(n)
    if n < 2:
        return mpf(2)
    x = mpf(n) * log(mpf(n))
    for _ in range(20):
        from sympy import mobius
        R_val = mpf(0)
        for k in range(1, 30):
            mu_k = mobius(k)
            if mu_k != 0:
                R_val += mpf(mu_k) / k * li(x ** (mpf(1)/k))
        R_prime = 1 / log(x)  # Approximate derivative
        x = x + (mpf(n) - R_val) / R_prime
        if abs(mpf(n) - R_val) < mpf('1e-10'):
            break
    return float(x)

def compute_delta_sequence(n_start, n_end):
    """Compute delta(n) = p(n) - R^{-1}(n) for a range of n."""
    deltas = []
    for n in range(n_start, n_end + 1):
        p_n = int(prime(n))
        r_inv = R_inv_approx(n)
        delta = p_n - r_inv
        deltas.append((n, p_n, r_inv, delta))
    return deltas

def analyze_delta_structure(deltas):
    """Analyze the structure of delta(n) for compressibility."""
    delta_vals = np.array([float(d[3]) for d in deltas])
    n_vals = np.array([d[0] for d in deltas])

    # 1. Basic statistics
    print(f"  Range: n={n_vals[0]}..{n_vals[-1]}")
    print(f"  delta mean={np.mean(delta_vals):.4f}, std={np.std(delta_vals):.4f}")
    print(f"  delta range: [{np.min(delta_vals):.4f}, {np.max(delta_vals):.4f}]")

    # 2. First differences (related to prime gaps)
    diffs = np.diff(delta_vals)
    print(f"  delta' mean={np.mean(diffs):.4f}, std={np.std(diffs):.4f}")

    # 3. FFT analysis - is delta sparse in frequency domain?
    fft_vals = np.fft.fft(delta_vals)
    magnitudes = np.abs(fft_vals)
    sorted_mags = np.sort(magnitudes)[::-1]
    total_energy = np.sum(magnitudes**2)

    cumulative = np.cumsum(sorted_mags**2) / total_energy
    k_90 = np.searchsorted(cumulative, 0.90) + 1
    k_95 = np.searchsorted(cumulative, 0.95) + 1
    k_99 = np.searchsorted(cumulative, 0.99) + 1
    N = len(delta_vals)

    print(f"  FFT sparsity: 90%@{k_90}/{N}, 95%@{k_95}/{N}, 99%@{k_99}/{N}")

    # 4. Wavelet-like analysis via Haar coefficients
    # Haar transform at multiple scales
    haar_coeffs = []
    vals = delta_vals.copy()
    level = 0
    while len(vals) > 1:
        n_pairs = len(vals) // 2
        details = [(vals[2*i+1] - vals[2*i]) / np.sqrt(2) for i in range(n_pairs)]
        averages = [(vals[2*i] + vals[2*i+1]) / np.sqrt(2) for i in range(n_pairs)]
        haar_coeffs.extend([(level, abs(d)) for d in details])
        vals = np.array(averages)
        level += 1

    haar_mags = sorted([c[1] for c in haar_coeffs], reverse=True)
    haar_total = sum(m**2 for m in haar_mags)
    haar_cum = np.cumsum(np.array(haar_mags)**2) / haar_total if haar_total > 0 else [0]
    haar_k_90 = np.searchsorted(haar_cum, 0.90) + 1
    print(f"  Haar sparsity: 90%@{haar_k_90}/{len(haar_mags)} coefficients")

    return delta_vals

def binary_search_with_approximate_oracle(n_target, tol=0.5):
    """
    Binary search for p(n) using approximate pi(x).

    At each step, compute R(x_mid) as approximation to pi(x_mid).
    Track the error and see when the approximation fails to guide the search.
    """
    p_n = int(prime(n_target))

    # Binary search bounds
    lo = max(2, int(n_target * math.log(n_target) * 0.8))
    hi = int(n_target * (math.log(n_target) + math.log(math.log(n_target))) * 1.2)

    steps = 0
    oracle_correct = 0
    oracle_total = 0

    while lo < hi:
        mid = (lo + hi) // 2
        steps += 1

        # Exact oracle
        exact_pi = int(primepi(mid))

        # Approximate oracle via R(mid)
        from sympy import mobius
        R_val = mpf(0)
        for k in range(1, 30):
            mu_k = mobius(k)
            if mu_k != 0:
                R_val += mpf(mu_k) / k * li(mpf(mid) ** (mpf(1)/k))
        approx_pi = round(float(R_val))

        oracle_total += 1
        exact_direction = exact_pi >= n_target
        approx_direction = approx_pi >= n_target

        if exact_direction == approx_direction:
            oracle_correct += 1

        # Use exact for actual search
        if exact_pi >= n_target:
            hi = mid
        else:
            lo = mid + 1

    accuracy = oracle_correct / oracle_total if oracle_total > 0 else 0
    return {
        'n': n_target,
        'p_n': p_n,
        'found': lo,
        'correct': lo == p_n,
        'steps': steps,
        'oracle_accuracy': accuracy,
        'oracle_correct': oracle_correct,
        'oracle_total': oracle_total,
    }

def test_delta_autocorrelation(deltas, max_lag=50):
    """Test if delta(n) has significant autocorrelation structure."""
    delta_vals = np.array([float(d[3]) for d in deltas])
    delta_vals = delta_vals - np.mean(delta_vals)

    if np.std(delta_vals) == 0:
        return []

    N = len(delta_vals)
    autocorr = []
    for lag in range(1, min(max_lag + 1, N)):
        c = np.corrcoef(delta_vals[:-lag], delta_vals[lag:])[0, 1]
        autocorr.append((lag, c))

    return autocorr


if __name__ == '__main__':
    print("=" * 80)
    print("PROPOSAL 3: Compressed Sensing on Delta(n)")
    print("=" * 80)

    # Use powers of 2 for clean Haar analysis
    print("\n--- Part A: Delta structure analysis (n=1..256) ---")
    deltas_small = compute_delta_sequence(1, 256)
    delta_vals = analyze_delta_structure(deltas_small)

    print("\n--- Part B: Delta structure analysis (n=500..756) ---")
    deltas_mid = compute_delta_sequence(500, 756)
    analyze_delta_structure(deltas_mid)

    print("\n--- Part C: Autocorrelation of delta(n) ---")
    autocorr = test_delta_autocorrelation(deltas_small, max_lag=30)
    significant = [(lag, c) for lag, c in autocorr if abs(c) > 0.1]
    print(f"  Significant autocorrelations (|r|>0.1): {len(significant)}/{len(autocorr)}")
    for lag, c in significant[:10]:
        print(f"    lag={lag:3d}: r={c:+.4f}")

    print("\n--- Part D: Binary search with approximate oracle ---")
    for n in [100, 500, 1000, 2000, 5000]:
        result = binary_search_with_approximate_oracle(n)
        print(f"  n={n:5d}: found={result['found']}, correct={result['correct']}, "
              f"steps={result['steps']}, oracle_accuracy={result['oracle_accuracy']:.3f}")

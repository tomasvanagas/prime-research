"""
Convergence Acceleration for the Riemann Explicit Formula
==========================================================
Tests multiple techniques to accelerate the conditionally convergent series:
  pi(x) = R(x) - sum_{k=1}^{N} 2*Re(R(x^{rho_k}))

Techniques:
  1. Windowed summation (Hanning, Gaussian, Lanczos)
  2. Cesaro summation
  3. Richardson extrapolation
  4. Euler/van Wijngaarden transform
  5. Shanks transformation (Wynn epsilon algorithm)
"""

import mpmath
import math
import time
import os
import sys

mpmath.mp.dps = 30

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# Import from existing module
from riemann_explicit import R_function, R_at_rho, MOBIUS


def load_zeros(count):
    """Load precomputed zeros, trying largest file available."""
    for try_count in [1000, 500, 300, 200]:
        fname = os.path.join(SCRIPT_DIR, f"zeta_zeros_{try_count}.txt")
        if os.path.exists(fname):
            with open(fname) as f:
                zeros = [line.strip() for line in f if line.strip()]
            if len(zeros) >= count:
                return zeros[:count]
    raise FileNotFoundError(f"Cannot find enough zeros (need {count})")


def compute_individual_terms(x, zeros):
    """
    Compute each term t_k = 2*Re(R(x^{rho_k})) for all zeros.
    Returns list of floats.
    """
    x_mpf = mpmath.mpf(x)
    terms = []
    for gs in zeros:
        gamma = mpmath.mpf(gs)
        rho = mpmath.mpc(0.5, gamma)
        val = 2 * float(mpmath.re(R_at_rho(x_mpf, rho)))
        terms.append(val)
    return terms


def plain_sum(Rx, terms, N):
    """Plain truncated sum using first N terms."""
    return Rx - sum(terms[:N])


# ============================================================
# 1. Windowed summation
# ============================================================

def window_hanning(terms, N):
    """Hanning window: w(t) = 0.5*(1 + cos(pi*t)) for t = k/N."""
    result = 0.0
    for k in range(N):
        t = (k + 1) / N  # t in (0, 1]
        w = 0.5 * (1 + math.cos(math.pi * t))
        result += w * terms[k]
    return result


def window_gaussian(terms, N, alpha=3.0):
    """Gaussian window: w(t) = exp(-alpha*t^2)."""
    result = 0.0
    for k in range(N):
        t = (k + 1) / N
        w = math.exp(-alpha * t * t)
        result += w * terms[k]
    return result


def window_lanczos(terms, N):
    """Lanczos sigma factor: w(t) = sinc(t) = sin(pi*t)/(pi*t)."""
    result = 0.0
    for k in range(N):
        t = (k + 1) / N
        if t < 1e-15:
            w = 1.0
        else:
            w = math.sin(math.pi * t) / (math.pi * t)
        result += w * terms[k]
    return result


def windowed_sum(Rx, terms, N, window_func, **kwargs):
    """Apply windowed summation."""
    return Rx - window_func(terms, N, **kwargs)


# ============================================================
# 2. Cesaro summation
# ============================================================

def cesaro_sum(Rx, terms, N):
    """
    Cesaro mean: C_N = (1/N) * sum_{k=1}^{N} S_k
    where S_k = sum_{j=1}^{k} t_j
    Equivalent to weighting term j by (N - j + 1) / N (Fejer kernel).
    """
    # More efficient: term j gets weight (N - j) / N  (0-indexed)
    weighted = 0.0
    for j in range(N):
        weight = (N - j) / N
        weighted += weight * terms[j]
    return Rx - weighted


# ============================================================
# 3. Richardson extrapolation
# ============================================================

def richardson_extrapolation(Rx, terms, N1, N2, alpha=1):
    """
    Given S(N1) and S(N2) with error ~ C/N^alpha:
      S_extrap = (N2^a * S(N2) - N1^a * S(N1)) / (N2^a - N1^a)
    """
    S1 = plain_sum(Rx, terms, N1)
    S2 = plain_sum(Rx, terms, N2)
    n1a = N1 ** alpha
    n2a = N2 ** alpha
    return (n2a * S2 - n1a * S1) / (n2a - n1a)


def richardson_3point(Rx, terms, N1, N2, N3, alpha=1):
    """Three-point Richardson using N1, N2, N3."""
    # First eliminate leading error between (N1,N2) and (N2,N3)
    S12 = richardson_extrapolation(Rx, terms, N1, N2, alpha)
    S23 = richardson_extrapolation(Rx, terms, N2, N3, alpha)
    # Then extrapolate again with alpha+1
    # Treat S12 as "at N2" and S23 as "at N3" with error ~ C/N^(alpha+1)
    n2a = N2 ** (alpha + 1)
    n3a = N3 ** (alpha + 1)
    return (n3a * S23 - n2a * S12) / (n3a - n2a)


# ============================================================
# 4. Euler / van Wijngaarden acceleration
# ============================================================

def euler_transform(terms, N):
    """
    Euler transform for oscillating series.
    Compute forward differences of the tail terms and apply binomial weighting.
    Returns the accelerated sum of terms[0:N].
    """
    # Build partial sums sequence
    partial = []
    s = 0.0
    for k in range(N):
        s += terms[k]
        partial.append(s)

    # Apply Euler transform to the partial sums sequence
    # E_n = sum_{k=0}^{n} C(n,k) * S_k / 2^(n+1)
    # Use the last chunk of partial sums
    m = min(N, 20)  # use last m partial sums for transform
    seq = partial[N - m:]

    result = 0.0
    for k in range(m):
        binom = math.comb(m - 1, k)
        result += binom * seq[k]
    result /= 2 ** (m - 1)
    return result


def van_wijngaarden(Rx, terms, N):
    """Apply Euler transform to accelerate the zero-sum."""
    accelerated_sum = euler_transform(terms, N)
    return Rx - accelerated_sum


# ============================================================
# 5. Shanks / Wynn epsilon algorithm
# ============================================================

def wynn_epsilon(partial_sums):
    """
    Wynn's epsilon algorithm (Shanks transformation).
    Given a sequence of partial sums, returns the best estimate.

    The epsilon table:
      eps_{-1}(n) = 0
      eps_0(n) = S_n
      eps_{k+1}(n) = eps_{k-1}(n+1) + 1/(eps_k(n+1) - eps_k(n))

    The even-indexed columns eps_{2k}(0) are the Shanks transforms.
    """
    n = len(partial_sums)
    if n < 3:
        return partial_sums[-1]

    # Build epsilon table
    # eps[k][j] where k is column, j is row
    eps = [[0.0] * n for _ in range(n + 1)]

    # Column 0 = partial sums
    for j in range(n):
        eps[0][j] = partial_sums[j]

    # Fill columns
    max_col = n - 1
    for k in range(1, max_col + 1):
        for j in range(n - k):
            diff = eps[k - 1][j + 1] - eps[k - 1][j]
            if abs(diff) < 1e-50:
                # Avoid division by zero; keep previous
                eps[k][j] = eps[k - 1][j + 1]
            else:
                if k == 1:
                    eps[k][j] = 1.0 / diff
                else:
                    eps[k][j] = eps[k - 2][j + 1] + 1.0 / diff

    # Return the highest even-column entry at row 0
    best = partial_sums[-1]
    for k in range(0, max_col + 1, 2):
        if k < len(eps) and len(eps[k]) > 0:
            best = eps[k][0]
    return best


def shanks_sum(Rx, terms, N):
    """Apply Wynn epsilon (Shanks) to partial sums of the zero-sum."""
    partial = []
    s = 0.0
    for k in range(N):
        s += terms[k]
        partial.append(s)

    accelerated = wynn_epsilon(partial)
    return Rx - accelerated


# ============================================================
# Combined: Shanks on windowed partial sums
# ============================================================

def shanks_on_windowed(Rx, terms, N, window_func):
    """Apply Shanks to a sequence of windowed partial sums at varying N."""
    # Build sequence of windowed estimates for N/4, N/3, N/2, 2N/3, 3N/4, N
    fractions = [0.25, 0.33, 0.4, 0.5, 0.6, 0.67, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0]
    estimates = []
    for frac in fractions:
        n = max(5, int(frac * N))
        val = windowed_sum(Rx, terms, n, window_func)
        estimates.append(val)

    if len(estimates) < 3:
        return estimates[-1]
    return wynn_epsilon(estimates)


# ============================================================
# Combined: Richardson on windowed sums
# ============================================================

def richardson_windowed(Rx, terms, N, window_func, alpha=1):
    """Richardson extrapolation applied to windowed sums at N/2 and N."""
    N1 = N // 2
    N2 = N
    S1 = windowed_sum(Rx, terms, N1, window_func)
    S2 = windowed_sum(Rx, terms, N2, window_func)
    n1a = N1 ** alpha
    n2a = N2 ** alpha
    return (n2a * S2 - n1a * S1) / (n2a - n1a)


# ============================================================
# Run all experiments
# ============================================================

EXACT = {100: 25, 1000: 168, 10000: 1229, 100000: 9592, 1000000: 78498}
X_VALUES = [100, 1000, 10000, 100000, 1000000]
N_VALUES = [50, 100, 200, 300, 500]


# ============================================================
# 6. Tail correction: analytic estimate of truncated tail
# ============================================================

def tail_correction(x, T_last, num_zeros_used):
    """
    Estimate the tail of the zero sum beyond height T_last.

    The density of zeros at height T is ~ (1/(2*pi)) * log(T/(2*pi)).
    Each zero rho = 1/2 + i*gamma contributes ~ 2*Re(R(x^rho)).

    For large gamma, R(x^rho) ~ li(x^rho) ~ x^rho / (rho * ln(x)).
    So 2*Re(R(x^rho)) ~ 2*x^{1/2} * cos(gamma*ln(x)) / (|rho| * ln(x)).

    The truncated tail sum from T_last to infinity can be estimated using
    the integral of the density times the term magnitude.

    A simpler approach: use the explicit formula remainder.
    The truncation error at height T is approximately:
      err ~ (x^{1/2} / pi) * sum of 1/gamma for gamma > T_last
         ~ (x^{1/2} / pi) * log(T_last) / T_last  (using Euler-Maclaurin)

    But empirically, the sign and magnitude depend on where we cut.
    """
    ln_x = math.log(x)
    sqrt_x = math.sqrt(x)

    # Empirical: estimate contribution of zeros beyond T_last
    # Using the integral approximation of the oscillatory sum
    # The "DC component" of the tail can be approximated
    # via the argument principle / counting function N(T)

    # More practical approach: use the last few terms to extrapolate the tail
    return 0.0  # placeholder - the analytic estimate is complex


def smoothed_explicit_formula(x, zeros_str, N):
    """
    Use the smoothed explicit formula with a compactly-supported test function.

    Instead of the sharp cutoff pi(x) which requires infinitely many zeros,
    use a smoothed version that converges faster.

    The key idea: convolve pi(x) with a bump function of width delta.
    This makes the zero sum converge like exp(-gamma*delta) instead of 1/gamma.

    Specifically, use the integrated Chebyshev function:
      psi_1(x) = integral_0^x psi(t) dt = x^2/2 - sum_rho x^{rho+1}/(rho*(rho+1)) - ...

    The sum over zeros converges as 1/gamma^2 instead of 1/gamma.
    Then recover pi(x) by differencing.
    """
    # This is complex to implement properly; placeholder
    return None


def tail_integral_correction(x, terms, N, zeros_str):
    """
    Estimate the tail of the sum beyond the N-th zero using
    the Riemann-von Mangoldt formula for zero density.

    The N(T) zeros up to height T satisfy N(T) ~ T/(2*pi) * log(T/(2*pi)) - T/(2*pi).

    For large gamma_k, the term 2*Re(R(x^{rho_k})) oscillates with decreasing envelope.
    The *average* magnitude of the k-th term is approximately:
      |t_k| ~ 2 * x^{1/2} / (gamma_k * ln(x))

    But the SUM of oscillating terms from N+1 to infinity has a much smaller net value.
    The systematic offset comes from the non-oscillatory part of each R(x^rho) evaluation.

    A key correction: the -1/(ln(x)) term from R(x) should include the
    contribution from the trivial zeros rho = -2, -4, -6, ...
    """
    ln_x = math.log(x)
    sqrt_x = math.sqrt(x)

    # Correction from trivial zeros: sum_{n=1}^{inf} R(x^{-2n})
    # For large x: R(x^{-2n}) ~ li(x^{-2n}) which is negligible for x >> 1
    # Actually for x^{-2n} < 1, li(x^{-2n}) is complex

    # The standard "other terms" in the explicit formula:
    # -ln(2) + integral_x^inf dt / (t*(t^2-1)*ln(t))
    # This integral is tiny for x >= 2.

    # Correction: the integral term
    if x > 2:
        # integral from x to inf of 1/(t*(t^2-1)*ln(t)) dt
        # For large x, this is approximately 1/(x^2 * ln(x))
        integral_correction = 0.0
        t = float(x)
        # Numerical integration (simple)
        dt = 0.1
        for _ in range(1000):
            if t > 1e10:
                break
            integral_correction += dt / (t * (t * t - 1) * math.log(t))
            t += dt
            dt *= 1.01
    else:
        integral_correction = 0.0

    return integral_correction


def polynomial_extrapolation(Rx, terms, N_values_and_results):
    """
    Given plain_sum results at multiple N values, fit a polynomial in 1/N
    and extrapolate to 1/N -> 0 (i.e., N -> infinity).
    """
    import numpy as np

    Ns = [n for n, _ in N_values_and_results]
    vals = [v for _, v in N_values_and_results]

    # Fit: S(N) = S_inf + a1/N + a2/N^2 + ...
    # Use 1/N as variable, extrapolate to 0
    xs = [1.0 / n for n in Ns]
    coeffs = np.polyfit(xs, vals, min(len(xs) - 1, 3))
    # Value at x=0 is the last coefficient
    return coeffs[-1]


def multi_N_extrapolation(Rx, terms, N_max):
    """Compute plain sums at several N values and extrapolate to N->inf."""
    try:
        import numpy as np
    except ImportError:
        return plain_sum(Rx, terms, N_max)

    # Use several sub-values
    fractions = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    pairs = []
    for f in fractions:
        n = max(10, int(f * N_max))
        s = plain_sum(Rx, terms, n)
        pairs.append((n, s))

    return polynomial_extrapolation(Rx, terms, pairs)


def aitken_delta2(Rx, terms, N):
    """
    Aitken's delta-squared process applied to the sequence of partial sums
    at N-2, N-1, N.
    If S_n -> L with error ~ C*r^n, then:
      L ~ S_n - (S_{n+1} - S_n)^2 / (S_{n+2} - 2*S_{n+1} + S_n)
    Apply this repeatedly to the partial sum sequence.
    """
    # Build partial sums of zero terms
    partial = []
    s = 0.0
    for k in range(N):
        s += terms[k]
        partial.append(Rx - s)

    # Apply Aitken repeatedly
    seq = partial
    for _ in range(min(5, len(seq) // 3)):
        new_seq = []
        for i in range(len(seq) - 2):
            denom = seq[i + 2] - 2 * seq[i + 1] + seq[i]
            if abs(denom) < 1e-30:
                new_seq.append(seq[i + 2])
            else:
                new_seq.append(seq[i] - (seq[i + 1] - seq[i]) ** 2 / denom)
        if len(new_seq) < 3:
            break
        seq = new_seq

    return seq[-1] if seq else partial[-1]


def run_experiments():
    """Run all convergence acceleration experiments."""

    max_zeros_needed = 500
    zeros_str = load_zeros(max_zeros_needed)

    results = {}  # results[(technique, x, N)] = (computed, rounded, exact_match)

    for x_val in X_VALUES:
        exact = EXACT[x_val]
        print(f"\n{'='*80}")
        print(f"  x = {x_val:,}    exact pi(x) = {exact}")
        print(f"{'='*80}")

        t0 = time.time()
        Rx = float(R_function(x_val))
        print(f"  R(x) = {Rx:.6f}  (computed in {time.time()-t0:.2f}s)")

        t0 = time.time()
        all_terms = compute_individual_terms(x_val, zeros_str[:max_zeros_needed])
        print(f"  Computed {max_zeros_needed} zero terms in {time.time()-t0:.2f}s")

        # --- Header ---
        header = f"  {'Technique':<35s}"
        for N in N_VALUES:
            header += f"  N={N:>3d}"
        print(header)
        print("  " + "-" * (35 + 8 * len(N_VALUES)))

        techniques = [
            ("Plain (no accel)", lambda Rx, terms, N: plain_sum(Rx, terms, N)),
            ("Hanning window", lambda Rx, terms, N: windowed_sum(Rx, terms, N, window_hanning)),
            ("Gaussian window (a=3)", lambda Rx, terms, N: windowed_sum(Rx, terms, N, window_gaussian, alpha=3.0)),
            ("Gaussian window (a=6)", lambda Rx, terms, N: windowed_sum(Rx, terms, N, window_gaussian, alpha=6.0)),
            ("Lanczos sigma", lambda Rx, terms, N: windowed_sum(Rx, terms, N, window_lanczos)),
            ("Cesaro mean", lambda Rx, terms, N: cesaro_sum(Rx, terms, N)),
            ("Richardson (a=1, N/2,N)", lambda Rx, terms, N: richardson_extrapolation(Rx, terms, N // 2, N, alpha=1)),
            ("Richardson (a=2, N/2,N)", lambda Rx, terms, N: richardson_extrapolation(Rx, terms, N // 2, N, alpha=2)),
            ("Richardson 3pt (a=1)", lambda Rx, terms, N: richardson_3point(Rx, terms, N // 4, N // 2, N, alpha=1)),
            ("Euler/van Wijngaarden", lambda Rx, terms, N: van_wijngaarden(Rx, terms, N)),
            ("Shanks (Wynn epsilon)", lambda Rx, terms, N: shanks_sum(Rx, terms, N)),
            ("Shanks+Hanning", lambda Rx, terms, N: shanks_on_windowed(Rx, terms, N, window_hanning)),
            ("Shanks+Gaussian(3)", lambda Rx, terms, N: shanks_on_windowed(Rx, terms, N, window_gaussian)),
            ("Richardson+Hanning(a=1)", lambda Rx, terms, N: richardson_windowed(Rx, terms, N, window_hanning, alpha=1)),
            ("Aitken delta^2", lambda Rx, terms, N: aitken_delta2(Rx, terms, N)),
            ("Poly extrap (1/N)", lambda Rx, terms, N: multi_N_extrapolation(Rx, terms, N)),
        ]

        for name, func in techniques:
            row = f"  {name:<35s}"
            for N in N_VALUES:
                try:
                    val = func(Rx, all_terms, N)
                    rounded = round(val)
                    err = val - exact
                    match = (rounded == exact)
                    results[(name, x_val, N)] = (val, rounded, match)
                    if match:
                        row += f"  {'OK':>6s}"
                    else:
                        row += f"  {err:>+6.1f}"
                except Exception as e:
                    row += f"  {'ERR':>6s}"
            print(row)

    return results


def summary_analysis(results):
    """Analyze which techniques give best results."""
    print("\n" + "=" * 80)
    print("SUMMARY: Which techniques make rounded result exact?")
    print("=" * 80)

    # For each technique, count how many (x, N) pairs give exact results
    from collections import defaultdict
    tech_scores = defaultdict(lambda: {"exact": 0, "total": 0})

    for (tech, x_val, N), (val, rounded, match) in results.items():
        tech_scores[tech]["total"] += 1
        if match:
            tech_scores[tech]["exact"] += 1

    print(f"\n  {'Technique':<35s}  {'Exact/Total':>12s}  {'Rate':>6s}")
    print("  " + "-" * 56)
    for tech in sorted(tech_scores, key=lambda t: -tech_scores[t]["exact"]):
        s = tech_scores[tech]
        rate = s["exact"] / s["total"] * 100
        print(f"  {tech:<35s}  {s['exact']:>4d}/{s['total']:<4d}     {rate:>5.1f}%")

    # Key question: does accelerated-300 beat plain-500?
    print("\n" + "=" * 80)
    print("KEY QUESTION: Does accelerated-300-zeros beat plain-500-zeros?")
    print("=" * 80)

    for x_val in X_VALUES:
        exact = EXACT[x_val]
        plain_500 = results.get(("Plain (no accel)", x_val, 500))
        if plain_500 is None:
            continue
        p500_err = abs(plain_500[0] - exact)

        print(f"\n  x = {x_val:>10,}  (exact = {exact})")
        print(f"    Plain 500 zeros: error = {p500_err:.4f}, exact = {plain_500[2]}")

        for tech in sorted(results):
            name, xv, N = tech
            if xv != x_val or N != 300 or name == "Plain (no accel)":
                continue
            val, rounded, match = results[tech]
            err = abs(val - exact)
            beats = err < p500_err
            print(f"    {name:<35s} 300z: err={err:.4f} {'BEATS' if beats else 'loses'} plain-500 (exact={match})")

    # The big question: pi(10^6) with 500 zeros?
    print("\n" + "=" * 80)
    print("THE BIG QUESTION: Can we get exact pi(10^6) with 500 zeros?")
    print("=" * 80)

    for tech in sorted(results):
        name, x_val, N = tech
        if x_val != 1000000 or N != 500:
            continue
        val, rounded, match = results[tech]
        err = abs(val - EXACT[1000000])
        status = "EXACT!" if match else f"off by {rounded - EXACT[1000000]:+d}"
        print(f"  {name:<35s}  N=500: {val:.4f}  rounded={rounded}  {status}  (err={err:.2f})")


def main():
    print("Convergence Acceleration for Riemann Explicit Formula")
    print("=" * 80)

    results = run_experiments()
    summary_analysis(results)


if __name__ == "__main__":
    main()

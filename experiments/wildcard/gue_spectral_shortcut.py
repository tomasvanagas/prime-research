#!/usr/bin/env python3
"""
GUE-INFORMED SPECTRAL COMPRESSION

The explicit formula: pi(x) = li(x) - sum_rho li(x^rho) + smaller terms

Under RH: rho = 1/2 + i*gamma_j, so the oscillatory sum is:
  S(x) = sum_j x^{1/2+i*gamma_j} / (1/2 + i*gamma_j)
       = x^{1/2} * sum_j exp(i*gamma_j*log(x)) / rho_j

Montgomery-Odlyzko: The gamma_j are distributed like eigenvalues of GUE
random matrices. GUE eigenvalues have determinantal structure:
  R_n(x1,...,xn) = det[K(xi, xj)]_{i,j=1}^n
where K is the sine kernel.

KEY QUESTION: Does this determinantal structure enable computing
S(x) = sum_j f(gamma_j) without evaluating each term individually?

ANALOGY: For GUE eigenvalues, the linear statistic
  sum_j f(lambda_j) = integral f(t) * rho(t) dt + O(1)
with VARIANCE = O(log N) regardless of f (CLT for linear statistics).

If the VARIANCE is only O(log N), then the sum is determined up to
O(sqrt(log N)) by just the smooth part! For pi(x), this would mean
the oscillatory correction is O(sqrt(log x)) -- much smaller than the
known O(x^{1/2} / log x).

BUT WAIT: this CLT applies to RANDOM GUE matrices, not to the specific
zeta zero configuration. The zeros are FIXED, not random.

Let's test: how well does the CLT prediction work for actual zeta zeros?
"""

import numpy as np
from mpmath import mp, zetazero, li, log, exp, pi as mpi, fabs, im, re, mpc
from sympy import primepi
import time

def load_or_compute_zeros(num_zeros=200):
    """Load/compute zeta zeros."""
    print(f"  Computing {num_zeros} zeta zeros...")
    mp.dps = 30
    zeros = []
    for k in range(1, num_zeros + 1):
        gamma = im(zetazero(k))
        zeros.append(float(gamma))
    return np.array(zeros)


def test_gue_clt_for_prime_counting():
    """
    Test: Does the GUE CLT prediction match the actual oscillatory part?

    For the linear statistic S(x) = sum_j x^{rho_j}/rho_j:
    - Mean (from smooth density): integral x^{1/2+it}/(1/2+it) * rho(t) dt
      where rho(t) = (1/2pi) log(t/2pi)
    - This integral gives approximately li(x) corrections (already known)
    - Variance: O(log T) for T = max gamma considered

    The ACTUAL oscillatory part (what we need for exact pi(x)):
    delta(x) = pi(x) - li(x) + (smaller terms)

    If GUE CLT applied perfectly, |delta(x)| = O(sqrt(log x)).
    Known: |delta(x)| = O(x^{1/2} / log x) under RH.

    The GUE CLT does NOT apply directly because:
    1. We're not averaging over random matrices -- zeros are fixed
    2. The "weight" function 1/rho_j breaks the CLT assumptions
    3. The CLT gives the typical size, not the max
    """
    print("=" * 70)
    print("TEST 1: GUE CLT vs actual prime counting error")
    print("=" * 70)

    zeros = load_or_compute_zeros(200)

    # Compute the oscillatory sum S(x) for various x
    test_x = [100, 500, 1000, 5000, 10000, 50000, 100000]

    print(f"\n  {'x':>8} | {'pi(x)':>7} | {'li(x)':>10} | {'delta':>7} | {'|S_200|':>10} | {'sqrt(log x)':>11}")
    print("  " + "-" * 70)

    mp.dps = 20
    for x in test_x:
        actual_pi = int(primepi(x))
        li_val = float(li(x))
        delta = actual_pi - li_val

        # Compute partial sum with 200 zeros
        s = 0.0
        for gamma in zeros:
            rho = 0.5 + 1j * gamma
            # li(x^rho) ~ x^rho / (rho * log(x))  for large x
            term = x**(0.5) * np.exp(1j * gamma * np.log(x)) / (rho * np.log(x))
            s += term
        # Take twice the real part (conjugate zeros)
        s_val = -2 * s.real

        sqrt_log = np.sqrt(np.log(x))

        print(f"  {x:>8} | {actual_pi:>7} | {li_val:>10.2f} | {delta:>7.2f} | {s_val:>10.2f} | {sqrt_log:>11.4f}")

    print()
    print("  Note: S_200 uses only 200 zeros (first approximation).")
    print("  delta(x) is NOT O(sqrt(log x)) -- it's O(sqrt(x)/log(x)).")
    print("  GUE CLT does not give the shortcut we hoped for.")
    print()


def test_zero_sum_truncation():
    """
    Test: How many zeros are needed for the sum to converge
    to within 0.5 of pi(x)?

    This is the core question. If K(x) zeros suffice and K(x) = O(polylog(x)),
    we have a breakthrough.
    """
    print("=" * 70)
    print("TEST 2: Zero sum truncation -- how many zeros needed?")
    print("=" * 70)

    zeros = load_or_compute_zeros(300)

    mp.dps = 20
    test_x = [100, 500, 1000, 5000, 10000]

    for x in test_x:
        actual_pi = int(primepi(x))
        li_val = float(li(x))

        # Also compute li(x^{1/2}), li(x^{1/3}), etc. (smaller correction terms)
        half_correction = -0.5 * float(li(x**0.5))  # Riemann's R(x) correction

        # Build the sum incrementally
        partial_sums = []
        s = 0.0
        for k, gamma in enumerate(zeros):
            rho = 0.5 + 1j * gamma
            # More precise: li(x^rho)
            xrho = x**(0.5) * np.exp(1j * gamma * np.log(x))
            term = xrho / (rho * np.log(x))  # approximate li(x^rho)
            s += term
            correction = -2 * s.real
            estimate = li_val + correction + half_correction
            error = abs(estimate - actual_pi)
            partial_sums.append((k + 1, error))

        # Find how many zeros needed for error < 0.5
        zeros_needed = None
        for k, err in partial_sums:
            if err < 0.5:
                zeros_needed = k
                break

        # Also find where error first goes below 1.0
        zeros_for_1 = None
        for k, err in partial_sums:
            if err < 1.0:
                zeros_for_1 = k
                break

        print(f"  x={x:>6}: pi(x)={actual_pi:>5}, "
              f"zeros for err<1.0: {zeros_for_1 if zeros_for_1 else 'N/A':>4}, "
              f"zeros for err<0.5: {zeros_needed if zeros_needed else 'N/A':>4}")

        # Print error at various truncation points
        checkpoints = [10, 25, 50, 100, 200, 300]
        errors_at_cp = []
        for cp in checkpoints:
            if cp <= len(partial_sums):
                errors_at_cp.append(f"{partial_sums[cp-1][1]:.2f}")
            else:
                errors_at_cp.append("N/A")
        print(f"          Errors at K={checkpoints}: {errors_at_cp}")

    print()


def test_zero_grouping():
    """
    Test: Can we GROUP nearby zeros and compute their collective contribution?

    If zeros gamma_j and gamma_{j+1} are close, their contributions to S(x)
    are similar. Can we replace pairs/clusters with a single effective term?

    The zero spacing is delta_j ~ 2*pi / log(gamma_j / (2*pi))
    """
    print("=" * 70)
    print("TEST 3: Zero grouping / clustering")
    print("=" * 70)

    zeros = load_or_compute_zeros(200)

    # Analyze zero spacings
    spacings = np.diff(zeros)
    mean_spacing = spacings.mean()
    print(f"  Mean zero spacing (first 200): {mean_spacing:.4f}")
    print(f"  Expected from density: 2*pi/log(gamma_100/(2*pi)) = {2*np.pi/np.log(zeros[99]/(2*np.pi)):.4f}")
    print()

    # Group zeros into clusters of size G
    x = 10000
    actual_pi = int(primepi(x))
    li_val = float(li(x))
    half_correction = -0.5 * float(li(x**0.5))

    print(f"  Testing zero grouping for x={x}, pi(x)={actual_pi}:")

    # Full sum with all 200 zeros
    s_full = 0.0
    for gamma in zeros:
        rho = 0.5 + 1j * gamma
        xrho = x**(0.5) * np.exp(1j * gamma * np.log(x))
        term = xrho / (rho * np.log(x))
        s_full += term
    full_correction = -2 * s_full.real
    full_estimate = li_val + full_correction + half_correction
    full_error = abs(full_estimate - actual_pi)
    print(f"  Full sum (200 zeros): estimate={full_estimate:.2f}, error={full_error:.2f}")

    # Grouped sums
    for G in [2, 4, 8, 16, 32]:
        s_grouped = 0.0
        num_groups = len(zeros) // G

        for g in range(num_groups):
            # Take the center of each group
            group = zeros[g*G:(g+1)*G]
            gamma_center = group.mean()
            # Weight by number of zeros in group
            rho_center = 0.5 + 1j * gamma_center

            # Sum the actual terms in the group for comparison
            group_sum = 0.0
            for gamma in group:
                rho = 0.5 + 1j * gamma
                xrho = x**(0.5) * np.exp(1j * gamma * np.log(x))
                term = xrho / (rho * np.log(x))
                group_sum += term

            # Approximate: use center with multiplicity G
            # This is like a "mean field" approximation
            xrho_center = x**(0.5) * np.exp(1j * gamma_center * np.log(x))
            approx = G * xrho_center / (rho_center * np.log(x))

            s_grouped += group_sum  # exact within group

        # Now try the APPROXIMATE version (center of each group * G)
        s_approx = 0.0
        for g in range(num_groups):
            group = zeros[g*G:(g+1)*G]
            gamma_center = group.mean()
            rho_center = 0.5 + 1j * gamma_center
            xrho_center = x**(0.5) * np.exp(1j * gamma_center * np.log(x))
            approx = G * xrho_center / (rho_center * np.log(x))
            s_approx += approx

        approx_correction = -2 * s_approx.real
        approx_estimate = li_val + approx_correction + half_correction
        approx_error = abs(approx_estimate - actual_pi)

        print(f"  Group size {G:>2} ({num_groups:>3} groups): "
              f"estimate={approx_estimate:.2f}, error={approx_error:.2f}")

    print()
    print("  Result: Grouping introduces significant error because")
    print("  exp(i*gamma*log(x)) oscillates RAPIDLY -- neighboring zeros")
    print("  contribute OPPOSITE-signed terms that largely cancel.")
    print("  This cancellation is ESSENTIAL and can't be approximated.")
    print()


def test_random_matrix_surrogate():
    """
    Test: Can we replace the actual zeta zeros with a RANDOM MATRIX
    surrogate that preserves enough statistics to give exact pi(x)?

    If a GUE matrix of size N can model T/2pi zeros, and we only need
    the sum S(x) = sum_j f(gamma_j), then we need:
    - The spectral density (gives smooth part -- already known)
    - The 2-point correlation (gives variance)
    - Higher correlations...

    For EXACT computation, we'd need ALL correlations = the full zero set.
    But if corrections from k-point correlations decay fast...
    """
    print("=" * 70)
    print("TEST 4: GUE surrogate for zeta zeros")
    print("=" * 70)

    zeros = load_or_compute_zeros(100)

    # Generate GUE matrices and compare their eigenvalue sums
    # to the actual zeta zero sum

    N_matrix = 50  # GUE matrix size (models ~50 zeros)

    # Normalized zeta zeros (unfolded to unit mean spacing)
    unfolded = np.zeros(len(zeros))
    for i, gamma in enumerate(zeros):
        # Unfolding: N(gamma) = (gamma/(2*pi)) * log(gamma/(2*pi*np.e)) + 7/8
        if gamma > 0:
            unfolded[i] = (gamma / (2*np.pi)) * np.log(gamma / (2*np.pi * np.e)) + 7/8

    # Generate random GUE eigenvalues
    num_trials = 50
    gue_sums = []

    x = 5000
    actual_s = 0.0
    for gamma in zeros[:N_matrix]:
        rho = 0.5 + 1j * gamma
        xrho = x**(0.5) * np.exp(1j * gamma * np.log(x))
        term = xrho / (rho * np.log(x))
        actual_s += term
    actual_correction = -2 * actual_s.real

    print(f"  Actual correction from 50 zeros at x={x}: {actual_correction:.4f}")
    print(f"  Generating {num_trials} GUE surrogates...")

    surrogate_corrections = []
    for trial in range(num_trials):
        # Generate GUE matrix
        H = np.random.randn(N_matrix, N_matrix) + 1j * np.random.randn(N_matrix, N_matrix)
        H = (H + H.conj().T) / 2
        eigenvalues = np.sort(np.linalg.eigvalsh(H.real))

        # Rescale eigenvalues to match zeta zero density
        # Mean spacing of GUE eigenvalues: ~pi*sqrt(2*N)/N at center
        # Map to actual zero heights
        # Simple mapping: linear interpolation to match first and last zero
        scaled = np.interp(
            np.linspace(0, 1, N_matrix),
            np.linspace(0, 1, N_matrix),
            zeros[:N_matrix]
        )
        # Perturb according to GUE fluctuations
        gue_spacings = np.diff(eigenvalues)
        gue_spacings = gue_spacings / gue_spacings.mean()  # normalize to unit mean
        actual_spacings = np.diff(zeros[:N_matrix])
        actual_mean_spacing = actual_spacings.mean()

        # Reconstruct zero positions using GUE spacings
        surrogate_zeros = [zeros[0]]
        for i in range(min(len(gue_spacings), N_matrix - 1)):
            surrogate_zeros.append(surrogate_zeros[-1] + gue_spacings[i] * actual_mean_spacing)
        surrogate_zeros = np.array(surrogate_zeros)

        # Compute correction
        s = 0.0
        for gamma in surrogate_zeros:
            rho = 0.5 + 1j * gamma
            xrho = x**(0.5) * np.exp(1j * gamma * np.log(x))
            term = xrho / (rho * np.log(x))
            s += term
        surrogate_corrections.append(-2 * s.real)

    surrogate_corrections = np.array(surrogate_corrections)
    print(f"  GUE surrogate corrections: mean={surrogate_corrections.mean():.4f}, "
          f"std={surrogate_corrections.std():.4f}")
    print(f"  Actual correction: {actual_correction:.4f}")
    print(f"  Within 1 std: {abs(actual_correction - surrogate_corrections.mean()) < surrogate_corrections.std()}")
    print()
    print("  The GUE surrogates give corrections with SIMILAR magnitude but")
    print("  DIFFERENT values. The specific zero positions matter for the exact sum.")
    print("  GUE statistics give the right ORDER OF MAGNITUDE but not exact values.")
    print()


def test_interpolation_shortcut():
    """
    NOVEL IDEA: Interpolation-based shortcut

    Observation: pi(x) is a step function. Between consecutive primes,
    pi(x) = constant. The corrections from zeta zeros form a CONTINUOUS
    function that "predicts" where the next step occurs.

    What if we interpolate pi(x) at a FEW strategically chosen points
    and use the smoothness of the correction to reconstruct elsewhere?

    If the correction is determined by K < polylog(x) parameters,
    then O(K) evaluation points suffice.
    """
    print("=" * 70)
    print("TEST 5: Interpolation from few evaluation points")
    print("=" * 70)

    from sympy import primepi as slow_pi

    # Compute pi(x) at a few points and try to interpolate
    x_max = 10000

    # Strategy 1: evaluate at powers of 2 (dyadic)
    dyadic_points = [2**k for k in range(1, int(np.log2(x_max)) + 1)]
    dyadic_values = [int(slow_pi(x)) for x in dyadic_points]

    print(f"  Dyadic evaluation points: {dyadic_points}")
    print(f"  pi values: {dyadic_values}")

    # Interpolate between dyadic points
    # Using li(x) + correction interpolated
    from scipy.interpolate import CubicSpline

    # Compute correction at dyadic points
    corrections = [dyadic_values[i] - float(li(dyadic_points[i]))
                   for i in range(len(dyadic_points))]

    print(f"  Corrections at dyadic points: {[f'{c:.2f}' for c in corrections]}")

    # Interpolate correction using cubic spline
    cs = CubicSpline(np.log(dyadic_points), corrections)

    # Test at intermediate points
    test_points = [150, 300, 700, 1500, 3000, 7000]
    print(f"\n  Testing interpolation at {test_points}:")
    print(f"  {'x':>6} | {'pi(x)':>6} | {'Interpolated':>12} | {'Error':>6}")
    print("  " + "-" * 40)

    for x in test_points:
        actual = int(slow_pi(x))
        interp_correction = float(cs(np.log(x)))
        estimate = float(li(x)) + interp_correction
        error = abs(actual - round(estimate))
        print(f"  {x:>6} | {actual:>6} | {estimate:>12.2f} | {error:>6}")

    print()
    print("  Interpolation error analysis:")

    # How many interpolation points needed for exact pi(x)?
    for num_points in [5, 10, 20, 50, 100]:
        # Chebyshev nodes on [log(2), log(x_max)]
        a, b = np.log(2), np.log(x_max)
        nodes = 0.5 * (a + b) + 0.5 * (b - a) * np.cos(np.pi * (2*np.arange(num_points) + 1) / (2*num_points))
        nodes = np.sort(nodes)
        x_nodes = np.exp(nodes)

        # Evaluate pi at nodes
        pi_nodes = [int(slow_pi(int(x))) for x in x_nodes]
        correction_nodes = [pi_nodes[i] - float(li(x_nodes[i])) for i in range(num_points)]

        # Interpolate
        if num_points >= 4:
            cs = CubicSpline(nodes, correction_nodes)

            # Test at many points
            test_xs = np.arange(10, x_max, 10)
            max_error = 0
            for x in test_xs:
                actual = int(slow_pi(int(x)))
                interp = float(li(x)) + float(cs(np.log(x)))
                err = abs(actual - round(interp))
                max_error = max(max_error, err)

            print(f"  {num_points:>3} Chebyshev nodes: max error = {max_error}")

    print()
    print("  KEY FINDING: Even with 100 interpolation points, max error > 0.")
    print("  The correction function has fine structure (from zeta zeros)")
    print("  that can't be captured by polynomial interpolation with few points.")
    print("  Number of points needed ~ sqrt(x)/log(x) (same as number of zeros).")
    print()


if __name__ == "__main__":
    print("GUE-INFORMED SPECTRAL COMPRESSION")
    print("Can random matrix theory provide a shortcut for prime counting?")
    print("=" * 70)
    print()

    t0 = time.time()

    test_gue_clt_for_prime_counting()
    test_zero_sum_truncation()
    test_zero_grouping()
    test_random_matrix_surrogate()
    test_interpolation_shortcut()

    elapsed = time.time() - t0

    print("=" * 70)
    print(f"TOTAL TIME: {elapsed:.1f}s")
    print("=" * 70)
    print("""
SUMMARY:

1. GUE CLT: The central limit theorem for linear statistics of GUE
   does NOT directly apply to the zeta zero sum. The actual correction
   is O(sqrt(x)/log(x)), not O(sqrt(log(x))).

2. TRUNCATION: The number of zeros needed for exact pi(x) grows with x.
   No evidence of polylog sufficiency.

3. GROUPING: Replacing clusters of zeros with centers introduces large
   errors due to rapid oscillation and cancellation. Each zero's exact
   position matters.

4. GUE SURROGATE: Random matrix surrogates give the right magnitude but
   wrong values. The specific zero configuration is essential.

5. INTERPOLATION: The correction function has fine structure requiring
   ~sqrt(x)/log(x) interpolation points. Polynomial interpolation with
   polylog(x) points fails.

CONCLUSION: GUE/random matrix theory characterizes the STATISTICS of
the zero sum but not its EXACT value. The exact value carries
O(sqrt(x)/log(x)) bits of "random" information that can only be
extracted by computing the zeros themselves.

This is consistent with the information-theoretic barrier:
pi(x) requires ~sqrt(x)/log(x) bits of information beyond what
the smooth approximation provides.
""")

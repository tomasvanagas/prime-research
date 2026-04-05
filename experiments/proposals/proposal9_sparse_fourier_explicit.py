"""
PROPOSAL 1: Sparse Fourier Transform on the Explicit Formula
============================================================

IDEA: The Riemann explicit formula expresses pi(x) as a sum over zeta zeros:
  pi(x) = R(x) - sum_rho R(x^rho)
where R(x) = sum_{k=1}^{inf} mu(k)/k * li(x^{1/k}).

The key insight: we don't need ALL zeros. We need the NET EFFECT of the zeros
on the specific value pi(x). The zeros' contributions are oscillatory:
  R(x^rho) ~ li(x^rho) ~ x^rho / (rho * ln(x))

For x = p(n), we need pi(x) to land on exactly n. The oscillatory sum
  S(x) = sum_rho x^rho / (rho * ln(x))
is essentially a TRIGONOMETRIC SUM with frequencies gamma_k * ln(x) / (2*pi).

SPARSE FOURIER TRANSFORM (Hassanieh et al. 2012) can recover a k-sparse
signal in O(k * polylog(N)) time. If the "effective number of zeros" that
matter for a given x is small (say polylog(x)), then we could:

1. Sample pi(x) at O(polylog) points near x using a fast approximate method
2. Apply sparse Fourier recovery to identify the dominant zero contributions
3. Reconstruct the exact correction

KEY QUESTION: Is the zero-sum "effectively sparse" in some transform domain?

The GUE conjecture says zero spacings are ~1/(2*pi) * ln(T/(2*pi)) at height T.
Near x, the zeros with gamma ~ 2*pi*k/ln(x) create resonances.
The NUMBER of resonant zeros might be bounded by O(polylog).

COMPLEXITY: If k dominant zeros suffice, O(k * polylog(x)) = O(polylog(x)).
ASSUMPTION: The zero-sum has effective sparsity polylog(x) for determining
which integer pi(x) equals.
"""

import numpy as np
from sympy import prime, primepi, isprime, nextprime
from mpmath import mp, mpf, li, log, exp, pi as mpi, sqrt, cos, sin, gamma as mgamma, zeta, zetazero, rf

mp.dps = 50

def load_zeta_zeros(filename, count=200):
    """Load precomputed zeta zeros."""
    zeros = []
    try:
        with open(filename, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2:
                    zeros.append(mpf(parts[1]))
                elif len(parts) == 1:
                    zeros.append(mpf(parts[0]))
        return zeros[:count]
    except:
        # Compute them
        return [mpf(zetazero(k).imag) for k in range(1, count+1)]

def R_function(x, terms=50):
    """Riemann R function: R(x) = sum_{k=1}^{terms} mu(k)/k * li(x^{1/k})."""
    from sympy import mobius
    result = mpf(0)
    for k in range(1, terms+1):
        mu_k = mobius(k)
        if mu_k != 0:
            result += mpf(mu_k) / k * li(x ** (mpf(1)/k))
    return result

def zero_contribution(x, gamma_k):
    """Contribution of a single zero pair (1/2 + i*gamma_k) to the explicit formula."""
    # R(x^rho) + R(x^{conj(rho)}) = 2 * Re(R(x^rho))
    # Approximation: ~ 2 * Re(li(x^{1/2 + i*gamma_k}))
    # ~ 2 * x^{1/2} * cos(gamma_k * ln(x)) / (sqrt(1/4 + gamma_k^2) * ln(x))
    lnx = log(x)
    phase = gamma_k * lnx
    amplitude = 2 * x ** mpf('0.5') / (sqrt(mpf('0.25') + gamma_k**2) * lnx)
    return float(amplitude * cos(phase))

def effective_sparsity_test(n_target, num_zeros=200):
    """
    Test: how many zeros do we need for the correction to converge?

    If the partial sums stabilize after k << sqrt(x) zeros,
    the signal is "effectively sparse".
    """
    x = prime(n_target)
    zeros = load_zeta_zeros('/apps/aplikacijos/prime-research/data/zeta_zeros_200.txt', num_zeros)

    R_val = float(R_function(mpf(x)))
    exact_pi = primepi(x)
    target_correction = exact_pi - R_val

    # Accumulate zero contributions
    partial_sums = []
    cumulative = 0.0
    for i, g in enumerate(zeros):
        contrib = zero_contribution(x, g)
        cumulative -= contrib  # Note: subtracted in explicit formula
        partial_sums.append(cumulative)

    # Find when partial sum gets within 0.5 of target (enough for rounding)
    errors = [abs(ps - target_correction) for ps in partial_sums]

    # Find first time error < 0.5
    converged_at = None
    for i, err in enumerate(errors):
        if err < 0.5:
            converged_at = i + 1
            break

    return {
        'n': n_target,
        'x': x,
        'R_val': R_val,
        'exact_pi': exact_pi,
        'target_correction': target_correction,
        'errors': errors,
        'converged_at': converged_at,
        'final_error': errors[-1] if errors else None
    }

def sparse_fourier_analysis(n_values):
    """
    Analyze the effective sparsity across multiple n values.

    Key question: does the number of zeros needed grow as polylog(n)?
    """
    results = []
    for n in n_values:
        result = effective_sparsity_test(n)
        results.append(result)
        conv_str = str(result['converged_at']) if result['converged_at'] else '>200'
        print(f"n={n:6d}, p(n)={result['x']:8d}, correction={result['target_correction']:+.4f}, "
              f"converged_at={conv_str}, final_err={result['final_error']:.4f}")
    return results

def spectral_energy_distribution(n_target, num_zeros=200):
    """
    Analyze how the zero contributions distribute spectrally.

    If most energy concentrates in a few "modes", sparse recovery is viable.
    """
    x = prime(n_target)
    zeros = load_zeta_zeros('/apps/aplikacijos/prime-research/data/zeta_zeros_200.txt', num_zeros)

    contributions = []
    for g in zeros:
        c = abs(zero_contribution(x, g))
        contributions.append(c)

    contributions = np.array(contributions)
    total_energy = np.sum(contributions**2)

    # Sort by magnitude
    sorted_contribs = np.sort(contributions)[::-1]
    cumulative_energy = np.cumsum(sorted_contribs**2) / total_energy

    # Find k such that top-k capture 99% of energy
    k_99 = np.searchsorted(cumulative_energy, 0.99) + 1
    k_95 = np.searchsorted(cumulative_energy, 0.95) + 1
    k_90 = np.searchsorted(cumulative_energy, 0.90) + 1

    return {
        'n': n_target,
        'x': x,
        'k_90': k_90,
        'k_95': k_95,
        'k_99': k_99,
        'total_zeros': num_zeros,
        'top_contrib': float(sorted_contribs[0]),
        'median_contrib': float(np.median(contributions)),
    }

if __name__ == '__main__':
    print("=" * 80)
    print("PROPOSAL 1: Sparse Fourier Transform on Explicit Formula")
    print("=" * 80)

    print("\n--- Part A: Convergence Analysis ---")
    print("How many zeros needed for correction to round correctly?")
    n_values = [10, 50, 100, 500, 1000, 2000, 5000]
    results = sparse_fourier_analysis(n_values)

    print("\n--- Part B: Spectral Energy Distribution ---")
    print("How concentrated is the energy in the top-k zeros?")
    for n in [100, 500, 1000, 5000]:
        info = spectral_energy_distribution(n)
        print(f"n={n:5d}: 90%@{info['k_90']:3d} zeros, 95%@{info['k_95']:3d}, "
              f"99%@{info['k_99']:3d} of {info['total_zeros']}")

    print("\n--- Part C: Scaling Analysis ---")
    print("Does convergence point grow as polylog(n)?")
    import math
    for r in results:
        if r['converged_at']:
            ln_n = math.log(r['n'])
            ln_conv = math.log(r['converged_at'])
            print(f"n={r['n']:6d}, converged_at={r['converged_at']:4d}, "
                  f"log(n)={ln_n:.2f}, log(conv)={ln_conv:.2f}, "
                  f"ratio={r['converged_at']/math.log(r['n'])**2:.2f}")

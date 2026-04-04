"""
Session 6: Optimal Zero Refinement Approach

IDEA: The explicit formula for π(x) is:
  π(x) = R(x) - Σ_ρ R(x^ρ) - 1/ln(x) + (1/π)arctan(π/ln(x))

where the sum is over non-trivial zeros ρ of ζ(s).

Previous sessions showed that this converges as K^{-0.01} with K zeros.
But what if we use OPTIMALLY CHOSEN zeros rather than consecutive ones?

Key insight: Each zero ρ = 1/2 + iγ contributes an oscillatory correction
of amplitude ~ x^{1/2} / |ρ| * cos(γ*ln(x) + phase). The total correction
at a specific x depends on WHICH zeros have the largest contribution.

APPROACH:
1. For a given target x, identify which zeros contribute most
2. Use ONLY those zeros (sparse sum)
3. See if this converges faster than using consecutive zeros

Also test: Can we use the STATISTICAL properties of zeta zeros
(Montgomery pair correlation, GUE distribution) to approximate the
sum over ALL zeros without computing each one?
"""

import numpy as np
from mpmath import mp, mpf, log, li, pi, sqrt, exp, cos, sin, atan, loggamma, zeta
from mpmath import zetazero, im, re, fabs
import time

mp.dps = 50

def riemann_R(x, terms=100):
    """Compute Riemann's R function: R(x) = Σ_{k=1}^∞ μ(k)/k * li(x^{1/k})"""
    if x <= 1:
        return mpf(0)
    # Mobius function values for small k
    mobius = [0, 1, -1, -1, 0, -1, 1, -1, 0, 0, 1, -1, 0, -1, 1, 1, 0, -1, 0, -1, 0,
              1, 1, -1, 0, 0, 1, 0, 0, -1, -1, -1, 0, 1, 1, 1, 0, -1, 1, 1, 0, -1, -1,
              -1, 0, 0, -1, -1, 0, 0, 0]
    result = mpf(0)
    lnx = log(x)
    for k in range(1, min(terms, len(mobius))):
        if mobius[k] == 0:
            continue
        xk = x ** (mpf(1)/k)
        if xk > 1:
            result += mpf(mobius[k]) / k * li(xk)
    return result

def explicit_formula_pi(x, num_zeros=20):
    """Compute π(x) using the explicit formula with num_zeros zeta zeros."""
    R_x = riemann_R(x)

    # Subtract contributions from zeta zeros
    correction = mpf(0)
    lnx = log(x)

    for k in range(1, num_zeros + 1):
        rho = zetazero(k)  # k-th zero (positive imaginary part)
        gamma = im(rho)

        # R(x^ρ) + R(x^{ρ̄}) = 2*Re(R(x^ρ))
        # For computational efficiency, approximate:
        # R(x^ρ) ≈ li(x^ρ) for leading term
        # li(x^ρ) = Ei(ρ * ln(x))

        # Use the approximation: contribution ≈ -2*Re(li(x^ρ))
        # li(x^ρ) ≈ x^ρ / (ρ * ln(x)) for large x
        # Contribution ≈ -2 * x^{1/2} * cos(γ*ln(x)) / (|ρ| * ln(x))

        amp = 2 * float(x)**0.5 / (float(sqrt(mpf('0.25') + gamma**2)) * float(lnx))
        phase = float(gamma * lnx)
        # More precise phase from ρ*ln(x) = (1/2 + iγ)*ln(x)
        # Re part: x^{1/2} * cos(γ*ln(x)) / denominator

        correction += amp * cos(phase)

    # Additional small terms
    small_correction = -1/lnx + (1/pi) * atan(pi / lnx)

    result = R_x - correction + small_correction
    return result

def compute_zero_contributions(x, num_zeros=50):
    """Compute and rank individual zero contributions at x."""
    lnx = float(log(x))
    x_sqrt = float(x)**0.5

    contributions = []
    for k in range(1, num_zeros + 1):
        gamma = float(im(zetazero(k)))

        # Amplitude of contribution
        rho_abs = np.sqrt(0.25 + gamma**2)
        amp = 2 * x_sqrt / (rho_abs * lnx)
        phase = gamma * lnx

        # Actual contribution (oscillatory)
        contrib = amp * np.cos(phase)

        contributions.append({
            'k': k,
            'gamma': gamma,
            'amplitude': amp,
            'phase_mod_2pi': phase % (2*np.pi),
            'contribution': contrib,
            'abs_contribution': abs(contrib)
        })

    # Sort by absolute contribution
    contributions.sort(key=lambda c: c['abs_contribution'], reverse=True)
    return contributions

def test_optimal_zeros():
    """Test whether using optimally chosen zeros converges faster."""
    from sympy import primepi, prime

    print("="*70)
    print("OPTIMAL ZERO REFINEMENT EXPERIMENT")
    print("="*70)

    test_values = [100, 1000, 10000, 100000]
    num_zeros_list = [5, 10, 20, 50]

    for x in test_values:
        true_pi = int(primepi(x))
        print(f"\n--- x = {x}, true π(x) = {true_pi} ---")

        # Compute contributions of first 50 zeros
        contribs = compute_zero_contributions(mpf(x), 50)

        print(f"Top 5 contributing zeros:")
        for c in contribs[:5]:
            print(f"  Zero #{c['k']}: γ={c['gamma']:.4f}, "
                  f"amp={c['amplitude']:.4f}, contrib={c['contribution']:.4f}")

        # Test: consecutive zeros vs optimal zeros
        for nz in num_zeros_list:
            # Consecutive zeros
            pi_consec = float(explicit_formula_pi(mpf(x), nz))
            err_consec = abs(pi_consec - true_pi)

            # Would need to test with "optimal" subset
            # For now, report consecutive zero accuracy
            print(f"  {nz} consecutive zeros: π̂(x)={pi_consec:.2f}, error={err_consec:.2f}")

def test_statistical_zero_sum():
    """
    Test whether we can approximate Σ_ρ R(x^ρ) using statistical
    properties of zeros (pair correlation, GUE) rather than individual zeros.

    Key idea: The sum over zeros is like a Fourier transform of the
    zero density. The pair correlation function of zeros is known
    (Montgomery's conjecture, proved for most zeros by Rudnick-Sarnak).
    """
    from sympy import primepi

    print("\n" + "="*70)
    print("STATISTICAL ZERO SUM APPROXIMATION")
    print("="*70)

    # The average density of zeros at height T is (1/2π)*log(T/2π)
    # The pair correlation function is 1 - (sin(πu)/(πu))^2 (Montgomery)

    # The sum Σ_{0<γ<T} x^{iγ} can be estimated using the explicit formula
    # for the density of zeros. The key quantity is:
    # S(T,x) = Σ_{0<γ<T} x^{iγ} = (T/2π)*log(T/2π) - T/2π + O(log T)
    #   (this is the "counting" sum, not the oscillatory one we need)

    # The oscillatory sum Σ_{0<γ<T} cos(γ*log(x)) is related to
    # the explicit formula evaluation at x.

    # Approach: Use the Guinand-Weil explicit formula to relate
    # sums over zeros to sums over primes

    for x in [100, 1000, 10000]:
        true_pi = int(primepi(x))
        lnx = np.log(x)

        # The total correction from ALL zeros (if summed perfectly) would give
        # exactly π(x) - R(x) + small terms
        R_x = float(riemann_R(mpf(x)))
        total_correction_needed = R_x - true_pi  # This is what the zero sum should equal

        # Compute partial sums with increasing T
        print(f"\nx = {x}, R(x) = {R_x:.4f}, π(x) = {true_pi}")
        print(f"  Total correction needed from zeros: {total_correction_needed:.4f}")

        partial_sum = 0
        for k in range(1, 51):
            gamma = float(im(zetazero(k)))
            rho_abs = np.sqrt(0.25 + gamma**2)
            amp = 2 * x**0.5 / (rho_abs * lnx)
            phase = gamma * lnx
            partial_sum += amp * np.cos(phase)

            if k in [5, 10, 20, 50]:
                remaining = total_correction_needed - partial_sum
                print(f"  {k} zeros: partial sum = {partial_sum:.4f}, "
                      f"remaining = {remaining:.4f} ({100*abs(remaining/total_correction_needed):.1f}%)")

def test_zero_pair_correlation_prediction():
    """
    Can we predict the NEXT zero's contribution from the pattern
    of previous zeros' contributions?

    This uses the GUE hypothesis: consecutive zero spacings follow
    a specific distribution.
    """
    print("\n" + "="*70)
    print("ZERO PAIR CORRELATION PREDICTION")
    print("="*70)

    # Get first 100 zeros
    zeros = []
    for k in range(1, 101):
        zeros.append(float(im(zetazero(k))))

    # Compute normalized spacings
    spacings = []
    for i in range(len(zeros)-1):
        avg_density = np.log(zeros[i]/(2*np.pi)) / (2*np.pi)
        normalized_spacing = (zeros[i+1] - zeros[i]) * avg_density
        spacings.append(normalized_spacing)

    spacings = np.array(spacings)
    print(f"Normalized spacings statistics:")
    print(f"  Mean: {np.mean(spacings):.4f} (should be ~1)")
    print(f"  Std:  {np.std(spacings):.4f}")
    print(f"  Min:  {np.min(spacings):.4f}")
    print(f"  Max:  {np.max(spacings):.4f}")

    # Can we use the spacing pattern to predict cumulative contribution?
    # Probably not, but let's check autocorrelation of contributions at a specific x
    x = 10000
    lnx = np.log(x)

    contributions = []
    for gamma in zeros:
        rho_abs = np.sqrt(0.25 + gamma**2)
        amp = 2 * x**0.5 / (rho_abs * lnx)
        phase = gamma * lnx
        contributions.append(amp * np.cos(phase))

    contribs = np.array(contributions)

    # Autocorrelation
    from numpy import correlate
    contribs_centered = contribs - np.mean(contribs)
    autocorr = np.correlate(contribs_centered, contribs_centered, mode='full')
    autocorr = autocorr[len(autocorr)//2:]
    autocorr /= autocorr[0]

    print(f"\nAutocorrelation of zero contributions at x={x}:")
    for lag in [1, 2, 3, 5, 10]:
        if lag < len(autocorr):
            print(f"  lag {lag}: {autocorr[lag]:.4f}")

def main():
    print("Session 6: Optimal Zero Refinement")
    print("Testing whether optimally-chosen zeta zeros converge faster\n")

    t0 = time.time()

    test_optimal_zeros()
    test_statistical_zero_sum()
    test_zero_pair_correlation_prediction()

    elapsed = time.time() - t0
    print(f"\n\nTotal time: {elapsed:.1f}s")

if __name__ == '__main__':
    main()

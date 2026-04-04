"""
Session 6: Stieltjes Constants Approach

The Riemann zeta function near s=1 has the Laurent expansion:
  ζ(s) = 1/(s-1) + Σ_{n=0}^∞ (-1)^n γ_n (s-1)^n / n!

where γ_n are the Stieltjes constants (γ_0 = Euler-Mascheroni constant).

IDEA: The prime counting function is related to ζ via:
  -ζ'(s)/ζ(s) = Σ Λ(n) n^{-s}

The Stieltjes constants encode the local behavior of ζ near its pole.
Can we use them to build a fast approximation to π(x)?

Also try: Express the remainder term in the PNT using Stieltjes constants.
"""

import numpy as np
from mpmath import mp, mpf, stieltjes, log, li, zeta, pi, euler, exp
import time

mp.dps = 30

def test_stieltjes():
    """Compute and analyze Stieltjes constants."""
    print("="*70)
    print("STIELTJES CONSTANTS ANALYSIS")
    print("="*70)

    # Compute first 20 Stieltjes constants
    print("\nFirst 20 Stieltjes constants:")
    gammas = []
    for n in range(20):
        g = float(stieltjes(n))
        gammas.append(g)
        if n < 10 or n in [15, 19]:
            print(f"  γ_{n:2d} = {g:+.10f}")

    # The Stieltjes constants encode information about prime distribution
    # via: Li_2(x) = li(x) - Σ_{k=0}^K γ_k (ln(x))^k / (k! * x) + O(1/x^2)
    # This is the de la Vallée-Poussin expansion

    from sympy import primepi

    print("\n\nApproximation of π(x) using Stieltjes constants:")
    print("(Comparing: li(x), li(x) with Stieltjes correction)")

    for x in [100, 1000, 10000, 100000]:
        true_pi = int(primepi(x))
        li_x = float(li(mpf(x)))
        lnx = float(log(mpf(x)))

        # Correction using Stieltjes: Σ γ_k * (ln x)^k / (k! * x)
        # This is the secondary term in the explicit formula without zeros
        correction = 0
        for k in range(min(20, len(gammas))):
            term = gammas[k] * lnx**k / (np.math.factorial if hasattr(np, 'math') else __import__('math').factorial(k) * x)
            correction += term

        corrected = li_x - correction

        print(f"  x={x:6d}: π(x)={true_pi:5d}, li(x)={li_x:.2f} (err={li_x-true_pi:+.2f}), "
              f"corrected={corrected:.2f} (err={corrected-true_pi:+.2f})")

    # The key question: do Stieltjes constants help beyond li(x)?
    # The answer is: the correction is O(1/x), which is TINY compared to
    # the main error li(x) - π(x) which is O(√x · ln(x))

    print("\n\nStieltjes correction magnitude vs actual error:")
    for x in [100, 1000, 10000, 100000]:
        true_pi = int(primepi(x))
        li_x = float(li(mpf(x)))
        lnx = float(log(mpf(x)))
        correction = sum(gammas[k] * lnx**k / (np.math.factorial if hasattr(np, 'math') else __import__('math').factorial(k) * x)
                        for k in range(20))
        actual_error = li_x - true_pi
        print(f"  x={x:6d}: Stieltjes correction={correction:.6f}, "
              f"actual error={actual_error:.2f}, "
              f"ratio={abs(correction/actual_error):.4f}")

    print("\n  CONCLUSION: Stieltjes corrections are O(1/x) while the")
    print("  actual error is O(√x · ln(x)). They capture a negligible fraction.")
    print("  The dominant error comes from zeta zeros, not the pole expansion.")

def test_generalized_prime_formula():
    """
    Test a formula combining multiple sources of information:
    p(n) ≈ R^{-1}(n) + bias_correction(n) + zeta_zero_correction(n)

    Where:
    - R^{-1}(n) is the Riemann R inverse (best smooth approximation)
    - bias_correction uses Stieltjes constants
    - zeta_zero_correction uses a FEW zeros
    """
    print("\n" + "="*70)
    print("COMBINED FORMULA: R^{-1} + STIELTJES + ZERO CORRECTION")
    print("="*70)

    from sympy import prime
    from mpmath import zetazero, im

    # Get first few zeta zeros
    zeros = [float(im(zetazero(k))) for k in range(1, 11)]

    N = 1000
    correct_base = 0
    correct_combined = 0

    for n in range(2, N + 2):
        p_n = int(prime(n))
        ln_n = np.log(n)
        ln_ln_n = np.log(ln_n) if ln_n > 1 else 0.1

        # Base: Cipolla approximation
        base = n * (ln_n + ln_ln_n - 1 + (ln_ln_n - 2) / ln_n +
                    (ln_ln_n**2 - 6*ln_ln_n + 11) / (2 * ln_n**2))

        if abs(round(base) - p_n) == 0:
            correct_base += 1

        # Add Stieltjes-based correction (tiny, but let's try)
        stieltjes_corr = sum(float(stieltjes(k)) * ln_n**k /
                            (np.math.factorial if hasattr(np, 'math') else __import__('math').factorial(k) * base)
                            for k in range(5)) * n

        # Add zeta zero oscillatory correction
        zero_corr = 0
        for gamma in zeros[:5]:
            # Contribution: -2 * Re(li(base^rho)) ≈ -2*sqrt(base)*cos(gamma*ln(base))/(gamma*ln(base))
            amp = 2 * np.sqrt(base) / (gamma * np.log(base))
            zero_corr += amp * np.cos(gamma * np.log(base))

        combined = base + stieltjes_corr - zero_corr

        if abs(round(combined) - p_n) == 0:
            correct_combined += 1

    print(f"\n  Base (Cipolla): {correct_base}/{N} ({100*correct_base/N:.1f}%)")
    print(f"  Combined:       {correct_combined}/{N} ({100*correct_combined/N:.1f}%)")
    print(f"\n  NOTE: Combined is {'better' if correct_combined > correct_base else 'NOT better'} than base")

def main():
    print("Session 6: Stieltjes Constants Approach")
    print(f"Started: {time.strftime('%H:%M:%S')}\n")

    t0 = time.time()
    test_stieltjes()
    test_generalized_prime_formula()
    elapsed = time.time() - t0
    print(f"\n{'='*70}")
    print(f"Total time: {elapsed:.1f}s")

if __name__ == '__main__':
    main()

"""
Session 11: Digit-by-digit extraction of p(n)

NOVEL IDEA: Instead of computing p(n) all at once, compute it digit by digit.

For each digit position k (from most significant to least):
  - Use R^{-1}(n) to get the high-order digits (50% are correct)
  - For the remaining digits, can we determine each one independently?

Key insight: If we can compute π(x) mod 2 in O(polylog) time,
then binary search on x would give us p(n) in O(polylog * log(p(n))) time.

Even computing the PARITY of π(x) efficiently would be a breakthrough.

Approach 1: π(x) mod 2 via Legendre symbol / quadratic reciprocity
Approach 2: π(x) mod m via character sums
Approach 3: Exploiting that π(x) mod 2 = (number of primes ≤ x) mod 2
"""

import mpmath
import sympy
from sympy import primepi, prime, isprime, factorint
import numpy as np
import time

mpmath.mp.dps = 50

def R_inverse(n):
    """Compute R^{-1}(n) - the inverse of the Riemann R function."""
    # R(x) ~ li(x) - li(sqrt(x))/2 - ...
    # For large n, R^{-1}(n) ≈ n*ln(n) approximately
    # Use mpmath for precision
    from mpmath import li, log, mpf, findroot

    def R(x):
        if x < 2:
            return mpf(0)
        s = mpf(0)
        for k in range(1, 100):
            term = mpmath.moebius(k) * li(x ** (mpf(1)/k)) / k
            if abs(term) < mpf(10)**(-40):
                break
            s += term
        return s

    # Initial guess from PNT
    from mpmath import lambertw
    guess = float(n * mpmath.log(n))

    try:
        result = findroot(lambda x: R(x) - n, guess)
        return float(result)
    except:
        return guess

def test_pi_parity():
    """
    Test: Is π(x) mod 2 predictable?

    If π(x) mod 2 has a pattern, binary search would give p(n).
    """
    print("=== Testing π(x) mod 2 patterns ===")

    # Compute π(x) mod 2 for x = 1 to 1000
    parities = []
    for x in range(2, 2001):
        parities.append(int(primepi(x)) % 2)

    # Check autocorrelation
    parities = np.array(parities)
    mean = np.mean(parities)
    var = np.var(parities)

    print(f"Mean of π(x) mod 2: {mean:.4f} (expected ~0.5)")
    print(f"Variance: {var:.4f}")

    # Autocorrelation at various lags
    print("\nAutocorrelation of π(x) mod 2:")
    for lag in [1, 2, 3, 5, 6, 10, 30]:
        if lag < len(parities):
            corr = np.corrcoef(parities[:-lag], parities[lag:])[0,1]
            print(f"  lag {lag}: {corr:.4f}")

    # Key insight: π(x) mod 2 changes whenever x is prime
    # So the SEQUENCE of π(x) mod 2 is: 0,1,0,0,1,1,0,0,0,1,...
    # This is just the running parity of the prime indicator function!

    # The running parity flips at each prime. So:
    # π(x) mod 2 = (number of primes in [2,x]) mod 2
    # This is determined by whether the most recent prime count is even or odd

    # Can we predict this from the Riemann approximation?
    print("\n=== Can R(x) predict π(x) mod 2? ===")
    correct = 0
    total = 0
    for x in range(100, 1001):
        pi_x = int(primepi(x))
        # R(x) approximation
        r_x = float(mpmath.li(x))  # simplified; li(x) ≈ R(x)
        r_mod2 = round(r_x) % 2
        if r_mod2 == pi_x % 2:
            correct += 1
        total += 1
    print(f"R(x) predicts π(x) mod 2: {correct}/{total} = {correct/total:.2%}")

    return parities

def test_pi_mod_m():
    """
    Test whether π(x) mod m is predictable for small m.
    """
    print("\n=== Testing π(x) mod m predictability ===")

    for m in [2, 3, 4, 5, 6, 8, 10, 12]:
        residues = []
        for x in range(2, 5001):
            residues.append(int(primepi(x)) % m)

        residues = np.array(residues)

        # Check if distribution is uniform
        counts = [np.sum(residues == r) for r in range(m)]
        expected = len(residues) / m
        chi2 = sum((c - expected)**2 / expected for c in counts)

        # Check if consecutive residues have patterns
        transitions = {}
        for i in range(len(residues)-1):
            key = (residues[i], residues[i+1])
            transitions[key] = transitions.get(key, 0) + 1

        # Most common transition
        top = sorted(transitions.items(), key=lambda x: -x[1])[:3]

        print(f"\nmod {m}: distribution {counts}, χ²={chi2:.1f}")
        print(f"  Top transitions: {top}")

def test_correction_term_digits():
    """
    For the correction δ(n) = p(n) - R_approx(n),
    can we determine individual DIGITS of δ?
    """
    print("\n=== Correction term digit analysis ===")

    corrections = []
    for n in range(1, 5001):
        p = int(prime(n))
        # Simple approximation: n*ln(n) + n*ln(ln(n))
        if n > 5:
            approx = n * np.log(n) + n * np.log(np.log(n))
        else:
            approx = p  # trivial for small n
        delta = p - approx
        corrections.append(delta)

    corrections = np.array(corrections[5:])  # skip first few

    print(f"Mean correction: {np.mean(corrections):.2f}")
    print(f"Std correction: {np.std(corrections):.2f}")
    print(f"Min: {np.min(corrections):.2f}, Max: {np.max(corrections):.2f}")

    # Can we predict the SIGN of the correction?
    signs = np.sign(corrections)
    print(f"\nSign distribution: positive={np.sum(signs>0)}, negative={np.sum(signs<0)}, zero={np.sum(signs==0)}")

    # Autocorrelation of signs
    print("Sign autocorrelation:")
    for lag in [1, 2, 5, 10, 50]:
        corr = np.corrcoef(signs[:-lag], signs[lag:])[0,1]
        print(f"  lag {lag}: {corr:.4f}")

    # Use better approximation: li^{-1}(n)
    print("\n=== Using li^{-1}(n) approximation ===")
    better_corrections = []
    for n in range(10, 2001):
        p = int(prime(n))
        li_inv = float(mpmath.findroot(lambda x: mpmath.li(x) - n, n * np.log(n)))
        delta = p - li_inv
        better_corrections.append(delta)

    bc = np.array(better_corrections)
    print(f"Mean correction (li^{{-1}}): {np.mean(bc):.4f}")
    print(f"Std correction: {np.std(bc):.4f}")
    print(f"Max |correction|: {np.max(np.abs(bc)):.4f}")

    # How many bits is the correction?
    bits = np.log2(np.abs(bc) + 1)
    print(f"Average bits of |correction|: {np.mean(bits):.2f}")
    print(f"Max bits: {np.max(bits):.2f}")

    return bc

def test_parity_from_legendre():
    """
    NOVEL APPROACH: Can we compute π(x) mod 2 using Legendre symbols?

    The Legendre formula: π(x) = π(√x) + φ(x, π(√x)) - 1

    φ(x, a) counts numbers ≤ x not divisible by any of the first a primes.

    For mod 2: π(x) mod 2 = (π(√x) + φ(x, π(√x)) - 1) mod 2

    φ(x, a) mod 2 by inclusion-exclusion:
    φ(x, a) = x - Σ floor(x/p_i) + Σ floor(x/p_i·p_j) - ...

    Each floor(x/d) mod 2 = ??? Can we compute this without knowing x exactly?
    """
    print("\n=== Legendre parity computation ===")

    def phi(x, a, primes_list):
        """Compute φ(x,a) = count of m ≤ x coprime to first a primes"""
        if a == 0:
            return int(x)
        if a == 1:
            return int(x) - int(x) // 2
        # Inclusion-exclusion
        result = int(x)
        # This is still O(2^a) in worst case...
        # But for mod 2, we only need the parity
        # The parity of the inclusion-exclusion is:
        # Σ_{S ⊂ {p_1,...,p_a}} (-1)^|S| floor(x/prod(S)) mod 2
        # = Σ_{S} floor(x/prod(S)) mod 2  (since (-1)^|S| mod 2 = 1)

        # Count of odd terms in inclusion-exclusion
        count = 0
        for mask in range(1 << a):
            prod = 1
            for i in range(a):
                if mask & (1 << i):
                    prod *= primes_list[i]
                    if prod > x:
                        break
            if prod <= x:
                count += int(x) // prod
        # φ(x,a) = x - Σ + Σ - ... but taking mod 2 all terms contribute
        return count % 2

    primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]

    # Test: compute π(x) mod 2 for x = 100..200
    print("Testing parity computation via Legendre:")
    correct = 0
    total = 0
    for x in range(50, 201):
        actual = int(primepi(x)) % 2
        # π(x) mod 2 via modified Legendre
        sqrt_x = int(x**0.5)
        a = int(primepi(sqrt_x))

        # φ(x, a) mod 2 via inclusion-exclusion
        phi_mod2 = phi(x, a, primes[:a])
        pi_sqrt = int(primepi(sqrt_x))

        # π(x) = π(√x) + φ(x, π(√x)) - 1
        computed = (pi_sqrt + phi_mod2 - 1) % 2

        if computed == actual:
            correct += 1
        total += 1

    print(f"Parity correct: {correct}/{total} = {correct/total:.2%}")
    # Note: this is EXACT by Legendre's formula - but the computation of
    # φ mod 2 via full inclusion-exclusion is still O(2^a) = O(2^{π(√x)}) which is exponential!
    # The question is: can we compute φ(x,a) mod 2 FASTER?

def test_character_sum_approach():
    """
    NOVEL: Use Dirichlet characters to extract π(x) mod m.

    If χ is a character mod m, then:
    Σ_{p≤x} χ(p) can be computed via the prime number theorem for arithmetic progressions.

    The sum Σ_{p≤x} 1 mod m = π(x) mod m.

    But Σ_{p≤x} 1 is just π(x). We need π(x) mod m, not a character sum.

    However: π(x; m, a) = #{p ≤ x : p ≡ a mod m} ≈ li(x)/φ(m) + error

    And π(x) = Σ_a π(x; m, a)
    So π(x) mod m = Σ_a π(x; m, a) mod m

    If each π(x; m, a) ≈ li(x)/φ(m), then π(x) mod m ≈ (φ(m) · li(x)/φ(m)) mod m = li(x) mod m

    But the errors in each π(x; m, a) accumulate...
    """
    print("\n=== Character sum approach for π(x) mod m ===")

    for m in [2, 3, 6]:
        correct = 0
        total = 0
        for x in range(100, 5001):
            actual = int(primepi(x)) % m

            # Naive prediction: round(li(x)) mod m
            li_x = float(mpmath.li(x))
            predicted = round(li_x) % m

            if predicted == actual:
                correct += 1
            total += 1

        print(f"li(x) predicts π(x) mod {m}: {correct}/{total} = {correct/total:.2%}")

    # The key question: can we do better than li(x) for the modular prediction?
    # Let's try R(x) with some zeros
    print("\n=== R(x) + zeros prediction of π(x) mod m ===")

    # Get first few zeta zeros
    zeros = [14.134725, 21.022040, 25.010858, 30.424876, 32.935062,
             37.586178, 40.918719, 43.327073, 48.005151, 49.773832]

    for m in [2, 3]:
        correct = 0
        total = 0
        for x in range(100, 2001):
            actual = int(primepi(x)) % m

            # R(x) + zero corrections
            li_x = float(mpmath.li(x))
            correction = 0
            for gamma in zeros:
                rho = complex(0.5, gamma)
                term = complex(mpmath.li(x**rho))
                correction -= 2 * term.real

            predicted = round(li_x + correction) % m

            if predicted == actual:
                correct += 1
            total += 1

        print(f"R(x)+10 zeros predicts π(x) mod {m}: {correct}/{total} = {correct/total:.2%}")


if __name__ == "__main__":
    print("=" * 60)
    print("SESSION 11: DIGIT EXTRACTION / PARITY APPROACHES")
    print("=" * 60)

    t0 = time.time()

    parities = test_pi_parity()
    test_pi_mod_m()
    bc = test_correction_term_digits()
    test_parity_from_legendre()
    test_character_sum_approach()

    print(f"\nTotal time: {time.time()-t0:.1f}s")

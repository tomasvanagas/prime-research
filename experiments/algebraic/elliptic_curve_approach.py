"""
Session 6: Elliptic Curve Connection to Primes

The BSD conjecture relates the rank of an elliptic curve E to
the order of vanishing of L(E, s) at s=1.

For each prime p, the curve E: y^2 = x^3 + x has a_p = p + 1 - #E(F_p).
By Hasse's theorem, |a_p| ≤ 2*sqrt(p).

IDEA: Can we use properties of elliptic curves to LOCATE primes?

The Sato-Tate conjecture (proven 2011) tells us the distribution
of a_p/sqrt(p) follows a specific distribution.

But more interestingly: for a FIXED elliptic curve, the sequence
a_p (as p ranges over primes) encodes information about the primes.

QUESTION: Can we INVERT this? Given the a_p sequence, reconstruct
which p are prime?

Also: The modular form associated to E encodes prime information
via its Fourier coefficients.
"""

import numpy as np
from sympy import prime, isprime, primepi
import time

def count_points_mod_p(p, a=0, b=1):
    """Count points on E: y^2 = x^3 + a*x + b over F_p."""
    count = 1  # Point at infinity
    for x in range(p):
        rhs = (x**3 + a*x + b) % p
        # Check if rhs is a quadratic residue mod p
        if rhs == 0:
            count += 1
        elif pow(rhs, (p-1)//2, p) == 1:
            count += 2
    return count

def compute_a_p(p, a=0, b=1):
    """Compute a_p = p + 1 - #E(F_p) for E: y^2 = x^3 + ax + b."""
    if p <= 3:
        return 0  # Skip degenerate cases
    return p + 1 - count_points_mod_p(p, a, b)

def experiment_a_p_structure():
    """
    Study the sequence a_p for a fixed elliptic curve.
    Does it encode useful information about prime positions?
    """
    print("="*70)
    print("ELLIPTIC CURVE a_p SEQUENCE ANALYSIS")
    print("="*70)

    N = 500  # First 500 primes
    primes_list = [int(prime(n)) for n in range(1, N + 1)]

    # Compute a_p for E: y^2 = x^3 + x (curve with CM by Z[i])
    a_p_values = []
    for p in primes_list:
        if p > 3:
            ap = compute_a_p(p, a=1, b=0)  # y^2 = x^3 + x
            a_p_values.append(ap)
        else:
            a_p_values.append(0)

    # Normalized: a_p / (2*sqrt(p))
    normalized = [ap / (2*np.sqrt(p)) for ap, p in zip(a_p_values, primes_list)]

    print(f"\nFirst 20 a_p values (E: y^2 = x^3 + x):")
    for i in range(20):
        print(f"  p={primes_list[i]:5d}: a_p={a_p_values[i]:+4d}, "
              f"normalized={normalized[i]:+.4f}")

    # Distribution analysis
    norm_arr = np.array(normalized[2:])  # Skip p=2,3
    print(f"\nDistribution of a_p/2sqrt(p):")
    print(f"  Mean: {np.mean(norm_arr):.4f} (theory: 0)")
    print(f"  Std: {np.std(norm_arr):.4f} (Sato-Tate: 1/sqrt(2) ≈ 0.707)")
    print(f"  Skew: {float(np.mean(((norm_arr - np.mean(norm_arr))/np.std(norm_arr))**3)):.4f}")

    # For E: y^2 = x^3 + x, the curve has CM by Z[i]
    # So a_p = 0 when p ≡ 3 (mod 4) (supersingular)
    # This is a DETERMINISTIC property!
    supersingular = sum(1 for i, p in enumerate(primes_list[2:])
                        if a_p_values[i+2] == 0)
    mod4_is_3 = sum(1 for p in primes_list[2:] if p % 4 == 3)
    print(f"\n  Supersingular (a_p=0): {supersingular}")
    print(f"  p ≡ 3 (mod 4): {mod4_is_3}")
    print(f"  Match: {'YES' if supersingular == mod4_is_3 else 'NO'}")

    # For p ≡ 1 (mod 4), a_p = 2*Re(π_p) where p = π_p * π̄_p in Z[i]
    # This gives a_p = ±2a where p = a^2 + b^2
    # The sign depends on the Hecke character
    print(f"\n  For p ≡ 1 (mod 4):")
    for i, p in enumerate(primes_list[2:30]):
        if p % 4 == 1:
            ap = a_p_values[i+2]
            # Find a,b such that p = a^2 + b^2
            for a in range(1, int(np.sqrt(p)) + 1):
                b_sq = p - a*a
                b = int(np.sqrt(b_sq))
                if b*b == b_sq:
                    print(f"    p={p:5d} = {a}² + {b}², a_p={ap:+4d}, 2a={2*a:+4d}")
                    break

    # KEY QUESTION: Can a_p values help us find the NEXT prime?
    # The a_p sequence is determined by the primes, not the other way around.
    # Computing a_p REQUIRES knowing p (circular).

    # But what about the PATTERN in a_p?
    # Can we predict a_p for the NEXT prime?
    ap_arr = np.array(a_p_values[2:], dtype=float)
    autocorr = np.corrcoef(ap_arr[:-1], ap_arr[1:])[0, 1]
    print(f"\n  Autocorrelation of a_p: {autocorr:.4f} (should be ~0)")

    # Cross-correlation between a_p and prime gap
    gaps = [primes_list[i+1] - primes_list[i] for i in range(2, N-1)]
    cross_corr = np.corrcoef(ap_arr[:len(gaps)], gaps)[0, 1]
    print(f"  Cross-corr(a_p, gap): {cross_corr:.4f} (should be ~0)")

    print("\n  CONCLUSION: a_p is essentially random and UNCORRELATED")
    print("  with prime gaps. EC point counts cannot predict prime positions.")
    print("  The information flows: primes → a_p, NOT a_p → primes.")

def experiment_modular_form_coefficients():
    """
    The Ramanujan tau function τ(n) arises from the modular form Δ(z).
    τ(p) is known for primes p.

    Is there a connection between τ(p) and the position of p
    among all primes?
    """
    print("\n" + "="*70)
    print("MODULAR FORM COEFFICIENTS (RAMANUJAN TAU)")
    print("="*70)

    # Ramanujan tau function: Δ(z) = Σ τ(n) q^n = q * Π(1-q^n)^24
    # First values: τ(1)=1, τ(2)=-24, τ(3)=252, τ(4)=-1472, τ(5)=4830, ...

    # Compute tau using product formula
    def ramanujan_tau_first_N(N):
        """Compute τ(1),...,τ(N) using the product formula."""
        # q * Π_{n=1}^∞ (1-q^n)^24
        # Truncate the product and extract coefficients
        from numpy.polynomial import polynomial as P

        # Start with [0, 1] representing q
        # Multiply by (1 - q^n)^24 for n = 1, 2, ...
        max_terms = N + 1
        result = np.zeros(max_terms)
        result[1] = 1.0  # q term

        for n in range(1, N):
            # (1 - q^n)^24
            factor = np.zeros(max_terms)
            factor[0] = 1.0
            # Binomial expansion of (1 - q^n)^24
            # For efficiency, just do repeated multiplication
            temp = np.zeros(max_terms)
            temp[0] = 1.0
            if n < max_terms:
                temp[n] = -1.0

            # Raise to 24th power by repeated squaring
            power = np.zeros(max_terms)
            power[0] = 1.0
            base = temp.copy()

            for bit in range(5):  # 24 = 11000 in binary
                if (24 >> bit) & 1:
                    new_power = np.convolve(power, base)[:max_terms]
                    power = new_power
                base = np.convolve(base, base)[:max_terms]

            result = np.convolve(result, power)[:max_terms]

        return result

    N = 200
    tau = ramanujan_tau_first_N(N)

    # Show tau at prime indices
    primes_small = [int(prime(n)) for n in range(1, 50)]
    print(f"\nRamanujan τ(p) for first primes:")
    for p in primes_small[:15]:
        if p < N:
            print(f"  τ({p:3d}) = {tau[p]:+12.0f}")

    # Known values for verification
    print(f"\n  Verification: τ(2) should be -24: got {tau[2]:.0f}")
    print(f"  τ(3) should be 252: got {tau[3]:.0f}")
    print(f"  τ(5) should be 4830: got {tau[5]:.0f}")
    print(f"  τ(7) should be -16744: got {tau[7]:.0f}")

    # Correlation with prime index
    prime_indices = []
    tau_at_primes = []
    for i, p in enumerate(primes_small):
        if p < N:
            prime_indices.append(i + 1)
            tau_at_primes.append(tau[p])

    corr = np.corrcoef(prime_indices, tau_at_primes)[0, 1]
    print(f"\n  Correlation of τ(p) with prime index: {corr:.4f}")
    print(f"  (Would need strong correlation to be useful; ~0 means useless)")

    # τ(p) satisfies |τ(p)| ≤ 2*p^{11/2} (Deligne's theorem)
    # Normalized: τ(p)/p^{11/2} follows Sato-Tate distribution
    normalized_tau = [tau[p] / (p**5.5) for p in primes_small if p < N]
    print(f"\n  τ(p)/p^{{11/2}} distribution: mean={np.mean(normalized_tau):.4f}, "
          f"std={np.std(normalized_tau):.4f}")

    print("\n  CONCLUSION: τ(p) is essentially random (Sato-Tate),")
    print("  uncorrelated with prime index. Modular form coefficients")
    print("  cannot help locate the nth prime.")

def main():
    print("Session 6: Elliptic Curve / Modular Form Approach")
    print(f"Started: {time.strftime('%H:%M:%S')}\n")

    t0 = time.time()
    experiment_a_p_structure()
    experiment_modular_form_coefficients()
    elapsed = time.time() - t0
    print(f"\n{'='*70}")
    print(f"Total time: {elapsed:.1f}s")

if __name__ == '__main__':
    main()

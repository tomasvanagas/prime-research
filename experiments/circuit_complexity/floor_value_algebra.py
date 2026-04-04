"""
Floor Value Algebraic Structure
Session 12 - April 2026

Key observation: The set V = {floor(x/k) : k=1,...,sqrt(x)} has size O(sqrt(x)).
The Lucy DP computes a function of V by iterating over primes.

Novel question: does V have ALGEBRAIC CLOSURE properties that could be exploited?

Specifically:
1. Is V closed under floor-division? (i.e., if v in V and p prime, is floor(v/p) in V?)
   YES - this is what makes the Lucy DP work.

2. Does V form a LATTICE under divisibility?
3. Can the Meissel-Lehmer formula be expressed as a polynomial in the values of V?
4. Is there a RECURRENCE relation among values of V that avoids the sieve?

If we could compute pi(x) as a LINEAR COMBINATION of {v-1 : v in V}
(which are the initial values of the Lucy DP), this would be a
non-sieve approach.
"""

import math
import sympy
from collections import defaultdict

def get_floor_values(x):
    vals = set()
    k = 1
    while k * k <= x:
        vals.add(x // k)
        vals.add(k)
        k += 1
    return sorted(vals)

def check_closure(x):
    """Check closure of V under floor-division by primes."""
    V = set(get_floor_values(x))
    sqrtx = int(math.isqrt(x))
    primes = list(sympy.primerange(2, sqrtx + 1))

    closed = True
    for v in sorted(V):
        for p in primes:
            vp = v // p
            if vp > 0 and vp not in V:
                closed = False
                # print(f"  NOT closed: floor({v}/{p}) = {vp} not in V")
                break
        if not closed:
            break

    print(f"x={x}: V closed under floor-division by primes? {closed}")
    return closed

def analyze_linear_dependence(x):
    """
    After the Lucy DP, S[v] for each v is determined.
    Can we express pi(x) = S[x] as a LINEAR combination of the
    initial values {v-1 : v in V}?

    If so, what are the coefficients?

    S[v, final] = sum_{w in V} c_{v,w} * (w - 1)  + sum of affine terms

    Finding the coefficients c_{x,w} would give pi(x) = sum c_{x,w} * (w-1) + const.
    """
    V = get_floor_values(x)
    sqrtx = int(math.isqrt(x))
    n = len(V)

    # Index mapping
    v_to_idx = {v: i for i, v in enumerate(V)}

    # Initialize: S[v] = v - 1 (identity + offset)
    # Represent S[v] as a vector of coefficients over the initial values
    # S[v, initial] = 1 * (v - 1)
    # coeffs[v] = dict mapping w -> coefficient of (w-1) in expression of S[v]
    # Plus a constant term

    coeffs = {v: defaultdict(float) for v in V}
    const = {v: 0.0 for v in V}

    for v in V:
        coeffs[v][v] = 1.0  # S[v] = 1*(v-1), so coeff of (v-1) is 1

    # Apply Lucy DP
    for p in range(2, sqrtx + 1):
        if coeffs[p][p] * ((p-1) - 1) + const[p] > coeffs[p-1][p-1] * ((p-1-1) - 1) + const[p-1]:
            # p is prime in the current sieve state
            # S[v, p] = S[v, p-1] - S[floor(v/p), p-1] + S[p-1, p-1]
            # In coefficient form:
            # coeffs[v] -= coeffs[floor(v/p)] + coeffs[p-1] ... wait this is wrong

            # Actually: S[v] -= (S[v//p] - S[p-1])
            # So new_coeffs[v][w] = coeffs[v][w] - coeffs[v//p][w] + coeffs[p-1][w]
            # And new_const[v] = const[v] - const[v//p] + const[p-1]

            # Need to check if p is prime by comparing S[p] > S[p-1]
            pass

    # Simpler approach: just compute S directly and check if pi(x) is a
    # specific linear combination of V values

    S = {}
    for v in V:
        S[v] = v - 1

    # Record initial values
    initial = {v: v - 1 for v in V}

    for p in range(2, sqrtx + 1):
        if S[p] > S[p - 1]:  # p is prime
            for v in reversed(V):
                if v < p * p:
                    break
                vp = v // p
                S[v] -= S[vp] - S[p - 1]

    pi_x = S[x]
    print(f"\nx = {x}, pi(x) = {pi_x}")

    # Can we express pi(x) as a simple function of the floor values?
    # E.g., pi(x) = sum_{v in V} c_v * v for some coefficients c_v?

    # Check: is pi(x) a Mobius-type sum over V?
    # pi(x) = sum_{v in V} mu_V(x, v) * (v - 1) ?
    # where mu_V is some Mobius function on V

    # Actually, the Lucy DP IS such a sum. The final S[x] is obtained
    # by a linear transformation of the initial values (v-1).
    # The transformation matrix is determined by the primes up to sqrt(x).

    # Let me compute this explicitly for small x.
    if n <= 50:
        import numpy as np

        # Perturbation method: set initial S[w] = (w-1) + epsilon * delta_{w,w0}
        # for each w0 in V. The change in S[x] gives the coefficient of (w0-1).

        base_S = {}
        for v in V:
            base_S[v] = v - 1
        for p in range(2, sqrtx + 1):
            if base_S[p] > base_S[p - 1]:
                for v in reversed(V):
                    if v < p * p:
                        break
                    base_S[v] -= base_S[v // p] - base_S[p - 1]

        # For each w0, perturb and see effect on S[x]
        epsilon = 0.001
        linear_coeffs = {}
        for w0 in V:
            perturbed_S = {}
            for v in V:
                perturbed_S[v] = (v - 1) + (epsilon if v == w0 else 0)

            for p in range(2, sqrtx + 1):
                if perturbed_S[p] > perturbed_S[p - 1]:
                    for v in reversed(V):
                        if v < p * p:
                            break
                        perturbed_S[v] -= perturbed_S[v // p] - perturbed_S[p - 1]

            coeff = (perturbed_S[x] - base_S[x]) / epsilon
            if abs(coeff) > 1e-6:
                linear_coeffs[w0] = round(coeff, 2)

        print(f"  Non-zero coefficients for pi({x}) as linear combination of initial values:")
        for w, c in sorted(linear_coeffs.items()):
            if abs(c) > 0.01:
                print(f"    coeff[{w}] = {c}")

        # Verify
        pi_check = sum(c * (w - 1) for w, c in linear_coeffs.items())
        print(f"  Verification: sum = {pi_check:.2f}, pi({x}) = {pi_x}")

        # How many non-zero coefficients?
        nonzero = sum(1 for c in linear_coeffs.values() if abs(c) > 0.01)
        print(f"  Number of non-zero coefficients: {nonzero} out of {n} floor values")

def analyze_coefficient_structure():
    """
    For various x, compute the coefficients of pi(x) as a linear combination
    of initial floor values, and look for patterns.
    """
    print("\n=== Coefficient Structure Analysis ===")

    for x in [100, 500, 1000, 5000]:
        V = get_floor_values(x)
        sqrtx = int(math.isqrt(x))
        n = len(V)

        # Compute base Lucy DP
        S_base = {}
        for v in V:
            S_base[v] = v - 1
        for p in range(2, sqrtx + 1):
            if S_base[p] > S_base[p - 1]:
                for v in reversed(V):
                    if v < p * p:
                        break
                    S_base[v] -= S_base[v // p] - S_base[p - 1]

        pi_x = S_base[x]

        # Perturbation to find coefficients
        epsilon = 0.0001
        coeffs = {}
        for w0 in V:
            S_pert = {}
            for v in V:
                S_pert[v] = (v - 1) + (epsilon if v == w0 else 0)
            for p in range(2, sqrtx + 1):
                if S_pert[p] > S_pert[p - 1]:
                    for v in reversed(V):
                        if v < p * p:
                            break
                        S_pert[v] -= S_pert[v // p] - S_pert[p - 1]
            c = (S_pert[x] - pi_x) / epsilon
            if abs(c) > 0.001:
                coeffs[w0] = round(c, 4)

        nonzero = sum(1 for c in coeffs.values() if abs(c) > 0.001)
        positive = sum(1 for c in coeffs.values() if c > 0.001)
        negative = sum(1 for c in coeffs.values() if c < -0.001)

        # Coefficient magnitude distribution
        magnitudes = sorted([abs(c) for c in coeffs.values()], reverse=True)

        print(f"\nx = {x}, pi(x) = {pi_x}, |V| = {n}")
        print(f"  Non-zero coefficients: {nonzero} ({nonzero/n*100:.1f}% of V)")
        print(f"  Positive: {positive}, Negative: {negative}")
        if magnitudes:
            print(f"  Largest |coeff|: {magnitudes[0]:.4f}")
            print(f"  Smallest |coeff|: {magnitudes[-1]:.4f}")
            print(f"  Sum of absolute coefficients: {sum(magnitudes):.4f}")

        # Top 10 coefficients
        top = sorted(coeffs.items(), key=lambda x: abs(x[1]), reverse=True)[:10]
        print(f"  Top 10 coefficients:")
        for w, c in top:
            print(f"    v={w}: c={c:+.4f} (initial value = {w-1})")

def main():
    # Closure check
    for x in [100, 1000, 10000]:
        check_closure(x)

    # Linear dependence analysis
    for x in [100, 500, 1000]:
        analyze_linear_dependence(x)

    # Coefficient structure
    analyze_coefficient_structure()

    print("\n" + "=" * 70)
    print("CONCLUSIONS")
    print("=" * 70)
    print("""
The Lucy DP computes pi(x) as a LINEAR transformation of the initial
values {v-1 : v in V}. The transformation matrix encodes all the
sieving operations.

Key question: does this linear transformation have LOW RANK or other
structure that could be exploited for faster computation?

If the transformation has rank r << |V|, then pi(x) depends on only
r independent linear combinations of the initial values, which could
potentially be computed in O(r * polylog) time.

If the transformation matrix is FULL RANK, then all |V| = O(sqrt(x))
initial values are needed, confirming the sieve barrier.
""")

if __name__ == "__main__":
    main()

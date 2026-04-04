#!/usr/bin/env python3
"""
Session 13: Can pi(x) be expressed as a low-degree integer polynomial
in the bits of x?

Over GF(2), the prime indicator has degree Theta(N). But over Z (integers),
the representation might be different because integer arithmetic is more
powerful than GF(2) arithmetic.

Key insight: the prime indicator function maps {0,1}^N -> {0,1}.
Over Z, any such function can be written as a MULTILINEAR polynomial:
  f(x_1,...,x_N) = sum_{S} c_S * prod_{i in S} x_i

where c_S are INTEGER (or rational) coefficients.

The multilinear representation over Z is UNIQUE for Boolean functions.
The DEGREE is the max |S| with c_S != 0.

For the prime indicator: what is this degree over Z?
And what is the magnitude of the coefficients?

If degree = O(polylog(N)): then pi(x) = sum of poly(N)^d evaluations,
potentially computable in NC.
If degree = Theta(N): same barrier as GF(2).
"""

import numpy as np
from sympy import isprime
from itertools import combinations

def compute_fourier_over_Z(f_table, N):
    """
    Compute the multilinear representation of f over Z.

    For f: {0,1}^N -> Z, the unique multilinear representation is:
    f(x) = sum_{S ⊆ [N]} c_S * prod_{i ∈ S} x_i

    where c_S = sum_{T ⊆ S} (-1)^{|S|-|T|} f(T)
    (Mobius inversion on the Boolean lattice)

    This is the same as the ANF over GF(2) but with integer coefficients.
    """
    c = np.zeros(2**N, dtype=np.float64)

    for S_idx in range(2**N):
        S_bits = [i for i in range(N) if S_idx & (1 << i)]
        val = 0.0
        # Sum over all subsets T of S
        for mask in range(2**len(S_bits) + 1):
            if mask >= 2**len(S_bits):
                break
            T_idx = 0
            T_size = 0
            for j, bit in enumerate(S_bits):
                if mask & (1 << j):
                    T_idx |= (1 << bit)
                    T_size += 1
            sign = (-1) ** (len(S_bits) - T_size)
            val += sign * f_table[T_idx]
        c[S_idx] = val

    return c

def experiment_integer_degree():
    """What is the degree of pi(x) as an integer polynomial in bits of x?"""
    print("="*70)
    print("EXPERIMENT: Integer Multilinear Representation of Prime Indicator")
    print("="*70)

    for N in range(3, 17):
        # Truth table
        f = np.array([1 if isprime(n) else 0 for n in range(2**N)], dtype=np.float64)

        # Compute integer multilinear coefficients
        c = compute_fourier_over_Z(f, N)

        # Analyze
        max_deg = 0
        nonzero = 0
        max_coeff = 0
        deg_counts = {}
        coeff_by_deg = {}

        for idx in range(2**N):
            if abs(c[idx]) > 0.5:  # Integer, so threshold at 0.5
                deg = bin(idx).count('1')
                max_deg = max(max_deg, deg)
                nonzero += 1
                max_coeff = max(max_coeff, abs(c[idx]))
                deg_counts[deg] = deg_counts.get(deg, 0) + 1
                coeff_by_deg.setdefault(deg, []).append(abs(c[idx]))

        print(f"\nN = {N}: degree = {max_deg}/{N}, nonzero = {nonzero}/{2**N} ({100*nonzero/2**N:.1f}%)")
        print(f"  Max |coeff| = {max_coeff:.0f}")

        # Show coefficient statistics by degree
        for d in sorted(coeff_by_deg.keys()):
            coeffs = coeff_by_deg[d]
            print(f"  deg {d}: {len(coeffs)} terms, max|c|={max(coeffs):.0f}, "
                  f"mean|c|={np.mean(coeffs):.1f}")

    print("\n" + "="*70)
    print("KEY QUESTION: Over Z, is the degree the same as over GF(2)?")
    print("If coefficients at high degree are all EVEN, the GF(2) degree")
    print("could be lower than the Z degree, but this seems unlikely for primes.")
    print("="*70)

def experiment_pi_as_polynomial():
    """
    Instead of the prime INDICATOR, look at pi(x) itself.

    pi(x) = sum_{n=2}^{x} 1_P(n)

    If we write x in binary as x = sum x_i * 2^i, then pi(x) is a function
    of x_1, ..., x_N. It's NOT a multilinear polynomial (because x ranges
    over all 2^N values, not just {0,1}^N in the usual sense).

    But we CAN study the INCREMENTS:
    Delta_i pi(x) = pi(x XOR 2^i) - pi(x) (flip bit i)

    If pi(x) had low Z-degree d as a function of the bits:
    Delta_{i1} Delta_{i2} ... Delta_{id+1} pi(x) = 0 for all x.
    """
    print("\n" + "="*70)
    print("EXPERIMENT: Finite Differences of pi(x)")
    print("="*70)

    for N in [8, 10, 12, 14]:
        # Compute pi for all N-bit numbers
        pi_table = np.array([int(sum(1 for k in range(2, n+1) if isprime(k)))
                            if n >= 2 else 0 for n in range(2**N)])

        # First-order differences: Delta_i pi(x) = pi(x ^ 2^i) - pi(x)
        # For each bit position, compute the max absolute difference
        max_first = np.zeros(N)
        for i in range(N):
            mask = 1 << i
            for x in range(2**N):
                diff = abs(pi_table[x ^ mask] - pi_table[x])
                max_first[i] = max(max_first[i], diff)

        print(f"\nN = {N}:")
        print(f"  Max |Delta_i pi(x)| for each bit position:")
        for i in range(N):
            print(f"    bit {i} (weight 2^{i}={2**i}): max diff = {int(max_first[i])}")

        # Second-order: Delta_i Delta_j pi(x) for i != j
        # Sample a few pairs
        max_second = 0
        for i in range(min(N, 8)):
            for j in range(i+1, min(N, 8)):
                mask_i = 1 << i
                mask_j = 1 << j
                for x in range(min(2**N, 10000)):
                    diff2 = (pi_table[x ^ mask_i ^ mask_j]
                            - pi_table[x ^ mask_i]
                            - pi_table[x ^ mask_j]
                            + pi_table[x])
                    max_second = max(max_second, abs(diff2))

        print(f"  Max |Delta_i Delta_j pi(x)|: {max_second}")

        # High-order differences for small bit subset
        if N <= 12:
            # Check: are (N-1)th order differences zero? (would mean degree < N)
            # Use subset {0, 1, ..., k-1} of bit positions
            for order in range(1, min(N, 8)):
                max_diff = 0
                bits = list(range(order))
                for x in range(min(2**N, 5000)):
                    # Compute order-th difference over these bits
                    val = 0
                    for mask_val in range(2**order):
                        flip = 0
                        sign_exp = 0
                        for b in range(order):
                            if mask_val & (1 << b):
                                flip |= (1 << bits[b])
                                sign_exp += 1
                        sign = (-1)**sign_exp
                        val += sign * pi_table[x ^ flip]
                    max_diff = max(max_diff, abs(val))
                print(f"  Order-{order} diff (bits 0..{order-1}): max = {max_diff}")

        # KEY: the max difference at order k tells us about degree
        # If max = 0 at order d+1, then degree ≤ d

if __name__ == "__main__":
    experiment_integer_degree()
    experiment_pi_as_polynomial()

    print("\n" + "="*70)
    print("IMPLICATIONS")
    print("="*70)
    print("""
If the integer multilinear degree is also Theta(N):
- The prime indicator is fundamentally "high-degree" regardless of field
- No low-degree polynomial (over ANY ring) can compute primality
- This strengthens the barrier: even arithmetic circuits need Omega(N) depth

If the integer coefficients grow exponentially:
- The function has large "influence" at high degrees
- No efficient evaluation via partial sums or truncation
""")

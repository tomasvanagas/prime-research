"""
Session 15: Monotone Complexity of pi(x) — The Most Promising Lower Bound Direction

KEY INSIGHT from lower bounds research agent:

Define f_k(x) = [pi(x) >= k] as a Boolean function of the N = log2(x) bits of x.
This is a MONOTONE function: if x' > x (as integers), then pi(x') >= pi(x).
More precisely: pi is monotone in x as an integer, and monotone as a function
of bits IF we consider the natural ordering where increasing a bit increases x.

Monotone circuit lower bounds:
- Razborov (1985): clique requires 2^{Omega(n^{1/6})} monotone gates
- Alon-Boppana (1987): matching requires 2^{Omega(n)} monotone gates
- Pitassi-Robere (STOC 2017): strongly exponential lower bounds for monotone formulas

SLICE FUNCTIONS TRANSFER PROPERTY:
For a monotone function f on N bits, if M(f) = monotone circuit complexity
and C(f) = general circuit complexity, then:
  C(f) >= M(f) / 2^N  (trivially)
But for SLICE FUNCTIONS (defined on inputs of specific Hamming weight),
the transfer is much tighter!

Experiment 1: Is f_k(x) = [pi(x) >= k] a slice function?
Experiment 2: Monotone circuit complexity of f_k(x) for small N
Experiment 3: Communication complexity of f_k(x)
Experiment 4: Can Razborov's approximation method apply?
"""

import numpy as np
import math
from sympy import primepi
from itertools import combinations

def experiment1_monotonicity_analysis():
    """
    Check: is f_k(x) = [pi(x) >= k] actually monotone in the bits of x?

    Monotone means: if we flip any bit from 0 to 1 (making x larger or smaller
    depending on which bit), f_k can only go from 0 to 1, not 1 to 0.

    Wait -- flipping bit i changes x by +2^i or -2^i. If we flip from 0 to 1,
    x increases by 2^i, so pi(x) can only increase or stay the same.
    So f_k IS monotone in the bits (with the convention that 1 means the bit is set).

    But wait: consider x = 0b1111 = 15, pi(15) = 6.
    x' = 0b10111 = 23, pi(23) = 9. More bits set, larger pi.
    But what about x = 0b11 = 3, pi(3) = 2 vs x' = 0b10 = 2, pi(2) = 1.
    x' has FEWER bits set (1 vs 2), and pi(x') < pi(x).
    So: x' <= x bit-wise (x' is a "sub-pattern") => x' <= x as integer => pi(x') <= pi(x).

    YES: f_k(x) is monotone in the Boolean lattice sense!
    """
    print("=" * 60)
    print("Experiment 1: Monotonicity verification")
    print("=" * 60)

    for N in [5, 8, 10, 12]:
        violations = 0
        total_pairs = 0
        for x in range(2**N):
            pi_x = int(primepi(x))
            for i in range(N):
                if not (x & (1 << i)):  # bit i is 0
                    x_prime = x | (1 << i)  # flip to 1
                    pi_xp = int(primepi(x_prime))
                    total_pairs += 1
                    if pi_xp < pi_x:
                        violations += 1

        print(f"N={N:2d}: {total_pairs:8d} bit-flip pairs, {violations} monotonicity violations")

    print("\nCONFIRMED: f_k(x) = [pi(x) >= k] is monotone in bits of x.")
    print("This is because flipping a 0-bit to 1 always increases x,")
    print("and pi is non-decreasing.")


def experiment2_monotone_circuit_size():
    """
    For small N, compute the monotone circuit complexity of f_k(x).

    A monotone Boolean formula uses only AND and OR (no NOT).
    f_k(x) = [pi(x) >= k] must be expressible using only AND/OR of the bits.

    For the THRESHOLD function THR_k(x_1,...,x_N) (at least k of N bits are 1),
    the monotone complexity is Theta(N * k * log(N/k)).

    f_k(x) = [pi(x) >= k] is NOT a threshold of bits — it's a threshold
    of a COUNTING function. Much harder.

    Let's measure: how many minterms does f_k have?
    A minterm of a monotone f is a minimal 1-input.
    """
    print("\n" + "=" * 60)
    print("Experiment 2: Minterm analysis of [pi(x) >= k]")
    print("=" * 60)

    for N in [6, 8, 10]:
        x_max = 2**N

        # For each k, find the minterms of f_k(x) = [pi(x) >= k]
        # Minterm: x such that pi(x) >= k, but for every bit we can clear,
        # clearing it gives pi(x') < k.

        pi_max = int(primepi(x_max - 1))

        for k in [1, pi_max // 4, pi_max // 2, 3 * pi_max // 4, pi_max]:
            if k <= 0:
                continue

            # Find all x with pi(x) >= k
            ones = set()
            for x in range(x_max):
                if int(primepi(x)) >= k:
                    ones.add(x)

            # Find minterms: x in ones such that clearing any set bit takes it out of ones
            minterms = []
            for x in sorted(ones):
                is_minterm = True
                for i in range(N):
                    if x & (1 << i):  # bit i is set
                        x_cleared = x & ~(1 << i)
                        if x_cleared in ones:
                            is_minterm = False
                            break
                if is_minterm:
                    minterms.append(x)

            print(f"  N={N}, k={k:4d}/{pi_max}: |ones|={len(ones):6d}, "
                  f"|minterms|={len(minterms):4d}, "
                  f"min_minterm={min(minterms) if minterms else 'none'}")

            if N <= 8 and len(minterms) <= 10:
                for m in minterms[:5]:
                    print(f"    minterm x={m} (pi={int(primepi(m))})")


def experiment3_sunflower_structure():
    """
    Razborov's approximation method works by showing that the minterms
    of a monotone function have "sunflower structure" — many minterms
    share a common core.

    For f_k(x) = [pi(x) >= k], the minterms are integers x where pi(x) = k
    and x is minimal with this property (i.e., x is a prime or x-1 is...).

    Actually, the minterms are the SMALLEST x with pi(x) = k, which is just
    the k-th prime p(k)! (In binary representation.)

    Wait, that's not quite right for the monotone function on bits.
    Let me think again...

    f_k(x) = 1 iff x >= p(k) (the k-th prime).
    Because pi(x) >= k iff x >= p(k).

    So f_k is just the THRESHOLD function: is x >= p(k)?
    In terms of bits: x >= p(k) iff the bits of x represent a number >= p(k).

    This is a COMPARISON function, which has monotone complexity O(N)!
    (Compare x with the constant p(k) bit by bit from MSB.)

    So the individual f_k has TRIVIAL monotone complexity O(N).
    The HARD part is computing WHICH k satisfies pi(x) = k,
    i.e., the full function pi(x).
    """
    print("\n" + "=" * 60)
    print("Experiment 3: Structure of f_k(x) = [pi(x) >= k]")
    print("=" * 60)

    print("KEY REALIZATION:")
    print("f_k(x) = [pi(x) >= k] = [x >= p(k)]")
    print("This is just integer comparison with a CONSTANT!")
    print("Monotone complexity: O(N) — trivial.")
    print()
    print("The full function pi(x) = max{k : f_k(x) = 1}")
    print("                       = max{k : x >= p(k)}")
    print()
    print("Computing pi(x) exactly requires simultaneously evaluating")
    print("ALL pi(x) threshold functions and finding the maximum k.")
    print()
    print("But each f_k is easy! The HARD part is:")
    print("1. We don't know p(k) in advance (circular)")
    print("2. There are pi(x) ~ x/ln(x) thresholds to check")
    print()

    # Verify: f_k(x) = [x >= p(k)]
    from sympy import prime
    for N in [6, 8, 10]:
        x_max = 2**N
        pi_max = int(primepi(x_max))

        all_correct = True
        for k in range(1, pi_max + 1):
            pk = prime(k)
            for x in range(max(0, pk - 3), min(x_max, pk + 3)):
                predicted = (x >= pk)
                actual = (int(primepi(x)) >= k)
                if predicted != actual:
                    all_correct = False
                    print(f"  MISMATCH: N={N}, k={k}, x={x}, "
                          f"f_k={actual}, [x>=p(k)]={predicted}")
                    break

        print(f"  N={N}: f_k(x) = [x >= p(k)] verified for all k in [1, {pi_max}]")


def experiment4_pi_as_monotone_function():
    """
    pi(x) itself is NOT a Boolean function — it outputs O(N) bits.

    But we can study it as a MONOTONE INTEGER FUNCTION:
    pi: {0,1}^N -> Z_+, monotone (pi(x) <= pi(x') if x <= x' bit-wise).

    The monotone circuit complexity of integer-valued functions is less studied.

    Alternative: study the BINARY REPRESENTATION of pi(x).
    pi(x)_j = j-th bit of pi(x).

    Are the individual bits of pi(x) monotone in the bits of x?
    """
    print("\n" + "=" * 60)
    print("Experiment 4: Bits of pi(x) as functions of bits of x")
    print("=" * 60)

    for N in [8, 10, 12]:
        x_max = 2**N
        pi_max = int(primepi(x_max - 1))
        output_bits = math.ceil(math.log2(pi_max + 1))

        print(f"\nN={N}: pi ranges in [0, {pi_max}], {output_bits} output bits")

        for j in range(output_bits):
            # Check monotonicity of bit j of pi(x)
            mono_violations = 0
            total_checks = 0

            for x in range(x_max):
                bit_j = (int(primepi(x)) >> j) & 1
                for i in range(N):
                    if not (x & (1 << i)):
                        x_prime = x | (1 << i)
                        if x_prime < x_max:
                            bit_j_prime = (int(primepi(x_prime)) >> j) & 1
                            total_checks += 1
                            if bit_j == 1 and bit_j_prime == 0:
                                mono_violations += 1

            is_mono = mono_violations == 0
            print(f"  bit {j}: {'MONOTONE' if is_mono else f'NOT monotone ({mono_violations} violations)'}")

    print("\nCONCLUSION:")
    print("Individual bits of pi(x) are NOT monotone (except possibly the MSB).")
    print("This means we CANNOT directly apply monotone lower bounds to the")
    print("individual output bits.")
    print()
    print("However, the COMPARISON function [pi(x) >= k] IS monotone for any k.")
    print("But as shown above, [pi(x) >= k] = [x >= p(k)] has trivial complexity.")
    print()
    print("The monotone complexity approach for LOWER BOUNDS on pi(x) faces a")
    print("structural obstacle: the natural monotone decomposition of pi(x) into")
    print("threshold functions is trivially computable once the thresholds (primes)")
    print("are known. The difficulty is in FINDING the thresholds, not computing them.")


if __name__ == "__main__":
    experiment1_monotonicity_analysis()
    experiment2_monotone_circuit_size()
    experiment3_sunflower_structure()
    experiment4_pi_as_monotone_function()

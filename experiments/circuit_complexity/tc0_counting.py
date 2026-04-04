"""
Session 15: #TC^0 Counting — The Most Promising Remaining Framework

KEY IDEA: If PRIMES is in TC^0 (e.g., BPSW is unconditionally correct),
then pi(x) = #{n <= x : C(n) = 1} where C is a TC^0 circuit.

The question becomes: can we COUNT satisfying inputs of a specific TC^0
circuit efficiently?

Relevant facts:
- #P-hard in general (even for AC^0 circuits — Toda's theorem)
- BUT specific circuits may have structure that allows efficient counting
- The BPSW circuit has VERY specific structure: scalar powering + 2x2 MPOW + GCD
- Threshold circuits with polylog depth can be "unrolled" into polynomial-size
  arithmetic formulas

Experiment 1: Structure of the BPSW TC^0 circuit
Experiment 2: Can we count satisfying inputs of simple threshold circuits?
Experiment 3: Connection to #TC^0 counting literature
Experiment 4: Exploiting BPSW structure for batch counting
"""

import numpy as np
import math
from sympy import primepi, isprime, primerange, nextprime
from collections import defaultdict

def experiment1_bpsw_circuit_structure():
    """
    BPSW = MR(2) AND Strong Lucas(D) where D is the first Jacobi -1.

    MR(2): compute 2^(n-1) mod n, check it's 1 or -1 at each squaring step.
    Strong Lucas: compute U_d, V_d mod n via 2x2 matrix powering.
    Jacobi: compute (D/n) via GCD-like recursion.

    All in TC^0 via:
    - Scalar pow (a^b mod m): TC^0 [Allender 1999, HAB 2002]
    - 2x2 matrix pow: TC^0 [Mereghetti-Palano 2000]
    - Jacobi/GCD: TC^0 [HAB 2002]

    The circuit C(n) has:
    - Input: n (N bits)
    - Output: 1 iff BPSW says "prime"
    - Size: poly(N) gates
    - Depth: O(1) (constant depth with threshold gates)

    For COUNTING: pi(x) = sum_{n=2}^{x} C(n) = sum of 2^N evaluations
    of a constant-depth threshold circuit.

    Key question: does the arithmetic structure of BPSW allow
    the sum to be computed without evaluating each C(n)?
    """
    print("=" * 60)
    print("Experiment 1: BPSW Circuit Structure Analysis")
    print("=" * 60)

    print("""
BPSW Circuit Components:
========================

1. Miller-Rabin base 2:
   - Write n-1 = 2^s * d (d odd)
   - Compute a_0 = 2^d mod n
   - Compute a_1 = a_0^2 mod n, ..., a_s = a_{s-1}^2 mod n
   - Accept if a_0 = 1, or any a_i = n-1

   TC^0 implementation: O(N^2) multiplications mod n, each in TC^0.
   Entire chain: iterated multiplication = TC^0 (HAB 2002).

2. Strong Lucas test (Selfridge parameters):
   - Find first D in {5, -7, 9, -11, ...} with Jacobi(D,n) = -1
   - Set P=1, Q=(1-D)/4
   - Compute U_d, V_d mod n via Lucas chain (2x2 matrix powering)
   - Accept if U_d = 0 mod n, or V_{2^j * d} = 0 for some 0 <= j < s

   TC^0 implementation: 2x2 matrix powering is in TC^0 (Mereghetti-Palano).
   Jacobi symbol: TC^0 via GCD (HAB 2002).

3. Combined: BPSW(n) = MR2(n) AND StrongLucas(n)

COUNTING CHALLENGE:
==================
pi(x) = sum_{n=2}^{x} BPSW(n)

For a TC^0 circuit C of size s and depth d:
- #inputs satisfying C is in #TC^0
- #TC^0 ⊆ #P (trivially)
- Is #TC^0 ⊆ NC? UNKNOWN — this is essentially our question!

The circuit C(n) acts on the BITS of n. Each bit is an independent variable.
The circuit computes modular arithmetic on n.

KEY STRUCTURAL OBSERVATION:
The BPSW circuit has the form C(n) = f(2^{n-1} mod n, U_d(n) mod n, ...)
where the arguments depend on n in a highly nonlinear way (modular exponentiation).

The function f is simple (comparison with 1, n-1, 0).
The HARD PART is the modular exponentiation, which creates a VERY nonlinear
dependency on the bits of n.
""")

    # Verify BPSW structure for small n
    print("\nBSPW components for small n:")
    print(f"{'n':>6} {'MR2':>4} {'prime':>6} {'2^(n-1)%n':>12}")
    for n in range(2, 50):
        if n < 2:
            continue
        mr2 = pow(2, n-1, n) == 1 or pow(2, n-1, n) == n-1
        # More careful: write n-1 = 2^s * d
        d = n - 1
        s = 0
        while d % 2 == 0:
            d //= 2
            s += 1
        x = pow(2, d, n)
        mr2_strong = False
        if x == 1 or x == n - 1:
            mr2_strong = True
        else:
            for r in range(1, s):
                x = pow(x, 2, n)
                if x == n - 1:
                    mr2_strong = True
                    break

        p = isprime(n)
        if mr2_strong != p and n < 50:
            print(f"{n:6d} {str(mr2_strong):>4} {str(bool(p)):>6} {pow(2, n-1, n):12d}  <- DISCREPANCY")


def experiment2_counting_threshold_circuits():
    """
    Can we count satisfying inputs of simple threshold circuits?

    Consider the simplest primality-related threshold function:
    "Is n odd AND n > 1?" — This is a trivial threshold function.
    Count = floor(x/2) for x >= 2.

    A slightly harder one: "Is n ≡ 1 or 5 mod 6?" (not div by 2 or 3)
    Count = floor(x/3) + floor((x+1)/6) approximately.

    The pattern: for any LINEAR threshold function of bits, counting is
    easy because it's a lattice point count in a halfspace.

    For MODULAR functions (a^b mod n): counting is hard because the
    function is cryptographically pseudorandom.

    Experiment: measure the correlation between consecutive BPSW outputs
    to see if there's exploitable structure in the counting.
    """
    print("\n" + "=" * 60)
    print("Experiment 2: Counting threshold circuit structure")
    print("=" * 60)

    # For counting linear threshold functions: trivial
    # sum_{n=1}^{x} [w^T bits(n) >= theta]
    # = number of integers in {1,...,x} in a halfspace of bit representation
    # This is a lattice point count in a polytope — Barvinok gives poly(N) for fixed dim.
    # But the "dimension" is N (number of bits), so Barvinok doesn't help directly.

    # However, for the specific halfspace defined by a single threshold gate:
    # sum_{n=0}^{x} [n >= k] = max(0, x - k + 1)
    # This is trivially O(1).

    # For conjunction of two thresholds:
    # sum [n >= a AND n <= b] = max(0, min(x, b) - a + 1)
    # Still O(1).

    # For MODULAR predicates: sum_{n=1}^{x} [n mod m == r] = floor((x-r)/m) + 1
    # O(1) as well!

    # The hard part: MODULAR EXPONENTIATION introduces extreme nonlinearity.
    # pow(2, n-1, n) is NOT a simple function of n's bits.

    # Let's measure autocorrelation of the BPSW output
    x_max = 10000
    results = [int(isprime(n)) for n in range(2, x_max + 1)]

    # Autocorrelation at lag k
    print("\nAutocorrelation of prime indicator (lag analysis):")
    mean = np.mean(results)
    var = np.var(results)
    for lag in [1, 2, 3, 4, 5, 6, 10, 30]:
        if lag >= len(results):
            break
        corr = np.corrcoef(results[:-lag], results[lag:])[0, 1]
        print(f"  lag {lag:3d}: correlation = {corr:+.6f}")

    # Joint distribution of consecutive primality indicators
    print("\nJoint distribution of (is_prime(n), is_prime(n+2)):")
    joint = defaultdict(int)
    for i in range(len(results) - 2):
        joint[(results[i], results[i+2])] += 1
    for key in sorted(joint.keys()):
        print(f"  {key}: {joint[key]:6d} ({joint[key]/len(results)*100:.2f}%)")

    # The point: prime indicators are essentially uncorrelated
    # This means the counting problem has NO exploitable short-range structure


def experiment3_modexp_batch_structure():
    """
    Can modular exponentiation be "batched" across consecutive n?

    For MR base 2: we compute 2^{n-1} mod n for n = 2, 3, ..., x.
    Is there ANY shared structure between these computations?

    Key observation: 2^{n-1} mod n depends on the FACTORIZATION of n
    (by CRT and Fermat's little theorem).

    For prime p: 2^{p-1} mod p = 1 (Fermat)
    For composite n = p*q: 2^{n-1} mod n depends on 2^{n-1} mod p and mod q
      = 2^{(n-1) mod (p-1)} mod p * ... (CRT)

    The key: (n-1) mod (p-1) depends on n mod p and p itself.
    For CONSECUTIVE n: n mod p cycles with period p.
    So 2^{n-1} mod p has period p*(p-1) (or a divisor).

    Could we use Fourier analysis over these periods?
    """
    print("\n" + "=" * 60)
    print("Experiment 3: Batch modular exponentiation structure")
    print("=" * 60)

    # For each small prime p, compute the period of 2^{n-1} mod p
    print("Period of 2^{n-1} mod p for small primes p:")
    for p in [3, 5, 7, 11, 13, 17, 19, 23, 29, 31]:
        values = [pow(2, n-1, p) for n in range(1, 200)]
        # Find period
        for period in range(1, 200):
            is_period = True
            for j in range(period, min(period + 50, len(values))):
                if values[j] != values[j - period]:
                    is_period = False
                    break
            if is_period:
                break
        else:
            period = "none found"

        print(f"  p={p:3d}: period = {period}, "
              f"first values: {values[:min(10, int(period) if isinstance(period, int) else 10)]}")

    # For composite n, the situation is more complex
    # 2^{n-1} mod n for consecutive n:
    print("\n2^{n-1} mod n for n = 2..50:")
    print(f"{'n':>4} {'2^(n-1)%n':>10} {'prime':>6} {'Fermat=1':>9}")
    for n in range(2, 51):
        val = pow(2, n-1, n)
        p = isprime(n)
        print(f"{n:4d} {val:10d} {str(bool(p)):>6} {str(val == 1):>9}")


def experiment4_tc0_counting_theory():
    """
    Theoretical analysis of #TC^0 counting.

    Key question: Is #TC^0 ⊆ FP (function P)?

    Known results:
    - #AC^0 (counting satisfying assignments of AC^0 circuits) is in FP
      for bounded-depth AND/OR circuits (Razborov-type arguments)
    - #TC^0 is NOT known to be in FP
    - Allender (1999): certain counting problems in TC^0 reduce to
      iterated multiplication, which is in #L ⊆ NC^2

    For our specific circuit (BPSW):
    - BPSW is a depth-O(1) TC^0 circuit of polynomial size
    - Its structure: threshold gates → modular arithmetic → comparison
    - The modular arithmetic creates the nonlinearity

    Could the COUNTING of BPSW be easier than general #TC^0?
    The BPSW circuit has the specific form:
    C(n) = [2^{d(n)} ≡ 1 or n-1 mod n] AND [U_{d(n)} ≡ 0 mod n]
    where d(n) = (n-1) / 2^{v_2(n-1)}

    The modular exponentiation 2^d mod n is NOT a polynomial in the bits of n.
    It's a TOWER of multiplications — reducible to TC^0 by HAB but with
    massive cancellation in the reduction.
    """
    print("\n" + "=" * 60)
    print("Experiment 4: #TC^0 Counting Theory")
    print("=" * 60)

    print("""
THEORETICAL ANALYSIS:
====================

Hierarchy of counting complexity:

#AC^0 ⊆ #TC^0 ⊆ #NC^1 ⊆ #L ⊆ #P

Key results:
1. #AC^0 ⊆ TC^0 (Allender-Gore 1993): counting AC^0 circuits can be done
   by TC^0 circuits. So the count ITSELF is in TC^0.

2. For TC^0 circuits: the count of satisfying assignments is NOT known
   to be in any low class. Even counting threshold gate outputs is hard.

3. The specific structure of BPSW:
   - BPSW(n) depends on pow(2, n-1, n) — modular exponentiation
   - pow(a, b, m) is in TC^0 when a, b, m are all given as input
   - But for COUNTING, we need sum_{n=2}^x [pow(2,n-1,n) passes test]
   - Here a=2 is FIXED, b=n-1 depends on n, and m=n also depends on n
   - The coupling between b and m (both functions of n) is the source of difficulty

4. FIXED-BASE modular exponentiation:
   - f(n) = 2^{n-1} mod n as a function of n
   - This is related to the Euler quotient q(n) = (2^{phi(n)} - 1) / n
   - The distribution of f(n) values is NOT understood

5. TOTIENT APPROACH:
   - By Fermat: 2^{n-1} mod n = 2^{(n-1) mod phi(n)} mod n * (something)
   - But phi(n) requires factoring n!
   - For primes: phi(p) = p-1, so 2^{p-1} mod p = 1 (trivially)
   - The point: the identity holds ONLY for primes, not for all n

CONCLUSION:
The #TC^0 counting approach for BPSW faces the same fundamental barrier:
the modular exponentiation 2^{n-1} mod n couples the exponent and modulus
through n's factorization. Counting how many n pass this test requires
understanding the distribution of Fermat residues, which is equivalent
to understanding the distribution of primes.

The approach is NOT circular (it doesn't need primes to count primes),
but it IS equivalent in difficulty to the original problem.

REFINED QUESTION: Is #TC^0 ⊆ NC? If YES, then pi(x) ∈ NC (since BPSW ∈ TC^0).
This is a COMPLEXITY THEORY question, independent of number theory.
""")

    # Empirical: distribution of 2^{n-1} mod n
    print("\nDistribution of Fermat residues 2^{n-1} mod n:")
    residue_counts = defaultdict(int)
    for n in range(2, 10001):
        r = pow(2, n-1, n)
        residue_counts[r] += 1

    # Most common residues
    sorted_residues = sorted(residue_counts.items(), key=lambda x: -x[1])
    print(f"{'residue':>10} {'count':>8} {'fraction':>10}")
    for r, c in sorted_residues[:20]:
        print(f"{r:10d} {c:8d} {c/9999:10.6f}")

    print(f"\nTotal distinct residues: {len(residue_counts)}")
    print(f"Residue = 1 (Fermat passes): {residue_counts[1]} "
          f"({residue_counts[1]/9999*100:.2f}%)")
    print(f"Number of primes in [2, 10000]: {int(primepi(10000))}")
    print(f"Ratio (Fermat passes) / primes: {residue_counts[1] / int(primepi(10000)):.4f}")


if __name__ == "__main__":
    experiment1_bpsw_circuit_structure()
    experiment2_counting_threshold_circuits()
    experiment3_modexp_batch_structure()
    experiment4_tc0_counting_theory()

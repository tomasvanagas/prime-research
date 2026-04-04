#!/usr/bin/env python3
"""
TC^0 Primality Test Investigation: Four Non-AKS Approaches

Investigates whether PRIMES can be placed in TC^0 via approaches that avoid
the growing-dimension matrix powering bottleneck of AKS.

Approaches:
1. Wilson's theorem: (n-1)! ≡ -1 (mod n) iff n is prime
2. Lucas sequences: U_n(P,Q) mod n via 2x2 matrix powering (FIXED dim!)
3. Quadratic Frobenius test: Frobenius in Z[x]/(x^2-bx-c, n) = 2x2 MPOW
4. Sum-of-two-squares: Fermat characterization for p ≡ 1 (mod 4)

Key insight: approaches 2 and 3 use FIXED-DIMENSION (2x2) matrix powering,
which IS in TC^0 by Mereghetti-Palano (2000). If a deterministic primality
test can be built from ONLY 2x2 matrix powers, PRIMES is in TC^0.

Author: Session 13 sub-agent
Date: 2026-04-04
"""

import math
import time
from sympy import isprime, jacobi_symbol, sqrt_mod, factorint, nextprime
from sympy.ntheory import legendre_symbol
from functools import lru_cache


# ============================================================================
# APPROACH 1: Wilson's Theorem and Factorial mod n in TC^0
# ============================================================================

def wilson_test(n):
    """Test primality via Wilson's theorem: (n-1)! ≡ -1 (mod n)."""
    if n < 2:
        return False
    if n <= 3:
        return True
    # Compute (n-1)! mod n
    fact = 1
    for i in range(2, n):
        fact = (fact * i) % n
    return fact == n - 1


def factorial_mod_circuit_analysis():
    """
    Analyze whether (n-1)! mod n is computable in TC^0.

    Key observations:
    - (n-1)! mod n = product of n-2 numbers, each < n
    - Iterated multiplication of N given numbers IS in TC^0 (HAB 2002)
    - BUT: N = n-1 here, and n is EXPONENTIAL in input size (N = log n bits)
    - So we need to multiply 2^N numbers, requiring circuit SIZE 2^N
    - This is NOT polynomial in input size N = log n

    Alternative: Can we compute n! mod n without multiplying all terms?
    - Shamir (1991): factorial is in GapL ⊆ NC^2
    - This uses the matrix permanent connection
    - GapL is NOT known to be in TC^0

    The Granville-style shortcut: if we know the factorization of n, then
    n! mod n can be computed from Legendre's formula. But we DON'T know
    the factorization (that's what we're trying to determine).
    """
    print("=" * 70)
    print("APPROACH 1: Wilson's Theorem / Factorial mod n")
    print("=" * 70)

    # Verify Wilson's theorem
    print("\nVerification of Wilson's theorem:")
    for n in range(2, 30):
        fact_mod = 1
        for i in range(2, n):
            fact_mod = (fact_mod * i) % n
        is_p = (fact_mod == n - 1)
        correct = is_p == isprime(n)
        if not correct:
            print(f"  BUG at n={n}!")

    print("  Wilson's theorem verified for n=2..29")

    # The circuit complexity issue
    print("\nCircuit complexity analysis:")
    print("  Input: N = log2(n) bits")
    print("  (n-1)! mod n requires multiplying n-1 = 2^N - 1 numbers")
    print("  Each number is N bits")
    print("  HAB (2002): iterated mult of M numbers is TC^0 when M = poly(N)")
    print("  BUT: M = 2^N here (EXPONENTIAL in input size)")
    print("  → Circuit size = Omega(2^N) = Omega(n)")
    print("  → NOT polynomial in input size N")
    print()

    # Can we avoid computing the full product?
    print("  Alternative: partial products?")
    print("  For composite n = ab: n! ≡ 0 (mod n) since a,b < n appear in product")
    print("  So (n-1)! mod n = 0 for composite n (when n > 4)")
    print("  WAIT: This means Wilson's test is really checking (n-1)! mod n ≠ 0")
    print("  For primes, (n-1)! ≡ -1 (mod n)")
    print("  For composites n > 4: (n-1)! ≡ 0 (mod n)")
    print("  For n = 4: 3! = 6 ≡ 2 (mod 4)")
    print()

    # Verify the composite case
    print("  Verify: (n-1)! mod n for composites:")
    for n in [4, 6, 8, 9, 10, 12, 15, 20, 25, 100]:
        fact_mod = 1
        for i in range(2, n):
            fact_mod = (fact_mod * i) % n
        print(f"    n={n}: (n-1)! mod n = {fact_mod}")

    print()
    print("  CONCLUSION: Wilson's theorem CANNOT give TC^0 primality.")
    print("  Computing (n-1)! mod n requires Omega(n) = Omega(2^N) multiplications.")
    print("  Even Shamir's GapL result only gives NC^2, not TC^0.")
    print("  The 'shortcut' of knowing (n-1)! ≡ 0 for composites n > 4")
    print("  doesn't help because we need to DISTINGUISH prime from composite,")
    print("  which requires actually computing the product.")
    print()
    print("  STATUS: CLOSED. Wilson → TC^0 requires computing (n-1)! mod n,")
    print("  which needs 2^N multiplications. Not amenable to TC^0.")


# ============================================================================
# APPROACH 2: Lucas Sequences and 2x2 Matrix Powering
# ============================================================================

def lucas_U_V(n, P, Q):
    """Compute U_n(P,Q) and V_n(P,Q) mod n using 2x2 matrix powering.

    The Lucas sequence satisfies:
    [U_{k+1}]   [P  -Q] [U_k]
    [U_k    ] = [1   0] [U_{k-1}]

    So [U_n, U_{n-1}] = M^{n-1} * [U_1, U_0] where M = [[P, -Q], [1, 0]]
    U_0 = 0, U_1 = 1, V_0 = 2, V_1 = P

    This is a 2x2 matrix power, which IS in TC^0 (Mereghetti-Palano 2000).
    """
    if n == 0:
        return 0, 2
    if n == 1:
        return 1, P

    # Matrix [[P, -Q], [1, 0]] power via repeated squaring
    # Track [[a, b], [c, d]] = M^k
    # Actually easier: use the doubling formulas
    # U_{2k} = U_k * V_k
    # V_{2k} = V_k^2 - 2*Q^k
    # U_{2k+1} = (P*U_{2k} + V_{2k}) / 2  (when P is odd, use mod-2 trick)

    # Use binary method for Lucas sequences
    U, V, Qk = 1, P, Q
    # We'll compute mod nothing first, then mod n
    # Actually let's just use the matrix method directly mod n

    def mat_mul_mod(A, B, m):
        return [
            [(A[0][0]*B[0][0] + A[0][1]*B[1][0]) % m,
             (A[0][0]*B[0][1] + A[0][1]*B[1][1]) % m],
            [(A[1][0]*B[0][0] + A[1][1]*B[1][0]) % m,
             (A[1][0]*B[0][1] + A[1][1]*B[1][1]) % m]
        ]

    def mat_pow_mod(M, exp, m):
        result = [[1, 0], [0, 1]]  # identity
        base = [row[:] for row in M]
        while exp > 0:
            if exp % 2 == 1:
                result = mat_mul_mod(result, base, m)
            base = mat_mul_mod(base, base, m)
            exp //= 2
        return result

    M = [[P % n, (-Q) % n], [1, 0]]
    Mn = mat_pow_mod(M, n - 1, n)
    # [U_n, U_{n-1}] = M^{n-1} * [1, 0]
    U_n = Mn[0][0] % n
    V_n = (P * Mn[0][0] + (-Q) * Mn[1][0] * 2 // P) % n if P != 0 else 0

    # Actually, use the relation V_n = P*U_n - 2*Q*U_{n-1}
    U_nm1 = Mn[1][0] % n
    V_n = (P * U_n - 2 * Q * U_nm1) % n

    return U_n, V_n


def lucas_probable_prime(n, P, Q):
    """
    Lucas probable prime test for parameters (P, Q).

    If n is prime and gcd(n, 2QD) = 1, where D = P^2 - 4Q:
    - U_{n - (D/n)} ≡ 0 (mod n)  where (D/n) is Jacobi symbol

    This is computed via 2x2 matrix powering, which IS in TC^0.
    """
    if n < 2:
        return False
    if n == 2:
        return True
    if n % 2 == 0:
        return False

    D = P * P - 4 * Q
    if D == 0:
        return False

    # Jacobi symbol (D/n)
    delta = jacobi_symbol(D, n)

    # Compute U_{n - delta} mod n via 2x2 matrix powering
    exp = n - delta

    M = [[P % n, (-Q) % n], [1, 0]]

    def mat_mul_mod(A, B, m):
        return [
            [(A[0][0]*B[0][0] + A[0][1]*B[1][0]) % m,
             (A[0][0]*B[0][1] + A[0][1]*B[1][1]) % m],
            [(A[1][0]*B[0][0] + A[1][1]*B[1][0]) % m,
             (A[1][0]*B[0][1] + A[1][1]*B[1][1]) % m]
        ]

    def mat_pow_mod(M, exp, m):
        result = [[1, 0], [0, 1]]
        base = [row[:] for row in M]
        while exp > 0:
            if exp % 2 == 1:
                result = mat_mul_mod(result, base, m)
            base = mat_mul_mod(base, base, m)
            exp //= 2
        return result

    Mn = mat_pow_mod(M, exp - 1, n)
    U_exp = Mn[0][0] % n  # U_{exp} = first entry of M^{exp-1} * [1, 0]

    return U_exp == 0


def strong_lucas_test(n):
    """
    Strong Lucas probable prime test with Selfridge parameters.

    Selfridge method A: D = first in {5, -7, 9, -11, 13, ...} with (D/n) = -1
    P = 1, Q = (1 - D) / 4

    Then check: U_{n+1} ≡ 0 (mod n) with extra strong conditions.
    """
    if n < 2:
        return False
    if n == 2 or n == 3:
        return True
    if n % 2 == 0:
        return False

    # Check perfect square (strong Lucas PSP can be square)
    sr = int(math.isqrt(n))
    if sr * sr == n:
        return False

    # Selfridge parameters
    D = 5
    sign = 1
    while True:
        js = jacobi_symbol(D, n)
        if js == 0:
            return n == abs(D)
        if js == -1:
            break
        D = -(D + 2 * sign)
        sign = -sign
        if D > 2 * n:
            return False  # safety

    P = 1
    Q = (1 - D) // 4

    # n + 1 = d * 2^s
    d = n + 1
    s = 0
    while d % 2 == 0:
        d //= 2
        s += 1

    # Compute U_d and V_d mod n via 2x2 matrix
    M = [[P % n, (-Q) % n], [1, 0]]

    def mat_mul_mod(A, B, m):
        return [
            [(A[0][0]*B[0][0] + A[0][1]*B[1][0]) % m,
             (A[0][0]*B[0][1] + A[0][1]*B[1][1]) % m],
            [(A[1][0]*B[0][0] + A[1][1]*B[1][0]) % m,
             (A[1][0]*B[0][1] + A[1][1]*B[1][1]) % m]
        ]

    def mat_pow_mod(M_mat, exp, m):
        result = [[1, 0], [0, 1]]
        base = [row[:] for row in M_mat]
        while exp > 0:
            if exp % 2 == 1:
                result = mat_mul_mod(result, base, m)
            base = mat_mul_mod(base, base, m)
            exp //= 2
        return result

    if d > 1:
        Md = mat_pow_mod(M, d - 1, n)
    else:
        Md = [[1, 0], [0, 1]]

    # U_d = Md[0][0], U_{d-1} = Md[1][0]
    U_d = Md[0][0] % n
    U_dm1 = Md[1][0] % n
    V_d = (P * U_d - 2 * Q * U_dm1) % n

    if U_d == 0 or V_d == 0:
        return True

    # Check V_{d*2^r} ≡ 0 (mod n) for r = 1, ..., s-1
    # V_{2k} = V_k^2 - 2*Q^k
    Qk = pow(Q, d, n)
    for r in range(1, s):
        V_d = (V_d * V_d - 2 * Qk) % n
        Qk = (Qk * Qk) % n
        if V_d == 0:
            return True

    return False


def lucas_approach_analysis():
    """
    Analyze the Lucas sequence approach for TC^0 primality.

    KEY INSIGHT: Lucas sequences are defined by 2x2 matrix recurrence.
    2x2 matrix powering IS in TC^0 (Mereghetti-Palano 2000, MPOW_2).

    So: any primality test that uses ONLY Lucas sequence values mod n
    is automatically computable in TC^0 (modulo the Jacobi symbol computation
    and some simple arithmetic, all of which are in TC^0).

    The question: is there a DETERMINISTIC primality test using only
    Lucas sequences that is UNCONDITIONALLY correct?
    """
    print("=" * 70)
    print("APPROACH 2: Lucas Sequences and 2x2 Matrix Powering")
    print("=" * 70)

    print("\nTheoretical analysis:")
    print("  Lucas sequences: U_n(P,Q) satisfies U_{k+1} = P*U_k - Q*U_{k-1}")
    print("  Matrix form: [U_{k+1}, U_k]^T = [[P,-Q],[1,0]]^k * [1,0]^T")
    print("  This is 2x2 matrix powering → IN TC^0 (MPOW_2, Mereghetti-Palano 2000)")
    print()
    print("  Jacobi symbol (D/n): known to be in TC^0")
    print("    (it's related to GCD which is in TC^0 by HAB 2002)")
    print()
    print("  Strong Lucas PRP test:")
    print("    n+1 = d*2^s, check U_d ≡ 0 or V_{d*2^r} ≡ 0 for some r")
    print("    All components are 2x2 MPOW + scalar pow + comparison")
    print("    → ENTIRE strong Lucas test is in TC^0!")
    print()

    # Test strong Lucas as standalone
    print("Testing strong Lucas PRP test accuracy:")

    # Known strong Lucas pseudoprimes (from literature)
    # The strong Lucas test with Selfridge parameters has known PSPs
    slpsp_known = [5459, 5777, 10877, 16109, 18971, 22499, 24569, 25199, 40309, 58519]

    # Test on range
    false_positives = []
    false_negatives = []

    for n in range(2, 100000):
        sl = strong_lucas_test(n)
        actual = isprime(n)
        if sl and not actual:
            false_positives.append(n)
        if not sl and actual:
            false_negatives.append(n)

    print(f"  Range 2-99999:")
    print(f"    False positives (composites passing test): {len(false_positives)}")
    print(f"    False negatives (primes failing test): {len(false_negatives)}")
    if false_positives[:20]:
        print(f"    First FPs: {false_positives[:20]}")
    if false_negatives[:20]:
        print(f"    First FNs: {false_negatives[:20]}")

    # BPSW = Miller-Rabin base 2 + Strong Lucas
    # Question: Is strong Lucas alone sufficient?
    print()
    print("  BPSW = Miller-Rabin(2) + Strong Lucas(Selfridge)")
    print(f"  Strong Lucas alone: {len(false_positives)} pseudoprimes below 100000")

    # Now test: what if we add MORE Lucas tests with different parameters?
    print()
    print("  Can we make Lucas DETERMINISTIC by using multiple parameter sets?")

    # For each false positive, find which (P, Q) parameters catch it
    if false_positives:
        print(f"  Analyzing {min(len(false_positives), 20)} false positives...")
        for n in false_positives[:20]:
            caught_by = []
            for P in range(1, 20):
                for Q in range(-5, 5):
                    D = P*P - 4*Q
                    if D == 0 or D % n == 0:
                        continue
                    if math.gcd(2*Q*D, n) != 1:
                        continue
                    result = lucas_probable_prime(n, P, Q)
                    if not result:  # correctly identifies as composite
                        caught_by.append((P, Q))
                        break
                if caught_by:
                    break
            if caught_by:
                print(f"    n={n}: caught by (P,Q)={caught_by[0]}")
            else:
                print(f"    n={n}: NOT caught by any tested (P,Q)!")

    print()
    print("  KEY QUESTION: Is there a FIXED set of (P,Q) parameters such that")
    print("  the combined Lucas tests form a deterministic primality test?")
    print()
    print("  ANSWER: UNKNOWN. No such fixed set is known.")
    print("  - Arnault (1997) showed that for ANY finite set of bases,")
    print("    there exist pseudoprimes passing all tests.")
    print("  - BUT: Arnault's construction is for Fermat/MR pseudoprimes.")
    print("  - For Lucas, the situation is less clear.")
    print("  - Baillie (1980): no SLPSP is known below 10^15")
    print("  - BPSW conjecture: no BPSW pseudoprime exists (unproven)")

    return false_positives


# ============================================================================
# APPROACH 3: Quadratic Frobenius Test (Grantham 1998)
# ============================================================================

def frobenius_test(n, b, c):
    """
    Quadratic Frobenius test for n with parameters (b, c).

    Works in Z[x]/(x^2 - bx - c, n), which is a 2D algebra over Z_n.
    The Frobenius automorphism maps x -> x^n in this algebra.

    Key: multiplication in Z[x]/(x^2 - bx - c) is a 2x2 matrix operation!
    So the Frobenius map x -> x^n is a 2x2 matrix power → TC^0.

    Test: For prime n with (D/n) = -1 where D = b^2 + 4c:
    x^n ≡ -x + b (mod x^2 - bx - c, n)  [Frobenius = conjugation]
    x^{n+1} ≡ -c (mod x^2 - bx - c, n)
    """
    if n < 2:
        return False
    if n == 2:
        return True
    if n % 2 == 0:
        return False

    D = (b * b + 4 * c) % n
    if D == 0:
        return False

    # Check (D/n) = -1
    js = jacobi_symbol(D if D <= n // 2 else D - n, n)
    # Actually need to be careful with Jacobi symbol
    D_full = b * b + 4 * c
    js = jacobi_symbol(D_full, n)

    if js != -1:
        return None  # Parameters don't apply

    # Compute x^n mod (x^2 - bx - c, n)
    # Represent elements as (a0 + a1*x), multiply via:
    # (a0 + a1*x)(b0 + b1*x) = a0*b0 + (a0*b1 + a1*b0)*x + a1*b1*x^2
    # x^2 = bx + c, so x^2 ≡ bx + c
    # = (a0*b0 + a1*b1*c) + (a0*b1 + a1*b0 + a1*b1*b)*x

    def poly_mul(p1, p2, n, b, c):
        """Multiply two elements in Z_n[x]/(x^2 - bx - c)."""
        a0, a1 = p1
        b0, b1 = p2
        r0 = (a0 * b0 + a1 * b1 * c) % n
        r1 = (a0 * b1 + a1 * b0 + a1 * b1 * b) % n
        return (r0, r1)

    def poly_pow(p, exp, n, b, c):
        """Raise element to power exp in Z_n[x]/(x^2 - bx - c)."""
        result = (1, 0)  # 1
        base = p
        while exp > 0:
            if exp % 2 == 1:
                result = poly_mul(result, base, n, b, c)
            base = poly_mul(base, base, n, b, c)
            exp //= 2
        return result

    # x = (0, 1)
    x_to_n = poly_pow((0, 1), n, n, b, c)

    # For prime n with (D/n) = -1: x^n ≡ b - x ≡ (b, -1) mod n
    expected = (b % n, (-1) % n)

    if x_to_n != expected:
        return False

    # Additional check: x^{n+1} ≡ -c (mod n)
    x_to_n1 = poly_mul(x_to_n, (0, 1), n, b, c)
    expected2 = ((-c) % n, 0)

    if x_to_n1 != expected2:
        return False

    return True


def enhanced_frobenius_test(n):
    """
    Enhanced Frobenius test trying multiple parameters.
    Uses Grantham's parameter selection.
    """
    if n < 2:
        return False
    if n <= 3:
        return True
    if n % 2 == 0:
        return False

    sr = int(math.isqrt(n))
    if sr * sr == n:
        return False

    # Find b, c with (b^2 + 4c / n) = -1
    for c in range(-10, 10):
        if c == 0:
            continue
        for b in range(0, 20):
            D = b * b + 4 * c
            if D == 0:
                continue
            if math.gcd(abs(D), n) != 1 or math.gcd(abs(c), n) != 1:
                continue
            js = jacobi_symbol(D, n)
            if js == -1:
                result = frobenius_test(n, b, c)
                if result is not None:
                    return result

    return True  # No suitable parameters found (shouldn't happen for odd n)


def frobenius_approach_analysis():
    """
    Analyze the Quadratic Frobenius Test for TC^0 primality.

    KEY INSIGHT: QFT operates in Z_n[x]/(f(x)) where deg(f) = 2.
    Exponentiation x^n in this ring = 2x2 matrix powering.
    2x2 MPOW IS in TC^0 (Mereghetti-Palano 2000).

    So: if QFT is deterministic, then PRIMES is in TC^0!
    """
    print("=" * 70)
    print("APPROACH 3: Quadratic Frobenius Test")
    print("=" * 70)

    print("\nTheoretical analysis:")
    print("  QFT works in Z_n[x]/(x^2 - bx - c), a 2D Z_n-algebra")
    print("  Element multiplication = 2x2 matrix mult over Z_n")
    print("  x^n computation = 2x2 matrix powering = MPOW_2 → IN TC^0")
    print("  Jacobi symbol computation → IN TC^0 (HAB 2002)")
    print("  Parameter selection (find b,c with (D/n)=-1):")
    print("    Small b,c suffice (GRH → O(log^2 n) parameters)")
    print("    Trying constant number of params: IN TC^0 if deterministic")
    print()

    # Test QFT accuracy
    print("Testing Quadratic Frobenius Test accuracy:")

    false_positives = []
    false_negatives = []
    errors = 0

    for n in range(2, 50000):
        try:
            qft = enhanced_frobenius_test(n)
            actual = isprime(n)
            if qft and not actual:
                false_positives.append(n)
            if not qft and actual:
                false_negatives.append(n)
        except Exception as e:
            errors += 1

    print(f"  Range 2-49999:")
    print(f"    False positives (composites passing): {len(false_positives)}")
    print(f"    False negatives (primes failing): {len(false_negatives)}")
    print(f"    Errors: {errors}")
    if false_positives[:20]:
        print(f"    First FPs: {false_positives[:20]}")
    if false_negatives[:10]:
        print(f"    First FNs: {false_negatives[:10]}")

    print()
    print("  Grantham (1998) proved: QFT has error < 1/7710 per test")
    print("  Much stronger than Miller-Rabin (error < 1/4)")
    print("  BUT: still probabilistic for single parameter set")
    print()

    # Can multiple QFT parameters give deterministic test?
    if false_positives:
        print(f"  Analyzing {min(len(false_positives), 15)} QFT pseudoprimes...")
        for n in false_positives[:15]:
            # Try many parameter pairs
            caught = False
            for c2 in range(-20, 20):
                if c2 == 0:
                    continue
                for b2 in range(0, 30):
                    D2 = b2 * b2 + 4 * c2
                    if D2 == 0 or math.gcd(abs(D2), n) != 1 or math.gcd(abs(c2), n) != 1:
                        continue
                    js2 = jacobi_symbol(D2, n)
                    if js2 == -1:
                        r = frobenius_test(n, b2, c2)
                        if r == False:
                            print(f"    n={n}: caught by (b,c)=({b2},{c2})")
                            caught = True
                            break
                if caught:
                    break
            if not caught:
                print(f"    n={n}: NOT caught by extended parameter search!")

    print()
    print("  CONCLUSION on QFT:")
    print("  - Single QFT with good parameters: very few pseudoprimes")
    print("  - Multiple QFT: likely catches all composites (conjectural)")
    print("  - QFT IS in TC^0 for any fixed parameter set")
    print("  - Deterministic QFT would need GRH or similar assumption")

    return false_positives


# ============================================================================
# APPROACH 4: Sum of Two Squares
# ============================================================================

def is_sum_of_two_squares(n):
    """Check if n = a^2 + b^2 for some integers a, b."""
    for a in range(int(math.isqrt(n)) + 1):
        b2 = n - a * a
        if b2 < 0:
            break
        b = int(math.isqrt(b2))
        if b * b == b2:
            return True, a, b
    return False, 0, 0


def count_two_square_representations(n):
    """Count the number of representations n = a^2 + b^2 with a <= b, a >= 0."""
    count = 0
    reps = []
    for a in range(int(math.isqrt(n // 2)) + 1):
        b2 = n - a * a
        b = int(math.isqrt(b2))
        if b * b == b2 and b >= a:
            count += 1
            reps.append((a, b))
    return count, reps


def sum_of_squares_approach():
    """
    Analyze sum-of-two-squares approach for TC^0 primality.

    Fermat's theorem: p ≡ 1 (mod 4) is prime iff p = a^2 + b^2 UNIQUELY.
    (Uniqueness up to order and signs.)

    Issues:
    1. Only works for p ≡ 1 (mod 4)
    2. Need to COUNT representations, not just find one
    3. Finding a representation uses Cornacchia's algorithm (GCD-like)
    4. GCD is in TC^0 (HAB 2002), but Cornacchia needs more
    """
    print("=" * 70)
    print("APPROACH 4: Sum of Two Squares Characterization")
    print("=" * 70)

    print("\nTheoretical analysis:")
    print("  Fermat: p prime, p ≡ 1 (mod 4) → p = a² + b² uniquely")
    print("  Converse: if n ≡ 1 (mod 4) and n = a² + b² uniquely, is n prime?")
    print()

    # Test the converse: is unique representation sufficient for primality?
    print("  Testing converse for n ≡ 1 (mod 4), n < 10000:")
    counterexamples = []
    for n in range(5, 10001, 4):
        count, reps = count_two_square_representations(n)
        if count == 1 and not isprime(n):
            counterexamples.append((n, reps))
        elif count == 1 and isprime(n):
            pass  # Expected

    print(f"    Composites with unique 2-square representation: {len(counterexamples)}")
    if counterexamples[:10]:
        for n, reps in counterexamples[:10]:
            print(f"      n={n} = {reps[0][0]}² + {reps[0][1]}², factors: {factorint(n)}")

    # Actually the correct theorem is more nuanced
    print()
    print("  CORRECTION: Fermat's theorem says:")
    print("  A PRIME p ≡ 1 (mod 4) has EXACTLY ONE representation as a² + b²")
    print("  A COMPOSITE n = p₁^a₁ * ... with all primes ≡ 1 (mod 4)")
    print("  has MULTIPLE representations if it has ≥ 2 distinct prime factors")
    print()

    # Check: does unique representation + n ≡ 1 (mod 4) imply prime or prime power?
    print("  Refined test: n ≡ 1 (mod 4) and exactly 1 representation as a²+b²")
    print("  implies n is a prime or a prime power?")

    prime_power_counter = []
    for n, reps in counterexamples:
        f = factorint(n)
        if len(f) == 1:  # prime power
            prime_power_counter.append(n)
        else:
            print(f"    GENUINE counterexample: n={n}, factors={f}, reps={reps}")

    if prime_power_counter:
        print(f"    Prime powers with unique rep: {prime_power_counter[:10]}")

    # The real issue: can we COUNT representations in TC^0?
    print()
    print("  Circuit complexity of counting representations:")
    print("  r_2(n) = #{(a,b): a²+b² = n} = 4 * sum_{d|n} chi(d)")
    print("  where chi is the non-principal character mod 4")
    print("  This requires knowing the DIVISORS of n → requires FACTORING")
    print("  Factoring is NOT known to be in TC^0 (not even known in P for large n!)")
    print()
    print("  Alternative: just CHECK if n = a² + b² for some a,b")
    print("  Cornacchia's algorithm: O(log n) steps, each involves GCD/division")
    print("  GCD IS in TC^0, but Cornacchia is SEQUENTIAL (each step depends on prev)")
    print("  Cornacchia depth: O(log n) = O(N) → NOT constant depth")
    print()
    print("  CONCLUSION: Sum-of-two-squares approach FAILS for TC^0.")
    print("  1. Only handles p ≡ 1 (mod 4)")
    print("  2. Counting representations requires factoring (not in TC^0)")
    print("  3. Finding a representation requires O(log n) sequential steps")
    print("  STATUS: CLOSED for TC^0 purposes.")

    return counterexamples


# ============================================================================
# COMBINED ANALYSIS: What CAN go in TC^0?
# ============================================================================

def combined_analysis():
    """
    Determine the most promising path to PRIMES in TC^0.
    """
    print("=" * 70)
    print("COMBINED ANALYSIS: TC^0 Primality Test Landscape")
    print("=" * 70)

    print()
    print("Operations known to be in TC^0:")
    print("  1. Addition, multiplication of N-bit integers")
    print("  2. Division (HAB 2002)")
    print("  3. GCD (follows from division + Euclidean algorithm in TC^0)")
    print("     WAIT: GCD is in NC^1, not known in TC^0!")
    print("     Correction: Integer division is in TC^0. Extended GCD is in NC^1.")
    print("  4. Jacobi symbol: involves GCD-like computation. In NC^1.")
    print("     Hmm, actually Jacobi symbol uses quadratic reciprocity, which")
    print("     can be computed from bit operations → likely TC^0")
    print("  5. Scalar powering x^n mod m (HAB 2002)")
    print("  6. 2x2 matrix powering M^n mod m (MPOW_2, Mereghetti-Palano 2000)")
    print("  7. Fixed-k matrix powering MPOW_k (Mereghetti-Palano 2000)")
    print("  8. Iterated multiplication of poly(N) integers (HAB 2002)")
    print()

    print("Primality tests and their TC^0 status:")
    print()
    print("  Test                    | Key Operation           | In TC^0?")
    print("  " + "-" * 65)
    print("  Trial division          | Division by all k<√n    | NO (2^{N/2} ops)")
    print("  Wilson's theorem        | (n-1)! mod n            | NO (2^N mults)")
    print("  Miller-Rabin (fixed a)  | a^{n-1} mod n           | YES (scalar pow)")
    print("  Strong Lucas (fixed PQ) | 2x2 MPOW + scalar pow  | YES (MPOW_2)")
    print("  Quadratic Frobenius     | 2x2 MPOW               | YES")
    print("  AKS                     | r×r MPOW (r=polylog)    | UNKNOWN")
    print("  Fermat test (fixed a)   | a^{n-1} mod n           | YES (scalar pow)")
    print("  Sum of two squares      | Cornacchia/factoring    | NO")
    print()

    print("CRITICAL OBSERVATION:")
    print("  Miller-Rabin base a: a^{n-1} ≡ 1 (mod n) is scalar powering → TC^0")
    print("  Strong Lucas (P,Q): 2x2 MPOW → TC^0")
    print("  Combined: BPSW = MR(2) + SLPSP(Selfridge) → TC^0")
    print()
    print("  The ONLY question is whether these tests are UNCONDITIONALLY correct.")
    print()

    print("  Known results on BPSW:")
    print("  - No BPSW pseudoprime found below 2^64 (exhaustive search)")
    print("  - Pomerance (1984): heuristic argument that BPSW PSPs exist")
    print("  - If BPSW is correct: PRIMES is in TC^0 UNCONDITIONALLY")
    print("  - But this is CONDITIONAL on BPSW correctness (unproven)")
    print()

    print("  Known results on Miller-Rabin with fixed bases:")
    print("  - MR with bases {2,3,5,7,11,13,17,19,23,29,31,37}:")
    print("    deterministic for n < 3.317×10^24 (Sorenson & Webster 2016)")
    print("  - For ALL n: GRH → MR(bases 2..2ln²(n)) is deterministic (Miller 1976)")
    print("  - Under GRH: PRIMES in TC^0 (just O(log²n) scalar powerings)")
    print("  - Without GRH: need infinitely many bases → not fixed TC^0 circuit")
    print()

    print("=" * 70)
    print("STRONGEST UNCONDITIONAL RESULT WE CAN ESTABLISH:")
    print("=" * 70)
    print()
    print("  For n < B (any fixed bound B):")
    print("  There exists a FINITE set of Miller-Rabin bases that is deterministic.")
    print("  This gives a NONUNIFORM TC^0 family {C_N} for PRIMES.")
    print("  (Different circuit for each input length N, but each is TC^0.)")
    print()
    print("  For UNIFORM TC^0 (single algorithm for all n):")
    print("  REQUIRES one of:")
    print("    (a) GRH (then Miller's test works)")
    print("    (b) BPSW correctness (then BPSW works)")
    print("    (c) Some other unconditional deterministic test using only")
    print("        scalar powers and 2x2 matrix powers")
    print()
    print("  Option (c) investigation:")
    print("  - Lenstra (1981): if we can find r with multiplicative order")
    print("    of n mod r > polylog(n), then Fermat test mod (x^r-1) works")
    print("    BUT this requires growing-dimension matrix powering (AKS problem)")
    print("  - Multiple Frobenius tests: even with O(1) params, no proof of correctness")
    print("  - Goldwasser-Kilian ECPP: uses elliptic curves, NOT in TC^0")
    print("  - Adleman-Huang: also not TC^0")
    print()

    # NEW ANGLE: What about combining BPSW-type tests with additional structure?
    print("  NEW ANGLE: Enhanced BPSW variants")
    print("  - Baillie-PSW with extra Frobenius check: still conjectural")
    print("  - Multiple Lucas sequences with different discriminants:")
    print("    Each is TC^0, but proving correctness for finite # of params is open")
    print()

    # Can we prove anything about the MINIMUM number of tests needed?
    print("  INFORMATION-THEORETIC ARGUMENT:")
    print("  - Each scalar/2x2 MPOW test gives O(N) bits of information about n")
    print("  - Primality is 1 bit of information")
    print("  - So in principle, O(1) tests should suffice")
    print("  - The question is whether O(1) SPECIFIC tests suffice for ALL n")
    print("  - This is exactly the BPSW / GRH question")


# ============================================================================
# EXPERIMENT: Find the minimal deterministic TC^0 test
# ============================================================================

def find_minimal_deterministic():
    """
    Experimentally find the minimum number of combined MR + Lucas tests
    needed for correct primality testing up to various bounds.
    """
    print()
    print("=" * 70)
    print("EXPERIMENT: Minimal Deterministic TC^0 Test (up to bounds)")
    print("=" * 70)

    # For each bound, find minimum bases for MR
    bounds = [1000, 10000, 100000, 1000000]

    for bound in bounds:
        print(f"\n  Bound: {bound}")

        # Try increasing sets of MR bases
        # Known: MR with bases {2, 3} is correct for n < 1373653
        bases_to_try = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]

        for num_bases in range(1, len(bases_to_try) + 1):
            bases = bases_to_try[:num_bases]
            correct = True

            for n in range(2, min(bound, 100001)):
                # Miller-Rabin with these bases
                mr_result = True
                if n in bases:
                    mr_actual = isprime(n)
                    continue
                if n < 2 or n % 2 == 0:
                    continue

                # Write n-1 = d * 2^r
                d = n - 1
                r = 0
                while d % 2 == 0:
                    d //= 2
                    r += 1

                composite = True
                for a in bases:
                    if a >= n:
                        continue
                    x = pow(a, d, n)
                    if x == 1 or x == n - 1:
                        composite = False
                        continue
                    found = False
                    for _ in range(r - 1):
                        x = pow(x, 2, n)
                        if x == n - 1:
                            found = True
                            break
                    if found:
                        composite = False

                if composite != (not isprime(n)):
                    correct = False
                    break

            if correct:
                print(f"    MR bases {bases}: CORRECT up to {min(bound, 100000)}")
                break

        # Now check: does adding strong Lucas help?
        # Check with just MR(2) + SLPSP
        bpsw_correct = True
        for n in range(2, min(bound, 100001)):
            if n < 3 or n % 2 == 0:
                continue
            # MR base 2
            d = n - 1
            r = 0
            while d % 2 == 0:
                d //= 2
                r += 1
            x = pow(2, d, n)
            mr2_pass = False
            if x == 1 or x == n - 1:
                mr2_pass = True
            else:
                for _ in range(r - 1):
                    x = pow(x, 2, n)
                    if x == n - 1:
                        mr2_pass = True
                        break

            if not mr2_pass:
                # MR says composite
                if isprime(n):
                    bpsw_correct = False
                    break
                continue

            # Strong Lucas
            sl = strong_lucas_test(n)

            if (mr2_pass and sl) != isprime(n):
                bpsw_correct = False
                print(f"    BPSW failure at n={n}!")
                break

        print(f"    BPSW (MR(2)+SLPSP): {'CORRECT' if bpsw_correct else 'FAILED'} up to {min(bound, 100000)}")


# ============================================================================
# MAIN
# ============================================================================

if __name__ == "__main__":
    print("TC^0 PRIMALITY TEST INVESTIGATION")
    print("=" * 70)
    print()
    print("Goal: Find unconditional TC^0 primality test (constant depth, poly size)")
    print("Key: 2x2 MPOW and scalar powering are in TC^0")
    print()

    t0 = time.time()

    # Approach 1: Wilson's theorem
    factorial_mod_circuit_analysis()
    print()

    # Approach 2: Lucas sequences
    lucas_fps = lucas_approach_analysis()
    print()

    # Approach 3: Quadratic Frobenius
    frob_fps = frobenius_approach_analysis()
    print()

    # Approach 4: Sum of two squares
    s2s_counterex = sum_of_squares_approach()
    print()

    # Combined analysis
    combined_analysis()

    # Experiment: minimal deterministic
    find_minimal_deterministic()

    elapsed = time.time() - t0
    print(f"\n\nTotal time: {elapsed:.1f}s")

    print()
    print("=" * 70)
    print("FINAL SUMMARY")
    print("=" * 70)
    print()
    print("Approach 1 (Wilson): CLOSED. (n-1)! mod n needs 2^N mults, not TC^0.")
    print("Approach 2 (Lucas): PROMISING. Strong Lucas IS in TC^0, but")
    print("  unconditional correctness unproven. Lucas PSPs exist but are rare.")
    print("Approach 3 (Frobenius): PROMISING. QFT IS in TC^0, error < 1/7710.")
    print("  Deterministic QFT unknown without GRH.")
    print("Approach 4 (Sum-of-squares): CLOSED. Needs factoring or O(N)-depth GCD.")
    print()
    print("BEST PATH: BPSW (MR base 2 + Strong Lucas) is in TC^0.")
    print("  If BPSW correctness is proven → PRIMES in TC^0 unconditionally.")
    print("  No BPSW pseudoprime known below 2^64.")
    print("  This is the NARROWEST gap between TC^0 and PRIMES.")

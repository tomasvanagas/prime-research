"""
Session 5: Arithmetic Term Construction for p(n)

BACKGROUND: Prunescu & Shunia (2024) showed a FIXED-LENGTH arithmetic term
for p(n) exists using only +, -, *, //, ^. But their formula involves
superexponential intermediate values.

OUR APPROACH: Try to construct a MORE EFFICIENT arithmetic term by:

1. Using the Willans-type construction with optimized primality testing
2. Using floor functions and modular arithmetic to encode sieve logic
3. Trying to encode Lucy DP as a fixed-depth arithmetic expression

KEY INSIGHT: The issue with Willans/Wilson-based formulas is that
(n-1)! grows as n^n. But what if we use MODULAR arithmetic to keep
numbers small? If we compute (n-1)! mod n, we only need numbers up to n^2.

APPROACH A: Optimized Willans formula
  p(n) = 1 + sum_{k=1}^{2n*ln(n)} [pi(k) < n]
  where pi(k) = sum_{j=2}^{k} [j is prime]
  and [j is prime] = [(j-1)! mod j == j-1]  (Wilson's theorem)

  Total operations: O(n * ln(n)) multiplications, each involving numbers up to k.
  For p(10^100): k ~ 10^102, so each (k-1)! mod k needs ~10^102 multiplications
  of numbers up to 10^102. Total: ~10^204 multiplications. INFEASIBLE.

APPROACH B: Miller-Rabin as arithmetic term
  MR test for n with witness a:
  Write n-1 = 2^s * d
  Compute a^d mod n, then square s times
  If a^d ≡ 1 (mod n) or a^{2^r * d} ≡ -1 (mod n) for some r, PROBABLY prime

  For deterministic MR: witnesses {2,3,5,7,11,13} suffice for n < 3.3×10^24
  For n < 10^{102}: need ~O(ln(n)^2) witnesses under GRH

  Can this be expressed as a fixed-length arithmetic term?
  YES! Each witness test is: pow(a, d, n) which is O(log(n)) multiplications mod n.

  But encoding "find the n-th number passing MR" requires a loop/sum.

APPROACH C: Novel — R^{-1}(n) as arithmetic term + rounding
  R^{-1}(n) ≈ n * ln(n) * (1 + correction terms)
  If we could express this with floor/ceiling and get error < gap/2,
  then p(n) = nearest_prime(R^{-1}(n)).

  But error ~ √p(n) >> gap ~ ln(p(n)).

APPROACH D: Encode the sieve of Eratosthenes as matrix operations
  The sieve marks composites. Can matrix power encode sieve steps?
  A^n gives the state after n sieve steps.
  If A is small and fixed, A^n is computable in O(log n) multiplications.

Let's test what's actually achievable.
"""

import time
import sympy
from sympy import prime as sympy_prime, isprime, factorial, primepi
from mpmath import mp, mpf, log, exp, floor as mpfloor, ceil as mpceil, li
from mpmath import nstr, fabs
import numpy as np

mp.dps = 50


def willans_formula(n):
    """
    Willans' formula (1964):
    p(n) = 1 + sum_{k=1}^{2^n} floor(n / (1 + pi(k)))^{1/n}

    This is exact but requires computing pi(k) for ALL k up to 2^n.
    Completely impractical for large n.

    Simplified version using Wilson's theorem:
    pi(k) = sum_{j=2}^{k} floor(cos^2(pi * (j-1)!/j))  [Willans]

    Let's just test correctness for tiny n.
    """
    count = 0
    for k in range(2, 1000):
        if isprime(k):
            count += 1
        if count == n:
            return k
    return -1


def wilson_primality(k):
    """Test if k is prime using Wilson's theorem: (k-1)! ≡ -1 (mod k)"""
    if k < 2:
        return False
    if k == 2:
        return True
    fact_mod = 1
    for j in range(2, k):
        fact_mod = (fact_mod * j) % k
    return fact_mod == k - 1


def miller_rabin_det(n, witnesses=None):
    """
    Deterministic Miller-Rabin.
    For n < 3,317,044,064,679,887,385,961,981:
    witnesses = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]
    """
    if n < 2:
        return False
    if n < 4:
        return True
    if n % 2 == 0:
        return False

    if witnesses is None:
        witnesses = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]

    d = n - 1
    s = 0
    while d % 2 == 0:
        d //= 2
        s += 1

    for a in witnesses:
        if a >= n:
            continue
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(s - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                break
        else:
            return False
    return True


def approach_d_sieve_matrix():
    """
    APPROACH D: Can we encode the sieve as a matrix?

    Consider the sequence s(k) = 1 if k is not sieved, 0 otherwise.
    After sieving by prime p: s(k) = 0 for k divisible by p.

    This is equivalent to: s(k) *= (1 - floor(k/p) + floor((k-1)/p))
    which is 0 when p|k and 1 otherwise.

    So: "k is prime" iff k >= 2 and product_{p <= sqrt(k)} (1 - [p|k]) = 1
    where [p|k] = 1 iff p divides k.

    [p|k] = floor(k/p) - floor((k-1)/p) = k//p - (k-1)//p

    This IS an arithmetic term! But it requires knowing primes up to sqrt(k)
    to sieve — CIRCULAR.

    UNLESS we sieve with ALL integers 2..sqrt(k), not just primes.
    Product_{d=2}^{sqrt(k)} (1 - [d|k]) = 0 for ALL composite k (since
    composites have at least one factor d in [2, sqrt(k)]).
    But it's also 0 for primes p when we test d=p! So this doesn't work.

    Correct test: k is prime iff the smallest divisor > 1 of k is k itself.
    min_divisor(k) = min{d >= 2 : d | k}
    k is prime iff min_divisor(k) = k

    Can min_divisor be an arithmetic term?
    min_divisor(k) = sum_{d=2}^{k} d * product_{j=2}^{d-1} ceil((k mod j) / k)
    (each term is d times "no j<d divides k")

    But this is O(k) operations per test, and we need O(n*ln(n)) tests.
    """
    print("=" * 80)
    print("APPROACH D: SIEVE AS ARITHMETIC OPERATIONS")
    print("=" * 80)

    # Test the arithmetic primality test
    def is_prime_arithmetic(k):
        """Purely arithmetic primality test using only +,-,*,//,mod"""
        if k < 2:
            return 0
        for d in range(2, int(k**0.5) + 1):
            if k % d == 0:
                return 0
        return 1

    # This is trial division — O(sqrt(k)) per test
    # To find p(n): need to test ~n*ln(n) candidates
    # Total: O(n * ln(n) * sqrt(n*ln(n))) = O(n^{3/2} * ln^{3/2}(n))

    # For p(10^100): O(10^{150} * 230^{3/2}) ≈ O(10^{154})
    # COMPLETELY INFEASIBLE.

    print("\nTrial division as arithmetic term:")
    t0 = time.time()
    count = 0
    k = 2
    target = 10000
    while count < target:
        if is_prime_arithmetic(k):
            count += 1
            if count == target:
                result = k
        k += 1
    t1 = time.time()
    print(f"  p({target}) = {result} (actual: {sympy_prime(target)})")
    print(f"  Time: {t1-t0:.3f}s")
    print(f"  Correct: {result == sympy_prime(target)}")


def approach_e_floor_sum():
    """
    APPROACH E: Floor function identities for prime counting.

    There exist identities like:
    pi(x) = sum_{k=2}^{x} floor(2/k * sum_{j=1}^{k-1} cos(2*pi*j/k))

    But the cosine makes this non-arithmetic.

    Pure arithmetic version:
    pi(x) = -1 + sum_{k=2}^{x} (1 - floor(sum_{d=2}^{sqrt(k)} (1 - ceil(k%d / k)) / (k-1)))

    This is just counting primes by trial division in disguise.

    Better: use Legendre's identity
    pi(x) = pi(sqrt(x)) + x - 1 - sum_{p<=sqrt(x)} (floor(x/p) - pi(p) + 1) + ...

    But Legendre/Meissel is already what Lucy DP implements!
    The question is: can it be expressed as a FIXED-DEPTH arithmetic expression?

    NO — because the recursion depth depends on x. Specifically:
    - Meissel: depth depends on pi(x^{1/3}) ~ x^{1/3}/ln(x)
    - Lucy DP: needs O(sqrt(x)) distinct values
    - Both require VARIABLE-depth computation

    A FIXED-DEPTH expression cannot adapt to the input size.
    Unless... we use exponentiation to encode variable-depth loops.

    KEY: A^n (matrix power) encodes n iterations of a linear map.
    Can we encode the sieve/DP as a linear map?
    """
    print("\n" + "=" * 80)
    print("APPROACH E: FIXED-DEPTH ARITHMETIC EXPRESSIONS")
    print("=" * 80)

    # Test: can pi(x) be computed by a fixed-depth expression with large intermediates?

    # Willans-type: pi(n) = sum_{k=2}^{n} floor(((k-1)! mod k) / (k-1))
    # This works because (k-1)! mod k = k-1 if k is prime, 0 otherwise
    def pi_willans(n):
        count = 0
        for k in range(2, n + 1):
            fact_mod = 1
            for j in range(2, k):
                fact_mod = (fact_mod * j) % k
            if fact_mod == k - 1:
                count += 1
        return count

    # Test
    for x in [10, 20, 50, 100]:
        t0 = time.time()
        result = pi_willans(x)
        t1 = time.time()
        print(f"  pi({x}) = {result} (actual: {primepi(x)}, time: {t1-t0:.4f}s)")


def approach_f_encoding_trick():
    """
    APPROACH F: Encode p(n) in a single large number.

    IDEA: The "prime constant" C = sum_{k=1}^{inf} 2^{-p(k)} encodes ALL primes.
    Binary: 0.0110101000101000101000100...
    The 1-bits are at positions p(1)=2, p(2)=3, p(3)=5, etc.

    If we KNEW C to sufficient precision, we could extract p(n) by finding
    the n-th 1-bit. But C is not computable in finite time to arbitrary precision.

    Alternative: Use a TRUNCATED version.
    C_N = sum_{k=1}^{N} 2^{-p(k)} (known for small N)

    This encodes the first N primes. But computing C_N requires knowing
    the first N primes — circular.

    HOWEVER: What if we could compute C_N from a FORMULA?

    C = sum_p 2^{-p} = sum_{n=2}^{inf} [n is prime] * 2^{-n}
    = sum_{n=2}^{inf} floor(((n-1)! mod n) / (n-1)) * 2^{-n}  [Wilson]

    This is an exact expression but requires computing for ALL n.

    BETTER IDEA: Mills' constant approach.
    Mills (1947): There exists A > 1 such that floor(A^{3^n}) is prime for all n.
    A ≈ 1.30637788386308069046...

    But computing A requires knowing the primes — circular again.

    UNLESS: A has a fast-converging series that doesn't require knowing primes.
    Does it? Let's check.
    """
    print("\n" + "=" * 80)
    print("APPROACH F: ENCODING TRICKS")
    print("=" * 80)

    # Mills' constant
    # A^3 ≈ 2.229..., A^9 ≈ 11.077..., A^27 ≈ 1361.000..., etc.
    # floor(A^3) = 2, floor(A^9) = 11, floor(A^27) = 1361
    # These are primes!

    # But this only gives a sparse subsequence of primes (towers of 3)
    # p(1), p(?), p(??), ... growing super-exponentially

    # What about the "prime-generating polynomial" approach?
    # Euler: n^2 + n + 41 is prime for n = 0, 1, ..., 39
    # But fails at n = 40 (40^2 + 40 + 41 = 41^2)

    # No single polynomial can generate ALL primes (proved by Goldbach)

    # HOWEVER: Can we use DIFFERENT polynomials for different ranges of n?
    # p(n) = f_1(n) for n in range 1
    # p(n) = f_2(n) for n in range 2
    # etc.

    # This requires knowing the breakpoints — which requires knowing the primes.

    # Let's test how well simple polynomials fit primes locally
    print("\nLocal polynomial fitting of p(n):")

    for start_n in [100, 1000, 10000]:
        ns = list(range(start_n, start_n + 50))
        ps = [sympy_prime(n) for n in ns]

        # Fit polynomial
        for deg in [1, 2, 3, 5]:
            coeffs = np.polyfit(ns, ps, deg)
            pred = np.polyval(coeffs, ns)
            errors = np.array(ps) - pred
            exact = sum(1 for e in errors if abs(e) < 0.5)

            # Test extrapolation
            test_ns = list(range(start_n + 50, start_n + 60))
            test_ps = [sympy_prime(n) for n in test_ns]
            test_pred = np.polyval(coeffs, test_ns)
            test_errors = np.array(test_ps) - test_pred
            test_exact = sum(1 for e in test_errors if abs(e) < 0.5)

            if deg <= 3 or exact > 0:
                print(f"  n~{start_n}, deg={deg}: {exact}/50 train exact, "
                      f"{test_exact}/10 test exact, max_err={max(abs(errors)):.1f}")


def approach_g_recursive_halving():
    """
    APPROACH G: Recursive doubling for p(n).

    Idea: p(2n) ≈ 2 * p(n) * (ln(2n)/ln(n)) ≈ 2 * p(n) * (1 + ln(2)/ln(n))

    More precisely: p(2n) = p(n) + sum_{k=n+1}^{2n} g(k) where g(k) = gap

    If we could estimate this sum of gaps, we could compute p(2n) from p(n).
    Then p(n) = T(n, p(n)) where T builds up from known base cases.

    This is a divide-and-conquer: p(n) from p(n/2), which from p(n/4), etc.
    Depth: O(log n). At each level, we need the sum of ~n/2^k gaps.

    Average sum of n gaps = p(n) - 2 ≈ n*ln(n)
    Variance of sum ≈ n * Var(gap) ≈ n * (ln(n))^2
    Std ≈ sqrt(n) * ln(n)

    So relative error ≈ sqrt(n) * ln(n) / (n * ln(n)) = 1/sqrt(n)

    For n = 10^100: relative error ≈ 10^{-50}
    Absolute error ≈ 10^{-50} * 10^{102} = 10^{52}
    Gap ≈ ln(10^{102}) ≈ 235
    Error/gap ≈ 10^{52}/235 ≈ 10^{50}

    STILL TOO LARGE. The variance of gap sums is unavoidable.
    """
    print("\n" + "=" * 80)
    print("APPROACH G: RECURSIVE HALVING")
    print("=" * 80)

    # Test: compute p(2n) from p(n) + average gap correction
    errors = []
    for n in [100, 200, 500, 1000, 2000, 5000, 10000]:
        pn = sympy_prime(n)
        p2n = sympy_prime(2*n)

        # Estimate p(2n) from p(n)
        # Average gap around p(n) is ln(p(n))
        # Number of primes from p(n) to p(2n) is n
        # So p(2n) ≈ p(n) + n * ln(p(2n))
        # Better: p(2n) ≈ p(n) + integral from n to 2n of ln(p(t)) dt
        # ≈ p(n) + n * (ln(2n) + ln(ln(2n)))

        from mpmath import mpf, log
        est = pn + int(float(mpf(n) * (log(mpf(2*n)) + log(log(mpf(2*n))))))
        error = est - p2n
        rel_error = abs(error) / p2n * 100
        gap = sympy_prime(2*n + 1) - p2n

        print(f"  n={n}: p(n)={pn}, p(2n)={p2n}")
        print(f"    Estimate: {est}, error: {error:+d}, |error|/gap: {abs(error)/gap:.1f}")
        errors.append(abs(error) / p2n)

    print(f"\n  Relative errors: {[f'{e:.6f}' for e in errors]}")
    print(f"  Trend: {'decreasing' if errors[-1] < errors[0] else 'not decreasing'}")


def approach_h_number_theory_identities():
    """
    APPROACH H: Use number theory identities to build p(n) from simpler quantities.

    Identity 1: sum_{p<=x} 1 = pi(x)
    Identity 2: sum_{p<=x} ln(p) = theta(x) ~ x (PNT)
    Identity 3: sum_{p<=x} 1/p = ln(ln(x)) + M (Mertens)
    Identity 4: prod_{p<=x} (1-1/p) = e^{-gamma}/ln(x) * (1 + O(1/ln(x))) (Mertens)

    Can we INVERT any of these to get p(n)?

    From identity 2: theta(p(n)) = sum_{k=1}^{n} ln(p(k))
    So p(n) = exp(theta(p(n)) - theta(p(n-1))) ... circular.

    From identity 3: sum_{k=1}^{n} 1/p(k) = ln(ln(p(n))) + M + o(1)
    So if we knew sum_{k=1}^{n} 1/p(k), we could estimate p(n).
    But computing this sum requires knowing the primes.

    NOVEL IDEA: What about INVERSE problems?
    Given that theta(x) = sum_{p<=x} ln(p), we know theta(p(n)) very precisely.
    theta(p(n)) = p(n) - correction (by PNT, theta(x) ~ x)

    More precisely: theta(p(n)) = p(n) + O(p(n)/ln(p(n)))
    Hmm, this is weaker than what we already have.

    Actually: theta(x) = x - sum_rho x^rho/rho - ln(2pi) - (1/2)ln(1-x^{-2})
    (explicit formula for theta)

    This has the SAME zeta zero sum as pi(x). No improvement.
    """
    print("\n" + "=" * 80)
    print("APPROACH H: NUMBER THEORY IDENTITY INVERSIONS")
    print("=" * 80)

    # Test: can we use theta(x) to improve on R^{-1}?
    # theta(x) = sum_{p<=x} ln(p)
    # psi(x) = sum_{p^k<=x} ln(p)  (includes prime powers)

    # Compute theta and psi for test values
    for x in [100, 1000, 10000, 100000]:
        primes_up_to_x = list(sympy.primerange(2, x + 1))
        theta = sum(np.log(p) for p in primes_up_to_x)
        psi = 0
        for p in primes_up_to_x:
            pk = p
            while pk <= x:
                psi += np.log(p)
                pk *= p

        pi_x = len(primes_up_to_x)

        print(f"\nx={x}: pi(x)={pi_x}")
        print(f"  theta(x)/x = {theta/x:.6f} (should -> 1)")
        print(f"  psi(x)/x = {psi/x:.6f} (should -> 1)")
        print(f"  (psi-theta)/sqrt(x) = {(psi-theta)/np.sqrt(x):.6f}")

    # Key finding: theta(x)/x -> 1 but the ERROR theta(x) - x is
    # of the same order as pi(x) - li(x), i.e., O(sqrt(x) * ln(x)).
    # Using theta gives NO IMPROVEMENT over using pi and li.


if __name__ == "__main__":
    approach_d_sieve_matrix()
    approach_e_floor_sum()
    approach_f_encoding_trick()
    approach_g_recursive_halving()
    approach_h_number_theory_identities()

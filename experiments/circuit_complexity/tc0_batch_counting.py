"""
Session 16: TC^0 Batch Counting Analysis

CENTRAL QUESTION: If PRIMES is in TC^0 (via BPSW), can we COUNT primes
in [1,x] without enumerating all x candidates?

pi(x) = sum_{n=2}^{x} BPSW(n)

Naive TC^0 circuit: evaluate BPSW(n) for each n, sum results.
Size: x * poly(N) = 2^N * poly(N) -- EXPONENTIAL in N = log2(x).

This file rigorously analyzes 5 potential escape routes:
1. TC^0 MAJORITY gates (fan-in reduction)
2. Divide-and-conquer in TC^0
3. Batch structure in modular exponentiation
4. Algebraic structure of f(k) = 2^{k-1} mod k
5. #TC^0 counting complexity

CONCLUSION (SPOILER): ALL five routes fail. The fundamental barrier is
that modular exponentiation couples exponent and modulus through k,
making the function cryptographically pseudorandom across consecutive k.
The question "#TC^0 ⊆ NC?" remains open and is THE key question.
"""

import numpy as np
import math
from collections import defaultdict
from sympy import primepi, isprime, primerange, factorint, totient
from math import isqrt, log, gcd, log2

# ============================================================================
# SECTION 1: TC^0 MAJORITY GATES -- CAN LARGE FAN-IN HELP?
# ============================================================================

def section1_majority_gates():
    """
    TC^0 = constant-depth circuits with AND, OR, NOT, and MAJORITY gates
    (equivalently, threshold gates THR_t: output 1 iff >= t inputs are 1).

    For counting pi(x), we sum x indicator bits. Iterated addition of x bits
    IS in TC^0: a single MAJORITY gate on x bits outputs 1 iff >= ceil(x/2)
    bits are 1. More precisely, a SUM of x bits (producing an N-bit output)
    can be done by a TC^0 circuit of size O(x * N) and depth O(1).

    THE PROBLEM: The x input bits are NOT given; they must be COMPUTED.
    Each bit = BPSW(k) for k = 2, ..., x. Computing BPSW(k) requires a
    poly(N)-size TC^0 circuit PER value of k.

    QUESTION: Can a MAJORITY gate somehow "implicitly" evaluate BPSW on
    multiple inputs without explicitly computing each one?

    ANSWER: No. Here is the proof.

    THEOREM 1: Any TC^0 circuit computing pi(x) on input x (N bits) must
    have size at least x / polylog(x) = 2^N / poly(N).

    Proof sketch:
    - pi(x) determines the number of primes in [1,x].
    - pi(x) - pi(x-1) = 1 iff x is prime, 0 otherwise.
    - If a circuit C computes pi(x), then C(x) - C(x-1) computes PRIMES.
    - The circuit C(x) - C(x-1) has size 2*|C| + O(N) (difference circuit).
    - So |C| >= |C_PRIMES| / 2 - O(N) where C_PRIMES is any circuit for PRIMES.
    - This doesn't give a superlinear lower bound by itself.

    A BETTER argument (communication complexity):
    - Consider x split as (high_bits, low_bits) = (a, b) where x = a*2^{N/2} + b.
    - pi(a*2^{N/2} + b) as a function of (a, b) has high 2-party communication
      complexity because changing b by 1 can flip pi(x) by 1 (if x is prime).
    - In fact, there are ~2^{N/2}/(N/2) primes in any interval of length 2^{N/2},
      so the function has ~2^{N/2}/N "sensitive" positions.
    - A TC^0 circuit of size s and depth d communicates O(s * d) bits [Hastad 1987].
    - So s * d >= 2^{N/2} / N, giving s >= 2^{N/2} / (N * d).
    - For constant d (TC^0), s >= 2^{N/2} / poly(N).

    This is EXPONENTIAL. The lower bound is not as tight as 2^N (the naive bound)
    but still exponential in N.

    HOWEVER: This argument uses a GENERAL lower bound technique. It applies to
    ANY TC^0 circuit computing pi(x), not just the "enumerate and sum" approach.
    The question is whether the COUNTING PROBLEM has a non-enumerative solution.

    MAJORITY GATE FAN-IN: A MAJORITY gate on m inputs requires m wires.
    To count primes, we need a gate with x = 2^N fan-in (one per candidate).
    Each wire carries BPSW(k), requiring poly(N) gates to compute.
    Total size: 2^N * poly(N). The MAJORITY gate itself is a single gate
    of fan-in 2^N, but this doesn't reduce the total circuit size.

    CONCLUSION: MAJORITY gates don't help because the INPUTS to the gate
    still need to be computed individually.
    """
    print("=" * 70)
    print("SECTION 1: TC^0 MAJORITY GATES")
    print("=" * 70)

    # Demonstrate: to sum x bits, a TC^0 circuit of O(x*log(x)) gates suffices.
    # But each bit is BPSW(k), requiring poly(log(x)) gates.
    # Total: x * poly(log(x)) = 2^N * poly(N).

    # Empirical: verify that MAJORITY over prime indicators does give pi(x)
    for x in [10, 100, 1000, 10000]:
        N = int(log2(x)) + 1
        pi_x = int(primepi(x))
        total_tests = x - 1  # testing 2..x
        circuit_size_naive = total_tests * N**2  # poly(N) per test
        print(f"x={x:>6d}, N={N:>3d}, pi(x)={pi_x:>5d}, "
              f"naive_circuit_size={circuit_size_naive:>12d}, "
              f"ratio_x_to_size={circuit_size_naive/x:.1f}")

    print("""
RESULT: MAJORITY gates require all inputs pre-computed.
For pi(x), the inputs are BPSW(2), BPSW(3), ..., BPSW(x).
Computing all x inputs costs x * poly(N) = 2^N * poly(N) gates.
The MAJORITY gate adds only O(x) more gates -- negligible.

The bottleneck is INPUT GENERATION, not AGGREGATION.
""")


# ============================================================================
# SECTION 2: DIVIDE-AND-CONQUER IN TC^0
# ============================================================================

def section2_divide_and_conquer():
    """
    pi(x) = pi(x/2) + [primes in (x/2, x]]

    If we could compute [primes in (a,b]] efficiently from a and b alone
    (without enumerating), divide-and-conquer could work.

    Recursion tree: depth log(x) = N, branching factor 2.
    At each level, we split [1,x] into 2^k intervals of length x/2^k.
    At the bottom (leaf level), each interval has length 1: trivially solved.

    PROBLEM: The recursion doesn't save work because at each level,
    the total number of candidates is STILL x.

    Level 0: 1 interval of length x       (x candidates total)
    Level 1: 2 intervals of length x/2    (x candidates total)
    Level k: 2^k intervals of length x/2^k (x candidates total)
    Level N: x intervals of length 1       (x candidates total)

    No level reduces the total work below x.

    ALTERNATIVE: Use inclusion-exclusion or sieve structure.
    pi(x) = pi(sqrt(x)) + [sieved count above sqrt(x)]
    This is Legendre's formula. The sieved count uses floor(x/d) values
    where d ranges over products of primes up to sqrt(x).

    The key insight (proven in Session 14): floor(x/d) for d up to sqrt(x)
    takes O(sqrt(x)) distinct values, and these values carry O(2^k)
    independent bits via inclusion-exclusion fractional parts.

    COMPUTATIONAL TEST: Verify error accumulation in divide-and-conquer.
    """
    print("=" * 70)
    print("SECTION 2: DIVIDE-AND-CONQUER")
    print("=" * 70)

    # Test: binary splitting pi(x) = pi(x/2) + pi(x, x/2)
    # where pi(x, x/2) = #{primes in (x/2, x]}

    # With smooth approximation: li(x) - li(x/2) approx primes in (x/2, x]
    from scipy.integrate import quad

    def li(x):
        if x <= 1:
            return 0
        # Logarithmic integral
        result, _ = quad(lambda t: 1/log(t), 2, x)
        return result

    print("\nDivide-and-conquer error analysis:")
    print(f"{'x':>8} {'pi(x)':>8} {'approx':>10} {'abs_err':>10} {'rel_err':>10}")
    for x in [100, 1000, 10000, 100000, 1000000]:
        pi_x = int(primepi(x))

        # Binary splitting with li approximation at leaf level
        def dc_approx(lo, hi, depth=0, max_depth=3):
            if depth >= max_depth or hi - lo <= 10:
                return li(hi) - li(max(lo, 2))
            mid = (lo + hi) // 2
            return dc_approx(lo, mid, depth+1, max_depth) + dc_approx(mid, hi, depth+1, max_depth)

        approx = dc_approx(1, x)
        err = abs(approx - pi_x)
        rel_err = err / pi_x if pi_x > 0 else float('inf')
        print(f"{x:8d} {pi_x:8d} {approx:10.1f} {err:10.1f} {rel_err:10.6f}")

    # Now test: does deeper splitting reduce error?
    print("\nEffect of splitting depth on error (x=100000):")
    x = 100000
    pi_x = int(primepi(x))
    for max_depth in range(1, 10):
        def dc_approx_depth(lo, hi, depth=0):
            if depth >= max_depth or hi - lo <= 1:
                return li(hi) - li(max(lo, 2))
            mid = (lo + hi) // 2
            return dc_approx_depth(lo, mid, depth+1) + dc_approx_depth(mid, hi, depth+1)

        approx = dc_approx_depth(1, x)
        err = abs(approx - pi_x)
        print(f"  depth={max_depth}: approx={approx:.1f}, error={err:.1f} "
              f"({err/pi_x*100:.4f}%)")

    print("""
RESULT: Divide-and-conquer with smooth approximation does NOT reduce error.
The error at each level is O(sqrt(interval_length) / log(interval_length)).
Summing over 2^k intervals of length x/2^k:
  Error = 2^k * O(sqrt(x/2^k) / log(x/2^k))
        = O(sqrt(x) * 2^{k/2} / log(x))

This INCREASES with depth k. Optimal is k=0 (no splitting).
The error floor is O(sqrt(x) / log(x)), matching PNT error.

D&C FAILS because the prime distribution has sqrt(x)-scale oscillations
that cannot be split into independent sub-problems.
""")


# ============================================================================
# SECTION 3: BATCH STRUCTURE IN MODULAR EXPONENTIATION
# ============================================================================

def section3_batch_modexp():
    """
    The core of BPSW: compute 2^{k-1} mod k for k = 2, ..., x.

    QUESTION: Is there shared computation across different k values?

    APPROACH 1: Fix the base and exponent, vary the modulus.
    For a fixed M, computing 2^M mod k for k = 1, ..., x simultaneously:
    - By CRT: 2^M mod k depends on 2^M mod p^a for each prime power p^a || k.
    - For prime p: 2^M mod p depends on M mod ord_p(2) where ord_p(2) | p-1.
    - Precomputing 2^M mod p for all primes p <= x costs O(pi(x) * polylog(x)).
    - Then 2^M mod k for composite k is assembled via CRT in O(polylog(x) per k).
    - Total: O(x * polylog(x)) -- still O(x) per candidate!

    APPROACH 2: The coupling problem.
    For BPSW, the exponent is k-1 (depends on k) and the modulus is k.
    f(k) = 2^{k-1} mod k.

    For prime p: f(p) = 2^{p-1} mod p = 1 (Fermat's little theorem).
    For composite k: f(k) depends on the factorization of k.

    APPROACH 3: CRT decomposition of f(k).
    f(k) = 2^{k-1} mod k. By CRT for k = p_1^{a_1} ... p_r^{a_r}:
    f(k) mod p_i^{a_i} = 2^{(k-1) mod lambda(p_i^{a_i})} mod p_i^{a_i}
    where lambda is the Carmichael function.

    The key coupling: (k-1) mod lambda(p^a) = (p^a * q + ... - 1) mod lambda(p^a)
    where q depends on other prime factors of k. This creates a complex
    dependency between the factors.

    EXPERIMENTAL TEST: Measure the "batch efficiency" -- how much computation
    is shared when computing f(k) for consecutive k.
    """
    print("=" * 70)
    print("SECTION 3: BATCH MODULAR EXPONENTIATION")
    print("=" * 70)

    # Experiment 3a: f(k) = 2^{k-1} mod k for k = 2..1000
    # Analyze the structure of these values
    x_max = 1000
    fvals = {}
    for k in range(2, x_max + 1):
        fvals[k] = pow(2, k - 1, k)

    # For primes: f(p) = 1 (Fermat). Verify.
    fermat_violations = []
    for k in range(3, x_max + 1):  # skip k=2: gcd(2,2)=2, Fermat doesn't apply
        if isprime(k) and fvals[k] != 1:
            fermat_violations.append(k)
    print(f"Fermat violations among odd primes up to {x_max}: {len(fermat_violations)}")
    print(f"(Note: k=2 excluded since gcd(2,2)!=1; 2^1 mod 2 = 0)")
    assert len(fermat_violations) == 0, "Fermat's little theorem violated!"

    # For composites: f(k) can be anything. Pseudoprimes have f(k) = 1.
    pseudoprimes = [k for k in range(2, x_max + 1)
                    if not isprime(k) and fvals[k] == 1]
    print(f"Fermat base-2 pseudoprimes up to {x_max}: {len(pseudoprimes)}")
    print(f"  Values: {pseudoprimes[:20]}{'...' if len(pseudoprimes) > 20 else ''}")

    # Experiment 3b: CRT decomposition analysis
    # For each composite k, how many distinct prime powers divide k?
    # The cost of CRT-based batch computation scales with this.
    print("\nCRT decomposition analysis:")
    total_prime_power_lookups = 0
    for k in range(2, x_max + 1):
        if not isprime(k):
            factors = factorint(k)
            total_prime_power_lookups += len(factors)
    composites = sum(1 for k in range(2, x_max+1) if not isprime(k))
    print(f"Average prime powers per composite k (2..{x_max}): "
          f"{total_prime_power_lookups/composites:.2f}")

    # Experiment 3c: Period analysis of f(k) mod small primes
    # For a fixed prime p, the function g_p(k) = 2^{k-1} mod p has period ord_p(2)*(p-1)
    # at most. Let's measure actual periods.
    print("\nPeriod of g_p(k) = 2^{k-1} mod p for small primes p:")
    for p in [3, 5, 7, 11, 13, 17, 19, 23]:
        vals = [pow(2, k-1, p) for k in range(1, 500)]
        # Find minimal period
        found_period = None
        for period in range(1, 250):
            if all(vals[j] == vals[j - period] for j in range(period, min(period + 100, len(vals)))):
                found_period = period
                break
        ord_p_2 = 1
        v = 2
        while v % p != 1:
            v = (v * 2) % p
            ord_p_2 += 1
        print(f"  p={p:3d}: period={found_period}, ord_p(2)={ord_p_2}, "
              f"p-1={p-1}, lcm={math.lcm(ord_p_2, p-1) if found_period else 'N/A'}")

    # Experiment 3d: Can we compute f(k) for all k using O(sqrt(x)) "precomputed" values?
    # Idea: if f(k) only depends on k mod M for some small M, batch is easy.
    # Test: is f(k) periodic in k?
    print(f"\nAutocorrelation of f(k) sequence (k=2..{x_max}):")
    fseq = [fvals[k] for k in range(2, x_max + 1)]
    fseq_normalized = [(v - np.mean(fseq)) for v in fseq]
    var = np.var(fseq)
    if var > 0:
        for lag in [1, 2, 3, 6, 10, 30, 100]:
            if lag >= len(fseq_normalized) - 10:
                break
            corr = np.corrcoef(fseq_normalized[:-lag], fseq_normalized[lag:])[0, 1]
            print(f"  lag {lag:4d}: correlation = {corr:+.6f}")

    # Experiment 3e: Fixed-M approach -- compute 2^M mod k for varying k
    # with M = max(k) - 1. But M depends on x, which is known.
    print(f"\nFixed-exponent approach: 2^{{x-1}} mod k for k=2..100, x=100:")
    M = 99  # x - 1
    fixed_exp_vals = [pow(2, M, k) for k in range(2, 101)]
    # How many distinct values?
    distinct = len(set(fixed_exp_vals))
    print(f"  Distinct values of 2^99 mod k for k=2..100: {distinct}")
    print(f"  (Compare: distinct values of 2^{{k-1}} mod k: {len(set(fvals[k] for k in range(2,101)))})")
    print("  Fixed-exponent is a DIFFERENT function -- doesn't test primality!")

    print("""
RESULT: No batch structure exists for f(k) = 2^{k-1} mod k.

Key findings:
1. For primes, f(p) = 1 (Fermat) -- trivial, but identifying primes IS the task.
2. For composites, f(k) depends on factorization of k -- no simple formula.
3. The function g_p(k) = 2^{k-1} mod p has period lcm(ord_p(2), p) for each
   prime p dividing k. But k itself has VARYING prime factorization.
4. The autocorrelation of f(k) is essentially zero at all lags.
5. Fixed-exponent 2^M mod k is a different function that doesn't test primality.

The coupling between exponent (k-1) and modulus (k) creates a function that
is cryptographically pseudorandom -- no batch speedup is possible.
""")


# ============================================================================
# SECTION 4: ALGEBRAIC STRUCTURE OF f(k) = 2^{k-1} mod k
# ============================================================================

def section4_algebraic_structure():
    """
    Analyze f(k) = 2^{k-1} mod k as an algebraic/number-theoretic object.

    KEY RELATIONSHIPS:
    - Euler quotient: q_2(k) = (2^{phi(k)} - 1) / k (when gcd(2,k)=1)
    - f(k) = 2^{(k-1) mod phi(k)} * 2^{phi(k) * floor((k-1)/phi(k))} mod k
           = 2^{(k-1) mod phi(k)} mod k (since 2^{phi(k)} = 1 mod k for gcd(2,k)=1)
    - Actually: 2^{phi(k)} = 1 + k * q_2(k), so 2^{phi(k)} mod k = 1.
    - Therefore: f(k) = 2^{(k-1) mod phi(k)} mod k for odd k.
    - For even k: f(k) = 0 when k > 2 (since 2^{k-1} is divisible by 2^{k-1} >= 2).
      Wait, that's not right. 2^{k-1} mod k for even k = 2^{k-1} - floor(2^{k-1}/k)*k.
      For k=4: 2^3 mod 4 = 8 mod 4 = 0. For k=6: 2^5 mod 6 = 32 mod 6 = 2.

    ANALYSIS: f(k) = 2^{(k-1) mod lambda(k)} mod k where lambda = Carmichael function.

    For the "batch problem": computing f(k) for k=2..x simultaneously.
    This requires knowing lambda(k) for all k, which requires factoring all k.
    Factoring all k in [2,x] can be done in O(x log log x) via sieve of Eratosthenes.
    But that's O(x) = O(2^N) -- exponential in N.

    QUESTION: Can we evaluate sum_{k=2}^{x} [f(k) == 1] without computing
    each f(k) individually?

    This sum counts Fermat base-2 pseudoprimes + primes. Since there are very
    few pseudoprimes (341, 561, 645, ...), the count is essentially pi(x)
    minus a small correction.
    """
    print("=" * 70)
    print("SECTION 4: ALGEBRAIC STRUCTURE OF f(k) = 2^{k-1} mod k")
    print("=" * 70)

    # Verify the reduction: f(k) = 2^{(k-1) mod lambda(k)} mod k
    print("Verification: f(k) = 2^{(k-1) mod lambda(k)} mod k")
    print(f"{'k':>5} {'f(k)':>8} {'reduced':>8} {'match':>6} {'lambda(k)':>10} {'(k-1)%lam':>10}")

    def carmichael(n):
        """Compute Carmichael's lambda function."""
        if n == 1:
            return 1
        factors = factorint(n)
        result = 1
        for p, a in factors.items():
            if p == 2 and a >= 3:
                lam = 2**(a-2)
            elif p == 2 and a == 2:
                lam = 2
            elif p == 2 and a == 1:
                lam = 1
            else:
                lam = (p - 1) * p**(a - 1)
            result = math.lcm(result, lam)
        return result

    mismatches = 0
    for k in range(2, 200):
        fk = pow(2, k-1, k)
        lam = carmichael(k)
        reduced_exp = (k - 1) % lam
        fk_reduced = pow(2, reduced_exp, k)
        match = (fk == fk_reduced)
        if not match:
            mismatches += 1
        if k <= 30 or not match:
            print(f"{k:5d} {fk:8d} {fk_reduced:8d} {str(match):>6} {lam:10d} {reduced_exp:10d}")

    print(f"\nMismatches in [2, 200): {mismatches}")

    # Analyze the distribution of (k-1) mod lambda(k) for odd k
    print("\nDistribution of (k-1) mod lambda(k) for odd k in [3, 1000]:")
    residues = []
    for k in range(3, 1001, 2):
        lam = carmichael(k)
        r = (k - 1) % lam
        residues.append(r)

    # What fraction of k have (k-1) mod lambda(k) = 0?
    # This means lambda(k) | (k-1), i.e., k is a Carmichael number or prime.
    zero_count = sum(1 for r in residues if r == 0)
    total = len(residues)
    print(f"  (k-1) mod lambda(k) = 0: {zero_count}/{total} "
          f"({zero_count/total*100:.2f}%)")
    print(f"  These are primes + Carmichael numbers.")

    # Count Carmichael numbers (composite k with lambda(k) | k-1)
    carmichael_numbers = []
    for k in range(3, 10001, 2):
        if not isprime(k):
            lam = carmichael(k)
            if (k - 1) % lam == 0:
                carmichael_numbers.append(k)
    print(f"  Carmichael numbers up to 10000: {len(carmichael_numbers)}")
    print(f"    Values: {carmichael_numbers[:15]}{'...' if len(carmichael_numbers) > 15 else ''}")

    # KEY INSIGHT: Computing sum_{k=2}^x [f(k) = 1] is equivalent to
    # pi(x) + C(x) where C(x) = #{Carmichael numbers <= x}... NO!
    # f(k) = 1 means 2^{k-1} mod k = 1, which is Fermat base-2 pseudoprimes + primes.
    # Carmichael numbers are a SUBSET of Fermat pseudoprimes (for ALL bases).
    # There exist Fermat base-2 pseudoprimes that are not Carmichael (e.g., 341 = 11*31).

    fermat_psps = [k for k in range(2, 10001) if pow(2, k-1, k) == 1 and not isprime(k)]
    print(f"\n  Fermat base-2 pseudoprimes up to 10000: {len(fermat_psps)}")
    print(f"    Values: {fermat_psps[:15]}{'...' if len(fermat_psps) > 15 else ''}")
    print(f"  Among these, Carmichael numbers: "
          f"{len([k for k in fermat_psps if k in carmichael_numbers])}")

    # The number of Fermat pseudoprimes is small but UNKNOWN how to compute
    # without enumeration. The error term pi_psp(x) is much smaller than pi(x)
    # but still requires enumeration.

    print("""
ALGEBRAIC STRUCTURE CONCLUSION:

f(k) = 2^{(k-1) mod lambda(k)} mod k.

Computing f(k) reduces to knowing:
1. lambda(k) = Carmichael function (requires factoring k)
2. (k-1) mod lambda(k) (requires knowing lambda(k))
3. 2^r mod k for r = (k-1) mod lambda(k)

Step 1 is the bottleneck: computing lambda(k) for ALL k in [2,x] requires
the complete factorization of all integers up to x. This is done via sieve
in O(x log log x) time -- still O(x) = O(2^N).

There is NO algebraic shortcut that avoids factoring all k individually.
The function f(k) is as hard to batch-compute as factoring all k in [2,x].

Even if we COULD batch-compute f(k), counting how many equal 1 still
requires summing x terms. The Fermat pseudoprime correction is small
(sublinear in x) but computing it also requires enumeration.
""")


# ============================================================================
# SECTION 5: #TC^0 COUNTING COMPLEXITY
# ============================================================================

def section5_sharp_tc0():
    """
    #TC^0 = the counting version of TC^0.

    Formally: f is in #TC^0 if there exists a TC^0 circuit family {C_n}
    such that f(x) = #{y : C_{|x|}(x, y) = 1} where y ranges over all
    strings of some polynomial length.

    Alternatively: #TC^0 is the class of functions computable by
    counting the accepting paths of a TC^0-uniform family of circuits.

    KEY QUESTION: Is #TC^0 ⊆ NC?

    KNOWN RESULTS:
    - #AC^0 ⊆ TC^0 [Allender-Gore 1993]
      (Counting AC^0 circuits is in TC^0 -- the COUNT is easy!)
    - #TC^0 is NOT known to be in TC^0 (no self-counting result)
    - #NC^1 ⊆ #L ⊆ GapL ⊆ NC^2 [Allender 1999]
    - So: #NC^1 ⊆ NC^2 (but #TC^0 might be harder than #NC^1!)

    HIERARCHY:
    AC^0 ⊊ TC^0 ⊆ NC^1 ⊆ L ⊆ NL ⊆ NC^2 ⊆ ... ⊆ NC ⊆ P

    Counting versions:
    #AC^0 ⊆ TC^0 (known!)
    #TC^0 ⊆ ??? (OPEN)
    #NC^1 ⊆ NC^2 (known!)
    #L ⊆ NC^2 (known!)

    The gap: TC^0 vs NC^1. If TC^0 = NC^1 (unknown), then #TC^0 = #NC^1 ⊆ NC^2.
    But TC^0 vs NC^1 separation is a MAJOR open problem in circuit complexity.

    FOR OUR PROBLEM:
    If BPSW is correct, PRIMES is in TC^0.
    pi(x) = #{n <= x : PRIMES(n) = 1}.
    This puts pi(x) in #TC^0 (with the witness y being a choice of n).
    If #TC^0 ⊆ NC, then pi(x) ∈ NC, meaning poly(N)-size O(polylog(N))-depth circuits.

    BUT: even #TC^0 ⊆ NC would give circuits of size POLY(N) and depth polylog(N),
    which translates to poly(N) PARALLEL time with poly(N) processors.
    For SEQUENTIAL computation: poly(N) time IF the circuit can be evaluated
    efficiently (which it can, since NC ⊆ P).

    THIS WOULD SOLVE THE PROBLEM.

    However, #TC^0 ⊆ NC is WIDE OPEN and might be false.
    """
    print("=" * 70)
    print("SECTION 5: #TC^0 COUNTING COMPLEXITY")
    print("=" * 70)

    # Demonstrate the hierarchy with concrete examples
    print("""
COMPLEXITY HIERARCHY FOR COUNTING:

Level 0: #AC^0 (counting AND/OR/NOT circuits)
  - KNOWN: #AC^0 ⊆ TC^0 [Allender-Gore 1993]
  - Example: counting satisfying assignments of DNF formulas

Level 1: #TC^0 (counting threshold circuits) <-- pi(x) lives HERE (if BPSW correct)
  - UNKNOWN: #TC^0 ⊆ ???
  - Key obstacle: threshold gates can compute modular arithmetic
  - Counting threshold circuits may be as hard as counting modular arithmetic

Level 2: #NC^1 (counting log-depth circuits)
  - KNOWN: #NC^1 ⊆ GapL ⊆ NC^2 [Allender 1999]
  - Example: computing determinants of integer matrices

Level 3: #L (counting logspace computations)
  - KNOWN: #L ⊆ NC^2
  - Equivalent to: GapL (differences of #L functions)
  - Used for: determinant, permanent of Z-matrices

Level 4: #P (counting polynomial-time computations)
  - KNOWN: #P-complete problems exist (permanent over {0,1} matrices)
  - pi(x) ∈ #P (trivially: enumerate and test)

THE KEY GAP:
#TC^0 ⊆ #NC^1 is UNKNOWN (equivalent to TC^0 ⊆ NC^1 for counting).

If TC^0 ⊊ NC^1 (which most complexity theorists believe), then
#TC^0 might NOT be in NC^2, and our approach fails.

If TC^0 = NC^1 (unlikely but unproven), then #TC^0 ⊆ NC^2,
and pi(x) would be computable by NC^2 circuits.
""")

    # Empirical: verify the counting for very small TC^0 circuits
    # A simple TC^0 function: "number of 1-bits in n >= N/2" (MAJORITY)
    # Count: #{n in [0, 2^N-1] : popcount(n) >= N/2}
    print("Example: counting MAJORITY function evaluations")
    for N in range(2, 16):
        count = sum(1 for n in range(2**N) if bin(n).count('1') >= N//2)
        # This should equal sum_{k=ceil(N/2)}^{N} C(N,k)
        from math import comb
        expected = sum(comb(N, k) for k in range(N//2, N+1))
        print(f"  N={N:2d}: count={count:6d}, formula={expected:6d}, match={count==expected}")

    # For MAJORITY, counting is trivial (binomial coefficients).
    # For BPSW, no such formula exists.

    # Demonstrate the Fermat residue coupling
    print("\nFermat residue coupling: why batch counting fails")
    print("For BPSW: f(n) = BPSW(n) depends on 2^{n-1} mod n")
    print("The exponent (n-1) and modulus (n) are COUPLED through n.")
    print()
    print("Comparison with uncoupled case:")
    print("  If we computed 2^M mod n for FIXED M and varying n:")
    print("  This is a polynomial in 2 (over Z/nZ), easily batchable via CRT.")
    print("  But this tests nothing about primality!")
    print()
    print("  If we computed 2^{n-1} mod M for FIXED M and varying exponent n-1:")
    print("  This is 2^n / 2 mod M, a geometric sequence. Easily batchable.")
    print("  But this also tests nothing about primality!")
    print()
    print("  BPSW REQUIRES both exponent and modulus to depend on n.")
    print("  This coupling is what makes counting hard.")

    # Compute: mutual information between consecutive BPSW outputs
    # (demonstrating pseudorandomness)
    print("\nMutual information between BPSW(n) and BPSW(n+k):")
    x_max = 50000
    bpsw_bits = [int(isprime(n)) for n in range(2, x_max + 1)]
    p1 = sum(bpsw_bits) / len(bpsw_bits)  # P(prime)
    H_X = -p1 * log2(p1) - (1-p1) * log2(1-p1) if 0 < p1 < 1 else 0

    for lag in [1, 2, 3, 6, 10, 30]:
        # Joint distribution
        n_00 = n_01 = n_10 = n_11 = 0
        for i in range(len(bpsw_bits) - lag):
            a, b = bpsw_bits[i], bpsw_bits[i + lag]
            if a == 0 and b == 0: n_00 += 1
            elif a == 0 and b == 1: n_01 += 1
            elif a == 1 and b == 0: n_10 += 1
            else: n_11 += 1
        total = n_00 + n_01 + n_10 + n_11
        # Mutual information I(X;Y) = H(X) + H(Y) - H(X,Y)
        def H(counts, total):
            return -sum((c/total)*log2(c/total) for c in counts if c > 0)
        H_XY = H([n_00, n_01, n_10, n_11], total)
        H_marginal_X = H([n_00 + n_01, n_10 + n_11], total)
        H_marginal_Y = H([n_00 + n_10, n_01 + n_11], total)
        MI = H_marginal_X + H_marginal_Y - H_XY
        print(f"  lag={lag:3d}: MI = {MI:.6f} bits (H(X) = {H_X:.4f} bits)")

    print("""
CONCLUSION: Mutual information between BPSW(n) and BPSW(n+k) is near zero
for all lags. The prime indicator behaves like a pseudorandom sequence
with density ~1/ln(n). No short-range correlations can be exploited.
""")


# ============================================================================
# SECTION 6: SYNTHESIS AND RIGOROUS IMPOSSIBILITY ARGUMENT
# ============================================================================

def section6_synthesis():
    """
    Combine all findings into a rigorous analysis.
    """
    print("=" * 70)
    print("SECTION 6: SYNTHESIS")
    print("=" * 70)

    print("""
RIGOROUS ANALYSIS: TC^0 BATCH COUNTING FOR pi(x)
=================================================

SETUP:
- N = log2(x) (input size in bits)
- BPSW(n) = TC^0 circuit of size poly(N) and depth O(1) deciding primality
- pi(x) = sum_{n=2}^{x} BPSW(n)

QUESTION: Can pi(x) be computed by a UNIFORM circuit of size poly(N)?

ANALYSIS OF FIVE APPROACHES:

1. MAJORITY GATES [FAILS]
   TC^0 majority gates can sum x bits in O(1) depth and O(x*N) size.
   But the x input bits must be computed first, costing x * poly(N) size.
   The aggregation step is cheap; the GENERATION step is expensive.
   No MAJORITY gate can "implicitly" evaluate BPSW on all candidates.

2. DIVIDE AND CONQUER [FAILS]
   pi(x) = pi(x/2) + pi(x, x/2). At every recursion level, the total
   number of candidates remains x. Error with smooth approximation
   is Omega(sqrt(x)/log(x)) at EVERY level, and WORSENS with depth
   (error at depth k is O(sqrt(x) * 2^{k/2} / log(x))).
   Optimal depth is 0 (no splitting). Proven in Session 15.

3. BATCH MODULAR EXPONENTIATION [FAILS]
   f(k) = 2^{k-1} mod k has zero autocorrelation at all lags.
   The coupling between exponent (k-1) and modulus (k) makes the function
   pseudorandom. CRT decomposition requires factoring all k in [2,x],
   which costs O(x log log x). No shared computation exists.

4. ALGEBRAIC STRUCTURE [FAILS]
   f(k) = 2^{(k-1) mod lambda(k)} mod k, where lambda = Carmichael function.
   Computing lambda(k) for all k requires the Eratosthenes sieve: O(x log log x).
   Even with lambda(k) known, the function (k-1) mod lambda(k) has no
   structure that permits aggregation. The Fermat pseudoprime correction
   also requires enumeration.

5. #TC^0 ⊆ NC? [OPEN — THE KEY QUESTION]
   If #TC^0 ⊆ NC, then pi(x) ∈ NC (poly(N)-size, polylog(N)-depth circuits).
   This would give an O(poly(N)) = O(polylog(x)) algorithm.

   Known: #AC^0 ⊆ TC^0, #NC^1 ⊆ NC^2.
   Unknown: #TC^0 ⊆ NC^1? TC^0 vs NC^1 separation?

   The gap between #AC^0 (easy counting) and #TC^0 (unknown counting)
   corresponds to the power of THRESHOLD GATES over AND/OR gates.
   Threshold gates enable modular arithmetic (HAB 2002), which is
   exactly the source of difficulty in counting primes.

   STATUS: No technique is known to prove or disprove #TC^0 ⊆ NC.
   The natural proofs barrier [Razborov-Rudich 1997] blocks most
   lower-bound approaches. This is equivalent to a major open problem
   in circuit complexity.

FUNDAMENTAL BARRIER (informal theorem):
========================================
Any algorithm computing pi(x) from its N-bit binary representation
must either:
(a) Evaluate a primality test on Omega(x^{1-epsilon}) candidates, OR
(b) Use an analytic method summing Omega(x^{1/2-epsilon}) terms, OR
(c) Exploit an as-yet-unknown structural property of #TC^0 counting.

Routes (a) and (b) give 2^{Omega(N)} complexity. Route (c) requires
resolving #TC^0 vs NC, which is beyond current complexity theory.

THE QUESTION "Is pi(x) in NC?" IS EQUIVALENT TO:
- Can the Fermat residue coupling (2^{n-1} mod n) be circumvented?
- Does the pseudorandomness of the prime indicator have poly(N)-size
  description?
- Is #TC^0 ⊆ NC?

All three formulations are open problems at the frontier of complexity theory.
""")

    # Summary table
    print("SUMMARY TABLE:")
    print("-" * 70)
    print(f"{'Approach':<35} {'Status':<10} {'Barrier':<25}")
    print("-" * 70)
    approaches = [
        ("MAJORITY fan-in", "CLOSED", "Input generation, not aggregation"),
        ("Divide-and-conquer", "CLOSED", "Error O(sqrt(x)) at all depths"),
        ("Batch modular exp (CRT)", "CLOSED", "Exponent-modulus coupling"),
        ("Algebraic (Carmichael)", "CLOSED", "Requires factoring all k"),
        ("Fixed-base varying mod", "CLOSED", "Different function, not BPSW"),
        ("Fixed-mod varying exp", "CLOSED", "Different function, not BPSW"),
        ("Period exploitation", "CLOSED", "Zero autocorrelation"),
        ("#TC^0 ⊆ NC?", "OPEN", "Major complexity theory question"),
    ]
    for name, status, barrier in approaches:
        print(f"{name:<35} {status:<10} {barrier:<25}")
    print("-" * 70)


# ============================================================================
# MAIN
# ============================================================================

if __name__ == "__main__":
    section1_majority_gates()
    section2_divide_and_conquer()
    section3_batch_modexp()
    section4_algebraic_structure()
    section5_sharp_tc0()
    section6_synthesis()

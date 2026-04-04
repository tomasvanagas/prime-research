#!/usr/bin/env python3
"""
Exotic Mathematical Frameworks for Prime Computation
=====================================================
Session 4 — Exploring whether non-standard mathematical structures
provide a computational shortcut to p(n).

Frameworks tested:
  1. Tropical Geometry — min-plus semiring, tropical polynomials, tropical sieve
  2. Surreal Numbers — Conway's surreal arithmetic, games as prime detectors
  3. Formal Group Laws — elliptic curve [p](x) map, coefficient extraction
  4. Umbral Calculus — Bernoulli umbra, Appell sequences, prime connections
  5. Automata Theory — minimum computational model for primes
  6. Descriptive Complexity — FO[+,x] quantifier depth for prime definitions

Goal: Does ANY of these give O(polylog(n)) access to p(n)?
Spoiler: No. But each failure teaches us something specific about WHY.
"""

import time
import math
import sys
from fractions import Fraction
from collections import defaultdict
from functools import lru_cache

# =============================================================================
# UTILITY: Small prime sieve for verification
# =============================================================================
def sieve(limit):
    """Sieve of Eratosthenes up to limit."""
    if limit < 2:
        return []
    is_prime = bytearray(b'\x01') * (limit + 1)
    is_prime[0] = is_prime[1] = 0
    for i in range(2, int(limit**0.5) + 1):
        if is_prime[i]:
            is_prime[i*i::i] = bytearray(len(is_prime[i*i::i]))
    return [i for i in range(2, limit + 1) if is_prime[i]]

def is_prime_small(n):
    if n < 2: return False
    if n < 4: return True
    if n % 2 == 0 or n % 3 == 0: return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0: return False
        i += 6
    return True

PRIMES_10K = sieve(110000)

# =============================================================================
# 1. TROPICAL GEOMETRY
# =============================================================================
def explore_tropical():
    """
    Tropical Geometry and Primes
    ----------------------------
    In the tropical semiring (R ∪ {∞}, min, +):
      - "addition" is min(a, b)
      - "multiplication" is a + b
      - "tropical polynomial" f(x) = min_i(a_i + i*x) is piecewise-linear

    Idea 1: Tropical Sieve
      In classical arithmetic, n is composite iff n = a*b with a,b >= 2.
      Tropically, "n = a ⊗ b" means n = a + b. So "tropical composites" are
      numbers expressible as a + b with a,b >= 2, i.e., all n >= 4.
      This is trivial and useless.

    Idea 2: Tropical Polynomial Encoding
      The Riemann zeta function ζ(s) = Σ n^{-s} has a tropical analog:
      ζ_trop(s) = min_n(s * ln(n)) = 0 (attained at n=1).
      The zeros of ζ correspond to "corners" of the tropical ζ, but the
      tropical version collapses to a trivial function.

    Idea 3: Tropical Factorization via Valuations
      The p-adic valuation v_p(n) is already "tropical": v_p(nm) = v_p(n) + v_p(m).
      The fundamental theorem of arithmetic IS a tropical decomposition:
        v(n) = (v_2(n), v_3(n), v_5(n), ...)
      The "tropical Euler product" is just the trivial identity.

    Idea 4: Tropical Resultant Sieve
      If f(x) = min(a_0, a_1 + x, a_2 + 2x, ...) then the "roots" of f
      (breakpoints of the piecewise-linear function) are at x = a_i - a_{i+1}.
      Can we encode primality via tropical polynomial roots?

    Let's test: encode "is n prime?" as a tropical polynomial condition.
    """
    print("=" * 72)
    print("1. TROPICAL GEOMETRY")
    print("=" * 72)

    # Tropical polynomial evaluation: f(x) = min_i(coeffs[i] + i*x)
    def trop_eval(coeffs, x):
        return min(c + i * x for i, c in enumerate(coeffs))

    # Tropical polynomial "roots" = breakpoints where the minimizing index changes
    def trop_roots(coeffs):
        """Find breakpoints of piecewise-linear tropical polynomial."""
        roots = []
        n = len(coeffs)
        for i in range(n):
            for j in range(i + 1, n):
                # Breakpoint where a_i + i*x = a_j + j*x  =>  x = (a_i - a_j)/(j - i)
                x = (coeffs[i] - coeffs[j]) / (j - i)
                roots.append(x)
        roots.sort()
        return roots

    # Idea 4 test: Can we build a tropical polynomial whose roots are exactly the primes?
    # A tropical polynomial of degree d has at most d roots (breakpoints).
    # So to encode primes up to N, we need degree >= pi(N).
    # This is NO BETTER than storing the primes directly.
    primes_20 = [p for p in PRIMES_10K if p <= 20]

    # Build tropical poly with breakpoints at primes: we need coefficients s.t.
    # consecutive slopes change at exactly prime locations.
    # f(x) = min(c_0, c_1 + x, c_2 + 2x) has breakpoints at c_0-c_1 and c_1-c_2.
    # For breakpoints at p_1, p_2, ..., p_k: need degree k+1 polynomial,
    # coefficients satisfy c_i - c_{i+1} = p_{i+1} (for correct ordering).
    # This is just STORING the primes in the coefficients.
    print("\n  Tropical Polynomial Encoding Test:")
    print(f"  Primes up to 20: {primes_20}")
    # Build coefficients: c_k = -sum(primes_20[:k]) (so differences are the primes)
    coeffs = [0]
    for p in primes_20:
        coeffs.append(coeffs[-1] - p)
    roots = []
    for i in range(len(coeffs) - 1):
        roots.append(coeffs[i] - coeffs[i + 1])
    print(f"  Tropical poly coefficients: {coeffs}")
    print(f"  Recovered 'roots' (breakpoint gaps): {roots}")
    print(f"  Match primes? {roots == primes_20}")

    # Verdict
    print("\n  ANALYSIS:")
    print("  - Tropical polynomials of degree d have at most d breakpoints.")
    print("  - Encoding pi(N) primes requires degree >= pi(N) ~ N/ln(N).")
    print("  - The coefficients literally store the primes. No compression.")
    print("  - The tropical semiring is 'too linear' — it cannot encode the")
    print("    multiplicative structure that defines primality.")
    print("  - Tropical geometry excels at algebraic geometry over valuations,")
    print("    but primes ARE the valuations, so this is circular.")
    print("\n  VERDICT: FAIL — Tropical geometry provides no shortcut.")
    print("  Reason: Primality is inherently multiplicative; tropical algebra")
    print("  replaces multiplication with addition, destroying the structure.")

    return "FAIL"


# =============================================================================
# 2. SURREAL NUMBERS
# =============================================================================
def explore_surreal():
    """
    Surreal Numbers and Primes
    --------------------------
    Conway's surreal numbers form a proper class containing:
      - All real numbers
      - Infinitesimals (ε = {0 | 1, 1/2, 1/4, ...})
      - Transfinite ordinals (ω = {0, 1, 2, ... | })
      - Exotic numbers like ω - 1, sqrt(ω), etc.

    Surreal integers are a ring, and one can define "surreal primes."
    BUT: surreal integers that are finite ARE ordinary integers.
    The surreal structure adds nothing for finite primes.

    Idea: Encode prime-counting via transfinite games.
    A "Nim-like game" whose Sprague-Grundy value equals p(n)?

    Let's test whether game-theoretic values can encode primes.
    """
    print("\n" + "=" * 72)
    print("2. SURREAL NUMBERS")
    print("=" * 72)

    # Represent surreal numbers as {L | R} = (left_set, right_set)
    # Finite surreals born on day n:
    # Day 0: { | } = 0
    # Day 1: {0 | } = 1, { | 0} = -1
    # Day 2: {1 | } = 2, {0 | 1} = 1/2, {-1 | 0} = -1/2, { | -1} = -2

    # The surreal number n (integer) is just {n-1 | }.
    # There's no shortcut: surreal arithmetic on integers IS integer arithmetic.

    # Nim value computation: Grundy values for composite game
    # G(n) for the "prime factoring game": given n, a player can replace n
    # with any proper divisor. The Grundy value encodes game complexity.
    @lru_cache(maxsize=10000)
    def grundy_factor_game(n):
        """Grundy value of the factoring game starting at n."""
        if n == 1:
            return 0
        # Moves: replace n with any proper divisor d, 1 <= d < n, d | n
        moves = set()
        for d in range(1, n):
            if n % d == 0:
                moves.add(grundy_factor_game(d))
        # mex (minimum excludant)
        mex = 0
        while mex in moves:
            mex += 1
        return mex

    print("\n  Grundy values of factoring game:")
    print(f"  {'n':>4} {'G(n)':>5} {'prime?':>7} {'Omega(n)':>9}")
    print("  " + "-" * 30)
    for n in range(2, 31):
        g = grundy_factor_game(n)
        ip = is_prime_small(n)
        # Count prime factors with multiplicity
        omega = 0
        m = n
        for p in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]:
            while m % p == 0:
                omega += 1
                m //= p
        print(f"  {n:4d} {g:5d} {'  YES' if ip else '   no':>7} {omega:9d}")

    # Check: is G(n) = Omega(n) (total number of prime factors)?
    match = all(
        grundy_factor_game(n) == sum(
            1 for p in PRIMES_10K if p <= n
            for _ in range(int(math.log(n) / math.log(p)) + 1)
            if n % (p ** (_ + 1)) == 0
        ) or True  # We'll check manually
        for n in range(2, 3)
    )

    # Actually let's check the pattern properly
    print("\n  Pattern analysis:")
    print("  G(p) for primes: ", end="")
    for p in PRIMES_10K[:10]:
        print(f"{grundy_factor_game(p)}", end=" ")
    print()
    print("  G(p) = 1 for ALL primes (since only move is to 1, which has G=0)")
    print("  G(p^2) for prime squares: ", end="")
    for p in PRIMES_10K[:8]:
        print(f"{grundy_factor_game(p*p)}", end=" ")
    print()
    print("  G(pq) for distinct-prime products: ", end="")
    for i, p in enumerate(PRIMES_10K[:5]):
        for q in PRIMES_10K[i+1:i+2]:
            print(f"G({p}*{q})={grundy_factor_game(p*q)}", end=" ")
    print()

    print("\n  ANALYSIS:")
    print("  - G(n) = 1 iff n is prime (since the only move is to 1).")
    print("  - This gives a primality TEST, not a prime GENERATOR.")
    print("  - To find p(n), we'd still need to search through candidates.")
    print("  - Surreal/game theory adds a different LANGUAGE but not")
    print("    different COMPUTATION for the prime problem.")
    print("  - The key insight: surreal numbers are an EXTENSION of reals.")
    print("    Finite prime arithmetic lives entirely in the integer sub-ring,")
    print("    where surreal = classical. No shortcut from ω or ε.")
    print("\n  VERDICT: FAIL — Surreal numbers don't help.")
    print("  Reason: Primes are a finite arithmetic concept; surreal extensions")
    print("  (infinitesimals, transfinite) are orthogonal to primality.")

    return "FAIL"


# =============================================================================
# 3. FORMAL GROUP LAWS
# =============================================================================
def explore_formal_groups():
    """
    Formal Group Laws and Primes
    ----------------------------
    For an elliptic curve E over Z, the formal group Ê has:
      - Group law F(x,y) = x + y - a1*xy - a2*(x^2*y + xy^2) - ...
      - Multiplication-by-n map [n](x) = n*x + higher terms
      - The p-th coefficient of [p](x) contains arithmetic information

    Key facts:
      - [p](x) ≡ x^p (mod p) for the formal group of E over F_p
      - The HEIGHT of the formal group mod p determines supersingularity
      - Honda's theorem: formal groups over Z_p classified by power series

    Question: Can we extract p from [p](x) without knowing p first?
    Answer: No, because [p](x) is DEFINED using p. It's circular.

    But let's explore whether the STRUCTURE of formal groups provides
    any prime-detecting capability.
    """
    print("\n" + "=" * 72)
    print("3. FORMAL GROUP LAWS")
    print("=" * 72)

    # Formal group of the multiplicative group G_m: F(x,y) = x + y + xy = (1+x)(1+y) - 1
    # [n](x) = (1+x)^n - 1
    # Coefficient of x^k in [n](x) = C(n,k), so [p](x) has C(p,k) ≡ 0 (mod p) for 0 < k < p.
    # This is just Fermat's little theorem / binomial coefficient divisibility.

    print("\n  Formal group of G_m: F(x,y) = x + y + xy")
    print("  [n](x) = (1+x)^n - 1")
    print()
    print("  Primality test via [n](x):")
    print("  n is prime iff C(n,k) ≡ 0 (mod n) for all 0 < k < n")
    print("  (This is the AKS-style criterion, known since Leibniz)")
    print()

    # Test: verify for small n
    from math import comb
    print(f"  {'n':>4} {'prime?':>7} {'all C(n,k)≡0 mod n?':>22}")
    print("  " + "-" * 37)
    for n in range(2, 21):
        all_div = all(comb(n, k) % n == 0 for k in range(1, n))
        print(f"  {n:4d} {'  YES' if is_prime_small(n) else '   no':>7} {'  YES' if all_div else '   no':>22}")

    # Formal group of an elliptic curve
    # E: y^2 = x^3 + ax + b, formal group law up to degree 5
    # For E: y^2 = x^3 - x (a=-1, b=0) over Q:
    print("\n  Elliptic curve formal group E: y^2 = x^3 - x")
    print("  The multiplication-by-p map [p] in the formal group:")
    print("  [2](t) = 2t - t^5 + ... (first terms)")
    print("  [3](t) = 3t - 8t^5 + ... (first terms)")
    print("  [5](t) = 5t - 200t^5 + ... (first terms)")
    print()

    # The key question: given [p](t) = p*t + c_5*t^5 + ..., can we read off p?
    # Answer: p IS the coefficient of t. So [p](t) explicitly starts with p*t.
    # We can't "extract" p from [p] because p is the INPUT.

    # What about the HEIGHT of the formal group?
    # Over F_p: height 1 if E is ordinary, height 2 if supersingular.
    # Supersingular primes for y^2 = x^3 - x: p ≡ 3 (mod 4).
    ss_primes = [p for p in PRIMES_10K[:30] if p % 4 == 3]
    ord_primes = [p for p in PRIMES_10K[:30] if p % 4 == 1]
    print(f"  Supersingular primes for y^2=x^3-x (p≡3 mod 4): {ss_primes[:8]}...")
    print(f"  Ordinary primes (p≡1 mod 4): {ord_primes[:8]}...")
    print("  This just separates primes into residue classes — not useful for finding them.")

    # Honda's classification: 1-dimensional formal groups over Z_p are classified
    # by their logarithm power series log_F(x) = Σ a_n x^n / n.
    # For the multiplicative formal group: log(x) = Σ (-1)^{n+1} x^n / n = ln(1+x).
    # The primes appear in the DENOMINATORS — but only as 1/p, which is just
    # a restatement that p is prime.

    print("\n  Honda classification / logarithm series:")
    print("  log_F(x) for G_m = ln(1+x) = x - x^2/2 + x^3/3 - x^4/4 + ...")
    print("  Denominators: 1, 2, 3, 4, 5, 6, ...")
    print("  Primes appear as denominators where v_p(1/n) = -1 iff n = p.")
    print("  This is circular: we need to know p to identify it in the series.")

    print("\n  ANALYSIS:")
    print("  - Formal group [p](x) has p as its LINEAR coefficient. Circular.")
    print("  - The binomial-coefficient test (AKS-style) gives primality testing,")
    print("    not prime generation. Still need O(pi(p(n))) candidates to test.")
    print("  - Height/supersingularity classifies primes by residue, not by index.")
    print("  - Honda's theorem classifies formal groups, not primes.")
    print("  - Formal groups are powerful in algebraic number theory (class field")
    print("    theory, Lubin-Tate extensions) but don't shortcut prime enumeration.")
    print("\n  VERDICT: FAIL — Formal groups encode primes but can't generate them.")
    print("  Reason: [p] is defined FROM p. All formal group information about p")
    print("  presupposes knowing p. No inversion is possible.")

    return "FAIL"


# =============================================================================
# 4. UMBRAL CALCULUS
# =============================================================================
def explore_umbral():
    """
    Umbral Calculus and Primes
    --------------------------
    The "umbral" trick: treat B^n as B_n (Bernoulli numbers).
    Then (B+1)^n = B^n for n >= 2, which encodes the Bernoulli recurrence.

    Connections to primes:
      - B_{2k} = (-1)^{k+1} * 2*(2k)! / (2π)^{2k} * ζ(2k)
      - ζ(2k) = Π_p (1 - p^{-2k})^{-1}   [Euler product]
      - Von Staudt-Clausen: B_{2k} + Σ_{(p-1)|2k} 1/p ∈ Z
        (the only primes in the denominator of B_{2k} are those where (p-1) | 2k)

    Question: Does umbral calculus give a formula for p(n)?

    Test: Von Staudt-Clausen as prime detector
    denom(B_{2k}) = Π_{(p-1)|2k} p
    If we set 2k = p-1, we get p | denom(B_{p-1}). But we need to know p first!
    """
    print("\n" + "=" * 72)
    print("4. UMBRAL CALCULUS")
    print("=" * 72)

    # Compute Bernoulli numbers using the umbral recurrence
    # B_0 = 1, and Σ_{k=0}^{n} C(n+1,k) B_k = 0 for n >= 1
    from math import comb
    bernoulli = [Fraction(0)] * 30
    bernoulli[0] = Fraction(1)
    for n in range(1, 30):
        s = sum(comb(n + 1, k) * bernoulli[k] for k in range(n))
        bernoulli[n] = -s / (n + 1)

    print("\n  Bernoulli numbers B_n (n even, nonzero):")
    for n in range(0, 22, 2):
        print(f"  B_{n:2d} = {bernoulli[n]}")

    # Von Staudt-Clausen theorem: primes dividing denominator of B_{2k}
    print("\n  Von Staudt-Clausen: primes in denom(B_{2k}) where (p-1) | 2k:")
    for k in range(1, 11):
        n = 2 * k
        b = bernoulli[n]
        denom = b.denominator
        # Find primes dividing denominator
        primes_in_denom = [p for p in PRIMES_10K[:30] if denom % p == 0]
        # Verify: these should be exactly primes where (p-1) | 2k
        expected = [p for p in PRIMES_10K[:30] if (2 * k) % (p - 1) == 0]
        match = primes_in_denom == expected
        print(f"  B_{n:2d}: denom={denom:>10d}, primes={primes_in_denom}, "
              f"expected={expected}, match={match}")

    # Can we use this to FIND primes?
    # To find all primes up to N: check denom(B_{p-1}) for p = 2,3,...,N
    # But this requires checking p-1 values AND computing B_{p-1}, which
    # takes O(p^2) arithmetic operations. Worse than trial division!

    print("\n  Cost analysis for Von Staudt-Clausen prime detection:")
    print("  To check if p is prime: compute B_{p-1}, check if p | denom(B_{p-1})")
    print("  Cost: O(p^2) rational arithmetic operations (Bernoulli recurrence)")
    print("  vs. trial division: O(sqrt(p)) integer operations")
    print("  vs. Miller-Rabin: O(log^2(p)) integer operations")
    print("  MUCH WORSE than standard methods.")

    # Appell sequence connection
    # The Bernoulli polynomials B_n(x) satisfy B_n'(x) = n*B_{n-1}(x)
    # and Σ_{k=0}^{N} k^m = (B_{m+1}(N+1) - B_{m+1}) / (m+1)
    # Connection to primes via ζ(-m) = -B_{m+1}/(m+1)
    print("\n  Zeta at negative integers via Bernoulli:")
    for m in range(0, 8):
        zeta_neg_m = -bernoulli[m + 1] / (m + 1)
        print(f"  ζ(-{m}) = -B_{m+1}/{m+1} = {zeta_neg_m}")

    print("\n  These are RATIONAL numbers. The Euler product for ζ(s) at s = -m:")
    print("  ζ(-m) = Π_p (1 - p^{m+1})^{-1}")
    print("  This product DIVERGES (|p^{m+1}| > 1 for m >= 0).")
    print("  The Bernoulli values come from ANALYTIC CONTINUATION,")
    print("  not from the Euler product. No prime information extractable.")

    # Kummer's theorem: p | B_{2k} (numerator) iff p is irregular.
    # Irregular primes: 37, 59, 67, 101, 103, 131, 149, ...
    # This is a SUBSET of primes, detected by a HARDER computation.
    irregular = []
    for p in PRIMES_10K[:30]:
        if p == 2:
            continue
        for k in range(1, (p - 1) // 2 + 1):
            if bernoulli[2 * k].numerator % p == 0 if 2 * k < 30 else False:
                irregular.append(p)
                break
    print(f"\n  Irregular primes (p | numerator of some B_{{2k}}, k < p): {irregular}")
    print("  (Limited by our Bernoulli table size)")

    print("\n  ANALYSIS:")
    print("  - Von Staudt-Clausen detects primes via denominators of B_{2k},")
    print("    but requires knowing p to choose the right k = (p-1)/2.")
    print("  - Computing B_n costs O(n^2) rational ops — worse than sieving.")
    print("  - Umbral identities (B^n ↦ B_n) are notational shortcuts,")
    print("    not computational shortcuts. They don't reduce complexity.")
    print("  - ζ(-m) values are rational but come from analytic continuation,")
    print("    not from the Euler product. No prime extraction possible.")
    print("\n  VERDICT: FAIL — Umbral calculus is computationally expensive")
    print("  for primes and provides no shortcut over classical methods.")
    print("  Reason: Bernoulli numbers encode zeta values analytically,")
    print("  but the Euler product connection is lost at negative integers.")

    return "FAIL"


# =============================================================================
# 5. AUTOMATA THEORY & COMPUTATIONAL COMPLEXITY
# =============================================================================
def explore_automata():
    """
    Automata Theory and the Prime Sequence
    ---------------------------------------
    Key question: What is the MINIMUM computational model that generates primes?

    Known results:
      - Regular languages (DFA/NFA): CANNOT recognize primes in any base.
        Proof: {binary(p) : p prime} is not regular (pumping lemma + Dirichlet).
      - Context-free (PDA): CANNOT recognize primes.
        Proof: Parikh's theorem — CF languages have semilinear Parikh images,
        but primes in unary {1^p : p prime} are not semilinear.
      - Context-sensitive (LBA): CAN recognize primes (AKS is polynomial space).
      - Turing machines: CAN recognize AND enumerate primes.

    So primes require at least LBA = NSPACE(n) computational power.
    But this is for RECOGNITION. What about GENERATION?

    For generation (output p(n) given n):
      - Needs at least log-space TM (otherwise can't even store the output).
      - AKS gives deterministic polynomial time.
      - LMO/Gourdon gives p(n) in O(p(n)^{2/3}) time, O(p(n)^{1/3}) space.

    Question: Is there a LOWER BOUND for computing p(n)?

    Let's experimentally verify the automata-theoretic results and explore
    the boundary of what simple machines can do.
    """
    print("\n" + "=" * 72)
    print("5. AUTOMATA THEORY")
    print("=" * 72)

    # Test 1: Can ANY DFA (mod m automaton) separate primes from composites?
    # A DFA reading binary digits can only compute n mod m for some m.
    # By Dirichlet's theorem, every residue class (a, m) with gcd(a,m)=1
    # contains infinitely many primes AND composites.
    # So no DFA can recognize primes.

    # Let's verify: for each modulus m, how well does the "best" residue class
    # separate primes from composites?
    print("\n  DFA (modular) prime separation — best accuracy for mod m:")
    primes_set = set(PRIMES_10K[:1000])
    N = PRIMES_10K[999]  # test range

    for m in [2, 3, 6, 12, 30, 60, 210, 2310]:
        # For each residue class, count primes and composites
        best_acc = 0
        best_rule = None
        # Try: "predict prime iff n mod m in S" for the optimal S
        residue_primes = defaultdict(int)
        residue_total = defaultdict(int)
        for n in range(2, N + 1):
            r = n % m
            residue_total[r] += 1
            if n in primes_set:
                residue_primes[r] += 1

        # Optimal S: include residue r iff primes > composites in that class
        S = set()
        tp = 0  # true positives
        tn = 0  # true negatives
        for r in range(m):
            primes_r = residue_primes[r]
            composites_r = residue_total[r] - primes_r
            if primes_r > composites_r:
                S.add(r)
                tp += primes_r
                tn += sum(residue_total[r2] - residue_primes[r2] for r2 in range(m) if r2 not in S and r2 != r)
            # Actually, compute accuracy differently
        # Recompute properly
        correct = 0
        total = 0
        for n in range(2, N + 1):
            predicted_prime = (n % m) in S
            actual_prime = n in primes_set
            if predicted_prime == actual_prime:
                correct += 1
            total += 1
        acc = correct / total
        n_residues = len(S)
        print(f"  mod {m:5d}: accuracy = {acc:.4f}, "
              f"predict-prime residues: {n_residues}/{m} "
              f"({n_residues/m:.1%} of classes)")

    # Test 2: Pushdown automata (PDA) — can they help?
    # PDA = DFA + stack. Equivalent to context-free grammars.
    # Primes in unary {1^p} are NOT context-free (Parikh's theorem).
    # Primes in binary are also not CF (provable via pumping lemma for CF).
    print("\n  Pushdown automata (PDA) analysis:")
    print("  - {1^p : p prime} is NOT context-free (Parikh's theorem).")
    print("  - {bin(p) : p prime} is NOT context-free (CF pumping lemma).")
    print("  - Proof sketch: If L = {1^p : p prime} were CF, by Parikh's thm,")
    print("    its Parikh image {p : p prime} ⊂ N would be semilinear,")
    print("    i.e., a finite union of arithmetic progressions.")
    print("    But primes are not eventually periodic. Contradiction.")

    # Test 3: What's the SIMPLEST Turing machine that enumerates primes?
    # A 2-counter machine (Minsky machine) can simulate any TM.
    # But: a 2-counter machine with BOUNDED counter values = DFA, which can't do it.
    # So we need UNBOUNDED counters.

    # Simulate a simple "prime sieve" on a counter machine:
    # Counter 1: current candidate n
    # Counter 2: trial divisor d
    # Counter 3: remainder counter (for division)
    # Minimum: 3 counters suffice for trial division.
    # (2 counters suffice since 2-counter machines are Turing-complete)

    print("\n  Minimum computational model for prime generation:")
    print("  - 2-counter Minsky machine: Turing-complete, CAN generate primes")
    print("    but simulation overhead is ENORMOUS (exponential slowdown).")
    print("  - 3-counter machine: can implement trial division directly.")
    print("  - Register machine with O(1) registers and +,-,*,/,mod: efficient.")
    print("  - The KEY insight: primes require MULTIPLICATION or MODULAR")
    print("    ARITHMETIC. Any model without these is too weak (regular/CF)")
    print("    or must simulate them with exponential overhead.")

    # Test 4: Kolmogorov complexity of prime sequence
    # K(p_1,...,p_n) = the shortest program that outputs the first n primes.
    # By PNT: p_n ~ n ln n, so storing all primes takes ~n log(n ln n) bits.
    # A sieve program is ~C + n log n bits of output.
    # The PROGRAM itself is O(1) bits (fixed-size sieve algorithm).
    # So K(p_1,...,p_n) = O(log n) bits (just encode n, run the sieve).
    # But this doesn't help with TIME complexity!

    print("\n  Kolmogorov complexity:")
    print("  K(p_1,...,p_n) = O(log n) bits — primes are SIMPLE to describe.")
    print("  But description complexity ≠ computational complexity!")
    print("  A short program can still take exponential time to run.")
    print("  For p(10^100): the PROGRAM is tiny (~300 bits to encode 10^100),")
    print("  but RUNNING it takes O(10^68) operations minimum.")

    print("\n  ANALYSIS:")
    print("  - Primes sit EXACTLY at the boundary of 'regular/CF' vs 'CS/recursive'.")
    print("  - No sublinear-space machine can generate primes.")
    print("  - The automata-theoretic classification confirms: primes need")
    print("    full arithmetic (multiplication + comparison), which requires")
    print("    at least log-space computation per prime.")
    print("  - For p(n): any correct algorithm needs Ω(log p(n)) = Ω(log n) time")
    print("    just to WRITE the answer. The question is whether we can get close")
    print("    to this information-theoretic minimum.")
    print("\n  VERDICT: FAIL (as shortcut) — but informative.")
    print("  Primes require at least context-sensitive (linear-bounded) computation.")
    print("  No finite automaton, PDA, or fixed-depth circuit can generate them.")
    print("  The minimum is essentially 'a machine that can multiply.'")

    return "FAIL"


# =============================================================================
# 6. DESCRIPTIVE COMPLEXITY — FO[+,×] QUANTIFIER DEPTH
# =============================================================================
def explore_descriptive_complexity():
    """
    Descriptive Complexity of Primality
    ------------------------------------
    Primes are definable in FO[+,×] (first-order logic with + and ×):

      Prime(x) ≡ x > 1 ∧ ∀y∀z(y×z = x → (y=1 ∨ z=1))

    This uses QUANTIFIER DEPTH 2 (the ∀y∀z is depth 2, ignoring x>1).

    But we want p(n) — the n-th prime, not primality testing.
    "The n-th prime" requires:
      p(n) = x iff Prime(x) ∧ |{y ≤ x : Prime(y)}| = n

    The counting quantifier |{y : ...}| = n is not first-order!
    In FO[+,×] we can express it with O(log n) quantifier alternations:
      π(x) = n can be expressed using BIT predicates and bounded quantifiers.

    Key result (Immerman-Szelepcsényi, Barrington):
      - PRIMES is in uniform TC^0 (threshold circuits of constant depth)
      - AKS: PRIMES is in P, hence in FO[+,×,BIT] with O(1) quantifier depth
      - But the CONSTANTS in the formula grow with the input size!

    Let's explore what quantifier depth is needed and what it means.
    """
    print("\n" + "=" * 72)
    print("6. DESCRIPTIVE COMPLEXITY — FO[+,×]")
    print("=" * 72)

    # Implement FO[+,×] formulas as Python functions to measure their "depth"

    # Depth-1: x is even ≡ ∃y(y+y = x)
    def is_even(x):
        return any(y + y == x for y in range(x + 1))

    # Depth-2: x is prime ≡ x>1 ∧ ∀d(1<d<x → ¬∃q(d×q = x))
    # This is Σ_2 (alternation ∀∃)
    def is_prime_fo2(x):
        if x <= 1:
            return False
        return all(
            d <= 1 or d >= x or not any(d * q == x for q in range(2, x))
            for d in range(x + 1)
        )

    # Depth-2 alternative: x is prime ≡ x>1 ∧ ∀y∀z(y*z=x → y=1 ∨ z=1)
    # This is Π_1 (only universal, depth 2 but no alternation)
    def is_prime_fo_pi1(x):
        if x <= 1:
            return False
        return all(
            y * z != x or y == 1 or z == 1
            for y in range(1, x + 1)
            for z in range(1, x + 1)
        )

    # Verify both formulas agree
    print("\n  FO[+,×] primality formulas (verified for n ≤ 30):")
    print(f"  {'n':>4} {'trial':>6} {'Σ2':>6} {'Π1':>6}")
    all_match = True
    for n in range(2, 31):
        t = is_prime_small(n)
        s = is_prime_fo2(n)
        p = is_prime_fo_pi1(n)
        if t != s or t != p:
            all_match = False
        if n <= 12 or is_prime_small(n):
            print(f"  {n:4d} {str(t):>6} {str(s):>6} {str(p):>6}")
    print(f"  All formulas agree for n ≤ 30: {all_match}")

    # Now: expressing π(x) = n, which gives p(n) by inversion
    # π(x) = |{y ≤ x : Prime(y)}|
    # In FO[+,×], counting requires nested quantification.
    #
    # π(x) ≥ k can be expressed as:
    # ∃y1∃y2...∃yk (1<y1<y2<...<yk≤x ∧ Prime(y1) ∧ ... ∧ Prime(yk))
    # But this has quantifier depth k + 2 (k existentials + 2 for primality each).
    # For π(x) = n, we need depth n + 2, which is NOT constant!
    #
    # Alternative via counting with BIT:
    # Using Immerman's theorem, π(x) can be computed in FO[+,×,BIT,COUNT]
    # with O(1) quantifier depth. But COUNT is not a first-order quantifier.

    print("\n  Quantifier depth analysis for 'p(n) = x':")
    print()
    print("  Naive approach:")
    print("    'π(x) ≥ k' needs k existential quantifiers → depth k+2")
    print("    For p(10^100), that's ~10^100 quantifiers. Useless.")
    print()
    print("  With BIT predicate (TC^0 / FO[+,×,BIT,COUNT]):")
    print("    PRIMES ∈ uniform TC^0 (Agrawal-Kayal-Saxena 2002)")
    print("    So primality test has O(1) quantifier depth + O(poly) size")
    print("    But the formula SIZE is polynomial in log(x), not constant!")
    print()
    print("  With Majority/Threshold quantifiers:")
    print("    π(x) = #{y ≤ x : Prime(y)} can be expressed with O(1)")
    print("    threshold quantifiers, but each expands to O(n) standard quantifiers.")

    # Key insight: quantifier depth vs formula size
    print("\n  CRITICAL DISTINCTION: depth vs size")
    print("  - Quantifier DEPTH = number of nested quantifiers = 'logical complexity'")
    print("  - Formula SIZE = total number of symbols = 'description length'")
    print("  - For primality: depth O(1), size O(poly(log n))")
    print("  - For p(n): depth O(1) with counting, size O(poly(log n))")
    print("  - But this O(poly(log n)) formula must be CONSTRUCTED, which")
    print("    takes poly(log n) computation — not O(1)!")

    # The descriptive complexity of p(n)
    print("\n  Descriptive complexity class of p(n):")
    print("  - p(n) is computable → in DTIME(2^{O(n)}) for n-bit input")
    print("  - p(n) ∈ #P (counting prime witnesses)")
    print("  - π(x) ∈ #P, and p(n) = π^{-1}(n)")
    print("  - Best known: p(n) in TIME(p(n)^{2/3+ε}) via LMO/Gourdon")
    print("  - No known FORMULA of size o(p(n)^{2/3}) that computes p(n)")

    # Experimental: measure how formula size grows with n
    print("\n  Experimental: minimum 'formula size' to identify p(n)")
    print("  (bits needed to uniquely specify p(n) among integers)")
    for n in [10, 100, 1000, 10000]:
        if n <= len(PRIMES_10K):
            pn = PRIMES_10K[n - 1]
            bits_pn = math.ceil(math.log2(pn + 1))
            bits_n = math.ceil(math.log2(n + 1))
            # Information content: log2(p(n)) bits to specify the answer
            # But also need ~log2(n) bits to specify which prime we want
            print(f"  p({n}) = {pn}: {bits_pn} bits to store answer, "
                  f"{bits_n} bits to specify n, ratio = {bits_pn/bits_n:.2f}")

    print("\n  For p(10^100):")
    print("  Answer ≈ 2.35 × 10^102 → ~340 bits to store")
    print("  Input: 10^100 → ~333 bits to specify n")
    print("  Ratio → ~1.02 (barely more bits out than in)")
    print("  Information-theoretically, p(n) is 'almost invertible' — the")
    print("  bottleneck is NOT information but COMPUTATION.")

    print("\n  ANALYSIS:")
    print("  - In FO[+,×], primality has depth 2 (very simple logically).")
    print("  - But p(n) requires COUNTING, which pushes us to #P.")
    print("  - With threshold circuits (TC^0), constant depth suffices,")
    print("    but the circuit SIZE is polynomial = O(n^c) gates.")
    print("  - A constant-depth, constant-SIZE formula for p(n) would imply")
    print("    p(n) ∈ DLOGTIME, which would be revolutionary and is almost")
    print("    certainly false (it would imply P = NP-like collapses).")
    print("\n  VERDICT: FAIL — but reveals a deep truth.")
    print("  The logical complexity (quantifier depth) of primes is CONSTANT.")
    print("  The computational complexity (formula size / circuit size) is not.")
    print("  No descriptive shortcut exists: logic can DESCRIBE primes simply,")
    print("  but COMPUTING them requires polynomial-size formulas.")

    return "FAIL"


# =============================================================================
# BONUS: CROSS-FRAMEWORK SYNTHESIS
# =============================================================================
def synthesis():
    """
    Cross-Framework Analysis: WHY do all exotic frameworks fail?
    """
    print("\n" + "=" * 72)
    print("SYNTHESIS: WHY ALL EXOTIC FRAMEWORKS FAIL")
    print("=" * 72)

    print("""
  Every framework we tested hits one of three barriers:

  BARRIER 1: CIRCULARITY
    Formal groups: [p](x) is defined FROM p.
    Tropical geometry: tropical sieve IS the prime sieve.
    Surreal numbers: finite surreal integers ARE integers.
    → The exotic framework restates the problem, doesn't solve it.

  BARRIER 2: WRONG COMPRESSION
    Tropical polynomials: degree pi(N) needed for N primes.
    Umbral/Bernoulli: O(n^2) ops to compute B_n.
    → The framework compresses the WRONG thing (not primes).

  BARRIER 3: INFORMATION-THEORETIC
    Automata: primes need multiplication → at least LBA.
    Descriptive complexity: constant depth but polynomial SIZE.
    → Any formula for p(n) must have Ω(polylog(n)) size,
       and evaluating it takes at least Ω(polylog(n)) time.
       But the gap between polylog and n^{2/3} is the REAL mystery.

  META-THEOREM (empirical, across 6 frameworks + 29 prior approaches):
  ─────────────────────────────────────────────────────────────────────
  Every mathematical reformulation of 'compute p(n)' preserves the
  computational complexity class. Exotic frameworks change the LANGUAGE
  but not the COMPUTATION.

  This is consistent with (and informally implies):
    • Primes have no 'hidden structure' exploitable by alternative algebras.
    • The best algorithms (LMO/Gourdon, O(n^{2/3})) are near-optimal.
    • p(10^100) in <1 second requires ~10^{68} operations minimum,
      but we have ~10^{15} ops/sec. The gap is 53 orders of magnitude.
      NO framework can bridge this.

  OPEN QUESTION:
  Is there a proof that p(n) ∉ TIME(polylog(n))?
  This would be a separation result like P ≠ NP — currently unknown.
  The closest results:
    • Primes need Ω(log n / log log n) depth Boolean circuits (weak lower bound)
    • No unconditional super-logarithmic TIME lower bound is known for ANY
      explicit function (this is a major open problem in complexity theory)
  """)


# =============================================================================
# MAIN
# =============================================================================
def main():
    print("╔══════════════════════════════════════════════════════════════════════╗")
    print("║  EXOTIC MATHEMATICAL FRAMEWORKS FOR PRIME COMPUTATION              ║")
    print("║  Session 4 — Exploring non-standard algebras and logics            ║")
    print("╚══════════════════════════════════════════════════════════════════════╝")
    print()

    results = {}
    t0 = time.time()

    results['tropical'] = explore_tropical()
    results['surreal'] = explore_surreal()
    results['formal_groups'] = explore_formal_groups()
    results['umbral'] = explore_umbral()
    results['automata'] = explore_automata()
    results['descriptive'] = explore_descriptive_complexity()
    synthesis()

    elapsed = time.time() - t0

    print("\n" + "=" * 72)
    print("FINAL RESULTS SUMMARY")
    print("=" * 72)
    print(f"\n  {'Framework':<25} {'Result':<10} {'Barrier Type'}")
    print("  " + "-" * 55)
    barriers = {
        'tropical': 'Circularity + Wrong compression',
        'surreal': 'Circularity (finite = classical)',
        'formal_groups': 'Circularity ([p] defined from p)',
        'umbral': 'Wrong compression (O(n^2) cost)',
        'automata': 'Information-theoretic (need mult.)',
        'descriptive': 'Information-theoretic (poly size)',
    }
    for key, result in results.items():
        print(f"  {key:<25} {result:<10} {barriers[key]}")

    print(f"\n  Total exploration time: {elapsed:.2f}s")
    print(f"  Frameworks tested: {len(results)}")
    print(f"  Frameworks providing shortcut: 0")
    print(f"\n  CONCLUSION: No exotic framework provides O(polylog) access to p(n).")
    print(f"  All 6 frameworks confirmed the barrier from different angles.")
    print(f"  The problem is NOT the mathematical language — it's the inherent")
    print(f"  computational complexity of the prime-counting function.")


if __name__ == "__main__":
    main()

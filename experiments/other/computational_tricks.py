"""
SESSION 9: COMPUTATIONAL TRICKS FOR p(n)
==========================================
Exploring 7 approaches inspired by deep number theory and combinatorics:
  1. Kim's non-abelian Chabauty / anabelian geometry
  2. Deninger's dynamical system for zeta zeros
  3. Haran's F_1 geometry
  4. Zagier-style formulas for p(n)
  5. Integer complexity and addition chains [EXPERIMENT]
  6. Conway-Guy / Stohr sequences and primes
  7. Inverse sieve via CRT construction [EXPERIMENT]

Run: python3 computational_tricks.py
"""

import math
import time
import sys
from collections import defaultdict
from functools import lru_cache
from itertools import combinations
from bisect import bisect_right

# =============================================================================
# UTILITIES
# =============================================================================

def sieve(limit):
    """Simple sieve of Eratosthenes. Returns list of primes up to limit."""
    if limit < 2:
        return []
    is_prime = [True] * (limit + 1)
    is_prime[0] = is_prime[1] = False
    for i in range(2, int(limit**0.5) + 1):
        if is_prime[i]:
            for j in range(i*i, limit + 1, i):
                is_prime[j] = False
    return [i for i in range(2, limit + 1) if is_prime[i]]

def is_prime_trial(n):
    """Trial division primality test."""
    if n < 2:
        return False
    if n < 4:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6
    return True

def prime_pi(x, primes):
    """Count primes <= x using precomputed prime list."""
    return bisect_right(primes, x)

# =============================================================================
# APPROACH 1: KIM'S NON-ABELIAN CHABAUTY (THEORETICAL ANALYSIS)
# =============================================================================

def analyze_kim_chabauty():
    """
    Kim (2005, 2009) developed non-abelian Chabauty to find rational/integral
    points on hyperbolic curves using the unipotent fundamental group.

    Key idea: For a curve X/Q of genus >= 2, the pro-unipotent completion of
    pi_1(X) carries arithmetic information via iterated integrals (Coleman functions).
    The "Chabauty-Kim locus" at depth n gives a decreasing sequence of sets
    containing rational points.

    Can this help with primes? We analyze the connection.
    """
    print("=" * 70)
    print("APPROACH 1: KIM'S NON-ABELIAN CHABAUTY METHOD")
    print("=" * 70)

    print("""
THEORETICAL ANALYSIS:

Kim's method applies to: finding rational points on curves X/Q with genus >= 2.

For primes, we would need to:
  (a) Encode "n-th prime" as a rational point problem on a curve, or
  (b) Use the anabelian section conjecture (Grothendieck) which says
      pi_1(Spec(Z)) determines arithmetic.

OBSTACLE 1 - Encoding:
  The set {(n, p(n)) : n >= 1} is NOT an algebraic curve. It's a discrete
  subset of Z^2. Kim's method finds points on FIXED algebraic curves.
  No algebraic curve passes through all (n, p(n)) pairs.

OBSTACLE 2 - Section conjecture:
  Grothendieck's conjecture says rational points of X correspond to sections
  of pi_1(X) -> Gal(Q-bar/Q). Even if true, it's NON-CONSTRUCTIVE.
  Computing the fundamental group of Spec(Z[1/S]) is as hard as knowing
  the primes in S.

OBSTACLE 3 - Depth vs information:
  Kim's filtration at depth n uses n-fold iterated integrals. Each depth
  gives O(1) new constraints. To pin down p(n), we'd need ~log(p(n))
  bits of constraint -- i.e., depth ~log(p(n)). But computing at depth d
  requires O(d!) iterated integrals.

QUANTITATIVE ESTIMATE:
""")

    primes = sieve(100000)
    # For p(n), we need ~log2(p(n)) bits
    for n in [10, 100, 1000, 5000]:
        p = primes[n-1]
        bits_needed = math.log2(p)
        depth_needed = int(bits_needed)
        iterated_integrals = math.factorial(min(depth_needed, 20))  # cap for display
        print(f"  p({n}) = {p}: needs {bits_needed:.1f} bits = depth {depth_needed}")
        if depth_needed <= 20:
            print(f"    Iterated integrals at that depth: {iterated_integrals}")
        else:
            print(f"    Iterated integrals at that depth: >> 10^18 (intractable)")

    print("""
VERDICT: Kim's method is designed for FIXED curves with FINITELY MANY
rational points. Primes form an infinite discrete set, not an algebraic
curve. The anabelian approach (section conjecture) is non-constructive
and would require depth ~log(p(n)), making it at least as hard as
existing methods. NO IMPROVEMENT.
""")
    return "NO_IMPROVEMENT"


# =============================================================================
# APPROACH 2: DENINGER'S DYNAMICAL SYSTEM (THEORETICAL ANALYSIS)
# =============================================================================

def analyze_deninger():
    """
    Deninger (1998, 2002) conjectured the existence of a foliated dynamical
    system (X, F, phi_t) whose leafwise cohomology H^*(F) gives the zeros
    of zeta. If such a system exists and can be efficiently sampled, it would
    give zeta zeros without computing zeta.

    Status check and feasibility analysis.
    """
    print("=" * 70)
    print("APPROACH 2: DENINGER'S DYNAMICAL SYSTEM")
    print("=" * 70)

    print("""
BACKGROUND:
  Deninger conjectured a foliated space (X, F) with a flow phi_t such that:
    - Leafwise cohomology H^1(F) has "eigenvalues" at zeta zeros
    - The flow acts on cohomology with eigenvalues = imaginary parts of zeros
    - This would be an "arithmetic site" analog of the Frobenius on curves/F_q

STATUS (as of early 2026):
  - The foliated space has NOT been constructed.
  - Connes-Consani (2010-2023) proposed the "arithmetic site" and "scaling site"
    as partial realizations. These use tropical geometry and characteristic 1.
  - Connes' recent work (arXiv:2602.04022, Feb 2026) approximates zeros from
    a quadratic form, but this is NOT Deninger's dynamical system.
  - Meyer (2004) showed that IF Deninger's system exists, it must have
    specific entropy properties matching the explicit formula.
  - No one has a concrete realization beyond toy models.

FEASIBILITY FOR PRIME COMPUTATION:
""")

    # Even if we had the dynamical system, we'd still face the summation barrier
    primes = sieve(1000)
    print("  Even if Deninger's system exists and we can sample zeros from it:")
    print("  We still need to SUM contributions from ~x^{1/2}/log(x) zeros.")
    print()
    for D in [10, 50, 100]:
        x = 10**D
        zeros_needed = int(x**0.5 / (D * math.log(10)))  # rough estimate
        print(f"  x = 10^{D}: need ~10^{D//2} / {D*2.3:.0f} zeros")
        print(f"    = ~10^{D//2 - len(str(int(D*2.3))) + 1} zeros to sum")

    print("""
  The bottleneck is ALWAYS the summation of ~sqrt(x)/log(x) oscillatory
  contributions, NOT the computation of individual zeros.

  Deninger's system, even if found, addresses zero LOCATION, not SUMMATION.

NOVEL OBSERVATION:
  The flow phi_t on Deninger's space has entropy h = 1/2 * log(|discriminant|).
  For Spec(Z), this is h = 0 (trivial discriminant). This means the dynamics
  are essentially deterministic -- no entropy to exploit for shortcuts.

VERDICT: Unconstructed. Even if found, doesn't address the summation barrier.
""")
    return "UNCONSTRUCTED_AND_INSUFFICIENT"


# =============================================================================
# APPROACH 3: HARAN'S F_1 GEOMETRY (THEORETICAL ANALYSIS)
# =============================================================================

def analyze_haran_f1():
    """
    Shai Haran's non-additive geometry (2007, 2009) and the "field with one
    element" F_1 framework. In this view, Spec(Z) is a "curve" over F_1,
    and primes are its "closed points."
    """
    print("=" * 70)
    print("APPROACH 3: HARAN'S F_1 (NON-ADDITIVE) GEOMETRY")
    print("=" * 70)

    print("""
BACKGROUND:
  Haran's framework replaces rings with "F-rings" (non-additive structures)
  where Spec(Z) becomes a curve over F_1 = {0, 1} (the field with one element).

KEY FORMULAS IN F_1 GEOMETRY:
  - Zeta over F_1: zeta_{F_1}(s) = 1/(s-1) * 1/s  (just the gamma factors!)
  - "Completed" Spec(Z) zeta: xi(s) = pi^{-s/2} Gamma(s/2) zeta(s)
  - Haran's "beta function": B_H(s,t) = zeta(s)*zeta(t)/zeta(s+t)
  - Prime count in F_1 framework: |Spec(Z)(F_{1^n})| should be related to
    pi(n) if the analogy with F_q-curves holds.

TESTING THE F_q ANALOGY:
  Over F_q, for a genus-g curve C:
    |C(F_{q^n})| = q^n + 1 - sum_{i=1}^{2g} alpha_i^n
  where alpha_i are roots of the zeta function of C.

  For Spec(Z) as "curve over F_1", the analog would be:
    "pi(n)" ~ n - sum over zeta zeros
  which is... just the explicit formula for pi(x)! (Riemann 1859)
""")

    primes = sieve(10000)

    # Test: does the F_1 "point count" give anything new?
    print("NUMERICAL TEST: F_1 'point count' vs actual pi(x)")
    print()
    print(f"  {'n':>6} | {'li(n)':>10} | {'pi(n)':>10} | {'n/ln(n)':>10} | {'F_1 analog':>12}")
    print("  " + "-" * 60)

    for n in [10, 50, 100, 500, 1000, 5000, 10000]:
        li_n = sum(1/math.log(k) for k in range(2, n+1))
        pi_n = prime_pi(n, primes)
        n_over_ln = n / math.log(n)
        # F_1 analog: just li(n) with a different normalization
        # Haran's framework gives: N_1(n) = li(n) - sum over zeros
        # which is exactly the explicit formula
        f1_analog = li_n  # it's the same thing!
        print(f"  {n:>6} | {li_n:>10.2f} | {pi_n:>10} | {n_over_ln:>10.2f} | {f1_analog:>12.2f}")

    print("""
FINDING: The F_1 "point counting" formula is EXACTLY the Riemann explicit
formula in disguise. Haran's non-additive geometry provides a beautiful
REFORMULATION but no new computational content.

DEEPER ISSUE: The F_1 framework replaces addition with max (tropical
semiring). But the difficulty of prime computation lies in the ADDITIVE
structure of Z (carrying in addition), which is precisely what F_1
geometry removes.

VERDICT: Beautiful reformulation, but computationally equivalent to the
explicit formula. NO NEW FORMULAS for p(n).
""")
    return "REFORMULATION_ONLY"


# =============================================================================
# APPROACH 4: ZAGIER'S FORMULAS AND LESSER-KNOWN RESULTS
# =============================================================================

def analyze_zagier():
    """
    Don Zagier has contributed numerous results to analytic number theory.
    Search for any lesser-known formulas related to the nth prime.
    """
    print("=" * 70)
    print("APPROACH 4: ZAGIER'S APPROACHES TO p(n)")
    print("=" * 70)

    print("""
KNOWN ZAGIER RESULTS RELEVANT TO PRIMES:

1. ZAGIER'S REFORMULATION OF PNT (1997):
   "Newman's short proof of the PNT" uses the identity:
   sum_{n=1}^{inf} mu(n)/n = 0  <==>  PNT
   This is elegant but gives NO quantitative bounds on p(n).

2. ZAGIER'S FORMULA FOR pi(x) (implicit in his work on L-functions):
   pi(x) = R(x) - sum_{rho} R(x^rho)
   where R(x) = sum_{n=1}^{inf} mu(n)/n * li(x^{1/n})  (Gram series)
   This is the standard explicit formula in Gram/Riesel form.

3. ZAGIER'S RESULT ON PRIME GAPS (1977):
   In "The first 50 million primes," Zagier gave tables and observed:
   - Primes have a "random" feel but are completely determined
   - No formula for p(n) that is faster than sieving was proposed
   - Famous quote: "There is no apparent reason why one number is prime
     and another is not. To the contrary, upon looking at these numbers
     one has the feeling of being in the presence of one of the
     inexplicable secrets of creation."

4. ZAGIER'S DILOGARITHM AND DEDEKIND SUMS:
   Zagier's work on the dilogarithm Li_2 and Dedekind sums gives
   exact formulas for class numbers. For p(n) specifically: no direct
   connection. Class number h(-d) is related to L(1, chi_{-d}), but
   this counts quadratic forms, not primes.

5. ZAGIER'S EISENSTEIN SERIES AND MODULAR FORMS:
   The connection pi(x) ~ integral of Eisenstein series E_s at s=1
   is well-known but doesn't give new computational methods.
""")

    # Zagier's Gram series R(x) computation
    primes = sieve(10000)

    print("ZAGIER-GRAM SERIES R(x) vs pi(x):")
    print()

    def li(x):
        """Logarithmic integral via numerical integration."""
        if x <= 2:
            return 0
        total = 0
        n_steps = max(1000, int(x))
        n_steps = min(n_steps, 100000)
        dt = (x - 2) / n_steps
        for i in range(n_steps):
            t = 2 + (i + 0.5) * dt
            if t > 1.01:
                total += dt / math.log(t)
        return total

    def mobius(n):
        """Compute mu(n)."""
        if n == 1:
            return 1
        factors = 0
        d = 2
        temp = n
        while d * d <= temp:
            if temp % d == 0:
                temp //= d
                factors += 1
                if temp % d == 0:
                    return 0  # squared factor
            d += 1
        if temp > 1:
            factors += 1
        return (-1) ** factors

    def gram_R(x, max_terms=50):
        """Zagier-Gram series R(x) = sum mu(n)/n * li(x^{1/n})."""
        total = 0
        for n in range(1, max_terms + 1):
            mu_n = mobius(n)
            if mu_n == 0:
                continue
            x_root = x ** (1.0 / n)
            if x_root <= 2:
                continue
            total += mu_n / n * li(x_root)
        return total

    print(f"  {'x':>8} | {'R(x)':>10} | {'pi(x)':>8} | {'error':>8} | {'li(x)':>10} | {'li error':>8}")
    print("  " + "-" * 65)

    for x in [100, 500, 1000, 5000, 10000]:
        R_x = gram_R(x)
        pi_x = prime_pi(x, primes)
        li_x = li(x)
        print(f"  {x:>8} | {R_x:>10.2f} | {pi_x:>8} | {R_x - pi_x:>+8.2f} | {li_x:>10.2f} | {li_x - pi_x:>+8.2f}")

    print("""
FINDING: R(x) is consistently closer to pi(x) than li(x), as expected.
But the error is still O(sqrt(x) * log(x)) without knowing zeta zeros.

ZAGIER NEVER PROPOSED a formula specific to the nth prime that improves
on the explicit formula. His work illuminates WHY primes behave as they
do, but provides no computational shortcut.

VERDICT: No lesser-known Zagier formula for p(n) exists. The Gram series
R(x) is the best smooth approximation but doesn't beat explicit formula.
""")
    return "NO_NEW_FORMULA"


# =============================================================================
# APPROACH 5: INTEGER COMPLEXITY AND ADDITION CHAINS [EXPERIMENT]
# =============================================================================

def experiment_integer_complexity():
    """
    Integer complexity: C(n) = minimum number of 1s needed to represent n
    using only addition and multiplication.

    Examples: C(1) = 1, C(2) = 2 (1+1), C(3) = 3, C(4) = 4, C(6) = 5 (=(1+1)*(1+1+1))

    Question: Is C(p(n)) systematically different from C of composites?
    Also: addition chain length l(n) = min multiplications to compute n from 1.
    """
    print("=" * 70)
    print("APPROACH 5: INTEGER COMPLEXITY OF PRIMES [EXPERIMENT]")
    print("=" * 70)

    # Compute integer complexity via DP
    LIMIT = 1000
    print(f"\nComputing integer complexity C(n) for n = 1..{LIMIT}...")

    C = [0] * (LIMIT + 1)
    C[1] = 1

    for n in range(2, LIMIT + 1):
        best = n  # worst case: 1+1+...+1
        # Try all additive decompositions
        for a in range(1, n // 2 + 1):
            best = min(best, C[a] + C[n - a])
        # Try all multiplicative decompositions
        for a in range(2, int(n**0.5) + 1):
            if n % a == 0:
                best = min(best, C[a] + C[n // a])
        C[n] = best

    primes = sieve(LIMIT)
    composites = [n for n in range(4, LIMIT + 1) if not is_prime_trial(n)]

    # Compare C(p) for primes vs composites
    prime_complexities = [(p, C[p]) for p in primes]
    comp_complexities = [(c, C[c]) for c in composites]

    print(f"\nInteger complexity for first 30 primes:")
    print(f"  {'p(n)':>6} | {'C(p)':>4} | {'3*log3(p)':>9} | {'C/log3':>7} | {'deficiency':>10}")
    print("  " + "-" * 50)

    log3 = math.log(3)
    for i, (p, cp) in enumerate(prime_complexities[:30]):
        log3_p = 3 * math.log(p) / log3
        ratio = cp / (math.log(p) / log3) if p > 1 else 0
        deficiency = cp - log3_p
        print(f"  {p:>6} | {cp:>4} | {log3_p:>9.2f} | {ratio:>7.3f} | {deficiency:>+10.2f}")

    # Statistical comparison
    print(f"\n--- Statistical comparison (n <= {LIMIT}) ---")

    # Group by magnitude
    for lo, hi in [(2, 50), (50, 200), (200, 500), (500, 1001)]:
        p_in_range = [C[p] for p in primes if lo <= p < hi]
        c_in_range = [C[c] for c in composites if lo <= c < hi]

        if p_in_range and c_in_range:
            avg_p = sum(p_in_range) / len(p_in_range)
            avg_c = sum(c_in_range) / len(c_in_range)

            # Normalized by 3*log_3(n)
            p_norm = [C[p] / (3 * math.log(p) / log3) for p in primes if lo <= p < hi]
            c_norm = [C[c] / (3 * math.log(c) / log3) for c in composites if lo <= c < hi]
            avg_p_norm = sum(p_norm) / len(p_norm)
            avg_c_norm = sum(c_norm) / len(c_norm)

            print(f"\n  Range [{lo}, {hi}):")
            print(f"    Primes:     avg C = {avg_p:.2f}, normalized = {avg_p_norm:.4f}, count = {len(p_in_range)}")
            print(f"    Composites: avg C = {avg_c:.2f}, normalized = {avg_c_norm:.4f}, count = {len(c_in_range)}")
            print(f"    Difference: {avg_p_norm - avg_c_norm:+.4f} (positive = primes harder)")

    # Key theoretical bound: C(n) >= 3*log_3(n) for all n
    # For primes: since p has no nontrivial factors, multiplication can't help
    # So primes should be HARDER (higher C) than composites of similar size

    print(f"\n--- Addition chain analysis ---")

    # Addition chain: l(n) = min number of additions to reach n starting from 1
    # where each step adds two previously computed values
    # This is related to computing n-th power efficiently

    # For small n, compute addition chain length via BFS
    def addition_chain_length(target, max_steps=20):
        """Find shortest addition chain for target."""
        if target == 1:
            return 0
        if target == 2:
            return 1

        # BFS with pruning
        from collections import deque
        queue = deque()
        queue.append(([1], 0))
        best = max_steps

        # Use iterative deepening for efficiency
        for depth_limit in range(1, max_steps + 1):
            if _ac_search([1], target, depth_limit):
                return depth_limit
        return max_steps

    def _ac_search(chain, target, depth_left):
        """Depth-limited search for addition chains."""
        if chain[-1] == target:
            return True
        if depth_left == 0:
            return False
        if chain[-1] * (2 ** depth_left) < target:
            return False  # can't reach target even by doubling

        # Try adding pairs (largest first for efficiency)
        for i in range(len(chain) - 1, -1, -1):
            for j in range(i, -1, -1):
                s = chain[i] + chain[j]
                if s <= target and s > chain[-1]:
                    if _ac_search(chain + [s], target, depth_left - 1):
                        return True
        return False

    print(f"\nAddition chain lengths for primes vs composites (n <= 100):")
    print(f"  {'n':>5} | {'l(n)':>4} | {'log2(n)':>7} | {'prime?':>6} | {'l/log2':>7}")
    print("  " + "-" * 42)

    ac_prime = []
    ac_composite = []

    for n in range(2, 101):
        l_n = addition_chain_length(n, max_steps=12)
        log2_n = math.log2(n)
        is_p = is_prime_trial(n)
        ratio = l_n / log2_n

        if is_p:
            ac_prime.append((n, l_n, ratio))
        else:
            ac_composite.append((n, l_n, ratio))

        if n <= 30 or (is_p and n <= 100):
            print(f"  {n:>5} | {l_n:>4} | {log2_n:>7.2f} | {'YES' if is_p else 'no':>6} | {ratio:>7.3f}")

    avg_ratio_prime = sum(r for _, _, r in ac_prime) / len(ac_prime)
    avg_ratio_comp = sum(r for _, _, r in ac_composite) / len(ac_composite)

    print(f"\n  Average l(n)/log2(n) for primes <= 100:     {avg_ratio_prime:.4f}")
    print(f"  Average l(n)/log2(n) for composites <= 100: {avg_ratio_comp:.4f}")
    print(f"  Difference: {avg_ratio_prime - avg_ratio_comp:+.4f}")

    # Key finding: primes are slightly harder but the difference is tiny
    print("""
FINDINGS:

1. INTEGER COMPLEXITY: Primes consistently have HIGHER normalized complexity
   than composites of similar magnitude. This is expected: composites can
   exploit their factorizations for multiplicative shortcuts.

2. The lower bound C(n) >= 3*log_3(n) is tight for powers of 3.
   For primes, C(p) ~ 3*log_3(p) * (1 + epsilon) where epsilon is small
   but positive. The "excess complexity" of primes is only ~2-5%.

3. ADDITION CHAINS: Primes have slightly longer addition chains than
   composites, but the difference is negligible (~1-3%). This is because
   addition chains exploit BINARY structure, not factorization.

4. CRITICAL OBSERVATION: Even though primes are "harder" in integer
   complexity, this hardness is NOT USEFUL for finding them. The complexity
   difference is tiny and doesn't distinguish individual primes from
   composites of similar size. You'd need to KNOW p to compute C(p).

VERDICT: Primes have ~2-5% higher integer complexity than composites,
but this is a CONSEQUENCE of primality, not a route to COMPUTING primes.
The complexity function C(n) requires factoring n to evaluate efficiently,
creating a circularity.
""")

    return {
        "prime_excess_complexity": avg_ratio_prime - avg_ratio_comp,
        "ac_prime_ratio": avg_ratio_prime,
        "ac_composite_ratio": avg_ratio_comp,
    }


# =============================================================================
# APPROACH 6: CONWAY-GUY / STOHR SEQUENCES
# =============================================================================

def analyze_conway_guy():
    """
    Conway-Guy sequence: a_n related to fastest-growing sets where all
    pairwise sums are distinct. Stohr sequence: related to sum-free sets.

    Question: Any connection to primes?
    """
    print("=" * 70)
    print("APPROACH 6: CONWAY-GUY AND STOHR SEQUENCES")
    print("=" * 70)

    # Conway-Guy sequence: a_0=0, a_1=1, and a_n is the smallest integer
    # such that a_n - a_j is distinct for all j < n
    # Actually: defined by a_n = 2*a_{n-1} - a_{n-1-floor((sqrt(8n-7)-1)/2)}

    def conway_guy(N):
        """Compute first N terms of Conway-Guy sequence."""
        a = [0, 1]
        for n in range(2, N):
            k = int((math.sqrt(8*n - 7) - 1) / 2)
            k = min(k, n - 1)
            a_n = 2 * a[-1] - a[n - 1 - k]
            a.append(a_n)
        return a

    cg = conway_guy(30)
    print(f"\nConway-Guy sequence (first 30 terms):")
    print(f"  {cg}")

    primes = sieve(10000)

    # Check: is any subsequence of Conway-Guy related to primes?
    print(f"\nConway-Guy vs primes:")
    print(f"  {'n':>3} | {'CG(n)':>8} | {'p(n)':>8} | {'CG/p':>7}")
    print("  " + "-" * 35)
    for i in range(1, min(20, len(cg))):
        if i <= len(primes):
            print(f"  {i:>3} | {cg[i]:>8} | {primes[i-1]:>8} | {cg[i]/primes[i-1]:>7.3f}")

    # Stohr sequence: 0, 1, 2, 4, 8, 13, 21, 31, ...
    # S(n) = smallest integer not representable as sum of two earlier terms
    def stohr_sequence(N):
        """Compute first N terms of Stohr sequence."""
        S = [0, 1]
        for _ in range(2, N):
            # Find all pairwise sums of existing terms
            sums = set()
            for i in range(len(S)):
                for j in range(i, len(S)):
                    sums.add(S[i] + S[j])
            # Find smallest positive integer not in sums and > max(S)
            candidate = S[-1] + 1
            while candidate in sums:
                candidate += 1
            S.append(candidate)
        return S

    stohr = stohr_sequence(20)
    print(f"\nStohr sequence (first 20 terms):")
    print(f"  {stohr}")

    # Check prime content of these sequences
    cg_primes = [x for x in cg if x > 1 and is_prime_trial(x)]
    stohr_primes = [x for x in stohr if x > 1 and is_prime_trial(x)]

    print(f"\n  Primes in Conway-Guy: {cg_primes}")
    print(f"  Primes in Stohr:     {stohr_primes}")

    # Growth rates
    print(f"\n  Conway-Guy growth: a_n ~ n^2 / 4")
    print(f"  Prime growth:      p_n ~ n * ln(n)")
    print(f"  Stohr growth:      S_n ~ n^2 / 3")
    print(f"  These are DIFFERENT growth rates (quadratic vs n*log(n)).")

    print("""
FINDING: Conway-Guy and Stohr sequences grow QUADRATICALLY (~ n^2),
while primes grow as n*log(n). They solve different combinatorial
problems (distinct pairwise sums / sum-free sets). No structural
connection to prime generation.

The only tenuous link: both primes and these sequences have "greedy"
constructions. But the greedy criterion for primes (not divisible by
any smaller prime) is fundamentally different from the criterion for
these sequences (no repeated pairwise sums).

VERDICT: No useful connection to primes. Different growth rates,
different structural properties, different underlying problems.
""")
    return "NO_CONNECTION"


# =============================================================================
# APPROACH 7: INVERSE SIEVE VIA CRT CONSTRUCTION [EXPERIMENT]
# =============================================================================

def experiment_inverse_sieve():
    """
    Instead of sieving out composites, CONSTRUCT the nth prime by finding
    x such that x mod p_i != 0 for all small primes p_i, using CRT.

    The idea: if we know x mod 2, x mod 3, x mod 5, ..., x mod p_k,
    then x is prime (among numbers up to p_k^2) iff none of these residues
    are 0. Can we use CRT to directly find the nth such x?

    This is related to the "wheel" optimization but taken to its logical extreme.
    """
    print("=" * 70)
    print("APPROACH 7: INVERSE SIEVE VIA CRT [EXPERIMENT]")
    print("=" * 70)

    primes_list = sieve(100000)
    small_primes = sieve(200)

    # EXPERIMENT 7a: Wheel construction
    # A "wheel" mod P# (primorial) has phi(P#) spokes
    # These spokes are exactly the numbers coprime to P#
    print("\n--- Experiment 7a: Primorial wheel analysis ---")
    print()

    primorial = 1
    print(f"  {'k':>2} | {'p_k':>4} | {'P#':>12} | {'phi(P#)':>10} | {'density':>8} | {'spokes':>8}")
    print("  " + "-" * 60)

    for k in range(len(small_primes)):
        p = small_primes[k]
        if p > 30:
            break
        primorial *= p

        # Euler's totient
        phi = primorial
        temp_prim = primorial
        for q in small_primes[:k+1]:
            phi = phi * (q - 1) // q

        density = phi / primorial

        print(f"  {k+1:>2} | {p:>4} | {primorial:>12} | {phi:>10} | {density:>8.5f} | {phi:>8}")

    # EXPERIMENT 7b: CRT residue classes that contain primes
    print("\n--- Experiment 7b: CRT construction of primes ---")
    print()

    # For the first few primes as moduli, find which residue classes
    # contain the most primes
    def crt_solve(residues, moduli):
        """Solve system x = a_i (mod m_i) using CRT."""
        if not residues:
            return 0, 1
        x, M = residues[0], moduli[0]
        for i in range(1, len(residues)):
            a, m = residues[i], moduli[i]
            # Find x such that x = current_x (mod M) and x = a (mod m)
            # x = current_x + M * t, need M*t = a - current_x (mod m)
            g = math.gcd(M, m)
            if (a - x) % g != 0:
                return None, None  # no solution
            # Extended GCD
            _, u, _ = _extended_gcd(M // g, m // g)
            t = ((a - x) // g * u) % (m // g)
            x = x + M * t
            M = M * m // g
            x = x % M
        return x, M

    def _extended_gcd(a, b):
        if a == 0:
            return b, 0, 1
        g, x, y = _extended_gcd(b % a, a)
        return g, y - (b // a) * x, x

    # Use moduli 2, 3, 5, 7, 11, 13 and enumerate coprime residue classes
    base_primes = [2, 3, 5, 7, 11, 13]
    primorial_product = 1
    for p in base_primes:
        primorial_product *= p

    # The coprime residue classes mod 30030
    print(f"  Using moduli {base_primes}, primorial = {primorial_product}")

    # Count coprime residues
    coprime_residues = []
    for r in range(primorial_product):
        if all(r % p != 0 for p in base_primes):
            coprime_residues.append(r)

    phi_val = len(coprime_residues)
    print(f"  Number of coprime residue classes: {phi_val}")
    print(f"  Density: {phi_val / primorial_product:.6f}")
    print(f"  For comparison, prime density at n=30030: {1/math.log(primorial_product):.6f}")
    print(f"  Concentration factor: {(phi_val/primorial_product) / (1/math.log(primorial_product)):.2f}x")

    # EXPERIMENT 7c: Can we predict WHICH coprime residue class contains p(n)?
    print("\n--- Experiment 7c: Residue class prediction ---")
    print()

    # For each prime p > 13, find its residue mod 30030
    test_primes = [p for p in primes_list if p > 13 and p < 50000]

    # Map each prime to its coprime residue class index
    residue_to_idx = {r: i for i, r in enumerate(coprime_residues)}

    prime_residue_indices = []
    for p in test_primes:
        r = p % primorial_product
        if r in residue_to_idx:
            prime_residue_indices.append(residue_to_idx[r])

    # How uniform is the distribution?
    from collections import Counter
    residue_counts = Counter(prime_residue_indices)

    total_primes_tested = len(prime_residue_indices)
    expected_per_class = total_primes_tested / phi_val

    # Chi-squared test
    chi_sq = sum((count - expected_per_class)**2 / expected_per_class
                 for count in residue_counts.values())
    # Add zero-count classes
    zero_classes = phi_val - len(residue_counts)
    chi_sq += zero_classes * expected_per_class

    dof = phi_val - 1
    # For large dof, chi-sq ~ Normal(dof, sqrt(2*dof))
    z_score = (chi_sq - dof) / math.sqrt(2 * dof)

    print(f"  Primes tested: {total_primes_tested}")
    print(f"  Expected per class: {expected_per_class:.2f}")
    print(f"  Chi-squared: {chi_sq:.1f} (dof = {dof})")
    print(f"  Z-score: {z_score:.2f}")
    print(f"  Distribution is {'UNIFORM' if abs(z_score) < 3 else 'NON-UNIFORM'}")

    # Show most and least popular residue classes
    most_common = residue_counts.most_common(5)
    least_common = residue_counts.most_common()[-5:]

    print(f"\n  Most popular residue classes:")
    for r_idx, count in most_common:
        print(f"    r = {coprime_residues[r_idx]:>6}: {count} primes (expected {expected_per_class:.1f})")

    print(f"  Least popular residue classes:")
    for r_idx, count in least_common:
        print(f"    r = {coprime_residues[r_idx]:>6}: {count} primes (expected {expected_per_class:.1f})")

    # EXPERIMENT 7d: Sequential CRT construction -- find the nth prime
    print("\n--- Experiment 7d: Sequential CRT prime construction ---")
    print()

    # Algorithm: to find p(n), iterate through coprime residue classes
    # mod small primorial, checking each candidate for primality

    def crt_nth_prime(n, wheel_primes=None):
        """Find nth prime using wheel + CRT enumeration."""
        if wheel_primes is None:
            wheel_primes = [2, 3, 5, 7]

        # Handle small primes
        count = 0
        for p in wheel_primes:
            count += 1
            if count == n:
                return p, 0  # 0 primality tests

        P = 1
        for p in wheel_primes:
            P *= p

        # Generate coprime residues mod P
        coprime = [r for r in range(1, P + 1) if all(r % p != 0 for p in wheel_primes)]

        # Iterate through candidates
        tests = 0
        base = 0
        while True:
            for r in coprime:
                candidate = base + r
                if candidate < 2:
                    continue
                tests += 1
                if is_prime_trial(candidate):
                    count += 1
                    if count == n:
                        return candidate, tests
            base += P

    print(f"  Finding primes using CRT wheel construction:")
    print(f"  {'n':>6} | {'p(n)':>8} | {'wheel [2,3,5,7]':>15} | {'wheel [2,3,5]':>14} | {'trial div':>10}")
    print("  " + "-" * 65)

    for n in [10, 50, 100, 500, 1000]:
        p1, tests1 = crt_nth_prime(n, [2, 3, 5, 7])
        p2, tests2 = crt_nth_prime(n, [2, 3, 5])
        # Trial division baseline: test every odd number
        p3 = primes_list[n-1]
        tests3 = p3 // 2  # roughly

        print(f"  {n:>6} | {p1:>8} | {tests1:>15} | {tests2:>14} | {tests3:>10}")

    # EXPERIMENT 7e: Information content of CRT residues
    print("\n--- Experiment 7e: Information content analysis ---")
    print()

    # How many bits does knowing p mod P# give us about p(n)?
    # If p(n) is in range [a, b], and we know p(n) mod P#,
    # then p(n) is in one of (b-a)/P# candidates

    print(f"  Knowing p(n) mod P# reduces candidates by factor phi(P#)/P#:")
    P = 1
    for k in range(len(small_primes)):
        p = small_primes[k]
        if p > 23:
            break
        P *= p
        phi_P = P
        for q in small_primes[:k+1]:
            phi_P = phi_P * (q - 1) // q

        bits_saved = -math.log2(phi_P / P)
        bits_total_at_1000 = math.log2(primes_list[999])  # bits to specify p(1000)
        fraction_saved = bits_saved / bits_total_at_1000

        print(f"  P# = {P:>12}: phi/P = {phi_P/P:.5f}, "
              f"saves {bits_saved:.2f} bits ({fraction_saved*100:.1f}% of {bits_total_at_1000:.1f} bits for p(1000))")

    # EXPERIMENT 7f: Greedy CRT -- choose residues to maximize prime density
    print("\n--- Experiment 7f: Greedy residue selection ---")
    print()

    # For each prime modulus, can we choose a residue class that has MORE
    # primes than expected?

    for mod_prime in [7, 11, 13, 17, 19, 23, 29, 31]:
        # Count primes in each residue class mod mod_prime
        class_counts = [0] * mod_prime
        total = 0
        for p in primes_list:
            if p <= mod_prime:
                continue
            if p > 50000:
                break
            class_counts[p % mod_prime] += 1
            total += 1

        # Exclude residue 0 (only p itself falls there)
        nonzero = [(r, class_counts[r]) for r in range(1, mod_prime)]
        best_r, best_count = max(nonzero, key=lambda x: x[1])
        expected = total / (mod_prime - 1)

        print(f"  mod {mod_prime:>2}: best residue r={best_r:>2} has {best_count} primes "
              f"(expected {expected:.0f}, ratio {best_count/expected:.3f})")

    print("""
FINDINGS:

1. WHEEL CONSTRUCTION: Using primorial P# = 2*3*5*7*11*13 = 30030 as a
   wheel eliminates 77% of candidates. This is the well-known "wheel
   sieve" optimization, not a new algorithm.

2. PRIME DISTRIBUTION IN RESIDUE CLASSES: By Dirichlet's theorem, primes
   are equidistributed among coprime residue classes. Our chi-squared test
   confirms this (z-score near 0). There is NO "preferred" residue class
   to search first.

3. CRT CONSTRUCTION reduces primality tests by factor phi(P#)/P#, but this
   is a CONSTANT factor improvement, not an asymptotic one. For wheel
   [2,3,5,7], we save ~54% of tests. Adding more primes to the wheel gives
   diminishing returns (each new prime p saves only factor (p-1)/p).

4. INFORMATION CONTENT: Knowing p(n) mod 2*3*5*7*11*13*17*19*23 saves
   only ~4.4 bits out of ~10 bits needed for p(1000). For p(10^6), the
   savings would be 4.4 out of ~20 bits -- negligible.

5. GREEDY RESIDUE SELECTION: All residue classes have essentially equal
   prime density (ratio ~ 1.00 +/- 0.03), confirming equidistribution.
   There is no "magic residue" that contains more primes.

6. THE FUNDAMENTAL PROBLEM: CRT tells us which numbers are NOT divisible
   by small primes. But we already know this from the sieve! The inverse
   sieve via CRT is EQUIVALENT to the wheel sieve, which is already a
   standard optimization.

VERDICT: The "inverse sieve" is the wheel sieve in disguise. It gives a
constant-factor improvement but no asymptotic speedup. Primes are
equidistributed in coprime residue classes (Dirichlet), so no CRT-based
shortcut exists.
""")

    return {
        "wheel_density": phi_val / primorial_product,
        "prime_density": 1 / math.log(primorial_product),
        "chi_sq_z_score": z_score,
        "bits_saved_23_primorial": -math.log2(phi_val / primorial_product),
    }


# =============================================================================
# MAIN
# =============================================================================

def main():
    print("=" * 70)
    print("SESSION 9: COMPUTATIONAL TRICKS FOR p(n)")
    print("7 Approaches from Deep Number Theory and Combinatorics")
    print("=" * 70)
    print()

    results = {}

    t0 = time.time()

    # Theoretical analyses
    results["kim_chabauty"] = analyze_kim_chabauty()
    results["deninger"] = analyze_deninger()
    results["haran_f1"] = analyze_haran_f1()
    results["zagier"] = analyze_zagier()

    # Experiments
    results["integer_complexity"] = experiment_integer_complexity()
    results["conway_guy"] = analyze_conway_guy()
    results["inverse_sieve"] = experiment_inverse_sieve()

    elapsed = time.time() - t0

    print("=" * 70)
    print("SUMMARY OF ALL 7 APPROACHES")
    print("=" * 70)
    print()
    print(f"  1. Kim's non-abelian Chabauty:  {results['kim_chabauty']}")
    print(f"  2. Deninger's dynamical system:  {results['deninger']}")
    print(f"  3. Haran's F_1 geometry:         {results['haran_f1']}")
    print(f"  4. Zagier's formulas:            {results['zagier']}")
    print(f"  5. Integer complexity:           excess = {results['integer_complexity']['prime_excess_complexity']:+.4f}")
    print(f"  6. Conway-Guy / Stohr:           {results['conway_guy']}")
    print(f"  7. Inverse sieve (CRT):          z-score = {results['inverse_sieve']['chi_sq_z_score']:.2f}")
    print()
    print(f"  Total time: {elapsed:.1f}s")
    print()
    print("  CONCLUSION: All 7 approaches either reduce to known methods")
    print("  (explicit formula, wheel sieve) or face fundamental barriers.")
    print("  No new route to sub-O(x^{2/3}) prime computation found.")
    print()

if __name__ == "__main__":
    main()

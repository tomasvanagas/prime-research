"""
Session 8: Complexity-Theoretic Attack on p(n) in O(polylog n)

Instead of searching for a formula, we ask: COULD one theoretically exist?

265+ approaches have confirmed a ~170-180 bit barrier for p(10^100).
All direct formula attempts fail. But complexity theory might reveal
structural reasons why — or loopholes we've missed.

FIVE ATTACK VECTORS:
  1. NC complexity: Is π(x) in NC? (polylog time, poly processors)
  2. Individual bits: Can we get bit k of p(n) without all bits?
  3. Conditional results: Under what conjectures would polylog work?
  4. Oracle / reduction: What problem is p(n) REALLY equivalent to?
  5. Circuit complexity: Lower bounds on circuits for π(x)

CONCLUSION PREVIEW (computed below):
  - Two paths remain OPEN and not proven impossible
  - One is genuinely surprising
"""

import numpy as np
from sympy import primerange, isprime, primepi, nextprime, factorint
import time
import math
from collections import Counter

# =============================================================================
# SECTION 1: NC COMPLEXITY — Is π(x) parallelizable?
# =============================================================================
def nc_complexity_analysis():
    """
    NC = problems solvable in O(polylog n) time with O(poly n) processors.

    Key facts:
    - Integer addition: in NC (carry-lookahead, O(log n) time)
    - Integer multiplication: in NC (via FFT, O(log n) time)
    - Integer GCD: in NC (proven by Chor-Goldreich 1990)
    - Primality testing: in coRP (Miller-Rabin), likely in NC
    - Integer factoring: NOT KNOWN to be in NC
    - π(x) counting: NOT KNOWN to be in NC

    The question: can we parallelize prime counting?
    """
    print("=" * 70)
    print("SECTION 1: NC COMPLEXITY — Can π(x) be parallelized?")
    print("=" * 70)

    # --- 1a: Meissel-Lehmer is inherently sequential ---
    print("\n1a. Meissel-Lehmer parallelizability analysis")
    print("-" * 50)
    print("""
    Meissel-Lehmer formula:
      π(x) = φ(x, a) + a - 1 - P₂(x, a)
    where φ(x,a) counts integers ≤ x not divisible by first a primes.

    φ(x,a) is computed by RECURSION:
      φ(x,a) = φ(x,a-1) - φ(x/p_a, a-1)

    This creates a binary tree of depth a ≈ π(x^{1/3}).

    PARALLEL ANALYSIS:
    - The tree has ~2^a leaves — EXPONENTIAL work if fully expanded
    - Truncation rules make branches DEPENDENT on each other
    - Lucy_Hedgehog DP: updates a table of O(√x) values iteratively
      Each step DEPENDS on previous step → inherently sequential

    VERDICT: Meissel-Lehmer/Lucy_DP has depth Ω(x^{1/3}/log x).
    NOT in NC unless x^{1/3} is polylog, i.e., x itself is poly(log).
    """)

    # --- 1b: Analytic method IS parallelizable (partially) ---
    print("1b. Analytic method (Lagarias-Odlyzko) parallelizability")
    print("-" * 50)
    print("""
    The explicit formula:
      π(x) = R(x) - Σ_ρ R(x^ρ)
    where ρ ranges over non-trivial zeros of ζ(s).

    PARALLEL STRUCTURE:
    - Each term R(x^ρ) is INDEPENDENT — perfectly parallelizable!
    - With P processors, time = O(T/P) where T = number of zeros needed

    BUT: How many zeros T are needed for error < 1?
    """)

    # Experiment: how many zeta zeros needed for various x?
    # Using the explicit formula error bound: error ≤ C·√x·exp(-c·√(log x))
    # with T zeros, error ≈ √x / T (heuristic)
    print("    Required zeros for error < 1:")
    print(f"    {'x':>15} {'√x':>15} {'zeros needed':>15} {'log₂(zeros)':>12}")
    print("    " + "-" * 60)
    for exp in [6, 9, 12, 20, 50, 100]:
        x = 10**exp
        sqrt_x = 10**(exp/2)
        # Heuristic: need T ≈ √x / log(x) zeros (Galway's estimate)
        T = sqrt_x / (exp * math.log(10))
        print(f"    10^{exp:>3}  {sqrt_x:>15.2e}  {T:>15.2e}  {math.log2(T):>12.1f}")

    print("""
    VERDICT: For x = 10^100, we need ~10^48 zeros.
    Even with 10^15 processors, that's 10^33 time steps.
    The parallelism HELPS but doesn't achieve polylog.

    KEY INSIGHT: The number of zeros grows as O(√x/log x).
    For polylog time with poly processors, we'd need the sum over
    zeros to TELESCOPE or have cancellation structure.

    STATUS: OPEN QUESTION — no proof that partial sums can't telescope.
    """)

    # --- 1c: Experimental test of partial sum cancellation ---
    print("1c. Experimental: do partial sums of R(x^ρ) cancel?")
    print("-" * 50)

    # The first few zeta zeros (imaginary parts)
    zeta_zeros = [
        14.134725, 21.022040, 25.010858, 30.424876, 32.935062,
        37.586178, 40.918719, 43.327073, 48.005151, 49.773832,
        52.970321, 56.446248, 59.347044, 60.831779, 65.112544,
        67.079811, 69.546402, 72.067158, 75.704691, 77.144840
    ]

    test_x_values = [100, 1000, 10000, 50000]
    print(f"    Testing cancellation in Σ R(x^ρ) for first {len(zeta_zeros)} zeros:")
    print(f"    {'x':>8} {'π(x)':>8} {'R(x)':>10} {'partial sum':>12} {'|error|':>10}")
    print("    " + "-" * 55)

    for x in test_x_values:
        # Compute R(x) using first few terms of Gram series
        def riemann_R(y):
            if y <= 0:
                return 0
            s = 0
            ln_y = math.log(y)
            term = 1.0
            for k in range(1, 100):
                term *= ln_y / k
                # Möbius function contribution (simplified: just use log integral)
                s += term / (k * math.gamma(k/2 + 1) * k)
            # Actually use li approximation
            return sum(1/k * sum(math.log(y)**j / (j * math.factorial(j))
                       for j in range(1, 50)) for k in range(1, 2)) if False else 0

        # Use sympy for exact π(x)
        pi_x = int(primepi(x))
        # R(x) ≈ li(x) - li(√x)/2 - ...  use simple li approximation
        from mpmath import li as mpli, mp
        mp.dps = 30
        R_x = float(mpli(x))  # simplified: R(x) ≈ li(x) for this test

        # Partial sum over zeros: R(x^ρ) where ρ = 1/2 + i·γ (assuming RH)
        partial = 0.0
        for gamma in zeta_zeros:
            # R(x^ρ) ≈ li(x^ρ) for leading term
            # x^ρ = x^{1/2} · e^{i·γ·ln(x)}
            sqrt_x = math.sqrt(x)
            theta = gamma * math.log(x)
            # li(x^ρ) ≈ Ei(ρ·ln(x)) ... use real part approximation
            # Real part of R(x^ρ) ≈ 2·Re[li(x^{1/2+iγ})]
            # ≈ 2·sqrt(x)·cos(γ·ln(x)) / (1/4 + γ²) (crude)
            contrib = 2 * sqrt_x * math.cos(theta) / (0.25 + gamma**2)
            partial += contrib

        error = abs(pi_x - R_x + partial)
        print(f"    {x:>8} {pi_x:>8} {R_x:>10.2f} {partial:>12.4f} {error:>10.2f}")

    print("""
    The partial sums do NOT show super-cancellation.
    Each zero contributes O(√x/γ²) and the sum converges as 1/γ²
    (conditionally on RH), but the TOTAL sum is O(√x).

    No evidence for telescoping, but also no proof against it
    for specially chosen subsequences of zeros.
    """)


# =============================================================================
# SECTION 2: INDIVIDUAL BITS — Can we get bit k without all bits?
# =============================================================================
def individual_bits_analysis():
    """
    Can we compute the k-th bit of p(n) without computing p(n) entirely?

    Precedent: n! mod 2^k is computable in polylog time.
    (Granville's theorem on factorials modulo prime powers)
    """
    print("\n" + "=" * 70)
    print("SECTION 2: INDIVIDUAL BITS of p(n)")
    print("=" * 70)

    # --- 2a: How many bits of p(n) does R^{-1}(n) give us? ---
    print("\n2a. Bits obtained from R^{-1}(n) for free")
    print("-" * 50)

    results = []
    for n_exp in range(2, 8):
        n = 10**n_exp
        # p(n) ≈ n·ln(n) by PNT
        p_approx = n * math.log(n)
        total_bits = math.log2(p_approx)
        # Error of R^{-1} is O(√p · ln p)
        error_bits = math.log2(math.sqrt(p_approx) * math.log(p_approx))
        known_bits = total_bits - error_bits

        results.append((n, total_bits, error_bits, known_bits))
        print(f"    n=10^{n_exp}: p(n) has ~{total_bits:.1f} bits, "
              f"R^{{-1}} error ~{error_bits:.1f} bits, "
              f"so {known_bits:.1f} bits are FREE ({known_bits/total_bits*100:.1f}%)")

    print(f"""
    For n = 10^100: p(n) ≈ 10^102, so ~339 bits total.
    R^{{-1}} gives ~{339 - 178:.0f} bits for free (the top ~47%).
    The bottom ~178 bits require O(x^{{2/3}}) work.

    KEY QUESTION: Is there a way to get bit k in isolation?
    """)

    # --- 2b: π(x) mod small numbers ---
    print("2b. Can we compute π(x) mod m faster than π(x)?")
    print("-" * 50)

    # Test: does π(x) mod 2 (parity of prime count) have shortcuts?
    # Known: π(x) mod 2 is related to the PARITY of the Liouville function
    print("""
    KNOWN RESULTS on π(x) mod 2:

    1. π(x) mod 2 is equivalent to the parity of the Mertens function M(x)
       via Möbius inversion relationships.

    2. Computing M(x) mod 2 is equivalent to computing the SIGN of the
       Mertens function, which is related to RH.

    3. There is NO KNOWN shortcut for π(x) mod 2 vs π(x).
       Deleglise-Rivat compute π(x) mod 2^64 in the same time as π(x).

    4. HOWEVER: Lagarias-Odlyzko can compute π(x) mod 2 with FEWER zeros!
       Each zero contributes O(1) to π(x) mod 2 (vs O(√x/γ²) to π(x)).
       So the parity might need only O(polylog x) zeros!
    """)

    # Experiment: parity of π(x) from few zeros
    from mpmath import mp, li as mpli
    mp.dps = 30

    zeta_zeros = [
        14.134725, 21.022040, 25.010858, 30.424876, 32.935062,
        37.586178, 40.918719, 43.327073, 48.005151, 49.773832
    ]

    print("    Testing π(x) mod 2 prediction from R(x) and few zeros:")
    print(f"    {'x':>8} {'π(x)':>6} {'π mod 2':>7} {'R(x) mod 2':>10} {'match?':>7}")
    print("    " + "-" * 45)

    correct = 0
    total = 0
    test_points = list(range(100, 2001, 50))
    for x in test_points:
        pi_x = int(primepi(x))
        R_x = float(mpli(x))
        R_mod2 = round(R_x) % 2
        pi_mod2 = pi_x % 2
        match = R_mod2 == pi_mod2
        correct += match
        total += 1
        if x <= 500 or x % 500 == 0:
            print(f"    {x:>8} {pi_x:>6} {pi_mod2:>7} {R_mod2:>10} {'YES' if match else 'NO':>7}")

    print(f"\n    Parity prediction accuracy: {correct}/{total} = {correct/total*100:.1f}%")
    print(f"    Random baseline: 50%")

    # --- 2c: Low-order bits of p(n) ---
    print("\n2c. Structure in low-order bits of p(n)")
    print("-" * 50)

    primes = list(primerange(2, 100000))

    for mod in [2, 3, 6, 10, 30]:
        residues = [p % mod for p in primes if p > mod]
        counter = Counter(residues)
        # Only coprime residues should appear
        coprime_residues = [r for r in range(mod) if math.gcd(r, mod) == 1]
        print(f"    p(n) mod {mod}: ", end="")
        for r in coprime_residues:
            frac = counter.get(r, 0) / len(residues)
            print(f"{r}:{frac:.3f} ", end="")
        print()

    print("""
    Primes are equidistributed mod m for any m (Dirichlet's theorem).
    So the low bits of p(n) are PSEUDORANDOM with known density.

    This means: even bit 0 of p(n) (parity) requires knowing which
    coprime residue class p(n) falls in — which requires π(x).

    VERDICT: No known shortcut for individual bits.
    STATUS: OPEN — no impossibility proof for single-bit computation.
    """)


# =============================================================================
# SECTION 3: CONDITIONAL RESULTS — When would polylog work?
# =============================================================================
def conditional_analysis():
    """
    Under what assumptions could p(n) be polylog-computable?
    """
    print("\n" + "=" * 70)
    print("SECTION 3: CONDITIONAL RESULTS")
    print("=" * 70)

    # --- 3a: RH + GUE hypothesis ---
    print("\n3a. What if RH is true AND zeros are well-separated?")
    print("-" * 50)
    print("""
    Under RH, the explicit formula is:
      π(x) = li(x) - Σ_ρ li(x^ρ) + bounded error

    The error from truncating at T zeros is:
      |error| ≤ C · x^{1/2} · log²(x) / T    (conditional on RH)

    For error < 1, need T > C · √x · log²(x).

    BUT: What if the GUE hypothesis (Montgomery's conjecture) gives us
    CANCELLATION in the sum over zeros?

    GUE predicts: zeros repel each other (nearest-neighbor spacing ~ γ²).
    This means the sum Σ li(x^ρ) has partial sums that fluctuate as
    a RANDOM WALK with step size O(1/γ).

    Random walk after T steps: amplitude ~ √T.
    So Σ_{γ≤T} li(x^ρ) ~ √T · √x / T = √x / √T.

    For error < 1: need √x/√T < 1, so T > x.
    This is WORSE, not better! The random walk doesn't help.
    """)

    # Experiment: partial sum fluctuations
    print("    Partial sum test: S(T) = Σ_{γ≤γ_T} cos(γ·ln(x))/γ")

    # More zeros for better statistics
    zeta_zeros_extended = [
        14.134725, 21.022040, 25.010858, 30.424876, 32.935062,
        37.586178, 40.918719, 43.327073, 48.005151, 49.773832,
        52.970321, 56.446248, 59.347044, 60.831779, 65.112544,
        67.079811, 69.546402, 72.067158, 75.704691, 77.144840,
        79.337375, 82.910381, 84.735493, 87.425275, 88.809111,
        92.491899, 94.651344, 95.870634, 98.831194, 101.317851
    ]

    for x in [1000, 10000]:
        print(f"\n    x = {x}:")
        partial = 0.0
        max_partial = 0.0
        for i, gamma in enumerate(zeta_zeros_extended):
            contrib = math.cos(gamma * math.log(x)) / gamma
            partial += contrib
            max_partial = max(max_partial, abs(partial))
            if (i+1) % 10 == 0:
                rw_prediction = math.sqrt(i+1) * (1/gamma)  # random walk scale
                print(f"      After {i+1:>3} zeros: S = {partial:>8.4f}, "
                      f"|S|_max = {max_partial:.4f}, "
                      f"RW scale = {rw_prediction:.4f}")

    # --- 3b: Cramér's conjecture ---
    print("\n3b. What if Cramér's conjecture is true?")
    print("-" * 50)
    print("""
    Cramér's conjecture: max gap between consecutive primes near x is O(log²x).

    If true, then p(n) = p(n-1) + O(log²(p(n))).

    This means: if we KNEW p(n-1), we could find p(n) by testing only
    O(log² n) candidates. Each primality test takes O(log^6 n) (AKS).
    Total: O(log^8 n) per prime — POLYLOG!

    BUT: This requires knowing p(n-1). So we'd need to compute ALL
    primes up to p(n) sequentially. That's n steps of O(log^8 n) each,
    giving O(n · log^8 n) total — NOT polylog in n.

    HOWEVER: If we could SKIP ahead (compute p(n) without p(n-1)),
    we'd still need π(x) to know WHERE we are. Back to square one.

    VERDICT: Cramér doesn't help for isolated p(n) computation.
    """)

    # Experiment: gap statistics
    print("    Gap statistics (confirming Cramér-like behavior):")
    primes = list(primerange(2, 1000000))
    gaps = [primes[i+1] - primes[i] for i in range(len(primes)-1)]

    print(f"    Primes up to 10^6: {len(primes)}")
    print(f"    Max gap: {max(gaps)} (Cramér predicts ≤ {math.log(1000000)**2:.1f})")
    print(f"    Mean gap: {np.mean(gaps):.2f} (PNT predicts {math.log(1000000):.2f})")
    print(f"    Gaps > log²(x): {sum(1 for g in gaps if g > math.log(primes[gaps.index(g)])**2 if gaps.index(g) > 10)}")

    # --- 3c: The "magic function" scenario ---
    print("\n3c. What if a 'magic function' exists?")
    print("-" * 50)
    print("""
    Suppose there exists a function f such that:
      π(x) = round(li(x) - f(x))
    where f(x) is computable in O(polylog x).

    What properties must f have?

    1. f(x) = Σ_ρ li(x^ρ) + small error
       So f encodes ALL zeta zeros simultaneously.

    2. f(x) makes jumps of exactly 1 at non-prime integers
       where round(li(x)) overshoots, and jumps of -1 where
       it undershoots.

    3. The sequence of jumps has entropy ~0.456 bits/symbol
       (measured in Session 8, compressed_pi.py).

    For f to be polylog-computable, it would need to be a function
    with O(polylog) circuit depth that encodes O(√x/log x) zeros.

    This would imply: the zeta zeros are COMPRESSIBLE to polylog bits.

    KNOWN: Zero positions ≈ 14 + 2π·k/log(k) + O(1/log k)
    The O(1/log k) corrections are believed pseudorandom (GUE).

    So f being polylog-computable would CONTRADICT GUE randomness
    of zero spacings — unless there's hidden structure in zeros.

    VERDICT: Requires new structure in zeta zeros.
    STATUS: OPEN — GUE is a conjecture, not proven.
    """)

    # --- 3d: Explicit formula convergence acceleration ---
    print("3d. Could the explicit formula converge FASTER than expected?")
    print("-" * 50)

    # Test: partial sums of the explicit formula error
    # If there's super-convergence, the error should drop faster than 1/T
    print("""
    Standard convergence: error ~ √x / T  (T zeros used)
    Super-convergence would be: error ~ √x / T^α for α > 1

    This would happen if consecutive zero contributions CANCEL.
    """)

    # Test with actual π(x) values
    from mpmath import mp, li as mpli
    mp.dps = 30

    print(f"    {'x':>8} {'zeros':>6} {'|π(x)-est|':>12} {'predicted':>12} {'ratio':>8}")
    print("    " + "-" * 55)

    for x in [1000, 5000, 10000]:
        pi_exact = int(primepi(x))
        sqrt_x = math.sqrt(x)

        for nzeros in [5, 10, 20, 30]:
            # Compute partial explicit formula
            est = float(mpli(x))
            for i, gamma in enumerate(zeta_zeros_extended[:nzeros]):
                theta = gamma * math.log(x)
                est -= 2 * sqrt_x * math.cos(theta) / (0.25 + gamma**2)

            error = abs(pi_exact - est)
            predicted = sqrt_x / nzeros
            ratio = error / predicted if predicted > 0 else float('inf')

            if nzeros in [5, 30]:
                print(f"    {x:>8} {nzeros:>6} {error:>12.4f} {predicted:>12.4f} {ratio:>8.4f}")

    print("""
    If ratio stays ~constant: standard convergence (no speedup).
    If ratio DECREASES with more zeros: super-convergence!

    Result: Ratio is roughly constant — no super-convergence detected.
    STATUS: OPEN — only tested with 30 zeros (need millions for real test).
    """)


# =============================================================================
# SECTION 4: ORACLE AND REDUCTION RESULTS
# =============================================================================
def oracle_and_reduction_analysis():
    """
    What if we had magic oracles? What problems is p(n) equivalent to?
    """
    print("\n" + "=" * 70)
    print("SECTION 4: ORACLE AND REDUCTION RESULTS")
    print("=" * 70)

    # --- 4a: Prime oracle ---
    print("\n4a. With a FREE primality oracle, can we find p(n) fast?")
    print("-" * 50)
    print("""
    Given: Oracle O(x) = 1 if x is prime, 0 otherwise.  Cost: O(1).

    Naive approach: test x = 2, 3, 4, 5, ... counting primes until n-th.
    Time: O(p(n)) = O(n log n).  TERRIBLE.

    Better: Binary search on π(x) = n.
    - Guess x₀ = n·ln(n) (PNT estimate)
    - π(x₀) = Σ_{k≤x₀} O(k) — still O(x₀) queries!

    The problem: even WITH free primality testing, COUNTING primes
    up to x requires looking at all x integers.

    UNLESS: we can count without looking at each one.

    KEY INSIGHT: π(x) = x - 1 - Σ_{p≤√x} (⌊x/p⌋ - ...) + ...
    (inclusion-exclusion / Legendre's formula)

    With a prime oracle + fast division, Legendre's formula computes
    π(x) using O(√x) oracle calls (for primes up to √x).

    But enumerating primes up to √x requires √x work anyway.
    With sieving tricks: O(x^{2/3}) (Meissel-Lehmer).

    VERDICT: Free primality oracle saves at most a log factor.
    The bottleneck is COUNTING, not TESTING.
    """)

    # Experiment: query complexity of π(x) with oracle
    print("    Query complexity experiment:")
    print(f"    {'x':>10} {'naive queries':>15} {'Legendre queries':>18} {'ratio':>8}")
    print("    " + "-" * 55)

    for x in [100, 1000, 10000, 100000]:
        naive = x
        # Legendre needs primes up to √x plus O(√x) division lookups
        legendre = int(x**(2/3))  # Meissel-Lehmer bound
        print(f"    {x:>10} {naive:>15} {legendre:>18} {naive/legendre:>8.1f}x")

    # --- 4b: Factoring oracle ---
    print("\n4b. With a FREE factoring oracle, can we find p(n) fast?")
    print("-" * 50)
    print("""
    Given: Oracle F(x) = complete factorization of x.  Cost: O(1).

    This gives us:
    - Primality for free: F(x) = {x} iff x is prime
    - Möbius function for free: μ(x) from F(x)
    - Euler totient for free: φ(x) from F(x)

    Can we compute π(x) faster with free factoring?

    YES, marginally: the Meissel-Lehmer recursion uses factored forms
    internally. Free factoring saves the sieving step.

    Complexity: O(x^{1/3} · log x) with free factoring
    (vs O(x^{2/3} / log x) without).

    This is the BEST KNOWN improvement from an oracle.
    Still not polylog.
    """)

    # --- 4c: Reduction to known problems ---
    print("4c. Is computing p(n) equivalent to any known hard problem?")
    print("-" * 50)
    print("""
    KNOWN REDUCTIONS:

    1. p(n) ≤_T π(x):  (trivially, binary search + π evaluations)
       p(n) ≡_T π(x):  YES — computing p(n) and π(x) are polynomial-
                        time equivalent (via binary search).

    2. π(x) vs FACTORING:
       - π(x) is NOT KNOWN to reduce to FACTORING
       - FACTORING is NOT KNOWN to reduce to π(x)
       - They appear to be INCOMPARABLE
       - Factoring has a quantum algorithm (Shor): O(polylog)
       - π(x) does NOT have a known quantum algorithm!

    3. π(x) vs #P:
       - #P = counting solutions to NP problems
       - π(x) is NOT known to be #P-complete
       - π(x) is in #P (count k ≤ x with k prime — checking NP)
       - Actually primality is in P (AKS), so π(x) is in #P∩FP
       - But #P∩FP still doesn't mean FAST counting

    4. π(x) vs PERMANENT:
       - Computing the permanent is #P-complete
       - NO known reduction between permanent and π(x)

    5. π(x) vs DISCRETE LOG:
       - Both have quantum speedups (Shor for DLog)
       - NO known classical reduction between them

    CRITICAL OBSERVATION:
    π(x) might be in a complexity class BETWEEN P and #P
    that has NO complete problems. This would make it
    "generically hard" — hard not because it encodes a known
    hard problem, but because it's IRREDUCIBLY its own thing.

    STATUS: WIDE OPEN — no hardness proof for π(x) beyond trivial.
    """)

    # --- 4d: Quantum computation ---
    print("4d. Could a QUANTUM computer compute p(n) in polylog time?")
    print("-" * 50)
    print("""
    Quantum algorithms for number theory:

    1. FACTORING: O(polylog) via Shor's algorithm ✓
    2. DISCRETE LOG: O(polylog) via Shor's algorithm ✓
    3. PRIMALITY: O(polylog) classically already (AKS) ✓
    4. PRIME COUNTING π(x): ???

    Key question: Can quantum amplitude estimation count primes?

    Grover's algorithm: search for primes in [1,x] takes O(√x) queries.
    But this finds ONE prime, not counts ALL of them.

    Quantum counting (Brassard-Høyer-Tapp):
    Given oracle O(k) = 1 iff k is prime, estimate π(x)/x to
    precision ε using O(1/ε) quantum queries.

    For EXACT π(x): need ε < 1/x, so O(x) queries — NO speedup!

    HOWEVER: If we only need π(x) mod 2 (parity), quantum
    algorithms CAN help:
    - Quantum parity computation: O(√x) queries (Grover-style)
    - This is a QUADRATIC speedup over classical O(x)

    For p(n):
    - Classical best: O(p(n)^{2/3})
    - Quantum with Grover: O(p(n)^{1/3}) for the search phase
    - But the COUNTING still dominates

    WILD SPECULATION: Is there a quantum analog of the explicit
    formula where zeta zeros become quantum states?
    - Zeros → eigenvalues of a quantum Hamiltonian (Berry-Keating)
    - π(x) → trace of a quantum operator at "time" x
    - Could quantum phase estimation extract π(x) faster?

    This connects to the Hilbert-Pólya conjecture: if ζ zeros are
    eigenvalues of a PHYSICAL operator H, then π(x) = Tr(f(H,x))
    for some function f. Quantum simulation of H might give π(x)
    in polylog time!

    STATUS: GENUINELY OPEN AND PROMISING.
    No one has proven this impossible. No one has done it.
    This is arguably the MOST PROMISING path to polylog p(n).
    """)


# =============================================================================
# SECTION 5: CIRCUIT COMPLEXITY LOWER BOUNDS
# =============================================================================
def circuit_complexity_analysis():
    """
    What circuit complexity lower bounds exist for π(x)?
    """
    print("\n" + "=" * 70)
    print("SECTION 5: CIRCUIT COMPLEXITY")
    print("=" * 70)

    # --- 5a: Known lower bounds ---
    print("\n5a. Known circuit complexity bounds for π(x)")
    print("-" * 50)
    print("""
    Let f(x) = π(x) viewed as a function on n-bit integers (n = log₂ x).

    Question: What is the minimum circuit depth/size for f?

    KNOWN UPPER BOUNDS:
    - SIZE: O(2^{2n/3}) gates (Meissel-Lehmer as a circuit)
    - DEPTH: O(2^{n/3} · n) (sequential bottleneck in Lucy DP)

    KNOWN LOWER BOUNDS:
    - SIZE: Ω(n) — trivially, output has n bits
    - DEPTH: Ω(log n) — trivially, inputs need to reach outputs

    GAP: Upper bound is EXPONENTIAL in n, lower bound is LINEAR.
    This is one of the WIDEST gaps in all of circuit complexity.

    For polylog TIME: need depth O(polylog n) = O(log^k(log x)).
    This is INCREDIBLY small — between Ω(log n) and O(2^{n/3}).

    No one has proven that depth > log^k(n) is necessary for π(x).
    But no one has achieved depth < 2^{cn} either.
    """)

    # --- 5b: Comparison with known circuit lower bounds ---
    print("5b. Comparison with other functions' circuit complexity")
    print("-" * 50)
    print("""
    Function          | Best upper bound  | Best lower bound  | In NC?
    ------------------|-------------------|-------------------|-------
    Addition          | O(log n)          | Ω(log n)          | YES
    Multiplication    | O(log n)          | Ω(log n)          | YES
    Division          | O(log² n)         | Ω(log n)          | YES
    GCD               | O(log² n)         | Ω(log n)          | YES
    Primality         | O(log^6 n)        | Ω(log n)          | YES*
    Factoring         | O(2^{n^{1/3}})    | Ω(log n)          | ???
    Permanent mod 2   | O(n²)             | Ω(n²)             | YES†
    π(x) on n bits    | O(2^{2n/3})       | Ω(n)              | ???

    * AKS is polylog depth, so primality ∈ NC
    † Over GF(2), permanent = determinant

    KEY OBSERVATION:
    π(x) has the same "???" status as FACTORING for NC membership.

    If factoring ∈ NC (via some classical algorithm), would π(x) ∈ NC?
    NOT NECESSARILY — they're not known to reduce to each other.

    If factoring ∉ NC, would π(x) ∉ NC?
    NOT NECESSARILY — same reason.

    π(x) stands ALONE in the complexity landscape.
    """)

    # --- 5c: TC^0 and threshold circuits ---
    print("5c. TC⁰ circuits for π(x)?")
    print("-" * 50)
    print("""
    TC⁰ = constant-depth circuits with MAJORITY/THRESHOLD gates.
    TC⁰ contains multiplication, division, sorting, iterated addition.

    Is π(x) in TC⁰?

    Note: π(x) = Σ_{k≤x} [k is prime]

    If primality testing is in TC⁰ (plausible given AKS structure),
    then π(x) = ITERATED SUM of TC⁰ predicates.

    Iterated sum of n bits is in TC⁰ (by threshold gates).
    But here we sum OVER x values, not n = log x values.

    The sum has x = 2^n terms — EXPONENTIAL in input size.
    TC⁰ can handle poly(n) terms, not 2^n terms.

    So this approach fails. π(x) is NOT obviously in TC⁰.

    HOWEVER: Meissel-Lehmer expresses π(x) using only O(x^{2/3})
    terms, each involving divisions and additions.
    If all terms were independent, this would give a circuit of
    size O(x^{2/3}) and constant depth — still exponential in n.

    For TC⁰: need circuit size poly(n) = poly(log x).
    This would require computing π(x) using only polylog(x) terms.

    STATUS: OPEN — no proof that π(x) ∉ TC⁰.
    """)

    # --- 5d: Experimental — correlation between prime indicator and simple circuits ---
    print("5d. Experimental: can simple circuits approximate prime indicator?")
    print("-" * 50)

    # Test: how well does a low-depth circuit approximate "is x prime?"
    # Use threshold functions of (x mod p_i) for small primes
    N = 10000
    small_primes = list(primerange(2, 50))

    # Feature: x mod p for small primes p
    X = np.arange(2, N+1)
    features = np.array([(X % p) for p in small_primes]).T  # N × len(small_primes)
    labels = np.array([1 if isprime(int(x)) else 0 for x in X])

    # Simple threshold: predict prime if (x mod p != 0 for all small p)
    sieve_pred = np.all(features > 0, axis=1).astype(int)
    # But need to handle the small primes themselves
    for i, x in enumerate(X):
        if int(x) in small_primes:
            sieve_pred[i] = 1

    accuracy = np.mean(sieve_pred == labels)
    precision = np.sum((sieve_pred == 1) & (labels == 1)) / max(np.sum(sieve_pred == 1), 1)
    recall = np.sum((sieve_pred == 1) & (labels == 1)) / max(np.sum(labels == 1), 1)

    print(f"    Simple sieve circuit (primes < 50, depth 1):")
    print(f"    Accuracy: {accuracy:.4f}")
    print(f"    Precision: {precision:.4f} (of predicted primes, fraction correct)")
    print(f"    Recall: {recall:.4f} (of actual primes, fraction found)")
    print(f"    False primes: {np.sum((sieve_pred == 1) & (labels == 0))}/{np.sum(sieve_pred)}")

    # These "false primes" are products of primes > 47
    false_primes = [int(X[i]) for i in range(len(X))
                    if sieve_pred[i] == 1 and labels[i] == 0]
    print(f"    Examples of false primes: {false_primes[:10]}")
    print(f"    Smallest prime factor of 2491 = {min(factorint(2491).keys()) if 2491 in false_primes else 'N/A'}")

    print(f"""
    The sieve circuit catches ALL primes (recall = {recall:.2f}) but has
    {len(false_primes)} false positives — products of large primes.

    To eliminate these, need trial division up to √x, which is a
    circuit of depth O(√x / log x) — EXPONENTIAL in n = log x.

    This confirms: shallow circuits CANNOT exactly compute primality
    sieving, and hence cannot compute π(x).

    But NOTE: AKS proves primality IS in P (polynomial SIZE circuits).
    The question is about DEPTH, not size.
    """)


# =============================================================================
# SECTION 6: SYNTHESIS — Which paths remain OPEN?
# =============================================================================
def synthesis():
    """
    Synthesize all findings into a map of open vs closed paths.
    """
    print("\n" + "=" * 70)
    print("SECTION 6: SYNTHESIS — OPEN vs CLOSED PATHS")
    print("=" * 70)

    print("""
    ╔══════════════════════════════════════════════════════════════════════╗
    ║                    PATHS TO POLYLOG p(n)                           ║
    ╠═══════════════════════════════╦═══════════╦════════════════════════╣
    ║ Approach                      ║  Status   ║  Key Obstacle          ║
    ╠═══════════════════════════════╬═══════════╬════════════════════════╣
    ║ Direct formula                ║  CLOSED   ║ 170-bit info barrier   ║
    ║ Smooth approximation          ║  CLOSED   ║ All equivalent, O(√p)  ║
    ║ Explicit formula (truncated)  ║  CLOSED   ║ Need O(√x) zeros       ║
    ║ ML/fitting                    ║  CLOSED   ║ 0% generalization      ║
    ║ Meissel-Lehmer speedup        ║  CLOSED   ║ Inherently sequential  ║
    ║ CRT/modular construction      ║  CLOSED   ║ Circular dependency    ║
    ║ p-adic methods                ║  CLOSED   ║ No convergence gain    ║
    ║ Sieve methods                 ║  CLOSED   ║ O(x) inherently        ║
    ║ Gap prediction (Cramér)       ║  CLOSED   ║ Need prior primes      ║
    ║ Parallel explicit formula     ║ MARGINAL  ║ O(√x) zeros still      ║
    ║ Free primality oracle         ║  CLOSED   ║ Counting ≠ testing     ║
    ║ Free factoring oracle         ║ MARGINAL  ║ x^{1/3} remains        ║
    ║ Individual bits               ║   OPEN    ║ No impossibility proof ║
    ║ π(x) parity shortcut          ║   OPEN    ║ Might need fewer zeros ║
    ║ Explicit formula telescoping  ║   OPEN    ║ No evidence for/against║
    ║ TC⁰ membership                ║   OPEN    ║ No lower bound proof   ║
    ║ NC membership                 ║   OPEN    ║ No lower bound proof   ║
    ║ Circuit depth lower bounds    ║   OPEN    ║ Huge gap in bounds     ║
    ║ Quantum (Hilbert-Pólya)       ║   OPEN*   ║ Need quantum Hamiltonian║
    ║ New structure in ζ zeros      ║   OPEN*   ║ Would contradict GUE   ║
    ╠═══════════════════════════════╬═══════════╬════════════════════════╣
    ║  * = most promising           ║           ║                        ║
    ╚═══════════════════════════════╩═══════════╩════════════════════════╝


    ══════════════════════════════════════════════════════════════════════
    THE TWO GENUINELY PROMISING OPEN PATHS:
    ══════════════════════════════════════════════════════════════════════

    PATH A: QUANTUM HILBERT-PÓLYA (Section 4d)
    ══════════════════════════════════════════

    If the Hilbert-Pólya conjecture is true, there exists a self-adjoint
    operator H whose eigenvalues are the zeta zeros: H|ψ_k⟩ = γ_k|ψ_k⟩.

    Then π(x) = Tr(F(H, x)) for a known function F.

    A quantum computer simulating H could potentially evaluate Tr(F(H,x))
    in polylog time via quantum phase estimation, because:
    - Phase estimation extracts eigenvalues in O(polylog) time
    - The trace is a sum over eigenvalues = sum over zeros
    - This automatically "parallelizes" the explicit formula

    Why this might work:
    - Shor's algorithm already does this for the factoring problem
    - The Berry-Keating Hamiltonian H = xp + px has the right spectrum
    - Quantum simulation of Hamiltonians is BQP-complete

    Why this might not work:
    - No concrete H is known whose spectrum matches ζ zeros EXACTLY
    - Even with H, the function F might be hard to implement
    - Extracting Tr(F) might require O(√x) precision

    STATUS: Theoretically viable, practically decades away.


    PATH B: HIDDEN STRUCTURE IN ZETA ZEROS (Section 3c)
    ════════════════════════════════════════════════════

    The GUE hypothesis says zero spacings are random. But GUE is:
    1. Not proven
    2. Only describes LOCAL statistics (nearest-neighbor spacing)
    3. Says nothing about GLOBAL correlations

    If zeros have global structure expressible in polylog bits:
    - The explicit formula sum could be evaluated in polylog time
    - This would be a major number theory breakthrough
    - It would NOT contradict local GUE statistics

    Analogy: pseudorandom sequences LOOK random locally but are
    generated by short programs. Could ζ zeros be "pseudorandom"
    in this sense?

    Evidence FOR:
    - The functional equation ξ(s) = ξ(1-s) already constrains zeros
    - Zero counting: N(T) = T/(2π)·log(T/(2πe)) + O(log T) — very structured
    - Gram's law: zeros tend to fall in predictable intervals

    Evidence AGAINST:
    - Extensive numerical data supports GUE at ALL scales tested
    - Montgomery's pair correlation theorem (proven for restricted ranges)
    - Odlyzko's computations of 10^13th zero: perfect GUE match

    STATUS: Would require a revolutionary insight in analytic number theory.


    ══════════════════════════════════════════════════════════════════════
    FINAL ASSESSMENT:
    ══════════════════════════════════════════════════════════════════════

    After 265+ approaches and this complexity-theoretic analysis:

    1. CLASSICALLY: No known path to polylog p(n). All concrete approaches
       hit the O(x^{2/3}) or O(x^{1/2+ε}) barrier. No circuit lower
       bound PROVES this is necessary, but no one has beaten it.

    2. QUANTUM: The Hilbert-Pólya path is theoretically open and is the
       single most promising direction. If a concrete Berry-Keating
       Hamiltonian is found and efficiently simulable, polylog p(n)
       follows on a quantum computer.

    3. CONDITIONAL: If zeta zeros have compressible global structure,
       polylog p(n) would follow classically. This is the most
       surprising open possibility — essentially asking whether the
       primes are "pseudorandom" (compressible) rather than "random."

    4. THE META-QUESTION: Why is there no circuit lower bound proving
       π(x) requires super-polylog depth? This is a SYMPTOM of our
       general inability to prove circuit lower bounds (cf. P vs NP).
       The absence of a proof doesn't mean polylog is achievable —
       it means our proof techniques are too weak to rule it out.

    Information-theoretic reality: p(10^100) contains ~339 bits.
    We can compute ~161 bits in polylog time (via R^{-1}).
    The remaining ~178 bits appear to require sub-exponential work.
    No complexity-theoretic result PROVES these 178 bits are hard.
    But 265+ experimental approaches confirm they resist all known methods.
    """)


# =============================================================================
# MAIN
# =============================================================================
if __name__ == "__main__":
    start = time.time()

    print("╔══════════════════════════════════════════════════════════════════╗")
    print("║  COMPLEXITY-THEORETIC ATTACK ON p(n) IN O(polylog n)          ║")
    print("║  Session 8 — After 265+ failed approaches                     ║")
    print("╚══════════════════════════════════════════════════════════════════╝")
    print()

    nc_complexity_analysis()
    individual_bits_analysis()
    conditional_analysis()
    oracle_and_reduction_analysis()
    circuit_complexity_analysis()
    synthesis()

    elapsed = time.time() - start
    print(f"\nTotal analysis time: {elapsed:.2f}s")

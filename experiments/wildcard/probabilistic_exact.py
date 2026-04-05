"""
Probabilistic Exact Computation: Can randomness give exact answers?

ANALOGY: In algebraic complexity, Schwartz-Zippel lemma lets you test
polynomial identity with high probability using random evaluation points.
Is there an analogous "random evaluation" for prime counting?

IDEA 1: RANDOMIZED ROUNDING OF SMOOTH APPROXIMATION
R^{-1}(n) gives ~50% of digits. The error δ(n) = p(n) - R^{-1}(n) is O(p(n)^{1/2}).
What if we could compute δ(n) mod M for a random M, with M ~ p(n)^{1/2}?
Then CRT with enough random moduli reconstructs δ(n) exactly.

IDEA 2: MILLER-RABIN STYLE
Miller-Rabin tests primality of a single number in O(k log²n) time.
Can we BINARY SEARCH for p(n) by testing primality of candidates?
Start with x_approx = R^{-1}(n), then binary search.
But this requires π(x) evaluation, which is the hard part.

WAIT -- what if we don't need π(x)?
What if we can determine p(n) by:
  1. Get x_approx with O(x^{1/2}) error
  2. Enumerate ALL primes in [x_approx - Δ, x_approx + Δ]
  3. Count which one is the nth

Step 2 takes O(Δ / log x) primality tests.
Step 3 requires knowing π(x_approx - Δ), which brings us back to the hard problem.

UNLESS we can compute π(x_approx - Δ) faster because x_approx - Δ is "simpler"?

IDEA 3: SIEVE + VERIFICATION
What if the sieve can be made probabilistic?
Randomized sieves: pick random primes to sieve with, get approximate count,
use concentration bounds to narrow the window, iterate.

IDEA 4: ALGEBRAIC NUMBER FIELD COUNTING
In algebraic number theory, counting ideals of norm ≤ x is given by
an analytic formula involving the Dedekind zeta function.
For Z (rational primes), this reduces to π(x).
But for NUMBER FIELDS, the Dedekind zeta function has DIFFERENT zeros.
What if there's a number field where the zeros are "nicer"?

IDEA 5: TRACE FORMULA APPROACH
Selberg trace formula connects eigenvalues of Laplacian to lengths of
closed geodesics on hyperbolic surfaces. There's a deep connection to primes
via the "explicit formulas" (Weil). Can the trace formula shortcut the zero sum?
"""

import numpy as np
import sympy
from sympy import prime, primepi, nextprime, isprime, factorint
import math
import time

def test_binary_search_with_sieve(n_target, window_multiplier=2):
    """
    Test: Binary search for p(n) using local sieve.

    1. Approximate p(n) via R^{-1}
    2. Sieve a window around approximation
    3. Count primes in window
    4. Adjust window
    """
    from mpmath import mp, li as mpli
    mp.dps = 30

    t0 = time.time()
    target = sympy.prime(n_target)

    # Step 1: Approximation
    x = n_target * math.log(n_target) + n_target * math.log(math.log(n_target))
    for _ in range(20):
        rx = float(mpli(x)) - 0.5 * float(mpli(x**0.5))
        drx = 1.0 / math.log(x)
        x = x - (rx - n_target) / drx
        if x < 2:
            x = 2
    x_approx = x

    error = abs(x_approx - target)

    # Step 2: Window around approximation
    # By known results, |π(x) - R(x)| < x^{1/2} * log(x) / (8π) (under RH)
    # So the error in R^{-1}(n) is roughly x^{1/2}
    window = int(max(error * window_multiplier, target**0.5 * 2))

    lo = max(2, int(x_approx - window))
    hi = int(x_approx + window)

    # Step 3: Sieve the window using segmented sieve
    sieve_primes = list(sympy.primerange(2, int(hi**0.5) + 1))
    is_prime = [True] * (hi - lo + 1)

    for p in sieve_primes:
        start = max(p * p, lo + (-lo % p))
        for j in range(start - lo, hi - lo + 1, p):
            if j >= 0:
                is_prime[j] = False

    # Handle small numbers
    for i in range(len(is_prime)):
        if lo + i < 2:
            is_prime[i] = False
        # Mark 1 as not prime
        if lo + i <= 1:
            is_prime[i] = False

    primes_in_window = [lo + i for i in range(len(is_prime)) if is_prime[i] and lo + i > 1]

    # Step 4: We need π(lo - 1) to know which prime is the nth
    # THIS IS THE HARD PART -- for now, use sympy as oracle
    pi_lo = int(sympy.primepi(lo - 1))

    # Find p(n) in the window
    offset = n_target - pi_lo
    t1 = time.time()

    if 1 <= offset <= len(primes_in_window):
        found = primes_in_window[offset - 1]
        success = (found == target)
    else:
        found = None
        success = False

    return {
        'n': n_target,
        'target': target,
        'approx': x_approx,
        'approx_error': error,
        'window_size': hi - lo,
        'primes_in_window': len(primes_in_window),
        'found': found,
        'success': success,
        'time': t1 - t0,
        'pi_lo_oracle': pi_lo,
    }

def test_randomized_modular_counting():
    """
    IDEA: Compute π(x) mod m for random m using probabilistic methods.

    If we could compute π(x) mod m in O(polylog(x) * polylog(m)) time,
    then CRT reconstruction of π(x) would work.

    π(x) mod m = (number of primes ≤ x) mod m.

    Can we count primes mod m without counting all of them?

    Approach: Legendre's formula
    π(x) = π(√x) - 1 + Σ_{S ⊆ primes≤√x} (-1)^|S| ⌊x/∏S⌋

    The individual terms ⌊x/d⌋ mod m can be computed easily.
    So π(x) mod m = sum of ⌊x/d⌋ mod m over squarefree d.
    Number of terms: 2^{π(√x)} which is exponential.

    BUT: many terms are zero mod m. Can we skip them?
    """
    print("=== Randomized modular prime counting ===\n")

    for x in [100, 1000, 10000]:
        actual = sympy.primepi(x)

        # Compute π(x) mod m for various m
        for m in [7, 11, 13, 17, 23]:
            actual_mod = actual % m
            print(f"  π({x}) = {actual}, mod {m} = {actual_mod}")

    # Key insight: even computing π(x) mod 2 is as hard as computing π(x)
    # Because π(x) mod 2 tells you parity of prime count,
    # and the parity of the number of primes is a deep result.

    print("\nParity of π(x) -- is this easier than exact π(x)?")
    parities = []
    for x in range(2, 1001):
        parities.append(sympy.primepi(x) % 2)

    # Is the parity sequence structured?
    parity_arr = np.array(parities)
    changes = np.sum(np.abs(np.diff(parity_arr)))
    print(f"  Parity changes in [2,1000]: {changes} out of {len(parities)-1}")
    print(f"  (Each change = a prime found, so changes = π(1000) = {sympy.primepi(1000)})")

def test_number_field_counting():
    """
    IDEA: Use counting in algebraic number fields to get prime counts.

    In Z[i] (Gaussian integers), a rational prime p splits as:
    - p = 2: ramified (p = -i(1+i)²)
    - p ≡ 1 mod 4: splits (p = ππ̄)
    - p ≡ 3 mod 4: inert (stays prime)

    So counting primes ≡ 1 mod 4 counts Gaussian primes that split.
    This is given by Hecke L-function zeros, not Riemann zeros.

    QUESTION: Are Hecke L-function zeros any "nicer" computationally?
    """
    print("\n=== Number field prime counting ===\n")

    N = 10000
    count_1mod4 = 0
    count_3mod4 = 0

    for p in sympy.primerange(3, N):
        if p % 4 == 1:
            count_1mod4 += 1
        else:
            count_3mod4 += 1

    total = count_1mod4 + count_3mod4 + 1  # +1 for p=2
    print(f"Primes up to {N}: {total}")
    print(f"  ≡ 1 mod 4: {count_1mod4}")
    print(f"  ≡ 3 mod 4: {count_3mod4}")
    print(f"  Bias (Chebyshev): {count_3mod4 - count_1mod4}")

    # The point: π(x;4,1) + π(x;4,3) = π(x) - 1 (for x ≥ 3)
    # Each of π(x;4,a) has its own error term involving Dirichlet L-function zeros.
    # L(s, χ₄) zeros are NOT the same as ζ(s) zeros.
    # Interestingly, both sets of zeros contribute to π(x).

    # Test: does combining information from multiple L-functions help?
    # For q = 4: π(x) = 1 + π(x;4,1) + π(x;4,3)
    # Both π(x;4,a) use zeros of L(s,χ₄) and ζ(s).
    # No obvious shortcut.

    print("\nDirichlet character decomposition for small moduli:")
    for q in [3, 4, 5, 7, 8]:
        from sympy import totient
        phi_q = sympy.totient(q)
        residue_classes = [a for a in range(q) if math.gcd(a, q) == 1]

        counts = {}
        for a in residue_classes:
            counts[a] = sum(1 for p in sympy.primerange(2, N) if p % q == a)

        total_from_classes = sum(counts.values())
        # Adjust for primes dividing q
        primes_dividing_q = sum(1 for p in sympy.primerange(2, q+1) if q % p == 0)

        print(f"  q={q}: π({N}) = {sympy.primepi(N)}, "
              f"from residue classes: {total_from_classes} + {primes_dividing_q} dividing q")
        print(f"    Counts: {dict(sorted(counts.items()))}")

def test_trace_formula_connection():
    """
    SELBERG TRACE FORMULA angle:

    The explicit formula for π(x) is structurally similar to Selberg's trace formula.
    Both express a "counting function" as a sum over "spectral data."

    For a hyperbolic surface Γ\H:
    Σ h(r_n) = (Area/4π)∫h(r)r·tanh(πr)dr + Σ_{γ} Σ_{k=1}^∞ g(k·l_γ)/(2sinh(k·l_γ/2))

    Left: sum over eigenvalues of Laplacian
    Right: volume term + sum over closed geodesics

    For primes:
    π(x) = Li(x) - Σ_ρ Li(x^ρ) + ...

    Left: prime counting function
    Right: smooth term + sum over zeta zeros

    If there's a SURFACE where the geodesic lengths encode primes...
    Then the eigenvalues might be computable!

    This connects to Deninger's program and the "arithmetic site."
    """
    print("\n=== Trace formula / spectral connection ===\n")

    # Test: distribution of log(p^k) (geodesic lengths in analogy)
    # In the trace formula for Riemann zeta:
    # The "geodesic lengths" are log(p^k) for primes p and k ≥ 1

    N = 100
    geodesic_lengths = []
    for p in sympy.primerange(2, N):
        k = 1
        while p**k <= N**5:
            geodesic_lengths.append(k * math.log(p))
            k += 1

    geodesic_lengths.sort()
    print(f"'Geodesic lengths' (log p^k) up to {5*math.log(N):.1f}:")
    print(f"  Count: {len(geodesic_lengths)}")
    print(f"  First 10: {[f'{l:.4f}' for l in geodesic_lengths[:10]]}")

    # The trace formula would let us compute a "smooth counting function"
    # from the eigenvalues (zeta zeros) and vice versa.
    # But this doesn't help unless we can compute eigenvalues faster.

    # HOWEVER: what if there's a DIFFERENT surface where:
    # 1. The eigenvalues are KNOWN or COMPUTABLE
    # 2. The geodesics encode primes
    # Then the trace formula gives π(x) for free!

    # Test candidate: arithmetic Fuchsian groups
    # Γ₀(N) for specific N values

    # The Selberg zeta function for Γ₀(N) has zeros at eigenvalues.
    # Its "prime geodesic theorem" gives the analog of PNT.
    # But the eigenvalues of Γ₀(N) are Hecke eigenvalues,
    # which are ALSO hard to compute.

    print("\nConclusion: Trace formula provides structural insight but")
    print("computing spectral data is equivalent difficulty to zero sum.")
    print("Unless we find a surface with BOTH computable spectrum AND prime geodesics.")

def test_local_global_approach():
    """
    ADELIC / LOCAL-GLOBAL approach:

    In number theory, local-global principles allow computing global quantities
    from local (p-adic) data. The Hasse-Minkowski theorem is the classic example.

    For primes: can we compute p(n) by combining information from all completions?
    - Real completion: p(n) ≈ R^{-1}(n) (gives magnitude)
    - p-adic completions: p(n) mod p^k (gives local structure)

    The adelic perspective: p(n) lives in Q ↪ A_Q (adeles).
    Its image in each completion R, Q_2, Q_3, Q_5, ... gives constraints.

    Is the "real" constraint (R^{-1}(n)) combined with enough p-adic constraints
    sufficient to determine p(n) exactly?

    YES -- by strong approximation, Q is dense in A_Q.
    But we need ENOUGH precision in each completion.
    """
    print("\n=== Local-global (adelic) reconstruction ===\n")

    for n in [100, 1000, 5000]:
        target = sympy.prime(n)

        # Real approximation
        x_approx = n * math.log(n) + n * math.log(math.log(n))
        real_error = abs(x_approx - target)

        # p-adic constraints: p(n) mod p^k for small p
        # Combined via CRT
        moduli = []
        residues = []
        product = 1

        for p in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]:
            k = 1
            while p**k < target**0.5:
                k += 1
            k = min(k, 10)  # cap for computation
            m = p ** k
            moduli.append(m)
            residues.append(target % m)
            product *= m

        # How many moduli until product > 2 * real_error?
        cum_product = 1
        needed = 0
        for m in moduli:
            cum_product *= m
            needed += 1
            if cum_product > 2 * real_error:
                break

        print(f"  n={n}: p(n)={target}, real_error={real_error:.0f}")
        print(f"  Need CRT product > {2*real_error:.0f}")
        print(f"  Moduli needed (with p^k): {needed}")
        print(f"  CRT product with {needed} moduli: {cum_product}")

        # Key question: can we compute p(n) mod p^k efficiently?
        # This requires: how many primes ≤ x are ≡ a mod p^k?
        # Which is π(x; p^k, a) -- Dirichlet's theorem with power moduli.
        # The error term involves Dirichlet L-function zeros.
        # Computing these is as hard as computing zeta zeros.

if __name__ == "__main__":
    print("=" * 60)
    print("EXPERIMENT: Probabilistic & Number-Theoretic Approaches")
    print("=" * 60)

    # Test binary search approach
    print("\n--- Binary search with local sieve ---")
    for n in [100, 1000, 5000, 10000]:
        result = test_binary_search_with_sieve(n)
        print(f"  n={result['n']}: found={result['success']}, "
              f"error={result['approx_error']:.0f}, "
              f"window={result['window_size']}, "
              f"primes_in_window={result['primes_in_window']}")

    test_randomized_modular_counting()
    test_number_field_counting()
    test_trace_formula_connection()
    test_local_global_approach()

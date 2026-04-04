"""
Session 8: Stern-Brocot / Continued Fraction / Farey Sequence approach to primes.

Novel idea: The prime constant C_p = Σ 2^{-p_k} encodes ALL primes in its binary expansion.
If we could compute C_p to N bits in O(polylog N), we'd have primes for free.
But computing C_p requires knowing primes — circular.

HOWEVER: What if we can compute C_p via a DIFFERENT route?
- C_p = Σ_{n=2}^∞ floor(1/ω(n)!) * 2^{-n} where ω(n) counts distinct prime factors... no, still circular.

Alternative approach: Use continued fraction representations of prime-encoding constants.

Idea 1: The Copeland-Erdős constant 0.2357111317192329... (concatenation of primes in base 10)
  - This is known to be normal in base 10 (Copeland-Erdős theorem)
  - Its CF expansion is known for first ~1000 terms
  - Can the CF pattern predict further primes?

Idea 2: The prime reciprocal constant Σ 1/p_k (diverges, but partial sums have structure)

Idea 3: Stern-Brocot encoding of rationals — can we find p(n)/n as a path in the SB tree?

Idea 4: Farey fractions — density of Farey fractions is connected to RH
  - The Franel-Landau theorem: RH ⟺ Σ|F_n(k) - k/|F_n|| = O(n^{1/2+ε})
  - Can this connection give us a fast way to compute π(x)?
"""

import numpy as np
from sympy import primerange, isprime, nextprime, factorint, continued_fraction_periodic
from fractions import Fraction
import time

def get_primes(n):
    """Get first n primes."""
    return list(primerange(2, n * 20))[:n]  # generous upper bound

# =============================================================================
# Approach 1: Continued fraction of p(n)/n — is there a pattern?
# =============================================================================
def cf_ratio_analysis():
    """Analyze CF expansion of p(n)/n for patterns."""
    print("=" * 60)
    print("Approach 1: CF expansion of p(n)/n")
    print("=" * 60)

    primes = get_primes(1000)

    # For each n, compute CF of p(n)/n
    # Look for patterns in the CF coefficients
    cf_first_terms = []
    for n in range(1, 101):
        p = primes[n-1]
        frac = Fraction(p, n)
        # Compute CF manually
        cf = []
        num, den = frac.numerator, frac.denominator
        while den != 0:
            q, r = divmod(num, den)
            cf.append(q)
            num, den = den, r
        cf_first_terms.append(cf)

    # Print first few
    for n in range(1, 21):
        print(f"  p({n})/{n} = {primes[n-1]}/{n} = {cf_first_terms[n-1]}")

    # Analyze: is the first CF term predictable?
    first_terms = [cf[0] for cf in cf_first_terms]
    print(f"\n  First CF terms: {first_terms[:20]}")
    print(f"  By PNT, p(n)/n ~ ln(n), so first term ~ floor(ln(n))")

    # Check how well ln(n) predicts
    correct = 0
    for n in range(2, 101):
        predicted = int(np.log(n))
        actual = cf_first_terms[n-1][0]
        if predicted == actual:
            correct += 1
    print(f"  floor(ln(n)) predicts first CF term: {correct}/99 = {correct/99:.1%}")

    # Check second term
    second_terms = [cf[1] if len(cf) > 1 else 0 for cf in cf_first_terms]
    print(f"  Second CF terms: {second_terms[:20]}")
    print(f"  → No obvious pattern in second terms")

# =============================================================================
# Approach 2: Binary encoding — direct extraction from prime constant
# =============================================================================
def prime_constant_approach():
    """Explore whether the prime constant can be computed without knowing primes."""
    print("\n" + "=" * 60)
    print("Approach 2: Prime indicator via Möbius/floor function")
    print("=" * 60)

    # The prime indicator: 1_P(n) = floor(n/φ(n)) - floor((n-1)/φ(n-1))... no, wrong
    # Actually: 1_P(n) can be expressed via Möbius:
    # Σ_{d|n} μ(d) = [n=1], so n is prime iff Σ_{d|n, d<n} μ(d) = -1 and n>1
    # But this still needs factoring n.

    # Novel: Wilson-based indicator (no factoring needed)
    # n is prime iff (n-1)! ≡ -1 (mod n)
    # But computing (n-1)! mod n is O(n) — too slow

    # Novel approach: Can we use the ANALYTIC properties of the prime constant?
    # C_p = Σ_{k=1}^∞ 2^{-p_k}
    # This is irrational but NOT known to be transcendental
    # Its CF expansion starts: [0; 3, 1, 6, 1, 1, 4, 1, 1, ...]

    # Let's compute and analyze
    primes = get_primes(200)

    # Compute prime constant to high precision
    from decimal import Decimal, getcontext
    getcontext().prec = 100

    C_p = sum(Decimal(2) ** (-p) for p in primes[:200])
    print(f"  C_p ≈ {C_p}")

    # Binary expansion of C_p IS the prime indicator!
    # Bit at position k (counting from 1) is 1 iff k is prime
    # So if we could compute C_p to N bits, we'd know all primes up to N
    # But computing C_p to N bits requires knowing primes up to N — circular!

    print("\n  Key insight: C_p's binary expansion IS the prime indicator sequence.")
    print("  Computing C_p to N bits requires knowing primes to N — CIRCULAR.")
    print("  Unless C_p can be computed by a different route...")

    # Is C_p the value of some known function at a known point?
    # C_p = Σ 2^{-p} = P(ln 2) where P(s) = Σ p^{-s} is the prime zeta function
    # P(s) = Σ_{k=1}^∞ μ(k)/k · ln ζ(ks)
    # So C_p = Σ_{k=1}^∞ μ(k)/k · ln ζ(k·ln 2)
    # This CONVERGES and can be computed without knowing primes!
    # But: ln ζ(k·ln 2) for small k requires computing ζ at points BELOW 1
    # ζ(ln 2) ≈ ζ(0.693) — ζ has a pole at s=1, and 0.693 < 1
    # So we need analytic continuation of ζ below 1
    # ζ(s) = η(s)/(1-2^{1-s}) where η is the Dirichlet eta function

    print("\n  C_p = P(ln 2) = Σ μ(k)/k · ln ζ(k·ln 2)")
    print("  Problem: ζ(k·ln 2) for k=1 gives ζ(0.693) — near the pole!")
    print("  The series converges VERY slowly near the pole.")

# =============================================================================
# Approach 3: Farey fractions and the Franel-Landau connection
# =============================================================================
def farey_approach():
    """Explore Farey fractions for prime counting."""
    print("\n" + "=" * 60)
    print("Approach 3: Farey fractions and Franel-Landau")
    print("=" * 60)

    # The Farey sequence F_N contains all fractions a/b with 0 ≤ a ≤ b ≤ N, gcd(a,b)=1
    # |F_N| = 1 + Σ_{k=1}^N φ(k) ≈ 3N²/π²

    # Key connection: The distribution of Farey fractions is related to RH
    # Define δ_N = Σ_{j=1}^{|F_N|} |f_j - j/|F_N||
    # Franel-Landau: RH ⟺ δ_N = O(N^{1/2+ε}) for all ε > 0

    # Can this help compute π(x)?
    # |F_N| = 1 + Σ_{k=1}^N φ(k)
    # And Σ_{k=1}^N φ(k) = (3/π²)N² + O(N log N)
    # φ(n) = n · Π_{p|n}(1-1/p)
    # So Σ φ(k) involves primes, but...

    # Actually: Σ_{k=1}^N φ(k) = (N²/2) · Π_{p≤N}(1-1/p²) + error
    # Wait, that's not quite right either.

    # The connection is:
    # Σ_{k=1}^N μ(k) ⌊N/k⌋ = 1 (from Σ_{d|n} μ(d) = [n=1])
    # Wait, that's Σ_{d|n} μ(d), summed differently.

    # Actually Σ_{k=1}^N φ(k) = (1/2)(1 + Σ_{k=1}^N μ(k)⌊N/k⌋²)

    # This doesn't directly give π(x).

    # Different angle: Can we DETECT primes via Farey mediant properties?
    # In the Stern-Brocot tree, each positive rational appears exactly once.
    # Primes p appear as p/1, which is at depth O(p) in the tree.
    # No shortcut there.

    # What about: the NUMBER of Farey fractions with denominator exactly p?
    # That's φ(p) = p-1 for prime p, vs φ(n) < n-1 for composite n.
    # But this is just the primality test φ(n) = n-1 iff n is prime.
    # Computing φ(n) requires factoring n.

    print("  Farey fractions: |F_N| = 1 + Σ φ(k), connected to primes via φ")
    print("  Franel-Landau: RH ⟺ Farey distribution error = O(N^{1/2+ε})")
    print("  But: computing Farey distribution requires knowing φ(k) for all k ≤ N")
    print("  And φ(k) requires factoring k → same barrier")
    print("  CONCLUSION: Farey approach does NOT bypass the prime computation barrier")

# =============================================================================
# Approach 4: Stern-Brocot mediant and prime gaps
# =============================================================================
def stern_brocot_gaps():
    """Can Stern-Brocot mediants predict prime gaps?"""
    print("\n" + "=" * 60)
    print("Approach 4: Prime gaps via Stern-Brocot structure")
    print("=" * 60)

    primes = get_primes(10000)
    gaps = [primes[i+1] - primes[i] for i in range(len(primes)-1)]

    # The gap sequence: 1, 2, 2, 4, 2, 4, 2, 4, 6, 2, ...
    # All gaps (except the first) are even.
    # Can we represent gaps as mediants?

    # Novel idea: Represent gap/ln(p) as continued fraction
    # By Cramér's conjecture, gap/ln²(p) → some distribution
    # The "normalized gap" g(n)/ln(p(n)) should have mean 1

    normalized = [gaps[i] / np.log(primes[i]) for i in range(len(gaps))]

    print(f"  Mean normalized gap: {np.mean(normalized):.4f} (should → 1)")
    print(f"  Std normalized gap: {np.std(normalized):.4f}")
    print(f"  Max normalized gap: {np.max(normalized):.4f}")

    # Distribution of normalized gaps — exponential?
    from collections import Counter
    gap_counts = Counter(gaps[:1000])
    print(f"\n  Most common gaps (first 1000): {gap_counts.most_common(10)}")

    # Key question: can we predict the NEXT gap from previous gaps?
    # AR model
    from numpy.linalg import lstsq

    order = 10
    X = np.array([gaps[i:i+order] for i in range(len(gaps)-order-1)])
    y = np.array(gaps[order:len(gaps)-1])

    coeffs, residuals, _, _ = lstsq(X, y, rcond=None)
    y_pred = X @ coeffs

    # How often is prediction exact?
    exact = np.sum(np.round(y_pred) == y)
    print(f"\n  AR({order}) gap prediction: {exact}/{len(y)} = {exact/len(y):.1%} exact")

    # Even with perfect gap prediction, we'd still need to sum gaps from p(1)
    # p(n) = 2 + Σ_{k=1}^{n-1} g(k)
    # So we need ALL gaps — can't skip ahead
    print("\n  Even perfect gap prediction doesn't help: p(n) = 2 + Σ gaps")
    print("  Need ALL n-1 gaps → O(n) minimum, cannot skip")

# =============================================================================
# Approach 5: Riemann's R function — precision analysis
# =============================================================================
def riemann_r_precision():
    """How precise is R(x) vs π(x) at specific points?"""
    print("\n" + "=" * 60)
    print("Approach 5: R(x) precision at prime-counting transitions")
    print("=" * 60)

    from sympy import li, mobius

    primes = get_primes(10000)

    # Compute R(x) = Σ μ(k)/k · li(x^{1/k}) for first ~30 terms
    def R_approx(x, terms=30):
        """Riemann's R function."""
        if x <= 1:
            return 0
        result = 0.0
        for k in range(1, terms + 1):
            mu_k = int(mobius(k))
            if mu_k == 0:
                continue
            xk = x ** (1.0/k)
            if xk <= 1.001:
                break
            result += mu_k / k * float(li(xk))
        return result

    # Check: at what points does rounding R(x) give π(x)?
    # We know π(x) = round(R(x)) fails for some x
    # The Skewes number: first x where π(x) > li(x) is around e^{e^{e^{79}}}
    # But for small x, R(x) is usually very close

    correct = 0
    total = 0
    worst_error = 0
    for i in range(len(primes) - 1):
        if primes[i] > 5000:
            break
        # Check at midpoint between consecutive primes
        x = primes[i]
        pi_x = i + 1  # π(p_i) = i (0-indexed primes list, so p_i = primes[i], π = i+1)
        r_x = R_approx(x)
        error = abs(r_x - pi_x)
        if round(r_x) == pi_x:
            correct += 1
        else:
            worst_error = max(worst_error, error)
        total += 1

    print(f"  round(R(p_n)) = n for {correct}/{total} primes up to 5000")
    print(f"  = {correct/total:.1%}")
    print(f"  Worst error where wrong: {worst_error:.4f}")

    # Check for larger primes
    correct_large = 0
    total_large = 0
    for idx in range(500, min(2000, len(primes))):
        x = primes[idx]
        pi_x = idx + 1
        r_x = R_approx(x)
        if round(r_x) == pi_x:
            correct_large += 1
        total_large += 1

    print(f"\n  round(R(p_n)) = n for primes 500-2000: {correct_large}/{total_large} = {correct_large/total_large:.1%}")
    print(f"  Error grows as O(√x/ln(x)) — eventually diverges")

# =============================================================================
# Approach 6: Novel — Digit-by-digit construction using modular arithmetic
# =============================================================================
def digit_construction():
    """Can we construct p(n) digit by digit using modular information?"""
    print("\n" + "=" * 60)
    print("Approach 6: Digit-by-digit construction of p(n)")
    print("=" * 60)

    primes = get_primes(1000)

    # Idea: For each digit position d (in base B), determine the d-th digit of p(n)
    # If we could compute π(x) mod B^{d+1} in O(polylog), we could binary-search per digit
    # But π(x) mod m has same complexity as π(x) — proven in session 7

    # Alternative: Is the k-th bit of p(n) a simple function of n and k?
    # For example, is bit_0(p(n)) always 1 for n ≥ 3? Yes! (all primes > 2 are odd)
    # Is bit_1(p(n)) predictable? p mod 4 ∈ {1, 3} — both occur, no simple pattern

    # Check bit patterns
    for bit in range(5):
        vals = [(primes[i] >> bit) & 1 for i in range(2, 100)]
        ones = sum(vals)
        print(f"  Bit {bit} of p(n) for n=3..100: {ones}/98 are 1 ({ones/98:.1%})")

    # Check p(n) mod small numbers
    for m in [3, 5, 7, 11, 13]:
        residues = [primes[i] % m for i in range(m, 200)]
        from collections import Counter
        dist = Counter(residues)
        # For p > m, p mod m is never 0 and roughly uniform over {1,...,m-1}\{0}
        print(f"  p(n) mod {m}: {dict(sorted(dist.items()))}")

    print("\n  p(n) mod m is equidistributed over non-zero residues (Dirichlet's theorem)")
    print("  No predictable pattern in individual digits/bits")
    print("  CONCLUSION: Digit-by-digit construction requires knowing π(x) per digit")

# =============================================================================
# Main
# =============================================================================
if __name__ == "__main__":
    print("SESSION 8: Stern-Brocot / CF / Novel Encoding Approaches")
    print("=" * 60)

    t0 = time.time()

    cf_ratio_analysis()
    prime_constant_approach()
    farey_approach()
    stern_brocot_gaps()
    riemann_r_precision()
    digit_construction()

    print(f"\n\nTotal time: {time.time() - t0:.1f}s")

    print("\n" + "=" * 60)
    print("SESSION 8 CONCLUSIONS (Stern-Brocot/CF/Novel Encoding)")
    print("=" * 60)
    print("""
1. CF of p(n)/n: First term ~ ln(n) (trivial), higher terms unpredictable
2. Prime constant C_p: Binary expansion IS prime indicator, but circular to compute
3. Farey fractions: Connected to RH but require φ(k) → factoring → same barrier
4. Stern-Brocot / gap prediction: AR(10) gives ~17% exact, need ALL gaps (O(n))
5. R(x) precision: 94-99% correct for small x, but error grows as O(√x/ln x)
6. Digit construction: Each digit of p(n) is equidistributed, no shortcut

ALL SIX APPROACHES FAIL. The barrier is information-theoretic.
""")

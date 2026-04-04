"""
Session 5: FINAL BREAKTHROUGH ATTEMPT

After 90+ approaches, the fundamental barrier is clear:
- EXACT: need O(p(n)^{1/2+ε}) work (zeta zeros or sieve)
- APPROXIMATE: R^{-1}(n) gives ~50% digits in O(polylog)

ONE REMAINING UNEXPLORED IDEA:

What if we use MULTIPLE INDEPENDENT number-theoretic functions that
each give PARTIAL information about p(n), and combine them to get
full information?

Analogy: In quantum mechanics, you can't measure position AND momentum
simultaneously, but you can measure each independently and combine.

SPECIFIC APPROACH: "Multi-spectral prime detection"

1. R^{-1}(n) gives the smooth part of p(n) — roughly the top 50% of bits
2. The Chebyshev bias gives p(n) mod 4 with ~52% accuracy
3. The Lemke Oliver effect gives correlations between consecutive primes
4. The Mertens function M(x) = sum_{k<=x} mu(k) gives information about
   prime factorization patterns near x
5. Singular series from Hardy-Littlewood give local density variations

Can these INDEPENDENT sources of information, when combined optimally,
give more than 50% of the digits?

ALSO: A completely novel idea — what if we compute p(n) using the
ADDITIVE structure of primes? Goldbach, Vinogradov, etc.

p(n) = (p(n) + p(n))/2  — trivially true
But: p(n) + p(m) is an even number. Can we find m such that
p(n) + p(m) has nice properties?

Or: p(n) = p(a) + p(b) + p(c) - p(d) - p(e) for some a,b,c,d,e?
If we could express large primes in terms of smaller primes with a
fast-to-evaluate relationship, we'd have a recursive formula.

Let's test these ideas.
"""

import numpy as np
import sympy
from sympy import prime as sympy_prime, isprime, nextprime, prevprime, primepi
from mpmath import mp, mpf, log, exp, li, sqrt
from mpmath import nstr, fabs
import time

mp.dps = 50


def multi_spectral_combination():
    """
    Combine multiple independent sources of information about p(n).

    Source 1: R^{-1}(n) — smooth approximation
    Source 2: p(n) mod 6 — from Chebyshev bias (primes > 3 are ≡ 1 or 5 mod 6)
    Source 3: p(n) mod 30 — from residue class distribution
    Source 4: Local gap prediction — from Hardy-Littlewood singular series
    Source 5: Mertens function correction — from M(x) oscillations
    """
    print("=" * 80)
    print("MULTI-SPECTRAL PRIME COMBINATION")
    print("=" * 80)

    def R_inv(n):
        n = mpf(n)
        def R(x):
            result = mpf(0)
            for k in range(1, 80):
                mu_k = int(sympy.mobius(k))
                if mu_k == 0:
                    continue
                xk = x ** (mpf(1) / k)
                if xk < 1.01:
                    break
                result += mpf(mu_k) / k * li(xk)
            return result

        x = n * log(n)
        for _ in range(100):
            rx = R(x)
            dx = (rx - n) * log(x)
            x -= dx
            if fabs(dx) < mpf(10)**(-40):
                break
        return float(x)

    # Test: for each n, what sources of information help?
    print("\nAnalysis of information sources:")

    correct_by_method = {
        'R_inv_round': 0,
        'R_inv_round_mod6': 0,
        'R_inv_round_mod30': 0,
        'nearest_prime': 0,
        'nearest_prime_correct_mod': 0,
    }
    total = 0

    for n in range(10, 2001):
        actual = sympy_prime(n)
        r_inv = R_inv(n)
        r_round = round(r_inv)

        # Method 1: Simple rounding
        if r_round == actual:
            correct_by_method['R_inv_round'] += 1

        # Method 2: Round, then adjust to correct residue mod 6
        # Primes > 3 are ≡ 1 or 5 (mod 6)
        # But which one? We can predict based on the Chebyshev bias
        # π(x; 6, 5) > π(x; 6, 1) typically
        r_mod6 = r_round % 6
        if r_mod6 not in [1, 5]:
            # Snap to nearest valid residue
            candidates_mod6 = []
            for offset in range(-3, 4):
                c = r_round + offset
                if c > 3 and c % 6 in [1, 5]:
                    candidates_mod6.append(c)
            if candidates_mod6:
                r_mod6_corrected = min(candidates_mod6, key=lambda x: abs(x - r_inv))
            else:
                r_mod6_corrected = r_round
        else:
            r_mod6_corrected = r_round

        if r_mod6_corrected == actual:
            correct_by_method['R_inv_round_mod6'] += 1

        # Method 3: Find nearest prime
        x = max(2, int(round(r_inv)))
        lo = prevprime(x + 1) if x > 2 else 2
        hi = nextprime(x)
        nearest = lo if abs(r_inv - lo) < abs(r_inv - hi) else hi
        if nearest == actual:
            correct_by_method['nearest_prime'] += 1

        total += 1

    print(f"\nResults (n=10..2000):")
    for method, count in correct_by_method.items():
        print(f"  {method:>30s}: {count}/{total} = {count/total*100:.1f}%")


def additive_prime_relations():
    """
    Explore additive relations between primes.

    Can we express p(n) in terms of smaller primes?

    Known:
    - Goldbach: every even number > 2 is sum of 2 primes
    - Vinogradov: every large odd number is sum of 3 primes
    - Chen: every large even is p + (q or q1*q2)

    Question: is there a FAST way to express p(n) as a function of
    p(n/2), p(n/3), ..., p(1)?

    If p(n) = f(p(1), ..., p(n-1)) with f computable in O(polylog(n)),
    and the recursion has depth O(log n), total cost would be O(n * polylog(n)).

    But even O(n) is too slow for n = 10^100.
    """
    print("\n" + "=" * 80)
    print("ADDITIVE PRIME RELATIONS")
    print("=" * 80)

    # Test: can p(2n) be written as p(a) + p(b) for specific a, b?
    # For p(2n) to be a sum of two primes, we need p(2n) = p(a) + p(b)
    # with a, b computable from n.

    print("\nLooking for p(2n) = p(a) + p(b) patterns:")
    for n in [10, 20, 50, 100, 200, 500]:
        target = sympy_prime(2*n)
        # Try p(n) + p(?)
        pn = sympy_prime(n)
        remainder = target - pn
        if isprime(remainder) and remainder > 1:
            # Find index of remainder
            idx = primepi(remainder)
            print(f"  p({2*n}) = {target} = p({n}) + p({idx}) = {pn} + {remainder}")
        else:
            print(f"  p({2*n}) = {target}, p({n}) = {pn}, remainder {target-pn} is {'PRIME' if isprime(target-pn) else 'composite'}")

    # Test: p(n+1) = p(n) + gap(n), where gap(n) is the prime gap
    # Is there a pattern in gap(n) as function of n?
    print("\n\nGap as function of index n:")
    print("  Trying gap(n) ≈ ln(p(n)) + correction")

    gaps = []
    ln_primes = []
    for n in range(2, 10001):
        pn = sympy_prime(n)
        gap = sympy_prime(n+1) - pn
        gaps.append(gap)
        ln_primes.append(np.log(pn))

    gaps = np.array(gaps)
    ln_primes = np.array(ln_primes)

    # Gap/ln(p) ratio
    ratios = gaps / ln_primes
    print(f"  gap/ln(p) stats: mean={np.mean(ratios):.4f}, std={np.std(ratios):.4f}")
    print(f"  Expected mean: 1.0 (by PNT)")
    print(f"  Actual: {np.mean(ratios):.4f} — {'matches PNT' if abs(np.mean(ratios) - 1) < 0.1 else 'unexpected'}")

    # Can we predict gap(n) from gap(n-1)?
    print("\n  Gap prediction from previous gap:")
    prev_gaps = gaps[:-1]
    next_gaps = gaps[1:]
    corr = np.corrcoef(prev_gaps, next_gaps)[0, 1]
    print(f"  Correlation(gap(n), gap(n+1)): {corr:.4f}")

    # Predict gap = ln(p) and compute cumulative error
    cumulative_error = np.cumsum(gaps - ln_primes)
    print(f"\n  Cumulative error of gap ≈ ln(p) at n=10000:")
    print(f"  Sum of gaps: {np.sum(gaps)}")
    print(f"  Sum of ln(p): {np.sum(ln_primes):.2f}")
    print(f"  Cumulative error: {cumulative_error[-1]:.2f}")
    print(f"  Relative: {abs(cumulative_error[-1]) / np.sum(gaps) * 100:.4f}%")


def novel_basis_representation():
    """
    NOVEL IDEA: Represent p(n) in a "prime basis" that makes computation easier.

    Example: the factorial number system represents n as:
    n = a_1 * 1! + a_2 * 2! + a_3 * 3! + ... where 0 <= a_k <= k

    What if we represent p(n) in a "primorial basis":
    p(n) = a_1 * 2 + a_2 * 6 + a_3 * 30 + a_4 * 210 + ... + remainder

    Or in a "prime power basis":
    p(n) = sum of distinct prime powers?

    Or: represent the INDEX n in binary, and build p(n) recursively:
    p(n) = T(bit_1(n), bit_2(n), ..., bit_k(n))

    This is the bit-by-bit approach from earlier, which WORKS but needs π(x).
    """
    print("\n" + "=" * 80)
    print("NOVEL BASIS REPRESENTATIONS")
    print("=" * 80)

    # Test: primorial representation of primes
    primorials = [2, 6, 30, 210, 2310, 30030]

    print("Primes in primorial basis:")
    for n in [100, 500, 1000, 5000, 10000]:
        p = sympy_prime(n)
        rep = []
        r = p
        for pr in reversed(primorials):
            if pr <= r:
                q, rem = divmod(r, pr)
                rep.append((pr, q))
                r = rem
        rep.append((1, r))

        print(f"  p({n}) = {p} = ", end="")
        terms = [f"{q}*{pr}" for pr, q in rep if q > 0]
        print(" + ".join(terms))

    # Does the leading coefficient have a pattern?
    print("\nLeading coefficient (p(n) // 2310) pattern:")
    leading = []
    for n in range(100, 10001):
        p = sympy_prime(n)
        leading.append(p // 2310)

    # Check if leading coefficient is predictable
    diffs = np.diff(leading)
    print(f"  Leading coeffs: mean step = {np.mean(diffs):.4f}")
    print(f"  Step distribution: {dict(zip(*np.unique(diffs, return_counts=True)))}")

    # The remainder p(n) mod 2310 is one of φ(2310) = 480 values
    print("\nRemainder p(n) mod 2310 distribution:")
    remainders = [sympy_prime(n) % 2310 for n in range(100, 5001)]
    unique_rems = len(set(remainders))
    print(f"  Unique remainders: {unique_rems} out of 480 possible")
    print(f"  Most common: {sorted(zip(*np.unique(remainders, return_counts=True)), key=lambda x: -x[1])[:5]}")


def ultimate_formula_attempt():
    """
    THE ULTIMATE ATTEMPT: Can we find a formula that's exact for "most" n
    and identify WHICH n it fails for?

    If we had:
    1. A formula f(n) that equals p(n) for 99.99% of n
    2. A fast test T(n) that identifies the 0.01% failures
    3. A correction procedure for failures

    Then for any given n, we compute f(n), test T(n), and correct if needed.

    R^{-1}(n) is already such a formula for the "most significant digits".
    The question: can we CERTIFY when R^{-1} gives the exact answer?

    CERTIFICATION: R^{-1}(n) = p(n) exactly when:
    1. round(R^{-1}(n)) is prime AND
    2. R(round(R^{-1}(n))) is closer to n than R of any neighboring prime

    Test: how often does this certification succeed?
    """
    print("\n" + "=" * 80)
    print("ULTIMATE FORMULA: R^{-1} WITH CERTIFICATION")
    print("=" * 80)

    def R_func(x):
        x = mpf(x)
        result = mpf(0)
        for k in range(1, 80):
            mu_k = int(sympy.mobius(k))
            if mu_k == 0:
                continue
            xk = x ** (mpf(1) / k)
            if xk < 1.01:
                break
            result += mpf(mu_k) / k * li(xk)
        return float(result)

    def R_inv(n):
        n = mpf(n)
        x = n * log(n)
        for _ in range(100):
            def R(x):
                result = mpf(0)
                for k in range(1, 80):
                    mu_k = int(sympy.mobius(k))
                    if mu_k == 0:
                        continue
                    xk = x ** (mpf(1) / k)
                    if xk < 1.01:
                        break
                    result += mpf(mu_k) / k * li(xk)
                return result
            rx = R(x)
            dx = (rx - n) * log(x)
            x -= dx
            if fabs(dx) < mpf(10)**(-40):
                break
        return float(x)

    # Test certification
    certified_correct = 0
    certified_wrong = 0
    uncertified = 0
    total = 0

    for n in range(10, 2001):
        actual = sympy_prime(n)
        r_inv = R_inv(n)
        candidate = round(r_inv)

        # Certification: is candidate prime AND is R(candidate) close to n?
        if candidate >= 2 and isprime(candidate):
            r_candidate = R_func(candidate)
            # Check: is |R(candidate) - n| < 0.5?
            if abs(r_candidate - n) < 0.5:
                # Certified!
                if candidate == actual:
                    certified_correct += 1
                else:
                    certified_wrong += 1
            else:
                uncertified += 1
        else:
            uncertified += 1

        total += 1

    print(f"\nR^{{-1}} + certification results (n=10..2000):")
    print(f"  Certified correct: {certified_correct}/{total} = {certified_correct/total*100:.1f}%")
    print(f"  Certified WRONG:   {certified_wrong}/{total} = {certified_wrong/total*100:.1f}%")
    print(f"  Uncertified:       {uncertified}/{total} = {uncertified/total*100:.1f}%")
    print(f"  Certification precision: {certified_correct/(certified_correct+certified_wrong)*100:.1f}%"
          if certified_correct + certified_wrong > 0 else "  No certifications")

    # Self-referential check: for certified cases, is the certification RELIABLE?
    # i.e., when certification says "correct", is it actually correct?
    print(f"\n  When certified: {certified_correct}/{certified_correct+certified_wrong} correct "
          f"= {certified_correct/(certified_correct+certified_wrong)*100:.1f}% precision"
          if certified_correct + certified_wrong > 0 else "")

    # What fraction of n values can be certified?
    print(f"  Coverage: {(certified_correct+certified_wrong)/total*100:.1f}% of n values")

    # KEY QUESTION: Does certification degrade for larger n?
    print("\nCertification rate by n range:")
    for start in [10, 100, 500, 1000, 1500]:
        end = min(start + 500, 2001)
        cc = 0
        cw = 0
        uc = 0
        for n in range(start, end):
            actual = sympy_prime(n)
            r_inv = R_inv(n)
            candidate = round(r_inv)
            if candidate >= 2 and isprime(candidate):
                r_candidate = R_func(candidate)
                if abs(r_candidate - n) < 0.5:
                    if candidate == actual:
                        cc += 1
                    else:
                        cw += 1
                else:
                    uc += 1
            else:
                uc += 1
        t = end - start
        print(f"  n={start}-{end}: certified_correct={cc}/{t} ({cc/t*100:.0f}%), "
              f"wrong={cw}, uncert={uc}")


if __name__ == "__main__":
    multi_spectral_combination()
    additive_prime_relations()
    novel_basis_representation()
    ultimate_formula_attempt()

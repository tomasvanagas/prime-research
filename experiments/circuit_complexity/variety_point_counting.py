"""
Session 14: Can pi(x) be expressed as the number of F_q-rational points
on an algebraic variety?

Background: For a variety V over F_q, the Weil conjectures (proved by
Deligne) give:
  #V(F_q) = sum_i (-1)^i Tr(Frob_q | H^i_c(V))

If we could find a family V_x of varieties over F_q (for fixed q, say q=2)
such that:
  #V_x(F_2) = pi(x)

and V_x has a poly(N)-bit description with dim = O(poly(N)), then
counting points on V_x would potentially be efficient.

For smooth projective curves of genus g, #C(F_q) = q + 1 - sum_{i=1}^{2g} alpha_i
where |alpha_i| = sqrt(q). So #C(F_q) ranges from q + 1 - 2g*sqrt(q) to
q + 1 + 2g*sqrt(q). For F_2: 3 - 2g*sqrt(2) to 3 + 2g*sqrt(2).
For pi(x) ~ x/ln(x), we'd need g ~ pi(x) / (2*sqrt(2)) ~ x / (2*sqrt(2)*ln(x)).
This is exponential in N — genus too large.

But HIGHER-DIMENSIONAL varieties can have more points. For a d-dimensional
variety over F_q: #V(F_q) ~ q^d. For F_2: #V(F_2) ~ 2^d.
To get pi(x) ~ 2^N / N, we need d ~ N - log N. This is polynomial in N!

So in principle, a d-dimensional variety with d ≈ N could have ~2^N / N points,
which matches pi(x).

The question: can we CONSTRUCT such a variety whose point count equals pi(x) exactly?
"""

import numpy as np
from sympy import primepi, isprime, primerange
from itertools import product as cartesian_product
from math import gcd, isqrt

def affine_variety_point_count():
    """
    For an affine variety V ⊂ F_2^d defined by polynomial equations,
    #V(F_2) = #{x ∈ F_2^d : f_1(x) = ... = f_m(x) = 0}.

    Over F_2, x_i^2 = x_i, so every polynomial reduces to multilinear.
    #V(F_2) is the number of 0/1 solutions to a system of linear equations
    over F_2 (if the polynomials are linearized).

    Actually, over F_2, any polynomial is equivalent to its multilinear form.
    So #V(F_2) = #{x ∈ F_2^d : f_i(x) = 0 for all i} where each f_i is
    a multilinear polynomial over F_2.

    Can we find f_1, ..., f_m in d variables such that #V(F_2) = pi(x)?

    This is equivalent to finding a Boolean function f: {0,1}^d → {0,1}
    with exactly pi(x) zeros (or ones). Trivially possible but circular.

    The NON-TRIVIAL question: can the variety V_x be described by poly(N)
    equations of degree poly(N) in poly(N) variables, with COEFFICIENTS
    efficiently computable from x?
    """
    print("=" * 70)
    print("Affine variety point counting over F_2")
    print("=" * 70)

    # For F_2^d, the total number of points is 2^d.
    # We want #V(F_2) = pi(x).
    # If V is defined by one linear equation, #V = 2^{d-1}.
    # If V is defined by k linearly independent equations, #V = 2^{d-k}.
    # So for pi(x) = 2^{d-k}, we need d - k = log2(pi(x)).
    # This only works if pi(x) is a power of 2.

    # For non-powers-of-2: need NONLINEAR equations.
    # The number of F_2-points of a degree-2 (quadratic) variety
    # is not restricted to powers of 2.

    # Example: V: x1*x2 = 0 in F_2^2 has solutions (0,0), (0,1), (1,0) = 3 points.
    # V: x1*x2 + x1 = 0 has solutions x1=0 (2 points) and x1=1,x2=1 (1 point) = 3.

    # More generally: a single quadratic form Q(x) = 0 in F_2^d has
    # 2^{d-1} + 2^{(d-1)/2 - 1} or 2^{d-1} - 2^{(d-1)/2 - 1} solutions
    # depending on the rank and type of Q.

    # For a SYSTEM of quadratic forms, the point count can be arbitrary.

    # The question reduces to: can we construct a poly(N)-size system of
    # polynomial equations over F_2 whose solution count = pi(x)?

    # This is EQUIVALENT to: is pi(x) computable by a poly(N)-size
    # arithmetic circuit over F_2? (The variety encodes the circuit.)

    # But over F_2, we showed (Session 13) that the ANF degree of the
    # prime indicator is Theta(N) and sparsity is 50%. This means
    # the polynomial DESCRIBING primes is essentially random over F_2.

    # However, we don't need the prime indicator — we need the COUNT pi(x).
    # These are different problems!

    print("Over F_2: variety V_x with #V_x(F_2) = pi(x) is equivalent to")
    print("finding a GF(2) circuit of poly(N) size computing the COUNT pi(x).")
    print("This is the NC question over GF(2).")
    print()

    # Let's try a different field: F_p for various primes p.
    # Over F_p, #V(F_p) can be any integer in [0, p^d].
    # For p = 2: #V ∈ [0, 2^d], need pi(x) ≈ 2^N / N.
    # For p = 3: #V ∈ [0, 3^d], need pi(x) ≈ 3^d for some d ~ N/1.58.

    for x in [100, 1000]:
        pi_x = int(primepi(x))
        N = x.bit_length()
        for q in [2, 3, 5, 7]:
            # Minimum dimension needed: d such that q^d >= pi(x)
            import math
            d_min = math.ceil(math.log(pi_x) / math.log(q)) if pi_x > 0 else 0
            print(f"x={x}, pi(x)={pi_x}, F_{q}: min dimension d={d_min}, N={N}")


def hypersurface_approach():
    """
    A hypersurface V: f(x_1, ..., x_d) = 0 in F_q^d has
    #V(F_q) = q^{d-1} + O(q^{(d-1)/2}) for non-singular hypersurfaces.

    For a specific count c, we need:
    q^{d-1} - c has absolute value at most (d * q^{(d-1)/2})
    (by Weil-Deligne bounds for degree-d hypersurfaces).

    Can we find a SPECIFIC polynomial f_x of degree D in d variables over F_q
    such that #{solutions in F_q^d} = pi(x)?

    The Chevalley-Warning theorem: if deg(f) < d, then #V(F_q) ≡ 0 (mod q).
    So if q does not divide pi(x), we need deg(f) >= d.

    More importantly: can the COEFFICIENTS of f_x be computed efficiently from x?
    """
    print("\n" + "=" * 70)
    print("Hypersurface point counting approach")
    print("=" * 70)

    # For small cases: find explicit polynomials
    # Over F_2, a polynomial f in d variables is multilinear.
    # A multilinear polynomial of degree d in d variables:
    # f = sum_{S ⊆ [d]} c_S * prod_{i in S} x_i where c_S ∈ F_2

    # The number of zeros of f in F_2^d = 2^d - wt(f) where
    # wt(f) is the number of x with f(x) = 1.

    # So #{f(x) = 0} = 2^d - #{f(x) = 1}.
    # We need #{f(x) = 0} = pi(x), so #{f(x) = 1} = 2^d - pi(x).

    # For x = 10, pi(10) = 4. We need #{f(x) = 1} = 2^d - 4.
    # For d = 3: #{f=1} = 4, #{f=0} = 4. So f has weight 4.
    # Any multilinear function on 3 variables with exactly 4 ones.
    # Example: f = x1 XOR x2 has weight 4 on {0,1}^3 (NO, that's wrong)
    # f = x1 on F_2^3 has weight 4 (when x1=1: 4 points). #{f=0} = 4.
    # So f = x_1 gives #{f=0} = 4 = pi(10). But this depends on d, not x!

    # For x = 100, pi(100) = 25. Need #{f=0} = 25. So #{f=1} = 2^d - 25.
    # For d = 5: 2^5 = 32, need #{f=1} = 7. A function on 5 vars with weight 7.
    # For d = 6: 2^6 = 64, need #{f=1} = 39.

    # The PROBLEM: constructing f from x requires knowing pi(x). Circular!

    # UNLESS there's a FORMULA for the coefficients of f_x that's
    # computable in poly(N) time without knowing pi(x).

    # Over larger fields: more room. Over F_3^d:
    # A polynomial f ∈ F_3[x_1,...,x_d] has #{f=0 in F_3^d} values.
    # By Chevalley-Warning: if deg(f) < d, then #{f=0} ≡ 0 (mod 3).

    print("Hypersurface approach: constructing f_x requires knowing pi(x).")
    print("Over F_q: Chevalley-Warning gives congruences #{V} ≡ 0 (mod q)")
    print("when deg < dim. This constrains but doesn't determine #{V}.")
    print()
    print("The fundamental issue: the COEFFICIENTS of f_x encode pi(x),")
    print("so computing them is as hard as computing pi(x) itself.")
    print()

    # Check: what congruences does pi(x) satisfy?
    print("Congruence structure of pi(x):")
    for q in [2, 3, 5, 7]:
        residues = [int(primepi(x)) % q for x in range(2, 1001)]
        from collections import Counter
        counts = Counter(residues)
        print(f"  pi(x) mod {q} for x=2..1000: {dict(sorted(counts.items()))}")

    # The distribution of pi(x) mod q is essentially uniform (from the
    # prime indicator's entropy of 0.537 bits per step).


def trace_of_frobenius():
    """
    For an elliptic curve E over F_p, #E(F_p) = p + 1 - a_p where
    a_p = Tr(Frob_p) satisfies |a_p| <= 2*sqrt(p).

    The sequence a_p (over varying p) encodes the L-function of E.
    But we want #E(F_p) for FIXED E and varying p... that gives a_p
    which relates to modular forms.

    What if we fix p (say p=2) and VARY the curve E?
    #E(F_2) ∈ {1, 2, 3, 4, 5} for elliptic curves over F_2.
    These are small numbers — not useful for encoding pi(x).

    For HIGHER-GENUS curves C over F_2:
    #C(F_2) can range from 0 to 2g + 1 + 2*floor(2*sqrt(2)*g)
    by the Hasse-Weil bound.

    For abelian varieties of dimension g over F_q:
    #A(F_q) can be any value in [a*q^g, b*q^g] for certain a, b > 0.

    The point: higher-dimensional varieties give more flexibility.
    """
    print("\n" + "=" * 70)
    print("Trace of Frobenius / variety approach")
    print("=" * 70)

    # For a genus-g curve over F_2:
    # #C(F_2) = 2 + 1 - sum_{i=1}^{2g} alpha_i where |alpha_i| = sqrt(2)
    # So #C(F_2) = 3 - sum, and sum ∈ [-2g*sqrt(2), 2g*sqrt(2)]
    # So #C(F_2) ∈ [3 - 2g*sqrt(2), 3 + 2g*sqrt(2)]

    # For pi(x) ≈ x / ln(x), we need:
    # 3 + 2g*sqrt(2) >= pi(x) → g >= (pi(x) - 3) / (2*sqrt(2))

    import math
    for x in [100, 1000, 10**6, 10**10, 10**100]:
        if x <= 10**7:
            pi_x = int(primepi(x))
        else:
            pi_x = int(x / math.log(x))
        N = math.ceil(math.log2(x))
        g_min = math.ceil((pi_x - 3) / (2 * math.sqrt(2)))
        print(f"x=10^{math.log10(x):.0f}: pi(x)≈{pi_x:.2e}, genus g≈{g_min:.2e}, N={N}")
        print(f"  Genus / N = {g_min/N:.1f} — {'poly(N)' if g_min <= N**3 else 'NOT poly(N)'}")

    # For x = 10^100: pi(x) ≈ 4.3 * 10^97, g ≈ 1.5 * 10^97.
    # N = 333. g/N = 4.5 * 10^94. NOT polynomial in N.

    # So the curve approach requires genus exponential in N.
    # What about higher-dimensional varieties?

    print("\nFor d-dimensional abelian varieties over F_2:")
    print("  #A(F_2) ~ 2^d, so need d ≈ log2(pi(x)) ≈ N - log(N)")
    print("  d is POLYNOMIAL in N!")
    print("  But constructing such a variety efficiently from x is the challenge.")
    print()

    # The KNOWN efficient point-counting algorithms for varieties:
    # - Schoof's algorithm for elliptic curves: O(log^8 p) over F_p
    # - Kedlaya's algorithm for hyperelliptic curves: O(g^3 * p^{1/2} * polylog)
    # - Point counting on general varieties: exponential in dimension

    # For a d-dimensional variety with d ≈ N, Kedlaya-type algorithms
    # would take time polynomial in 2^d ≈ pi(x) — which is our current best!

    print("Kedlaya point counting on d-dimensional variety: O(d^3 * p^d * polylog)")
    print("For d ≈ N over F_2: O(N^3 * 2^N * polylog) ≈ O(x * polylog)")
    print("This is EQUIVALENT to the sieve method — no improvement!")
    print()
    print("CONCLUSION: Algebraic geometry point counting does NOT provide")
    print("a faster algorithm because:")
    print("1. Low-dimensional varieties can't encode pi(x) (too few points)")
    print("2. High-dimensional varieties have slow point-counting algorithms")
    print("3. The Weil-Deligne formula encodes the SAME information as the")
    print("   Riemann explicit formula (Frobenius eigenvalues = zeta zeros)")


if __name__ == '__main__':
    affine_variety_point_count()
    hypersurface_approach()
    trace_of_frobenius()

#!/usr/bin/env python3
"""
A7 (extension) — degree-d components of f_chi_P^{(n)} and support-hypergraph
rigidity.

The polynomial f_chi_P^{(n)} has monomials of mixed degrees (1, 2, ..., d_max).
Standard GCT compares HOMOGENEOUS polynomials.  Decomposing by degree:

    f_chi_P^{(n)}_d := sum_{|S|=d, val(S) prime} prod_{i in S} x_i

Each f_chi_P^{(n)}_d is a homogeneous multi-affine degree-d polynomial in n
variables, and falls in the GL_n-irrep Sym^d(C^n)*.

We compute Stab Lie alg dim for each degree component, plus the
support-hypergraph rigidity test:
    "does there exist any choice of coefficients (c_S)_{S in support} such
     that the polynomial sum_S c_S prod_{i in S} x_i has non-trivial Stab?"

If NO, the support pattern is "Lie-rigid" — no f with that support can have
non-trivial Lie stabilizer.  If YES, the χ_P-coefficient choice (all 1's)
might be a non-rigid case, in which case the rigidity is a structural fact
about prime-encoding hypergraphs (rather than a no-op statement).
"""

import itertools
import random

import sympy as sp


def sieve_primes(M):
    if M < 2:
        return set()
    s = [True] * (M + 1)
    s[0] = s[1] = False
    for i in range(2, int(M**0.5) + 1):
        if s[i]:
            for j in range(i * i, M + 1, i):
                s[j] = False
    return {i for i, x in enumerate(s) if x}


def chi_p_support_by_degree(n):
    primes = sieve_primes(2**n)
    by_deg = {}
    for k in range(n + 1):
        for S in itertools.combinations(range(1, n + 1), k):
            if sum(2**(i - 1) for i in S) in primes:
                by_deg.setdefault(k, []).append(S)
    return by_deg


def stab_lie_dim_from_support(coeffs, n, x_syms=None):
    if x_syms is None:
        x_syms = sp.symbols(f'x_1:{n + 1}')
    f = sp.Integer(0)
    for S, c in coeffs.items():
        term = sp.Integer(c)
        for i in S:
            term *= x_syms[i - 1]
        f += term
    grad = [sp.diff(f, xi) for xi in x_syms]
    poly_x_coeffs = {}
    for i in range(n):
        gi = sp.expand(grad[i])
        if gi == 0:
            continue
        gi_poly = sp.Poly(gi, *x_syms)
        for monom, coef in gi_poly.terms():
            for j in range(n):
                new_mon = list(monom)
                new_mon[j] += 1
                new_mon = tuple(new_mon)
                if new_mon not in poly_x_coeffs:
                    poly_x_coeffs[new_mon] = {}
                key = (i, j)
                poly_x_coeffs[new_mon][key] = poly_x_coeffs[new_mon].get(key, sp.Rational(0)) + sp.Rational(coef)
    rows = []
    for monom, ad in poly_x_coeffs.items():
        row = [sp.Rational(0)] * (n * n)
        for (i, j), c in ad.items():
            if c != 0:
                row[i * n + j] = c
        if any(v != 0 for v in row):
            rows.append(row)
    if not rows:
        return n * n
    return n * n - sp.Matrix(rows).rank()


def support_hypergraph_rigidity(support, n, n_trials=200, seed=0):
    """Sample n_trials random coefficient vectors on the support; report
    distribution of Stab dim."""
    rng = random.Random(seed)
    dims = []
    for trial in range(n_trials):
        c = {S: rng.randint(-10, 10) or 1 for S in support}
        dims.append(stab_lie_dim_from_support(c, n))
    from collections import Counter
    return Counter(dims)


def main():
    print("=" * 78)
    print("Degree-components and support rigidity for f_chi_P^{(n)}")
    print("=" * 78)
    for n in [4, 5, 6]:
        by_deg = chi_p_support_by_degree(n)
        print(f"\n--- n = {n} ---")
        print(f"  primes <= 2^{n}: {sieve_primes(2**n)}")
        print(f"  degree decomposition of support:")
        for d in sorted(by_deg):
            print(f"    deg {d}: {by_deg[d]}  ({len(by_deg[d])} monos)")

        # Stab dim for each degree-d component (just χ_P coefs, all 1)
        x_syms = sp.symbols(f'x_1:{n + 1}')
        for d in sorted(by_deg):
            coeffs = {S: 1 for S in by_deg[d]}
            try:
                dim = stab_lie_dim_from_support(coeffs, n, x_syms)
                print(f"    deg {d}: Stab Lie dim of f_d (chi_P) = {dim}")
            except Exception as e:
                print(f"    deg {d}: ERROR {e}")

        # Full polynomial Stab — known to be 0 from prior result
        # Support rigidity test (= over coefficient space, always 0?)
        full_support = sum(by_deg.values(), [])
        print(f"\n  Support-rigidity test (random integer coefs in [-10,10]):")
        n_trials = 100 if n <= 5 else 50
        counter = support_hypergraph_rigidity(full_support, n, n_trials=n_trials,
                                               seed=42 + n)
        for k in sorted(counter):
            print(f"    Stab dim = {k:>2}: {counter[k]} / {n_trials}")
        if list(counter.keys()) == [0]:
            print(f"    => Support is LIE-RIGID (no coefficient choice gives non-trivial Stab)")
        else:
            print(f"    => Support has SOME coefficient choices with non-trivial Stab "
                  f"({sum(v for k, v in counter.items() if k > 0)} / {n_trials})")


if __name__ == "__main__":
    main()

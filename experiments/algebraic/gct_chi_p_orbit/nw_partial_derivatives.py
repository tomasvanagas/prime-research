#!/usr/bin/env python3
"""
A7 (extension) — Nisan-Wigderson partial-derivative spaces of f_chi_P^{(n)}.

For a polynomial f in x_1,...,x_n, the k-th order partial-derivative space is
    PD_k(f) := span{ d^k f / dx_{i_1}...dx_{i_k} : 1 <= i_1 < ... < i_k <= n },
treated as polynomials in x_1,...,x_n.  dim PD_k(f) is a Nisan-Wigderson
invariant: a strong lower bound on the formula size of f scales as
    L(f) >= dim PD_k(f) / dim PD_k(f_generic)
when the latter ratio is < 1.  For multi-affine polynomials, dim PD_k is
combinatorially bounded by (n choose k) and equals the size of the support
of f restricted to subsets containing the differentiation index set, modulo
linear dependencies.

The KEY question (FAL-NW): does dim PD_k(f_chi_P^{(n)}) deviate from
matched-support random multi-affine polynomials?

If dim PD_k(f_chi_P) < dim PD_k(matched random) at some k, then f_chi_P has
HIDDEN linear dependencies among its partial derivatives — a Nisan-Wigderson
formula lower bound applies.  If dim is matched, NW does not give a lower
bound.
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


def chi_p_support_indices(n):
    """Subsets S of {1..n} with val(S) prime (val(S) = sum 2^{i-1})."""
    primes = sieve_primes(2**n)
    return tuple(
        S for k in range(n + 1)
        for S in itertools.combinations(range(1, n + 1), k)
        if sum(2**(i - 1) for i in S) in primes
    )


def pd_space_dim_multilinear(coeffs, n, k):
    """
    Given a multilinear polynomial as {S -> coef} (S is a frozenset/tuple of
    variable indices in 1..n), compute dim of k-th order partial-derivative
    space.

    For multilinear f = sum_S c_S prod_{i in S} x_i:
        d^k f / dx_T = sum_{S: T subset S} c_S prod_{i in S\T} x_i
                     = sum_{S: T subset S} c_S * monomial(S\T)
    treating only |T| = k.

    Build matrix M[T][U] = coefficient of monomial prod_{i in U} x_i in
    d^k f / dx_T:
        M[T][U] = c_{T cup U} if (T cup U) in support and T disjoint U,
                  0 otherwise.
    Then dim PD_k = rank(M).

    Args:
      coeffs: dict {tuple S -> int/Rational coefficient}
      n: number of variables
      k: differentiation order

    Returns: rank (int)
    """
    # Rows: T subsets of size k.  Columns: U subsets of size up to (n-k).
    # We only need cols corresponding to monomials that appear in d^k f/dx_T
    # for some T.
    if k > n:
        return 0
    # Collect all monomials U that arise
    rows = list(itertools.combinations(range(1, n + 1), k))
    # All possible U: subsets of {1..n} of size up to n-k that disjoint from
    # some T such that T ∪ U is in support
    support_set = set(coeffs.keys())
    # All "tail" monomials: U = S \ T for S in support, T subset S, |T|=k
    cols_set = set()
    for S in support_set:
        if len(S) < k:
            continue
        for T in itertools.combinations(S, k):
            U = tuple(sorted(set(S) - set(T)))
            cols_set.add(U)
    cols = sorted(cols_set)
    col_index = {U: i for i, U in enumerate(cols)}
    if not cols:
        return 0
    # Build the matrix
    M = []
    for T in rows:
        row = [sp.Rational(0)] * len(cols)
        Tset = set(T)
        for S in support_set:
            if not Tset.issubset(set(S)):
                continue
            U = tuple(sorted(set(S) - Tset))
            row[col_index[U]] += sp.Rational(coeffs[S])
        M.append(row)
    return sp.Matrix(M).rank()


def chi_p_coeffs(n):
    """coeffs dict: subset -> 1 for each subset whose val is prime."""
    primes = sieve_primes(2**n)
    return {
        S: 1
        for k in range(n + 1)
        for S in itertools.combinations(range(1, n + 1), k)
        if sum(2**(i - 1) for i in S) in primes
    }


def random_matched_support(coeffs, rng):
    """Same support as coeffs, random integer coefs from 1..7."""
    return {S: rng.randint(1, 7) for S in coeffs}


def pds_test(n, n_trials=100, seed=0):
    """Compute dim PD_k for f_chi_P^{(n)} and compare to baseline."""
    f_coeffs = chi_p_coeffs(n)
    print(f"\nn = {n}, support size = {len(f_coeffs)}")
    print(f"  support = {sorted(f_coeffs.keys())}")
    rng = random.Random(seed)

    for k in range(n + 1):
        d_chi = pd_space_dim_multilinear(f_coeffs, n, k)
        baseline = []
        for _ in range(n_trials):
            g = random_matched_support(f_coeffs, rng)
            d_g = pd_space_dim_multilinear(g, n, k)
            baseline.append(d_g)
        mu = sum(baseline) / n_trials
        sd = (sum((x - mu)**2 for x in baseline) / n_trials)**0.5
        if sd > 0:
            z = (d_chi - mu) / sd
        else:
            z = float('inf') if d_chi != mu else 0.0
        mn, mx = min(baseline), max(baseline)
        z_str = f"{z:+.2f}" if z != float('inf') else "  inf"
        print(f"  k={k}: dim PD_k(chi_P) = {d_chi:>3}, "
              f"baseline mean={mu:5.2f} std={sd:.2f} range=[{mn},{mx}], z={z_str}")


def main():
    print("=" * 78)
    print("Nisan-Wigderson partial-derivative space dim of f_chi_P^{(n)}")
    print("=" * 78)
    for n in [4, 5, 6, 7]:
        pds_test(n, n_trials=100 if n <= 6 else 30, seed=42 + n)


if __name__ == "__main__":
    main()

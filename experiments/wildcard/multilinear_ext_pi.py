"""
Multilinear extension (MLE) of the prime indicator chi_P over the boolean
hypercube, viewed as an n-bit integer x in [0, 2^n).

Idea (sumcheck-style): for f: {0,1}^n -> F_q, the unique multilinear
polynomial F: F_q^n -> F_q that agrees with f on the hypercube has
coefficients given by a Walsh-Hadamard-style transform. If chi_P viewed
this way has SPARSE multilinear coefficients (poly(n) nonzeros), then
F can be evaluated in poly(n) = O(polylog N) and so pi(x) = sum over
prefix-bounded x via sumcheck.

This experiment: compute MLE coefficients of chi_P on n-bit integers
for n = 4..14 (N = 16..16384). Measure number of nonzero coefficients,
their effective sparsity, and the L1/L2 norms.

If sparsity scales like 2^(c*n) with c near 1, the MLE is "dense" and
gives no advantage. If sparsity scales like poly(n), this is a genuine
shortcut - sumcheck would yield O(polylog N) evaluation.

Note: the MLE is unique and equal to the Walsh-Hadamard transform on
{0,1}-valued coordinates after a change of variable. We compute the
"Mobius" / zeta transform over the boolean lattice and count nonzeros.
"""

import numpy as np
from sympy import isprime


def chi_P_vector(n):
    """Return chi_P over 0..2^n-1 as a numpy array of int8."""
    N = 1 << n
    out = np.zeros(N, dtype=np.int64)
    for x in range(N):
        if isprime(x):
            out[x] = 1
    return out


def boolean_mobius(f):
    """In-place subset-sum (zeta) transform over boolean lattice; this
    converts the function-on-hypercube to its MLE coefficient vector
    (in the multilinear basis prod_{i in S} x_i).
    Length must be 2^n.
    Reference: standard subset-convolution machinery; this transform
    of a {0,1} function gives the unique multilinear polynomial that
    matches it on the cube via interpolation in Lagrange basis.
    """
    g = f.copy().astype(object)  # exact int arithmetic
    N = len(g)
    n = int(np.log2(N))
    assert (1 << n) == N
    # Mobius transform: g[T] = sum_{S subset T} (-1)^{|T-S|} f[S]
    # Implemented as iterated 1D Mobius across each bit.
    for bit in range(n):
        step = 1 << bit
        for i in range(0, N, step << 1):
            for j in range(step):
                # subset transform with sign for Mobius
                a = g[i + j]
                b = g[i + j + step]
                g[i + j] = a
                g[i + j + step] = b - a
    return g


def main():
    print("=== Multilinear extension sparsity of chi_P on {0,1}^n ===\n")
    print(f"  n   N     popcount  nonzero_MLE  L1_norm   max_|coef|   nz/N    nz/n^c?")

    for n in range(4, 15):
        N = 1 << n
        f = chi_P_vector(n)
        popc = int(f.sum())  # = pi(N-1) approx
        c = boolean_mobius(f)
        nz = int(sum(1 for v in c if v != 0))
        l1 = int(sum(abs(int(v)) for v in c))
        m = int(max(abs(int(v)) for v in c))
        ratio = nz / N
        # check polynomial-in-n sparsity: log nz / log n
        from math import log
        sparsity_exponent = log(max(nz, 1)) / log(max(n, 2))
        print(f"  {n:2d}  {N:5d}  {popc:5d}     {nz:6d}      {l1:8d}  {m:10d}  {ratio:6.4f}  exp={sparsity_exponent:5.2f}")

    print("\n" + "=" * 70)
    print("Verdict heuristic:")
    print("  - If 'nz/N' converges to ~1 -> MLE is dense, no shortcut.")
    print("  - If 'nz' grows polynomially in n (exp << log2 N) -> sparse, polylog eval.")
    print("  - If 'max_|coef|' grows like 2^n -> norm explodes; sumcheck exact arithmetic")
    print("    needs O(N) bits anyway, killing the polylog advantage.")
    print("=" * 70)


if __name__ == "__main__":
    main()

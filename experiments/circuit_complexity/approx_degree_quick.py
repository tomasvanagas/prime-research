#!/usr/bin/env python3
"""
Quick approximate degree computation for the prime indicator.
Uses LP to find minimum-degree polynomial approximating chi_P on {0,1}^N.
"""
import numpy as np
from sympy import isprime
from scipy.optimize import linprog
from itertools import combinations
import time

def truth_table_prime(N):
    """Prime indicator on {0,1}^N."""
    size = 1 << N
    table = np.zeros(size, dtype=np.float64)
    for x in range(size):
        table[x] = 1.0 if (x >= 2 and isprime(x)) else 0.0
    return table

def monomial_matrix(N, max_degree):
    """Build matrix of monomial evaluations on {0,1}^N up to given degree."""
    size = 1 << N
    # Enumerate all subsets S of [N] with |S| <= max_degree
    columns = []
    for d in range(max_degree + 1):
        for S in combinations(range(N), d):
            col = np.ones(size, dtype=np.float64)
            for i in S:
                for x in range(size):
                    if not ((x >> i) & 1):
                        col[x] = 0.0
            columns.append(col)
    return np.column_stack(columns) if columns else np.ones((size, 1))

def approx_degree_for_threshold(f, N, epsilon):
    """Find minimum degree d such that there exists polynomial p with |p(x)-f(x)| <= epsilon for all x."""
    size = 1 << N

    for d in range(N + 1):
        M = monomial_matrix(N, d)
        n_coeffs = M.shape[1]

        # LP: minimize 0 subject to M*c - f <= epsilon, f - M*c <= epsilon
        # Variables: c (n_coeffs)
        # Constraints: M*c <= f + epsilon, -M*c <= -f + epsilon
        # i.e., M*c <= f + eps AND M*c >= f - eps

        c_obj = np.zeros(n_coeffs)  # feasibility problem
        A_ub = np.vstack([M, -M])
        b_ub = np.concatenate([f + epsilon, -f + epsilon])

        result = linprog(c_obj, A_ub=A_ub, b_ub=b_ub, method='highs',
                        options={'presolve': True})

        if result.success:
            return d

    return N  # worst case

def main():
    print("APPROXIMATE DEGREE OF PRIME INDICATOR")
    print("=" * 60)

    for N in range(4, 13):
        t0 = time.time()
        f = truth_table_prime(N)
        size = 1 << N
        n_primes = int(np.sum(f))
        density = n_primes / size

        print(f"\nN={N}: {size} inputs, {n_primes} primes (density {density:.3f})")

        # Find approximate degree for various epsilon
        for eps in [0.49, 0.4, 0.3, 0.2, 0.1, 0.01]:
            d = approx_degree_for_threshold(f, N, eps)
            elapsed = time.time() - t0
            print(f"  eps={eps:.2f}: adeg={d} ({elapsed:.1f}s)")
            if elapsed > 120:
                print("  (timeout, skipping smaller eps)")
                break

    print("\n" + "=" * 60)
    print("SCALING ANALYSIS")
    print("=" * 60)
    print("\nKey: epsilon=0.49 is sufficient for rounding (determines if prime or not)")
    print("If adeg(chi_P, 0.49) = O(polylog N), quantum query complexity is polylog")
    print("If adeg(chi_P, 0.49) = Theta(N), then Omega(N) quantum queries needed")

if __name__ == '__main__':
    main()

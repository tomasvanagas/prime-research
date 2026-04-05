#!/usr/bin/env python3
"""
Approximate degree of prime indicator - small N only.
Focus on N=4..10 where LP is tractable.
"""
import numpy as np
from sympy import isprime
from scipy.optimize import linprog
from itertools import combinations
import time
import sys

def truth_table_prime(N):
    size = 1 << N
    return np.array([1.0 if (x >= 2 and isprime(x)) else 0.0 for x in range(size)])

def monomial_matrix(N, max_degree):
    size = 1 << N
    columns = []
    for d in range(max_degree + 1):
        for S in combinations(range(N), d):
            col = np.ones(size)
            for i in S:
                bits = np.array([(x >> i) & 1 for x in range(size)], dtype=np.float64)
                col *= bits
            columns.append(col)
    return np.column_stack(columns)

def check_feasibility(f, M, epsilon):
    """Check if degree-d polynomial can epsilon-approximate f."""
    n = M.shape[1]
    c_obj = np.zeros(n)
    A_ub = np.vstack([M, -M])
    b_ub = np.concatenate([f + epsilon, -f + epsilon])
    result = linprog(c_obj, A_ub=A_ub, b_ub=b_ub, method='highs',
                    options={'presolve': True, 'time_limit': 30})
    return result.success

def find_min_epsilon(f, M):
    """Find minimum epsilon for given monomial matrix via LP."""
    n = M.shape[1]
    size = len(f)
    # Variables: [coeffs..., epsilon]
    # Minimize epsilon
    # Subject to: M*c - f <= epsilon, f - M*c <= epsilon, epsilon >= 0
    c_obj = np.zeros(n + 1)
    c_obj[-1] = 1.0  # minimize epsilon

    # M*c - epsilon*1 <= f  =>  [M, -1] * [c; eps] <= f
    # -M*c - epsilon*1 <= -f  =>  [-M, -1] * [c; eps] <= -f
    ones_col = -np.ones((size, 1))
    A_ub = np.vstack([
        np.hstack([M, ones_col]),
        np.hstack([-M, ones_col])
    ])
    b_ub = np.concatenate([f, -f])

    bounds = [(None, None)] * n + [(0, None)]  # epsilon >= 0

    result = linprog(c_obj, A_ub=A_ub, b_ub=b_ub, bounds=bounds,
                    method='highs', options={'presolve': True, 'time_limit': 60})

    if result.success:
        return result.x[-1]
    return None

def main():
    print("APPROXIMATE DEGREE OF PRIME INDICATOR chi_P(x)")
    print("=" * 70)
    print("chi_P(x) = 1 if x is prime, 0 otherwise")
    print("For each N and degree d, find min epsilon such that")
    print("|p(x) - chi_P(x)| <= epsilon for all x in {0,1}^N")
    print()

    results = {}

    for N in range(4, 11):
        f = truth_table_prime(N)
        size = 1 << N
        n_primes = int(np.sum(f))

        print(f"N={N}: {size} inputs, {n_primes} primes")

        results[N] = {}

        for d in range(0, N + 1):
            t0 = time.time()
            M = monomial_matrix(N, d)
            n_coeffs = M.shape[1]

            eps = find_min_epsilon(f, M)
            elapsed = time.time() - t0

            if eps is not None:
                results[N][d] = eps
                marker = ""
                if eps < 0.5:
                    marker = " <-- SUFFICIENT FOR ROUNDING"
                elif eps < 0.01:
                    marker = " <-- NEAR EXACT"
                print(f"  d={d:2d} ({n_coeffs:5d} monomials): eps={eps:.6f} [{elapsed:.1f}s]{marker}")
            else:
                print(f"  d={d:2d} ({n_coeffs:5d} monomials): LP FAILED [{elapsed:.1f}s]")

            if elapsed > 30 and d < N:
                print(f"  (skipping higher degrees for N={N})")
                break

            if eps is not None and eps < 1e-10:
                print(f"  (exact at d={d}, skipping higher)")
                break

        print()

    # Summary
    print("=" * 70)
    print("SUMMARY: Minimum degree for epsilon < 0.49 (rounding threshold)")
    print("=" * 70)
    for N in sorted(results.keys()):
        for d in sorted(results[N].keys()):
            if results[N][d] < 0.49:
                print(f"  N={N}: adeg(0.49) = {d}")
                break
        else:
            print(f"  N={N}: adeg(0.49) > {max(results[N].keys()) if results[N] else '?'}")

    print()
    print("SCALING:")
    print("If adeg(chi_P, 0.49) = Theta(N) → Omega(N) quantum queries (no speedup)")
    print("If adeg(chi_P, 0.49) = O(sqrt(N)) → quantum speedup possible")
    print("If adeg(chi_P, 0.49) = O(polylog(N)) → polylog quantum queries!")

if __name__ == '__main__':
    main()

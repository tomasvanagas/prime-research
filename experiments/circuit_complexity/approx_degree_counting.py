#!/usr/bin/env python3
"""
Approximate degree of pi(x) mod 2 — the counting parity function.
Session 12 proved pi(x) mod 2 is as hard as pi(x).
Compare with approximate degree of chi_P (the indicator).
"""
import numpy as np
from sympy import isprime
from scipy.optimize import linprog
from itertools import combinations
import time

def pi_mod2_table(N):
    """pi(x) mod 2 for all x in [0, 2^N)."""
    size = 1 << N
    table = np.zeros(size, dtype=np.float64)
    count = 0
    for x in range(size):
        if x >= 2 and isprime(x):
            count += 1
        table[x] = float(count % 2)
    return table

def prime_indicator_table(N):
    """chi_P(x) for all x in [0, 2^N)."""
    size = 1 << N
    return np.array([1.0 if (x >= 2 and isprime(x)) else 0.0 for x in range(size)])

def monomial_matrix(N, d):
    size = 1 << N
    columns = []
    for deg in range(d + 1):
        for S in combinations(range(N), deg):
            col = np.ones(size)
            for i in S:
                col *= np.array([(x >> i) & 1 for x in range(size)], dtype=np.float64)
            columns.append(col)
    return np.column_stack(columns)

def find_min_epsilon(f, M):
    n = M.shape[1]
    size = len(f)
    c_obj = np.zeros(n + 1)
    c_obj[-1] = 1.0
    ones_col = -np.ones((size, 1))
    A_ub = np.vstack([np.hstack([M, ones_col]), np.hstack([-M, ones_col])])
    b_ub = np.concatenate([f, -f])
    bounds = [(None, None)] * n + [(0, None)]
    result = linprog(c_obj, A_ub=A_ub, b_ub=b_ub, bounds=bounds,
                    method='highs', options={'presolve': True, 'time_limit': 60})
    return result.x[-1] if result.success else None

def main():
    print("APPROXIMATE DEGREE: chi_P vs pi(x) mod 2")
    print("=" * 70)

    for N in range(4, 11):
        f_ind = prime_indicator_table(N)
        f_par = pi_mod2_table(N)
        size = 1 << N

        print(f"\nN={N}: {size} inputs")
        print(f"  chi_P density: {np.mean(f_ind):.3f}")
        print(f"  pi mod 2 density: {np.mean(f_par):.3f}")
        print(f"  {'deg':>4} {'eps(chi_P)':>12} {'eps(pi mod 2)':>14}")
        print(f"  {'-'*34}")

        adeg_ind = None
        adeg_par = None

        for d in range(0, N + 1):
            t0 = time.time()
            M = monomial_matrix(N, d)

            eps_ind = find_min_epsilon(f_ind, M)
            eps_par = find_min_epsilon(f_par, M)
            elapsed = time.time() - t0

            if eps_ind is not None and eps_par is not None:
                mark_i = " *" if eps_ind < 0.49 and adeg_ind is None else ""
                mark_p = " *" if eps_par < 0.49 and adeg_par is None else ""
                if eps_ind < 0.49 and adeg_ind is None:
                    adeg_ind = d
                if eps_par < 0.49 and adeg_par is None:
                    adeg_par = d
                print(f"  {d:4d} {eps_ind:12.6f}{mark_i} {eps_par:14.6f}{mark_p}")

            if elapsed > 30:
                break
            if (eps_ind is not None and eps_ind < 1e-10 and
                eps_par is not None and eps_par < 1e-10):
                break

        print(f"\n  adeg(chi_P, 0.49) = {adeg_ind}")
        print(f"  adeg(pi mod 2, 0.49) = {adeg_par}")
        if adeg_ind and adeg_par:
            print(f"  Ratio: {adeg_par/adeg_ind:.2f}")

    print("\n" + "=" * 70)
    print("KEY QUESTION: Is pi(x) mod 2 harder than chi_P in approximate degree?")
    print("If adeg(pi mod 2) >> adeg(chi_P), the COUNTING step adds difficulty.")
    print("If adeg(pi mod 2) ≈ adeg(chi_P), counting is not the bottleneck.")

if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""Extract and inspect the actual annihilator polynomial of chi_P at low degree.

For each N, find the explicit smallest-degree polynomial g(x_0, ..., x_{N-1}) over F_2
such that g * chi_P = 0 (or g * (1 + chi_P) = 0). Inspect monomial support.
"""
import numpy as np
from sympy import isprime
from itertools import combinations
from algebraic_immunity_chi_p import (
    chi_P, build_truth_table, evaluate_monomials, monomials_up_to_degree, gauss_elim_F2
)


def find_annihilator(f_tt, N, target="f", max_d=None):
    """Find nontrivial annihilator g of f (or 1+f), of minimum degree.

    Returns (degree, c) where c is the coefficient vector indexed by
    monomials_up_to_degree(N, degree).
    """
    if max_d is None:
        max_d = N
    if target == "f":
        sub = np.nonzero(f_tt)[0]
    else:
        sub = np.nonzero(1 - f_tt)[0]
    for d in range(0, max_d + 1):
        masks = monomials_up_to_degree(N, d)
        V = evaluate_monomials(N, masks)
        A = V[sub].copy()
        # Find nullspace of A over F_2
        # Reduce columns to RREF
        n_rows, n_cols = A.shape
        # We want a c with A·c = 0.
        # Standard: row-reduce A, then read off free columns.
        B = A.copy()
        # Row reduction tracking pivot columns
        pivot_col = []
        pivot_row = 0
        for col in range(n_cols):
            if pivot_row >= n_rows:
                break
            idx = None
            for r in range(pivot_row, n_rows):
                if B[r, col] == 1:
                    idx = r
                    break
            if idx is None:
                continue
            if idx != pivot_row:
                B[[pivot_row, idx]] = B[[idx, pivot_row]]
            for r in range(n_rows):
                if r != pivot_row and B[r, col] == 1:
                    B[r] ^= B[pivot_row]
            pivot_col.append(col)
            pivot_row += 1
        rank = len(pivot_col)
        if rank < n_cols:
            # Nontrivial null space exists.
            free_cols = [c for c in range(n_cols) if c not in pivot_col]
            # Pick the first free column to construct kernel vector.
            free_c = free_cols[0]
            # Solve for pivot columns: B[i, pivot_col[i]] = 1, want c such that
            # B@c = 0. Set c[free_c] = 1, c[other free] = 0.
            c = np.zeros(n_cols, dtype=np.int8)
            c[free_c] = 1
            # For each pivot row i: B[i, pivot_col[i]] * c[pivot_col[i]] +
            # sum over j != pivot_col[i] of B[i, j] * c[j] = 0
            # → c[pivot_col[i]] = sum over free cols of B[i, free] * c[free]
            for i, pc in enumerate(pivot_col):
                # Bit-sum of B[i, free_c] (since other free cols have c=0)
                c[pc] = int(B[i, free_c])
            # Verify
            assert np.all((A @ c) % 2 == 0), "Kernel verify failed"
            return d, masks, c, target
    return None


def describe_annihilator(masks, c, N, max_show=12):
    """Pretty-print the annihilator polynomial."""
    terms = []
    for j, m in enumerate(masks):
        if c[j] == 1:
            if m == 0:
                terms.append("1")
            else:
                bits = [str(i) for i in range(N) if (m >> i) & 1]
                terms.append("x_{" + ",".join(bits) + "}")
    if not terms:
        return "0"
    if len(terms) > max_show:
        return " + ".join(terms[:max_show]) + f" + ... ({len(terms)} total terms)"
    return " + ".join(terms)


def verify_annihilator(masks, c, N, f_tt, target):
    """Verify g * f == 0 (or g * (1+f) == 0) over the full input space."""
    rng = 2 ** N
    g_eval = np.zeros(rng, dtype=np.int8)
    xs = np.arange(rng, dtype=np.int64)
    for j, m in enumerate(masks):
        if c[j] == 1:
            g_eval = (g_eval + ((xs & m) == m).astype(np.int8)) % 2
    if target == "f":
        prod = (g_eval * f_tt) % 2
    else:
        prod = (g_eval * (1 - f_tt)) % 2
    return np.all(prod == 0), int(np.sum(g_eval))


def main():
    for N in [5, 6, 7, 8, 9, 10, 11, 12, 13]:
        print(f"\n{'=' * 60}")
        print(f"N = {N}")
        print(f"{'=' * 60}")
        tt = build_truth_table(chi_P, N)
        n_primes = int(np.sum(tt))
        print(f"  primes in [0, {2**N}): {n_primes}, density = {n_primes / 2**N:.4f}")

        # Try ann of f first
        for target in ["f", "1+f"]:
            res = find_annihilator(tt, N, target=target, max_d=N)
            if res is None:
                print(f"  No annihilator of {target} found.")
                continue
            d, masks, c, _ = res
            ok, weight = verify_annihilator(masks, c, N, tt, target)
            terms_str = describe_annihilator(masks, c, N, max_show=20)
            n_terms = int(np.sum(c))
            print(f"  Annihilator of {target}: deg = {d}, {n_terms} monomials, "
                  f"weight |g^{-1}(1)| = {weight} / {2**N}")
            print(f"    g(x) = {terms_str}")
            print(f"    verified g * {'f' if target == 'f' else '(1+f)'} == 0: {ok}")
        # Now find the SMALLEST degree across both
        best = None
        for target in ["f", "1+f"]:
            res = find_annihilator(tt, N, target=target, max_d=N)
            if res is None:
                continue
            d = res[0]
            if best is None or d < best[0]:
                best = (d, target, res)
        if best is not None:
            d, target, _ = best
            print(f"  AI(chi_P, N={N}) = {d}  (annihilator of {target})")


if __name__ == "__main__":
    main()

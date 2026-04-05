#!/usr/bin/env python3
"""
NOF Discrepancy of the Prime Indicator Function chi_P.

k=3 party Number-On-Forehead model, balanced partition into N/3 bits each.
For N = 6, 9, 12:
  - Build tensor T[a][b][c] = (-1)^{chi_P(x)} where x = a + b*2^{N/3} + c*2^{2N/3}
  - Compute spectral norm (max singular value of unfoldings) -> upper bound on disc
  - Compute exact discrepancy for N=6 by full enumeration of cylinder intersections
  - Estimate discrepancy for N=9, 12 by random sampling + greedy optimization
  - Compute correlation with low-degree F_2 polynomials
  - Compare all measures with random functions of same density

This computation has never been done before. Previous work only computed
mode-unfolding ranks.
"""

import numpy as np
from itertools import product, combinations
import time
import sys
from functools import reduce

# ============================================================
# UTILITIES
# ============================================================

def sieve(n):
    """Sieve of Eratosthenes, returns boolean list."""
    if n < 2:
        return [False] * (n + 1)
    is_prime = [True] * (n + 1)
    is_prime[0] = is_prime[1] = False
    for i in range(2, int(n**0.5) + 1):
        if is_prime[i]:
            for j in range(i * i, n + 1, i):
                is_prime[j] = False
    return is_prime


def build_tensor(N):
    """
    Build the 3-party NOF tensor for N-bit inputs.
    T[a][b][c] = -1 if (a + b*2^{N/3} + c*2^{2N/3}) is prime, else +1.
    """
    k = 3
    n3 = N // k
    dim = 2 ** n3
    max_val = 2 ** N - 1

    primes = sieve(max_val)

    T = np.zeros((dim, dim, dim), dtype=np.int8)
    for a in range(dim):
        for b in range(dim):
            for c in range(dim):
                x = a + b * (2 ** n3) + c * (2 ** (2 * n3))
                if primes[x]:
                    T[a][b][c] = -1  # prime
                else:
                    T[a][b][c] = +1  # composite/0/1

    return T, primes


def tensor_stats(T, N):
    """Print basic tensor statistics."""
    n3 = N // 3
    dim = 2 ** n3
    total = dim ** 3
    n_prime = np.sum(T == -1)
    n_comp = np.sum(T == +1)
    density = n_prime / total
    print(f"  Dimension: {dim}x{dim}x{dim} = {total} entries")
    print(f"  Primes (T=-1): {n_prime}, Composites (T=+1): {n_comp}")
    print(f"  Prime density: {density:.6f}")
    print(f"  Bias (mean of T): {np.mean(T):.6f}")
    return density


# ============================================================
# 1. SPECTRAL NORM (upper bound on discrepancy)
# ============================================================

def spectral_norm_unfoldings(T):
    """
    Compute spectral norm (max singular value) of each mode-unfolding.
    For a 3-tensor T of shape (d1, d2, d3):
      Mode-0 unfolding: reshape to (d1, d2*d3) -> max SVD
      Mode-1 unfolding: reshape to (d2, d1*d3)
      Mode-2 unfolding: reshape to (d3, d1*d2)
    The spectral norm gives an upper bound on discrepancy.
    """
    d1, d2, d3 = T.shape
    total = d1 * d2 * d3

    norms = []
    for mode in range(3):
        if mode == 0:
            M = T.reshape(d1, d2 * d3).astype(np.float64)
        elif mode == 1:
            M = T.transpose(1, 0, 2).reshape(d2, d1 * d3).astype(np.float64)
        else:
            M = T.transpose(2, 0, 1).reshape(d3, d1 * d2).astype(np.float64)

        sv = np.linalg.svd(M, compute_uv=False)
        norms.append(sv[0])

    return norms


# ============================================================
# 2. EXACT DISCREPANCY (N=6 only: 4x4x4, feasible)
# ============================================================

def exact_discrepancy(T):
    """
    Enumerate all cylinder intersections C = A x B x C where
    A ⊆ [d1], B ⊆ [d2], C ⊆ [d3], all non-empty.
    disc(C) = |sum_{(a,b,c) in A x B x C} T[a,b,c]| / total
    Return maximum discrepancy.
    """
    d1, d2, d3 = T.shape
    total = d1 * d2 * d3

    # Precompute partial sums for efficiency
    # For each subset A of [d1], compute sum over a in A of T[a, :, :]
    # Then for each B, C, the sum is sum over b in B, c in C of partial[b, c]

    # Actually for 4x4x4, brute force over all 2^4-1 subsets per dim is fine
    # (2^4 - 1)^3 = 15^3 = 3375 combinations

    best_disc = 0.0
    best_sets = None

    # Generate all non-empty subsets of [d]
    def all_subsets(d):
        subs = []
        for mask in range(1, 2**d):
            s = [i for i in range(d) if mask & (1 << i)]
            subs.append(s)
        return subs

    subsA = all_subsets(d1)
    subsB = all_subsets(d2)
    subsC = all_subsets(d3)

    T_float = T.astype(np.float64)

    for A in subsA:
        T_A = T_float[A, :, :]  # shape (|A|, d2, d3)
        sum_A = T_A.sum(axis=0)  # shape (d2, d3)
        for B in subsB:
            sum_AB = sum_A[B, :].sum(axis=0)  # shape (d3,)
            for C in subsC:
                s = sum_AB[C].sum()
                disc = abs(s) / total
                if disc > best_disc:
                    best_disc = disc
                    best_sets = (A, B, C)

    return best_disc, best_sets


# ============================================================
# 3. ESTIMATED DISCREPANCY (random sampling + greedy)
# ============================================================

def random_cylinder_discrepancy(T, n_samples=1000000):
    """
    Estimate discrepancy by random sampling of cylinder intersections.
    Each sample: pick each A, B, C by including each element independently with prob 1/2.
    """
    d1, d2, d3 = T.shape
    total = d1 * d2 * d3
    T_float = T.astype(np.float64)

    best_disc = 0.0

    # Vectorized batch sampling
    batch = min(10000, n_samples)
    for start in range(0, n_samples, batch):
        actual_batch = min(batch, n_samples - start)

        # Random masks for each dimension
        maskA = np.random.randint(0, 2, size=(actual_batch, d1)).astype(np.float64)
        maskB = np.random.randint(0, 2, size=(actual_batch, d2)).astype(np.float64)
        maskC = np.random.randint(0, 2, size=(actual_batch, d3)).astype(np.float64)

        # Ensure non-empty: if all zeros, flip a random bit
        for masks, d in [(maskA, d1), (maskB, d2), (maskC, d3)]:
            empty = masks.sum(axis=1) == 0
            if empty.any():
                idx = np.where(empty)[0]
                flip = np.random.randint(0, d, size=len(idx))
                masks[idx, flip] = 1.0

        # Compute sums: sum_{a in A, b in B, c in C} T[a,b,c]
        # Step 1: contract with A: partial[s,j,k] = sum_i maskA[s,i] * T[i,j,k]
        partial = np.einsum('si,ijk->sjk', maskA, T_float)  # (batch, d2, d3)
        # Step 2: contract with B: partial2[s,k] = sum_j maskB[s,j] * partial[s,j,k]
        partial = np.einsum('sj,sjk->sk', maskB, partial)  # (batch, d3)
        # Step 3: contract with C
        sums = np.einsum('sk,sk->s', maskC, partial)  # (batch,)

        discs = np.abs(sums) / total
        batch_best = discs.max()
        if batch_best > best_disc:
            best_disc = batch_best

    return best_disc


def greedy_discrepancy(T, n_restarts=1000):
    """
    Greedy optimization: start with random A, B, C, then iteratively
    flip elements to maximize |sum|.
    """
    d1, d2, d3 = T.shape
    total = d1 * d2 * d3
    T_float = T.astype(np.float64)

    best_disc = 0.0
    best_sets = None

    for _ in range(n_restarts):
        # Random initial sets (each element included with prob 1/2)
        A = np.random.randint(0, 2, size=d1).astype(np.float64)
        B = np.random.randint(0, 2, size=d2).astype(np.float64)
        C = np.random.randint(0, 2, size=d3).astype(np.float64)

        # Ensure non-empty
        if A.sum() == 0: A[np.random.randint(d1)] = 1.0
        if B.sum() == 0: B[np.random.randint(d2)] = 1.0
        if C.sum() == 0: C[np.random.randint(d3)] = 1.0

        improved = True
        while improved:
            improved = False

            # Current sum
            cur_sum = np.einsum('i,j,k,ijk->', A, B, C, T_float)

            # Try flipping each element of A
            for i in range(d1):
                A[i] = 1.0 - A[i]
                if A.sum() == 0:
                    A[i] = 1.0 - A[i]
                    continue
                new_sum = np.einsum('i,j,k,ijk->', A, B, C, T_float)
                if abs(new_sum) > abs(cur_sum):
                    cur_sum = new_sum
                    improved = True
                else:
                    A[i] = 1.0 - A[i]

            # Try flipping each element of B
            for j in range(d2):
                B[j] = 1.0 - B[j]
                if B.sum() == 0:
                    B[j] = 1.0 - B[j]
                    continue
                new_sum = np.einsum('i,j,k,ijk->', A, B, C, T_float)
                if abs(new_sum) > abs(cur_sum):
                    cur_sum = new_sum
                    improved = True
                else:
                    B[j] = 1.0 - B[j]

            # Try flipping each element of C
            for k in range(d3):
                C[k] = 1.0 - C[k]
                if C.sum() == 0:
                    C[k] = 1.0 - C[k]
                    continue
                new_sum = np.einsum('i,j,k,ijk->', A, B, C, T_float)
                if abs(new_sum) > abs(cur_sum):
                    cur_sum = new_sum
                    improved = True
                else:
                    C[k] = 1.0 - C[k]

        disc = abs(cur_sum) / total
        if disc > best_disc:
            best_disc = disc
            best_sets = (A.copy(), B.copy(), C.copy())

    return best_disc, best_sets


# ============================================================
# 4. CORRELATION WITH LOW-DEGREE F_2 POLYNOMIALS
# ============================================================

def generate_f2_monomials(n_vars, max_degree):
    """
    Generate all monomials of degree <= max_degree over n_vars variables.
    Each monomial is a frozenset of variable indices.
    Returns list of frozensets.
    """
    monos = [frozenset()]  # degree 0: constant 1
    for d in range(1, max_degree + 1):
        for combo in combinations(range(n_vars), d):
            monos.append(frozenset(combo))
    return monos


def eval_monomial(mono, x_bits):
    """Evaluate a monomial (product of variables) on a bit vector."""
    result = 1
    for idx in mono:
        result &= x_bits[idx]
    return result


def correlation_with_f2_poly(T, N, max_degree=5):
    """
    For each degree d = 1, ..., max_degree:
      Find max over all degree-d polynomials p over F_2 of
      |E[(-1)^{chi_P(x) + p(x)}]|

    We use the fact that for F_2 polynomials, the correlation equals
    the maximum magnitude of the Walsh-Hadamard coefficients restricted
    to monomials of degree <= d.

    Specifically: E[(-1)^{f(x) + p(x)}] = sum of f-hat(S) * (-1)^{p includes S}
    The maximum over all degree-d p is achieved by aligning signs with
    the largest Walsh coefficients at degree <= d.

    Actually, simpler: the max correlation with degree-d F_2 poly equals
    sum_{|S| <= d} |f-hat(S)| where f-hat are Walsh coefficients.
    No wait -- the max over linear functions is max |f-hat(S)| for |S|=1.
    For general degree d, max correlation = max over degree-d p of |<f, (-1)^p>|.

    For exact computation: since p(x) = sum of monomials in some set M (each of degree <= d),
    (-1)^{p(x)} is NOT simply a product of (-1)^{monomials}. So we need to
    work in ±1 domain directly.

    Approach: Compute Walsh-Hadamard transform of f (in ±1 encoding).
    f-hat(S) = (1/2^N) sum_x f(x) (-1)^{<S,x>}
    The correlation with the degree-1 linear function L_S(x) = <S,x> mod 2 is exactly f-hat(S).
    The max correlation with degree-d polynomials requires more work.

    For feasible N (6, 9, 12), we can compute full Walsh spectrum.
    """
    n3 = N // 3
    dim = 2 ** n3
    total = 2 ** N

    # Flatten tensor to function f: {0,1}^N -> {-1,+1}
    f = np.zeros(total, dtype=np.float64)
    for a in range(dim):
        for b in range(dim):
            for c in range(dim):
                x = a + b * (2 ** n3) + c * (2 ** (2 * n3))
                f[x] = T[a, b, c]

    # Compute full Walsh-Hadamard transform
    # f-hat(S) = (1/2^N) sum_x f(x) (-1)^{popcount(S & x)}
    # Use fast Walsh-Hadamard transform
    fhat = f.copy()
    for i in range(N):
        half = 1 << i
        for j in range(0, total, 2 * half):
            for k in range(half):
                u = fhat[j + k]
                v = fhat[j + k + half]
                fhat[j + k] = u + v
                fhat[j + k + half] = u - v
    fhat /= total  # normalize

    # Now fhat[S] = Walsh coefficient for subset S (S encoded as bitmask)
    # |fhat[S]| for S with popcount(S) = d gives correlation with degree-1 chars

    # For degree-d correlation: max over all degree-d F_2 polynomials
    # A degree-d F_2 polynomial p(x) = XOR of monomials, each of degree <= d.
    # (-1)^{p(x)} is a product of (-1)^{monomial} values... but XOR makes this complex.
    #
    # Key insight: (-1)^{p(x)} = product over monomials m in support(p) of (-1)^{m(x)}
    # where m(x) = product of variables in m.
    # This is NOT the same as a single Walsh character.
    #
    # For degree 1: p(x) = sum of x_i (mod 2) = XOR of subset.
    # (-1)^{p(x)} = product (-1)^{x_i} = (-1)^{sum x_i} = (-1)^{<S,x>} = Walsh char.
    # So max degree-1 correlation = max |fhat[S]| over |S| = any (since <S,x> = XOR of bits in S).
    # Actually for degree 1, p(x) = c_0 + c_1 x_1 + ... + c_N x_N,
    # and (-1)^{p(x)} = (-1)^{c_0} * prod (-1)^{c_i x_i} = ±1 * Walsh char S = {i: c_i=1}.
    # So max degree-1 correlation = max |fhat[S]| over all S. This is the spectral norm in
    # a different sense.

    # For degree d > 1: we need to consider that (-1)^{p(x)} for degree-d p
    # is a multilinear function of degree d in the {-1,+1} variables y_i = (-1)^{x_i}.
    # The correlation E[f(x) * (-1)^{p(x)}] involves products of y_i's up to degree d.
    # This is exactly the sum of Walsh coefficients weighted by the p-structure.

    # The MAX correlation with ANY degree-d F_2 polynomial is bounded by:
    # sum_{|S| <= d} |fhat[S]|  (L1 norm of low-degree Walsh spectrum)
    # But the exact max is harder. For our purposes, we compute:
    # 1) max |fhat[S]| over |S| <= d (= correlation with single Walsh char of weight <= d)
    # 2) L1 norm of degree-d Walsh spectrum
    # 3) For small d and N, try to optimize directly

    results = {}

    # Organize Walsh coefficients by degree
    for d in range(min(N, max_degree) + 1):
        coeffs_at_d = []
        for S in range(total):
            if bin(S).count('1') == d:
                coeffs_at_d.append((S, fhat[S]))

        max_coeff = max(abs(c) for _, c in coeffs_at_d) if coeffs_at_d else 0
        l1_norm = sum(abs(c) for _, c in coeffs_at_d)
        l2_norm = np.sqrt(sum(c**2 for _, c in coeffs_at_d))

        results[d] = {
            'max_abs_coeff': max_coeff,
            'l1_norm': l1_norm,
            'l2_norm': l2_norm,
            'n_coeffs': len(coeffs_at_d),
        }

    # Cumulative: max correlation with degree-d poly (approximated by max single char)
    # and L1 bound
    cumulative = {}
    for d in range(min(N, max_degree) + 1):
        cum_max = max(results[dd]['max_abs_coeff'] for dd in range(d + 1))
        cum_l1 = sum(results[dd]['l1_norm'] for dd in range(d + 1))
        cumulative[d] = {
            'max_single_char': cum_max,
            'l1_bound': min(cum_l1, 1.0),  # correlation is at most 1
        }

    return results, cumulative, fhat


# ============================================================
# 5. RANDOM COMPARISON
# ============================================================

def build_random_tensor(dim, density):
    """Build a random ±1 tensor with given density of -1 entries."""
    T = np.ones((dim, dim, dim), dtype=np.int8)
    mask = np.random.random((dim, dim, dim)) < density
    T[mask] = -1
    return T


# ============================================================
# MAIN
# ============================================================

def main():
    np.random.seed(42)

    print("=" * 80)
    print("NOF DISCREPANCY OF THE PRIME INDICATOR FUNCTION chi_P")
    print("k=3 parties, balanced partition into N/3 bits each")
    print("=" * 80)

    N_values = [6, 9, 12]
    n_random = 10  # number of random comparisons

    all_results = {}

    for N in N_values:
        n3 = N // 3
        dim = 2 ** n3
        total = dim ** 3

        print(f"\n{'=' * 80}")
        print(f"N = {N}  (dim = {dim}^3 = {total}, range [0, {2**N - 1}])")
        print(f"{'=' * 80}")

        t0 = time.time()
        T, primes = build_tensor(N)
        density = tensor_stats(T, N)
        print(f"  Build time: {time.time() - t0:.3f}s")

        results = {'N': N, 'dim': dim, 'density': density}

        # --- Spectral Norm ---
        print(f"\n--- Spectral Norm (upper bound on discrepancy) ---")
        norms = spectral_norm_unfoldings(T)
        for mode, norm in enumerate(norms):
            disc_bound = norm / total
            print(f"  Mode-{mode} unfolding: spectral norm = {norm:.4f}, "
                  f"disc bound = {disc_bound:.6f}")
        results['spectral_norms'] = norms
        results['spectral_disc_bounds'] = [n / total for n in norms]

        # --- Exact Discrepancy (only for N=6) ---
        if N == 6:
            print(f"\n--- Exact Discrepancy (full enumeration) ---")
            t0 = time.time()
            exact_disc, exact_sets = exact_discrepancy(T)
            print(f"  Exact disc_3(chi_P) = {exact_disc:.6f}")
            print(f"  Achieved by A={exact_sets[0]}, B={exact_sets[1]}, C={exact_sets[2]}")

            # Show the actual values
            A, B, C = exact_sets
            print(f"  Cylinder sum = ", end="")
            s = sum(T[a, b, c] for a in A for b in B for c in C)
            print(f"{s} (out of {total})")
            print(f"  Enumeration time: {time.time() - t0:.3f}s")
            results['exact_disc'] = exact_disc

        # --- Estimated Discrepancy ---
        print(f"\n--- Estimated Discrepancy ---")

        # Random sampling
        n_samples = min(1000000, max(100000, total * 100))
        print(f"  Random sampling ({n_samples} samples)...")
        t0 = time.time()
        rand_disc = random_cylinder_discrepancy(T, n_samples=n_samples)
        print(f"    Best random disc = {rand_disc:.6f} ({time.time() - t0:.1f}s)")
        results['random_sample_disc'] = rand_disc

        # Greedy optimization
        n_restarts = min(5000, max(500, 1000))
        print(f"  Greedy optimization ({n_restarts} restarts)...")
        t0 = time.time()
        greedy_disc, greedy_sets = greedy_discrepancy(T, n_restarts=n_restarts)
        print(f"    Best greedy disc = {greedy_disc:.6f} ({time.time() - t0:.1f}s)")
        if greedy_sets is not None:
            A_idx = np.where(greedy_sets[0])[0]
            B_idx = np.where(greedy_sets[1])[0]
            C_idx = np.where(greedy_sets[2])[0]
            print(f"    Sets: A={list(A_idx)}, B={list(B_idx)}, C={list(C_idx)}")
        results['greedy_disc'] = greedy_disc

        # --- Walsh Spectrum / F_2 Polynomial Correlation ---
        max_deg = min(N // 3, 5)
        print(f"\n--- Correlation with degree-d F_2 polynomials (max_deg={max_deg}) ---")
        t0 = time.time()
        walsh_by_deg, cumulative, fhat = correlation_with_f2_poly(T, N, max_degree=max_deg)

        print(f"  {'Degree':>6} | {'#Coeffs':>8} | {'Max |f-hat|':>12} | "
              f"{'L1 norm':>10} | {'L2 norm':>10} | {'Cum Max':>10} | {'Cum L1':>10}")
        print(f"  {'-'*6}-+-{'-'*8}-+-{'-'*12}-+-{'-'*10}-+-{'-'*10}-+-{'-'*10}-+-{'-'*10}")
        for d in range(max_deg + 1):
            wd = walsh_by_deg[d]
            cd = cumulative[d]
            print(f"  {d:>6} | {wd['n_coeffs']:>8} | {wd['max_abs_coeff']:>12.6f} | "
                  f"{wd['l1_norm']:>10.6f} | {wd['l2_norm']:>10.6f} | "
                  f"{cd['max_single_char']:>10.6f} | {cd['l1_bound']:>10.6f}")
        print(f"  Walsh computation time: {time.time() - t0:.3f}s")
        results['walsh'] = walsh_by_deg
        results['walsh_cumulative'] = cumulative

        # --- Top Walsh coefficients ---
        print(f"\n  Top 10 Walsh coefficients by magnitude:")
        sorted_idx = np.argsort(np.abs(fhat))[::-1]
        for rank, S in enumerate(sorted_idx[:10]):
            bits = bin(S)[2:].zfill(N)
            weight = bin(S).count('1')
            print(f"    #{rank+1}: S={bits} (wt={weight}), f-hat(S) = {fhat[S]:+.6f}")

        # --- Comparison with Random Functions ---
        print(f"\n--- Comparison with {n_random} random functions of same density ---")
        random_spectral = []
        random_discs = []
        random_greedy = []
        random_walsh_max = {d: [] for d in range(max_deg + 1)}

        for trial in range(n_random):
            Tr = build_random_tensor(dim, density)

            # Spectral norm
            rn = spectral_norm_unfoldings(Tr)
            random_spectral.append([n / total for n in rn])

            # Random sampling disc
            rd = random_cylinder_discrepancy(Tr, n_samples=min(100000, n_samples // 10))
            random_discs.append(rd)

            # Greedy disc
            gd, _ = greedy_discrepancy(Tr, n_restarts=min(200, n_restarts // 5))
            random_greedy.append(gd)

            # Walsh
            rw, rc, _ = correlation_with_f2_poly(Tr, N, max_degree=max_deg)
            for d in range(max_deg + 1):
                random_walsh_max[d].append(rw[d]['max_abs_coeff'])

        print(f"\n  {'Measure':>30} | {'chi_P':>10} | {'Rand Mean':>10} | {'Rand Std':>10} | {'Ratio':>8}")
        print(f"  {'-'*30}-+-{'-'*10}-+-{'-'*10}-+-{'-'*10}-+-{'-'*8}")

        # Spectral disc bound (use min across modes as tightest bound)
        chi_spec = min(results['spectral_disc_bounds'])
        rand_spec_vals = [min(r) for r in random_spectral]
        print(f"  {'Spectral disc bound (min mode)':>30} | {chi_spec:>10.6f} | "
              f"{np.mean(rand_spec_vals):>10.6f} | {np.std(rand_spec_vals):>10.6f} | "
              f"{chi_spec / np.mean(rand_spec_vals):>8.3f}")

        # Greedy disc
        chi_greedy = results['greedy_disc']
        print(f"  {'Greedy disc':>30} | {chi_greedy:>10.6f} | "
              f"{np.mean(random_greedy):>10.6f} | {np.std(random_greedy):>10.6f} | "
              f"{chi_greedy / np.mean(random_greedy):>8.3f}")

        # Random sample disc
        chi_rand = results['random_sample_disc']
        print(f"  {'Random sample disc':>30} | {chi_rand:>10.6f} | "
              f"{np.mean(random_discs):>10.6f} | {np.std(random_discs):>10.6f} | "
              f"{chi_rand / np.mean(random_discs):>8.3f}")

        # Walsh by degree
        for d in range(max_deg + 1):
            chi_w = walsh_by_deg[d]['max_abs_coeff']
            rw_vals = random_walsh_max[d]
            label = f'Max |Walsh coeff| deg {d}'
            ratio = chi_w / np.mean(rw_vals) if np.mean(rw_vals) > 0 else float('inf')
            print(f"  {label:>30} | {chi_w:>10.6f} | "
                  f"{np.mean(rw_vals):>10.6f} | {np.std(rw_vals):>10.6f} | "
                  f"{ratio:>8.3f}")

        all_results[N] = results
        print()

    # ============================================================
    # SUMMARY
    # ============================================================
    print("\n" + "=" * 80)
    print("SUMMARY: DISCREPANCY OF chi_P vs RANDOM")
    print("=" * 80)

    print(f"\n{'N':>4} | {'Exact Disc':>12} | {'Greedy Disc':>12} | "
          f"{'Spec Bound':>12} | {'Rand Greedy':>12} | {'Ratio G/R':>10}")
    print(f"{'-'*4}-+-{'-'*12}-+-{'-'*12}-+-{'-'*12}-+-{'-'*12}-+-{'-'*10}")

    for N in N_values:
        r = all_results[N]
        exact = f"{r.get('exact_disc', float('nan')):.6f}" if 'exact_disc' in r else "N/A"
        greedy = f"{r['greedy_disc']:.6f}"
        spec = f"{min(r['spectral_disc_bounds']):.6f}"
        print(f"{N:>4} | {exact:>12} | {greedy:>12} | {spec:>12} | "
              f"{'(see above)':>12} | {'(see above)':>10}")

    print(f"""
INTERPRETATION:
- Discrepancy < 2^{{-k}} would imply Omega(k) communication complexity.
- For k=3, disc < 1/8 = 0.125 would give Omega(3) = non-trivial lower bound.
- If chi_P has LOWER discrepancy than random functions of same density,
  this suggests additional structure that makes communication harder.
- If chi_P has HIGHER discrepancy, the prime indicator might be "easier"
  to approximate by cylinder intersections (consistent with being in ACC^0/TC^0).
- Walsh correlation with low-degree F_2 polynomials measures distance from AC^0[2].
""")


if __name__ == '__main__':
    main()

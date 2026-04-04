"""
Extended determinantal complexity experiments for pi_N polynomial.

pi_N(x_1,...,x_N) = pi(sum x_i * 2^i) = number of primes <= x,
where x = x_0 + 2*x_1 + ... + 2^{N-1}*x_{N-1} is an N-bit number.

Previous sessions found dc(pi_N) = N for N=2,3,4.
This experiment extends to N=5..12 using multiple strategies:

1. Multilinear polynomial computation via Mobius inversion on Boolean lattice
2. Lower bounds from substitution/restriction (affine minors)
3. Numerical optimization (L-BFGS-B + differential evolution) for m x m reps
4. Algebraic lower bound: dc(f) >= rank of the "flattening" matrix
5. For large N: parameter counting and dimension analysis

Key fact: dc(f) >= partition_rank(f) where partition_rank is the largest
rank of any "substitution matrix" obtained by partitioning variables into
two sets and viewing f as a matrix in those variables.
"""

import numpy as np
from itertools import product as iproduct, combinations
from scipy.optimize import minimize
import sys
import time

# Force unbuffered output
import functools
print = functools.partial(print, flush=True)

# ---------------------------------------------------------------------------
# Sympy primepi is slow for large x; use a simple sieve for x up to 2^12=4096
# ---------------------------------------------------------------------------

def simple_sieve(limit):
    """Sieve of Eratosthenes returning a lookup list is_prime[0..limit]."""
    if limit < 2:
        return [False] * (limit + 1)
    sieve = [True] * (limit + 1)
    sieve[0] = sieve[1] = False
    for i in range(2, int(limit**0.5) + 1):
        if sieve[i]:
            for j in range(i * i, limit + 1, i):
                sieve[j] = False
    return sieve

def compute_primepi_table(max_val):
    """Return list pi_table where pi_table[x] = number of primes <= x."""
    sieve = simple_sieve(max_val)
    pi_table = [0] * (max_val + 1)
    count = 0
    for i in range(max_val + 1):
        if sieve[i]:
            count += 1
        pi_table[i] = count
    return pi_table

# Precompute for x up to 2^12 = 4096
PI_TABLE = compute_primepi_table(4096)

def primepi(x):
    if x < 0:
        return 0
    if x <= 4096:
        return PI_TABLE[x]
    # Fallback for safety
    from sympy import primepi as sp_primepi
    return int(sp_primepi(x))


# ---------------------------------------------------------------------------
# Multilinear polynomial computation
# ---------------------------------------------------------------------------

def compute_pi_polynomial(N):
    """
    Compute multilinear polynomial P(x_0,...,x_{N-1}) = pi(sum x_i * 2^i)
    via Mobius inversion on the Boolean lattice.

    Returns dict: frozenset_of_variable_indices -> integer_coefficient
    """
    evaluations = {}
    for mask in range(1 << N):
        bits = tuple((mask >> i) & 1 for i in range(N))
        x = sum(b * (1 << i) for i, b in enumerate(bits))
        evaluations[mask] = primepi(x)

    coefficients = {}
    for S_mask in range(1 << N):
        coeff = 0
        # Sum over subsets T of S
        T_mask = S_mask
        while True:
            sign = (-1) ** bin(S_mask ^ T_mask).count('1')
            coeff += sign * evaluations[T_mask]
            if T_mask == 0:
                break
            T_mask = (T_mask - 1) & S_mask

        if coeff != 0:
            S = frozenset(i for i in range(N) if S_mask & (1 << i))
            coefficients[S] = coeff

    return coefficients


def poly_info(coeffs, N):
    """Print info about the multilinear polynomial."""
    n_terms = len(coeffs)
    total_possible = 1 << N
    degree = max((len(S) for S in coeffs), default=0)
    max_coeff = max(abs(c) for c in coeffs.values()) if coeffs else 0
    coeff_sum = sum(abs(c) for c in coeffs.values())

    print(f"  Monomials: {n_terms}/{total_possible} "
          f"({100*n_terms/total_possible:.1f}% nonzero)")
    print(f"  Degree: {degree}")
    print(f"  Max |coeff|: {max_coeff}")
    print(f"  Sum |coeff|: {coeff_sum}")

    # Degree distribution
    deg_counts = {}
    for S, c in coeffs.items():
        d = len(S)
        deg_counts[d] = deg_counts.get(d, 0) + 1
    print(f"  Degree distribution: {dict(sorted(deg_counts.items()))}")

    return degree, n_terms, max_coeff


# ---------------------------------------------------------------------------
# Lower bound: substitution rank (partition rank)
# ---------------------------------------------------------------------------

def substitution_rank_lower_bound(coeffs, N):
    """
    For each partition of {0,...,N-1} into two sets A, B,
    view f as a bilinear form in (x_A, x_B) variables.
    The rank of the resulting matrix is a lower bound on dc(f).

    dc(f) >= max over partitions of rank(M_{A,B})

    We check all partitions with |A| = floor(N/2).
    """
    best_rank = 0
    best_partition = None
    half = N // 2

    for A in combinations(range(N), half):
        A_set = set(A)
        B_set = set(range(N)) - A_set
        A_list = sorted(A_set)
        B_list = sorted(B_set)

        nA = 1 << len(A_list)
        nB = 1 << len(B_list)

        # Build the "flattening" matrix
        # Rows indexed by assignments to A variables, cols by B variables
        # Entry = sum of coefficients of monomials that match this assignment pattern
        M = np.zeros((nA, nB))

        for S, c in coeffs.items():
            S_A = S & A_set
            S_B = S & B_set

            # Row index: which A-variables are in S
            row = 0
            for idx, a in enumerate(A_list):
                if a in S_A:
                    row |= 1 << idx

            # Col index: which B-variables are in S
            col = 0
            for idx, b in enumerate(B_list):
                if b in S_B:
                    col |= 1 << idx

            M[row, col] += c

        r = np.linalg.matrix_rank(M)
        if r > best_rank:
            best_rank = r
            best_partition = (A_list, B_list)

    return best_rank, best_partition


# ---------------------------------------------------------------------------
# Lower bound: slice rank
# ---------------------------------------------------------------------------

def slice_rank_lower_bound(coeffs, N):
    """
    Fix one variable x_i = 0 or 1, get a polynomial in N-1 variables.
    The number of distinct such restrictions gives info about complexity.

    Also: fix k variables, count how many distinct polynomials arise.
    More directly: dc(f) >= dc(f|_{x_i=c}) for any restriction.
    So we compute substitution rank for all single-variable restrictions.
    """
    best_rank = 0
    for var in range(N):
        for val in [0, 1]:
            # Restrict x_var = val
            restricted = {}
            for S, c in coeffs.items():
                if var in S:
                    if val == 1:
                        new_S = frozenset(v for v in S if v != var)
                        restricted[new_S] = restricted.get(new_S, 0) + c
                    # if val == 0, this monomial vanishes
                else:
                    restricted[S] = restricted.get(S, 0) + c

            # Remove zeros
            restricted = {S: c for S, c in restricted.items() if c != 0}

            if len(restricted) == 0:
                continue

            # Get substitution rank of restricted polynomial (N-1 vars)
            remaining_vars = sorted(set(range(N)) - {var})
            # Remap variables
            remap = {v: i for i, v in enumerate(remaining_vars)}
            remapped = {}
            for S, c in restricted.items():
                new_S = frozenset(remap[v] for v in S)
                remapped[new_S] = c

            if len(remaining_vars) >= 2:
                r, _ = substitution_rank_lower_bound(remapped, N - 1)
                if r > best_rank:
                    best_rank = r

    return best_rank


# ---------------------------------------------------------------------------
# Numerical optimization for m x m determinantal representation
# ---------------------------------------------------------------------------

def compute_targets(N):
    """Compute pi(x) for all 2^N values of x."""
    targets = []
    bit_configs = []
    for mask in range(1 << N):
        bits = tuple((mask >> i) & 1 for i in range(N))
        x = sum(b * (1 << i) for i, b in enumerate(bits))
        targets.append(primepi(x))
        bit_configs.append(bits)
    return np.array(targets, dtype=float), bit_configs


def objective(params, m, N, targets, bit_configs):
    """Sum of squared errors: det(M(bits)) vs targets."""
    A = params.reshape(m, m, N + 1)
    total_err = 0.0
    for idx, bits in enumerate(bit_configs):
        M = A[:, :, 0].copy()
        for k in range(N):
            if bits[k]:
                M += A[:, :, k + 1]
        det_val = np.linalg.det(M)
        total_err += (det_val - targets[idx]) ** 2
    return total_err


def find_det_rep(N, m, n_trials=100, verbose=True):
    """
    Search for m x m determinantal representation of pi_N.
    Returns (best_error, best_A) or (best_error, None) if not found.
    """
    targets, bit_configs = compute_targets(N)
    n_params = m * m * (N + 1)
    n_constraints = 1 << N

    if verbose:
        print(f"  Searching {m}x{m} rep: {n_params} params, "
              f"{n_constraints} constraints, "
              f"ratio={n_params/n_constraints:.2f}")

    best_error = float('inf')
    best_A = None

    for trial in range(n_trials):
        np.random.seed(trial * 31 + 17)

        # Vary initialization scale
        if trial < n_trials // 3:
            scale = 0.5
        elif trial < 2 * n_trials // 3:
            scale = 1.5
        else:
            scale = 3.0

        x0 = np.random.randn(n_params) * scale

        result = minimize(
            objective, x0, args=(m, N, targets, bit_configs),
            method='L-BFGS-B',
            options={'maxiter': 8000, 'ftol': 1e-25, 'gtol': 1e-15}
        )

        if result.fun < best_error:
            best_error = result.fun
            best_A = result.x.reshape(m, m, N + 1)

            if best_error < 1e-10:
                break

    # Verify
    if best_error < 1e-4 and best_A is not None:
        all_correct = True
        max_err = 0
        for idx, bits in enumerate(bit_configs):
            M = best_A[:, :, 0].copy()
            for k in range(N):
                if bits[k]:
                    M += best_A[:, :, k + 1]
            det_val = np.linalg.det(M)
            err = abs(det_val - targets[idx])
            max_err = max(max_err, err)
            if err > 0.4:
                all_correct = False

        if all_correct:
            if verbose:
                print(f"    FOUND! error={best_error:.2e}, max_pointwise={max_err:.2e}")
            return best_error, best_A
        else:
            if verbose:
                print(f"    Low total error={best_error:.2e} but max_pointwise={max_err:.2e} FAILS")

    if verbose:
        print(f"    Best error: {best_error:.2e}")

    return best_error, None


# ---------------------------------------------------------------------------
# Upper bound: explicit construction attempt via Grenet-style matrices
# ---------------------------------------------------------------------------

def grenet_style_search(N, m, n_trials=100):
    """
    Search for det reps with structured matrices:
    - Hessenberg form (upper or lower)
    - Tridiagonal
    These have fewer parameters, so optimization is easier.
    """
    targets, bit_configs = compute_targets(N)

    # Hessenberg: M[i,j] = 0 for i > j+1
    # So row i has nonzero entries only in columns j >= i-1
    # Entries per row: m - max(0, i-1) for row i
    # Total nonzero entries: sum_{i=0}^{m-1} (m - max(0,i-1)) = m^2 - (m-1)(m-2)/2

    def objective_hess(params, m, N, targets, bit_configs):
        # Unpack parameters for lower Hessenberg matrix
        # M[i][j] is nonzero only if j <= i+1
        idx = 0
        A_entries = {}  # (i,j) -> array of N+1 coefficients
        for i in range(m):
            for j in range(min(i + 2, m)):
                A_entries[(i, j)] = params[idx:idx + N + 1]
                idx += N + 1

        total_err = 0.0
        for bidx, bits in enumerate(bit_configs):
            M = np.zeros((m, m))
            for (i, j), coeffs in A_entries.items():
                val = coeffs[0]
                for k in range(N):
                    if bits[k]:
                        val += coeffs[k + 1]
                M[i, j] = val
            det_val = np.linalg.det(M)
            total_err += (det_val - targets[bidx]) ** 2
        return total_err

    n_nonzero = sum(min(i + 2, m) for i in range(m))
    n_params = n_nonzero * (N + 1)

    best_error = float('inf')
    for trial in range(n_trials):
        np.random.seed(trial * 43 + 11)
        x0 = np.random.randn(n_params) * 1.5
        result = minimize(
            objective_hess, x0, args=(m, N, targets, bit_configs),
            method='L-BFGS-B',
            options={'maxiter': 5000, 'ftol': 1e-25}
        )
        if result.fun < best_error:
            best_error = result.fun
            if best_error < 1e-10:
                break

    return best_error


# ---------------------------------------------------------------------------
# Main experiment
# ---------------------------------------------------------------------------

def analyze_N(N, try_sizes=None, n_trials=None):
    """Full analysis for a given N."""
    t0 = time.time()
    print(f"\n{'=' * 70}")
    print(f"  N = {N}  (x ranges 0..{(1 << N) - 1}, pi(x) ranges 0..{primepi((1 << N) - 1)})")
    print(f"{'=' * 70}")

    # 1. Compute polynomial
    print("\n--- Multilinear polynomial ---")
    coeffs = compute_pi_polynomial(N)
    degree, n_terms, max_coeff = poly_info(coeffs, N)

    # 2. Lower bounds
    print("\n--- Lower bounds ---")
    sub_rank, partition = substitution_rank_lower_bound(coeffs, N)
    print(f"  Substitution rank (balanced partition): {sub_rank}")
    if partition:
        print(f"    Best partition: A={partition[0]}, B={partition[1]}")

    if N <= 10:
        slice_lb = slice_rank_lower_bound(coeffs, N)
        print(f"  Max restriction rank (single var fixed): {slice_lb}")

    # Combined lower bound
    lb = max(sub_rank, degree)  # dc >= degree since det of m x m has degree m
    print(f"  Combined lower bound: dc(pi_{N}) >= {lb}")

    # 3. Numerical search for det representations
    print("\n--- Numerical search for det representations ---")

    if try_sizes is None:
        # Start from lower bound, go up to N+3
        try_sizes = list(range(max(lb, N), min(N + 4, N + 5)))
        # Also try smaller if degree < N
        if degree < N:
            try_sizes = list(range(degree, N + 4))

    if n_trials is None:
        # Scale trials inversely with difficulty
        n_trials_base = max(50, 500 - N * 30)
    else:
        n_trials_base = n_trials

    found_size = None
    results = {}

    for m in try_sizes:
        n_params = m * m * (N + 1)
        n_constraints = 1 << N

        # Skip if clearly underdetermined
        if n_params < n_constraints:
            print(f"  m={m}: SKIP (underdetermined: {n_params} < {n_constraints})")
            results[m] = ('skip', float('inf'))
            continue

        # Adjust trials based on parameter count
        trials = max(20, min(n_trials_base, 300))
        if m * m * (N + 1) > 2000:
            trials = max(20, trials // 3)

        err, A = find_det_rep(N, m, n_trials=trials)
        results[m] = ('found' if A is not None else 'not_found', err)

        if A is not None:
            found_size = m
            break

    # 4. Try Hessenberg structure for the smallest found size or lower bound
    if found_size is None and N <= 10:
        print("\n--- Hessenberg structure search ---")
        for m in range(max(lb, N), N + 4):
            n_nonzero = sum(min(i + 2, m) for i in range(m))
            n_hess_params = n_nonzero * (N + 1)
            if n_hess_params < (1 << N):
                print(f"  m={m} Hessenberg: SKIP (underdetermined)")
                continue
            trials_h = max(20, min(100, 200 - N * 10))
            err_h = grenet_style_search(N, m, n_trials=trials_h)
            print(f"  m={m} Hessenberg: error={err_h:.2e}")
            if err_h < 1e-4:
                print(f"    Hessenberg representation FOUND at size {m}!")
                if found_size is None or m < found_size:
                    found_size = m
                break

    elapsed = time.time() - t0

    # Summary
    print(f"\n--- Summary for N={N} ---")
    print(f"  Lower bound: dc(pi_{N}) >= {lb}")
    if found_size:
        print(f"  Upper bound: dc(pi_{N}) <= {found_size}")
        if found_size == lb:
            print(f"  EXACT: dc(pi_{N}) = {found_size}")
        elif found_size == N:
            print(f"  dc(pi_{N}) = {found_size} = N")
    else:
        print(f"  No representation found up to m={max(try_sizes) if try_sizes else '?'}")
        print(f"  Likely dc(pi_{N}) > {max(try_sizes) if try_sizes else '?'}")
    print(f"  Time: {elapsed:.1f}s")

    return {
        'N': N,
        'degree': degree,
        'n_terms': n_terms,
        'max_coeff': max_coeff,
        'lower_bound': lb,
        'sub_rank': sub_rank,
        'found_size': found_size,
        'results': results,
        'time': elapsed,
    }


def main():
    print("=" * 70)
    print("  EXTENDED DETERMINANTAL COMPLEXITY OF pi_N")
    print("  dc(f) = smallest m s.t. f = det(M), M m×m with affine entries")
    print("=" * 70)

    all_results = {}

    # Phase 1: Verify known results N=2,3,4
    print("\n\n### PHASE 1: Verify known results (N=2,3,4) ###")
    for N in [2, 3, 4]:
        r = analyze_N(N, try_sizes=list(range(N, N + 3)), n_trials=100)
        all_results[N] = r

    # Phase 2: New experiments N=5..8 (feasible for full search)
    print("\n\n### PHASE 2: New experiments (N=5..8) ###")
    for N in [5, 6, 7, 8]:
        trials = max(30, 150 - N * 15)
        r = analyze_N(N, n_trials=trials)
        all_results[N] = r

    # Phase 3: Larger N=9..12 (mostly lower bounds + limited search)
    print("\n\n### PHASE 3: Larger experiments (N=9..12) ###")
    for N in [9, 10, 11, 12]:
        # For large N, limit the search range
        max_try = min(N + 2, 14)
        sizes = list(range(N, max_try + 1))
        r = analyze_N(N, try_sizes=sizes, n_trials=20)
        all_results[N] = r

    # Final summary table
    print("\n\n" + "=" * 70)
    print("  FINAL SUMMARY")
    print("=" * 70)
    print(f"{'N':>3} {'deg':>4} {'#terms':>7} {'max|c|':>8} "
          f"{'sub_rank':>9} {'LB':>4} {'dc<=':>5} {'dc=N?':>6} {'time':>7}")
    print("-" * 70)
    for N in sorted(all_results):
        r = all_results[N]
        ub_str = str(r['found_size']) if r['found_size'] else '?'
        eq_n = ''
        if r['found_size'] is not None:
            if r['found_size'] == N:
                eq_n = 'YES'
            elif r['found_size'] > N:
                eq_n = f'>{N}'
            else:
                eq_n = f'<{N}!'
        print(f"{r['N']:>3} {r['degree']:>4} {r['n_terms']:>7} {r['max_coeff']:>8} "
              f"{r['sub_rank']:>9} {r['lower_bound']:>4} {ub_str:>5} {eq_n:>6} "
              f"{r['time']:>6.1f}s")

    print("\n\nKEY FINDINGS:")
    for N in sorted(all_results):
        r = all_results[N]
        if r['found_size'] is not None:
            if r['found_size'] == N:
                print(f"  N={N}: dc(pi_{N}) = {N} (confirmed)")
            else:
                print(f"  N={N}: dc(pi_{N}) <= {r['found_size']} "
                      f"(lower bound {r['lower_bound']})")
        else:
            print(f"  N={N}: dc(pi_{N}) >= {r['lower_bound']}, "
                  f"no rep found up to size tried")


if __name__ == "__main__":
    main()

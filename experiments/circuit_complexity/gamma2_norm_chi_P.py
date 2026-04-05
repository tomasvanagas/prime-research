#!/usr/bin/env python3
"""
Gamma-2 norm (factorization norm) and sign-rank of the communication matrix
of the prime indicator function chi_P.

For N-bit inputs, split x = (a, b) where a = top N/2 bits, b = bottom N/2 bits.
Matrix M[a][b] = (-1)^{chi_P(x)} where x = a * 2^{N/2} + b.
  +1 if composite/0/1, -1 if prime.

Computes for N = 4, 6, 8, 10, 12 (matrix sizes 4x4 to 64x64):
  1. Rank over R, rank over F_2
  2. Full SVD, nuclear/spectral/Frobenius norms
  3. Gamma-2 norm via SDP (cvxpy)
  4. Sign-rank bounds (Forster lower bound, explicit search for small N)
  5. Comparison with random ±1 matrices of same size and bias
"""

import numpy as np
from sympy import isprime
import scipy.linalg as la
import time
import warnings
warnings.filterwarnings('ignore')

try:
    import cvxpy as cp
    HAS_CVXPY = True
except ImportError:
    HAS_CVXPY = False
    print("WARNING: cvxpy not available. Gamma-2 will use approximations only.")


def build_prime_matrix(N):
    """Build the ±1 communication matrix for N-bit inputs."""
    half = N // 2
    rows = 2 ** half  # number of values for top half
    cols = 2 ** half  # number of values for bottom half
    M = np.ones((rows, cols), dtype=float)
    prime_count = 0
    for a in range(rows):
        for b in range(cols):
            x = a * (2 ** half) + b
            if isprime(x):
                M[a][b] = -1.0
                prime_count += 1
    return M, prime_count


def matrix_rank_F2(M_pm1):
    """Rank over F_2 of the {0,1} version (where -1 -> 1, +1 -> 0)."""
    rows, cols = M_pm1.shape
    # Convert: -1 -> 1, +1 -> 0
    B = ((1 - M_pm1) / 2).astype(int) % 2
    # Gaussian elimination over F_2
    B = B.copy()
    rank = 0
    for col in range(cols):
        # Find pivot
        pivot = None
        for row in range(rank, rows):
            if B[row, col] == 1:
                pivot = row
                break
        if pivot is None:
            continue
        # Swap
        B[[rank, pivot]] = B[[pivot, rank]]
        # Eliminate
        for row in range(rows):
            if row != rank and B[row, col] == 1:
                B[row] = (B[row] + B[rank]) % 2
        rank += 1
    return rank


def compute_svd_properties(M):
    """Compute SVD and derived norms."""
    U, s, Vt = la.svd(M, full_matrices=False)
    rank_R = np.sum(s > 1e-10)
    nuclear = np.sum(s)
    spectral = s[0]
    frobenius = np.sqrt(np.sum(s**2))
    return {
        'singular_values': s,
        'rank_R': int(rank_R),
        'nuclear_norm': nuclear,
        'spectral_norm': spectral,
        'frobenius_norm': frobenius,
    }


def gamma2_norm_sdp(M):
    """
    Compute gamma-2 norm via SDP using cvxpy.

    gamma_2(M) = min_{M=UV^T} (max_i ||u_i||_2) * (max_j ||v_j||_2)

    Equivalent SDP (Linial-Shraibman, Lee-Shraibman):
    gamma_2(M) = min sqrt(alpha * beta) such that:
        X is m x m PSD, Y is n x n PSD,
        [[X, M], [M^T, Y]] >> 0,
        X[i,i] <= alpha for all i,
        Y[j,j] <= beta for all j.

    Setting alpha = beta = t (by symmetry optimization):
    gamma_2(M) = min t such that:
        X is m x m PSD, Y is n x n PSD,  (FULL PSD, not diagonal!)
        [[X, M], [M^T, Y]] >> 0,
        X[i,i] <= t for all i,
        Y[j,j] <= t for all j.

    Key: X and Y are FULL PSD matrices with bounded diagonal entries,
    NOT diagonal matrices. Using diagonal matrices gives the trivially
    loose bound gamma_2 = spectral norm.
    """
    if not HAS_CVXPY:
        return None

    m, n = M.shape

    # Variables: full PSD matrices X (m x m) and Y (n x n)
    X = cp.Variable((m, m), symmetric=True)
    Y = cp.Variable((n, n), symmetric=True)
    t = cp.Variable(nonneg=True)

    # Build the block matrix [[X, M], [M^T, Y]]
    block = cp.bmat([
        [X, M],
        [M.T, Y]
    ])

    constraints = [
        block >> 0,  # PSD
        X >> 0,      # X is PSD
        Y >> 0,      # Y is PSD
    ]

    # Diagonal entries bounded by t
    for i in range(m):
        constraints.append(X[i, i] <= t)
    for j in range(n):
        constraints.append(Y[j, j] <= t)

    prob = cp.Problem(cp.Minimize(t), constraints)
    try:
        prob.solve(solver=cp.SCS, verbose=False, max_iters=20000, eps=1e-8)
        if prob.status in ['optimal', 'optimal_inaccurate']:
            return float(t.value)
        else:
            # Try another solver
            prob.solve(solver=cp.CLARABEL, verbose=False)
            if prob.status in ['optimal', 'optimal_inaccurate']:
                return float(t.value)
            return None
    except Exception as ex:
        print(f"  SDP solver error: {ex}")
        return None


def gamma2_bounds_from_svd(M, svd_props):
    """
    Bounds on gamma_2 from SVD properties.

    For an m x n ±1 matrix M:
    - Lower bound: gamma_2(M) >= ||M||_F^2 / (m * ||M||_2) = n / ||M||_2
      (since ||M||_F^2 = m*n for ±1 matrices)
      This comes from: for any factorization M = UV^T,
      ||M||_F^2 <= (max||u_i||^2 * m) * (max||v_j||^2 * ... ) -- not tight.

      Better: gamma_2(M) >= nuclear_norm / min(m,n)
      (Linial-Shraibman: gamma_2(M) * min(m,n) >= ||M||_*)

    - Upper bound from SVD: M = U diag(s) V^T = (U diag(sqrt(s))) (V diag(sqrt(s)))^T
      gamma_2 <= max_i ||row_i(U diag(sqrt(s)))|| * max_j ||row_j(V diag(sqrt(s)))||
    """
    m, n = M.shape
    s = svd_props['singular_values']
    spectral = svd_props['spectral_norm']
    nuclear = svd_props['nuclear_norm']

    # Lower bound: nuclear / min(m,n)
    lower = nuclear / min(m, n)

    # Upper bound from SVD factorization
    U_full, s_full, Vt_full = la.svd(M, full_matrices=False)
    sqrt_s = np.sqrt(s_full)
    A = U_full * sqrt_s[np.newaxis, :]   # m x k, row i = U[i,:] * sqrt(s)
    B = Vt_full.T * sqrt_s[np.newaxis, :] # n x k, row j = V[j,:] * sqrt(s)
    max_row_A = np.max(np.sqrt(np.sum(A**2, axis=1)))
    max_row_B = np.max(np.sqrt(np.sum(B**2, axis=1)))
    upper = max_row_A * max_row_B

    return lower, upper


def forster_sign_rank_bound(M):
    """
    Forster's bound: For an m x n ±1 matrix M with sign-rank r,
    ||M||_2 >= sqrt(m*n) / sqrt(r)

    So: r >= m*n / ||M||_2^2
    """
    m, n = M.shape
    spectral = la.svd(M, compute_uv=False)[0]
    lower_bound = (m * n) / (spectral ** 2)
    return max(1, int(np.ceil(lower_bound)))


def search_sign_rank(M, max_rank_to_try=None):
    """
    For small matrices, try to find the sign-rank by searching for
    rank-r matrices with the same sign pattern.

    sign-rank(M) = min r such that there exists rank-r matrix S with sign(S) = sign(M).

    Use nuclear norm minimization as a heuristic for rank minimization.
    """
    if not HAS_CVXPY:
        return None, None

    m, n = M.shape
    rank_M = int(np.linalg.matrix_rank(M))

    if max_rank_to_try is None:
        max_rank_to_try = rank_M

    # Binary search: for each candidate rank r, check feasibility
    # Use nuclear norm minimization with sign constraints

    # First, try nuclear norm minimization with sign constraints
    S = cp.Variable((m, n))
    constraints = []
    for i in range(m):
        for j in range(n):
            if M[i, j] > 0:
                constraints.append(S[i, j] >= 1.0)
            else:
                constraints.append(S[i, j] <= -1.0)

    prob = cp.Problem(cp.Minimize(cp.normNuc(S)), constraints)
    try:
        prob.solve(solver=cp.SCS, verbose=False, max_iters=20000)
        if prob.status in ['optimal', 'optimal_inaccurate']:
            S_val = S.value
            # Check the rank of the solution
            sv = la.svd(S_val, compute_uv=False)
            effective_rank = np.sum(sv > 1e-4 * sv[0])
            nuclear_val = np.sum(sv)
            return int(effective_rank), nuclear_val
    except Exception as ex:
        print(f"  Sign-rank search error: {ex}")

    return None, None


def generate_random_matrix_same_bias(m, n, neg_fraction, rng):
    """Generate a random ±1 matrix with approximately the same fraction of -1s."""
    M = np.ones((m, n))
    total = m * n
    num_neg = int(round(neg_fraction * total))
    indices = rng.choice(total, size=num_neg, replace=False)
    M.flat[indices] = -1.0
    return M


def analyze_matrix(M, label, compute_sign_rank=True):
    """Full analysis of a ±1 matrix."""
    m, n = M.shape
    results = {}

    # Basic properties
    num_neg = np.sum(M < 0)
    total = m * n
    results['size'] = (m, n)
    results['neg_fraction'] = num_neg / total
    results['num_primes'] = int(num_neg)

    # SVD properties
    svd = compute_svd_properties(M)
    results.update(svd)

    # F2 rank
    results['rank_F2'] = matrix_rank_F2(M)

    # Gamma-2 norm via SDP
    t0 = time.time()
    g2 = gamma2_norm_sdp(M)
    t_g2 = time.time() - t0
    results['gamma2'] = g2
    results['gamma2_time'] = t_g2

    # Gamma-2 bounds from SVD
    g2_lower, g2_upper = gamma2_bounds_from_svd(M, svd)
    results['gamma2_lower_spectral'] = g2_lower
    results['gamma2_upper_svd'] = g2_upper

    # Forster sign-rank lower bound
    results['sign_rank_forster_lb'] = forster_sign_rank_bound(M)

    # Sign-rank search (only for small matrices)
    if compute_sign_rank and m * n <= 1024:  # up to 32x32
        sr, nuc_min = search_sign_rank(M)
        results['sign_rank_ub_nuc'] = sr
        results['sign_rank_nuc_val'] = nuc_min
    else:
        results['sign_rank_ub_nuc'] = None
        results['sign_rank_nuc_val'] = None

    return results


def print_results(label, results):
    """Pretty-print analysis results."""
    m, n = results['size']
    print(f"\n{'='*70}")
    print(f"  {label}  (size {m}x{n})")
    print(f"{'='*70}")
    print(f"  Entries: {m*n}, Negatives (primes): {results['num_primes']} ({results['neg_fraction']:.3f})")
    print(f"  Rank over R:  {results['rank_R']}")
    print(f"  Rank over F2: {results['rank_F2']}")
    print()

    sv = results['singular_values']
    print(f"  Singular values ({len(sv)} total):")
    for i, s in enumerate(sv[:min(10, len(sv))]):
        print(f"    sigma_{i+1} = {s:.6f}")
    if len(sv) > 10:
        print(f"    ... ({len(sv)-10} more)")
    print()

    print(f"  Nuclear norm   ||M||_* = {results['nuclear_norm']:.6f}")
    print(f"  Spectral norm  ||M||_2 = {results['spectral_norm']:.6f}")
    print(f"  Frobenius norm ||M||_F = {results['frobenius_norm']:.6f}")
    print(f"  Ratio nuclear/spectral = {results['nuclear_norm']/results['spectral_norm']:.4f}")
    print(f"  Ratio nuclear/frobenius = {results['nuclear_norm']/results['frobenius_norm']:.4f}")
    print()

    if results['gamma2'] is not None:
        print(f"  GAMMA-2 NORM (SDP):     {results['gamma2']:.6f}")
    print(f"  Gamma-2 lower (nuc/dim):  {results['gamma2_lower_spectral']:.6f}")
    print(f"  Gamma-2 upper (SVD):      {results['gamma2_upper_svd']:.6f}")
    if results['gamma2'] is not None:
        print(f"  Gamma-2 / spectral norm:  {results['gamma2']/results['spectral_norm']:.6f}")
        print(f"  log2(gamma2):             {np.log2(results['gamma2']):.4f}")
    print()

    print(f"  Forster sign-rank LB: {results['sign_rank_forster_lb']}")
    print(f"  Actual rank:          {results['rank_R']}")
    if results['sign_rank_ub_nuc'] is not None:
        print(f"  Sign-rank UB (nuc-min): {results['sign_rank_ub_nuc']}")
        print(f"  Nuclear norm of sign-rank minimizer: {results['sign_rank_nuc_val']:.4f}")
    print()


def main():
    print("=" * 70)
    print("  GAMMA-2 NORM & SIGN-RANK OF PRIME INDICATOR COMMUNICATION MATRIX")
    print("=" * 70)

    N_values = [4, 6, 8, 10, 12]
    prime_results = {}
    random_results = {}

    rng = np.random.default_rng(42)
    NUM_RANDOM = 10

    for N in N_values:
        half = N // 2
        dim = 2 ** half
        print(f"\n\n{'#'*70}")
        print(f"  N = {N} bits, matrix size = {dim} x {dim}")
        print(f"{'#'*70}")

        # Build prime matrix
        t0 = time.time()
        M, pc = build_prime_matrix(N)
        t_build = time.time() - t0
        print(f"  Built prime matrix in {t_build:.3f}s, {pc} primes in range [0, {2**N})")

        # Analyze prime matrix
        do_sign_rank = (dim <= 32)
        res = analyze_matrix(M, f"Prime indicator N={N}", compute_sign_rank=do_sign_rank)
        prime_results[N] = res
        print_results(f"PRIME INDICATOR chi_P, N={N}", res)

        # Random comparison
        neg_frac = res['neg_fraction']
        print(f"\n  --- Random ±1 matrices (same size {dim}x{dim}, bias={neg_frac:.3f}) ---")
        rand_g2_list = []
        rand_spectral_list = []
        rand_nuclear_list = []
        rand_rank_list = []
        rand_sign_rank_lb_list = []
        rand_frobenius_list = []
        rand_rank_f2_list = []

        for trial in range(NUM_RANDOM):
            Mr = generate_random_matrix_same_bias(dim, dim, neg_frac, rng)
            rr = analyze_matrix(Mr, f"Random #{trial+1}", compute_sign_rank=False)
            rand_g2_list.append(rr['gamma2'])
            rand_spectral_list.append(rr['spectral_norm'])
            rand_nuclear_list.append(rr['nuclear_norm'])
            rand_rank_list.append(rr['rank_R'])
            rand_sign_rank_lb_list.append(rr['sign_rank_forster_lb'])
            rand_frobenius_list.append(rr['frobenius_norm'])
            rand_rank_f2_list.append(rr['rank_F2'])

        random_results[N] = {
            'gamma2': rand_g2_list,
            'spectral': rand_spectral_list,
            'nuclear': rand_nuclear_list,
            'rank_R': rand_rank_list,
            'rank_F2': rand_rank_f2_list,
            'sign_rank_lb': rand_sign_rank_lb_list,
            'frobenius': rand_frobenius_list,
        }

        # Print comparison
        print(f"\n  {'Measure':<25} {'Prime':>12} {'Rand mean':>12} {'Rand std':>10} {'Ratio P/R':>10}")
        print(f"  {'-'*69}")

        def compare(name, pval, rvals):
            rvals_clean = [v for v in rvals if v is not None]
            if not rvals_clean:
                print(f"  {name:<25} {pval:>12.4f} {'N/A':>12} {'N/A':>10} {'N/A':>10}")
                return
            rmean = np.mean(rvals_clean)
            rstd = np.std(rvals_clean)
            ratio = pval / rmean if rmean > 0 else float('inf')
            print(f"  {name:<25} {pval:>12.4f} {rmean:>12.4f} {rstd:>10.4f} {ratio:>10.4f}")

        compare("Rank (R)", res['rank_R'], rand_rank_list)
        compare("Rank (F2)", res['rank_F2'], rand_rank_f2_list)
        compare("Spectral norm", res['spectral_norm'], rand_spectral_list)
        compare("Nuclear norm", res['nuclear_norm'], rand_nuclear_list)
        compare("Frobenius norm", res['frobenius_norm'], rand_frobenius_list)
        if res['gamma2'] is not None:
            compare("GAMMA-2 NORM", res['gamma2'], rand_g2_list)
        compare("Forster sign-rank LB", res['sign_rank_forster_lb'], rand_sign_rank_lb_list)

    # Summary table
    print(f"\n\n{'#'*70}")
    print(f"  SUMMARY: SCALING OF MEASURES WITH N")
    print(f"{'#'*70}")
    print(f"\n  {'N':>3} {'dim':>5} {'rank_R':>7} {'rank_F2':>8} {'spectral':>10} {'nuclear':>10} "
          f"{'gamma2':>10} {'log2(g2)':>10} {'Forster':>8} {'sign_rk_ub':>11}")
    print(f"  {'-'*93}")

    for N in N_values:
        r = prime_results[N]
        dim = 2 ** (N // 2)
        g2 = r['gamma2']
        g2_str = f"{g2:.4f}" if g2 is not None else "N/A"
        lg2_str = f"{np.log2(g2):.4f}" if g2 is not None else "N/A"
        sr_ub = r.get('sign_rank_ub_nuc')
        sr_str = str(sr_ub) if sr_ub is not None else "N/A"
        print(f"  {N:>3} {dim:>5} {r['rank_R']:>7} {r['rank_F2']:>8} "
              f"{r['spectral_norm']:>10.4f} {r['nuclear_norm']:>10.4f} "
              f"{g2_str:>10} {lg2_str:>10} {r['sign_rank_forster_lb']:>8} {sr_str:>11}")

    # Scaling analysis
    print(f"\n\n  SCALING ANALYSIS:")
    g2_vals = [(N, prime_results[N]['gamma2']) for N in N_values if prime_results[N]['gamma2'] is not None]
    if len(g2_vals) >= 2:
        for i in range(1, len(g2_vals)):
            N_prev, g2_prev = g2_vals[i-1]
            N_curr, g2_curr = g2_vals[i]
            if g2_prev > 0:
                ratio = g2_curr / g2_prev
                dim_ratio = 2 ** ((N_curr - N_prev) // 2)
                print(f"  N={N_prev}->{N_curr}: gamma2 ratio = {ratio:.4f}, "
                      f"dim ratio = {dim_ratio}, "
                      f"gamma2 growth per dim doubling = {ratio/dim_ratio:.4f}")

    # Compare gamma2/dim scaling
    print(f"\n  gamma2 / dim ratios:")
    for N in N_values:
        g2 = prime_results[N]['gamma2']
        dim = 2 ** (N // 2)
        if g2 is not None:
            print(f"    N={N:>2}: gamma2/dim = {g2/dim:.6f}, gamma2/sqrt(dim) = {g2/np.sqrt(dim):.6f}, "
                  f"gamma2/log(dim) = {g2/np.log2(dim):.6f}")

    # Prime vs Random comparison summary
    print(f"\n\n  PRIME vs RANDOM GAMMA-2 COMPARISON:")
    print(f"  {'N':>3} {'Prime g2':>10} {'Rand g2 mean':>13} {'Rand g2 std':>12} {'P/R ratio':>10}")
    print(f"  {'-'*50}")
    for N in N_values:
        g2_p = prime_results[N]['gamma2']
        g2_r = random_results[N]['gamma2']
        g2_r_clean = [v for v in g2_r if v is not None]
        if g2_p is not None and g2_r_clean:
            rmean = np.mean(g2_r_clean)
            rstd = np.std(g2_r_clean)
            print(f"  {N:>3} {g2_p:>10.4f} {rmean:>13.4f} {rstd:>12.4f} {g2_p/rmean:>10.4f}")

    print(f"\n\n  KEY QUESTIONS:")
    print(f"  1. Does gamma_2 grow polynomially or exponentially in N?")
    g2_vals_only = [prime_results[N]['gamma2'] for N in N_values if prime_results[N]['gamma2'] is not None]
    if len(g2_vals_only) >= 3:
        # Fit log(gamma2) vs N
        Ns = [N for N in N_values if prime_results[N]['gamma2'] is not None]
        log_g2 = [np.log2(g) for g in g2_vals_only]
        coeffs = np.polyfit(Ns, log_g2, 1)
        print(f"     log2(gamma2) ~ {coeffs[0]:.4f} * N + {coeffs[1]:.4f}")
        print(f"     => gamma2 ~ 2^({coeffs[0]:.4f} * N)")
        if coeffs[0] > 0.4:
            print(f"     EXPONENTIAL growth (base ~{2**coeffs[0]:.4f} per bit)")
        else:
            print(f"     Sub-exponential / polynomial growth")

    for N in N_values:
        r = prime_results[N]
        if r['sign_rank_ub_nuc'] is not None and r['sign_rank_ub_nuc'] < r['rank_R']:
            print(f"  2. N={N}: sign-rank ({r['sign_rank_ub_nuc']}) < rank ({r['rank_R']}) => HIDDEN LOW-RANK STRUCTURE")
        elif r['sign_rank_ub_nuc'] is not None:
            print(f"  2. N={N}: sign-rank UB ({r['sign_rank_ub_nuc']}) >= rank ({r['rank_R']}) => no hidden structure detected")

    print(f"\n  Done.")


if __name__ == '__main__':
    main()

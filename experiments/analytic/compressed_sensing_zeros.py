"""
Session 7: Compressed Sensing and Sparse Recovery for Zeta Zero Sums

The explicit formula: pi(x) = R(x) - sum_rho R(x^rho) + small terms

Previous sessions showed:
- Need ~sqrt(x) zeros for exact results (infeasible for large x)
- Raw summation converges as K^{-0.01} (barely converges)
- Smoothing kernels help somewhat but don't change asymptotics

THIS SESSION explores whether SPARSE/ADAPTIVE zero selection can do better:

1. For a given x, which zeros contribute most? (cos(gamma*ln(x))/gamma)
2. Can we select O(polylog(x)) zeros that capture the correction to O(1)?
3. Matching pursuit / greedy zero selection
4. Random sampling with importance weights
5. Resonance analysis (Gram-like points)
6. Prony/ESPRIT signal recovery from sparse samples of psi(x)

Key insight to test: The zero sum is NOT a generic trigonometric sum.
The zeros have GUE statistics (Montgomery-Odlyzko). Does this structure
allow compression?
"""

import numpy as np
from mpmath import mp, mpf, li, log, sqrt, pi as mpi, zetazero, im, re, exp, cos, sin, ei
import time
import json
from sympy import primepi, prime as sympy_prime, isprime, nextprime

mp.dps = 40

# =============================================================================
# INFRASTRUCTURE
# =============================================================================

def get_zeta_zeros(K):
    """Compute first K zeta zero imaginary parts."""
    print(f"  Computing {K} zeta zeros...")
    t0 = time.time()
    zeros = []
    for k in range(1, K + 1):
        z = zetazero(k)
        zeros.append(float(im(z)))
    print(f"  Done in {time.time()-t0:.1f}s")
    return zeros

def R_function(x):
    """Riemann's R function."""
    x = mpf(x)
    mobius = {1:1, 2:-1, 3:-1, 5:-1, 6:1, 7:-1, 10:1, 11:-1, 13:-1,
              14:1, 15:1, 17:-1, 19:-1, 21:1, 22:1, 23:-1, 26:1, 29:-1, 30:-1}
    result = mpf(0)
    for k in range(1, 50):
        from sympy import mobius as sympy_mobius
        mu_k = mobius.get(k, int(sympy_mobius(k)))
        if mu_k == 0:
            continue
        xk = x ** (mpf(1)/k)
        if xk < 1.01:
            break
        result += mpf(mu_k) / k * li(xk)
    return result

def zero_contribution(gamma, x):
    """
    Compute the contribution of zero rho = 1/2 + i*gamma to the explicit formula.

    Contribution = -2 * Re(li(x^rho)) where li(z) = Ei(rho * ln(x))

    For large x, this is approximately:
    -2 * sqrt(x) * cos(gamma*ln(x) - arg(rho)) / (|rho| * ln(x))
    """
    x_mp = mpf(x)
    lnx = log(x_mp)
    rho = mp.mpc(0.5, gamma)
    # Exact: li(x^rho) = Ei(rho * ln(x))
    val = ei(rho * lnx)
    return -2 * float(re(val))

def zero_contribution_approx(gamma, x):
    """
    Fast approximate contribution of zero at gamma for value x.

    li(x^rho) ~ x^rho / (rho * ln(x))
    Re part ~ sqrt(x) * [0.5*cos(t) + gamma*sin(t)] / ((0.25+gamma^2)*ln(x))
    where t = gamma * ln(x)
    """
    lnx = np.log(x)
    sqrtx = np.sqrt(x)
    t = gamma * lnx
    denom = (0.25 + gamma**2) * lnx
    real_part = sqrtx * (0.5 * np.cos(t) + gamma * np.sin(t)) / denom
    return -2 * real_part

def pi_with_zeros(x, gamma_list, weights=None):
    """
    Compute pi(x) using R(x) + weighted sum of zero contributions.

    If weights is None, all weights are 1.
    """
    result = float(R_function(mpf(x)))
    for i, gamma in enumerate(gamma_list):
        w = 1.0 if weights is None else weights[i]
        result += w * zero_contribution_approx(gamma, x)
    return result

def pi_with_zeros_exact(x, gamma_list, weights=None):
    """Same but with exact zero contributions (slower)."""
    result = float(R_function(mpf(x)))
    for i, gamma in enumerate(gamma_list):
        w = 1.0 if weights is None else weights[i]
        result += w * zero_contribution(gamma, x)
    return result


# =============================================================================
# EXPERIMENT 1: SPARSE ZERO SELECTION - Which zeros matter most for given x?
# =============================================================================

def experiment1_sparse_selection(zeros, test_cases):
    """
    For each test x, rank zeros by |contribution| and select top K.
    Compare: top-K selection vs first-K (sequential) vs random-K.
    """
    print("\n" + "="*75)
    print("EXPERIMENT 1: SPARSE ZERO SELECTION (contribution-ranked)")
    print("="*75)

    results = {}

    for name, x, true_pi in test_cases:
        print(f"\n--- x = {x} (pi(x) = {true_pi}) ---")

        # Compute all contributions
        contribs = []
        for gamma in zeros:
            c = zero_contribution_approx(gamma, x)
            contribs.append(c)
        contribs = np.array(contribs)

        # Base value (no zeros)
        R_x = float(R_function(mpf(x)))
        base_error = abs(R_x - true_pi)
        print(f"  R(x) = {R_x:.4f}, base error = {base_error:.4f}")

        # Rank by absolute contribution
        ranked_indices = np.argsort(-np.abs(contribs))

        # Show top contributors
        print(f"  Top 10 contributors (gamma, |contribution|):")
        for i in range(min(10, len(ranked_indices))):
            idx = ranked_indices[i]
            print(f"    #{i+1}: gamma_{idx+1} = {zeros[idx]:.4f}, "
                  f"contrib = {contribs[idx]:+.6f}, "
                  f"|contrib| = {abs(contribs[idx]):.6f}")

        # Compare selection strategies for various K
        K_values = [1, 2, 3, 5, 8, 10, 15, 20, 30, 50, 75, 100]
        K_values = [k for k in K_values if k <= len(zeros)]

        print(f"\n  {'K':>4} | {'First-K err':>12} | {'Top-K err':>12} | {'Random-K err':>12} | {'Top-K indices':>20}")
        print("  " + "-"*75)

        case_results = []
        for K in K_values:
            # Strategy 1: First K zeros (sequential)
            first_k_sum = R_x + np.sum(contribs[:K])
            first_k_err = abs(first_k_sum - true_pi)

            # Strategy 2: Top K by absolute contribution
            top_k_idx = ranked_indices[:K]
            top_k_sum = R_x + np.sum(contribs[top_k_idx])
            top_k_err = abs(top_k_sum - true_pi)

            # Strategy 3: Random K (average over 20 trials)
            rng = np.random.RandomState(42)
            rand_errs = []
            for _ in range(20):
                rand_idx = rng.choice(len(zeros), K, replace=False)
                rand_sum = R_x + np.sum(contribs[rand_idx])
                rand_errs.append(abs(rand_sum - true_pi))
            rand_k_err = np.mean(rand_errs)

            top_indices_str = str(sorted(top_k_idx[:5].tolist()))
            print(f"  {K:4d} | {first_k_err:12.4f} | {top_k_err:12.4f} | "
                  f"{rand_k_err:12.4f} | {top_indices_str}")

            case_results.append({
                'K': K,
                'first_k_err': first_k_err,
                'top_k_err': top_k_err,
                'random_k_err': rand_k_err
            })

        results[name] = case_results

    return results


# =============================================================================
# EXPERIMENT 2: MATCHING PURSUIT (greedy zero selection)
# =============================================================================

def experiment2_matching_pursuit(zeros, test_cases):
    """
    Greedy algorithm: iteratively select the zero that reduces the residual most.
    This is matching pursuit adapted to the zero sum problem.
    """
    print("\n" + "="*75)
    print("EXPERIMENT 2: MATCHING PURSUIT (greedy residual reduction)")
    print("="*75)

    results = {}

    for name, x, true_pi in test_cases:
        print(f"\n--- x = {x} (pi(x) = {true_pi}) ---")

        R_x = float(R_function(mpf(x)))

        # Compute contribution of each zero
        contribs = np.array([zero_contribution_approx(g, x) for g in zeros])

        # Target: we want sum of selected contributions to equal (true_pi - R_x)
        target = true_pi - R_x
        residual = target
        selected = []
        selected_set = set()

        print(f"  Target correction = {target:.4f}")
        print(f"  {'Step':>4} | {'Selected gamma':>14} | {'Contrib':>10} | {'Residual':>12} | {'|Error|':>10}")
        print("  " + "-"*65)

        K_max = min(50, len(zeros))
        case_results = []

        for step in range(K_max):
            # Find zero that reduces |residual| most
            best_idx = None
            best_new_resid = abs(residual)

            for i in range(len(zeros)):
                if i in selected_set:
                    continue
                new_resid = abs(residual - contribs[i])
                if new_resid < best_new_resid:
                    best_new_resid = new_resid
                    best_idx = i

            if best_idx is None:
                break

            selected.append(best_idx)
            selected_set.add(best_idx)
            residual -= contribs[best_idx]
            error = abs(residual)

            if step < 20 or (step + 1) % 10 == 0:
                print(f"  {step+1:4d} | gamma_{best_idx+1:>4d} = {zeros[best_idx]:8.3f} | "
                      f"{contribs[best_idx]:+10.4f} | {residual:+12.4f} | {error:10.4f}")

            case_results.append({
                'step': step + 1,
                'zero_index': best_idx + 1,
                'gamma': zeros[best_idx],
                'error': error
            })

            # Early stop
            if error < 0.5:
                print(f"  *** Error < 0.5 reached at step {step+1}! ***")
                break

        results[name] = case_results

        # Compare: how many sequential zeros needed for same accuracy?
        mp_error_at_20 = case_results[min(19, len(case_results)-1)]['error'] if len(case_results) >= 1 else float('inf')
        seq_errors = [abs(true_pi - R_x - np.sum(contribs[:k])) for k in range(1, len(zeros)+1)]
        for k, e in enumerate(seq_errors, 1):
            if e <= mp_error_at_20:
                print(f"  Sequential needs K={k} zeros for same accuracy as MP with 20")
                break
        else:
            print(f"  Sequential never reaches MP-20 accuracy ({mp_error_at_20:.4f})")

    return results


# =============================================================================
# EXPERIMENT 3: IMPORTANCE SAMPLING with GUE statistics
# =============================================================================

def experiment3_importance_sampling(zeros, test_cases):
    """
    Random sampling of zeros with importance weights.

    Key idea: contribution of zero gamma is ~1/gamma, so sample with
    probability proportional to 1/gamma (or |contribution|).
    Use importance weighting to get unbiased estimate.
    """
    print("\n" + "="*75)
    print("EXPERIMENT 3: IMPORTANCE SAMPLING OF ZEROS")
    print("="*75)

    results = {}
    N = len(zeros)

    for name, x, true_pi in test_cases:
        print(f"\n--- x = {x} (pi(x) = {true_pi}) ---")

        R_x = float(R_function(mpf(x)))
        contribs = np.array([zero_contribution_approx(g, x) for g in zeros])
        full_correction = np.sum(contribs)
        full_estimate = R_x + full_correction
        full_error = abs(full_estimate - true_pi)

        print(f"  Full sum ({N} zeros): est = {full_estimate:.4f}, err = {full_error:.4f}")

        # Importance distribution: p_k proportional to 1/gamma_k
        gammas = np.array(zeros)
        p_importance = 1.0 / gammas
        p_importance /= p_importance.sum()

        # Uniform distribution for comparison
        p_uniform = np.ones(N) / N

        # Contribution-weighted distribution
        abs_contribs = np.abs(contribs)
        p_contrib = abs_contribs / abs_contribs.sum() if abs_contribs.sum() > 0 else p_uniform

        K_values = [5, 10, 20, 50, 100]
        K_values = [k for k in K_values if k <= N]
        n_trials = 100

        print(f"\n  {'K':>4} | {'Uniform RMSE':>12} | {'1/gamma RMSE':>12} | {'|c| RMSE':>12} | {'Uniform bias':>12}")
        print("  " + "-"*65)

        case_results = []
        for K in K_values:
            rng = np.random.RandomState(123)

            # Uniform sampling (rescale by N/K)
            uniform_estimates = []
            for _ in range(n_trials):
                idx = rng.choice(N, K, replace=False)
                est = R_x + np.sum(contribs[idx]) * (N / K)
                uniform_estimates.append(est)

            # Importance sampling with 1/gamma
            importance_estimates_1g = []
            for _ in range(n_trials):
                idx = rng.choice(N, K, replace=False, p=p_importance)
                # IS correction: weight each sample by 1/(K * p_k)
                weights = 1.0 / (K * p_importance[idx])
                est = R_x + np.sum(contribs[idx] * weights)
                importance_estimates_1g.append(est)

            # Importance sampling with |contribution|
            importance_estimates_c = []
            for _ in range(n_trials):
                idx = rng.choice(N, K, replace=False, p=p_contrib)
                weights = 1.0 / (K * p_contrib[idx])
                est = R_x + np.sum(contribs[idx] * weights)
                importance_estimates_c.append(est)

            u_rmse = np.sqrt(np.mean([(e - true_pi)**2 for e in uniform_estimates]))
            i1_rmse = np.sqrt(np.mean([(e - true_pi)**2 for e in importance_estimates_1g]))
            ic_rmse = np.sqrt(np.mean([(e - true_pi)**2 for e in importance_estimates_c]))
            u_bias = np.mean(uniform_estimates) - true_pi

            print(f"  {K:4d} | {u_rmse:12.4f} | {i1_rmse:12.4f} | {ic_rmse:12.4f} | {u_bias:+12.4f}")

            case_results.append({
                'K': K,
                'uniform_rmse': u_rmse,
                'importance_1g_rmse': i1_rmse,
                'importance_contrib_rmse': ic_rmse,
                'uniform_bias': u_bias
            })

        results[name] = case_results

    return results


# =============================================================================
# EXPERIMENT 4: OPTIMAL K-ZERO SELECTION (optimization)
# =============================================================================

def experiment4_optimal_selection(zeros, test_cases):
    """
    Given budget K, find the K zeros that minimize |pi(x) - estimate|.

    This is a combinatorial optimization problem. We use:
    - Brute force for small K
    - Greedy (matching pursuit from exp 2) as heuristic
    - Contribution-based ranking (exp 1)
    """
    print("\n" + "="*75)
    print("EXPERIMENT 4: OPTIMAL K-ZERO SUBSET (combinatorial)")
    print("="*75)

    from itertools import combinations

    results = {}

    for name, x, true_pi in test_cases:
        print(f"\n--- x = {x} (pi(x) = {true_pi}) ---")

        R_x = float(R_function(mpf(x)))
        contribs = np.array([zero_contribution_approx(g, x) for g in zeros])
        target = true_pi - R_x

        N = min(30, len(zeros))  # Limit for combinatorial search
        contribs_small = contribs[:N]

        case_results = []

        for K in [1, 2, 3, 4, 5]:
            if K > N:
                break

            best_err = float('inf')
            best_subset = None
            n_combos = 0

            # Brute force for small K and N
            if N <= 30 and K <= 5:
                for combo in combinations(range(N), K):
                    s = sum(contribs_small[i] for i in combo)
                    err = abs(target - s)
                    if err < best_err:
                        best_err = err
                        best_subset = combo
                    n_combos += 1

            # Also compute sequential and top-K for comparison
            seq_err = abs(target - np.sum(contribs[:K]))
            ranked = np.argsort(-np.abs(contribs[:N]))
            topk_err = abs(target - np.sum(contribs[ranked[:K]]))

            gamma_strs = [f"g{i+1}={zeros[i]:.1f}" for i in best_subset] if best_subset else ["?"]
            print(f"  K={K}: optimal_err={best_err:.6f} ({','.join(gamma_strs[:3])}), "
                  f"sequential_err={seq_err:.6f}, topK_err={topk_err:.6f} "
                  f"[searched {n_combos} combos]")

            case_results.append({
                'K': K,
                'optimal_error': best_err,
                'sequential_error': seq_err,
                'topk_error': topk_err,
                'optimal_subset': list(best_subset) if best_subset else [],
                'combos_searched': n_combos
            })

        results[name] = case_results

    return results


# =============================================================================
# EXPERIMENT 5: RESONANCE PHENOMENA (Gram-like points)
# =============================================================================

def experiment5_resonances(zeros, test_cases_extended):
    """
    Are there special values of x where the zero sum simplifies?

    Gram points: theta(g_n) = n*pi where theta is Riemann-Siegel theta.
    At these points, Z(g_n) = (-1)^n zeta(1/2 + i*g_n), which is real.

    For the zero sum, examine constructive/destructive interference.
    The contribution is cos(gamma_k * ln(x)). If ln(x) = 2*pi*m/gamma_1
    for integer m, then zero 1 contributes maximally.

    But for ALL zeros to interfere constructively, we'd need gamma_k * ln(x)
    near multiples of 2*pi for all k simultaneously.
    The gamma_k are incommensurate, so this is impossible (equidistribution).

    BUT: Maybe for some x, MOST zeros nearly cancel, leaving a sparse sum?
    """
    print("\n" + "="*75)
    print("EXPERIMENT 5: RESONANCE ANALYSIS")
    print("="*75)

    results = {}

    # For a range of x values, compute the "effective number of zeros"
    # = (sum of contribs)^2 / sum of contribs^2 (participation ratio)

    gammas = np.array(zeros[:100])

    x_range = np.logspace(2, 6, 200)

    print("\n  Participation ratio: how many zeros effectively contribute")
    print(f"  (1 = single zero dominates, N = all contribute equally)\n")

    participation_ratios = []
    total_corrections = []

    for x in x_range:
        lnx = np.log(x)
        sqrtx = np.sqrt(x)

        # Approximate contributions (magnitude only)
        amplitudes = []
        for g in gammas:
            t = g * lnx
            denom = (0.25 + g**2) * lnx
            amp = sqrtx * np.sqrt(0.25 * np.cos(t)**2 + g**2 * np.sin(t)**2 +
                                   g * np.cos(t) * np.sin(t)) / denom
            amplitudes.append(amp)

        # Actually use the signed contributions
        contribs = np.array([zero_contribution_approx(g, x) for g in gammas])
        abs_contribs = np.abs(contribs)

        if abs_contribs.sum() > 0:
            # Participation ratio
            p2 = np.sum(abs_contribs)**2
            p4 = np.sum(abs_contribs**2)
            pr = p2 / p4 if p4 > 0 else 0
            participation_ratios.append(pr)
            total_corrections.append(np.sum(contribs))
        else:
            participation_ratios.append(0)
            total_corrections.append(0)

    pr_arr = np.array(participation_ratios)
    tc_arr = np.array(total_corrections)

    # Find x values where participation ratio is smallest (most sparse)
    sorted_idx = np.argsort(pr_arr)

    print(f"  Participation ratio statistics (100 zeros):")
    print(f"    Mean: {np.mean(pr_arr):.2f}")
    print(f"    Std:  {np.std(pr_arr):.2f}")
    print(f"    Min:  {np.min(pr_arr):.2f} at x = {x_range[sorted_idx[0]]:.0f}")
    print(f"    Max:  {np.max(pr_arr):.2f} at x = {x_range[sorted_idx[-1]]:.0f}")

    # Show most "sparse" points
    print(f"\n  Most sparse x values (lowest participation ratio):")
    for i in range(min(10, len(sorted_idx))):
        idx = sorted_idx[i]
        x = x_range[idx]
        actual_pi = int(primepi(int(x)))
        R_x = float(R_function(mpf(int(x))))
        corr = tc_arr[idx]
        est = R_x + corr
        err = abs(est - actual_pi)
        print(f"    x={x:>10.0f}: PR={pr_arr[idx]:.2f}, "
              f"correction={corr:+.2f}, est={est:.2f}, pi(x)={actual_pi}, err={err:.2f}")

    # Phase analysis: at what fraction of zeros is cos(gamma*lnx) near 0?
    print(f"\n  Phase cancellation analysis:")
    for x in [100, 1000, 10000, 100000]:
        lnx = np.log(x)
        phases = (gammas * lnx) % (2 * np.pi)
        # Count near 0 or pi (cos near +/-1) vs near pi/2 (cos near 0)
        cos_vals = np.cos(phases)
        near_zero = np.sum(np.abs(cos_vals) < 0.1)
        near_one = np.sum(np.abs(cos_vals) > 0.9)
        print(f"    x={x:>6d}: |cos| < 0.1: {near_zero}/100, |cos| > 0.9: {near_one}/100, "
              f"mean |cos|: {np.mean(np.abs(cos_vals)):.4f}")

    results['participation_ratio_stats'] = {
        'mean': float(np.mean(pr_arr)),
        'std': float(np.std(pr_arr)),
        'min': float(np.min(pr_arr)),
        'max': float(np.max(pr_arr))
    }

    return results


# =============================================================================
# EXPERIMENT 6: PRONY / ESPRIT SIGNAL RECOVERY
# =============================================================================

def experiment6_prony_esprit(test_cases):
    """
    psi(x) - x = -sum_rho x^rho/rho is a sum of complex exponentials.

    On a log scale: let u = ln(x), then
    psi(e^u) - e^u = -sum_k e^{(1/2 + i*gamma_k)*u} / (1/2 + i*gamma_k)
                    = -e^{u/2} * sum_k e^{i*gamma_k*u} / (1/2 + i*gamma_k)

    After dividing by e^{u/2}, the signal is:
    S(u) = sum_k c_k * e^{i*gamma_k*u}

    This is EXACTLY the form amenable to Prony's method / ESPRIT / MUSIC.

    Can we recover gamma_k from samples of S(u)?
    """
    print("\n" + "="*75)
    print("EXPERIMENT 6: PRONY/ESPRIT SIGNAL RECOVERY FROM psi(x)")
    print("="*75)

    # Compute psi(x) exactly for small x using sympy
    from sympy import factorint

    def chebyshev_psi(x):
        """Exact Chebyshev psi function."""
        x = int(x)
        result = 0.0
        for p in range(2, x + 1):
            if not isprime(p):
                continue
            pk = p
            while pk <= x:
                result += np.log(p)
                pk *= p
        return result

    # Sample S(u) at uniform points
    print("\n  Computing Chebyshev psi(x) samples...")

    # Use u = ln(x), sample at uniform u spacing
    u_min, u_max = np.log(100), np.log(10000)
    N_samples = 200
    u_vals = np.linspace(u_min, u_max, N_samples)

    S_vals = []
    for u in u_vals:
        x = np.exp(u)
        psi_x = chebyshev_psi(int(x))
        S_u = (psi_x - x) / np.sqrt(x)  # Remove sqrt(x) factor
        S_vals.append(S_u)

    S_vals = np.array(S_vals)

    print(f"  Computed {N_samples} samples of S(u) = (psi(e^u) - e^u)/e^(u/2)")
    print(f"  u range: [{u_min:.2f}, {u_max:.2f}]")
    print(f"  S range: [{S_vals.min():.4f}, {S_vals.max():.4f}]")

    # Prony's method: build Hankel matrix and extract frequencies
    print("\n  --- Matrix Pencil Method ---")

    M = N_samples // 2  # Pencil parameter

    # Build Hankel-like matrices Y0 and Y1
    L = N_samples - M
    Y0 = np.zeros((L, M), dtype=complex)
    Y1 = np.zeros((L, M), dtype=complex)

    for i in range(L):
        for j in range(M):
            Y0[i, j] = S_vals[i + j]
            Y1[i, j] = S_vals[i + j + 1] if i + j + 1 < N_samples else 0

    # SVD of Y0
    U, sigma, Vh = np.linalg.svd(Y0, full_matrices=False)

    # Singular value spectrum - how many significant components?
    print(f"\n  Singular value spectrum (first 20):")
    for i in range(min(20, len(sigma))):
        print(f"    sigma_{i+1:>3d} = {sigma[i]:.6f}")

    # Effective rank
    thresh = sigma[0] * 1e-3
    rank = np.sum(sigma > thresh)
    print(f"\n  Effective rank (threshold {thresh:.4f}): {rank}")

    # Try to recover frequencies using truncated SVD + generalized eigenvalues
    for K in [5, 10, 20]:
        if K > rank:
            K = rank

        # Truncated matrices
        U_k = U[:, :K]
        S_k = np.diag(sigma[:K])
        Vh_k = Vh[:K, :]

        # Matrix pencil: find eigenvalues of (U_k^H Y1 Vh_k^H S_k^{-1})
        try:
            S_k_inv = np.diag(1.0 / sigma[:K])
            Z = U_k.conj().T @ Y1 @ Vh_k.conj().T @ S_k_inv
            eigvals = np.linalg.eigvals(Z)

            # Convert eigenvalues to frequencies
            # z_k = exp(i * gamma_k * du) where du is the u-spacing
            du = u_vals[1] - u_vals[0]
            freqs = np.angle(eigvals) / du

            # Sort by magnitude (most significant)
            freq_mags = np.abs(freqs)
            sorted_idx = np.argsort(freq_mags)
            freqs_sorted = freqs[sorted_idx]

            # Compare with known zeta zeros
            known_zeros = [14.134725, 21.022040, 25.010858, 30.424876, 32.935062,
                          37.586178, 40.918719, 43.327073, 48.005151, 49.773832]

            print(f"\n  K={K} components: recovered frequencies vs known zeta zeros:")
            matched = 0
            for kz in known_zeros[:min(K, 10)]:
                dists = np.abs(np.abs(freqs) - kz)
                best_match = np.min(dists)
                best_freq = np.abs(freqs)[np.argmin(dists)]
                match_str = "MATCH" if best_match < 1.0 else "miss"
                if best_match < 1.0:
                    matched += 1
                if kz < 55:
                    print(f"    gamma = {kz:.3f} -> nearest recovered = {best_freq:.3f} "
                          f"(diff = {best_match:.3f}) {match_str}")
            print(f"    Matched: {matched}/{min(K, len(known_zeros))}")

        except Exception as e:
            print(f"  K={K}: failed - {e}")

    return {'n_samples': N_samples, 'effective_rank': int(rank)}


# =============================================================================
# EXPERIMENT 7: ADAPTIVE ZERO SELECTION WITH CROSS-VALIDATION
# =============================================================================

def experiment7_adaptive_crossval(zeros, test_cases):
    """
    Select zeros adaptively using multiple x values.

    Key question: if we find good zeros for x1, do they also work for x2?
    (If yes, there's a universal sparse set; if no, the problem is x-dependent.)
    """
    print("\n" + "="*75)
    print("EXPERIMENT 7: CROSS-VALIDATION OF ZERO SELECTION")
    print("="*75)

    # Use matching pursuit on training x values, test on held-out values
    all_x = [100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000]
    all_pi = [int(primepi(x)) for x in all_x]

    train_idx = [0, 2, 4, 6, 8]  # 100, 500, 2000, 10000, 50000
    test_idx = [1, 3, 5, 7, 9]   # 200, 1000, 5000, 20000, 100000

    print(f"  Train: {[all_x[i] for i in train_idx]}")
    print(f"  Test:  {[all_x[i] for i in test_idx]}")

    # For each train x, do matching pursuit to get order of zero importance
    zero_ranks = {}  # zero_index -> list of ranks from each training x

    for ti in train_idx:
        x = all_x[ti]
        R_x = float(R_function(mpf(x)))
        contribs = np.array([zero_contribution_approx(g, x) for g in zeros])
        target = all_pi[ti] - R_x

        # Matching pursuit order
        residual = target
        selected = []
        remaining = set(range(len(zeros)))

        for step in range(min(50, len(zeros))):
            best_idx = None
            best_resid = abs(residual)
            for i in remaining:
                new_resid = abs(residual - contribs[i])
                if new_resid < best_resid:
                    best_resid = new_resid
                    best_idx = i
            if best_idx is None:
                break
            selected.append(best_idx)
            remaining.remove(best_idx)
            residual -= contribs[best_idx]

        for rank, idx in enumerate(selected):
            if idx not in zero_ranks:
                zero_ranks[idx] = []
            zero_ranks[idx].append(rank)

    # Average rank across training x values
    avg_ranks = {}
    for idx, ranks in zero_ranks.items():
        avg_ranks[idx] = np.mean(ranks)

    # Sort by average rank (lower = more universally important)
    universal_order = sorted(avg_ranks.keys(), key=lambda i: avg_ranks[i])

    print(f"\n  Most universally important zeros (by avg MP rank):")
    for i in range(min(20, len(universal_order))):
        idx = universal_order[i]
        n_appearances = len(zero_ranks[idx])
        print(f"    Zero #{idx+1} (gamma={zeros[idx]:.3f}): avg_rank={avg_ranks[idx]:.1f}, "
              f"appears in {n_appearances}/{len(train_idx)} training sets")

    # Now test: use universal order vs sequential vs per-x optimal
    print(f"\n  Cross-validation: error on TEST set")
    print(f"  {'K':>4} | {'Sequential':>12} | {'Universal MP':>12} | {'Per-x optimal':>12}")
    print("  " + "-"*50)

    for K in [5, 10, 20, 30, 50]:
        if K > len(zeros):
            break

        seq_errs = []
        univ_errs = []
        perx_errs = []

        for ti in test_idx:
            x = all_x[ti]
            R_x = float(R_function(mpf(x)))
            contribs = np.array([zero_contribution_approx(g, x) for g in zeros])
            true = all_pi[ti]

            # Sequential
            seq_est = R_x + np.sum(contribs[:K])
            seq_errs.append(abs(seq_est - true))

            # Universal MP order
            univ_idx = universal_order[:K]
            univ_est = R_x + np.sum(contribs[univ_idx])
            univ_errs.append(abs(univ_est - true))

            # Per-x top-K (oracle - knows the answer)
            ranked = np.argsort(-np.abs(contribs))
            perx_est = R_x + np.sum(contribs[ranked[:K]])
            perx_errs.append(abs(perx_est - true))

        print(f"  {K:4d} | {np.mean(seq_errs):12.4f} | {np.mean(univ_errs):12.4f} | {np.mean(perx_errs):12.4f}")

    # KEY QUESTION: Does the universal order generalize?
    print(f"\n  Sequential order: 1, 2, 3, ... (by gamma magnitude)")
    print(f"  Universal MP order: {universal_order[:10]}")
    sequential_order = list(range(len(zeros)))

    # Correlation between sequential and universal
    so_ranks = {i: i for i in range(len(zeros))}
    uo_ranks = {idx: rank for rank, idx in enumerate(universal_order)}
    common = set(so_ranks.keys()) & set(uo_ranks.keys())
    if common:
        so_vals = [so_ranks[i] for i in sorted(common)]
        uo_vals = [uo_ranks[i] for i in sorted(common)]
        corr = np.corrcoef(so_vals, uo_vals)[0, 1]
        print(f"  Rank correlation (sequential vs universal): {corr:.4f}")

    return {'universal_order': universal_order[:50]}


# =============================================================================
# EXPERIMENT 8: ERROR SCALING ANALYSIS
# =============================================================================

def experiment8_error_scaling(zeros):
    """
    Key question: How does the minimum achievable error scale with K zeros?

    If error ~ K^{-alpha}, what is alpha for each strategy?
    Is there any strategy with alpha > 1 (super-linear convergence)?
    """
    print("\n" + "="*75)
    print("EXPERIMENT 8: ERROR SCALING WITH NUMBER OF ZEROS")
    print("="*75)

    x_values = [1000, 10000, 100000]

    for x in x_values:
        true = int(primepi(x))
        R_x = float(R_function(mpf(x)))
        contribs = np.array([zero_contribution_approx(g, x) for g in zeros])
        target = true - R_x

        print(f"\n--- x = {x}, pi(x) = {true}, target correction = {target:.4f} ---")

        K_values = list(range(1, min(len(zeros), 101)))

        # Strategy errors
        seq_errors = []
        mp_errors = []
        topk_errors = []

        # Sequential
        cumsum = np.cumsum(contribs)
        for K in K_values:
            seq_errors.append(abs(target - cumsum[K-1]))

        # Top-K by contribution
        ranked = np.argsort(-np.abs(contribs))
        for K in K_values:
            s = np.sum(contribs[ranked[:K]])
            topk_errors.append(abs(target - s))

        # Matching pursuit
        residual = target
        selected_set = set()
        for K in K_values:
            best_idx = None
            best_resid = abs(residual)
            for i in range(len(zeros)):
                if i in selected_set:
                    continue
                nr = abs(residual - contribs[i])
                if nr < best_resid:
                    best_resid = nr
                    best_idx = i
            if best_idx is not None:
                selected_set.add(best_idx)
                residual -= contribs[best_idx]
            mp_errors.append(abs(residual))

        seq_errors = np.array(seq_errors)
        mp_errors = np.array(mp_errors)
        topk_errors = np.array(topk_errors)

        # Fit power law: error ~ K^{-alpha}
        K_arr = np.array(K_values, dtype=float)

        def fit_alpha(errors, label):
            # Only fit where error > 0.01 to avoid log(0)
            mask = errors > 0.01
            if np.sum(mask) > 5:
                log_K = np.log(K_arr[mask])
                log_err = np.log(errors[mask])
                try:
                    alpha, c = np.polyfit(log_K, log_err, 1)
                    return -alpha  # negative because we want positive alpha for decay
                except:
                    return 0.0
            return 0.0

        alpha_seq = fit_alpha(seq_errors, "sequential")
        alpha_mp = fit_alpha(mp_errors, "matching_pursuit")
        alpha_topk = fit_alpha(topk_errors, "top_k")

        print(f"  Scaling exponents (error ~ K^{{-alpha}}):")
        print(f"    Sequential:      alpha = {alpha_seq:.4f}")
        print(f"    Matching pursuit: alpha = {alpha_mp:.4f}")
        print(f"    Top-K:           alpha = {alpha_topk:.4f}")

        # How many zeros for error < 0.5?
        for label, errors in [("Sequential", seq_errors), ("MP", mp_errors), ("Top-K", topk_errors)]:
            achieved = np.where(errors < 0.5)[0]
            if len(achieved) > 0:
                print(f"    {label}: error < 0.5 first at K = {achieved[0]+1}")
            else:
                print(f"    {label}: error never < 0.5 with {len(K_values)} zeros (min err = {errors.min():.4f})")

    return {}


# =============================================================================
# EXPERIMENT 9: CESARO/RIESZ SMOOTHING + SPARSE SELECTION COMBINED
# =============================================================================

def experiment9_smoothing_plus_sparse(zeros, test_cases):
    """
    Combine smoothing kernels (which reduce Gibbs-like oscillation from
    truncation) with sparse selection. Does smoothing + sparsity help?
    """
    print("\n" + "="*75)
    print("EXPERIMENT 9: SMOOTHING KERNELS + SPARSE SELECTION")
    print("="*75)

    kernels = {
        'none': lambda k, K: 1.0,
        'cesaro': lambda k, K: 1.0 - k / (K + 1),
        'riesz': lambda k, K: 1.0 - (k / (K + 1))**2,
        'gaussian': lambda k, K: np.exp(-2 * (k / (K + 1))**2),
        'hann': lambda k, K: 0.5 * (1 + np.cos(np.pi * k / (K + 1))),
    }

    for name, x, true_pi in test_cases:
        print(f"\n--- x = {x} (pi(x) = {true_pi}) ---")
        R_x = float(R_function(mpf(x)))
        contribs = np.array([zero_contribution_approx(g, x) for g in zeros])

        print(f"  {'Kernel':>10} | {'K=10 err':>10} | {'K=20 err':>10} | {'K=50 err':>10} | {'K=100 err':>10}")
        print("  " + "-"*60)

        for kname, kfunc in kernels.items():
            errs = []
            for K in [10, 20, 50, 100]:
                if K > len(zeros):
                    errs.append(float('nan'))
                    continue
                weighted_sum = sum(kfunc(k, K) * contribs[k] for k in range(K))
                est = R_x + weighted_sum
                errs.append(abs(est - true_pi))
            print(f"  {kname:>10} | {errs[0]:10.4f} | {errs[1]:10.4f} | {errs[2]:10.4f} | {errs[3]:10.4f}")

    return {}


# =============================================================================
# MAIN
# =============================================================================

def main():
    print("="*75)
    print("SESSION 7: COMPRESSED SENSING AND SPARSE RECOVERY FOR ZETA ZEROS")
    print("="*75)
    print(f"Started: {time.strftime('%Y-%m-%d %H:%M:%S')}")

    overall_start = time.time()
    all_results = {}

    # Precompute zeros
    zeros = get_zeta_zeros(150)

    # Test cases
    test_cases = [
        ("x=1000", 1000, int(primepi(1000))),
        ("x=10000", 10000, int(primepi(10000))),
        ("x=100000", 100000, int(primepi(100000))),
    ]

    test_cases_small = [
        ("x=1000", 1000, int(primepi(1000))),
        ("x=10000", 10000, int(primepi(10000))),
    ]

    print(f"\n  Test cases:")
    for name, x, pi_x in test_cases:
        print(f"    {name}: pi({x}) = {pi_x}")

    # Run experiments
    print("\n")

    r1 = experiment1_sparse_selection(zeros, test_cases)
    all_results['exp1_sparse_selection'] = r1

    r2 = experiment2_matching_pursuit(zeros, test_cases)
    all_results['exp2_matching_pursuit'] = r2

    r3 = experiment3_importance_sampling(zeros, test_cases)
    all_results['exp3_importance_sampling'] = r3

    r4 = experiment4_optimal_selection(zeros, test_cases_small)
    all_results['exp4_optimal_selection'] = r4

    r5 = experiment5_resonances(zeros, test_cases)
    all_results['exp5_resonances'] = r5

    r6 = experiment6_prony_esprit(test_cases)
    all_results['exp6_prony_esprit'] = r6

    r7 = experiment7_adaptive_crossval(zeros, test_cases)
    all_results['exp7_crossval'] = r7

    r8 = experiment8_error_scaling(zeros)
    all_results['exp8_scaling'] = r8

    r9 = experiment9_smoothing_plus_sparse(zeros, test_cases)
    all_results['exp9_smoothing_sparse'] = r9

    elapsed = time.time() - overall_start
    print(f"\n{'='*75}")
    print(f"ALL EXPERIMENTS COMPLETE in {elapsed:.1f}s")
    print(f"{'='*75}")

if __name__ == "__main__":
    main()

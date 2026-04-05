"""
Spectral Compression of Zeta Zero Sum
======================================

The explicit formula for pi(x):
    pi(x) = R(x) - sum_rho R(x^rho) - 1/ln(2) + integral term

where the sum runs over ALL nontrivial zeros rho = 1/2 + i*gamma of zeta(s).

For exact pi(x), we need the zero sum accurate to +/- 0.5.
Standard approach: ~x^{1/2}/log(x) zeros needed => far too many.

THIS EXPERIMENT tests whether the sum is compressible:
  1. Top-k zeros by magnitude of contribution
  2. Random subsampling
  3. Hierarchical grouping of nearby zeros
  4. Low-rank approximation of "zero matrix" Z[i,j] = x_i^{rho_j}/rho_j
  5. Scaling analysis: how does the number of required zeros grow with x?

Key question: does K(x) (zeros needed for |error| < 0.5) grow sub-linearly in x?
"""

import numpy as np
from mpmath import mp, mpf, log as mplog, exp as mpexp, li
import sympy
import time
import os

# ---------------------------------------------------------------------------
# Load precomputed zeros from data files
# ---------------------------------------------------------------------------

DATA_DIR = "/apps/aplikacijos/prime-research/data"

def load_zeros(num_zeros=500):
    """Load imaginary parts of zeta zeros from data files."""
    fname = os.path.join(DATA_DIR, f"zeta_zeros_{num_zeros}.txt")
    if not os.path.exists(fname):
        # Fall back to largest available
        for n in [1000, 500, 300, 200]:
            fname = os.path.join(DATA_DIR, f"zeta_zeros_{n}.txt")
            if os.path.exists(fname):
                break
    gammas = []
    with open(fname) as f:
        for line in f:
            line = line.strip()
            if line:
                gammas.append(float(line))
    return np.array(gammas[:num_zeros])


# ---------------------------------------------------------------------------
# Core: compute the zero-sum contribution for a single x
# ---------------------------------------------------------------------------

def zero_contributions(x, gammas):
    """
    Compute individual zero contributions to the explicit formula.

    Each conjugate pair (rho, rho_bar) contributes:
        c_k = 2 * Re[ x^rho / (rho * log(x)) ]
            = 2 * x^{1/2} * cos(gamma_k * log(x)) / (|rho_k| * log(x))

    Returns array of contributions (one per zero).
    """
    log_x = np.log(x)
    sqrt_x = x ** 0.5
    phases = gammas * log_x
    rho_abs = np.sqrt(0.25 + gammas**2)
    contribs = 2.0 * sqrt_x * np.cos(phases) / (rho_abs * log_x)
    return contribs


def R_approx(x):
    """Gram series R(x) approx: li(x) - 0.5*li(sqrt(x))."""
    mp.dps = 30
    x = float(x)
    return float(li(mpf(x)) - mpf('0.5') * li(mpf(x) ** mpf('0.5')))


def full_zero_sum(x, gammas):
    """Sum all zero contributions."""
    return np.sum(zero_contributions(x, gammas))


# ---------------------------------------------------------------------------
# Experiment 1: Top-k by magnitude
# ---------------------------------------------------------------------------

def experiment_topk(x_values, gammas, pi_values):
    """
    For each x, rank zeros by |contribution| and keep only top-k.
    Measure error vs k.
    """
    print("=" * 70)
    print("EXPERIMENT 1: Top-k zeros by magnitude of contribution")
    print("=" * 70)

    for x, pi_x in zip(x_values, pi_values):
        R_x = R_approx(x)
        target_correction = R_x - pi_x  # this is what the zero sum should equal
        contribs = zero_contributions(x, gammas)
        full_sum = np.sum(contribs)

        # Sort by |contribution|
        order = np.argsort(-np.abs(contribs))
        sorted_contribs = contribs[order]

        print(f"\n  x = {x:>8}, pi(x) = {pi_x}, R(x) = {R_x:.4f}, "
              f"target correction = {target_correction:.4f}")
        print(f"  Full sum ({len(gammas)} zeros) = {full_sum:.4f}, "
              f"residual = {abs(R_x - full_sum - pi_x):.4f}")

        # Find minimum k for |error| < 0.5
        cumsum = np.cumsum(sorted_contribs)
        found = False
        for k_check in [5, 10, 20, 50, 100, 200, len(gammas)]:
            k_check = min(k_check, len(gammas))
            error = abs(R_x - cumsum[k_check - 1] - pi_x)
            print(f"    top-{k_check:>3} zeros: sum = {cumsum[k_check-1]:>10.4f}, "
                  f"|error| = {error:.4f}"
                  + (" ***" if error < 0.5 and not found else ""))
            if error < 0.5 and not found:
                found = True


# ---------------------------------------------------------------------------
# Experiment 2: Random subsampling
# ---------------------------------------------------------------------------

def experiment_random_subsample(x_values, gammas, pi_values, n_trials=20):
    """
    Uniformly subsample k zeros, scale the sum by N/k.
    Measure error statistics.
    """
    print("\n" + "=" * 70)
    print("EXPERIMENT 2: Random subsampling with rescaling")
    print("=" * 70)

    N = len(gammas)
    for x, pi_x in zip(x_values, pi_values):
        R_x = R_approx(x)
        contribs = zero_contributions(x, gammas)
        target_correction = R_x - pi_x

        print(f"\n  x = {x}, pi(x) = {pi_x}, target correction = {target_correction:.4f}")

        for k in [10, 25, 50, 100, 200]:
            if k > N:
                continue
            errors = []
            for _ in range(n_trials):
                idx = np.random.choice(N, size=k, replace=False)
                # Scale by N/k to estimate full sum
                approx_sum = np.sum(contribs[idx]) * (N / k)
                errors.append(abs(R_x - approx_sum - pi_x))
            errors = np.array(errors)
            print(f"    k={k:>3}: mean |error| = {errors.mean():.4f}, "
                  f"median = {np.median(errors):.4f}, "
                  f"min = {errors.min():.4f}, "
                  f"P(|e|<0.5) = {np.mean(errors < 0.5):.0%}")


# ---------------------------------------------------------------------------
# Experiment 3: Hierarchical grouping
# ---------------------------------------------------------------------------

def experiment_hierarchical(x_values, gammas, pi_values):
    """
    Group consecutive zeros into blocks, represent each block by its sum.
    Test whether block sums are predictable from the block's gamma range.
    """
    print("\n" + "=" * 70)
    print("EXPERIMENT 3: Hierarchical grouping of zeros")
    print("=" * 70)

    N = len(gammas)
    for x, pi_x in zip(x_values, pi_values):
        R_x = R_approx(x)
        contribs = zero_contributions(x, gammas)
        full_sum = np.sum(contribs)
        target_correction = R_x - pi_x

        print(f"\n  x = {x}, pi(x) = {pi_x}")

        for block_size in [2, 5, 10, 20, 50]:
            n_blocks = N // block_size
            block_sums = []
            for b in range(n_blocks):
                block = contribs[b * block_size:(b + 1) * block_size]
                block_sums.append(np.sum(block))
            block_sums = np.array(block_sums)

            total_from_blocks = np.sum(block_sums)
            # Measure cancellation within blocks
            mean_block_magnitude = np.mean(np.abs(block_sums))
            mean_individual_magnitude = np.mean(np.abs(contribs[:n_blocks * block_size]))
            cancellation_ratio = mean_block_magnitude / (block_size * mean_individual_magnitude)

            # Can we approximate by keeping only largest blocks?
            block_order = np.argsort(-np.abs(block_sums))
            top_half = np.sum(block_sums[block_order[:n_blocks // 2]])
            error_half = abs(R_x - top_half - pi_x)

            print(f"    block_size={block_size:>2}: {n_blocks:>3} blocks, "
                  f"cancellation={cancellation_ratio:.4f}, "
                  f"top-half error={error_half:.4f}")


# ---------------------------------------------------------------------------
# Experiment 4: Low-rank approximation of zero matrix
# ---------------------------------------------------------------------------

def experiment_lowrank(x_values, gammas, pi_values):
    """
    Build the matrix Z[i,j] = contribution of zero j to x_i.
    Check its singular values -- if it has low effective rank,
    we can compress the computation.
    """
    print("\n" + "=" * 70)
    print("EXPERIMENT 4: Low-rank approximation of zero matrix")
    print("=" * 70)

    # Build matrix: rows = x values, columns = zeros
    x_dense = np.arange(
        int(x_values[0]),
        min(int(x_values[-1]) + 1, int(x_values[0]) + 500)
    )
    # Limit to manageable size
    if len(x_dense) > 500:
        x_dense = np.linspace(x_values[0], x_values[-1], 500)

    max_zeros = min(len(gammas), 200)
    g = gammas[:max_zeros]

    print(f"\n  Matrix dimensions: {len(x_dense)} x-values x {max_zeros} zeros")

    # Build matrix
    Z = np.zeros((len(x_dense), max_zeros))
    for i, x in enumerate(x_dense):
        Z[i, :] = zero_contributions(x, g)

    # SVD
    U, S, Vt = np.linalg.svd(Z, full_matrices=False)

    print(f"  Top 20 singular values:")
    for i in range(min(20, len(S))):
        print(f"    sigma_{i+1:>2} = {S[i]:.6f}")

    # Effective rank at various thresholds
    total_energy = np.sum(S**2)
    cumulative = np.cumsum(S**2) / total_energy

    print(f"\n  Effective rank (fraction of energy captured):")
    for frac in [0.5, 0.8, 0.9, 0.95, 0.99, 0.999]:
        rank = np.searchsorted(cumulative, frac) + 1
        print(f"    {frac*100:.1f}%: rank {rank}/{len(S)}")

    # Reconstruction error with rank-k approximation
    # For each target x in x_values, compute error using rank-k approx
    print(f"\n  Reconstruction error (for pi(x) rounding) using rank-k approx:")

    # Find indices of x_values in x_dense (approximately)
    test_indices = []
    for xv in x_values:
        idx = np.argmin(np.abs(x_dense - xv))
        test_indices.append(idx)

    for rank in [1, 2, 5, 10, 20, 50, 100]:
        if rank > len(S):
            continue
        Z_approx = U[:, :rank] @ np.diag(S[:rank]) @ Vt[:rank, :]
        max_err = 0
        for idx, (xv, pi_x) in zip(test_indices, zip(x_values, pi_values)):
            row_err = np.abs(Z[idx, :].sum() - Z_approx[idx, :].sum())
            max_err = max(max_err, row_err)
        # Also overall Frobenius
        frob_err = np.linalg.norm(Z - Z_approx, 'fro') / np.linalg.norm(Z, 'fro')
        print(f"    rank {rank:>3}: max sum-error = {max_err:.6f}, "
              f"relative Frobenius = {frob_err:.6f}")


# ---------------------------------------------------------------------------
# Experiment 5: Scaling analysis -- K(x) vs x
# ---------------------------------------------------------------------------

def experiment_scaling(gammas):
    """
    For increasing x, determine the minimum K such that using the top-K
    zeros (by contribution magnitude) gives |error| < 0.5.

    Key question: does K(x) grow as O(x^alpha) with alpha < 0.5?
    """
    print("\n" + "=" * 70)
    print("EXPERIMENT 5: Scaling -- how many zeros needed vs x?")
    print("=" * 70)

    # We can only test modest x since we have 500 zeros
    # For larger x the full 500-zero sum itself may not be accurate
    test_x = [50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000]

    print(f"\n  {'x':>8}  {'pi(x)':>8}  {'K_topk':>7}  {'K_seq':>6}  "
          f"{'K_topk/sqrt(x)':>14}  {'K_topk/x^0.33':>14}  {'full_err':>10}")

    results = []
    for x in test_x:
        pi_x = int(sympy.primepi(x))
        R_x = R_approx(x)
        target = R_x - pi_x  # What the zero sum should equal

        contribs = zero_contributions(x, gammas)
        full_sum = np.sum(contribs)
        full_err = abs(R_x - full_sum - pi_x)

        # Method A: top-k by magnitude
        order = np.argsort(-np.abs(contribs))
        cumsum_topk = np.cumsum(contribs[order])
        K_topk = len(gammas)  # default: all needed
        for k in range(1, len(gammas) + 1):
            if abs(R_x - cumsum_topk[k - 1] - pi_x) < 0.5:
                K_topk = k
                break

        # Method B: sequential (first K zeros)
        cumsum_seq = np.cumsum(contribs)
        K_seq = len(gammas)
        for k in range(1, len(gammas) + 1):
            if abs(R_x - cumsum_seq[k - 1] - pi_x) < 0.5:
                K_seq = k
                break

        ratio_sqrt = K_topk / (x ** 0.5) if K_topk < len(gammas) else float('nan')
        ratio_cube = K_topk / (x ** 0.333) if K_topk < len(gammas) else float('nan')

        marker = "" if K_topk < len(gammas) else " (all)"
        print(f"  {x:>8}  {pi_x:>8}  {K_topk:>7}{marker:>5}  {K_seq:>6}  "
              f"{ratio_sqrt:>14.4f}  {ratio_cube:>14.4f}  {full_err:>10.4f}")

        results.append((x, K_topk, K_seq, full_err))

    # Fit power law K_topk ~ C * x^alpha
    valid = [(x, k) for x, k, _, _ in results if k < len(gammas)]
    if len(valid) >= 3:
        log_x = np.log([v[0] for v in valid])
        log_k = np.log([v[1] for v in valid])
        alpha, log_C = np.polyfit(log_x, log_k, 1)
        C = np.exp(log_C)
        print(f"\n  Power law fit: K(x) ~ {C:.2f} * x^{alpha:.4f}")
        print(f"  (For reference: x^0.5 exponent would be 0.5, polylog would be ~0)")
    else:
        print(f"\n  Not enough data points with K < {len(gammas)} for power-law fit.")

    return results


# ---------------------------------------------------------------------------
# Experiment 6: GUE-motivated fast summation
# ---------------------------------------------------------------------------

def experiment_gue_structure(gammas):
    """
    GUE statistics imply specific pair correlations between zeros.
    Test if these correlations create exploitable structure in partial sums.
    """
    print("\n" + "=" * 70)
    print("EXPERIMENT 6: GUE statistics and partial sum behavior")
    print("=" * 70)

    # Normalized spacings
    spacings = np.diff(gammas)
    mean_spacing = np.mean(spacings)
    normalized_spacings = spacings / mean_spacing

    print(f"\n  Zero spacing statistics (N={len(gammas)} zeros):")
    print(f"    Mean spacing: {mean_spacing:.4f}")
    print(f"    Std of normalized spacings: {np.std(normalized_spacings):.4f}")
    print(f"    Min normalized spacing: {np.min(normalized_spacings):.4f}")
    print(f"    (GUE prediction for std: ~0.42, min repulsion near 0)")

    # For a few x values, check if partial sum variance matches random-walk prediction
    print(f"\n  Partial sum variance vs random walk prediction:")
    for x in [100, 1000, 10000]:
        contribs = zero_contributions(x, gammas)
        partial_sums = np.cumsum(contribs)

        # For random walk: Var(S_K) = K * Var(c_k)
        # With GUE correlations, it should differ
        var_contribs = np.var(contribs)
        for K in [20, 50, 100, 200]:
            if K > len(gammas):
                continue
            actual_var_S_K = partial_sums[K - 1] ** 2  # single realization
            rw_prediction = K * var_contribs
            print(f"    x={x:>6}, K={K:>3}: "
                  f"|S_K|^2 = {actual_var_S_K:.4f}, "
                  f"K*Var(c) = {rw_prediction:.4f}, "
                  f"ratio = {actual_var_S_K / rw_prediction:.4f}" if rw_prediction > 0
                  else f"    x={x:>6}, K={K:>3}: degenerate")

    # Test: pair correlation of contributions
    print(f"\n  Pair correlation of contributions at x=10000:")
    contribs = zero_contributions(10000, gammas)
    autocorr_lags = [1, 2, 3, 5, 10, 20]
    centered = contribs - np.mean(contribs)
    var_c = np.var(contribs)
    for lag in autocorr_lags:
        if lag < len(contribs):
            corr = np.mean(centered[:-lag] * centered[lag:]) / var_c
            print(f"    lag {lag:>2}: autocorrelation = {corr:.6f}")


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

def summarize(scaling_results, gammas):
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)

    print(f"""
  Zeros available: {len(gammas)}
  Gamma range: [{gammas[0]:.2f}, {gammas[-1]:.2f}]

  KEY FINDINGS:
  """)

    valid = [(x, k) for x, k, _, _ in scaling_results if k < len(gammas)]
    if len(valid) >= 3:
        log_x = np.log([v[0] for v in valid])
        log_k = np.log([v[1] for v in valid])
        alpha, log_C = np.polyfit(log_x, log_k, 1)
        if alpha > 0.4:
            verdict = "NEGATIVE: K(x) grows roughly as x^{:.2f}, close to sqrt(x) barrier.".format(alpha)
        elif alpha > 0.1:
            verdict = "INTERESTING: K(x) grows as x^{:.2f}, slower than sqrt(x) but not polylog.".format(alpha)
        else:
            verdict = "PROMISING: K(x) grows very slowly (x^{:.2f}), possibly polylog!".format(alpha)
        print(f"  Scaling exponent: K(x) ~ x^{alpha:.4f}")
        print(f"  Verdict: {verdict}")
    else:
        print("  Insufficient data for scaling analysis (all x needed all zeros).")
        print("  Verdict: NEGATIVE -- even 500 zeros insufficient for moderate x.")

    print(f"""
  INTERPRETATION:
  - If exponent ~0.5: confirms the sqrt(x) barrier, no compression possible
  - If exponent <0.5: partial compression exists, but probably not enough
  - Only exponent ~0 (polylog) would indicate a breakthrough
  - Low-rank analysis shows whether the zero matrix is compressible across x values
  """)


# ===========================================================================
# Main
# ===========================================================================

if __name__ == "__main__":
    np.random.seed(42)

    print("=" * 70)
    print("SPECTRAL COMPRESSION OF ZETA ZERO SUM")
    print("Testing whether the explicit formula zero sum is compressible")
    print("=" * 70)

    # Load zeros
    gammas = load_zeros(500)
    print(f"\nLoaded {len(gammas)} zeta zeros, gamma in [{gammas[0]:.4f}, {gammas[-1]:.4f}]")

    # Test x values for experiments 1-4
    x_values = np.array([100, 500, 1000, 5000, 10000])
    pi_values = np.array([int(sympy.primepi(int(x))) for x in x_values])

    print(f"\nTest x values: {x_values.tolist()}")
    print(f"Corresponding pi(x): {pi_values.tolist()}")

    t0 = time.time()

    experiment_topk(x_values, gammas, pi_values)
    experiment_random_subsample(x_values, gammas, pi_values)
    experiment_hierarchical(x_values, gammas, pi_values)
    experiment_lowrank(x_values, gammas, pi_values)
    scaling_results = experiment_scaling(gammas)
    experiment_gue_structure(gammas)
    summarize(scaling_results, gammas)

    elapsed = time.time() - t0
    print(f"\nTotal runtime: {elapsed:.1f}s")

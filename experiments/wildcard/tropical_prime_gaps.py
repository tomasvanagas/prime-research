"""
Tropical / Min-Plus Structure in Prime Gaps: Comprehensive Test

HYPOTHESIS: If prime gaps g(n) = p(n+1) - p(n) have hidden low-dimensional
structure in some transformed space, then:
  1. Hankel matrices of gap sequences would have rapidly decaying singular values
  2. Delay embeddings would reveal a low-dimensional attractor
  3. Short recurrences would predict gaps accurately
  4. Min-plus matrix exponentiation could fast-forward partial sums

If ANY of these hold, we could compute p(n) = 2 + sum_{k=1}^{n-1} g(k)
via repeated squaring in O(d^3 log n) operations.

PRIOR WORK (from CLOSED_PATHS.md):
  - Tropical/min-plus: 7+ tests, all FAIL. Min-plus is optimization not counting.
  - Gap recurrence (linear, AR(1..50)): R^2 < 0 for all orders.
  - Dynamical system on gaps: ~90% random, MI(g_n;g_{n+1})~0.3 bits, no attractor.
  - Hankel matrices: full rank in prior tests.

THIS EXPERIMENT adds: SVD spectrum analysis at multiple scales, normalized gap
analysis, polynomial recurrence fitting, explicit min-plus matrix test, and
rigorous dimension estimation via correlation dimension.

Usage:
  python tropical_prime_gaps.py [--nprimes N] [--quick]
"""

import numpy as np
import time
import argparse
import sys
from scipy import linalg
from scipy.spatial.distance import pdist
from sympy import prime, nextprime

# ---------------------------------------------------------------------------
# 0. Generate prime gaps
# ---------------------------------------------------------------------------

def generate_gaps(n_primes):
    """Generate first n_primes primes and their gaps."""
    t0 = time.time()
    primes = []
    p = 2
    for _ in range(n_primes):
        primes.append(p)
        p = nextprime(p)
    primes = np.array(primes, dtype=np.int64)
    gaps = np.diff(primes)
    elapsed = time.time() - t0
    print(f"Generated {n_primes} primes in {elapsed:.2f}s")
    print(f"  Range: p(1)={primes[0]} to p({n_primes})={primes[-1]}")
    print(f"  Gap stats: mean={gaps.mean():.2f}, std={gaps.std():.2f}, "
          f"min={gaps.min()}, max={gaps.max()}")
    return primes, gaps

# ---------------------------------------------------------------------------
# 1. Hankel matrix SVD analysis
# ---------------------------------------------------------------------------

def hankel_svd_analysis(gaps, window_sizes=None):
    """
    Build Hankel matrices from the gap sequence and analyze SVD spectrum.
    A rapidly decaying spectrum would indicate low-rank (compressible) structure.
    """
    print("\n" + "="*70)
    print("EXPERIMENT 1: Hankel Matrix SVD Analysis")
    print("="*70)

    if window_sizes is None:
        window_sizes = [50, 100, 200, 500]

    results = {}
    for w in window_sizes:
        if 2 * w > len(gaps):
            print(f"  Window {w}: skipped (need {2*w} gaps, have {len(gaps)})")
            continue

        # Build Hankel matrix: H[i,j] = gaps[i+j]
        n_rows = w
        n_cols = w
        H = np.zeros((n_rows, n_cols))
        for i in range(n_rows):
            for j in range(n_cols):
                H[i, j] = gaps[i + j]

        svs = linalg.svdvals(H)
        svs_norm = svs / svs[0]  # normalize by largest

        # Effective rank: number of SVs > 1% of max
        rank_1pct = np.sum(svs_norm > 0.01)
        rank_10pct = np.sum(svs_norm > 0.10)
        # Energy concentration: fraction of Frobenius norm in top k
        total_energy = np.sum(svs**2)
        energy_top5 = np.sum(svs[:5]**2) / total_energy
        energy_top10 = np.sum(svs[:10]**2) / total_energy
        energy_top20 = np.sum(svs[:20]**2) / total_energy

        # Compare to random: generate IID sequence with same mean/std
        rng = np.random.default_rng(42)
        rand_gaps = rng.choice(gaps[:2*w], size=2*w, replace=True)
        H_rand = np.zeros((n_rows, n_cols))
        for i in range(n_rows):
            for j in range(n_cols):
                H_rand[i, j] = rand_gaps[i + j]
        svs_rand = linalg.svdvals(H_rand)
        svs_rand_norm = svs_rand / svs_rand[0]
        rand_rank_1pct = np.sum(svs_rand_norm > 0.01)
        rand_energy_top10 = np.sum(svs_rand[:10]**2) / np.sum(svs_rand**2)

        results[w] = {
            'svs': svs_norm,
            'rank_1pct': rank_1pct,
            'rank_10pct': rank_10pct,
            'energy_top5': energy_top5,
            'energy_top10': energy_top10,
            'energy_top20': energy_top20,
            'rand_rank_1pct': rand_rank_1pct,
            'rand_energy_top10': rand_energy_top10,
        }

        print(f"\n  Window size w={w}:")
        print(f"    Effective rank (SV > 1% of max):  {rank_1pct}/{w} "
              f"(random: {rand_rank_1pct}/{w})")
        print(f"    Effective rank (SV > 10% of max): {rank_10pct}/{w}")
        print(f"    Energy in top 5 SVs:  {energy_top5:.4f} "
              f"(random: {np.sum(svs_rand[:5]**2)/np.sum(svs_rand**2):.4f})")
        print(f"    Energy in top 10 SVs: {energy_top10:.4f} "
              f"(random: {rand_energy_top10:.4f})")
        print(f"    Energy in top 20 SVs: {energy_top20:.4f}")
        print(f"    SV decay (first 10): {svs_norm[:10].round(4)}")

    return results

# ---------------------------------------------------------------------------
# 2. Delay embedding / Takens theorem
# ---------------------------------------------------------------------------

def delay_embedding_analysis(gaps, max_dim=20, n_points=5000):
    """
    Embed the gap sequence in d-dimensional space and estimate
    the intrinsic dimension of the embedded manifold.
    Uses correlation dimension (Grassberger-Procaccia).
    """
    print("\n" + "="*70)
    print("EXPERIMENT 2: Delay Embedding & Dimension Estimation")
    print("="*70)

    n_use = min(n_points, len(gaps) - max_dim)
    dims_to_test = [2, 3, 5, 7, 10, 15, 20]

    results = {}

    for d in dims_to_test:
        if d > max_dim or n_use < 100:
            continue

        # Build delay embedding vectors
        X = np.zeros((n_use, d))
        for i in range(n_use):
            X[i] = gaps[i:i+d]

        # Correlation dimension via Grassberger-Procaccia
        # Use a subsample for speed
        n_sub = min(2000, n_use)
        idx = np.random.default_rng(42).choice(n_use, n_sub, replace=False)
        X_sub = X[idx]

        dists = pdist(X_sub, metric='chebyshev')
        dists = dists[dists > 0]

        if len(dists) == 0:
            continue

        # Compute correlation integral at multiple radii
        log_eps = np.linspace(np.log(np.percentile(dists, 1)),
                              np.log(np.percentile(dists, 50)), 20)
        eps_vals = np.exp(log_eps)
        C_vals = np.array([np.mean(dists < eps) for eps in eps_vals])
        C_vals = C_vals[C_vals > 0]

        if len(C_vals) < 5:
            print(f"  d={d}: insufficient data for dimension estimate")
            continue

        log_C = np.log(C_vals[:len(eps_vals)])
        log_e = log_eps[:len(C_vals)]

        # Linear fit in the scaling region (middle 60%)
        n_pts = len(log_C)
        i_lo = n_pts // 5
        i_hi = 4 * n_pts // 5
        if i_hi - i_lo < 3:
            i_lo, i_hi = 0, n_pts

        coeffs = np.polyfit(log_e[i_lo:i_hi], log_C[i_lo:i_hi], 1)
        corr_dim = coeffs[0]

        # Compare to what random IID would give (should equal embedding dim)
        rng = np.random.default_rng(123)
        X_rand = np.zeros((n_sub, d))
        rand_gaps = rng.choice(gaps[:n_use+d], size=n_use+d, replace=True)
        for i in range(n_sub):
            ii = idx[i] if idx[i] + d <= len(rand_gaps) else i
            X_rand[i] = rand_gaps[ii:ii+d]

        dists_rand = pdist(X_rand, metric='chebyshev')
        dists_rand = dists_rand[dists_rand > 0]
        C_rand = np.array([np.mean(dists_rand < eps) for eps in eps_vals])
        C_rand = C_rand[C_rand > 0]
        if len(C_rand) >= 5:
            log_Cr = np.log(C_rand[:len(eps_vals)])
            log_er = log_eps[:len(C_rand)]
            n_r = len(log_Cr)
            i_lo_r, i_hi_r = n_r // 5, 4 * n_r // 5
            if i_hi_r - i_lo_r < 3:
                i_lo_r, i_hi_r = 0, n_r
            rand_dim = np.polyfit(log_er[i_lo_r:i_hi_r], log_Cr[i_lo_r:i_hi_r], 1)[0]
        else:
            rand_dim = float('nan')

        results[d] = {
            'corr_dim': corr_dim,
            'rand_dim': rand_dim,
        }

        print(f"  Embedding dim d={d}: correlation dim = {corr_dim:.2f} "
              f"(random baseline: {rand_dim:.2f})")

    # Check if correlation dimension saturates (would indicate finite attractor)
    if len(results) >= 3:
        dims = sorted(results.keys())
        cd_vals = [results[d]['corr_dim'] for d in dims]
        print(f"\n  Correlation dimension vs embedding dimension:")
        print(f"    Embedding dims: {dims}")
        print(f"    Corr dims:      {[f'{v:.2f}' for v in cd_vals]}")
        if cd_vals[-1] > cd_vals[0] * 1.5:
            print(f"    VERDICT: Correlation dim GROWS with embedding dim => NO low-dim attractor")
        else:
            print(f"    VERDICT: Correlation dim SATURATES => possible low-dim attractor at ~{cd_vals[-1]:.1f}")

    return results

# ---------------------------------------------------------------------------
# 3. Recurrence fitting (linear + polynomial)
# ---------------------------------------------------------------------------

def recurrence_fitting(gaps, max_order=20, n_train=None):
    """
    Fit g(n) = f(g(n-1), ..., g(n-k)) for various k.
    Test linear (AR) and quadratic models.
    """
    print("\n" + "="*70)
    print("EXPERIMENT 3: Recurrence Fitting")
    print("="*70)

    if n_train is None:
        n_train = min(50000, len(gaps) - 1000)
    n_test = min(5000, len(gaps) - n_train - max_order)

    if n_test < 100:
        print("  Not enough data for train/test split")
        return {}

    mean_gap = np.mean(gaps[:n_train])
    baseline_mse = np.mean((gaps[max_order:n_train] - mean_gap)**2)
    baseline_mae = np.mean(np.abs(gaps[max_order:n_train] - mean_gap))

    print(f"  Training on {n_train} gaps, testing on {n_test}")
    print(f"  Baseline (predict mean={mean_gap:.2f}): MSE={baseline_mse:.2f}, MAE={baseline_mae:.2f}")

    results = {}

    orders_to_test = [1, 2, 3, 5, 10, 15, 20]

    for k in orders_to_test:
        if k > max_order:
            continue

        # Build training data
        X_train = np.zeros((n_train - k, k))
        y_train = np.zeros(n_train - k)
        for i in range(k, n_train):
            X_train[i-k] = gaps[i-k:i][::-1]  # g(n-1), g(n-2), ..., g(n-k)
            y_train[i-k] = gaps[i]

        X_test = np.zeros((n_test, k))
        y_test = np.zeros(n_test)
        for i in range(n_test):
            idx = n_train + i
            X_test[i] = gaps[idx-k:idx][::-1]
            y_test[i] = gaps[idx]

        # Linear fit
        try:
            coeffs, residuals, rank, sv = np.linalg.lstsq(
                np.column_stack([X_train, np.ones(len(X_train))]),
                y_train, rcond=None
            )
            y_pred_train = X_train @ coeffs[:-1] + coeffs[-1]
            y_pred_test = X_test @ coeffs[:-1] + coeffs[-1]
            mse_train = np.mean((y_train - y_pred_train)**2)
            mse_test = np.mean((y_test - y_pred_test)**2)
            r2_train = 1 - mse_train / baseline_mse
            r2_test = 1 - mse_test / np.mean((y_test - mean_gap)**2)
            mae_test = np.mean(np.abs(y_test - y_pred_test))
        except Exception as e:
            r2_train = r2_test = float('nan')
            mae_test = float('nan')

        # Quadratic fit (only for small k to avoid combinatorial explosion)
        r2_quad_test = float('nan')
        mae_quad_test = float('nan')
        if k <= 5:
            try:
                # Add pairwise products
                from itertools import combinations_with_replacement
                pairs = list(combinations_with_replacement(range(k), 2))
                X_quad_train = np.column_stack([
                    X_train,
                    np.array([X_train[:, i] * X_train[:, j] for i, j in pairs]).T,
                    np.ones(len(X_train))
                ])
                X_quad_test = np.column_stack([
                    X_test,
                    np.array([X_test[:, i] * X_test[:, j] for i, j in pairs]).T,
                    np.ones(len(X_test))
                ])
                coeffs_q, _, _, _ = np.linalg.lstsq(X_quad_train, y_train, rcond=None)
                y_pred_q = X_quad_test @ coeffs_q
                mse_q = np.mean((y_test - y_pred_q)**2)
                r2_quad_test = 1 - mse_q / np.mean((y_test - mean_gap)**2)
                mae_quad_test = np.mean(np.abs(y_test - y_pred_q))
            except Exception:
                pass

        results[k] = {
            'r2_train': r2_train,
            'r2_test': r2_test,
            'mae_test': mae_test,
            'r2_quad_test': r2_quad_test,
            'mae_quad_test': mae_quad_test,
        }

        quad_str = f", Quad R²={r2_quad_test:.4f}" if not np.isnan(r2_quad_test) else ""
        print(f"  Order k={k:2d}: Linear R²(train)={r2_train:.4f}, "
              f"R²(test)={r2_test:.4f}, MAE={mae_test:.2f}{quad_str}")

    return results

# ---------------------------------------------------------------------------
# 4. Normalized gap analysis (g(n)/log(p(n)))
# ---------------------------------------------------------------------------

def normalized_gap_analysis(primes, gaps, window_sizes=None):
    """
    Analyze g(n)/log(p(n)) which should be ~1 on average by Cramer's model.
    Check if normalization reveals hidden structure.
    """
    print("\n" + "="*70)
    print("EXPERIMENT 4: Normalized Gap Analysis (g/log p)")
    print("="*70)

    log_p = np.log(primes[:-1].astype(float))
    norm_gaps = gaps / log_p

    print(f"  Normalized gap stats: mean={norm_gaps.mean():.4f}, "
          f"std={norm_gaps.std():.4f}")

    # Autocorrelation of normalized gaps
    ng = norm_gaps - norm_gaps.mean()
    var = np.var(ng)
    autocorr = []
    for lag in range(1, 21):
        c = np.mean(ng[:-lag] * ng[lag:]) / var
        autocorr.append(c)

    print(f"  Autocorrelation (lags 1-20):")
    print(f"    {[f'{a:.4f}' for a in autocorr]}")
    max_ac = max(abs(a) for a in autocorr)
    print(f"    Max |autocorr| = {max_ac:.4f}")

    # Hankel SVD on normalized gaps
    if window_sizes is None:
        window_sizes = [50, 100, 200]

    for w in window_sizes:
        if 2 * w > len(norm_gaps):
            continue
        H = np.zeros((w, w))
        for i in range(w):
            for j in range(w):
                H[i, j] = norm_gaps[i + j]
        svs = linalg.svdvals(H)
        svs_n = svs / svs[0]
        rank_1pct = np.sum(svs_n > 0.01)
        energy_top10 = np.sum(svs[:10]**2) / np.sum(svs**2)
        print(f"  Hankel SVD (normalized, w={w}): "
              f"rank(1%)={rank_1pct}/{w}, energy(top10)={energy_top10:.4f}")

    # Mutual information between consecutive normalized gaps (binned estimate)
    n_bins = 30
    hist_2d, _, _ = np.histogram2d(norm_gaps[:-1], norm_gaps[1:], bins=n_bins)
    p_xy = hist_2d / hist_2d.sum()
    p_x = p_xy.sum(axis=1)
    p_y = p_xy.sum(axis=0)
    mi = 0
    for i in range(n_bins):
        for j in range(n_bins):
            if p_xy[i, j] > 0 and p_x[i] > 0 and p_y[j] > 0:
                mi += p_xy[i, j] * np.log2(p_xy[i, j] / (p_x[i] * p_y[j]))
    print(f"  MI(ng_n; ng_{'{n+1}'}) = {mi:.4f} bits")

    # Entropy of individual normalized gaps (binned)
    h_x = -np.sum(p_x[p_x > 0] * np.log2(p_x[p_x > 0]))
    print(f"  H(ng_n) = {h_x:.4f} bits (binned, {n_bins} bins)")
    print(f"  MI / H = {mi/h_x:.4f} (fraction of entropy explained)")

    return {
        'mean': float(norm_gaps.mean()),
        'std': float(norm_gaps.std()),
        'max_autocorr': float(max_ac),
        'mi_bits': float(mi),
        'entropy_bits': float(h_x),
        'mi_fraction': float(mi / h_x),
    }

# ---------------------------------------------------------------------------
# 5. Min-plus matrix approach
# ---------------------------------------------------------------------------

def minplus_test(gaps, test_dims=None):
    """
    Test if partial sums S(n) = sum_{k=0}^{n-1} g(k) can be expressed via
    min-plus matrix exponentiation.

    In the min-plus semiring, (a ⊕ b) = min(a,b), (a ⊗ b) = a + b.
    Matrix 'multiplication' uses these operations.

    For this to work for SUMS (not mins), we'd need the PLUS-TIMES semiring
    which is just ordinary arithmetic -- no shortcut there.

    But we test: can the GAP SEQUENCE be generated by a min-plus linear
    recurrence of small dimension d?
    g(n) = min_{j} (A_j + g(n-j)) for some constants A_j
    """
    print("\n" + "="*70)
    print("EXPERIMENT 5: Min-Plus / Tropical Structure")
    print("="*70)

    if test_dims is None:
        test_dims = [2, 3, 5, 10]

    print("  Testing min-plus linear recurrence: g(n) = min_j(A_j + g(n-j))")
    print("  (This would allow O(d^3 log n) computation via repeated squaring)")

    n_test = min(10000, len(gaps))
    g = gaps[:n_test].astype(float)

    for d in test_dims:
        # Try to find A_1,...,A_d such that g(n) ≈ min_j(A_j + g(n-j))
        # Optimal A_j = median of (g(n) - g(n-j)) over training data
        n_train = n_test // 2
        A = np.zeros(d)
        for j in range(1, d + 1):
            diffs = g[j:n_train] - g[:n_train - j]
            A[j-1] = np.median(diffs)

        # Predict on second half
        correct = 0
        total = 0
        abs_errors = []
        for n in range(n_train, n_test):
            candidates = []
            for j in range(1, d + 1):
                if n - j >= 0:
                    candidates.append(A[j-1] + g[n-j])
            if candidates:
                pred = min(candidates)
                err = abs(pred - g[n])
                abs_errors.append(err)
                if err < 0.5:
                    correct += 1
                total += 1

        if total > 0:
            accuracy = correct / total
            mean_err = np.mean(abs_errors)
            print(f"  d={d:2d}: accuracy={accuracy:.4f}, mean|error|={mean_err:.2f}, "
                  f"A={A.round(2)}")
        else:
            print(f"  d={d:2d}: no predictions possible")

    # Also test max-plus version
    print("\n  Testing max-plus: g(n) = max_j(B_j + g(n-j))")
    for d in [5, 10]:
        n_train = n_test // 2
        B = np.zeros(d)
        for j in range(1, d + 1):
            diffs = g[j:n_train] - g[:n_train - j]
            B[j-1] = np.median(diffs)

        correct = 0
        total = 0
        abs_errors = []
        for n in range(n_train, n_test):
            candidates = []
            for j in range(1, d + 1):
                if n - j >= 0:
                    candidates.append(B[j-1] + g[n-j])
            if candidates:
                pred = max(candidates)
                err = abs(pred - g[n])
                abs_errors.append(err)
                if err < 0.5:
                    correct += 1
                total += 1

        if total > 0:
            accuracy = correct / total
            mean_err = np.mean(abs_errors)
            print(f"  d={d:2d}: accuracy={accuracy:.4f}, mean|error|={mean_err:.2f}")

    # KEY THEORETICAL POINT
    print("\n  THEORETICAL NOTE:")
    print("  Min-plus recurrence computes g(n) = min over paths, not sums.")
    print("  p(n) = 2 + SUM of gaps requires the (+, x) semiring = ordinary arithmetic.")
    print("  Even if gaps had min-plus structure, recovering p(n) needs the SUM,")
    print("  which destroys the min-plus shortcut.")

    return {}

# ---------------------------------------------------------------------------
# 6. Summary / Verdict
# ---------------------------------------------------------------------------

def compute_verdict(hankel_results, embed_results, recurrence_results, norm_results):
    """Synthesize all results into a clear verdict."""
    print("\n" + "="*70)
    print("OVERALL VERDICT")
    print("="*70)

    issues = []

    # Check Hankel rank
    if hankel_results:
        worst_w = max(hankel_results.keys())
        r = hankel_results[worst_w]
        if r['rank_1pct'] > worst_w * 0.5:
            issues.append(f"Hankel rank {r['rank_1pct']}/{worst_w} (>50% = high rank)")
        ratio = r['rank_1pct'] / r['rand_rank_1pct'] if r['rand_rank_1pct'] > 0 else 1.0
        print(f"  Hankel rank (w={worst_w}): {r['rank_1pct']}/{worst_w} "
              f"(random: {r['rand_rank_1pct']}/{worst_w}, ratio: {ratio:.2f})")

    # Check embedding dimension
    if embed_results:
        dims = sorted(embed_results.keys())
        if len(dims) >= 2:
            cd_first = embed_results[dims[0]]['corr_dim']
            cd_last = embed_results[dims[-1]]['corr_dim']
            if cd_last > cd_first * 1.3:
                issues.append(f"Correlation dim grows: {cd_first:.1f} -> {cd_last:.1f} (no attractor)")
            print(f"  Correlation dimension: {cd_first:.1f} (d={dims[0]}) -> "
                  f"{cd_last:.1f} (d={dims[-1]})")

    # Check recurrence R^2
    if recurrence_results:
        best_r2 = max(r['r2_test'] for r in recurrence_results.values()
                      if not np.isnan(r['r2_test']))
        best_quad = max((r['r2_quad_test'] for r in recurrence_results.values()
                        if not np.isnan(r['r2_quad_test'])), default=float('nan'))
        if best_r2 < 0.05:
            issues.append(f"Best linear R²={best_r2:.4f} (< 5%)")
        print(f"  Best linear R²(test): {best_r2:.4f}")
        if not np.isnan(best_quad):
            print(f"  Best quadratic R²(test): {best_quad:.4f}")

    # Check MI fraction
    if norm_results:
        mi_frac = norm_results['mi_fraction']
        if mi_frac < 0.1:
            issues.append(f"MI fraction = {mi_frac:.4f} (< 10% of entropy)")
        print(f"  MI(g_n; g_{{n+1}}) / H(g_n) = {mi_frac:.4f}")

    print()
    if len(issues) >= 3:
        print("  VERDICT: NO LOW-DIMENSIONAL STRUCTURE FOUND.")
        print("  Prime gaps are effectively incompressible for fast-forwarding.")
        print("  This confirms prior results (CLOSED_PATHS entries 472, 501, 562).")
    elif len(issues) >= 1:
        print(f"  VERDICT: MARGINAL -- {len(issues)} of 4 tests negative.")
        print(f"  Issues: {'; '.join(issues)}")
    else:
        print("  VERDICT: POSSIBLE STRUCTURE DETECTED -- investigate further!")

    print(f"\n  Failure reasons ({len(issues)}/4 tests negative):")
    for iss in issues:
        print(f"    - {iss}")

    return issues

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Tropical structure in prime gaps")
    parser.add_argument('--nprimes', type=int, default=100000,
                        help='Number of primes to generate (default: 100000)')
    parser.add_argument('--quick', action='store_true',
                        help='Quick mode: fewer primes, smaller windows')
    args = parser.parse_args()

    n_primes = args.nprimes
    if args.quick:
        n_primes = min(n_primes, 20000)

    print("Tropical / Min-Plus Structure in Prime Gaps")
    print("=" * 70)

    primes, gaps = generate_gaps(n_primes)

    # 1. Hankel SVD
    if args.quick:
        hankel_results = hankel_svd_analysis(gaps, window_sizes=[50, 100])
    else:
        hankel_results = hankel_svd_analysis(gaps, window_sizes=[50, 100, 200, 500])

    # 2. Delay embedding
    embed_results = delay_embedding_analysis(
        gaps, max_dim=20,
        n_points=min(10000, len(gaps) - 20))

    # 3. Recurrence fitting
    recurrence_results = recurrence_fitting(gaps, max_order=20)

    # 4. Normalized gap analysis
    norm_results = normalized_gap_analysis(primes, gaps)

    # 5. Min-plus test
    minplus_test(gaps)

    # 6. Verdict
    issues = compute_verdict(hankel_results, embed_results,
                             recurrence_results, norm_results)

    print("\nDone.")
    return len(issues)


if __name__ == '__main__':
    main()

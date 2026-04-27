"""
ATTACK_VECTORS §C2 — Higher-order arithmetic corrections (Conrey-Snaith).

Test orders 4, 5, 6 of zeta-zero n-correlations at N=8000.
Three complementary probes:

  (A) k-th nearest-neighbor spacing distributions P_k(s) for k=0..5.
      P_k probes the (k+2)-point correlation directly (P_5 sees orders
      up to 7-point clusters). Strong-statistics 1D histograms.

  (B) Higher cumulants kappa_n(L) of zero counts in disjoint windows
      for n=2..6 at L in {1, 2, 4, 8, 16, 32, 64, 128}. kappa_n is
      sensitive to n-point correlations; for a GUE DPP the closed form
      is

         log det(I + (e^t - 1) K_L) = sum_n kappa_n(L) t^n / n!

      with sine-kernel K(x, y) = sin(pi (x-y)) / (pi (x-y)).

  (C) Equally-spaced n-tuple density R_n(0, s, 2s, ..., (n-1)s) for
      n=4, 5, 6, swept over s in [0.5, 5]. GUE prediction:
          R_n(s_1,...,s_n) = det[K(s_i - s_j)]_{i,j}.

For each test:
  - Empirical from the N=8000 unfolded zeta zeros.
  - GUE Monte Carlo: eigenvalues of a single N=8000 GUE matrix
    (Hermitian-Wigner; central bulk extracted, unfolded by semicircle).
  - Gap-shuffled null: K=20 seeds preserving the gap-marginal
    (kills higher-order correlations beyond 1-spacing).

Z-scores report whether empirical zeta deviates from either control.
"""
from pathlib import Path
import json
import time
import numpy as np

ROOT = Path('/apps/aplikacijos/prime-research')
HERE = Path(__file__).parent
ZEROS_PATH = ROOT / 'data' / 'zeta_zeros_8000.txt'
RESULTS_JSON = HERE / 'n_correlations_4_5_6_results.json'
RESULTS_MD = HERE / 'n_correlations_4_5_6_results.md'

RNG_SEED = 20260427


# ---------------------------------------------------------------------------
# Loading and unfolding
# ---------------------------------------------------------------------------
def load_zeros(path):
    return np.array([float(line.strip()) for line in open(path) if line.strip()],
                    dtype=np.float64)


def unfold_riemann(gammas):
    """Riemann-von Mangoldt smooth-counting unfolding to unit mean spacing."""
    return (gammas / (2 * np.pi)) * np.log(gammas / (2 * np.pi * np.e)) + 7 / 8


# ---------------------------------------------------------------------------
# GUE Monte Carlo eigenvalues (semicircle unfolded)
# ---------------------------------------------------------------------------
def gue_eigenvalues(N, seed):
    """Sample N GUE eigenvalues, semicircle-unfold to unit mean spacing.
    Use the central 60% of the spectrum to avoid edge effects.
    """
    rng = np.random.default_rng(seed)
    A = (rng.standard_normal((N, N)) + 1j * rng.standard_normal((N, N))) / np.sqrt(2)
    H = (A + A.conj().T) / np.sqrt(2)
    H /= np.sqrt(2 * N)  # so spectral radius -> 1 (Wigner semicircle on [-1, 1])
    evs = np.linalg.eigvalsh(H)
    # Semicircle density rho(x) = (2/pi) sqrt(1 - x^2) on [-1, 1]
    # Unfolded coordinate u(x) = N * integral_{-1}^x rho(t) dt
    # = N * (1/pi) * (x sqrt(1-x^2) + arcsin(x) + pi/2)
    x = np.clip(evs, -0.999999, 0.999999)
    u = (N / np.pi) * (x * np.sqrt(1 - x ** 2) + np.arcsin(x) + np.pi / 2)
    # Take central 60% to avoid edge behaviour (where unfolding is least accurate)
    lo = int(0.20 * N)
    hi = int(0.80 * N)
    u_central = np.sort(u[lo:hi])
    return u_central


def gue_pool_eigenvalues(N_target, batch=2000, seed=RNG_SEED):
    """Build a pool of >= N_target central GUE eigenvalues by batching
    multiple smaller GUE matrices. Each batch contributes ~0.6*batch
    central eigenvalues. Cheaper than one giant diagonalization.
    """
    pool = []
    s = seed
    n_central_per_batch = int(0.6 * batch)
    n_batches = (N_target + n_central_per_batch - 1) // n_central_per_batch
    for _ in range(n_batches):
        u = gue_eigenvalues(batch, s)
        pool.append(u)
        s += 1
    # We can't simply concatenate: each batch has its own [start, end] interval.
    # For statistics, use each batch independently and aggregate.
    return pool  # list of arrays


# ---------------------------------------------------------------------------
# Gap-shuffled null
# ---------------------------------------------------------------------------
def gap_shuffled(unfolded, seed):
    rng = np.random.default_rng(seed)
    gaps = np.diff(unfolded).copy()
    rng.shuffle(gaps)
    out = np.empty_like(unfolded)
    out[0] = unfolded[0]
    out[1:] = unfolded[0] + np.cumsum(gaps)
    return out


# ---------------------------------------------------------------------------
# Test (A): k-th nearest neighbor spacing distributions P_k(s)
# ---------------------------------------------------------------------------
def kth_spacings(unfolded, k):
    """Distribution of (u_{n+k+1} - u_n) for n=0..N-k-2.
    For k=0 = consecutive spacing; k=1 = 2nd-nearest; etc.
    """
    return unfolded[k + 1:] - unfolded[: -(k + 1)]


def histogram(values, bins, weight=None):
    H, _ = np.histogram(values, bins=bins, density=True)
    return H


def test_A_kth_spacings(zeta_unfolded, gue_pool, n_null_seeds=20, k_max=5):
    """Compare zeta P_k(s) to GUE Monte Carlo and gap-shuffled null."""
    s_max = (k_max + 1) * 3  # plenty of headroom — k-th spacing mean ~k+1
    bin_w = 0.10
    bins = np.arange(0, s_max + bin_w, bin_w)
    centers = (bins[:-1] + bins[1:]) / 2

    out = {'centers': centers.tolist(), 'k_max': k_max, 'bin_width': bin_w}

    for k in range(k_max + 1):
        # Empirical zeta
        s_zeta = kth_spacings(zeta_unfolded, k)
        H_zeta = histogram(s_zeta, bins)

        # GUE pool: aggregate over batches
        s_gue_all = []
        for u in gue_pool:
            s_gue_all.append(kth_spacings(u, k))
        s_gue = np.concatenate(s_gue_all)
        H_gue = histogram(s_gue, bins)

        # Gap-shuffled null: K=n_null_seeds
        H_null_runs = []
        for seed in range(n_null_seeds):
            u_null = gap_shuffled(zeta_unfolded, seed=RNG_SEED + seed)
            s_null = kth_spacings(u_null, k)
            H_null_runs.append(histogram(s_null, bins))
        H_null_runs = np.array(H_null_runs)
        H_null_mean = H_null_runs.mean(axis=0)
        H_null_std = H_null_runs.std(axis=0)

        # Z-score per bin: zeta vs null
        # Restrict to bins with H_null_std > 0 and finite-mass region
        with np.errstate(divide='ignore', invalid='ignore'):
            z_vs_null = np.where(H_null_std > 1e-9,
                                 (H_zeta - H_null_mean) / np.maximum(H_null_std, 1e-9),
                                 np.nan)
        # Restrict to high-mass bins
        mass_mask = H_null_mean > 1e-3
        z_vs_null_masked = z_vs_null.copy()
        z_vs_null_masked[~mass_mask] = np.nan

        # KS-type RMS of differences
        rms_zeta_vs_gue = float(np.sqrt(np.mean((H_zeta - H_gue) ** 2)))
        rms_zeta_vs_null = float(np.sqrt(np.mean((H_zeta - H_null_mean) ** 2)))

        # Max |z| in masked region
        abs_z = np.abs(z_vs_null_masked)
        if np.all(np.isnan(abs_z)):
            max_z = float('nan')
            argmax_s = float('nan')
        else:
            max_z = float(np.nanmax(abs_z))
            argmax_bin = int(np.nanargmax(abs_z))
            argmax_s = float(centers[argmax_bin])

        out[f'k={k}'] = {
            'H_zeta': H_zeta.tolist(),
            'H_gue': H_gue.tolist(),
            'H_null_mean': H_null_mean.tolist(),
            'H_null_std': H_null_std.tolist(),
            'z_vs_null': z_vs_null.tolist(),
            'rms_zeta_vs_gue': rms_zeta_vs_gue,
            'rms_zeta_vs_null': rms_zeta_vs_null,
            'max_abs_z': max_z,
            'argmax_z_s': argmax_s,
            'n_zeta_samples': len(s_zeta),
            'n_gue_samples': int(len(s_gue)),
        }
        print(f"  P_{k}: n_zeta={len(s_zeta)}, "
              f"rms vs GUE = {rms_zeta_vs_gue:.4f}, "
              f"rms vs null = {rms_zeta_vs_null:.4f}, "
              f"max|z| vs null = {max_z:.2f} at s={argmax_s:.2f}",
              flush=True)
    return out


# ---------------------------------------------------------------------------
# Test (B): Higher cumulants kappa_n(L)
# ---------------------------------------------------------------------------
def cumulants_in_disjoint_windows(unfolded, L, n_max=6):
    """Empirical cumulants kappa_n for n=1..n_max of the count of
    unfolded zeros in disjoint windows of length L.
    """
    span = unfolded[-1] - unfolded[0]
    n_win = int(span / L)
    if n_win < 30:
        return {n: np.nan for n in range(1, n_max + 1)}, n_win
    starts = unfolded[0] + np.arange(n_win) * L
    idx_lo = np.searchsorted(unfolded, starts, side='left')
    idx_hi = np.searchsorted(unfolded, starts + L, side='left')
    counts = (idx_hi - idx_lo).astype(np.float64)
    mu = counts.mean()
    cents = counts - mu
    # Central moments
    m2 = (cents ** 2).mean()
    m3 = (cents ** 3).mean()
    m4 = (cents ** 4).mean()
    m5 = (cents ** 5).mean()
    m6 = (cents ** 6).mean()
    # Cumulants from central moments
    k1 = mu
    k2 = m2
    k3 = m3
    k4 = m4 - 3 * m2 ** 2
    k5 = m5 - 10 * m3 * m2
    k6 = m6 - 15 * m4 * m2 - 10 * m3 ** 2 + 30 * m2 ** 3
    out = {1: k1, 2: k2, 3: k3, 4: k4, 5: k5, 6: k6}
    return out, n_win


def gue_cumulants_via_pool(gue_pool, Ls, n_max=6):
    """For each L, aggregate cumulants over all GUE batches by treating
    each batch as a separate sample of disjoint windows.
    """
    rec = {L: {n: [] for n in range(1, n_max + 1)} for L in Ls}
    n_wins_total = {L: 0 for L in Ls}
    for u in gue_pool:
        for L in Ls:
            cums, n_win = cumulants_in_disjoint_windows(u, L, n_max=n_max)
            for n in range(1, n_max + 1):
                if not np.isnan(cums[n]):
                    rec[L][n].append((n_win, cums[n]))
            n_wins_total[L] += n_win
    out = {}
    for L in Ls:
        cell = {}
        for n in range(1, n_max + 1):
            samples = rec[L][n]
            if samples:
                # weighted mean by n_win (proxy for sample count)
                weights = np.array([s[0] for s in samples])
                vals = np.array([s[1] for s in samples])
                wm = float(np.sum(weights * vals) / np.sum(weights))
                # standard deviation across batches (unweighted)
                std = float(np.std(vals)) if len(vals) > 1 else 0.0
                cell[n] = wm
                cell[f'{n}_std'] = std
                cell[f'{n}_n'] = len(vals)
            else:
                cell[n] = np.nan
                cell[f'{n}_std'] = np.nan
                cell[f'{n}_n'] = 0
        cell['n_wins_total'] = n_wins_total[L]
        out[L] = cell
    return out


def test_B_cumulants(zeta_unfolded, gue_pool, n_null_seeds=20, n_max=6):
    """Compare empirical zeta kappas to GUE MC and gap-shuffled null."""
    Ls = [1, 2, 4, 8, 16, 32, 64, 128]
    out = {'Ls': Ls, 'n_max': n_max}

    # Empirical zeta
    zeta_cums = {}
    for L in Ls:
        cums, n_win = cumulants_in_disjoint_windows(zeta_unfolded, L, n_max=n_max)
        zeta_cums[L] = {**{int(k): float(v) for k, v in cums.items()},
                        'n_win': n_win}
    out['zeta'] = zeta_cums

    # GUE MC
    out['gue'] = gue_cumulants_via_pool(gue_pool, Ls, n_max=n_max)
    # Convert keys to str for JSON
    out['gue'] = {str(L): {str(k): (float(v) if not np.isnan(v) else None)
                            for k, v in cell.items()}
                   for L, cell in out['gue'].items()}

    # Gap-shuffled null
    null_runs = {L: {n: [] for n in range(1, n_max + 1)} for L in Ls}
    for seed in range(n_null_seeds):
        u_null = gap_shuffled(zeta_unfolded, seed=RNG_SEED + 100 + seed)
        for L in Ls:
            cums, _ = cumulants_in_disjoint_windows(u_null, L, n_max=n_max)
            for n in range(1, n_max + 1):
                if not np.isnan(cums[n]):
                    null_runs[L][n].append(cums[n])
    null_summary = {}
    for L in Ls:
        cell = {}
        for n in range(1, n_max + 1):
            arr = np.array(null_runs[L][n])
            if arr.size > 0:
                cell[n] = float(arr.mean())
                cell[f'{n}_std'] = float(arr.std())
            else:
                cell[n] = None
                cell[f'{n}_std'] = None
        null_summary[L] = cell
    out['null'] = {str(L): {str(k): v for k, v in cell.items()}
                   for L, cell in null_summary.items()}

    # Z-scores: zeta vs null
    print("  cumulants vs gap-shuffled null (z-scores):", flush=True)
    z_table = {}
    for L in Ls:
        row = {}
        for n in range(1, n_max + 1):
            null_mean = null_summary[L][n]
            null_std = null_summary[L].get(f'{n}_std')
            zeta_v = zeta_cums[L].get(n)
            if (null_mean is not None and null_std and null_std > 1e-9
                    and zeta_v is not None and not np.isnan(zeta_v)):
                z = (zeta_v - null_mean) / null_std
            else:
                z = float('nan')
            row[n] = z
        z_table[L] = row
        print(f"    L={L:>3}: " + " ".join(
            f"k_{n}_z={row[n]:+5.2f}" if not np.isnan(row[n]) else f"k_{n}_z=  nan"
            for n in range(2, n_max + 1)), flush=True)
    out['z_vs_null'] = {str(L): {str(k): (float(v) if not np.isnan(v) else None)
                                  for k, v in row.items()}
                        for L, row in z_table.items()}

    return out


# ---------------------------------------------------------------------------
# Test (C): Equally-spaced n-tuple density R_n(0, s, 2s, ..., (n-1)s)
# ---------------------------------------------------------------------------
def K_sine(t):
    """Sine kernel K(t) = sin(pi t) / (pi t), K(0) = 1."""
    out = np.ones_like(t, dtype=float)
    nz = np.abs(t) > 1e-14
    pt = np.pi * t[nz]
    out[nz] = np.sin(pt) / pt
    return out


def gue_R_n_equal_spaced(s, n):
    """GUE prediction R_n(0, s, 2s, ..., (n-1)s) = det[K(s_i - s_j)].
    The matrix M[i, j] = K((j - i) * s) is Toeplitz; size n x n.
    """
    arr = np.arange(n)
    diffs = arr[:, None] - arr[None, :]  # (j - i) takes values -n+1..n-1
    K_vals = K_sine(diffs.astype(float) * s)
    return float(np.linalg.det(K_vals))


def empirical_R_n_equal_spaced(unfolded, s, n, tol):
    """For each ref zero z_i, count n-tuples (i, j_1, j_2, ..., j_{n-1})
    with z_{j_k} - z_i in [k*s - tol/2, k*s + tol/2] for k=1..n-1.

    Density estimator: R_n(0, s, ..., (n-1)s) ~ count / (N_ref * tol^{n-1})
    where N_ref = number of reference zeros with full window coverage.
    """
    # For each k=1..n-1, find the set of zero indices j such that
    # z_j - z_i in [k*s - tol/2, k*s + tol/2]. Since the indexing is
    # implicitly ordered, use binary search.
    N = len(unfolded)
    half = tol / 2
    counts = 0
    n_ref = 0
    last_window_end = unfolded[-1]
    for i in range(N):
        z_i = unfolded[i]
        if z_i + (n - 1) * s + half > last_window_end:
            break
        n_ref += 1
        # For each k=1..n-1, get the slice of zeros in the window
        slices = []
        for k in range(1, n):
            lo = z_i + k * s - half
            hi = z_i + k * s + half
            j_lo = np.searchsorted(unfolded, lo, side='left')
            j_hi = np.searchsorted(unfolded, hi, side='left')
            slices.append(j_hi - j_lo)
        # Number of n-tuples = product of slice counts
        if all(c > 0 for c in slices):
            counts += int(np.prod(slices))
    if n_ref == 0:
        return float('nan'), 0
    density = counts / (n_ref * tol ** (n - 1))
    return float(density), n_ref


def test_C_equal_spaced_R_n(zeta_unfolded, gue_pool, n_null_seeds=20):
    """Sweep s in [0.5, 5.0] for n in {4, 5, 6}. tol scaled with s
    to keep window count ~ similar."""
    s_grid = np.array([0.5, 0.7, 1.0, 1.3, 1.6, 2.0, 2.5, 3.0, 4.0, 5.0])
    out = {'s_grid': s_grid.tolist()}

    for n in [4, 5, 6]:
        rows = []
        # tol large enough to give reasonable counts at large n
        tol = 0.20
        gue_pred = np.array([gue_R_n_equal_spaced(s, n) for s in s_grid])
        emp_zeta = np.array([empirical_R_n_equal_spaced(zeta_unfolded, s, n, tol)[0]
                             for s in s_grid])
        # GUE empirical from pool (concat-like aggregation)
        emp_gue = []
        for s in s_grid:
            density_sum = 0.0
            ref_sum = 0
            for u in gue_pool:
                d, nref = empirical_R_n_equal_spaced(u, s, n, tol)
                if not np.isnan(d):
                    density_sum += d * nref
                    ref_sum += nref
            emp_gue.append(density_sum / ref_sum if ref_sum > 0 else float('nan'))
        emp_gue = np.array(emp_gue)
        # Gap-shuffled null
        null_runs = []
        for seed in range(n_null_seeds):
            u_null = gap_shuffled(zeta_unfolded, seed=RNG_SEED + 200 + seed)
            row = [empirical_R_n_equal_spaced(u_null, s, n, tol)[0] for s in s_grid]
            null_runs.append(row)
        null_runs = np.array(null_runs)
        null_mean = null_runs.mean(axis=0)
        null_std = null_runs.std(axis=0)
        with np.errstate(divide='ignore', invalid='ignore'):
            z_vs_null = (emp_zeta - null_mean) / np.maximum(null_std, 1e-9)
        max_abs_z = float(np.nanmax(np.abs(z_vs_null)))
        argmax_idx = int(np.nanargmax(np.abs(z_vs_null)))
        out[f'n={n}'] = {
            'tol': tol,
            'gue_pred': gue_pred.tolist(),
            'emp_zeta': emp_zeta.tolist(),
            'emp_gue': emp_gue.tolist(),
            'null_mean': null_mean.tolist(),
            'null_std': null_std.tolist(),
            'z_vs_null': z_vs_null.tolist(),
            'max_abs_z': max_abs_z,
            'argmax_s': float(s_grid[argmax_idx]),
        }
        print(f"  R_{n} equal-spaced sweep:", flush=True)
        for i, s in enumerate(s_grid):
            print(f"    s={s:.2f} GUE={gue_pred[i]:+.3f}  "
                  f"emp_zeta={emp_zeta[i]:.3f}  emp_gue={emp_gue[i]:.3f}  "
                  f"null={null_mean[i]:.3f}±{null_std[i]:.3f}  "
                  f"z={z_vs_null[i]:+.2f}", flush=True)
        print(f"    max|z|={max_abs_z:.2f} at s={s_grid[argmax_idx]:.2f}", flush=True)
    return out


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    t0 = time.time()
    print(f"Loading {ZEROS_PATH}...", flush=True)
    gammas = load_zeros(ZEROS_PATH)
    print(f"  N={len(gammas)} zeros, gamma range [{gammas[0]:.4f}, {gammas[-1]:.4f}]",
          flush=True)
    u_zeta = unfold_riemann(gammas)
    print(f"  unfolded mean spacing = {np.diff(u_zeta).mean():.6f}",
          flush=True)

    print("\nGenerating GUE pool (batch=2000, target~8000 central evs)...", flush=True)
    t1 = time.time()
    gue_pool = gue_pool_eigenvalues(N_target=8000, batch=2000, seed=RNG_SEED)
    n_pool = sum(len(u) for u in gue_pool)
    print(f"  Pool: {len(gue_pool)} batches, {n_pool} total central GUE eigenvalues "
          f"({time.time() - t1:.1f}s)", flush=True)

    results = {
        'meta': {
            'N_zeta': len(gammas),
            'gamma_min': float(gammas[0]),
            'gamma_max': float(gammas[-1]),
            'unfolded_mean_spacing': float(np.diff(u_zeta).mean()),
            'gue_pool_batches': len(gue_pool),
            'gue_pool_size': int(n_pool),
            'rng_seed': RNG_SEED,
        }
    }

    print("\n=== TEST (A): k-th nearest-neighbor spacing distributions ===",
          flush=True)
    t1 = time.time()
    results['test_A'] = test_A_kth_spacings(u_zeta, gue_pool,
                                             n_null_seeds=20, k_max=5)
    print(f"  Test A done ({time.time() - t1:.1f}s)", flush=True)

    print("\n=== TEST (B): higher cumulants kappa_n(L) for n=2..6 ===",
          flush=True)
    t1 = time.time()
    results['test_B'] = test_B_cumulants(u_zeta, gue_pool,
                                          n_null_seeds=20, n_max=6)
    print(f"  Test B done ({time.time() - t1:.1f}s)", flush=True)

    print("\n=== TEST (C): equally-spaced R_n for n=4, 5, 6 ===", flush=True)
    t1 = time.time()
    results['test_C'] = test_C_equal_spaced_R_n(u_zeta, gue_pool,
                                                 n_null_seeds=20)
    print(f"  Test C done ({time.time() - t1:.1f}s)", flush=True)

    results['meta']['total_runtime_s'] = float(time.time() - t0)
    RESULTS_JSON.write_text(json.dumps(results, indent=2))
    print(f"\nWROTE {RESULTS_JSON}", flush=True)
    print(f"TOTAL runtime: {time.time() - t0:.1f}s", flush=True)


if __name__ == '__main__':
    main()

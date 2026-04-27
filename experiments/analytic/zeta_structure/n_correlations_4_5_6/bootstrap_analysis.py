"""
Follow-up analysis for §C2 wild swing: the gap-shuffled null is NOT a
proper GUE-comparator (it destroys rigidity, conflating GUE-rigidity
with arithmetic structure). Two further checks:

  1. Bootstrap empirical zeta R_n by resampling reference zeros — get
     empirical zeta's OWN sampling variance.

  2. Pure-GUE z-scores: for each statistic, compute it on K=20
     independent GUE batches of size N=2000 (central 60% kept ~ 1200
     evs each), and use cross-batch std as the GUE noise. Compare
     empirical zeta against this distribution.
"""
from pathlib import Path
import json
import time
import numpy as np

ROOT = Path('/apps/aplikacijos/prime-research')
HERE = Path(__file__).parent
ZEROS_PATH = ROOT / 'data' / 'zeta_zeros_8000.txt'
RESULTS_JSON = HERE / 'bootstrap_analysis_results.json'

RNG_SEED = 20260427_2


def load_zeros(path):
    return np.array([float(line.strip()) for line in open(path) if line.strip()],
                    dtype=np.float64)


def unfold_riemann(gammas):
    return (gammas / (2 * np.pi)) * np.log(gammas / (2 * np.pi * np.e)) + 7 / 8


def K_sine(t):
    out = np.ones_like(t, dtype=float)
    nz = np.abs(t) > 1e-14
    pt = np.pi * t[nz]
    out[nz] = np.sin(pt) / pt
    return out


def gue_R_n_equal_spaced(s, n):
    arr = np.arange(n)
    diffs = arr[:, None] - arr[None, :]
    K_vals = K_sine(diffs.astype(float) * s)
    return float(np.linalg.det(K_vals))


def empirical_R_n_equal_spaced(unfolded, s, n, tol):
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
        slices = []
        for k in range(1, n):
            lo = z_i + k * s - half
            hi = z_i + k * s + half
            j_lo = np.searchsorted(unfolded, lo, side='left')
            j_hi = np.searchsorted(unfolded, hi, side='left')
            slices.append(j_hi - j_lo)
        if all(c > 0 for c in slices):
            counts += int(np.prod(slices))
    if n_ref == 0:
        return float('nan'), 0
    density = counts / (n_ref * tol ** (n - 1))
    return float(density), n_ref


def gue_eigenvalues(N, seed):
    rng = np.random.default_rng(seed)
    A = (rng.standard_normal((N, N)) + 1j * rng.standard_normal((N, N))) / np.sqrt(2)
    H = (A + A.conj().T) / np.sqrt(2)
    H /= np.sqrt(2 * N)
    evs = np.linalg.eigvalsh(H)
    x = np.clip(evs, -0.999999, 0.999999)
    u = (N / np.pi) * (x * np.sqrt(1 - x ** 2) + np.arcsin(x) + np.pi / 2)
    lo = int(0.20 * N)
    hi = int(0.80 * N)
    u_central = np.sort(u[lo:hi])
    return u_central


def kth_spacings(unfolded, k):
    return unfolded[k + 1:] - unfolded[: -(k + 1)]


def cumulants_in_disjoint_windows(unfolded, L, n_max=6):
    span = unfolded[-1] - unfolded[0]
    n_win = int(span / L)
    if n_win < 30:
        return {n: float('nan') for n in range(1, n_max + 1)}, n_win
    starts = unfolded[0] + np.arange(n_win) * L
    idx_lo = np.searchsorted(unfolded, starts, side='left')
    idx_hi = np.searchsorted(unfolded, starts + L, side='left')
    counts = (idx_hi - idx_lo).astype(np.float64)
    mu = counts.mean()
    cents = counts - mu
    m2 = (cents ** 2).mean()
    m3 = (cents ** 3).mean()
    m4 = (cents ** 4).mean()
    m5 = (cents ** 5).mean()
    m6 = (cents ** 6).mean()
    return {1: mu, 2: m2, 3: m3,
            4: m4 - 3 * m2 ** 2,
            5: m5 - 10 * m3 * m2,
            6: m6 - 15 * m4 * m2 - 10 * m3 ** 2 + 30 * m2 ** 3}, n_win


def main():
    t0 = time.time()
    gammas = load_zeros(ZEROS_PATH)
    u_zeta = unfold_riemann(gammas)
    print(f"N_zeta = {len(u_zeta)}", flush=True)

    # Generate K=20 GUE batches of N=2000 (central 60% kept = 1200 each)
    K_BATCHES = 20
    BATCH_SIZE = 2000
    print(f"Generating K={K_BATCHES} GUE batches of N={BATCH_SIZE}...", flush=True)
    t1 = time.time()
    gue_batches = []
    for b in range(K_BATCHES):
        gue_batches.append(gue_eigenvalues(BATCH_SIZE, RNG_SEED + b))
    print(f"  Pool: {K_BATCHES} batches × ~{int(0.6*BATCH_SIZE)} central eigenvalues "
          f"({time.time() - t1:.1f}s)", flush=True)

    # ALSO subdivide empirical zeta into K=20 sub-windows of 400 zeros each
    # (matched size to GUE batches) for matched-noise comparison
    K_ZETA_CHUNKS = 20
    chunk_size = len(u_zeta) // K_ZETA_CHUNKS
    zeta_chunks = []
    for c in range(K_ZETA_CHUNKS):
        chunk = u_zeta[c * chunk_size:(c + 1) * chunk_size]
        # Re-unfold the chunk to local mean spacing 1
        # (its mean spacing already ~1 since each chunk is small)
        zeta_chunks.append(chunk)

    # ---------------------------------------------------------------
    # (1) R_n equal-spaced: bootstrap the EMPIRICAL ZETA + GUE batches
    # ---------------------------------------------------------------
    s_grid = np.array([0.5, 0.7, 1.0, 1.3, 1.6, 2.0, 2.5, 3.0, 4.0, 5.0])
    out = {'s_grid': s_grid.tolist(), 'K_batches': K_BATCHES}

    print("\n--- R_n equal-spaced: full empirical zeta (8000) vs K=20 GUE (1200 each)"
          " + 20 zeta sub-chunks (400 each) ---", flush=True)
    for n in [4, 5, 6]:
        tol = 0.20
        gue_pred = np.array([gue_R_n_equal_spaced(s, n) for s in s_grid])
        # empirical zeta on full 8000
        emp_zeta_full = np.array([empirical_R_n_equal_spaced(u_zeta, s, n, tol)[0]
                                  for s in s_grid])
        # GUE batches: per-batch values
        gue_per_batch = np.zeros((K_BATCHES, len(s_grid)))
        for b, u in enumerate(gue_batches):
            for i, s in enumerate(s_grid):
                d, _ = empirical_R_n_equal_spaced(u, s, n, tol)
                gue_per_batch[b, i] = d if not np.isnan(d) else 0.0
        gue_mean = gue_per_batch.mean(axis=0)
        gue_std = gue_per_batch.std(axis=0)
        # Empirical zeta sub-chunks (matched-size noise estimator)
        zeta_chunk_per_chunk = np.zeros((K_ZETA_CHUNKS, len(s_grid)))
        for c, chunk in enumerate(zeta_chunks):
            for i, s in enumerate(s_grid):
                d, _ = empirical_R_n_equal_spaced(chunk, s, n, tol)
                zeta_chunk_per_chunk[c, i] = d if not np.isnan(d) else 0.0
        zeta_chunk_mean = zeta_chunk_per_chunk.mean(axis=0)
        zeta_chunk_std = zeta_chunk_per_chunk.std(axis=0)
        # Z-scores: empirical zeta full vs GUE batch distribution
        with np.errstate(divide='ignore', invalid='ignore'):
            z_zeta_full_vs_gue = np.where(gue_std > 1e-9,
                                           (emp_zeta_full - gue_mean) / np.maximum(gue_std, 1e-9),
                                           np.nan)
        # Z-scores: zeta chunk mean vs GUE batch (matched scale)
        with np.errstate(divide='ignore', invalid='ignore'):
            z_zeta_chunk_vs_gue = np.where(gue_std > 1e-9,
                                            (zeta_chunk_mean - gue_mean) / np.maximum(gue_std, 1e-9),
                                            np.nan)
        # Z-scores: GUE expected R = R_n^{GUE_pred} from determinant; zeta vs theoretical
        # Standard error of empirical zeta_full estimate ~ zeta_chunk_std / sqrt(K_ZETA_CHUNKS)
        # (since 8000 zeros = 20 chunks of 400 zeros pooled)
        zeta_full_se = zeta_chunk_std / np.sqrt(K_ZETA_CHUNKS)
        with np.errstate(divide='ignore', invalid='ignore'):
            z_zeta_vs_theory = np.where(zeta_full_se > 1e-9,
                                         (emp_zeta_full - gue_pred) / np.maximum(zeta_full_se, 1e-9),
                                         np.nan)

        print(f"\nn={n}:", flush=True)
        for i, s in enumerate(s_grid):
            print(f"  s={s:.2f} | GUE_pred={gue_pred[i]:+.3f} | "
                  f"emp_zeta_full={emp_zeta_full[i]:+.3f} | "
                  f"GUE_batch_mean={gue_mean[i]:+.3f}±{gue_std[i]:.3f} | "
                  f"zeta_chunk_mean={zeta_chunk_mean[i]:+.3f}±{zeta_chunk_std[i]:.3f} | "
                  f"z_full_vs_GUE={z_zeta_full_vs_gue[i]:+5.2f} | "
                  f"z_full_vs_theory={z_zeta_vs_theory[i]:+5.2f}",
                  flush=True)
        max_z_full_gue = float(np.nanmax(np.abs(z_zeta_full_vs_gue)))
        max_z_full_theory = float(np.nanmax(np.abs(z_zeta_vs_theory)))
        argmax_idx_gue = int(np.nanargmax(np.abs(z_zeta_full_vs_gue)))
        argmax_idx_th = int(np.nanargmax(np.abs(z_zeta_vs_theory)))
        print(f"  max|z_full_vs_GUE_batch|={max_z_full_gue:.2f} at s={s_grid[argmax_idx_gue]:.2f}",
              flush=True)
        print(f"  max|z_full_vs_theory|={max_z_full_theory:.2f} at s={s_grid[argmax_idx_th]:.2f}",
              flush=True)

        out[f'n={n}'] = {
            'tol': tol,
            'gue_pred': gue_pred.tolist(),
            'emp_zeta_full': emp_zeta_full.tolist(),
            'gue_mean': gue_mean.tolist(),
            'gue_std': gue_std.tolist(),
            'zeta_chunk_mean': zeta_chunk_mean.tolist(),
            'zeta_chunk_std': zeta_chunk_std.tolist(),
            'zeta_full_se': zeta_full_se.tolist(),
            'z_zeta_full_vs_gue_batch': z_zeta_full_vs_gue.tolist(),
            'z_zeta_full_vs_theory': z_zeta_vs_theory.tolist(),
            'max_z_vs_gue_batch': max_z_full_gue,
            'argmax_s_vs_gue_batch': float(s_grid[argmax_idx_gue]),
            'max_z_vs_theory': max_z_full_theory,
            'argmax_s_vs_theory': float(s_grid[argmax_idx_th]),
        }

    # ---------------------------------------------------------------
    # (2) Cumulants kappa_n(L): empirical zeta vs K=20 GUE batches
    # ---------------------------------------------------------------
    print("\n--- Cumulants kappa_n(L): empirical zeta vs 20 GUE batches ---",
          flush=True)
    Ls = [1, 2, 4, 8, 16, 32, 64]  # L=128 has too few windows per batch
    gue_cums_per_batch = {L: {n: [] for n in range(2, 7)} for L in Ls}
    for u in gue_batches:
        for L in Ls:
            cums, _ = cumulants_in_disjoint_windows(u, L, n_max=6)
            for n in range(2, 7):
                if not np.isnan(cums[n]):
                    gue_cums_per_batch[L][n].append(cums[n])
    # Empirical zeta cumulants
    zeta_cums = {}
    for L in Ls:
        cums, _ = cumulants_in_disjoint_windows(u_zeta, L, n_max=6)
        zeta_cums[L] = {n: float(cums[n]) for n in range(2, 7)}
    # Print z-scores
    out_cums = {}
    for L in Ls:
        row = {}
        print(f"  L={L:>3}: ", end='', flush=True)
        for n in range(2, 7):
            arr = np.array(gue_cums_per_batch[L][n])
            if arr.size > 1:
                gm = arr.mean()
                gs = arr.std()
                zv = (zeta_cums[L][n] - gm) / max(gs, 1e-9)
                row[n] = {
                    'zeta': zeta_cums[L][n],
                    'gue_mean': float(gm),
                    'gue_std': float(gs),
                    'z': float(zv),
                    'gue_n_batches': int(arr.size),
                }
                print(f"k_{n}: zeta={zeta_cums[L][n]:+.3f} GUE={gm:+.3f}±{gs:.3f} "
                      f"z={zv:+.2f}  ", end='', flush=True)
            else:
                row[n] = None
                print(f"k_{n}: NaN  ", end='', flush=True)
        print(flush=True)
        out_cums[L] = row
    out['cumulants_vs_gue_batches'] = {str(k): v for k, v in out_cums.items()}

    out['runtime_s'] = float(time.time() - t0)
    RESULTS_JSON.write_text(json.dumps(out, indent=2))
    print(f"\nWROTE {RESULTS_JSON}", flush=True)
    print(f"TOTAL runtime: {time.time() - t0:.1f}s", flush=True)


if __name__ == '__main__':
    main()

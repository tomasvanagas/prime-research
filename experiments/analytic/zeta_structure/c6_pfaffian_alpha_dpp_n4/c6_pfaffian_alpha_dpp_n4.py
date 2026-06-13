"""
ATTACK_VECTORS §C6 — Pfaffian / alpha-determinantal point process structure of
ζ zeros at order n=4.

Cross-domain ingredient: Pfaffian point processes (Borodin 2009 arXiv:0911.1153
§2.4-2.6, Soshnikov 2000 arXiv:math/0002099 §3, Vere-Jones 1997 alpha-permanent
generalisation).

The Vere-Jones alpha-determinantal interpolation:
    R_n^{alpha}(s_1, ..., s_n) := det_{alpha}[ K_sine(s_i - s_j) ]_{n x n}
where det_alpha is the alpha-determinant
    det_alpha(A) = sum_{sigma in S_n} alpha^{n - c(sigma)} prod_i A[i, sigma(i)]
with c(sigma) the number of cycles in sigma.

Special values:
    alpha = -1  -->  standard determinant (DPP, GUE)
    alpha =  0  -->  product of diagonal (Poisson independence)
    alpha = +1  -->  permanent (PPP, boson)

Pfaffian processes are the matrix-valued generalisation, native to GOE/GSE
(beta = 1, 4) ensembles.

QUESTION: Are the empirical n=4 correlations of ζ zeros best fit by:
  (i)   the canonical sine-kernel determinantal model (alpha = -1, GUE, beta=2)
  (ii)  a Pfaffian process (GOE beta=1 or GSE beta=4)
  (iii) an alpha-DPP with alpha != -1

If (ii) or (iii) at >= 5 sigma over (i), we have a NEW structural fact about
ζ-zero arithmetic content beyond E7.1 / E3.13 / S123.

If (i) is the unique best fit at all alpha in [-1.2, -0.8] and dominates GOE/GSE
to >= 5 sigma, we strengthen E7.1 from "DPP-typical at orders 4-6" to
"DPP-typical AND not Pfaffian-typical at order 4" — a strictly stronger
structural negative.

METHOD:
1. Load 8000 unfolded ζ zeros.
2. Generate matched-density Monte-Carlo ensembles for GUE, GOE, GSE (each
   yielding ~8000 central eigenvalues, semicircle-unfolded).
3. Sample 32 randomly-spaced 4-tuple offsets (s_1, s_2, s_3) with
   s_1 < s_2 < s_3 in [0.5, 6.0], all spacings >= 0.4 to avoid kernel
   degeneracies.
4. For each tuple, compute empirical R_4 on each of (zeros, GUE-MC, GOE-MC,
   GSE-MC) via window-counting.
5. For each tuple, compute theoretical R_4^{det} = det[K_sine(s_i - s_j)]
   (the sine-kernel DPP / GUE prediction).
6. Compute residuals between (zeros, GUE-MC, GOE-MC, GSE-MC) and the GUE
   theoretical prediction. The model that gives the smallest residual on
   ζ zeros is the best fit.
7. Compute alpha-DPP best-fit by 1-D scan over alpha in [-1.2, -0.8] using
   the sine-kernel as base.
8. Bootstrap z-scores by partitioning each ensemble into K=20 sub-batches.
"""
from __future__ import annotations
from pathlib import Path
import json
import itertools
import time
import numpy as np

ROOT = Path('/apps/aplikacijos/prime-research')
HERE = Path(__file__).parent
ZEROS_PATH = ROOT / 'data' / 'zeta_zeros_8000.txt'
RESULTS_JSON = HERE / 'c6_pfaffian_alpha_dpp_n4_results.json'
RESULTS_MD = HERE / 'c6_pfaffian_alpha_dpp_n4_results.md'

RNG_SEED = 20260501
N_TUPLES = 96
N_BOOTSTRAP_BATCHES = 20
N_ALPHA_BOOTSTRAP = 200

# Window half-width for empirical R_4 counting (in units of mean spacing).
TOL = 0.20


def load_zeros(path):
    return np.array(
        [float(line.strip()) for line in open(path) if line.strip()],
        dtype=np.float64,
    )


def unfold_riemann(gammas):
    return (gammas / (2 * np.pi)) * np.log(gammas / (2 * np.pi * np.e)) + 7 / 8


# ---------------------------------------------------------------------------
# Random matrix ensembles
# ---------------------------------------------------------------------------
def gue_eigs(N, seed):
    rng = np.random.default_rng(seed)
    A = (rng.standard_normal((N, N)) + 1j * rng.standard_normal((N, N))) / np.sqrt(2)
    H = (A + A.conj().T) / np.sqrt(2)
    H /= np.sqrt(2 * N)
    return np.linalg.eigvalsh(H)


def goe_eigs(N, seed):
    rng = np.random.default_rng(seed)
    A = rng.standard_normal((N, N))
    H = (A + A.T) / 2
    H /= np.sqrt(2 * N)
    return np.linalg.eigvalsh(H)


def gse_eigs(N, seed):
    """Quaternion self-dual matrix; 2N complex matrix with sym/antisym structure.
    Practically: generate a 2N x 2N complex matrix that's self-dual (its
    quaternion conjugate equals itself) and Hermitian; eigenvalues come in
    Kramers degenerate pairs — take one of each pair."""
    rng = np.random.default_rng(seed)
    # Build a 2N x 2N matrix M = [[A, B], [-B*, A*]] with A Hermitian and B
    # antisymmetric. Then M is Hermitian self-dual (quaternion-real).
    A = (rng.standard_normal((N, N)) + 1j * rng.standard_normal((N, N))) / np.sqrt(2)
    A = (A + A.conj().T) / 2
    B = (rng.standard_normal((N, N)) + 1j * rng.standard_normal((N, N))) / np.sqrt(2)
    B = (B - B.T) / 2  # antisymmetric
    M = np.zeros((2 * N, 2 * N), dtype=complex)
    M[:N, :N] = A
    M[N:, N:] = A.conj()
    M[:N, N:] = B
    M[N:, :N] = -B.conj()
    # Normalise so that the spectral radius -> 1 on Wigner semicircle on [-1, 1]
    M /= np.sqrt(2 * (2 * N))  # scale matches GUE on the doubled matrix
    evs = np.linalg.eigvalsh(M)
    # Kramers degeneracy: every eigenvalue appears twice (within numerical noise)
    # Take every-other eigenvalue from sorted spectrum.
    return evs[::2]


def semicircle_unfold(evs, N_actual):
    """Semicircle unfolding to unit mean spacing on the central 60% bulk.

    Robust version: detect spectral support [a, b] empirically (rather than
    assuming [-1, 1] which depends on matrix normalisation), apply Wigner-
    semicircle CDF unfolding, then rescale central section to enforce unit
    mean spacing exactly. This corrects for the matrix-normalisation choice.
    """
    evs = np.sort(evs)
    # Empirical support boundary: use 1st/99th percentile padded a bit
    lo_e, hi_e = np.percentile(evs, [0.5, 99.5])
    # Wigner-semicircle on [a, b] with center c = (a+b)/2 and radius R = (b-a)/2
    c = (lo_e + hi_e) / 2
    R = (hi_e - lo_e) / 2
    x = np.clip((evs - c) / R, -0.999999, 0.999999)
    u = (N_actual / np.pi) * (x * np.sqrt(1 - x ** 2) + np.arcsin(x) + np.pi / 2)
    lo = int(0.20 * len(u))
    hi = int(0.80 * len(u))
    u_central = np.sort(u[lo:hi])
    # Enforce unit mean spacing exactly (corrects any residual scale drift)
    sp = np.mean(np.diff(u_central))
    if sp > 0:
        u_central = u_central / sp
    return u_central


def build_ensemble_pool(ens, N_target, batch=2000, seed=RNG_SEED):
    """Build a pool of ensemble eigenvalues batched into multiple matrices.
    Each batch yields ~0.6*batch central eigenvalues."""
    pool = []
    s = seed
    n_central_per_batch = int(0.6 * batch)
    n_batches = (N_target + n_central_per_batch - 1) // n_central_per_batch
    for _ in range(n_batches):
        if ens == 'GUE':
            evs = gue_eigs(batch, s)
            u = semicircle_unfold(evs, N_actual=batch)
        elif ens == 'GOE':
            evs = goe_eigs(batch, s)
            u = semicircle_unfold(evs, N_actual=batch)
        elif ens == 'GSE':
            # For GSE we generate a batch x batch self-dual matrix, eigs come in pairs
            # so we get batch eigenvalues but only batch/2 distinct ones; effective
            # density doubles. Use batch//2 for the "matrix dim" so we get ~0.6*batch
            # central eigenvalues at unit density.
            evs = gse_eigs(batch // 2, s)
            u = semicircle_unfold(evs, N_actual=batch // 2)
        else:
            raise ValueError(ens)
        pool.append(u)
        s += 1
    return pool


# ---------------------------------------------------------------------------
# Sine kernel + alpha-determinant
# ---------------------------------------------------------------------------
def K_sine(t):
    if abs(t) < 1e-14:
        return 1.0
    return np.sin(np.pi * t) / (np.pi * t)


def K_sine_matrix(s_offsets):
    n = len(s_offsets)
    K = np.empty((n, n))
    for i in range(n):
        for j in range(n):
            K[i, j] = K_sine(s_offsets[i] - s_offsets[j])
    return K


def num_cycles(perm):
    n = len(perm)
    seen = [False] * n
    c = 0
    for i in range(n):
        if seen[i]:
            continue
        c += 1
        j = i
        while not seen[j]:
            seen[j] = True
            j = perm[j]
    return c


def alpha_det(A, alpha):
    """alpha-determinant: sum over permutations sigma of
    alpha^{n - cycles(sigma)} * prod_i A[i, sigma(i)].
    For alpha = -1: standard determinant.
    For alpha = +1: permanent.
    For alpha = 0: product of diagonal entries.
    """
    n = A.shape[0]
    total = 0.0
    for sigma in itertools.permutations(range(n)):
        c = num_cycles(sigma)
        prod = 1.0
        for i in range(n):
            prod *= A[i, sigma[i]]
        total += (alpha ** (n - c)) * prod
    return total


# ---------------------------------------------------------------------------
# Empirical R_4 estimator
# ---------------------------------------------------------------------------
def empirical_R4(unfolded, s_offsets, tol):
    """Empirical density R_4(0, s_1, s_2, s_3) where s_offsets = (s_1, s_2, s_3).

    For each reference zero z_i, count the number of triples (j_1, j_2, j_3)
    with z_{j_k} - z_i in [s_k - tol/2, s_k + tol/2] for k=1,2,3.

    R_4 ~ (count) / (n_ref * tol^3).

    n_ref = number of i such that the window [z_i, z_i + s_3 + tol/2] is within
    the ensemble.
    """
    s_arr = np.array(sorted(s_offsets), dtype=float)
    half = tol / 2
    last_end = unfolded[-1]
    counts = 0
    n_ref = 0
    N = len(unfolded)
    for i in range(N):
        z_i = unfolded[i]
        if z_i + s_arr[-1] + half > last_end:
            break
        n_ref += 1
        slices = []
        for sk in s_arr:
            lo = z_i + sk - half
            hi = z_i + sk + half
            j_lo = np.searchsorted(unfolded, lo, side='left')
            j_hi = np.searchsorted(unfolded, hi, side='left')
            slices.append(j_hi - j_lo)
        if all(c > 0 for c in slices):
            counts += int(np.prod(slices))
    if n_ref == 0:
        return float('nan'), 0
    density = counts / (n_ref * tol ** 3)
    return float(density), n_ref


def empirical_R4_pool(pool, s_offsets, tol):
    """Aggregate R_4 across pool of unfolded ensembles, weighted by n_ref."""
    weight_sum = 0
    val_sum = 0.0
    for u in pool:
        d, nref = empirical_R4(u, s_offsets, tol)
        if not np.isnan(d):
            val_sum += d * nref
            weight_sum += nref
    if weight_sum == 0:
        return float('nan')
    return val_sum / weight_sum


def empirical_R4_pool_per_batch(pool, s_offsets, tol):
    """Return list of per-batch R_4 measurements (for bootstrap)."""
    out = []
    for u in pool:
        d, nref = empirical_R4(u, s_offsets, tol)
        if not np.isnan(d):
            out.append(d)
    return out


# ---------------------------------------------------------------------------
# Tuple sampling
# ---------------------------------------------------------------------------
def sample_tuples(n_tuples, seed, s_min=0.5, s_max=6.0, gap_min=0.4):
    """Sample n_tuples 3-tuples (s_1, s_2, s_3) with all pairwise differences
    >= gap_min and 0 < s_1 < s_2 < s_3 in [s_min, s_max]. Also include s_1
    > gap_min (gap from origin).
    """
    rng = np.random.default_rng(seed)
    out = []
    while len(out) < n_tuples:
        ss = sorted(rng.uniform(s_min, s_max, size=3))
        if ss[0] >= gap_min and ss[1] - ss[0] >= gap_min and ss[2] - ss[1] >= gap_min:
            out.append(tuple(round(x, 3) for x in ss))
    return out


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    t0 = time.time()
    print("Loading zeros...", flush=True)
    gammas = load_zeros(ZEROS_PATH)
    u_zeta = unfold_riemann(gammas)
    print(f"  N={len(gammas)}, mean spacing = {np.diff(u_zeta).mean():.6f}",
          flush=True)

    # Use ensemble pool 4x larger than ζ for tighter control noise floors
    print("\nGenerating GUE pool...", flush=True)
    t1 = time.time()
    gue_pool = build_ensemble_pool('GUE', N_target=32000, batch=2000, seed=RNG_SEED)
    print(f"  {len(gue_pool)} batches, {sum(len(u) for u in gue_pool)} eigs "
          f"({time.time() - t1:.1f}s)", flush=True)

    print("\nGenerating GOE pool...", flush=True)
    t1 = time.time()
    goe_pool = build_ensemble_pool('GOE', N_target=32000, batch=2000,
                                    seed=RNG_SEED + 1000)
    print(f"  {len(goe_pool)} batches, {sum(len(u) for u in goe_pool)} eigs "
          f"({time.time() - t1:.1f}s)", flush=True)

    print("\nGenerating GSE pool...", flush=True)
    t1 = time.time()
    gse_pool = build_ensemble_pool('GSE', N_target=32000, batch=2000,
                                    seed=RNG_SEED + 2000)
    print(f"  {len(gse_pool)} batches, {sum(len(u) for u in gse_pool)} eigs "
          f"({time.time() - t1:.1f}s)", flush=True)

    # Sample 4-tuple offsets
    tuples = sample_tuples(N_TUPLES, seed=RNG_SEED + 5000)
    print(f"\nSampled {len(tuples)} 4-tuple offsets", flush=True)

    # For each tuple compute:
    # - R_4^GUE_theory (sine-kernel determinant)
    # - empirical R_4 on (zeros, GUE-MC, GOE-MC, GSE-MC)
    # - alpha-DPP fit residual

    rows = []
    for k, ss in enumerate(tuples):
        s_offsets = (0.0,) + ss  # (0, s_1, s_2, s_3)

        # Theoretical sine-kernel determinant
        K = K_sine_matrix(s_offsets)
        R4_det_GUE = float(np.linalg.det(K))

        # Empirical R_4
        R4_zeta, nref_zeta = empirical_R4(u_zeta, s_offsets[1:], TOL)
        R4_GUE_emp = empirical_R4_pool(gue_pool, s_offsets[1:], TOL)
        R4_GOE_emp = empirical_R4_pool(goe_pool, s_offsets[1:], TOL)
        R4_GSE_emp = empirical_R4_pool(gse_pool, s_offsets[1:], TOL)

        # Per-batch for bootstrap
        R4_zeta_batches = empirical_R4_pool_per_batch(
            np.array_split(u_zeta, N_BOOTSTRAP_BATCHES), s_offsets[1:], TOL,
        )

        rows.append({
            'tuple': list(ss),
            'R4_det_GUE_theory': R4_det_GUE,
            'R4_zeta_emp': R4_zeta,
            'nref_zeta': nref_zeta,
            'R4_GUE_emp': R4_GUE_emp,
            'R4_GOE_emp': R4_GOE_emp,
            'R4_GSE_emp': R4_GSE_emp,
            'R4_zeta_batches': [float(x) for x in R4_zeta_batches],
        })
        print(f"  tuple {k+1}/{len(tuples)} ss={ss}: "
              f"theory={R4_det_GUE:.4f}  zeta={R4_zeta:.4f}  "
              f"GUE_emp={R4_GUE_emp:.4f}  GOE_emp={R4_GOE_emp:.4f}  "
              f"GSE_emp={R4_GSE_emp:.4f}", flush=True)

    # Aggregate analysis
    R4_det = np.array([r['R4_det_GUE_theory'] for r in rows])
    R4_zeta = np.array([r['R4_zeta_emp'] for r in rows])
    R4_GUE = np.array([r['R4_GUE_emp'] for r in rows])
    R4_GOE = np.array([r['R4_GOE_emp'] for r in rows])
    R4_GSE = np.array([r['R4_GSE_emp'] for r in rows])

    # Per-tuple residual against THEORETICAL sine-kernel determinant
    res_zeta_vs_theory = R4_zeta - R4_det
    res_GUE_vs_theory = R4_GUE - R4_det
    res_GOE_vs_theory = R4_GOE - R4_det
    res_GSE_vs_theory = R4_GSE - R4_det

    # Direct empirical comparison: zeta vs each ensemble
    res_zeta_vs_GUE = R4_zeta - R4_GUE
    res_zeta_vs_GOE = R4_zeta - R4_GOE
    res_zeta_vs_GSE = R4_zeta - R4_GSE

    # L2 norms of residuals (lower = better fit)
    L2_zeta_vs_theory = float(np.sqrt(np.mean(res_zeta_vs_theory ** 2)))
    L2_GUE_vs_theory = float(np.sqrt(np.mean(res_GUE_vs_theory ** 2)))
    L2_GOE_vs_theory = float(np.sqrt(np.mean(res_GOE_vs_theory ** 2)))
    L2_GSE_vs_theory = float(np.sqrt(np.mean(res_GSE_vs_theory ** 2)))

    L2_zeta_vs_GUE = float(np.sqrt(np.mean(res_zeta_vs_GUE ** 2)))
    L2_zeta_vs_GOE = float(np.sqrt(np.mean(res_zeta_vs_GOE ** 2)))
    L2_zeta_vs_GSE = float(np.sqrt(np.mean(res_zeta_vs_GSE ** 2)))

    print(f"\n--- RMS residuals vs theoretical sine-kernel det[K] ---", flush=True)
    print(f"  zeta vs det:  {L2_zeta_vs_theory:.5f}", flush=True)
    print(f"  GUE  vs det:  {L2_GUE_vs_theory:.5f}", flush=True)
    print(f"  GOE  vs det:  {L2_GOE_vs_theory:.5f}", flush=True)
    print(f"  GSE  vs det:  {L2_GSE_vs_theory:.5f}", flush=True)

    print(f"\n--- RMS residuals zeta vs each empirical ensemble ---", flush=True)
    print(f"  zeta vs GUE_emp: {L2_zeta_vs_GUE:.5f}", flush=True)
    print(f"  zeta vs GOE_emp: {L2_zeta_vs_GOE:.5f}", flush=True)
    print(f"  zeta vs GSE_emp: {L2_zeta_vs_GSE:.5f}", flush=True)

    # alpha-DPP best fit + bootstrap CI on alpha*
    print(f"\n--- alpha-DPP best fit on zeta R_4 ---", flush=True)
    # Fine grid to get tight alpha* point estimate
    alphas = np.linspace(-1.5, -0.5, 41)
    # Pre-compute alpha-det predictions for each tuple at every alpha
    preds_grid = np.empty((len(rows), len(alphas)))
    for i, r in enumerate(rows):
        s_offsets = (0.0,) + tuple(r['tuple'])
        K = K_sine_matrix(s_offsets)
        for j, alpha in enumerate(alphas):
            preds_grid[i, j] = alpha_det(K, alpha)
    losses = np.sqrt(((R4_zeta[:, None] - preds_grid) ** 2).mean(axis=0))
    for j in range(len(alphas)):
        print(f"  alpha={alphas[j]:+.3f}: L2 = {losses[j]:.5f}", flush=True)

    best_idx = int(np.argmin(losses))
    alpha_star = float(alphas[best_idx])
    L2_star = float(losses[best_idx])
    L2_at_minus1 = float(losses[np.argmin(np.abs(alphas - (-1.0)))])
    print(f"\n  alpha* (best fit) = {alpha_star:.3f}, L2* = {L2_star:.5f}", flush=True)
    print(f"  L2 at alpha=-1 (GUE/DPP) = {L2_at_minus1:.5f}", flush=True)
    print(f"  delta L2 (alpha* vs -1) = {L2_at_minus1 - L2_star:.5f}", flush=True)

    # Bootstrap CI on alpha* by resampling tuples (32-of-32 with replacement)
    print(f"\n--- bootstrap alpha* CI ({N_ALPHA_BOOTSTRAP} resamples) ---",
          flush=True)
    rng_boot = np.random.default_rng(RNG_SEED + 9999)
    alpha_boot = []
    for b in range(N_ALPHA_BOOTSTRAP):
        idx = rng_boot.choice(len(rows), size=len(rows), replace=True)
        l_b = np.sqrt(((R4_zeta[idx, None] - preds_grid[idx, :]) ** 2).mean(axis=0))
        alpha_boot.append(float(alphas[int(np.argmin(l_b))]))
    alpha_boot = np.array(alpha_boot)
    alpha_p2 = float(np.percentile(alpha_boot, 2.5))
    alpha_p97 = float(np.percentile(alpha_boot, 97.5))
    alpha_p16 = float(np.percentile(alpha_boot, 16))
    alpha_p84 = float(np.percentile(alpha_boot, 84))
    alpha_mean = float(alpha_boot.mean())
    print(f"  bootstrap alpha* mean = {alpha_mean:.3f}", flush=True)
    print(f"  bootstrap alpha* 95% CI = [{alpha_p2:.3f}, {alpha_p97:.3f}]",
          flush=True)
    print(f"  bootstrap alpha* 68% CI = [{alpha_p16:.3f}, {alpha_p84:.3f}]",
          flush=True)
    rejects_minus1 = (alpha_p2 > -1.0) or (alpha_p97 < -1.0)
    print(f"  alpha=-1 rejected at 95% CI: {rejects_minus1}", flush=True)

    # Per-tuple R_4 measurement on each individual ensemble batch — gives an
    # empirical noise distribution PER TUPLE for ensembles of ~2000 eigenvalues
    # (each batch). The matched-zero ensemble has 8000 zeros so its noise is
    # ~sqrt(N_batch / N_zeta) = sqrt(2000/8000) = 0.5 the noise of one batch.
    print(f"\n--- per-tuple z-score discrimination ---", flush=True)
    # For each tuple compute std across GUE batches; this is the matched
    # noise model per tuple (since GUE-MC batches have ~2000 eigenvalues each
    # and ζ has 8000, scale by 1/2 to convert to ζ-equivalent SE)
    sq_factor = float(np.sqrt(2000.0 / 8000.0))  # 0.5
    per_tuple_GUE_std = []
    per_tuple_GOE_std = []
    per_tuple_GSE_std = []
    for k, ss in enumerate(tuples):
        s_offsets = (0.0,) + ss
        per_batch_GUE = [empirical_R4(u, s_offsets[1:], TOL)[0] for u in gue_pool]
        per_batch_GOE = [empirical_R4(u, s_offsets[1:], TOL)[0] for u in goe_pool]
        per_batch_GSE = [empirical_R4(u, s_offsets[1:], TOL)[0] for u in gse_pool]
        # Use the std of per-batch as the per-batch noise estimate.
        per_tuple_GUE_std.append(float(np.std([x for x in per_batch_GUE if not np.isnan(x)])))
        per_tuple_GOE_std.append(float(np.std([x for x in per_batch_GOE if not np.isnan(x)])))
        per_tuple_GSE_std.append(float(np.std([x for x in per_batch_GSE if not np.isnan(x)])))
    per_tuple_GUE_std = np.array(per_tuple_GUE_std)
    per_tuple_GOE_std = np.array(per_tuple_GOE_std)
    per_tuple_GSE_std = np.array(per_tuple_GSE_std)
    # zeta SE per tuple ≈ GUE batch std × sq_factor
    zeta_SE_per_tuple = per_tuple_GUE_std * sq_factor
    # Mask out tuples where SE ≈ 0
    mask = zeta_SE_per_tuple > 1e-4
    print(f"  zeta SE per tuple (matched-noise model): "
          f"mean = {zeta_SE_per_tuple[mask].mean():.4f}", flush=True)

    # Per-tuple z-score for each model: (R4_emp - R4_pred) / SE
    z_per_tuple_GUE_pred = (R4_zeta - R4_GUE) / np.maximum(zeta_SE_per_tuple, 1e-6)
    z_per_tuple_GOE_pred = (R4_zeta - R4_GOE) / np.maximum(zeta_SE_per_tuple, 1e-6)
    z_per_tuple_GSE_pred = (R4_zeta - R4_GSE) / np.maximum(zeta_SE_per_tuple, 1e-6)
    z_per_tuple_det = (R4_zeta - R4_det) / np.maximum(zeta_SE_per_tuple, 1e-6)

    # Aggregate: chi-square sum and rms z
    chi2_GUE_pred = float(np.nansum(z_per_tuple_GUE_pred[mask] ** 2))
    chi2_GOE_pred = float(np.nansum(z_per_tuple_GOE_pred[mask] ** 2))
    chi2_GSE_pred = float(np.nansum(z_per_tuple_GSE_pred[mask] ** 2))
    chi2_det = float(np.nansum(z_per_tuple_det[mask] ** 2))
    dof = int(mask.sum())
    print(f"  chi^2 per tuple summed (dof={dof}):", flush=True)
    print(f"    zeta vs GUE-MC: {chi2_GUE_pred:.1f}  (chi^2/dof = {chi2_GUE_pred/dof:.2f})", flush=True)
    print(f"    zeta vs GOE-MC: {chi2_GOE_pred:.1f}  (chi^2/dof = {chi2_GOE_pred/dof:.2f})", flush=True)
    print(f"    zeta vs GSE-MC: {chi2_GSE_pred:.1f}  (chi^2/dof = {chi2_GSE_pred/dof:.2f})", flush=True)
    print(f"    zeta vs det:    {chi2_det:.1f}  (chi^2/dof = {chi2_det/dof:.2f})", flush=True)
    # If chi2/dof ~ 1: model matches data within noise
    # If chi2/dof >> 1: model rejected

    # Pairwise z-scores from L2 RMS / overall noise
    rms_zeta_SE = float(np.sqrt(np.mean(zeta_SE_per_tuple[mask] ** 2)))

    def z_score(rms_residual):
        if rms_zeta_SE < 1e-12:
            return float('nan')
        return rms_residual / rms_zeta_SE

    z_zeta_vs_GUE_emp = z_score(L2_zeta_vs_GUE)
    z_zeta_vs_GOE_emp = z_score(L2_zeta_vs_GOE)
    z_zeta_vs_GSE_emp = z_score(L2_zeta_vs_GSE)
    z_zeta_vs_theory = z_score(L2_zeta_vs_theory)
    z_GUE_vs_theory = z_score(L2_GUE_vs_theory)
    z_GOE_vs_theory = z_score(L2_GOE_vs_theory)
    z_GSE_vs_theory = z_score(L2_GSE_vs_theory)

    print(f"\n--- L2 RMS / RMS-zeta-SE z-scores ---", flush=True)
    print(f"  zeta SE RMS = {rms_zeta_SE:.5f}", flush=True)
    print(f"  z(zeta vs theory)   = {z_zeta_vs_theory:.2f}", flush=True)
    print(f"  z(GUE-MC vs theory) = {z_GUE_vs_theory:.2f}", flush=True)
    print(f"  z(GOE-MC vs theory) = {z_GOE_vs_theory:.2f}", flush=True)
    print(f"  z(GSE-MC vs theory) = {z_GSE_vs_theory:.2f}", flush=True)
    print(f"  z(zeta vs GUE_emp)  = {z_zeta_vs_GUE_emp:.2f}", flush=True)
    print(f"  z(zeta vs GOE_emp)  = {z_zeta_vs_GOE_emp:.2f}", flush=True)
    print(f"  z(zeta vs GSE_emp)  = {z_zeta_vs_GSE_emp:.2f}", flush=True)

    # Discriminating power: how cleanly does the test reject GOE/GSE?
    discrimination_GOE = (L2_GOE_vs_theory - L2_GUE_vs_theory) / rms_zeta_SE
    discrimination_GSE = (L2_GSE_vs_theory - L2_GUE_vs_theory) / rms_zeta_SE
    print(f"\n--- discrimination of Pfaffian (GOE/GSE) vs DPP (GUE) ---",
          flush=True)
    print(f"  Δ(L2_GOE - L2_GUE) / SE = {discrimination_GOE:.2f}σ",
          flush=True)
    print(f"  Δ(L2_GSE - L2_GUE) / SE = {discrimination_GSE:.2f}σ",
          flush=True)

    out = {
        'meta': {
            'N_zeta': len(gammas),
            'gamma_min': float(gammas[0]),
            'gamma_max': float(gammas[-1]),
            'mean_spacing': float(np.diff(u_zeta).mean()),
            'GUE_pool_size': int(sum(len(u) for u in gue_pool)),
            'GOE_pool_size': int(sum(len(u) for u in goe_pool)),
            'GSE_pool_size': int(sum(len(u) for u in gse_pool)),
            'n_tuples': N_TUPLES,
            'tol': TOL,
            'rng_seed': RNG_SEED,
            'runtime_s': float(time.time() - t0),
        },
        'rows': rows,
        'L2_residuals_vs_theory': {
            'zeta': L2_zeta_vs_theory,
            'GUE': L2_GUE_vs_theory,
            'GOE': L2_GOE_vs_theory,
            'GSE': L2_GSE_vs_theory,
        },
        'L2_residuals_zeta_vs_ensemble': {
            'GUE_emp': L2_zeta_vs_GUE,
            'GOE_emp': L2_zeta_vs_GOE,
            'GSE_emp': L2_zeta_vs_GSE,
        },
        'alpha_dpp_scan': {
            'alphas': alphas.tolist(),
            'L2_losses': losses.tolist(),
            'alpha_star': alpha_star,
            'L2_star': L2_star,
            'L2_at_alpha_minus1': L2_at_minus1,
            'alpha_boot_mean': alpha_mean,
            'alpha_boot_p2': alpha_p2,
            'alpha_boot_p97': alpha_p97,
            'alpha_boot_p16': alpha_p16,
            'alpha_boot_p84': alpha_p84,
            'alpha_minus1_rejected_at_95pct': bool(rejects_minus1),
        },
        'bootstrap': {
            'rms_zeta_SE': rms_zeta_SE,
            'z_zeta_vs_theory': z_zeta_vs_theory,
            'z_GUE_vs_theory': z_GUE_vs_theory,
            'z_GOE_vs_theory': z_GOE_vs_theory,
            'z_GSE_vs_theory': z_GSE_vs_theory,
            'z_zeta_vs_GUE_emp': z_zeta_vs_GUE_emp,
            'z_zeta_vs_GOE_emp': z_zeta_vs_GOE_emp,
            'z_zeta_vs_GSE_emp': z_zeta_vs_GSE_emp,
            'discrimination_GOE_sigma': float(discrimination_GOE),
            'discrimination_GSE_sigma': float(discrimination_GSE),
            'chi2_GUE_pred': chi2_GUE_pred,
            'chi2_GOE_pred': chi2_GOE_pred,
            'chi2_GSE_pred': chi2_GSE_pred,
            'chi2_det': chi2_det,
            'dof': dof,
        },
    }
    RESULTS_JSON.write_text(json.dumps(out, indent=2))
    print(f"\nWROTE {RESULTS_JSON}", flush=True)
    print(f"TOTAL runtime: {time.time() - t0:.1f}s", flush=True)


if __name__ == '__main__':
    main()

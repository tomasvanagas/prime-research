"""
C6 BIAS CONTROL — critical falsifier for the alpha* = -1.10 finding.

The main C6 experiment found alpha* = -1.10 with bootstrap CI [-1.15, -1.05]
on the 8000-zero ensemble. The CI excludes alpha = -1, which would indicate
ζ-zero R_4 deviates from sine-kernel DPP toward an alpha-DPP with alpha < -1.

CRITICAL CONTROL: does the same procedure applied to a MATCHED-SIZE
8000-eigenvalue GUE Monte-Carlo ensemble also yield alpha* < -1?
- If yes: the observed shift is finite-sample / finite-tol estimator bias,
  NOT a real ζ-zero structural fact. alpha* = -1.10 is a measurement
  artifact and the C6 attack closes B-grade structural negative.
- If no: the shift on ζ is genuine and exceeds what GUE-MC produces.
  alpha* = -1.10 is structural content. C6 closes A-grade.

PROCEDURE: generate K independent 8000-eigenvalue GUE-MC ensembles, run
the alpha-DPP fit on each (using the SAME 96 tuples as the main run),
report the distribution of alpha*.
"""
from __future__ import annotations
from pathlib import Path
import json
import time
import sys
import numpy as np

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
from c6_pfaffian_alpha_dpp_n4 import (  # noqa: E402
    gue_eigs, semicircle_unfold,
    K_sine_matrix, alpha_det,
    empirical_R4, sample_tuples,
)

ROOT = Path('/apps/aplikacijos/prime-research')
RESULTS_JSON = HERE / 'c6_alpha_bias_control_results.json'

RNG_SEED = 20260501
N_TUPLES = 96
TOL = 0.20

# K independent 8000-eigenvalue GUE ensembles. Each requires ~600 central
# eigenvalues from a 1000x1000 GUE matrix to keep cost low. To get 8000
# central eigenvalues, use ~14 batches of 1000 each, but this is "one
# ensemble" — concatenated as if we had 8000 ζ zeros from one realization.
N_TRIALS = 30
N_EIGS_PER_TRIAL = 8000


def build_one_gue_ensemble(N_target, seed, batch=2000):
    """Build a single GUE 'ensemble' of N_target eigenvalues by stacking
    contiguous unfolded blocks. Important: each block must be unfolded
    independently then offset to be contiguous (mimicking what we'd get
    from a single 8000-zero ζ-stretch)."""
    pool = []
    s = seed
    n_central_per_batch = int(0.6 * batch)
    n_batches = (N_target + n_central_per_batch - 1) // n_central_per_batch
    offset = 0.0
    for _ in range(n_batches):
        evs = gue_eigs(batch, s)
        u = semicircle_unfold(evs, N_actual=batch)
        u = u - u[0] + offset + 1.0  # shift contiguous
        pool.append(u)
        offset = u[-1]
        s += 1
    full = np.concatenate(pool)
    return full[:N_target]


def fit_alpha_star(R4_emp, preds_grid, alphas):
    """Find alpha that minimises L2(R4_emp - alpha-det predictions)."""
    losses = np.sqrt(((R4_emp[:, None] - preds_grid) ** 2).mean(axis=0))
    best_idx = int(np.argmin(losses))
    return float(alphas[best_idx]), float(losses[best_idx]), losses.tolist()


def main():
    t0 = time.time()
    print(f"Generating {N_TRIALS} independent 8000-eigenvalue GUE ensembles "
          f"and computing per-trial alpha*.", flush=True)
    print(f"This is the BIAS CONTROL for the main C6 experiment's "
          f"alpha* = -1.10 finding.", flush=True)

    # Same tuples as main run
    tuples = sample_tuples(N_TUPLES, seed=RNG_SEED + 5000)
    print(f"Using same {len(tuples)} 4-tuple offsets as main run.", flush=True)

    alphas = np.linspace(-1.5, -0.5, 41)
    # Pre-compute alpha-det predictions for each tuple at every alpha
    print("Pre-computing alpha-det theoretical predictions...", flush=True)
    preds_grid = np.empty((len(tuples), len(alphas)))
    for i, ss in enumerate(tuples):
        s_offsets = (0.0,) + ss
        K = K_sine_matrix(s_offsets)
        for j, alpha in enumerate(alphas):
            preds_grid[i, j] = alpha_det(K, alpha)

    alpha_star_results = []
    L2_star_results = []
    for trial in range(N_TRIALS):
        seed_t = RNG_SEED + 100000 + trial * 1000
        u_one = build_one_gue_ensemble(N_EIGS_PER_TRIAL, seed=seed_t)
        # Compute empirical R_4 for each tuple
        R4_emp = np.empty(len(tuples))
        for i, ss in enumerate(tuples):
            d, _ = empirical_R4(u_one, ss, TOL)
            R4_emp[i] = d
        a_star, L2_star, losses = fit_alpha_star(R4_emp, preds_grid, alphas)
        alpha_star_results.append(a_star)
        L2_star_results.append(L2_star)
        print(f"  trial {trial+1}/{N_TRIALS}: alpha* = {a_star:+.3f}, "
              f"L2* = {L2_star:.4f}", flush=True)

    alpha_star_results = np.array(alpha_star_results)
    L2_star_results = np.array(L2_star_results)

    print(f"\n--- Distribution of alpha* over {N_TRIALS} GUE-MC trials ---",
          flush=True)
    print(f"  mean   = {alpha_star_results.mean():.3f}", flush=True)
    print(f"  median = {np.median(alpha_star_results):.3f}", flush=True)
    print(f"  std    = {alpha_star_results.std():.3f}", flush=True)
    print(f"  min    = {alpha_star_results.min():.3f}", flush=True)
    print(f"  max    = {alpha_star_results.max():.3f}", flush=True)
    p2, p25, p75, p97 = np.percentile(alpha_star_results, [2.5, 25, 75, 97.5])
    print(f"  95% CI = [{p2:.3f}, {p97:.3f}]", flush=True)
    print(f"  IQR    = [{p25:.3f}, {p75:.3f}]", flush=True)

    # Compare to ζ result
    zeta_alpha_star = -1.10
    print(f"\n--- COMPARISON to zeta finding ---", flush=True)
    print(f"  zeta alpha* = {zeta_alpha_star:.3f}", flush=True)
    print(f"  GUE-MC alpha* mean = {alpha_star_results.mean():.3f}", flush=True)
    print(f"  delta (zeta - GUE-MC mean) = "
          f"{zeta_alpha_star - alpha_star_results.mean():+.3f}", flush=True)
    z_zeta = (zeta_alpha_star - alpha_star_results.mean()) / max(alpha_star_results.std(), 1e-6)
    print(f"  zeta alpha* z-score relative to GUE-MC distribution: {z_zeta:.2f}σ",
          flush=True)
    n_below = int(np.sum(alpha_star_results <= zeta_alpha_star))
    print(f"  GUE-MC trials with alpha* <= -1.10: {n_below}/{N_TRIALS} "
          f"= {100 * n_below / N_TRIALS:.1f}%", flush=True)

    out = {
        'meta': {
            'n_trials': N_TRIALS,
            'n_eigs_per_trial': N_EIGS_PER_TRIAL,
            'n_tuples': len(tuples),
            'tol': TOL,
            'rng_seed': RNG_SEED,
            'runtime_s': float(time.time() - t0),
        },
        'alpha_star_distribution': {
            'values': alpha_star_results.tolist(),
            'L2_values': L2_star_results.tolist(),
            'mean': float(alpha_star_results.mean()),
            'median': float(np.median(alpha_star_results)),
            'std': float(alpha_star_results.std()),
            'min': float(alpha_star_results.min()),
            'max': float(alpha_star_results.max()),
            'p2_5': float(p2),
            'p97_5': float(p97),
            'p25': float(p25),
            'p75': float(p75),
        },
        'comparison_to_zeta': {
            'zeta_alpha_star': zeta_alpha_star,
            'gue_mc_mean': float(alpha_star_results.mean()),
            'z_score_zeta_vs_gue_distribution': float(z_zeta),
            'gue_trials_at_or_below_zeta_alpha_pct':
                float(100 * n_below / N_TRIALS),
        },
    }
    RESULTS_JSON.write_text(json.dumps(out, indent=2))
    print(f"\nWROTE {RESULTS_JSON}", flush=True)
    print(f"TOTAL runtime: {time.time() - t0:.1f}s", flush=True)


if __name__ == '__main__':
    main()

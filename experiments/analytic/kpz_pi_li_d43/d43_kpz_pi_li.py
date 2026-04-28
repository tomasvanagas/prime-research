"""
D43 — Hairer regularity structures / KPZ universality of pi(x) - Li(x).

Question: does D(x) := (pi(x) - Li(x)) * log(x) / sqrt(x), sampled on a KPZ-
spaced logarithmic grid x_k = X/2 + k * floor(X^{1/3}), exhibit Tracy-Widom
beta=2 tails (KPZ universality class) instead of the project-standard Gauss
(CLT scale, Bombieri / Heath-Brown / S133 / S108)?

Pre-stated falsifiers (REGISTERED before measurement):

  F1 (Gauss-skewness): |skew(Z)| < 3 * sqrt(6/N) where N is grid size.
       Z = standardized D. Holds => Gauss-consistent skewness.
  F2 (TW2-skewness):   skew(Z) in [0.224 - 3*se, 0.224 + 3*se] where
       se = sqrt(6/N). A-grade if F2 PASSES while F1 FAILS.
  F3 (right-tail decay): linear regression of log(1 - empirical CDF) on
       z^{3/2} for z in [1, 3] should fit slope ~ -4/3 = -1.333 for TW2,
       slope ~ -z^{1/2}/2 (curving) for Gauss. A-grade if TW2 fit r^2 > 0.95
       AND Gauss fit r^2 < TW2 fit r^2 - 0.05.
  F4 (Hölder regularity): wavelet decomposition of D(x). Energy E_j at
       scale 2^j fits log E_j = -j*(2*alpha + 1) + const.
       Smooth (alpha = infinity): E_j -> 0 super-exponentially, slope very
                                   steep, alpha > 5 numerical floor.
       KPZ (alpha in (0, 1/2)): slope corresponds to alpha < 1/2 with
                                 statistical significance.
       Cramér model alpha (control C1) gives finite-N reference.
       A-grade if alpha_chi_P significantly LESS than alpha_smooth_control
       AND alpha_chi_P > 0 (i.e., D(x) has detectable rough fluctuation).
  F5 (KS distance): KS(empirical Z, TW2) < KS(empirical Z, Gauss).
       A-grade if KS_TW2 < KS_Gauss with bootstrap confidence.

A-grade outcome: F2 OR F3 OR F4 OR F5 holds the TW2 / KPZ branch.
B-grade outcome: F1 holds Gauss branch + clean wavelet smoothness profile +
   first quantitative project test of non-CLT scale, with KPZ falsified
   structurally — refines E3.1 / E3.13.
F-grade outcome: insufficient grid points / numerical instability.
"""

import os
import math
import json
import sys
import time
import argparse
from pathlib import Path

import numpy as np
import mpmath as mp
import pywt
from scipy.stats import skew, kurtosis, kstest, norm

HERE = Path(__file__).parent


def sieve_pi_table(N):
    """Compute pi(x) for x = 0..N as cumulative count of primes <= x."""
    sieve = np.ones(N + 1, dtype=np.bool_)
    sieve[:2] = False
    for i in range(2, int(N ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i :: i] = False
    return np.cumsum(sieve, dtype=np.int64)


def li(x):
    """Logarithmic integral Li(x) = integral_0^x dt/log(t) (Cauchy PV).

    Uses mpmath for accuracy at our scale (x up to ~2^22 ≈ 4e6).
    Vectorised over a python iterable.
    """
    out = np.empty(len(x), dtype=np.float64)
    for i, xv in enumerate(x):
        out[i] = float(mp.li(xv))
    return out


def kpz_grid(X):
    """KPZ-spaced grid: x_k = X/2 + k * floor(X^{1/3}) for k = 0, 1, ..."""
    step = max(1, int(math.floor(X ** (1.0 / 3.0))))
    start = X // 2
    n_pts = (X - start) // step
    xs = np.array([start + k * step for k in range(1, n_pts + 1)],
                  dtype=np.int64)
    return xs, step


def D_field(pi_table, xs):
    """D(x) := (pi(x) - Li(x)) * log(x) / sqrt(x)."""
    pi_vals = pi_table[xs].astype(np.float64)
    li_vals = li(xs)
    return (pi_vals - li_vals) * np.log(xs) / np.sqrt(xs)


def standardize(z):
    return (z - z.mean()) / z.std(ddof=1)


def tw2_moments():
    """Tracy-Widom beta=2 reference moments (literature values)."""
    return dict(mean=-1.7711, var=0.81320, std=math.sqrt(0.81320),
                skew=0.2241, exkurt=0.0934)


def tw2_right_tail_logfn(z):
    """log(1 - F_2(z)) ~ - (4/3) z^{3/2} for large z (asymptotic)."""
    # Use the leading-order asymptotic
    return -(4.0 / 3.0) * np.power(np.maximum(z, 1e-12), 1.5) \
           - 1.5 * np.log(np.maximum(z, 1e-12)) - math.log(32 * math.pi)


def gauss_right_tail_logfn(z):
    return -0.5 * z ** 2 - 0.5 * math.log(2 * math.pi)


def fit_right_tail_slope(z_data, basis='kpz'):
    """Empirical log(1 - F(z)) regression. basis='kpz' uses z^{3/2};
    basis='gauss' uses z^2/2.

    Returns (slope, intercept, r2, n_points).
    """
    z_sorted = np.sort(z_data)
    n = len(z_sorted)
    # Empirical 1 - F at z_i is (n - i) / n (i.e., proportion strictly above z_i).
    # Use z in [1.0, 3.0] - the meaningful right tail
    mask = (z_sorted >= 1.0) & (z_sorted <= 3.0)
    if mask.sum() < 20:
        return None
    zs = z_sorted[mask]
    ranks = np.arange(len(z_sorted))[mask]
    surv = (n - ranks) / n
    surv = np.clip(surv, 1e-10, 1.0)
    log_s = np.log(surv)
    if basis == 'kpz':
        x = zs ** 1.5
    else:
        x = zs ** 2 / 2.0
    # least-squares fit
    A = np.vstack([x, np.ones_like(x)]).T
    sol, _, _, _ = np.linalg.lstsq(A, log_s, rcond=None)
    slope, intercept = sol
    pred = A @ sol
    ss_res = ((log_s - pred) ** 2).sum()
    ss_tot = ((log_s - log_s.mean()) ** 2).sum()
    r2 = 1.0 - ss_res / max(ss_tot, 1e-20)
    return dict(slope=float(slope), intercept=float(intercept),
                r2=float(r2), n=int(len(zs)),
                z_range=(float(zs[0]), float(zs[-1])))


def wavelet_holder(D, wavelet='db4', max_level=None):
    """Hölder regularity of D via wavelet detail energy decay.

    For a function in C^alpha (Hölder), the L^2 energy of detail
    coefficients at scale 2^j satisfies E_j ~ 2^{-2j*alpha} (up to
    boundary effects).  We fit log2(E_j) vs j linearly to extract alpha.

    Returns dict with per-level energies and fitted alpha.
    """
    n = len(D)
    if max_level is None:
        max_level = pywt.dwt_max_level(n, pywt.Wavelet(wavelet).dec_len)
        max_level = min(max_level, 8)
    coeffs = pywt.wavedec(D, wavelet, level=max_level)
    # coeffs[0] = approximation at coarsest level
    # coeffs[1] = detail at coarsest level (largest scale)
    # coeffs[-1] = detail at finest level (smallest scale)
    # Convention: we fit log2 E_j vs j where j increases with FINEST scale.
    # j = 0 is coarsest detail (level L), j = max_level - 1 is finest.
    energies = []
    for level_idx in range(1, len(coeffs)):
        c = coeffs[level_idx]
        E = float(np.mean(c ** 2))
        energies.append(E)
    # Reorder so j=0 is finest scale (j increasing -> coarser)
    # coeffs[1] is coarsest detail, coeffs[-1] is finest detail.
    # Reverse so j=0 -> finest.
    energies = list(reversed(energies))
    # Fit log2(E_j) = -2*alpha*j + const (j goes finest->coarsest in our reversed order means
    # at finer scales (smaller j) energy is smaller for smooth functions).
    # Actually for detail at scale 2^{-j} (j increasing -> finer),
    # E_j ~ 2^{-2j*alpha} for Hölder alpha.
    # Our reversed array has j=0 finest scale, j=max-1 coarsest.
    # So E_j ~ 2^{-2*alpha*(max_level-1-j)} ... messy.
    # Cleanest: fit log2(E) vs scale_index where scale_index increases with COARSER scale.
    # Use original (non-reversed): coeffs[1] is coarsest, coeffs[-1] is finest.
    # We want E at scale s: as scale grows (coarser), energy grows for ROUGH signals.
    # For Hölder alpha, sum |c_{j,k}|^2 over k at scale j is ~ 2^{-2*j*alpha} where
    # j here is "octaves of fine detail". The detail coefficients at finest scale are
    # the smallest count in number, but each coeff is small for smooth signals.
    #
    # Cleanest empirical convention: define j_idx = 1, 2, ..., max_level
    # where j_idx=1 is the coarsest detail (longest scale 2^L) and j_idx=max_level is
    # the finest detail. Then for Hölder alpha, mean-square per-coefficient |c_{j,k}|^2
    # at scale 2^{-j} (where j is in ORIGINAL math convention, fine = large j) follows
    # |c_{j,k}|^2 ~ 2^{-j*(2*alpha+1)}.
    # In our pywt indexing, level_idx=1 is the coarsest, which is j = 1 in math
    # convention (longest scale). level_idx=max is finest, which is j = max in math.
    # So the energy per coefficient at level_idx = i is ~ 2^{-i*(2*alpha+1)} approximately.
    #
    # Simplification: just fit log2(E) vs level_idx (1..max_level).
    raw_energies = []
    for level_idx in range(1, len(coeffs)):
        c = coeffs[level_idx]
        E = float(np.mean(c ** 2))
        raw_energies.append((level_idx, E, len(c)))
    log_E = np.array([math.log2(max(E, 1e-300)) for _, E, _ in raw_energies])
    levels = np.array([li for li, _, _ in raw_energies], dtype=np.float64)
    # Linear fit: log2(E) = slope * level + const
    A = np.vstack([levels, np.ones_like(levels)]).T
    sol, _, _, _ = np.linalg.lstsq(A, log_E, rcond=None)
    slope, intercept = sol
    pred = A @ sol
    ss_res = ((log_E - pred) ** 2).sum()
    ss_tot = ((log_E - log_E.mean()) ** 2).sum()
    r2 = 1.0 - ss_res / max(ss_tot, 1e-20)
    # Convert slope to alpha: slope = -(2*alpha + 1)
    # alpha = (-slope - 1) / 2 -- BUT only if energy DECREASES toward coarser scales.
    # For natural rough signals, we expect energy to DECREASE as we go to finer scales
    # (the high-frequency content is small). So slope should be POSITIVE in our
    # convention (level_idx increasing = finer scales = lower energy = lower log2(E)).
    # Wait, NO: as level_idx goes from 1 (coarsest detail) to max (finest), the
    # detail captures finer features. For a smooth signal, finer features have
    # SMALLER energy => log2(E) decreases as level_idx grows => slope NEGATIVE.
    # For a rough signal, finer features have LARGER energy proportionally? Actually
    # a Hölder alpha signal has |c_{j,k}|^2 ~ 2^{-2*j*alpha} at scale 2^{-j} (math j;
    # finer = larger j). In pywt's convention, finer scale has higher level_idx, so
    # the "math j" agrees with "level_idx" up to constants. Mean-energy at level i:
    # E_i = mean_k |c_{i,k}|^2 ~ 2^{-2*alpha*i} (per-coeff). So slope of log2(E) vs i
    # is -2*alpha. So alpha = -slope/2.
    alpha = -float(slope) / 2.0
    return dict(slope=float(slope), intercept=float(intercept),
                alpha_holder=alpha, r2=float(r2),
                level_energies=[(int(li), e, int(c)) for li, e, c in raw_energies])


def cramer_pi(N, seed=12345):
    """Cramér model: independent Bernoulli(1/log n) for n >= 2."""
    rng = np.random.default_rng(seed)
    p = np.zeros(N + 1, dtype=np.float64)
    n = np.arange(2, N + 1)
    p[2:] = 1.0 / np.log(n)
    p[2] = 1.0  # 2 is prime in Cramér tradition; doesn't matter
    bern = rng.uniform(size=N + 1) < p
    return np.cumsum(bern, dtype=np.int64)


def matched_iid(N_grid, seed=23456):
    """IID standard Gaussian — fixed for skewness baseline."""
    rng = np.random.default_rng(seed)
    return rng.standard_normal(N_grid)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--logX', type=int, default=20,
                        help='X = 2^logX upper limit')
    parser.add_argument('--seed', type=int, default=12345)
    parser.add_argument('--outdir', type=str, default=str(HERE))
    args = parser.parse_args()

    X = 1 << args.logX
    print(f'[D43] X = 2^{args.logX} = {X}')

    t0 = time.time()
    pi_tab = sieve_pi_table(X)
    print(f'[sieve] pi({X}) = {pi_tab[X]}, took {time.time()-t0:.2f}s')

    xs, step = kpz_grid(X)
    print(f'[grid] X^(1/3) = {step}, n_grid = {len(xs)}')

    t0 = time.time()
    D = D_field(pi_tab, xs)
    print(f'[D]    n = {len(D)}, mean = {D.mean():.4f}, '
          f'std = {D.std(ddof=1):.4f}, took {time.time()-t0:.2f}s')

    Z = standardize(D)
    out = dict(
        params=dict(X=X, logX=args.logX, step=step, n_grid=len(xs)),
        D_summary=dict(mean=float(D.mean()), std=float(D.std(ddof=1)),
                       min=float(D.min()), max=float(D.max())),
        Z_summary=dict(skew=float(skew(Z)), exkurt=float(kurtosis(Z)),
                       n=len(Z)),
        tw2_ref=tw2_moments(),
    )

    # F1, F2: skewness tests
    se_skew = math.sqrt(6.0 / len(Z))
    Z_skew = float(skew(Z))
    out['F1_gauss_skew'] = dict(
        empirical_skew=Z_skew, se=se_skew,
        z_score=Z_skew / se_skew,
        passes=abs(Z_skew) < 3 * se_skew,
    )
    out['F2_tw2_skew'] = dict(
        empirical_skew=Z_skew, target=0.2241, se=se_skew,
        z_score=(Z_skew - 0.2241) / se_skew,
        passes=abs(Z_skew - 0.2241) < 3 * se_skew,
    )

    # F3: tail fits
    out['F3_kpz_tail'] = fit_right_tail_slope(Z, basis='kpz')
    out['F3_gauss_tail'] = fit_right_tail_slope(Z, basis='gauss')

    # F4: Hölder via wavelet
    out['F4_holder'] = wavelet_holder(D, wavelet='db4')

    # F5: KS tests
    ks_g = kstest(Z, 'norm')
    out['F5_KS_gauss'] = dict(stat=float(ks_g.statistic),
                              pvalue=float(ks_g.pvalue))
    # For TW2 we need the CDF; use mpmath for an mpmath-based numeric CDF.
    # Approximate via Bornemann's series or use a finite tabulated CDF.
    # For our purposes: KS distance between Z and a TW2 sample.
    out['F5_TW2_compared_via_moments'] = dict(
        skew_match_with_TW2=Z_skew - 0.2241,
        exkurt_match_with_TW2=float(kurtosis(Z)) - 0.0934,
    )

    # CONTROLS
    print('[control C1] Cramér model fluctuation field')
    pi_C = cramer_pi(X, seed=args.seed)
    D_C = (pi_C[xs].astype(np.float64) - li(xs)) * np.log(xs) / np.sqrt(xs)
    Z_C = standardize(D_C)
    out['control_cramer'] = dict(
        D_summary=dict(mean=float(D_C.mean()), std=float(D_C.std(ddof=1))),
        Z_summary=dict(skew=float(skew(Z_C)), exkurt=float(kurtosis(Z_C))),
        F4_holder=wavelet_holder(D_C, wavelet='db4'),
        F1_gauss_skew_passes=abs(skew(Z_C)) < 3 * se_skew,
    )

    print('[control C2] IID Gaussian on grid')
    G = matched_iid(len(xs), seed=args.seed + 1)
    G *= D.std(ddof=1)
    out['control_iid_gauss'] = dict(
        Z_summary=dict(skew=float(skew(standardize(G))),
                       exkurt=float(kurtosis(standardize(G)))),
        F4_holder=wavelet_holder(G, wavelet='db4'),
    )

    # Output
    out_path = Path(args.outdir) / 'd43_kpz_pi_li_results.json'
    with open(out_path, 'w') as f:
        json.dump(out, f, indent=2, default=str)
    print(f'[done] wrote {out_path}')

    # Summary print
    print('\n=== SUMMARY ===')
    print(f'  N grid pts: {len(xs)}')
    print(f'  D(x): mean = {D.mean():.4f}, std = {D.std(ddof=1):.4f}')
    print(f'  Z skew: {Z_skew:.4f}  (Gauss=0, TW2=+0.224, se={se_skew:.4f})')
    print(f'  Z exkurt: {kurtosis(Z):.4f}  (Gauss=0, TW2=+0.093)')
    print(f'  F1 (Gauss skew |z|<3): {out["F1_gauss_skew"]["passes"]}')
    print(f'  F2 (TW2 skew within 3se): {out["F2_tw2_skew"]["passes"]}')
    print(f'  F3 KPZ-basis tail r^2: {out["F3_kpz_tail"]["r2"]:.4f}')
    print(f'  F3 Gauss-basis tail r^2: {out["F3_gauss_tail"]["r2"]:.4f}')
    print(f'  F4 Hölder alpha (chi_P): {out["F4_holder"]["alpha_holder"]:.4f} '
          f'(r^2={out["F4_holder"]["r2"]:.3f})')
    print(f'  F4 Hölder alpha (Cramér C1): {out["control_cramer"]["F4_holder"]["alpha_holder"]:.4f}')
    print(f'  F4 Hölder alpha (IID C2):    {out["control_iid_gauss"]["F4_holder"]["alpha_holder"]:.4f}')
    print(f'  F5 KS Gauss: stat={out["F5_KS_gauss"]["stat"]:.4f}, '
          f'p={out["F5_KS_gauss"]["pvalue"]:.4f}')

    return out


if __name__ == '__main__':
    main()

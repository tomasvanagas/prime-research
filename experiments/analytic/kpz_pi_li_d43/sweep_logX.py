"""Sweep over multiple X values to test moment stability of D(x) on KPZ grid."""
import json
import math
import sys
import time
from pathlib import Path

import numpy as np
import mpmath as mp
import pywt
from scipy.stats import skew, kurtosis, kstest

sys.path.insert(0, str(Path(__file__).parent))
from d43_kpz_pi_li import (
    sieve_pi_table, kpz_grid, D_field, standardize, wavelet_holder,
    cramer_pi, fit_right_tail_slope,
)


def detrended(D, win=51):
    """Subtract a moving-average trend (window 51 points, ~5*X^{1/3}-spacing)."""
    if win % 2 == 0:
        win += 1
    half = win // 2
    pad = np.pad(D, half, mode='edge')
    kernel = np.ones(win) / win
    trend = np.convolve(pad, kernel, mode='valid')
    return D - trend


def main():
    # Largest sieve: build once for X_max = 2^24
    LOGX_MAX = 24
    XMAX = 1 << LOGX_MAX
    print(f'[sieve] building pi table for X_max = 2^{LOGX_MAX} = {XMAX}')
    t0 = time.time()
    pi_tab = sieve_pi_table(XMAX)
    print(f'[sieve] done in {time.time()-t0:.2f}s, pi(X_max) = {pi_tab[XMAX]}')

    # Build Cramér once at X_max for control
    print(f'[cramer] building Cramér model at X_max')
    pi_C = cramer_pi(XMAX, seed=42)

    rows = []
    for logX in [18, 19, 20, 21, 22, 23, 24]:
        X = 1 << logX
        xs, step = kpz_grid(X)
        # Use mpmath li only for the full grid points
        D = (pi_tab[xs].astype(np.float64) - np.array(
                 [float(mp.li(int(xv))) for xv in xs]
             )) * np.log(xs) / np.sqrt(xs)
        Z = standardize(D)
        Dd = detrended(D, win=51)
        Zd = standardize(Dd)
        # Cramér control
        D_C = (pi_C[xs].astype(np.float64) - np.array(
                   [float(mp.li(int(xv))) for xv in xs]
               )) * np.log(xs) / np.sqrt(xs)
        Z_C = standardize(D_C)
        Dd_C = detrended(D_C, win=51)
        Zd_C = standardize(Dd_C)

        # Wavelet Hölder
        H_chi = wavelet_holder(D, wavelet='db4')
        H_det = wavelet_holder(Dd, wavelet='db4')
        H_C = wavelet_holder(D_C, wavelet='db4')

        row = dict(
            logX=logX,
            X=X,
            step=step,
            n_grid=len(xs),
            D_mean=float(D.mean()),
            D_std=float(D.std(ddof=1)),
            # Whole-window standardised moments
            skew_Z=float(skew(Z)),
            exkurt_Z=float(kurtosis(Z)),
            # Detrended (moving-avg residuals) moments — kills the Skewes bias drift
            skew_Zd=float(skew(Zd)),
            exkurt_Zd=float(kurtosis(Zd)),
            # Hölder
            alpha_chi=H_chi['alpha_holder'],
            alpha_chi_r2=H_chi['r2'],
            alpha_chi_detrend=H_det['alpha_holder'],
            # Cramér control
            cramer_skew_Z=float(skew(Z_C)),
            cramer_exkurt_Z=float(kurtosis(Z_C)),
            cramer_skew_Zd=float(skew(Zd_C)),
            cramer_exkurt_Zd=float(kurtosis(Zd_C)),
            cramer_alpha=H_C['alpha_holder'],
            # Standard error for skew
            se_skew=math.sqrt(6.0 / len(xs)),
        )
        rows.append(row)
        print(f'logX={logX:2d} n={len(xs):6d} '
              f'skew_Z={row["skew_Z"]:+.4f} skew_Zd={row["skew_Zd"]:+.4f} '
              f'exkurt_Zd={row["exkurt_Zd"]:+.4f} '
              f'alpha_chi={row["alpha_chi"]:.4f} '
              f'C_skew_Zd={row["cramer_skew_Zd"]:+.4f} '
              f'C_exkurt_Zd={row["cramer_exkurt_Zd"]:+.4f}')

    out = dict(rows=rows, tw2=dict(skew=0.2241, exkurt=0.0934))
    with open(Path(__file__).parent / 'sweep_logX_results.json', 'w') as f:
        json.dump(out, f, indent=2, default=str)
    print('[done] wrote sweep_logX_results.json')


if __name__ == '__main__':
    main()

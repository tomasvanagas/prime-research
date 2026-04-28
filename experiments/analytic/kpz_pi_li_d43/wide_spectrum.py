"""Wide-range spectral & wavelet diagnostics on D(x).

Uses x in [10^4, 10^7] (u-span ~6.9, ~15 oscillations of γ_1=14.13) so the
FFT resolves the dominant zeta-zero frequencies. Compares to Cramér model.
"""
import math
import json
import sys
from pathlib import Path

import numpy as np
import mpmath as mp
import pywt
from scipy.stats import skew, kurtosis

sys.path.insert(0, str(Path(__file__).parent))
from d43_kpz_pi_li import sieve_pi_table, cramer_pi, wavelet_holder
from spectral_signature import ZETA_ZEROS_GAMMA


def detrended(D, win=51):
    if win % 2 == 0:
        win += 1
    half = win // 2
    pad = np.pad(D, half, mode='edge')
    kernel = np.ones(win) / win
    trend = np.convolve(pad, kernel, mode='valid')
    return D - trend


def main():
    XMAX = 10_000_000
    XMIN = 10_000
    print(f'[wide] X range = [{XMIN:.0e}, {XMAX:.0e}]')
    pi_tab = sieve_pi_table(XMAX)
    pi_C = cramer_pi(XMAX, seed=1729)

    # KPZ-spaced grid over the full range. Step = X^{1/3} for X = XMAX.
    step = max(1, int(math.floor(XMAX ** (1.0 / 3.0))))
    n_pts = (XMAX - XMIN) // step
    xs = np.array([XMIN + k * step for k in range(1, n_pts + 1)],
                  dtype=np.int64)
    print(f'[wide] step = {step}, n_grid = {len(xs)}')

    li_vals = np.array([float(mp.li(int(xv))) for xv in xs])
    pi_vals = pi_tab[xs].astype(np.float64)
    D = (pi_vals - li_vals) * np.log(xs) / np.sqrt(xs)
    print(f'[chi_P]  mean(D) = {D.mean():.4f}, std(D) = {D.std(ddof=1):.4f}')
    pi_C_vals = pi_C[xs].astype(np.float64)
    D_C = (pi_C_vals - li_vals) * np.log(xs) / np.sqrt(xs)
    print(f'[Cramér] mean(D_C) = {D_C.mean():.4f}, std(D_C) = {D_C.std(ddof=1):.4f}')

    # Spectrum on a uniform u-grid
    u = np.log(xs.astype(np.float64))
    n_u = len(u)
    u_uni = np.linspace(u[0], u[-1], n_u)
    D_uni = np.interp(u_uni, u, D - D.mean())
    DC_uni = np.interp(u_uni, u, D_C - D_C.mean())
    F = np.fft.rfft(D_uni)
    F_C = np.fft.rfft(DC_uni)
    freqs = np.fft.rfftfreq(n_u, d=(u_uni[1] - u_uni[0]))
    gamma = 2 * math.pi * freqs
    P = np.abs(F) ** 2 / n_u
    P_C = np.abs(F_C) ** 2 / n_u
    res_gamma = gamma[1] - gamma[0]
    print(f'[spec] u-span = {u[-1] - u[0]:.4f}, gamma resolution = {res_gamma:.4f}')

    # Find peaks near zeta-zero gammas (use halfwidth = 1.5 * res_gamma)
    halfwidth = 1.5 * res_gamma
    median_P = float(np.median(P[gamma > 1.0]))
    median_P_C = float(np.median(P_C[gamma > 1.0]))
    peaks_chi, peaks_C = [], []
    for g in ZETA_ZEROS_GAMMA:
        mask = np.abs(gamma - g) < halfwidth
        if mask.sum() == 0:
            peaks_chi.append({'gamma_target': g, 'ratio': 0.0})
            peaks_C.append({'gamma_target': g, 'ratio': 0.0})
            continue
        peaks_chi.append({'gamma_target': g,
                          'gamma_peak': float(gamma[mask][np.argmax(P[mask])]),
                          'P_peak': float(P[mask].max()),
                          'ratio': float(P[mask].max() / median_P)})
        peaks_C.append({'gamma_target': g,
                        'gamma_peak': float(gamma[mask][np.argmax(P_C[mask])]),
                        'P_peak': float(P_C[mask].max()),
                        'ratio': float(P_C[mask].max() / median_P_C)})

    print(f'\nFirst 12 zeta-zero peaks: peak/baseline ratio')
    print(f'   {"γ":>8s}  {"chi_P":>10s}  {"Cramér":>10s}')
    for pc, pcc in zip(peaks_chi[:12], peaks_C[:12]):
        print(f'   {pc["gamma_target"]:>8.4f}  {pc["ratio"]:>10.2f}  {pcc["ratio"]:>10.2f}')

    # Hölder regularity: chi_P (raw + detrended), Cramér (raw + detrended)
    H_chi = wavelet_holder(D, wavelet='db4')
    H_C = wavelet_holder(D_C, wavelet='db4')
    Dd = detrended(D)
    DdC = detrended(D_C)
    H_chi_d = wavelet_holder(Dd, wavelet='db4')
    H_C_d = wavelet_holder(DdC, wavelet='db4')
    print(f'\nHölder regularity:')
    print(f'  chi_P  raw      α = {H_chi["alpha_holder"]:.4f}  r² = {H_chi["r2"]:.4f}')
    print(f'  chi_P  detrend  α = {H_chi_d["alpha_holder"]:.4f}  r² = {H_chi_d["r2"]:.4f}')
    print(f'  Cramér raw      α = {H_C["alpha_holder"]:.4f}  r² = {H_C["r2"]:.4f}')
    print(f'  Cramér detrend  α = {H_C_d["alpha_holder"]:.4f}  r² = {H_C_d["r2"]:.4f}')

    Z = (D - D.mean()) / D.std(ddof=1)
    Z_C = (D_C - D_C.mean()) / D_C.std(ddof=1)
    Zd = (Dd - Dd.mean()) / Dd.std(ddof=1)
    ZdC = (DdC - DdC.mean()) / DdC.std(ddof=1)
    print(f'\nMoments:')
    print(f'  chi_P  raw       skew = {skew(Z):+.4f}  exkurt = {kurtosis(Z):+.4f}')
    print(f'  chi_P  detrend   skew = {skew(Zd):+.4f}  exkurt = {kurtosis(Zd):+.4f}')
    print(f'  Cramér raw       skew = {skew(Z_C):+.4f}  exkurt = {kurtosis(Z_C):+.4f}')
    print(f'  Cramér detrend   skew = {skew(ZdC):+.4f}  exkurt = {kurtosis(ZdC):+.4f}')
    print(f'  TW2              skew = +0.2241 exkurt = +0.0934')
    print(f'  Gauss            skew = +0.0000 exkurt = +0.0000')
    se_skew_chi = math.sqrt(6.0 / len(Z))

    # Save everything
    out = dict(
        params=dict(XMIN=XMIN, XMAX=XMAX, step=step, n_grid=len(xs)),
        u_span=float(u[-1] - u[0]),
        gamma_resolution=float(res_gamma),
        D_summary=dict(chi_P=dict(mean=float(D.mean()), std=float(D.std(ddof=1))),
                       cramer=dict(mean=float(D_C.mean()),
                                   std=float(D_C.std(ddof=1)))),
        zeta_zero_peaks_chi_P=peaks_chi,
        zeta_zero_peaks_cramer=peaks_C,
        spectrum_baseline_P_chi=median_P,
        spectrum_baseline_P_cramer=median_P_C,
        holder=dict(
            chi_P_raw=H_chi['alpha_holder'],
            chi_P_detrend=H_chi_d['alpha_holder'],
            cramer_raw=H_C['alpha_holder'],
            cramer_detrend=H_C_d['alpha_holder'],
        ),
        moments=dict(
            chi_P_raw=dict(skew=float(skew(Z)), exkurt=float(kurtosis(Z))),
            chi_P_detrend=dict(skew=float(skew(Zd)), exkurt=float(kurtosis(Zd))),
            cramer_raw=dict(skew=float(skew(Z_C)), exkurt=float(kurtosis(Z_C))),
            cramer_detrend=dict(skew=float(skew(ZdC)), exkurt=float(kurtosis(ZdC))),
            tw2_target=dict(skew=0.2241, exkurt=0.0934),
            se_skew=se_skew_chi,
        ),
    )
    out_path = Path(__file__).parent / 'wide_spectrum_results.json'
    with open(out_path, 'w') as f:
        json.dump(out, f, indent=2, default=str)
    print(f'\n[done] wrote {out_path}')


if __name__ == '__main__':
    main()

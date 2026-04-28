"""Compute the spectral / Fourier signature of D(x) to expose its deterministic
almost-periodic structure (zeta-zero oscillations) vs a true stochastic field.

This is the KEY test that distinguishes KPZ-class fields (broadband white-noise-
like spectrum + roughness) from explicit-formula-driven D (sharp zeta-zero peaks
at γ_k log x frequencies + smooth trend).
"""
import json
import math
import sys
import time
from pathlib import Path

import numpy as np
import mpmath as mp
import pywt
from scipy.stats import skew, kurtosis

sys.path.insert(0, str(Path(__file__).parent))
from d43_kpz_pi_li import sieve_pi_table, kpz_grid, cramer_pi


def compute_D(pi_table, xs):
    li_vals = np.array([float(mp.li(int(xv))) for xv in xs])
    return (pi_table[xs].astype(np.float64) - li_vals) * np.log(xs) / np.sqrt(xs)


def power_spectrum(D, log_xs):
    """Power spectrum of D as a function of u = log(x).

    Resample D onto a uniform u-grid (interpolation) so we can FFT in
    log-space; the explicit formula predicts dominant Fourier components at
    frequencies γ_k (the imaginary parts of nontrivial ζ-zeros).
    """
    # Resample on a uniform u-grid
    u = log_xs
    n_u = len(u)
    u_uniform = np.linspace(u[0], u[-1], n_u)
    D_uniform = np.interp(u_uniform, u, D)
    D_uniform = D_uniform - D_uniform.mean()
    # FFT
    F = np.fft.rfft(D_uniform)
    freqs = np.fft.rfftfreq(n_u, d=(u_uniform[1] - u_uniform[0]))
    # Convert frequency to "γ" units: u-space frequency f corresponds to γ = 2π f
    gamma = 2 * math.pi * freqs
    P = np.abs(F) ** 2 / n_u
    return gamma, P, F


# First few non-trivial zeta zeros (Odlyzko table truncated; well-known values)
ZETA_ZEROS_GAMMA = [
    14.134725, 21.022040, 25.010858, 30.424876, 32.935062, 37.586178,
    40.918719, 43.327073, 48.005151, 49.773832, 52.970321, 56.446248,
    59.347044, 60.831779, 65.112544, 67.079811, 69.546402, 72.067158,
    75.704691, 77.144840,
]


def find_peaks_near(gamma_axis, P, target_gammas, halfwidth=0.5):
    """For each target γ_k, return max P in window |γ - γ_k| < halfwidth and
    median P over the broader spectrum (off-peak baseline)."""
    out = []
    median_P = float(np.median(P[gamma_axis > 1.0]))
    for g in target_gammas:
        mask = np.abs(gamma_axis - g) < halfwidth
        if mask.sum() == 0:
            out.append(dict(gamma_target=g, gamma_peak=None,
                            P_peak=0.0, P_baseline=median_P, ratio=0.0))
            continue
        peak = float(P[mask].max())
        peak_idx = np.argmax(P * mask)
        peak_gamma = float(gamma_axis[peak_idx])
        out.append(dict(gamma_target=g, gamma_peak=peak_gamma,
                        P_peak=peak, P_baseline=median_P,
                        ratio=peak / max(median_P, 1e-30)))
    return out, median_P


def main():
    LOGX = 24
    X = 1 << LOGX
    print(f'[D43-spec] X = 2^{LOGX} = {X}')
    pi_tab = sieve_pi_table(X)
    pi_C = cramer_pi(X, seed=99)
    xs, step = kpz_grid(X)
    log_xs = np.log(xs.astype(np.float64))

    D = compute_D(pi_tab, xs)
    D_C = compute_D(pi_C, xs)
    print(f'[chi_P]  D mean = {D.mean():.4f}, std = {D.std(ddof=1):.4f}')
    print(f'[Cramér] D mean = {D_C.mean():.4f}, std = {D_C.std(ddof=1):.4f}')

    gamma_chi, P_chi, _ = power_spectrum(D, log_xs)
    gamma_C, P_C, _ = power_spectrum(D_C, log_xs)

    peaks_chi, base_chi = find_peaks_near(gamma_chi, P_chi, ZETA_ZEROS_GAMMA,
                                          halfwidth=0.5)
    peaks_C, base_C = find_peaks_near(gamma_C, P_C, ZETA_ZEROS_GAMMA,
                                      halfwidth=0.5)
    print(f'\n[chi_P spectrum] median (γ>1) baseline P = {base_chi:.4e}')
    print(f'[Cramér spectrum] median (γ>1) baseline P = {base_C:.4e}')
    print(f'\nFirst 10 zeta-zero peak ratios (peak / baseline):')
    print(f'   {"γ":>8s}  {"chi_P ratio":>12s}  {"Cramér ratio":>12s}')
    for pchi, pC in zip(peaks_chi[:10], peaks_C[:10]):
        print(f'   {pchi["gamma_target"]:>8.4f}  '
              f'{pchi["ratio"]:>12.2f}  {pC["ratio"]:>12.2f}')

    # Total "spike energy" at zeta-zero locations vs total spectrum energy
    # mask out the DC / very low frequency
    valid = gamma_chi > 1.0
    total_chi = float(P_chi[valid].sum())
    total_C = float(P_C[valid].sum())
    spike_chi = sum(p['P_peak'] for p in peaks_chi[:20] if p['P_peak'])
    spike_C = sum(p['P_peak'] for p in peaks_C[:20] if p['P_peak'])
    print(f'\n[chi_P]  spike (top 20 zeros) / total = {spike_chi/total_chi:.4f}')
    print(f'[Cramér] spike (top 20 zeros) / total = {spike_C/total_C:.4f}')

    out = dict(
        params=dict(logX=LOGX, X=X, step=step, n_grid=len(xs)),
        D_summary=dict(chi_P=dict(mean=float(D.mean()), std=float(D.std(ddof=1))),
                       cramer=dict(mean=float(D_C.mean()),
                                   std=float(D_C.std(ddof=1)))),
        baseline_P_chi=base_chi,
        baseline_P_cramer=base_C,
        zeta_zero_peaks_chi_P=peaks_chi,
        zeta_zero_peaks_cramer=peaks_C,
        spike_fraction_chi_P=spike_chi / total_chi,
        spike_fraction_cramer=spike_C / total_C,
        total_spectrum_energy_chi=total_chi,
        total_spectrum_energy_cramer=total_C,
    )
    out_path = Path(__file__).parent / 'spectral_signature_results.json'
    with open(out_path, 'w') as f:
        json.dump(out, f, indent=2, default=str)
    print(f'\n[done] wrote {out_path}')


if __name__ == '__main__':
    main()

# sweep_logX — multi-X sweep for D43

Auxiliary script for the D43 KPZ-universality experiment. Sweeps the
KPZ-grid measurement at logX ∈ {18, 19, 20, 21, 22, 23, 24} and
records: whole-window standardised skew/exkurt of D, detrended skew/
exkurt (moving-avg window 51 KPZ-grid points), wavelet Hölder α (raw
+ detrended), and Cramér-control values at each X.

Purpose: verify that the wavelet Hölder α(D) ≈ 0.85 measurement is
**stable across logX**, ruling out a single-window artifact.

**Key result** (full table in `../d43_kpz_pi_li_results.md`):
α(D, raw) ranges in [0.810, 0.867] across the 7 logX values, with
linear-fit r² > 0.998 in all cases. The detrended skew(Z_d) ranges in
[-0.062, +0.031] — Gauss-consistent. KPZ rejected uniformly.

See `sweep_logX_results.json` for raw numbers and
`d43_kpz_pi_li_results.md` for the consolidated narrative.

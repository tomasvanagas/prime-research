# spectral_signature — short-window FFT diagnostic for D43

Auxiliary diagnostic for the D43 KPZ experiment. Runs at logX = 24,
KPZ-grid window = [X/2, X], u-span = log(2) ≈ 0.693 (small).

**Purpose**: attempt to detect zeta-zero peaks γ_k in the FFT of D(x).

**Result**: at this window the u-span (0.693) is too narrow — FFT
gamma resolution is ~9, comparable to gamma_1=14.13. Peaks not
clearly resolved (peak/median ratios all 0.0 with halfwidth 0.5).

The wider-range run `wide_spectrum.py` extends u-span to 6.89 and
DOES resolve γ_1, γ_2, γ_3 cleanly. See `wide_spectrum_results.md`
and `../d43_kpz_pi_li_results.md` for the consolidated spectral
discussion.

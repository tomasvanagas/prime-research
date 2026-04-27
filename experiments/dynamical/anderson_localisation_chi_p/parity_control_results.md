# parity_control.py — Results

Auxiliary script for the §C4 Anderson localisation attack. See the
primary writeup `anderson_localisation_chi_p_results.md` (Result 2) for
the full discussion.

**Summary:** parity-matched controls (random odd-subset of size pi(N)-1
plus V(2)=1, and chi_P shuffled within odd indices) reduce the naive
88 sigma deviation to 23-33 sigma, leaving a residual peak at
E = 1.088 ~ -2 cos(2 pi / 3) = +1, the **mod-3 resonance**. This
isolates the parity confound and exposes the next-layer structure
(mod-3 / mod-6 sieving) that motivates the W-trick cascade in
`wtrick_control.py`.

Output JSON: `parity_N{N}_s{seeds}_e{energies}_a{amp}.json`.

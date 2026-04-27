# Sub-window Robustness Test

Auxiliary diagnostic script for ATTACK_VECTORS.md §C5 wild swing
(S108). Aggregated results live in `stein_wasserstein_pi_results.md`.

## What this script does

`test_low_zero_robustness.py` partitions `[10^6, 10^7]` into 10
overlapping sub-windows of width 0.5 in log10. For each sub-window:

1. Computes empirical `W_1(D̂_emp window, fitted Gaussian)`.
2. Computes theoretical `W_1(D_th(50 zeros) window, fitted Gaussian)`.
3. Reports skewness and excess kurtosis.

A separate "random-phase null" replaces the actual Riemann γ_k phases
with uniform random phases (50 trials).

## Key output

Sub-window (10 partitions of [10^6, 10^7]):
- Empirical `W_1`: mean = 0.01257, std = 0.00400 (range 0.0075–0.0211)
- Theoretical `W_1`: mean = 0.01785, std varies similarly
- **Pearson correlation r = 0.906** (n = 10 windows)

Skew/kurt patterns: empirical `skew ∈ [−0.08, +0.25]`,
`kurt ∈ [−0.62, −0.27]` — sub-Gaussian across all sub-windows.

Random-phase null (50 trials, same `1/|ρ|` weights):
- W_1 mean = 0.01158, std = 0.00310 across 50 trials
- Empirical W_1 = 0.00829 → z = -1.06 vs random-phase null

## Interpretation

The strong correlation (r = 0.906) of `W_1` magnitudes across
sub-windows between empirical and explicit-formula prediction is
the cleanest *quantitative* test of the structural explanation.
The sub-windows have different low-zero-mode interference patterns
(due to different `log x` ranges), and the empirical W_1 *follows*
the predicted variation.

The random-phase null indicates that the actual zero phases give
results within 1σ of a generic random-phase low-zero sum — D̂'s
behaviour is not distinguishable from the explicit-formula
prediction with arbitrary phases.

## Falsifier

If r < 0.5 across 10+ sub-windows, the structural matching
collapses and the plateau `c(X)` becomes a coincidental empirical
constant rather than an explicit-formula consequence.

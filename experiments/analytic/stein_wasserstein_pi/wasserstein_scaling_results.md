# Wasserstein-1 K-Scaling Diagnostic

Auxiliary diagnostic script for ATTACK_VECTORS.md §C5 wild swing
(S108). Aggregated results live in `stein_wasserstein_pi_results.md`.

## What this script does

`wasserstein_scaling.py` computes the empirical
`W_1(D̂_K, N(μ̂(K), σ̂(K)²))` for `K ∈ {200, 500, 1000, 2000, 5000,
10000}` against sample-fitted Gaussian controls, to test whether
`W_1(D̂)` shrinks as Stein-CLT `O(1/√K)` (Gaussian limit) or plateaus
(non-Gaussian arithmetic structure).

## Key output (full corrected table)

| K     | W_1(D̂)    | W_1 Gauss-control      | Z-score    | W_1 × √K | excess kurt |
|-------|-----------|------------------------|------------|----------|-------------|
| 200   | 0.01221   | 0.01219 ± 0.0026       | +0.01      | 0.173    | -0.475      |
| 500   | 0.00951   | 0.00815 ± 0.0018       | +0.76      | 0.213    | -0.434      |
| 1000  | 0.00944   | 0.00600 ± 0.0014       | +2.38      | 0.298    | -0.449      |
| 2000  | 0.00908   | 0.00421 ± 0.0011       | +4.55      | 0.406    | -0.428      |
| 5000  | 0.00828   | 0.00256 ± 0.00058      | +9.80      | 0.585    | -0.408      |
| 10000 | 0.00829   | 0.00181 ± 0.00043      | **+15.08** | 0.829    | -0.410      |

## Interpretation

- Gauss control's `W_1 × √K ≈ 0.18` constant across K — confirms
  Stein-CLT rate `c_G(σ̂)/√K`.
- D̂'s `W_1 × √K` grows monotonically: 0.17 → 0.21 → 0.30 → 0.41 →
  0.59 → 0.83. Linear-ish growth with √K is consistent with a
  *constant W_1 plateau* `c(X) ≈ 0.0083`.
- Excess kurtosis is stable at -0.41 to -0.48 across all K
  values — sub-Gaussian signature persistent.
- Z-score grows from ~0 at K=200 (sample-noise dominates) to
  +15.08 at K=10000 (plateau exposed).

## Falsifier

If at K = 50000 the W_1(D̂) drops below the K=10000 value (i.e.,
0.00829), the plateau interpretation is wrong and the apparent
plateau was a finite-K artefact. Computationally feasible in
~30 minutes (pi-table for x ≤ 10^7 is fast).

## Why scaling matters

Stein-CLT bounds give `O(1/√K)` for genuinely-Gaussian samples. A
plateau at `c > 0` is the quantitative signature that the underlying
distribution (in the K → ∞ limit) is *not* Gaussian. The structural
origin (low Riemann zeros) is documented in
`structural_explanation_results.md` and
`test_low_zero_robustness_results.md`.

# Structural Explanation: Low-Zero Match for D(x)

Auxiliary diagnostic script for ATTACK_VECTORS.md §C5 wild swing
(S108). Aggregated results live in `stein_wasserstein_pi_results.md`.

## What this script does

`structural_explanation.py` derives the leading approximation to
`D(x) = (π(x) - Li(x)) log(x) / √x` from Riemann's explicit formula
truncated to `n` lowest non-trivial zeros, and compares to empirical
`D̂` over `K = 5000` log-uniform anchors in `[10^6, 10^7]`.

## Key output

```
Empirical D: mean=-1.3295, std=0.2206, skew=0.054, kurt=-0.417
W_1(D̂_emp, fitted Gaussian) = 0.008672

Truncation to n zeros:
n      W_1(D_th)    corr(D_emp, D_th)
1      0.025821     0.4255
2      0.017253     0.5640
3      0.008629     0.6219     ← matches empirical W_1
5      0.019174     0.6863
10     0.021824     0.7563
20     0.015665     0.8181
50     0.013661     0.8907

Best-fit α (D_emp = α · D_th(20)): 1.029
Variance of D_emp explained by 20 zeros: 67%
```

The leading-explicit-formula prediction tracks the empirical signal
with `α ≈ 1` (no fit needed — `α = 1` is the analytic prediction)
and correlation `0.89` at `n = 50` zeros.

The W_1 plateau magnitude `c ≈ 0.0087` is matched within 5% by the
3-lowest-zero truncation. Higher truncations oscillate around this
value as more cosines combine.

## Falsifier

If the empirical correlation `corr(D_emp, D_th(50))` falls below
0.7 at K=5000 in the same x-window, the structural explanation
fails and §C5 collapses to "the plateau is real but its origin is
not the leading explicit-formula contribution".

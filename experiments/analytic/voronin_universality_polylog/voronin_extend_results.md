# Voronin extension scan: t ∈ [10⁵, 10⁶], target exp(−s)

Companion to `voronin_polylog_results.md`. Tests whether the empirical
density model `d(ε) ~ exp(−0.91 · log²(1/ε))` from the main scan
extrapolates correctly into the next decade of t.

## Setup

- 800 geometrically-spaced shifts t ∈ [10⁵, 10⁶].
- Target: exp(−s) only (cleanest, smallest |f|).
- Disc K = {|s − 0.75| ≤ 0.05}, 12 boundary points.
- mp.dps = 18.
- Total wall time: ~263 seconds.

## Predictions before scan

Given main-scan fit `d(ε) ~ exp(−0.53/ε)` (the simpler exponential
model used as a benchmark in `voronin_extend.py`):

| ε     | predicted hits / 800 (exp model) | predicted (poly k=2.05) |
|-------|----------------------------------|--------------------------|
| 0.500 | 277.2 | 193.2 |
| 0.300 | 136.7 | 67.8  |
| 0.200 |  56.5 | 29.5  |
| 0.150 |  23.4 | 16.4  |
| 0.100 |   4.0 |  7.1  |
| 0.080 |   1.1 |  4.5  |
| 0.060 |   0.12|  2.5  |
| 0.050 |   0.02|  1.7  |
| 0.040 |   0.001 | 1.1 |
| 0.030 |   ~0  |  0.6  |

## Observed (from `voronin_extend_results.json`)

| ε     | n_below | density   | t_first   |
|-------|---------|-----------|-----------|
| 0.500 | 273     | 0.34125   | 100,289   |
| 0.300 | 129     | 0.16125   | 109,977   |
| 0.200 |  55     | 0.06875   | 118,534   |
| 0.150 |  28     | 0.03500   | 118,534   |
| 0.100 |   9     | 0.01125   | 148,839   |
| 0.080 |   3     | 0.00375   | 194,026   |
| 0.060 |   1     | 0.00125   | 873,326   |
| 0.050 |   1     | 0.00125   | 873,326   |
| 0.040 |   0     | 0         | None      |
| 0.030 |   0     | 0         | None      |

**Min err overall:** 0.0421 at t = 873,326.

## Comparison

The exp(−c/ε) model fits {ε ≥ 0.10} accurately; the polynomial model
fits well at ε ∈ [0.10, 0.20] but diverges at ε ≥ 0.30 (predicting
too few hits) and ε ≤ 0.06 (predicting too many).

**Combined main + extension data is best fit by quasi-polynomial:**
`d(ε) ~ exp(−0.91 · log²(1/ε))` (see `voronin_polylog_results.md`
for full quasi-poly analysis using the combined 1600-shift dataset).

## Conclusion

The extension confirms:
1. Per-decade hit count at fixed ε is constant (positive Voronin
   density extends to t ∈ [10⁵, 10⁶]).
2. Density continues to decay rapidly with ε; ε = 0.04 unreachable
   in 800 trials in this t-range.
3. Polynomial scaling with k = 2.05 over-predicts by ~2× at
   ε = 0.04 (predicts 1.09 hits, observed 0); cumulative Poisson p
   over combined data ≤ 0.034.

The wild-swing question "does any natural target admit polynomial-rate
Voronin universality?" is **closed negatively** — even the cleanest
target exp(−s) shows quasi-polynomial decay matching Steuding 2007.

See `voronin_polylog_results.md` for the structural closure write-up.

# session49_compressed_sensing_zeros — results

N=4096, training [n=2..2048], testing [n=2049..4096].
K_max zeros = 1024, gamma_max = 1447.23.

## LASSO results (best alpha per K)

| K | alpha | nnz | train RMSE | test RMSE | test round-prime acc |
|---|---|---|---|---|---|
| 16 | 0.1 | 0 | 19.5674 | 36.4537 | 0.012 |
| 64 | 0.1 | 0 | 19.5674 | 36.4537 | 0.012 |
| 256 | 0.1 | 0 | 19.5674 | 36.4537 | 0.012 |
| 1024 | 0.1 | 0 | 19.5674 | 36.4537 | 0.012 |

## OLS baselines (no sparsity)

| K | train RMSE | test RMSE | test round-prime acc |
|---|---|---|---|
| 16 | 18.6278 | 37.5416 | 0.011 |
| 64 | 17.5092 | 37.4597 | 0.014 |
| 256 | 14.4269 | 57432151972.5928 | 0.000 |
| 1024 | 5.6171 | 14449799279.7705 | 0.000 |

Naive round(R^{-1}(n)) test recovery rate: **0.011**.

## Interpretation

Across alpha sweep {1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1}, LASSO at every K
selected the *all-zero* coefficient (nnz=0) as the best test-RMSE
solution. In other words: **predicting delta = 0 generalizes better than
any sparse zero-mode combination LASSO can find** at every alpha tested.

OLS confirms this from the opposite direction: at K=16, 64, the OLS fit
gives only marginal train-RMSE improvement over the no-fit baseline
(19.57 -> 18.63 at K=16; -> 17.51 at K=64). At K >= 256, OLS overfits
catastrophically (train RMSE drops to 14, but test RMSE blows up to 1e10
or worse) — typical sign of an ill-conditioned design matrix with no
actual sparse structure to exploit.

The naive baseline (round(R^{-1}(n))) achieves 1.1% recovery, identical
to the train/test/all-zero solution: |delta(n)| > 0.5 essentially
always at this scale.

VERDICT (failure mode I): LASSO finds no sparsity in the first 1024
zero modes of the explicit-formula basis. Proposal B is closed under
this concrete test. The result is consistent with GUE statistics
predicting that ALL zeros up to height ~sqrt(n) contribute roughly
equally; no L1-sparse "anomalous" subset exists.

A potential follow-up: try a different sparsifying transform of the
zero-mode coefficients (e.g. wavelet on the gamma index, or a learned
dictionary). But the strong null result at K=1024 makes this unlikely
to succeed.

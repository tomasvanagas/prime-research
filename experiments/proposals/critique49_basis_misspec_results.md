# critique49_basis_misspec — results

Critic test: did Proposal C's `cos(gamma log n)/sqrt(n)` basis
(prime index n) hide a positive result that the **correct** basis
`sqrt(p_n) * cos(gamma log p_n)` (prime value) would surface?

N=4096, train [n=2..2048], test [n=2049..4096].

| K | orig acc | correct acc | orig test RMSE | correct test RMSE |
|---|---|---|---|---|
| 32 | 0.015 | 0.014 | 37.395 | 25.325 |
| 128 | 0.011 | 0.017 | 37.387 | 23.655 |
| 512 | 0.012 | 0.000 | 38.122 | 1670.741 |
| 2000 | 0.011 | 0.000 | 38.274 | 795.423 |

Naive round(R^{-1}) baseline: **0.011**.

## Verdict

Both bases plateau near naive (0.017 vs 0.015); misspecification did NOT hide a positive result. Proposal C closed.

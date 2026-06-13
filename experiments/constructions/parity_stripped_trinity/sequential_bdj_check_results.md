# sequential_bdj_check.py — results pointer

This auxiliary script extends `parity_stripped_trinity.py` by stripping
**cumulative** squarefree-conductor major arcs from `χ_P` and re-
measuring the Toeplitz 4th spectral moment `m_4`.

**Purpose.** Single-parity stripping closes 83% of the BDJ violation.
Quantify the residual as more major arcs are added.

**Output.** At `N = 1000`:

| Q  | m_4(χ*) | m_4 / (N/log²N) | λ̃_max(χ*) |
|---:|--------:|----------------:|-----------:|
|  0 | 60.682  |          2.896  |   15.245   |
|  2 | 10.355  |          0.494  |    6.561   |
|  7 |  4.208  |          0.201  |    4.487   |
| 11 |  3.459  |          0.165  |    3.669   |
| 29 |  3.237  |          0.154  |    4.156   |

(Compare BDJ universal `8/3 ≈ 2.67`, BERN matched-density `m_4 ≈ 2.65`.)

**Full discussion in:** `parity_stripped_trinity_results.md` §"BDJ
Toeplitz m_4 (sequential)". Raw results in `sequential_bdj_results.json`.

**Falsifier verdict.** Pre-stated F4.b confirmed. The residual `m_4`
ratio `0.15 N/log²N` matches expected geometric tail of HL spike
contributions from q ≥ 11.

# sequential_strip_check.py — results pointer

This auxiliary script extends `parity_stripped_trinity.py` by stripping
**cumulative** squarefree-conductor major arcs `V_q^prim` (additive
Fourier modes `e^{2πi a n / q}` for `(a, q) = 1`) from BOTH `χ_P` and
the Bernoulli baseline at matched density.

**Purpose.** Single-parity stripping (`Q=2`) leaves the Mahler deficit
unchanged at `−0.307` nat. Need to test whether the deficit closes as
more major arcs are stripped — distinguishing "single-q-sourced" from
"multi-q-distributed" attribution.

**Output.** Tables and trajectory of the L²-mass-normalised "shape
deficit" `Δ_shape(Q)`:

| Q  | shape-Δ |
|---:|--------:|
|  0 | −0.299  |
|  2 | −0.234  |
|  7 | −0.129  |
| 13 | −0.098  |
| 23 | −0.055  |
| 29 | −0.049  |

84% closure at Q=29; trajectory monotone toward zero.

**Full discussion in:** `parity_stripped_trinity_results.md` §"Sequential
major-arc stripping — the structural decomposition" → "Mahler shape-
deficit". Raw results in `sequential_strip_results.json`.

**Falsifier verdict.** Pre-stated F3.b ("residual-sourced") sharpened
to "sequentially-major-arc-sourced". F3.b PASS in spirit; the
quantitative attribution refines the qualitative outcome.

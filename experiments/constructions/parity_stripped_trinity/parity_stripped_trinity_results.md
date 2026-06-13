# Parity-stripped trinity — results

**Session:** S239 (paradigm-shift mode, 2026-04-30).
**Edges composed:** E2.20 (Mahler deficit `−0.307` nat for `χ_P`),
E2.21 (Newman L^∞ peak at `z = −1` from q=2 parity major arc),
E2.31 (Toeplitz `m_4 ≈ 2.95 N/log²N`, attributed in S204 to parity
spike), E1.6 + E2.10 (parity bisection / `L mod 2 = x mod 2`).
**Self-grade:** B (refines E2.20 + E2.31; new structural hierarchy
of major-arc concentration across three previously-disjoint statistics).
**No external technique imported.**

## Construction recap

**Parity-strip operator.** Given `χ_P : {1, ..., N} → {0, 1}`, define
`α_2(N) := (1/N) Σ χ_P(n)(−1)^n`; set `χ̃_P := χ_P − α_2(−1)^n`.
By construction `f̃_N(−1) = 0` (q=2 parity major arc nulled).

**Sequential major-arc-strip operator.** For squarefree `Q`-list
`{1, 2, 3, 5, 6, 7, ...}`, project out the additive-Fourier subspaces
`V_q^prim = span{e^{2πi a n/q} : (a, q) = 1}` for each `q ≤ Q`.
At Q=∞, all squarefree-conductor major arcs are removed.

The Mahler measure `log m(f) := (1/2π)∫log|f(e^{iθ})| dθ` is
estimated via Jensen-FFT at grid M = 4N. The Toeplitz `m_4` follows
the S204 protocol (standardised entries; symmetric Toeplitz; eigenvalues
scaled by `1/√N`).

## Trinity verdicts at parity-only stripping (Q={1, 2})

### F1 — sanity (parity major arc nulled)

`|f̃_N(−1)|` measured by direct dot-product:

| N        | `|f_N(−1)|` | `|f̃_N(−1)|` |
|---------:|------------:|-------------:|
| `2^{12}` | 562.00      | `3.66e−12`   |
| `2^{14}` | 1898.00     | `1.81e−10`   |
| `2^{16}` | 6540.00     | `1.47e−09`   |
| `2^{18}` | 22998.00    | `1.81e−08`   |

**F1 PASS.** `|f̃(−1)|` < `10⁻⁷` (numerical floor).

### F2 — L^∞ profile after stripping

Argmax of `|f̃_N(e^{2πi k/M})|` on M=4N FFT grid: stays at `k=0`
(z=1, the DC value `f_N(1) = π(N)`) for all N. The pre-stated F2
"L^∞ shifts to q ≥ 3" was naive: the trivial DC peak at `q=1` (= `z=1`)
has magnitude `π(N)` and is invariant under parity stripping (since
`Σ_{n=1}^{N} (−1)^n = 0` for even N). Effectively E2.21's content
is the *non-trivial* peak structure (q ≥ 2).

The non-trivial second-largest peak after parity stripping
empirically lands near `e^{2πi/3}` (q=3), magnitude ≈ `π(N)/2`,
matching the HL prediction `μ²(3)/φ(3) · π(N) = (1/2)π(N)`.

### F3 — Mahler deficit fate (parity-only stripping)

Per N, with 50 Bernoulli matched-density replicates each (chi
stripped vs BERN-stripped baseline):

| N        | log m(f_N) | log m(f̃_N) | BERN     | BERN-strip | Δ(χ_P) | Δ(χ̃_P) |
|---------:|-----------:|------------:|---------:|-----------:|-------:|--------:|
| `2^{12}` |    +2.5207 |    +2.4686  | 2.806    | 2.766      | −0.286 | −0.298  |
| `2^{14}` |    +3.1277 |    +3.1137  | 3.427    | 3.416      | −0.299 | −0.303  |
| `2^{16}` |    +3.7462 |    +3.7420  | 4.052    | 4.050      | −0.306 | −0.308  |
| `2^{18}` |    +4.3788 |    +4.3775  | 4.688    | 4.687      | −0.309 | −0.309  |

The deficit `Δ_∞ ≈ −0.307` (S134 / E2.20) is **PRESERVED** after
parity stripping with both `χ` and BERN stripped consistently. At
N=2^{18} the difference `Δ(χ̃) − Δ(χ)` is below the BERN noise floor.

**Naive reading:** "Mahler deficit is residual-sourced (NOT
parity-sourced)." This was the pre-stated F3.b interpretation.
However, the proper measurement requires **also normalising for the
L²-mass change induced by stripping** — see §"Sequential strip"
below. The naive deficit-vs-stripped-BERN measurement holds, but
admits a sharper interpretation.

### F4 — BDJ Toeplitz `m_4` after parity stripping

20 Bernoulli replicates per cell:

| N    | m_4(χ_P) | m_4(χ̃_P) | reduction | BERN       | BDJ universal |
|-----:|---------:|----------:|----------:|-----------:|---------------|
|  500 |   38.35  |    6.871  |   −82.1%  |  2.65±0.28 | 8/3 ≈ 2.67    |
| 1000 |   60.68  |   10.355  |   −82.9%  |  2.61±0.21 | 8/3 ≈ 2.67    |
| 1500 |   77.91  |   13.096  |   −83.2%  |  2.73±0.17 | 8/3 ≈ 2.67    |

**F4 OUTCOME:** parity stripping closes ~83% of the BDJ violation.
`m_4(χ̃_P) ≈ 0.50 · N/log²N` (vs S204's `m_4(χ_P) ≈ 2.95 · N/log²N`).
The constant 0.50 is consistent with cumulative q ≥ 3 major-arc
contributions.

## Sequential major-arc stripping — the structural decomposition

The principal new content is the **hierarchical attribution** of
the three statistics across cumulative squarefree-q strips.

### Mahler shape-deficit (`χ_P` and BERN both stripped, L²-normalised)

The "shape-deficit" `Δ_shape(Q) := (log m(f^{(Q)}_χ) − log ‖f^{(Q)}_χ‖_2)
− (log m(f^{(Q)}_BERN) − log ‖f^{(Q)}_BERN‖_2)` controls for the
L²-mass loss induced by stripping. Without this normalisation the
"raw" deficit drifts mechanically as `‖f^{(Q)}‖_2` decreases.

`N = 2^{14}`, 20 BERN replicates:

| Q    | qs stripped       | log m(χ*) | log m(BERN*) | raw Δ    | shape Δ  | ‖χ*‖₂ |
|-----:|-------------------|----------:|-------------:|---------:|---------:|------:|
|   0  | none              |   3.1277  |   3.4265     | −0.2987  | −0.2987  | 43.59 |
|   1  | {1}               |   3.1136  |   3.4140     | −0.3003  | −0.3003  | 40.98 |
|   2  | {1, 2}            |   3.0996  |   3.4034     | −0.3038  | −0.2337  | 38.21 |
|   3  | {1, 2, 3}         |   3.1168  |   3.4241     | −0.3073  | −0.1981  | 36.74 |
|   5  | {1, 2, 3, 5}      |   3.1129  |   3.4239     | −0.3110  | −0.1812  | 35.99 |
|   7  | {1, 2, 3, 5, 6, 7}|   3.1047  |   3.4233     | −0.3186  | −0.1292  | 33.90 |
|  11  | {... up to 11}    |   3.0965  |   3.4222     | −0.3257  | −0.1022  | 32.75 |
|  13  | {... up to 13}    |   3.0916  |   3.4214     | −0.3298  | −0.0981  | 32.47 |
|  17  | {... up to 17}    |   3.0785  |   3.4192     | −0.3407  | −0.0718  | 31.25 |
|  23  | {... up to 23}    |   3.0593  |   3.4145     | −0.3552  | −0.0551  | 30.24 |
|  29  | {... up to 29}    |   3.0486  |   3.4114     | −0.3628  | −0.0494  | 29.80 |

**The shape-deficit collapses from −0.299 (Q=0) to −0.049 (Q=29) —
83% closure in Q ≤ 29.** Trajectory monotone, asymptotically → 0.
The q=2 parity arc closes only **22%** of the deficit
`(0.299 − 0.234)/0.299`; the remainder is distributed across q ≥ 3.

This **falsifies the implicit S134-attribution of `Δ_∞` to a single
mechanism**: the deficit is sequentially-major-arc-sourced and
*broadly distributed*, not parity-dominated. (Cf. S204 attributing
the BDJ `m_4` violation to parity dominantly — the opposite extreme
of the hierarchy below.)

### BDJ Toeplitz m_4 (sequential)

`N = 1000`:

| Q    | qs stripped       | m_4(χ*) | m_4 / (N/log²N) | λ̃_max(χ*) |
|-----:|-------------------|--------:|----------------:|-----------:|
|   0  | none              | 60.682  |          2.896  |   15.245   |
|   1  | {1} (DC only)     | 60.682  |          2.896  |   15.245   |
|   2  | {1, 2}            | 10.355  |          0.494  |    6.561   |
|   3  | {1, 2, 3}         |  8.614  |          0.411  |    6.991   |
|   5  | up to 5           |  8.976  |          0.428  |    7.199   |
|   7  | up to 7           |  4.208  |          0.201  |    4.487   |
|  11  | up to 11          |  3.459  |          0.165  |    3.669   |
|  13  | up to 13          |  3.524  |          0.168  |    3.710   |
|  17  | up to 17          |  3.079  |          0.147  |    3.655   |
|  23  | up to 23          |  3.089  |          0.147  |    3.967   |
|  29  | up to 29          |  3.237  |          0.154  |    4.156   |

The DC strip alone (Q=1) leaves m_4 unchanged (DC contributes only
to the rank-1 mean shift, not to the spike). Parity strip (Q=2)
closes 83% of m_4 excess; further stripping closes another ~10%.
Asymptotic m_4 / (N/log²N) ≈ 0.15 at Q=29 vs the 0 needed for BDJ
universality. Sub-leading constant declines slowly with Q — pattern
consistent with progressive HL spike removal.

### Hierarchy of major-arc concentration

| Statistic                  | q=2 attribution | q ≤ 7 attribution | q ≤ 29 attribution |
|----------------------------|----------------:|------------------:|-------------------:|
| Newman L^∞ peak (E2.21)    |          100%   |             100%  |              100%  |
| BDJ m_4 (E2.31)            |           83%   |              93%  |               95%  |
| Mahler shape-deficit (E2.20)|           22%   |              57%  |               84%  |

(Percentages = fraction of the "structural deviation from random"
attributable to stripping squarefree major arcs ≤ Q.)

The three statistics, formerly read as "parity major-arc
fingerprints" of E2.20+E2.21+E2.31, are NOT parity-unified.
They form a **hierarchy of major-arc concentration** ordered by
their q=2 attribution fraction. L^∞ peak is parity-pure
(by definition), BDJ m_4 is parity-dominant (~83%), Mahler deficit
is broadly-distributed across many q.

## What edges this work refines

**E2.20 (Mahler deficit `−0.307`):** SHARPENED. The deficit is
sequentially-major-arc-sourced; q=2 contributes only ~22%. Open: is
the Q→∞ limit `Δ_shape(Q) → 0`? Empirical extrapolation suggests
yes; theoretical anchor would be S168's squarefree-q spike formula
(E(q,N) = μ²(q)·(π(N)−r(q))² / (φ(q) N)) summed over all sqf q.

**E2.31 (BDJ m_4 ≈ 2.95 N/log²N):** SHARPENED. Parity contribution
is 83% (matches S204's "89% from largest eigenvalue" within
rank-1-vs-parity-vector overlap). Residual ~17% from q ≥ 3
collectively.

**E2.21 (Newman L^∞):** unchanged content; this work uses the q=2
peak as the primary subtractand to test the other two statistics.

## What this does NOT do

- Does not produce a new EDGES.md entry. The findings refine three
  existing edges; refinement annotations are appended in-place per
  CLAUDE.md "When you discover a new edge" rules (does not meet the
  bar — no new theorem, no new identity not derivable from existing
  HL machinery).
- Does not produce algorithmic content. The hierarchy is descriptive
  / structural, and the major-arc-sourced view is implicit in S168
  (which the project already had). The novelty is empirically
  separating the three statistics by q-attribution profile.
- Does not unify with E1.5 (h_2 entropy of π(X)/X mod m). E1.5 lives
  on a different scale (entropy of small-modulus projections) and
  composes with these via S168 only at q ≤ N^{0.185}.

## Algorithmic implication

None directly. The hierarchy refines the project's internal model
of which arithmetic statistics carry which fraction of HL singular-
series content. Useful for choosing the right baseline when
designing a new pseudorandomness measure: "a measure dominated by
q=2 content (e.g., Newman L^∞ or BDJ m_4) does NOT separate from
random under parity stripping; a measure distributed across many
q (e.g., Mahler deficit) requires multi-arc stripping to test."

## Open questions / next session candidates

1. **Sharpen the Mahler shape-deficit asymptotic.** Predict
   `Δ_shape(Q) ≈ −Σ_{q sqfree, q > Q} c_q · log(...)` analytically
   from a Jensen-integral expansion in HL major arcs.
2. **Test on Liouville / Möbius (E2.20 cross-domain comparison).**
   For `f^λ_N(z) := Σ λ(n) z^n`, does the same hierarchy hold?
   E2.18 says λ Anderson-Lyapunov is BDJ-universal already (no
   spike); the Mahler shape-deficit test is the right diagnostic.
3. **Sequential strip as algorithmic primitive.** Each `V_q^prim`
   projection costs `O(φ(q) N)`; cumulative up to Q=29 cost
   `O(N · Σ φ(q)) = O(N log²Q)`. Negligible for N up to 2^20.
   Possible application: pre-conditioner for spectral approximation
   of `χ_P`-related operators where the parity peak dominates
   condition number.

## File manifest

- `definition.md` — formal construction + pre-stated falsifiers.
- `parity_stripped_trinity.py` — main FFT + BDJ measurement script.
- `parity_stripped_trinity_results.json` — raw FFT + BDJ block results.
- `sequential_strip_check.py` — multi-Q Mahler stripping with
   matched-Q BERN baseline.
- `sequential_strip_results.json` — raw sequential Mahler results.
- `sequential_bdj_check.py` — multi-Q BDJ m_4 stripping.
- `sequential_bdj_results.json` — raw sequential BDJ m_4 results.
- `parity_stripped_trinity_results.md` — this file.
- `run.log`, `sequential_strip.log`, `sequential_bdj.log` — stdouts.

# Session 239 — parity-stripped trinity (PARADIGM-SHIFT mode)

**Date:** 2026-04-30.
**Mode:** PARADIGM-SHIFT (no cross-domain imports permitted).
**Self-grade:** **B**.
**Edges refined:** **E2.20** (Mahler deficit `−0.307`),
**E2.31** (BDJ Toeplitz `m_4 ≈ 2.95 N/log²N`), **E2.21** (Newman L^∞
parity peak — used as construction primitive). **No new EDGES entry.**
**Predecessors:** S134 (E2.20 build), S138 (E2.21 build), S204 (E2.31
build), S55 / S70 (E1.6 bisection), S168/S190 (squarefree-q spike
formula).

## Mode discipline

* **No `WebFetch` / `WebSearch`.** Confirmed.
* **No new cross-domain technique.** The construction (orthogonal
  projection along `(−1)^n`; sequential additive-Fourier projection
  onto `V_q^prim` for sqfree q) is elementary linear algebra, already
  USED in the project (Möbius inversion of `[gcd(n,W)=1]` from S219;
  additive-Fourier subspace decomposition from S168). Jensen-FFT for
  Mahler measure and Toeplitz-spectrum framework are USED-mode
  techniques in CROSS_DOMAIN_TECHNIQUES.md.
* **No new `ATTACK_VECTORS.md` entries.**
* **Only existing project content read** (EDGES.md, CLOSED_PATHS.md,
  selected `experiments/constructions/`, recent paradigm-shift sessions
  S191/S200/S219, OPEN_POSITIVE_TARGETS.md).

## Target picked (paradigm-shift composition)

Three "parity-fingerprint" edges (E2.20 Mahler deficit `−0.307`,
E2.21 Newman L^∞ peak at `z=−1`, E2.31 BDJ `m_4 ≈ 2.95 N/log²N`)
have been read implicitly in the project as a **unified
parity-major-arc fingerprint** of `χ_P`. S204 explicitly attributes
~89% of E2.31 to the parity rank-1 spike. E2.21 is parity-defined.
E2.20's mechanism is "TBD; HL singular series the most plausible
candidate" (S134 open question).

The paradigm-shift composition: **subtract the parity major arc from
`χ_P` and re-measure all three statistics simultaneously**. If the
parity attribution of E2.20 is similar to E2.31's, all three should
collapse to baseline together. If different, the three edges are
NOT structurally unified by q=2 alone.

## What I built

`experiments/constructions/parity_stripped_trinity/`:
* `definition.md` — formal construction + 4 pre-stated falsifiers.
* `parity_stripped_trinity.py` — main Mahler + L^∞ + Toeplitz `m_4`
  measurement with matched-density Bernoulli baselines.
* `parity_stripped_trinity_results.md` — F1-F4 verdict tables +
  hierarchy decomposition.
* `parity_stripped_trinity_results.json` — raw FFT+BDJ measurements.
* `sequential_strip_check.py` — extension to cumulative Q-stripping
  for Mahler shape-deficit, with matched-Q BERN baseline.
* `sequential_bdj_check.py` — extension to cumulative Q-stripping
  for BDJ Toeplitz `m_4`.
* `sequential_strip_results.json`, `sequential_bdj_results.json`,
  `run.log`, `sequential_strip.log`, `sequential_bdj.log`.

## The three findings

### Finding 1 — Parity attribution of the trinity is hierarchical, NOT uniform

`χ_P` parity correlation `α_2(N) = (1/N)Σ χ_P(n)(−1)^n = −(π(N)−2)/N`
defines `χ̃_P := χ_P − α_2·(−1)^n`. Then `f̃_N(−1) = 0` exactly.
Single-strip Mahler / L^∞ / BDJ measurements at `N ∈ {2^{12}, 2^{14},
2^{16}, 2^{18}}` (FFT block) and `N ∈ {500, 1000, 1500}` (BDJ block):

| Statistic                  | q=2 attribution |
|----------------------------|----------------:|
| Newman L^∞ peak (E2.21)    |          100%   |
| BDJ m_4 (E2.31)            |           83%   |
| Mahler shape-deficit (E2.20)|           22%   |

The Mahler deficit at N = 2^{18} is **identical** in `χ_P` and `χ̃_P`
(both at `−0.309 nat` against matched-stripped-BERN baseline) — the
parity strip removes `< 1%` of the deficit at this N. By contrast,
BDJ `m_4` drops from 60.7 → 10.4 (83% reduction) at N=1000.

### Finding 2 — Sequential cumulative stripping of squarefree-conductor major arcs

At `N = 2^{14}`, sequentially strip `V_q^prim` for sqfree q ≤ Q
(primitive additive-Fourier characters mod q) from BOTH `χ_P` AND
the BERN baseline. Mahler shape-deficit (L²-mass-normalised):

| Q  | shape-Δ | closure of `Δ_∞` |
|---:|--------:|-----------------:|
|  0 | −0.299  |        0%        |
|  2 | −0.234  |       22%        |
|  7 | −0.129  |       57%        |
| 13 | −0.098  |       67%        |
| 23 | −0.055  |       82%        |
| 29 | −0.049  |       84%        |

Trajectory monotone, asymptotically heading toward zero. The Mahler
deficit `−0.307` is **distributed across many sqfree-conductor major
arcs**, NOT concentrated at q=2.

BDJ `m_4` sequential strip at `N = 1000`:

| Q  | m_4(χ*) | m_4 / (N/log²N) |
|---:|--------:|----------------:|
|  0 | 60.682  |          2.896  |
|  2 | 10.355  |          0.494  |
|  7 |  4.208  |          0.201  |
| 11 |  3.459  |          0.165  |
| 29 |  3.237  |          0.154  |

(Compare BDJ universal `8/3 ≈ 2.67` and BERN matched-density `m_4 ≈
2.65`.) Q=2 strip closes 83% of the violation; Q ≤ 29 closes 95%.

### Finding 3 — The E2.20+E2.21+E2.31 trinity is NOT parity-unified

The naive "all three are parity-major-arc fingerprints" reading is
**falsified**. The three statistics span a hierarchy of major-arc
concentration:

* L^∞ peak: parity-pure (definitional).
* BDJ m_4: parity-dominant (83%); residual ~17% from q ≥ 3.
* Mahler shape-deficit: parity-minor (22%); ~78% from q ≥ 3
  collectively, with a slow geometric tail.

Pre-stated F3.b ("Mahler deficit residual-sourced") is the actual
outcome, but the stronger claim is the **inversion** of parity
attribution between E2.20 and E2.31: 22% vs 83% — opposite extremes
of the hierarchy.

## Mechanism interpretation (without new technique import)

The Toeplitz `m_4` measures `‖f_N‖_4^4 / N²` (L^4 norm of `f_N`
on the unit circle, modulo Szegő). Major-arc contributions to L^4
add as `(μ²(q)/φ(q))² · π(N)²` per sqfree q, summing as Σ μ²(q)/φ(q)²
(Mertens / Selberg-Delange type series, dominated by small q —
`q=2` weight `1`, `q=3` weight `1/4`, `q=5` weight `1/16`, ...
geometric tail). **Parity (q=2) dominates ≈ 83%.**

The Mahler measure is the GEOMETRIC mean of `|f_N(e^{iθ})|`, not a
power-sum. Major-arc contributions add **logarithmically**: each q
contributes `log(μ²(q)/φ(q))` weighted by the local measure of its
arc-neighborhood. Local measures of arc neighborhoods at q scale as
`1/q` (φ(q) primitive characters × `O(1/(qN))` width each). The
weighted log-sum is **distributed evenly** across many q because
the log-flat geometry of Mahler integration suppresses peaks. The
22% q=2 vs 78% rest split is the natural log-vs-power-sum behavior
on the same HL-major-arc data.

This is a *structural* explanation (no external import). The exact
22% vs 83% inversion is empirically verified, not derived; the
predicted asymptotic at large Q would require analytical continuation
of the Jensen-integral expansion in HL major arcs (open question).

## Pre-stated F1-F4 verdicts

| F# | Falsifier                                          | Verdict | Outcome |
|----|----------------------------------------------------|---------|---------|
| F1 | `|f̃_N(−1)| ≤ 10⁻¹⁰` for parity strip             | PASS    | Numerical floor `< 10⁻⁷`. |
| F2 | L^∞ shifts to q ≥ 3 after stripping               | PARTIAL | DC peak (q=1) supersedes; non-trivial peak at q=3 with `π(N)/2` magnitude (HL prediction). |
| F3 | F3.b (residual-sourced Mahler) vs F3.a / F3.c     | F3.b    | Sharper: NOT residual-sourced — *sequentially major-arc-sourced*, broadly distributed. |
| F4 | F4.a (BDJ → 8/3) vs F4.b (residual divergence)    | F4.b    | 83% closure at Q=2; residual ~17% from q ≥ 3 collectively. |

## What edges this work refines

- **E2.20** (Mahler deficit): *S239 sequential-strip refinement
  inline*. The deficit is sequentially major-arc-sourced; q=2
  contributes 22%, not the implied "dominant" share. The mechanism
  candidate "HL major arcs" (S134 open question) is now **empirically
  supported** at `Q ≤ 29` closure of 84% of the deficit.

- **E2.31** (BDJ m_4): *S239 sequential-strip refinement inline*.
  Confirms S204's parity attribution (~83%); residual ~17% comes
  from q ≥ 3 collectively, with `m_4 / (N/log²N) ≈ 0.15` at Q=29
  (consistent with geometric-tail HL contribution).

- **E2.21** (Newman L^∞): unchanged content; serves as the primary
  subtractand.

## What this does NOT do

- Does not produce algorithmic content for π(x). The hierarchy is
  descriptive / structural; no new computational primitive.
- Does not write a new EDGES entry. The findings refine three
  existing edges in place per CLAUDE.md "When you discover a new
  edge" criteria — no new theorem / identity not derivable from
  existing HL machinery.
- Does not contradict S168 / S190's squarefree-q spike formula. The
  S168 formula `E(q, N) = μ²(q)·(π(N)−r(q))² / (φ(q) N)` predicts
  the L²-energy of `χ_P` projected on `V_q^prim`. My sequential
  strip measures the residual `χ_P − P_{V_q^prim}χ_P`'s Mahler
  measure, which is a different statistic. The 22% vs 83% inversion
  is consistent with S168 because Mahler is a log-flat integral
  while m_4 is an L^4 power-sum.

## Self-evaluation (CLAUDE.md template)

1. **What did I produce that was not in the project before?**
   The hierarchy decomposition table (Newman 100% / BDJ 83% /
   Mahler 22% q=2 attribution) and the sequential-strip Mahler
   shape-deficit closure trajectory (`Δ_shape(Q) → 0` empirically
   monotone). The construction (`χ̃_P` parity-strip + cumulative
   `V_q^prim` strip) and the inversion observation (22% vs 83%)
   are new project artifacts.

2. **What edges did my work compose or cite?**
   E2.20, E2.21, E2.31 (composed); E1.6 + E2.10 (parity bisection
   identity used to derive `α_2 = −(π(N)−2)/N`); S168 / S190
   (squarefree-q spike formula, cited as the structural anchor for
   the residual q ≥ 3 contributions); E7.21 (Pólya-Carlson — flagged
   as alternative mechanism candidate but not invoked).

3. **Did my session produce only duplicate closures?**
   No. The hierarchy decomposition and the 22%/83% inversion were
   not previously stated in the project. The S239 finding falsifies
   an implicit (but never explicit) "all three are parity-fingerprint"
   reading and produces a new structural decomposition.

4. **Next-action for the next agent.**
   Compute the asymptotic of `Δ_shape(Q→∞)` analytically: predict
   the per-q contribution as a Jensen-integral over the q-th
   major-arc neighborhood. If the prediction matches the empirical
   `0.299 → 0.049` trajectory across Q ∈ {0..29}, E2.20 reduces
   structurally to E2.13 (HL singular series) and the deficit is no
   longer an independent edge — it's a Jensen-Mahler corollary of
   the HL imprint. Filed in `parity_stripped_trinity_results.md` as
   "Open question 1".

## Grade defense

**B** is appropriate because:
- The session produced new structural content (hierarchy
  decomposition + parity-attribution inversion) not previously stated.
- It refines TWO existing edges (E2.20, E2.31) inline with concrete
  empirical evidence.
- It does NOT produce a new theorem, new identity, or algorithmic
  content (those are the A-grade thresholds).
- The construction is elementary (orthogonal projection); novelty
  is in the EXPERIMENT, not the technique.
- Pre-stated falsifiers F3 and F4 produced sharper outcomes than
  predicted (F3.b sharpened to "sequentially major-arc-sourced";
  F4.b sharpened with quantitative residual constant 0.15).

Inflated A would require: a closed-form prediction for the per-q
Mahler contribution matching empirical 22%/+57%/+27%/+/+15% etc.
This was not produced; left as the next-session open question.

# Session 124 — §D2.a.1 PH on continuous-marginal-matched baseline B3

**Mode:** novelty (single-session, B-grade target).
**Self-grade: B** (refinement of E2.17 with three-way decomposition;
pre-stated F_a holds on absolute thresholds; partial violation on
relative thresholds is the new content).

## Target

`NOVELTY_CHALLENGES.md` §D2.a.1 (proposed S117): "PH on the residual
marginal-distribution component. Construct a continuous baseline
B3 = IID samples from a KDE-fitted continuous distribution matched
to the W-tricked marginal. z(P_W; B3) should land near zero on T0
AND T1, isolating the marginal-distribution effect from any residual
serial structure. Predicted: B3 absorbs the |z| ≈ 7-12 vs B1
deficit entirely. Cost: 1 session."

## What this session produced

### One-line outcome

Continuous-quantile baseline B3 (inverse-transform sampling on the
linearly-interpolated empirical CDF of the W-tricked Cramér-
normalised gaps; Devroye 1986 §II.2.1) yields |z(P_W; B3)| ≤ 0.65
on T0 AND T1 across d ∈ {2, 3, 4} — the predicted full absorption
holds. The S117 marginal-component decomposes further into a
**marginal-envelope** sub-component (dominant) and a **discreteness**
sub-component (small, ~1-3σ on T0).

### New mathematical content not in the project before

1. **B3 baseline construction**: a probabilistic null-model primitive
   that has not been used in the project before — IID inverse-
   transform sampling on the linearly-interpolated empirical CDF.
   This complements the existing B1 (parametric IID Exp(1)) and B2
   (permutation = exact discrete marginal) baselines. CROSS_DOMAIN
   import: Devroye 1986 §II.2.1.

2. **Quantitative three-way decomposition of E2.17 PH deficit**:
   refines S117's two-way (serial + marginal) into three components
   (envelope + discreteness + serial-residual), with envelope
   identified as dominant by ~7-9σ on T0. The discreteness component
   is the smallest novel finding — a real but small (~1-3σ on T0)
   sub-effect that the existing B2 baseline cannot disentangle from
   serial residual.

3. **Cancellation observation**: the discreteness and serial-residual
   sub-components partially null on (PRIMES_W vs B3), explaining why
   z(P_W; B3) ≈ 0 even though both individual sub-components measure
   1-3σ. This is a precise structural observation about how the
   discrete-grid PH lattice effect (B2 > B3 by ~5 on T0) combines
   with serial correlation (PRIMES_W < B2 by ~5 on T0).

4. **EDGES.md E2.17 updated inline** with the three-component
   decomposition statement.

### Edges composed / cited

- **Composes** E2.13 (Gowers W=210 uniformity) + E2.17 (PH deficit)
  + Devroye 1986 inverse-transform sampling as a structural
  refinement of the singular-series fingerprint family.
- **Refines E2.17** inline (no new edge — sub-decomposition of an
  existing component).
- **Cross-references** E2.14, E2.15, E2.16, E2.20 — sister
  W-trick-erases-everything HL fingerprints.

### Files changed / added

```
experiments/topological/persistent_homology_w_trick_marginal_b3/
  persistent_homology_w_trick_marginal_b3.py
  persistent_homology_w_trick_marginal_b3_results.md
  b3_pilot_b1.json, b3_main.json, b3_d2.json, b3_d4.json
  b3_pilot_b1.log, b3_main.log

EDGES.md                                  — refined E2.17
NOVELTY_CHALLENGES.md                     — §D2.a.1 marked CLOSED;
                                            §D2.a.1.i and .ii successors
status/CLOSED_PATHS.md                    — new row at session 124
status/SESSION_INSIGHTS.md                — new entry session 124
archive/sessions/session124_d2a1_ph_marginal_b3.md  (this file)
```

## Cross-domain technique

**Devroye 1986** "Non-Uniform Random Variate Generation" Springer,
Chapter II §II.2.1 (Inversion Method). Standard text on probabilistic
sampling. Used to construct the linearly-interpolated empirical CDF
sampler. Not augmenting CROSS_DOMAIN_TECHNIQUES.md (the technique is
elementary probability, not a frontier import).

## CLAUDE.md self-evaluation (4 questions)

### Q1. What did I produce that was not in the project before this session?

(a) The B3 baseline (continuous inverse-transform sampler from the
empirical CDF of W-tricked Cramér-normalised gaps) and its PH
characterisation across d ∈ {2, 3, 4}, M=1000, x ≈ 10⁶, three
residues b ∈ {1, 11, 13} pooled, K=20.

(b) A three-way decomposition of E2.17 that refines S117's two-way
decomposition: (envelope ~7-9σ) + (discreteness ~1-3σ) + (serial-
residual ~1-2σ) on T0. The cancellation of (ii) and (iii) on the
(PRIMES_W vs B3) comparison is a new structural observation.

(c) The EDGES.md E2.17 entry has been refined inline with the
three-component statement and the new B3 z-scores. The mathematical
content was not in the project before this session.

(d) Empirical observation that B3 has uniformly larger std (~8 on T0
for d=3) than B2 (~3) — a structural consequence of continuous
resampling vs discrete permutation, not previously documented.

### Q2. What edges did my work compose or cite?

- Composes E2.13 (Gowers W=210 W-trick) + E2.17 (PH deficit) +
  Devroye 1986 inverse-transform sampling.
- Refines E2.17 inline.
- Cross-references E2.14 / E2.15 / E2.16 / E2.20 as sister
  fingerprints in the W-trick-erases-everything HL family.

### Q3. If my session produced only duplicate closures, why?

N/A — the session produced refined structural content (the three-way
decomposition with quantitative attributions of each sub-component)
that did not exist in the project before. The B3 baseline is novel
to this project. The pre-stated F_a holds on absolute thresholds;
the partial violation on relative thresholds (T0 d=2 Δ=2.89, d=3
Δ=1.94) IS the new content (it surfaces the discreteness sub-
component as a real but small effect).

### Q4. What is the next-action for the next agent?

The two successors are written in NOVELTY_CHALLENGES §D2.a.1.i and
§D2.a.1.ii:

- **§D2.a.1.i** Pure-discrete IID baseline B4 (sample from empirical
  PMF, no interpolation) — isolates discreteness sub-component
  directly by comparing B4 (discrete IID) vs B2 (discrete
  permutation); identical marginal but different serial structure.
  Predicted: z(P_W; B4) within 0.5σ of z(P_W; B2) if S117/S124
  serial-residual is truly noise. Cost: 1 session.

- **§D2.a.1.ii** Sliding-bandwidth KDE B5(σ) — vary σ from 0
  (≈ B2 discrete) to large (≈ near-Exp(1) smooth) and trace
  z(P_W; B5(σ)). Predicted: sigmoidal crossover at σ ≈ 0.16
  (grid-spacing / 2) where the discreteness component switches
  off. Would produce a "PH bandwidth-scan curve" parallel to
  E2.13's Gowers W-scan. Cost: 1 session.

Either is a B-grade single-session refinement.

## Honest failure assessment

The session did NOT produce A-grade content. The three-way
decomposition is a structural refinement of an existing edge; the
B3 absorption was the predicted outcome from S117. The unexpected
content was the partial relative-threshold violation (B3 absorbs
even better than B2 on T0), which led to identifying the
discreteness sub-component — a small but consistent ~1-3σ effect.
This is B-grade: substantive refinement, not novelty leap.

The discreteness sub-component is at the edge of statistical
significance (1.7-4.8σ across (d, b) cells on T0). The successor
challenge §D2.a.1.i (B4 baseline) is required to confirm
discreteness as its own distinct component beyond the cancellation
with serial-residual.

No A-grade upgrade was warranted: the session followed the planned
B-grade target, which produced a clean refinement.

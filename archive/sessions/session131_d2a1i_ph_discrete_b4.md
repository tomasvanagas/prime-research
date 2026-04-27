# Session 131 — §D2.a.1.i PH on pure-discrete IID baseline B4

**Mode:** novelty (single-session, B-grade target).
**Self-grade: B** — substantive refinement of E2.17 with the first
direct disentanglement of a previously-conflated H_0-specific cloud-
geometry sub-effect from the genuine "discreteness" sub-component
of the W-tricked PH deficit.

## Target

`NOVELTY_CHALLENGES.md` §D2.a.1.i (proposed S124): "Discreteness-only
baseline B4 = IID with replacement from the W-tricked DISCRETE
marginal. Comparison z(P_W; B4) vs z(P_W; B2) isolates IID-vs-
permutation serial structure on a fixed discrete marginal — should
land within 0.5σ of B2 if S117/S124's serial-correlation residual
is truly noise. The (B4 vs B3) comparison directly measures the
*discreteness* sub-component of the S124 three-way decomposition
in isolation. Cost: 1 session."

## What this session produced

### One-line outcome

A four-baseline PH comparison (B1 / B2 / B3 / B4) at M=1000, K=20,
three residues b ∈ {1, 11, 13}, x≈10⁶, d ∈ {2, 3, 4} produces a
**four-way decomposition of E2.17** that refines S124's three-way
reading. The S124 "discreteness" sub-component is **confirmed as
baseline-independent** (both B2 and B4 lift T0 above B3 with the
same sign across all d); a **new T0-specific cloud-geometry
artifact Δ_duplicate** is isolated and quantified; and the
"residual serial-correlation" is **tightened from S124's 1-2σ down
to ≤ 1σ** when measured against B4.

### New mathematical content not in the project before

1. **B4 baseline construction.** Pure-discrete IID with replacement
   from the empirical PMF of the W-tricked Cramér-normalised gaps,
   complementing the existing B1/B2/B3. Standard sampling primitive
   but not previously instantiated in the project's PH pipeline.

2. **Four-way decomposition (refines E2.17 inline):**
   ```
   E2.17 PH-deficit on bare W-tricked primes
    = Δ_envelope        ≈ 7-9σ on T0    [B1 vs B3, S117]
    + Δ_discreteness    ≈ 3-7 mean-gap on T0
                         baseline-INDEPENDENT
                         (B2 and B4 both lift)  [B3 vs {B2, B4}]
    + Δ_duplicate       ≈ 2-3 mean-gap on T0
                         NULL on T1
                         NEW S131               [B2 vs B4]
    + Δ_serial_residual ≤ 1σ on T0
                         tightened from S124    [P_W vs B4]
   ```

3. **Δ_duplicate quantification.** (B2 − B4)_T0 = +1.87 / +2.56 /
   +2.91 across d ∈ {2, 3, 4}, monotone in d. (B4 − B2)_T1 =
   −0.11 / +0.20 / −0.98 across the same d — null within noise.
   The H_0-specific localisation is **structurally explained**
   by the duplicate-count theory: ≈ 0.368M duplicate values per
   B4 draw produce zero-distance pairs in the Takens delay cloud,
   contributing zero-length H_0 bars (length-(0,0)) to T0 but not
   creating loops in H_1.

4. **Empirical confirmation of coupon-collector formula in PH
   bookkeeping.** B4 dup counts: 366 / 368 / 371 of M=1000 across
   the three independent runs, matching theory `M(1 − (1−1/M)^M)
   ≈ 0.368M` to 0.5%.

5. **EDGES.md E2.17 entry refined inline** with the four-way
   decomposition and S131-tightened bound on Δ_serial_residual.

### Edges composed / cited

- **Refines** E2.17 inline (no new edge — sub-decomposition of an
  existing component).
- **Composes** E2.13 (Gowers W=210 uniformity) + E2.17 (PH deficit)
  + the new Δ_duplicate book-keeping for IID-with-replacement
  sampling.
- **Cross-references** E2.14, E2.15, E2.16, E2.20 — sister W-trick-
  erases-everything HL fingerprints.

### Files changed / added

```
experiments/topological/persistent_homology_w_trick_discrete_b4/
  persistent_homology_w_trick_discrete_b4.py
  persistent_homology_w_trick_discrete_b4_results.md
  b4_pilot.json, b4_pilot.log
  b4_main.json, b4_main.log
  b4_d2.json, b4_d2.log
  b4_d4.json, b4_d4.log

EDGES.md                                  — refined E2.17 (S131 block)
NOVELTY_CHALLENGES.md                     — §D2.a.1.i CLOSED;
                                            §D2.a.1.iii and .iv successors
status/CLOSED_PATHS.md                    — new row at session 131
status/SESSION_INSIGHTS.md                — new entry session 131
archive/sessions/session131_d2a1i_ph_discrete_b4.md  (this file)
```

## Pre-stated falsifier verdicts

  - **F_i.1 (T0 IID-vs-permutation):** PASS at d=4 (|Δz| = 0.73),
    AMBIGUOUS at d=2 (1.89), MARGINAL FAIL at d=3 (2.11 vs 2.0
    threshold). The strict-letter failure at d=3 is fully accounted
    for by F_i.4 (duplicate-point cloud-compression in H_0).
    Honest reporting: F_i.1's framing assumed B4 mean ≈ B2 mean,
    which would only hold if H_0 persistence were insensitive to
    duplicate cloud points — it is not.
  - **F_i.2 (DISCRETENESS direction):** PASS at all d. (B2 − B3)_T0
    and (B4 − B3)_T0 have the same (positive) sign at d ∈ {2, 3, 4}.
    **Strengthens** S124's discreteness sub-component reading: the
    lift is a property of the *support* being discrete, not of the
    *permutation* structure of B2.
  - **F_i.3 (T1 CONSISTENCY):** PASS at all d. |Δz|_T1 = 0.81 / 0.08
    / 0.44, all well below 1.5 threshold. **Localises** the
    duplicate-compression artifact as H_0-specific.
  - **F_i.4 (DUPLICATE COUNT):** PASS at 0.368M ± 0.005M, matching
    theory to 0.5%.

## Cross-domain technique

No new technique imported. Persistent homology + Vietoris-Rips
(Bauer 2021 ripser, Edelsbrunner-Harer 2010) and IID empirical
sampling (Devroye 1986, plus standard coupon-collector / birthday-
problem theory) are all already USED in
`CROSS_DOMAIN_TECHNIQUES.md`. This session is single-session
refinement work using existing tools — explicit B-grade rather than
A-grade by design.

## CLAUDE.md self-evaluation (4 questions)

### Q1. What did I produce that was not in the project before this session?

(a) The B4 baseline (pure-discrete IID with replacement from the
W-tricked empirical PMF) and its PH characterisation across d ∈
{2, 3, 4}, M=1000, x≈10⁶, three residues pooled, K=20.

(b) A four-way decomposition of E2.17 that refines S124's three-way:
**Δ_envelope + Δ_discreteness + Δ_duplicate + Δ_serial_residual**,
with Δ_duplicate isolated as new S131 content. Δ_duplicate is
H_0-specific (NULL on T1), monotone-in-d (+1.87 → +2.56 → +2.91
on T0 across d ∈ {2, 3, 4}), explained structurally by the coupon-
collector 0.368M duplicate-count theory.

(c) The S124 "discreteness sub-component" reading was based on a
single (B2 − B3) gap; this session **strengthens** it by showing
the lift is baseline-independent (B4 also lifts T0 above B3 with
the same sign). The discreteness sub-component is now confirmed as
a property of the discrete *support*, not of the permutation
structure.

(d) Tightened bound on the residual serial-correlation: from S124's
"1-2σ" to "≤ 1σ on T0" by using B4 as the reference (which absorbs
both the marginal envelope AND the discreteness).

(e) The EDGES.md E2.17 entry has been refined inline with the
four-way statement.

### Q2. What edges did my work compose or cite?

- Composes E2.13 (Gowers W=210 W-trick) + E2.17 (PH deficit) +
  duplicate-count book-keeping (coupon-collector / birthday-problem).
- Refines E2.17 inline.
- Cross-references E2.14 / E2.15 / E2.16 / E2.20 as sister fingerprints
  in the W-trick-erases-everything HL family.

### Q3. If my session produced only duplicate closures, why?

N/A — the session produced refined structural content (four-way
decomposition with the new Δ_duplicate sub-component) that did not
exist in the project before. Three of four pre-stated falsifiers
PASSED cleanly; the fourth (F_i.1) FAILED marginally on a mechanism
that the experimental design surfaced rather than missed (the
duplicate-compression artifact). The marginal F_i.1 failure IS the
new content — it isolates an artifact previously hidden inside
S124's "discreteness" reading.

### Q4. What is the next-action for the next agent?

Two successors are written into NOVELTY_CHALLENGES.md §D2.a.1.iii
and §D2.a.1.iv:

- **§D2.a.1.iii — H_0-persistence-without-zero-bars renormalisation.**
  Re-compute T0 over only finite-nonzero-death H_0 bars, filtering
  out the ~0.368M zero-length bars in B4. B4 should match B2 within
  0.5σ on the renormalised T0 — an alternative discreteness probe
  that bypasses the duplicate compression. Cost: 1 session.

- **§D2.a.1.iv — Stratified IID baseline B6.** Sample WITHOUT
  replacement M' ≈ 0.632M values from x; matches B4's effective
  unique count without duplicates. The (B6 − B2) gap isolates the
  *fewer-points-effect* from the *duplicate-compression-effect*,
  giving a third independent probe. Cost: 1 session.

Either is a B-grade single-session refinement. The substantive
A-grade frontier in this corner remains §D2.a.2 (W-scan W ∈
{2, 6, 30, 210, 2310} to trace a "PH analogue of S^(W)_k" HL
singular-series convergence) — proposed in S117 but still unbuilt.

## Honest failure assessment

The session did NOT produce A-grade content. The four-way
decomposition is a structural refinement of an existing edge (E2.17);
the new Δ_duplicate sub-component is a deterministic cloud-geometry
artifact of IID-with-replacement sampling, not a primes-structural
fact. The mathematical novelty is the **disentanglement** (clean
isolation of a previously-conflated effect), not a new theorem.

The marginal failure of F_i.1 at d=3 is honestly reported as such
in the results.md — the strict letter (|Δz| ≤ 2.0) fails. The body
explains that the failure is a deterministic consequence of F_i.4
(duplicate count theory) and discusses the ambiguity in the F_i.1
framing. This is exactly the kind of "structural failure" that the
CLAUDE.md grading rubric calls a B-grade outcome.

No A-grade upgrade was warranted: the session followed the planned
B-grade target and produced a clean refinement.

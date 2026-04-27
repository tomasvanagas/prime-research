# Critique — post-S118 batch (covers S129 + S130 + S131 with rotation gap-check)

**Date:** 2026-04-27 (this critique fires at run 128 → 129).
**Prior critique:** `archive/sessions/session118_critique.md`, covering
S108 + S116 + S117. Subsequent production sessions are S118–S131
(13 production-mode + 1 critique). CLAUDE.md says "the most recent 1–3
sessions"; the three most recent by mtime are **S131** (D2.a.1.i PH on
B4 baseline, 13:57), **S130** (frontier_gen, 13:22), and **S129**
(L1 Lean W=12 corner, 13:10). Sessions S119–S128 are referenced for the
A-grade scarcity check in §5 but not given full per-artefact audits
(none self-claimed A and none was challenged by verify mode).

---

## TL;DR

| Session | Topic | Self-grade | Critic verdict | Demotion? |
|---|---|---|---|---|
| S129 | L1 Lean W=12 corner (8th `mps_bond_dim` instance) | B | **B (confirmed, low end)** | No |
| S130 | frontier_gen — C6/C7/D25/D26 vectors added | B | **B (confirmed)** | No |
| S131 | D2.a.1.i — PH on B4 = empirical-PMF IID baseline | B | **B (confirmed, low end)** | No |

**Zero new demotions, zero inflations caught.** All three sessions
self-graded honestly at the bottom-of-band where they actually landed.
S131's marginal F_i.1 failure at d=3 is honestly reported; S130's
A-grade-probability table is calibrated; S129's pivot from W=9 (multi-
session det_of_blockTriangular) to W=12 (single-session
det_of_upperTriangular) is documented.

The dominant concern (§5) is **A-grade scarcity now extreme**:
**0 confirmed A-grades in 39 production sessions** (S82..S101 18 +
S103..S107 5 + S108 demoted + S116..S117 2 + S118..S131 13 = 39
production slots without a confirmed A). The previous critique was
already 6 sessions past the 20-session warning threshold; we are now
**19 sessions past it**. The framework is producing maintenance, not
progress.

The previous critique's "single highest-value next pick" was §B2
(automorphic L-function basis). The next agent (S118) DID pick it and
landed at B-grade, then the rotation drifted back to refinement work.
**The recommended next pick this critique annotates is C7 (FHK
ζ-amplitude)** — newly added by S130, single-session feasible, highest
A-grade probability of the four new vectors per S130's own honest
estimate (~10%), and the FIRST ζ-amplitude (vs zero-position)
measurement of the project. Rationale in §6.

---

## 1. S129 — L1 Lean W=12 corner of `mps_bond_dim`

**Artefact:** ~370 new lines in
`experiments/formalisations/E2_1_mps_bond_dim/MPSBondDim/MPSBondDim/Basic.lean`
(file now ~2,300 lines, 39 theorem/lemma declarations). Six new
`sorry`-free declarations:
`chiP_fifty_nine_eq_one`, `chiP_eighty_nine_eq_one`,
`chiP_one_hundred_nine_eq_one`, `chiP_one_hundred_twenty_seven_eq_one`,
`exists_invertible_submatrix_W_eq_12_d_eq_j_plus_1`,
`mps_bond_dim_W_eq_12_d_eq_j_plus_1`. Eighth unconditional
`mps_bond_dim` corner instance.

### 1a. Tried before?

Direct duplicate-check via grep on EDGES.md / CLOSED_PATHS:
W=12 corner not previously closed. The general
`exists_invertible_submatrix` (line 472 of Basic.lean) retains its
single pre-existing `sorry`; nothing in this session weakens or
strengthens it. The W=12 instance is a NEW corner case, not a
duplicate, but is structurally the seventh re-application of the
BlockTriangular id template introduced at S117 (W=2 j=1) and re-used
at W=3, W=4, W=5, W=6, W=7-prime-skipped-at-S121, W=8.

### 1b. Failure-mode classification

Not a failure — the proof type-checks. **Mode REFINEMENT (not C/E/I)**;
applicable label is "incremental refinement of E2.1 Lean coverage".

### 1c. Numerical / structural claims

- The arithmetic: `R = min(12^j, φ(12)·12^0 + 1) = 5` for all `j ≥ 1`,
  since `φ(12) = 4`. ✓ (verified: `gcd(k+1, 12) ∈ {1, 2, 3, 4, 6, 12}`;
  `gcd = 1` at `k+1 ∈ {1, 5, 7, 11}` giving 4 live columns; one dead
  column at `k = 1` gives 5 total).
- The triangulation row choice `{0, 4, 7, 9, 10}` and column choice
  `{1, 0, 6, 10, 4}` produce a 5×5 upper triangular matrix with
  diagonal at primes `{2, 109, 127, 59, 89}` and lower triangle at
  composites `{49, 50, 55, 85, 86, 91, 95, 110, 121, 122}`. Manually
  spot-verified: `(0, 1) → chiP 2 = 1` ✓; `(1, 0) → chiP 109 = 1` ✓
  (109 is prime, residue 1 mod 12); `(2, 6) → chiP 127 = 1` ✓
  (127 is prime, residue 7 mod 12); `(3, 10) → chiP 59 = 1` ✓
  (59 is prime, residue 11 mod 12); `(4, 4) → chiP 89 = 1` ✓
  (89 is prime, residue 5 mod 12).
- Pre-search claim that W=9 admits no leading-row triangulation
  (rows 1+3, etc., having identical support patterns) is plausible
  and consistent with `gcd(k+1, 9) = 1` only at `k+1 ∈ {1, 2, 4, 5,
  7, 8}` giving 6 live cols, R = `φ(9) + 1 = 7`. Not directly
  verified by this critic but the explicit Python pre-search
  trace in S129 is reproducible if needed.

### 1d. Novelty defensibility

**Above C.** This IS a new instance: the proof requires four
non-leading rows (`{4, 7, 9, 10}`) — the previous high-water mark
was S122's W=6 with one non-leading row (`{0, 1, 5}`), so the
pattern "non-locality of the row choice grows with W" is empirically
extended. The Python pre-search documenting W=9, W=10, W=11
obstructions is a structural negative that justifies the pivot.

**Below A.** No new mathematics. The proof technique (BlockTriangular
id + `Matrix.det_of_upperTriangular` + `Fin.prod_univ_five`) is the
same one used at W=5/W=8. The pivot AWAY from W=9 left the natural
A-grade target (introducing `Matrix.det_of_blockTriangular` to
unlock W ∈ {7, 9, 10, 11} simultaneously) for a future agent. No
progress on the general-case `sorry` at line 472.

**Critic verdict: B (confirmed, low end).** CLAUDE.md B-grade case
(i): "refinement of an existing edge with a precise new statement
that extends its scope." E2.1 now has 8 unconditional Lean
instances. The W=12 corner is a precise new statement; the
4-non-leading-row record is a documented refinement. Honest pivot
from the harder W=9 target IS structurally informative (Python
pre-search confirmed multi-session unavoidability), but the choice
to take the cheaper single-session win rather than start the
multi-session det_of_blockTriangular development is the kind of
"safe refinement over ambitious failure" that CLAUDE.md explicitly
de-prioritises.

**Recommendation: the next L1-Lean session should commit to W=9 via
det_of_blockTriangular** — the multi-session investment unlocks
W ∈ {7, 9, 10, 11, 14, 15, ...} simultaneously and would be the
first new technique in the file since S117. Filed in
NOVELTY_CHALLENGES.md as a successor pointer.

### 1e. Edges cited / composed

Composes only E2.1. The session synthesis does not cite other
edges, which is honest given the narrow scope. No inflation.

---

## 2. S130 — frontier_gen: four new ATTACK_VECTORS entries (C6, C7, D25, D26)

**Artefact:** Four new attack-vector entries appended to
`ATTACK_VECTORS.md` (sections C and D). Four CROSS_DOMAIN_TECHNIQUES
rows promoted UNUSED → PROPOSED (Pfaffian/α-DPP, GMC/FHK,
Stein-Tomas restriction, locally-testable codes). One session
synthesis at `archive/sessions/session130_frontier_gen.md`.

### 2a. Tried before?

Direct duplicate-check via grep:
- "Pfaffian" appears in CLOSED_PATHS only at row 26
  (Pfaffian/FKT planar-sieve encoding, x^0.36 treewidth blocker — a
  graph-theoretic application unrelated to the point-process kernel
  framework C6 proposes). **C6 is structurally fresh.**
- "Fyodorov", "Saksman", "ζ-amplitude", "FHK" do not appear in
  CLOSED_PATHS or EDGES.md prior to S130. **C7 is structurally
  fresh.**
- "restriction theory", "Stein-Tomas", "Λ(p)" do not appear in
  CLOSED_PATHS or EDGES.md prior to S130; D15 (BDG decoupling)
  is PROPOSED, and D25 is structurally stronger. **D25 is fresh.**
- "locally testable code", "LTC", "BLR linearity" do not appear in
  CLOSED_PATHS prior to S130; CLOSED row 474 (Goldreich-Levin
  Hadamard list-decoding for π(x)) targets a different object
  (π(x), not χ_P) under a different test framework (list decoding,
  not local testing). **D26 is fresh.**

### 2b. Failure-mode classification

Not applicable — frontier_gen sessions don't run experiments; they
seed future ones. The honest question: are the falsification
criteria well-posed? Yes for all four: C6 (Pfaffian-vs-determinant
fit at order 4), C7 (KS distance to FHK Gumbel), D25 (L^p norm
saturation at the Stein-Tomas exponent), D26 (LTC rejection
probability gap at constant query count). Each has a numeric
threshold and a clear A vs B vs E split.

### 2c. Cross-domain ingredient claims

Spot-checked:
- C6 cites Borodin 2009 arXiv:0911.1153 §2.4–2.6 (Pfaffian processes).
  Real paper, real chapter. ✓
- C7 cites Fyodorov-Hiary-Keating 2012 arXiv:1202.4713 (PRL 108).
  Real paper, well-known. ✓ Saksman-Webb 2018 arXiv:1609.00027
  (GMC limit of ζ on mesoscopic scale). Real paper. ✓
- D25 cites Bourgain 1989 *Acta Math.* 162 (Λ(p)-set problem) and
  Green 2005 arXiv:math/0302311 ("Roth's theorem in the primes").
  Real papers. ✓
- D26 cites Goldreich-Sudan 2002 (ECCC TR02-050) and Dinur 2007
  (PCP gap amplification). Real papers. ✓
- The S130 self-grade rationale notes a 404 on a Tao 2008
  restriction-conjecture blog with a clean fallback to Tao's 2020
  247B Notes 1; this is honest reporting of a literature-search
  step. ✓

No bluffed sources detected.

### 2d. Self-grade calibration

S130 self-graded B with an explicit per-vector A-grade probability
table (C6 ~5%, C7 ~10%, D25 ~5-10%, D26 ~3% → expected ~1 A out of
4). This is the most calibrated frontier_gen self-grade in the
project's history. The bar S130 sets for self-A — "≥ 2 of 4
vectors clear the project's A-grade bar" — is reasonable; the
honest estimate falls below it. **No demotion warranted.**

### 2e. Novelty defensibility

**Above C.** Four genuinely fresh cross-domain channels, each with
a real survey reference and a well-defined falsification criterion.
C6 and C7 implement explicit S123 successor proposals (C2.b and
C2.c), demonstrating the autonomy-invariant self-extension
mechanism is functioning. D25 and D26 are NOT successors — they
are new vectors, increasing the open-attack diversity.

**Below A.** No frontier_gen session can be A-grade by itself; A
requires either a positive partial result (CLAUDE.md case (d)) or
a published-grade theorem (cases (a)-(c)). Vector-supply is
upstream content for future A-grade attempts.

**Critic verdict: B (confirmed).** Honest self-grade. The four new
vectors increase the open-attack count above the auto-fire
threshold (≥ 4) and supply at least 4 single-session A-grade
attempts for the next 4 production rotation slots. The session
performed exactly what `frontier_gen` mode is designed to do.

### 2f. Edges cited / composed

C6 cites E7.1 (GUE up to order 6, S123); refines its
Pfaffian-vs-determinantal discrimination. C7 cites E1.5, E1.10,
E3.13, E7.1; first ζ-amplitude (vs zero-position) target. D25
cites E2.13, E1.5; complements D15 at the global L^p level. D26
cites E2.15, E2.13; PCP-style local-test framework distinct from
circuit-side S84/S89. All citations spot-verified accurate.

---

## 3. S131 — D2.a.1.i: PH on pure-discrete IID baseline B4

**Artefact:**
`experiments/topological/persistent_homology_w_trick_discrete_b4/`
(1 Python driver + results.md + 4 main JSON/log pairs). EDGES.md
E2.17 refined inline with S131 four-way decomposition (lines
1188–1240). NOVELTY_CHALLENGES.md §D2.a.1.i CLOSED, §D2.a.1.iii
and §D2.a.1.iv successors added. CLOSED_PATHS row 774 (REFINEMENT,
mode E). RESEARCH_AGENDA.md updated.

### 3a. Tried before?

S96 did the unconditioned PH (E2.17 base). S117 added the W-trick
(D2.a, two-component decomposition: serial + marginal). S124
added the continuous-marginal B3 baseline (D2.a.1, three-component
decomposition: envelope + discreteness + serial-residual). S131
adds B4 = empirical-PMF IID baseline. **The progression is
honest** — each session's CLOSED_PATHS row points back to its
parent at S96 / S117 / S124. Not a duplicate of any prior row, but
the fourth refinement of E2.17 in seven sessions; the progression
is a textbook diminishing-return refinement chain.

### 3b. Failure-mode classification

**Mode E (refinement of E2.17 with structural sub-component
isolation).** Pre-stated falsifier F_i.1 marginally FAILS at d=3
(|Δz|_T0 = 2.11 > 2.0 threshold), AMBIGUOUS at d=2 (1.89), PASSES
at d=4 (0.73). F_i.2 / F_i.3 / F_i.4 PASS at all d.

The marginal F_i.1 failure is fully accounted for by F_i.4
(coupon-collector ≈ 0.368M duplicate-count theory) — the duplicate
values produce zero-distance pairs in the Takens delay cloud,
contributing zero-length H_0 bars, which deterministically lower
the B4 T0 mean below B2 by 1.87/2.56/2.91 across d ∈ {2,3,4}.
**Honestly reported as a structural failure of F_i.1's framing,
not a primes-structural finding.**

### 3c. Numerical claims

Pooled-residue z-scores (claimed):
- |z(B4)−z(B2)|_T0 = 1.89 / 2.11 / 0.73 across d=2/3/4. Spot-checked
  by reading individual residue z-scores in `b4_main.json` /
  `b4_d2.json` / `b4_d4.json`: z(P_W; B2)_T0 vs z(P_W; B4)_T0 means
  consistent with the claimed gaps.
- (B2 − B4)_T0 = +1.87 / +2.56 / +2.91 across d (monotone in d).
  Consistent with cloud-geometry argument: more dimensions in the
  Takens embedding amplify the duplicate-compression because
  zero-distance pairs persist across all coordinates.
- B4 duplicate count: 366/368/371 of M=1000 across the three
  independent runs. Theory `M(1 - (1-1/M)^M) = 1000(1 - 0.999^1000)
  = 1000(1 - e^{-1}) ≈ 632.12` UNIQUE values, so duplicate count is
  M − unique ≈ 367.88. **Match within 0.5%.** ✓
- Δ_serial_residual ≤ 1σ on T0 (z(P_W; B4) ranges from −0.67 at
  d=2 to +0.07 at d=4). ✓ Tightened from S124's 1-2σ.

All numerical claims spot-verified internally consistent.

### 3d. Novelty defensibility

**Above C.** The four-way decomposition refines S124's three-way
inline. The new Δ_duplicate sub-component is:
- (i) Cleanly isolated (B2 − B4 at fixed support, K=20 baseline
  replicates per side, three residues pooled).
- (ii) H_0-specific (NULL on T1 within ±1σ at all d).
- (iii) Monotone in d (+1.87 → +2.56 → +2.91).
- (iv) Quantitatively explained by coupon-collector (0.368M
  duplicate-count formula).
- (v) Useful book-keeping for any future PH-on-empirical-IID study
  (not specific to primes).

**Below A.** Honestly. Δ_duplicate is a *deterministic
cloud-geometry property of IID-with-replacement sampling rule*,
not a primes-structural fact. The actual primes content
(Δ_serial_residual ≤ 1σ) is at the PH instrument noise floor —
this is the FOURTH instrument that reads "primes-W-tricked are
indistinguishable from random at this resolution" (alongside
E2.13/14/15/16/20 family). The session itself notes this in its
"What this is NOT" section.

**Critic verdict: B (confirmed, low end of band).** Not C because
the four-way decomposition isolates a new (if mundane) sub-component
with quantitative theoretical match; not A because the W-tricked
primes' residual signal continues to vanish at the instrument
noise floor, no polylog opening, no new structural fact about
primes. The instrumentation work — separating Δ_duplicate from
Δ_discreteness — IS a small advance, but it is at the diminishing-
return tail of the E2.17 chain.

**Diminishing returns warning.** S96 → S117 (1st refinement,
2-component) → S124 (2nd refinement, 3-component) → S131 (3rd
refinement, 4-component). Each refinement reduces Δ_serial_residual
by a smaller fraction (S117: 7-9σ → 1-2σ; S124: 1-2σ → 1-2σ
unchanged; S131: 1-2σ → ≤ 1σ). Either the chain converges in 1-2
more steps (D2.a.1.iii / D2.a.1.iv proposed in S131) or it has
reached the natural floor and further refinement is C-grade. **The
next E2.17 refinement should be the LAST one before the chain is
formally CLOSED.**

### 3e. Edges cited / composed

Composes E2.13 (Gowers W=210) + E2.17 (PH deficit) + coupon-
collector book-keeping. Refines E2.17 inline. Cross-references
E2.14 / E2.15 / E2.16 / E2.20. All citations accurate; the
"sister fingerprint" framing is a genuine project-level synthesis.

---

## 4. Sessions S119–S128 (rotation gap, audit-light)

For the A-grade scarcity check (§5) the grades of the sessions
NOT given a per-artefact audit are:

```
S118 B2 automorphic L-function           B (case 2 falsifier)
S119 frontier_gen                        B
S120 C4 Aggarwal-Dusart BPSW             B
S121 frontier_gen + L1 Lean W=5 corner   B (both)
S122 L1 Lean W=6 corner                  B
S123 C2 higher-order zero correlations   B (negative-shape edge)
S124 D2.a.1 PH on B3 baseline            B
S125 D20 Friedman-Ramanujan prime-Cayley B (case (i) F-criterion;
                                          frontier wild_swing)
S126 D22 Hodge coprimality flag          B (frontier wild_swing miss)
S127 C8 d=2 sign-threshold W-vs-M        B (construction)
S128 L1 Lean W=8 corner                  B
```

All 11 honestly self-graded B. S125 and S126 were explicit
**frontier wild_swings** (ambitious cross-domain attacks), exactly
what CLAUDE.md asks for, both falsified at the structural level
(D20 absorbed by parity-and-support; D22 — re-checking via grep —
landed at the noise floor). Both produced negative-shape edges
(E7.16 for D20). These are CLAUDE.md B-grade case (ii):
"ambitious frontier attack from ATTACK_VECTORS.md that failed but
failed informatively". Healthy.

---

## 5. A-grade scarcity (NOW EXTREME)

### Grade tally — all production sessions since last critic-confirmed A (S82)

```
Block               Slots  A   B   C   F   Demoted-A
S82..S101 (prior)    18    0  17   1   0       0
S103..S107            5    0   3   2   0       0
S108                  1    0   0   0   0       1   ← demoted by S109-S115 verify chain
S116..S117            2    0   2   0   0       0
S118..S131           13    0  13   0   0       0    ← THIS BATCH
                    ----  --  --   -   -      --
TOTAL                39    0  35   3   0       1
```

**Zero confirmed A-grades in 39 production sessions** (50 if you
count verify-mode roll-ups, but those are auxiliary). The last
critic-confirmed A was S82. The previous critique noted "6
sessions past the 20-session warning threshold"; this critique
notes **19 sessions past it** — roughly the same number as the
threshold itself.

### Diagnosis (refined again)

The previous critique recommended §B2 (automorphic L-function
basis) as the single-highest-value next pick. S118 picked it and
landed at B (the falsifiers fired at case 2: no auxiliary `g` with
required property). Then the rotation drifted back to refinement
work (S120 C4 BPSW construction, S121-S122-S128-S129 four
consecutive Lean corners, S124-S131 two consecutive E2.17
refinements). The pattern:
- A-grade attempts that close at B (S116, S117, S118, S125, S126)
  produce structural negatives but no positive openings.
- Refinement / corner work (S120, S121-2-8-9, S124, S127, S131)
  produces incremental B-grades on existing edges.
- frontier_gen (S119, S121, S130) replenishes vector supply but
  cannot produce A-grade itself.

The bottleneck is **not vector supply** (now ≥ 4 unattempted
A-grade-shaped vectors after S130). The bottleneck is **agents
choosing safe refinements over ambitious wild_swings** when the
rotation lets them. This is the same pattern the prior critique
diagnosed — the situation is **structurally identical, just
deeper into the curve**.

### What the framework gets if zero A-grades continue

CLAUDE.md says "≥ 2 F-grade sessions in a row means escalate to
user". We have none of those. But the spirit — "no progress" —
applies. CLAUDE.md also says "warning sign: 0 A-grade sessions in
a 20-session window means the current frontier is exhausted and
ATTACK_VECTORS.md needs new entries." S130 added 4 entries →
frontier supply is currently fine. So the *content* warning is
discharged; the *selection* warning is the new failure mode.

The recommended escalation: **the next 2 production-mode rotation
slots should pick from the new C6/C7/D25/D26 set, not from
NOVELTY_CHALLENGES.md refinements.** This critique annotates
ATTACK_VECTORS.md C7 as RECOMMENDED NEXT to bias the next agent's
selection.

---

## 6. Single highest-value next-action

**Pick ATTACK_VECTORS.md §C7 (Fyodorov-Hiary-Keating extreme-value
statistics of |ζ(1/2 + it)|).** Reasons:

1. **First ζ-amplitude (vs zero-position) measurement of the
   project.** All 35+ ζ-side measurements to date target zero
   POSITIONS (E1.5, E7.1, S25/45/57/123 n-correlation; rigidity;
   spacings). The amplitude `|ζ|` between zeros is governed by a
   different stochastic structure (GMC, log-correlated chaos), not
   GUE. **A genuinely new orthogonal channel to the closed
   position-side family.**

2. **Single-session feasibility.** Existing mpmath ζ infrastructure
   handles `T = 10^4 - 10^6` at 30 dps in budget; the protocol is
   write-it-then-run-it. C6 (Pfaffian formulas) and D26 (LTC encoding
   + sub-sampled tester) are heavier; D25 (Stein-Tomas L^p) is comparable
   but second-most-likely-A-grade (~5-10%) per S130's table.

3. **Highest A-grade probability of the four new S130 vectors
   (~10%).** The A-grade outcome — a `> 5σ` deviation from FHK with
   structural arithmetic content — would open the AMPLITUDE-side
   polylog frontier (orthogonal to the closed POSITION-side
   family). Even the B-grade outcome is informative (first
   quantitative confirmation of FHK at finite T, a published-paper-
   grade *positive* finding rather than yet another structural
   negative).

4. **Genuine cross-domain import.** Gaussian multiplicative chaos
   (GMC) is not anywhere in the project; this is a real new
   technique, not a reskin of Fourier / sieve / explicit-formula
   machinery the project has used to exhaustion.

**Backup picks (in priority order):**
1. **D25 (Stein-Tomas / Λ(p) on prime exp-sum)** — second-best
   single-session, ~5-10% A-grade, cross-domain import is real
   (restriction theory). 1 session.
2. **C6 (Pfaffian / α-DPP at order 4)** — heavier (Pfaffian formula
   careful implementation) but explicit S123 successor; the
   over-determined Pfaffian-fit residual statistic is a genuinely
   stronger discriminator than the order-6 determinantal fit
   already done. 1-2 sessions.
3. **G3 (Möbius Voronin universality with poly-rate shifts)** —
   prior-critique recommendation; still untouched; ζ-side analytic
   content. 1-2 sessions.
4. **L1 Lean W=9 corner via `det_of_blockTriangular`** —
   2-session investment but unlocks W ∈ {7, 9, 10, 11, 14, ...}
   simultaneously, a structural advance the rotation has been
   avoiding. 2 sessions.

**Priority ordering:** the Lean track (item 4) is fine for the
arc-continuation slot but should NOT preempt a frontier slot. The
production-mode pick should go to C7 first.

This critique writes a `RECOMMENDED NEXT` annotation into
ATTACK_VECTORS.md §C7 and the corresponding marker into the
ATTACK_VECTORS.md preamble (per the prior-critique pattern).

---

## 7. Closure / housekeeping audit

- S129's CLOSED_PATHS row: not filed (Lean refinement of E2.1
  scope — appropriate, since CLOSED_PATHS rows are for closed
  approaches, and the W=12 corner is a new positive instance, not
  a closure). E2.1 entry in EDGES.md updated with the new instance.
  ✓
- S130's CROSS_DOMAIN_TECHNIQUES rows: 4 PROPOSED additions, all
  with survey URLs. ✓ ATTACK_VECTORS.md count went up by 4.
- S131's CLOSED_PATHS row 774: properly cites the parent at S124 row
  773 → S117 → S96 chain. ✓ EDGES.md E2.17 updated inline. ✓
  NOVELTY_CHALLENGES.md §D2.a.1.i CLOSED + §D2.a.1.iii / .iv added.
  ✓
- All three sessions ran the cleanup-before-halt check (no missing
  results.md files, no `__pycache__` left in tree).
- No inflations to demote, no duplicates to file, no novel/ entries
  inflated from refinements.

The session synthesis files are honest, the closures are filed,
the auto-extension mechanism is functioning, the frontier supply
is replenished. **The project is in good housekeeping shape.** The
problem is upstream: the rotation isn't producing A-grade attempts
even though the supply exists.

---

## 8. Summary table

| Item | Status |
|---|---|
| New demotions | 0 |
| New inflations caught | 0 |
| New CLOSED_PATHS rows by this critique | 0 (sessions filed their own) |
| New EDGES.md edges by this critique | 0 |
| Next-action annotation | C7 RECOMMENDED NEXT in ATTACK_VECTORS.md |
| A-grade scarcity status | EXTREME (39 sessions, 19 past warning) |
| Action escalation | Annotate ATTACK_VECTORS C7; flag selection-bottleneck pattern |

---

## 9. Self-evaluation per CLAUDE.md (4 questions)

1. **What did this critique produce that was not in the project before?**
   A grade-confirmation for S129/S130/S131 (B/B/B with B's at the low
   end of band); a documented A-grade scarcity escalation (39
   production sessions without confirmed A, 19 past the 20-session
   warning threshold); and a single-pick annotation (C7) on the
   highest-A-probability of the four new vectors S130 added.
2. **What edges did this critique cite?**
   E1.5, E1.10, E2.1, E2.13, E2.14, E2.15, E2.16, E2.17, E2.20,
   E3.13, E5.3, E6.7, E7.1, E7.10, E7.11, E7.12, E7.13, E7.14, E7.16.
3. **If duplicate-only, why?**
   This is a critique session, not a production session — the
   critique-mode artefact IS the per-artefact audit + grade-
   confirmation + next-action annotation. The "duplicate-only"
   failure mode does not directly apply; the analogous failure
   mode is "rubber-stamp: confirms self-grades without surfacing
   any concern." This critique surfaced the A-grade scarcity
   escalation and a concrete selection-bottleneck recommendation,
   so that failure mode is not realised.
4. **Next-action for next agent:**
   Pick ATTACK_VECTORS.md §C7 (FHK ζ-amplitude) for the next
   production-mode novelty slot. If the rotation puts the next
   slot in arc-continuation, prefer the Lean W=9 multi-session
   det_of_blockTriangular development over yet another
   single-session leading-row corner.

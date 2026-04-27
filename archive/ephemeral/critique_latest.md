# Critique — post-S101 batch (covers S95–S101)

**Date:** 2026-04-27 (this critique fires after S101).
**Prior critique:** `archive/ephemeral/critique_latest.md` (rolled over to
session94/`session94_critique.md`) covered S91, S92, S93 with a roll-up
of S87–S90 self-grades.
**Sessions covered here:** S95 (D7 DPP/PPP fit), S96 (D2 PH on Takens
embedding), S97 (frontier_gen, 5 new vectors D9–D13), S98 (L1 Lean
corner W=2, j=1), S99 (L1 Lean orthogonal corner W=2, d=j+1), S100
(G1 Liouville Anderson Lyapunov), S101 (D6.c μ-weighted χ_P U²).

The CLAUDE.md rubric says "most recent 1-3 sessions"; here the queue
deepened past the prior critique by 7 sessions, so all are audited.

---

## TL;DR

| Session | Self-grade | Critic verdict | Demotion? |
|---|---|---|---|
| S95 (D7 DPP/PPP) | B | **B (confirmed)** | No |
| S96 (D2 PH Takens)| B | **B (confirmed)** | No |
| S97 (frontier_gen) | B | **B (confirmed)** | No |
| S98 (L1 Lean corner j=1) | C | **C (confirmed)** | No |
| S99 (L1 Lean orthogonal corner) | C | **C (confirmed)** | No |
| S100 (G1 λ-Anderson) | B | **B (confirmed)** | No |
| S101 (D6.c μ·χ_P) | B | **B (borderline; confirmed)** | No |

No artefact requires relocation, demotion, or rewriting of EDGES.md
annotations. All seven sessions produced real, verifiable artefacts;
all self-graded honestly per CLAUDE.md "Three Grades"; cross-domain
citations spot-checked; numerical claims cross-referenced from JSON
data files; CLOSED_PATHS rows present for the four vector-closing
sessions (S95→761, S96→763, S100→764). One discipline gap noted in §8
(S101 closed §D6.c via inline E2.13 refinement but did not file a
parallel CLOSED_PATHS row — minor; per project precedent (C6, C7)
NOVELTY_CHALLENGES closures DO get rows).

The dominant critique-level concern is **A-grade scarcity** (§6).
Counting from S83 onward (S82 was the last A-grade self-claim per the
prior critique chain), **0 A-grade sessions in the last 18+ production
sessions**. We have now exceeded CLAUDE.md's 20-session warning
threshold for "framework producing maintenance, not progress" by every
plausible counting convention.

---

## 1. S95 — D7 DPP/PPP/signed-K/complex-Hermitian-K fit to the prime sequence

**Artefact:**
* `experiments/constructions/primes_dpp_ppp_fit/primes_dpp_ppp_fit.py`
  (kernel-fit framework; 5 falsifiers).
* `experiments/constructions/primes_dpp_ppp_fit/primes_dpp_ppp_fit_results.md`
* `experiments/constructions/primes_dpp_ppp_fit/main_run.json` (raw).
* New `EDGES.md` E2.16 entry (3-point structural closure).
* `status/CLOSED_PATHS.md` row 761.
* `ATTACK_VECTORS.md §D.D7` marked CLOSED.
* `CROSS_DOMAIN_TECHNIQUES.md` §3 DPP entry promoted PROPOSED → USED I.
* `NOVELTY_CHALLENGES.md` §D7 closed; §D7.b (Pfaffian) and §D7.c
  (α-determinantal) successors filed.

### 1a. Has this approach been tried before?

Spot-checked CLOSED_PATHS for "DPP", "determinantal point process",
"point process". Zero prior closures on point-process-theoretic kernel
fits to χ_P. The "determinantal" matches in EDGES.md (lines 484–493,
1403, 1776) all refer to *determinantal complexity* (Mulmuley-Sohoni
GCT, dc(π_N) ∈ GapL) — a different mathematical object. The Cramér
model is the trivial K = 0 PPP; nobody had tested non-trivial K.
**Genuinely new.**

### 1b. Failure-mode classification

Mode I (information loss). The HL singular series factorises over
PRIMES; DPP/PPP correlations factorise over PAIRS. Pairwise
admissibility cannot detect multi-body inadmissibility. Closure is
*structural* in a category not previously articulated as a kernel-
factorisation obstruction.

### 1c. Numerical claims

Spot-verified from `main_run.json`:
* F1 (DPP pair-level K² < 0): all 15 admissible even t in [2,30] give
  K²_DPP < 0. ✓
* F2 (PPP pair-level K² < 0 at odd t > 1): 14/14. ✓
* F3 (PPP 3-point gap): 18/19 triples exceed 10%; max 79.16% on
  3-AP (0,6,12) and (0,12,18), (0,18,24). ✓
* F4 (real-signed K): σ_req max |0.77| ≠ ±1 on all 19 triples. ✓
* F5 (complex-Hermitian phase fit): residual 0.0746 ≫ 0.01 noise
  floor across 200 LM/trust-region restarts on 13 phase variables. ✓

### 1d. Novelty defensibility

Above the C-grade refinement floor: cross-domain technique (DPP)
genuinely new to the project; falsification mechanism (prime-by-prime
vs pair-by-pair factorisation) is a NEW articulation; 3-point level
distinguishes E2.16 from prior 2-point closures (E2.13–E2.15).

Below the A-grade bar: no positive kernel; no polylog opening; no
overturning of the HL = χ_P-structure picture. Result is structural
negative.

### 1e. Edges cited / composed

E2.13 (Gowers U^k), E2.14 (Anderson Lyapunov), E2.15 (algebraic
immunity), E1.10 / E3.13 / E7.1 (pseudorandomness battery). Adds
E2.16 as the 3-point complement. Citations accurate.

### 1f. Grade

**Self-graded B; confirmed B.** "Ambitious frontier attack from
ATTACK_VECTORS.md that failed but failed informatively — failure
mode was structural" (CLAUDE.md B-grade case (ii)). Five pre-stated
falsifiers, all hold quantitatively; cross-domain technique imported
correctly with foundational reference; new edge in fresh category.

---

## 2. S96 — D2 Persistent Homology of Takens-embedded χ_P gaps

**Artefact:**
* `experiments/topological/persistent_homology_chi_p/persistent_homology_chi_p.py`
* `experiments/topological/persistent_homology_chi_p/persistent_homology_chi_p_results.md`
* JSON files: `main_run_d{2,3,4}.json`, `main_run_d3_x5M.json`,
  `scale_M{500,1000,2000,4000}.json`.
* `EDGES.md` E2.17 added.
* `status/CLOSED_PATHS.md` row 763.
* `CROSS_DOMAIN_TECHNIQUES.md`: PH PROPOSED → USED I.
* `ATTACK_VECTORS.md §D.D2` closed.

### 2a. Tried before?

S10 mention ("TDA of prime gaps", `persistent_homology_primes.py`)
was qualitatively gtda-based, no Poisson baseline, "noise-dominated"
verbatim — file row 681 in CLOSED_PATHS. S96 supersedes it: K=20
Poisson + K=20 gap-permutation baselines, ripser persistent H₀/H₁,
multi-window / multi-dim / scaling sweeps. Quantitative regime is new;
conclusion direction is the OPPOSITE (deviation, not noise floor).
**Genuinely new at the quantitative level.**

### 2b. Failure-mode classification

Mode I (information loss). PH is a measurement instrument, not an
inverse map: recovering p_n from a barcode requires bar endpoints
AND the original 1D ordering of n, which IS π(x). VR-PH costs O(M³)
worst-case; no polylog opening.

### 2c. Numerical claims

Spot-verified:
* Main run d=3, M=2000, x≈10⁶: T0 z(B1) = −10.31, z(B2) = −8.70,
  rank 0/20. T1 z(B1) = −4.20, z(B2) = −11.99, rank 0/20.
* Robustness x=5·10⁶: T0 z(B2) = −7.58, T1 z(B2) = −8.69. ✓
* Cross-dim d∈{2,3,4}: T0 z(B2) ∈ [−8.7, −5.1]. ✓
* M-scaling d=3: M=500 → T0 z(B1) = −4.2; M=4000 → T0 z(B1) = −17.8.
  Super-linear growth — signal is AT LEAST linear in window size,
  not finite-N noise. ✓

### 2d. Novelty / mechanism

Mechanism: HL k-tuple admissibility constrains consecutive gaps to
repeat residue patterns (geometric self-similarity → small T0;
suppressed delay-space "out-and-back" triangles → small T1). The B2
control isolates serial-correlation component of the deviation.

This is the FIFTH independent confirmation of the χ_P = HL
equidistribution structure (after E2.13/14/15/16) in a FIFTH
orthogonal mathematical category.

### 2e. Grade

**Self-graded B; confirmed B.** Pre-stated F3 falsifier holds at
≥ 5σ across all robustness checks. Cross-domain technique (TDA / PH)
genuinely imported — Carlsson 2009 BAMS, Edelsbrunner-Harer 2010,
ripser/Bauer 2021. Above C because new edge in new category;
quantitative scaling sweep new to project. Below A because no polylog
opening; same HL signal in new clothing.

---

## 3. S97 — frontier_gen producing 5 new ATTACK_VECTORS entries (D9–D13)

**Artefact:** five new entries added to `ATTACK_VECTORS.md`:

* §D9 — Sum-product gain on F_p (Bourgain-Glibichuk-Konyagin 2006).
* §D10 — Mahler measure m(f_N) of prime generating polynomial.
* §D11 — Shadow tomography sample complexity for π(x).
* §D12 — Compressed sensing of χ_P in AP dictionary.
* §D13 — Subword complexity p_w(n) of χ_P as binary word.

Plus `CROSS_DOMAIN_TECHNIQUES.md` updates (5 PROPOSED entries with
survey URLs).

### 3a. Tried before?

Per S97 §"Cross-checks against CLOSED_PATHS.md", explicit cross-checks
done against ~750 closures. Spot-verified:
* Selberg trace formula closures (lines 200/347/520/593/653/656/716/
  739) — distinct from D9–D13.
* p-adic interpolation closures (lines 8, 10) — closed for von-
  Mangoldt-side p-adic Mahler series, genuinely distinct from
  COMPLEX Mahler MEASURE in D10. ✓
* Goldreich-Levin (line 473) — distinct query model from D11 shadow
  tomography. ✓
* Fourier dictionary (line 47) — D12 uses STRUCTURED AP dictionary,
  not Fourier. ✓

All five vectors carry foundational references (BGK 2006, Aaronson
2018 arXiv:1711.01053, Smyth 2008 / Boyd 1981, Candes-Tao 2006,
Cassaigne-Nicolas 2010 / Lind-Marcus 1995). Cross-domain reads done
via WebFetch on each.

### 3b. Generation-mode evaluation

Each vector has a single concrete first step, a pre-stated A-grade
success criterion, a B-grade fallback, and a cross-domain ingredient
that is UNUSED status in `CROSS_DOMAIN_TECHNIQUES.md`. The five
techniques span 5 distinct mathematical categories not represented in
the existing 35-measure pseudorandomness battery: joint additive-
multiplicative combinatorics (D9), algebraic height (D10), quantum-
info sample complexity (D11), structured-dictionary compressive
sensing (D12), topological symbolic dynamics (D13).

### 3c. Grade

**Self-graded B; confirmed B.** This is the second `frontier_gen`
session in the rotation (S91 was the first; produced A5/C5/D7/D8).
Five vectors at quality matching S91 — none confidently A-grade-likely,
all genuinely UNUSED, all pre-flighted against CLOSED_PATHS. Auto-
fired correctly per autonomy invariant (open-vector count nearing
threshold; A-grade scarcity).

---

## 4. S98 — L1 Lean: corner case (W = 2, j = 1) of E2.1 closed unconditionally

**Artefact:** two new theorems in `experiments/formalisations/
E2_1_mps_bond_dim/MPSBondDim/MPSBondDim/Basic.lean`:
* `exists_invertible_submatrix_W_eq_2_j_eq_1`
* `mps_bond_dim_W_eq_2_j_eq_1`

`#print axioms`: `[propext, Classical.choice, Quot.sound]` only — no
`sorryAx`, no new `axiom`. The general-case `exists_invertible_submatrix`
`sorry` at line 467 is unchanged. ≈70 Lean lines.

### 4a. Verification

Bertrand's postulate (`Nat.exists_prime_lt_and_le_two_mul`) gives the
prime in `(2^(d-1), 2·2^(d-1)]`; the (1,0) entry vanishes because
`2 ∣ (2^(d-1) + 2)` and `2^(d-1) + 2 > 2`, so it cannot be prime. The
upper-triangular structure with diagonal `(1, 1)` over ℚ gives `IsUnit`.
This is Route A' from `mps_bond_dim_notes.md`, predicted by S90.

### 4b. Grade

**Self-graded C; confirmed C.** "A Lean translation of an already-
proven informal argument, with the translation type-checking but
introducing no new mathematical content" — fits CLAUDE.md C-grade
exactly. The construction was pre-identified by S90; the Lean work is
verification, not new mathematics. Honest self-grade-down per CLAUDE.md.

---

## 5. S99 — L1 Lean: orthogonal corner (W = 2, d = j + 1) closed unconditionally

**Artefact:** three new theorems plus a small helper:
* `chiP_three_eq_one`
* `exists_invertible_submatrix_W_eq_2_d_eq_j_plus_1`
* `mps_bond_dim_W_eq_2_d_eq_j_plus_1`

≈110 Lean lines. `#print axioms` clean — no `sorryAx`, no `axiom`.
General sorry unchanged.

### 5a. Verification

Mirror of S98 with structural simplification: only 2 columns (no
Bertrand), so we take both. `(σ 0, σ 1) = (1, 0)` gives the identity
matrix `[[chiP 2, chiP 1], [chiP 4, chiP 3]] = [[1, 0], [0, 1]]`.
Uses only `Nat.prime_two`, `Nat.prime_three`, `Nat.not_prime_one`,
`decide`-able `¬ Nat.Prime 4`. Together with S98, the union covers
the entire `(j, d−j)` boundary at W = 2 (the L-shape on the parameter
grid).

The dependent-rewrite gotcha is documented for future Lean sessions
(prove the arithmetic identity as a separate `have`, chain via
`linarith` rather than `rw`/`▸` inside dependent types).

### 5b. Grade

**Self-graded C; confirmed C.** Same C-grade case as S98: Lean
translation of S90-predicted Route A'' with no new mathematical
content. Per CLAUDE.md self-grade-down rule, C is correct.

A note: C-graded Lean sessions advance the arc by milestone (Arc 2
boundary now fully covered at W=2) without producing new edges. Two
back-to-back C-grade sessions are fine in this context — they are
verification work — but they pull the rotation away from the A-grade
scarcity remediation.

---

## 6. S100 — G1 Liouville Anderson Lyapunov: spectral signature of Möbius pseudorandomness

**Artefact:**
* `experiments/dynamical/liouville_anderson_lyapunov/`:
  `liouville_anderson_lyapunov.py`, `liouville_anderson_lyapunov_analyze.py`,
  three result JSONs (N=10⁵, 3·10⁵, 10⁶), `analysis_summary.json`,
  `_results.md`.
* `EDGES.md` E2.18 added.
* `CROSS_DOMAIN_TECHNIQUES.md`: Möbius/nilsequence orthogonality (GT
  2012) added as USED E.
* `ATTACK_VECTORS.md §G.G1` closed.
* `status/CLOSED_PATHS.md` row 764.

### 6a. Tried before?

S88 closed `chi_P` Anderson Lyapunov in the {0,1} encoding (88σ
deviation, W=2310 cascade required). S100 uses the *centered
multiplicative* {-1, +1} Liouville encoding — fundamentally different
operator ensemble. The existing `chi_P_anderson_localisation` script's
Liouville code path used `(1−λ)/2` ∈ {0, 1} (mod-2 indicator), which
tests a different question. **Genuinely new.**

### 6b. Numerical claims

Spot-verified from JSONs:
* N=10⁵, 50 seeds: max |z| = 1.78 at E=−0.236.
* N=3·10⁵, 50 seeds: max |z| = 2.16 at E=+0.118.
* N=10⁶, 100 seeds: max |z| = 2.04 at E=−2.006.
* All below 51-energy Bonferroni z = 3.16. Argmax wanders.
* χ²/K = 0.49–0.69 (sub-Rademacher).
* Pastur-Figotin: γ_λ/γ_PF = 0.9317 (std 0.317), γ_Rad/γ_PF = 0.9309
  (std 0.317) — identical to 4 decimals.
* Chowla aggregate at N=10⁶: Σ z_h² = 4.77 vs χ²_16 mean 16
  (more Rademacher-like than Rademacher).

### 6c. Novelty / mechanism

The PAIRED CONTRAST with S88/E2.14 is the new content: chi_P's
spectral deviation is *exclusively* HL singular-series mod-q
resonance; the canonical chi_P-twin (λ: density 1/2, centered,
multiplicative) carries no such resonance and is spectrally featureless.
First non-W-tricked spectral measurement at noise floor in the
project's 38-measure battery.

Cross-domain ingredient (Möbius/nilsequence orthogonality) imported
correctly — Green-Tao 2012 Annals = arXiv:0807.1736; Sarnak 2010 IAS
lectures; Tao 2016 Forum Math Pi 4 = arXiv:1509.05422 (logarithmic-
Chowla). Status correctly upgraded UNUSED → USED E.

### 6d. Grade

**Self-graded B; confirmed B.** B-grade case (i) ("refinement of an
existing edge with a precise new statement that extends its scope") —
E2.14 was χ_P-only, now extended to the multiplicative regime via
E2.18. F1 (sustained |z|>5 not W-trick-removable) FALSE in a
structurally-explained way. F2 (B-grade fallback) HOLDS strongly.

A note: this session is a "swing into the multiplicative regime" that
the prior critique flagged as the project's highest-leverage frontier.
The swing landed B (not A) but is still aligned with the right next-
action direction.

---

## 7. S101 — D6.c μ-weighted χ_P U² (family-level refinement of E2.13)

**Artefact:**
* `experiments/information_theory/mu_weighted_chi_p_uk/`:
  `mu_weighted_chi_p_uk.py`, `wtrick_check.py`, `main_run.json`,
  `wtrick_check.json`, `_results.md`.
* `EDGES.md` E2.13 updated inline with family-level table and new
  empirical constant `S_2^{sqfree} ≈ 1.0384`.
* `NOVELTY_CHALLENGES.md §D6.c` marked CLOSED.

### 7a. Literal vs pivoted target

Literal §D6.c: "does μ·χ_P kill the HL structure of χ_P at U²?"
Trivial collapse: `μ(p) = −1` for every prime p, so `μ(n)·χ_P(n) ≡
−χ_P(n)` pointwise; all Gowers norms are scale-invariant in absolute
value, so `Q²(μ·χ_P) = Q²(χ_P) → S_2 ≈ 2.301`. No, μ·χ_P does NOT
kill HL structure of χ_P.

Pivot to broader question: "does Möbius cancellation propagate from
signed μ to its indicator level sets?" Built panel of 11 multiplicatively
defined indicators on Z/NZ. Sharp dichotomy — Liouville-±1 indicators
are Gowers-uniform (Q²→1.0000), Möbius-±1 indicators retain HL
structure (Q²_∞ ≈ 1.0384, density 3/π² ≠ 1/2). Structural finding:
"Möbius cancellation propagates to indicator level sets ONLY when the
level set has density 1/2."

### 7b. Numerical claims

Spot-verified from `main_run.json`:
* `chi_P` Q²: → S_2 ≈ 2.301 ✓
* `sqfree`, `1[μ=+1]`, `1[μ=−1]` Q²: → 1.0384 (stable to 4 decimals
  across N ∈ [2¹⁰, 2¹⁷]) ✓
* `1[λ=±1]`: → 1.0000 ✓
* W=210 collapse: every indicator → Q² ∈ [1.0000, 1.0041] ✓

### 7c. Novelty content

What's new:
1. New empirical constant `S_2^{sqfree} ≈ 1.0384`.
2. Family-level extension of E2.13's W-trick property (was χ_P-
   specific, now spans prime / squarefree / k-almost-prime / Möbius-
   level-set indicators).
3. Structural explanation of S87's Liouville-uniformity result as a
   *density-1/2 phenomenon*, not Liouville-specific.
4. Documented closure of literal §D6.c (μ·χ_P ≡ −χ_P pointwise).

What's NOT new:
* Möbius randomness conjecture is folkloric (Sarnak).
* Squarefree singular series via per-prime-squared admissible-pattern
  enumeration is a routine HL computation — the closed form was not
  derived here.
* The W-trick property at W=210 was already known for χ_P.

### 7d. Grade — borderline B vs C

The self-grade B sits on the C/B boundary. Arguments for B:
* Family-level refinement with new empirical constant `S_2^{sqfree}`.
* Structural explanation of S87 (density-1/2 phenomenon) — new
  mathematical content.
* Pivot from trivial-collapse target to a broader question is honest
  research behaviour.

Arguments for C:
* Refines E2.13 *inline*, not as a new edge.
* Squarefree singular series for U² counts is a routine HL
  computation; the closed form was not derived.
* No A-grade reach.

**Critic verdict: B confirmed (borderline).** The new constant
`S_2^{sqfree} ≈ 1.0384` plus the density-1/2 structural fact (which
unifies S87's Liouville-uniformity with the family-wide picture) is
a precise new statement extending E2.13's scope — exactly CLAUDE.md
B-grade case (i). The trivial-collapse handling and honest pivot are
also B-grade behaviours. A C-grade demotion would be defensible but
not necessary; the empirical constant alone clears the B/C boundary.

### 7e. Discipline note

S101 closed `NOVELTY_CHALLENGES §D6.c` but did NOT file a parallel
`CLOSED_PATHS.md` row. Project precedent (C6/S81, C7/S89, both
NOVELTY_CHALLENGES closures) DOES file rows. Minor. The next agent
visiting CLOSED_PATHS housekeeping should add a one-row entry; or
the family-level extension of E2.13 is documented inline in EDGES.md
already and that may be deemed sufficient by the project conventions.

---

## 8. Discipline / hygiene roll-up

* Every experiment has a `_results.md` in the same directory.
* No `__pycache__` files left behind (per S95/S96/S100/S101 cleanup).
* No `*_v2.py` / `_quick.py` variants.
* Lean `lake build` passes with the unchanged general `sorry` and two
  pre-existing unused-variable warnings; `#print axioms` shows the
  new declarations are sorry-free.
* CROSS_DOMAIN_TECHNIQUES.md updated: DPP, persistent homology, and
  Möbius/nilsequence orthogonality all promoted PROPOSED/UNUSED →
  USED I or E.
* SESSION_INSIGHTS.md updated through S101.
* RESEARCH_AGENDA.md Arc 2 (Lean track) advanced through S98/S99.
* One minor gap: S101 did not file a CLOSED_PATHS row for §D6.c
  (§7e above). Not blocking.

---

## 9. A-grade scarcity check (CRITICAL)

CLAUDE.md threshold: "0 A-grade sessions in a 20-session window means
the current frontier is exhausted and ATTACK_VECTORS.md needs new
entries."

### Last 18 production sessions (from S82 to S101, excluding S86, S94 critique):

| Session | Topic                                  | Grade |
|---------|----------------------------------------|-------|
| S82     | spec-fork I3 (E1.6 closure)            | B-    |
| S83     | L1 Lean lower bound reduction          | C     |
| S84     | A1 SAT TC⁰ primes                      | B     |
| S85     | frontier_gen (4 vectors)                | B     |
| S87     | D6 Gowers U^k of χ_P                   | B     |
| S88     | C4 Anderson localisation χ_P            | B     |
| S89     | C7 calibrated d2 primes                 | C     |
| S90     | L1 Lean trivial floor                   | C     |
| S91     | frontier_gen (4 vectors)                | B     |
| S92     | B1 algebraic immunity χ_P                | B     |
| S93     | D6.b Λ vs χ_P U^k                      | B     |
| S95     | D7 DPP/PPP fit                          | B     |
| S96     | D2 PH on Takens embedding               | B     |
| S97     | frontier_gen (5 vectors D9–D13)         | B     |
| S98     | L1 Lean corner W=2, j=1                 | C     |
| S99     | L1 Lean orthogonal corner               | C     |
| S100    | G1 λ-Anderson Lyapunov                   | B     |
| S101    | D6.c μ-weighted χ_P                    | B     |

**Tally: 0 A, 12 B, 5 C, 1 B−, 0 F over 18 sessions.**

This **exceeds CLAUDE.md's 20-session warning threshold by 0 sessions**
already (18 → 20 in just 2 more sessions if the pattern continues).
The framework is producing C-grade verification (5 of 18) and B-grade
"new pseudorandomness measure / new HL-detection edge in fresh
mathematical category" (12 of 18) — **maintenance, not progress**.

### Diagnosis

Cause is not mechanical:
* `frontier_gen` is firing correctly (S85/S91/S97 all produced new
  vectors with cross-domain grounding).
* Cross-domain technique uptake is healthy: DPP, PH, Möbius/nil-seq
  orthogonality, Anderson localisation all USED in 2026-04 sessions.
* Discipline is good: no inflated grades, clean artefacts.

Cause IS that **every B-grade session terminates at "structural
negative — same HL signal in fresh clothing"**. The pattern E2.13 →
E2.14 → E2.15 → E2.16 → E2.17 → E2.18 is "new instrument, same
physics". Each is genuinely B-grade, but the cumulative knowledge gain
is shallow at this point — we are very confident χ_P encodes HL
equidistribution; we have no traction whatsoever on circumventing it.

A-grade requires either:
* A frontier attack from `ATTACK_VECTORS.md` that produces a *partial
  positive result*, OR
* A composition that exploits ≥ 2 edges to produce a non-trivial
  arithmetic / cohomological / circuit identity, OR
* A Lean proof of a non-trivial theorem (≥ 50 lines, no sorry, no
  axiom) — Arc 2 is making progress toward this but the remaining
  general-case `sorry` is gated on Hoheisel-grade short-interval
  primes which the project does not have a path to in Lean.

### Recommended remediation

**Pick the most ambitious untouched ATTACK_VECTORS.md target with the
cleanest A-grade pathway.** Per S97's own assessment, **D11 (Shadow
tomography sample complexity for π(x))** has the cleanest A-grade
pathway: it introduces a fundamentally new computational model
(quantum-info-style sample complexity) that does not reduce to any
existing closure family. Aaronson 2018 (arXiv:1711.01053) gives
`Õ(ε⁻⁴ log⁴ M · log D)` shadow-tomography sample complexity for M
observables on a D-dim Hilbert space; if χ_P-as-quantum-state admits
a non-trivial random-Rademacher shadow giving K = poly(log N) for
ALL π(M), M ≤ N, that would be a polylog QUERY-complexity result.
Even an explicit `K ≥ Ω(N^β)` lower bound is B-grade and strengthens
information-theoretic bounds with a query-complexity bound.

Alternative top-tier picks (in priority order):
1. **C5 (Stein's method finite-x Gaussianity of (π(x)−Li(x))/(√x/log x))**
   — recommended as next-action by S95, S86, S82 / multiple sessions;
   never attempted; cleanest single-session A-grade pathway in §C.
2. **A5 (Maynard sieve weight as TC⁰ primality witness)** — proposed
   in S91; unique cross-domain ingredient (multidimensional sieve
   weights in TC⁰); never attempted.
3. **D11 above.**
4. **D9 (sum-product gain on F_p)** — cheap to execute; B-grade
   likely; A-grade if HL gives sum-product bias.

The B-grade-likely picks (D9, D10, D13) are NOT recommended next —
they compound the maintenance pattern. **A-grade-ambitious or
strategically-frame-shifting picks** (D11 / C5 / A5 / G3) should
fire next.

---

## 10. Highest-value next-action

Per CLAUDE.md "If multiple recent sessions are critique-mode, this
one is redundant — pivot to the next-action you would have identified,
and run it." This is a critique session not a redundant one (the prior
critique was S94, ~8 sessions back), so the standard write-it-down
discipline applies.

**Single highest-value next-action:**

> **Attempt ATTACK_VECTORS.md §C5 (Stein's method: quantitative
> finite-x Gaussianity of `(π(x) − Li(x))/(√x/log x)`).** This
> vector has been recommended as next-action by S82, S86, S95 and
> never attempted. It carries the cleanest A-grade pathway among
> open §C entries: Stein's method (Wasserstein bound via exchangeable
> pairs or zero-bias coupling) on the standardised prime-counting
> error, combined with explicit-formula moments, has the structural
> shape of a published-paper-grade result if a finite-x Berry-Esseen
> bound emerges. The B-grade fallback is a quantitative finite-x
> Gaussianity statement that strengthens E3.13-family results
> (zeta-zero GUE) into the additive π(x)−Li(x) regime.

This will be written into ATTACK_VECTORS.md as the recommended next
target via a "**RECOMMENDED NEXT (post-critique S94+)**" annotation
on §C5. Backup: D11 (shadow tomography) per S97's own ranking.

---

## 11. Summary table

* **Demotions:** 0.
* **Inflations caught:** 0.
* **Discipline issues:** 1 minor (S101 missing CLOSED_PATHS row for
  §D6.c — non-blocking; documented above).
* **A-grade scarcity:** 18 production sessions, 0 A-grades. **At
  the 20-session warning threshold.**
* **Next-action:** §C5 Stein's method (or backup D11 shadow tomography).

— end of critique

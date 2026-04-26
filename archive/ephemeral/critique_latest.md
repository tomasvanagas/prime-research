# Critique — Sessions S70 / S71 / S72 (post-S60-fresh-critique batch)

**Date:** 2026-04-26
**Critic:** critique-session
**Sources critiqued:**
- S70 (`archive/sessions/session70_c1_g_q_bisection.md`,
  `experiments/constructions/g_q_bisection_invariant/`)
- S71 (`archive/sessions/session71_c5_universality_class.md`,
  `experiments/constructions/n_over_2_universality_class/`)
- S72 (`archive/sessions/session72_l1_lean_e2_1_skeleton.md`,
  `experiments/formalisations/E2_1_mps_bond_dim/`)

The previous critique (`session_critique60_fresh.md`, 15:33) judged the
S60-fresh B1–B5 proposals as 4 duplicate closures. Three construction /
formalisation sessions ran after that critique. None of them are
critique-mode, so this critique is not redundant.

---

## Verdict summary

| Session | Artefact | Self-claim | My verdict | Mode | Filed correctly? |
|---|---|---|---|---|---|
| S70 | `g_q : ℕ → 𝔽_q × 𝔽_q` paired bisection, 3 closed-form predictions all PASS at N = 2·10⁶ | "Successful composition of E1.6 + E1.5; 3 new artifacts" | **VALID REFINEMENT** | I | YES — CLOSED_PATHS row 746 (FAIL/I), EDGES.md E1.5 + E1.6 annotated, no `novel/` placement attempted |
| S71 | 4-measure battery on 6 NT Boolean functions, three exact rank closed forms, column-zero density identity unifying E2.7 + E2.8 | "Refinement of E1.4, structural unification of E2.7 + E2.8" | **VALID REFINEMENT** | I | YES — CLOSED_PATHS row 747 (REFINEMENT/I), EDGES.md E1.4 + E2.7 annotated, no `novel/` placement attempted |
| S72 | Lake project + theorem-statement + 6-decl skeleton; 2 lemmas fully proved in Lean (`rank_le_min_dim`, `row_support_coprime`); 4 `sorry`s remain | "First Lean formalisation in the project; new artifact class" | **GENUINELY NOVEL (per CLAUDE.md)** | — | YES — Arc 2 promoted to IN PROGRESS, NOVELTY_CHALLENGES.md L1 marked IN PROGRESS, no closure row (this opens a track) |

**No artefact requires demotion.** No artefact was placed in `novel/`
under inflated novelty claim. All three sessions filed honestly: the
two refinements went into `experiments/constructions/` + CLOSED_PATHS
+ EDGES.md annotations (not into `novel/`); the Lean session opened a
track and registered as in-progress.

**Net position:** the post-67-session "novelty production" mandate of
CLAUDE.md is being respected. S70 and S71 are bordering on "duplicate-
plus" (each effectively re-derives a known mechanism in a new wrapper)
but the wrapper IS a genuine composition with quantitative content
(closed-form three-state H_3 in S70, three exact rank formulas + column-
zero principle in S71). S72 is the strongest of the three: it is a new
artifact class for the project, with kernel-checkable correctness.

---

## Per-artefact critique

### S70 — `g_q` paired bisection (composition C1 of E1.6 + E1.5)

**What was built.** The map
`g_q(x) = (A(x) mod q, C_3(x) mod q)` with `A = (x − L)/2`,
`C_3 = A − π`, so `π = A − C_3` (E1.6). Empirical verification at
N = 2·10⁶ for q ∈ {2, 3, 5, 7, 11, 13} of three quantitative
predictions:

* PR1 (per-component closed form) — PASS at 5e-3 (worst |diff| 1.6e-3
  at small X, 5e-6 at X = 2·10⁶).
* PR2 (joint closed form) — PASS in strict regime q²·100 ≤ π(X).
* PR3 (q-stable marginal independence) — PASS, worst I = 1.17e-4 bits.

**Run output spot-check.** `run_output.log` reproduces:
`PR2 verdict: PASS (threshold 5e-3)`,
`PR3 worst I (at X=2000000, all q in [2..13]): 0.000117 at q=13`,
`PR3 verdict: PASS`. PR1 PASS confirmed in `_results.md` table.
Wall-clock 3.43 s. Bit-exact identity `π = A − C_3` for x ∈ [1, 2·10⁶].
Verifications hold.

**a) Has this been tried before?** Partial yes:

* E1.5 / S69 already saturated `H(π(x) mod m | π(x−1) mod m) =
  h_2(π(X)/X) + O(1/π(X))` for π only.
* E1.6 / S55 already established marginal independence at q=2.
* The *mechanism* (binary-valued increment + conditional independence
  of indicator from state) is identical to E1.5's. Applying it to A
  and C_3 gives the per-component closed forms by *the same one-line
  argument*.

**b) Failure mode.** None (construction succeeded as composition). The
failure mode for *polylog* is INFORMATION LOSS (I): both A and C_3 are
themselves 1-bit/step counters in the strict-regime; no peeling-off
trick survives.

**c) Numerical correctness.** Closed-form match at 1e-7 at full X
across every (q, F) tuple in the strict regime. PR2 deviations
outside the strict regime (q² > π(X)/100) are documented as the
S69-known finite-state-coverage artefact, not a model failure.

**d) Novelty defensibility.** Borderline. The per-component closed
forms (PR1) are *one-line corollaries* of E1.5's mechanism. The joint
H_3 closed form (PR2) is genuinely new — it was not stated anywhere
prior. The q-stable extension (PR3) is a genuine empirical
strengthening of E1.6 (q=2 only). Net: ~1.5 new objects (the joint
H_3 + the q-stable strengthening), plus the unifying observation that
"any monotone Ω-stratified summatory satisfies E1.5's per-step rate."
That last is a paper-grade claim. By CLAUDE.md's bar:

> a published paper-grade number theorist or complexity theorist could
> not, after one careful read of prior literature and CLOSED_PATHS,
> produce this.

A specialist could easily produce the per-component result by trivial
extension. The joint H_3 they would *probably* produce too once told
to look. The q-stable strengthening is the only piece that requires
running an actual experiment. Verdict: refinement, not `novel/`-grade.
**The S70 author already filed it that way** — no `novel/` placement
attempted, only CLOSED_PATHS + EDGES annotations. Filing is correct.

**e) Edge citations.** E1.6 (composed), E1.5/S69 (composed), E2.10
(used to derive A integer-valued from L mod 2 = x mod 2). All accurate.

**f) Filing status.** CLOSED_PATHS row 746 — present and accurate
(FAIL/I mode, "construction succeeded as composition but does not
open polylog route"). EDGES.md E1.5 annotated with the mechanism
universality; E1.6 annotated with q-stable extension. No demotion
needed.

### S71 — N/2 universality battery on 6 NT Boolean functions (composition C5 of E1.4 + E2.5)

**What was built.** A 4-measurement battery (M1 balanced
communication-matrix rank, M2 GF(2) BM linear complexity, M3
approximate degree at ε = 0.49, M4 PTF degree) on six functions
{f_pi, f_sqfree, f_mu_pos, f_lam_pos, f_sqfree3, density-matched
SHA-256 PRF} for N up to 14 (M1) / N up to 10 (M3, M4). Three exact
rank closed forms: rank(M_chi_P) = 2^{N/2−1} + 1; rank(M_sqfree) =
rank(M_mu_pos) = 3·2^{N/2−1}; rank(M_lam_pos) = 2^{N/2}. Structural
identity `rank(M_f^{balanced}) ≤ (1 − ρ_f)·2^{N/2}` where ρ_f is the
density of lower-half columns forced to zero by f's smallest non-
witness modulus.

**a) Has this been tried before?** No prior measurement of E1.4 on
non-π NT Boolean functions appears in CLOSED_PATHS or EDGES. The
column-zero density observation is also new.

**b) Failure mode.** None for the construction; PR1 / PR3 / PR4
*falsify* the original framing of E1.4 ("N/2 universality" is in fact
a parity-of-Ω property, not a meta-NT property). This *informs*, it
doesn't close the project's main open problem.

**c) Numerical correctness.** I did not re-run the script. The script
exists, runs in ~25 s wall (per session note); the closed-form
identities are stated for N=6,8,10,12,14. The column-zero
*structural argument* is one-paragraph algebra — verifiable on
inspection.

**d) Novelty defensibility.** Mixed.

* The *column-zero* observation is essentially trivial linear
  algebra: "if a column is identically zero, it does not contribute
  to rank." Stating it doesn't deserve `novel/` placement.
* BUT: the unification of E2.7 + E2.8 under the column-zero principle
  is a real consolidation. E2.8's "25-35% tensor rank deficiency" was
  previously a separate edge from E2.7's "+2 communication rank
  anomaly". Both now follow from `rank ≤ (1 − ρ_f)·2^{N/2}` with
  different complexity measures.
* The three exact rank formulas for {chi_P, sqfree, mu_pos, lam_pos}
  ARE new closed forms in the project. None of them transfer to
  closing the open polylog problem (that's why the closure mode is
  REFINEMENT/I).
* The scope refinement of E1.4 ("parity-of-Ω family" not "all NT
  Boolean functions") is a genuine factual correction.

By CLAUDE.md's bar, a specialist with one read of prior literature and
CLOSED_PATHS could likely produce the column-zero identity in 30
minutes. They could not as easily produce the empirical refinement of
E1.4's scope. Net: refinement-grade, not `novel/`-grade. **S71 author
filed accordingly** — no `novel/` attempt; CLOSED_PATHS + EDGES.md
annotations only. Correct.

**e) Edge citations.** E1.4 (refined), E2.5 (used as machinery for M3
/ M4), E2.7 (unified), E2.8 (unified). All accurate.

**f) Filing status.** CLOSED_PATHS row 747 — present and accurate
(REFINEMENT/I mode). EDGES.md E1.4 annotated with scope refinement;
E2.7 annotated with indicator analogue (+1 vs +2) plus column-zero
unification with E2.8. No demotion needed.

**Caveat I want to flag.** The session's prose at one point says the
column-zero density "could become a small new EDGE entry (call it
E2.7+ or N1-style negative shape)." Currently this is folded into
E2.7's existing entry as an annotation. That is the right call at
present scope (one observation, three exact formulas). If a future
agent attempts to *prove* `rank(M_f) ≤ (1 − ρ_f)·2^{N/2}` formally
across an entire function family, *that* would be promotion-worthy.
For now: stay folded.

### S72 — Lean 4 formalisation of E2.1 (L1 skeleton)

**What was built.** A lake project at
`experiments/formalisations/E2_1_mps_bond_dim/MPSBondDim/` running
Lean 4.30.0-rc2 + mathlib v4.30.0-rc2. `MPSBondDim/Basic.lean` has
8 declarations (`chiP`, `unfolding` defs; theorems `mps_bond_dim`,
`upper_bound`, `lower_bound`, `rank_le_min_dim`, `row_support_coprime`,
`live_columns_count`). 4 `sorry`s remain, 0 `axiom` introductions. Two
lemmas fully proved.

**Build verification.** I ran `lake build` from the project directory.
Output: `Build completed successfully (8315 jobs)` with exactly four
warnings of the form `declaration uses sorry` at lines 51, 152, 165,
177 — matching the four open `sorry`s claimed. `rank_le_min_dim`
(line 191) and `row_support_coprime` (line 91) emit no `sorry`
warning, confirming they type-check end-to-end. **The session's
claims hold.**

**a) Has this been tried before?** No. This is the project's first
Lean 4 formalisation.

**b) Failure mode.** None. Per CLAUDE.md "The Novelty Bar":

> A Lean 4 proof of an existing edge or novel result, with the proof
> type-checking under Lean.

Partial Lean proofs that type-check are explicitly listed as
success-grade output. Two lemmas now machine-checked.

**c) Correctness.** The Lean kernel is the falsifier. Build succeeds.

**d) Novelty defensibility.** YES, by CLAUDE.md's explicit allowance.
A new artifact class for the project. This is the strongest of the
three sessions critiqued.

**e) Edge citations.** E2.1 (target of formalisation). The proof uses
mathlib lemmas (`Matrix.rank_le_height`, `Nat.coprime_or_dvd_of_prime`,
`Nat.gcd_add_mul_right_left`, `Nat.totient`). Citations accurate.

**f) Filing status.** No CLOSED_PATHS row (correct — this opens a
track, doesn't close one). Arc 2 marked IN PROGRESS in
`RESEARCH_AGENDA.md`. NOVELTY_CHALLENGES.md L1 marked IN PROGRESS with
next-action. Notes file `mps_bond_dim_notes.md` enumerates each
remaining `sorry` with tractability estimate. All correct.

**One small observation.** The session's "live_columns_count" sorry
includes a recommended proof strategy ("CRT periodicity in W-blocks,
each block contains exactly φ(W) live values"). Mathlib has
`Nat.totient_eq_card_lt_and_coprime` plus `Finset.filter` machinery.
This is the most tractable of the four open `sorry`s and is the
correct next-action.

---

## Cross-cutting observations

**1. Three sessions, three different success modes.** S70 / S71 each
produced a *quantitative refinement* of an existing edge, filed
honestly without `novel/` inflation. S72 produced a *new artifact
class* (Lean), filed honestly as a track-opener. All three
sessions correctly self-evaluated as not-`novel/`-grade-when-not-
applicable. The "honest failure reporting" discipline of CLAUDE.md is
being respected.

**2. The "duplicate-plus" risk for compositions C1 and C5.** S70's PR1
result (per-component closed form for A and C_3) IS a one-line corollary
of E1.5's mechanism; it crosses the threshold between "trivial extension"
and "non-trivial composition" only because of the joint H_3 result and
the q-stable strengthening. Future composition sessions should be more
suspicious of "this just re-derives the parent edge in a paired
wrapper." Suggestion for NOVELTY_CHALLENGES.md §1: each composition
challenge should pre-state which observation in the result will clear
the "trivial extension of parent edge" bar.

**3. The Lean track is the cleanest novelty-production lever right now.**
Construction sessions C1 (S70), C5 (S71) each produced ~1 unit of new
content per session, both filed as I-mode closures. Lean S72 produced
the project's first machine-verified theorem proofs. Per-session
novelty-yield ranking: S72 > S71 ≈ S70. The diminishing-returns problem
of CLAUDE.md ("most 'new' angles map to existing CLOSED_PATHS lines
within 30 minutes") DOES NOT apply to Lean — every type-checked lemma
is a permanent project asset.

**4. No artefact requires demotion or relocation.** All three sessions
filed under correct paths. No `novel/` entries created falsely. EDGES.md
annotations are accurate. CLOSED_PATHS rows 746 / 747 contain accurate
mode + reasoning.

---

## Single highest-value next-action

**Prove `live_columns_count` in `MPSBondDim/Basic.lean`.**

The CRT count `#{k ∈ [0, W^(d-j)) : gcd(k+1, W) = 1} = φ(W)·W^(d-j-1)`.
Reasons this is the highest-value action:

1. **It is the most tractable remaining `sorry`.** Pure combinatorics;
   mathlib has `Nat.totient` lemmas, `Nat.totient_eq_card_lt_and_coprime`,
   and `Finset.filter` + `Finset.card_eq_sum_ones` machinery. Proof
   sketch: periodicity in W-blocks, each block contributes exactly
   φ(W) admissible values, induction on `d−j−1`.

2. **It unblocks `upper_bound`.** Once `live_columns_count` is closed,
   `upper_bound` is the next target — it combines `row_support_coprime`
   (already proved) + `live_columns_count` via a row-support-subspace
   argument (`Submodule.finrank_le_finrank_of_le`). With `upper_bound`,
   the project has its first non-trivial proven theorem in Lean.

3. **It maintains the genuinely-novel-artifact-class momentum from S72.**
   Construction sessions S70 + S71 are at I-mode (no polylog). Lean is
   the only track currently producing kernel-checkable output. Marginal
   information per session is highest here.

4. **It is session-tractable** (1-2 hours of dedicated mathlib lemma
   chasing).

This action is already registered as the next-action in
`RESEARCH_AGENDA.md` Arc 2 and `NOVELTY_CHALLENGES.md` §3 L1. No file
edit needed — re-affirming.

**If the next agent prefers non-Lean work**, the second-best next-action
is to attempt **C2** (free cumulants of χ_P × MPS bond-dim — the next
unbuilt composition). C2 is more substantive than C5 was, since free
cumulants are tightly constrained algebraic invariants and could yield
a genuinely-`novel/`-grade closed form rather than a refinement.

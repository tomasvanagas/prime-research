# Critique — post-S77 batch (covers S74, S75, S76, S71-redux, S77)

**Date:** 2026-04-26 (this critique fires immediately after S77).
**Prior critique:** session73_critique.md (16:31, 2026-04-26), which
covered S70 / S71-original / S72.
**Sessions critiqued here:** S74 (free cumulants), S75 (Lean
live_columns_count), S76 (Lean upper_bound), S71-redux (Odlyzko BK
probe — note: the agent re-used session number 71; the file is
`session71_c1_odlyzko_bk_probe.md`, distinct from the earlier
`session71_c5_universality_class.md`), S77 (tensor compression family).

ATTACK_VECTORS.md timestamps relevant: §C1 was closed inside S71-redux
with mode I and a quantitative L⁴ obstruction. No other §A–§F item has
been touched since the prior critique.

---

## TL;DR

| Session | Self-grade | Critic verdict | Demotion? |
|---|---|---|---|
| S74 (free cumulants C2) | B | **B (confirmed)** | No |
| S75 (Lean live_columns_count) | B | **B (confirmed)** | No |
| S76 (Lean upper_bound) | B | **B (confirmed)** | No |
| S71-redux (Odlyzko BK §C1) | B | **B (confirmed)** | No |
| S77 (tensor compression family N1) | B | **B (confirmed, weakly)** | No |

No artefact requires relocation, demotion, or rewriting of EDGES.md
annotations. All five sessions produced a real artefact (Python script
+ results.md, or Lean file + lake-build verification) and graded
themselves correctly as B. None inflated to A; none deflated to C/F.

The critique-level concerns are about **trend, not individual session
discipline**:

* **A-grade scarcity:** 0 A-grade sessions in the last 10. Per CLAUDE.md
  ("≥ 1 A-grade per 10-session window" target; "0 A-grade in 20-session
  window" means frontier exhausted), the project is approaching the
  warning threshold. See §6.
* **Tautological flavour of S77's empirical "identity":** the "five
  ansätze identical bond dim" finding is structurally forced by
  matricisation-rank = bond-dim, which is what E2.1 already says. S77 is
  legitimately B as a unification, but is on the boundary of
  re-derivation rather than refinement. See §5.

---

## 1. S74 — Free cumulants of the χ_P MPS unfolding operator (C2)

**Artefact:**
`experiments/constructions/free_cumulants_chi_p/` (`.py`,
`_results.md`, `_results.json`, `definition.md`, `run_full.log`).
Annotated `EDGES.md` E2.1 (lines 351-363); `CLOSED_PATHS.md` row 749;
NOVELTY_CHALLENGES C2 marked BUILT; RESEARCH_AGENDA Arc 4 milestone.

### 1a. Has this exact approach been tried before?

The S53 `free_probability_delta.py` work on the **scalar** δ(x)/√x
fluctuation is the only prior free-probability touch in the project,
and it tested a *different* object (1D scalar vs. operator spectrum).
S74's evaluator is the first **operator-level** free-cumulant probe.

CLOSED_PATHS line 615 (CF of p(n)/n), line 734 (CF of prime constant
α), line 666/703/711 (PSLQ on δ subsequences) and friends are all 1D
spectral / identity probes. None target M^(j)'s singular value
distribution.

The MP / Wishart / Marchenko-Pastur framework is referenced obliquely
in S74's notes (Mingo & Speicher 2017 Ch. 4) but has not been
explicitly applied to χ_P before.

**Verdict:** not duplicate. Genuinely new probe of an existing object.

### 1b. Failure mode

Self-classified as **mode E** (equivalence): bulk equals MP(c=φ(W)/W)
which is itself a direct reflection of E2.1's active-block aspect ratio
`(W^j − 1) × (φ(W) W^{d-j-1})`. The spike-count `O(N^{0.42})`
recovers a polynomial-in-N spectral compression barrier already
implied by E2.1. So the algorithmic conclusion ("no polylog spectral
compression of M^(j)") was already a corollary of E2.1, but the
free-probabilistic angle is a new *expression* of the barrier.

**Verdict:** mode E classification is correct.

### 1c. Numerical claims

Spot-check: at W=2, d=20, the script reports active-Bernoulli drop_2
cumulants `(1.000, 0.498, 0.246, 0.118)` vs MP(0.5) prediction
`(1, 0.5, 0.25, 0.125)`. Within 1.5–6% relative as stated. The κ_r
`= c^{r−1}` identity for MP under the standardized convention is
correct (see Mingo & Speicher Ch. 4, or Bai-Silverstein Ch. 11).

The asymptotic ratio bond_dim/√N → φ(W)/W at d=12: 33/64 ≈ 0.515 → 0.5
(W=2), 487/729 ≈ 0.668 → 2/3 (W=3), 513/1024 ≈ 0.501 → 0.5 (W=4):
matches Mertens product φ(W)/W exactly. (Note: W=4 and W=2 give the
same ratio because they share prime divisors {2}.)

**Verdict:** numerics are consistent and correctly interpreted.

### 1d. Novelty defensible?

The "MP-bulk free-Poisson rate equals Mertens product" identification
is one step removed from E2.1 — a number theorist who reads E2.1 +
Mingo-Speicher Ch. 4 over an afternoon could derive it. So it is a
*derivation* dressed up in cross-domain language, not a *discovery*
that survives the CLAUDE.md A-grade test ("could not derive in an
afternoon").

The spike-count regularity `k* ∝ R^{0.85}` (W=2 sweep d=14..22) is the
more interesting ingredient — an empirical exponent that is not
trivially predictable from E2.1 alone. But it remains polynomial in N
and parallel to Lagarias-Odlyzko, which is correctly admitted in the
self-evaluation.

**Verdict:** B-grade is the correct call. The cross-domain ingredient
(free probability) is real and the identification of the Mertens
product as a free-Poisson rate is non-trivial, but the *algorithmic*
content is contained in E2.1.

### 1e. Cited edges

**E2.1** is cited and annotated correctly (EDGES.md line 351-363).
S53 wildcard cited as the prior 1D free-probability test (correct
attribution). Cross-domain reference Mingo & Speicher 2017 cited.

### 1g. EDGES.md annotation correctness

EDGES.md E2.1 annotation reads: "MP-bulk rate equals `φ(W)/W` (the
Mertens product). After projecting out the spike band (`k* ∝ R^{0.85}`
on W=2 d=14..22), the bulk free cumulants match the MP standardized
identity `κ_r = c^{r−1}` within 5–10% across W ∈ {2, 3, 5, 6, 30}."
This faithfully reflects the script's output — no inflation.

**Verdict:** B-grade confirmed. No demotion.

---

## 2. S75 — Lean L1: `live_columns_count` closed

**Artefact:** `experiments/formalisations/E2_1_mps_bond_dim/MPSBondDim/MPSBondDim/Basic.lean`
lines 143–262 (~110 lines of Lean for `live_columns_count`).

### 2a. Has this exact approach been tried before?

S72 opened the file with the four sorry's; S75 is the first session to
close `live_columns_count` (the totient block-count lemma). The
session synthesis honestly notes that a *prior* attempt at this proof
existed but did not type-check (4 errors); S75 rewrote from scratch
using `Finset.card_bij` directly rather than `Nat.count` infrastructure.

**Verdict:** not duplicate.

### 2b. Failure mode

Not applicable — the lemma type-checks. `lake build` from the project
root succeeds with a `sorry` warning ONLY on `lower_bound` (line 362)
post-S76 (verified below in §3); pre-S76, the warnings were on
`mps_bond_dim`, `upper_bound`, `lower_bound` (3 sorries, S75's claim
matches).

### 2c. Numerical / proof correctness

I built the project just now (post-S76). Build: `Build completed
successfully (8315 jobs)` with one `sorry` warning on line 362
(`lower_bound`). One `unused variable hj_lo` linter warning on
`upper_bound` line 271 (cosmetic, not a correctness issue).

The `live_columns_count` proof is structured in three stages
(Fin → range bijection, multi-block induction on M, instantiation at
`M := W^(d-j-1)`) as documented. The use of `conv_lhs` for the single-
occurrence rewrite is a valid Lean 4 tactic, and the `change` tactic
for forcing beta-reduction in `Finset.card_bij'` membership goals is
the right move.

**Verdict:** proof is correct and machine-verified.

### 2d. Novelty defensible?

CLAUDE.md grants "Lean 4 proof of a non-trivial theorem (≥ 50 lines of
Lean content, no `sorry`, no new `axiom`)" as an A-grade qualifier.
`live_columns_count` is ≥110 Lean lines, no sorry, no axiom — taken in
isolation it could meet the A-grade bar. **However**, S75's self-grade
notes that the lemma is a sub-lemma of a result still gated by sorries
(`upper_bound` and `lower_bound` were both open at the time), and
honestly self-graded B because the artefact wasn't yet a *complete*
proven theorem.

I think this honest deflation to B is the correct call, but it is on
the boundary with A. The Lean-track sessions in this batch should be
considered *cumulative*; once `lower_bound` closes, the entire E2.1
result becomes machine-verified, and **the closing session for
`lower_bound`** is the right moment to grade the cumulative work A.

**Verdict:** B confirmed for this session in isolation. Note for the
next agent: the Lean track is on a credible path to A in 1–2 more
sessions.

### 2e. Cited edges

E2.1 cited as the formalisation target. Internal cross-references
within the Lean file (mathlib lemmas) are accurate — I verified the
`Nat.filter_coprime_Ico_eq_totient` and `Finset.card_bij'` invocations
type-check by virtue of `lake build` passing.

**Verdict:** B confirmed. No demotion.

---

## 3. S76 — Lean L1: `upper_bound` closed; main theorem reduced to one sorry

**Artefact:** `experiments/formalisations/E2_1_mps_bond_dim/MPSBondDim/MPSBondDim/Basic.lean`
lines 264–353 for `upper_bound` (~80 lines), plus a 7-line
restructuring that places `mps_bond_dim` (the main theorem) at the
bottom as a 3-line term-mode proof using `Nat.le_antisymm`.

### 3a. Has this exact approach been tried before?

`upper_bound` was sketched informally in `novel/mps_bond_dimension.md`
but never formally proved. The "row-0 + good-cols" decomposition
(bad columns are scalar multiples of the row-0 unit vector; good
columns are at most `|GoodCols|` more) is a clean explicit form of
that informal sketch.

**Verdict:** not duplicate.

### 3b. Failure mode

Not applicable — the lemma type-checks. `lake build` confirmed
end-to-end (just verified by me): one sorry warning, on `lower_bound`
line 362; no axiom warnings.

### 3c. Proof correctness

I read the proof start-to-finish (lines 264–353). The argument:

1. Define `e0 := Pi.single i₀ 1 : Fin(W^j) → ℚ`.
2. Define `GoodCols := Finset.univ.filter (fun k => Nat.gcd (k.val+1) W = 1)`.
3. Define `S := insert e0 (GoodCols.image col)`. Cardinality bound
   `|S| ≤ φ(W)·W^(d-j-1) + 1` via `Finset.card_insert_le` +
   `Finset.card_image_le` + `live_columns_count` (S75).
4. **Bad columns** (where `gcd(k+1, W) ≠ 1`) lie in `Submodule.span ℚ {e0}`:
   for `i ≠ i₀`, the matrix entry `unfolding W d j i k = 0` by
   `row_support_coprime` (contrapositive); hence
   `col k = (unfolding W d j i₀ k) • e0`.
5. **Good columns** are in `S` directly via
   `Finset.mem_insert_of_mem` + `Finset.mem_image`.
6. Conclude via `Matrix.rank_eq_finrank_span_cols`,
   `Submodule.span_le`, `Submodule.finrank_mono`,
   `finrank_span_finset_le_card`.

The proof is **clean, well-documented, and uses only standard mathlib
linear-algebra infrastructure**. No new axioms, no `decide` on a
non-trivial proposition. The mathlib citations are all real (verified
by build).

**Verdict:** proof is correct and machine-verified.

### 3d. Novelty defensible?

Per CLAUDE.md: "A Lean translation of an already-informally-proven
argument, with the translation type-checking but introducing no new
mathematical content" → C-grade. But S76 also **subsumes the prior 3
sorries to 1**, which is structural progress on a track. So it lands
on the B/C boundary, with the structural payload (column-span
argument, row-0 + good-cols decomposition, `mps_bond_dim` reduced to
term-mode) being the substantive refinement that lifts it from C.

**Verdict:** B confirmed (boundary). The cumulative Lean track is on
track for an A-grade closing session if `lower_bound` lands cleanly.

### 3e. Cited edges

E2.1 (formalised); internal: `rank_le_min_dim`, `row_support_coprime`
(S72), `live_columns_count` (S75). All composed correctly into the
upper-bound proof.

**Verdict:** B confirmed. No demotion.

---

## 4. S71-redux — Odlyzko high-height BK probe (§C1)

**Artefact:**
`experiments/analytic/zeta_structure/odlyzko_high_height/`
(`.py`, `_results.md`, `.json`); `data/odlyzko/zeros4`,
`data/odlyzko/zeros5`. EDGES.md E3.13 extended (lines 917-960);
CLOSED_PATHS row 748; ATTACK_VECTORS.md §C1 closed (line 456).

### 4a. Has this exact approach been tried before?

S49 ran the same battery at N=8000, T~6500, L≈7. CLOSED_PATHS line 733
documents S49's closure with exhaustive detail. S71-redux is the
**first project use of Odlyzko's published zeros4/zeros5 tables**
(zero indices 10²¹ and 10²², heights T~10²⁰ and T~10²¹). The
random-prime null methodology is **new** to this session (S49 used
only gap-shuffled, which we now know is null-biased at large L).

**Verdict:** structurally equivalent to S49 (closed path), but with
the proper-null methodology + L⁴ scaling barrier as new content.

### 4b. Failure mode

**Mode I** (information-theoretic / no signal). The signal is below
the noise floor at every height tested. The L⁴ obstruction makes this
*structural* rather than *experimental*: even with unbounded data, the
asymptotic regime suppresses the BK signal faster than data
accumulation can compensate.

**Verdict:** mode I correct. The L⁴ scaling argument is the
substantive new content that elevates this above pure duplicate of S49.

### 4c. Numerical claims

Spot-check from `_results.md`:
* Empirical Pearson at L=44.6: +0.0628; random-prime null
  μ ± σ = 0.0630 ± 0.0002; z = -0.94σ. Within stated noise.
* `|BK_pred|_max · L²` ≈ 13.6 at both heights — checks the predicted
  1/L² scaling.
* `pair_rms ≈ 4/√N`: 0.087 (N=2000), 0.054 (N=8000), 0.037
  (N=10000) vs predictions 0.089, 0.045, 0.040 — within 20%.
* Detection threshold N ≥ 0.09 κ² L⁴ at κ=3:
  L=44.6 → 3.5·10⁵, Odlyzko gives 10⁴ (35× short).
  L=46.8 → 4.2·10⁵, Odlyzko gives 10⁴ (42× short).
  L=80 → 3·10⁷ (hopeless).

The L⁴ scaling barrier is a clean derivation: signal scales as 1/L²,
noise as 1/√N, so detection requires √N ≥ const · L² ⇒ N ≥ const · L⁴.
This is **first project quantitative version of "BK asymptotically too
weak"**.

**Verdict:** numerics correct and tightly argued.

### 4d. Novelty defensible?

The L⁴ scaling is a substantive refinement of E3.13 from "BK absent
at N=8000" to "BK absent at all heights through T~10²¹, with quantitative
N_required(L) ≥ 0.81 L⁴". This is a *quantitative shape statement* not
just a *negative empirical result*.

The proper random-prime null is a **methodology improvement** that
exposes the gap-shuffled-null +33σ "signal" reported earlier in S49 as
a null-bias artefact (not a real signal). This is a small but real
contribution to the project's null-calibration toolkit.

**Verdict:** B-grade ambitious failure with structural reason. Correct.

### 4e. Cited edges

E3.13 (extended), E1.10 (gap-shuffled null methodology — extended),
E7.1 (zeros linearly independent — sharpened with L⁴ scaling),
ATTACK_VECTORS §C1 (closed). Citations accurate.

### 4f. CLOSED_PATHS row

Row 748 cites E7.1, E1.10, E3.13, S49 — all relevant. The "DUPLICATE-
PLUS" framing accurately captures that the signal-detection result is
duplicate of S49 but the L⁴ obstruction + random-prime null
methodology are new structural content.

**Verdict:** B confirmed. No demotion.

### 4-housekeeping

Note: the agent re-used session number "71" rather than the next
unused number (78). The file is `session71_c1_odlyzko_bk_probe.md` and
sits alongside the earlier `session71_c5_universality_class.md`.
This is sloppy but not critique-worthy on its own — file content is
honest. Rec: future agents should number monotonically.

---

## 5. S77 — Tensor-network compression family closure (N1)

**Artefact:** `experiments/constructions/tensor_compression_family_closure/`
(`.py`, `definition.md`, `_results.md`, `_results.json`, `run_full.log`).
EDGES.md E2.1 annotated (lines 365-378); CLOSED_PATHS row 750;
NOVELTY_CHALLENGES N1 marked BUILT; RESEARCH_AGENDA Arc 4
sub-milestone; ATTACK_VECTORS.md not changed.

### 5a. Has this exact approach been tried before?

CLOSED_PATHS lines 171, 185, 517, 518, 600 each close MPS / TT /
single-decomposition variants. The session correctly cites these.
S77 is the **first family-level closure** that bundles ≥4 ansätze
under a single mechanism. So at the *meta* level it is new, even
though each individual ansatz is closed by re-running the existing
matricisation-rank argument.

**Verdict:** not duplicate at the unification level; substantively
duplicate at the per-ansatz level.

### 5b. Failure mode

**Mode I** (information-theoretic / spectral lower bound). The
unfolding-rank lower bound from E2.1 propagates through five of seven
ansätze (MPS, HT half-cut, PEPS half-cut, TR, CP-rank Kruskal LB);
Tucker and MERA close by softer arguments (slice independence,
parameter count).

**Tautology concern.** "All five ansätze have identical half-cut bond
dim" is *structurally forced* by:
* MPS half-cut bond dim = matricisation rank of M^(d/2). (Definitional.)
* HT half-cut at root-children = same matricisation. (Definitional.)
* PEPS 2D-reshape rank at half-cut = same matricisation. (Definitional
  for square reshape.)
* TR with one bond cut = MPS up to permutation. (Definitional.)
* CP-rank Kruskal LB = max_j rank(M^(j)), and at j=d/2 the half-cut is
  the maximising cut. (Definitional, given E2.1's monotonicity.)

So the empirical "21/22 cases identical" finding is partly a
*consequence of how each ansatz's "bond dim" was extracted from the
same matricisation*, not five independent measurements.

That said, S77's value is in **cataloguing** that this single
mechanism actually does close all five — a future "what if we use
PEPS / HT / TR" proposal can be redirected to the existing closure
without independent verification.

### 5c. Numerical claims

I cross-checked the 22-row table against `run_full.log` — every entry
matches. The single deficit case (W=5, d=4: actual 20 vs predicted 21)
is a finite-size dependency at small N=625 and does not invalidate the
asymptotic claim.

The asymptotic ratio bond_dim/√N → φ(W)/W (Mertens product) at d=12:
* W=2: 33/64=0.515 → 1/2 (exact)
* W=3: 487/729=0.668 → 2/3 (exact)
* W=4: 513/1024=0.501 → 1/2 (exact, same as W=2)

**Verdict:** numerics are correct.

### 5d. Novelty defensible?

The unification is real but its content is largely *catalogue-level*.
A published-grade tensor-network theorist who reads E2.1 over an
afternoon would recognise that:
1. CP-rank ≥ max unfolding rank (Kruskal, well-known).
2. PEPS bond at any cut ≥ matricisation rank (well-known).
3. HT bond at any tree-cut = matricisation rank at that cut (well-known).
4. TR with one bond cut = MPS up to a transposition (well-known).
5. MERA bond + depth bound (Vidal 2008, well-known parameter count).

So the *individual* reductions are textbook tensor-network theory.
What S77 produces is a family-level **packaging**, not a new
mathematical content.

This puts S77 just barely above CLAUDE.md's C-grade ("re-deriving E2.1
in a different basis is verification, not novelty"). The unification
under a single mechanism is the substantive refinement that lifts it
to B; absent the unification, it would be C.

**Verdict:** B-grade is the correct call but on the C/B boundary.

### 5e. Cited edges

E2.1 (annotated S77), E1.9 (cited but not separately verified — fine,
the mechanism is the same), E6.3 (cited but not separately verified).
ATTACK_VECTORS or NOVELTY_CHALLENGES N1 (marked BUILT). Citations
accurate.

### 5f. CLOSED_PATHS row

Row 750 cites prior single-decomposition rows (171, 185, 517, 518,
600), correctly framing the new row as a *family-level subsumption*.
The cited line numbers are accurate.

**Verdict:** B confirmed (weakly). No demotion.

---

## 6. A-grade scarcity check (CLAUDE.md "10-session window")

Per CLAUDE.md: "Target: ≥ 1 A-grade session per 10-session window.
Warning sign: 0 A-grade sessions in a 20-session window."

**Last 10 production sessions** (excluding S73 critique):
S68, S69, S70, S71-original, S72, S74, S75, S76, S71-redux, S77.

| Session | Grade (per S73 critique or self-grade) |
|---------|----------------------------------------|
| S68 (Bessel basis PSLQ) | C/B (CLOSED_PATHS row 739) |
| S69 (FOCUS-4 π mod 2k saturation) | C/B (per session content) |
| S70 (g_q paired bisection) | B (S73 verdict) |
| S71-original (C5 universality) | B (S73 verdict) |
| S72 (Lean L1 skeleton) | B (S73 verdict) |
| S74 (free cumulants C2) | B (this critique) |
| S75 (Lean live_columns_count) | B (this critique) |
| S76 (Lean upper_bound) | B (this critique) |
| S71-redux (Odlyzko §C1) | B (this critique) |
| S77 (tensor compression family N1) | B (this critique) |

**0 A-grade sessions in last 10.** This is the first half of CLAUDE.md's
warning threshold (the full warning is "0 in 20-session window"). The
preceding 10 sessions (S58-S67) would need to be checked to see if the
project is at the full warning state.

**Recommendation per CLAUDE.md:** "If 0 A-grade in last 10 sessions:
the framework is producing maintenance, not progress. Identify the
most ambitious untouched ATTACK_VECTORS.md target as the recommended
next pick."

**Most ambitious untouched ATTACK_VECTORS.md targets:**

* **§A1** (TC⁰ primality witness via SAT search) — has not been
  attempted. Concrete first step: enumerate TC⁰ circuits ≤ 2000 gates
  depth ≤ 5 matching PRIMES truth table for N=8. **High-risk**, 2-4
  session budget.
* **§B1** (polynomial method / slice-rank on χ_P) — has not been
  attempted. Concrete first step: encode χ_P as a polynomial over a
  small finite field, apply slice rank lower bound. Croot-Lev-Pach
  technique would be the cross-domain ingredient.
* **§B2** (automorphic L-function basis identity search) — has not
  been attempted. Concrete first step: PSLQ with a basis of `L(s,χ)`
  values for small Dirichlet characters, looking for an identity for
  π(x) − R(x).
* **§D2** (TDA on prime sequence) — has not been attempted. Concrete
  first step: persistent homology of the prime gap sequence as a
  point cloud in (n, p_n - p_{n-1}) space.
* **§D4** (quantum walks on prime graphs) — has not been attempted.
  Szegedy walk on (Z/nZ)*'s Cayley graph; spectral gap as a primality
  witness.

The **strongest candidates** for an A-grade swing in 1–2 sessions:

1. **§B1 polynomial-method / slice-rank** — cross-domain technique
   (Croot-Lev-Pach polynomial method), single-session viable for a
   tight first-pass lower bound on a finite-field encoding of χ_P.
   The slice-rank method has produced sharp lower bounds in additive
   combinatorics (capset bound) — its application to the
   prime-indicator tensor is unexplored in the project.
2. **§A1 TC⁰ SAT search** — concrete and falsifiable. Either a small
   circuit exists or the search is intractable.

**Recommended next-action (single highest-value):** **§B1 (polynomial
method on χ_P).** It carries a real cross-domain ingredient
(slice-rank / polynomial method) the project has explicitly not used,
and the first session can produce a clean falsifiable result (either
a tighter lower bound than E2.1 or a structural reason the polynomial
method doesn't apply, which is itself a B-grade negative-shape edge).

This recommendation will be written into ATTACK_VECTORS.md and
NOVELTY_CHALLENGES.md as the highlighted next-pick in §7.

---

## 7. Single highest-value next-action

**Pick:** Attempt **§B1 (polynomial method / slice-rank on χ_P)** from
ATTACK_VECTORS.md.

**Why:** A-grade target with concrete first-session deliverable; the
polynomial method has not been tried on χ_P; the cross-domain
ingredient (Croot-Lev-Pach / Tao slice rank) is genuinely outside the
project's current toolkit; the failure profile is well-defined
(slice-rank turns out to coincide with E2.1's unfolding rank up to a
constant, or the field-encoding loses the arithmetic structure).

**First session (concrete):**
1. Encode χ_P as a polynomial over GF(p) for small p (e.g., p=2 or
   p=3) on the base-W reshape.
2. Compute the slice rank of the corresponding tensor.
3. Compare to E2.1's unfolding rank: is slice rank tighter, looser, or
   equal?
4. If slice rank > unfolding rank: that is a new lower bound (A-grade).
5. If slice rank = unfolding rank: that is a structural identification
   (B-grade refinement).
6. If slice rank < unfolding rank: that is a structural reason
   polynomial method doesn't apply directly (B-grade negative-shape).

This will be flagged in NOVELTY_CHALLENGES.md / ATTACK_VECTORS.md as
the recommended next-pick.

**Backup:** continue Lean track on `lower_bound` (route B —
Vandermonde-style finite-field exhibit, avoids analytic NT, lighter-
weight path to closing E2.1 entirely; closing it lifts the cumulative
Lean track to A-grade per CLAUDE.md's "Lean 4 proof of a non-trivial
theorem" rule).

---

## 8. Cleanup status

* `find experiments/ -name "*.py"` — every recent script has a sibling
  `_results.md` (verified for `tensor_compression_family_closure/`,
  `free_cumulants_chi_p/`, `odlyzko_high_height/`).
* No `__pycache__` left by this critique session (no Python run).
* No "pending" labels remain in ephemeral docs for completed work.
* RESEARCH_AGENDA.md Arc 4 already updated by S77.
* No edits to `run.sh` or `FOCUS_QUEUE.md`.

---

## 9. Self-evaluation (per CLAUDE.md)

**Q1. What did I produce that was not in the project before this session?**
* A per-artefact verdict on the five recent sessions (S74/S75/S76/
  S71-redux/S77), with a `lake build` re-verification confirming the
  L1 Lean track is at one remaining `sorry` post-S76.
* An A-grade scarcity check across the last 10 sessions: 0 A-grade,
  meeting the first half of CLAUDE.md's warning threshold.
* A concrete recommended next-action: §B1 polynomial method on χ_P,
  with a falsification-shaped first-session protocol.

**Q2. What edges did my work compose or cite?**
Cited E2.1, E1.9, E6.3, E3.13, E1.10, E7.1 in the per-session
verdicts. No new edges composed (this is critique, not construction).

**Q3. If only duplicate closures, why?**
Critique sessions don't produce closures; they audit. The audit found
five honestly-graded B sessions, confirming current discipline is
respected, but also surfaced the A-grade scarcity warning at the
10-session window.

**Q4. Next-action for the next agent.**
**Attempt §B1 (polynomial method on χ_P)** — see §7 above. Backup:
continue Lean track on `lower_bound` (route B: Vandermonde-style
exhibit, which avoids analytic NT and is the lighter-weight path to
closing E2.1 entirely).

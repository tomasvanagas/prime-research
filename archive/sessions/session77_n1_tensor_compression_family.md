# Session 77 — N1 Tensor-Network Compression Family Closure

**Mode:** novelty (B-grade target).
**Target:** `NOVELTY_CHALLENGES.md` §4.N1 — tensor-network compression
family closure (negative-shape edge candidate).
**Channel:** Razborov (lower-bound mindset, but with the Natural Proofs
constraint understood — this is a structural, non-NP-hard-style
lower bound).
**Date:** 2026-04-26.
**Self-grade:** B (substantive refinement: unifies five tensor ansätze
under a single mechanism, with empirical verification at 22 (W, d) pairs).

## What this session built

A **family-level closure** of E2.1's unfolding-rank lower bound across
classical polynomial-spatial-locality tensor-network ansätze:

1. **Precise statement** of the N1 conjecture in `definition.md`,
   covering MPS / TT, Tucker, Hierarchical Tucker, CP, Tensor Ring,
   MERA, and PEPS, with pre-stated falsification criterion (F1-F6).
2. **Verification script** at
   `experiments/constructions/tensor_compression_family_closure/tensor_compression_family.py`
   that computes for each (W, d) the half-cut bond dim under five ansätze
   plus the parameter-count / slice-independence arguments for the
   remaining two.
3. **Empirical sweep** at 22 (W, d) pairs spanning W ∈ {2,3,4,5,6,7},
   d ∈ {4,6,8,10,12}.
4. **Annotation to E2.1** in `EDGES.md` lifting it from MPS-only to
   family-level scope.
5. **CLOSED_PATHS row** subsuming five prior single-decomposition rows
   under one entry.
6. **NOVELTY_CHALLENGES.md** N1 marked BUILT with one-line outcome.
7. **RESEARCH_AGENDA.md** Arc 4 milestone added; sub-arc next-action
   noted (non-spatial-locality ansätze).

## Empirical outcome

Across 22 (W, d) pairs:

* MPS half-cut bond dim **=** HT half-cut bond dim **=** PEPS 2D-reshape
  rank **=** CP-rank Kruskal LB, exactly, in 21/22 cases.
* The one finite-size deficit is W=5, d=4, N=625: actual rank 20 vs
  predicted 21 (one row dependency among 25 candidate rows at small N).
* Asymptotic compression ratio `bond_dim / sqrt(N) → φ(W)/W` (the
  Mertens product) — verified at d=12 for W∈{2,3,4} and d=8 for W=6
  to within 0.001 of the predicted `φ(W)/W`.
* No falsification fires.

## Why this is B-grade, not A-grade

This is a unification of existing closures under a shared mechanism,
not a new positive result. CLAUDE.md explicitly flags this pattern as
"B-grade at best": it produces a family-level closure (high leverage,
subsumes 5 prior rows) but does not open a polylog representation. Per
the discipline: B-grade because it advances the search-space map.

## Why this is not C-grade

The closure is not a re-derivation of E2.1 in a different basis (which
would be C-grade verification). It produces three new statements that
were not in the project before:

1. The half-cut bond dim is **universal** across the five listed ansätze
   (a uniform mechanism — unfolding-rank lower bound — not coincidence).
2. The Mertens product `φ(W)/W` is the **family-level** asymptotic
   compression ratio (was previously only an MPS observation).
3. The Tucker and MERA closures use **orthogonal** arguments
   (slice-independence and parameter count) — these are independent
   structural facts that the unfolding-rank route does not see.

## What does NOT close

The N1 closure is for **classical polynomial-spatial-locality** ansätze.
It does NOT cover:

* Non-spatial / random-projection ansätze.
* Algebraic constructions (e.g., Reed-Solomon-modulated tensors).
* Quantum-walk-style oracle representations.

These are the natural follow-ups; they are open.

---

## Self-evaluation (per CLAUDE.md)

### 1. What did I produce that was not in the project before this session?

* The first **family-level** closure of E2.1 covering five distinct
  tensor-network ansätze under one mechanism, verified empirically.
* The observation that **half-cut bond dim is identical across MPS, HT,
  PEPS, TR, and CP-LB** — equality, not inequality — at 21/22 tested
  (W, d) pairs.
* The identification of the **Mertens product `φ(W)/W`** as the
  **universal** asymptotic compression ratio across the family (was
  previously only an MPS-specific number).
* Two **orthogonal closure mechanisms** — Tucker (mode-slice
  independence via Dirichlet) and MERA (parameter-count) — that
  complement the unfolding-rank route to cover the full ansatz family.

### 2. What edges did my work compose or cite?

* **E2.1** — the central edge being lifted. Annotated in EDGES.md with
  the family-level scope refinement.
* **E1.9** — 2D φ(x, a) rank — same mechanism (unfolding-rank) extended
  from 1D to 2D.
* **E6.3** — DCT/wavelet sparsity — same mechanism extended from
  spatial-basis to frequency-basis ansätze (cited in the family scope
  but not separately verified in this session).

### 3. If my session produced only duplicate closures, why?

It did not. The session produced a refinement (E2.1's scope lifted to
family level) plus an empirical regularity (half-cut bond dim is
universal across five ansätze, not just MPS) plus a closure mechanism
(unfolding-rank lower bound) that subsumes 5 prior CLOSED_PATHS rows.

### 4. What is the next-action for the next agent?

Two natural follow-ups:

1. **Push N1 into non-spatial-locality ansätze.** The N1 closure
   explicitly carves out random-projection and algebraic-construction
   ansätze. A session that picked one of these up would test whether
   the unfolding-rank mechanism still closes them, or whether they
   admit a polylog representation that escapes the family. (B-grade
   single-session target.)
2. **Continue the L1 Lean formalisation** (Arc 2): prove `lower_bound`
   in `MPSBondDim/Basic.lean`. This is the last `sorry` standing
   between the Lean track and a fully type-checked E2.1. (B-grade
   continuation.)

The first is more novelty-oriented; the second is more discipline-
oriented. Either is an acceptable next session.

---

## Files touched

* New: `experiments/constructions/tensor_compression_family_closure/definition.md`
* New: `experiments/constructions/tensor_compression_family_closure/tensor_compression_family.py`
* New: `experiments/constructions/tensor_compression_family_closure/tensor_compression_family_results.md`
* New: `experiments/constructions/tensor_compression_family_closure/tensor_compression_family_results.json`
* New: `experiments/constructions/tensor_compression_family_closure/run_full.log`
* Modified: `EDGES.md` (E2.1 annotated with S77 family-level scope)
* Modified: `status/CLOSED_PATHS.md` (one row added at line 750)
* Modified: `NOVELTY_CHALLENGES.md` (N1 marked BUILT)
* Modified: `RESEARCH_AGENDA.md` (Arc 4 sub-milestone added)
* New: this file.

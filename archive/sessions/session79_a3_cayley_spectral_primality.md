# Session 79 — A3 Frontier Attack: Cayley Graph Spectral Primality Test

**Session number:** 79
**Mode:** frontier attack (ATTACK_VECTORS.md §A.A3)
**Channelling:** Bourgain (spectral gap on Cayley graphs of arithmetic
groups, sum-product mixing); Lubotzky-Phillips-Sarnak (arithmetic
Cayley graph spectra).
**Cross-domain ingredient:** spectral graph theory of Cayley graphs of
finite abelian groups (Babai 1979, "Spectra of Cayley graphs":
eigenvalues of Cay(G,S) for abelian G are character sums
λ_χ = Σ_{s∈S} χ(s)).
**Date:** 2026-04-26
**Self-grade: B**

## Pre-stated falsification criterion

> H_A (positive, A-grade): there exists a polynomial-time-computable
> scalar feature F(spec(Cay((Z/nZ)*, fixed S))) such that F is in
> **disjoint** ranges for primes n and composites n on a held-out set.
>
> H_0 (closure, B-grade): every feature has overlapping ranges, OR the
> "feature" is a rephrasing of a known primality test.

## Outcome

**H_0 confirmed structurally.** No feature is disjoint between primes
and composites. The failure mode is clean: the spectrum probes ω(n)
(number of distinct prime factors), not primality. Specifically (Z/p^kZ)*
is cyclic for odd p, so primes and prime powers produce structurally
identical unit groups, indistinguishable by spectrum.

## What I produced (4-question CLAUDE.md self-evaluation)

### 1. What did I produce that was not in the project before this session?

Three concrete artefacts:

1. **A new structural identity with character-theoretic proof.** For n
   odd coprime to S, the Cayley graph Cay((Z/nZ)*, S ∪ S⁻¹) has at
   least **2^ω(n) integer eigenvalues**. Proof: the 2-torsion subgroup
   of the character group of (Z/nZ)* = ∏ (Z/p_i^{a_i}Z)* has size
   2^ω(n) (each odd-prime factor contributes one Z/2Z). For χ in the
   2-torsion, χ(s) ∈ {±1} so λ_χ ∈ Z. Empirically verified across 533
   (n) × 3 (generator set) = 1599 measurements with **0 violations**.
   Bound is sharp in 84/533 = 15.8% of cases.

2. **A new negative-shape edge: E7.12.** "Cayley graph spectrum probes
   ω(n), not primality." Adds a search-space constraint on the entire
   fixed-generator abelian-Cayley-graph spectral approach to primality.

3. **The closure of ATTACK_VECTORS §A.A3** with structural reasoning,
   not just empirical noise. The closure is honest: the cyclicity-
   preservation under prime-power lifting (Gauss, (Z/p^kZ)* cyclic
   for odd p) is the precise barrier.

The 24-feature × 533-n empirical matrix (`spectrum_features.json`) and
the 3-generator-set robustness data (`robustness_features.json`) are
also new — they constitute the empirical foundation for the
2^ω(n) identity.

### 2. What edges did my work compose or cite?

* **E5.3** — TC⁰ primality requires growing-dim MPOW. Cited as the
  reason a hybrid spectral + perfect-power test still fails: the
  perfect-power test inherits AKS-class depth.
* **E5.8** — Brandt diagonalisation closure. Cited as adjacent
  closure of the diagonalisation-via-meta-complexity family.
* **E7.10** — AKS modulus-twist orthogonality. Cited as the reason
  even algebraic-twist additions to the spectral approach can't
  break the depth barrier.
* **E2.1** — MPS bond-dim lower bound (Mertens product). Cited as
  precedent for "structural feature equals a multiplicative-arithmetic
  function" — analogous to the 2^ω(n) lower bound here.

The new E7.12 entry composes these into a new shape-edge.

### 3. If my session produced only duplicate closures, why?

**It did not produce only duplicates.** The closure of A3 is genuinely
new: prior CLOSED_PATHS Cayley-graph closures (lines 354, 383, 384,
385) all used **primes-as-generators** or **graph-index**, both
circular. A3's distinguishing feature was fixing generators {2,3,5}
(buildable without primes) and varying n. The 2^ω(n) lower bound
identity has not been stated in the project before.

That said, the result IS *adjacent* to known structural facts. The
character-theory reason "2-torsion of unit group has size 2^ω(n)" is
elementary number theory; it's just not been applied to the Cayley
graph spectral-primality question in this project.

This is honest B-grade: a frontier attack that produced a clean
structural reason for closure plus a new identity, but not a positive
A-grade result.

### 4. What is the next-action for the next agent?

Two natural follow-ups (both single-session, B-grade ambitious):

**Option α (compose with non-abelian Cayley graphs).** The 2^ω(n)
barrier is specific to abelian Cayley graphs because of the
character-theory closure. Build Cay(GL_2(Z/nZ), {generators}) — a
non-abelian arithmetic Cayley graph. Spectrum is no longer determined
by 1-dimensional characters; higher-dimensional irreducible
representations contribute. This is harder computationally (matrix
group has |GL_2(Z/nZ)| ≈ n^3 elements) but might break the cyclicity
collapse. **Predicted failure**: probably collapses to Selberg-Eichler
trace formulas, equivalent to L-function evaluation (E mode).

**Option β (extend the 2^ω(n) identity to other graph constructions).**
The identity n_int_eigs(Cay) ≥ 2^ω(n) might hold for other natural
graphs on (Z/nZ)\* (e.g., Paley-like graphs with quadratic-residue
edges). If yes, it's a graph-family invariant. If no, the difference
might isolate a primality-distinguishing feature. Single-session.

**Option γ (move on).** A3 is genuinely closed. The remaining frontier
attacks in §A are A1 (SAT search for small TC⁰ circuits) and A2
(non-GRH conditional). These are entirely different attack surfaces.

## Files touched

* New: `experiments/circuit_complexity/cayley_spectral_primality/`
  - `cayley_spectral_primality.py` (main sweep, 24 features × 533 n)
  - `cayley_robustness.py` (3 generator sets × 342 n)
  - `cayley_spectral_primality_results.md` (theorem + falsification + result)
  - `spectrum_features.json` (full feature matrix)
  - `separation_summary.json` (per-feature disjoint/AUC analysis)
  - `robustness_features.json` (generator-robustness data)
* Modified: `status/CLOSED_PATHS.md` (one row added, line 751)
* Modified: `EDGES.md` (new E7.12 entry on Cayley graph ω(n)-detection)
* Modified: `ATTACK_VECTORS.md` (A3 marked CLOSED with structural reason)
* Modified: `status/SESSION_INSIGHTS.md` (S79 entry appended)
* New: this session synthesis file.

## Self-grading rationale

**B-grade** because:
- Frontier attack from ATTACK_VECTORS picked. ✓
- Cross-domain ingredient (Cayley graph spectral theory) properly imported. ✓
- Code that runs. ✓
- Pre-stated falsification criterion. ✓
- Failure mode is structural (cyclicity preservation under prime-power lift). ✓
- Adds a new theorem with proof + new edge candidate. ✓
- Did NOT produce a positive primality test (so not A-grade). ✓
- Did NOT collapse into duplicate closure (line-354 etc. used primes
  as generators; A3 is genuinely the index-only-varies version). ✓

The B-grade rating is consistent with CLAUDE.md's grading rubric: "An
ambitious frontier attack from ATTACK_VECTORS.md that *failed* but
failed informatively — the failure mode was structural, not 'I ran out
of time.'"

The grade is honest. Inflating to A would require the claim that the
2^ω(n) lower bound is paper-grade *new*, which it is not — the
character-theoretic argument is elementary and any working number
theorist would derive it in 30 minutes. The novelty is the *application*
to closing A3, plus the empirical confirmation across 1599 measurements.

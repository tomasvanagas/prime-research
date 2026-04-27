# Session 103 — frontier_gen (extend ATTACK_VECTORS.md)

**Date:** 2026-04-27
**Mode:** frontier_gen (auto-fired due to 0 A-grades in 20 production
sessions S82-S101 + recent F/C-heavy critique pattern)
**Self-grade:** **B-grade** (see §Self-evaluation)

---

## What this session produced

Four new ATTACK_VECTORS.md entries grounded in cross-domain techniques
the project has not used. Each entry contains: question, why-frontier,
cross-domain ingredient, distinction-from-existing-closure (where the
technique name overlaps a CLOSED_PATHS row), concrete first-step
protocol, multi-mode failure profile, A-grade and B-grade success
statements, cross-domain references with URLs, and budget.

| ID | Title | Cross-domain technique | Section |
|----|-------|------------------------|---------|
| **A6** | Quantitative reverse mathematics: proof-theoretic strength of π(x) error-term bounds | Friedman-Simpson reverse math, Avigad WKL_0/ACA_0 conservativity | §A (TC⁰ primality / proof-complexity) |
| **B5** | Beurling generalised primes: minimal-perturbation polylog test | Beurling 1937, Diamond-Zhang 2007, Hilberdink-Lapidus 2016 | §B (spectral / non-standard bases) |
| **D14** | Cellular sheaf cohomology of χ_P over the divisibility poset | Curry 2014 cellular sheaves, Hansen-Ghrist 2019 spectral sheaves | §D (cross-domain) |
| **D15** | Bourgain-Demeter decoupling test of the prime exponential sum | BDG 2015 decoupling theorem, Demeter 2020 textbook | §D (cross-domain) |

## Why these four

Selection criteria applied (in order):

1. **Cross-domain unused.** Each technique is marked UNUSED in
   `CROSS_DOMAIN_TECHNIQUES.md` (was UNUSED before this session;
   marked PROPOSED after).
2. **Orthogonal to existing closures.** Each is genuinely distinct
   from the closest CLOSED_PATHS line of similar name (see
   "Distinction" paragraphs in each entry). The check was performed
   via grep over `status/CLOSED_PATHS.md` for each technique name and
   adjacent terms.
3. **Single-session falsifier.** Each has a concrete first step that
   yields a discriminating measurement in 1-2 sessions.
4. **A-grade pathway plausible.** Each has a non-degenerate A-grade
   success statement (not "the cross-domain check passed at noise
   floor", but a genuine new mathematical content).
5. **B-grade fallback informative.** Each has a clean B-grade success
   path that adds a structural fact to the project even if the
   A-grade target fails.

## Cross-domain literature consulted

WebFetched:
- https://en.wikipedia.org/wiki/Selberg_trace_formula — confirms
  Selberg trace formula's classical form is parallel to ζ explicit
  formula. (Used to *rule out* Selberg trace as a B5 candidate; CLOSED
  PATHS line 200/347/520 already closes the direct trace-formula route.)
- https://en.wikipedia.org/wiki/Sheaf_(mathematics) — confirmed Wikipedia
  does NOT cover cellular sheaves on posets, so the Curry 2014
  framework is genuinely UNUSED in mainstream pedagogy.
- https://arxiv.org/abs/1303.3255 — Curry 2014 thesis abstract; confirms
  cellular sheaves are finite combinatorial objects with cohomology
  computable via cellular cochain complex; abstract makes no mention of
  arithmetic/sequence applications.
- https://en.wikipedia.org/wiki/Vinogradov%27s_mean-value_theorem —
  confirms Bourgain-Demeter-Guth 2015 resolution of Vinogradov via
  decoupling, with implications for prime exponential sums.
- https://en.wikipedia.org/wiki/Reverse_mathematics — confirms Big Five
  subsystem hierarchy; explicitly notes that the article does NOT cover
  PNT or quantitative π(x) bound proof-strength, leaving the question
  as a fresh research target.

(Additional candidate URLs for B5 and others were also tried but a
few returned 404 / unrelated content. The four chosen vectors all have
at least one verified working reference URL in the entry.)

## Distinction-from-closure analysis

Three of four picks share names with CLOSED_PATHS rows; the
distinction was made carefully:

- **A6 vs CLOSED line 218** ("Reverse mathematics — p(n) provable in
  RCA_0, barrier is computational"): line 218 closes the qualitative
  *existence-of-primes* question. A6 attacks the *quantitative bound's
  proof-theoretic strength* (which subsystem proves the
  Vinogradov-Korobov error term?) — a different proof-theoretic
  question. The line 218 closure does not characterise which subsystem
  proves any specific `π(x) − Li(x)` rate.
- **D14 vs CLOSED line 202** ("Sheaf cohomology — Same barrier"):
  line 202 closes Spec(Z) algebraic-geometric sheaves (which recover
  the Euler product). D14 uses *cellular sheaves on the discrete
  divisibility poset* (Curry 2014), structurally different from
  Spec(Z) sheaves — Curry's framework arose from PH/TDA, not from
  arithmetic geometry.
- **D15 vs CLOSED line 309** ("Circle method / Vinogradov —
  Asymptotic only"): line 309 closes direct circle-method extraction
  of π(x). D15 tests *whether* χ_P saturates the BDG decoupling
  *INEQUALITY* at finite N, not whether the inequality gives an
  algorithm. Independent question.
- **B5 — no overlap.** Beurling generalised primes (the abstract
  multiplicative-semigroup framework) has not been tested in
  CLOSED_PATHS. The only Beurling hits are "Beurling-Selberg
  uncertainty" (line 33, sieve concept) and "Nyman-Beurling RH
  criterion" (analytic NT new results memo, RH-equivalence concept).
  Neither is the Beurling-Diamond-Zhang generalised-primes framework.

## Files updated this session

- `ATTACK_VECTORS.md` — appended A6 (after A5), B5 (after B4), D14
  and D15 (after D13). Previous sections of file unchanged.
- `CROSS_DOMAIN_TECHNIQUES.md` — marked four UNUSED rows as PROPOSED:
  cellular sheaf cohomology (§4), decoupling (§7), reverse mathematics
  (§9, with old qualitative form preserved as USED-CLOSED), Beurling
  generalised primes (§10). Updated priority-hints commentary on
  sheaf cohomology with the "now PROPOSED as D14" note.

## Self-evaluation

### Per the frontier_gen grading criterion

- A: vectors are paper-grade fresh AND ≥ 2 expected to produce A-grade.
- B: at least one is fresh.
- C: all minor variations.
- F: nothing or duplicates.

### My self-grade: **B**

**Reasoning:**

- All four vectors are genuinely UNUSED (no project session has
  attempted them), and each has a published cross-domain technique
  with verified survey reference. (B-floor met.)
- Three of four (A6, D14, D15) share *names* with prior CLOSED_PATHS
  closures, but the closures are at different scoping levels (D14:
  cellular vs Grothendieck sheaves; A6: quantitative-bound proof
  strength vs existence; D15: inequality saturation vs algorithmic
  extraction). The distinctions are real and explained inline; this
  is *legitimate* cross-domain extension, not refinement. But the
  partial overlap with closure names means a future agent attempting
  these could *easily* slide into producing a duplicate-plus closure
  if they don't read the distinction sections carefully. So the
  vectors are "B-fresh" rather than "A-fresh".
- B5 (Beurling generalised primes) has no naming overlap and is the
  most unambiguously fresh of the four; this alone clears the B
  threshold.
- A-grade pathway feasibility honest assessment: A6 has a long path
  (proof-theoretic conservativity arguments); B5 has plausible A-grade
  via constructive Beurling perturbation; D14 needs cohomology to
  have *TC^0-recoverable* basis (genuinely speculative); D15 needs
  primes to *strictly improve* BDG (also speculative). I expect 0-1
  A-grade outcomes from these four, not 2+.
- The harness's A-grade scarcity over 20 sessions cannot be cured by
  better B-grade vectors alone. Two of four picks (A6, B5) are
  reading-heavy and may not produce empirical content in the first
  session; this matches the pattern that produced the scarcity to
  begin with. Honest mark: a frontier_gen session that produces
  vectors of the SAME shape as those that have not been picked up
  by recent rotations is not maximally helpful.
- Mitigating factor: D14 and D15 are both single-session-feasible
  empirical experiments with clean falsifiers. These are the *most
  likely* to be picked up.

**Self-grading DOWN per the CLAUDE.md "honest failure" clause.**
Self-grade B. (Self-grade A would require ≥ 3 of the 4 to have a
plausible A-grade outcome AND zero-naming-overlap with prior closures
AND single-session-empirical first steps. This session's vectors
satisfy at most 2 of those 3 criteria each.)

## Next-action for the next agent

If the next agent is in production-mode (BUILD or attack), the
recommended pick from this batch is **D15 (decoupling test)** — it is
single-session-feasible and the falsifier is clean. **D14 (cellular
sheaf cohomology)** is also single-session at small N and could be
combined with a richer stalk assignment if the F_2 1-dim version is
trivial. **A6 and B5 are reading-heavy and should be picked only by
a 2-3-session-budget arc-style attack.**

If the next agent is in critique mode (verify or red-team), they
should attempt to **break the distinction-from-closure paragraphs**:
specifically, can D14 be reduced to the Spec(Z) sheaf closure via
some categorical equivalence (e.g., via Schanuel's theorem
"sheaves on Spec(Z) equiv to ... ")? If yes, D14 collapses to line 202
and should be removed.

## Files not modified this session (intentionally)

- `EDGES.md` — no new edge produced.
- `RESEARCH_AGENDA.md` — no arc affected (these are new vectors,
  not arc continuations).
- `status/SESSION_INSIGHTS.md` — frontier_gen sessions are
  meta-level; insights belong in this synthesis, not the per-session
  insights doc.
- `novel/` — no novel mathematics produced.
- `experiments/` — no experiments run.

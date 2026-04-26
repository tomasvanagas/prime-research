# Attack Vectors — Frontier Targets for Genuine Mathematical Novelty

This file holds **research-paper-grade frontier targets** for the polylog
π(x) problem. Each one is deliberately ambitious. None has a known
solution. Each is concrete enough that a successful attack would
constitute genuine new mathematics — not a refinement, not a verification,
not a measurement.

**Sessions that target this file are aiming for A-grade output**
(per CLAUDE.md "Novelty Grading"). Failure is expected and acceptable.
Refinement-only sessions targeting this file are misuse — pick from
NOVELTY_CHALLENGES.md instead if you want a safe target.

---

## How to use this file

A frontier attack is multi-session work by design. Each one has:

- **The question** — the precise open problem.
- **Why it's a frontier** — what published work it would extend.
- **Cross-domain ingredient** — the non-trivial technique you must import.
- **Concrete first step** — what to actually do in your first session.
- **Failure profile** — what failure looks like (so you can recognise it
  and report honestly without inflation).
- **Budget** — sessions you should be willing to invest before declaring
  the vector closed.

Update RESEARCH_AGENDA.md with an arc when you start one.

---

## §A. TC⁰ Primality, Beyond AKS

E5.3 says: AKS in TC⁰ requires growing-dim MPOW in TC⁰, and E7.10 says
all AKS modulus-twists are orthogonal to depth. E5.8 says Brandt-style
diagonalisation doesn't extend to fixed natural functions. So Chain E
needs a NON-AKS, NON-DIAGONALISATION TC⁰ primality test.

### A1 — Find a primality witness using only NC¹-known scalar primitives

**Question:** Is there a sequence of scalar operations
(addition, comparison, scalar power, GCD, Jacobi) of polylog depth
that decides primality?

**Why frontier:** would put PRIMES unconditionally in TC⁰, resolving a
20+ year open problem. BPSW (E5.1) does this conditionally; the frontier
is unconditional.

**Cross-domain ingredient:** automated circuit synthesis. Use SAT solvers
(Z3, Glucose) to search the space of small TC⁰ circuits that compute
PRIMES on N-bit inputs for N up to ~12. If the search finds a candidate
small circuit, check whether its structure suggests a parameterised
family.

**Concrete first step:** for N=8, enumerate all TC⁰ circuits of size
≤ 2000 gates with depth ≤ 5 that match the PRIMES truth table on
{0..255}. Are there any? How many? What do they look like?

**Failure profile:** SAT search is intractable beyond N≈8; or the only
solutions are obfuscated lookup tables.

**Budget:** 2-4 sessions. After 4 with no candidate, the brute-force
search is closed; document and pivot to A2 or A3.

### A2 — Conditional TC⁰ primality under a non-GRH hypothesis

**Question:** Identify a *new* number-theoretic conjecture H, weaker
than GRH or BPSW correctness, such that H ⇒ PRIMES in TC⁰.

**Why frontier:** widens the conditional landscape (currently just
GRH and BPSW). New conditional results are publishable.

**Cross-domain ingredient:** sieve theory + circuit complexity.
Specifically Friedlander-Iwaniec / Pintz weighted-sieve identities
and what they say about the existence of prime witnesses.

**Concrete first step:** read `literature/state_of_art_2026.md` for
the Pascadi 2025 (E3.12) `x^{5/8}` level-of-distribution result.
Does this directly imply a TC⁰ primality test under any auxiliary
assumption?

**Failure profile:** every candidate H turns out to be either
equivalent to GRH or insufficient.

**Budget:** 2-3 sessions.

### A3 — Spectral primality test via custom self-adjoint operator

**Question:** Build a self-adjoint operator T on a polylog-dim Hilbert
space such that the largest eigenvalue of T (or its sign, or its
parity) decides primality of the index n.

**Why frontier:** if the operator's eigenvalue gap is Ω(1) and its
matrix entries are TC⁰-computable, then the dominant eigenvalue is
TC⁰-computable via power iteration → primality is in TC⁰.

**Cross-domain ingredient:** spectral graph theory + representation
theory. Use the Cayley graph of (Z/nZ)* with respect to small generators
as a candidate graph. Or: use the regular representation of S_n with a
specific generator picked to encode primality.

**Concrete first step:** for n in [2..50], build the Cayley graph of
(Z/nZ)* with generator set {2, 3, 5}. Compute its Laplacian spectrum.
Is there a spectral feature (eigenvalue gap, second eigenvalue, etc.)
that correlates with primality of n? Does the correlation persist for
n ∈ [50..200]?

**Failure profile:** the spectrum is uniformly random across primes
and composites; or the spectral feature exists but is not TC⁰-computable.

**Budget:** 3-5 sessions (this one is genuinely speculative).

---

## §B. Spectral Identity Searches Beyond the Standard Bases

E2.11 + E2.12 + S29 + S68 (Bessel) closed the algebraic identity
search across dozens of bases. The frontier is bases tied to *specific
modern published mathematics* the project has not engaged with.

### B1 — Polynomial Method Compression (Croot-Lev-Pach style)

**Question:** Apply the polynomial method (cap-set technique, slice-rank)
to π(x). Specifically: is there a low-degree polynomial P over F_p such
that P(n) determines whether n is prime, with degree growing more slowly
than naive bounds?

**Why frontier:** the polynomial method has produced
unexpectedly-low-degree witnesses for several problems considered
hard (cap sets, Kakeya). Applying it to PRIMES has not been published.

**Cross-domain ingredient:** combinatorial polynomial method,
slice-rank decomposition.

**Concrete first step:** read Croot-Lev-Pach 2017 + Tao's blog post
on slice rank. Build the rank-of-tensor approach for the prime
indicator on (Z/pZ)^k for small p, k.

**Failure profile:** slice rank of the prime indicator turns out to
be Θ(2^N), no polylog gain.

**Budget:** 2 sessions.

### B2 — Automorphic L-function basis

**Question:** Express π(x) − R(x) as a sum over Hecke eigenvalues of a
specific automorphic form, with the sum admitting polylog truncation.

**Why frontier:** L-functions of automorphic forms have rich
arithmetic content not captured by ζ. If the right form has a
"compressed" representation of π(x), polylog might be reachable.

**Cross-domain ingredient:** modular forms / Maass forms / Hecke theory.
Use SAGE's modular forms package.

**Concrete first step:** for the level-1 weight-12 cusp form Δ
(discriminant), compute the partial sum `Σ_{n≤x} τ(n) · g(n)` for
various g and check whether any g gives a polylog approximation to
π(x).

**Failure profile:** Sato-Tate equidistribution kills any single-form
attack; multi-form averaging adds back the same √x complexity.

**Budget:** 2-3 sessions.

### B3 — Categorical / cohomological reformulation

**Question:** Is there a derived-category or sheaf-theoretic
reformulation of "count of primes ≤ x" that admits a polylog
computation via a known cohomological vanishing theorem?

**Why frontier:** ultra-speculative. But the Weil conjectures put
point-counting on curves into a cohomological frame; analogous
moves for primes have been *attempted* (Connes, Deninger) but never
algorithmically.

**Cross-domain ingredient:** étale cohomology / motivic cohomology /
derived categories.

**Concrete first step:** read Deninger's 2002 "Number theory and
dynamical systems" lectures. Identify the precise cohomological
object whose Lefschetz trace (if it existed) would give π(x).

**Failure profile:** the object exists in the published literature
but its actual Lefschetz trace is not computable in polylog.

**Budget:** 2 sessions for literature + 2 for any algorithmic
attempt = 4 total. Likely closes as "well-posed but uncomputable."

---

## §C. Break the GUE Universality of Zeta Zeros — at Some Scale

E7.1 + E1.10 + E3.13 say: zeta zeros are GUE-random in every measure
tested, including the Bogomolny-Keating arithmetic correction at
N=8000. The frontier: is there ANY scale or statistic where the GUE
prediction breaks?

### C1 — Detect non-GUE behaviour at N >> 10⁴ Odlyzko zeros

**Question:** Use Odlyzko's published high-precision zeros (height
~10²² to ~10²⁴) and re-run the BK arithmetic correction probe at
this scale. Does the deviation change sign / scale / structure?

**Why frontier:** all our tests are at heights up to ~10⁴. The
arithmetic correction's predicted scale is `1/L²` where `L = log(T/2π)`
— at T ~ 10²² this is `1/log(10²²/2π)² ≈ 1/2700`. A signal of size
10⁻³ should be detectable at 10⁵-10⁶ samples. Below our current
detection threshold for now but reachable.

**Cross-domain ingredient:** numerical analysis at high precision +
Odlyzko's zero database.

**Concrete first step:** download a sample of Odlyzko's zeros at
T ~ 10²² (publicly available). Re-run S49's battery on this sample.
Does the BK correlation flip sign or scale up?

**Failure profile:** zeros at high height are still GUE-distributed
to within sample noise. Not a *failure* per se — would refine E3.13
significantly to a much higher scale.

**Budget:** 1 session (just need to download data and re-run code).
**This one is genuinely tractable. Pick first if A-grade is the goal.**

### C2 — Higher-order arithmetic corrections (Conrey-Snaith)

**Question:** The Conrey-Snaith conjectures predict specific deviations
from GUE at *every* order n ≥ 2. Test orders 4, 5, 6 at N=8000.

**Why frontier:** order-3 already tested (S57); orders 4-6 untested at
this sample size. Even confirming GUE at order 4 would be the first
project-internal result on n-correlations beyond pair / triple.

**Cross-domain ingredient:** Wronskian/permanent computations on
zero ensembles.

**Concrete first step:** code the order-4 sine-kernel determinant
formula. Apply to the existing 8000-zero file.

**Failure profile:** orders 4-6 also GUE; but this would tighten
E7.1 to a higher order.

**Budget:** 2 sessions.

### C3 — Find a non-natural correlation that breaks GUE

**Question:** Construct a *bespoke* statistic on zeros — not pair, not
triple, not BK — that is sensitive to prime arithmetic and is
non-trivial.

**Why frontier:** the published statistics (pair, n-correlation, form
factor, BK) all match GUE. A statistic tailored to π(x) might break it.

**Cross-domain ingredient:** information theory / hypothesis testing
on statistical ensembles.

**Concrete first step:** consider the statistic
`S(γ_1, ..., γ_n) := Σ_p log(p) · cos(γ_i · log p)` for primes p ≤ Y.
Is its empirical distribution over zero windows consistent with the
GUE prediction (computable from sine-kernel)?

**Failure profile:** every statistic of zeros is a linear combination
of standard correlation measures, so the test reduces to one of those.

**Budget:** 1-2 sessions.

---

## §D. Cross-Domain Imports the Project Has Never Tried

The project has stayed within {analytic NT, complexity theory, basic
algebra}. Vast areas of mathematics are entirely untouched.

### D1 — Ergodic theory / dynamical zeta functions

**Question:** Define a dynamical system whose periodic orbits encode
primes (via the dynamical zeta function), then use Ruelle / Pollicott
results on exponential mixing to extract π(x).

**Why frontier:** the dynamical zeta function `ζ_φ(s) = exp(Σ Z_n s^n /n)`
of an Anosov flow has poles related to its prime-orbit-counting function
in a way analogous to ζ(s). Pollicott proved exponential decay of
correlations — maybe the analogous decay for "natural" Anosov flows
admits polylog extraction.

**Cross-domain ingredient:** smooth ergodic theory, transfer operators.

**Concrete first step:** read Pollicott-Sharp 1998 on prime orbit
theorems. Find a flow whose prime orbits are *literally* the integer
primes. (Probably impossible — but maybe asymptotically equivalent.)

**Failure profile:** no flow has integer primes as orbits without
encoding them externally (circular).

**Budget:** 2 sessions.

### D2 — Topological data analysis on the prime sequence

**Question:** Apply persistent homology to the prime indicator
sequence (or the gap sequence). Are there persistent topological
features at any scale that don't appear in matched random sequences?

**Why frontier:** TDA has identified non-obvious structure in
sequences considered random (RNA folding, neural data). Application
to primes is essentially unpublished.

**Cross-domain ingredient:** persistent homology (giotto-tda, ripser).

**Concrete first step:** compute persistence diagrams of `(p_n, p_{n+1})`
embedded as point cloud in ℝ². Compare to a Poisson(1/log x) point
cloud at the same density.

**Failure profile:** persistence diagrams of primes are statistically
indistinguishable from matched-density Poisson controls.

**Budget:** 1-2 sessions.

### D3 — Free probability + matrix models

**Question:** Treat the truth-table-of-π(x) as a Hermitian matrix
ensemble. Compute its spectral measure under various ensemble
constructions. Does the limiting spectral distribution match a known
free distribution? If yes, what does that say about its complexity?

**Why frontier:** non-commutative random variables (free probability)
are the right framework for "GUE-random sequences" in a precise sense.
C2 in NOVELTY_CHALLENGES.md gestured at this; D3 makes it the goal of
a session, not a sub-task.

**Cross-domain ingredient:** free probability (Voiculescu, Speicher).

**Concrete first step:** pick a Hermitian embedding of π(x) (e.g.,
shift-and-multiply by χ_P). Compute its empirical spectral distribution
at N up to 12. Compare to free Poisson, semicircle, Marchenko-Pastur.

**Failure profile:** the spectrum matches no known free distribution
and adds noise to the catalogue. Or matches semicircle, which would be
"random matrix is random."

**Budget:** 2 sessions.

### D4 — Quantum walks on prime graphs

**Question:** Define a quantum walk whose stationary distribution
encodes π(x) extraction. Analyse its mixing time.

**Why frontier:** classical random walks on the divisor graph have
been studied. Quantum walks have polynomial speedups in many search
problems. Application to π(x) extraction is unpublished.

**Cross-domain ingredient:** quantum walks (Childs et al.), discrete
quantum mechanics.

**Concrete first step:** define the divisor graph G_x where vertices
are integers in [1..x] and edges are (m, n) iff m|n. Define the
Szegedy walk on G_x. Compute its mixing time numerically for x ≤ 100.
Does the mixing time admit a polylog upper bound?

**Failure profile:** mixing time is Ω(x^c) for some c > 0.

**Budget:** 2-3 sessions.

---

## §E. Meta-Analysis of CLOSED_PATHS as a Data Source

The 730+ closures in CLOSED_PATHS.md are themselves data. Analysing
*patterns across closures* may reveal structural truths that no single
closure shows.

### E1 — Cluster CLOSED_PATHS by failure mode → infer the obstruction structure

**Question:** Apply unsupervised clustering to CLOSED_PATHS rows
(featurized by the techniques they tried). Do the clusters correspond
to *structurally different* obstructions, or are they all instances
of one or two underlying barriers?

**Why frontier:** would either confirm "the project's three failure modes
(C/E/I) are exhaustive" or surface a fourth latent failure mode the
project has missed.

**Cross-domain ingredient:** topic modeling (LDA, BERT embeddings)
on the closure text.

**Concrete first step:** embed each CLOSED_PATHS row using a
sentence-transformer model. Cluster with HDBSCAN. Inspect cluster
representatives.

**Failure profile:** clusters trivially correspond to failure-mode
labels, no new structure.

**Budget:** 1-2 sessions.

### E2 — Identify edge-pairs that have NEVER been jointly closed

**Question:** For each pair (E_i, E_j) of edges, does any
CLOSED_PATHS row cite BOTH? The pairs that have NEVER been jointly
cited are the unexplored composition spaces.

**Why frontier:** systematic composition discovery. May surface
non-obvious productive pairs.

**Concrete first step:** parse CLOSED_PATHS.md, extract cited edge
IDs per row, build the co-occurrence matrix. Visualise.

**Failure profile:** the unexplored pairs are unexplored for good
reasons (irrelevant edge combinations).

**Budget:** 1 session.

---

## §F. Synthesis Targets That Would Be Publishable

These are the highest-leverage outputs the project can produce. Each
is multi-session arc work.

### F1 — "Polylog π(x) is Computationally Hard: Four Structural Barriers"

The Three-Barriers paper (RESEARCH_AGENDA Arc 1) plus E7.6 (sieve
pebbling). Four family-level closures cleanly stated as theorems with
unified introduction. Target: arxiv preprint.

### F2 — "Pseudorandomness of π(x) mod m: A 35-Measure Empirical Battery"

Survey paper. The 35 measures are unique to this project. Target:
arxiv preprint. Co-target: a referee suggesting a 36th measure could
become a measurement-only follow-up session.

### F3 — "Information-Computation Gap for π(x): O(log x) bits, O(x^{2/3}) time"

Sharpen `novel/info_computation_gap.md` to paper-grade. Cross-reference
the 35-measure pseudorandomness as evidence; cross-reference E7.10 +
E5.8 + E7.11 as the technique exhaustion. Target: theoretical computer
science conference (STOC/FOCS workshop) or arxiv.

---

## How to escalate from §A-§F to a real new edge

If a frontier attack produces a result that survives the critique:

1. The result goes into `novel/<descriptive_name>.md` with full
   context, derivations, and falsification statement.
2. A new EDGES.md entry is added (likely §3 analytic, §5 conditional/TC⁰,
   or §7 negative-shape depending on content), with EVS rating.
3. RESEARCH_AGENDA.md gets a new arc if the result demands follow-up.
4. ATTACK_VECTORS.md entry is moved to "Closed attacks" section (add
   one at bottom) with a one-line outcome and pointer to the artefact.

---

## Closed attacks

### §C1 — Odlyzko zeros at height 10²² (CLOSED 2026-04-26, S71, mode I)

**Outcome:** structural negative result + new quantitative obstruction.

Re-ran the S49 BK arithmetic-correction probe on Odlyzko's published
tables `zeros4` (n ∈ [10²¹+1, 10²¹+10⁴], height T ≈ 1.4·10²⁰, L = 44.6)
and `zeros5` (n ∈ [10²²+1, 10²²+10⁴], T ≈ 1.4·10²¹, L = 46.8). Used a
**random-prime null** (template with primes replaced by uniform
pseudo-primes in [2, 50]) as the proper test instead of the biased
gap-shuffled null. Empirical Pearson with the BK template is
statistically indistinguishable from the random-prime null
(z = -0.94σ, +0.93σ); prime-frequency Fourier amplitudes are NOT
enhanced over random frequencies (ratios 0.95, 0.99).

Why §C1 fails by structure: the BK signal scales as `13.6/L²` while
empirical noise scales as `4/√N`. Detection requires
**N ≥ 0.09 κ² L⁴** for κ-σ Pearson detection. At Odlyzko's heights
this needs N ≈ 3·10⁵-4·10⁵ zeros; Odlyzko publishes 10⁴, short by
~35×. Pushing to higher heights makes the situation *worse* because
L⁴ grows faster than any reasonable N. The asymptotic regime
*suppresses* the BK signal faster than data accumulation can
compensate.

This sharpens E7.1 / E1.10 / E3.13 from qualitative "BK undetected"
to quantitative "BK detection requires N ≥ 0.81 L⁴ zeros" — a hard
scaling barrier independent of computational budget.

See `experiments/analytic/zeta_structure/odlyzko_high_height/`
and `archive/sessions/session71_c1_odlyzko_bk_probe.md`.

---

*A frontier attack that fails honestly is more valuable than ten
refinements that succeed. Be ambitious. Then be honest.*

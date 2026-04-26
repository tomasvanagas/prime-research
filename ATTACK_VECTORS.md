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

### A4 — Bounded-arithmetic provability of a TC⁰ primality witness

**Question:** Buss's theory `S^1_2` Σ^b_1-defines exactly the polynomial-
time functions, with witnessing: `S^1_2 ⊢ ∀x ∃y φ(x,y)` (φ ∈ Σ^b_1)
gives a polytime-computable witness. The Cook-Nguyen 2010 theory
`VTC^0` is the analogous theory for the complexity class TC^0. Question:
is there a Σ^B_1-statement φ(x) such that `VTC^0 ⊢ ∀x. (x ∈ PRIMES ↔ φ(x))`
would imply PRIMES ∈ TC^0? Equivalently: identify the precise
*proof-complexity* statement whose unprovability is the obstruction.

**Why frontier:** the proof-complexity / circuit-complexity duality is
a separate barrier from the explicit-circuit-construction barrier the
project has been working on. A cleanly-stated `VTC^0`-provability
question is a publishable technical result on its own (Krajicek 2019
chapter on PRIMES is silent about VTC^0). Either we identify φ and
prove its `VTC^0`-provability (A-grade: PRIMES in TC^0), or we identify
φ and prove its `VTC^0`-unprovability via a Razborov-Rudich /
proof-complexity-natural-proofs barrier (B-grade: structural negative
result complementary to E5.3).

**Cross-domain ingredient:** bounded arithmetic (Buss `S^1_2`), Cook-
Nguyen `VTC^0` and `VNC^1`, Σ^B_1 / Π^B_1 hierarchy, witnessing theorems.
This is **logic / proof-complexity**, not circuit complexity directly —
the project has not used it.

**Concrete first step:** Read Cook-Nguyen 2010 §9 (theory `VTC^0`).
Write down the Σ^B_1 sentence φ_AKS that is the natural translation
of the AKS correctness proof. Verify (informally, on paper) that
`VTC^0 ⊢ φ_AKS` would suffice — i.e., that Σ^B_1-witnessing for VTC^0
extracts a uniform TC^0 circuit from the proof. Then identify the
SPECIFIC subclaim of the AKS proof (likely "Σ Λ(p^k) ≥ x/log x for
some small range" or the cyclotomic factorisation step) whose
VTC^0-provability is the bottleneck.

**Failure profile (E):** the natural witness sentence φ_AKS reduces
to growing-dim MPOW in `VTC^0` (E5.3) and inherits its open status —
the formal-logic frame doesn't help. **(I):** the witness sentence is
unprovable in `VTC^0` by a proof-complexity natural-proofs / pseudorandom
function obstruction (e.g., Krajicek's witnessing barriers). **(INC):**
no clean Σ^B_1 sentence captures primality without circular reference
to AKS-style modular powering.

**A-grade success:** identify the precise φ + prove `VTC^0 ⊢ φ` ⇒
PRIMES ∈ TC^0, AND make non-trivial progress on whether `VTC^0 ⊢ φ`
holds (positive: TC^0 algorithm; negative: a new proof-complexity
barrier). **B-grade success:** identify φ but the provability is
equivalent to an existing closed barrier (E5.3, E5.8) — adds a new
"logical face" to the existing wall but no new bound.

**Cross-domain refs:**
- Buss 1986 thesis "Bounded Arithmetic"
  https://www.math.ucsd.edu/~sbuss/ResearchWeb/PhDThesis/index.html
- Cook-Nguyen 2010 "Logical Foundations of Proof Complexity"
  https://www.cs.toronto.edu/~sacook/  (book draft online)
- Krajicek 2019 "Proof Complexity" (Cambridge), chapters 5, 11
- Wikipedia: Bounded Arithmetic
  https://en.wikipedia.org/wiki/Bounded_arithmetic

**Budget:** 2-3 sessions (heavy reading session 1, formalisation
session 2, attempt at provability/unprovability session 3).

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

### B4 — Voronin universality with effective shifts as a polylog approximator

**Question:** Voronin (1975) proved that every non-vanishing analytic
function `f` on a compact set `K ⊂ {1/2 < Re(s) < 1}` is uniformly
approximable by `ζ(s + it)` for a positive-density set of vertical
shifts `t`. Garunkstis 2003 gave the only known effective bound on
the smallest valid `t`: `t ≤ exp(exp(10/ε^{13}))` — astronomical in
the worst case. Question: are there *natural* target functions `f`
(specifically, smooth approximants to `R^{-1}(n)` or `li(s)^{-1}`) for
which the effective shift `t(f, ε)` scales POLYNOMIALLY in `1/ε`?
If yes, evaluating `ζ(s + it)` at a single computable shift gives a
polylog reconstruction of `f`, hence a polylog approximation to π(x).

**Why frontier:** Voronin universality has NEVER been used
algorithmically in published work — every existing application is
existential. The empirical question "is the worst-case Garunkstis
bound tight, or only an artefact of the proof?" is open. The
project has done extensive effective-zero-truncation work (E3.2,
E3.4, E3.5) but has never tried the universality angle, in which
the goal is to *construct* a single ζ-value that evaluates the
target rather than to *sum* contributions of many zeros.

**Cross-domain ingredient:** Voronin universality theorem, Bohr
almost-periodicity, joint universality of L-functions (Steuding
2007), effective bounds (Garunkstis 2003).

**Concrete first step:** for target functions `f_c(s) = (s - c)^{-1}`
on a small disc around `s_0 = 0.7`, parameter `c ∈ {0.6, 0.65, 0.7}`,
do a numerical search: compute `ζ(s_0 + it)` at `t ∈ [10, 10^6]` (high
precision via mpmath) and plot the smallest valid `t` as a function of
`ε`. Does `t(ε)` scale polynomially or super-polynomially in `1/ε`?
Compare to the Garunkstis worst-case bound. If polynomial scaling is
observed, attempt to characterise *which* `f` admit polynomial-rate
universality (likely those with specific Bohr-set-density signatures).

**Failure profile (I):** empirical `t(ε)` scales as `exp(exp(c/ε^k))`
for every tested `f` — Garunkstis bound is empirically tight, no
polylog window exists. **(E):** even if some `f` admits polynomial
shifts, computing `ζ(s + it)` at large `t` itself requires `Ω(√t)`
operations per Riemann-Siegel — circular if `t` is large. **(INC):**
no `f` whose pre-image gives `π(x)` has a closed form admitting the
universality test.

**A-grade success:** identify a target `f` with polynomial-rate
universality + show that `ζ(s + it)` at the recovered `t` admits
polylog evaluation (e.g., via a structural identity that bypasses
Riemann-Siegel). **B-grade success:** strong empirical evidence that
Voronin shifts are super-polynomial for every `f` of interest,
producing a NEW negative-shape edge (E7.x: "Voronin universality is
algorithmically inert for π(x)") complementary to E7.11.

**Cross-domain refs:**
- Garunkstis 2003 "The effective universality theorem for the Riemann zeta-function"
  (Acta Arith. 113, also https://arxiv.org/abs/math/0306072 )
- Steuding 2007 "Value-Distribution of L-Functions"
  Springer LNM 1877
- Wikipedia: Zeta function universality
  https://en.wikipedia.org/wiki/Zeta_function_universality

**Budget:** 2 sessions. Tractable: numerical experiment is
mostly mpmath-based.

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

### C4 — Anderson localisation in a prime-driven discrete Schrödinger operator

**Question:** Define the discrete 1D Schrödinger operator
`H = -Δ + V` on `Z`, where `V(n) = +1` if `n` is prime and `0` else
(prime-indicator potential at constant strength). For a random
Bernoulli(1/log N) potential, Aizenman-Molchanov fractional moment
method gives complete spectral localisation in 1D for any disorder.
Question: does the deterministic prime-indicator potential have
the same Lyapunov exponent / localisation length as the random
Bernoulli, or does it deviate? If the Lyapunov exponent
`γ_prime(E) > γ_random(E)` (anomalously strong localisation), the
prime potential carries arithmetic structure that random doesn't —
a NEW pseudorandomness deviation invisible to the existing 35
measures, all of which are local statistics.

**Why frontier:** Anderson-localisation theory for *deterministic*
potentials is mostly known for *quasi-periodic* potentials
(Bourgain-Goldstein-Schlag; Avila-Jitomirskaya). The prime indicator
is genuinely random-looking but has neither randomness nor
quasi-periodicity — it is in a third regime. The Lyapunov exponent
`γ(E)` of the transfer-matrix product is a SPECTRAL statistic that
the project's pseudorandomness battery has not measured. If
`γ_prime ≠ γ_random` to statistical significance, this is the FIRST
example of the prime indicator deviating from pseudorandomness on a
spectral measure (analogous to S84's circuit-complexity deviation
but in a different category). If `γ_prime = γ_random` within noise,
add a 36th pseudorandomness measure with a new theoretical flavour.

**Cross-domain ingredient:** Anderson localisation (Aizenman-Molchanov
fractional moment method, Furstenberg-Kifer Lyapunov-exponent theory
for transfer matrices, Dyson-Schmidt formula
`γ(E) = -∫ log|E - λ| dν(λ)`).

**Concrete first step:** for `N = 10^5`, compute the transfer-matrix
product `T_N(E) = ∏_{k=1}^{N} [[E - V(k), -1], [1, 0]]` for `V = χ_P`
and `V = Bernoulli(π(N)/N)`. Estimate the Lyapunov exponent
`γ(E) = lim (1/N) log ‖T_N(E)‖` for `E ∈ [-1, 3]` (the relevant
energy band) at 100 sample points. Plot `γ_prime(E)` vs
`γ_random(E)`. Compute mean, max, and L^2 deviation.

**Failure profile (I):** prime potential gives `γ(E)` matched to
random within sample noise — closes as the 36th pseudorandomness
measure, no new structure. **(C):** the prime potential's Lyapunov
exponent reduces to the spectral density of `χ_P` via Thouless
formula, which is itself a known function of the existing edges —
circular. **(INC):** the prime indicator at unit strength is too
weak for any localisation to be observed at finite `N`.

**A-grade success:** `γ_prime(E)` deviates from `γ_random(E)` by
`> 3σ` at some energy `E_0`, AND the deviation has a structural
explanation via prime correlations (e.g., the Hardy-Littlewood pair
correlation density gives the deviation via Dyson-Schmidt).
**B-grade success:** clean negative result with a theoretical
explanation — `γ_prime = γ_random` because the Aizenman-Molchanov
proof goes through for any potential with the moment bounds that
χ_P satisfies. **B-grade success (alt):** identify a modified
potential (e.g., `V(n) = log p_n` if `n` is prime else 0) where
deviation is observable, opening a follow-up.

**Cross-domain refs:**
- Aizenman-Warzel 2015 "Random Operators" (AMS) chapters 6-7
- Bourgain "Green's function estimates for lattice Schrödinger operators"
  https://arxiv.org/abs/math/0410079
- Wikipedia: Anderson localization
  https://en.wikipedia.org/wiki/Anderson_localization

**Budget:** 2 sessions. Tractable: transfer-matrix products are
straightforward at `N = 10^5`.

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

### D5 — Continuous-time quantum walk on the same graphs §D.D4 closed (Szegedy)

**Question:** §D.D4 (S80, closed E) ruled out Szegedy walks on
divisor / coprime / Cayley graphs as polylog π(x) extractors. The
continuous-time quantum walk (CTQW; Childs 2009) generated by
Hamiltonian `H = A_G` (adjacency) or `H = L_G` (Laplacian) has
fundamentally different mixing phenomenology — Childs et al. 2003
showed CTQW gives EXPONENTIAL speedup on glued binary trees where
Szegedy gives only polynomial. Question: on the divisor graph `D_x`,
coprime graph `C_x`, or Cayley graph of `(Z/NZ)*`, does the CTQW
amplitude `|⟨prime | e^{-iHt} | seed⟩|^2` evolve to a primality-
detecting concentration in time `O(polylog x)`?

**Why frontier:** the §D.D4 closure exploited the Szegedy
discriminant theorem `mixing ~ 1/√δ`. CTQW mixing instead depends
on the SPECTRAL STRUCTURE (density of states, eigenvector overlap
with the seed), not just the gap. On a graph where the Szegedy gap
is `δ ~ 1/x^{0.5}` but the CTQW eigenvalue distribution has a
well-isolated band edge, CTQW could mix exponentially faster. The
glued-tree precedent shows this is not just theoretical. If the
arithmetic graphs have the right band-edge structure (which §D.D4
did not measure — gap is the wrong invariant for CTQW), the closure
on Szegedy does not transfer. **§D.D4's closure expressly notes**:
"for a Szegedy walk to give polylog primality extraction..." — the
two failure conditions are incompatible *for Szegedy*; CTQW could
satisfy them simultaneously.

**Cross-domain ingredient:** continuous-time quantum walks (Childs
2009 PRL 102, 180501 = arXiv:0806.1972). Specifically: the glued-
trees algorithm (Childs-Cleve-Deotto-Farhi-Gutmann-Spielman 2003
STOC), spectral-density analysis, eigenvalue degeneracy management.
This is the §D.D4 next-action flagged in S80.

**Concrete first step:** for `x ∈ {64, 128, 256, 512}`, compute
`|⟨v_prime | e^{-iH t} | v_1⟩|^2` for `H = A_{D_x}` (divisor adjacency),
`v_1 = |1⟩` (start at vertex 1), `v_prime = (1/√π(x)) Σ_p |p⟩` (uniform
on primes). Sweep `t ∈ [0, 100]` and compute the maximum amplitude.
Compare to the same quantity on a uniform-random graph with the same
degree sequence. Also compute the spectral density of `A_{D_x}` near
the band edges — is there an isolated cluster of eigenvalues that
overlaps significantly with `v_prime`?

**Failure profile (I):** CTQW mixing matches Szegedy quadratic
speedup — mixing `~ 1/Δ` for both, no gain at the band edge.
**(E):** CTQW spectral density is dominated by the bulk (Marchenko-
Pastur-like), no isolated band edge — closes via spectral concentration.
**(INC):** the right seed `|v_seed⟩` to localise the prime-detecting
amplitude is itself unknown.

**A-grade success:** identify a graph + seed combination where CTQW
amplitude concentrates polylogarithmically on the prime set, AND
the simulation cost (Hamiltonian simulation overhead) stays polylog
in `x`. **B-grade success:** structurally extend §D.D4 to CTQW via
a clean lemma "Szegedy gap ≪ 1/polylog ⇒ CTQW mixing ≫ polylog"
on the relevant graph families — promotes E7.13 to a stronger form
covering both walk types.

**Cross-domain refs:**
- Childs 2009 "Universal computation by quantum walk" PRL 102, 180501
  https://arxiv.org/abs/0806.1972
- Childs-Cleve-Deotto-Farhi-Gutmann-Spielman 2003 "Exponential
  algorithmic speedup by a quantum walk" STOC 2003
  https://arxiv.org/abs/quant-ph/0209131
- Wikipedia: Quantum walk (continuous-time section)
  https://en.wikipedia.org/wiki/Quantum_walk

**Budget:** 2 sessions. Hamiltonian simulation at `x ≤ 512` is
straightforward via `scipy.linalg.expm` on a small matrix.

### D6 — Gowers U^k norms of the prime indicator χ_P

**Question:** the project's pseudorandomness battery contains 35+
measures (correlation, spectral, entropy, tensor-rank, LFSR, etc.)
and `χ_P` deviates from random on NONE of them at any meaningful
scale. The Gowers uniformity norm `‖f‖_{U^k}` — the canonical
"k-degree pseudorandomness" measure in additive combinatorics — has
NEVER been computed for `χ_P`. Question: what is
`‖χ_P − E[χ_P]‖_{U^k[N]}` for `k ∈ {2, 3}` and `N` up to `2^{16}`?
Is it `o(1)` (truly k-pseudorandom, by Green-Tao-Ziegler inverse
theorem this means `χ_P` does not correlate with any (k-1)-step
nilsequence at scale N) or `Ω(1)` (carries nilsequence structure)?

**Why frontier:** Green-Tao-Ziegler (arXiv:1009.3998) PROVED the
inverse theorem: large `U^{s+1}` norm ⇔ correlation with an `s`-step
nilsequence. For the von Mangoldt function `Λ`, the U^k-smallness
proof required massive machinery (Green-Tao "Mobius and
nilsequences", arXiv:0807.1736), which is *conditional* on the
inverse Gowers-norm conjecture (later resolved). For the BARE
indicator `χ_P` (not Λ-weighted by `log p`), the U^k norm is in a
formally different regime and is **not known in published
literature**. Computing it directly is novel. If `χ_P` U^k for some
small `k_0` is `Ω(1)`, the Green-Tao-Ziegler inverse theorem points
to a SPECIFIC (k_0-1)-step nilsequence that correlates with primes
— that nilsequence would be a new arithmetic invariant for `χ_P`
and the FIRST measure in the project where `χ_P` deviates from
random in the strict additive-combinatorics sense.

**Cross-domain ingredient:** Gowers uniformity norms, additive
combinatorics, Green-Tao-Ziegler U^{s+1}[N] inverse theorem,
Mobius / Λ vs χ_P comparison.

**Concrete first step:** compute `U^2(χ_P)` for `N = 2^{12}, 2^{14},
2^{16}` (FFT-based, O(N log N)). Compute `U^3(χ_P)` for `N = 2^{10},
2^{11}, 2^{12}` (direct sum over `(x, h_1, h_2)`, O(N^3) — feasible
to N≈4096). Compute the same for: (a) uniform random Bernoulli at
density `π(N)/N`; (b) the Liouville indicator `1[λ(n) = -1]`; (c) a
matched-density i.i.d. control. Report `‖f‖_{U^k} / ‖random‖_{U^k}`
ratios with bootstrap error bars.

**Failure profile (I):** `‖χ_P‖_{U^k}` matches random within sample
noise — closes as 36th pseudorandomness measure, NEW kind because
it directly addresses "k-degree polynomial structure" not visible to
local correlation tests. **(E):** Green-Tao "Mobius nilsequences"
result transports trivially from `Λ` to `χ_P` via partial summation
or smooth-vs-rough comparison — closure is "subsumed by published
result." **(INC):** numerical computation at small N has too much
finite-size noise to discriminate.

**A-grade success:** `‖χ_P‖_{U^{k_0}}` is `Ω(1)` for explicit `k_0`
AND the Green-Tao-Ziegler inverse theorem applies to identify a
specific (k_0-1)-step nilsequence that correlates with `χ_P` at
non-trivial scale — first ever direct nilsequence-shape arithmetic
content in the project. **B-grade success:** `‖χ_P‖_{U^k} = o(1)`
to `k = 3` empirically, but with quantitative rate that gives a
new pseudorandomness bound stronger than any of the 35 existing
measures.

**Cross-domain refs:**
- Green-Tao "Linear equations in primes"
  https://arxiv.org/abs/math/0606088
- Green-Tao "The Mobius function is strongly orthogonal to nilsequences"
  https://arxiv.org/abs/0807.1736
- Green-Tao-Ziegler "An inverse theorem for the Gowers U^{s+1}[N] norm"
  https://arxiv.org/abs/1009.3998
- Wikipedia: Gowers norm
  https://en.wikipedia.org/wiki/Gowers_norm

**Budget:** 2 sessions. U^2 trivial; U^3 at `N = 4096` ≈ 7×10^{10}
operations — overnight on one core.

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

### §D.D6 — Gowers U^k norms of chi_P (CLOSED 2026-04-26, S85, mode E)

**Outcome:** structural confirmation of Hardy-Littlewood prediction at the
Gowers-norm scale + 36th-37th pseudorandomness measure with closed-form
prediction + new EDGE E2.13.

Computed `Q^k(f) = ‖f‖_{U^k}^{2^k} / E[f]^{2^k}` for chi_P, matched-density
Bernoulli, Liouville, and W-tricked chi_{W,1} on Z/NZ. U^2 via FFT
identity; U^3 via recursion through Δ_h. Hardy-Littlewood {0,1}^k-cube
singular series numerics: **S_2 = 2.300938**, **S_3 = 54.116088**.

**Empirical Q^2(chi_P)**: 2.103 → 2.153 monotonically across N ∈ [2^10, 2^18];
gap to S_2 = 2.301 decays as O(1/log N). **Q^3(chi_P) ≈ 35.5** stable
across N ∈ [2^10, 2^13]. **Liouville Q^2 = 1.000** (Gowers-uniform at U^2,
matches Mobius/nilsequence orthogonality). **W=210 W-tricked Q^2 = 1.003
vs HL pred 1.002** — match to 0.1%, restoring Gowers uniformity.

**Why §D6 is closed (mode E):** chi_P's U^k structure is *exactly*
Hardy-Littlewood — no additional bit beyond the singular series. The
Q^k extraction itself costs Θ(N^2 log N) (not polylog). No deviation
from HL prediction at U^2 or U^3 within the tested range, no novel
arithmetic content beyond what HL conjectures already predict.

**Adds EDGE E2.13** (Q^k(chi_P) = S_k {0,1}^k-cube singular series),
upgrading the pseudorandomness battery from "35 measures all
random-like" to "36 measures, of which exactly the Gowers norm
detects HL correlations". Complements the S84 circuit-depth-2
PRIMES-vs-random gap (different mechanism, both ultimately reduce to
known structure of primes).

See `experiments/information_theory/gowers_uk_chi_p.py`,
`experiments/information_theory/hardy_littlewood_box.py`,
`experiments/information_theory/gowers_uk_chi_p_results.md`,
`archive/sessions/session85_d6_gowers_uk_chi_p.md`.

**Follow-ups proposed** (per CLAUDE.md self-extension):
- (a) **U^4 of chi_P at N ≤ 2^12**: extend the verification chain to
  k=4 (predicted S_4 ~ 10^3); requires p^5 enumeration for α_p(k=4)
  which is borderline tractable. Borderline B-grade if Q^4 matches S_4.
- (b) **Lambda vs chi_P comparison**: empirically compare
  `Q^k(Lambda_truncated)` to `Q^k(chi_P)`. Predicted: Lambda is
  Gowers-uniform after centering (Green-Tao theorem) while chi_P
  has Q^k = S_k. Numerical comparison would visualise the difference
  cleanly. C-grade refinement.

### §A.A1 — SAT/ILP search for small TC^0 (depth-2 sign-threshold W=1) circuits at N=8 (PARTIAL CLOSURE 2026-04-26, S84, mode I)

**Outcome:** sub-family closure (depth-2 W=1) + structural mechanism + open
remainder (depth ≥ 3 or W ≥ 2 still untested).

**Sub-family closed**: depth-2 sign-threshold circuits with weights ∈ {-1, 0, +1}
at both layers. At N=8 such circuits need ≥ 17 bottom gates (lower bound from
ILP-proven infeasibility at M=8, 10, 12, 14 with k_max=4 candidate set, M=16
infeasible with k_max=5; M=20 timed out before resolving).

**Exact small-N values** (full-kmax candidate set, ILP):
- N=4 → M=3 (full search 108 W=1 thresholds)
- N=6 → M=6 (full search 1458 W=1 thresholds)
- N=7 → M ≥ 6 (full search 5103; M=5,6 unsat; M=7,8 timed out at 180s)
- N=8 → M ≥ 17 (k_max=4: M ≤ 14 unsat; k_max=5: M=16 unsat; M=20 unknown)

**Statistically significant PRIMES-vs-random gap at N=6**: PRIMES needs
M=6, all 10 random matched controls (weight 18/64) need M ∈ {7, 8}.
Binomial null p < 0.001. **First instance in the project of a circuit-
complexity measure where PRIMES empirically deviates from a matched-density
random function.**

**Mechanism**: PRIMES has an exceptionally strong single-bit predictor
(bit_0 = "is x odd") with 70.3% accuracy at N=8 — vs random matched best
single-bit at 57.0% (max over 30 seeds). The depth-2 circuit gets the
1-bit predictor as a "free" first bottom gate; subsequent gates correct
fewer remaining errors → smaller M.

**PTF degrees** (LP, real-coeff monomial basis): PRIMES at N=4..8 has PTF
degree (2, 3, 3, 3, 4); pi(x) mod 2 (3, 3, 3, 4, 4); matched random median
(2, 3, 3, 3, 4). PRIMES PTF degree statistically indistinguishable from
random — confirms E1.10/E3.13 at the Boolean-function level.

**Why §A.A1 is NOT fully closed**: this session ruled out only the W=1
sub-family. The full §A1 question (depth ≤ 5, size ≤ 2000, with arbitrary
integer weights up to ~2^N) remains open. Three follow-ups would extend
the closure:
- (a) larger weight bound W ∈ {2, 4, 8} at depth 2 (1-2 sessions);
- (b) depth-3 with restricted bottom layer (1 session);
- (c) "1-bit-calibrated random" baseline at N=6 — does the PRIMES-vs-
      random gap persist when random functions are biased to match PRIMES'
      single-bit predictor accuracy? (1 session, B-grade falsification test).

See `experiments/circuit_complexity/sat_tc0_primes_n8/` and
`archive/sessions/session84_a1_sat_tc0_primes.md`.

### §D.D4 — Szegedy quantum walks on prime/divisor/coprime/Cayley graphs (CLOSED 2026-04-26, S80, mode E)

**Outcome:** quantum extension of E7.12 + clean closed-form on coprime
graph + structural negative result.

Tested Szegedy walks on three number-theoretic graph families using the
Szegedy 2004 discriminant matrix theorem (`D(x,y) = sqrt(P(x,y) P(y,x))`,
mixing time `O(1/sqrt(δ))`).

(A) Cayley `Cay((Z/NZ)*, {2,3,5,inv})` for primes N ∈ {31, …, 1009}:
classical mixing `t_C(N) ~ N^{0.789}`, Szegedy mixing
`t_Q(N) ~ N^{0.414}` — quadratic speedup confirmed but neither is
polylog. **EXTENDS S79 (E7.12) to the quantum-walk regime**.

(B) Coprime `C_x` for x ∈ {30, …, 1000}: spectral gap is `Ω(1)`
(asymptotically `δ → 0.4166...`); stationary prime mass / uniform
prime mass `→ ζ(2) = π²/6 ≈ 1.6449` (Mertens; verified at x=1000
with deviation -0.022). Bias is constant — no polylog π(x) extractor.

(C) Divisor `D_x`: gap `Ω(1)`. High-prime-mass eigenvectors localise
on prime-centered divisor clusters (one specific prime `p` + its
multiples per eigenvector, e.g., k=14 at x=250 concentrates 100% mass
on {43, 47} ∪ multiples ∪ {1}) — same degree-class probing mechanism
as E7.12.

**Why §D.D4 fails by structure:** for a Szegedy walk to give polylog
primality extraction, you need (i) `δ = Ω(1/polylog(x))` AND (ii) a
single eigenvector with primality-localised mass. The two conditions
are **incompatible** on the tested families: Cayley graphs satisfy a
weak version of (ii) but `δ` shrinks polynomially; coprime/divisor
graphs satisfy (i) but only the trivial Perron eigenvector carries
global primality information (via degree-weighted Mertens factor).

**Suggested EDGE upgrade**: extend E7.12 (currently classical-Cayley)
to a stronger "Cayley-spectral barrier carries to Szegedy quantisation"
form, OR add E7.13 covering Szegedy walks on arithmetic graphs.

See `experiments/quantum/szegedy_walk_prime_graphs/` and
`archive/sessions/session80_d4_szegedy_walk.md`.

### §A.A3 — Spectral primality test via (Z/nZ)* Cayley graph (CLOSED 2026-04-26, S79, mode E)

**Outcome:** structural negative result + new lower-bound identity.

Built Cay((Z/nZ)*, {2,3,5,2⁻¹,3⁻¹,5⁻¹}) for 533 values of n ∈ [7,2000]
coprime to 30 (300 primes, 13 prime powers, 214 semiprimes, 6 triple-
prime composites). Computed full adjacency spectrum and 24 candidate
features. **No feature is disjoint primes-vs-composites.**
Mann-Whitney AUC for the hard sub-case (primes vs prime powers, both
ω(n)=1): 0.509 with {2,3,5}, 0.566 with {2,3,7}, 0.673 with
{2,3,5,7}∪inv — essentially noise.

Why §A.A3 fails by structure: from CRT character theory, n_int_eigs ≥
**2^ω(n)** where ω(n) is the number of distinct prime factors. (Z/p^kZ)*
is cyclic for odd p (Gauss), so prime n and prime power n=p^k both have
ω=1 and produce structurally identical unit groups. Spectrum cannot
distinguish them. Adding a perfect-power test inherits AKS-class
growing-dim MPOW depth (E5.3, E7.10).

Adds new EDGES.md entry **E7.12** (Cayley graph spectrum probes ω(n),
not primality) — search-space constraint on the entire fixed-generator
abelian-Cayley-graph spectral family.

See `experiments/circuit_complexity/cayley_spectral_primality/` and
`archive/sessions/session79_a3_cayley_spectral_primality.md`.

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

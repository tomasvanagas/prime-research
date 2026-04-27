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

### A5 — Maynard multidimensional sieve weight as a TC⁰ primality witness  [CLOSED S116, see "Closed attacks"]

**Question:** Maynard's 2015 multidimensional sieve produces an explicit
positive weight `w(n) = (Σ_{d_i | n+h_i} λ_{d_1,...,d_k})^2` (with
`λ_{d_1,...,d_k} = F(log d_1 / log R, ..., log d_k / log R) · μ(d_1)…μ(d_k)`
for an optimised polynomial `F` on the simplex `Σ x_i ≤ 1`) such that
`Σ_{N ≤ n ≤ 2N} w(n) χ_P(n+h_i) > 0` for some `i ∈ [k]`, proving bounded
gaps unconditionally. Question: is the sieve weight `w(n)` itself
TC⁰-computable, AND does a calibrated threshold `w(n) > τ*` correlate
with primality of `n+h_i*` strongly enough to give a polylog-time
primality witness?

**Why frontier:** the Maynard sieve is the **most refined explicit
prime-detection machinery** in modern analytic number theory. By
construction `w(n) χ_P(n+h_i) ≥ 0` and the sum is provably positive —
it IS a primality detector at the *aggregate* level. The CRITICAL
question, never asked in the published literature, is whether the
weight itself, evaluated at a SINGLE n, has a structured TC⁰-realisable
representation. Pascadi 2025 (E3.12) extended Maynard-Tao to `x^{5/8}`
level of distribution, raising the urgency. If `w` is TC⁰-computable
AND `w(n) > τ*` separates primes from composites with high
precision-recall, this is the first known TC⁰ primality test outside
the AKS family and would resolve PRIMES ∈ TC⁰ unconditionally.

**Cross-domain ingredient:** Maynard 2015 multidimensional GPY sieve,
specifically the optimal weight construction (Maynard's Lemma 4.1) and
its companion Selberg-sieve divisor enumeration. The TC⁰-side is
divisor enumeration up to `R = N^{θ/2}` and arithmetic on rational
coefficients of bounded denominator.

**Concrete first step:** for `k = 3`, `R = N^{0.1}`, fix
`H = {0, 2, 6}` (admissible 3-tuple). Use Maynard's optimal `F*` for
k=3 (a degree-2 symmetric polynomial maximising the eigenvalue ratio
of the M-matrix from §6.1). Tabulate `w(n)` for
`n ∈ [10^4, 10^4 + 5000]`. Compare distributions of `w(n)` for `n`
such that `(n, n+2, n+6)` contains ≥ 1 prime vs `n` such that none of
`{n, n+2, n+6}` is prime. Compute Mann-Whitney AUC and threshold `τ*`
maximising (precision × recall). Count operations to evaluate `w` at
one `n`: divisors of `n+h_i` up to `R = n^{0.1}` times polylog overhead
per divisor — is this `polylog(n)` (unlikely except for divisor-bounded
n) or `n^{0.05}` (typical) or `n^{0.1}` (worst case)?

**Failure profile (E):** evaluating `w(n)` requires enumerating
divisors of `n` up to `R = n^{0.1}`, costing `Ω(n^{0.05})` per `n` via
Eratosthenes — sub-poly but not polylog. Closes as "Maynard sieve is a
sub-polynomial primality witness, not polylog" — sharpens E6.7
(sieve-pebbling) by showing Maynard does NOT escape the divisor-
enumeration barrier. **(I):** `w(n)` distributions overlap strongly
between primes and non-primes — `w` does not separate at single-n
resolution; only its window average over [N, 2N] is non-zero
(Cauchy-Schwarz averaging). **(INC):** `F*` for k > 3 requires
solving a multivariate quadratic eigenvalue problem itself; the
weight is implicit, not closed-form.

**A-grade success:** `w(n)` is TC⁰-computable AND `w(n) > τ*` has
> 0.95 precision and recall for primality on a held-out test set of
size ≥ 10^6 — gives PRIMES ∈ TC⁰ unconditionally; resolves the
project's only open problem.

**B-grade success:** explicit lower bound `op(w) ≥ Ω(d(n))` showing
the Selberg-sieve sum cannot escape divisor enumeration — produces a
new negative-shape edge "Maynard weight requires `Ω(d(n))` ops per
single-n evaluation," refining E6.7 by closing the **most refined
known prime-detection sieve** to the same `x^{0.034}`-style barrier.

**Cross-domain refs:**
- Maynard 2015 "Small gaps between primes" Annals of Math 181, 383
  https://arxiv.org/abs/1311.4600
- Polymath8b 2014 "Variants of the Selberg sieve, and bounded
  intervals containing many primes" arXiv:1407.4897
- Pascadi 2025 (E3.12, see `literature/state_of_art_2026.md`)
- Wikipedia: Maynard's theorem (gap conjectures)
  https://en.wikipedia.org/wiki/Yitang_Zhang#Bounded_gaps_between_primes

**Budget:** 2-3 sessions. Session 1: implement Maynard weight
evaluator + AUC test at k=3. Session 2: depth analysis of the
evaluator (TC⁰ feasibility). Session 3 (only if first two yield
A-grade signal): scale to k=5 and N = 10^6 — and at that point the
project should escalate to verify mode.

### A6 — Quantitative reverse mathematics: proof-theoretic strength of π(x) error-term bounds

**Question:** the Friedman-Simpson reverse-mathematics program classifies
mathematical theorems by the weakest subsystem of second-order arithmetic
in which they are provable: `RCA_0 ⊂ WKL_0 ⊂ ACA_0 ⊂ ATR_0 ⊂ Π^1_1-CA_0`.
The PNT itself is provably in `WKL_0` (Simpson textbook §III.1, also
Avigad-Friedman 1998). However, **the proof-theoretic strength of
QUANTITATIVE bounds on `π(x) − Li(x)` is not characterised in the
published reverse-math literature**. Specifically: what is the weakest
subsystem proving the de la Vallée Poussin error term
`|π(x) − Li(x)| ≤ x · exp(-c · √log x)` for some explicit `c > 0`?
Is it `WKL_0` (giving polynomial-time witness extraction by a Friedman-
Sieg conservativity result over `Sigma^0_1`-witnessing) or does it
require `ACA_0` (giving only computable but not polynomial-time
witness)? What about the Korobov-Vinogradov bound
`exp(-c (log x)^{3/5} (log log x)^{-1/5})`, which uses sharper
zero-free-region estimates?

**Why frontier:** the Friedman-Sieg bounded-arithmetic conservativity
results convert proof-theoretic strength into computational strength:
- a `Π^0_2` statement provable in `WKL_0` admits a recursive witness.
- a `Π^0_2` statement provable in `RCA_0` admits a *primitive recursive*
  witness — and combined with `BΣ^0_1`-bounded induction, a polytime
  witness in many cases.
- a statement requiring `ACA_0` typically only has `Δ^0_2` witness, no
  polytime structure.
**No published work characterises the proof-theoretic strength of
quantitative π(x) bounds.** Friedman 1976 noted PNT in WKL_0 without
addressing the error term. If the Vinogradov-Korobov bound is provable
in `WKL_0`, then by the conservativity theorem there is a polytime
algorithm that on input `x` produces a witness `n` such that the bound
holds — this would not be polylog π(x) per se, but it would be a
*proof-theoretic upper bound* on the witness complexity of the bound,
the FIRST such characterisation. If Korobov-Vinogradov is *unprovable*
in `WKL_0`, that's a B-grade structural barrier complementary to E5.3.

**Cross-domain ingredient:** Reverse mathematics (Friedman-Simpson),
specifically the analytic-NT subprograms in Simpson 2009 *Subsystems
of Second-Order Arithmetic* (Cambridge) and Friedman-Sieg 1989. The
key tools: (i) `WKL_0`'s connection to Π^0_1 conservativity over
`PRA`; (ii) Avigad's reverse-math analysis of Dirichlet's theorem
(Avigad-Towsner 2012); (iii) the Simpson-Tanaka-Yamazaki analysis of
the analytic continuation of ζ. UNUSED in the project (§9 row).

**Distinction from line 218 closure:** `CLOSED_PATHS.md` line 218
("Reverse mathematics — p(n) provable in RCA_0, barrier is computational")
closed the QUALITATIVE existence-of-primes question. A6 attacks the
QUANTITATIVE bound's proof strength, which is a fundamentally different
proof-theoretic question. The line 218 closure does not characterise
which subsystem proves any specific `π(x) − Li(x)` rate.

**Concrete first step:** read Simpson 2009 Chapter III (mathematics
in `WKL_0`) and Avigad-Friedman 1998 "Combining decision procedures
for the reals". Identify the precise step in the de la Vallée Poussin
proof of PNT-with-error-term where ACA_0-strength reasoning is invoked
(typically: existence of zero-free regions via complex analysis, or
extraction of the implied constant `c` in the error term). Determine
whether that step admits a `WKL_0`-conservative reformulation. If yes,
extract a uniform polytime witness function for the bound. If no,
identify the precise instance where ACA_0 is necessary and characterise
the obstruction.

**Failure profile (E):** the Vinogradov-Korobov error-term proof
*essentially* uses ACA_0-strength complex analysis (via residue
calculus on the zero-free region), and no `WKL_0` reformulation is
possible — closes as "quantitative π(x) bounds inherit ACA_0
proof strength, complementing the explicit-formula barrier."
**B-grade.** **(I):** the bound is provable in `WKL_0` but the
implied polytime witness is not polylog (e.g., O(x^{1/2}) or
x^{2/3}) — adds proof-theoretic confirmation of existing time bounds
but no new complexity insight. **(C):** the proof-strength
characterisation reduces to a `Π^0_1` statement equivalent to
unprovable-in-PA — circular.

**A-grade success:** identify a quantitative π(x) bound provable in
`WKL_0` that gives a NEW polytime witness extraction strictly improving
existing computational time bounds, OR prove that the Korobov-Vinogradov
bound is unprovable in `WKL_0` via a Friedman-style conservativity
argument — first proof-theoretic obstruction to polylog π(x).

**B-grade success:** complete proof-strength characterisation of a
specific quantitative π(x) bound, with explicit reduction to a known
subsystem — first reverse-math result on quantitative prime
counting, complementary to A4 (VTC^0 circuit-complexity proof
strength) but in the `WKL_0`/`ACA_0` hierarchy.

**Cross-domain refs:**
- Simpson 2009 *Subsystems of Second-Order Arithmetic* Cambridge
  (Chapters II-V cover RCA_0, WKL_0, ACA_0)
- Avigad-Friedman 1998 "Combining decision procedures for the reals"
  https://www.andrew.cmu.edu/user/avigad/Papers/dpr.pdf
- Friedman-Sieg 1989 "Inductive definitions, constructive ordinals,
  and normal derivations" Lecture Notes in Math 1320
- Avigad-Towsner 2012 "Functional interpretation and inductive
  definitions" J. Symbolic Logic 77, 357
- Wikipedia: Reverse mathematics
  https://en.wikipedia.org/wiki/Reverse_mathematics

**Budget:** 2-3 sessions. Session 1: literature read (Simpson Ch III,
Avigad-Friedman). Session 2: identify the ACA_0-invoking step in
the standard PNT-with-error proof. Session 3: attempt either
WKL_0 reformulation or non-conservativity argument.

---

## §B. Spectral Identity Searches Beyond the Standard Bases

E2.11 + E2.12 + S29 + S68 (Bessel) closed the algebraic identity
search across dozens of bases. The frontier is bases tied to *specific
modern published mathematics* the project has not engaged with.

### B1 — Polynomial Method Compression (Croot-Lev-Pach style)  [CLOSED S92, see "Closed attacks"]

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

**Outcome (S92, mode E):** algebraic immunity AI_F_2(chi_P) = 2 for all
N in [4,13] via the explicit annihilator g(x) = (1+x_0)(1+x_1)
encoding the trivial mod-4 sieve fact. F_p multilinear ANF degree
near-saturated; slice rank brackets non-informative (p=2) or match
random (p=3, k>=3). W-trick at W >= 6 fully removes the AI deviation.
Adds **EDGE E2.15**, third independent confirmation of "chi_P
structure = HL equidistribution mod q" alongside E2.13 (Gowers, S85)
and E2.14 (Anderson, S88). See `experiments/algebraic/algebraic_immunity_chi_p/`
and Closed attacks below.

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

### B5 — Beurling generalised primes: minimal-perturbation polylog test

**Question:** A Beurling generalised prime system is an abstract
multiplicative semigroup `(P, ·)` of positive reals `1 < p_1 ≤ p_2 ≤ …`
together with the integers it generates `N(P) = {p_{i_1} p_{i_2} …}`,
yielding counting functions `π_P(x) = #{p_i ≤ x}` and `N_P(x) =
#{n ∈ N(P) : n ≤ x}`. Beurling 1937 proved that `N_P(x) = Cx + O(x/log^γ x)`
for `γ > 3/2` ⇒ PNT analog `π_P(x) ~ Li(x)`. The CONVERSE problem
(Diamond-Zhang 2007 *Trans. AMS* 359, Hilberdink-Lapidus 2016
*Funct. Approx.* 54): for any prescribed regularity rate on `N_P(x)`,
how close can `π_P` be to integer primes? Specifically: **is there a
Beurling system `P` with (i) `π_P(x)` is `polylog(x)`-evaluable AND
(ii) symmetric difference `|{p ∈ P : p ≤ x} △ {p prime : p ≤ x}| =
o(π(x)/log x)`?** If yes, the integer-prime sequence is in a polylog-
evaluable Beurling neighbourhood — strong evidence that *information*
in π(x) compresses to polylog up to negligibly many primes.

**Why frontier:** Beurling generalised primes are a **one-parameter
deformation family** of the integer primes — perturbations parameterised
by a regularity index `γ` on the integer-counting function `N_P(x)`.
The published literature (Diamond-Zhang, Hilberdink-Lapidus, Pintz)
characterises the analytic regime where Beurling π_P inherits PNT-style
asymptotics, but the **algorithmic / polylog question has never been
posed**: which Beurling systems near the integers admit polylog π_P
extractors? This lifts the problem to a moduli-space question — if the
integer primes are an isolated point in the polylog-evaluable Beurling
locus, the compressibility is "robust"; if they are *not* in the
polylog locus but the locus is dense, the integer-prime hardness has
a *transcendental* character. The frontier separates two cases:
(a) all Beurling systems sufficiently close to integer primes inherit
the same polylog hardness (rigidity result); (b) most close-by Beurling
systems are polylog-evaluable, and integer primes are an exceptional
non-polylog point (singularity result).

**Cross-domain ingredient:** Beurling generalised primes (Beurling 1937
Acta Math. 68; Hilberdink 2009 J. Number Theory 129; Hilberdink-
Lapidus 2016 *Functiones et Approximatio* 54). Specifically: the
**Diamond-Zhang correspondence** between abstract sieves and Beurling
systems, and the **Hilberdink quasi-Beurling spaces** characterising
deformation of integer primes via continuous regularity index `γ`.
UNUSED in §10. Distinct from the Nyman-Beurling RH criterion (a
different concept; closed in `experiments/analytic/analytic_nt_new_results.md`)
and from the Beurling-Selberg uncertainty (CLOSED line 33 — also a
different concept, sieve-theoretic).

**Concrete first step:** for parameter ε ∈ {0.1, 0.05, 0.01}, construct
the Beurling system `P_ε := {p ∈ primes : p ≤ N} ∪ {q : q = p · (1 + δ_p),
|δ_p| < ε, p prime, p ≤ N}` — a perturbed prime system where each
prime is shifted by at most ε in log-scale. Verify: (a) Beurling
counting `N_{P_ε}(x)` admits the regularity required for PNT;
(b) compute `π_{P_ε}(x)` directly for `N = 10^6` via sieve over
log-shifted primes; (c) test whether **the perturbation enables a
polylog evaluator** — does the log-shifted system `P_ε` have any
algebraic structure (rational shifts? geometric shifts?) that gives
its π_P a polylog formula? Specifically, try `P_ε = {p^{1+ε_p}}` for
small `ε_p` (Hilberdink's "scale-perturbed" family) and ask if any
Hilberdink-class system with deformation parameter `γ < 3/2` has a
polylog π_P extractor in finite N. Compare empirical `π_{P_ε}(x)` to
true `π(x)` and to the singular series prediction.

**Failure profile (E):** every Beurling perturbation that PRESERVES
the polylog-evaluable property reduces to a deterministic-shift
parametric family (e.g., `p_i → p_i + c · log p_i`) which is
trivially polylog-evaluable but is at integer-distance > π(x)/log x
from the true primes — closes as "polylog Beurling locus does not
intersect any o(π/log)-neighbourhood of integer primes." **B-grade**:
new structural rigidity result. **(I):** every Beurling system
passing the regularity test has π_P matching π up to lower-order
terms, but the polylog property is not inherited from analytical
regularity — closes as "Beurling regularity ≠ algorithmic
compressibility." **(INC):** the deformation-space topology is
insufficiently characterised in the published literature to make
the question precise.

**A-grade success:** construct a Beurling system `P_*` with
`|π_{P_*}(x) − π(x)| = o(π(x)/log x)` AND a polylog-evaluable
`π_{P_*}` formula — first known polylog perturbative neighbourhood of
the integer primes; would either give a polylog π(x) algorithm via
δ-correction OR a structural argument why integer primes are
*algorithmically isolated* in Beurling deformation space.

**B-grade success:** clean rigidity theorem "the Beurling polylog
locus has positive distance from integer primes in any reasonable
metric on prime-counting functions" — first reverse-direction result
on Beurling-polylog correspondence, complementing E7.10 (AKS
modulus-twist orthogonality) by showing the *orthogonal* deformation
direction (counting-regularity) also doesn't give polylog.

**Cross-domain refs:**
- Beurling 1937 "Analyse de la loi asymptotique de la distribution
  des nombres premiers généralisés" Acta Math. 68
- Diamond-Zhang 2007 "A PNT equivalence for Beurling numbers"
  Trans. AMS 359
- Hilberdink 2009 "An arithmetical mapping and applications to Ω-results
  for the Riemann zeta function" J. Number Theory 129
- Hilberdink-Lapidus 2016 "Beurling zeta functions, generalised primes,
  and fractal membranes" *Functiones et Approximatio* 54
  https://projecteuclid.org/euclid.facm/1428589969
- Wikipedia: Generalised primes
  https://en.wikipedia.org/wiki/Generalized_Riemann_hypothesis (mentions
  Beurling-zeta in passing); see also https://en.wikipedia.org/wiki/Riemann_hypothesis

**Budget:** 2-3 sessions. Session 1: read Hilberdink-Lapidus 2016 +
Diamond-Zhang 2007 + characterise the regularity-vs-polylog dichotomy
on paper. Session 2: numerical experiments with parametrised Beurling
families (log-shift, scale-perturbation). Session 3: attempt either
A-grade construction or rigidity proof.

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

### C4 — Anderson localisation in a prime-driven discrete Schrödinger operator [CLOSED S88, see "Closed attacks"]

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

### C5 — Stein's method: quantitative finite-x Gaussianity of `(π(x) - Li(x))/(√x/log x)`  [CLOSED S108, see "Closed attacks"]

**S108 outcome (provisional A-grade pending verify):** plateau detected
with `W_1(D̂_K, N(μ̂, σ̂²)) ≈ 0.0083` at K=10000 over `x ∈ [10^6, 10^7]`,
not shrinking as `1/√K`; structurally explained (correlation 0.906
across sub-windows) by the explicit-formula low-zero contribution.
Closure mode: **E** (reduces to E1.5 explicit-formula, joins GUE-sieve-
circuit closure family). See "Closed attacks" below.

**Question:** Under RH, the renormalised deviation
`D(x) := (π(x) - Li(x)) · log x / √x` is conjectured to converge
weakly to a centered (non-degenerate) distribution as `x → ∞`. Stein's
method provides EXPLICIT QUANTITATIVE bounds on the Wasserstein-1
distance `W_1(D̂_X, N(0, σ̂²))` for any FIXED finite `X`, where `D̂_X`
is the empirical distribution of `D(x)` sampled over an interval —
distinct from the asymptotic CLT statements (Selberg, Hejhal). The
exchangeable-pair Stein construction extracts an explicit Berry-Esseen
rate from a SINGLE differential operator on test functions. Question:
does empirical `W_1(D̂_X, N(0, σ̂²))` scale as the Stein-CLT prediction
`O(1/√K)` (genuine Gaussianity) or plateau at a positive constant
(quantitative non-Gaussian arithmetic structure invisible to the
asymptotic Selberg framework)?

**Why frontier:** Selberg's CLT for `log|ζ(1/2 + it)|` is well-known,
but **the analogous QUANTITATIVE finite-x Gaussianity statement for
`π(x) - Li(x)` has never been proved**. Existing unconditional bounds
(Pintz 1980; Korevaar 2002) are pointwise-discrepancy bounds, not
Wasserstein-shape bounds against a Gaussian target. **Stein's method
has never been applied to π(x) deviations** despite being the
canonical machinery for this class of problem. The exchangeable-pair
Wasserstein bound — which the project already has the data and tools
to compute — is a single explicit operator computation. If the
empirical Wasserstein plateaus at `Ω(1)`, this is the FIRST
quantitative non-Gaussianity result for `π(x) - Li(x)`, orthogonal to
the 35+ existing pseudorandomness measures (which test mod-m, Fourier,
spectral, etc., not the *direct* central-limit shape).

**Cross-domain ingredient:** Stein's method (Stein 1986; Chen-Goldstein-
Shao 2011 *Normal Approximation by Stein's Method* Springer; Ross 2011
*Probability Surveys* 8, 210). Specifically the exchangeable-pair
Wasserstein bound `W_1(W, Z) ≤ E|1 - (1/2λ) E[(W' - W)^2 | W]| + …`
applied to the increments `W = π(x + Δ) - π(x) - Li(x+Δ) + Li(x)` for
fixed window `Δ = √x` and varying anchor `x`.

**Concrete first step:** for `x ∈ [10^6, 10^7]` at `K = 1000`
log-uniform anchors `x_k`, compute the precise `π(x_k)` (sympy /
isprime over [x_k - 1000, x_k + 1000] aggregated; or precomputed sieve
table) and `D(x_k) = (π(x_k) - Li(x_k)) log(x_k) / √(x_k)`. Compute
empirical `W_1(D, N(0, σ̂²))` via sorted-CDF integration. Compare
against (a) the Stein-CLT theoretical `O(1/√K)` rate; (b) empirical
`W_1` for a matched-variance Gaussian sample of size K; (c) empirical
`W_1` for a Bernoulli-discrepancy null at the same density. If
`W_1(D) ≫ 1/√K` and ≫ Gaussian-control, the deviation has higher-than-
Gaussian arithmetic structure — characterise via Stein's perturbation
lemma (which Stein operator perturbation explains the gap).

**Failure profile (E):** `W_1(D, N(0, σ²))` matches the `1/√K`
Gaussian baseline within sample noise → confirms quantitative
finite-x Gaussianity and adds a 38th pseudorandomness measure of the
strongest possible type (Wasserstein at finite x, not asymptotic CLT).
**B-grade.** **(I):** the empirical Stein-Chen bound matches exactly
the unconditional Korevaar 2002 quantitative bound — reduces to a
known result. **(INC):** sample size `K = 1000` anchors is too small
to discriminate Wasserstein deviations of size `< 0.05` against
sample noise (`1/√K ≈ 0.032`).

**A-grade success:** `W_1(D, N(0, σ²)) ≥ c > 0` even as `K → ∞`, AND
the gap is structurally explained by a Stein operator perturbation
that ties to a specific zeta-zero contribution (e.g., Bessel-type
oscillation). FIRST quantitative non-Gaussianity result for
`π(x) - Li(x)` at finite x, opening a new structural angle on the
GUE-sieve closure family.

**B-grade success:** explicit Stein-Chen bound `W_1 ≤ c log(x)/√x`
with constant `c < 1` — quantitatively sharper than known
unconditional bounds (Pintz 1980's `O(√x exp(-c√(log x)))` is a
discrepancy bound, not a Wasserstein-shape bound).

**Cross-domain refs:**
- Ross 2011 "Fundamentals of Stein's method" *Probability Surveys* 8, 210
  https://projecteuclid.org/euclid.ps/1314825876
- Chen-Goldstein-Shao 2011 *Normal Approximation by Stein's Method* Springer
- Barbour-Chen 2005 *An Introduction to Stein's Method* (Singapore Univ Press)
- Stein 1986 *Approximate Computation of Expectations* (IMS Lecture Notes 7)
- Wikipedia: Stein's method
  https://en.wikipedia.org/wiki/Stein%27s_method

**Budget:** 1 session. Wasserstein computation is `O(K log K)`;
`π(x)` to 10^7 is precomputable in seconds via segmented sieve.
Outcome decidable in a single session.

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

### D2 — Topological data analysis on the prime sequence  [CLOSED S96, see "Closed attacks"]

**Question:** Apply persistent homology to the prime indicator
sequence (or the gap sequence). Are there persistent topological
features at any scale that don't appear in matched random sequences?

**Why frontier:** TDA has identified non-obvious structure in
sequences considered random (RNA folding, neural data). Application
to primes is essentially unpublished.

**Cross-domain ingredient:** persistent homology (Carlsson 2009 BAMS,
ripser, Edelsbrunner-Harer 2010).

**Closure (S96, mode I, B-grade negative-shape edge E2.17):** Takens-
embedded (delay 1, dim d ∈ {2, 3, 4}) Cramér-normalised prime gaps
on a sliding window of M=2000 consecutive gaps near p ≈ 10^6, run
through ripser Vietoris-Rips persistent H_0 / H_1 with thresh=4. The
PRE-REGISTERED F3 falsifier holds (PRIMES > 3σ below BOTH IID Exp(1)
AND gap-permuted controls): T0 (total H_0 persistence) z(B1) =
−10.31, z(B2) = −8.70 (rank 0/20 in K=20 bootstrap each); T1 (total
H_1 persistence) z(B1) = −4.20, z(B2) = −11.99 (rank 0/20 each).
Robust at second window (x ≈ 5·10^6), all d ∈ {2, 3, 4}, and across
M ∈ {500..4000} with z growing super-linearly. Mechanism: HL k-tuple
admissibility constrains consecutive gaps to repeat residue patterns
more often than random, creating tighter clusters (smaller T0) and
fewer random-loop H_1 cycles (smaller T1) in delay-embedding space.
**Fifth orthogonal HL-detection category** after E2.13 (Gowers),
E2.14 (Anderson), E2.15 (alg. immunity), E2.16 (DPP failure) — first
project measurement in algebraic-topological / metric-space geometry.
NOT an algorithmic opening: VR-PH costs O(M^3); no polylog evaluator
suggested. See `experiments/topological/persistent_homology_chi_p/`,
`archive/sessions/session96_d2_persistent_homology_chi_p.md`.

**Successor challenges (proposed in S96):** (D2.a) PH of W=210
W-tricked gap sequence — does the W-trick erase the T0/T1 deficit
(would link to E2.13)? (D2.b) Persistence-image vector classifier
on PRIMES vs B1 vs B2 windows — interpretable spectrum. (D2.c)
Sliding-window {χ_P indicator} embedding instead of gap embedding,
PH on Hamming distance.

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

### D7 — Determinantal point process (DPP) fit to the integer prime sequence  [CLOSED S95, see "Closed attacks"]

**Question:** A determinantal point process is a random simple point
configuration whose `k`-point correlation function equals
`ρ_k(x_1, …, x_k) = det[K(x_i, x_j)]_{i,j ≤ k}` for some self-adjoint,
positive trace-class kernel `K`. The GUE eigenvalues are a DPP with
the sine kernel; spanning trees are DPPs; non-intersecting random walks
are DPPs. Question: empirically, is the integer prime sequence
`{p_n}` consistent with being a DPP for some translation-invariant
kernel `K(n - m)`? Specifically: do the empirical 3-point
correlations of `χ_P` match the determinantal prediction
`R_3(t_1, t_2) = det[[K(0), K(t_1), K(t_2)], [K(t_1), K(0), K(t_2-t_1)],
[K(t_2), K(t_2-t_1), K(0)]]` for `K(t)` recovered from the empirical
pair correlation? If yes, the kernel `K` is a NEW arithmetic object
distinct from {ζ-zeros, Hardy-Littlewood singular series} and gives
`π(N) = N · K(0)` as a candidate polylog formula.

**Why frontier:** the Cramér model treats primes as a Poisson process
(the trivial DPP `K = 0`), and Cramér's predictions DEVIATE from
observation (Maier 1985 disproof). Hardy-Littlewood pair correlations
indicate non-Poisson REPULSION at small distances — exactly the
qualitative signature of a DPP with NON-ZERO kernel. **No published
work fits a DPP kernel to the integer prime sequence** despite the
direct analogy to GUE eigenvalues (which DO form a DPP, sine kernel).
The frontier question: is there a translation-invariant kernel
`K(n - m)` such that the empirical `k`-point correlations of `χ_P`
match the Pfaffian/determinant identities for `k = 2, 3, 4`? If yes,
`K` is a new arithmetic invariant orthogonal to the singular series.
If no, the DEVIATION at the first failing `(k, t)` tuple is itself a
new structural fact (anti-DPP).

**Cross-domain ingredient:** Determinantal point processes (Hough-
Krishnapur-Peres-Virag 2009 textbook *Zeros of Gaussian Analytic
Functions and Determinantal Point Processes* AMS ULect 51; Soshnikov
2000 "Determinantal random point fields" Russ. Math. Surv. 55, 923).
The single technique: the **test of determinantality** via over-
determined `k`-point correlation identities — for a translation-
invariant putative kernel `K(t)` recovered from `R_2`, the Pfaffian
prediction at `k = 3, 4` is over-determined and either matches or
falsifies DPP structure.

**Concrete first step:** for `N = 10^7`, compute the empirical
pair-correlation `R_2(t) = (1/N) Σ_{n ≤ N} χ_P(n) χ_P(n+t)` for
`t ∈ [1, 200]`. Solve for the candidate kernel via
`K(t)^2 = (1/log N)^2 - R_2(t)` (DPP identity, treating primes as a
density-`1/log N` DPP). Then **test**: predict `R_3(t_1, t_2)` from
`K(t)` via the 3×3 determinantal identity above and compare to
empirical `R_3(t_1, t_2) = (1/N) Σ_n χ_P(n) χ_P(n+t_1) χ_P(n+t_2)`
for admissible `(t_1, t_2)` (e.g., `(2, 6)`, `(2, 8)`, `(6, 12)`,
`(4, 10)` — H-L admissible 3-tuples). Compute the relative error
`|R_3^{emp} - R_3^{DPP}| / R_3^{HL}` (HL = Hardy-Littlewood
prediction). If error `< 0.05` across all tested triples, primes are
empirically consistent with a DPP at this scale. If error `> 0.5` at
any single triple, primes are NOT a DPP.

**Failure profile (E):** the empirical `K(t)` reduces to a known
function of the Hardy-Littlewood pair-correlation density `S(t)`
(itself derivable from the singular series E2.13 / E2.14). The DPP
frame adds no new arithmetic content — closes as "DPP kernel of
primes IS HL singular series in a different parametrisation."
**(I):** primes are NOT a DPP — empirical `R_3` deviates from the
determinantal prediction by `> 5σ` at some `(t_1, t_2)`. **B-grade
quantitative anti-DPP result** — new structural fact about prime
correlations. **(INC):** at `N = 10^7` the third cumulant signal-to-
noise is `< 1`, so the test is inconclusive without `N ≥ 10^9`.

**A-grade success:** primes are consistent with a DPP with a kernel
`K(t)` admitting a closed form OR a polylog evaluator — gives a NEW
polylog representation of `π(N) = N · K(0) + correction`. Would be a
fundamentally new arithmetic structure on χ_P.

**B-grade success:** primes are NOT a DPP — quantitative breakdown of
the determinantal identity at a specific `(k, t)` tuple, producing a
new structural edge "primes are anti-determinantal beyond pair
correlation." Complements E2.13 / E2.14 with a *cumulant-shape* fact.

**Cross-domain refs:**
- Hough-Krishnapur-Peres-Virag 2009 *Zeros of Gaussian Analytic
  Functions and Determinantal Point Processes* (AMS ULect 51)
  https://www.ams.org/books/ulect/051/
- Soshnikov 2000 "Determinantal random point fields" Russ. Math. Surv. 55, 923
  https://arxiv.org/abs/math/0002099
- Wikipedia: Determinantal point process
  https://en.wikipedia.org/wiki/Determinantal_point_process

**Budget:** 1-2 sessions. Pair correlation already known; the new
computation is the 3-point DPP test (`O(N · t_max^2) ≈ 10^{11}` ops),
feasible overnight on one core.

### D8 — Differentiable architecture search (DARTS) for π(x) circuits

**Question:** §A.A1 (S84) used SAT/ILP to enumerate small TC⁰
circuits computing PRIMES at N=8, finding that depth-2 sign-threshold
W=1 circuits need ≥ 17 gates AND that PRIMES has a statistically
significant 1-bit-predictor advantage over matched random functions
(p < 0.001 at N=6). The SAT approach is provably intractable beyond
N ≈ 10. Differentiable Architecture Search (Liu-Simonyan-Yang 2019
ICLR) replaces discrete circuit enumeration with a continuous
relaxation: each gate is a softmax-mixture of `{AND, OR, XOR, MAJ_k,
ID}` operations, weights are trained by gradient descent on the
`(n, χ_P(n))` regression loss, then the architecture is discretised
to the highest-weight operations. DARTS reaches search spaces of size
10^4 gates that SAT cannot. Question: does DARTS, applied to a
fixed-depth threshold-circuit search space at N=12, 14, 16, find
low-loss circuits computing PRIMES that EXTRAPOLATE to held-out
N ∈ {20, 24}? If yes, the discovered circuits are CANDIDATE evidence
for PRIMES ∈ TC⁰; their structure can be inspected and potentially
formalised.

**Why frontier:** **no published work applies neural architecture
search to primality**. The S84 SAT approach is intrinsically
size-limited (exponential blow-up at N ≥ 10). DARTS is the canonical
bridge between discrete circuit synthesis and continuous optimisation,
achieving 10^4-gate searches in hours where SAT would take centuries.
Two outcomes are both informative: (i) DARTS finds a low-loss
generalising architecture → CANDIDATE TC⁰ circuit, structural
inspection could yield a proof; (ii) DARTS provably fails (loss
plateaus, no extrapolation) → empirical *circuit-search lower bound*
complementary to S84, in a much larger search space.

**Cross-domain ingredient:** Differentiable architecture search
(Liu-Simonyan-Yang 2019 "DARTS: Differentiable Architecture Search"
ICLR; Bender et al. 2018 "Understanding and Simplifying One-Shot
Architecture Search" ICML; JAX automatic differentiation Bradbury et
al. 2018). The single technique: continuous relaxation of discrete
circuit choice via softmax `α` over a fixed gate library, bilevel
optimisation (architecture α and gate weights w trained alternately),
then arg-max discretisation. Applied to a TC⁰-shaped search space
(depth `d ≤ 4`, threshold gates with bounded weights, fan-in ≤ N).

**Concrete first step:** implement a DARTS-style search in JAX over
a circuit search space of depth 3, bottom-layer fan-in `N=12`,
mid/top layer fan-in 16, with gate library
`{AND, OR, XOR, MAJ_3, MAJ_5, MAJ_7, ID}`. Train on `(n, χ_P(n))` for
`n ∈ [0, 2^{12})` for 1000 epochs (Adam, lr=10^{-3}, batch=256).
Examine: (i) final training loss vs the best feasible loss at this
size; (ii) generalisation to held-out `n ∈ [2^{12}, 2^{12} + 1000)`;
(iii) discretised circuit's behaviour at `N ∈ {14, 16, 20}` — does
the same architecture template still fit? Compare to a control task
(uniformly random Boolean function of matched density 18/64) — does
DARTS achieve LOWER loss on PRIMES than on random? If yes, this is a
circuit-complexity analogue of S84's "PRIMES has a 1-bit predictor
advantage" but at N=12 rather than N=6 and in a vastly larger search
paradigm.

**Failure profile (E):** DARTS converges to a high-loss plateau on
PRIMES (no low-depth circuit exists at this size and search budget),
closing as "DARTS-confirmed depth-3 hardness, circuit-complexity
analogue of S84 at higher N." **(I):** DARTS finds low-loss circuits
that overfit — fail to extrapolate to `N + δ`. The discovered
circuits are lookup tables in disguise (matrix factorisation of the
truth table). **(INC):** the DARTS bilevel optimisation is unstable
on PRIMES; loss is dominated by gradient noise, no signal at all.

**A-grade success:** DARTS finds a depth ≤ 4, size ≤ 10^4 gate
circuit whose discretised form computes PRIMES with > 99% accuracy
on held-out N=20 — PRIMES ∈ TC⁰ strongly suggested empirically;
structural inspection could yield a proof.

**B-grade success:** DARTS achieves provably-lower loss on PRIMES
than on density-matched random Boolean functions (statistically
significant at p < 0.001 over ≥ 30 random seeds) — extends the S84
empirical PRIMES-vs-random circuit-complexity gap to a NEW circuit-
search method at higher N. Complements E5.x / S84's depth-2
threshold finding with a depth-3+ DARTS finding.

**Cross-domain refs:**
- Liu-Simonyan-Yang 2019 "DARTS: Differentiable Architecture Search"
  ICLR 2019 https://arxiv.org/abs/1806.09055
- Bender et al. 2018 "Understanding and Simplifying One-Shot
  Architecture Search" ICML 2018
- Baydin et al. 2018 "Automatic Differentiation in Machine Learning:
  a Survey" JMLR 18, 1-43
- Wikipedia: Differentiable programming
  https://en.wikipedia.org/wiki/Differentiable_programming

**Budget:** 2 sessions. Session 1: implement DARTS search at N=12,
depth 3, train 1000 epochs, evaluate generalisation. Session 2:
scale to N=16, ablate gate library, attempt extrapolation tests at
N ∈ {20, 24}.

### D9 — Sum-product gain on the prime set in F_p (Bourgain-Glibichuk-Konyagin)

**Question:** for `A_prime = {p_i mod p : p_i prime, p_i < p} ⊂ F_p`,
the Bourgain-Glibichuk-Konyagin (2006) sum-product theorem gives
`max(|A+A|, |A·A|) ≥ |A|^{15/14} / (log|A|)^{2/7}` for any `A` with
`|A| ≤ p^{7/13}(log p)^{-4/13}`. The bound is for ARBITRARY `A`. The
question: is the *empirical* sum-product gain
`g(A) := max(|A+A|, |A·A|) / |A|^{15/14}` strictly LARGER for the prime
set than for matched-cardinality random subsets of `F_p`? If yes,
primes carry an additive-multiplicative joint structure invisible to
the existing 35+ pseudorandomness measures (all of which are local /
1-dimensional). If matched, primes are BGK-saturating in the strongest
combinatorial sense.

**Why frontier:** sum-product is THE canonical test of "no joint
structure under both addition and multiplication." Primes ARE the
multiplicative atoms but their additive structure is constrained only
by Hardy-Littlewood singular series. **The prime set in F_p has never
been directly tested as a sum-product problem**; published BGK
applications all target generic subsets, multiplicative subgroups, or
random-set baselines (Garaev 2008 is closest, but tests *random*
A ⊂ F_p, not A_prime). If `g(prime) ≫ g(random)`, the gain has to
come from the HL mod-q residue distribution → links additive
combinatorics to the singular series at the *joint* level (orthogonal
to E2.13 Gowers, which is purely additive). If `g(prime) ≈ g(random)`,
this is the FIRST sum-product-class pseudorandomness measure of χ_P,
and the strongest such test (joint-structure invariant).

**Cross-domain ingredient:** Bourgain-Glibichuk-Konyagin sum-product
theorem in F_p (Bourgain-Konyagin 2003 + Bourgain-Glibichuk-Konyagin
2006), Tao-Vu 2006 *Additive Combinatorics* textbook, additive
energy / multiplicative energy machinery. UNUSED in the project (§7
"Sum-product theorems" row).

**Concrete first step:** for `p ∈ {1009, 10007, 100003, 1000003}`,
let `A_prime := {q mod p : q prime, q < p}` (cardinality `~ π(p)`).
Direct computation: `|A_prime + A_prime|` (mod p), `|A_prime · A_prime|`
(mod p), each via Python set union over `O(|A|^2)` pairs (`p ≤ 10^6`
gives `|A|^2 ≤ 6×10^9`, feasible on one core overnight with memory-
optimised hashing). Compute the gain `g(A_prime) =
max(|A+A|, |A·A|) / |A_prime|^{15/14}`. Compare against:
(a) 100 random subsets `B` of `F_p` of matched cardinality
`|B| = |A_prime|`, drawn uniformly without replacement;
(b) the BGK theoretical lower bound `|A|^{15/14}/(log|A|)^{2/7}`;
(c) the *multiplicative subgroup* baseline (e.g., quadratic residues)
of matched cardinality. Bootstrap z-score:
`z = (g(prime) - mean(g_random)) / std(g_random)`. Is `|z| > 3` at
all 4 primes?

**Failure profile (E):** `g(prime) = g(random)` within sample noise
— closes as 38th pseudorandomness measure (sum-product gain). **B-grade**.
**(I):** `g(prime) < g(random)` — primes' singular-series structure
*reduces* sum-product gain (they're more compressible than random
under the joint operation). Would be a NEW structural fact and a
B-grade negative-shape edge "primes are sub-BGK in joint additive-
multiplicative gain." **(INC):** for `p ≥ 10^7`, direct enumeration
of |A+A|, |A·A| is `O(|A|^2) ≈ 10^{11}` ops — borderline; needs
GPU acceleration or sketch-based estimation.

**A-grade success:** `g(prime) > g(random)` by a factor that GROWS
with `p`, AND the gain has a closed-form interpretation via HL
mod-q residue distribution (specifically: the singular series
predicts the under-representation of certain `(a+b) mod p` residues
when both `a, b` are prime, giving a `|A+A|` reduction). Would be
the FIRST sum-product result for the prime set and would link
additive combinatorics directly to HL.

**B-grade success:** clean negative result — `g(prime) - g(random)`
is `o(1)` as `p → ∞` (joint additive-multiplicative pseudorandomness
of primes). Adds a 38th pseudorandomness measure of a structurally
distinct type.

**Cross-domain refs:**
- Bourgain-Glibichuk-Konyagin 2006 "Estimates for the number of sums
  and products and for exponential sums in fields of prime order"
  J. London Math. Soc.
  https://londmathsoc.onlinelibrary.wiley.com/doi/abs/10.1112/S0024610706022721
- Garaev 2008 "An explicit sum-product estimate in F_p"
  https://arxiv.org/abs/math/0702780
- Iosevich 2009 "Sum-product phenomena in F_p: a brief introduction"
  https://arxiv.org/abs/0904.2075
- Tao-Vu 2006 *Additive Combinatorics* Cambridge Univ Press

**Budget:** 1-2 sessions. `p ≤ 10^6` direct computation tractable in
one session.

### D10 — Mahler measure of the prime generating polynomial f_N(z) = Σ_{n≤N} χ_P(n) z^n

**Question:** define the polynomial
`f_N(z) := Σ_{n=1}^{N} χ_P(n) · z^n` with 0/1 coefficients. Its Mahler
measure
`m(f_N) = exp(∫_0^1 log|f_N(e^{2π i θ})| dθ) = |a_d| · ∏_i max(1, |α_i|)`
(Jensen / Boyd) is a multiplicative height invariant tied to
Lehmer's conjecture and the logarithmic Weil height. Question: does
`m(f_N)` deviate from random matched-density 0/1 polynomials? Three
candidate scaling regimes: (i) `m(f_N) = O((log N)^c)` (near-cyclotomic
factorisation, *full algebraic compressibility*); (ii) `m(f_N) = N^α`
for `α ∈ (0, 1)` distinguishable from random-matched-density baseline
(intermediate height structure); (iii) `m(f_N) ~ √N` Lehmer-typical
for non-cyclotomic random 0/1 polynomials (no compression).

**Why frontier:** Mahler measure connects to logarithmic Weil heights,
Boyd's range-of-Mahler's-measure conjecture, Lehmer's conjecture on
the smallest non-cyclotomic measure, and entropy-rate of polynomial
spectral content. **No published work computes `m(f_N)` for `f_N` the
prime generating polynomial.** Cyclotomic polynomials have `m = 1`
(Kronecker's theorem); if `f_N(z)` factors near-cyclotomically over
`Q[z]`, primes are radically compressible and π(N) is computable in
polylog via roots-of-unity sampling. Smyth 2008 (height survey)
explicitly notes that Mahler measures of arithmetic-indicator
polynomials are an unsolved height question.

**Cross-domain ingredient:** Mahler measure / Lehmer conjecture /
Weil height. Boyd 1981 / Smyth 2008 survey. UNUSED in the project
(§2 / §10 algebraic-height — proposed entry).

**Concrete first step:** for `N ∈ {2^{10}, 2^{12}, 2^{14}, 2^{16}}`,
build the coefficient vector `(χ_P(n))_{n=1}^{N}` (length `N`,
support size π(N)). Estimate `m(f_N)` via Jensen-FFT integral:
`log m(f_N) ≈ (1/M) Σ_{k=0}^{M-1} log |f_N(e^{2π i k/M})|`
with `M = 2^{18}` sample points (FFT cost `O(M log M)`).
Numerical-precision care: `f_N` evaluated via mpmath at 50 dps.
For small `N ≤ 2^{10}` also factor `f_N(z)` over `Q[z]` via
sympy and inspect cyclotomic factors (count of cyclotomic
factors and their indices). Compare `m(f_N)` to:
(a) ensemble of 100 random 0/1 polynomials at matched density
`π(N)/N`; (b) cyclotomic-only baseline (any cyclotomic factors
of `f_N` go to `m = 1`); (c) HL-singular-series-weighted random
polynomial (correlations match HL pair correlation; this is the
**most stringent random baseline**). Plot `log m(f_N) / log N` vs
`log N`. Is the slope `1/2` (random Lehmer-typical), or strictly
less?

**Failure profile (E):** `m(f_N) ~ N^{1/2}` matching density-matched
random — primes are Lehmer-typical, no algebraic-height structure
beyond density. Closes as 38th pseudorandomness measure. **B-grade.**
**(I):** Jensen-FFT precision insufficient at `N ≥ 2^{14}` to
discriminate sub-linear from linear scaling — addresses with mpmath.
**(INC):** `f_N` has very few cyclotomic factors (since the 0/1
indicator is irreducible-typical) — but a *partial* cyclotomic share
matters for `m`.

**A-grade success:** `m(f_N) = O((log N)^c)` for some `c > 0` —
near-cyclotomic factorisation gives an immediate polylog evaluator
for π(N) via roots-of-unity sampling at the cyclotomic factors.
Would be the FIRST polylog representation of χ_P.

**B-grade success:** `m(f_N) = N^α` with `α` statistically
distinguishable from random-baseline `α_random ≈ 1/2` (deviation
> 3σ across 100 controls). FIRST project measure detecting
*algebraic-height* structure of χ_P; complements the existing
35-measure pseudorandomness battery with a new mathematical category
(transcendence / algebraic height).

**Cross-domain refs:**
- Smyth 2008 "The Mahler measure of algebraic numbers: a survey"
  in *Number Theory and Polynomials* CUP
  https://homepages.ed.ac.uk/cjsmyth/papers/107.pdf
- Boyd 1981 "Speculations concerning the range of Mahler's measure"
  Canad. Math. Bull. 24
- Lehmer 1933 "Factorization of certain cyclotomic functions" Ann.
  Math. 34
- Wikipedia: Mahler measure
  https://en.wikipedia.org/wiki/Mahler_measure

**Budget:** 1 session. Jensen-FFT at `N=2^{16}`, `M=2^{18}` is
< 1 minute on one core; the 100-control ensemble is ~1 hour.

### D11 — Shadow tomography sample complexity for π(x) extraction

**Question:** Aaronson 2018 ("Shadow Tomography of Quantum States",
STOC 2018, arXiv:1711.01053) proved that for an n-qubit state `ρ`,
`M` observables can be estimated to ε accuracy using
`Õ(ε^{-4} log^4 M · log D)` copies — exponentially better than naive.
The CLASSICAL analogue (Huang-Kueng-Preskill 2020 "Predicting many
properties of a quantum system from very few measurements" Nature
Physics 16, 1050) gives `poly(log M)` for local observables.
Question: under the *random-shadow query model* — where each query is
`Σ_n σ(n) χ_P(n)` for a random Rademacher mask `σ : [1, N] → {-1, +1}`
— how many queries suffice to estimate `π(M)` for ALL `M ∈ [1, N]`
simultaneously to within ±1? If `K = poly(log N)`, this is a polylog
QUERY complexity for the full π profile (orthogonal to direct
evaluation; new model of computation).

**Why frontier:** the project has E1.5 (Shannon entropy lower bound:
`H(π(X) mod m) ≥ log m`) as an *information-theoretic* bound, and
explicit-formula closures as *time-complexity* bounds. **Shadow
tomography gives a SAMPLE / QUERY complexity bound** — a third
independent computational complexity axis the project has never
explored. If the random-shadow query model gives `K = poly(log N)`
for full-profile estimation, this is a polylog *query* complexity for
π — not polylog evaluation but polylog ORACLE complexity, a direct
complement to the project's information-computation gap framework
(F3 paper). If the random-shadow lower bound is `K = Ω(N^c)` for any
unbiased estimator, this strengthens the project's existing bounds
with a NEW model. Either outcome is a publishable result orthogonal
to existing closures.

**Cross-domain ingredient:** Shadow tomography (Aaronson 2018 STOC),
classical shadows protocol (Huang-Kueng-Preskill 2020 Nature Physics
16, 1050 = arXiv:2002.08953), random projections / Johnson-
Lindenstrauss. UNUSED in §8.

**Concrete first step:** for `N = 2^{20}`, treat
`x = (χ_P(n))_{n=1}^{N} ∈ {0, 1}^N` as the unknown signal. Define
the family of *cumulative-window* observables
`O_M := (1/M) · 1_{[1, M]}^T x = π(M)/M`, for `M ∈ {2^k : k ≤ 20}`.
Sample `K ∈ {100, 1000, 10000}` random Rademacher masks `σ_j : [1, N]
→ {-1, +1}` i.i.d.; for each, query the linear functional
`y_j := Σ_n σ_j(n) χ_P(n)` (one "shadow"). Use the Huang-Kueng-
Preskill linear unbiased estimator
`π̂(M) = N · Σ_j y_j · ⟨σ_j, 1_{[1, M]}⟩ / (K · M · N)` (rescaling the
random-Rademacher-shadow inversion formula). Compute the empirical
L^∞ error `max_M |π̂(M) - π(M)|`. Plot error vs `K`. Is the scaling
`K^{-1/2}` (CLT-like, no special structure) or sub-`K^{-1/2}` (better
than naive sample-mean)? Compare to direct-sieve baseline `O(N)` ops.

**Failure profile (E):** classical-shadow scaling on cumulative-window
observables is exactly equivalent to direct sample-mean estimation of
`π(M)` — `K = Ω(M)` for ε = 1 — no gain over sieve. Closes as
"shadow protocol on global observables = sample mean." **(I):**
shadow protocol gives unbiased estimators with constant sample
complexity per fixed M, but `K` blows up with the number of
simultaneously-estimated M's via union bound — `K = Ω(log(log N))`
trivially, no polylog tightness. **(INC):** the random-Rademacher
shadow protocol is suboptimal for global indicators; need a SMARTER
projection ensemble (e.g., random tensor-product Walsh-Hadamard).

**A-grade success:** classical-shadow protocol with `K = O(log^c N)`
queries to χ_P estimates `π(M)` for ALL `M ≤ N` to ε = 1
simultaneously. POLYLOG QUERY complexity for π — new computational-
model result; would be a *publishable* result on its own.

**B-grade success:** explicit query lower bound `K ≥ Ω(N^β)` for any
unbiased estimator under the random-shadow query model. Strengthens
the project's information-theoretic bounds with a new query-
complexity bound complementary to E1.5 / F3.

**Cross-domain refs:**
- Aaronson 2018 "Shadow Tomography of Quantum States" STOC 2018
  https://arxiv.org/abs/1711.01053
- Huang-Kueng-Preskill 2020 "Predicting many properties of a
  quantum system from very few measurements" Nature Physics 16, 1050
  https://arxiv.org/abs/2002.08953
- Wikipedia: Quantum tomography (Shadow tomography section)
  https://en.wikipedia.org/wiki/Quantum_tomography

**Budget:** 1-2 sessions. Random-shadow simulation at `N = 2^{20}`,
`K = 10^4` is `O(KN) = 10^{10}` ops — borderline tractable on one
core; GPU optional.

### D12 — Compressed sensing of χ_P in structured arithmetic-progression dictionaries

**Question:** compressed sensing (Candes-Tao 2006) tells us a signal
`x ∈ R^N` that is `K`-sparse in dictionary `D` with restricted
isometry property (RIP) is recoverable from `O(K log(N/K))` random
linear measurements via L1-minimization. Question: is `χ_P ∈ {0,1}^N`
`K`-sparse in a STRUCTURED arithmetic-progression dictionary
`D = {1_{aZ + b}|_{[1, N]} : a ∈ S, 0 ≤ b < a}` for some scale set
`S ⊂ [1, √N]`, with `K` polylog in `N`? If yes, then π(M) for `M ≤ N`
is computable in `O(K log N)` operations from `K` dictionary atoms —
a polylog evaluator at finite `N`. The Cramér model and HL singular
series both predict `K = Θ(π(N))` (full-rank, no sparsity). The
empirical question: at what `N` does the OMP-greedy decomposition
saturate, and is `K` strictly sub-linear?

**Why frontier:** compressed sensing has NEVER been applied to χ_P
in a STRUCTURED arithmetic dictionary. Naive Fourier dictionary
tested previously (CLOSED_PATHS line 47: "82% of Fourier modes
needed") — but the AP dictionary is structurally distinct: its
atoms encode CRT residue information, the natural language for
primes via HL singular series. The empirical question: does
OMP-greedy give `K = poly(log N)` for ε ≤ 0.1 reconstruction error
in the AP dictionary, or does it match `K = Θ(π(N))` (Cramér-typical)?

**Cross-domain ingredient:** compressed sensing / RIP / OMP / L1-min
recovery (Candes-Tao 2006; Foucart-Rauhut 2013 textbook). UNUSED in
§6 ("Compressed sensing on χ_P" row).

**Concrete first step:** at `N = 2^{14}`, build the AP dictionary
`D = {δ_{aZ + b}|_{[0, N)} : a ∈ {2, 3, 5, 7, 11, 13, 30, 210, 2310,
30030}, 0 ≤ b < a}` (≈ 33000 atoms; HL-relevant sieve moduli).
Solve OMP greedy decomposition: at each step, pick the atom
`d ∈ D` maximising `|⟨d, r⟩| / ‖d‖` where `r` is the current
residual, then update `r ← r - ⟨d, r⟩ d / ‖d‖²`. Stop when
`‖r‖²/‖χ_P‖² < ε` for `ε ∈ {0.5, 0.2, 0.1, 0.05}`. Plot `K_ε`
(sparsity needed) vs `ε` and vs `log N`. Compare to:
(a) ensemble of 100 random Bernoulli vectors at matched density
`π(N)/N`; (b) Liouville indicator `1[λ(n) = -1]`; (c) Walsh-Hadamard
dictionary baseline (orthonormal). Does `K_ε = O((log N)^c)` for any
ε, or `K_ε = Θ(π(N))`?

**Failure profile (E):** `K_ε = Θ(π(N))` for ε ≤ 0.1 — AP dictionary
cannot approximate χ_P with poly-log sparsity. Reduces to Cramér-
model duplicate. Closes as "AP-dictionary L1-min has full effective
rank for χ_P". **B-grade if quantitative.**
**(I):** AP dictionary highly coherent (overlapping atoms with
shared support) — OMP unstable, but L1-min relaxation may give
different `K_ε`. Test both.
**(INC):** at small `N` finite-size noise dominates — need
`N ≥ 2^{18}` for clean signal; OMP at `N = 2^{18}` with `~10^5`
atoms is `O(N · |D| · K)` ≈ `10^{11}` ops, borderline.

**A-grade success:** `K_ε = poly(log N)` for some structured
dictionary `D` (AP or HL-weighted-AP or W-trick-shifted-AP from
Green-Tao). Provides a polylog REPRESENTATION of χ_P at finite `N`,
hence a polylog π(M) evaluator for M ≤ N. POLYLOG π(x) at the
project's target scale.

**B-grade success:** prove rigorously `K_ε ≥ Ω(π(N))` for the AP
dictionary family — closing this dictionary as not-compressible
and adding a new EDGE.

**Cross-domain refs:**
- Candes-Tao 2006 "Decoding by linear programming" IEEE Trans. Inf.
  Theory 51(12)
- Foucart-Rauhut 2013 *A Mathematical Introduction to Compressive
  Sensing* Birkhäuser
- Wikipedia: Restricted isometry property
  https://en.wikipedia.org/wiki/Restricted_isometry_property
- Wikipedia: Compressed sensing
  https://en.wikipedia.org/wiki/Compressed_sensing

**Budget:** 1-2 sessions. `N = 2^{14}` OMP feasible in minutes;
scale to `2^{18}` if signal observed.

### D13 — Subword complexity / topological entropy of χ_P as binary infinite word  [CLOSED S104, see "Closed attacks"]

**Closure (S104, mode E, B-grade negative-shape edge E2.19):** at
N = 5·10^6, n_max = 22, K = 20, on five chi_P-derived streams (RAW /
ODD / W6 / W30 / W210), the cascade max |z_perm| in matched-density
Bernoulli + permutation baselines is RAW 132.7 → ODD 277.1 → W6
120.5 → W30 24.8 → W210 8.4 — clean monotone reduction by ~1.5
orders of magnitude across W ∈ {1, 6, 30, 210}. Pre-registered F3
(PRIMES > 3σ from BOTH B1 AND B2 at every n in [n_lo, n_hi] with
n_hi ≥ 18) HOLDS for ODD/W6/W30; W=210 partially erases (residual
8.4σ at n=12, sign-flips at n=17). Effective topological entropy
`h_eff(n) = log_2 p_w(n) / n` of W=210 stream matches Bernoulli
matched-density to ≤ 0.001 across n ∈ [1, 22]. **Mechanism:** subword
complexity counts distinct prime-position configurations in sliding
length-n windows; chi_P configurations are constrained by mod-p
admissibility for p ≤ n; W = primorial(p_k) sieves out those
constraints. Same HL k-tuple-admissibility engine that drives
E2.13/E2.14/E2.16/E2.17. **Adds EDGE E2.19**; **38th pseudorandomness
measure**; **SIXTH orthogonal HL-detection family** (first in
symbolic-dynamics / factor-complexity category). NOT an algorithmic
opening: `p_w(n)` evaluation is `O(L log L)`. See
`experiments/dynamical/subword_complexity_chi_p/`,
`archive/sessions/session104_d13_subword_complexity.md`.

**Successor challenges (proposed in S104):** (D13.a) Scale to W=2310
with N ≥ 5·10^7 to test whether the W=210 residual (8.4σ at n=12)
collapses to ≤ 3σ or persists — single session. (D13.b) Subword
complexity of `lambda(n)` binarised as `1[lambda = -1]` — predicted
no W-trick needed (cf. E2.18); extends multiplicative-regime
pattern. (D13.c) Joint p_w distribution across W ∈ {6, 30, 210,
2310} as a single Mahalanobis statistic.

**Question (original):** the prime indicator infinite word
`w = χ_P(2) χ_P(3) χ_P(4) ... ∈ {0, 1}^∞` has SUBWORD COMPLEXITY
`p_w(n) := #{distinct length-n factors of w}`. For Bernoulli words at
density `ρ`, `p_w(n) ~ 2^n` and topological entropy `h_w = log 2`.
For Sturmian words, `p_w(n) = n + 1` and `h_w = 0`. For automatic
sequences (Thue-Morse), `p_w(n) = O(n)`. Question: what is `p_w(n)`
for `w = χ_P` at `n ∈ [1, 30]` empirically? Specifically, is `h_w`
strictly less than `log 2` (sub-Bernoulli entropy → genuine arithmetic
structure), or is `p_w(n) / 2^n → 1` (Bernoulli-typical)?

**Why frontier:** subword complexity is THE canonical combinatorial
invariant of an infinite word (Lind-Marcus 1995 textbook *Symbolic
Dynamics and Coding*; Cassaigne-Nicolas 2010 "Factor complexity"
survey). The Morse-Hedlund theorem (1938-40) classifies sequences by
`p_w(n)` growth: ultimately periodic ⇔ `p_w(n) ≤ n` for some `n`;
aperiodic ⇒ `p_w(n) ≥ n+1`. **Subword complexity of χ_P has not been
computed in the published literature** — closest known result is for
the Liouville characteristic word (which is *normal*). The FRACTRAN
closure (CLOSED_PATHS line 179) is about a different object (FRACTRAN
as a Turing-machine encoding of primes); subword complexity of χ_P
itself is fresh in topological / symbolic dynamics. Either outcome is
informative: (a) `p_w(n) = 2^n - o(2^n)` matching Bernoulli — 38th
pseudorandomness measure in a NEW category; (b) `p_w(n) ≤ 2^n · (1 -
δ)` with `δ > 0` measurable at `n ≤ 30` — FIRST topological-dynamics
deviation of χ_P from random.

**Cross-domain ingredient:** subword complexity / topological entropy
of binary words. Lind-Marcus 1995; Cassaigne-Nicolas 2010 in
*Combinatorics, Automata and Number Theory* (Cambridge); Morse-
Hedlund 1938-40. UNUSED in §5 ("Symbolic dynamics on prime sequences"
row).

**Concrete first step:** for `N = 10^7`, build the binary string
`w_N = (χ_P(k))_{k=2}^{N}` (length `N-1`). For `n ∈ {1, 2, ..., 30}`,
count distinct length-`n` factors of `w_N` via rolling-hash sliding
window: maintain a Python set of length-`n` integer-encodings (each
window encoded as `Σ_{j=0}^{n-1} w_N[i+j] · 2^j ≤ 2^{30}`); insert at
each window position; final set size is `p_{w_N}(n)`. Memory cost
`O(N · 4 bytes) = 40 MB`; time cost `O(N · n) = 3×10^8` ops per `n`.
Compare to:
(a) 100 IID Bernoulli words at matched density `π(N)/N ≈ 0.0577`:
theoretical `p_Bern(n) ≈ Σ_k C(n,k) (ρ)^k(1-ρ)^{n-k} · (1 -
exp(-(N-n) · prob))` — saturates at `p(n) ≤ 2^n` and at `p(n) ≤ N - n`;
(b) HL-singular-series-corrected Bernoulli (n-tuple-density
prediction); (c) Liouville characteristic word (a normal word —
known `h = log 2`). Plot `log p_{w_N}(n)` vs `n`. Compute empirical
slope `h_emp = (log p(30) - log p(20))/10` and compare to `log 2`.
Is `h_emp` statistically distinguishable from `log 2` at `> 3σ`?

**Failure profile (E):** `p_{w_N}(n)` matches Bernoulli within sample
variance — closes as 38th pseudorandomness measure (subword
complexity), NEW category (topological dynamics / symbolic dynamics).
**B-grade.**
**(I):** for `n` sufficiently large (`n ≥ log_2 N ≈ 23`), all
length-`n` windows are distinct (since `N - n < 2^n` for the
relevant range), so `p_w(n) = N - n + 1` is saturated TRIVIALLY
and carries no signal. Real signal lives in `n ∈ [10, log_2 N - 2]`.
**(INC):** finite-`N` corrections at `n ≈ log_2 N` make the empirical
entropy estimate biased; need careful Lempel-Ziv-style entropy
estimator with finite-size correction.

**A-grade success:** `p_{w_N}(n)` provably bounded by a sub-
exponential function (e.g., `2^n · exp(-c√n)` for explicit `c > 0`)
AND the rate matches an HL-prediction-corrected formula — gives an
EXPLICIT factor-complexity formula for χ_P, linking symbolic
dynamics directly to HL singular series. Would be a publishable
combinatorial-on-words result.

**B-grade success:** `p_{w_N}(n)` at `n ∈ [10, 25]` statistically
distinguishable from Bernoulli matched-density (deviation > 3σ
across 100 controls). FIRST topological-dynamics deviation of
χ_P from random; complements existing 35+ pseudorandomness measures
with a new mathematical category.

**Cross-domain refs:**
- Lind-Marcus 1995 *An Introduction to Symbolic Dynamics and Coding*
  Cambridge Univ Press
- Cassaigne-Nicolas 2010 "Factor complexity" in *Combinatorics,
  Automata and Number Theory* Cambridge Encyclopedia of Math 135
- Morse-Hedlund 1938 "Symbolic dynamics" Amer. J. Math. 60
- Wikipedia: Complexity function
  https://en.wikipedia.org/wiki/Complexity_function

**Budget:** 1 session. Sliding-window factor counting at `N=10^7`,
`n ≤ 30` is `O(N · 30) ≈ 3×10^8` ops in optimised Python — ~5
minutes. 100-control ensemble: ~10 hours; or use bootstrap
resampling.

### D14 — Cellular sheaf cohomology of χ_P over the divisibility poset

**Question:** Curry's 2014 thesis *Sheaves, Cosheaves and Applications*
(Penn State; arXiv:1303.3255) develops **cellular sheaves** as finite,
algorithmically-tractable sheaves of vector spaces parametrised over
a cell complex (or, equivalently, over a poset via the order complex).
For the divisibility poset `(N_{≤N}, |)` — finite truncation of
`(N, |)` — define the cellular sheaf `F` with stalks
`F(n) := F_2 · χ_P(n)` (1-dim if `n` is prime, 0-dim otherwise) and
restriction maps `F(n) → F(m)` for `n | m` given by 0 (since composites
"lose" the prime stalk). Question: what is the cellular cohomology
`H^k(N_{≤N}, F)` for `k = 0, 1`? Does it carry HL-singular-series
structure (closure path), is it trivial (38th pseudorandomness measure
in a NEW topological-combinatorial category), or does it admit a
TC^0-recoverable basis (A-grade — gives polylog primality witness via
cohomology class evaluation)?

**Why frontier:** **no published work computes cellular sheaf
cohomology of an arithmetic indicator over the divisibility poset.**
The persistent-homology closure (S96, E2.17) used Vietoris-Rips
filtration on a *metric* delay-embedding — a metric-space topological
invariant. Cellular sheaf cohomology is **structurally orthogonal**:
it uses the *order/divisibility* structure of integers, not metric
distances. Curry's framework gives:
(a) explicit cellular cochain complex `C^k = ⊕_σ F(σ)` over k-dim cells
    with combinatorial differential;
(b) computable `H^k = ker(d^k)/im(d^{k-1})` via Smith normal form
    over `F_2` (polynomial-time in `|N_{≤N}|`);
(c) a Verdier dual `D F` whose dimensions track support-cohomology
    of χ_P — distinct from H^k(F).
This is a **categorical / order-based topological probe** entirely
absent from the project's pseudorandomness battery. After PH closure
(metric-topological), cellular sheaf cohomology is the natural
successor in topological methods, flagged in CROSS_DOMAIN_TECHNIQUES.md
as an UNUSED orthogonal route (§4 row, marked "candidate").

**Distinction from line 202 closure:** `CLOSED_PATHS.md` line 202
("Sheaf cohomology — FAIL E — Same barrier") closes sheaves on
`Spec(Z)` (i.e., algebraic-geometric sheaves; line 198 = topos theory
of Spec(Z) recovers Euler product). D14 uses **cellular sheaves on
the discrete divisibility poset** (Curry 2014 framework, finite
combinatorial vector-space sheaves), which is a categorically
different object from Spec(Z) sheaves — Curry's framework explicitly
arose from PH/TDA, not from arithmetic geometry. The two objects
share a name only.

**Cross-domain ingredient:** Cellular sheaves on posets (Curry 2014
thesis https://arxiv.org/abs/1303.3255), Hansen-Ghrist 2019 *J. Appl.
Comput. Topol.* 3 spectral sheaf theory, Robinson 2014 *Topological
Signal Processing* (Springer). The single technique: build the
order-complex of the divisibility poset truncated to `[1, N]`, attach
a stalk `F(n) = F_2^{χ_P(n)}` to each vertex, set restriction maps
`r_{n,m} = 0` for `n|m, n ≠ m, n composite`, and `r_{1,p} = id` for
`p` prime. Compute `H^0` (global sections — should pick up isolated
prime "tops" without composite covers) and `H^1` (1-cocycles /
1-coboundaries on prime-pair edges).

**Concrete first step:** for `N ∈ {64, 256, 1024, 4096}`, build the
divisibility poset via `pos = {(n, m) : n | m, n ≤ m ≤ N}`, the order
complex `K = Hasse diagram`, the cellular sheaf with the stalk
assignment above. Compute `dim H^0(K, F)` and `dim H^1(K, F)` over
`F_2` via Smith normal form (sympy `Matrix.rref` over `F_2` or
explicit boundary-matrix kernel/image). Compare to:
(a) cellular sheaf with stalks `F'(n) = F_2 · Bernoulli(π(N)/N)`
    (matched-density random) over 100 random seeds;
(b) `F''(n) = F_2 · 1[lambda(n) = -1]` (Liouville-indicator stalks);
(c) `F'''(n) = F_2 · χ_{[N/2, N]}(n)` (top-half indicator — known
    trivially structured).
Plot `dim H^k(F) / dim H^k(F')` as a function of N. Does the prime
sheaf cohomology grow at the same rate as Bernoulli (closure: trivial
combinatorial topology) or strictly slower/faster (deviation)?

**Failure profile (E):** `dim H^k(F)` matches Bernoulli matched-density
within sample noise — closes as 38th pseudorandomness measure
(cellular-sheaf cohomology), NEW mathematical category (categorical
topology / cellular cohomology) orthogonal to the metric-topological
PH closure. **B-grade.** **(I):** cellular sheaf cohomology of any
0-stalk-on-composite assignment trivialises (because composites all
have 0-stalks, the differential collapses) — closes as "the chosen
sheaf structure is too rigid to detect arithmetic." Mitigation: use
`F(n) = F_2^{Ω(n)}` (with stalk dimension `Ω(n)`, the prime-factor
count) — richer; orthogonal sub-attack. **(INC):** Smith normal form
over `F_2` at `N = 4096` is `O(|cells|^3)` — `|cells|` grows like
`Σ_{n ≤ N} d(n) ≈ N log N`, so cost is `O((N log N)^3)`; tractable
to N=4096 (`~10^{12}` ops, ~1 day) but not 65536.

**A-grade success:** `H^1(K, F)` is non-trivial AND its dimension has
a closed-form expression in `π(N)` AND a TC^0-recoverable basis (i.e.,
each cohomology class has a polylog-evaluable representative cocycle).
Then *primality of n* is decidable via cohomology class evaluation:
"n is prime iff the dual basis vector e^*_n has non-zero pairing with
some explicit cohomology class." Would be PRIMES ∈ TC^0 unconditionally,
plus a structurally NEW characterisation of primality (categorical
rather than analytic).

**B-grade success:** explicit dimension formula for `dim H^k(K, F)` in
terms of arithmetic functions of N (e.g., `dim H^1 = π(N) − ω(N) + …`)
that distinguishes the prime sheaf from Bernoulli matched-density at
> 3σ. FIRST cellular-sheaf invariant of χ_P; complements PH closure
(E2.17) with a categorical-topological measure.

**Cross-domain refs:**
- Curry 2014 *Sheaves, Cosheaves and Applications* PhD thesis Penn State
  https://arxiv.org/abs/1303.3255
- Hansen-Ghrist 2019 "Toward a spectral theory of cellular sheaves"
  J. Appl. Comput. Topol. 3
- Robinson 2014 *Topological Signal Processing* Springer
- Wikipedia: Sheaf cohomology
  https://en.wikipedia.org/wiki/Sheaf_cohomology

**Budget:** 1-2 sessions. Session 1: implement cellular sheaf and
boundary matrices at N = 256, 1024; compare to Bernoulli baseline.
Session 2 (only if signal observed at session 1): scale to N = 4096
and try richer stalk assignments (`Ω(n)`-dim, multiplicative
characters).

### D15 — Bourgain-Demeter decoupling test of the prime exponential sum

**Question:** the Bourgain-Demeter-Guth (2015) `l^2` decoupling
theorem proves: for any partition of a curved hypersurface (moment
curve, paraboloid, sphere) into caps `θ_j` of size `δ`, and any
function `f` Fourier-supported on the hypersurface, we have
`||f||_{L^p(B_R)} ≤ K_{p, ε}(δ) · (Σ_j ||f_θ_j||_{L^p(B_R)}^2)^{1/2}`
with `K_{p, ε}(δ) = δ^{-ε}` for `p ≤ 2(d+2)/d` (sharp). The decoupling
constant `K(δ)` quantifies "how independent" the cap-supported pieces
are. Question: define the prime exponential sum
`f_N(x) := Σ_{p prime, p ≤ N} e^{2π i p x / N}` (a finite sum of
exponentials supported on the prime set inside [1, N]). Empirically:
**does f_N saturate the BDG decoupling inequality** — i.e., is the
empirical decoupling constant `K_emp(N) := ||f_N||_{L^p}/(Σ_j ||f_N
restricted to cap_j||_{L^p}^2)^{1/2}` consistent with the BDG bound,
or does it deviate? A statistically significant deviation would be
a NEW arithmetic invariant of χ_P at the *Fourier-restriction* level.

**Why frontier:** **no published work directly tests χ_P's saturation
of the BDG decoupling inequality**, despite it being THE canonical
modern measure of "joint cap-decoupling structure" of a discrete
set. Bourgain-Demeter-Guth 2015 used decoupling to resolve Vinogradov's
mean value conjecture (a statement about EXTREMAL sums of e(α p^k));
the inequality is sharp on extremal configurations. The prime set is
density-`1/log N`, with mod-q residue concentration (HL singular
series). If primes saturate the decoupling inequality with constant
strictly *better* than worst-case `δ^{-ε}` — say `(log N)^c` — that
gives a **tighter Fourier-restriction bound for prime sums** than
existing Vinogradov-Korobov estimates. If primes saturate exactly the
worst-case bound, primes are decoupling-extremal; both outcomes are
informative. The ε-stability of the BDG bound makes it discriminating
even at finite N.

**Distinction from line 309 closure (Circle method / Vinogradov):**
that closure rules out direct circle-method extraction of π(x). D15
asks the *INEQUALITY-saturation* question: how tight is BDG's bound
on prime sums empirically? This is a structural test of χ_P, not an
algorithm. The two questions are independent.

**Cross-domain ingredient:** Bourgain-Demeter-Guth 2015 decoupling
theorem (arXiv:1512.02384, Annals of Math 2015), Demeter 2020
*Fourier Restriction, Decoupling and Applications* (Cambridge),
Vinogradov mean value resolution. UNUSED in §7 ("Decoupling"
candidate row).

**Concrete first step:** for `N ∈ {2^{12}, 2^{14}, 2^{16}}`, build
`f_N(x) = (1/√π(N)) Σ_{p ≤ N} e^{2π i p x / N}`. Partition the unit
circle into `M = 2^{j}` caps `[k/M, (k+1)/M]` for `j ∈ {4, 6, 8}` (so
cap size `δ = 1/M`). For each cap, compute `f_N|_{cap_j}` (the
restriction of `f_N` to frequencies in cap j) by Fourier-band-pass,
then evaluate `||f_N|_{cap_j}||_{L^4(unit interval)}` via FFT
(`L^4 = (∫ |g|^4)^{1/4}`). Compute the empirical decoupling constant
`K_emp(N, M) := ||f_N||_{L^4}/(Σ_j ||f_N|_{cap_j}||_{L^4}^2)^{1/2}`.
Compare to the BDG bound `K_BDG(δ) = δ^{-ε}` and to:
(a) random-density-matched Bernoulli sample at density `1/log N`;
(b) Liouville-supported `f_N^{lambda}(x) := Σ_n lambda(n)/√N · e^{2πinx/N}`
    (centred multiplicative; should be decoupling-extremal per Möbius
    orthogonality);
(c) arithmetic-progression-supported `f^{AP}_N` (known: AP supports
    saturate a *trivial* decoupling bound).
Plot `log K_emp / log(1/δ)` (effective ε exponent) vs `log N`. Is
the prime exponential sum's effective-ε significantly *smaller* than
random matched-density (better-than-BDG bound for primes) or larger
(decoupling-extremal)?

**Failure profile (E):** `K_emp(prime) = K_emp(random)` within sample
noise — primes saturate BDG at the random level; closes as 38th
pseudorandomness measure, **NEW kind** (Fourier-restriction `L^p`
geometry). **B-grade**. **(I):** `K_emp(prime) > K_emp(random)` —
primes are *worse* than random under decoupling, indicating cap
correlations from HL singular series; structurally interesting but
likely reduces to Vinogradov-Korobov bounds via duality. **(INC):**
ε-discrimination at finite N is borderline; need `N ≥ 2^{18}` for
clean signal at `δ = 2^{-8}`.

**A-grade success:** `K_emp(prime) ≤ (log N)^c` for some explicit
`c > 0`, strictly better than random matched-density `(log N)^{c+1}`
or worse — gives a NEW Fourier-restriction bound for prime
exponential sums tighter than Vinogradov-Korobov, which directly
constrains the explicit-formula error term in `π(x) − Li(x)`.
**Quantitatively new analytic-NT result.**

**B-grade success:** `K_emp(prime)` deviates from random matched-
density by a quantity reproducing the HL singular-series structure
(specifically: cap pairs `(j, k)` with `(jk) mod q` concentrated for
small q have anomalous joint contribution) — first **decoupling-
saturation** characterisation of χ_P. Adds an EDGE in the Fourier-
restriction category (structurally distinct from Gowers norms which
test additive-shift uniformity).

**Cross-domain refs:**
- Bourgain-Demeter-Guth 2015 "Proof of the main conjecture in
  Vinogradov's mean value theorem for degrees higher than three"
  Annals of Math 184 (2016), 633
  https://arxiv.org/abs/1512.02384
- Bourgain-Demeter 2015 "The proof of the l^2 decoupling conjecture"
  Annals of Math 182 (2015), 351
  https://arxiv.org/abs/1403.5335
- Demeter 2020 *Fourier Restriction, Decoupling and Applications*
  Cambridge Studies in Advanced Mathematics 184
- Wikipedia: Vinogradov's mean-value theorem
  https://en.wikipedia.org/wiki/Vinogradov%27s_mean-value_theorem

**Budget:** 1-2 sessions. FFT-based cap restriction at `N = 2^{16}`,
`M = 2^8` is `O(N log N · M) = 10^{11}` ops — borderline single-core
overnight. Single-session feasibility for `N = 2^{14}, M = 2^6`.

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

## §G. The Multiplicative Regime — Beyond the W-Trick Wall

**Why this section exists:** S87 (Gowers `U^k` of χ_P), S88 (Anderson
Lyapunov of χ_P-Schrödinger), and S89 (parity-calibrated depth-2
threshold circuits) jointly established a triple confirmation that
**every detectable deviation of χ_P from random is exactly the Hardy-
Littlewood singular series mod q, fully removed by the Green-Tao
W-trick**. The additive / spectral / circuit-complexity regime on the
*bare prime indicator* is now a closed family. Future productive
attacks must move to the **multiplicative regime**: Liouville `λ(n)`,
Möbius `μ(n)`, von Mangoldt `Λ(n)` truncated, etc. These are predicted
to be Gowers-uniform *without* W-tricking (Green-Tao 2012 Möbius/
nilsequence orthogonality theorem). Any deviation found here is fresh
territory the W-trick framework cannot explain.

The next-action recommendation from S88 explicitly flagged this pivot.
This section is its concrete realisation.

### G1 — Liouville Anderson Lyapunov: spectral signature of Möbius pseudorandomness [CLOSED S100, see "Closed attacks"]

**Question:** repeat S88's Schrödinger experiment with potential
`V(n) = λ(n) ∈ {-1, +1}` (Liouville function, completely multiplicative
on primes). Is `γ_λ(E) - γ_Rademacher(E)` zero within sample noise
across all energies AND without any W-trick?

**Why frontier:** the W-trick story for χ_P (S88, E2.14) is driven by
small-prime residue-class structure. Liouville is *centered*
(`E[λ] = 0` asymptotically) and is the canonical Möbius-nilsequence-
orthogonal sequence. If `γ_λ - γ_Rademacher = 0` to machine precision
*at all energies* — that's the spectral analogue of GT Möbius
orthogonality, and the project would have *direct empirical
confirmation* that Liouville carries no detectable arithmetic structure
visible to Anderson-Lyapunov machinery.

If, however, `γ_λ ≠ γ_Rademacher` at some energy (without W-tricking)
— that's an A-grade result: the first project measure showing
Liouville/Möbius deviation that the W-trick framework cannot explain.

**Cross-domain ingredient:** Anderson localisation theory (already
imported S88), Möbius/nilsequence orthogonality (Green-Tao
arXiv:0807.1736).

**Concrete first step:** modify
`experiments/dynamical/anderson_localisation_chi_p/anderson_localisation_chi_p.py`
to take potential `V(n) = lambda(n)` (computable in O(N log log N)
via the Liouville sieve). Run at N = 10^5 with 50 seeds of matched
random `±1` signs (Rademacher) and 51 energies. Compute pointwise
z-score sweep.

**Failure profile (E mode):** `γ_λ ≈ γ_Rademacher` across all energies
to within sample noise. Adds a 38th pseudorandomness measure
*confirming* GT Möbius orthogonality at the spectral level. **B-grade.**

**A-grade success:** any sustained `|z| > 5` deviation at any energy,
NOT W-tricked away by sieving λ on coprime-to-W subsets. Would falsify
the "everything reduces to HL via W-trick" picture.

**Budget:** 1 session (existing infrastructure).

### G2 — Liouville Gowers U^k norms: testing the Green-Tao orthogonality theorem at scale

**Question:** parallel to §D6 (S87, E2.13), compute `‖λ‖_{U^k}` for
k = 2, 3 at large N. Green-Tao predicts `‖λ‖_{U^k} → 0` polynomially
in `1/log N`. Empirical *rate* of decay has not been measured for the
bare Liouville on Z/NZ.

**Why frontier:** combined with G1, gives two independent measures of
"is there ANY structure in λ that the additive/spectral W-trick
framework misses". Empirical decay-rate of `‖λ‖_{U^k}` would also
give a numerical handle on the size of the implied error term in
GT's theorem.

**Cross-domain ingredient:** Gowers norms (already imported S85),
GT-Ziegler U^{s+1} inverse theorem (arXiv:1009.3998), Möbius-
nilsequence orthogonality (arXiv:0807.1736).

**Concrete first step:** reuse `experiments/information_theory/gowers_uk_chi_p.py`
with `f(n) = λ(n)` instead of `χ_P(n)`. Note we want `λ(n) ∈ {-1, +1}`
(with mean 0), NOT the indicator `1[λ = -1]` (which S87 already
showed is Gowers-uniform at U^2). Run for N ∈ {2^14, 2^16, 2^18, 2^20}
at U^2; N up to 2^14 at U^3.

**Failure profile (E mode):** `Q^k(λ) → 0` polynomially in `1/log N`,
matching Green-Tao prediction. **B-grade**: 38th-39th pseudorandomness
measure with closed-form prediction.

**A-grade success:** `Q^k(λ) = Ω(1)` at any tested k, OR sub-polynomial
decay rate, OR a stable energy where empirical sample mean differs
from GT prediction by more than `O(1/sqrt(N))`. Any of these would be
a deviation from Möbius orthogonality at the empirical level.

**Budget:** 1 session.

### G3 — Möbius Voronin universality: algorithmic content of `ζ(s)·1/ζ(s)`

**Question:** Voronin universality (S85, §B4) says ζ(s) approximates
every analytic non-vanishing function on a fixed disk. The MULTIPLICATIVE
inverse `1/ζ(s) = Σ μ(n)/n^s` is the Dirichlet series for Möbius. Does
`1/ζ(s)` (or some natural Möbius-side variant) admit *effective*
universality with **polynomial-rate** shifts, OR is it provably
super-poly?

**Why frontier:** Voronin universality of ζ has been numerically
explored, but the universality of `1/ζ` and other Möbius-side
zeta-related functions has not. If `1/ζ` admits poly-rate
universality, it gives a concrete polylog-time approximator for any
analytic function — including potentially `π(x)` itself via the
explicit formula `π(x) ~ Li(x) + Σ Li(x^ρ) · μ_ρ`.

If `1/ζ` is provably super-poly universal (no poly-rate shifts), it
would extend Garunkstis 2003's effective bound to multiplicative
zeta-data — a B-grade structural negative.

**Cross-domain ingredient:** zeta universality (Voronin 1975 / Steuding
2007 LNM 1877 / Garunkstis 2003 arXiv:math/0306072), Möbius-
nilsequence orthogonality (Green-Tao arXiv:0807.1736), L-function
effective approximation (Bombieri-Hejhal 1995 / Hejhal-Rackner 1992).

**Concrete first step:** numerical search. For target `f(s) = e^s` on
a small disk, find shifts `t_ε` such that `|1/ζ(s + i·t_ε) - f(s)| < ε`
on the disk. Compare empirical `t_ε` scaling to: (a) `1/ζ`'s analogue
of Garunkstis's bound for ζ, (b) polynomial in `1/ε`, (c) plain
`exp(1/ε)`.

**Failure profile (I mode):** `t_ε` super-polynomial in `1/ε`,
matching the heuristic that universality of `1/ζ` is no faster than
ζ's. **B-grade**: extends the negative-shape edge from B4 to the
multiplicative side.

**A-grade success:** `t_ε = poly(1/ε)` for some natural target `f`
on `1/ζ`. Would constitute the first known effective polynomial-rate
universality for any zeta-related function and be a candidate polylog
approximator route.

**Budget:** 2 sessions (Voronin numerics are heavy; the existing B4
infrastructure helps but the multiplicative side is fresh).

### Why §G is the project's current highest-leverage frontier

**Empirical case:** the additive/spectral/circuit-complexity regime
on χ_P is now a triple-confirmed closed family (E2.13 Gowers + E2.14
Anderson + E1.10/E3.13 prior + S89 calibrated controls). Every new
attack on χ_P pseudorandomness collapses to "yes this is Hardy-
Littlewood mod q, removed by the W-trick." The first three §G vectors
are designed to be *not W-tricked*: they target multiplicative objects
(λ, μ, 1/ζ) that the Green-Tao machinery predicts are uniform
WITHOUT W-tricking. Any deviation is genuinely new.

**Theoretical case:** the polylog π(x) explicit formula is
`π(x) = Li(x) - Σ_ρ Li(x^ρ) + small`, where `μ` controls the inverse
of ζ and the secondary terms. If μ-side has more structure than
predicted, the polylog frontier reopens via direct partial-sum
construction. If μ-side is Gowers-uniform exactly as predicted, the
polylog problem reduces cleanly to "compute the bulk of `Σ_ρ Li(x^ρ)`,"
a smaller and structurally-clean target.

**Sequencing:** G2 (Liouville Gowers) is the lowest-cost first step —
existing infrastructure, single script change. G1 (Liouville Anderson)
is also single-session. Both should fire in the next 2-3 cross-domain
or wild-swing slots. G3 (Möbius Voronin) is heavier and worth a
dedicated wild_swing.

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

### §A.A5 — Maynard 2015 multidimensional sieve weight as TC⁰ primality witness (CLOSED 2026-04-27, S116, mode E, B-grade)

**Outcome:** Maynard sieve weight `w(n) = (Σ μ(d_1)…μ(d_k) F(log d_i/log R))^2`
is **NOT** a polylog primality witness. Two independent obstructions:

1. **Aggregate-not-pointwise.** Maynard's main theorem says
   `Σ_{N≤n<2N} w(n) χ_P(n+h_i) ≥ c(k) Σ w(n)` for some `i ∈ [k]`.
   This empirically holds: aggregate ratio = 0.094–0.212 across θ.
   But pointwise AUC restricted to odd `n` (the only n that can
   contain primes given H = {0, 2, 6}) stays in **[0.66, 0.69]**
   across θ ∈ {0.10, …, 0.40} and four F choices including
   Maynard's symmetric polynomial. Best F1 = 0.62.

2. **Divisor-enumeration cost.** Mean coprime simplex tuple count
   per single-n evaluation scales as N^α for α ∈ {0.10, 0.11, 0.12}
   at θ ∈ {0.20, 0.30, 0.40} on N ∈ [10^4, 10^6]. Power-law,
   sub-θ but positive — **not polylog**. Listing squarefree divisors
   of `n + h_i` up to `R = N^θ` reduces to growing-dim modular
   powering / divisor enumeration — exactly E5.3.

**Mechanism:** Maynard's eigenvalue ratio `M_k/I_k` (= 1.515 for k=3
with optimal SDP F) controls *aggregate* weight efficiency, not
*pointwise* primality contrast. Maynard-symmetric F (with 2-parameter
optimization toward k=3 optimum) does no better than vanilla GPY
F = (1−Σx_i)². Pre-stated A-grade target (AUC_oddN ≥ 0.95 AND mean
ops = polylog) is empirically falsified at all θ tested.

**What this rules out:** Maynard sieve is the most refined explicit
prime-detection machinery in modern analytic NT — the candidate that
came closest to "an unconditional TC⁰ primality test outside the
AKS family" in the project's enumeration. With this closure, the
sieve-route attack family is closed alongside AKS-modulus (E7.10),
Brandt MKtP (E5.8), and convergence-acceleration (E7.11) — these
form the four major families systematically exhausted on the only
open problem (status/OPEN_PROBLEMS.md polylog π(x)).

Adds **EDGE E7.14**: Maynard sieve weight is not a TC⁰ primality
witness (aggregate-not-pointwise + divisor-enumeration sub-poly).
Cites E6.7, E5.3, E7.10. CLOSED_PATHS row at session 116. See
`experiments/sieve/maynard_weight_pointwise/`.

**Successor proposal (different cross-domain technique):**
- **A5.a:** spectral sieve — instead of Maynard's *positivity-via-
  Selberg-square* approach, try a *spectral-projector approach*:
  define a Hilbert-space-valued sieve weight via the orthogonal
  projector onto the prime indicator's level set in `L^2(Z/NZ)`,
  using **Sarnak-Xue-style trace-class operators on automorphic
  forms** (cross-domain: representation theory / spectral theory of
  L-functions). Different obstruction class than Maynard's
  divisor-enumeration cost.

### §C.C5 — Stein's method finite-x Wasserstein test of D(x)=(π(x)-Li(x))log(x)/√x (CLOSED 2026-04-27, S108, mode E, A-grade-provisional)

**Outcome:** the §C5 question — does empirical `W_1(D̂_K, N(μ̂, σ̂²))`
scale as Stein-CLT `O(1/√K)` or plateau at `c > 0`? — is answered:
**plateau detected**, structurally explained.

Protocol: K log-uniform anchors `x_k ∈ [10^6, 10^7]`, `D(x) :=
(π(x) − Li(x)) · log(x) / √x`. W_1 computed in closed form via sorted-
CDF sub-bin trapezoidal integration; controls drawn from N(μ̂, σ̂²)
with sample-fitted parameters per trial.

| K     | W_1(D̂_K, N(μ̂, σ̂²)) | W_1 Gaussian (sample-fit) | Z-score |
|-------|----------------------|----------------------------|---------|
| 1000  | 0.0087               | ~0.0019 (1/√K rate)        | ≥ 3σ    |
| 5000  | 0.0087               | ~0.00088                   | ≥ 5σ    |
| 10000 | 0.0083               | ~0.00062                   | ≥ 5σ    |

Empirical D̂'s W_1 plateaus at `c ≈ 0.0083 (–0.0087)`; Gaussian-control
W_1 shrinks as `c_G/√K` as predicted by Stein-CLT. The plateau is
NOT consistent with a centred Gaussian limit at finite x.

**Structural explanation (this is what makes the closure mode E):**
the plateau is matched within ≤ 5% by the contribution of the lowest
~3 non-trivial Riemann zeros via the explicit formula

```
D(x) = -2 Σ_k cos(γ_k log x − arctan(2 γ_k)) / √(¼ + γ_k²)
       − log(2)·log(x)/√x  +  small remainder.
```

Across 10 sub-windows of `[10^6, 10^7]` of width 0.5 in log10:
- Empirical `W_1(D̂_emp window)` ranges 0.0075 – 0.0211
- Theoretical `W_1(D_th(50 zeros) window)` ranges 0.0141 – 0.0267
- **Pearson correlation r = 0.906** (n=10 windows)

A *random-phase* control (same γ_k weights, uniform random phases
replacing the actual zero-phase structure) gives `W_1 = 0.0116 ± 0.003`
across 50 trials; empirical D̂ is z = −1.06 vs this null — *not*
distinguishable from a generic random-phase low-zero sum.

The negative-excess-kurtosis signature `kurt(D̂) = −0.41 (95% CI
[−0.46, −0.36] at K=10000)` confirms sub-Gaussianity sourced from
arcsine-distributed individual cosine modes.

**Why §C5 is closed (mode E — explicit-formula reduction):** the
finite-x Wasserstein deviation reduces *quantitatively* to the
contribution of low Riemann zeros via the explicit formula. No
arithmetic structure orthogonal to the explicit formula is detected.
The metric (Wasserstein-1) is new for π(x) - Li(x) at finite x —
this is a quantitative refinement of E1.5 / E7.1 — but it joins the
GUE-sieve-circuit closure family rather than opening a new
bit-extraction angle.

**Why this is A-grade provisional (despite closure mode E):**
- (a) the cross-domain technique (Stein's method) had never been
  applied to π(x) - Li(x) in this project;
- (b) this is the first quantitative finite-x Wasserstein-shape
  bound for `π(x) - Li(x)` — orthogonal to the asymptotic Selberg /
  Hejhal CLT and to the pointwise-discrepancy bounds of Pintz /
  Korevaar;
- (c) the structural matching (correlation 0.906 across sub-windows
  with the explicit-formula prediction) is non-trivial and concrete,
  not just curve-fitting (the constant `α = 1.029 ≈ 1` is the
  analytic prediction, not a fit parameter);
- (d) the result satisfies §C5's verbatim A-grade success criterion.

If the verify-mode adversarial check finds a flaw (e.g., the
plateau collapses at higher K or higher x, or the structural
correlation is an artefact of how zeros are weighted in the explicit
formula), this gets demoted to B-grade (refinement of E1.5 with
explicit Wasserstein bound).

**Successor entries to consider:**
- C5b: same machinery, x ∈ [10^7, 10^8] — does `c(X)` shrink
  monotonically? (Initial test at K=1000: c ≈ 0.0067, smaller than
  c(10^6) ≈ 0.0087 — suggests `c(X) → 0` as `X → ∞`, consistent
  with asymptotic Hejhal CLT).
- C5c: replace D(x) with `(π(x) mod 2) · log(x)`-style discretised
  versions and test Wasserstein plateau on those (does the discretised
  version inherit the same low-zero-driven plateau?).

Pointer: `experiments/analytic/stein_wasserstein_pi/`,
`novel/finite_x_wasserstein_plateau.md`.

### §D.D13 — Subword complexity / topological entropy of χ_P (CLOSED 2026-04-27, S104, mode E)

**Outcome:** structural negative + new EDGE **E2.19** + 38th
pseudorandomness measure + sixth orthogonal HL-detection family
(symbolic-dynamics / factor-complexity / topological-entropy).

Protocol: at N = 5 · 10^6, n_max = 22, K = 20, five chi_P-derived
binary streams (RAW = chi_P(2..N); ODD = chi_P(2k+1); W{q} = Green-
Tao W-trick at primorial q ∈ {6, 30, 210}, residue 1). For each
stream and each n ∈ [1, 22], compute subword complexity p_w(n) :=
#{distinct length-n factors of w} via vectorised rolling encoding
(numpy uint64; cost O(nL); memory O(L)). Compare to K=20 random
PERMUTATIONS of the stream (B2: preserves density and 1-marginal
exactly, kills all >1-gram structure) and K=20 iid Bernoulli
matched-density samples (B1).

**Pre-registered F3 falsifier** (PRIMES > 3σ from BOTH B1 AND B2 at
every n in [n_lo, n_hi] with n_hi ≥ 18, on at least one of {ODD, W6,
W30}) **HOLDS** at z >> 3σ for ODD, W6, W30. W=210 erases to
|z_perm| ≤ 8.4σ at the worst n.

**Cascade table** (max |z_perm| over n ∈ [1, 22]):

| Stream | density | L         | max\|z_perm\| | at n | p_chi(22)/p_perm(22) |
|--------|--------:|----------:|--------------:|-----:|---------------------:|
| RAW    | 0.0697 | 4 999 999 | 132.7 | 22 | 0.018 (98 % deficit) |
| ODD    | 0.1394 | 2 499 999 | 277.1 |  7 | 0.216 (78 % deficit) |
| W6     | 0.2090 |   833 334 | 120.5 |  8 | 0.806 (19 % deficit) |
| W30    | 0.2611 |   166 667 |  24.8 | 17 | 0.994 (≈ noise)      |
| W210   | 0.3053 |    23 810 |   8.4 | 12 | 1.011 (≈ noise)      |

Effective topological entropy `h_eff(n) = log_2 p_w(n) / n` of the
W=210 stream agrees with Bernoulli matched-density to ≤ 0.001 across
n ∈ [1, 22] (saturation regime at n=22: log_2(L)/22 = 0.661, h_chi =
0.6581, h_perm = 0.6574).

**Mechanism:** subword complexity counts distinct prime-position
configurations inside a sliding length-n window; chi_P configurations
are arithmetically constrained by mod-p admissibility for every prime
p ≤ n (only one prime-multiple of p can lie in the window unless the
window contains p itself). Green-Tao W-trick at W = primorial(p_k)
sieves out integers not coprime to p_1 ... p_k, removing exactly
those constraints. Same Hardy-Littlewood k-tuple-admissibility engine
that drives E2.13 Gowers, E2.14 Anderson, E2.16 DPP, E2.17 PH.

**Why §D.D13 is closed (mode E):** the subword complexity at fixed
window length is `O(L log L)` (rolling encode + sort) — no polylog
opening. The W-trick fingerprint encodes only HL singular-series
information already known via E2.13. New mathematical category
(symbolic dynamics / factor complexity), same physics.

CLOSED_PATHS line 181 ("Symbolic dynamics — near-random block
complexity, S4") is now **promoted** from informal placeholder to
a precise quantitative measurement with finite-N value, K=20
baselines, and full W-trick cascade.

See `experiments/dynamical/subword_complexity_chi_p/` and
`archive/sessions/session104_d13_subword_complexity.md`.

**Follow-ups proposed** (per CLAUDE.md self-extension; one successor
with a DIFFERENT cross-domain technique):

- (D13.a) **Scale to W = 2310** (next primorial). Requires N ≥
  5 · 10^7 to keep stream length L ≥ 2 · 10^4. Tests whether the
  W=210 residual (|z_perm| = 8.4σ at n=12) collapses to ≤ 3σ or
  persists as genuine higher-mod structure.
  *Cross-domain technique*: same (factor complexity).
  *Why useful*: tightens the W-trick fingerprint convergence.
  Single session.
- (D13.b) **Subword complexity of `1[lambda(n) = -1]`** binarised
  Liouville — predicted no W-trick needed (cf. E2.18, where Anderson
  Lyapunov of Liouville matches Rademacher without W-trick).
  *Cross-domain technique*: same (factor complexity) + Möbius
  orthogonality framework.
  *Why useful*: extends the E2.18-style multiplicative-regime
  pattern to a sixth measure category. Single session.
- (D13.c) **Lempel-Ziv complexity LZ78 / LZ77 of chi_P** — distinct
  but related cross-domain measure. Compresses the FULL infinite
  word into a parsing tree; LZ77 length / N is an estimator of the
  measure-theoretic entropy `h_μ`. *Cross-domain technique*: data
  compression / Kolmogorov complexity surrogate (Ziv-Lempel 1977
  IEEE Trans. Inf. Theory 23). UNUSED in §6 of
  CROSS_DOMAIN_TECHNIQUES.md as a proper *quantitative* edge —
  only mentioned generically. Recommended successor introducing a
  new technique to the registry. Single session.

### §G.G1 — Liouville Anderson Lyapunov: spectral signature of Möbius pseudorandomness (CLOSED 2026-04-27, S100, mode E)

**Outcome:** structural negative + new EDGE E2.18 + first
**non-W-tricked** spectral measurement at noise floor + first
project use of cross-domain technique "Möbius/nilsequence
orthogonality" (Green-Tao 2012, Sarnak 2010, Tao 2016 logarithmic-
Chowla).

Protocol: discrete 1D Schrödinger operator
`H psi(n) = -psi(n+1) - psi(n-1) + V(n) psi(n)` on Z, transfer
matrix `T_n(E) = [[V(n) - E, -1], [1, 0]] in SL(2, R)`, Lyapunov
exponent `gamma(E) = lim_N (1/N) log ||T_N ... T_1||`. **Centered
multiplicative encoding** `V(n) = lambda(n) ∈ {-1, +1}` (mean → 0
by PNT for Liouville, variance = 1 exactly), Pastur-Figotin
prediction in-band `gamma_PF(E) = 1/(8 sin^2 k)`. Compared against
i.i.d. Rademacher ±1 baseline. 51 energies in [-2.95, 2.95]. Three
sample sizes: N ∈ {10^5 (50 seeds), 3·10^5 (50 seeds), 10^6 (100
seeds)}.

**Pre-registered F1 falsifier (sustained `|z| > 5` not removed by
multiplicative W-trick) FALSE.** Empirical results:

| N        | seeds | max |z| | argmax E  | χ²/K  | L²-rank λ |
|----------|-------|---------|-----------|-------|-----------|
| 10^5     | 50    | 1.78    | -0.236    | 0.63  | 21 / 50   |
| 3·10^5   | 50    | 2.16    | +0.118    | 0.49  |  7 / 50   |
| 10^6     | 100   | 2.04    | -2.006    | 0.69  | 41 / 100  |

`max |z|` flat in N (well below 51-energy Bonferroni threshold
3.16). Argmax energy WANDERS — statistical, not arithmetic. χ²/K <
1 throughout (lambda is sub-Rademacher on the spectral chi^2
aggregate). Pastur-Figotin agreement: γ_λ/γ_PF = 0.9317 (std 0.32),
identical to γ_Rademacher/γ_PF = 0.9309 to 4 decimals.

**Independent off-spectral check.** Two-point Chowla aggregate at
N=10^6, h=1..16: Σ z_h² = **4.77** vs Rademacher χ²_16 mean 16, std
√32 ≈ 5.7 (empirical p ≈ 0.997). Lambda is *more Rademacher-like
than Rademacher* on this independent test. All 16 individual
|z_h| ≤ 1.11. Consistent with Tao 2016 logarithmic-Chowla.

**Stark contrast with E2.14 (chi_P, S88) at the same N=10^5 grid:**
chi_P max |z| = 88.5σ at LOCKED E=0.108 (parity resonance) and a
secondary E=1.088 (mod-3 resonance), required W=2310 sieve to
reduce to ~4σ. **Lambda needs no W-trick at all.** The contrast IS
the new content: chi_P's spectral deviation is exclusively HL
singular-series mod-q resonance, and the canonical multiplicative-
regime companion (lambda) carries no such resonance.

**Why §G.G1 is closed (mode E — structural match to Möbius/Sarnak
orthogonality + Tao logarithmic-Chowla):** the experiment was
designed to be the spectral analogue of GT's Möbius/nilsequence
orthogonality theorem (arXiv:0807.1736). Lambda matching Rademacher
to within seed noise IS that orthogonality at the spectral level.
No deviation found ⇒ the intended polylog opening (lambda-side
explicit-formula partial sum) is closed.

**Adds EDGE E2.18** (gamma_lambda(E) ≈ gamma_Rademacher(E) within
seed noise across [-3, 3], flat in N up to 10^6, no W-trick
needed). FIRST non-W-tricked spectral measure at noise floor in the
project's 38-measure pseudorandomness battery.

See `experiments/dynamical/liouville_anderson_lyapunov/` and
`archive/sessions/session100_g1_liouville_anderson_lyapunov.md`.

**Follow-ups proposed** (per CLAUDE.md self-extension; one
successor with a DIFFERENT cross-domain technique):

- (G1.a) **Möbius `μ(n) ∈ {-1, 0, +1}` instead of lambda.** μ has
  density 6/π² ≈ 0.61 (squarefree) and removes prime-power signs.
  Predicted: μ-Anderson Lyapunov also matches Rademacher at scale.
  *Cross-domain technique*: same (Anderson + Möbius orthogonality).
  *Why useful*: confirms the spectral closure is robust to the
  squarefree-thinning step, not an artefact of lambda's specific
  parity-of-Ω structure. Single session.
- (G1.b) **Liouville × characters: `V(n) = λ(n) χ(n)` for χ a
  non-trivial Dirichlet character mod q.** Tests whether multiplying
  by a small-mod-q twist breaks the Möbius orthogonality picture.
  *Cross-domain technique*: extends to twisted-multiplicative
  setting (S56 character-twisted Liouville closures suggest this
  is also at noise floor — direct check needed). Single session.
- (G1.c) **`V(n) = log p_n` along primes** (different from G1: V
  supported on primes only, but with weight log p_n). Distinct from
  chi_P (uniform weight) and from lambda (full-support). Cross-
  domain ingredient: Pastur-Figotin generalisation to non-binary
  potentials. *Why useful*: introduces a NEW cross-domain technique
  (heavy-tail random-Schrödinger spectra; cf. Bourgain-Goldstein-
  Schlag 2002 *Annals* 154 for log-bounded potentials). Single
  session, but the ingredient genuinely extends the project's tool
  set. Recommended successor.

### §D.D2 — Topological data analysis (persistent homology) on prime gap sequence (CLOSED 2026-04-27, S96, mode I)

**Outcome:** structural negative + new EDGE E2.17 + fifth orthogonal
HL-detection category (algebraic-topological / metric-space geometry,
complementing E2.13 Gowers, E2.14 Anderson, E2.15 alg. immunity, E2.16
DPP failure).

Protocol: window of M consecutive Cramér-normalised gaps
`x_n = (p_{n+1} - p_n) / log(p_n)` near `p ≈ 10^6`, Takens delay
embedding `y_n = (x_n, ..., x_{n+d-1})` at `tau = 1`, Vietoris-Rips
filtration via ripser (`thresh = 4`), summaries on persistence
diagram: `T0` = total finite H_0 persistence, `T1` = total H_1
persistence, `L0/L1` = max persistence in each. Two baselines, K=20
each: B1 = IID Exp(1) (Poisson process); B2 = random permutation of
empirical window (preserves gap marginal, kills serial correlation).

**Pre-registered F3 falsifier (PRIMES > 3σ from BOTH B1 AND B2 on at
least one summary) HOLDS** at d=3, M=2000, x ≈ 10^6:

| Statistic | PRIMES | B1 mean | B2 mean | z(B1) | z(B2) |
|-----------|--------|---------|---------|-------|-------|
| **T0**    | 243.34 | 349.32 ± 10.28 | 277.43 ± 3.92 | **−10.31** | **−8.70** |
| **T1**    |  37.24 |  45.09 ± 1.87  |  56.09 ± 1.57 |  **−4.20** | **−11.99** |
| L0        |   1.77 |   2.39 ± 0.64  |   1.68 ± 0.27 |   −0.96    |   +0.35    |
| L1        |   0.37 |   0.45 ± 0.08  |   0.37 ± 0.04 |   −0.93    |   +0.05    |

PRIMES T0, T1 rank 0/20 in BOTH baseline distributions.

Robustness:
- Different window (x ≈ 5·10^6, same M, d): T0 z(B2) = −7.58,
  T1 z(B2) = −8.69. Effect persists.
- Different embedding dimension (d ∈ {2, 3, 4}): T0 z(B2) ∈ [−8.7,
  −5.1], T1 z(B2) ∈ [−12.0, −3.6]. Effect dimension-stable.
- M-scaling at d=3: T0 z(B1) = −4.2 (M=500), −6.0 (1000), −10.3
  (2000), −17.8 (4000). Z-scores grow super-linearly — signal is at
  least linear in window size, not finite-N artifact.

**Mechanism:** HL k-tuple admissibility constrains consecutive gaps
to repeat residue patterns more often than random, creating
geometric self-similarity in the delay-embedding cloud (small T0 =
clusters merge faster) and suppressing random "out-and-back" loops
in delay space (small T1). The B2 control preserves the gap MARGINAL
but destroys serial correlation — the deficit relative to B2
isolates the *serial-correlation* component of the deviation.

**Why §D.D2 is closed (mode I — structural negative on the
algorithmic question):** PH is a measurement instrument, not an
algorithm. VR-PH costs O(M^2) distance + O(M^3) worst-case
persistence; no closed-form polylog evaluator suggested. The
deviation is REAL and large, but has no algorithmic implication
for π(x).

**Adds EDGE E2.17** (PH of Takens-embedded normalised prime gaps
deviates ≥ 5σ from both Poisson and gap-permuted controls). FIFTH
INDEPENDENT CONFIRMATION of HL equidistribution structure in chi_P
in a NEW mathematical category (algebraic-topological / metric-space
geometry).

See `experiments/topological/persistent_homology_chi_p/` and
`archive/sessions/session96_d2_persistent_homology_chi_p.md`.

**Follow-ups proposed** (per CLAUDE.md self-extension):
- (D2.a) PH of W=210 W-tricked normalised gaps. Predicted: T0/T1
  deficit reduces toward Poisson if HL serial structure is the only
  mechanism (linking E2.17 to E2.13).
- (D2.b) Persistence-image vector classifier on PRIMES vs B1 vs B2
  windows for interpretable spectral discrimination axis.
- (D2.c) Sliding-window χ_P indicator embedding (Hamming/Manhattan
  distance) instead of gap embedding — different signal channel.

### §D.D7 — Determinantal point process (DPP) fit to the integer prime sequence (CLOSED 2026-04-26, S95, mode I)

**Outcome:** structural negative + new EDGE E2.16 + first 3-point
structural confirmation that chi_P deviation IS the HL singular
series (complementing pair-level closures E2.13 / E2.14 / E2.15).

Tested four progressively-flexible kernel hypotheses:

(1) **Real DPP at pair level**: `K^2_DPP(t) := rho^2 - R_2(t) >= 0`
required. Empirically `K^2_DPP(t) < 0` for ALL 15 admissible even t in
[2, 30] at N = 10^7 (HL gives `S(0, t) > 1`, hence `R_2 > rho^2`).
Pair-level DPP infeasible.

(2) **Real PPP at pair level**: `K^2_PPP(t) := R_2(t) - rho^2 >= 0`.
Empirically `K^2_PPP(t) < 0` for ALL 14 odd t > 1 (`R_2 ≈ 0`,
`rho^2 ≈ 4.4e-3`). Pair-level PPP infeasible. **Sign of `R_2 - rho^2`
flips between odd/even — no real translation-invariant kernel.**

(3) **Real-signed PPP on all-even sub-lattice**: PPP feasible at pair
level, test 3-point identity:
`R_3 = perm[K]` with `K(t) = ±rho sqrt(S(0,t)-1)`. Required sign cross-
term `sigma_req ∈ (-0.541, +0.769)`, NEVER `±1` (closest 0.769) at
any of 19 admissible all-even triples up to t_2 ≤ 26. Real signed K
of right magnitudes fails universally.

(4) **Complex Hermitian K on all-even sub-lattice**: phases free.
Least-squares fit of `phi : t -> [0, 2π)` over 13 distinct offsets,
200 random starts (LM + trust-region). Best max residual 0.0746,
≫ 0.01 sample-noise floor. **No global complex Hermitian phase
assignment matches HL.**

**Mechanism**: HL `S(0, t_1, ..., t_k)` factorises over PRIMES, each
factor `alpha_p` depending on `nu_p({0, t_1, ..., t_k})` =
#distinct residues mod p. DPP/PPP correlations factorise over PAIRS.
Pair admissibility ≠ triple admissibility: e.g., `(0, 4, 14)` is
pair-admissible mod 3 but `nu_3 = 3` saturates, giving `R_3^HL = 0`
while PPP predicts `1.02e-3`. Pure structural failure.

| Falsifier | Hypothesis | Result | Magnitude |
|-----------|------------|--------|-----------|
| F1 | DPP pair: K²_DPP < 0 at every even t | HOLDS | 15/15 |
| F2 | PPP pair: K²_PPP < 0 at every odd t > 1 | HOLDS | 14/14 |
| F3 | PPP 3pt overshoots HL by > 10% | HOLDS | 18/19 (max 79.16%) |
| F4 | sigma_req never ±1 (real signed K fail) | HOLDS | 19/19 (max \|σ\|=0.77) |
| F5 | No complex Hermitian phase fits globally | HOLDS | best res 0.0746 |

**Why §D.D7 is closed (mode I — structural negative):** primes are
quantitatively NOT a translation-invariant DPP / PPP / signed-real-K
/ complex-Hermitian-K point process. The HL prime-by-prime
factorisation cannot be reduced to pair correlations because
admissibility is a multi-body property. Same structural reason as
E2.13 (Gowers `U^k` -> HL singular series), E2.14 (Anderson Lyapunov
captured by W-trick), E2.15 (algebraic immunity = mod-4 fact), now
extended to the 3-point level in the random-matrix-theoretic
category.

**Adds EDGE E2.16** (primes are not a translation-invariant DPP/PPP/
signed-K/complex-Hermitian-K point process). PROVIDES FIRST
3-point STRUCTURAL CONFIRMATION (the prior three E2.13-E2.15 are all
2-point statements) that chi_P structure = HL equidistribution mod q,
in a fourth orthogonal mathematical category (point-process theory).

See `experiments/constructions/primes_dpp_ppp_fit/` and
`archive/sessions/session95_d7_dpp_ppp_fit.md`.

**Follow-ups proposed** (per CLAUDE.md self-extension):
- (a) **Pfaffian / matrix-valued kernel fit**: a Pfaffian point
  process (matrix-valued K with anti-symmetric block) admits richer
  3- and 4-point identities than DPP / PPP. Test whether W-tricked
  `chi_{210, 1}` (where pair correlations approach 1) is a Pfaffian
  process. Predicted: same structural failure (multi-body
  admissibility), but the falsification mechanism would be different.
- (b) **alpha-determinantal generalisation** (Vere-Jones 1997):
  `R_2(t) = rho^2 - alpha(t) K(t)^2` with `alpha(t) ∈ R` allowed
  offset-varying. The sign-flip across odd/even is then captured by
  `alpha(odd) = +1, alpha(even) = -1` (since `R_2(odd) < rho^2`
  requires +alpha and `R_2(even) > rho^2` requires -alpha). The
  question: do the 3-point alpha-identities then match HL? Predicted:
  still fails (same prime-vs-pair factorisation obstruction).
- (c) **A-grade reach**: derive a *non-pair-factorisable* K-like
  object whose multi-point statistics reproduce HL. Goal: a new
  arithmetic invariant of `chi_P` that factorises over primes (as
  HL does), not over pairs. If successful, this would be the first
  project object capturing the prime-by-prime structure. Speculative;
  budget 2-3 sessions.

### §C.C4 — Anderson localisation Lyapunov exponent of chi_P-driven Schrödinger operator (CLOSED 2026-04-26, S88, mode E)

**Outcome:** structural confirmation of Hardy-Littlewood / W-trick at the
Anderson-localisation Lyapunov-exponent scale + new EDGE E2.14 + spectral
analogue of E2.13 (S85 Gowers norms).

Setup: discrete 1D Schrödinger `H psi(n) = -psi(n+1) - psi(n-1) + V(n)psi(n)`,
transfer matrix `T_n = [[V(n) - E, -1], [1, 0]] in SL(2, R)`,
`gamma(E) = lim (1/N) log ||T_N ... T_1||`. Pastur-Figotin gives
`gamma(E) ~ rho(1-rho)/(8 sin^2 k)` for sparse `V in {0,1}` at density
`rho` inside band `E = -2 cos k in (-2, 2)`. Vectorised numerical
estimator with periodic L^2 renormalisation, O(N) per energy.

**Naive Bernoulli baseline at N=10^5, 50 seeds, 51 energies**:
`max |z|(gamma_prime - gamma_bern) = 88.5 sigma` at E=0.108 (k=pi/2,
parity resonance). Z-score grows as sqrt(N): 12.3 at N=10^4 -> 88.5
at N=10^5 — real bias, not finite-size noise.

**Confounder:** chi_P concentrated on odd indices (every prime > 2 is
odd). Parity-matched controls C1 (random odd-subset of size pi(N)-1
plus V(2)=1) and C3 (chi_P shuffled within odd indices) reduce to
33 sigma at E = 1.088 ~ -2 cos(2pi/3) = +1, the **mod-3 resonance**
(free transfer matrix rotates 2pi/3 per step, maximally couples to
mod-3-periodic potentials; chi_P has mod-3 structure).

**W-trick cascade** (random potential supported on
{n : gcd(n, W) = 1} plus delta at small primes p|W):

| W    | sieved        | max |z| (N=10^5)  | max |z| (N=2*10^5) |
|------|---------------|-------------------|---------------------|
| -    | -             | 88.5              | -                   |
| 2    | parity        | 32.7              | -                   |
| 6    | mod 2, 3      | 8.95              | 11.93               |
| 30   | mod 2, 3, 5   | 5.66              | 6.29                |
| 210  | mod 2..7      | 3.99              | 6.07                |
| 2310 | mod 2..11     | -                 | 3.96                |

Cascade ends at residual ~ 4 sigma (over 31 energies, borderline noise
after Bonferroni: `0.05 / 31 -> z_thresh ~ 3.16`).

**Why §C.C4 is closed (mode E):** chi_P's Anderson-Lyapunov deviation
from random IS the spectral signature of small-modulus residue-class
structure of primes, fully captured by the Green-Tao W-trick. Same
structural reason as E2.13 (Gowers `U^k` -> Hardy-Littlewood singular
series). gamma(E) extraction itself costs Theta(N) per energy, not
polylog. No deviation beyond HL prediction at any tested W within
sample noise.

**Adds EDGE E2.14** (Anderson Lyapunov gamma(E) of chi_P-Schrödinger
deviates from random by amount captured by W-trick). PROVIDES SECOND
INDEPENDENT CONFIRMATION (after E2.13/S85) that chi_P structure = HL
equidistribution mod q, in a category orthogonal to additive
combinatorics (spectral / transfer-matrix Lyapunov).

See `experiments/dynamical/anderson_localisation_chi_p/` and
`archive/sessions/session88_c4_anderson_localisation_chi_p.md`.

**Follow-ups proposed** (per CLAUDE.md self-extension):
- (a) **Modified-amplitude potential `V(n) = log(p_n)` if n prime else 0**:
  the original §C4 success-alt suggested log-weighted potential might
  give larger localisation. Pastur-Figotin scales as `Var(V)`, so a
  log-weighted potential boosts gamma by `(<log p>^2 / 1^2) ~ (log N)^2`
  amplification — could either (i) make the W-trick cascade easier to
  measure or (ii) reveal a deviation that the constant-amplitude version
  missed. C-grade (refinement) or B-grade (if deviation persists).
- (b) **Anderson localisation of `Lambda(n)` (von Mangoldt) potential**:
  Lambda is Gowers-uniform after centering (Green-Tao) so the cascade
  prediction is `gamma_Lambda - gamma_random -> 0` with NO W-trick needed.
  Direct empirical verification of the spectral analogue of Green-Tao's
  Lambda-vs-Mu nilsequence orthogonality. C-grade.

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

### §B.B1 — Slice rank / Croot-Lev-Pach polynomial method on chi_P (CLOSED 2026-04-26, S92, mode E)

**Outcome:** structural confirmation of Hardy-Littlewood / W-trick at the
algebraic-immunity level + new EDGE E2.15. THIRD INDEPENDENT
CONFIRMATION (after E2.13 Gowers, E2.14 Anderson) that chi_P structure
equals HL equidistribution mod q.

Three structurally distinct polynomial-method invariants tested:

(1) **Algebraic immunity AI(chi_P) over F_2** — direct LP
(Gauss elimination over F_2 on the `2^N x sum_i C(N,i)`
monomial-evaluation matrix), N=4..13. Plus annihilator extraction via
nullspace recovery.

(2) **F_p multilinear ANF degree** of chi_P viewed as `(F_p)^k -> F_p`
via base-p expansion, for `(p, k)` in `{(3,2..5), (5,2..3), (7,2..3),
(11,2)}`.

(3) **Slice-rank brackets** via Tao 2016 inequality
`slice_rank(T) >= max_axis rank(flatten_axis(T))` (LB) plus greedy
heaviest-slice peeling (UB), over F_2.

**Empirical AI(chi_P, N) raw** (no W-trick):

```
  N | rho_chi  | AI_chi_P |  AI_lam+ | AI_mu!=0 |  AI_random_mean
  4 |   0.3750 |        2 |        2 |        1 |            1.88
  6 |   0.2812 |        2 |        3 |        2 |            2.25
  8 |   0.2109 |        2 |        4 |        2 |            3.00 (std 0)
 11 |   0.1509 |        2 |        5 |        2 |            4.00 (std 0)
 12 |   0.1377 |        2 |        6 |        2 |            4.00 (std 0)
 13 |   0.1255 |        2 |        6 |        2 |              --
```

**AI(chi_P, N) = 2 for all N in [4, 13].** Random matched-density grows
as `Theta(log_2(1/rho))`, reaching 4 at N=11 with zero std.

**Mechanism (annihilator extraction, the SAME polynomial for all
N >= 5):**

```
g(x) = 1 + x_0 + x_1 + x_{0,1} = (1 + x_0)(1 + x_1)
```

`g(n) = 1` iff `n ≡ 0 mod 4`. The annihilation works because no prime
> 2 is divisible by 4 (and chi_P(2) gets killed by `bit_1(2) = 1`).
**AI(chi_P) = 2 is exactly the polynomial-method encoding of the
trivial mod-4 sieve fact.**

`Mobius!=0` (squarefree) inherits the same annihilator (n ≡ 0 mod 4
implies n is non-squarefree). Liouville+ does NOT inherit
(`lambda(4) = +1`, so 4 is in support).

**W-trick correction (`chi_P_{W,b}(n) = chi_P(W*n + b)`):**

```
    W   b   N |  rho_chi | AI_chi | AI_random_mean
    1   0   8 |   0.2109 |      2 |       3.00 (std 0)
    1   0  11 |   0.1509 |      2 |       4.00 (std 0)
    6   1   8 |   0.4531 |      4 |       4.00 (std 0)
    6   1  11 |   0.3535 |      5 |       5.00 (std 0)
   30   7  11 |   0.3823 |      5 |       5.00 (std 0)
   30  11  11 |   0.3789 |      5 |       5.00 (std 0)
```

**The deviation is fully removed by W >= 6**: AI(chi_P_{W,b}) tracks
AI(random matched-density) within 1 step at all tested N, and exact
match (zero std) at most cells.

**F_p multilinear ANF degree**: chi_P near-saturates max possible
degree `(p-1)*k` at all tested `(p, k)` other than `(3, 3)` and
`(5, 3)`, where chi_P drops by 1 below the random max — well within
"almost all coefficients populated" noise.

**Slice-rank brackets**: at p=2 non-informative (2-row flattenings
cap rank at 2 trivially, both chi_P and random); at p=3 chi_P matches
random for k>=3.

**Why §B.B1 is closed (mode E):** chi_P's polynomial-method deviation
from random IS the polynomial-method encoding of the trivial mod-4
sieve fact (no prime > 2 divisible by 4), fully captured by the
Green-Tao W-trick. Same structural reason as E2.13 (Gowers `U^k` ->
HL singular series, S85) and E2.14 (Anderson Lyapunov captured by
W-trick cascade, S88). AI(chi_P) extraction is exponential time in N
(LP over `2^N x sum_i C(N,i)` F_2 matrix), not polylog. No deviation
beyond mod-q at any tested W.

**Adds EDGE E2.15** (Algebraic immunity AI_F_2(chi_P) = 2 with
explicit annihilator `(1+x_0)(1+x_1)`, removed by W-trick at W >= 6).
PROVIDES THIRD INDEPENDENT CONFIRMATION (after E2.13 Gowers and
E2.14 Anderson) that chi_P structure = HL equidistribution mod q,
in a third orthogonal mathematical category (Boolean polynomial
method / algebraic cryptanalysis vs additive combinatorics vs
spectral / transfer-matrix Lyapunov).

See `experiments/algebraic/algebraic_immunity_chi_p/` and
`archive/sessions/session92_b1_algebraic_immunity_chi_p.md`.

**Follow-ups proposed** (per CLAUDE.md self-extension):
- (a) **Algebraic immunity of the multiplicative-regime functions**
  (Liouville `lambda(n) in {-1, +1}` centered, Mobius `mu(n) in {-1, 0, +1}`)
  per ATTACK_VECTORS §G — these are predicted to be Gowers-uniform
  WITHOUT W-tricking by Green-Tao Mobius/nilsequence orthogonality.
  AI(lambda_centered) and AI(mu) over F_3 (since codomain is 3-valued
  for lambda; over F_2 for the indicator
  `1[lambda = +1]`) untested. If AI tracks random WITHOUT W-tricking,
  fourth-leg confirmation in the multiplicative regime.
- (b) **Higher-degree structured annihilators**: extend the mod-4 fact
  to mod-q for `q in {6, 12, 30}`. Specifically: are there
  unexpectedly-low-degree degree-`O(1)` annihilators of `chi_P_{W,b}`
  beyond what the random matched-density baseline allows? At W=6, b=1
  the residual structure (n = p^2 for primes p > 5) might admit a
  small-degree annihilator over F_2 that random-Bernoulli-matched
  baselines lack. Concrete sub-falsification: AI(chi_P_{6,1}) <
  AI(random_{6,1}) by ≥ 1 at N >= 10 — would extend E2.15.

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

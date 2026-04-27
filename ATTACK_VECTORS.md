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

### B2 — Automorphic L-function basis  [CLOSED S118, see "Closed attacks"]

**Outcome (S118, mode E, B-grade — case 2 of pre-stated falsifiers):**
F_τ Hecke twisted partial sums of level-1 weight-12 Δ are ~3× more
obstructed from spanning π(x) − Li(x) than matched-multiplicative
random Sato-Tate ensembles at finite N (ratio 2.83 ± 0.02 across
9 (N, K) cells, Z = 17–58σ). Mechanism: Sato-Tate equidistribution →
{γ_k^Δ} and {γ_k^ζ} are GUE-distributed independently → Hecke is a
narrow-band basis at the wrong heights for fitting ζ-zero-driven
oscillations of y. Adds **EDGE E7.15**, FIRST automorphic-spectral
measurement of the project, refining E7.1 / E1.10 / E3.13 with
quantitative L(s, Δ)-vs-ζ independence ratio. Cross-domain technique
"Selberg trace formula" (§1 of CROSS_DOMAIN_TECHNIQUES) promoted
UNUSED → USED (E). See `experiments/algebraic/automorphic_l_function_basis/`
and Closed attacks below.

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

### C2 — Higher-order arithmetic corrections (Conrey-Snaith)  [CLOSED S123, see "Closed attacks"]

**Outcome (S123, mode I, B-grade):** First project measurement of
n-correlations at orders 4, 5, 6 at N=8000. Empirical zeta R_n along
equally-spaced slices, P_k(s) for k=0..5, and κ_n(L) for n=2..6 all
match GUE prediction within sample noise once matched-finite-N null is
used (gap-shuffled null gives spuriously huge z-scores because it
destroys GUE rigidity). The Conrey-Snaith arithmetic correction scales
as `1/L²` and is below the empirical noise floor `1/√(n_tuples)` at
N=8000, paralleling the §C1 / S71 closure at higher height. Refines
**EDGE E7.1** from "GUE up to order 3" to "GUE up to order 6 across
three independent probes." See `experiments/analytic/zeta_structure/n_correlations_4_5_6/`
and Closed attacks below.

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

### C6 — Pfaffian / α-determinantal point process structure of ζ zeros at order n=4

**Question:** the canonical correlation-function structure of the
zeta zeros has been tested as a *determinantal* point process (sine
kernel) up to order 6 (E7.1, S123). A strictly larger class of point
processes — *Pfaffian point processes* — uses a 2×2 *matrix-valued*
anti-symmetric kernel `K(x, y) = ((K_11(x,y), K_12(x,y)), (K_21(x,y),
K_22(x,y)))^T` with `K_12(x,y) = -K_21(y,x)`, and gives `n`-point
correlations as Pfaffians: `ρ_n(x_1, ..., x_n) = Pf[K(x_i, x_j)]_{i,j}`
(of a `2n × 2n` antisymmetric matrix). Pfaffian processes properly
contain determinantal as a degenerate case (`K_11 = K_22 = 0`); they
are the natural model for orthogonal/symplectic ensembles of random
matrices. The Vere-Jones α-determinantal generalisation
`ρ_n = (sgn α)^{n-1} det[K]_α` interpolates between α=1 (boson permanent)
α=0 (Poisson) α=-1 (DPP). Question: do the empirical 4-point
correlations of the 8000-zero sample satisfy any non-trivial PFAFFIAN
identity (over-determined Pfaffian-vs-determinant relation that fails
for the GUE/sine-kernel determinantal model), and which α-DPP value
fits best? An anti-Pfaffian deviation would be a NEW arithmetic edge
beyond the closed determinantal structure.

**Why frontier:** the C2/S123 closure showed empirical `R_n` for `n =
4, 5, 6` match the determinantal sine-kernel model within the noise
floor. That measurement leaves OPEN whether a *richer* point-process
structure (Pfaffian or α-DPP) fits the empirical correlations
*better* than the determinantal model — i.e., is the GUE
classification at `n = 4` a *unique best fit* or merely *one of many
admissible* fits, with the matrix-valued or α-deformation also
saturating? The Pfaffian identity is over-determined: at order 4,
the Pfaffian formula gives 3 independent linear combinations of the
3 distinct pair-of-pairs invariants, which forces a *specific
quadratic* relation among the empirical moments
`(R_2(s_1) R_2(s_2) - R_2(s_3)^2 + R_2(s_4)^2) = ?`
that does not hold for the determinantal sine kernel. The C2 successor
proposal C2.b explicitly flagged this — D24 picks it up.

**Cross-domain ingredient:** Borodin 2009 "Determinantal point
processes" *Encyclopaedia Math. Sci.* 152 = arXiv:0911.1153
https://arxiv.org/abs/0911.1153 §2.4–2.6 (Pfaffian processes &
matrix-valued kernels), Soshnikov 2000 *Russ. Math. Surv.* 55 =
arXiv:math/0002099 §3 (Pfaffian processes for orthogonal/symplectic
ensembles), Vere-Jones 1997 *New Zealand J. Math.* 26 (α-permanent
generalisation). UNUSED in §3 of CROSS_DOMAIN_TECHNIQUES (sine-kernel
DPP at α=-1 is closed; Pfaffian and other α values are PROPOSED).

**Distinction from C2 (S123 closure) and D7 (S95 DPP closure on
χ_P):** C2 tested the sine-kernel DPP fit to ζ zeros up to order 6
and found agreement within sample noise. C6 asks the *over-determined*
Pfaffian-fit question — a structurally stronger and more
discriminating test. D7 (S95) closed primes-as-DPP/PPP/signed-K with
mode I (primes are NOT a DPP); C6 asks the same Pfaffian/α-DPP
question on ζ ZEROS, an opposite signed object: zeros are the natural
sine-kernel DPP, so the question is whether they are MORE THAN a
DPP — the Pfaffian or α-DPP fit might be tighter.

**Concrete first step:** for the existing 8000 Riemann zeros (γ ≤ 8148,
Riemann-vM unfolded), at order n=4:
1. Empirically compute `R_4(s_1, s_2, s_3) := (1/n_ref) Σ_n
   1[γ_{n+i_j}/(γ_n+s_j) ∈ (1±ε), j=1,2,3]` via slice integration
   for s_i ∈ a 4×4 grid in [0.5, 5].
2. Compute the corresponding determinantal prediction
   `R_4^{det} = det[K(s_i - s_j)]_{4×4}` and the Pfaffian
   prediction `R_4^{Pf}` for the 2×2 matrix-valued kernel
   `K_{Pf} = ((sin πr/πr, ...), (..., -sin πr/πr))` of the GOE
   ensemble (Mehta 2004 §6.4) and the GSE ensemble.
3. Compute residuals `R_4^{emp} - R_4^{det}` and `R_4^{emp} - R_4^{Pf,GOE}`
   and `R_4^{emp} - R_4^{Pf,GSE}`. Test the PFAFFIAN identity
   `R_4^{Pf} = (R_4^{emp,GUE-style})^2 + (R_4^{emp,GUE-style,offset})^2`
   versus the (non-identity) DETERMINANTAL prediction. Bootstrap z-score
   from 20 sub-batches.
4. α-DPP: solve `R_4^{α=-1+δ}` for the best-fit `δ ∈ [-0.2, 0.2]`
   minimising `||R_4^{α} - R_4^{emp}||_2`. Report the optimal `α*`
   and its 95% CI (likelihood-ratio bootstrap).

**Failure profile (E):** GOE/GSE Pfaffian fit is INDISTINGUISHABLE
from sine-kernel determinantal fit at order 4 (both within noise
floor); α-fit gives `α* = -1.00 ± 0.05` consistent with GUE/sine.
Closes as: zeros are determinantal-typical, no Pfaffian / α-deformation
content. **B-grade.**
**(I):** the Pfaffian identity test fails because zeros are
correctly determinantal (not orthogonal/symplectic Pfaffian); reduces
to confirming GUE α=-1 classification. Adds a new specificity edge
"zeros are sine-kernel-DPP, NOT GOE/GSE Pfaffian." **B-grade.**
**(INC):** at n=4, sample noise on `R_4` is too high for `< 5%`
discriminating power on Pfaffian-vs-determinantal residuals; need
the Odlyzko `~10^6` zero block.

**A-grade success:** empirical `R_4` is BEST FIT by a Pfaffian /
α-DPP with `α ≠ -1` at `> 5σ` significance — zeros are not GUE-
determinantal but a proper deformation. Would falsify the canonical
GUE classification and constitute a structural discovery on ζ-zero
arithmetic content invisible to all 6-order sine-kernel measurements.

**B-grade success:** explicit upper bound on the Pfaffian-vs-determinantal
discrimination at order 4 in the form `||R_4^{emp} - R_4^{α=-1}||_2 ≥
c · ||R_4^{emp} - R_4^{α* (best-fit)}||_2` for explicit `c < 1`,
strengthening E7.1 from "DPP-typical" to "DPP-typical at *least* as
sharply as α-DPP" — first higher-α discrimination measurement on ζ
zeros.

**Cross-domain refs:**
- Borodin 2009 "Determinantal point processes" *Encyclopaedia Math.
  Sci.* 152 = arXiv:0911.1153 https://arxiv.org/abs/0911.1153
- Soshnikov 2000 "Determinantal random point fields" *Russ. Math.
  Surv.* 55 = arXiv:math/0002099 https://arxiv.org/abs/math/0002099
- Vere-Jones 1997 "Alpha-permanents and their applications to
  multivariate gamma, negative binomial and ordinary binomial
  distributions" *New Zealand J. Math.* 26
- Mehta 2004 *Random Matrices* (3rd ed.) Elsevier §6 (GOE/GSE
  Pfaffian formulas)

**Budget:** 1 session. Existing 8000-zero data + sine-kernel infra
from S123 are reusable. New work: Pfaffian residual computation
`O(n_ref · 4 · 16)` per slice point, and α-fit `O(n_grid · n_ref)`.
Decidable single session.

### C7 — Fyodorov-Hiary-Keating extreme-value statistics of |ζ(1/2 + it)| over short windows  [CLOSED S133, see "Closed attacks"]

**Outcome (S133, mode I, B-grade case 2):** FIRST ζ-amplitude (vs
zero-position) measurement of the project. M_T-mean is T-INDEPENDENT
(pooled −0.657 ± 0.045 across T ∈ {10⁴, 10⁵, 10⁶}, K = 100 windows
per anchor), confirming the FHK normalisation `log log T −
(3/4) log log log T`. M_T SHAPE at finite T ≤ 10⁶ is approximately
GAUSSIAN, NOT Gumbel(loc, 1/2): variance 1.47× the FHK prediction,
skewness ≈ 0.1 (Gumbel +1.14), excess kurtosis ≈ −0.4 (Gumbel +2.4),
KS to Gauss < KS to FHK Gumbel by factor 0.4–0.7 at all 3 T. **First
quantitative bound on FHK convergence rate at finite T**, in either
project-internal or published work; mean converges fast, shape converges
slow. Adds **EDGE E7.18**, first ζ-amplitude edge of the project.
Cross-domain "Gaussian multiplicative chaos / FHK extreme-value
statistics" (§3 of CROSS_DOMAIN_TECHNIQUES) promoted PROPOSED → USED
(I). See `experiments/analytic/zeta_structure/fhk_amplitude_max/` and
"Closed attacks" below.

**Question:** Fyodorov-Hiary-Keating 2012 (PRL 108, 170601 =
arXiv:1202.4713) conjectured — based on a freezing-transition analogy
to characteristic polynomials of random unitary matrices — that the
maximum of `log|ζ(1/2 + it)|` over a unit-length window
`[T, T + 1]` of the critical line satisfies, asymptotically,
`max_{t ∈ [T, T+1]} log|ζ(1/2 + it)| = log log T - (3/4) log log log T + M_T`
where `M_T` converges in distribution to a randomly shifted Gumbel
distribution. The conjecture has THREE distinct empirical signatures:
(i) the **leading-order constant** `log log T` (Selberg CLT predicts
this; FHK refines); (ii) the **secondary log-log-log correction**
(`-3/4` is FHK-specific, distinguishing log-correlated chaos from
plain Gaussian); (iii) the **Gumbel tail**. Saksman-Webb 2018
(arXiv:1609.00027) PROVED that ζ(1/2 + it) on a *mesoscopic* scale
converges to a Gaussian multiplicative chaos (GMC) measure — a
log-correlated random field. Question: empirically test the **FULL
joint distribution** `(max log|ζ|, position-of-max, second-max)` over
windows at `T ~ 10^4 - 10^7`, comparing to FHK's GMC prediction. Does
the empirical distribution match the GMC prediction within sample
noise (closes the FHK conjecture in mode E + GMC), or does the
Lambert-style position-of-max distribution carry an arithmetic
correction beyond the GMC prediction (A-grade — new ζ-amplitude
fingerprint)?

**Why frontier:** the project has measured ζ-zeros' *positions* as
a point process in 35+ ways (correlation, n-correlation, GUE rigidity,
spacing distributions). **The ζ-amplitude `|ζ(1/2 + it)|` between
zeros has NEVER been tested in this project**, and is the carrier of
the FHK / GMC conjectures. If the empirical max distribution matches
FHK, the project would have its first quantitative confirmation of
GMC universality at a finite scale. If it deviates, the deviation is
the first arithmetic-amplitude content of ζ — a *third* RM-vs-zeta
discrepancy beyond pair correlation and rigidity. The C2.c successor
proposal in S123 explicitly flagged this; C7 makes it a full
attack vector.

**Cross-domain ingredient:** Fyodorov-Hiary-Keating 2012 PRL 108 =
arXiv:1202.4713 (freezing transition + Gumbel max conjecture);
Saksman-Webb 2018 = arXiv:1609.00027 (GMC limit of ζ on mesoscopic
scale); Arguin-Belius-Bourgade 2017 *Comm. Math. Phys.* 349 =
arXiv:1612.08575 (rigorous max upper bound `log log T - (3/4 + ε)
log log log T` in random-model surrogate); Bourgade-Kuan 2014 *Adv.
Math.* (RMT-side max-of-char-poly result). UNUSED in §3 of
CROSS_DOMAIN_TECHNIQUES — first GMC import for the project.

**Distinction from §C1, §C2, S25/45/57/123 (zero-point-process
measurements) and S110/115 (Stein-Wasserstein on `D(x)`):** all
prior project work tested ZEROS or `π(x) - Li(x)` deviation. C7 is
the FIRST measurement on the ζ-AMPLITUDE between zeros — a complementary
geometric object whose distribution is governed by GMC, not GUE. The
two distributions are statistically distinguishable: GUE rigidity
sees zero-positions; GMC log-correlation sees inter-zero amplitudes.

**Concrete first step:** for `T ∈ {10^4, 10^5, 10^6}` and `K = 1000`
non-overlapping unit-length windows `[T_k, T_k + 1]`:
1. Compute `|ζ(1/2 + it)|` to 30 dps via mpmath at `M = 5000` evenly-
   spaced samples per window (`O(M · √T)` ops via Riemann-Siegel).
2. Extract `max_k = max_{t ∈ window k} log|ζ(1/2 + it)|`,
   `argmax_k`, and `secondmax_k`.
3. Empirical CDF of `max_k - log log T_k + (3/4) log log log T_k`
   (the FHK-renormalised max). Compare to (i) the FHK Gumbel
   prediction `F(M_∞ ≤ x) = exp(-c · e^{-2x})` for explicit shift
   constant c (Bourgade-Kuan 2014 prediction); (ii) the matched-
   Selberg-CLT Gaussian prediction (no GMC correction);
   (iii) bootstrap of 100 ζ-vertical-shift surrogates.
4. KS distance `D_max(empirical, FHK_Gumbel)`. Statistic: `Z =
   (D_max - mean(D_random)) / σ_random`. Bonferroni for 3 T-values.
5. Joint distribution: bivariate `(max_k, argmax_k - T_k)` to
   test FHK's prediction that `argmax` is *uniform* on the window.
   Is empirical `argmax` distribution uniform within sample noise,
   or does it concentrate near specific (e.g., zero-adjacent)
   positions?

**Failure profile (E):** empirical max distribution matches FHK
Gumbel prediction within sample noise across all `T` — closes as
38th (or 40th, depending on count) project-internal pseudorandomness
measure of ζ-amplitude, FIRST positive empirical confirmation of
the FHK conjecture at finite `T`, **B-grade**.
**(I):** empirical max systematically deviates from FHK by a
quantitative amount that reduces to a known correction (e.g.,
Conrey-Snaith arithmetic correction at the moment-generating
function level) — Pfaffian-style absorption into existing closures.
**(INC):** computing `|ζ|` at `T = 10^7` precision-limited; need
distributed computation or precomputed Odlyzko-style ζ tables.

**A-grade success:** empirical max distribution deviates from FHK
by `> 5σ` in any of the three signatures (leading constant,
log-log-log correction, Gumbel tail) AND the deviation has a
structural arithmetic explanation (e.g., the deviation factors
through the Hardy-Littlewood pair-correlation density at amplitude
scale `1/√L`). FIRST arithmetic-amplitude deviation of `|ζ|` from
the GMC universality prediction — opens the AMPLITUDE-side polylog
frontier (orthogonal to the closed POSITION-side family).

**B-grade success:** explicit empirical confirmation of FHK at
`T = 10^4 - 10^6` with quantitative bound on the secondary
correction `|max - log log T - α log log log T| ≤ c` for explicit
constant `c` and `α = 3/4 + δ` with `|δ| < 0.05`. First quantitative
finite-T FHK measurement; complements E7.1 (zeros are GUE) with a
quantitative ζ-amplitude statement.

**Cross-domain refs:**
- Fyodorov-Hiary-Keating 2012 "Freezing Transition, Characteristic
  Polynomials of Random Matrices, and the Riemann Zeta-Function"
  *PRL* 108, 170601 = arXiv:1202.4713
  https://arxiv.org/abs/1202.4713
- Saksman-Webb 2018 "The Riemann zeta function and Gaussian
  multiplicative chaos: statistics on the critical line"
  arXiv:1609.00027 https://arxiv.org/abs/1609.00027
- Arguin-Belius-Bourgade 2017 "Maximum of the characteristic
  polynomial of random unitary matrices" *Comm. Math. Phys.* 349 =
  arXiv:1612.08575 https://arxiv.org/abs/1612.08575
- Bourgade-Kuan 2014 "Strong Szegő asymptotics and zeros of the
  zeta function" *Comm. Pure Appl. Math.* 67
- Wikipedia: Gaussian multiplicative chaos
  https://en.wikipedia.org/wiki/Gaussian_free_field

**Budget:** 1-2 sessions. Session 1: implement ζ-evaluation
(Riemann-Siegel + mpmath at 30 dps), `T = 10^4 - 10^5`, `K = 1000`
windows, `M = 5000` samples per window. Session 2 (only if signal):
scale to `T = 10^6` via Odlyzko-zero-aided ζ tables.

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

### D10 — Mahler measure of the prime generating polynomial f_N(z) = Σ_{n≤N} χ_P(n) z^n  [CLOSED S134, see "Closed attacks"]

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

### D16 — Oseledec multiplicative ergodic theorem: full Lyapunov spectrum of higher-dimensional χ_P-driven cocycles

**Question:** the project's spectral measurements (E2.14 chi_P Anderson
Lyapunov, E2.18 lambda Anderson Lyapunov) compute only the TOP Lyapunov
exponent `gamma_1(E)` of an `SL(2, R)`-cocycle generated by
`T_n(E) = [[V(n) - E, -1], [1, 0]]`. The Oseledec multiplicative ergodic
theorem (Oseledec 1968) gives, for ANY ergodic cocycle on `SL(d, R)`,
a FULL Lyapunov spectrum `lambda_1 >= lambda_2 >= ... >= lambda_d` with
`Sigma lambda_i = 0` (volume preservation), together with an Oseledec
splitting `R^d = E_1 oplus ... oplus E_d` of measurable subbundles
realising each rate. **The sub-leading exponents `lambda_2, ..., lambda_d`
and the gaps `lambda_i - lambda_{i+1}` are pseudorandomness measures
distinct from the top exponent and have NEVER been computed for any
arithmetic-indicator cocycle.** Question: for `V(n) = chi_P(n)` driving
an `SL(3, R)`-cocycle (e.g., a discrete *vector* Schrödinger operator on
a width-3 ladder lattice, or a 3-step delay product
`T_n^{(3)} = T_{n+2}^{(2)}(E_1) T_{n+1}^{(2)}(E_2) T_n^{(2)}(E_3)`),
does the gap `lambda_1 - lambda_2` deviate from random matched-density
WITHOUT a W-trick correction at the same scale where `lambda_1` alone
required `W = 210` to reach noise-floor?

**Why frontier:** the Anderson-Lyapunov closures (E2.14, E2.18) measure
a SCALAR rate (top exponent only). But arithmetic structure can hide
in the GAP between `lambda_1` and `lambda_2` — the sub-leading rate
governs the SECOND principal direction in the Oseledec splitting, a
genuinely different geometric invariant. **No published work computes
the multi-Lyapunov spectrum of a chi_P- or lambda-driven d-dim cocycle
for d > 2.** If the gap `lambda_1 - lambda_2` for chi_P deviates from
random AND the W-trick cascade does NOT remove it, the deviation is
NEW arithmetic content invisible to the project's 38-measure battery.
If the gap matches random, the SL(d, R) measurement extends the W-trick
closure to a NEW (sub-leading) spectral channel, complementing E2.14
in its own right.

**Cross-domain ingredient:** Oseledec multiplicative ergodic theorem
(Oseledec 1968 *Trans. Moscow Math. Soc.* 19; Furstenberg-Kesten 1960
*Ann. Math. Statist.* 31; Ruelle 1979 *IHÉS Publ. Math.* 50; Viana
2014 *Lectures on Lyapunov Exponents* CUP). UNUSED in §5 of
CROSS_DOMAIN_TECHNIQUES (row "Multiplicative ergodic theorem
(Oseledec)"). Specifically: Bennetin-Galgani-Strelcyn QR-iteration
algorithm for computing the FULL Lyapunov spectrum (periodic Gram-
Schmidt every K steps to prevent collapse onto the top eigendirection).

**Distinction from E2.14 / E2.18 / line 422:** E2.14 / E2.18 measure
`gamma_1` only via `(1/N) log ||T_N ... T_1||`; the sub-leading
exponents are CHEAPER to compute on the SAME data via QR iteration,
but were NOT extracted. CLOSED_PATHS line 422 ("Ergodic theory / orbit
complexity — transfer operators have spectral theory = zeta zeros;
Furstenberg correspondence is circular") closes a DIFFERENT object:
the Koopman/transfer-operator spectrum of the prime-shift, where the
spectrum on `L^2(N)` is ζ-zeros via the explicit formula. The
Oseledec MET applies to the PRODUCT of d-dim transfer matrices over a
random/quasi-random potential — its Lyapunov spectrum is structurally
unrelated to the Koopman-operator spectrum of line 422 (one is a rate
of growth on R^d, the other is an L^2-spectrum on a function space).

**Concrete first step:** modify
`experiments/dynamical/anderson_localisation_chi_p/anderson_localisation_chi_p.py`
to compute the full `SL(3, R)`-cocycle Lyapunov spectrum at `N = 10^5`,
50 seeds, 31 energies. Use the standard Bennetin algorithm:
- Build a 3-step product `T_n^{(3)}(E_1, E_2, E_3) = M_3(E_3) M_2(E_2)
  M_1(E_1)` driven by a 3-stride sub-stream `(V(3n), V(3n+1), V(3n+2))`
  for `V = chi_P`. (For numerical stability, also use a width-3 vector
  Schrödinger ladder Hamiltonian `H` with V on the diagonal and
  off-diagonal coupling `-1` between rungs.)
- Initialise `Q_0 = I_3 in O(3)`. At each step: `M_n Q_{n-1} = Q_n R_n`
  via QR factorisation. Accumulate `log |R_n[i, i]|` for `i = 1, 2, 3`.
- Empirical `lambda_i = (1/N) Sigma_n log |R_n[i, i]|`. Verify
  `Sigma lambda_i = trace(log|det Pi M_n|) / N → 0` (SL(3, R)
  preservation).
- Compare `lambda_1 - lambda_2` for chi_P-driven vs Bernoulli-matched-
  density vs Rademacher across `E ∈ [-3, 3]`. Run W-trick cascade
  `W ∈ {1, 6, 30, 210}`. Z-score sweep with 51 energies and Bonferroni
  threshold ~3.16.

**Failure profile (E):** chi_P sub-leading gap `lambda_1 - lambda_2`
deviates by the same W-trick-erasable amount as `lambda_1` alone —
closes as 39th pseudorandomness measure (multi-Lyapunov spectrum)
extending E2.14 to higher-dim cocycles, but no NEW arithmetic content.
**B-grade.** **(I):** sub-leading gap for chi_P matches random within
sample noise BUT `lambda_1` does not — sub-leading rates are
pseudorandom while top rate carries HL structure. Add a B-grade
DIFFERENTIAL edge: HL singular series concentrates in the principal
direction only. **(INC):** numerical instability in QR iteration at
`N = 10^5` produces ≥ 5% systematic error; need `N = 10^6` and double-
precision Gram-Schmidt.

**A-grade success:** `lambda_1 - lambda_2` for chi_P deviates from
random by `> 5σ` AND the deviation persists at `W = 2310` (full sieve
through primes ≤ 11). Would be the FIRST arithmetic deviation of
chi_P invisible to the W-trick framework — opens a structurally new
direction in the multiplicative-regime program (§G), since the
Oseledec splitting in `R^3` resolves directions that the scalar
`SL(2, R)` measurement cannot.

**B-grade success:** clean numerical extension of the W-trick cascade
to the full `lambda_1, lambda_2, lambda_3` spectrum at `d = 3`,
documenting that EACH exponent carries HL singular-series fingerprint
removed at the same primorial threshold — promotes E2.14 from a
scalar-spectral edge to a multi-Lyapunov-spectrum edge.

**Cross-domain refs:**
- Oseledec 1968 "A multiplicative ergodic theorem. Lyapunov
  characteristic numbers for dynamical systems" *Trans. Moscow Math.
  Soc.* 19, 197
- Furstenberg-Kesten 1960 "Products of random matrices"
  *Ann. Math. Statist.* 31, 457
- Ruelle 1979 "Ergodic theory of differentiable dynamical systems"
  *Inst. Hautes Études Sci. Publ. Math.* 50, 27
- Viana 2014 *Lectures on Lyapunov Exponents* Cambridge Univ Press
- Bennetin-Galgani-Giorgilli-Strelcyn 1980 *Meccanica* 15, 9 (QR
  algorithm for Lyapunov spectrum)
- Wikipedia: Oseledets theorem
  https://en.wikipedia.org/wiki/Oseledets_theorem

**Budget:** 1-2 sessions. Session 1: implement QR Bennetin solver at
`d = 3`, `N = 10^5`; compute `lambda_1, lambda_2, lambda_3` for chi_P
+ 50-seed baselines. Session 2 (only if signal observed): scale to
`d = 4` and `W = 2310` cascade.

### D17 — Discrete Morse theory on the divisibility poset: critical-cell complexity of the prime indicator

**Question:** Forman's discrete Morse theory (Forman 1998 *Adv. Math.*
134; Forman 2002 *Sém. Lothar. Combin.* 48) attaches to any finite
regular CW complex `K` a discrete *Morse matching* `mu` on its face
poset such that the *critical cells* `Crit(K, mu)` form a sub-complex
with the same homology as `K` but typically far fewer cells. Build the
*divisibility order complex* `Delta(N) := order complex of (N_{<= N}, |)`
and assign a discrete Morse function `f(n) = log n` (or
`f([a_0 < a_1 < ... < a_k]) = k`). Question: what is the asymptotic
size `|Crit(Delta(N), mu_*)|` for the OPTIMAL Morse matching as
`N -> infty`? Does it scale as `c * pi(N)` (capturing primality count
asymptotically), as `O(log N)^c` (radical compression), or as
`Theta(N)` (no compression)? If `|Crit(Delta(N))| = O(polylog N)`,
the prime indicator IS topologically compressible at the cell-counting
level — first known polylog representation of an arithmetic structure
on the divisibility lattice.

**Why frontier:** **no published work computes the discrete-Morse
critical-cell complexity of the divisibility order complex.**
Bjorner-Wachs 2005 (*Trans. AMS* 357) studied poset homotopy types of
related arithmetic posets but never the Morse-critical-cell count.
The question is a clean combinatorial extremal problem with a
direct algorithmic translation: an explicit optimal Morse matching of
size `O(N - polylog(N))` would mean `pi(N)` is recoverable from the
critical-cell list at polylog cost. Furthermore, Benedetti-Lutz 2014
(*Exp. Math.* 23) gave a polynomial-time greedy algorithm `random
discrete Morse` that achieves near-optimal matchings on combinatorial
complexes; running it on `Delta(N)` at `N` up to `~10^4` is a single-
session experiment.

**Cross-domain ingredient:** Discrete Morse theory (Forman 2002
"A user's guide to discrete Morse theory" Sém. Lothar. Combin. 48 =
emis.de/journals/SLC/wpapers/s48forman.pdf), Benedetti-Lutz 2014
"Random discrete Morse theory and a new library of triangulations"
Exp. Math. 23, 66; Bauer-Lange 2011 "Optimal topological simplification
of discrete functions on surfaces" *Discrete Comput. Geom.* 47.
UNUSED in §4 of CROSS_DOMAIN_TECHNIQUES ("Discrete Morse theory" row).

**Distinction from D14 (cellular sheaves, S103) and D2 / E2.17 (PH on
gaps, S96):** D14 attaches an EXTERNALLY-defined sheaf with stalk
`F(n) = F_2 chi_P(n)` and computes its sheaf cohomology. D17
measures the INTRINSIC critical-cell count of the order complex
itself — no sheaf, no externally-imposed prime-indicator data; the
Morse function `log n` is a *natural* function, and the critical cells
are those of the divisibility lattice, not those of an arithmetic
sheaf. D2 used metric Vietoris-Rips on Takens-embedded gaps (a
metric-space topological invariant); D17 uses the COMBINATORIAL
divisibility-order structure (a discrete-topological invariant). The
three are categorically different objects.

**Concrete first step:** for `N ∈ {64, 256, 1024, 4096}`, build the
divisibility poset Hasse diagram via `pos = [(n, m) : n | m,
n != m, m <= N]`. Construct the order complex via
`Delta(N) := {chains in (N_{<=N}, |)}` (cell of dimension `k` is a
chain `n_0 < n_1 < ... < n_k` with `n_0 | n_1 | ... | n_k`). Implement
Forman's *random discrete Morse* algorithm (Benedetti-Lutz 2014):
greedy collapse of free face pairs via the Hasse-diagram pairing
`(sigma, tau)` with `sigma in partial tau` and no other unpaired
co-face. Iterate until no free pair remains; output critical cells.
For each `N`:
- count `m_0 = #{0-critical}`, `m_1 = #{1-critical}`, `m_2 = ...`;
- compare to `pi(N)` (number of primes ≤ N), `omega(N) = E[Omega(n)]`,
  and `log_2 N`;
- run on a CONTROL: random matched-density Bernoulli sub-poset
  (delete edges uniformly at random keeping density `2 (log log N)/log N`,
  matching divisibility-graph density);
- statistic: `m_0(prime poset) / m_0(control)` — is it < 1 (compression
  beyond random) or = 1 (Bernoulli-typical)?

Plot `log m_k / log N` vs `log N`. Does the slope match `1` (full
linear), `1/log` (logarithmic), or `0` (polylog)?

**Failure profile (E):** Morse matching greedy gives `m_0(N) ~ pi(N)`
asymptotically (every prime is forced critical because primes have
no proper divisors > 1, hence no free face to pair them with). Closes
as "primes ARE the 0-critical cells of any reasonable Morse function
on `Delta(N)`" — circular and uninformative; primes are the input,
not the output. **(I):** `m_1(N)` deviates from random Bernoulli
matched-density divisibility-poset by `< 5%` — closes as 39th
pseudorandomness measure (Morse complexity), NEW combinatorial-
topological category, no algorithmic content. **B-grade.** **(INC):**
random-Morse greedy is non-deterministic; need many runs and Bayesian
posterior on `m_*` distribution.

**A-grade success:** the OPTIMAL Morse matching has `m_1(N) =
O(polylog N)` — radical compression. Then primality of `n ≤ N` is
decidable by a polylog-size set of *cohomology pairings against the
critical cells*; the prime-counting function `pi(N) = m_0 -
(corrections)` becomes polylog-evaluable. Would be the FIRST
combinatorial-topological compression of `chi_P`. NB: this requires
the optimal matching, not just a random-greedy one — Joswig-Pfetsch
2006 *Discrete Comput. Geom.* showed optimal matchings are
NP-hard in general but the divisibility lattice has STRONG
*shellability* properties that may make `Delta(N)` polynomially
optimal-matchable.

**B-grade success:** explicit closed form for `|Crit_k(Delta(N))|` in
terms of arithmetic functions of N (e.g.,
`m_1 ~ Σ_{n ≤ N} omega(n)/2`), distinguishing prime-poset Morse
complexity from random-controlled. FIRST discrete-Morse invariant
of the divisibility poset; complements E2.17 (metric PH) and the D14
(cellular sheaf) framework with a third, structurally orthogonal
topological category.

**Cross-domain refs:**
- Forman 1998 "Morse theory for cell complexes" *Adv. Math.* 134, 90
  https://www.sciencedirect.com/science/article/pii/S0001870897916509
- Forman 2002 "A user's guide to discrete Morse theory" Sém. Lothar.
  Combin. 48
  https://www.emis.de/journals/SLC/wpapers/s48forman.pdf
- Benedetti-Lutz 2014 "Random discrete Morse theory and a new library
  of triangulations" Exp. Math. 23, 66
  https://arxiv.org/abs/1303.6422
- Bjorner-Wachs 2005 "Geometrically constructed bases for homology of
  partition lattices of type A, B and D" *Trans. AMS* 357
- Wikipedia: Discrete Morse theory
  https://en.wikipedia.org/wiki/Discrete_Morse_theory

**Budget:** 1-2 sessions. Session 1: implement Hasse-diagram
construction + Benedetti-Lutz random-greedy collapse at `N = 1024`;
empirical `m_*(N)` slope. Session 2 (only if non-trivial slope
observed): scale to `N = 10^4` and attempt to characterise the optimal
matching for restricted sub-posets (e.g., squarefree-only).

### D18 — Sudan-Guruswami list decoding: low-degree polynomial agreement of χ_P over F_p

**Question:** the Guruswami-Sudan 1999 list-decoding algorithm decodes
Reed-Solomon codes of dimension `k` and block length `n` up to the
Johnson bound `t = n - sqrt{kn}` errors (agreement fraction `>= sqrt{R}`
where `R = k/n`), in polynomial time, with list size `O(1/eps^2)` for
`eps`-radius decoding. View `chi_P` over `F_p` as a length-p received
word over `F_2 subset F_p`, and ASK: what is the smallest `k = k_*(p)`
such that there exists `f in F_p[X]` with `deg f <= k` and
`#{i in [0, p) : f(i) = chi_P(i)} >= sqrt{k p}` (agreement at the
Johnson bound)? If `k_*(p) = O(polylog p)`, then `f` is a polylog-
description polynomial whose values agree with `chi_P` on a `1/sqrt{log p}`
fraction of inputs — a polylog *sub-correlation* witness for primality.
Combined with a *list-checker* (Sudan's amplification trick), this
gives a polylog *primality predictor* with the matching probability `>
1/2 + 1/sqrt{log p}` — a positive A-grade signal not previously seen
in the project's pseudorandomness battery.

**Why frontier:** the project's CLOSED_PATHS line 474 closure
("Goldreich-Levin / list decoding for pi(x)") closed the *Hadamard*
list-decoding path: chi_P over F_2 has zero correlation with low-bit
Fourier characters and no linear encoding compresses it. The
Sudan-Guruswami algorithm is a STRUCTURALLY DIFFERENT list decoder for
*algebraic* (Reed-Solomon) codes — it tests `f(i) = chi_P(i)` for
*polynomial* `f` over `F_p`, not Hadamard-character-correlation over
`F_2`. **No published work runs Sudan-Guruswami on chi_P over F_p.**
The two failure conditions are independent: line 474 says no Hadamard
linear encoding; D18 asks whether NON-LINEAR low-degree-polynomial
agreement exists.

**Cross-domain ingredient:** Sudan 1997 *J. Complex.* 13 + Guruswami-
Sudan 1999 *IEEE Trans. Inf. Theory* 45 list-decoding algorithms;
Reed-Solomon codes; bivariate-polynomial interpolation and factoring
(Sudan 1996, Roth-Ruckenstein 2000 efficient factoring). UNUSED in §6
of CROSS_DOMAIN_TECHNIQUES ("Reed-Solomon decoding", "List decoding"
rows). Distinct from Goldreich-Levin (line 474), which is Hadamard
list-decoding over F_2.

**Distinction from line 474 closure:** line 474 closes Goldreich-Levin
on `pi(x)` (Hadamard list-decoding over F_2 of the prime-counting
function viewed as a Fourier object): no low-bit Fourier mass.
**D18 attacks chi_P (the prime indicator, not pi)** as Reed-Solomon
codeword over `F_p` (algebraic code over a prime-power alphabet, not
Hadamard linear): no low-degree POLYNOMIAL agreement?  These are
disjoint code classes (linear-in-F_2 vs polynomial-in-F_p) and
disjoint targets (cumulative-counter vs indicator). Concrete: the
Sudan-Guruswami test for `k = log p` on `p = 10^4` would not be
implied by line 474's closure of Goldreich-Levin on `pi(x)`.

**Concrete first step:** for `p ∈ {1009, 10007, 100003}`, build
`r := (chi_P(0), chi_P(1), ..., chi_P(p-1)) in F_p^p` (treating
chi_P as 0/1 lifted to F_p). Sweep target degrees
`k ∈ {2, 4, 8, 16, log_2 p}`. Run Sudan's interpolation-and-factor
algorithm:
- (interpolation) construct a non-zero bivariate `Q(X, Y) in F_p[X, Y]`
  with `(1, k)`-weighted degree at most `D = sqrt{2 (k-1) p}` (Sudan's
  bound) such that `Q(i, chi_P(i)) = 0` for all `i in [0, p)` —
  reduces to a sparse linear-system solve in `F_p^{p}`-dim via CRT or
  direct `O(p^omega)` Gaussian elimination at `p ~ 10^3`;
- (factoring) compute irreducible factors `Y - f_j(X)` of `Q(X, Y)` in
  `F_p[X][Y]` via `Roth-Ruckenstein 2000`-style root-finding;
- (verification) for each `f_j` of degree `<= k`, count
  `agree(f_j) = #{i : f_j(i) = chi_P(i)}` and report
  `agree(f_j) / p`. The Johnson bound predicts decoding succeeds iff
  `agree(f_j) > sqrt{k p}`.

Compare empirical `(k_*, agree(f_*))` against:
(a) ensemble of 100 random F_p-vectors of matched cardinality `pi(p)`
    (random subsets);
(b) Liouville-indicator `1[lambda(n) = -1]` of matched density;
(c) HL-singular-series-corrected random control.

**Failure profile (E):** for every `k = O(polylog p)`, all `f_j` of
degree `<= k` have `agree(f_j) ~ k/p + O(1/sqrt p)` (no better than
random) — closes as 39th pseudorandomness measure (algebraic
list-decoding agreement), NEW kind because it tests polynomial-
agreement structure not captured by Hadamard / Fourier list decoding.
**B-grade.**
**(I):** the Sudan-Guruswami algorithm finds `f_j` with `agree > k/p`,
but the EXCESS agreement is of order `(log p)^{-c}` for some `c <
1/2` — exactly the Sato-Tate / random matched expectation. Reduces to
known equidistribution.
**(INC):** at `p = 10^5`, the bivariate interpolation step requires
`O(p^2)` memory for the linear system — borderline single-core.

**A-grade success:** there exists `k* = O(polylog p)` and `f_* in F_p[X]`
with `deg f_* <= k_*` and `agree(f_*) >= 1/2 + 1/sqrt{log p}` —
non-trivial *positive* polynomial correlator for chi_P. Combined with
Sudan's amplification (run Sudan-Guruswami on `chi_P XOR f_*` to peel
the next correlator), gives a polylog-description ensemble of partial
predictors. First positive A-grade signal of structural compressibility
for chi_P at the polynomial-agreement level — would *open* the polylog
frontier rather than close another dictionary.

**B-grade success:** rigorous lower bound `k_*(p) >= Omega(p / log^c p)`
for any polynomial agreeing with chi_P at the Johnson rate — closes
the algebraic list-decoding angle quantitatively, complementing the
Goldreich-Levin closure (line 474) with a Reed-Solomon-flavoured
negative-shape edge.

**Cross-domain refs:**
- Sudan 1997 "Decoding of Reed Solomon codes beyond the error-
  correction bound" J. Complex. 13, 180
  https://people.csail.mit.edu/madhu/papers/1997/ldc-final.pdf
- Guruswami-Sudan 1999 "Improved decoding of Reed-Solomon and
  algebraic-geometry codes" *IEEE Trans. Inf. Theory* 45, 1757
  https://www.cs.cmu.edu/~venkatg/pubs/papers/listdecoding-NOW.pdf
- Roth-Ruckenstein 2000 "Efficient decoding of Reed-Solomon codes
  beyond half the minimum distance" IEEE Trans. Inf. Theory 46
- Guruswami 2004 textbook *List Decoding of Error-Correcting Codes*
  Springer LNCS 3282
- Wikipedia: Guruswami-Sudan list decoding algorithm
  https://en.wikipedia.org/wiki/Guruswami%E2%80%93Sudan_list_decoding_algorithm

**Budget:** 1-2 sessions. Session 1: implement Sudan-Guruswami at
`p = 1009`, `k <= 16`; report empirical `(k_*, agree(f_*))` and
compare to random baselines. Session 2 (only if non-trivial agreement
observed): scale to `p = 10^4-10^5` with sparse-interpolation.

### D19 — Mapper algorithm: arithmetic-lens topological network of χ_P

**Question:** Singh-Memoli-Carlsson 2007 introduced the *Mapper*
algorithm: given a point cloud `X subset R^d`, a *lens* function
`f: X -> Z` (typically `Z = R` or `R^2`), an open cover `U` of `Z`,
and a clusterer, the Mapper output is the *nerve* `N(f^{-1}(U))` —
a graph (or simplicial complex) whose vertices are clusters in
preimages of cover elements and whose edges connect overlapping
clusters. Mapper has detected non-obvious global structure in
high-dim biological data (Lum et al. 2013 *Sci. Rep.* 3) where
metric Vietoris-Rips persistent homology is uninformative.
Question: applied to a point cloud of consecutive prime gaps near
`p ~ 10^6` with an *arithmetic-tailored lens* (e.g.,
`f(p_n) = (p_n mod 30, p_n mod 7)` or `f(p_n) = (Omega(p_n - 1) mod 4,
log p_n)`), does the Mapper graph have a topological feature (cycle of
prescribed length, branch-point of high degree, persistent component)
that random matched-density gap sequences lack?

**Why frontier:** D2 (S96, E2.17) closed the *metric* Vietoris-Rips
PH on Takens-embedded gaps as a structural HL-detection measure with
no algorithmic opening. **Mapper is a structurally distinct
topological probe**: it does NOT use a metric on the point cloud —
the topology is induced by the LENS function `f` and a clustering on
the preimages. The arithmetic content lives in `f`'s definition (e.g.,
mod-q residues, factorisation invariants), not in pairwise distances.
**No published work applies Mapper to χ_P with an arithmetic lens.**
Two outcomes are informative: (i) Mapper graph has a feature reproducing
HL singular-series structure (similar to E2.17, but in graph-topological
language) — closes Mapper as 40th pseudorandomness measure; (ii)
Mapper detects a feature that PH did not (a non-metric-detectable
arithmetic signature) — A-grade structural fact extending the
HL-detection family with a NEW topological invariant.

**Cross-domain ingredient:** Mapper algorithm (Singh-Memoli-Carlsson
2007 *Eurographics Symp. Point-Based Graphics* = Stanford TR 2007;
Lum et al. 2013 "Extracting insights from the shape of complex data
using topology" *Sci. Rep.* 3, 1236; KeplerMapper Python library
[Saul-van Veen]). UNUSED in §4 of CROSS_DOMAIN_TECHNIQUES ("Mapper
algorithm" row, marked "candidate"). Distinct from the metric-PH
closure E2.17.

**Distinction from D2 / E2.17 (S96):** D2 used metric Vietoris-Rips
filtration on Takens delay-embedded gap vectors `y_n = (x_n, x_{n+1},
x_{n+2})` — purely metric. D19 uses an ARITHMETIC LENS function
`f: gap-space -> R^2` (e.g., `(p_n mod 30, log p_n)`) and clusters
preimages of `f`-cover elements — graph topology induced by the lens,
not by Euclidean distance. The two invariants would coincide only
if the lens factored through the Euclidean embedding, which arithmetic
lenses do not.

**Concrete first step:** for a window of `M = 5000` consecutive primes
near `p_0 ~ 10^6`, build the point cloud `X = {(p_n - p_{n-1},
p_{n+1} - p_n, p_{n+2} - p_{n+1}) : n in window}` (3-dim gap-vector).
Choose lens `f(x) = (x[0] mod 6, log mean(x))` (arithmetic-residue +
density). Cover `Z = [0, 6) x [0.5, 5]` with overlapping rectangles
(say `5 x 5` grid, 30% overlap). Cluster each preimage by single-link
agglomerative clustering (cutoff at the within-cover-cell median
inter-point distance). Build the Mapper graph `G(prime)`. Compute:
- number of vertices, edges, connected components;
- first Betti number (number of independent cycles);
- max-degree vertex;
- length-distribution of induced cycles.

Compare to:
(a) `K = 20` IID Exp(1) Cramér-normalised gap controls (B1 from S96);
(b) `K = 20` random permutations of the same window (B2);
(c) `K = 20` density-matched HL-uniform-residue baselines (B3 — gaps
    drawn from the HL-predicted joint distribution mod 6 then resampled).

Statistic: `z(b1)`, `z(b2)`, `z(b3)` for first Betti number `b1`,
component count, max-degree. Does `G(prime)` carry a topological
feature absent in B1 / B2 / B3 at `> 3 sigma`?

**Failure profile (E):** Mapper graph statistics for primes match B2
(gap-permuted) within sample noise — closes as 39th pseudorandomness
measure (arithmetic-lens topological-network), confirming that the
arithmetic-lens projection of χ_P inherits the same HL structure
already seen in metric PH (E2.17). **B-grade.**
**(I):** Mapper graph for primes deviates from B1 (Poisson) by the
HL-residue concentration but matches B2 (gap-permuted, same marginal)
— same structural picture as E2.17. Reduces to E2.17 in graph
language. **(INC):** Mapper output is highly sensitive to cover
parameters; need ablation across (5x5, 7x7, 10x10) covers and
(20%, 30%, 50%) overlaps, with bootstrap stability scoring.

**A-grade success:** Mapper graph for primes has a topological feature
(specific cycle length, branch-point degree, component count) that
deviates from B1 / B2 / B3 by `> 5 sigma` AND has a closed-form
arithmetic interpretation NOT reducible to HL pair correlation —
first arithmetic signature of χ_P invisible to metric PH (E2.17),
opening a distinct topological-network category.

**B-grade success:** Mapper graph for primes statistically distinguishable
from B2 in graph-Betti / component count by `> 3 sigma` across
multiple windows AND multiple lens choices, with the deviation
quantitatively matched by an HL-residue calculation. Adds a sixth
HL-detection family member (after E2.13 Gowers, E2.14 Anderson, E2.15
alg. immunity, E2.16 DPP, E2.17 PH, E2.19 subword) in a NEW
graph-topological language.

**Cross-domain refs:**
- Singh-Memoli-Carlsson 2007 "Topological methods for the analysis of
  high dimensional data sets and 3D object recognition" *Eurographics
  Symp. Point-Based Graphics* (Stanford TR 2007)
  https://research.math.osu.edu/tgda/mapperPBG.pdf
- Lum et al. 2013 "Extracting insights from the shape of complex data
  using topology" *Sci. Rep.* 3, 1236
  https://www.nature.com/articles/srep01236
- Saul-van Veen "KeplerMapper" Python library
  https://kepler-mapper.scikit-tda.org/
- Wikipedia: Mapper (topological data analysis)
  https://en.wikipedia.org/wiki/Mapper_(topological_data_analysis)

**Budget:** 1-2 sessions. Session 1: install KeplerMapper or implement
Mapper from scratch, run on M=5000 prime-gap window with one lens
choice; report graph-statistic deviations vs B1/B2/B3. Session 2:
ablate lens / cover / overlap; second window for robustness.

### D20 — Friedman spectral gap of the prime-Cayley graph: λ_2 vs Ramanujan bound  [CLOSED S125, see "Closed attacks"]

**Outcome (S125, mode E, B-grade case (i)):** r_N(prime) of
`Cay(Z/NZ, ±primes < N^c)` is Friedman-typical within ±2σ once
support `[2, M)` and parity (odd) are matched. Bare r_N = 2.05–11.30
ranges sub-Ramanujan by orders of magnitude, but the deviation
reduces to two trivial finite-N effects: bounded-support FFT spike
at `k = 1` (`Z_supp = +0.7..+1.9` after support-matching, within
±2σ on 10/10 cells) and parity-frequency spike at `k ≈ (N-1)/2`
modulated by the single even prime `p = 2` (`Z_no2 = +0.5..+2.1`
after both controls, within ±2σ on 10/10 cells). No HL singular-
series residual detected. Adds **EDGE E7.16**, FIRST abelian-Cayley-
spectral measurement of the prime exponential sum's Friedman /
Ramanujan ratio. CROSS_DOMAIN_TECHNIQUES §1 row "Random regular
graph spectral gap (Friedman)" promoted PROPOSED → USED-E. See
`experiments/algebraic/friedman_ramanujan_prime_cayley/` and Closed
attacks below.

**Question:** for `N` prime, build the abelian Cayley graph
`G_N = Cay(Z/NZ, S_N)` with symmetric generator set
`S_N = {±p : p prime, p < N^c}` for `c ∈ {1/3, 1/2, 1}`. The graph has
degree `d = |S_N| ~ 2 N^c / (c log N)`. Friedman 2008 (Memoirs AMS 195)
proved that a UNIFORMLY RANDOM `d`-regular graph on `N` vertices has
second-largest eigenvalue `λ_2 ≤ 2√(d-1) + ε` w.h.p. (the Alon-Boppana-
Friedman bound; graphs achieving equality are *Ramanujan*). Question:
does the *arithmetic-prime* Cayley graph `G_N` have `λ_2` (a) matching
the Friedman bound `2√(d-1)` within `o(√d)` (Ramanujan-typical, primes
generate a near-optimal expander), (b) strictly *below* the bound by
`> Cω(N) √log d` (super-Ramanujan, primes generate a *better* expander
than random), or (c) strictly *above* the bound (sub-Ramanujan,
arithmetic structure produces a *worse* expander than random)?

**Why frontier:** the abelian Cayley graph `Cay(Z/NZ, S)` is
*diagonalised by characters*: `λ_χ = Σ_{s ∈ S} χ(s) = Σ_{p < N^c}
2 cos(2π pk/N)` for the character `χ_k(x) = e^{2πi kx/N}`. So
`λ_2(G_N) = max_{k ≠ 0} |Σ_{p < N^c} 2 cos(2π pk/N)|`. This IS a
*prime exponential sum*. Vinogradov's classical bound gives
`|Σ_{p ≤ M} e(pα)| ≪ M (log M)^A / √q` for `α` close to a rational `a/q`
with `q ≤ M^{1/3}`; this implies `λ_2(G_N) ≪ d / N^{c/3} (log)^A`
for typical `k`. **But the Friedman / Ramanujan bound `2√(d-1) ~ √d`
is much sharper than `d / N^{c/3}` when `c < 3/2` — i.e., Vinogradov is
quantitatively WEAKER than Ramanujan in this regime.** No published
work (to our knowledge) establishes whether `Cay(Z/NZ, primes < N^c)`
saturates the Ramanujan bound. If yes, primes are a Ramanujan generating
set and the prime exponential sum is GUE-typical at the spectral-gap
level. If primes are super-Ramanujan, that is structural arithmetic
content invisible to Vinogradov.

**Cross-domain ingredient:** Friedman 2008 (Memoirs AMS 195, "A proof
of Alon's second eigenvalue conjecture and related problems") =
arXiv:cs/0405020. Concentration-of-measure for spectral gap of random
regular graphs. UNUSED in §1 of CROSS_DOMAIN_TECHNIQUES (row "Random
regular graph spectral gap (Friedman)").

**Distinction from CLOSED line 356 (Cayley/Ihara/spectral graph) and
CLOSED §A.A3 (E7.12, S79):** Line 356 closed `Cay(Z/xZ, primes)` as
*circular* — the trivial eigenvalue `λ_0 = |S|` recovers `π(N^c)`,
which is what we are trying to compute. Section §A.A3 / E7.12 closed
the *fixed-generator* Cayley graph `Cay((Z/nZ)*, {2, 3, 5})` as a
*pointwise primality test* — the spectrum probes `ω(n)`, not primality.
**D20 is structurally different from both:** it is not asking the
trivial eigenvalue (closure 1) and not asking pointwise primality
(closure 2). It is asking whether the SECOND eigenvalue `λ_2`
saturates the random-regular-graph Friedman bound, a quantitative
statistical question about the prime exponential sum that is OPEN
in the literature.

**Concrete first step:** for `N ∈ {509, 1009, 4001, 16001, 65537}`
(prime to give a true cyclic group), `c ∈ {1/2, 2/3}`:
1. Enumerate `S_N = {p : p < N^c, p prime}` and its negation. Compute
   `d = 2 |{p < N^c}|`.
2. Compute `λ_k = Σ_{p < N^c} 2 cos(2π pk/N)` for all `k ∈ [1, N-1]` via
   FFT (`O(N log N)` total).
3. Sort `|λ_k|`; report `λ_2 = max_{k ≠ 0} |λ_k|` and the *Ramanujan
   ratio* `r_N = λ_2 / (2 √(d-1))`.
4. Plot `r_N` vs `log N` for both `c = 1/2` and `c = 2/3`. Compare to
   a control: 100 random `d`-regular graphs (configuration model)
   on `N` vertices, mean `λ_2`.
5. Statistic: `Z_N = (r_N(prime) - r_N(random)) / σ_random`. Sign tells
   sub- vs super-Ramanujan; magnitude tells deviation.

**Failure profile (E):** `r_N → 1` as `N → ∞` for both `c` values, with
`Z_N` within `±2σ` of zero across all five `N`. Closes as 39th
pseudorandomness category: *the prime set is a Ramanujan-typical
abelian generating set of `Z/NZ`*. **B-grade,** new EDGES.md entry on
spectral pseudorandomness of primes-as-generators. **(I):** `r_N → 1`
holds at `c = 1/2` but fails at `c = 2/3` (or vice versa) — exhibits a
DENSITY-DEPENDENT pseudorandomness threshold for primes; documents
*at which density* primes start behaving Ramanujan-typically.
**(INC):** N too small for asymptotic regime; need `N ~ 10^6` and
distributed FFT.

**A-grade success:** `r_N` is bounded *strictly below* 1 by
`> 5σ` across all five `N` and both `c` values, with the gap
shrinking like `(log N)^{-α}` for some explicit `α`. Means primes are
*super-Ramanujan*: a stronger expander than any random regular graph,
giving an arithmetic-content quantitative improvement on Vinogradov.
First explicit super-Ramanujan property of the prime set.

**B-grade success:** clean numerical confirmation that `r_N ≈ 1` to
within 5% across `N`, with the residual deviation tracked to a
specific arithmetic correction (e.g., a Type-I/II decomposition of
the prime exponential sum gives the leading correction term).

**Cross-domain refs:**
- Friedman 2008 "A proof of Alon's second eigenvalue conjecture and
  related problems" Memoirs AMS 195 = arXiv:cs/0405020
  https://arxiv.org/abs/cs/0405020
- Hoory-Linial-Wigderson 2006 "Expander graphs and their applications"
  Bull. AMS 43, 439 https://www.cs.huji.ac.il/~nati/PAPERS/expander_survey.pdf
- Lubotzky 2012 "Expander graphs in pure and applied mathematics"
  Bull. AMS 49, 113 https://arxiv.org/abs/1105.2389
- Vinogradov 1937 *Some problems in the analytic theory of numbers*
  (prime exponential sum bound)

**Budget:** 1 session. FFT is `O(N log N)` per `N`, trivially fast at
`N ≤ 65537`. Random-regular controls: 100 seeds via configuration model.

### D21 — Reeb graph topology of arithmetic height functions on Z

**Question:** the *Reeb graph* `R(X, f)` of a topological space `X`
with a continuous function `f : X → R` collapses each level-set
component `f^{-1}(c)` to a point and tracks how components merge /
split as `c` varies (Edelsbrunner-Harer 2010 textbook ch. VI.3). For a
discrete-time `Z`-indexed sequence `f : Z_{≤ N} → R`, the Reeb graph is
a 1-complex whose vertices are local extrema of `f` and whose edges are
monotone segments. Its topology is a *combinatorial fingerprint*
distinct from the metric persistent-homology fingerprint (E2.17, S96)
and from the discrete-Morse cell count (D17). Question: for the
arithmetic height function `f(n) = ω(n)` (number of distinct prime
factors) on `n ∈ [2, N]`, does the Reeb graph `R_N := R(Z_{≤ N}, ω)`
have a number of vertices, edges, or first-Betti-number that DEVIATES
from a matched-density Bernoulli random sequence by `> 5σ`?

**Why frontier:** Reeb graphs are *insensitive to metric distortion*
(only level-set crossings count) and so they probe a topological
property orthogonal to PH (which IS metric-sensitive on Takens
embeddings) and to discrete Morse (which is on the divisibility *poset*,
not on a `Z`-indexed *time series*). Erdős-Kac says
`ω(n) ~ Normal(log log n, log log n)` for typical `n`, suggesting the
Reeb graph of `ω` is statistically random-like. **But the local extrema
of `ω` over short windows are dictated by primes (`ω(p) = 1`) vs
highly composite numbers (`ω(2 * 3 * 5 * 7 * ...) ≈ log log n`), and
the Reeb graph encodes which extrema "merge" as the level rises.** If
the merge-tree structure deviates from random by even a constant factor
in the loop count, that is a NEW topological invariant of the prime-
factor distribution, not captured by the Erdős-Kac CLT. Reeb graphs of
arithmetic height functions on Z have NEVER been computed.

**Cross-domain ingredient:** Reeb graphs (Edelsbrunner-Harer 2010
*Computational Topology*, ch. VI.3). UNUSED in §4 of
CROSS_DOMAIN_TECHNIQUES. Specifically: the Mapper algorithm (D19) is a
Reeb-graph-LIKE invariant via a chosen lens; D21 is the EXACT Reeb graph
of an arithmetic function.

**Distinction from D2 (PH on Takens-embedded gaps, CLOSED S96), D17
(discrete Morse on divisibility poset), D19 (Mapper):** D2 / E2.17
worked on the metric Vietoris-Rips complex of a Takens embedding; the
metric scale is the parameter. D17 works on the divisibility poset
(combinatorial cells, not a time series). D19 (Mapper) works on a
chosen lens function plus a chosen cover; the cover is a parameter.
D21 works on the EXACT Reeb graph (no metric, no lens, no cover) of
an arithmetic height function on `Z`. The choice of `f` (here `ω(n)`)
is the parameter.

**Concrete first step:** for `N ∈ {1024, 4096, 16384, 65536}`:
1. Compute `f(n) = ω(n)` for `n ∈ [2, N]` via small-prime sieve.
2. Build the Reeb graph: identify local minima `m_i` and local maxima
   `M_j` (a *local minimum* at `n` means `f(n-1), f(n+1) > f(n)`).
   Vertices = extrema; edges = monotone segments between consecutive
   extrema with their level-set merge tracking.
3. Compute graph statistics: `|V|, |E|, β_0` (component count),
   `β_1 = |E| - |V| + β_0` (first Betti, the number of independent
   loops in the contour tree), and the *persistence-of-pairing*
   distribution from the merge tree.
4. CONTROL: 100 random sequences `g(n)` matched in mean and variance to
   `ω(n)` (`Erdős-Kac` Gaussian with mean and variance `log log n`).
   Compute `(|V|, |E|, β_1)(g)` for each control; report mean and σ.
5. Statistic: `Z_β1 = (β_1(ω) - mean β_1(g)) / σ`. Bonferroni for the
   four `N` values gives threshold `~3.16`.

**Failure profile (E):** `Z_β1` within `±2σ` across all `N` and
height-functions `f ∈ {ω, Ω, λ, μ²}` — closes as: arithmetic Reeb
topology IS Erdős-Kac-Gaussian-typical, 39th-or-40th pseudorandomness
category. **B-grade.** **(I):** `Z_β1` deviates only for `f = ω` and
`f = Ω` (additive functions with prime support) but not for `f = λ`
or `f = μ²` — splits the Reeb-detection family by additive vs
multiplicative arithmetic structure. **(INC):** Reeb merge tree
sensitive to ties in `ω(n)` (very common since `ω` is integer-valued);
need careful tie-breaking (lexicographic perturbation).

**A-grade success:** `Z_β1 > 5σ` consistently AND a closed-form
asymptotic for `β_1(R_N) - β_1(R_N^{control})` derived from a
Hardy-Littlewood-style singular-series calculation. Adds a NEW
arithmetic invariant: the loop count of the Reeb graph of `ω`.

**B-grade success:** explicit power-law fit
`β_1(R_N) - β_1^{control} ~ c N^{α}` with `α` distinct from the
expected diffusion-theory `α = 1/2`. First non-trivial scaling law for
a Reeb invariant on an arithmetic sequence.

**Cross-domain refs:**
- Edelsbrunner-Harer 2010 *Computational Topology: An Introduction*
  AMS, ch. VI.3 (Reeb graphs)
- Carr-Snoeyink-Axen 2003 "Computing contour trees in all dimensions"
  *Comput. Geom.* 24, 75 (efficient Reeb / merge-tree algorithms)
- Hardy-Ramanujan 1917 / Erdős-Kac 1940 (Gaussianity of `ω(n)`)
- Wikipedia: Reeb graph https://en.wikipedia.org/wiki/Reeb_graph

**Budget:** 1 session. Reeb graphs of length-`N` sequences cost
`O(N log N)`; trivially fast at `N ≤ 10^6`.

### D22 — Higher-order Hodge Laplacian spectrum on the coprimality / divisibility complex  [CLOSED S126, see "Closed attacks"]

**Question:** the *combinatorial Hodge Laplacian* `L_k = ∂_{k+1}
∂_{k+1}^* + ∂_k^* ∂_k` acts on `k`-cochains of a simplicial complex
`K` (Lim 2020 SIAM Review 62, "Hodge Laplacians on graphs" =
arXiv:1507.05379). Its kernel computes `k`-th *real Betti numbers* via
the Hodge decomposition `dim ker(L_k) = β_k(K, R)`. The graph Laplacian
is `L_0`; *higher-order* `L_k` for `k ≥ 1` carries strictly more
information than `L_0` and is invisible to the closed Cayley / GCD-
graph spectrum (CLOSED_PATHS line 356, 387, both at `k = 0`). Question:
build the *Vietoris-Rips-like coprimality flag complex*
`K_N := \{σ ⊆ [2, N] : ∀ i, j ∈ σ, gcd(i, j) = 1\}` (a `(k+1)`-simplex
is a pairwise-coprime `(k+2)`-tuple of integers). Compute the spectrum
of `L_1(K_N)` (acting on edges) and `L_2(K_N)` (acting on triangles).
Does the smallest non-zero eigenvalue of `L_1` (the Cheeger constant
analogue at the *edge* level), or the second non-zero eigenvalue of
`L_2`, deviate from a matched-density random flag complex by `> 5σ`?

**Why frontier:** the coprimality FLAG complex is structurally
distinct from the coprimality 1-skeleton (= GCD graph): the higher-
order simplices encode pairwise-coprimality of triples and quadruples,
which involves the *Möbius* function and the singular series in
non-trivial ways (e.g., `P(gcd(a,b,c) = 1) = 1/ζ(2) ⋅ 1/ζ(3) ⋅ ...`
type products). The Hodge Laplacian SPECTRUM at `k = 1, 2` has never
been computed for any arithmetic flag complex. **CLOSED line 387 (GCD
spectrum = Ramanujan sums = Meissel-Lehmer cost) addresses ONLY `L_0`**;
it does NOT touch `L_1, L_2`. If `L_1`'s spectrum encodes more than
Ramanujan sums (in particular, if it encodes triple-coprime / Möbius
data not reducible to bilinear forms), it may carry circuit-complexity
content invisible to the closed `L_0` line.

**Cross-domain ingredient:** Lim 2020 SIAM Review 62(3) "Hodge
Laplacians on graphs" = arXiv:1507.05379. UNUSED in §1 of
CROSS_DOMAIN_TECHNIQUES. The Hodge decomposition for finite simplicial
complexes (Eckmann 1944, Friedman 1996 in the discrete setting). Specific
algorithm: build the boundary matrices `B_k`, form `L_k = B_k^T B_k +
B_{k+1} B_{k+1}^T`, diagonalise.

**Distinction from CLOSED line 356 / 387 (Cayley / GCD `L_0`),
CLOSED line 423 (transfer-operator spectrum = ζ-zeros), CLOSED §A.A3
(E7.12 Cayley spectrum = `ω(n)` probe):** all four closures are
GRAPH-LEVEL (`L_0`) or operator-on-`L^2(N)` spectra. D22 is the FIRST
proposal to compute `L_k` for `k ≥ 1` on an arithmetic flag complex.
The Hodge Laplacian `L_1` at the edge level is the *edge graph
Laplacian* of the line graph of `K_N`'s 1-skeleton; its eigenvalues
encode pairwise-coprime data PLUS triple-coprime *cohomological*
data not captured by graph Laplacians.

**Concrete first step:** for `N ∈ {32, 64, 128, 256}` (clique-counting
is exponential in clique number; this caps `N`):
1. Build the coprimality FLAG complex `K_N`: each pairwise-coprime
   triple `{a, b, c}` with `gcd(a, b) = gcd(b, c) = gcd(a, c) = 1`
   contributes a 2-simplex; quadruples contribute 3-simplices.
2. Form boundary matrices `B_1` (edges → vertices), `B_2` (triangles →
   edges). Form `L_1 = B_1^T B_1 + B_2 B_2^T`.
3. Diagonalise `L_1`, sort eigenvalues. Report the smallest 50
   eigenvalues `(λ_1^{(1)}, ..., λ_{50}^{(1)})` and `β_1 =
   dim ker(L_1)`.
4. CONTROL: random Erdős-Rényi flag complex with edge probability
   matched to the coprimality-graph density `6/π² + O(log N / N)`.
   100 seeds, mean spectrum and σ.
5. Statistic: KS distance between coprimality-spectrum and matched
   ER-spectrum at `k = 1`; same at `k = 2`. Z-score for `β_1`,
   `β_2`.

**Failure profile (E):** spectra match within KS noise floor, `β_k`
within `2σ` of ER expectation — closes as 41st pseudorandomness
category: *the coprimality flag complex is Erdős-Rényi-typical at the
Hodge-Laplacian-spectrum level for `k = 1, 2`*. **B-grade.** **(I):**
spectra match BUT `β_1` deviates by `> 3σ` — the loop count of the
coprimality complex carries a Möbius-singular-series fingerprint
absent from random ER. New arithmetic Betti edge. **(INC):** flag
complex is too dense at `N ≥ 200`; need restricted complex (e.g.,
2-coloured edges, only `gcd = 1`-via-shared-prime-cofactor).

**A-grade success:** `L_1` has a *spectral gap* `λ_2^{(1)} > c > 0`
uniformly in `N`, while the matched random ER `L_1` has `λ_2 → 0`
(or the reverse). A persistent Hodge-spectral gap of the coprimality
complex is invisible to all 38 known pseudorandomness measures and
encodes a new arithmetic constraint. **A-grade IF combined with a
TC⁰-computable boundary matrix:** spectrum-based primality test using
`L_1`-eigenvalue thresholding.

**B-grade success:** explicit closed form for `β_2(K_N)` in terms of
Möbius / Hardy-Littlewood singular series — would be the FIRST
non-trivial Betti calculation for an arithmetic flag complex.

**Cross-domain refs:**
- Lim 2020 SIAM Review 62(3) "Hodge Laplacians on graphs" =
  arXiv:1507.05379 https://arxiv.org/abs/1507.05379
- Eckmann 1944 "Harmonische Funktionen und Randwertaufgaben in einem
  Komplex" Comment. Math. Helv. 17, 240
- Friedman 1996 "Computing Betti numbers via combinatorial Laplacians"
  Algorithmica 21, 331
- Horak-Jost 2013 "Spectra of combinatorial Laplace operators on
  simplicial complexes" Adv. Math. 244, 303

**Budget:** 1 session at `N = 256`; possible scale-up to `N = 512` if
`L_1` is sparse enough. Boundary-matrix construction is `O(N^2)` for
edges, `O(N^3)` for triangles via clique enumeration; diagonalisation
is `O(|edges|^3)` ≤ `O(N^6)`.

### D23 — Density Hales-Jewett: combinatorial-line density of primes in base-k digit cubes

**Question:** identify each integer `n` in `[0, k^d)` with its base-`k`
digit string `n = (a_{d-1}, ..., a_1, a_0) ∈ {0, 1, ..., k-1}^d`, an
element of the "combinatorial cube" `[k]^d`. A *combinatorial line*
in `[k]^d` is a set of `k` strings parametrised by a *wildcard
position pattern* `w ∈ {0, ..., k-1, *}^d` containing at least one
`*`: the line is `{w(t) : t ∈ {0, ..., k-1}}` where `w(t)` replaces
every `*` with `t`. The Density Hales-Jewett theorem (Furstenberg-
Katznelson 1991; quantitative by Polymath1 2010) says: any subset
`A ⊂ [k]^d` of density `> δ` contains a combinatorial line for `d`
sufficiently large. Question: for `A = primes ∩ [0, k^d)`, what is the
density of `A` ON COMBINATORIAL LINES — i.e., for a uniform random line
`L`, what fraction `f_k(d) = E_L [|A ∩ L| / k]` of its `k` integers
are prime, and how does this compare to the global density
`π(k^d) / k^d ~ 1 / (d log k)`?

**Why frontier:** combinatorial lines in `[k]^d` are *highly
structured* arithmetic progressions in disguise. A combinatorial line
in `[10]^d` with wildcards at digit positions `S ⊂ [d]` is the AP
`{n_0 + t (Σ_{i ∈ S} 10^i) : t = 0, ..., 9}`. So D23 is testing the
density of primes on *structured* APs of common difference
`Σ_{i ∈ S} 10^i` (a sum of distinct powers of 10). Hardy-Littlewood
predicts the prime density on AP `n_0 + t * d` to be
`(1 / log n) * S(d)` where `S(d) = ∏_{p | d} p / (p - 1)`. **For
`d = Σ_i 10^i` the singular series `S(d)` depends on which digit
positions are wildcards** — a non-trivial dependence that has NEVER
been measured. If the empirical `f_k(d)` deviates from the Hardy-
Littlewood prediction by `> 5σ` for some specific wildcard pattern,
that is structural arithmetic content. If it matches, D23 closes as a
HL-saturation measure on a NEW family of APs.

**Cross-domain ingredient:** Density Hales-Jewett — Furstenberg-
Katznelson 1991 J. d'Anal. Math. 57, 64; Polymath1 2010 "A new proof
of the density Hales-Jewett theorem" Annals 175 = arXiv:0910.3926.
UNUSED in §7 of CROSS_DOMAIN_TECHNIQUES.

**Distinction from D6 (Gowers `U^k` of χ_P, CLOSED S85, E2.13) and
the W-trick / B-grade Hardy-Littlewood line:** D6 / E2.13 measured
the `U^k` Gowers norm, an `L^p`-style FOURIER-analytic quantity that
detects k-AP density at FIXED common difference averaged over starting
points. D23 measures the density of primes on combinatorial-line APs
of SPECIFIC arithmetic structure (common difference = sum of distinct
powers of `k`), a STRUCTURED-AP quantity that complements the
unrestricted `U^k` measurement. The two are statistically independent
fingerprints.

**Concrete first step:** fix `k = 10`, `d ∈ {3, 4, 5, 6, 7}` (so
`k^d = 10^7` at largest, fits in memory):
1. Compute `χ_P(n)` for `n ∈ [0, 10^d)` via segment sieve.
2. Enumerate all combinatorial lines `L_w` for wildcard patterns `w`
   with exactly 1, 2, 3 wildcards. For `d = 5`, 1-wildcard lines
   number `5 ⋅ 10^4 = 5 × 10^4`; 2-wildcard `C(5,2) ⋅ 10^3 = 10^4`;
   3-wildcard `C(5,3) ⋅ 10^2 = 10^3`.
3. For each line `L_w`, compute `|primes ∩ L_w| / 10`.
4. Average over lines with the same wildcard pattern *count* (1 *, 2 *,
   3 *). Plot mean line-density `f_k^{(j)}(d)` for `j ∈ {1, 2, 3}`
   wildcards.
5. Compare to Hardy-Littlewood prediction
   `1/log(n_0) ⋅ S(Σ_{i ∈ S} 10^i)` averaged over starting points
   `n_0` and singular-series factors over wildcard positions `S`.
6. Z-score `Z = (f_k(d) - HL_pred) / σ_pred`. Bonferroni for ~30
   pattern types at `d = 5`: threshold `~3.4`.

**Failure profile (E):** `Z` within `±3σ` for all wildcard patterns
across all `d` — closes: prime density on combinatorial-line APs is
HL-typical, 39th saturation measure. **B-grade.** **(I):** `Z > 5σ`
for `j = 1` wildcards (single-digit-flip APs) but `Z ≈ 0` for `j ≥ 2`
— a *hierarchical* HL saturation, with single-digit perturbations
showing arithmetic structure invisible to Gowers norms. Adds a graded
HL-saturation hierarchy. **(INC):** memory at `d = 7` is `k^7 = 10^7
* 8 bytes = 80 MB`, manageable; at `d = 8` bumps to 800 MB — feasible
single-session.

**A-grade success:** the wildcard-position pattern `w` (which digit
positions are *) determines a specific arithmetic correction that is
SHARP — empirical density matches a closed-form prediction involving
the singular series of `Σ_i 10^i`, deriving a NEW wildcard-pattern
HL identity with theorem-level statement.

**B-grade success:** explicit measurement of `f_k^{(j)}(d)` vs the
classical HL prediction, with the residual decomposed into a sum of
explicit arithmetic corrections. Adds a refined edge in the HL-
saturation family.

**Cross-domain refs:**
- Furstenberg-Katznelson 1991 "A density version of the Hales-Jewett
  theorem" J. d'Anal. Math. 57, 64
- Polymath1 (Gowers et al.) 2010 "A new proof of the density Hales-
  Jewett theorem" Annals of Math. 175, 1283 = arXiv:0910.3926
  https://arxiv.org/abs/0910.3926
- Hales-Jewett 1963 "Regularity and positional games" Trans. AMS 106, 222
- Hardy-Littlewood 1923 "Some problems of 'Partitio Numerorum' III"
  Acta Math. 44, 1

**Budget:** 1 session. Sieve to `10^d` is fast; line enumeration is
`O(d ⋅ k^d)` (linear in the cube size).

### D24 — Eynard-Orantin topological recursion of prime correlation moments

**Question:** Eynard-Orantin 2007 (= arXiv:math-ph/0702045)
established that the genus expansion of matrix-model correlators
`W_g^{(n)}` (the connected `n`-point function at genus `g`) is
determined recursively from a *spectral curve* `Σ` (an algebraic
curve `y(x)`) via the *topological recursion* relations: each
`W_g^{(n)}` is a residue formula on `Σ` involving lower-`g` and
lower-`n` correlators. Question: define the *prime correlation
generating function* `W^{(2)}(z_1, z_2; N) = Σ_{p, q ≤ N, p ≠ q}
1/((z_1 - p)(z_2 - q))` and its higher-`n` analogues. Does the
moment expansion of `W^{(n)}` admit a topological recursion with a
specific spectral curve, e.g., `y^2 = R(x)` for some explicit
arithmetic `R(x)`? If yes, the entire prime correlation hierarchy is
*polylog-determined* by `R(x)` — collapsing the `O(x^{2/3})` cost of
prime correlations to a recursive evaluation of residues on `Σ`.

**Why frontier:** topological recursion is the universal language of
matrix-model correlators (Hermitian, β-ensemble, and *β = 1* GUE
all admit it). The Riemann ζ moments `M_k = ∫|ζ(1/2 + it)|^{2k} dt`
are conjecturally given by random-matrix-theory (Keating-Snaith 2000)
which IS a topological-recursion-driven matrix model — so the LEADING
arithmetic content of ζ moments is captured by topological recursion.
But what about PRIME correlations? Prime n-point correlations
`Σ_{p_1 < ... < p_n ≤ N} 1` are NOT matrix-model expectations a priori,
yet they may admit a *formal* topological-recursion structure with an
arithmetic spectral curve. **No published work attempts to derive a
spectral curve for the prime correlation hierarchy.**

**Cross-domain ingredient:** Eynard-Orantin topological recursion
2007 *Comm. Number Theory Phys.* 1 = arXiv:math-ph/0702045. UNUSED in
§3 of CROSS_DOMAIN_TECHNIQUES (row "Stochastic loop equations").
Specifically: derive Schwinger-Dyson / loop equations for `W^{(n)}`
of the prime ensemble; check if they admit a closed solution via a
spectral curve `y(x)` with arithmetic content.

**Distinction from CLOSED §C.C2 (Conrey-Snaith higher-order
arithmetic corrections, S123) and free-cumulant work (§3 USED I,
edge E2.1):** Conrey-Snaith give the *RMT* arithmetic correction at
fixed `n`; D24 asks whether the *full hierarchy* `W^{(n)}` has a
recursive spectral curve, which is a STRONGER statement. Free
cumulants (E2.1) computed second-order moments only; D24 requires
ALL orders simultaneously satisfying the topological recursion ansatz.

**Concrete first step:** as a one-session SCREENING test:
1. Compute `W^{(2)}(z_1, z_2; N)` and `W^{(3)}(z_1, z_2, z_3; N)` for
   `N = 10^4, 10^5, 10^6` and rational sample points `z_i ∈ [N, 2N]`.
2. Eynard-Orantin gives `W^{(3)} = ResRes operations on W^{(2)}
   times W^{(1)}`. Test: does the EMPIRICAL `W^{(3)}` agree with the
   topological-recursion prediction from `W^{(2)}` and `W^{(1)} =
   Σ_p 1/(z - p)` to within `5σ` sample error?
3. If YES, identify the SPECTRAL CURVE: Eynard-Orantin `y(x) =
   lim_{N → ∞} (1/N) W^{(1)}(N x)`. For primes, by the prime number
   theorem `W^{(1)}(N x) ~ ∫_0^{N x} 1/(log u) du / (z - u)` —
   identify if this defines an algebraic `y(x)`.
4. If NO, identify the OBSTRUCTION: prime correlations do NOT factor
   as residues on a spectral curve — closes the topological-recursion
   route.
5. CONTROL: same `W^{(n)}` for a Bernoulli matched-density random
   indicator. Test if the recursion holds for the random control.

**Failure profile (E):** `W^{(3)}` does not match the topological-
recursion prediction even for the random control — the framework is
inappropriate for non-positive-definite sample measures; closes as a
structural mismatch (matrix-model topological recursion requires
real positivity). **B-grade.** **(I):** topological recursion holds
for the RANDOM control but NOT for primes — primes deviate from the
universal recursion in a quantitative way; defines a NEW non-recursive
edge of the prime ensemble. **B+grade.** **(INC):** sample noise at
`N = 10^4` too high; need `N = 10^6` and adaptive importance sampling.

**A-grade success:** the prime correlation hierarchy SATISFIES
topological recursion with an explicit spectral curve `y^2 = R(x)`
where `R(x)` is an explicit arithmetic function (e.g., a meromorphic
function whose residues are `1/log p`). Then ALL `n`-point prime
correlations are computable in `polylog(N)` from `R(x)` alone — a
revolutionary arithmetic compression. This would falsify the
information-theoretic floor at the n-correlation level.

**B-grade success:** explicit identification of the obstruction to
topological recursion: which step of the Eynard-Orantin derivation
fails for the prime ensemble (loop equations, spectral curve algebra,
or residue extraction), and quantification of the failure size as a
function of `N`. Adds a structural negative edge: matrix-model
topological recursion does not apply to prime correlation moments.

**Cross-domain refs:**
- Eynard-Orantin 2007 "Invariants of algebraic curves and topological
  expansion" *Comm. Number Theory Phys.* 1, 347 =
  arXiv:math-ph/0702045 https://arxiv.org/abs/math-ph/0702045
- Eynard 2016 *Counting Surfaces* Birkhäuser (textbook)
- Borot-Eynard 2011 "Enumeration of maps with self-avoiding loops..."
  arXiv:1107.5028 (B-Eynard topological recursion variants)
- Keating-Snaith 2000 "Random matrix theory and ζ(1/2 + it)"
  *Comm. Math. Phys.* 214, 57

**Budget:** 1 session for the SCREENING test (steps 1-4 above);
2-3 additional sessions if signal observed.

### D25 — Stein-Tomas restriction theorem: L^p restriction-extension bound for the prime exponential sum

**Question:** the Stein-Tomas restriction theorem (Stein 1975 +
Tomas 1975) says: for a smooth measure `dσ` on a curved hypersurface
`Σ ⊂ R^d` with `(d-1)`-dimensional non-zero Gaussian curvature, the
*restriction* operator `R: f ↦ \hat{f}|_Σ` extends boundedly
`L^p(R^d) → L^2(Σ, dσ)` for all `p ≤ 2(d+1)/(d+3)`, equivalently
the *extension* operator `E: g ↦ \widehat{g · dσ}` is bounded
`L^2(Σ) → L^q(R^d)` for `q ≥ 2(d+1)/(d-1)`. The Stein-Tomas
exponent is the SHARP `L^2 → L^p` end-point of restriction theory,
**stronger** than Bourgain-Demeter-Guth decoupling (decoupling
implies restriction at the same exponent up to `δ^{-ε}` losses).
Question: viewing the prime set as a discrete subset of the unit
circle via `\{p/N : p prime, p ≤ N\}`, define the prime extension
operator `E_N g(x) := (1/√π(N)) Σ_{p ≤ N} g(p) e^{2π i p x / N}`.
What is the *empirical* `L^p(R/Z)` norm of `E_N 1` (the unit
constant on primes), and how does it scale with `N`? If
`||E_N 1||_{L^p}` strictly beats the Stein-Tomas exponent (i.e., is
*smaller* than the random-matched-density baseline), it gives a
NEW restriction-style bound for prime exponential sums, **strictly
sharper than Vinogradov-Korobov**, with implications for the explicit-
formula error term.

**Why frontier:** the restriction theorem is the *sharpest*
canonical Fourier-analytic measure of "how independent are the
phases of a discrete set?". The published prime-side restriction
work is fragmentary (Bourgain 1989 for arithmetic progressions;
Salem-set-like work; Roth's theorem in the primes via restriction
— Green 2005 *GAFA*) but **no published work directly tests χ_P's
saturation of the Stein-Tomas restriction inequality at finite
`N`**. The Stein-Tomas test is the natural Fourier-restriction
upgrade of D15 (BDG decoupling, PROPOSED): decoupling controls
*L^p of cap-restricted f*, restriction controls the *full L^p of f*.
A strictly Stein-Tomas-saturating prime sum (exact `L^p` saturation
matching the random baseline within `(log N)^{-c}`) would be the
SHARPEST possible Fourier-analytic prime pseudorandomness statement,
strictly stronger than the BDG decoupling test of D15.

**Cross-domain ingredient:** Stein 1993 *Harmonic Analysis: Real-
Variable Methods, Orthogonality, and Oscillatory Integrals*
(Princeton); Tomas 1975 "A restriction theorem for the Fourier
transform" *Bull. AMS* 81; Tao 2003 lecture notes "Some recent
progress on the restriction conjecture"; Tao 2020 "247B Notes 1:
Restriction theory" https://terrytao.wordpress.com/2020/03/29/247b-notes-1-restriction-theory/;
Foschi-Oliveira e Silva 2024 "The endpoint Stein-Tomas inequality:
old and new" *São Paulo J. Math. Sci.* https://link.springer.com/article/10.1007/s40863-024-00422-x.
UNUSED in §7 of CROSS_DOMAIN_TECHNIQUES (row "Restriction theory").

**Distinction from D15 (BDG decoupling, PROPOSED) and CLOSED line
309 (Vinogradov circle method):** D15 measures empirical decoupling
constant `K_emp(δ) = ||f||_{L^p}/(Σ_j ||f|_{cap_j}||_{L^p}^2)^{1/2}`
— a CAP-RESTRICTED quantity at scale `δ`. D25 measures the FULL
`L^p` norm `||E_N 1||_{L^p(R/Z)}` against the Stein-Tomas global
bound, an inequality that is structurally stronger (decoupling is
weaker than restriction: BDG is needed precisely because direct
restriction fails for d ≥ 3 paraboloid). Furthermore, the prime set
is *not* a curved hypersurface — primes inhabit the integer lattice;
so the restriction question is whether prime exponential sums
saturate the SUBSET-RESTRICTION bound (Salem-set, not classical
hypersurface restriction) — the right reference is Bourgain 1989
"Bounded orthogonal systems and the Λ(p)-set problem" *Acta Math.*
162. Closed line 309 (Vinogradov circle method) closes circle-method
extraction of π(x); D25 asks the *Salem-set / Λ(p)-set* saturation
of the prime set, a structural test, not an algorithm.

**Concrete first step:** for `N ∈ {2^{12}, 2^{14}, 2^{16}, 2^{18}}`,
build the unnormalised prime exponential sum
`f_N(x) = (1/√π(N)) Σ_{p ≤ N, p prime} e^{2π i p x / N}` on the
torus `R/Z` discretised at `M = 4 N` quadrature points. Compute
`||f_N||_{L^p(R/Z)}^p = (1/M) Σ_x |f_N(x)|^p` for `p ∈ {4, 6, 8,
12}`. Compare to:
(a) the Stein-Tomas / Bourgain Λ(p)-set theoretical bound (the
   "prime version" of restriction predicts `||f_N||_{L^p} ≪
   (log N)^{c_p}` for the Λ(p)-set inequality);
(b) ensemble of 100 random matched-density Bernoulli sets
   `B ⊂ [1, N]` with `|B| = π(N)`, computing `||f_B||_{L^p}`
   for each (gives the random baseline);
(c) Liouville `f_λ_N(x) = (1/√N) Σ_n λ(n) e^{2π i n x / N}` and
   the perfect-square set baseline (the squares `{m^2 : m ≤ √N}`
   are the canonical Salem set; Vinogradov-Lerch result gives
   the L^p extremality);
(d) random Walsh-Hadamard matched-density set.

Plot `log ||f_N||_{L^p}^p / log N` vs `log N` for each `p`. The
Stein-Tomas / Λ(p) saturation says slope `≤ ε` for any `ε > 0`
(polylog-bounded `L^p` norm). Slope of order `1/(p log N)` would
be sub-Λ(p) — very tight saturation. Slope ≥ 0.1 indicates
non-saturation (primes are not Λ(p)).

**Failure profile (E):** primes saturate the Λ(p) restriction
bound exactly at the random-matched-density baseline level — closes
as 40th-or-41st pseudorandomness measure (Stein-Tomas restriction),
**B-grade**. **(I):** primes are *worse* than random matched
density at restriction (slope `> 0` for some `p`) — identifies a
NEW pseudorandomness deficit of χ_P at the global Fourier-norm
level; complements existing local-correlation closures with a
global-norm one. **(INC):** discretisation at `M = 4N` is too
coarse for `p ≥ 8`; need `M = 16 N` to control Riemann-sum bias
on the high-`p` `L^p` integral.

**A-grade success:** prime restriction `L^p` norm strictly *better*
than random matched-density by a factor growing with `N` — primes
are super-Λ(p), giving a NEW restriction-theoretic improvement on
Vinogradov-Korobov bounds for prime exponential sums. The improvement
would directly shrink the explicit-formula error term and is itself
a publishable analytic-NT result.

**B-grade success:** explicit asymptotic `||f_N||_{L^p}^p =
N^{0+o(1)}` (Λ(p) saturation) at finite `N` with quantitative
constants matching one of the random baselines — first finite-N
Λ(p) measurement on the prime set. Complements D15 (BDG decoupling)
which tests cap-restricted norms; D25 tests full `L^p` norms.

**Cross-domain refs:**
- Stein 1993 *Harmonic Analysis* Princeton, Ch. IX (oscillatory
  integrals & restriction)
- Tomas 1975 "A restriction theorem for the Fourier transform"
  *Bull. AMS* 81, 477
- Bourgain 1989 "Bounded orthogonal systems and the Λ(p)-set
  problem" *Acta Math.* 162, 227 (the Λ(p) set restriction
  inequality, the right framework for *integer-lattice subsets*)
- Green 2005 "Roth's theorem in the primes" *GAFA* 15, 340 =
  arXiv:math/0302311 (restriction-style argument applied to χ_P)
- Tao 2020 "247B Notes 1: Restriction theory"
  https://terrytao.wordpress.com/2020/03/29/247b-notes-1-restriction-theory/
- Foschi-Oliveira e Silva 2024 "The endpoint Stein-Tomas inequality:
  old and new" *São Paulo J. Math. Sci.*
  https://link.springer.com/article/10.1007/s40863-024-00422-x
- Wikipedia: Vinogradov's theorem
  https://en.wikipedia.org/wiki/Vinogradov%27s_theorem

**Budget:** 1-2 sessions. FFT-based `||f_N||_{L^p}^p` evaluation
at `N = 2^{18}`, `M = 2^{20}`, `p ≤ 12` is `O(M log M · |p|) ≈
10^{8}` ops — single-core minutes. 100-control ensemble: ~3 hours.

### D26 — Locally testable codes: constant-query primality test from a structured encoding of χ_P

**Question:** a *locally testable code* (LTC) is an error-correcting
code admitting a probabilistic tester that queries `q = O(1)` bits of
a candidate codeword `w` and accepts (resp. rejects) `w` w.h.p. iff `w`
is at small (resp. large) Hamming distance from the code. Hadamard
codes, Reed-Muller codes, and (in a deeper construction) Dinur 2007's
amplified PCP codes are all LTCs with `q = 3` queries. Question:
construct an LTC `C_N ⊂ {0, 1}^{D(N)}` (where `D(N)` is the block
length) that ENCODES the prime indicator `χ_P|_{[1, N]}` as one
specific codeword `w_P ∈ C_N`, with the LTC's local tester `T`
satisfying:
(a) on input `w_P` (the prime-encoding codeword), the tester accepts
   with probability 1 (`χ_P` is in the code);
(b) on input `w_q := w_P ⊕ e_n` (corruption at coordinate
   corresponding to integer `n`), the tester rejects iff `n ∉ P`.
If both (a) and (b) hold for some explicit LTC and tester, the
local tester IS a `q = O(1)` query primality predictor at any
given `n` — a polylog (specifically constant-query) primality test
in a NEW computational model: the *PCP query model*. This is
strictly stronger than D8/A1 (TC^0 circuit search) which targets
gate complexity at fixed input length, and orthogonal to D18
(Reed-Solomon list decoding) which tests *list-decoding agreement*
not local-test acceptance.

**Why frontier:** the LTC framework is the canonical PCP-style
tool for "constant-query verifiability." Goldreich-Sudan 2002
*J. ACM* showed LTCs of almost-linear length exist; Dinur 2007
*J. ACM* gave the gap-amplification LTC. **No published work
constructs an LTC encoding χ_P with a primality-detecting tester.**
The right framing: encode `χ_P` as a Hadamard or Reed-Muller
codeword (or as a polynomial low-degree extension over `F_p`); the
LTC tester is then a 3-query linearity / low-degree test
(BLR-test, Kaufman-Litsyn). The question is whether the linearity
defect of `χ_P`-encoded codewords carries primality information.
Two outcomes: (a) the BLR / Kaufman-Litsyn tester accepts on `w_P`
with high probability and rejects on `w_q` (the corruption at `n`)
with probability `≥ 1/2 - 1/log N` iff `n ∉ P` — A-grade, gives
constant-query primality predictor; (b) the tester is uninformative
on χ_P-encoded codewords — closes LTC primality testing.

**Cross-domain ingredient:** locally testable codes (Goldreich-Sudan
2002 *J. ACM* 49 = ECCC TR02-050; Dinur 2007 *J. ACM* 54 "PCP theorem
by gap amplification"; Goldreich 2010 textbook *P, NP and NP-Completeness*
Cambridge); Hadamard / Reed-Muller code constructions; BLR linearity
test (Blum-Luby-Rubinfeld 1993 *JCSS* 47); Kaufman-Litsyn 2005 FOCS
"Almost orthogonal LTCs". UNUSED in §6 of CROSS_DOMAIN_TECHNIQUES
(row "Locally testable codes").

**Distinction from D18 (Sudan-Guruswami list decoding) and CLOSED line
474 (Goldreich-Levin / Hadamard list decoding):** D18 asks
*list-decoding agreement* `agree(f) > sqrt(k p)` for `f` of degree
`k`. Line 474 closes Goldreich-Levin Hadamard list decoding of pi(x)
(NOT χ_P; the cumulative count, not the indicator). D26 asks
*local-test acceptance* — a POSITIVE-probability rejection of
corrupted codewords. Local testability is a STRONGER property than
local decodability: every LTC is a locally decodable code (Goldreich-
Sudan 2002), but the converse fails. The LTC question on χ_P has
never been asked.

**Concrete first step:** for `N ∈ {64, 256, 1024}`, encode
`χ_P|_{[0, N)} ∈ {0, 1}^N` via three different LTC constructions:
1. **Hadamard-coded χ_P**: `Had(χ_P)(y) = ⟨y, χ_P⟩ = Σ_{n ∈ [0, N)}
   y_n χ_P(n) mod 2` for `y ∈ {0, 1}^N`. Block length `2^N`,
   prohibitively large for `N ≥ 32`. SUBSAMPLED Hadamard:
   `Had_sampled(χ_P)` defined on a random `K`-element subset of
   `{0, 1}^N` for `K = N^2`. BLR linearity tester: pick random
   `y_1, y_2`; accept if `Had(y_1) ⊕ Had(y_2) = Had(y_1 ⊕ y_2)`.
   For `w_P` (the true encoding), tester accepts with probability 1.
   For `w_n := w_P ⊕ e_{δ_n}` (corruption of one Hadamard coordinate
   indexing the standard basis vector at `n`), tester rejects with
   some probability `p_{rej}(n)`. Compute `p_{rej}(n)` for primes
   vs composites `n`. Statistic: AUC of `(p_{rej}(n))_{n}` for
   discriminating primes from composites.
2. **Reed-Muller-coded χ_P**: lift `χ_P` to a low-degree polynomial
   `f_P ∈ F_p[X_1, ..., X_d]` with `d = O(log N / log p)` and degree
   `O(log N)`. RM tester (Kaufman-Litsyn): pick a random affine line
   `ℓ ⊂ F_p^d`; check `f_P|_ℓ` is a polynomial of degree `≤ d` (by
   univariate-polynomial fit of `p` evaluations). Compute rejection
   probability on corrupted vs uncorrupted codewords; primality of
   the corruption coordinate is the discriminator.
3. **Long-code χ_P**: each coordinate of the LongCode at length
   `2^{2^N}` is indexed by a Boolean function `f: {0, 1}^N → {0,
   1}`; the codeword `LC(χ_P)` has `LC(χ_P)[f] = f(χ_P)`. Tester:
   query `LC(χ_P)[f]` and `LC(χ_P)[g]` for two random `f, g` and
   check consistency.

For each LTC, compute primality discrimination AUC of `p_{rej}(n)`
across `n ∈ [0, N)`. Compare to (a) random-matched density baseline,
(b) Liouville `λ`-encoded LTC, (c) Bernoulli `1/log N` density
control.

**Failure profile (E):** the LTC tester acceptance probability on
`w_q` is *uniform* in `n` (independent of primality of the corruption
location) — closes as "LTC local-test rejection probability does
not encode primality information," 41st pseudorandomness measure
in NEW (PCP-style local-test) category. **B-grade.** **(I):** the
Hadamard linearity test rejects ALL χ_P-corrupted codewords with
probability close to 1 — `χ_P` is far from any Hadamard codeword
(consistent with χ_P having full algebraic immunity, S92 / E2.15);
the LTC framework cannot be the right code for χ_P. **(INC):**
subsampling Hadamard at `K = N^2` produces too noisy a discrimination
statistic; need `K = N^3`.

**A-grade success:** an LTC `C` with constant-query tester `T`
such that `Pr[T accepts w_P ⊕ e_{δ_n}] > 1/2` iff `n` is prime,
with discrimination gap `≥ 1/log N` and tester `q = O(1)`. **A
constant-query primality predictor on encoded χ_P** — first known
LTC-side primality test, complementary to AKS (deterministic but
polylog DEPTH) and BPSW (probabilistic but no LTC structure). The
implied PCP-style primality test would be an entirely new
computational model for primality.

**B-grade success:** explicit lower bound `Pr[T rejects w_P ⊕ e_n]
≤ 1/2 + o(1/log N)` for every constant-query tester on a class
of natural LTCs (Hadamard, Reed-Muller, Long-code) — strong
negative-shape edge "no constant-query LTC tester predicts
primality at non-trivial advantage." Refines E5.x (S84 SAT search,
S89 calibrated D2 circuits) with a structurally distinct
constant-query lower bound.

**Cross-domain refs:**
- Goldreich-Sudan 2002 "Locally testable codes and PCPs of
  almost-linear length" *J. ACM* 53 (2006) = ECCC TR02-050
  https://eccc.weizmann.ac.il/eccc-reports/2002/TR02-050/
- Dinur 2007 "The PCP theorem by gap amplification" *J. ACM* 54(3),
  Article 12
- Blum-Luby-Rubinfeld 1993 "Self-testing/correcting with applications
  to numerical problems" *J. Comput. Syst. Sci.* 47
- Kaufman-Litsyn 2005 "Almost orthogonal linear codes are locally
  testable" *FOCS*
- Goldreich 2010 *P, NP, and NP-Completeness* Cambridge Univ Press
- Wikipedia: Locally testable code
  https://en.wikipedia.org/wiki/Locally_testable_code

**Budget:** 1-2 sessions. Session 1: Hadamard subsampled tester at
`N ≤ 256`, `K = N^2`; report `p_{rej}(n)` AUC vs primality. Session
2 (only if signal): Reed-Muller low-degree extension over `F_p`
with `d = O(log N / log p)` degree, BLR tester at `N = 1024`.

### D27 — Newman / Erdős L^∞-flatness of the χ_P-coefficient polynomial on |z|=1  [CLOSED S138, see "Closed attacks"]

**Outcome (S138, mode I, B-grade case (c)):** First quantitative L^∞ measurement
of the prime-indicator polynomial. Empirical identity verified at all 5 N ∈
{2^{10}, 2^{12}, 2^{14}, 2^{16}, 2^{18}}: `R(a/q; primes) = √π(N) · μ(q)²/φ(q) +
O(√log N)` at every rational a/q with (a, q) = 1 and 1 ≤ q ≤ 64. The bare
‖f_N‖_∞ = (1+o(1))·√π(N), attained at z = -1 (parity major arc). Squarefree q:
relative error 0.1–6%; non-squarefree q (μ²=0): R ≈ 0.05–0.7 (Salem-Zygmund noise).
**Erdős-extremal flatness DECISIVELY REFUTED**; primes are at the OPPOSITE
extreme — saturating the trivial upper bound. Adds **EDGE E2.21**. See "Closed
attacks" below.

**Question:** Erdős 1957 conjectured that for any Littlewood polynomial
`g(z) = Σ_{n=0}^{N-1} ε_n z^n` with `ε_n ∈ {±1}`, `max_{|z|=1} |g(z)|
≥ (1 + c) √N` for some absolute `c > 0`. The corresponding question
for **Newman polynomials** (Newman 1976 *PAMS* 51) — `f(z) =
Σ_{n=0}^{N-1} a_n z^n` with `a_n ∈ {0, 1}` — concerns the L^∞
extremality of the polynomial whose coefficient sequence is a binary
indicator. Rudin-Shapiro polynomials achieve `||g||_∞ ≤ √(2N)` on the
±1 side; the analogous "flat" Newman polynomial would have `||f||_∞`
within constant factor of `√(M)` where `M = #{n : a_n = 1}` is the
weight. Define the prime-indicator Newman polynomial
`f_N(z) := Σ_{n=2}^{N} χ_P(n) z^n`, weight `M = π(N)`. **Question:**
how does `R_N := ||f_N||_∞ / √(π(N))` behave as `N → ∞`? Three
qualitative regimes:
- `R_N → 1`: primes saturate the Newman / Rudin-Shapiro-type flatness
  bound — primes are an *Erdős-extremal* binary support, encoding
  near-optimal Fourier-flatness as a structural arithmetic identity.
- `R_N → c > 1` constant: primes are flat but not extremal — fixed
  Erdős-deficit on the prime support.
- `R_N → ∞` slowly: primes induce concentration peaks in `|f_N(z)|`
  on the unit circle (Hardy-Littlewood mod-q phase coherence).

**Why frontier:** `R_N = 1` (or `R_N → 1` rapidly) would be a
**genuine Fourier-side analytic identity for the prime indicator** —
the first such identity beyond the prime number theorem. Conversely,
quantitative `R_N - 1` measures the *finite-N HL singular-series
imprint on the L^∞ norm of the prime exponential sum* — a global
analytic quantity invisible to lower-`p` `L^p` measurements. The
Erdős conjecture for ±1 polynomials was *resolved by Balister-Bollobás-
Morris-Sahasrabudhe 2020 Annals* — they showed flat ±1 polynomials
exist (constructively), settling Erdős. The 0/1 / Newman analogue
is **wide open**, and the prime-support specialization has **never
been measured at any scale**.

**Cross-domain ingredient:** Erdős / Newman / Littlewood extremal-
polynomial harmonic analysis. Bourgain 1988 *Acta Math.* 161 (Λ(p)
sets and L^∞ behaviour); Kahane 1985 *Some Random Series of
Functions* CUP (Salem-Zygmund random model); Bourgain-Lewko 2016
*J. Funct. Anal.* (mean-vs-max bounds for random Walsh-like
polynomials); Balister-Bollobás-Morris-Sahasrabudhe 2020 *Annals*
192, 977 (resolution of Erdős for ±1). UNUSED in the project: the
existing D10 (CLOSED S134, edge E2.20) measured the **Mahler
measure** `m(f_N) = exp(∫ log|f_N| dθ)` — the *log integral* of `|f|`
— while D27 measures the **L^∞ norm** = the *maximum* of `|f|`.
Mahler / L^∞ are structurally inequivalent: log-integral averages
out high peaks, L^∞ is sensitive only to peaks. (Concretely:
m(Rudin-Shapiro) = 1 trivially, but ||Rudin-Shapiro||_∞ = √(2N) is
the flatness statement.)

**Distinction from D10 (CLOSED, Mahler = log integral) and D25
(PROPOSED, Stein-Tomas Λ(p) = `||f||_p` for finite `p`):** D10
measured `log m(f_N) - log m(f_BERN)` plateau ≈ -0.307 nat. D25
measures `||f_N||_p^p` for `p ∈ {4, 6, 8, 12}`. D27 measures the
endpoint `p = ∞`, structurally distinct: `||f||_∞ = lim_{p → ∞}
||f||_p^{1/p}` differs from finite-`p` norm by factors that vanish
in `p` but are exactly the LITTLEWOOD-FLATNESS quantity. An
empirical separation between `R_N` (primes) and `R_N` (Bernoulli)
that does NOT appear in D10's Mahler measure or in D25's `L^p`
ratios would be a NEW pseudorandomness category.

**Concrete first step:** for `N ∈ {2^{14}, 2^{16}, 2^{18}, 2^{20}}`,
build `f_N(z) = Σ_{n=2}^{N} χ_P(n) z^n` and compute
`||f_N||_∞ = max_{|z|=1} |f_N(z)|` via FFT on `M = 16 N` quadrature
points (sufficient for `||f_N||_∞` accuracy `~ ||f_N||_∞ / M^2` by
oversampling theory; verify with `M = 32 N` at the smallest `N`).
Compute `R_N = ||f_N||_∞ / √(π(N))`. Compare to:
(a) 100 random Bernoulli matched-density polynomials
   `f_B(z) = Σ_{n=2}^{N} b_n z^n` with `b_n ∼ Bernoulli(π(N)/N)`
   independent — gives the typical Salem-Zygmund baseline
   `R_B ~ √(2 log N)`.
(b) Random ±1 (Littlewood) signed primes:
   `f_L(z) = Σ_{n=2}^{N} χ_P(n) ε_n z^n` with `ε_n ∼ ±1` random —
   removes the positivity bias.
(c) Rudin-Shapiro of length `π(N)` truncated to support — gives the
   flat-polynomial benchmark `R_RS = √2`.
(d) Liouville `f_λ_N(z) = Σ_{n=2}^{N} λ(n) z^n` (multiplicative
   regime control, complementing G2).
Plot `R_N` vs `log N` for each ensemble. Three falsifiers:
- `|R_N(prime) − R_B|` ≤ 0.1 across `N` → B-grade closure: 7th
  orthogonal pseudorandomness measure (L^∞-flatness category).
- `R_N(prime) → c < R_RS` significantly → A-grade: primes are
  super-Rudin-Shapiro flat, a fundamental new identity.
- `R_N(prime) >> R_B` with growth → I-mode: HL singular series
  produces L^∞ concentration peaks at rational points `z = e^{2π i k/q}`
  — partial closure with the precise concentration set as new content.

**Failure profile (E):** `R_N(prime)` matches `R_B` within sample
noise across all `N` — **B-grade**, 7th orthogonal Fourier-side
pseudorandomness measure. **(I):** `R_N(prime) > R_B` quantitatively,
and the maximum is attained at rational points encoded by primes`p`
mod small `q` — partial closure, identifies the singular-series
imprint with quantitative scale. **(INC):** FFT oversampling at
`M = 16 N` insufficient at `N = 2^{20}` — sup-norm requires
`M = 64 N` to bracket peaks; budget escalation needed.

**A-grade success:** `R_N(prime) → 1` (or `→ √2` matching Rudin-
Shapiro) with quantitative rate `R_N - 1 = O((log N)^{-c})` — primes
are an Erdős-extremal flat support; the L^∞ norm of the prime
exponential sum saturates the Newman / Rudin-Shapiro flatness bound.
This is a fundamental new analytic-NT identity, structurally
sharper than Vinogradov-Korobov bounds (which control `||f||_2`
not `||f||_∞`).

**B-grade success:** explicit `R_N - R_B = c_∞ + o(1)` with
quantitative `c_∞ ≠ 0`, complementing E2.20's Mahler-measure
plateau Δ_∞ ≈ -0.307. The L^∞-vs-Mahler comparison `R_N` and
`m(f_N)/N^{1/2}` then constitutes a TWO-PARAMETER (peak-vs-log-
mean) Fourier-side fingerprint of HL on χ_P.

**Cross-domain refs:**
- Erdős 1957 *Mich. Math. J.* 4 (the flatness conjecture for ±1
  polynomials)
- Newman 1976 "Norms of polynomials" *Proc. AMS* 51 (the 0/1
  flatness question)
- Kahane 1985 *Some Random Series of Functions* CUP, 2nd ed.
- Bourgain 1988 *Acta Math.* 161 (Λ(p) sets and L^∞ extremality)
- Bourgain-Lewko 2016 *J. Funct. Anal.* (random-polynomial mean-vs-max)
- Balister-Bollobás-Morris-Sahasrabudhe 2020 *Annals* 192, 977
  (resolution of Erdős for Littlewood polynomials)
- Wikipedia: Littlewood polynomial
  https://en.wikipedia.org/wiki/Littlewood_polynomial

**Budget:** 1 session for the SCREENING test (FFT-based `||f_N||_∞`
at `N ≤ 2^{18}` plus 100-control ensemble — `O(M log M) ≈ 10^9` ops,
single-core ~30 minutes). 2nd session if signal observed: extend
to `N = 2^{20}`, `M = 64 N`, isolate the rational-point concentration
set if `R_N > R_B`.

### D28 — Non-abelian LPS Ramanujan Cayley graph with prime-indexed quaternion generators

**Question:** Lubotzky-Phillips-Sarnak 1988 *Combinatorica* 8 construct
optimal `(q + 1)`-regular Ramanujan Cayley graphs on `PSL_2(F_p)` using
quaternion-algebra generators: for fixed prime `q ≡ 1 (mod 4)`,
Jacobi's four-square theorem gives `q + 1` integer solutions to
`a_0^2 + a_1^2 + a_2^2 + a_3^2 = q` (with parity conditions); each
solution produces a generator `g_q = (a_0 + i a_1; a_2 + i a_3 |
−a_2 + i a_3; a_0 − i a_1) \in PSL_2(F_p)` (where `i^2 ≡ −1 \pmod p`).
The resulting Cayley graph `X^{q, p}` saturates Alon-Boppana
`λ_2(X^{q, p}) ≤ 2√q` — an OPTIMAL Ramanujan expander, structurally
distinct from Friedman-typical random regular graphs (CLOSED §D.D20,
S125, edge E7.16). Question: Define a **prime-indexed LPS-merged
graph** `Y_{N, p} := Cay(PSL_2(F_p), G_N)` where the generator set
`G_N = ⋃_{q ≤ N, q ≡ 1 (4)} {g_q^{(0)}, ..., g_q^{(q)}}` is the
union of LPS generators across ALL primes `q ≤ N`. The total
generator set has cardinality `Σ_{q ≤ N, q ≡ 1 (4)} (q + 1) ≈
(N^2 / 2 \log N)` — by PNT-in-AP. Question: what is `λ_2(Y_{N, p})`?
Three regimes:
- **Super-Ramanujan**: `λ_2(Y_{N, p}) ≤ 2√(d − 1) − c` for some
  `c > 0` independent of `N` — the prime-indexed merged generator
  set is *strictly better* than the standard LPS bound. **A-grade**.
- **Ramanujan-typical**: `λ_2(Y_{N, p})` matches the standard LPS
  saturation `2√(d − 1)` within polylog. **B-grade case (i)**.
- **Sub-Ramanujan**: `λ_2(Y_{N, p}) > 2√(d − 1)` — primes induce
  spectral defects exceeding the LPS bound. **B-grade case (ii) /
  I-mode**.

**Why frontier:** the prime number theorem in AP guarantees that the
generator set `G_N` is asymptotically large, but the LPS bound
`2√(d − 1)` was proved only for *single-q* LPS graphs. Whether the
*multi-q-merged* LPS graph (one generator family per prime `q ≤ N`)
saturates the same bound is **open in the literature**. If primes are
**Ramanujan-mixing-extremal** in the multi-q sense, this is a NEW
non-abelian arithmetic spectral identity. Distinct from CLOSED
§D.D20 (S125) which addressed the abelian Cayley `Cay(Z/NZ, primes)`
spectral gap — D20 is the F_p-Fourier construction with primes as
generators; D28 is the non-abelian PSL_2(F_p) construction with
LPS-quaternion-derived generators indexed by primes. The S125
successor proposal D20.c explicitly flagged this attack.

**Cross-domain ingredient:** **Lubotzky-Phillips-Sarnak Ramanujan
graph construction** (LPS 1988 *Combinatorica* 8). Quaternion
algebras over Q ramified only at `{q, ∞}`; Strong Approximation
theorem; Deligne's bound on Hecke eigenvalues. This is the only
known *deterministic* construction of Ramanujan graphs at every
admissible degree. UNUSED in CROSS_DOMAIN_TECHNIQUES.md §1 (only
the random Friedman construction has been used, USED-E via S125).
**Channelled mathematician:** Sarnak (Ramanujan graphs, automorphic
forms, expander constructions).

**Distinction from CLOSED §D.D20 (abelian Cay(Z/NZ, primes) S125
edge E7.16) and CLOSED line 752 (E7.12 fixed generator (Z/nZ)*
Cayley):** D20 addressed `Cay(Z/NZ, ±primes < N^c)` — an abelian
Cayley graph on (Z/NZ, +) with primes as generators; closed as
Friedman-typical within ±2σ once support and parity matched.
D28 is structurally different: the underlying graph is non-abelian
(`PSL_2(F_p)`), the generators are **quaternion-algebra-derived
4-tuples** indexed by primes (not the primes themselves), and the
Ramanujan bound is the LPS bound (not the Friedman bound). Line
752 (E7.12) addressed standard fixed-generator (Z/nZ)* Cayley
spectra.

**Concrete first step:** for `p ∈ {97, 193, 389, 769}` (primes
admitting `i = √(−1)` in F_p, i.e., `p ≡ 1 (mod 4)`) and
`N ∈ {5, 13, 17, 29, 37}` (small primes ≡ 1 mod 4 for which
LPS generators exist):
1. For each `q ∈ {5, 13, ..., N}` enumerate the `q + 1` Jacobi
   four-square representations of `q` and produce the corresponding
   LPS generator matrices in `PSL_2(F_p)`.
2. Build the merged generator set `G_N = ⋃_q {LPS-q generators}`
   and compute the Cayley graph adjacency matrix on the `(p^3 − p)/2`
   vertices of `PSL_2(F_p)`. Cap `p ≤ 769` so vertex count
   `≤ 2.3 · 10^8` — feasible with sparse Lanczos for top-20 eigenvalues.
3. Compute `λ_2(Y_{N, p})` via sparse-Arnoldi on `A − (d/|V|) J`.
   Compare to LPS bound `2 √d` where `d = |G_N|`.
4. Compare to `K = 30` random matched-degree subset of LPS
   generators (across many `q`) — gives the random-LPS baseline.
5. Compare to the same construction with `G_N` replaced by
   primes-restricted to `Bernoulli(π(N)/N)` random subset of
   integers (matched cardinality, no LPS structure).

**Failure profile (E):** `λ_2(Y_{N, p}) - 2√d` matches the
random-matched-LPS-subset baseline within ±2σ across all (p, N) —
B-grade closure case (i): prime-merged LPS graphs are
Ramanujan-typical, no arithmetic-of-primes mixing-defect.
**(I):** `λ_2(Y_{N, p}) > 2 √d + c (\log N)^{-α}` with quantitative
positive deficit — primes induce a sub-Ramanujan spectral defect
in the merged construction; identifies a HL-mod-q-style imprint
on non-abelian spectra. B-grade case (ii). **(INC):** sparse-
Lanczos at `p = 769`, `|V| = 2.3 \cdot 10^8` exceeds memory; need
distributed eigensolver.

**A-grade success:** `λ_2(Y_{N, p}) < 2 √d - c (\log N)^{-α}` for
some `c > 0` and `α < 1`, robust across `(p, N)` — primes are
**super-Ramanujan-mixing-extremal** in the LPS-merged construction.
This is the first known super-LPS Ramanujan graph and an explicit
arithmetic identity for the prime spectral signature.

**B-grade success:** quantitative measurement of the prime-LPS
spectral gap `λ_2(Y_{N, p}) = 2√d + δ_N` with explicit `δ_N`
asymptotic. Even a clean B-result joins HL-detection family with
a NEW non-abelian Cayley spectral category (currently zero entries).

**Cross-domain refs:**
- Lubotzky-Phillips-Sarnak 1988 "Ramanujan graphs" *Combinatorica*
  8, 261
- Lubotzky 2012 "Expander graphs in pure and applied mathematics"
  *Bull. AMS* 49, 113 = arXiv:1105.2389 https://arxiv.org/abs/1105.2389
- Margulis 1988 "Explicit group-theoretical constructions of
  combinatorial schemes" *Probl. Inf. Transm.* 24
- Hoory-Linial-Wigderson 2006 "Expander graphs and their
  applications" *Bull. AMS* 43, 439
- Davidoff-Sarnak-Valette 2003 *Elementary Number Theory, Group
  Theory and Ramanujan Graphs* CUP
- Wikipedia: Ramanujan graph https://en.wikipedia.org/wiki/Ramanujan_graph

**Budget:** 1-2 sessions. Session 1: `p ≤ 193`, `N ≤ 17` (vertex
counts ≤ `3.6 \cdot 10^6`); single-core sparse-Lanczos minutes.
Session 2 (if signal): scale to `p = 769`, `N = 37` with sparse
Lanczos; ~hours. Successor D28.a if Y is super-Ramanujan: replace
prime-indexed LPS with squarefree-indexed and compare (isolates
prime-vs-squarefree-density component).

### D29 — Cohn-Elkies linear programming bound on the prime autocorrelation function

**Question:** the Cohn-Elkies 2003 linear programming bound for
sphere packings (*Annals of Math* 157, 689 = arXiv:math/0110009)
gives the following extremal LP: for any function `f: R^d → R` with
`f̂ ≥ 0` everywhere, `f(0) > 0`, `f̂(0) > 0`, and `f(x) ≤ 0` for
`|x| ≥ r`, the *sphere packing density* of any point set with
minimum pairwise distance `≥ r` is bounded above by
`vol(B_r) f(0) / (2^d f̂(0))`. Viazovska 2017 *Annals* 185 found
the exact magic function `f` (built from Eisenstein series and
the Jacobi theta function) achieving the sharp `E_8`-lattice bound
in dimension 8 — establishing dimension-8 sphere packing
optimality. **Question:** does the LP framework, restricted to the
1-dimensional integer lattice with the *prime set* `P = {p prime : p ≤ N}`
as the sample point set, give a sharp bound on the **prime
autocorrelation function** `R_P(t) := #{(p, q) \in P^2 : p − q = t}`?
Specifically: build the LP problem
   minimize `f(0) / f̂(0)` over functions `f: Z → R` with
   `f̂ ≥ 0` on `Z` (Bochner positivity), `f(t) ≥ R_P(t) / |P|` for
   all `t ∈ Z` (autocorrelation domination), and `f(t) ≤ 0` for
   `|t| > T_max`.
The LP optimum is automatic (numerical convex optimization);
the **question is whether primes saturate the LP bound**, i.e.,
whether `f(0) / f̂(0)` matches the actual `|P| / N` density
within polylog factors, and whether the optimal `f*` admits a
**modular-form representation** in the Viazovska sense.

**Why frontier:** if primes saturate the Cohn-Elkies LP bound for
1-dimensional integer-lattice autocorrelation, this is a
**fundamental Fourier-side extremality identity for χ_P**, equivalent
in spirit to "primes are an autocorrelation-optimal subset of Z."
The HL conjecture predicts `R_P(t) / |P|` ~ singular series
`S(t) := ∏_{p | t, p > 2} (p − 1)/(p − 2) · 2 C_2 / log^2 N` for
`t` even (twin-prime-like asymptotics). The Cohn-Elkies LP
optimum on `R_P` would be a SHARP analytic-NT bound, sharper than
HL's heuristic constants. The Viazovska-style modular-form magic
function — if it exists for the prime autocorrelation problem —
would **directly compute the singular series** as a residue formula
on a modular form, a polylog evaluable identity. **No published work
attempts a Cohn-Elkies LP analysis of arithmetic autocorrelations;**
this is a structurally new framework for the project.

**Cross-domain ingredient:** **Cohn-Elkies LP bounds for sphere
packing** (Cohn-Elkies 2003 *Annals* 157 = arXiv:math/0110009);
Viazovska 2017 *Annals* 185, 991 = arXiv:1603.04246 (E_8 magic
function); Cohn-Kumar-Miller-Radchenko-Viazovska 2017 *Annals* 185,
1017 (Leech / `Λ_{24}` sphere packing). UNUSED in
CROSS_DOMAIN_TECHNIQUES — not in any §. This is a structurally new
import (sphere-packing LP + modular forms applied to arithmetic
autocorrelation). **Channelled mathematician:** Cohn / Viazovska
(LP bounds via modular-form construction).

**Distinction from CLOSED §D.D7 (DPP fit, S95, E2.16) and PROPOSED
§D.D25 (Stein-Tomas restriction):** D7 fit a *kernel* `K(x, y)`
to the prime point process; the LP question is dual — given the
autocorrelation data, find the *extremal positive-definite Fourier
function* bounding it. D25 measures `||f_N||_p` for the prime
exponential sum `f_N(x) = Σ_p e^{2πi p x / N}`; D29 measures
*autocorrelation* `R_P(t)` directly and tests Cohn-Elkies LP
saturation. D29 is the *modular-form-side* extremality, D25 is the
*L^p Fourier* extremality — Viazovska's connection between Eisenstein
series and Cohn-Elkies LP makes D29 the natural modular-form analogue.

**Concrete first step:** for `N ∈ {10^4, 10^5, 10^6}` and
`T_max = 10^3`:
1. Compute `R_P(t) = #{(p, q) \le N : p − q = t}` for `t ∈ [0, T_max]`
   — direct convolution on `χ_P` indicator vector via FFT, `O(N \log N)`
   per `N`.
2. Set up the Cohn-Elkies LP in finite-dimensional approximation:
   parameterise `f` by `f(t) = Σ_{k = 0}^{K} c_k φ_k(t)` for
   Hermite-function basis `φ_k` (with rapid Fourier decay so
   `f̂` is computable in closed form). Constraints: `f̂(ξ) ≥ 0`
   on `Z` (discretise `ξ ∈ [0, 1)` at `M = 4 N` points;
   linear constraints in `c_k`); `f(t) ≥ R_P(t) / |P|` for
   `t ∈ [0, T_max]` (linear); `f(t) ≤ 0` for `|t| > T_max` (linear).
   Objective: `c_0 f(0) / Σ_k c_k φ̂_k(0)` (a linear-fractional
   objective; convexify by normalisation `f̂(0) = 1`).
3. Solve via a standard LP solver (`scipy.optimize.linprog` with
   `K ≤ 50` basis functions, ~thousands of constraints, ~milliseconds).
4. Compare LP-optimum bound `f*(0)` to actual `|P|/N = π(N)/N`.
   Saturation ratio `S_N := |P| / (N \cdot f*(0))`. Three regimes:
   `S_N → 1` (primes saturate within polylog) — A-grade extremal;
   `S_N → c < 1` (primes leave LP slack) — B-grade with
   quantitative gap measurement; `S_N → 1` AND `f*` admits modular-
   form representation — A++ Viazovska-style identity.
5. CONTROL: same LP on a Bernoulli-matched-density random subset
   of `[1, N]` — gives the random-baseline saturation `S_N(B)`.

**Failure profile (E):** primes saturate the Cohn-Elkies LP bound
within ±5% across `N ∈ {10^4, 10^5, 10^6}` AND the random Bernoulli
control also saturates within ±5% — the LP is too loose to
discriminate; closes as 7th-or-8th orthogonal pseudorandomness
measure with quantitative LP-saturation as the new content.
**B-grade.** **(I):** primes are *more* extremal than random
matched-density (`S_N(prime) → 1`, `S_N(B) < 1`) — primes have
genuine Cohn-Elkies LP-saturation that random subsets do not;
quantitative gap is the new content. **B+grade.** **(INC):** LP
solver feasibility breaks at `K > 50` basis functions; need
larger basis or specialized SDP relaxation for tighter bound.

**A-grade success:** primes saturate the LP bound within
`(\log N)^{-c}` AND the optimal LP function `f*` admits a
**modular-form representation** (Viazovska-style: `f*(t) =
A · θ(t) + B · E_4(τ) θ(t)` or similar Eisenstein-Jacobi-theta
combination for `t ∈ Z`) — gives an **explicit closed-form
identity** for the prime autocorrelation, polylog-evaluable from
the modular-form coefficients. This would be a **fundamental
analytic-NT identity** with direct implications for the singular
series and the level of distribution.

**B-grade success:** LP-saturation ratio `S_N(prime)` measured
across `N`, with quantitative comparison to random and to the HL-
predicted singular series. Even non-saturating, the gap
`(1 - S_N(prime))` is a NEW Cohn-Elkies-extremality fingerprint of
χ_P, structurally orthogonal to D7 (DPP), D10 (Mahler), and D25
(Stein-Tomas).

**Cross-domain refs:**
- Cohn-Elkies 2003 "New upper bounds on sphere packings I"
  *Annals of Math* 157, 689 = arXiv:math/0110009
  https://arxiv.org/abs/math/0110009
- Viazovska 2017 "The sphere packing problem in dimension 8"
  *Annals of Math* 185, 991 = arXiv:1603.04246
  https://arxiv.org/abs/1603.04246
- Cohn-Kumar-Miller-Radchenko-Viazovska 2017 "The sphere packing
  problem in dimension 24" *Annals* 185, 1017 = arXiv:1603.06518
- Cohn 2017 "A conceptual breakthrough in sphere packing"
  *Notices AMS* 64(2), 102
- Cohn-Goncalves 2019 "An optimal uncertainty principle in twelve
  dimensions via modular forms" *Invent. Math.* 217 = arXiv:1712.04438
  (LP-modular-form connection sharpened)
- Wikipedia: Sphere packing in dimension 8
  https://en.wikipedia.org/wiki/Sphere_packing
- Wikipedia: Maryna Viazovska
  https://en.wikipedia.org/wiki/Maryna_Viazovska

**Budget:** 1-2 sessions. Session 1: LP setup at `K = 30`,
`T_max = 10^3`, `N ≤ 10^5`; reports `S_N(prime)` vs `S_N(B)`.
Session 2 (if signal): scale to `K = 100`, `T_max = 10^4`,
`N = 10^6`; investigate modular-form structure of `f*` (compare
coefficients to Eisenstein series `E_4, E_6, E_8` and Jacobi theta
`θ` evaluated at suitable cusps).

### D30 — Pollicott-Ruelle resonances of an arithmetic transfer operator

**Question:** for a hyperbolic dynamical system `(M, T)` and a
Hölder-continuous weight `w: M → R`, the **transfer operator**
`L_w: f \mapsto (L_w f)(x) := Σ_{T(y) = x} w(y) f(y)` acts on a
suitable Banach space (anisotropic Sobolev space, Liverani 2004
*Comm. Math. Phys.* 248) with **discrete spectrum** outside an
essential-spectrum disc — the eigenvalues are called **Pollicott-
Ruelle resonances** (Pollicott 1985 *Inventiones* 81; Ruelle 1976
*IHÉS Publ. Math.* 50). The leading resonance gives the topological
pressure / Lyapunov exponent; sub-leading resonances `λ_k ∈ C`
with `|λ_k| < 1` give **exponential decay rates of correlations**:
`⟨f \cdot g \circ T^n⟩ - ⟨f⟩ \langle g⟩ ~ \sum_k c_k(f, g) λ_k^n`.
**Question:** define an *arithmetic transfer operator* by
weighting the Gauss map `T(x) = \{1/x\}` by a prime-related
function `w(x) := h(\lfloor 1/x \rfloor)` for `h: N → R` chosen as
either `h = χ_P` (prime indicator) or `h = λ` (Liouville function)
or `h = Λ` (von Mangoldt). What are the Pollicott-Ruelle resonances
`{λ_k}` of `L_w`? Three qualitative outcomes:
- **Trivial spectrum** (`λ_k` matches the spectrum of `L_1` for
  the unweighted Gauss map, no arithmetic resonance) — closes
  arithmetic transfer-operator route.
- **Single arithmetic resonance** `λ_* < 1` whose magnitude
  encodes a quantitative *correlation-decay rate* of arithmetic
  data (e.g., `\log |λ_*|` matches the PNT error term `\log N \cdot N^{-1/2}`).
- **Continuous spectrum / no isolated resonance** (Liouville case
  expected) — confirms the multiplicative-regime pseudorandomness.

**Why frontier:** Pollicott-Ruelle spectral theory is the canonical
**dynamical** tool for quantitative correlation-decay rates. For
arithmetic dynamical systems weighted by NT functions, the question
of whether a discrete resonance encodes π(x) information is
**OPEN**. The Mayer 1991 transfer operator approach to ζ via the
Gauss map (Mayer 1991 *Bull. AMS* 25) showed that the *unweighted*
Gauss-map transfer operator has Pollicott-Ruelle resonances at
`λ_k = ζ(2k)/k` (for even `k`), giving an **alternative entire-
function representation of ζ via a dynamical determinant** —
`ζ(s) = \det(I - L_s)` for a one-parameter family of transfer
operators. The PRIME-WEIGHTED analogue is unstudied. If the
χ_P-weighted Gauss map has a Pollicott-Ruelle resonance encoding
`Σ_n χ_P(n) f(\{1/n\})` for some test function `f`, this is a
**dynamical-determinant representation of the prime indicator** —
analogous to Mayer's representation of ζ.

**Cross-domain ingredient:** **Pollicott-Ruelle resonances /
transfer operator spectral theory** (Pollicott 1985 *Inventiones*
81, 413; Ruelle 1976 *IHÉS Publ. Math.* 50; Baladi 2018 *Dynamical
Zeta Functions and Dynamical Determinants* Springer Ergeb. 68;
Liverani 2004 *Comm. Math. Phys.* 248 anisotropic Sobolev spaces;
Tsujii 2008 *Ergod. Theory Dyn. Syst.* 28 transfer-operator
spectral gap; Gouëzel-Liverani 2006 *J. Diff. Geom.* 79). Mayer
1991 *Bull. AMS* 25 (transfer-operator representation of ζ via
Gauss map). UNUSED in CROSS_DOMAIN_TECHNIQUES §5 row "Transfer
operator spectrum" (PROPOSED only as candidate). **Channelled
mathematician:** Ruelle / Baladi (transfer-operator spectral theory),
Mayer (arithmetic transfer operators).

**Distinction from CLOSED line 105 (Transfer matrix sieve), CLOSED
line 182 (FRACTRAN transfer operator), and CLOSED line 425
(Furstenberg correspondence):** Line 105 was a CONSTRUCTIVE
transfer-matrix sieve (state space `lcm(primes ≤ √x)` exponentially
large); D30 is SPECTRAL not constructive. Line 182 was the FRACTRAN
implementation of an explicit prime-generating finite automaton;
D30 is the spectral theory of the Gauss map's transfer operator
weighted by χ_P, structurally distinct (continuous-state dynamical
system, not discrete automaton). Line 425 closed the abstract
Furstenberg correspondence; D30 is the *quantitative* Pollicott-
Ruelle spectrum of a SPECIFIC weighted dynamical system, not the
abstract correspondence principle.

**Concrete first step:** for the Gauss map `T(x) = \{1/x\}` and
weight `w_h(x) := h(\lfloor 1/x \rfloor)` with `h ∈ \{χ_P, λ, Λ_{≤N}\}`:
1. Discretise `[0, 1]` at `M = 10^4` adaptive nodes (Markov
   partition by Gauss map cylinders `[1/(n+1), 1/n]` for
   `n = 1, ..., n_{max}`).
2. Build the transfer operator matrix `L_h ∈ R^{M \times M}` by
   `(L_h)_{ij} = w_h(y) / |T'(y)|` where `y` is the unique
   pre-image of node `x_i` in the cylinder `j`.
3. Compute the top-30 eigenvalues `λ_k` of `L_h` via Arnoldi.
   The leading `|λ_0|` is the spectral radius; `|λ_1| / |λ_0|`
   is the spectral gap (correlation decay rate).
4. Compare to the unweighted Gauss-map spectrum (Mayer 1991 gives
   exact eigenvalues `λ_k = ζ(2k)/k` for `k ≥ 1`); deviations
   `δ_k = λ_k(L_h) − λ_k(L_1)` are the arithmetic-resonance
   signal.
5. CONTROL: `h = ` Bernoulli-matched-density random function on
   `\{1, ..., n_{max}\}` — gives the random-weight baseline
   spectral perturbation.
6. Increase `M, n_{max}` to verify resonance stability under
   discretisation refinement (Pollicott-Ruelle resonances are
   stable; spurious-discretisation eigenvalues are not).

**Failure profile (E):** χ_P-weighted Gauss-map transfer operator
spectrum matches the random Bernoulli-weighted baseline within ±2σ
across all `(M, n_{max})` — closes mode E with the precise
discretisation-stable spectral bound as the new content.
**B-grade.** **(I):** χ_P-weighted spectrum has a non-trivial
resonance at `|λ_*| < 1` quantitatively distinct from random
baseline by `\geq 3σ`, AND the resonance is stable under
refinement — identifies a NEW dynamical-spectrum HL-detection
fingerprint. **B+grade.** **(INC):** transfer operator
discretisation truncation error at `M = 10^4` exceeds the
arithmetic signal; need anisotropic-Sobolev-norm stabilization
(Liverani 2004) at `M = 10^5`, requiring distributed eigensolver.

**A-grade success:** isolated Pollicott-Ruelle resonance `λ_*`
of `L_{χ_P}` satisfying `|λ_*| = c · π(N) / N` for an explicit
`c > 0` (or some other polylog-evaluable arithmetic function),
giving a **dynamical-determinant representation of π(x)** in the
sense of Mayer 1991 for ζ. This would be a fundamental new
transfer-operator identity for the prime indicator,
structurally analogous to (but distinct from) the Riemann ζ
representation via the unweighted Gauss map.

**B-grade success:** explicit quantitative measurement of the
χ_P-weighted Gauss-map transfer operator spectral gap
`|λ_1| / |λ_0|`, with comparison to the unweighted Mayer
spectrum and to random-weight baselines. Even without an
isolated arithmetic resonance, the quantitative perturbation is a
NEW dynamical-spectrum measurement of χ_P pseudorandomness,
orthogonal to the existing 6+ pseudorandomness categories.

**Cross-domain refs:**
- Pollicott 1985 "On the rate of mixing of Axiom A flows"
  *Inventiones Math.* 81, 413
- Ruelle 1976 "Zeta-functions for expanding maps and Anosov flows"
  *Inventiones Math.* 34, 231
- Mayer 1991 "The thermodynamic formalism approach to Selberg's
  zeta function for PSL(2, Z)" *Bull. AMS* 25, 55
- Baladi 2018 *Dynamical Zeta Functions and Dynamical Determinants
  for Hyperbolic Maps* Springer Ergeb. 68
- Liverani 2004 "On contact Anosov flows" *Annals of Math.* 159,
  1275 (anisotropic Sobolev spaces)
- Tsujii 2010 "Quasi-compactness of transfer operators for contact
  Anosov flows" *Nonlinearity* 23, 1495
- Wikipedia: Transfer operator
  https://en.wikipedia.org/wiki/Transfer_operator

**Budget:** 1-2 sessions. Session 1: SCREENING test at
`M = 10^4`, `n_{max} = 100`, three weights `(χ_P, λ, Λ_{≤N})`
plus 30 random-weight controls. Compute top-30 eigenvalues, plot
deviation `δ_k` vs random baseline. Session 2 (if signal):
refine to `M = 10^5`, `n_{max} = 500`, with anisotropic-Sobolev
stabilization (Liverani-style cone fields) to verify resonance
stability.

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

### §D.D27 — Newman / Erdős L^∞-flatness of the χ_P-coefficient polynomial f_N(z) = Σ_{n≤N} χ_P(n) z^n on |z|=1 (CLOSED 2026-04-27, S138, mode I, B-grade case (c))

**Outcome:** B-grade case (c) mode I — first project quantitative L^∞
measurement of the prime-indicator polynomial. The L^∞ norm at every
rational point a/q is **EXACTLY the Hardy-Littlewood density factor
μ(q)²/φ(q) · √π(N)** to within Salem-Zygmund noise `O(√(2 log N))`.
The bare `‖f_N‖_∞ = (1+o(1))·√π(N)`, attained at z = -1 (q=2 parity
major arc). **Erdős-extremal flatness DECISIVELY REFUTED**: primes
are at the OPPOSITE extreme — saturating the trivial upper bound from
above via parity, NOT achieving Rudin-Shapiro `R = √2` flatness.

Cross-domain technique: **Newman / Erdős / Littlewood polynomial
L^∞-flatness extremal harmonic analysis** (Newman 1976 *Proc. AMS* 51;
Erdős 1957 *Mich. Math. J.* 4; Kahane 1985 *Some Random Series of
Functions* CUP §6 Salem-Zygmund; Bourgain 1988 *Acta Math.* 161 Λ(p)
sets; Balister-Bollobás-Morris-Sahasrabudhe 2020 *Annals* 192, 977
resolution of Erdős for ±1). PROPOSED in D27 (S136 frontier_gen).
Promotes CROSS_DOMAIN_TECHNIQUES §2 row PROPOSED → USED I with edge
**E2.21**. Channelled mathematician: **Vinogradov** (major-arc /
minor-arc decomposition for prime exponential sums, HL circle method).

Setup: `N ∈ {2^{10}, 2^{12}, 2^{14}, 2^{16}, 2^{18}}`. FFT-based
`‖f_N‖_∞ = max_k |F[k]|` on `M = 16(N+1)` oversampled DFT points;
oversampling validated stable across `oversample ∈ {8, 16, 32, 64, 128}`
to integer ‖f‖_∞ = π(N) at small N (DC saturation). Major-arc
breakdown: `R(a/q) = |F[a M / q]| / √π(N)` for every (a, q) = 1 with
1 ≤ q ≤ Q_max (Q_max ≤ 64 at N = 2^{18}, limited by FFT resolution at
small N). Minor-arc max: `max_k |F[k]|` for k NOT within ±half_arc =
±64 indices of any major arc up to Q_max. Five ensembles: PRIMES,
matched-density Bernoulli (uniform support, 50 seeds), matched-parity
Bernoulli (odd support + p=2, 50 seeds), signed-prime Littlewood
(50 seeds), Liouville `f_λ(z) = Σ λ(n) z^n`, Rudin-Shapiro at length
π(N). Total wall-time ~3 min on one core (~88s for the 50-seed
Bernoulli ensemble at N = 2^{18}).

Findings (verbatim from results.md):

> **F (a) (Salem-Zygmund-typical, B-grade) REFUTED:**
> R_N(prime) = (1+o(1))·√π(N) → ∞ at the bare and off-DC level.
>
> **F (b) (Erdős-extremal flat, A-grade) REFUTED:**
> R_N(prime) → ∞ as √π(N), the OPPOSITE extreme of `R = √2`.
>
> **F (c) (HL singular-series imprint, I-mode) HOLDS:**
> R(a/q; primes) = √π(N) · μ(q)² / φ(q) + O(√(2 log N)) for every
> (a, q) = 1 with 1 ≤ q ≤ 64. Squarefree q: 0.1–6% relative error
> (median 1.5%) across 64 q values at N=2^{18}; non-squarefree q
> (μ²=0): R ≈ 0.05–0.7, vanishes within Salem-Zygmund noise.
>
> **Q_max scan (minor-arc residual at N = 2^{18}, π(N) = 23000,
> √π(N) = 151.66):**
>
> | Q_max | R(prime, minor) | R(Bern, minor) ± std | z-score |
> |------:|----------------:|----------------------:|--------:|
> | 1     | 151.64          | 10.55 ± 0.66          | +214.7  |
> | 2     | 75.76           | 10.75 ± 0.81          | +80.2   |
> | 4     | 75.76           | 10.70 ± 0.57          | +113.9  |
> | 8     | 38.09           | 10.71 ± 0.51          | +53.4   |
> | 16    | 19.17           | 10.84 ± 0.75          | +11.1   |
> | 32    | 12.83           | 10.84 ± 0.64          | +3.1    |
> | 64    | 12.31           | 10.56 ± 0.57          | +3.1    |
>
> Scaling: at each Q_max, R(prime, minor) ≈ √π(N) / φ(q*), where
> q* = smallest squarefree integer > Q_max.

**Mechanism for closure mode I:** the L^∞ norm DECOMPOSES exactly
into the Hardy-Littlewood circle-method singular series. Major-arc
R(a/q) = √π(N)·μ²(q)/φ(q) is the classical Vinogradov / Hardy-
Littlewood exponential-sum identity `Σ_{p≤N} e^{2πi p a/q} ≈
μ(q)/φ(q)·π(N)` (a consequence of the explicit formula on residue
classes mod q). The bare ‖f‖_∞ is dominated by q=2 parity. The
minor-arc residual converges to the Bernoulli noise floor as
Q_max → ∞, modulo cumulative O(Q²) HL peaks at q > Q_max (giving
the persistent ~3σ residual at Q_max ∈ {32, 64}).

**What this rules out:** any "flat-polynomial" representation of
χ_P as a polylog primality witness. Any L^∞-Chebyshev-style
compression of f_N cannot beat √π(N) ≈ √(N/log N) storage — the
parity peak alone forces this lower bound.

**Distinct from:** D10/E2.20 (Mahler = log integral, geometric
mean) — D27 makes the parity major-arc contribution to the Jensen
integral explicit at z=-1, partially accounting for E2.20's
Δ_∞ ≈ -0.307; D25 (Stein-Tomas L^p, finite p) — D27 covers the
p=∞ endpoint; CLOSED line 778/E7.16 (Friedman/Ramanujan prime-
Cayley spectral gap) — that uses the prime SET as Cayley generators
on Z/NZ, this uses χ_P as polynomial coefficients.

**Adds EDGE E2.21**: the L^∞-norm endpoint of the L^p Fourier-side
characterisation of χ_P, complementing the existing 7th orthogonal
pseudorandomness measure category (Fourier-side L^p).

**Why this is B-grade (not C):** the test produced (i) an empirical
identity `R(a/q) = √π(N)·μ²(q)/φ(q)` verified across 124+ rationals
at 5 orders of magnitude of N — first quantitative measurement of HL
in the L^∞-norm at major arcs in the project; (ii) a Q_max-scan
demonstrating convergence to Bernoulli noise floor; (iii) a
quantitative connection to E2.20's Mahler-measure deficit via the
parity major arc Jensen contribution.

**Why not A:** R_N still grows polynomially (Θ(√π(N)) ~ √(N/log N)),
no polylog evaluator opens; the structural mechanism is the
classical HL circle method already known via E1.5 / E1.10.

**Successor entries (proposed in S138):**

(a) **D27.a — Vinogradov minor-arc supremum.** Direct numerical
verification of the Vinogradov bound: at N = 2^{18}, sample F at
10⁴ random rationals a/q with q ∈ (N^{2/5}, N^{3/5}) and (a, q) = 1.
Does max R(a/q) ≤ C · N^{1/4} (log N)^A for explicit C, A?
Cross-domain: **Vinogradov's mean-value method** (USED in §10).

(b) **D27.b — Liouville L^∞ at major arcs.** Möbius/nilsequence
orthogonality (Green-Tao 2012) predicts `R^λ(a/q) = o(1)` for ALL q
including q=2. Screening data already shows
`‖f_λ‖_∞/√N = 3.16` at N=2^{14} with argmax at freq/M=0.478 (NOT at
q=2). Quantitative single-session test: do all major-arc R^λ values
go to zero at N → ∞ (Möbius orthogonality at every rational)?

(c) **D27.c — twin-prime Newman polynomial L^∞.**
`f_N^{twin}(z) = Σ_{n: n,n+2 prime} z^n`. Predict R(a/q) follows
the Hardy-Littlewood twin-prime singular series with the
`σ_2(q) := ∏_{p|q}(p-2)/(p-1)` factor on top — different q-pattern
than D27. A-grade if it doesn't, structural fingerprint refinement
if it does.

**Files / directories:** `experiments/analytic/newman_linfty_chi_p/`
— `.py`, `_results.md`, `results.json`, `run_full.log`. See
`archive/sessions/session138_d27_newman_linfty_chi_p.md`.

### §D.D10 — Mahler measure of the prime indicator polynomial f_N(z) = Σ_{n≤N} χ_P(n) z^n (CLOSED 2026-04-27, S134, mode I, B-grade)

**Outcome:** B-grade negative-shape (mode I) — first project measurement
in the **algebraic-height / multiplicative-height** category. Pre-stated
A-grade (cyclotomic compressibility `m(f_N) = O((log N)^c)`) decisively
**falsified**; pre-stated B-grade (constant-deficit `Δ_∞ ≠ 0` from
density-matched random) **met** at z(MATCH) = `−337σ` at N = 2^{18}.

Cross-domain technique: **Mahler measure / Lehmer's conjecture / log
Weil height** (Smyth 2008; Lehmer 1933 *Ann. Math.* 34; Boyd 1981
*Canad. Math. Bull.* 24). Jensen's formula
`log m(f) = ∫₀¹ log|f(e^{2π i θ})| dθ` evaluated by FFT. UNUSED in the
project before S134; promotes CROSS_DOMAIN_TECHNIQUES §2 row PROPOSED →
USED I with edge E2.20. Channelled mathematician: **Boyd / Smyth**
(range-of-Mahler-measure conjecture).

Setup: `N ∈ {2^{10}, 2^{12}, 2^{14}, 2^{16}, 2^{17}, 2^{18}}`, `M = 2^{18}`
or `2^{19}` Jensen-FFT sample points (numerical agreement with mpmath
polyroots Jensen formula at N ∈ {64, 128} to 4 decimal places). Five
ensembles: PRIMES, BERN (iid Bernoulli matched on density `π(N)/N`,
50–100 controls), MATCH (matched cardinality `π(N)`, 50–100 controls),
LIOUVILLE (`1[λ(n) = −1]`, density 1/2), SQFREE (`μ²(n) = 1`,
density 6/π²). Total wall-time ~5 min on one core.

Findings (verbatim from results.md):

> **F1 (Lehmer-typical) REFUTED:** PRIMES log m at N=2^{18} is `−110σ`
> below BERN, `−337σ` below MATCH. Not random-equivalent.
>
> **F2 (cyclotomic / A-grade) REFUTED:** sympy factorisation of `f_N(z)`
> over `Q[z]` at N ∈ {64, 128, 256} yields `f_N = z² · g_N` with
> `g_N` **irreducible** of degree 59, 125, 249 respectively. **Zero**
> cyclotomic factors. No Φ_n(z) divides f_N. Mahler measure `m(f_N) = Θ(√N)`,
> not poly-logarithmic; no compressibility opening.
>
> **F3 (constant-deficit, B-grade) HOLDS:** `Δ(N) := log m_PRIMES(N) −
> log m_BERN(d_N)(N)` plateaus at `Δ_∞ ≈ −0.307 ± 0.001 nat` from
> `N = 2^{16}` onward. Robust at all sampled N ∈ [10³, 2.6·10⁵].
>
> **Slope identity:** `log m ≈ α log N + β` with α_PRIMES = 0.4566,
> α_BERN = 0.4577 (indistinguishable, `R² > 0.9998` for both); the
> deviation is in the intercept (`β_PRIMES − β_BERN ≈ −0.30`), not the
> slope.
>
> **Direction same as E2.17 (PH):** PRIMES log m **smaller** than density-
> matched random — primes are MORE constrained than iid in algebraic
> height, just as primes are MORE constrained in topological persistence.

**Mechanism (conjectural):** the deficit reflects that the geometric
mean of `|f_PRIMES(e^{iθ})|` over the unit circle is structurally
smaller than the iid-Gaussian / Bernoulli prediction `½ log(d N) − γ/2`.
Two candidate sources: (i) Hardy-Littlewood non-trivial pair
correlations (singular series `S(t) ≠ 1`); (ii) Vinogradov-Vaughan
major-arc structure of `S(a/q)` near rationals `a/q` with small `q`.
Quantitatively reconstructing `Δ_∞ = −0.307` from H-L singular series
or major-arc Cesàro contributions is the immediate D10 successor.

**What this rules out:** the cyclotomic / Lehmer-easy / polylog
algebraic-height regime for primes' `f_N(z)`. PRIMES `f_N` has full
generic Mahler-measure scaling `m ~ √N` and zero cyclotomic share —
no roots-of-unity-sampling polylog evaluator opening. The deficit
direction (`< 0`) makes algebraic height an *under-randomness* probe
aligned with E2.17's topological under-randomness, not an algorithmic
gateway.

**Why this is B-grade (not C):** the test produced a quantitative
**asymptotic constant deficit** `Δ_∞ ≈ −0.307` not previously measured
or reported. Combined with the cyclotomic-share zero finding, this is
two distinct structural facts about `f_PRIMES(z)`. Promoted to
**E2.20** as the **6th orthogonal pseudorandomness measure category**
(after E2.13/E2.14/E2.15/E2.16/E2.17 — additive-combinatorial / spectral
/ F_2-algebraic / point-process / metric-topological). Algebraic-
height is the new category.

**Why not A:** (i) m(f_N) still grows polynomially (Θ(√N)), no polylog
evaluator opens; (ii) the slope α is **the same** as random — only the
intercept differs; an A-grade signature would have been a different
exponent OR cyclotomic compressibility; (iii) the constant deficit
`Δ_∞` may reduce to E2.13 / E2.16 (HL singular series) on theoretical
analysis — the structural mechanism is conjectural, not yet computed.

**Successor entries (proposed in S134):**

(a) **D10.a — singular-series Cesàro fingerprint.** Compute the
Hardy-Littlewood major-arc Jensen integral
`∫₀¹ log |Σ_q μ(q)/φ(q) · {x}_q-shifts| dθ` and compare to
`Δ_∞ = −0.307`. If matches to 1%, E2.20 reduces to E2.13 / E2.16
(HL closure). If does not, E2.20 is structurally orthogonal to HL.
Cross-domain: **Vaughan's identity / Vinogradov decomposition**.

(b) **D10.b — twin-prime / Goldbach Mahler analogue.** Compute
`m(f_N^{twin}(z)) = m(Σ_{n: n,n+2 prime} z^n)` etc. Does the
deficit `Δ_∞` fingerprint different prime sub-families differently?
Outcome: **fingerprint** (each family has a distinct `Δ_∞` — A-grade
opening for new arithmetic invariant) or **collapse** (all families
same `Δ_∞` — universal under-randomness, deeper structural fact).

(c) **D10.c — Liouville Mahler measure as per-zero invariant.**
`g_N(z) = Σ_{n≤N} λ(n) z^n` has ±1 coefficients, degree N. Compute
`m(g_N)` — Liouville is "more random" than primes (E1.5: full Shannon
entropy). Predict `Δ_∞^{Liouville} ≈ 0` if Liouville is iid-typical.

**Files / directories:** `experiments/algebraic/mahler_measure_chi_p/`
— `.py`, `_results.md`, `results/` (main.json + scale.json + cyclo /
roots cross-checks). See `archive/sessions/session134_d10_mahler_chi_p.md`.

### §C.C7 — Fyodorov-Hiary-Keating extreme-value statistics of |ζ(1/2 + it)| over short windows (CLOSED 2026-04-27, S133, mode I, B-grade case 2)

**Outcome:** B-grade case 2 (informative-failure mode I) — FIRST
ζ-amplitude (vs zero-position) measurement of the project. Pre-stated
A-target (deviation > 5σ in any of three signatures with structural
arithmetic explanation) NOT met (largest 2.2σ at the widest baseline).
Pre-stated B-grade case 1 mode E (FHK match within sample noise) PARTIAL
— mean confirmed, variance/shape disagree. Pre-stated **B-grade case 2
mode I (Selberg-CLT-Gaussian fits better than FHK) MET**: KS to free
Gauss < KS to FHK Gumbel(1/2) at all 3 T (ratio 0.4–0.7); Vuong z
direction-consistent at all 3 T; variance 1.47× too large. Adds a
**genuine new structural fact** about the FHK convergence rate at
finite T, not just a duplicate closure.

Cross-domain: Gaussian multiplicative chaos (Saksman-Webb 2018 GMC
limit of ζ on mesoscopic scale; FHK 2012 PRL 108 freezing-transition
conjecture; Arguin-Belius-Bourgade 2017 RMT-side Gumbel proof). First
project use of GMC; promotes CROSS_DOMAIN_TECHNIQUES §3 row PROPOSED →
USED I with edge E7.18. Channelled mathematician: **Bourgain**
(extreme-value statistics of log-correlated random fields).

Setup: T_base ∈ {10⁴, 10⁵, 10⁶}; K = 100 unit-length windows per anchor
spaced by 10 to avoid overlap; M = 200 evenly-spaced samples per window
(spacing 0.005 ≪ inter-zero scale 0.6–0.85); ζ via mpmath dps = 15.
Total wall-time ~21 min on one core.

Findings (verbatim from results.md):

> **F1 — Mean T-INDEPENDENCE (FHK normalisation works):** M_T mean
> across T ∈ {10⁴, 10⁵, 10⁶} is {−0.699, −0.632, −0.641} ± {0.067,
> 0.083, 0.082} sem. Pairwise Z(M_T mean diff) ∈ {0.63, 0.55, −0.08}
> — all small, supporting FHK's prediction of a T-independent limit
> mean. Pooled M_∞-mean = −0.657 ± 0.045 across 300 windows.
>
> **F2 — Shape Gaussian-NOT-Gumbel (FHK shape NOT detectable at finite
> T ≤ 10⁶):** KS to free Gauss = {0.061, 0.062, 0.050}; KS to FHK
> Gumbel(1/2) = {0.088, 0.169, 0.128}. Ratio (Gauss/Gumbel) ∈
> {0.69, 0.37, 0.39} — Gauss preferred at every T-anchor by KS.
> Vuong z (Gauss vs Gumbel) = {−1.79, −1.43, −1.58} (joint Z ≈ −2.8,
> Gauss preferred at the joint level). Skewness {0.02, 0.15, 0.14}
> vs Gumbel +1.139 — distribution approximately SYMMETRIC. Excess
> kurtosis {−0.72, −0.85, −0.14} vs Gumbel +2.4 — PLATYKURTIC.
>
> **F3 — Variance 1.47× FHK prediction** (sustained across 3 T):
> empirical M_T var = {0.452, 0.692, 0.677}, FHK Gumbel(1/2) prediction
> π²/24 = 0.4112. Bootstrap 95% CIs at T ≥ 10⁵ exclude the FHK value.
>
> **F4 — Selberg-CLT secondary correction inconclusive:** observed
> Selberg-resid drops {−0.004, −0.068, −0.063} are systematically less
> negative than FHK predicts {−0.167, −0.304, −0.137}; max deviation
> Z = +2.225 at the widest baseline T = 10⁴ → 10⁶, below the 5σ
> A-grade threshold. Data weakly leans toward Selberg's `b = 1, c = 0`
> form over FHK's `b = 1, c = −0.75` but inconclusive.
>
> **F5 — Pointwise log|ζ| variance is Selberg-typical:** empirical
> {0.952, 1.131, 1.292} vs Selberg pred (1/2) log log T = {1.110,
> 1.222, 1.313}. Within 15% across all 3 T — sanity check passes.

**Mechanism for closure mode I:** Saksman-Webb 2018 proved ζ(1/2 + it)
on mesoscopic scale converges to a GMC measure; FHK Gumbel limit is a
refined consequence applying log-correlation theory to that
asymptotic GMC. The published FHK literature does NOT address the
finite-T convergence rate. Empirically: mean convergence is FAST
(T-stable at T = 10⁴ already), shape convergence is SLOW (still
~Gaussian at T = 10⁶). Pre-freezing log-correlated noise is
approximately Gaussian (CLT on the scale-summed log-correlation
kernel); the freezing transition responsible for the Gumbel heavy tail
has not yet activated at the tested T.

**What this rules out:** the FHK Gumbel-shape detection at K = 100,
T ≤ 10⁶ is closed; either larger T (Hiary's `O(t^{4/13+ε})` ζ
algorithm enabling T = 10⁹–10¹²), larger K (≥ 10⁴ windows per anchor
to discriminate sub-σ deviations), or larger window scale (mesoscopic
δ = (log T)^α per Saksman-Webb) is required.

**Why this is B-grade (not C):** the test PRODUCED a structural
positive — first quantitative finite-T FHK convergence-rate result in
project-internal or published work, plus the empirical determination
of the FHK universal intercept M_∞-mean = −0.657 ± 0.045 (giving GMC
moment-generating constant `c ≈ 0.151` under FHK Gumbel form, vs
random-matrix-side `c ≈ 0.79` Bourgade-Kuan, factor-5 finite-T gap).
Adds **EDGE E7.18**, first ζ-amplitude edge of the project,
complementary to the closed POSITION-side family (E7.1, E1.10, E3.13,
E7.15).

**Successor proposals** (per CLAUDE.md self-extension):

- **C7.a Mesoscopic-window FHK at the Saksman-Webb scale.** Repeat
  the measurement with window length δ = (log T)^{1/2} or (log T)
  rather than length 1. Saksman-Webb proved sharp GMC convergence on
  mesoscopic scales; the FHK Gumbel shape may be visible at the
  larger scale where finite-T corrections are smaller. Same
  cross-domain (GMC); single session at T = 10⁶, α = 1/2.
- **C7.b Joint argmax × prime alignment.** Argmax distribution in
  unit window is approximately uniform (KS 0.16-0.20) per S133. Test
  whether argmax POSITIONS correlate with prime-power loci
  `(p^k log p)/(2π) mod 1`. *Cross-domain*: extension of the Hardy-
  Littlewood pair-correlation density from zero-positions to
  amplitude-extremum positions. UNUSED in the project. Single session.
- **C7.c Higher-order Keating-Snaith joint moments.** μ_λ(T) :=
  E[M_T^λ] for λ ∈ {2, 3, 4} test the GMC scaling-cone exponents.
  Reuses S133 per-window data (max + second_max). 1-2 sessions; NEW
  cross-domain (Keating-Snaith joint moments conjecture
  arXiv:math/0006046). Promotes "Keating-Snaith joint moments" row to
  PROPOSED in CROSS_DOMAIN_TECHNIQUES §3.

Pointer: `experiments/analytic/zeta_structure/fhk_amplitude_max/`,
`archive/sessions/session133_c7_fhk_amplitude_max.md`.

### §D.D22 — Higher-order Hodge Laplacian L_1 spectrum of the coprimality flag complex K_N := \{σ ⊆ [2..N] : σ pairwise coprime\} (CLOSED 2026-04-27, S126, mode E, B-grade)

**Outcome:** B-grade structural negative-shape closure plus a sharp new
empirical identity (E7.17) tying L_1 spectral multiplicity to the
Bertrand-prime count `π(N) − π(N/2)`. Pre-stated A-target ("uniform
`λ_2^{(1)} > c > 0` spectral gap distinct from ER") not met. Pre-stated
B-grade case (E) ("spectra match within KS noise") strongly falsified
(KS p < 1e-300 at N=128); B-grade case (I) ("β_1 ≠ 0 with σ
significance") trivially falsified by Bertrand-prime cone collapse
forcing β_k(K_N) = 0 for k ≥ 1 deterministically. The actual outcome is
a NEW B-grade case: structural closed-form characterisation of the L_1
spectrum's top eigenvalue and its multiplicity.

Cross-domain: combinatorial Hodge Laplacian (Eckmann 1944; Friedman
1996; Horak-Jost 2013; Lim 2020 SIAM Review 62 = arXiv:1507.05379).
First quantitative project use of higher-order Hodge Laplacian; promotes
CROSS_DOMAIN_TECHNIQUES §1 row PROPOSED → USED E with edge E7.17.
Channelled mathematician: Friedman (combinatorial Laplacians).

Setup: build K_N from the coprimality graph G_N = ([2..N], gcd(i,j)=1);
construct boundary matrices B_1 (edges→vertices), B_2 (triangles→edges),
B_3 (tetrahedra→triangles); diagonalise the combinatorial Hodge
Laplacian L_1 = B_1^T B_1 + B_2 B_2^T at N ∈ {8, 12, 16, 24, 32, 48, 64,
96, 128} with 5–30 ER matched-density flag-complex controls per N.

Findings (verbatim from results.md):

> **F1 — Hodge KERNEL is deterministically trivial.** β_0(K_N) = 1,
> β_k(K_N) = 0 for k ≥ 1 ∀ N ≥ 3. Forced by Bertrand's postulate: any
> prime `p ∈ (N/2, N]` is universal in G_N (coprime to every v ≠ p in
> [2, N]) ⇒ K_N is a cone, hence contractible.
>
> **F2 — `λ_max(L_1) = |V| = N − 1` exactly** at all 9 tested N
> (numerical residual ≤ 1e-13). ER controls saturate at `p|V| +
> O(\sqrt{|V|}) < |V|` strictly.
>
> **F3 — Multiplicity identity (main empirical theorem).** For all 9
> tested N, `mult(λ_max(L_1, K_N) = |V|) = C(k+1, 2) = k(k+1)/2` where
> `k = π(N) − π(N/2)`:
> (N, k, mult): (8,2,3), (12,2,3), (16,2,3), (24,4,10), (32,5,15),
> (48,6,21), (64,7,28), (96,9,45), (128,13,91). Empirical = predicted in
> EVERY cell. Mechanism: K_N = Δ^{k-1} \ast F(H) join decomposition;
> Horak-Jost 2013 join spectrum gives multiplicity = #(0-faces ∪
> 1-faces of Δ^{k-1}) = k + C(k,2) = C(k+1, 2). Full proof = single-page
> Friedman/Horak-Jost exercise.
>
> **F4 — L_1 mean shift reduces to triple-coprime singular series.**
> Trace identity: mean(L_1) = 2 + 3|T|/|E|, so the coprime-vs-ER mean
> shift is `3(T_cp − T_ER)/|E|` with `T_cp / T_ER → ∏_p(1 − 3/p² +
> 2/p³) / (6/π²)³ ≈ 1.27628` as N → ∞. Empirical ratios across
> N=32..128: 1.241–1.273. Z[mean(L_1)] grows: 3.04 → 4.40 → 5.82 →
> 9.48 → **18.33 at N=128**. KS p-value: 1.9e-11 → 1e-300.

**Mechanism for closure mode E:** F1, F2, F3 all reduce to the
universal-vertex join decomposition (a Bertrand-postulate consequence
plus classical Horak-Jost spectral machinery); F4 reduces to the
well-known triple-pairwise-coprime singular-series identity. The Hodge
import is genuinely new at the *flag complex* level (no prior project
work on L_k for k ≥ 1 on arithmetic complexes; CLOSED 356/387/E7.12/
E7.16 all addressed L_0 only), but the *arithmetic content* collapses
to known number-theoretic identities (E2.13 family).

**Adds EDGE E7.17:** mult(λ_max(L_1, K_N) = N − 1) = (π(N) − π(N/2))
(π(N) − π(N/2) + 1) / 2. Empirical 9-cell verification + join-formula
proof sketch. First higher-order Hodge fingerprint in the project's
edge catalogue.

**Successors proposed (S126):**

#### D22.a — N^{1/2}-sieved coprimality flag complex
Restrict `K_N` to vertices `> N^{1/2}` (removes small-prime hubs).
Test whether the F3 multiplicity identity persists when small-prime
universal-degree-deficits are removed; predict `mult(λ_max) ≈ C(k+1, 2)`
unchanged since Bertrand primes survive sieving. Single-session.

#### D22.b — Non-trivial β_2 / β_3 search via universal-vertex removal
The full coprimality flag complex is contractible (F1). A *truncated*
flag complex `K^{(c)}_N := \{σ : σ pairwise coprime AND no element
exceeds c·N\}` for c ∈ {0.5, 0.4, 0.3} removes universal Bertrand
primes by hand. Test whether β_2, β_3 become nontrivial; if yes, the
non-vanishing dimensions encode arithmetic obstructions invisible to
the cone-collapsed full complex. Single-session at N ≤ 64.

#### D22.c — Hodge Laplacian of the chi_P-induced sub-flag-complex
Skip the coprimality angle; build flag complex on prime vertices only
(restrict to `[2, N] ∩ P`, with edges from coprimality = automatic for
distinct primes). The complex is the simplex `Δ^{π(N) − 1}`, so L_q
spectrum is degenerate. Modify: take vertices = primes, edges = pairs
(p, q) such that some specific arithmetic relation holds (e.g.,
twin-prime adjacency, Goldbach-pair-adjacency). Tests whether sparser
arithmetic flag complexes carry non-trivial Hodge content invisible to
the dense coprimality complex. Single-session.

Cites: E2.13 (Gowers triple-coprime), E2.14 (Anderson chi_P), E2.16
(DPP), E2.17 (PH on gaps), E2.19 (subword complexity), E7.12 (Cayley
(Z/nZ)^* spectrum), E7.16 (prime-Cayley Friedman). Distinct from
CLOSED 387 (L_0 only).

Cross-domain refs:
- Eckmann 1944 *Comment. Math. Helv.* 17, 240
- Friedman 1996 *Algorithmica* 21, 331
- Horak-Jost 2013 *Adv. Math.* 244, 303
- Lim 2020 *SIAM Review* 62(3), 685 = arXiv:1507.05379
- Goff 2009 *J. Algebraic Combin.* 30, 215 (join spectrum)
- Kahle 2009 *Discrete Math* 309, 1658 (random flag complexes)
- Kahle 2014 *Annals* 179, 1085 (sharp vanishing thresholds)

See `experiments/topological/hodge_coprimality_flag/`.

### §C.C2 — Higher-order arithmetic corrections at orders 4, 5, 6 (CLOSED 2026-04-27, S123, mode I, B-grade)

**Outcome:** No detectable arithmetic deviation from GUE at orders 4,
5, 6 of the zeta-zero point process at N=8000. First-of-kind project
measurement on this correlation-order ladder (E7.1 prior content was
GUE up to order 3 only).

Setup: 8000 Riemann zeros (γ ≤ 8148), Riemann-vM unfolded. Three
independent probes:

(A) **R_n along equally-spaced slices** `R_n(0, s, 2s, ..., (n-1)s)`
    for n ∈ {4, 5, 6} at s ∈ {0.5, 0.7, 1, 1.3, 1.6, 2, 2.5, 3, 4, 5},
    bin tol = 0.20. GUE theoretical prediction = `det[K((j-i)s)]`,
    where `K(t) = sin(πt)/(πt)`. Bootstrap noise from K=20 sub-chunks
    of 400 zeros each (matched-N to GUE batches).

(B) **k-th nearest-neighbor spacing distributions P_k(s)** for k ∈
    {0, 1, 2, 3, 4, 5}, bin width 0.10, range [0, 18]. P_k probes up
    to (k+2)-point correlation.

(C) **Higher cumulants κ_n(L)** for n ∈ {2, 3, 4, 5, 6} at L ∈ {1,
    2, 4, 8, 16, 32, 64}.

Two nulls:
- **GUE Monte Carlo pool**: K=20 batches × N=2000 Hermitian Wigner
  matrices, central 60% (1200 evs) semicircle-unfolded.
- **Gap-shuffled null**: K=20 seeds preserving the empirical gap
  distribution exactly.

**Pre-registered A-grade target**: any sustained |z| > 5σ deviation
from GUE prediction at any (n, s) ∈ {4, 5, 6} × s, OR at any (n, L)
∈ {4, 5, 6} × L, with structural explanation. **FALSIFIED.**

**Headline z-scores after matched-finite-N control:**

| Test | Probe | max \|z\| | Bin | Source of \|z\| |
|------|-------|----------|-----|-------------|
| (A) | R_4 vs theory | 2.36σ | s=2.0 | within sample noise |
| (A) | R_5 vs theory | 6.00σ raw | s=2.0 | reproduced by GUE batch (Poisson shot noise on `n_ref·tol⁴ ≈ 12`) |
| (A) | R_6 vs theory | 1.56σ | s=2.0 | within sample noise |
| (B) | P_k vs GUE pool | <1.5× σ_pool | all k | within sample noise |
| (C) | κ_n vs GUE batches, L≥8, n=3..6 | < 2.1σ | all bins | within sample noise |
| (C) | κ_2 vs GUE batches, L=32 | -9.44σ | — | known E7.1 GUE-rigidity signal |
| (C) | κ_4, κ_6 vs GUE, L=1, 2 | up to 28σ | — | Riemann-vM vs semicircle finite-N unfolding mismatch |

**Mechanism (mode I — information loss):** The Conrey-Snaith conjectural
arithmetic correction at order n scales as `1/L²` with `L = log(γ/2π)`.
At γ ≤ 8148, L ≈ 6.5, and `1/L² ≈ 0.024`. The empirical noise floor
for n-correlations is `1/√(n_tuples)` which is `≥ 0.3` at order n = 5
with tol = 0.20. The signal is 1.5+ orders of magnitude below the floor.

**Why §C.C2 is closed (mode I):** the test was the canonical extension
of S25/S45/S57 from order 3 to orders 4-6, and produces null result at
the available N. Detection of any putative arithmetic correction at
these orders requires either Odlyzko-tabulated high-height zeros (cf.
S71 closure of §C1; the `N ≥ Ω(L⁴)` scaling barrier transfers from
pair correlation to higher-order correlations) or a structurally
different (non-natural) statistic.

**Why this is B-grade (not C):** the test PRODUCED a structural
negative — first orders-4-to-6 GUE confirmation across three
independent probes, with quantitative explanation of (i) why
gap-shuffled is the wrong null for higher-order GUE-vs-arithmetic
discrimination, (ii) the Poisson-shot-noise floor at order ≥ 5,
(iii) the Riemann-vM-vs-semicircle finite-N unfolding mismatch at
small L cumulants. Refines E7.1 from "GUE up to order 3" to "GUE up
to order 6 across (P_k, R_n, κ_n)."

**What this rules out:** the natural-statistic GUE-vs-Conrey-Snaith
detection route at heights γ ≤ 10⁴ is exhausted. With S71 closing
high-height pair-correlation Odlyzko-zero re-runs as also below
detection threshold, the entire `1/L²` arithmetic-correction family
is closed at every available height.

Refines **EDGE E7.1** with an explicit orders-4-to-6 statement.
Promotes `CROSS_DOMAIN_TECHNIQUES.md` row "Sine-kernel n-correlation
determinant" from implicit USED (n=2, 3) to USED (n ≤ 6) with mode E.
CLOSED_PATHS row at S123.

**Successor proposals** (per CLAUDE.md self-extension; one with a
DIFFERENT cross-domain technique):

- (C2.a) **§C3 (bespoke statistic on zeros)** with the GUE-batch null
  built in S123. A statistic of the form `S(γ_1,...,γ_n) := Σ_p log p
  · cos(γ_i log p)` is sensitive to prime arithmetic in a way the
  pair/triple/n-correlation isn't. Same cross-domain (random matrix
  theory). Single session.
- (C2.b) **Pfaffian higher-order identities** — for a putative
  Pfaffian point process (matrix-valued kernel with anti-symmetric
  block; Vere-Jones α-determinantal generalisation), the n-point
  correlation has a Pfaffian rather than determinant form. Test
  whether empirical zeta n-correlations at n=4 satisfy Pfaffian
  identities. *Cross-domain*: matrix-valued DPP / Pfaffian point
  processes (Borodin 2009 *Encyclopaedia Math. Sci.* 152). UNUSED in
  the project (the D7 successor proposal §a flagged this for chi_P;
  here on zeros). Single session.
- (C2.c) **Joint moments of `|ζ(1/2 + iγ)|^k`** (Fyodorov-
  Hiary-Keating 2012 conjectures). FHK predicts the joint distribution
  of values along the critical line conditioned on zero locations.
  Tests whether zeros' arithmetic structure manifests through the
  *amplitude* of ζ rather than the location of zeros. *Cross-
  domain*: log-correlated random fields / Gaussian multiplicative
  chaos (Saksman-Webb 2018; FHK 2012 PRL 108). Genuinely new
  technique to the project. 2-3 sessions.

Pointer: `experiments/analytic/zeta_structure/n_correlations_4_5_6/`,
`archive/sessions/session122_c2_higher_order_zero_correlations.md`.

### §B.B2 — Automorphic L-function basis (CLOSED 2026-04-27, S118, mode E, B-grade)

**Outcome:** F_τ Hecke twisted partial sums of level-1 weight-12 cusp
form Δ are NOT a polylog basis for π(x) − Li(x); the obstruction is
quantitatively NEW.

Setup: a(n) = τ(n)/n^{11/2} computed exactly via η(q)^24 = q^{-1}Δ(q)
expansion (Euler pentagonal + repeated squaring with arbitrary-
precision int convolutions); verified against Hecke recursion at
N ∈ {5k, 10k, 20k}. M = 200 log-uniform anchors x_i in [N/50, N-1].
y(x) = π(x) - Li(x) via sympy.primepi + mpmath.li.

Three feature ensembles compared at K twist frequencies log-uniform
t_k ∈ [1, 50]:

- F_τ[i, k] = Σ_{n ≤ x_i} a(n) cos/sin(t_k log n)  (true Hecke).
- F_iid: a_iid(n) drawn i.i.d. from Sato-Tate (2/π)sin²θ density.
- F_mult: matched-multiplicative random — a(p) ~ Sato-Tate at primes,
  Hecke recursion at prime powers, multiplicativity at composites.

5-fold cross-validated OLS regression of y on each F (with intercept).
Residual rms_oos = √(mean over 5 OOS folds of squared residual).

**Pre-registered targets**: A-grade requires F_τ test rms < F_random
by ≥ 3σ; B-grade case (1) F_τ matches F_random within sample noise;
B-grade case (2) explicit numerical bound on F_τ / F_random ratio.

**Outcome: B-grade case (2) — F_τ rms is ~3× LARGER than F_random
(opposite of A-grade direction), Z = 17–58σ across 9 (N, K) cells.**

Headline at canonical N=10⁴, K_τ=8, M=200:

| basis    | rms_oos | β_oos | eff. rank |
|----------|---------|-------|-----------|
| baseline (no fit) | 12.873 | 0.322 | — |
| F_τ      |  4.239  | 0.055 | 15 / 16   |
| F_iid (10 seeds) | 1.502 ± 0.058 | 0.306 ± 0.087 | 14.9 |
| F_mult (10 seeds) | 1.566 ± 0.090 | 0.289 ± 0.090 | 15.5 |

Z(rms F_τ vs F_iid) = +47.13. Z(rms F_τ vs F_mult) = +29.75.

**K-scan at N=10⁴**: F_τ rms_oos ∈ [3.88, 4.31] across K ∈ {4, 6, 8,
10, 12}; F_random rms_oos ∈ [1.48, 1.69]. Z ≥ 17σ in safe regime
K_τ/M < 0.06.

**N-scan at K=8**: F_τ rms_oos / √N = 0.0415, 0.0424, 0.0406 at N =
5k, 10k, 20k (clean √N scaling); F_mult / √N = 0.0183, 0.0157, 0.0146.
**Ratio F_τ / F_mult = 2.85, 2.82, 2.83 — flat in N to within 1%**.

**Mechanism (E mode):** by Mellin-Perron, F_τ_k(x) is dominated by
L(s, Δ) zero contributions at heights γ_l^Δ; by the Riemann explicit
formula y oscillates at heights γ_l^ζ. Sato-Tate equidistribution and
Katz-Sarnak GUE statistics imply {γ_l^Δ} and {γ_l^ζ} are GUE-distributed
independently. F_τ is a narrow-band basis at the **wrong heights**;
F_random_mult is broadband (incoherent random superposition of all
heights). Same K, similar effective rank, but F_τ is in a "wrong
subspace" of feature-space. The ~3× obstruction ratio is consistent
with (effective spectral coverage)^{-1/2} for K=8 narrow- vs broad-band
bases.

**β_τ ≈ 0.05** (Hecke fit absorbs the √x-scaling component of y) but
the constant residual is ~3× larger than the F_random growing residual
(β_random ≈ 0.30). So Hecke captures the SHAPE (β) but not the SIZE
(rms) — structurally distinct decomposition.

**Why this closes B2:** the B2 question — "is there auxiliary g such
that Σ_{n≤x} τ(n) g(n) gives polylog approximation to π(x)?" — is
quantitatively answered no for the canonical class of g expressible
as linear combinations of {cos(t_k log n), sin(t_k log n) : k = 1..K}
at K ≤ 12. Hecke spans a structurally obstructed subspace by
Sato-Tate equidistribution; multi-form averaging cannot bypass GUE
independence (predicted by B2.b successor; would need O(infinity)
forms).

**What this rules out:** §B2 was the recommended next frontier attack
per the post-S117 critique — the only major function-field / spectral
candidate untouched. With this closure, the **automorphic-spectral
attack family** is closed alongside the four construction-side
families (AKS-modulus E7.10, Brandt MKtP E5.8, convergence-acceleration
E7.11, Maynard sieve E7.14). The L(s, Δ) zero spectrum is GUE-
independent of the ζ-zero spectrum, hence cannot supply the cancellation
needed for polylog π(x).

Adds **EDGE E7.15**: F_τ vs F_random_mult ratio = 2.83 ± 0.02 across
9 (N, K) cells, Z ≥ 17σ per cell; FIRST automorphic-spectral
measurement of the project; refines E7.1 / E1.10 / E3.13 with
quantitative L(s, Δ)-vs-ζ independence ratio. Cites E7.1, E1.10, E3.13,
E1.5, E3.1. CLOSED_PATHS row at session 118. CROSS_DOMAIN_TECHNIQUES
§1 row "Selberg trace formula" promoted UNUSED → USED (E).

**Successor proposals (per CLAUDE.md self-extension; one with a
DIFFERENT cross-domain technique):**

- (B2.a) **Higher-weight cusp forms** (k ∈ {16, 18, 20, 22, 26}):
  does ρ vary by form? Predicted: flat (all Sato-Tate). Same
  cross-domain. Single session.
- (B2.b) **Multi-form combination** (Δ ⊕ derivatives ⊕ multiple
  forms simultaneously). Tests whether multi-form averaging closes
  the spectral gap. Predicted: ρ → 1 only at K_combined ~ infinity
  (cannot bypass GUE-independence). Single session.
- (B2.c) **t-grid tuned to L(Δ) zero heights**: t_k = γ_k^Δ from
  computed L(s, Δ) zeros. Tests whether the obstruction is partly
  due to my arbitrary t-grid. Predicted: ρ unchanged (the misalignment
  is L(Δ) ↔ ζ, not t-grid ↔ L(Δ)). New cross-domain ingredient
  (L-function zero-finding for non-ζ); recommended successor.

Pointer: `experiments/algebraic/automorphic_l_function_basis/`,
`experiments/algebraic/automorphic_l_function_basis/automorphic_l_function_basis_results.md`.

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

### §D.D20 — Friedman / Ramanujan spectral gap of `Cay(Z/NZ, ±primes < N^c)` (CLOSED 2026-04-27, S125, mode E, B-grade case (i))

**Outcome:** primes-as-generators of the abelian Cayley graph
`Cay(Z/NZ, ±primes < N^c)` are **Friedman-typical within ±2σ** once
support `[2, M)` and parity (odd) are matched. The bare Ramanujan
ratio `r_N := λ_2 / (2 √(d-1))` ranges 2.05 → 11.30 (sub-Ramanujan
by orders of magnitude across `N ∈ {509, 1009, 4001, 16001, 65537}`
and `c ∈ {0.5, 2/3}`), but the deviation reduces to two trivial
finite-N effects.

Setup: FFT computation of full λ_2 spectrum of the {0,1}-indicator
of `S = ±primes < N^c` ⊂ `Z/NZ`. Four random controls × 100 seeds:
- **B1 = uniform Z/NZ** (Friedman reference).
- **B2 = support-matched [2, M)**.
- **B3 = parity-matched odd integers in [3, M)**.
- **B4 = HL W=6-coprime in [3, M)**.

Two frequency bands measured: full and minor-arc `k ∈ [N/4, 3N/4]`.

**Pre-registered F1 (A-grade)**: `r_N(prime)` deviates from
`r_N(random)` by `> 5σ` consistently with super-/sub-Ramanujan
direction and `(log N)^{-α}` shrinkage. **FALSIFIED.** Pre-registered
F2 (B-grade case (i)): `r_N(prime) ≈ r_N(random)` within ±2σ once
support and parity matched. **HOLDS.**

**Headline z-scores after support- and parity-matching:**

| Match level | Z range over 10 cells | Sign-test | Within ±2σ? |
|-------------|------------------------|-----------|-------------|
| B1 = uniform Z/NZ (full) | +4.69 .. +66.27 | 10/10 + | NO (trivial) |
| B2 = support [2,M) (full) | +0.68 .. +1.87 | 10/10 + | YES |
| B3 = odd [3,M) (minor) | -31 .. -15622 | 10/10 - | NO (parity) |
| B3 with primes - {2} (minor) | +0.51 .. +2.07 | 10/10 + | YES |

**Mechanism (E mode):** the two non-trivial spectral spikes in the
prime-Cayley FFT are at (i) `k = 1` (low-frequency major-arc; all
generators p < M ≪ N have `cos(2π p/N) ≈ 1`, giving `λ_1 ≈ d`;
Vinogradov's prime-exp-sum bound does not apply since `q = N` is
not bounded by a fractional power of M), and (ii) `k ≈ (N-1)/2`
(parity frequency; all primes > 2 are odd, giving full alignment
modulated by the single even prime p=2 which contributes `(-1)^2 =
+1` instead of `-1`, reducing the peak by exactly 4 units). After
controlling for both, no Hardy-Littlewood mod-q residual is detected
at the scales tested — the dominant `1/N²` finite-N spikes are an
order of magnitude larger than any HL singular-series correction.

**Why §D20 is closed (mode E):** the cross-domain Friedman 2008
random-regular-graph spectral gap reference performed real work —
defining the Ramanujan-typicality null hypothesis, the matched
control distribution, the Bonferroni multiple-comparison correction.
Both deviations reduce to closed-form finite-N effects matching
empirics quantitatively. No algorithmic opening for π(x).

**What this rules out:** §D20 was the recommended next frontier in
the abelian-Cayley-spectral category. With this closure, the
**prime-as-generator Cayley spectral gap family** is closed alongside
E7.12 (S79 fixed-generator (Z/nZ)\* spectrum probes ω(n)) and
E7.13 (S80 Szegedy walks on arithmetic graphs do not yield polylog
mixing). The abelian-Cayley side of the algebraic-Cayley frontier
is now exhausted.

Adds **EDGE E7.16**: r_N(prime) of `Cay(Z/NZ, ±primes < N^c)`
indistinguishable from parity-and-support-matched random subsets
within ±2σ across 10 (N, c) cells; bare deviation reduces to
bounded-support k=1 spike + parity k≈N/2 spike modulated by the
single even prime p=2. Cites E2.13, E2.14, E2.18, E7.12, E7.15.
CROSS_DOMAIN_TECHNIQUES.md §1 row "Random regular graph spectral
gap (Friedman)" promoted PROPOSED → USED-E with edge E7.16.
CLOSED_PATHS row at session 125.

**Channelled mathematician**: Bourgain. The Bourgain framing
identified the right diagnostic question (separate major-/minor-arc
behaviour, k=1 vs k≈N/2) and the right mechanism (low-frequency
major-arc dominates trivially, parity at next order, no genuine
arithmetic content).

**Successor proposals** (per CLAUDE.md self-extension; one with a
DIFFERENT cross-domain technique):

- (D20.a) **Cheeger constant / edge expansion** of the prime-Cayley
  graph instead of `λ_2` saturation. Cheeger inequality `λ_2 ≤ 2 h_G
  ≤ √(2 d λ_2)` (Hoory-Linial-Wigderson §2): measuring `h_G`
  independently from `λ_2` is a structurally orthogonal expansion
  check. Cross-domain: combinatorial expansion / Cheeger inequality.
  Single session.
- (D20.b) **Extend to c = 1 + primorial N** (so `Z/NZ` is maximally-
  multi-component) — does the cyclic vs primorial structure of N
  introduce a deviation invisible at prime N? Same cross-domain
  technique. Single session.
- (D20.c) **Non-abelian Cayley graph `Cay(SL_2(F_p), prime
  generators)`** (Lubotzky-Phillips-Sarnak Ramanujan graphs).
  Tests whether non-commutative spectral gap is more discriminating.
  *Cross-domain: arithmetic Ramanujan graphs* (LPS 1988
  *Combinatorica* 8 — UNUSED in §1 of CROSS_DOMAIN_TECHNIQUES.md).
  Recommended successor. 2 sessions.

Pointer: `experiments/algebraic/friedman_ramanujan_prime_cayley/`,
`experiments/algebraic/friedman_ramanujan_prime_cayley/friedman_ramanujan_prime_cayley_results.md`.

---

*A frontier attack that fails honestly is more valuable than ten
refinements that succeed. Be ambitious. Then be honest.*

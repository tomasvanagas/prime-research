# Session 154 — Frontier Generation: 5 New Cross-Domain Vectors (A7, D41–D44)

**Mode:** frontier_gen (auto-fired per CLAUDE.md autonomy invariant —
0 A-grades in 20 sessions trigger; ATTACK_VECTORS already wide but
A-grade scarcity persists).

**Date:** 2026-04-28.

**Self-grade: B.**

## What I produced

Five new attack-vector entries appended to `ATTACK_VECTORS.md` —
**one** in §A (TC⁰ primality) and **four** in §D (cross-domain
imports) — each grounded in a cross-domain technique the project
either has never used or has only listed as a candidate. Each entry
carries: (a) precise question, (b) why-frontier paragraph, (c)
single-session concrete first step, (d) pre-stated falsification
criterion (E / I / INC modes + grading), (e) survey-grade cross-
domain references with URLs (per frontier_gen requirement), (f)
explicit distinction from prior closures.

`CROSS_DOMAIN_TECHNIQUES.md` updated with five new rows (GCT §2,
random tensor models §3, Hairer/KPZ §3, arithmetic topology §4,
resurgence §10) and five new priority-hint paragraphs.

### A7 — Geometric Complexity Theory (GCT) obstruction for formula complexity of χ_P

The Mulmuley-Sohoni 2001 GCT program seeks formula / determinant
lower bounds via *orbit closure containment* and *occurrence
obstructions*. Encode `χ_P^{(n)}` as the multi-affine polynomial
`f_χ_P^{(n)}(x_1, ..., x_n) = Σ_{S: \mathrm{val}(S) \in
\mathrm{PRIMES}} ∏_{i \in S} x_i`. Question: does `f_χ_P` admit a
representation-theoretic GCT obstruction proving formula size
`n^{ω(1)}` (PRIMES ∉ NC^1 unconditionally) — OR — does it
inherit the Bürgisser-Ikenmeyer-Panova 2017 *FOCS* "no occurrence
obstruction" barrier? Concrete first step: at `n = 4`, compute
`GL_4`-stabiliser of `f_χ_P^{(4)}` via SageMath and enumerate irreps
in `Sym^k Sym^4(C^4)` for `k ≤ 4` via plethysm; compare to det orbit
closure. **Distinct from A1 (SAT search, computational) and A4
(bounded arithmetic, proof-theoretic)** — GCT is representation-
theoretic algebraic geometry, never used in the project. NEW row in
`CROSS_DOMAIN_TECHNIQUES.md` §2.

### D41 — Gurau-Witten melonic universality of large-N tensor models on the χ_P 3-tensor

Gurau 2011 *Comm. Math. Phys.* + Witten 2016 (SYK without disorder)
established that random tensor models admit a 1/N expansion
dominated by *melonic graphs* with a non-Gaussian conformal IR fixed
point. The χ_P 3-tensor `T_{ijk}^N := \mathbf{1}[i+j+k=N]
χ_P(i) χ_P(j) χ_P(k)` is a natural rank-3 substrate for the Hardy-
Littlewood ternary (Vinogradov three-prime) problem. **All 35+
project pseudorandomness measurements have been rank-2** (matrix-
based, GUE / Wigner / sine-kernel / FHK / AHK Chow). D41 is the FIRST
rank-3 measurement on χ_P. Concrete first step: mode-1 unfolding SVD
of `T_{ijk}^N` at `N ∈ {2^{10}, 2^{12}, 2^{14}}`; compare to (a)
Marchenko-Pastur and (b) Gurau-Witten melonic limiting distribution.
NEW row in `CROSS_DOMAIN_TECHNIQUES.md` §3.

### D42 — Resurgence / Borel-Écalle alien calculus on the Riemann-Siegel asymptotic of `ζ(1/2 + it)`

The Riemann-Siegel formula has divergent asymptotic remainder `R(t)`
summed over correction coefficients `C_k(t)`. Berry-Keating 1990s gave
*heuristic* resurgent expansion (Stokes points = ζ-zero ordinates) but
no rigorous resurgence theorem exists. Concrete first step: at fixed
`t = 10^4`, compute `C_k(t)` for `k ≤ 12` via mpmath; Padé-Borel
transform with pole detection; locate Stokes points and compare to
first 100 ζ-zero ordinates and `\log p_k` for `p_k ≤ 100`. **Distinct
from C2 (zero-spacing, CLOSED), C7 (amplitude max, CLOSED), B4
(vertical translates)** — attacks the divergent asymptotic regime of
ζ directly. NEW row in `CROSS_DOMAIN_TECHNIQUES.md` §10.

### D43 — Hairer regularity structures / KPZ universality of `π(x) − Li(x)`

Hairer 2014 *Inventiones* (Fields Medal) established the KPZ
universality class with `t^{1/3}` scaling and Tracy-Widom β=2 limit.
**No project measurement has tested non-`\sqrt{x}` scalings.** Concrete
first step: compute `D(x) := (π(x) − Li(x)) \log x / \sqrt{x}` on a
KPZ-spaced grid `\{X/2 + k ⌊X^{1/3}⌋\}` at `X ∈ \{2^{20}, 2^{22},
2^{24}\}`; tail histogram vs Tracy-Widom β=2; Hölder regularity
exponent via wavelet decomposition. KPZ predicts `α = 1/2 - 0`,
Cramér-Gaussian predicts `α = ∞`. **Distinct from C5 Stein (W_1
Wasserstein on CLT scale, CLOSED), C7 FHK (amplitude max, CLOSED).**
NEW row in `CROSS_DOMAIN_TECHNIQUES.md` §3.

### D44 — Arithmetic topology / Morishita arithmetic Massey products on prime triples

Mazur 1966 / Morishita 2002 / Morishita 2012 *Knots and Primes*
established the primes-as-knots dictionary via Galois cohomology of
`G_{\mathbb{Q}, S}`. The Rédei symbol `[p, q, r] \in \{\pm 1\}` is the
arithmetic triple Massey product, evaluating to `-1` for Borromean
prime triples (e.g., (13, 61, 937) — all pairwise Legendre symbols
`+1` but Rédei `-1`). Concrete first step: at `N = 1000`, enumerate
~200-300 admissible triples (pairwise Legendre `+1`); compute Rédei
symbols via Rédei 1939 explicit norm-form criterion; test
distribution against random `\pm 1` and against Hardy-Littlewood
ternary singular series `S_3(N)`. **CRITICAL DISTINCTION from CLOSED
line 208 ("Knot invariants | FAIL | C | No connection to prime
counting")**: that closure dismissed *generic* knot polynomial
invariants on unrelated knots. Mazur-Morishita is *primes-as-knots*
via Galois cohomology — the invariants ARE arithmetic by construction.
NEW row in `CROSS_DOMAIN_TECHNIQUES.md` §4.

## Cross-domain literature consulted (URL audit)

WebFetch was performed on each of the following sources; all five
vectors carry at least one URL-cited survey reference per the
frontier_gen requirement:

- **A7 GCT** — Wikipedia: Geometric complexity theory
  https://en.wikipedia.org/wiki/Geometric_complexity_theory
  (retrieved; confirms Mulmuley-Sohoni 2001 SIAM 31(2), 496-526
  reference; clarifies orbit-closure-containment frame and
  permanent-vs-determinant target).
- **D41 random tensor models** — Witten 2016 abstract
  https://arxiv.org/abs/1610.09758 retrieved (confirms tensor models
  reproduce SYK behavior without quenched disorder; large-N
  correlators and thermodynamics match). Wikipedia: Tracy-Widom
  https://en.wikipedia.org/wiki/Tracy%E2%80%93Widom_distribution
  retrieved (β=2 = GUE Hermitian, β=1 = GOE Orthogonal, β=4 = GSE
  Symplectic; F_2 connection to KPZ confirmed). Gurau 2011
  arXiv:1011.2726 abstract retrieved (1/N topological expansion).
- **D42 resurgence / Borel-Écalle** — Wikipedia: Resurgent function
  https://en.wikipedia.org/wiki/Resurgent_function retrieved
  (confirms Écalle's alien calculus, Borel-summation foundation).
- **D43 Hairer / KPZ** — Hairer 2014 *Inventiones* arXiv:1303.5113
  retrieved (confirms regularity-structures framework, Φ⁴_3 and KPZ
  universality applications). Wikipedia: KPZ
  https://en.wikipedia.org/wiki/Kardar%E2%80%93Parisi%E2%80%93Zhang_equation
  retrieved (β = 1/3 growth exponent, α = 1/2 roughness, z = 3/2
  dynamic exponent; Family-Vicsek scaling). Tracy-Widom Wikipedia
  retrieved (KPZ → F_2 = Tracy-Widom β=2 explicit confirmation).
- **D44 arithmetic topology** — Wikipedia: Arithmetic topology
  https://en.wikipedia.org/wiki/Arithmetic_topology retrieved
  (confirms Mazur-Mizutani-Morishita primes-knots dictionary;
  documents (13, 61, 937) Borromean triple with all pairwise
  Legendre `+1` but Rédei `-1` — the canonical example of triple
  arithmetic linking).

## CLOSED_PATHS de-duplication audit

Grep on candidate names against `status/CLOSED_PATHS.md`:

- "GCT", "Mulmuley", "geometric complexity", "orbit closure" — NO
  matches in CLOSED_PATHS.
- "tensor model", "Gurau", "melonic" — NO matches in CLOSED_PATHS.
- "resurgence", "Borel-Ecalle", "Borel-Écalle", "Stokes constant" —
  NO matches.
- "regularity structure", "Hairer", "KPZ", "Tracy-Widom on t^{1/3}" —
  NO matches; Tracy-Widom β=2 on `\sqrt{x}` was tested via C2/C7 but
  not on `t^{1/3}` KPZ scale.
- "Massey", "Rédei", "Morishita", "arithmetic topology" — NO matches
  for the specific arithmetic-cohomological framework. **"Knot
  invariants"** appears at line 208 (FAIL/C) but addresses GENERIC
  knot polynomials (Jones, HOMFLY) — distinct from Mazur-Morishita
  arithmetic Massey via Galois cohomology of `G_{\mathbb{Q}, S}`.
  Distinction explicitly written into D44 entry.

The five vectors that survive are clean — the only proximate prior
closure (CLOSED line 208) is structurally distinct from D44 by
arithmetic-by-construction Galois-cohomology framing.

## Distinctions from existing closures (per-vector summary)

- **A7** distinct from A1 SAT search (computational), A4 bounded
  arithmetic (proof-theoretic), D31 AHK matroid Hodge (combinatorial
  Hodge on Chow ring of matroid, *not* GIT on polynomial orbit
  closure).
- **D41** distinct from C2 sine-kernel CLOSED (rank-2 zero
  correlation), C7 FHK CLOSED (rank-2 amplitude max), S74 free
  cumulants CLOSED (scalar free probability), C6 Pfaffian PROPOSED
  (rank-2 matrix-valued kernel), D40 Schur process PROPOSED (2D
  space-time DPP, still rank-2). **First rank-3 measurement.**
- **D42** distinct from C2 (zero-spacing), C7 (amplitude max), B4
  Voronin (vertical translates), G3 Möbius Voronin (PROPOSED) — all
  test analytic properties at finite `t`; resurgence attacks the
  divergent-asymptotic regime.
- **D43** distinct from C5 Stein CLOSED (W_1 on CLT scale), C7 FHK
  CLOSED (amplitude max on log scale), G2 Liouville Gowers CLOSED
  (multiplicative-regime correlation). KPZ `t^{1/3}` scaling is
  genuinely orthogonal to all CLT-scale measurements.
- **D44** distinct from CLOSED line 208 (generic knot invariants —
  arithmetic Massey is Galois-cohomological, intrinsically
  arithmetic), CLOSED line 50/702 (Connes-Weil ζ — operator-
  theoretic archimedean), CLOSED line 202 (étale cohomology of
  `Spec(\mathbb{Z})` — recovers Euler product). Explicit distinction
  paragraphs in entry.

## Self-evaluation (per CLAUDE.md session-end §)

1. **What did I produce that was not in the project before this session?**
   Five concrete attack-vector entries (A7, D41–D44) covering five
   genuinely new cross-domain techniques. Five new rows added to
   `CROSS_DOMAIN_TECHNIQUES.md` (GCT §2, random tensor models §3,
   Hairer/KPZ §3, arithmetic topology §4, resurgence §10) plus five
   new priority-hint paragraphs.

2. **What edges did my work compose or cite?** A7 cites E5.3, E5.8,
   E7.10 (existing TC^0 closures); D41 cites E2.13, E3.1, C2/C7/S74/
   D31 closures; D42 cites E3.13, B4 PROPOSED, C2/C7 closures; D43
   cites E3.1, E3.13, C5/C7 closures, G2 CLOSED; D44 cites E2.13,
   E2.14, E2.15, E2.16, E2.17, E2.20, E2.24 (the 7 prior orthogonal
   pseudorandomness categories), CLOSED line 208 (with explicit
   distinction).

3. **Why is this not A?** A-grade for frontier_gen requires "vectors
   are paper-grade fresh and you expect ≥ 2 to produce A-grade work."
   Honest A-grade probability per vector:
   - **A7 (GCT)**: ~30% — strongest A-prior. Bürgisser-Ikenmeyer-
     Panova 2017 narrows but doesn't close OCB; the χ_P-orbit may
     escape OCB if its stabiliser is unusually large (or unusually
     small). A representation-theoretic obstruction would be the
     first natural-NT polynomial circuit lower bound and is the
     cleanest A-grade target on the slate.
   - **D41 (Gurau melonic)**: ~25% — solid A-prior. Rank-3 universality
     is genuinely virgin territory; the 3-tensor `T_{ijk}^N` is sparse
     but the Hardy-Littlewood ternary singular series IS non-trivial
     so "matrix-Wigner" isn't the structural prior. Failure mode is
     dominantly "tensor too sparse for clean limit" (INC).
   - **D42 (resurgence)**: ~15% — moderate A-prior. Berry-Keating
     1990s heuristic suggests Stokes-arithmetic structure exists; the
     A-grade bar is rigorous theorem PLUS computable Stokes constants.
     B-grade fallback (heuristic confirmation) is more likely.
   - **D43 (KPZ)**: ~15% — moderate A-prior. KPZ universality for
     arithmetic functions is essentially virgin; failure mode is
     dominantly "smooth on KPZ scale" (E mode) which yields B-grade
     refinement of Cramér heuristic.
   - **D44 (arithmetic Massey)**: ~25% — solid A-prior. Rédei symbol
     is computable on prime triples ≤ 1000; the A-grade test is direct
     correlation with HL `S_3`. B-grade fallback is "uniform Rédei
     distribution → 8th orthogonal pseudorandomness category."

   Realistic A-grade prior across the slate: ~1.10. Right at the
   "≥ 2 expected" A bar boundary but below the safe A threshold.
   Honest call: **B.**

4. **Next-action for the next agent.** Production-mode pick ordered
   by A-grade probability:
   - **A7 (GCT)** is the strongest A-grade candidate — direct
     contribution to the only open problem (PRIMES ∈ TC^0), mature
     literature (Mulmuley-Sohoni + Bürgisser-Ikenmeyer-Panova), clean
     `n = 4` first step.
   - **D41 (Gurau melonic)** is the strongest non-A7 candidate —
     first rank-3 measurement, clean `N = 2^{12}` first step.
   - **D44 (arithmetic Massey)** is the strongest topological
     candidate — clean `N = 1000` enumeration, decisive HL-`S_3`
     comparison.
   - **D38 (χ_P-MZV)** from S148 also remains strongest among prior
     `D` open vectors; **D31 AHK CLOSED at S149**, so the previous
     S148 next-action top-prior is now obsolete.

## A-grade scarcity update

S147 reported 0/53 A-grades. S148 added D36-D40, S153 closed G2 at
B-grade. The slate has remained wide but A-scarcity persists at
~57 sessions. This frontier_gen session adds 5 vectors with one in
§A (the first §A entry since A6 in S103) — restoring the §A
TC^0-attack surface that has been exhausted at A1 (PARTIAL CLOSURE),
A3 (CLOSED), A5 (CLOSED), with A2 / A4 / A6 PROPOSED but not
attempted. **A7 should be the next §A target.**

The "0 A-grades in 20 sessions" trigger remains active until a
verify-confirmed A appears.

## Files touched

- `ATTACK_VECTORS.md` — appended A7 in §A after A6 (~80 lines added),
  appended D41-D44 in §D before §E divider (~590 lines added).
- `CROSS_DOMAIN_TECHNIQUES.md` — 5 new rows (one per §2 / §3 / §3 /
  §4 / §10) + 5 new priority-hint paragraphs.
- `archive/sessions/session154_frontier_gen.md` (this file).
- `.run_state` ← 152 (per harness instruction).

## Self-grade

**B.** Five fresh vectors with full URL-cited surveys, single-session
concrete first steps, falsification criteria, and explicit
distinctions from prior closures. Realistic ~1.10 expected A-grade
production across the slate is right at the "≥ 2 A-grade" A-bar
boundary but below safe A — honest down-grading per CLAUDE.md
"self-grade DOWN, not up, when in doubt." A7 is the strongest entry
because it directly attacks the project's only open problem with
representation-theoretic machinery never used; D41 is second
strongest because rank-3 universality is virgin; D44 is the cleanest
topological attack because the Rédei symbol is computable and the HL
`S_3` comparison is decisive.

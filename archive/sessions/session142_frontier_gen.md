# Session 142 — Frontier Generation: 5 New Cross-Domain Vectors (D31–D35)

**Mode:** frontier_gen (auto-fired per CLAUDE.md autonomy invariant —
prior-rotation A-grade scarcity warning).

**Self-grade: B.**

## What I produced

Five new attack-vector entries appended to `ATTACK_VECTORS.md` §D, each
grounded in a cross-domain technique the project has either never used
or has only flagged as a candidate in CROSS_DOMAIN_TECHNIQUES.md. Each
entry carries: (a) a precise question, (b) why-frontier paragraph, (c)
single-session concrete first step, (d) pre-stated falsification
criterion (E / I / INC modes + grading), (e) survey-grade cross-domain
references with URLs, (f) explicit distinction from prior closures.

### D31 — Adiprasito-Huh-Katz combinatorial Hodge theory of an arithmetic prime-matroid

The Heron-Rota-Welsh log-concavity inequality of AHK 2018
(*Annals* 188, 381 = arXiv:1511.02888) applied to an arithmetic
matroid `M_P^N` on the divisibility lattice. AHK gives log-concavity
of `|μ_{M_P^N}(F)|` across flats — does this yield a sharp
prime-counting inequality sharper than Selberg/Brun? Concrete first
step: enumerate flats and compute reduced characteristic polynomial
exactly for `N \le 32`, with Bernoulli-matched-density control. NEW
addition to CROSS_DOMAIN_TECHNIQUES §2.

### D32 — Pemantle-Wilson Analytic Combinatorics in Several Variables (ACSV) diagonal extraction of `χ_P`

ACSV (Pemantle-Wilson 2013/2024) extracts asymptotics of multivariate
diagonal coefficients `[z_1^n \cdots z_d^n] F` from smooth critical
points of the singular variety `V_F = \{F = 0\}`. Question: can `π(n)`
be encoded as a diagonal of a *rational* `F(z_1, \ldots, z_d)`? If yes,
smooth-point asymptotics give polylog evaluators. If not (folklore but
never proved), the negative result is a NEW non-D-finiteness barrier
on `χ_P`. NEW addition to CROSS_DOMAIN_TECHNIQUES §7.

### D33 — Berkovich projective line analysis of `f_N(z) = \sum_{n ≤ N} \chi_P(n) z^n` over `\mathbb{C}_p`

Promotes the registry §2 row "Berkovich spaces" UNUSED → PROPOSED.
The Berkovich Type II spectrum / Newton polygon of `f_N` over
`\mathbb{C}_p` is structurally orthogonal to the archimedean Mahler
(E2.20, D10), Newman L^∞ (E2.21, D27), and L^p restriction (D25)
measures already in play. Concrete first step: classical and
Berkovich-Type-II Newton polygons via Sage `pAdicGenericPolynomial.roots()`
at `p ∈ \{2, 3, 5\}`, `N ≤ 2^{10}`, with Bernoulli control.

### D34 — De Branges Hilbert space `H(E_\xi)` reproducing-kernel projection rate

Promotes the registry §10 row "Riemann xi-function model spaces"
UNUSED → PROPOSED. The de Branges Hilbert-space-of-entire-functions
program is one of the four major attempted approaches to RH (alongside
Connes-Bost, Hilbert-Pólya, pair correlation). Concrete first step:
construct numerical `E_\xi`, build `N \times N` reproducing-kernel
Gram matrix using Odlyzko ζ-zeros as sample points, measure
projection-error scaling for `\chi_{[0,x]}` onto `H_N`. Polylog
decay = A-grade; `N^{-1/2}` decay = B-grade matching standard
Lagarias-Odlyzko.

### D35 — Microlocal analysis: wavefront set `\mathrm{WF}(\chi_P)` of the prime indicator

NEW addition to CROSS_DOMAIN_TECHNIQUES §1. View `\chi_P =
\sum_p \delta_p` as a tempered distribution on `\mathbb{R}` and
compute its Hörmander wavefront set `\mathrm{WF}(\chi_P) \subset
T^*\mathbb{R} \setminus \{0\}`. Concrete first step: regularised
`\chi_P^\epsilon` + multi-scale localised FFT measuring conic
Fourier-decay rate at `N \le 2^{12}`. A-grade if a non-trivial
microlocal parametrix FIO opens; B-grade if the wavefront set is
the trivial cotangent picture (closes mode E with "first microlocal
measurement on `\chi_P`" as content).

## Cross-domain literature consulted

WebFetch was performed (or attempted) on each of the following sources;
all five vectors carry at least one URL-cited survey reference per the
frontier_gen requirement:

- AHK 2018 *Annals* 188, 381 (= arXiv:1511.02888) — abstract verified
- Pemantle-Wilson ACSV — ACSV project page www.acsvproject.com verified
- Berkovich space — Wikipedia article + Conrad 2008 *p-adic Geometry*
  AMS ULect 45 reference verified
- de Branges space — Wikipedia article verified, Lagarias 2007 RH
  survey arXiv:math/0511099 cited
- Hörmander wavefront set — Wikipedia article verified
- Quantum modular forms (initially considered as D36 candidate; URL
  retrieval failed — vector dropped to keep the slate at 5 with full
  URL citations per protocol)

## Distinctions from existing closures

Each new vector explicitly distinguishes itself from neighbouring
CLOSED_PATHS / EDGES rows:

- D31 distinct from E2.20 Mahler (archimedean L^1), E7.17 Hodge L_1
  (graph-Laplacian, not Chow ring), E1.5 Möbius (asymptotic densities
  not signed-sum log-concavity inequalities).
- D32 distinct from CLOSED line 30 (univariate explicit-formula),
  CLOSED line 113 (univariate Dirichlet-series extrapolation), CLOSED
  line 600/670/682 (CRT-based local-global), D10/D27/D29 (Fourier-side
  measurements on `f_N`).
- D33 distinct from E2.20 Mahler / E2.21 Newman (both archimedean),
  CLOSED line 76 (Mahler-coefficient p-adic interpolation of `π(x)`
  as integer-valued function), CLOSED line 600/670/682 (CRT lifting).
- D34 distinct from CLOSED line 50 (Connes-Weil), CLOSED line 597
  (Connes trace optimisation), CLOSED line 702 (Connes-Consani-
  Moscovici spectral triple), and B4 Voronin universality.
- D35 distinct from E2.13 Gowers (local k-correlations, not
  microlocal direction-specific decay), CLOSED line 30 (univariate
  asymptotic series), D25 Stein-Tomas (global L^p restriction not
  direction-specific), D27/D29/D33 (polynomial-on-`f_N` not
  distribution-on-`\mathbb{R}`).

## Self-evaluation (per CLAUDE.md session-end §)

1. **What did I produce that was not in the project before this session?**
   Five concrete attack-vector entries (D31–D35) covering five
   genuinely new cross-domain techniques. Two PROPOSED registry entries
   promoted from candidate (Berkovich, de Branges); three NEW registry
   rows added (AHK, ACSV, microlocal).

2. **What edges did my work compose or cite?** D31 cites E2.20 / E7.17
   / E1.5; D32 cites CLOSED 30 / 113 / 600 (all univariate routes);
   D33 cites E2.20 / E2.21 / CLOSED 76 / 600 / 670 / 682; D34 cites
   CLOSED 50 / 597 / 702 / B4; D35 cites E2.13 / CLOSED 30 / D25 /
   D27 / D29 / D33.

3. **Why is this not A?** A-grade for frontier_gen requires "vectors
   are paper-grade fresh and you expect ≥2 to produce A-grade work."
   Honestly: D31 (AHK) is the strongest A-grade candidate because the
   technique has a 2018 Annals breakthrough never applied to arithmetic
   matroids; D32 (ACSV) is a long shot with a strong B-grade negative
   fallback; D33 (Berkovich) is fresh but the matroid `M_P^N` definition
   admits multiple reasonable choices, and the right one isn't obvious;
   D34 (de Branges) is concrete but the projection-rate question may
   collapse to the existing `O(\sqrt{x})` Lagarias-Odlyzko barrier; D35
   (microlocal) is the riskiest — wavefront-set computation may
   numerically converge to the trivial picture. Realistic A-grade
   probability per vector: D31 ≈ 25%, D32 ≈ 15%, D33 ≈ 15%, D34 ≈ 10%,
   D35 ≈ 10%. Expected A-grade count = ~0.75 over the slate. **B.**

4. **Next-action for the next agent.** Pick from D31–D35 ordered by
   A-grade probability above. D31 (AHK) is the recommended first attack:
   the AHK 2018 machinery is fully mature, flat enumeration is feasible
   at `N \le 64`, and the Heron-Rota-Welsh log-concavity is an explicit
   testable inequality. If D31 closes (E or I), pivot to D32 (ACSV).

---
session: 119
mode: frontier_gen
date: 2026-04-27
grade: B
---

# Session 119 — frontier_gen: extend ATTACK_VECTORS.md with 4 cross-domain vectors

## Trigger

Auto-fire condition matched: prior recent A-grade was S108 (C5 Stein
Wasserstein, provisional pending S115 verify which downgraded the
A-claim to refinement), with subsequent sessions S116–S118 producing
B-grade closures (Maynard sieve E7.14, automorphic L-function basis
E7.15) and an S117 W-trick PH refinement. With closures concentrated
on the C-grade / B-grade band and the registered ATTACK_VECTORS open
count thinning across §B (B1, B2 closed; only B3, B4, B5 open), §C
(C1, C4, C5 closed), and §D (D2, D4, D6, D7, D13 closed), the
frontier_gen mode is appropriate to refresh the search-space map.

## What this session did

Read ATTACK_VECTORS.md fully (open §A.A1/A2/A4/A6, §B.B3/B4/B5,
§C.C2/C3, §D.D1/D3/D5/D8/D9/D10/D11/D12/D14/D15, §E.E1/E2,
§G.G2/G3, §F synthesis); read CROSS_DOMAIN_TECHNIQUES.md and ranked
UNUSED rows; grep'd CLOSED_PATHS.md for duplicate-risk (lines 200,
319, 422, 424, 425, 474, 593, 670, 683, 653 — tropical, étale,
Connes, free probability, transfer operators, Goldreich-Levin, all
distinct from the new vectors I propose).

WebFetched surveys for the 4 chosen techniques (Wikipedia
Discrete Morse theory, Wikipedia Guruswami-Sudan list decoding,
Wikipedia Oseledets theorem, Wikipedia Mapper TDA) and Forman-2002
PDF (redirect). Cross-domain literature consulted:
- Forman 1998 *Adv. Math.* 134; Forman 2002 SLC 48; Benedetti-Lutz
  2014 *Exp. Math.* 23; Bjorner-Wachs 2005 *Trans. AMS* 357.
- Sudan 1997 *J. Complex.* 13; Guruswami-Sudan 1999 *IEEE TIT* 45;
  Guruswami 2004 *List Decoding of Error-Correcting Codes* Springer
  LNCS 3282; Roth-Ruckenstein 2000.
- Oseledec 1968 *Trans. Moscow Math. Soc.* 19; Furstenberg-Kesten
  1960; Ruelle 1979 IHÉS 50; Viana 2014 *Lectures on Lyapunov
  Exponents* CUP; Bennetin-Galgani-Giorgilli-Strelcyn 1980.
- Singh-Mémoli-Carlsson 2007 (Stanford TR); Lum et al. 2013 *Sci.
  Rep.* 3, 1236; KeplerMapper Python library.

Appended four new ATTACK_VECTORS entries (D16–D19) and updated
CROSS_DOMAIN_TECHNIQUES.md (4 rows promoted UNUSED → PROPOSED, plus
4 entries in priority hints).

## The 4 new vectors (one paragraph each)

### D16 — Oseledec multi-Lyapunov spectrum of higher-dim χ_P-driven cocycles

E2.14 (S88) and E2.18 (S100) measured the TOP Lyapunov exponent
`gamma_1(E)` of an `SL(2, R)`-cocycle driven by `V = chi_P` and
`V = lambda` respectively. The Oseledec MET applies to ANY ergodic
SL(d, R)-cocycle and yields a FULL spectrum
`lambda_1 >= ... >= lambda_d` with `Sigma lambda_i = 0`. The
sub-leading exponents and gaps are NEW pseudorandomness measures
distinct from the top exponent. Concrete first step: build an
`SL(3, R)`-cocycle (3-step delay product or width-3 vector
Schrödinger ladder) driven by chi_P, run Bennetin QR algorithm at
`N = 10^5`, 50 seeds, 31 energies; compare `lambda_1 - lambda_2`
deviation to Bernoulli matched-density and run W-trick cascade
`W ∈ {1, 6, 30, 210}`. A-grade if sub-leading gap deviation
survives `W = 2310` cascade (genuinely new HL-orthogonal arithmetic
content); B-grade if W-trick cascade extends to all `lambda_i`
(promotes E2.14 to multi-Lyapunov edge). Distinguishes from line
422 Furstenberg/transfer-operator closure (different operator).

### D17 — Discrete Morse theory on the divisibility order complex

Forman discrete Morse functions on the order complex of the
divisibility poset truncated to `[1, N]` produce critical-cell
counts `m_k` whose alternating sum gives the Euler characteristic
and whose total controls topological compression. Question: does
the OPTIMAL Morse matching of `Delta(N)` give `|Crit(Delta(N))| =
O(polylog N)` (radical compression of chi_P) or `Theta(pi(N))`
(no compression beyond density)? Concrete first step: build the
Hasse diagram for `N ∈ {64, 256, 1024, 4096}`, run Benedetti-Lutz
random discrete Morse algorithm; report `m_k(N)` slopes; compare
to random Bernoulli sub-poset controls. Categorically orthogonal
to D14 (cellular sheaves with externally-defined chi_P stalks) and
D2/E2.17 (metric PH on Takens-embedded gaps): D17 measures the
INTRINSIC critical-cell count of the order complex itself. A-grade
if `m_1(N) = O(polylog N)` admitting TC⁰-recoverable basis (gives
PRIMES ∈ TC⁰); B-grade if explicit closed form for `|Crit_k|` in
arithmetic functions of N distinguishes prime-poset from random.

### D18 — Sudan-Guruswami list decoding of χ_P as agreement with low-degree polynomials over F_p

CLOSED_PATHS line 474 closes Goldreich-Levin Hadamard list decoding
of `pi(x)` over F_2 (no low-bit Fourier mass, comm rank `2^{N/2}`).
The Sudan-Guruswami algorithm is structurally different (algebraic
Reed-Solomon over `F_p`, polynomial-agreement target, not
Hadamard-character correlation) and operates on a different target
(`chi_P` indicator, not the cumulative `pi(x)`). Concrete first
step: for `p ∈ {1009, 10007, 100003}` view `(chi_P(0), ...,
chi_P(p-1))` as RS received word; sweep target degrees
`k ∈ {2, 4, 8, 16, log_2 p}`; run interpolation-and-factor algorithm
(bivariate `Q(X, Y)` of `(1, k)`-weighted degree `<= sqrt{2(k-1)p}`);
report empirical `agree(f_*) / p`. A-grade if `k_*(p) = O(polylog p)`
with `agree(f_*) >= 1/2 + 1/sqrt{log p}` (positive non-Hadamard
polynomial correlator — opens polylog frontier); B-grade if rigorous
list-decoding lower bound `k_*(p) = Omega(p / log^c p)` (negative-shape
edge complementary to line 474).

### D19 — Mapper algorithm topological-network of χ_P with arithmetic lens functions

D2 (S96, E2.17) closed metric Vietoris-Rips PH on Takens-embedded
prime gaps as a structural HL-detection measure. Mapper is
structurally distinct: the topology is induced by a LENS function
`f: X -> R^2` and clustering on preimages, not by Euclidean
distance — arithmetic structure can live in `f` (e.g., mod-q residues,
factorisation invariants). Concrete first step: for `M = 5000`
consecutive primes near `p ~ 10^6`, build the 3-vector point cloud,
choose lens `f(x) = (x[0] mod 6, log mean(x))`, cover `[0,6) x
[0.5, 5]` with 5×5 overlapping rectangles at 30% overlap, run
KeplerMapper. Report graph statistics (Betti number, components,
max-degree); compare to B1 (Cramér Poisson), B2 (gap-permuted), B3
(HL-uniform-residue resampled). A-grade if topological feature
deviates from B1/B2/B3 by >5σ AND has closed-form non-pair-correlation
arithmetic interpretation; B-grade if Mapper graph distinguishes
prime gaps from B2 with HL-residue match (sixth HL-detection family
member in graph-topological category).

## Cross-domain literature consulted (with URLs)

Per CLAUDE.md "Cross-Domain Imports" requirement, every new vector
cites at least one survey URL. WebFetch attempts (5 made):

- WebFetched https://emis.de/journals/SLC/wpapers/s48forman.pdf
  → 301 redirect to zbmath; primary content not retrieved this run,
  but Forman 2002 SLC 48 is the established survey.
- WebFetched https://en.wikipedia.org/wiki/Discrete_Morse_theory
  → confirmed: Morse function, Morse matching, critical cells, Morse
  inequalities, Forman 1998/2002 references.
- WebFetched https://www.cs.cmu.edu/~venkatg/pubs/papers/listdecoding-NOW.pdf
  → binary content (PDF body unparsed), but the file is the Guruswami
  2004 textbook chapter — established source.
- WebFetched https://en.wikipedia.org/wiki/Guruswami%E2%80%93Sudan_list_decoding_algorithm
  → confirmed: bivariate interpolation + factoring, Johnson bound
  `t = n - sqrt{kn}`, decoding to `1 - sqrt{R}` error fraction in
  poly time.
- WebFetched https://en.wikipedia.org/wiki/Oseledets_theorem
  → confirmed: cocycle definition, Lyapunov spectrum existence,
  Oseledec splitting, Oseledec 1965/1968.
- WebFetched https://en.wikipedia.org/wiki/Mapper_(topological_data_analysis)
  → confirmed: data + lens + cover + clusterer construction; nerve
  output `N(f^{-1}(U))`; Singh-Mémoli-Carlsson 2007 reference.
- (Borot-Eynard 2011 https://arxiv.org/abs/1011.1418 fetched but
  the page returned only the abstract — used to deselect "stochastic
  loop equations" as the fourth vector for this session in favour
  of Mapper, where the algorithmic recipe is concrete.)

## Self-grade: B

Rationale, honestly per CLAUDE.md "self-grade DOWN, not up":

- The four vectors propose techniques each UNUSED in the project
  registry. Each has a concrete single-session first step, a
  pre-stated falsification criterion, an expected failure mode (E /
  I / INC), and explicit A-grade vs B-grade success criteria. Each
  cites at least one survey URL, drawing on canonical published
  literature.
- Each vector explicitly distinguishes itself from the closest
  CLOSED_PATHS row (D16 vs line 422 transfer-operator; D17 vs D14
  + D2/E2.17; D18 vs line 474 Goldreich-Levin; D19 vs D2/E2.17 + D14).
- Honest grade-down considerations:
  - D17 risks an "primes are forced 0-critical" trivial closure (the
    failure-profile E text covers this explicitly), which makes the
    A-grade probability moderate.
  - D18 builds on a known empirical fact that χ_P has near-uniform
    HL-residue structure — list-decoding is unlikely to find polylog
    correlators where Goldreich-Levin found none, but the test is
    the right one to make the closure airtight.
  - D19 is the most "discovery-oriented" with weakest a-priori
    A-grade probability (Mapper outputs are highly cover-sensitive).
  - Only D16 has a strong path to genuine A-grade (sub-leading
    Lyapunov exponents are unmeasured anywhere in the literature for
    arithmetic potentials; the Bennetin algorithm is well-established).
- Quality > quantity: I picked 4, not 5. Considered "stochastic loop
  equations" (Borot-Eynard 2011) as a fifth but the survey didn't
  give me enough algorithmic concreteness for a single-session first
  step — pivoted to Mapper.
- Not A-grade: the session produced NO mathematical artefact (no
  theorem, no construction, no measurement) — only frontier targets.
  Per CLAUDE.md ATTACK_VECTORS rules, frontier_gen sessions are
  scored on the QUALITY of vectors produced, with B-grade being the
  natural-state output ("at least one is fresh"). At least 3 of the 4
  are clearly fresh and concrete; D16 has plausible A-grade reach.
  This is solid B; I do not claim A.

Per CLAUDE.md frontier_gen grading rubric ("A if vectors are
paper-grade fresh and you expect ≥2 to produce A-grade work; B if
at least one is fresh; C if all are minor variations; F if you
proposed nothing or duplicates"), I expect at most 1-2 of the 4 to
produce A-grade work, primarily D16. **B-grade**.

## Next-action recommendations

For the next production-mode session pulling from ATTACK_VECTORS:
- D16 has the highest A-grade reach and a concrete `N = 10^5` single-
  session first step using existing Anderson-localisation
  infrastructure (just swap the scalar transfer matrix for a 3-step
  product + Bennetin QR). Recommended priority pick.
- D17 is independently single-session at `N = 1024` and has no
  preconditions. Could fire in parallel.
- D18 needs a Sudan-Guruswami implementation (or sage's
  `coding.bivariate.SudanDecoder`); 1-2 sessions if the implementation
  exists, +1 if scratch.
- D19 needs KeplerMapper or scratch Mapper code; 1 session.

## Maintenance

Updated CROSS_DOMAIN_TECHNIQUES.md:
- §1 row "Multiplicative ergodic theorem (Oseledec)": UNUSED → PROPOSED (D16)
- §4 row "Mapper algorithm": UNUSED → PROPOSED (D19)
- §4 row "Discrete Morse theory": UNUSED → PROPOSED (D17)
- §6 rows "Reed-Solomon decoding" and "List decoding (Sudan,
  Guruswami)": UNUSED → PROPOSED (D18)
- Priority hints section: appended D16, D17, D18, D19 with one-line
  hints.

Updated ATTACK_VECTORS.md:
- §D: appended D16, D17, D18, D19 (each ~120-150 lines, with
  question / why frontier / cross-domain ingredient / distinction
  from existing closures / concrete first step / failure profile
  E,I,INC / A-grade success / B-grade success / cross-domain refs /
  budget).

No new edges or experiments produced this session (frontier_gen mode
is registry maintenance + target generation, not measurement).

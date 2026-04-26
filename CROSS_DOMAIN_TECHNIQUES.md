# Cross-Domain Technique Registry

This file is the source pool for `frontier_gen` mode. It tracks
mathematical techniques outside {classical analytic NT, complexity
theory, basic algebra} that the project either has used (as evidence
for / against polylog π(x)) or has not yet imported.

When `frontier_gen` mode fires (because ATTACK_VECTORS is thin or the
A-grade scarcity warning triggers), the agent picks 3-5 UNUSED rows,
WebFetches a foundational survey, and proposes new ATTACK_VECTORS
entries grounded in those techniques.

When any session imports a new technique that is not on this list,
the session must append it here with a survey reference URL.

---

## Status legend
- **USED (mode E)**: imported and reduced to a known closure (E mode).
- **USED (mode I)**: imported and produced a structural negative result.
- **USED (mode A)**: imported and produced an A-grade positive result.
- **PARTIAL**: superficial use only — name borrowed, machinery not imported.
- **UNUSED**: never imported.
- **PROPOSED**: scheduled for an upcoming attack vector.

---

## §1. Spectral / Operator-Theoretic

| Technique | Status | Survey ref | Linked vector / edge |
|-----------|--------|------------|----------------------|
| Cayley graph spectra (abelian) | USED (E) | Babai 1979 | A3 / E7.12 |
| Szegedy quantum walk | USED (I) | Szegedy 2004 (quant-ph/0401053) | D4 / E7.13 |
| Selberg trace formula | UNUSED | Iwaniec, *Spectral Methods of Automorphic Forms* (AMS, 2002) | candidate |
| Random regular graph spectral gap (Friedman) | UNUSED | Friedman 2008 (Memoirs AMS 195) | candidate |
| Anderson localisation (random Schrödinger) | PROPOSED (C4) | Aizenman-Warzel 2015 textbook; Bourgain arXiv:math/0410079 | C4 (S85 frontier_gen) |
| Continuous-time quantum walks (CTQW) | PROPOSED (D5) | Childs 2009 PRL 102 = arXiv:0806.1972; Childs-Cleve-Deotto-Farhi-Gutmann-Spielman 2003 = arXiv:quant-ph/0209131 | D5 (S85 frontier_gen) |
| Hodge / Laplacian on simplicial complexes | UNUSED | Lim 2020 SIAM Review 62(3) | candidate |
| Connes operator-theoretic ζ | UNUSED | Connes 1999 Selecta Math 5 | F1 zone |

## §2. Algebraic / Geometric

| Technique | Status | Survey ref | Linked vector / edge |
|-----------|--------|------------|----------------------|
| Croot-Lev-Pach polynomial method | PROPOSED | Tao 2016 blog post + Croot-Lev-Pach 2017 Annals 185 | B1 |
| Slice rank | PROPOSED | Tao 2016 (same) | B1 |
| Tropical geometry of arithmetic varieties | UNUSED | Maclagan-Sturmfels 2015 textbook | candidate |
| Étale cohomology / Frobenius traces | UNUSED | Deligne, Weil II | B3 (sketch only) |
| Motivic integration | UNUSED | Cluckers-Loeser surveys | candidate |
| Berkovich spaces (non-archimedean) | UNUSED | Baker-Rumely 2010 textbook | candidate |
| Schubert calculus / Grassmannian | UNUSED | Eisenbud-Harris 2016 textbook | candidate |
| Adelic / idele class group | UNUSED | Tate's thesis | F1 zone |

## §3. Probabilistic / Random Matrix

| Technique | Status | Survey ref | Linked vector / edge |
|-----------|--------|------------|----------------------|
| GUE pair correlation | USED (E) | Montgomery 1973 | E3.1, E3.13 |
| Bogomolny-Keating arithmetic correction | USED (I, structural) | BK 1996 | E3.13, §C1 closed |
| Free probability / free cumulants | USED (I) | Mingo-Speicher 2017 | C2 / E2.1 |
| Marchenko-Pastur | USED (PARTIAL) | Mingo-Speicher 2017 | S82 spike-bulk split |
| Stein's method | UNUSED | Ross 2011 Probability Surveys 8 | candidate |
| Determinantal point processes | UNUSED | Hough et al. 2009 textbook | candidate |
| Concentration of measure | PARTIAL | Boucheron-Lugosi-Massart 2013 | mentioned only |
| Stochastic loop equations | UNUSED | Borot-Eynard 2011 | candidate |

## §4. Topological / Geometric

| Technique | Status | Survey ref | Linked vector / edge |
|-----------|--------|------------|----------------------|
| Persistent homology | PROPOSED | Carlsson 2009 BAMS 46 | D2 |
| Mapper algorithm | UNUSED | Singh-Mémoli-Carlsson 2007 | candidate |
| Reeb graphs | UNUSED | Edelsbrunner-Harer 2010 textbook | candidate |
| Sheaf cohomology of sequences | UNUSED | Curry 2014 thesis | candidate |
| Discrete Morse theory | UNUSED | Forman 2002 | candidate |

## §5. Dynamical / Ergodic

| Technique | Status | Survey ref | Linked vector / edge |
|-----------|--------|------------|----------------------|
| Furstenberg correspondence | UNUSED | Furstenberg 1981 textbook | D1 candidate |
| Transfer operator spectrum | UNUSED | Baladi 2018 textbook | D1 |
| Mahler-style measures and ergodic averages | UNUSED | Everest-Ward 1999 textbook | candidate |
| Symbolic dynamics on prime sequences | UNUSED | Lind-Marcus 1995 textbook | candidate |
| Dynamical zeta functions (Ruelle / Selberg) | UNUSED | Ruelle 1976 | D1 |
| Multiplicative ergodic theorem (Oseledec) | UNUSED | Viana 2014 textbook | candidate |

## §6. Information-Theoretic / Coding

| Technique | Status | Survey ref | Linked vector / edge |
|-----------|--------|------------|----------------------|
| h_2 entropy of π(X)/X mod m | USED (A) | this project | E1.5 |
| LZ complexity of π(X) mod m | USED (E) | standard | CLOSED_PATHS |
| Reed-Solomon decoding | UNUSED | MacWilliams-Sloane 1977 textbook | candidate |
| List decoding (Sudan, Guruswami) | UNUSED | Guruswami textbook | candidate |
| Locally testable codes | UNUSED | Goldreich textbook | candidate |
| Compressed sensing on χ_P | UNUSED | Foucart-Rauhut 2013 textbook | candidate |

## §7. Combinatorial / Additive

| Technique | Status | Survey ref | Linked vector / edge |
|-----------|--------|------------|----------------------|
| Sum-product theorems | UNUSED | Tao-Vu 2006 textbook | candidate |
| Decoupling (Bourgain-Demeter) | UNUSED | Demeter 2020 textbook | candidate |
| Restriction theory (Stein-Tomas) | UNUSED | Tao 2003 lecture notes | candidate |
| Gowers norms (U^k uniformity) | USED E (S85, edge E2.13) | Green-Tao arXiv:math/0606088; Green-Tao arXiv:0807.1736 (Mobius nilsequences); Green-Tao-Ziegler arXiv:1009.3998 (U^{s+1} inverse); Hardy-Littlewood 1923 | D6 (CLOSED: Q^2(chi_P)→2.30=S_2 box singular series at N=2^18; Q^3 stable at 35.5; W-trick restores Gowers uniformity at U^2 within 0.1%; 36th pseudorandomness measure with closed-form HL prediction) |
| Bounded gaps / GPY sieve | UNUSED | Maynard 2015 Annals | candidate |
| Density Hales-Jewett | UNUSED | Polymath1 2010 | candidate |

## §8. Computational / Algorithmic

| Technique | Status | Survey ref | Linked vector / edge |
|-----------|--------|------------|----------------------|
| Tensor train (TT) | USED (I) | Oseledets 2011 | E2.1 / N1 |
| MPS / Hierarchical Tucker | USED (I) | Hackbusch 2012 textbook | E2.1 / N1 |
| MERA | USED (I) | Vidal 2007 PRL | N1 |
| Holographic / shadow tomography | UNUSED | Aaronson 2018 | candidate |
| Differentiable programming for π(x) | UNUSED | Bradbury et al. JAX docs | candidate |
| Boolean SAT for circuit synthesis | PROPOSED | Knuth Vol 4F | A1 |

## §9. Logic / Foundations

| Technique | Status | Survey ref | Linked vector / edge |
|-----------|--------|------------|----------------------|
| Lean 4 + mathlib4 formalisation | USED (PARTIAL) | mathlib docs | L1 (in progress) |
| Reverse mathematics (Friedman, Simpson) | UNUSED | Simpson 2009 textbook | candidate |
| Bounded arithmetic (Buss S^i_2 / Cook-Nguyen VTC^0) | PROPOSED (A4) | Buss 1986 thesis; Cook-Nguyen 2010 "Logical Foundations of Proof Complexity"; Krajicek 1995 textbook | A4 (S85 frontier_gen) |
| Proof complexity lower bounds | UNUSED | Krajicek 2019 textbook | candidate |

## §10. Frontier / Speculative

| Technique | Status | Survey ref | Linked vector / edge |
|-----------|--------|------------|----------------------|
| Voronin universality of ζ | PROPOSED (B4) | Steuding 2007 LNM 1877; Garunkstis 2003 Acta Arith. 113 = arXiv:math/0306072 | B4 (S85 frontier_gen) |
| Beurling generalised primes | UNUSED | Hilberdink-Lapidus 2016 | candidate |
| Riemann xi-function model spaces | UNUSED | de Branges 2004 unpublished + Lagarias surveys | candidate |
| Newman-Pollard hypothesis | UNUSED | Newman 1976 PAMS 51 | candidate |

---

## Maintenance rules

1. When `frontier_gen` mode picks an UNUSED technique to propose, mark
   it PROPOSED with the new ATTACK_VECTORS entry ID.
2. When a session attempts a PROPOSED technique, mark it USED with the
   correct mode (E / I / A) and the resulting edge ID.
3. When any session imports a technique not in this file, append a row
   to the appropriate section with a survey reference. Do NOT inflate
   superficial uses to "USED (A)"; if the import was just a name, mark
   it PARTIAL.
4. Do not delete rows. Status changes should be visible in the file
   history. The "graveyard" of USED-and-failed techniques is itself
   useful information.

---

## Priority hints for `frontier_gen`

The following UNUSED techniques are flagged as **highest priority** by
prior-session analysis (cross-domain ingredients with concrete
single-session protocols):

- **Voronin universality of ζ** (§10) — every analytic function is
  approximable by ζ vertical translates. If π(x) admits a Voronin-style
  approximation, the partial-sum representation could collapse.
- **Continuous-time quantum walks** (§1) — flagged in S80's next-action.
  Different mixing phenomenology than Szegedy walks.
- **Stein's method** (§3) — a quantitative Gaussian-deviation framework
  for π(x) mod m fluctuations not yet tried.
- **Sheaf cohomology of sequences** (§4) — gives an obstruction-theoretic
  view of the prime sequence as a "topological" object.
- **Furstenberg correspondence** (§5) — links arithmetic recurrence to
  measure-preserving systems; a way to make the GUE statistics into an
  ergodic statement.
- **Bounded arithmetic S^i_2** (§9) — formal-logic-side complement to
  TC⁰ circuit lower bounds. If π(x) is provable in S^1_2, then it has
  polylog-time witness extraction.

When in doubt, pick from this list first.

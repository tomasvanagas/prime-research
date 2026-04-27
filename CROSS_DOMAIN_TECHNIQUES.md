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
| Anderson localisation (random Schrödinger) | USED E (S88 chi_P edge E2.14; S100 lambda edge E2.18) | Aizenman-Warzel 2015 (AMS GSM 168); Furstenberg-Kifer 1983 Israel J. Math. 46; Pastur-Figotin 1992; Bourgain arXiv:math/0410079 | C4 (CLOSED: gamma_prime - gamma_bern cascade 88σ→4σ under W-trick W∈{2,6,30,210,2310}; mod-3 resonance peak at E~+1 isolated; second independent confirmation after E2.13 that chi_P deviation IS spectral signature of HL equidistribution mod q); G1 (CLOSED: gamma_lambda - gamma_Rademacher matched within seed noise across N up to 10^6 *without W-trick*; Liouville is spectrally featureless, isolating chi_P's deviation as exclusively HL mod-q) |
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
| Mahler measure / Lehmer's conjecture / Weil height of generating polynomial | PROPOSED (D10) | Smyth 2008 "The Mahler measure of algebraic numbers: a survey" CUP = https://homepages.ed.ac.uk/cjsmyth/papers/107.pdf; Boyd 1981 Canad. Math. Bull. 24; Lehmer 1933 Ann. Math. 34 | D10 (S97 frontier_gen) |

## §3. Probabilistic / Random Matrix

| Technique | Status | Survey ref | Linked vector / edge |
|-----------|--------|------------|----------------------|
| GUE pair correlation | USED (E) | Montgomery 1973 | E3.1, E3.13 |
| Bogomolny-Keating arithmetic correction | USED (I, structural) | BK 1996 | E3.13, §C1 closed |
| Free probability / free cumulants | USED (I) | Mingo-Speicher 2017 | C2 / E2.1 |
| Marchenko-Pastur | USED (PARTIAL) | Mingo-Speicher 2017 | S82 spike-bulk split |
| Stein's method | PROPOSED (C5) | Ross 2011 Probability Surveys 8; Chen-Goldstein-Shao 2011; Stein 1986 IMS Lect Notes 7 | C5 (S91 frontier_gen) |
| Determinantal point processes (incl. PPP, signed-K, complex-Hermitian) | USED I (S95, edge E2.16) | Hough-Krishnapur-Peres-Virag 2009 AMS ULect 51; Soshnikov 2000 Russ. Math. Surv. 55 = arXiv:math/0002099; Vere-Jones 1997 (alpha-DPP) | D7 (CLOSED: K^2_DPP < 0 at all 15 admissible even t and K^2_PPP < 0 at all 14 odd t > 1; 3-point PPP overshoots HL by up to 79.16% on 3-AP triples; signed-real-K sigma_req in (-0.541, 0.769) never ±1; complex Hermitian phase LS fit best residual 0.0746) |
| Concentration of measure | PARTIAL | Boucheron-Lugosi-Massart 2013 | mentioned only |
| Stochastic loop equations | UNUSED | Borot-Eynard 2011 | candidate |
| Möbius / nilsequence orthogonality (Sarnak conjecture, GT 2012, Tao logarithmic-Chowla) | USED E (S100, edge E2.18) | **Green-Tao 2012 *Annals* 175 = arXiv:0807.1736** (Möbius nilsequence orthogonality); **Sarnak 2010 IAS lectures** (Möbius randomness conjecture); **Tao 2016 Forum Math Pi 4 e8 = arXiv:1509.05422** (logarithmic-Chowla for two-point); Tao-Teräväinen 2017 arXiv:1708.02610 | G1 (CLOSED: gamma_lambda matches Rademacher Lyapunov within 2.2σ across 51 energies and 3 orders of N; Σ Chowla-z² = 4.77 vs χ²_16 = 16; first non-W-tricked spectral measure at noise floor) |

## §4. Topological / Geometric

| Technique | Status | Survey ref | Linked vector / edge |
|-----------|--------|------------|----------------------|
| Persistent homology | USED I (S96, edge E2.17) | Carlsson 2009 BAMS 46(2); Edelsbrunner-Harer 2010 textbook; Bauer 2021 *J. Appl. Comput. Topol.* 5 = arXiv:1908.02518; Cohen-Steiner-Edelsbrunner-Harer 2007 (PH stability) | D2 (CLOSED: total H_0 and H_1 persistence of Takens-embedded normalised prime gaps deviate ≥ 5σ from BOTH IID Exp(1) AND gap-permuted control across d ∈ {2,3,4} and M ∈ {500..4000}; primes more clustered + fewer loops than Poisson; fifth orthogonal HL-detection category after Gowers/Anderson/AlgImmunity/DPP) |
| Mapper algorithm | UNUSED | Singh-Mémoli-Carlsson 2007 | candidate |
| Reeb graphs | UNUSED | Edelsbrunner-Harer 2010 textbook | candidate |
| Cellular sheaf cohomology over divisibility poset (Curry framework) | PROPOSED (D14) | Curry 2014 PhD thesis = arXiv:1303.3255; Hansen-Ghrist 2019 *J. Appl. Comput. Topol.* 3; Robinson 2014 *Topological Signal Processing* Springer | D14 (S103 frontier_gen) |
| Sheaf cohomology of sequences (general — non-cellular) | UNUSED | Curry 2014 thesis | candidate |
| Discrete Morse theory | UNUSED | Forman 2002 | candidate |

## §5. Dynamical / Ergodic

| Technique | Status | Survey ref | Linked vector / edge |
|-----------|--------|------------|----------------------|
| Furstenberg correspondence | UNUSED | Furstenberg 1981 textbook | D1 candidate |
| Transfer operator spectrum | UNUSED | Baladi 2018 textbook | D1 |
| Mahler-style measures and ergodic averages | UNUSED | Everest-Ward 1999 textbook | candidate |
| Symbolic dynamics on prime sequences (subword complexity / topological entropy of χ_P) | USED E (S104, edge E2.19) | Lind-Marcus 1995 *An Introduction to Symbolic Dynamics and Coding* CUP; Cassaigne-Nicolas 2010 "Factor complexity" in CANT vol. 135; Morse-Hedlund 1938 Amer. J. Math. 60 | D13 (CLOSED: max \|z_perm\| cascade RAW 132.7 → ODD 277.1 → W6 120.5 → W30 24.8 → W210 8.4; W=210 W-trick erases the deficit to 8.4σ at L = 2.4·10^4; 38th pseudorandomness measure, sixth orthogonal HL-detection family in symbolic-dynamics / factor-complexity / topological-entropy category; subword complexity of chi_P never published) |
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
| Compressed sensing on χ_P (RIP / OMP / L1-min in structured AP dictionary) | PROPOSED (D12) | Candes-Tao 2006 IEEE Trans. Inf. Theory 51(12); Foucart-Rauhut 2013 *A Mathematical Introduction to Compressive Sensing* Birkhäuser | D12 (S97 frontier_gen) |

## §7. Combinatorial / Additive

| Technique | Status | Survey ref | Linked vector / edge |
|-----------|--------|------------|----------------------|
| Sum-product theorems (Bourgain-Glibichuk-Konyagin in F_p) | PROPOSED (D9) | Bourgain-Glibichuk-Konyagin 2006 J. London Math. Soc. = https://londmathsoc.onlinelibrary.wiley.com/doi/abs/10.1112/S0024610706022721; Garaev 2008 arXiv:math/0702780; Iosevich 2009 arXiv:0904.2075; Tao-Vu 2006 *Additive Combinatorics* CUP | D9 (S97 frontier_gen) |
| Decoupling (Bourgain-Demeter-Guth) | PROPOSED (D15) | Bourgain-Demeter 2015 *Annals* 182 = arXiv:1403.5335; Bourgain-Demeter-Guth 2015 *Annals* 184 = arXiv:1512.02384; Demeter 2020 *Fourier Restriction, Decoupling and Applications* CUP | D15 (S103 frontier_gen) |
| Restriction theory (Stein-Tomas) | UNUSED | Tao 2003 lecture notes | candidate |
| Gowers norms (U^k uniformity) | USED E (S85, edge E2.13) | Green-Tao arXiv:math/0606088; Green-Tao arXiv:0807.1736 (Mobius nilsequences); Green-Tao-Ziegler arXiv:1009.3998 (U^{s+1} inverse); Hardy-Littlewood 1923 | D6 (CLOSED: Q^2(chi_P)→2.30=S_2 box singular series at N=2^18; Q^3 stable at 35.5; W-trick restores Gowers uniformity at U^2 within 0.1%; 36th pseudorandomness measure with closed-form HL prediction) |
| Bounded gaps / GPY sieve (Maynard multidim. weights) | USED E (S116, edge E7.14) | Maynard 2015 Annals 181 = arXiv:1311.4600; Polymath8b arXiv:1407.4897 | A5 (CLOSED: aggregate-not-pointwise + divisor-enumeration cost ⇒ not TC⁰; 4th sieve-route family closure complementing E7.10 / E5.8 / E7.11) |
| Density Hales-Jewett | UNUSED | Polymath1 2010 | candidate |

## §8. Computational / Algorithmic

| Technique | Status | Survey ref | Linked vector / edge |
|-----------|--------|------------|----------------------|
| Tensor train (TT) | USED (I) | Oseledets 2011 | E2.1 / N1 |
| MPS / Hierarchical Tucker | USED (I) | Hackbusch 2012 textbook | E2.1 / N1 |
| MERA | USED (I) | Vidal 2007 PRL | N1 |
| Holographic / shadow tomography (Aaronson; classical shadows Huang-Kueng-Preskill) | PROPOSED (D11) | Aaronson 2018 STOC = arXiv:1711.01053; Huang-Kueng-Preskill 2020 Nature Physics 16, 1050 = arXiv:2002.08953 | D11 (S97 frontier_gen) |
| Differentiable programming for π(x) (DARTS architecture search) | PROPOSED (D8) | Liu-Simonyan-Yang 2019 ICLR = arXiv:1806.09055; Bender et al. 2018 ICML; Baydin et al. 2018 JMLR 18 | D8 (S91 frontier_gen) |
| Boolean SAT for circuit synthesis | PROPOSED | Knuth Vol 4F | A1 |

## §9. Logic / Foundations

| Technique | Status | Survey ref | Linked vector / edge |
|-----------|--------|------------|----------------------|
| Lean 4 + mathlib4 formalisation | USED (PARTIAL) | mathlib docs | L1 (in progress) |
| Reverse mathematics — quantitative π(x) bounds (Friedman, Simpson, Avigad) | PROPOSED (A6) | Simpson 2009 *Subsystems of Second-Order Arithmetic* CUP; Avigad-Friedman 1998 "Combining decision procedures for the reals" https://www.andrew.cmu.edu/user/avigad/Papers/dpr.pdf; Avigad-Towsner 2012 *J. Symbolic Logic* 77; Friedman-Sieg 1989 LNM 1320 | A6 (S103 frontier_gen) |
| Reverse mathematics — qualitative existence (CLOSED line 218 CLOSED_PATHS) | USED (E) | Simpson 2009 textbook | CLOSED line 218 (S7) |
| Bounded arithmetic (Buss S^i_2 / Cook-Nguyen VTC^0) | PROPOSED (A4) | Buss 1986 thesis; Cook-Nguyen 2010 "Logical Foundations of Proof Complexity"; Krajicek 1995 textbook | A4 (S85 frontier_gen) |
| Proof complexity lower bounds | UNUSED | Krajicek 2019 textbook | candidate |

## §10. Frontier / Speculative

| Technique | Status | Survey ref | Linked vector / edge |
|-----------|--------|------------|----------------------|
| Voronin universality of ζ | PROPOSED (B4) | Steuding 2007 LNM 1877; Garunkstis 2003 Acta Arith. 113 = arXiv:math/0306072 | B4 (S85 frontier_gen) |
| Beurling generalised primes (perturbative deformation family) | PROPOSED (B5) | Beurling 1937 *Acta Math.* 68; Diamond-Zhang 2007 *Trans. AMS* 359; Hilberdink 2009 *J. Number Theory* 129; Hilberdink-Lapidus 2016 *Functiones et Approximatio* 54 = https://projecteuclid.org/euclid.facm/1428589969 | B5 (S103 frontier_gen) |
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
  for π(x) mod m fluctuations. **NOW PROPOSED as C5 (S91).**
- **Sheaf cohomology of sequences** (§4) — gives an obstruction-theoretic
  view of the prime sequence as a "topological" object. **Persistent
  homology now CLOSED S96 mode I, edge E2.17** — the metric-space
  topological line is exhausted within the project's HL-detection family;
  sheaf cohomology of sequences remains an UNUSED orthogonal route.
  **NOW PROPOSED as D14 (S103)** via the cellular-sheaf framework
  (Curry 2014), structurally distinct from the Spec(Z) sheaf line
  CLOSED in CLOSED_PATHS line 202.
- **Furstenberg correspondence** (§5) — links arithmetic recurrence to
  measure-preserving systems; a way to make the GUE statistics into an
  ergodic statement. (CAVEAT: classical form closed line 421 of
  CLOSED_PATHS — circular under spectral=zero-zero duality. Needs
  carefully restricted sub-formulation to avoid duplicate.)
- **Bounded arithmetic S^i_2** (§9) — formal-logic-side complement to
  TC⁰ circuit lower bounds. **NOW PROPOSED as A4 (S85).**
- **Determinantal point processes** (§3) — fit a kernel to the integer
  prime sequence; tests whether primes are a DPP. **CLOSED S95, mode I,
  edge E2.16: primes are NOT a translation-invariant DPP / PPP /
  signed-real-K / complex-Hermitian-K point process; HL prime-factor
  factorisation is genuinely non-pairwise.**
- **Maynard multidimensional sieve weight** (§7) — direct test whether
  the most refined explicit prime-detection sieve is TC⁰-computable at
  single-n granularity. **NOW PROPOSED as A5 (S91).**
- **Differentiable architecture search (DARTS)** (§8) — gradient-based
  bridge over circuit-synthesis SAT barrier from S84. **NOW PROPOSED
  as D8 (S91).**

When in doubt, pick from this list first.

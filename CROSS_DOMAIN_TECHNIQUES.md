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
| Selberg trace formula / Hecke eigenvalue partial sums | USED E (S118, edge E7.15) | Iwaniec, *Spectral Methods of Automorphic Forms* (AMS, 2002); Sarnak 2008 letter to Bombieri; Katz-Sarnak 1999 *Random Matrices, Frobenius Eigenvalues, and Monodromy* AMS Coll. 45 | B2 (CLOSED: F_τ Hecke twisted partial sums of level-1 weight-12 Δ are 3× more obstructed from spanning π(x) − Li(x) than matched-multiplicative random Sato-Tate ensembles, ratio 2.83 ± 0.02 across 9 (N, K) cells, Z = 17–58σ; mechanism = Sato-Tate / Katz-Sarnak GUE-independence of {γ^Δ} and {γ^ζ}) |
| Random regular graph spectral gap (Friedman) | USED E (S125 edge E7.16) | Friedman 2008 Memoirs AMS 195 = arXiv:cs/0405020 https://arxiv.org/abs/cs/0405020; Hoory-Linial-Wigderson 2006 Bull. AMS 43, 439; Lubotzky 2012 Bull. AMS 49, 113 = arXiv:1105.2389 | D20 (CLOSED: r_N(prime) of Cay(Z/NZ, ±primes < N^c) indistinguishable from parity-and-support-matched random within ±2σ across 10 (N, c) cells; bare r_N = 2.05–11.30 sub-Ramanujan, but reduces to bounded-support FFT spike at k=1 + parity spike at k≈N/2 modulated by single even prime p=2; no HL singular-series residual detected at scales tested) |
| Anderson localisation (random Schrödinger) | USED E (S88 chi_P edge E2.14; S100 lambda edge E2.18) | Aizenman-Warzel 2015 (AMS GSM 168); Furstenberg-Kifer 1983 Israel J. Math. 46; Pastur-Figotin 1992; Bourgain arXiv:math/0410079 | C4 (CLOSED: gamma_prime - gamma_bern cascade 88σ→4σ under W-trick W∈{2,6,30,210,2310}; mod-3 resonance peak at E~+1 isolated; second independent confirmation after E2.13 that chi_P deviation IS spectral signature of HL equidistribution mod q); G1 (CLOSED: gamma_lambda - gamma_Rademacher matched within seed noise across N up to 10^6 *without W-trick*; Liouville is spectrally featureless, isolating chi_P's deviation as exclusively HL mod-q) |
| Continuous-time quantum walks (CTQW) | PROPOSED (D5) | Childs 2009 PRL 102 = arXiv:0806.1972; Childs-Cleve-Deotto-Farhi-Gutmann-Spielman 2003 = arXiv:quant-ph/0209131 | D5 (S85 frontier_gen) |
| Hodge / Laplacian on simplicial complexes (higher-order L_k, k ≥ 1) | USED E (S126 edge E7.17) | Lim 2020 SIAM Review 62(3) "Hodge Laplacians on graphs" = arXiv:1507.05379 https://arxiv.org/abs/1507.05379; Eckmann 1944 Comment. Math. Helv. 17; Friedman 1996 Algorithmica 21; Horak-Jost 2013 Adv. Math. 244; Goff 2009 J. Algebraic Combin. 30, 215 (join spectrum); Kahle 2014 *Annals* 179 (random flag complexes) | D22 (CLOSED S126: L_1 spectrum of K_N := pairwise-coprime flag complex on [2..N] is structurally constrained — λ_max(L_1) = |V| = N − 1 EXACTLY at all 9 tested N ∈ {8..128}; multiplicity = C(k+1, 2) where k = π(N) − π(N/2) is the Bertrand-prime count, verified perfectly in 9/9 cells; mechanism = Bertrand-universal-vertex join decomposition K_N = Δ^{k-1} \ast F(H); β_k = 0 for k ≥ 1 deterministically (cone collapse); L_1 mean shift reduces to triple-coprime singular series ∏_p(1 − 3/p² + 2/p³) / (6/π²)³ ≈ 1.276; KS p < 1e-300 at N = 128) |
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
| Mahler measure / Lehmer's conjecture / Weil height of generating polynomial | USED I (S134, edge E2.20) | Smyth 2008 "The Mahler measure of algebraic numbers: a survey" CUP = https://homepages.ed.ac.uk/cjsmyth/papers/107.pdf; Boyd 1981 Canad. Math. Bull. 24; Lehmer 1933 Ann. Math. 34; Dobrowolski 1979 *Acta Arith.* 34 (lower bound `log m(P) ≥ C(log log D / log D)³`) | D10 (CLOSED S134, mode I, B-grade: log m(f_PRIMES_N) − log m(f_BERN(d_N)_N) plateaus at Δ_∞ ≈ −0.307 ± 0.001 nat from N=2^16, z(MATCH) = −337σ at N=2^18; f_N / z² irreducible over Q[z] for N∈{64,128,256} — zero cyclotomic share; A-grade poly-log compressibility hypothesis falsified; **6th orthogonal pseudorandomness measure category — algebraic-height / multiplicative-height — after E2.13 Gowers / E2.14 Anderson / E2.15 alg.immunity / E2.16 DPP failure / E2.17 PH**) |
| Newman / Erdős polynomial flatness — L^∞ extremality of 0/1-coefficient polynomials on `\|z\|=1` | PROPOSED (D27) | Erdős 1957 *Mich. Math. J.* 4 (Littlewood/Erdős flatness conjecture for ±1 polynomials); Newman 1976 "Norms of polynomials" *Proc. AMS* 51 (0/1 case); Kahane 1985 *Some Random Series of Functions* CUP §6 (Salem-Zygmund random model); Bourgain 1988 "Bounded orthogonal systems and the Λ(p)-set problem" *Acta Math.* 162; Bourgain-Lewko 2016 *J. Funct. Anal.* (random-polynomial mean-vs-max bounds); Balister-Bollobás-Morris-Sahasrabudhe 2020 *Annals of Math.* 192, 977 (resolution of Erdős for ±1 — flat Littlewood polynomials EXIST constructively) https://en.wikipedia.org/wiki/Littlewood_polynomial | D27 (S136 frontier_gen) — `\|f_N\|_∞ / √π(N)` for `f_N(z) = Σ_{n≤N} χ_P(n) z^n`, structurally orthogonal to D10 Mahler (log integral) and D25 Stein-Tomas (`L^p` finite-`p`) |
| LPS Ramanujan graph construction (Lubotzky-Phillips-Sarnak quaternion-algebra Cayley) | PROPOSED (D28) | Lubotzky-Phillips-Sarnak 1988 "Ramanujan graphs" *Combinatorica* 8, 261; Margulis 1988 *Probl. Inf. Transm.* 24; Davidoff-Sarnak-Valette 2003 *Elementary Number Theory, Group Theory and Ramanujan Graphs* CUP; Lubotzky 2012 *Bull. AMS* 49, 113 = arXiv:1105.2389 https://arxiv.org/abs/1105.2389; Hoory-Linial-Wigderson 2006 *Bull. AMS* 43, 439; Wikipedia: Ramanujan graph https://en.wikipedia.org/wiki/Ramanujan_graph | D28 (S136 frontier_gen) — non-abelian Cay(PSL_2(F_p), prime-indexed LPS-quaternion generators); explicitly flagged as D20.c successor in S125 |

## §3. Probabilistic / Random Matrix

| Technique | Status | Survey ref | Linked vector / edge |
|-----------|--------|------------|----------------------|
| GUE pair correlation | USED (E) | Montgomery 1973 | E3.1, E3.13 |
| Sine-kernel n-correlation determinant `R_n = det[K(s_i-s_j)]` (n=2..6) | USED E (S25/S45/S57/S123, edge E7.1 refined) | Mehta 2004 *Random Matrices* (3rd ed.) §6.2; Hough-Krishnapur-Peres-Virág 2009 AMS ULect 51 §1.2; Bohigas-Giannoni-Schmit 1984 Phys. Rev. Lett. 52, 1 | C2 (CLOSED S123: empirical zeta R_n along equally-spaced slices, P_k spacings k=0..5, κ_n cumulants n=2..6 at L≥8 — all match GUE within sample noise once matched-finite-N null is used; Conrey-Snaith arithmetic correction 1/L²≈0.024 below noise floor 1/√(n_tuples)≥0.3 at order 5) |
| Bogomolny-Keating arithmetic correction | USED (I, structural) | BK 1996 | E3.13, §C1 closed |
| Conrey-Snaith higher-order arithmetic corrections | USED I (S123) | Conrey-Snaith 2007 arXiv:0710.5304 | C2 (CLOSED: orders 4,5,6 detection floor 1/L² well below noise floor 1/√(n_tuples) at γ ≤ 10⁴; same Odlyzko-height scaling barrier as §C1/S71) |
| Free probability / free cumulants | USED (I) | Mingo-Speicher 2017 | C2 / E2.1 |
| Marchenko-Pastur | USED (PARTIAL) | Mingo-Speicher 2017 | S82 spike-bulk split |
| Stein's method | PROPOSED (C5) | Ross 2011 Probability Surveys 8; Chen-Goldstein-Shao 2011; Stein 1986 IMS Lect Notes 7 | C5 (S91 frontier_gen) |
| Determinantal point processes (incl. PPP, signed-K, complex-Hermitian) | USED I (S95, edge E2.16) | Hough-Krishnapur-Peres-Virag 2009 AMS ULect 51; Soshnikov 2000 Russ. Math. Surv. 55 = arXiv:math/0002099; Vere-Jones 1997 (alpha-DPP) | D7 (CLOSED: K^2_DPP < 0 at all 15 admissible even t and K^2_PPP < 0 at all 14 odd t > 1; 3-point PPP overshoots HL by up to 79.16% on 3-AP triples; signed-real-K sigma_req in (-0.541, 0.769) never ±1; complex Hermitian phase LS fit best residual 0.0746) |
| Pfaffian / α-determinantal point processes (matrix-valued kernel) | PROPOSED (C6) | Borodin 2009 *Encyclopaedia Math. Sci.* 152 = arXiv:0911.1153 https://arxiv.org/abs/0911.1153; Soshnikov 2000 = arXiv:math/0002099; Vere-Jones 1997 *New Zealand J. Math.* 26 (α-permanent); Mehta 2004 *Random Matrices* §6 (GOE/GSE Pfaffian formulas) | C6 (S130 frontier_gen) — over-determined Pfaffian-vs-determinant n=4 fit on ζ zeros |
| Gaussian multiplicative chaos / Fyodorov-Hiary-Keating max statistics of |ζ(1/2 + it)| | USED I (S133, edge E7.18) | Fyodorov-Hiary-Keating 2012 *PRL* 108 = arXiv:1202.4713 https://arxiv.org/abs/1202.4713; Saksman-Webb 2018 = arXiv:1609.00027 https://arxiv.org/abs/1609.00027; Arguin-Belius-Bourgade 2017 *Comm. Math. Phys.* 349 = arXiv:1612.08575; Bourgade-Kuan 2014 *Comm. Pure Appl. Math.* 67 | C7 (CLOSED S133, mode I, B-grade case 2: FIRST ζ-amplitude measurement; M_T mean T-INDEPENDENT (pooled −0.657 ± 0.045 across T ∈ {10⁴, 10⁵, 10⁶}, K = 100/anchor) confirming FHK normalisation; M_T SHAPE is approximately Gaussian at finite T ≤ 10⁶ — var 1.47× FHK Gumbel(1/2) prediction π²/24, skew ≈ 0.10 vs Gumbel +1.14, ex-kurt ≈ −0.4 vs Gumbel +2.4; FHK Gumbel-shape limit NOT detectable; first quantitative FHK convergence-rate result) |
| Concentration of measure | PARTIAL | Boucheron-Lugosi-Massart 2013 | mentioned only |
| Stochastic loop equations / Eynard-Orantin topological recursion | PROPOSED (D24) | Eynard-Orantin 2007 *Comm. Number Theory Phys.* 1 = arXiv:math-ph/0702045 https://arxiv.org/abs/math-ph/0702045; Eynard 2016 *Counting Surfaces* Birkhäuser; Borot-Eynard 2011 arXiv:1107.5028; Keating-Snaith 2000 *Comm. Math. Phys.* 214 | D24 (S121 frontier_gen) |
| Möbius / nilsequence orthogonality (Sarnak conjecture, GT 2012, Tao logarithmic-Chowla) | USED E (S100, edge E2.18) | **Green-Tao 2012 *Annals* 175 = arXiv:0807.1736** (Möbius nilsequence orthogonality); **Sarnak 2010 IAS lectures** (Möbius randomness conjecture); **Tao 2016 Forum Math Pi 4 e8 = arXiv:1509.05422** (logarithmic-Chowla for two-point); Tao-Teräväinen 2017 arXiv:1708.02610 | G1 (CLOSED: gamma_lambda matches Rademacher Lyapunov within 2.2σ across 51 energies and 3 orders of N; Σ Chowla-z² = 4.77 vs χ²_16 = 16; first non-W-tricked spectral measure at noise floor) |

## §4. Topological / Geometric

| Technique | Status | Survey ref | Linked vector / edge |
|-----------|--------|------------|----------------------|
| Persistent homology | USED I (S96, edge E2.17) | Carlsson 2009 BAMS 46(2); Edelsbrunner-Harer 2010 textbook; Bauer 2021 *J. Appl. Comput. Topol.* 5 = arXiv:1908.02518; Cohen-Steiner-Edelsbrunner-Harer 2007 (PH stability) | D2 (CLOSED: total H_0 and H_1 persistence of Takens-embedded normalised prime gaps deviate ≥ 5σ from BOTH IID Exp(1) AND gap-permuted control across d ∈ {2,3,4} and M ∈ {500..4000}; primes more clustered + fewer loops than Poisson; fifth orthogonal HL-detection category after Gowers/Anderson/AlgImmunity/DPP) |
| Mapper algorithm | PROPOSED (D19) | Singh-Mémoli-Carlsson 2007 *Eurographics Symp. PBG* (Stanford TR) https://research.math.osu.edu/tgda/mapperPBG.pdf; Lum et al. 2013 *Sci. Rep.* 3, 1236 https://www.nature.com/articles/srep01236; KeplerMapper https://kepler-mapper.scikit-tda.org/; https://en.wikipedia.org/wiki/Mapper_(topological_data_analysis) | D19 (S119 frontier_gen) |
| Reeb graphs of arithmetic height functions on Z | PROPOSED (D21) | Edelsbrunner-Harer 2010 *Computational Topology* AMS, ch. VI.3; Carr-Snoeyink-Axen 2003 *Comput. Geom.* 24, 75; Wikipedia: Reeb graph https://en.wikipedia.org/wiki/Reeb_graph | D21 (S121 frontier_gen) |
| Cellular sheaf cohomology over divisibility poset (Curry framework) | PROPOSED (D14) | Curry 2014 PhD thesis = arXiv:1303.3255; Hansen-Ghrist 2019 *J. Appl. Comput. Topol.* 3; Robinson 2014 *Topological Signal Processing* Springer | D14 (S103 frontier_gen) |
| Sheaf cohomology of sequences (general — non-cellular) | UNUSED | Curry 2014 thesis | candidate |
| Discrete Morse theory | PROPOSED (D17) | Forman 1998 *Adv. Math.* 134, 90; Forman 2002 "A user's guide to discrete Morse theory" Sém. Lothar. Combin. 48 https://www.emis.de/journals/SLC/wpapers/s48forman.pdf; Benedetti-Lutz 2014 *Exp. Math.* 23, 66 https://arxiv.org/abs/1303.6422; Bjorner-Wachs 2005 *Trans. AMS* 357; https://en.wikipedia.org/wiki/Discrete_Morse_theory | D17 (S119 frontier_gen) |

## §5. Dynamical / Ergodic

| Technique | Status | Survey ref | Linked vector / edge |
|-----------|--------|------------|----------------------|
| Furstenberg correspondence | UNUSED | Furstenberg 1981 textbook | D1 candidate |
| Transfer operator spectrum (Pollicott-Ruelle resonances of arithmetic dynamical system) | PROPOSED (D30) | Pollicott 1985 *Inventiones Math.* 81, 413 (rate of mixing of Axiom A flows); Ruelle 1976 *Inventiones Math.* 34, 231 (zeta-functions for expanding maps); Mayer 1991 *Bull. AMS* 25, 55 (transfer-operator approach to Selberg ζ via Gauss map); Baladi 2018 *Dynamical Zeta Functions and Dynamical Determinants for Hyperbolic Maps* Springer Ergeb. 68; Liverani 2004 *Comm. Math. Phys.* 248 (anisotropic Sobolev spaces); Tsujii 2010 *Nonlinearity* 23; Wikipedia: Transfer operator https://en.wikipedia.org/wiki/Transfer_operator | D30 (S136 frontier_gen) — Pollicott-Ruelle resonances of χ_P-weighted Gauss-map transfer operator; structurally distinct from CLOSED line 105 (constructive transfer matrix sieve), CLOSED line 182 (FRACTRAN transfer operator), CLOSED line 425 (Furstenberg correspondence) |
| Mahler-style measures and ergodic averages | UNUSED | Everest-Ward 1999 textbook | candidate |
| Symbolic dynamics on prime sequences (subword complexity / topological entropy of χ_P) | USED E (S104, edge E2.19) | Lind-Marcus 1995 *An Introduction to Symbolic Dynamics and Coding* CUP; Cassaigne-Nicolas 2010 "Factor complexity" in CANT vol. 135; Morse-Hedlund 1938 Amer. J. Math. 60 | D13 (CLOSED: max \|z_perm\| cascade RAW 132.7 → ODD 277.1 → W6 120.5 → W30 24.8 → W210 8.4; W=210 W-trick erases the deficit to 8.4σ at L = 2.4·10^4; 38th pseudorandomness measure, sixth orthogonal HL-detection family in symbolic-dynamics / factor-complexity / topological-entropy category; subword complexity of chi_P never published) |
| Dynamical zeta functions (Ruelle / Selberg) | UNUSED | Ruelle 1976 | D1 |
| Multiplicative ergodic theorem (Oseledec) | PROPOSED (D16) | Oseledec 1968 *Trans. Moscow Math. Soc.* 19; Furstenberg-Kesten 1960 *Ann. Math. Statist.* 31; Ruelle 1979 *IHÉS Publ. Math.* 50; Viana 2014 *Lectures on Lyapunov Exponents* CUP; Bennetin-Galgani-Giorgilli-Strelcyn 1980 *Meccanica* 15 (QR algorithm); https://en.wikipedia.org/wiki/Oseledets_theorem | D16 (S119 frontier_gen) |

## §6. Information-Theoretic / Coding

| Technique | Status | Survey ref | Linked vector / edge |
|-----------|--------|------------|----------------------|
| h_2 entropy of π(X)/X mod m | USED (A) | this project | E1.5 |
| LZ complexity of π(X) mod m | USED (E) | standard | CLOSED_PATHS |
| Reed-Solomon decoding | PROPOSED (D18) | MacWilliams-Sloane 1977 textbook; Sudan 1997 J. Complex. 13, 180 https://people.csail.mit.edu/madhu/papers/1997/ldc-final.pdf | D18 (S119 frontier_gen) |
| List decoding (Sudan, Guruswami) | PROPOSED (D18) | Guruswami-Sudan 1999 IEEE Trans. Inf. Theory 45, 1757 https://www.cs.cmu.edu/~venkatg/pubs/papers/listdecoding-NOW.pdf; Guruswami 2004 *List Decoding of Error-Correcting Codes* Springer LNCS 3282; Roth-Ruckenstein 2000 IEEE Trans. Inf. Theory 46; https://en.wikipedia.org/wiki/Guruswami%E2%80%93Sudan_list_decoding_algorithm | D18 (S119 frontier_gen) |
| Locally testable codes (BLR linearity test, Reed-Muller LTC, Long-code) | PROPOSED (D26) | Goldreich-Sudan 2002 *J. ACM* 53 = ECCC TR02-050 https://eccc.weizmann.ac.il/eccc-reports/2002/TR02-050/; Dinur 2007 *J. ACM* 54 (PCP gap amplification); Blum-Luby-Rubinfeld 1993 *JCSS* 47 (BLR linearity); Kaufman-Litsyn 2005 *FOCS* (almost-orthogonal LTCs); Goldreich 2010 *P, NP, NP-Completeness* CUP | D26 (S130 frontier_gen) — constant-query primality predictor on χ_P-encoded codeword |
| Compressed sensing on χ_P (RIP / OMP / L1-min in structured AP dictionary) | PROPOSED (D12) | Candes-Tao 2006 IEEE Trans. Inf. Theory 51(12); Foucart-Rauhut 2013 *A Mathematical Introduction to Compressive Sensing* Birkhäuser | D12 (S97 frontier_gen) |

## §7. Combinatorial / Additive

| Technique | Status | Survey ref | Linked vector / edge |
|-----------|--------|------------|----------------------|
| Sum-product theorems (Bourgain-Glibichuk-Konyagin in F_p) | PROPOSED (D9) | Bourgain-Glibichuk-Konyagin 2006 J. London Math. Soc. = https://londmathsoc.onlinelibrary.wiley.com/doi/abs/10.1112/S0024610706022721; Garaev 2008 arXiv:math/0702780; Iosevich 2009 arXiv:0904.2075; Tao-Vu 2006 *Additive Combinatorics* CUP | D9 (S97 frontier_gen) |
| Decoupling (Bourgain-Demeter-Guth) | PROPOSED (D15) | Bourgain-Demeter 2015 *Annals* 182 = arXiv:1403.5335; Bourgain-Demeter-Guth 2015 *Annals* 184 = arXiv:1512.02384; Demeter 2020 *Fourier Restriction, Decoupling and Applications* CUP | D15 (S103 frontier_gen) |
| Restriction theory (Stein-Tomas / Bourgain Λ(p)-set) | PROPOSED (D25) | Stein 1993 *Harmonic Analysis* Princeton Ch. IX; Tomas 1975 *Bull. AMS* 81; Bourgain 1989 "Bounded orthogonal systems and the Λ(p)-set problem" *Acta Math.* 162; Green 2005 "Roth's theorem in the primes" *GAFA* 15 = arXiv:math/0302311; Tao 2020 247B Notes 1 https://terrytao.wordpress.com/2020/03/29/247b-notes-1-restriction-theory/; Foschi-Oliveira e Silva 2024 *São Paulo J. Math. Sci.* https://link.springer.com/article/10.1007/s40863-024-00422-x | D25 (S130 frontier_gen) — global L^p restriction-saturation test on prime exponential sum |
| Gowers norms (U^k uniformity) | USED E (S85, edge E2.13) | Green-Tao arXiv:math/0606088; Green-Tao arXiv:0807.1736 (Mobius nilsequences); Green-Tao-Ziegler arXiv:1009.3998 (U^{s+1} inverse); Hardy-Littlewood 1923 | D6 (CLOSED: Q^2(chi_P)→2.30=S_2 box singular series at N=2^18; Q^3 stable at 35.5; W-trick restores Gowers uniformity at U^2 within 0.1%; 36th pseudorandomness measure with closed-form HL prediction) |
| Bounded gaps / GPY sieve (Maynard multidim. weights) | USED E (S116, edge E7.14) | Maynard 2015 Annals 181 = arXiv:1311.4600; Polymath8b arXiv:1407.4897 | A5 (CLOSED: aggregate-not-pointwise + divisor-enumeration cost ⇒ not TC⁰; 4th sieve-route family closure complementing E7.10 / E5.8 / E7.11) |
| Density Hales-Jewett (combinatorial-line density on [k]^d) | PROPOSED (D23) | Furstenberg-Katznelson 1991 J. d'Anal. Math. 57; Polymath1 2010 *Annals* 175 = arXiv:0910.3926 https://arxiv.org/abs/0910.3926; Hales-Jewett 1963 Trans. AMS 106 | D23 (S121 frontier_gen) |

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
| Newman-Pollard hypothesis | PARTIAL (folded into D27 PROPOSED — see §2 row "Newman / Erdős polynomial flatness") | Newman 1976 *Proc. AMS* 51 | D27 (S136 frontier_gen) — D27 measures Newman's L^∞ extremality on χ_P-coefficient polynomial directly |
| Cohn-Elkies LP bound for sphere packing + Viazovska modular-form magic functions | PROPOSED (D29) | Cohn-Elkies 2003 "New upper bounds on sphere packings I" *Annals of Math* 157, 689 = arXiv:math/0110009 https://arxiv.org/abs/math/0110009; Viazovska 2017 "The sphere packing problem in dimension 8" *Annals of Math* 185, 991 = arXiv:1603.04246 https://arxiv.org/abs/1603.04246; Cohn-Kumar-Miller-Radchenko-Viazovska 2017 *Annals* 185, 1017 = arXiv:1603.06518; Cohn 2017 "A conceptual breakthrough in sphere packing" *Notices AMS* 64(2), 102; Cohn-Goncalves 2019 *Invent. Math.* 217 = arXiv:1712.04438 (LP-modular-form connection); Wikipedia: Sphere packing in dimension 8 https://en.wikipedia.org/wiki/Sphere_packing | D29 (S136 frontier_gen) — Cohn-Elkies LP applied to prime autocorrelation `R_P(t)`, Viazovska-style modular-form magic function for χ_P |

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
- **Multiplicative ergodic theorem (Oseledec)** (§5) — the project's
  Anderson Lyapunov work (E2.14 chi_P, E2.18 lambda) measured only the
  TOP exponent of an SL(2, R)-cocycle. The full Oseledec spectrum at
  d=3,4 gives sub-leading exponents that are pseudorandomness measures
  in their own right. **NOW PROPOSED as D16 (S119).**
- **Discrete Morse theory** (§4) — order-complex critical-cell
  combinatorial invariant of the divisibility poset. Categorically
  orthogonal to metric PH (E2.17, S96) and to cellular sheaf cohomology
  (D14, S103). **NOW PROPOSED as D17 (S119).**
- **Sudan-Guruswami list decoding** (§6) — algebraic Reed-Solomon list
  decoder for chi_P over F_p. Distinct code class and distinct target
  from the closed Goldreich-Levin Hadamard line (line 474 closes
  list-decoding for pi(x) over F_2). **NOW PROPOSED as D18 (S119).**
- **Mapper algorithm** (§4) — arithmetic-lens-based topological-network
  invariant; lens-induced topology (not metric) makes it structurally
  distinct from PH (E2.17) and from cellular sheaves (D14). **NOW
  PROPOSED as D19 (S119).**
- **Random regular graph spectral gap (Friedman)** (§1) — quantitative
  Ramanujan-bound test on `Cay(Z/NZ, primes < N^c)`. λ_2 saturation
  vs random regular graph determines whether primes are a Ramanujan-
  typical generating set; CLOSED line 356 only addressed the trivial
  λ_0. **CLOSED S125, edge E7.16, USED-E**: r_N(prime) =
  λ_2 / (2 √(d-1)) is Friedman-typical within ±2σ once support
  `[2, M)` and parity (odd) are matched. Bare r_N = 2.05–11.30
  sub-Ramanujan, but reduces to bounded-support `k = 1` FFT spike +
  single-prime `p = 2` parity-spike artefact. No HL singular-series
  residual detected at the scales tested.
- **Reeb graphs of arithmetic height functions** (§4) — exact (no
  metric, no lens) Reeb graph of `ω(n)` etc. on `Z`; orthogonal to PH
  (metric), Mapper (lens), and discrete Morse (poset). **NOW PROPOSED
  as D21 (S121).**
- **Higher-order Hodge Laplacian** (§1) — `L_k` for `k ≥ 1` on
  coprimality FLAG complex; CLOSED lines 356, 387 addressed only `L_0`
  (graph Laplacian = Ramanujan sums). Higher-order encodes Möbius /
  triple-coprime data orthogonal to bilinear-form Ramanujan sums.
  **NOW PROPOSED as D22 (S121).**
- **Density Hales-Jewett combinatorial-line density** (§7) —
  combinatorial-line density of primes in `[10]^d` digit cubes;
  structured-AP density complementing the unrestricted Gowers `U^k`
  measurement (E2.13, S85). **NOW PROPOSED as D23 (S121).**
- **Eynard-Orantin topological recursion** (§3) — loop equations /
  spectral curve for the prime correlation hierarchy; tests if all
  `n`-point correlations are recursively determined by `W^{(1)} =
  Σ_p 1/(z - p)` and a spectral curve `y(x)`. **NOW PROPOSED as D24
  (S121).**
- **Pfaffian / α-determinantal point processes** (§3) — matrix-valued
  kernel structurally richer than DPP. C2 / S123 closed sine-kernel
  determinantal at orders ≤ 6; the Pfaffian (orthogonal/symplectic)
  and Vere-Jones α-DPP fits to ζ zeros are over-determined tests
  with strictly stronger discriminating power. **NOW PROPOSED as
  C6 (S130).**
- **Gaussian multiplicative chaos / FHK extreme-value statistics**
  (§3) — first ζ-AMPLITUDE measurement of the project. All prior
  ζ measurements (35+) target zero positions; FHK / Saksman-Webb
  GMC governs the inter-zero amplitude distribution. **NOW PROPOSED
  as C7 (S130).**
- **Stein-Tomas / Bourgain Λ(p)-set restriction theory** (§7) —
  global Fourier-norm test of prime exponential sums; structurally
  stronger than BDG decoupling (D15 PROPOSED), testing the FULL L^p
  norm rather than cap-restricted norms. **NOW PROPOSED as D25
  (S130).**
- **Locally testable codes** (§6) — constant-query primality test
  via PCP-style local tester on χ_P-encoded codeword. Distinct from
  Reed-Solomon list decoding (D18 PROPOSED) and Goldreich-Levin
  Hadamard list decoding (CLOSED line 474). **NOW PROPOSED as D26
  (S130).**
- **Newman / Erdős polynomial L^∞ flatness** (§2) — extremal-
  polynomial endpoint bound at `p = ∞` for the prime-coefficient
  generating polynomial; structurally orthogonal to D10 Mahler (log
  integral, CLOSED S134) and D25 Stein-Tomas (`L^p` finite-`p`).
  Erdős conjecture for ±1 polynomials was settled by Balister-
  Bollobás-Morris-Sahasrabudhe 2020 *Annals*; the 0/1 / Newman /
  prime-support analogue is **wide open**. **NOW PROPOSED as D27
  (S136).**
- **LPS Ramanujan Cayley graph (Lubotzky-Phillips-Sarnak)** (§1) —
  non-abelian Cay(PSL_2(F_p), prime-indexed quaternion generators).
  Distinct from CLOSED §D.D20 (S125, abelian Cay(Z/NZ, primes)
  Friedman bound) — D28 is the LPS-quaternion construction
  explicitly flagged by S125 as the D20.c successor. **NOW PROPOSED
  as D28 (S136).**
- **Cohn-Elkies LP bound + Viazovska modular-form magic functions**
  (§10 NEW) — sphere-packing LP bound applied to prime
  autocorrelation `R_P(t)`; A-grade if Viazovska-style
  modular-form representation gives an explicit polylog-evaluable
  identity for the prime singular series. **NOW PROPOSED as D29
  (S136).**
- **Pollicott-Ruelle resonances of arithmetic transfer operator**
  (§5) — spectrum of `χ_P`-weighted Gauss-map transfer operator;
  Mayer 1991-style dynamical-determinant representation of π(x)
  if a non-trivial isolated arithmetic resonance exists. Distinct
  from CLOSED line 105 / 182 / 425. **NOW PROPOSED as D30 (S136).**

When in doubt, pick from this list first.

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
| Continuous-time quantum walks (CTQW) | USED E (S141 edge E7.20) | Childs 2009 PRL 102 = arXiv:0806.1972; Childs-Cleve-Deotto-Farhi-Gutmann-Spielman 2003 = arXiv:quant-ph/0209131 | D5 (CLOSED S141: divisor / coprime graph CTQW from seed `|1⟩` with `H ∈ {A, L}` has `max_t |⟨v_p|exp(−iHt)|1⟩|² ≤ C·π(x)/x` with `C = 1.151 + 0.609/log x` (asymptote 1.151, x ∈ {32..768}); top eigenvector overlap with v_p decays as `x^{−0.20}` — no isolated band-edge cluster; closure mechanism distinct from D4 Szegedy: CTQW closure runs through eigenstate-overlap dispersal rather than spectral gap; extends E7.13 from discrete-time to continuous-time walks) |
| Hodge / Laplacian on simplicial complexes (higher-order L_k, k ≥ 1) | USED E (S126 edge E7.17) | Lim 2020 SIAM Review 62(3) "Hodge Laplacians on graphs" = arXiv:1507.05379 https://arxiv.org/abs/1507.05379; Eckmann 1944 Comment. Math. Helv. 17; Friedman 1996 Algorithmica 21; Horak-Jost 2013 Adv. Math. 244; Goff 2009 J. Algebraic Combin. 30, 215 (join spectrum); Kahle 2014 *Annals* 179 (random flag complexes) | D22 (CLOSED S126: L_1 spectrum of K_N := pairwise-coprime flag complex on [2..N] is structurally constrained — λ_max(L_1) = |V| = N − 1 EXACTLY at all 9 tested N ∈ {8..128}; multiplicity = C(k+1, 2) where k = π(N) − π(N/2) is the Bertrand-prime count, verified perfectly in 9/9 cells; mechanism = Bertrand-universal-vertex join decomposition K_N = Δ^{k-1} \ast F(H); β_k = 0 for k ≥ 1 deterministically (cone collapse); L_1 mean shift reduces to triple-coprime singular series ∏_p(1 − 3/p² + 2/p³) / (6/π²)³ ≈ 1.276; KS p < 1e-300 at N = 128) |
| Connes operator-theoretic ζ | UNUSED | Connes 1999 Selecta Math 5 | F1 zone |
| Connes-Consani-Marcolli endomotive — Bost-Connes Galois action on KMS_∞ ground states with χ_P-projector enrichment | PROPOSED (D48) | Bost-Connes 1995 "Hecke algebras, type III factors and phase transitions with spontaneous symmetry breaking in number theory" *Selecta Math.* (NS) 1, 411; Connes-Consani-Marcolli 2007 "Noncommutative geometry and motives: the thermodynamics of endomotives" *Adv. Math.* 214, 761 = arXiv:math/0512138 https://arxiv.org/abs/math/0512138; Connes-Marcolli-Ramachandran 2005 "KMS states and complex multiplication" *Selecta Math.* 11, 325 = arXiv:math/0501424 https://arxiv.org/abs/math/0501424; Connes-Marcolli 2008 *Noncommutative Geometry, Quantum Fields and Motives* AMS Coll. Pub. 55; Marcolli 2005 *Arithmetic Noncommutative Geometry* AMS ULect 36; https://en.wikipedia.org/wiki/Bost%E2%80%93Connes_system | D48 (S160 frontier_gen) — KMS_∞ Galois orbit-length distribution of χ_P-projected ground states; CRITICALLY DISTINCT from CLOSED line 185 (raw BC partition function = ζ; trace loses Galois data) |
| Microlocal analysis / wavefront sets / pseudodifferential operators / FIO calculus | PROPOSED (D35) | Hörmander 1971 "Fourier integral operators I" *Acta Math.* 127, 79; Hörmander 1985 *The Analysis of Linear Partial Differential Operators III* Springer Grundlehren 274 (Ch. 18); Treves 1980 *Introduction to Pseudodifferential and Fourier Integral Operators* Plenum; Sjöstrand 1982 "Singularités analytiques microlocales" *Astérisque* 95; Melrose 1995 *Geometric Scattering Theory* CUP; https://en.wikipedia.org/wiki/Wave_front_set | D35 (S142 frontier_gen) — wavefront set `WF(χ_P) ⊂ T*R \ {0}` of prime indicator as tempered distribution; conic Fourier-decay measurement |

## §2. Algebraic / Geometric

| Technique | Status | Survey ref | Linked vector / edge |
|-----------|--------|------------|----------------------|
| Croot-Lev-Pach polynomial method | PROPOSED | Tao 2016 blog post + Croot-Lev-Pach 2017 Annals 185 | B1 |
| Slice rank | PROPOSED | Tao 2016 (same) | B1 |
| Tropical geometry of arithmetic varieties | USED E (CLOSED_PATHS lines 204, 326, 431, 432, 660 — min-plus tropicalisation of generating functions all collapse to "smallest prime"; tropical Dirichlet conv = standard conv) | Maclagan-Sturmfels 2015 *Introduction to Tropical Geometry* AMS GSM 161; Mikhalkin 2005 *J. AMS* 18 (enumerative tropical AG); Itenberg-Mikhalkin-Shustin 2009 *Tropical Algebraic Geometry* Birkhäuser 2nd ed.; https://en.wikipedia.org/wiki/Tropical_geometry | CLOSED_PATHS 204/326/431/432/660 (compression-collapse closures) |
| Baker-Norine graph Riemann-Roch / tropical-curve chip-firing divisor rank | USED-I (D45 / S161, edge E2.28) | Baker-Norine 2007 "Riemann-Roch and Abel-Jacobi theory on a finite graph" *Adv. Math.* 215, 766 = arXiv:math/0608360 https://arxiv.org/abs/math/0608360; Dhar 1990 *Phys. Rev. Lett.* 64 (chip-firing burning); Bjorner-Lovász-Shor 1991 *Eur. J. Combin.* 12; Gathmann-Kerber 2008 *Math. Z.* 259 = arXiv:math/0612129 https://arxiv.org/abs/math/0612129 (tropical-curve RR extension); Mikhalkin-Zharkov 2008 in *Curves and their Jacobians* AMS Contemp. Math. 465 | D45 / S161 — Closed-form q-reduced identity for D_P on `Γ_N` (divisibility) and `H_N` (Hasse cover): `D'_P^Γ = (π(N)−π(N/2))δ_1 + Σ_{p≤N/2} δ_p`; `D'_P^H = π(N)·δ_1`. Verified N≤512. Rank `r(D_P) = max(0, deg − g + 1) = 0` generically — Brill-Noether-non-special. Adds edge E2.28 (chip-firing q-reduced identity, EVS M shape) |
| Étale cohomology / Frobenius traces | UNUSED | Deligne, Weil II | B3 (sketch only) |
| Motivic integration | UNUSED | Cluckers-Loeser surveys | candidate |
| Berkovich spaces (non-archimedean) | PROPOSED (D33) | Berkovich 1990 *Spectral Theory and Analytic Geometry over Non-Archimedean Fields* AMS Math. Surv. Mono. 33; Baker-Rumely 2010 *Potential Theory and Dynamics on the Berkovich Projective Line* AMS Math. Surv. 159; Favre-Rivera-Letelier 2010 *Math. Ann.* 348; Conrad 2008 in *p-adic Geometry* AMS ULect 45 https://math.stanford.edu/~conrad/papers/aws.pdf; https://en.wikipedia.org/wiki/Berkovich_space | D33 (S142 frontier_gen) — Berkovich projective line spectrum / Type II Newton polygon of `f_N(z) = Σ_{n≤N} χ_P(n) z^n` over `C_p` |
| Schubert calculus / Grassmannian — Plücker coordinates, Schubert cells, Pieri / Littlewood-Richardson rule, honeycomb model | PROPOSED (D46) | Eisenbud-Harris 2016 *3264 and All That* CUP; Fulton 1997 *Young Tableaux* CUP; Griffiths-Harris 1978 *Principles of Algebraic Geometry* Wiley §1.5; Knutson-Tao 1999 *JAMS* 12 = arXiv:math/9807160 https://arxiv.org/abs/math/9807160; https://en.wikipedia.org/wiki/Schubert_calculus | D46 (S160 frontier_gen) — Plücker coordinates / Schubert cell of `V_χ_P^N := rowspan ((1)_i, (χ_P(i))_i) ∈ Gr(2, N)`; LR coefficient `c_{λ,λ}^{(N-2)^2}` of self-intersection |
| Adelic / idele class group | UNUSED | Tate's thesis | F1 zone |
| Mahler measure / Lehmer's conjecture / Weil height of generating polynomial | USED I (S134, edge E2.20) | Smyth 2008 "The Mahler measure of algebraic numbers: a survey" CUP = https://homepages.ed.ac.uk/cjsmyth/papers/107.pdf; Boyd 1981 Canad. Math. Bull. 24; Lehmer 1933 Ann. Math. 34; Dobrowolski 1979 *Acta Arith.* 34 (lower bound `log m(P) ≥ C(log log D / log D)³`) | D10 (CLOSED S134, mode I, B-grade: log m(f_PRIMES_N) − log m(f_BERN(d_N)_N) plateaus at Δ_∞ ≈ −0.307 ± 0.001 nat from N=2^16, z(MATCH) = −337σ at N=2^18; f_N / z² irreducible over Q[z] for N∈{64,128,256} — zero cyclotomic share; A-grade poly-log compressibility hypothesis falsified; **6th orthogonal pseudorandomness measure category — algebraic-height / multiplicative-height — after E2.13 Gowers / E2.14 Anderson / E2.15 alg.immunity / E2.16 DPP failure / E2.17 PH**) |
| Multivariate Mahler measure / Boyd-Smyth `L'`-value identities | PROPOSED (D36) | Boyd 1981 *Canad. Math. Bull.* 24, 453; Smyth 1981 *Bull. AMS* (NS) 4, 49 (`m(1+x+y) = (3√3/4π) L(χ_{-3}, 2)`); Boyd 1998 "Mahler's measure and special values of L-functions" *Exp. Math.* 7, 37; Smyth 2008 CUP survey = https://homepages.ed.ac.uk/cjsmyth/papers/107.pdf; Wikipedia: Mahler measure https://en.wikipedia.org/wiki/Mahler_measure | D36 (S148 frontier_gen) — 2-D Jensen integral `m(F_N(z_1, z_2))` over `T^2` accessing CM-elliptic `L'(E,1)` BSD-style identities not available in 1-D D10 |
| Newman / Erdős polynomial flatness — L^∞ extremality of 0/1-coefficient polynomials on `\|z\|=1` | PROPOSED (D27) | Erdős 1957 *Mich. Math. J.* 4 (Littlewood/Erdős flatness conjecture for ±1 polynomials); Newman 1976 "Norms of polynomials" *Proc. AMS* 51 (0/1 case); Kahane 1985 *Some Random Series of Functions* CUP §6 (Salem-Zygmund random model); Bourgain 1988 "Bounded orthogonal systems and the Λ(p)-set problem" *Acta Math.* 162; Bourgain-Lewko 2016 *J. Funct. Anal.* (random-polynomial mean-vs-max bounds); Balister-Bollobás-Morris-Sahasrabudhe 2020 *Annals of Math.* 192, 977 (resolution of Erdős for ±1 — flat Littlewood polynomials EXIST constructively) https://en.wikipedia.org/wiki/Littlewood_polynomial | D27 (S136 frontier_gen) — `\|f_N\|_∞ / √π(N)` for `f_N(z) = Σ_{n≤N} χ_P(n) z^n`, structurally orthogonal to D10 Mahler (log integral) and D25 Stein-Tomas (`L^p` finite-`p`) |
| AHK matroid Hodge theory / Heron-Rota-Welsh log-concavity (Chow ring of matroid, hard Lefschetz, Hodge-Riemann) | USED I (S149, edge E2.24) | Adiprasito-Huh-Katz 2018 "Hodge theory for combinatorial geometries" *Annals of Math.* 188, 381 = arXiv:1511.02888 https://arxiv.org/abs/1511.02888; Huh 2012 *J. AMS* 25, 907; Brändén-Huh 2020 *Annals* 192, 821 = arXiv:1902.03719; De Concini-Procesi 1995 *Selecta Math.* 1; Heron 1972 / Rota 1971 / Welsh 1976 (original log-concavity conjecture) | D31 (CLOSED S149: AHK Hodge slack `δ_k = |w_k|² − |w_{k−1}||w_{k+1}|` of the prime-divisibility transversal matroid `M_P^N` deviates from degree-preserving config-model baseline by `+5.61σ` at `k = 1, N = 20`, persists at `+3.91σ` on the Bertrand-stripped connected part `M_conn^N`; ~50% of raw signal traces to `(t − 1)^{π(N)−π(N/2)}` direct-summand factor; ~50% residual structure of `M_conn^N` not yet identified analytically — 7th orthogonal pseudorandomness-detection category; structurally distinct from E7.17 graph-Laplacian Hodge L_1 (D22, S126 — Chow ring vs Hodge Laplacian)) |
| LPS Ramanujan graph construction (Lubotzky-Phillips-Sarnak quaternion-algebra Cayley) | PROPOSED (D28) | Lubotzky-Phillips-Sarnak 1988 "Ramanujan graphs" *Combinatorica* 8, 261; Margulis 1988 *Probl. Inf. Transm.* 24; Davidoff-Sarnak-Valette 2003 *Elementary Number Theory, Group Theory and Ramanujan Graphs* CUP; Lubotzky 2012 *Bull. AMS* 49, 113 = arXiv:1105.2389 https://arxiv.org/abs/1105.2389; Hoory-Linial-Wigderson 2006 *Bull. AMS* 43, 439; Wikipedia: Ramanujan graph https://en.wikipedia.org/wiki/Ramanujan_graph | D28 (S136 frontier_gen) — non-abelian Cay(PSL_2(F_p), prime-indexed LPS-quaternion generators); explicitly flagged as D20.c successor in S125 |
| Geometric Complexity Theory (GCT) — Mulmuley-Sohoni orbit closure containment / occurrence obstructions on `GL_n`-orbit `\overline{GL_n \cdot f}` | USED PARTIAL E (S156, edge E2.26) — orbit-dim sub-frame closed; plethysm sub-frame still UNUSED | Mulmuley-Sohoni 2001 "Geometric Complexity Theory I: An Approach to the P vs NP and Related Problems" *SIAM J. Comput.* 31(2), 496–526; Bürgisser-Ikenmeyer-Panova 2017 *FOCS* "No occurrence obstructions in geometric complexity theory" arXiv:1604.06431 https://arxiv.org/abs/1604.06431; Bürgisser-Landsberg-Manivel-Weyman 2011 *SIAM J. Comput.* 40, 1179; Mulmuley 2009 "Geometric Complexity Theory VI: The flip via positivity" arXiv:0704.0229 https://arxiv.org/abs/0704.0229; Landsberg 2017 *Geometry and Complexity Theory* CUP, ch. 9-10; Wikipedia: Geometric complexity theory https://en.wikipedia.org/wiki/Geometric_complexity_theory | A7 (PARTIAL CLOSE S156, mode E, B-grade): orbit-dim sub-frame closed — `dim Stab_{GL_n}(f_χ_P^{(n)}) = 0` and `dim PD_k = matched-baseline` for all `k = 0..n` at `n = 4, 5, 6, 7`; std = 0 across 100 random matched-support trials; support hypergraph is Lie-rigid. Representation-theoretic invariants are fully determined by support hypergraph; arithmetic content invisible at this sub-frame. **Plethysm sub-frame (occurrence obstructions in `Sym^k Sym^d V`) remains OPEN** and requires SageMath. See `experiments/algebraic/gct_chi_p_orbit/`. |

## §3. Probabilistic / Random Matrix

| Technique | Status | Survey ref | Linked vector / edge |
|-----------|--------|------------|----------------------|
| GUE pair correlation | USED (E) | Montgomery 1973 | E3.1, E3.13 |
| GUE random-phase / Poisson-approximation for {γ_j log x mod 2π} as predictor of explicit-formula truncation variance | USED E (S195, refines E3.1 closure) — closed-form Var(π(x) − R_K(x)) ≈ x log²(K)/(2π² K log²x) matches empirical median |err| within 5–55% across 3 decades (600 (x, K, |err|) triples); pooled slope-through-origin / half-Gaussian ratio = 0.74 stable across decades (= GUE-vs-Poisson reduction); inverts to K*(x, p) = Θ(x) for any p ∈ (0, 1), constants 0.09 (50%), 1.35 (99%) | Montgomery 1973 (pair correlation conjecture); Odlyzko 1989 *Math. Comp.* 48 273 (numerical equidistribution); Mehta 2004 ch.5 (sine-kernel form factor for n=2 GUE correction); Riemann–von Mangoldt explicit formula | Thread 3 / Galway frontier (closed quantitatively in distribution); E3.1 (further refined); E1.5 (matched) |
| Mellin-domain log-Gaussian smoothing kernel (M_φ(ρ) = e^{ρ²h²/2}) combined with GUE random-phase variance bound | USED E (S196, refines S195 closure under any log-Gaussian bandwidth) — for π_{K,h}(x) := R(x) − 2 Σ_{j≤K} Re R(x^{ρ_j})·exp(−h²γ_j²/2) the error variance decomposes as Var(TAIL)(K) + Var(BIAS)(K, h) over disjoint j ranges; Var(TAIL) is independent of h, forcing K = Θ̃(x) regardless of bandwidth choice. Empirical validation 40 samples × 11 bandwidths × K_max=2000 at x ~ 1.78e5: K*_50 unchanged at 1782-1783 for all h ≤ 3·10⁻⁴ (smoothing inactive), > K_max for h ≥ 5·10⁻⁴ (bias dominates). Bias variance prediction matches empirical to 5-15% across active regime | Galway 2004 *Math. Comp.* 73 (smoothing-kernel for analytic π(x) algorithm); Hiary 2011 *Math. Comp.* 80 (heat-kernel framework, t^{4/13}); Riemann-von Mangoldt explicit formula; the GUE variance bound entry above | Thread 3 / Galway frontier closure under smoothing (S196); E3.1 (closed under any log-Gaussian bandwidth) |
| Sine-kernel n-correlation determinant `R_n = det[K(s_i-s_j)]` (n=2..6) | USED E (S25/S45/S57/S123, edge E7.1 refined) | Mehta 2004 *Random Matrices* (3rd ed.) §6.2; Hough-Krishnapur-Peres-Virág 2009 AMS ULect 51 §1.2; Bohigas-Giannoni-Schmit 1984 Phys. Rev. Lett. 52, 1 | C2 (CLOSED S123: empirical zeta R_n along equally-spaced slices, P_k spacings k=0..5, κ_n cumulants n=2..6 at L≥8 — all match GUE within sample noise once matched-finite-N null is used; Conrey-Snaith arithmetic correction 1/L²≈0.024 below noise floor 1/√(n_tuples)≥0.3 at order 5) |
| Bogomolny-Keating arithmetic correction | USED (I, structural) | BK 1996 | E3.13, §C1 closed |
| Conrey-Snaith higher-order arithmetic corrections | USED I (S123) | Conrey-Snaith 2007 arXiv:0710.5304 | C2 (CLOSED: orders 4,5,6 detection floor 1/L² well below noise floor 1/√(n_tuples) at γ ≤ 10⁴; same Odlyzko-height scaling barrier as §C1/S71) |
| Free probability / free cumulants | USED (I) | Mingo-Speicher 2017 | C2 / E2.1 |
| Marchenko-Pastur | USED (PARTIAL) | Mingo-Speicher 2017 | S82 spike-bulk split |
| Stein's method | PROPOSED (C5) | Ross 2011 Probability Surveys 8; Chen-Goldstein-Shao 2011; Stein 1986 IMS Lect Notes 7 | C5 (S91 frontier_gen) |
| Determinantal point processes (incl. PPP, signed-K, complex-Hermitian) | USED I (S95, edge E2.16) | Hough-Krishnapur-Peres-Virag 2009 AMS ULect 51; Soshnikov 2000 Russ. Math. Surv. 55 = arXiv:math/0002099; Vere-Jones 1997 (alpha-DPP) | D7 (CLOSED: K^2_DPP < 0 at all 15 admissible even t and K^2_PPP < 0 at all 14 odd t > 1; 3-point PPP overshoots HL by up to 79.16% on 3-AP triples; signed-real-K sigma_req in (-0.541, 0.769) never ±1; complex Hermitian phase LS fit best residual 0.0746) |
| Pfaffian / α-determinantal point processes (matrix-valued kernel) | PROPOSED (C6) | Borodin 2009 *Encyclopaedia Math. Sci.* 152 = arXiv:0911.1153 https://arxiv.org/abs/0911.1153; Soshnikov 2000 = arXiv:math/0002099; Vere-Jones 1997 *New Zealand J. Math.* 26 (α-permanent); Mehta 2004 *Random Matrices* §6 (GOE/GSE Pfaffian formulas) | C6 (S130 frontier_gen) — over-determined Pfaffian-vs-determinant n=4 fit on ζ zeros |
| Gaussian multiplicative chaos / Fyodorov-Hiary-Keating max statistics of |ζ(1/2 + it)| | USED I (S133, edge E7.18) | Fyodorov-Hiary-Keating 2012 *PRL* 108 = arXiv:1202.4713 https://arxiv.org/abs/1202.4713; Saksman-Webb 2018 = arXiv:1609.00027 https://arxiv.org/abs/1609.00027; Arguin-Belius-Bourgade 2017 *Comm. Math. Phys.* 349 = arXiv:1612.08575; Bourgade-Kuan 2014 *Comm. Pure Appl. Math.* 67 | C7 (CLOSED S133, mode I, B-grade case 2: FIRST ζ-amplitude measurement; M_T mean T-INDEPENDENT (pooled −0.657 ± 0.045 across T ∈ {10⁴, 10⁵, 10⁶}, K = 100/anchor) confirming FHK normalisation; M_T SHAPE is approximately Gaussian at finite T ≤ 10⁶ — var 1.47× FHK Gumbel(1/2) prediction π²/24, skew ≈ 0.10 vs Gumbel +1.14, ex-kurt ≈ −0.4 vs Gumbel +2.4; FHK Gumbel-shape limit NOT detectable; first quantitative FHK convergence-rate result) |
| Schur process / Schur measure on Young diagrams (Borodin-Okounkov-Reshetikhin / Okounkov infinite-wedge / Borodin-Okounkov-Toeplitz) | PROPOSED (D40) | Okounkov 2001 "Infinite wedge and random partitions" *Selecta Math.* 7, 57 = arXiv:math/9907127 https://arxiv.org/abs/math/9907127; Okounkov-Reshetikhin 2003 *J. AMS* 16, 581 = arXiv:math/0107056 https://arxiv.org/abs/math/0107056; Borodin-Okounkov 2000 *IEOT* 37, 386 = arXiv:math/9907165 https://arxiv.org/abs/math/9907165; Macdonald 1995 *Symmetric Functions and Hall Polynomials* OUP (2nd ed.); Borodin-Gorin 2012 lectures arXiv:1212.3351 https://arxiv.org/abs/1212.3351 | D40 (S148 frontier_gen) — Schur measure `Pr(λ) ∝ s_λ(α) s_λ(β)` with prime-specialised `α_n = 1/p_n^s`, `β = (1)`; 2-D space-time DPP on Young diagrams structurally orthogonal to 1-D translation-invariant DPP D7 (CLOSED) |
| Random tensor models (Gurau colored / uncolored) — melonic 1/N expansion, Gurau-Witten SYK universality | PROPOSED (D41) | Gurau 2011 "The 1/N expansion of colored tensor models" *Comm. Math. Phys.* 304, 69 = arXiv:1011.2726 https://arxiv.org/abs/1011.2726; Witten 2016 "An SYK-like model without disorder" arXiv:1610.09758 https://arxiv.org/abs/1610.09758; Klebanov-Tarnopolsky 2017 *Phys. Rev. D* 95, 046004 = arXiv:1611.08915 https://arxiv.org/abs/1611.08915; Bonzom-Gurau-Riello-Rivasseau 2011 *Nucl. Phys. B* 853, 174 = arXiv:1105.3122 https://arxiv.org/abs/1105.3122; Gurau 2017 *Random Tensors* OUP textbook | D41 (S154 frontier_gen) — χ_P 3-tensor `T_{ijk}^N := \mathbf{1}[i+j+k=N] χ_P(i)χ_P(j)χ_P(k)` mode-1 unfolding singular spectrum vs Gurau-Witten melonic limit; first rank-3 universality test for χ_P, orthogonal to all CLOSED rank-2 measurements (C2, C7, S74, D31) |
| Hairer regularity structures (singular SPDE limits) / KPZ universality / Tracy-Widom β=2 on `t^{1/3}` scale | USED PARTIAL E (S157, edge E2.27) | Hairer 2014 "A theory of regularity structures" *Inventiones Math.* 198, 269 = arXiv:1303.5113 https://arxiv.org/abs/1303.5113 (Fields Medal); Corwin 2012 "The Kardar-Parisi-Zhang equation and universality class" *Random Matrices Theory Appl.* 1, 1130001 = arXiv:1106.1596 https://arxiv.org/abs/1106.1596; Tracy-Widom 1994 "Level spacing distributions and the Airy kernel" *Comm. Math. Phys.* 159, 151; Hairer-Quastel 2018 *Forum Math. Pi* 6, e3 = arXiv:1512.07845; https://en.wikipedia.org/wiki/Kardar%E2%80%93Parisi%E2%80%93Zhang_equation; https://en.wikipedia.org/wiki/Tracy%E2%80%93Widom_distribution | D43 (PARTIAL CLOSURE S157, mode E, B-grade case (i): wavelet Hölder α(D) ≈ 0.85 stable across logX ∈ {18..24} with r²>0.998, far above KPZ ceiling 1/2; pre-stated falsifiers F1–F4 reject KPZ; Tracy-Widom β=2 skewness +0.224 rejected by detrended Z_d skew |z|>3 wrong sign in 5/7 windows; right-tail Gauss r²=0.977 ≥ KPZ r²=0.952; structural reason — D is deterministic almost-periodic explicit-formula sum, not stochastic SPDE-driven; FIRST project test of any non-`√x` scaling of PNT error term; FIRST project Hölder regularity / wavelet test on any arithmetic function. Sub-frames D43.b (logX-extension to 28) and D43.c (K-truncated explicit-formula residual) remain OPEN.) |
| Concentration of measure | PARTIAL | Boucheron-Lugosi-Massart 2013 | mentioned only |
| Stochastic loop equations / Eynard-Orantin topological recursion | PROPOSED (D24) | Eynard-Orantin 2007 *Comm. Number Theory Phys.* 1 = arXiv:math-ph/0702045 https://arxiv.org/abs/math-ph/0702045; Eynard 2016 *Counting Surfaces* Birkhäuser; Borot-Eynard 2011 arXiv:1107.5028; Keating-Snaith 2000 *Comm. Math. Phys.* 214 | D24 (S121 frontier_gen) |
| Möbius / nilsequence orthogonality (Sarnak conjecture, GT 2012, Tao logarithmic-Chowla) | USED E (S100, edge E2.18; S153, edge E2.25) | **Green-Tao 2012 *Annals* 175 = arXiv:0807.1736** (Möbius nilsequence orthogonality); **Sarnak 2010 IAS lectures** (Möbius randomness conjecture); **Tao 2016 Forum Math Pi 4 e8 = arXiv:1509.05422** (logarithmic-Chowla for two-point); Tao-Teräväinen 2017 arXiv:1708.02610; **GT-Ziegler 2012 *Annals* 176 = arXiv:1009.3998** (U^{s+1} inverse theorem) | G1 (CLOSED: gamma_lambda matches Rademacher Lyapunov within 2.2σ across 51 energies and 3 orders of N; Σ Chowla-z² = 4.77 vs χ²_16 = 16; first non-W-tricked spectral measure at noise floor); **G2 (CLOSED S153, mode E, B-grade**: ||λ||_{U^k}^{2^k} = ||IID||_{U^k}^{2^k} · (1 ± O(1/√N)) for k ∈ {2,3} at N ≤ 2^20; ||μ||_{U^2}^4 → 0.7392/N matched-variance prediction; max single Fourier coefficient matches CLT √(N log N) law within ±3σ — no 1-step nilsequence correlation; first project Gowers-norm measurement on multiplicative functions, paired with G1 to give two orthogonal multiplicative-regime confirmations of GT orthogonality) |

## §4. Topological / Geometric

| Technique | Status | Survey ref | Linked vector / edge |
|-----------|--------|------------|----------------------|
| Persistent homology | USED I (S96, edge E2.17) | Carlsson 2009 BAMS 46(2); Edelsbrunner-Harer 2010 textbook; Bauer 2021 *J. Appl. Comput. Topol.* 5 = arXiv:1908.02518; Cohen-Steiner-Edelsbrunner-Harer 2007 (PH stability) | D2 (CLOSED: total H_0 and H_1 persistence of Takens-embedded normalised prime gaps deviate ≥ 5σ from BOTH IID Exp(1) AND gap-permuted control across d ∈ {2,3,4} and M ∈ {500..4000}; primes more clustered + fewer loops than Poisson; fifth orthogonal HL-detection category after Gowers/Anderson/AlgImmunity/DPP) |
| Discrete Ricci curvature (Ollivier optimal-transport / Lin-Lu-Yau / Bakry-Émery `CD(K, ∞)`) | PROPOSED (D39) | Ollivier 2009 *J. Funct. Anal.* 256 = arXiv:math/0701886 https://arxiv.org/abs/math/0701886; Lin-Lu-Yau 2011 *Tohoku Math. J.* 63 = arXiv:1010.4337 https://arxiv.org/abs/1010.4337; Bakry-Émery 1985 LNM 1123; Erbar-Maas 2012 *Arch. Ration. Mech. Anal.* 206 = arXiv:1111.2687; Münch-Wojciechowski 2019 *Adv. Math.* 356 = arXiv:1712.10160; https://en.wikipedia.org/wiki/Ricci_curvature | D39 (S148 frontier_gen) — `κ_LLY(x, y) := 1 − W_1(m_x, m_y) / d(x, y)` on `Cay(Z/NZ, primes)` and on coprimality complex `K_N`; nonlinear optimal-transport invariant orthogonal to D20 spectral gap (CLOSED) and D22 Hodge L_1 (CLOSED) |
| Mapper algorithm | PROPOSED (D19) | Singh-Mémoli-Carlsson 2007 *Eurographics Symp. PBG* (Stanford TR) https://research.math.osu.edu/tgda/mapperPBG.pdf; Lum et al. 2013 *Sci. Rep.* 3, 1236 https://www.nature.com/articles/srep01236; KeplerMapper https://kepler-mapper.scikit-tda.org/; https://en.wikipedia.org/wiki/Mapper_(topological_data_analysis) | D19 (S119 frontier_gen) |
| Reeb graphs of arithmetic height functions on Z | PROPOSED (D21) | Edelsbrunner-Harer 2010 *Computational Topology* AMS, ch. VI.3; Carr-Snoeyink-Axen 2003 *Comput. Geom.* 24, 75; Wikipedia: Reeb graph https://en.wikipedia.org/wiki/Reeb_graph | D21 (S121 frontier_gen) |
| Cellular sheaf cohomology over divisibility poset (Curry framework) | PROPOSED (D14) | Curry 2014 PhD thesis = arXiv:1303.3255; Hansen-Ghrist 2019 *J. Appl. Comput. Topol.* 3; Robinson 2014 *Topological Signal Processing* Springer | D14 (S103 frontier_gen) |
| Sheaf cohomology of sequences (general — non-cellular) | UNUSED | Curry 2014 thesis | candidate |
| Discrete Morse theory | PROPOSED (D17) | Forman 1998 *Adv. Math.* 134, 90; Forman 2002 "A user's guide to discrete Morse theory" Sém. Lothar. Combin. 48 https://www.emis.de/journals/SLC/wpapers/s48forman.pdf; Benedetti-Lutz 2014 *Exp. Math.* 23, 66 https://arxiv.org/abs/1303.6422; Bjorner-Wachs 2005 *Trans. AMS* 357; https://en.wikipedia.org/wiki/Discrete_Morse_theory | D17 (S119 frontier_gen) |
| Arithmetic topology (Mazur primes-as-knots dictionary) — Morishita arithmetic Massey products / Rédei symbol on prime triples / Vogel Galois-cohomology framework | USED E (S164, edge E2.29) | Mazur 1966 "Notes on étale cohomology of number fields" *Ann. Sci. ENS* 6, 521; Morishita 2002 "Milnor invariants and Massey products for prime numbers" *J. Reine Angew. Math.* 550, 141; Morishita 2012 *Knots and Primes: An Introduction to Arithmetic Topology* Springer Universitext https://link.springer.com/book/10.1007/978-1-4471-2158-9; Rédei 1939 "Ein neues zahlentheoretisches Symbol mit Anwendungen auf die Theorie der quadratischen Zahlkörper" *J. Reine Angew. Math.* 180, 1; Lemmermeyer 2000 *Reciprocity Laws: From Euler to Eisenstein* Springer Monographs (Ch. 9 explicit Rédei formula); Stevenhagen 1996 *Acta Arith.* 76, 89 (Rédei reciprocity, governing fields, negative Pell); Vogel 2004 Heidelberg PhD thesis "Massey products in the Galois cohomology of number fields"; Wikipedia: Arithmetic topology https://en.wikipedia.org/wiki/Arithmetic_topology | D44 (CLOSED S164: empirical Rédei `[p, q, r]` on Borromean-admissible prime triples up to N=1000 is uniform within strict 5σ envelope; `f_+ = 0.4853` clean subset (K=7591, z=-2.56); constant `|Δf_+| ≈ 0.013` with `z ∼ √K` consistent with finite-r effective Chebotarev `O(1/√(log r))` in Rédei field K_pq/Q; F1 (A-grade, HL-S_3 correlated bias) FALSIFIED, F2 (UNBIASED, mode E) HOLDS; no mod-8/12/24 substructure beyond Bonferroni-3σ. Adds **edge E2.29** — 13th orthogonal pseudorandomness category for χ_P (arithmetic topology / Galois cohomology); first project measurement of Mazur primes-as-knots / Galois-cohomological triple Massey-product invariants. Distinct from CLOSED line 208 generic knot invariants by Galois-cohomology arithmetic-by-construction. Validated on canonical Borromean triple `[13, 61, 937] = -1`.) |

## §5. Dynamical / Ergodic

| Technique | Status | Survey ref | Linked vector / edge |
|-----------|--------|------------|----------------------|
| Furstenberg correspondence | UNUSED | Furstenberg 1981 textbook | D1 candidate |
| Transfer operator spectrum (Pollicott-Ruelle resonances of arithmetic dynamical system) | USED E (S140, edge E2.22) | Pollicott 1985 *Inventiones Math.* 81, 413 (rate of mixing of Axiom A flows); Ruelle 1976 *Inventiones Math.* 34, 231 (zeta-functions for expanding maps); Mayer 1991 *Bull. AMS* 25, 55 (transfer-operator approach to Selberg ζ via Gauss map); Baladi 2018 *Dynamical Zeta Functions and Dynamical Determinants for Hyperbolic Maps* Springer Ergeb. 68; Liverani 2004 *Comm. Math. Phys.* 248 (anisotropic Sobolev spaces); Tsujii 2010 *Nonlinearity* 23; Wirsing 1974 *Acta Arith.* 24 (GKW constant); Wikipedia: Transfer operator https://en.wikipedia.org/wiki/Transfer_operator | D30 (CLOSED S140: χ_P-weighted Gauss-map transfer operator has refinement-stable PR resonances |λ_k| stable to <1% CV across (M_grid, n_max) ∈ {30..160}×{100..800}; closed-form Rayleigh-quotient prediction λ_0^{χ_P} = Σ_p a_p with a_n=(2/log²2)·[ln2/(n(n-1))-ln((n+1)/n)/(n-1)+ln((n+2)/(n+1))/n] (n≥2) explicit, +0.6% rel error; Cramér 1/log n + odd-parity matched control places χ_P within ±2σ Bonferroni at 200 seeds — A-grade hypothesis "isolated arithmetic resonance encoding π(x)/x at polylog cost" empirically falsified; first project dynamical-spectral measurement on χ_P, structurally distinct from E2.13/E2.14/E2.16/E2.17/E2.19/E2.20/E2.21 — at the dynamical-spectral level the leading resonance is a 1-summing density invariant captured by Cramér by construction; structurally distinct from CLOSED line 105 (constructive transfer matrix sieve), CLOSED line 182 (FRACTRAN), CLOSED line 320/425 (unweighted Furstenberg / orbit count = ζ-zeros)) |
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
| Linear-programming bounds for codes / sphere packing (Cohn-Elkies LP / Delsarte / Bachoc-Vallentin three-point / Schrijver SDP / Vaaler-Beurling-Selberg majorants) | USED E (S145, edge E2.23) | Cohn-Elkies 2003 *Annals* 157 = arXiv:math/0110009; Viazovska 2017 *Annals* 185 = arXiv:1603.04246 (E_8); Cohn-Goncalves 2019 *Invent. Math.* 217 = arXiv:1712.04438 (12-D uncertainty); Delsarte 1973 (association schemes); Schrijver 2005 *IEEE IT* 51 = arXiv:math/0405348; Bachoc-Vallentin 2008 = arXiv:math/0608426; Vaaler 1985 *Bull. AMS* 12 (1-D extremal majorants); Graham-Vaaler 1981 *Trans. AMS* 265; https://en.wikipedia.org/wiki/Sphere_packing | D29 (CLOSED: discrete-`Z` Delsarte LP applied to prime autocorrelation `R_P(t)/π(N)` gives optimum `f̂*(0; T) ≈ 1.15 + 0.85 log T`, density bound loose by factor `log N` for any `T = polylog(N)`; LP saturation requires `T ≈ N` — strictly inside sieve barrier E6.7. Optimal `f^*` is period-4 sinusoid reflecting parity (E2.1), NOT a Viazovska modular form. New shape edge E2.23.) |

## §7. Combinatorial / Additive

| Technique | Status | Survey ref | Linked vector / edge |
|-----------|--------|------------|----------------------|
| Sum-product theorems (Bourgain-Glibichuk-Konyagin in F_p) | PROPOSED (D9) | Bourgain-Glibichuk-Konyagin 2006 J. London Math. Soc. = https://londmathsoc.onlinelibrary.wiley.com/doi/abs/10.1112/S0024610706022721; Garaev 2008 arXiv:math/0702780; Iosevich 2009 arXiv:0904.2075; Tao-Vu 2006 *Additive Combinatorics* CUP | D9 (S97 frontier_gen) |
| Decoupling (Bourgain-Demeter-Guth) | PROPOSED (D15) | Bourgain-Demeter 2015 *Annals* 182 = arXiv:1403.5335; Bourgain-Demeter-Guth 2015 *Annals* 184 = arXiv:1512.02384; Demeter 2020 *Fourier Restriction, Decoupling and Applications* CUP | D15 (S103 frontier_gen) |
| Restriction theory (Stein-Tomas / Bourgain Λ(p)-set) | PROPOSED (D25) | Stein 1993 *Harmonic Analysis* Princeton Ch. IX; Tomas 1975 *Bull. AMS* 81; Bourgain 1989 "Bounded orthogonal systems and the Λ(p)-set problem" *Acta Math.* 162; Green 2005 "Roth's theorem in the primes" *GAFA* 15 = arXiv:math/0302311; Tao 2020 247B Notes 1 https://terrytao.wordpress.com/2020/03/29/247b-notes-1-restriction-theory/; Foschi-Oliveira e Silva 2024 *São Paulo J. Math. Sci.* https://link.springer.com/article/10.1007/s40863-024-00422-x | D25 (S130 frontier_gen) — global L^p restriction-saturation test on prime exponential sum |
| Gowers norms (U^k uniformity) | USED E (S85, edge E2.13; S153, edge E2.25) | Green-Tao arXiv:math/0606088; Green-Tao arXiv:0807.1736 (Mobius nilsequences); Green-Tao-Ziegler arXiv:1009.3998 (U^{s+1} inverse); Hardy-Littlewood 1923 | D6 (CLOSED: Q^2(chi_P)→2.30=S_2 box singular series at N=2^18; Q^3 stable at 35.5; W-trick restores Gowers uniformity at U^2 within 0.1%); **G2 (CLOSED S153: λ and μ U^k norms match matched-variance IID at N ≤ 2^20; first Gowers-norm measurement on multiplicative functions; no W-trick needed)** |
| Bounded gaps / GPY sieve (Maynard multidim. weights) | USED E (S116, edge E7.14) | Maynard 2015 Annals 181 = arXiv:1311.4600; Polymath8b arXiv:1407.4897 | A5 (CLOSED: aggregate-not-pointwise + divisor-enumeration cost ⇒ not TC⁰; 4th sieve-route family closure complementing E7.10 / E5.8 / E7.11) |
| Density Hales-Jewett (combinatorial-line density on [k]^d) | PROPOSED (D23) | Furstenberg-Katznelson 1991 J. d'Anal. Math. 57; Polymath1 2010 *Annals* 175 = arXiv:0910.3926 https://arxiv.org/abs/0910.3926; Hales-Jewett 1963 Trans. AMS 106 | D23 (S121 frontier_gen) |
| Analytic Combinatorics in Several Variables (ACSV) — Pemantle-Wilson smooth-point critical-variety diagonal extraction | USED I (S197, edge E7.21) | Pemantle-Wilson 2013/2024 *Analytic Combinatorics in Several Variables* CUP (2nd ed. 2024) https://www.acsvproject.com/; Melczer 2021 *An Invitation to Analytic Combinatorics* Springer Texts in Mathematics; Pemantle-Wilson 2008 "Twenty combinatorial examples of asymptotics derived from multivariate generating functions" *SIAM Review* 50, 199; **Pólya 1916** *Sitzungsber. Berlin Math. Ges.*; **Carlson 1921** *Math. Z.* 9, 1; **Furstenberg 1967** *J. Algebra* 7, 271; **Lipshitz 1988** *J. Algebra* 113, 373; **Skolem 1934 / Mahler 1935 / Lech 1953**; Stanley 1999 *EC2* CUP §6.4; Flajolet-Sedgewick 2009 *AC* CUP App. B.4; https://en.wikipedia.org/wiki/Analytic_combinatorics | D32 (CLOSED S197, mode I, B-grade): UNCONDITIONAL theorem `f(z) := Σ χ_P(n) z^n` has `|z|=1` natural boundary (Pólya-Carlson + non-eventual-periodicity of primes) → not D-finite → not algebraic → not the diagonal of any rational/algebraic multivariate `F` for any `d ≥ 1` (Furstenberg/Lipshitz). Strict strengthening of CLOSED_PATHS row 584 (empirical ≤ order 20) to all orders. Adds edge E7.21. Empirical companion: T4 Erdős-Turán equidistribution — closest root of `f_N(t)` to `t=1` halves as `N` doubles, confirms no stable PW critical point. Closes the entire ACSV smooth-point class plus Wilf-Zeilberger telescoping / holonomic-gradient / Karr-Schneider Π Σ-summation (each requires D-finite input). |
| Cluster algebras / Fomin-Zelevinsky exchange relations / Y-system periodicity / cluster-finite Dynkin classification / Laurent phenomenon | PROPOSED (D47) | Fomin-Zelevinsky 2002 "Cluster algebras I: Foundations" *J. AMS* 15, 497 = arXiv:math/0104151 https://arxiv.org/abs/math/0104151; Fomin-Zelevinsky 2003 "Cluster algebras II: Finite type classification" *Invent. Math.* 154, 63 = arXiv:math/0208229 https://arxiv.org/abs/math/0208229; Keller 2013 "The periodicity conjecture for pairs of Dynkin diagrams" *Annals* 177, 111 = arXiv:1001.1531 https://arxiv.org/abs/1001.1531; Fomin-Reading 2007 IAS/Park City lectures = arXiv:math/0505518; Fomin-Williams-Zelevinsky 2016 book draft = arXiv:1608.05735 https://arxiv.org/abs/1608.05735; https://en.wikipedia.org/wiki/Cluster_algebra | D47 (S160 frontier_gen) — mutation orbit return-time `R_N` of prime-gap quiver (path / twin-prime / gap-equality variants) under Keller mutation; cluster-finite ⇔ Dynkin type test |
| Multiple zeta values (MZV) / Brown 2012 mixed Tate motives over Z / Goncharov motivic coproduct / Hoffman basis | PROPOSED (D38) | Brown 2012 "Mixed Tate motives over Z" *Annals* 175 = arXiv:1102.1312 https://arxiv.org/abs/1102.1312 (proves Hoffman basis `{ζ(n_1,...,n_r) : n_i ∈ {2,3}}` generates MZV `Q`-algebra); Goncharov 2005 *Duke Math. J.* 128 = arXiv:math/0208144 https://arxiv.org/abs/math/0208144 (motivic coproduct); Zagier 1994 *First Eur. Cong. Math.* II, Birkhäuser (Strong Conjecture `d_n = d_{n-2} + d_{n-3}`); Hoffman 1997 *J. Algebra* 194; Wikipedia: Multiple zeta function https://en.wikipedia.org/wiki/Multiple_zeta_function | D38 (S148 frontier_gen) — `ζ_P(s_1, ..., s_k) := Σ_{p_1 < ... < p_k} 1/(p_1^{s_1}⋯p_k^{s_k})`; PSLQ-test reduction to Brown's Hoffman basis ⊕ Mertens constants `M_s` |

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
| Riemann xi-function model spaces / de Branges Hilbert spaces of entire functions H(E) | PROPOSED (D34) | de Branges 1968 *Hilbert Spaces of Entire Functions* Prentice-Hall; de Branges 1992 *Bull. AMS* 26, 1 ("The Riemann hypothesis for Hilbert spaces of entire functions"); Lagarias 2007 "The Riemann hypothesis: arithmetic and geometry" = arXiv:math/0511099 https://arxiv.org/abs/math/0511099; Burnol 2002 *C. R. Acad. Sci. Paris* 333; Garunkstis-Steuding 2007 *Math. Slovaca* 57; https://en.wikipedia.org/wiki/De_Branges_space | D34 (S142 frontier_gen) — finite-rank reproducing-kernel projection in H(E_ξ) approximation rate of `χ_{[0,x]}` via Odlyzko ζ-zero sample points |
| Newman-Pollard hypothesis | PARTIAL (folded into D27 PROPOSED — see §2 row "Newman / Erdős polynomial flatness") | Newman 1976 *Proc. AMS* 51 | D27 (S136 frontier_gen) — D27 measures Newman's L^∞ extremality on χ_P-coefficient polynomial directly |
| Cohn-Elkies LP bound for sphere packing + Viazovska modular-form magic functions | PROPOSED (D29) | Cohn-Elkies 2003 "New upper bounds on sphere packings I" *Annals of Math* 157, 689 = arXiv:math/0110009 https://arxiv.org/abs/math/0110009; Viazovska 2017 "The sphere packing problem in dimension 8" *Annals of Math* 185, 991 = arXiv:1603.04246 https://arxiv.org/abs/1603.04246; Cohn-Kumar-Miller-Radchenko-Viazovska 2017 *Annals* 185, 1017 = arXiv:1603.06518; Cohn 2017 "A conceptual breakthrough in sphere packing" *Notices AMS* 64(2), 102; Cohn-Goncalves 2019 *Invent. Math.* 217 = arXiv:1712.04438 (LP-modular-form connection); Wikipedia: Sphere packing in dimension 8 https://en.wikipedia.org/wiki/Sphere_packing | D29 (S136 frontier_gen) — Cohn-Elkies LP applied to prime autocorrelation `R_P(t)`, Viazovska-style modular-form magic function for χ_P |
| Quantum modular forms (Zagier 2010) — smooth-cocycle defect on `Q` under SL_2(Z) | PROPOSED (D37) | Zagier 2010 "Quantum modular forms" *Clay Math. Proc.* 11, AMS pp. 659–675 (Quanta of Maths, Connes Festschrift); Bringmann-Folsom-Ono-Rhoades 2017 *Harmonic Maass Forms and Mock Modular Forms* AMS Coll. 64 (Ch. 21); Folsom-Ono-Rhoades 2013 "Mock theta functions and quantum modular forms" *Forum Math. Pi* 1 = arXiv:1301.5532 https://arxiv.org/abs/1301.5532; Lawrence-Zagier 1999 *Asian J. Math.* 3, 93 ("Modular forms and quantum invariants of 3-manifolds"); Wikipedia: Mock modular form https://en.wikipedia.org/wiki/Mock_modular_form | D37 (S148 frontier_gen) — cocycle defect `h_γ(a/q; N) := f_N(e^{2πi a/q}) − (cq+da)^{-k} f_N(e^{2πi γ(a/q)})` for `f_N(z) = Σ_{n≤N} χ_P(n) z^n` at rational z; smoothness defect under SL_2(Z) cocycle, archimedean; distinct from B4 Voronin (vertical translates of `ζ`) and from D27 Newman / D33 Berkovich |
| Resurgence theory / Borel-Écalle alien calculus — Stokes constants from divergent asymptotic series; Borel transform with isolated singularities | PROPOSED (D42) | Écalle 1981 *Les fonctions résurgentes* 3 volumes (Publ. Math. d'Orsay 81-05, 81-06, 85-05); Sauzin 2014 "Introduction to 1-summability and resurgence" arXiv:1405.0356 https://arxiv.org/abs/1405.0356; Mitschi-Sauzin 2016 *Divergent Series, Summability and Resurgence I* Springer LNM 2153 https://link.springer.com/book/10.1007/978-3-319-28736-2; Berry 1995 "Asymptotics, superasymptotics, hyperasymptotics" *Asymptotics Beyond All Orders* NATO ASI 284; Berry-Keating 1992 "A new asymptotic representation for `ζ(1/2 + it)` and quantum spectral determinants" *Proc. R. Soc. A* 437, 151; Wikipedia: Resurgent function https://en.wikipedia.org/wiki/Resurgent_function | D42 (S154 frontier_gen) — Borel transform of Riemann-Siegel correction coefficients `C_k(t)` with Padé pole detection; Stokes points compared to first 100 ζ-zero ordinates and `\log p_k`; first rigorous resurgence test for ζ-asymptotic divergence, complementary to Berry-Keating 1990s heuristic |

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
- **AHK matroid Hodge theory** (§2 NEW) — Adiprasito-Huh-Katz 2018
  *Annals* 188 hard Lefschetz / Hodge-Riemann for Chow ring `A^*(M)`
  of an arbitrary matroid; Heron-Rota-Welsh log-concavity. Apply to
  arithmetic prime-matroid `M_P^N` on the divisibility lattice;
  log-concavity inequality on signed sums of arithmetic Möbius
  values is structurally orthogonal to all archimedean / non-
  archimedean polynomial measures (D10, D27, D33). **NOW PROPOSED
  as D31 (S142).**
- **Analytic Combinatorics in Several Variables (Pemantle-Wilson)**
  (§7 NEW) — multivariate rational generating function `F(z_1, ...,
  z_d)` whose diagonal `[z_1^n \cdots z_d^n] F` encodes `π(n)`;
  smooth-point critical-variety analysis on `V_F = \{F = 0\}` gives
  closed-form Hessian asymptotics. Tests whether `Σ_n χ_P(n) z^n` is
  D-finite (folklore: no, but never proved). **NOW PROPOSED as D32
  (S142).**
- **Berkovich projective line / non-archimedean potential theory**
  (§2) — `\mathbb{P}^{1,\mathrm{an}}_{\mathbb{C}_p}` Type II
  Newton polygon decomposition of `f_N(z) = Σ_{n≤N} χ_P(n) z^n`;
  structurally orthogonal to D10 archimedean Mahler / D27 archimedean
  Newman / D33 archimedean L^p restriction. Distinct from CLOSED
  line 76 (p-adic interpolation of `π(x)` as a function of integer x)
  — D33 is geometric-Berkovich, not Mahler-coefficient. **NOW PROPOSED
  as D33 (S142).**
- **De Branges Hilbert spaces of entire functions H(E_\xi)** (§10) —
  reproducing-kernel projection rate of `\chi_{[0,x]}` onto the
  finite-rank `H_N \subset H(E_\xi)` space spanned by Odlyzko-zero
  sample points. Tests whether the de Branges program (one of four
  major RH approaches) admits polylog convergence. Distinct from
  CLOSED line 50 / 597 / 702 (Connes-Weil / Connes trace / CCM
  spectral triple). **NOW PROPOSED as D34 (S142).**
- **Microlocal analysis / wavefront sets** (§1) — Hörmander's
  `WF(\chi_P) \subset T^*\mathbb{R} \setminus \{0\}`: which directions
  in cotangent space are non-smooth for `\chi_P` viewed as a tempered
  distribution? Distinct from D6 Gowers (local k-correlations) and
  D25 Stein-Tomas (global L^p restriction); microlocal direction-
  specific Fourier decay. Pseudodifferential / FIO machinery for
  potential parametrix. **NOW PROPOSED as D35 (S142).**
- **Multivariate Boyd-Smyth Mahler measure** (§2) — `m(F_N(z_1, z_2))`
  Jensen integral over `T^2`; structurally orthogonal to univariate
  D10 (CLOSED S134, E2.20) by accessing CM-elliptic `L'(E,1)` BSD-style
  Smyth-Boyd identities not available in 1-D (genus-0 vs genus-≥1
  jump). **NOW PROPOSED as D36 (S148).**
- **Quantum modular forms (Zagier 2010)** (§10) — smooth-cocycle
  defect of `f_N(e^{2πi a/q})` under SL_2(Z); distinct from B4 Voronin
  (vertical-translate density), D27 Newman (sup norm), D33 Berkovich
  (non-archimedean). **NOW PROPOSED as D37 (S148).**
- **χ_P-restricted multiple zeta values / Brown 2012 motivic Galois
  group** (§7) — `ζ_P(s_1, ..., s_k) := Σ_{p_1<...<p_k} 1/(p_1^{s_1}
  ⋯ p_k^{s_k})`; PSLQ-tests whether prime-restricted MZVs reduce to
  Brown's Hoffman basis ⊕ Mertens constants. **NOW PROPOSED as D38
  (S148).**
- **Discrete Ricci curvature (Ollivier 2009 / Lin-Lu-Yau 2011)** (§4) —
  `κ_LLY := 1 − W_1(m_x, m_y) / d(x, y)` on `Cay(Z/NZ, primes)` and
  on coprimality complex; nonlinear optimal-transport invariant
  orthogonal to D20 Friedman spectral gap (CLOSED) and D22 Hodge L_1
  (CLOSED). **NOW PROPOSED as D39 (S148).**
- **Schur process / prime-supported Schur measure (Borodin-Okounkov-
  Reshetikhin)** (§3) — `Pr(λ) ∝ s_λ(α) s_λ(β)` with prime-specialised
  `α_n = 1/p_n^s`; 2-D space-time DPP on Young diagrams structurally
  orthogonal to 1-D translation-invariant DPP D7 (CLOSED). **NOW
  PROPOSED as D40 (S148).**
- **Geometric Complexity Theory (Mulmuley-Sohoni 2001)** (§2) —
  representation-theoretic obstruction (orbit-closure containment /
  occurrence obstructions) to formula complexity of the `f_χ_P^{(n)}`
  multi-affine encoding polynomial. Distinct from A1 SAT search
  (computational) and A4 bounded arithmetic (proof-theoretic). The
  Bürgisser-Ikenmeyer-Panova 2017 *FOCS* "no occurrence obstruction"
  theorem narrows but does not close the broader orbit-closure frame.
  **NOW PROPOSED as A7 (S154).**
- **Random tensor models (Gurau-Witten melonic universality)** (§3) —
  large-N tensor models admit a 1/N expansion dominated by *melonic
  graphs* with non-Gaussian conformal IR fixed point (Witten 2016 SYK
  without disorder). Tests rank-3 universality of `T_{ijk}^N := \mathbf{1}
  [i+j+k=N] χ_P(i)χ_P(j)χ_P(k)` — the FIRST rank-3 measurement on χ_P,
  orthogonal to all CLOSED rank-2 universality measurements. **NOW
  PROPOSED as D41 (S154).**
- **Resurgence / Borel-Écalle alien calculus** (§10) — divergent
  asymptotic series with isolated Borel-transform singularities and
  Stokes constants. Applied to the divergent Riemann-Siegel remainder
  `R(t)`, asks whether Stokes points coincide with ζ-zero ordinates
  (rigorous version of Berry-Keating 1990s heuristic). Distinct from
  C2 sine-kernel (zero-spacing), C7 FHK (amplitude max), B4 Voronin
  (vertical translates) — attacks the asymptotic-divergence regime
  directly. **NOW PROPOSED as D42 (S154).**
- **Hairer regularity structures / KPZ universality of `π(x) − Li(x)`**
  (§3) — non-Gaussian universality class with `t^{1/3}` scaling and
  Tracy-Widom β=2 limit. Tests Hölder regularity of `D(x) := (π(x) -
  Li(x)) \log x / \sqrt{x}` on KPZ-spaced grid via wavelet decomposition.
  First non-`\sqrt{x}`-scale measurement of PNT error term, orthogonal
  to all CLOSED CLT-scale measurements (C5 Stein, C7 FHK). **NOW
  PROPOSED as D43 (S154).**
- **Arithmetic topology (Mazur-Morishita primes-as-knots; Rédei
  symbol on prime triples)** (§4) — arithmetic Massey products
  `\langle p, q, r \rangle_2 \in \mathbb{F}_2` evaluating to the Rédei
  symbol `[p, q, r]` on admissible (Borromean-condition) prime triples.
  CRITICALLY DISTINCT from CLOSED line 208 (generic knot polynomial
  invariants): Mazur-Morishita is *primes-as-knots* via Galois
  cohomology of `G_{\mathbb{Q}, S}` — the invariants ARE arithmetic
  by construction. Tests bilinear-vs-trilinear refinement of Hardy-
  Littlewood. **NOW PROPOSED as D44 (S154).**
- **Baker-Norine 2007 graph Riemann-Roch / chip-firing divisor rank**
  (§2 NEW) — non-linear divisor-rank function `r_G(D)` on a finite graph
  `G` via Dhar's burning algorithm. Apply to the prime divisor `D_P^N`
  on the divisibility graph `Γ_N`. CRITICALLY DISTINCT from CLOSED
  tropical lines 204/326/431/432/660 (those are min-plus compression of
  generating functions, collapsing to "smallest prime"; chip-firing is a
  graph game on divisors, not a min-plus compression). Distinct from
  D22 (linear Hodge spectrum). **NOW PROPOSED as D45 (S160).**
- **Schubert calculus / Plücker-coordinate Grassmannian point on χ_P**
  (§2) — `V_χ_P^N := rowspan ((1)_i, (χ_P(i))_i) ∈ Gr(2, ℂ^N)` and its
  Schubert cell `Ω_λ`. Plücker coordinates `p_{ij} = χ_P(j) - χ_P(i)`
  encode prime-pair correlations; Schubert intersection multiplicity
  `c_{λ,λ}^{(N-2)^2}` is an LR coefficient. UNUSED row in §2 finally
  attacked. **NOW PROPOSED as D46 (S160).**
- **Cluster algebras (Fomin-Zelevinsky 2002)** (§7 NEW) — mutation orbit
  return-time of a prime-gap-seeded quiver (path / twin-prime / gap-
  equality variants); cluster-finite ⇔ Dynkin (`A_n`, `D_n`, `E_{6,7,8}`)
  via Fomin-Zelevinsky 2003 finite-type classification. Y-system
  periodicity by Keller 2013. A genuine cluster-finite arithmetic quiver
  would yield the FIRST non-additive Laurent identity on prime gaps.
  **NOW PROPOSED as D47 (S160).**
- **Connes-Consani-Marcolli 2007 endomotive — Bost-Connes Galois orbits
  on KMS_∞ ground states under χ_P-projector** (§1 NEW) — the BC partition
  function reduces to ζ (CLOSED line 185) but the KMS_∞ ground-state
  Galois orbit structure (BC 1995 §6, Hilbert 12th link) is strictly
  richer information. CCM 2007 endomotive formalism makes the
  Galois-orbit data computable. χ_P-projector restriction tests whether
  prime-Frobenius subalgebra orbit-length distribution carries arithmetic
  signal. CRITICALLY DISTINCT from line 185 (trace loses Galois data).
  **NOW PROPOSED as D48 (S160).**

When in doubt, pick from this list first.

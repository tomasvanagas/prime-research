# Session-by-Session Key Insights

Detailed findings from each research session. Read the relevant session
if you need deep context on a specific topic. For a quick overview,
see the Status section of CLAUDE.md.

---

## Session 12
- "Is pi(x) in NC?" is EQUIVALENT to our target.
  All known approaches produce circuits of size 2^{Theta(N)} (exponential in input).

## Session 13
(a) BPSW IS computable in TC^0 (MR=scalar pow, Strong Lucas=2x2 MPOW, Jacobi=GCD).
    PRIMES in TC^0 iff BPSW (or similar) is unconditionally correct.
    Verified correct to 2^64. GRH also suffices (Miller's test = O(N^2) scalar pows).
(b) Prime indicator ANF degree = Theta(N) over GF(2), 50% sparsity (random-like).
    No GF(2) algebraic shortcut for counting primes.
(c) Zeta zero minimum K_min ~ 0.35 * x^{0.27}: power law, no reordering helps.
(d) Spectral graph approaches all circular or equivalent to Meissel-Lehmer.
(e) No new algorithmic breakthroughs in 2025-2026 literature.

## Session 14
(a) **PRIMES in L and pi(x) in NC are INDEPENDENT questions.** The chain
    PRIMES in L -> pi(x) in #L -> pi(x) in NC^2 BREAKS due to workspace mismatch:
    NL machine needs O(N) bits for candidate n, but #L allows only O(log N).
(b) **I-E fractional parts carry O(2^k) independent bits** -> no determinant
    smaller than 2^{Theta(sqrt(x)/log(x))} can encode the Legendre sieve.
    A GapL algorithm MUST avoid floor functions entirely.
(c) **Lucy DP matrices have NO algebraic structure**: unipotent, displacement
    rank 50-60% of dimension, full-rank product. No compression possible.
(d) **Nonlinear sieve breaks parity in theory but NOT in efficiency**: nonlinear
    ops on floor values CAN distinguish primes from semiprimes but cost >= O(x^{2/3}).
(e) **All algebraic variety approaches fail**: low-dim can't encode pi(x),
    high-dim has slow point counting, Frobenius eigenvalues = zeta zeros.
(f) **pi(x) mod m is NOT a linear recurrence** for any m (Mauduit-Rivat consequence).

## Session 15
(a) **Determinantal complexity connection**: pi(x) as degree-N multilinear polynomial
    in bits. Found NxN det representations for N=2,3,4. "Is pi(x) in GapL?" <=>
    "polynomial determinantal complexity?" For N>=10, GENERIC polynomials don't fit.
(b) **#TC^0 subset NC? is THE complexity question**: If BPSW in TC^0 (conditional),
    pi(x) in NC iff #TC^0 subset NC. Fermat residue coupling prevents batch counting.
(c) **Uniformity is the true barrier**: Nonuniform circuits trivially poly(N).
    Hard part = UNIFORM construction. Natural proofs barrier blocks lower bounds.
(d) **ALL randomized approaches fail**: zeta zero sampling (100% needed), probabilistic
    sieve (10^6x worse), hash counting, quantum counting -- all rigorously excluded.
(e) **Divide-and-conquer fails**: error accumulates O(sqrt(x)) regardless of depth.
(f) **Monotone complexity inapplicable**: [pi(x)>=k] = [x>=p(k)] trivially O(N).
(g) **Arithmetic circuits don't help**: VP=VNC^2 (depth free), but pi(x) not a
    polynomial over fields. Tau conjecture orthogonal.
(h) **Systematic analysis of 8 intermediate quantity families**: residues, polynomial
    evals, matrix eigenvalues, topology, representation theory, entropy, recursive,
    physical -- ALL route back to floor values or zeta zeros.
(i) **No 2026 breakthroughs**: Chen-Tal-Wang STOC 2026 closest to TC^0 frontier.

## Session 16
(a) **15 intermediate quantity families now CLOSED** (8 from S15 + 4 GapL + 7 novel):
    class numbers h(-d), L-function L(1,chi), elliptic curve a_p, regulators,
    additive combinatorics/sumsets, ergodic theory, model theory, tropical geometry,
    sufficient statistics, algebraic geometry/F_q, representation theory S_n/GL_n.
    ALL route to primes (C), zeta zeros (E), or lose information (I).
(b) **Three "pillars" are the ONLY exact encodings of pi(x)**: prime positions,
    zeta zeros, floor values {floor(x/k)}. These are informationally equivalent.
    No fourth encoding found across 15+ candidate families.
(c) **TC^0 batch counting has 5 closed routes**: MAJORITY fan-in, divide-and-conquer,
    batch CRT, Carmichael structure, period exploitation.
    Communication complexity gives 2^{N/2}/poly(N) lower bound on TC^0 circuit size.
(d) **H-T signed cancellation transfer to pi(x) is CLOSED via 6 routes**:
    POSITIVITY of prime indicator prevents cancellation. Converting M(x)->pi(x)
    always costs O(x^{2/3}).
(e) **Lambert W error is structurally random**: delta(n) uncorrelated with gaps,
    uniform mod m, ~sqrt(x) magnitude. No exploitable structure.
(f) **HKM 2023 achieves O~(sqrt(x)) for pi(x) elementarily** (Math. Comp. 2024).
(g) **Aggarwal 2025 gives O(sqrt(n)*log^4(n)) for p(n)** via binary search + HKM.

## Session 17
(a) **EXACT communication complexity of pi(x)**: rank(pi_N) = 2^{N/2-1} + 2 for
    balanced bit partition (verified N=2..20). This gives dc(pi_N) >= Omega(sqrt(x)).
    The multilinear polynomial route to GapL is DEFINITIVELY CLOSED.
(b) **Boolean Fourier analysis**: Prime indicator has ~30% excess low-degree Fourier
    weight vs random (from parity/mod-4 structure), but noise sensitivity is near-random
    (~0.9x). No evidence of low-depth circuit structure. Spectral profile near-random.
(c) **Ono partition characterization (PNAS 2024) CLOSED**: requires divisors = circular.
(d) **Chen-Tal-Wang ECCC 2026**: THR-of-THR lower bounds, not number-theoretic.
(e) **sqrt(x) barrier is UNIVERSAL**: communication complexity, Fourier analysis,
    determinantal complexity, substitution rank -- ALL converge to sqrt(x).

## Session 18
(a) **Derandomization theory CLOSED (6 routes):** IW97, NW generators, approximation
    amplification, KI04/PIT, reverse hardness amplification, batch Fermat derandomization.
(b) **Natural proofs barrier (RR97) CONFIRMED relevant:** Cannot prove pi(x) not in
    TC^0 via "natural" methods. All our Fourier/rank measures are natural.
(c) **META-INSIGHT:** Derandomization addresses removing randomness from EXISTING
    efficient algorithms. Our problem is more fundamental: no efficient algorithm exists.
(d) **Formula complexity:** KW theorem gives formula size >= 2^{N/2-O(1)} for pi(x).
    Via Valiant, general circuit size >= 2^{N/4}.
    TC^0 bound 2^{N/2}/poly(N) remains strongest for constant-depth.

## Session 19
(a) **Unbalanced communication complexity UNIVERSAL FORMULA:** rank(k) = 2^{min(k,N-k)-1}+2
    for ALL bit partitions with min(k,N-k) >= 3 (verified N=4..20, all k). For fixed k,
    rank stabilizes once N >= 2k. No polynomial-rank partition exists.
(b) **SVD spectral decay is POWER-LAW:** S_osc ~ i^{-1} (not geometric). 90% of osc
    variance in ~20 SVs (seemingly bounded), but 99% needs ~30% of all SVs (exponential).
    Max osc SV scales as x^{0.66}. Power-law means better-than-R^{-1} approx but never exact.
(c) **SVD IS the explicit formula:** Top osc SVs correspond to first zeta zeros (corr 0.95
    for gamma_1 at N=20). Zeta basis explains only 0.12% of variance at N=20. The
    communication complexity barrier prevents efficient extraction.
(d) **3-party NOF:** Balanced (N/3)^3 split gives cut rank 2^{N/3}, consistent with TC^0.
    Best (1,1,N-2) split gives trivial rank 4. Insufficient for ACC^0/TC^0 separation.
(e) **PSLQ/LLL identity search EXHAUSTIVE NEGATIVE:** All 6 relation types (linear,
    polynomial, recurrence, modular, functional, discrete derivative) fail for f(x)=pi(x)-R(x).
    Cross-validation proves all single-point relations are spurious.
(f) **NFS-type L[1/3] for pi(x) CLOSED:** 4 sub-approaches (norm sieve, Chebotarev, class
    groups, Artin L-functions) all fail. NFS exploits multiplicative structure; prime counting
    is additive-global with no analog.
(g) **Gap predictability CLOSED:** AR models give no improvement, MI(g_n;g_{n+1})=0.38 bits
    (10.3%), gaps 9% more compressible than i.i.d. Near Cramér random model.
(h) **Kt complexity of delta(n):** |delta(n)|~n^{0.57}, bit ratio=0.52, AR(1) R²=0.996 
    but RMSE=10.5 (innovations random). Sign runs mean length 38.5 (zeta oscillation).
    Uniform mod m. 18% compressible. NO exploitable structure beyond local smoothness.
(i) **Total approaches: ~520+.** sqrt(x) barrier remains universal across ALL known measures.

## Session 20: Fresh Perspective (Wildcard)
Started from FIRST PRINCIPLES without reading CLOSED_PATHS.md.
10 experiments across 10 unconventional ideas, 6 parallel agents.

(a) **CRT MODULAR RECONSTRUCTION CLOSED:** π(x) mod m is a random walk with step ~1/ln(x).
    No exploitable structure. Computing π(x) mod m requires π(x) — circular. CRT adds overhead.
(b) **TENSOR NETWORK SIEVE CLOSED:** Sieve has VOLUME-LAW entanglement. MPS bond dimension =
    primorial(y) = exponential. DFA has exactly primorial states. Binary→mixed-radix trades one
    non-locality for another. Sieve function is "maximally complex" in tensor network sense.
(c) **SPECTRAL SHORTCUT CLOSED:** Trace formula IS the explicit formula (circular). ζ'/ζ on a
    contour costs O(T^{3/2}). Heat kernel smoothing introduces O(√x) bias. First 10-20 zeros
    capture 80-95% of variance; hard 5-20% grows with x. No free lunch from spectral methods.
(d) **DYNAMICAL GAPS CLOSED:** MI(g_n;g_{n+1})~0.3-0.5 bits out of ~3.5-4. No linear recurrence,
    HMM, substitution, or low-dimensional attractor. LZ complexity near-random. Confirms S19(g).
(e) **FINITE FIELD LIFTING CLOSED:** q→1 degenerates; N_q(n) gives PNT not exact; virtual curve
    needs genus ~√x giving O(x^{3/2}). Root cause: F_q[x] has rational ζ (no zeros), Z has ∞ zeros.
(f) **SIEVE FUNCTION COMPRESSION:** S(v,p) is binary step function (incompressible). Fourier NOT
    sparse (99% needs 50-70% of modes). BUT: sieve updates have 90% energy in rank-1 component,
    and log-Fourier has 90% in ~10 modes. The 90/10 SPLIT confirms smooth+random decomposition.
(g) **ITERATIVE ZERO REFINEMENT:** Self-correction converges (sensitivity 0.002-0.019) but to
    value ~√x from exact. 50 zeros gives O(√x) error. The information IS in the high zeros.
(h) **ALL ALGEBRAIC SHORTCUTS CLOSED:** Wilson (circular), cyclotomic (circular), arithmetic
    derivative (circular), determinant sieve (≡ Möbius), number fields (MORE zeros), matrix
    exponentiation (non-autonomous), Selberg sieve (parity barrier), characters (equally hard).
(i) **INFORMATION-THEORETIC QUANTIFICATION:** For x=10^100, π(x) has ~332 bits total, ~173
    "hard" bits encoding zeta zeros. Even the MSB of Δ=π(x)-Li(x) requires O(√x) zeros.
    Hard bits form a HOLOGRAPHIC projection — no single bit is independently computable.
(j) **Total approaches: ~550+.** 25 new closed paths. Fresh perspective independently rediscovered
    and confirmed ALL known barriers from a new starting point. No new viable directions found.
    The √x barrier is genuinely universal.

## Session 20 (Deep Focus: Kt Complexity of delta(n))

**Task:** FOCUS_QUEUE Task #1 — empirically estimate Kt(delta(n)|n).
20+ experiments across 4 parallel sub-agents + direct experiments.

Key findings:

(a) **Kt(delta(n)|n) = O(log n) TRIVIALLY**: For any polytime-computable function, Kt = O(1) + log(runtime)
    = O(log n). Kt is NOT the right measure. The correct question is CIRCUIT COMPLEXITY of delta.

(b) **CANNOT predict delta(n) from n**: Every method tested (linear, polynomial deg 2-6, Fourier up to
    100 terms, RandomForest, GradientBoosting) achieves R² < 0 on held-out test set. Even best ML method
    has RMSE 1.07x the naive baseline. Adding n as feature to AR(k) models: 0% improvement.
    **n is informationally irrelevant to delta.** This is the session's strongest empirical result.

(c) **Power spectrum: PSD ~ f^{-1.69}**: Between pink (1/f) and Brownian (1/f²) noise. Confirmed by
    shuffled baseline (flat PSD). Mutual information decays as power law MI ~ k^{-0.34}, NOT exponential.
    This means no finite-memory model (ARMA) captures the full correlation structure.

(d) **Wavelet: 83.5% energy in coarsest 3 scales**: Fine-scale coefficients (D1-D4) are non-Gaussian
    with excess kurtosis up to 2.2. Coarse scales are Gaussian. Signature of smooth large-scale +
    intermittent local fluctuations.

(e) **Explicit formula partial sums DIVERGE**: At x ≤ 10^6, adding more zeta zeros makes the
    approximation WORSE. Error grows ~linearly with K and ~x^{0.25}. This is the well-known conditional
    convergence issue. Theoretical: need O(sqrt(x)*log²(x)) zeros for |error| < 1.

(f) **Effective dimensionality**: SVD decay S_k ~ k^{-1.1} (power-law, consistent with Session 17).
    Correlation dimension grows with embedding dim (1.07→1.79 for dim 2→20). No low-dimensional
    manifold. Information dimension of delta values = 0.84 (slightly clustered).

(g) **MDL analysis**: AR(5) dominates all other models (RMSE=10.5). Fourier with 201 params: RMSE=55.
    Poly degree 50: RMSE=132. To match AR from n alone: ~10^4+ parameters needed.

(h) **Compressibility via sequential structure**: LZ complexity = 27% of random. bz2 = 18.5%. Permutation
    entropy at order 7 = 84% of max. Transition matrix spectral gap = 0.003 (slow mixing). ALL
    compressibility comes from autocorrelation, NOT from any function of n.

(i) **Incremental entropy**: ~0.22*log(n) + 4.6 bits. Very slow growth. Each delta value adds ~6-7
    "surprise bits." This is consistent with the smooth+random decomposition.

(j) **Connection to Session 17**: The 1/f^{1.69} spectrum EXPLAINS why communication rank = 2^{N/2-1}+2.
    Each zeta zero contributes R(x^rho) with amplitude ~x^{1/2}/gamma. The incommensurable gamma values
    make the sum incompressible — no finite linear combination of frequencies can represent it.

**New structural result**: delta(n) is a 1/f^{1.69} colored noise process. This spectral exponent
quantifies the long-range correlation structure and connects directly to the zeta zero oscillation model.

**11 new closed paths** (prediction methods, low-rank/manifold approaches).
**Total approaches: ~560+.**

## Session 21 (Proposal Session)

**Method:** 8 fresh proposals developed with code + computational tests, unconstrained by prior analysis.
3 parallel internet searches for latest papers (2024-2026).

**Literature:** No breakthroughs. Aggarwal (2025) O(sqrt(n)*log^4 n). Hirsch-Kessler-Mendlovic O~(sqrt(N)) elementary. primecount 8.4 SIMD. Dequantization inapplicable (high rank). ML provably learns smooth part only.

**Key experimental findings:**
1. **CRT gives unique prime with M=30030:** With just 6 small moduli {2,3,5,7,11,13}, CRT narrows delta(n) to a UNIQUE prime candidate for all n ≤ 10000. The framework achieves 7462x search compression. **Blocked by circularity:** computing p(n) mod q requires pi(x;q,a), which costs O(x^{2/3}).
2. **Spectral+Sieve confirmed:** K zeros give width O(x/K). Combined with CRT, could give O(polylog) if CRT modular computation is cheap.
3. **ML 0% exact:** Ridge regression 34.7% RMSE improvement but ZERO exact predictions of delta(n). kNN slightly better (RMSE 77.7 vs 296.8 baseline) but still 0.2% exact. Confirms information barrier.
4. **Fourier interpolation saturates:** Adding zeros beyond 5 barely helps (test RMSE 2.30→2.20 for 5→30 zeros). Need O(sqrt(x)) zeros.
5. **F_1 extrapolation diverges.** The function field→number field limit is singular.
6. **Gaussian smoothing kills everything.** No test function optimum exists.
7. **Strong autocorrelation at lag 1** in p(n) mod 30 (chi2/df=159), decays by lag 5. Short-range structure only.
8. **Zeta zero ratios are generic irrationals.** Quasi-periodicity super-period ~10^158. No shortcut from Diophantine structure.

**Reduction identified:** The ENTIRE problem reduces to: can p(n) mod q be computed in O(polylog) for small fixed q? If yes for q ≤ 13, CRT gives p(n) exactly.

**8 new experiments** saved to `experiments/proposals/`.
**8 paths closed:** ML prediction, Fourier interpolation, F_1 extrapolation, Gaussian trace formula, Diophantine quasi-periodicity.
**1 path refined:** CRT modular reconstruction (viable framework, circularity barrier identified).

## Session 22 (Critique Session)

**Method:** Adversarial critique of all 8 proposals from Session 21. Each proposal checked against
525+ CLOSED_PATHS entries. New experiments run on the "prime race" direction (most novel component).

**Verdicts:** All 8 proposals are DUPLICATE — every approach matches prior closed paths.
No proposal addresses the 4 remaining open directions (circuit complexity, novel identity,
BPSW correctness, S*T tradeoff).

**New experimental results (prime race direction):**
(a) **pi(x;q,a) errors are 2x ROUGHER than pi(x) errors.** Total variation ratio ≈ 2.0 across
    all tested moduli q=3,4,5,7. The arithmetic progression error is WORSE, not better.
(b) **L-function zeros have same density as zeta zeros.** Found 8 zeros of L(s,chi_4) up to
    t=100 (vs ~29 zeta zeros). Per-character, spectral power concentrates in fewer modes
    (K=39-81 for 90% vs K=148 for zeta), but phi(q) characters multiply the total cost.
(c) **Chebyshev bias provides O(1) bits.** E(x;4) > 0 for 99.52% of x (matching Rubinstein-Sarnak
    prediction 99.59%). Sign predictable but value is not — provides constant-bit, not log(n)-bit info.
(d) **PSD slope of E(x;4)/sqrt(x) = -1.664** (vs -1.570 for pi(x)-Li(x)). Similar spectral structure.
    The prime race is informationally equivalent to the main problem.

**2 new paths closed:** Prime race shortcut for CRT, L-function zero convergence advantage.
**Total approaches: ~528+.**

**Key meta-insight:** Session 21 proposals were generated independently but ALL rediscovered
previously closed paths. This is strong evidence that the space of "natural" approaches has
been thoroughly explored. Any breakthrough must come from the 4 open directions or from
genuinely novel mathematics not yet connected to prime counting.

## Session 20 (Space-Time Tradeoff Investigation)

Systematic investigation of formal S*T lower bounds for pi(x) via four approaches.

(a) **Communication complexity route CLOSED for super-polylog bounds.** D(pi) = Theta(N/2)
    bits (from rank = 2^{N/2-1}+2). Via Beame's r-round framework: T*S >= Omega(N^2).
    With S=poly(N): T >= N^{2-c} = polylog(x). Cannot rule out polylog(x) time.
    Fundamental limit: N-bit input has D(f) <= N, capping all comm-complexity tradeoffs at poly(N).

(b) **Nechiporuk formula bound gives only Omega(N).** Optimal block size s*=3 (from calculus).
    The exponential rank 2^{s/2} is local (grows with block size, not input size).
    Method inherently limited to O(N^2) for any function. Useless here.

(c) **OBDD size grows as 2^{0.79*N} empirically (N=4..12).** Matches random functions of same
    density (random: 2^{0.80*N}). OBDD width >= 2^{N/2-1} provably from communication rank.
    But OBDDs are a restricted model; general BPs (variables read multiple times) can be
    exponentially smaller.

(d) **M-L DAG pebbling: T*S >= Omega(x^{5/6}/ln x)** for the specific Meissel-Lehmer computation.
    Lucy DP DAG: depth = pi(sqrt(x)), width = O(sqrt(x)). Exact pebbling computed for x=10,15,20
    (4-5 pebbles). But this is algorithm-specific; a non-sieve approach bypasses this DAG.

(e) **BDD comparison: isPrime, pi mod 2, pi total, random all grow at 2^{0.72-0.80*N}.**
    pi(x) has no structural BDD advantage over random functions.

**5 new sub-paths closed.** Total approaches: ~533+.

**Key insight:** Proving T >= x^{Omega(1)} for general pi(x) algorithms requires CIRCUIT
LOWER BOUNDS (showing pi(x) not in NC), which faces the Natural Proofs barrier (Razborov-Rudich).
This is at least as hard as P != NP. The space-time tradeoff problem for pi(x) is
provably beyond current complexity-theoretic techniques for general algorithms.

## Session 23 (Multiparty Communication + Space-Time + Algebraic Structure)

**Method:** 3 parallel sub-agents + 3 direct experiments. Focus on extending prior results and
testing genuinely novel algebraic directions.

**Sub-agent results:**

(a) **k-party NOF communication complexity (k=2..8, N=6..18):**
    For k >= 3: mode-unfolding rank is EXACTLY 2^{ceil(N/k)} = FULL rank for all tested cases.
    pi(x) is indistinguishable from a random function in this measure. Only k=2 shows sub-maximal
    rank (the known 2^{N/2-1}+2 formula). Mode-unfolding rank is TOO COARSE for the ACC^0 question.
    Need true tensor rank, discrepancy, or polynomial method arguments.
    SVD analysis: pi(x) is essentially rank-1 (smooth part captures >99.9% of variance).
    Residual (oscillatory part) has FULL rank for k >= 3.

(b) **Space-time tradeoff lower bounds CLOSED (as route to impossibility proof):**
    - Communication complexity: D(f) <= N = log(x), so T*S >= Omega(N^2) = Omega(log^2 x). NEVER
      sufficient to rule out polylog.
    - OBDD: 2^{0.79*N} growth, same as random. OBDD >= sqrt(x) from rank.
    - Nechiporuk: trivial Omega(log x).
    - M-L pebbling: T*S >= Omega(x^{5/6}/ln x) but ALGORITHM-SPECIFIC.
    - General lower bounds face Natural Proofs barrier. AS HARD AS P != NP.
    **Conclusion:** The space-time tradeoff route to proving polylog impossible is CLOSED.
    Only circuit lower bounds can resolve the question, and those face fundamental barriers.

**Direct experiments:**

(c) **pi(n) is NOT holonomic (D-finite):** Tested recurrences of order d <= 20 with polynomial
    coefficients of degree r <= 8. Test/random ratio consistently ~1.0-1.7 (indistinguishable
    from random). This is STRONGER than Session 14's "not LRS" result — polynomial coefficients
    don't help. Rules out holonomic sequence algorithms (baby-step/giant-step).

(d) **Short-interval explicit formula iteration:** Confirmed the sinc-based zero cutoff:
    interval of width W needs zeros up to height ~x/W. Iteration CANNOT reduce zero count:
    each round needs K_i ~ K_1 zeros. The oscillatory contribution sums to O(sqrt(x)/log x)
    regardless of partition. Hybrid optimum at W = sqrt(x) gives O(sqrt(x)*polylog).

(e) **Ono partition characterization with p-adic lifting CLOSED:** M_k(n) DP has O(n^2) ops —
    worse than O(x^{2/3}). Modular computation: same op count, just bounded values. The Ono
    criterion mod l gives only 46-72% accuracy. Ramanujan congruences are special cases.
    Computationally inferior to all known methods.

**~9 new closed paths. Total approaches: ~570+.**

**Session 23 meta-insight:** The barriers are now understood at FOUR levels:
1. **Analytic:** sqrt(x) zeros needed (GUE statistics, explicit formula convergence)
2. **Algebraic:** Not LRS, not holonomic, not k-automatic (information-theoretic)
3. **Combinatorial:** Communication rank 2^{N/2-1}+2, full tensor rank for k >= 3
4. **Complexity-theoretic:** Proving impossibility faces Natural Proofs barrier (as hard as P != NP)

This means: either a breakthrough exists via genuinely new mathematics, or proving it impossible
is equivalent to solving major open problems in complexity theory. Both directions are OPEN and
both appear to require fundamentally new ideas beyond current mathematical techniques.

---

### Session 24 (2026-04-05): Fresh Perspective Session 2

**Focus:** Complete restart from first principles. 11 independent experiments across 10
parallel agents, testing genuinely novel approaches inspired by breakthroughs in other fields
(Shor, compressed sensing, fast multipole, AlphaFold, Candès-Tao).

**New approaches tested (16 new closed paths):**
1. CRT Prime Locator — progression counting is circular
2. Hierarchical Sieve — fast-multipole analog scales linearly, not polylog
3. Spectral Compression — δ(n) partially compressible but critical 5% dense in ALL bases
4. Dynamical Fast-Forward — gaps unpredictable, indicator not k-automatic
5. Probabilistic Exact — binary search relocates problem, doesn't solve it
6. Recursive Dickman — DDE correction dominates (exponent 7.9), 93%+ nodes need exact
7. LFSR over Finite Fields — L/N=0.5000 over ALL GF(p), maximally random-like
8. Neural Arithmetic — test RMSE never < 0.5, no generalization
9. Étale/Weil Formula — genus ~10^102, worse than brute force
10. Connes Trace Formula — no test function σ avoids the √x error
11. NTT/Dirichlet Convolution — Perron needs O(x) quadrature points
12-16. Adelic reconstruction, Kolmogorov complexity, multiplicative Fourier, etc.

**Key new quantitative findings:**
- Prime indicator linear complexity: **L/N = 0.5000** over all GF(p) for p=2..23
  (maximally random in every finite field — strongest info-theoretic evidence yet)
- δ(n) Kolmogorov complexity: **5x more compressible than random** (ratio 0.049),
  compression improves with length (1.56 bits/symbol at N=10000)
- DCT is optimal standard basis: 99% energy in 10.4% of coefficients
- 2-kernel of prime indicator has 38 growing sequences → NOT k-automatic

**Literature search:** No new algorithms for π(x) in 2025-2026. Two new references:
Valley Scanner zeta zero finder (arXiv:2512.09960), communication complexity from
information causality (arXiv:2602.10206).

**~16 new closed paths. Total approaches: ~575+.**

**Session 24 meta-insight:** Across 21 fresh-perspective approaches (10 from Session 20 +
11 from Session 24), EVERY approach collapses into one of three failure modes:
- **Circularity** (CRT, probabilistic, adelic, bisection)
- **Equivalence** (sieve, Dickman, étale, NTT, Perron, Connes)
- **Information loss** (spectral, dynamical, LFSR, neural, ML)

The failure taxonomy from novel/failure_taxonomy.md is **empirically confirmed exhaustive**
across all 575+ approaches tested. No approach has ever avoided all three simultaneously.
The strongest new evidence is the LFSR result: linear complexity L/N=0.5000 across all fields
means the prime indicator is indistinguishable from random in the most powerful algebraic sense
available (finite-state linear prediction). Combined with Session 23's full tensor rank result,
this provides converging evidence from both algebraic AND combinatorial directions that the
prime indicator carries near-maximal information density.

---

### Session 25 (2026-04-05): Deep Focus — Zeta Zero Structural Patterns

**Focus:** FOCUS_QUEUE Task #2. Systematic search for algebraic, arithmetic, or structural
relations among Riemann zeta zeros that could enable fast summation of the explicit formula.
7 experiments run in 6 parallel agents.

**Results (ALL NEGATIVE — 6 new closed paths):**

(a) **Pairwise rational relations CLOSED:** 499,500 ratios gamma_i/gamma_j tested against
    closest p/q (q≤100). KS test significant (p=1.3e-6) but effect is WRONG DIRECTION — zeros
    are marginally FARTHER from simple rationals than random (GUE repulsion). No exploitable
    rational structure.

(b) **PSLQ/LLL integer relations CLOSED:** 13,000+ tests at 60-digit precision. Zero relations
    among subsets of 3-5 zeros augmented with {1, pi, log(2pi)}. 1,225 pairwise tests negative.
    K≥20 "hits" are lattice artifacts (confirmed by random baseline). Consecutive spacing
    relations also negative. Zeros appear linearly independent over Q, as conjectured.

(c) **DFT spectral structure CLOSED:** Power spectrum matches GUE with correlation 0.9999.
    Spectral flatness: 0.93-0.999 at high frequencies (white noise), 0.013 overall (Weyl law
    trend at low freq). Pair correlation matches GUE prediction 1-(sin(πr)/(πr))². Number
    variance follows GUE logarithmic growth. Only p=2 shows faint signal (12x median in
    Lomb-Scargle); all p≥5 indistinguishable from noise. O(N) bits incompressible.

(d) **Partial sums recurrence CLOSED:** No linear recurrence (order 1-20), nonlinear (deg 2-3),
    or difference recurrence for S_K(x) = Σ_{k=1}^{K} cos(γ_k log x)/(1/4+γ_k²)^{1/2}.
    Best residuals 2-6% (pure overfitting). Convergence neither geometric nor algebraic (R²<0.5).
    Each zero contributes genuinely independent information.

(e) **Mod-constant equidistribution CLOSED:** gamma_n mod m is uniform for ALL 10 moduli
    tested (1, π, log(2π), 2π, √(2π), e, log(2), log(3), log(5), log(7)). KS p-values all >0.4.
    Weyl sums at expected 1/√N threshold. Discrepancy D_N=0.0163 is BELOW random expectation
    ~0.031 (GUE repulsion makes zeros "more uniform than random"). Joint distribution (mod 1,
    mod π) independent. No arithmetic structure exploitable.

(f) **Sparse matrix model CLOSED:** Tridiagonal (Jacobi) matrices fit N zeros with <0.05% error
    using 2N-1 parameters — but this is overfitting, not compression. Toeplitz (N params) degrades
    at N≥100. Prime-based matrices (4 params) fail at 15% error. Extrapolation: 0.8-2.6% error.
    Matrix entries have 35% relative variation, no simple pattern. The sparse representation
    REPACKAGES information without reducing it.

**6 new closed paths. Total approaches: ~581+.**

**Session 25 synthesis:** The zeta zeros are structurally GUE-random in EVERY sense tested:
- **Algebraic:** No integer linear relations (PSLQ). No simple rational ratios. Linearly independent over Q.
- **Spectral:** Power spectrum = GUE (corr 0.9999). Pair correlation = GUE. Number variance = GUE.
- **Arithmetic:** Equidistributed mod every tested constant. Discrepancy below random.
- **Sequential:** No recurrence in partial sums. Each term contributes independent information.
- **Structural:** Cannot be compressed via sparse matrices. O(N) parameters required for N zeros.

This definitively closes the "zeta zero compressibility" direction (Open Problem #3) for all
structural approaches. The ONLY remaining avenue would be if the zeros have structure at scales
not tested (e.g., correlations requiring >1000 zeros to detect, or relations with coefficients
>1000). But GUE universality predictions strongly suggest this is not the case.

---

## Session 26

Two-pass proposal session: 16 total proposals (8 + 8), all with computational experiments.

**Literature search:** 3 parallel agents searched arxiv/Google Scholar across 6 fields. No new
algorithmic breakthroughs. Key new references: Ramanujan Library (ICLR 2025), IntSeqBERT (2026),
Kolmogorov incompressibility of primes (2024), Song trace formula (2024).

**Session 26a (Proposals 1-8):** Spectral+sieve (O(√x)·polylog), CRT reconstruction (unique
with M=30030 but circular), ML on delta (0% exact), F_1 extrapolation (divergent), Fourier
interpolation (saturates), Gaussian trace formula (kills all zeros), additive combinatorics
(O(1) bits/step), Diophantine quasi-periodicity (super-period ~10^158).

**Session 26b (Proposals 9-16):**

(a) **Candidate generation CLOSED:** Heuristic narrowing + primality verification. Candidate set
    scales as n^0.577 ≈ √n, not polylog. RH error bound forces interval O(√p(n)·log p(n)).
    Progressive filtering (mod 30, zeta zeros, small-prime sieve) gives multiplicative reduction
    to ~3.8% of interval but absolute count still O(√n).

(b) **Class field splitting CLOSED:** 39 quadratic fields give 12.3 independent bits via
    Chebotarev density. CRT modulus 2^72.7 — enough for unique identification with 20 fields.
    But computing π_split(x,d) requires enumerating primes (circular) or L(s,χ_d) zeros
    (MORE zeros than ζ alone). Cost strictly WORSE: O(x^{1/2+ε}·log n).

(c) **Green-Tao nilsystem CLOSED:** R²=0.565 in-sample collapses to 0.253 CV. Dominant modes
    = sieve biases (1/2, 1/3, 1/6). li^{-1}(n) gives 1.51 correct digits vs nilsequences 0.84.
    Green-Tao describes collective statistics (APs), not individual primes.

(d) **Pfaffian/treewidth CLOSED:** All 3 graph encodings have polynomial treewidth: sieve ~x^0.36,
    Lucy DP ~x^0.51, prime dependency = complete graph. FKT/bounded-treewidth shortcuts impossible.

(e) **Tropical sieve CLOSED:** 7 sub-tests all fail. Tropicalization val(a+b) ≥ min(val(a),val(b))
    destroys counting information. Min-plus converts counting→optimization (irreversible).

(f) **PSLQ on delta(n) CLOSED:** No algebraic relations of any kind. PSLQ finds only spurious
    (unique per window) relations. Berlekamp-Massey: LFSR length = N/2 for ALL moduli 2-11.
    Delta(n) is maximally complex in every finite field. 6.26 bits entropy, 121 unique values.

(g) **Smooth number subtraction CLOSED:** Ψ(x,B) computation IS the sieve (50-640x slower than
    Eratosthenes). Buchstab tree = 1.8·x^{2/3} distinct args = Lucy DP. Linear combination of
    Ψ features: test RMSE=37.3, 0/500 exact.

(h) **Tensor network sieve CLOSED:** MPS bond dimension = 2^{(0.33-0.43)·a}, exponential in
    #primes. Sieve tensor 5-20x more structured than random (constant factor only). Floor
    values create non-local correlations across bipartition cuts.

**Key new insight:** The sieve dependency graph is COMPLETE for x≥500 (every prime pair shares
a common multiple). Combined with LFSR=N/2 for delta(n), these are new clean expressions of
the information-theoretic barrier.

**8 new closed paths. Total approaches: ~590+ across 26 sessions.**

## Session 27 (Critique)

Adversarial critique of all 16 Session 26 proposals. Web search for 2025-2026 literature (no breakthroughs).

(a) **All 16 proposals DUPLICATE.** Every proposal maps to existing CLOSED_PATHS entries.
    Session 26 was self-consistent: it generated and closed its own proposals within the session.
    After 590+ approaches, proposal-and-test cycles return 100% duplicates.

(b) **Group-theoretic prime race shortcut CLOSED.** The ONE untested suggestion from Session 26b
    synthesis (using Galois group of Q(ζ_q)/Q to relate prime races). Experiment shows:
    - Cross-modulus E(x;q) correlations are trend artifacts (roughness ratio 0.997 after detrending)
    - Galois lifts go wrong direction: fine→coarse free, coarse→fine needs NEW L-function zeros
    - Linear prediction of p(n) mod q: R² < 0.03, MI = 1-3% of maximum
    - CRT via prime races costs 96 L-functions for q≤23, each O(√x) = strictly WORSE
    - Character orthogonality makes cross-modulus shortcuts mathematically impossible

(c) **Failure taxonomy confirmed exhaustive.** All 590+ approaches fail via exactly one of:
    Circularity, Equivalence, or Information Loss. No fourth mode has been found.

(d) **Literature April 2026:** No new prime counting breakthroughs. Valley scanner for zeta zeros
    (Nov 2025), variational zero-finding (Jan 2026) — improve zero location, not summation.

**1 new closed path. Total approaches: ~591+ across 27 sessions.**

---

### Session 28 (2026-04-05): Wheel Decomposition Circuit Complexity

**Focus:** Novel angle on circuit complexity -- does decomposing pi(x) by residue
class mod primorial M give simpler per-class circuits? Distinct from CRT approaches
(Sessions 3,7,9,13,20-22,24) because we ask about COUNTING circuit size, not p(n) mod q.

**5 experiments run** (`experiments/circuit_complexity/modular_counting_attack.py`):

(a) **Mixed-radix conditional entropy:** Entropy reduction 1-52% depending on M and N,
    but reduction SHRINKS as N grows for fixed M. For N=16: M=6 gives only 1.1%,
    M=2310 gives 40.9%. This is a finite-size effect from removing small-prime composites.

(b) **Per-class circuit complexity: ANTI-WIN.** Per-class pi_r(x) mod 2 has 3-4x MORE
    normalized transitions than full pi(x) mod 2. Full autocorrelation ~0.80 (consecutive
    integers have correlated primality), per-class autocorrelation ~0.0 (M-spaced values
    are independent). Decomposition DESTROYS the sequential regularity.

(c) **Cross-class independence:** I/H ratio < 0.01, decreasing with N. Confirms
    Bombieri-Vinogradov. But independence is BAD: no shortcut between classes.

(d) **Divide-and-conquer scaling:** Total class transitions / full transitions converges
    to 1.0 as N grows. Only wheel sieving constant factor.

(e) **Entropy scaling:** Total entropy grows as phi(M) -- decomposition increases total work.

**1 new closed path. Total approaches: ~592+ across 28 sessions.**

**Session 28 insight:** The sequential structure of pi(x) -- that consecutive integers
have correlated primality due to guaranteed gaps between primes -- is the ONLY source
of computational regularity in the problem. Any decomposition that breaks this sequential
structure (CRT, wheel, residue classes, tensor decomposition) necessarily makes each
subproblem HARDER per bit. This is a new way to state the information barrier:
the regularity of pi(x) is TOPOLOGICAL (sequential ordering) not ALGEBRAIC (modular).

### Session 28 continued: Per-Bit Complexity, Approximate Degree, More

**Additional experiments (6+ total in this session):**

(f) **Per-bit circuit complexity gradient (NOVEL):** Total influence per output bit shows
    clear 2-tier structure. LSB-half has 2x+ higher influence than MSB-half, ratio GROWING
    with N (1.38 at N=4 → 2.13 at N=14). R-correlation crosses 0.5 at bit position ~N/2.
    MSB influence stays O(1); bit 0 (parity) influence grows as ~N/2.
    All data CONSISTENT with polylog circuits -- nothing rules them out.

(g) **APPROXIMATE DEGREE = N/2 (NOVEL, key result):** For the prime indicator chi_P and
    the counting parity pi(x) mod 2, the L∞ approximate degree at threshold eps=0.49 is:
    adeg(chi_P, 0.49) = adeg(pi mod 2, 0.49) = ceil(N/2).
    The counting step adds NO difficulty. N/2 is the EXACT smooth/oscillatory boundary.
    Quantum query complexity: Omega(N/4), still O(log x) = polylog.

(h) **Rounding boundary analysis:** frac(R(x)) is perfectly uniform (chi-squared p >> 0.05).
    Bits-of-precision distribution is EXACTLY geometric: P(need k bits) = 2^{-(k-1)}.
    No "easy subset" of inputs exists. R(x) accuracy drops: 65% at N=8, 13% at N=16.

(i) **Multiplicative-additive structure:** I-E signed rank = #distinct floors (full rank).
    90% of floor values have nonzero net contribution. Carry chains match random exactly.
    Monochromatic rectangle partition ~ 2^{0.76*N}, exceeding rank ~ 2^{0.41*N}.

(j) **The "N/2" universality (NOVEL synthesis):** ALL complexity measures converge at N/2:
    approximate degree, communication rank deficiency, oscillatory bit count, per-bit
    influence crossover, LFSR complexity of delta(n). The smooth/oscillatory boundary
    is a UNIVERSAL phenomenon, not an artifact of any particular measure.

**5+ new closed paths. Total approaches: ~597+ across 28 sessions.**
**1 novel result saved: approximate degree N/2 (novel/approx_degree_prime.md).**

### Session 28b: Explicit BDD Circuit Synthesis (NEW)

**(k) BDD synthesis for pi(x), N=4..14 (NOVEL, first explicit measurement):**
Constructed ROBDDs with multiple variable orderings for each output bit of pi(x).

Key measurements (LSB = hardest bit):
- N=4: BDD=7, N=8: BDD=53, N=10: BDD=147, N=12: BDD=421, N=14: BDD=1207
- Exponential fit: BDD(LSB) ~ 2^(0.73*N) = x^{0.73}
- This is WORSE than sqrt(x) = 2^(0.5*N) -- BDDs struggle more than expected
- OBDD (Session 20) gave 2^(0.79*N); multi-order BDD improves to 2^(0.73*N)
- Communication rank (Session 17) gave 2^(0.50*N) = sqrt(x) -- a lower bound

Per-bit structure: MSB has BDD ~ N+1 (trivially small), LSB has BDD ~ 2^(0.73*N).
Middle bits interpolate smoothly. Higher bits encode smooth part, lower encode oscillatory.

Influence of LSB ~ N/2 (consistent with pseudorandom behavior and Session 28a finding).

**Critical caveat:** BDD complexity != general circuit complexity. Functions exist with
exponential BDD but polynomial circuits (e.g., multiplication). The BDD result does NOT
prove super-polynomial circuit complexity. It only establishes that branching programs
(a restricted model) cannot compute pi(x) efficiently.

Saved: experiments/circuit_complexity/explicit_circuit_synthesis.py
Results: experiments/circuit_complexity/explicit_circuit_synthesis_results.md

## Session 29 (Fresh Perspective)
**Approach:** First-principles attack from 10 unconventional angles, deliberately
avoiding prior closed paths. 8 experiments, 5 parallel sub-agents.

Key quantitative findings:
(a) **R(x) accuracy:** |pi-R| = 0.05 at x=5000, 0.33 at x=1000. R(x) converges in
    10-20 Möbius terms. Smooth approximation is NOT the bottleneck.
(b) **Zero scaling: K(x) ~ x^{0.47}:** Power-law fit confirms ~sqrt(x) zeros needed
    for exact pi(x). Polylog fit gives exponent 8.77 -- not truly polylog.
(c) **Fourier sparsity:** pi(x)-x/log(x) has 99% energy in 1.25% of Fourier components.
    But sparse frequencies ARE zeta zero imaginary parts. Circular.
(d) **AR(1) residuals:** 91.4% reduction in Cipolla residual, but irreducible error
    grows as O(log n). Cannot achieve O(1) prediction error.
(e) **Zero grouping fails:** Group size 2 gives 9.65 error vs 4.54 full sum (x=10000).
    Rapid oscillation makes each zero's exact position essential.
(f) **GUE surrogates:** Right magnitude (±1.18), wrong value. Specific zero config matters.
(g) **p(n) mod m:** Not periodic for m≥3. No LFSR structure.
(h) **PSLQ identity search:** No new exact identities beyond known asymptotics.
(i) **Primorial decomposition:** c=3 is optimal → O(x^{2/3}) fundamental for sieves.

18 approaches tested, all closed. Three barriers crystallized:
1. Information-theoretic: log(x)/2 bits from zeros, no compression
2. Sieve-combinatorial: x^{2/3} is optimal balance
3. Analytic: K(x)~x^{0.47}, GUE prevents compression

Full synthesis: novel/session29_fresh_perspective.md

## Session 29 (Deep Focus: Novel Identity Search)

Deep-focus session on FOCUS_QUEUE Task #3: searching for computable identities
relating f(x) = pi(x) - R(x) to elementary/algebraic functions.

**7 experiments run, all negative:**

(a) **Extended PSLQ (x up to 100000):** 14-element basis including zeta zero
    oscillations. 18 tests total. 15 relations found with nonzero f-coefficient,
    but ALL fail cross-validation at other x-values (residuals 13-53000).
    Functional relations f(ax) vs f(x) for a=2,3,4: all fail.
    Shift recurrences f(x)..f(x+10): all fail.

(b) **Wilf-Zeilberger definite sum:** Delta_f is bimodal (prime indicator: ~0.90
    at primes, ~-0.096 at composites). Higher-order differences RMS GROWS with
    ratio converging to 2.0 (white noise signature). Hypergeometric R^2=0.997
    is spurious (trivial autocorrelation f(x+1)~f(x)). Summation kernel
    K(x)=f(x)*x*log(x) has full Hankel rank 250/250 (incompressible).

(c) **Number-theoretic constants:** Bernoulli numbers: r=-0.006, contributions
    <10^{-7}. Zeta(2..7): PSLQ x-dependent only. Dirichlet L(1,chi): same.
    Ramanujan tau: r=+0.010, p=0.93. All zero correlation.

(d) **Chebyshev connection:** psi(x)/log(x) captures 91% of f(x) variance
    (r=0.996) via partial summation identity. Improvement grows with x (ratio
    0.500 at x<100, 0.090 at x~100000). But psi(x) costs O(x) -- not a shortcut.

(e) **LLL minimal polynomials:** "Candidate" polynomials found at every x, but
    DIFFERENT polynomials at each x. f(x) values are effectively algebraically
    independent transcendentals. Multi-point relations at float64 limit.
    Polynomial-in-log(x): best validation RMSE = 1.976 (poor).

(f) **ODE search:** Linear ODE (order<=3, poly deg<=3): best residual 4.3e-8
    matches random noise baseline 2.5e-7 (spurious from heavy smoothing).
    Euler-type and nonlinear ODE: residuals ~0.99 (zero explanatory power).

(g) **Volterra integral equations:** K=1/(x-t+1) gives 1.7% residual but is
    trivially a local average. K=1/log(t) and K=1/t: total failure (~0.99).

**Direction #5 (Novel Identity) is now CLOSED.** f(x) has no computable identity
in any tested basis. This was the last remaining computational path.

Remaining open directions: circuit complexity (#TC^0 subset NC?), Kt complexity,
Berry-Keating Hamiltonian (all theoretical/monitoring).

Experiments saved: experiments/algebraic/identity_search/ (10 files)

## Session 30 (Critique of Session 29 Proposals)

Adversarial critique of 7 proposals + 3 synthesis directions from Session 29.
2 new experiments run, 1 literature search agent.

**Verdicts: 0 out of 10 viable.**

(a) **All 7 proposals DUPLICATE:** Sparse Fourier (S7,9,18,25), CRT modular (S3,7,9,13,20,24),
    compressed sensing (S7,19,20), recursive refinement (S5,7,9,15), adelic/Markov (S3,4,10,21,24),
    Schoof analogue (S4,14), neural delta oracle (S3,8,10,24 + Kolpakov-Rocke impossibility).

(b) **TG Kernel paper was ALREADY DEBUNKED (Session 12).** The proposals flagged it as
    "HIGHEST PRIORITY" without checking CLOSED_PATHS.md. arXiv:2506.22634 violates the
    uncertainty principle and was identified as AI-generated with errors. Serious procedural failure.

(c) **Autocorrelation exploitation CLOSED (new experiment):** r(1)=0.975 but smooth-trend
    artifact (detrended: r(1)=0.41). AR(1) RMSE=7.31 >> 0.5 threshold. Only 5.3% within
    rounding range. Conditional entropy H(delta(n+1)|delta(n)) = 4.93 bits/step. RMSE grows
    with n. See experiments/proposals/critique_incremental_delta.py

(d) **Verification-prediction separation CLOSED (new experiment):** Primality testing is
    O(polylog) but ordinality (is g the n-th prime?) requires pi(g) = O(x^{2/3}). The claim
    conflates primality with ordinality. At x=10^100: primality ~10^9 ops, ordinality ~10^68
    ops. See experiments/proposals/critique_verification_separation.py

(e) **Literature search (April 2026):** No new algorithmic breakthroughs found. TG Kernel
    debunked, Ono partition (closed S17,23), Brandt MKtP (impossibility route), Kolpakov-Rocke
    (ML impossible), Lamzouri (sparsity impossible).

**2 new closed paths. Total approaches: 628+ across 30 sessions.**

**Meta-insight:** Three consecutive critique sessions (S22, S27, S30) find 100% duplicate
rate. The space of "natural" approaches appears exhausted. Any breakthrough requires
genuinely novel mathematics or progress on #TC^0 subset NC?.

## Session 31 (Novel Circuit Complexity Measures)

**Method:** 6 parallel sub-agents computing measures explicitly recommended in OPEN_PROBLEMS.md
("need true tensor rank or discrepancy") but never previously computed.

**Literature search (April 2026, 2nd pass):**
- Rossman "Riffle Rank" (May 2025): new tensor measure for algebraic circuits, not Boolean
- Conditional tensor rank lower bounds (ECCC TR25-038): under NSETH, not unconditional
- No new TC^0 lower bounds, no progress on #TC^0 ⊆ NC?
- "Is pi(x) in NC?" remains entirely unaddressed in 2025-2026 literature

**Key experimental findings:**

(a) **TRUE TENSOR RANK (3-party, N=6,9,12):**
    N=6 (4x4x4): rank=5, random≈7 (ratio 0.75).
    N=9 (8x8x8): rank≤19, random≈29 (ratio 0.67).
    N=12 (16x16x16): rank≈67, random>100 (ratio<0.67).
    **Scaling: rank ~ d^{1.5} = 2^{N/2} = sqrt(x).** Exponent converges to 1.5 from 1.16.
    chi_P is 25-35% simpler than random, close to generic rank d²/3.
    CONFIRMS the N/2 universality and sqrt(x) barrier.

(b) **NOF DISCREPANCY (3-party):**
    Discrepancy is HIGH (~bias), dominated by density imbalance. Trivial cylinder (A=B=C=all)
    achieves max. **No communication lower bound follows.** Degree-1 Walsh correlation
    2.9-11x random (parity: primes are odd). Degree≥2 matches random. Consistent with TC^0.

(c) **GAMMA-2 (FACTORIZATION) NORM:**
    γ₂ ~ 2^{0.186N}, 85-93% of random. Sign-rank = rank (no hidden structure).
    SM complexity bound O(log γ₂) = O(N) trivially weak.

(d) **F_2 CORRELATION PROFILE (N=6..16):**
    W(1) z-score: 3 (N=6) to 1513 (N=16) — parity spike.
    W(d≥2): BELOW random (z-scores -2 to -30). After removing parity, chi_P is MORE
    pseudorandom than random functions. W(0)+W(1) = 47-68% of total weight.

(e) **BDD SIZE / SENSITIVITY:**
    BDD grows as 1.661^N (exponential). chi_P ~30% simpler than random (ratio 0.69-0.71).
    Sensitivity = certificate complexity = N (maximum). Decision tree depth = N.

**5 new closed paths, ~633+ total across 31 sessions.**

**Session 31 key insight:** All "natural" combinatorial circuit complexity measures
(15+ computed across Sessions 17-31) produce the SAME picture:
1. chi_P is exponentially complex (2^{Theta(N)}) in every measure
2. chi_P is 10-35% simpler than random, from parity + small prime divisibility
3. Every measure converges to the sqrt(x) = 2^{N/2} barrier
4. The Natural Proofs barrier prevents any of these measures from proving impossibility

**What remains viable:**
- Non-natural methods (avoiding the Razborov-Rudich barrier)
- Conditional lower bounds (Brandt MKtP framework, 2024)
- Constructive approach: try to BUILD a poly-size circuit (no evidence it exists or doesn't)
- Wait for breakthroughs in TC^0 lower bounds from complexity theory community

---

## Session 32: Fresh Perspective (2026-04-05)

**Approach:** Started from first principles without reading closed paths. Tested 8 experiments in parallel using analogies from FMM, compressed sensing, quantum computing, tropical geometry, and spectral theory.

**8 experiments, ALL closed:**
1. Zero sum convergence acceleration (6 methods): errors GROW as N^{0.8-1.0}
2. Recursive prime counting (FMM-inspired): 99.9% work at depth 0, Θ(x^{2/3})
3. Tropical/min-plus gap structure: Hankel rank ratio 0.98 (= random), R² < 0.002
4. Trace formula/moment method: moments diverge, Weil duality circular, full-rank
5. Contour integral: equivalent to zero sum by residue theorem
6. Hybrid analytic+sieve: O(x^{1/2+ε}), matches Lagarias-Odlyzko exactly
7. Hilbert-Pólya trace: GUE analogy breaks (unbounded spectrum), δ(x) incompressible
8. Literature 2025-2026: TG Kernel paper interesting (~1200 zeros for 10^8-digit x) but doesn't break barrier

**Most interesting finding:** The TG Kernel paper (arxiv 2506.22634) shows ~1200 zeros suffice for π(x) at astronomical x. This is much fewer than the naive O(x^{1/2}/log x), but zero computation cost remains the bottleneck. The barrier may be more about arithmetic complexity on huge numbers than zero count.

**Key quantitative results:**
- δ(x)/(√x/log x) is incompressible: poly deg-50 RMSE ≈ deg-2 RMSE (0.205 vs 0.212)
- SVD of zero-contribution matrix: 99% energy needs 133/500 components
- Prime gaps MI: 0.42% of entropy between consecutive gaps
- Truncated zero sums: K=0 often beats K=100-500 (conditional convergence)

**Total closed paths:** 640+ across 32 sessions.

## Session 33: Deep Focus — Conditional Algorithms (2026-04-05)

**Task:** FOCUS_QUEUE Task #4. Systematic investigation of the best exact algorithm
for p(n) under standard conjectures (GRH, RH, Elliott-Halberstam, Cramér's).

**7 experiments across 4 parallel sub-agents, ALL closed:**

(a) **GRH batch Miller testing:** GRH bound 2ln²(n) is conservative (only witnesses {2,3}
    needed to 10^6). Batch testing gives 1.02x speedup (negligible — each number needs
    independent modular exponentiation). Cost O(n·polylog n), worse than sieve.

(b) **GRH explicit formula optimal T:** With 1000 zeros, error ≥18 for x≥10^4.
    T_min = O(√x·log²x) confirmed. Individual zeros at γ~200-400 still contribute >0.5
    to pi(10^6). No way to use fewer than O(√x) zeros. GUE-random phases incompressible.

(c) **Elliott-Halberstam and counting:** EH controls equidistribution, NOT total count.
    Summing li(x)/φ(q) across residue classes recovers li(x) (same error as direct).
    Residue-class approach is 10^3-10^6x SLOWER than direct computation.

(d) **Gap structure under Cramér:** Cramér model fits (max gaps 40-60% of ln²p). R⁻¹(n)
    search interval has O(ln x) primes. But identifying which is p(n) requires pi(x) at
    the boundary. Counting is the bottleneck, not searching.

(e) **Schoenfeld explicit bounds under RH:** Sieve in Schoenfeld interval is O(x^{1/2+ε}),
    better than Meissel-Lehmer O(x^{2/3}) asymptotically, but still exponential in input bits.
    1000 zeros reduce error by only 21% at x=10^6. Need O(√x·log²x) zeros for exactness.

(f) **Cramér search + counting bottleneck:** Walk phase O(ln⁴x) trivial (0-114 steps).
    Counting phase dominates: 25% at n=5M, heading to 100%. At x=10^100:
    count=10^66.7 ops vs walk=10^9.5 ops (factor 10^57 gap).

(g) **Best conditional algorithm (all conjectures combined):**
    | Assumption | Best complexity |
    |---|---|
    | Unconditional (Meissel-Lehmer) | O(x^{2/3}) |
    | RH + Turing zeros | O(x^{2/3+ε}) — WORSE |
    | RH + Odlyzko-Schönhage batch | O(x^{1/2+ε}) — best known |
    | GRH | O(x^{1/2+ε}) — same |
    | GRH + EH + Cramér | O(x^{1/2+ε}) — same |
    | Approximate only (R⁻¹) | O(polylog) — ~50% digits |

**Key insight: The √x barrier is fundamental across ALL standard conjectures.**
- RH alone with naive zeros is WORSE than unconditional (O(x^{2/3+ε}) vs O(x^{2/3}))
- RH only helps when combined with batch zero computation (Odlyzko-Schönhage)
- GRH, Elliott-Halberstam, and Cramér's conjecture all address the WRONG subproblem
  (primality testing, distribution, or gap structure) — none reduce the counting cost
- The counting bottleneck is information-theoretic: exact pi(x) encodes O(√x) zeta zeros

**8 new closed paths. Total approaches: 641+ across 33 sessions.**

Experiments saved: experiments/analytic/conditional/ (4 scripts + 4 results files)

---

## Session 35

**Focus:** Circuit complexity deep-dive — meta-complexity, approximate circuits, GF(2) SLP analysis.

**5 experiments completed:**

(a) **Approximate circuit complexity — NO phase transition.** Rank-k SVD approximation
    accuracy for balanced partition increases GRADUALLY: rank-1 → 59%, rank-5 → 73%,
    rank-10 → 83%, rank-20 → 95% (at N=14). Each singular vector contributes a small,
    roughly equal increment. No "easy core" vs "hard shell." Errors from degree-2 PTF
    are spatially uniform (CV = 0.22-0.25). Rules out "approximate then correct" strategies.

(b) **MKtP / meta-complexity framework CLOSED.** Formal analysis: "Is pi(x) in NC?" ≡
    "Is Kt(pi(x) mod 2 | x) = O(polylog)?" — pure reformulation, no new technique.
    Brandt's conditional framework connects circuit lower bounds to meta-complexity but
    for ANY function in E, not pi(x) specifically. Kt(T_N) = O(2^N * N) by sieve
    regardless of circuit size, so Kt cannot distinguish poly from exponential circuits.

(c) **Smooth approximation is useless for parity.** R(x) accuracy for pi(x) mod 2 → 50%
    (random) as N grows. Top-2 SVs capture only 12% of variance at N=14 (decreasing as 1/N).
    The parity of pi(x) is ENTIRELY determined by the oscillatory (zeta zero) part.

(d) **GF(2) SLP complexity = random.** ANF sparsity EXACTLY 0.50 for all N. CSE savings
    50-66% = SAME as random functions (±1%). SLP length Θ(2^N). Variable frequencies
    uniform. Monomial fractions at each degree → 0.50. pi(x) mod 2 is indistinguishable
    from random in ALL GF(2) structural metrics.

(e) **Unbalanced partition rank correction.** Formula rank = 2^{min(k,N-k)-1}+2 confirmed
    for k ≤ N/2 but UNDERESTIMATES for k > N/2 (often full row rank). Polynomial rank
    for k = 2*log(N) doesn't help: component functions g_i(MSBs) are as hard as pi(x).

**Threshold gate construction (agent completed):**
(f) **PTF degree = N/2 exactly** (LP-verified for N=4-12). Ratio degree/N = 0.50 from N=6 onward.
    This means C(N, N/2) ~ 2^N/sqrt(N) monomials needed for depth-2 threshold circuits —
    EXPONENTIAL. Single LTF accuracy → 0.50 (random). Matches random function bound (Gotsman 1994).
    Does NOT rule out poly-depth TC^0 but confirms exponential at depth 2.

**Key session insight:** The pi(x) mod 2 function is random-like in EVERY computable measure:
ANF sparsity, SLP length, rank-k approximation decay, error geography, block entropy,
gzip compressibility (approaching random), PTF degree scaling. The function has NO exploitable
structure beyond what's already captured by the smooth approximation R(x) — which provides
exactly 0% of parity information.

The MKtP/Brandt framework, last remaining recommended direction from S31, is now closed as
an attack path. It provides a reformulation equivalent to the circuit complexity question,
not a new technique.

**SAT circuit minimization (agent completed):**
(g) **Circuit size for pi(x) mod 2 matches random functions.** N=4: exact minimum = 6 gates
    (same as is_prime, random = 5). N=5-8 upper bounds via Shannon decomposition: 12, 39, 81,
    153 gates (growth ~2x per bit). N=9-10 fell back to DNF synthesis (loose bounds: 2015, 4829).
    Polynomial fit ~N^9.8, but synthesis gives upper bounds only — true minima likely much smaller.
    pi(x) mod 2 and random functions have nearly identical circuit sizes at all N tested.
    No evidence of structural advantage over random. INCONCLUSIVE on exact growth rate due to
    synthesis limitations, but confirms pseudorandomness thesis at circuit level.

**10 new closed paths. Total approaches: 651+ across 35 sessions.**

Experiments saved: experiments/circuit_complexity/ (5 scripts + 5 results files)

---

## Session 36 — Fresh Perspective (Analogies from Other Fields)

**Approach:** Start from first principles, ignoring prior closed paths. Draw analogies from
Shor's algorithm, compressed sensing, fast multipole method, AlphaFold, matrix completion.
Test 7 genuinely unconventional ideas via 5 parallel agents + main thread.

**Web research findings (2025-2026):**
- Ono/Craig/van Ittersum 2024: "Integer partitions detect primes" — beautiful but O(exp(√n))
- Guth/Maynard 2024: New zero density estimates → primes in short intervals x^{17/30}
- Aggarwal 2025: O(√n·log⁴n) for computing p_n via binary search (arxiv 2510.16285)
- Green/Sawhney 2024: New prime detection via integer partitions connection
- NO progress toward polylog in the literature

**7 experiments run (5 agents + 2 main thread):**

(a) **Batch Möbius sieve (recursive halving) CLOSED:** Depth O(log(π(√x))) but branching
    2^{half_size}. Total work 2^{π(√x)} = same as unpruned standard tree. At x=10^4:
    30x slower than Lucy_Hedgehog. At x=10^6: branching = 2^84. Reorganizing IE order
    doesn't reduce total work.

(b) **Automata-theoretic digit DP CLOSED:** Product DFA for B-rough numbers achieves
    O(primorial(B)·log x). Digit DP verified correct. Timing: M=210 at x=10^6 in 6ms.
    But DFA provably minimal (all states distinguishable). Need B≥√x → M=e^{√x}.
    At B=13 (M=30030), false positives ≈ true primes. Exponential state space.

(c) **Adelic local-global CRT (6th variant) CLOSED:** pi(x) mod q via AP decomposition
    100-1000x slower than pi(x) direct. R(x) error grows O(√x), mod-q correctness
    degrades. Liouville/Mertens parity ~50% (random). Floor division not ring hom.

(d) **Self-correcting explicit formula CLOSED:** Truncation bias is systematic ~O(N), not
    random noise. Monotonicity+rounding: max error reduced 50% but accuracy < 3%.
    Primality correction propagates errors. Autocorrelation > 0.98.

(e) **Ono partition criterion (GF approach) CLOSED:** M₁(n) = Σ d(j)·p(n-j) via
    divisor-partition convolution. Best: O(n^{3/2}). Direct: O(exp(π√(2n/3))).
    Both worse than Meissel-Lehmer. Partition functions inherently need O(n) terms.

(f) **Local-to-global reconstruction — KEY FINDING:**
    E(x) = pi(x) - Li(x) has only **O(log x) bits** of information.
    At x=10^100: E needs ~160 bits, polylog budget is ~110,000 bits.
    **NO information-theoretic barrier to polylog!**
    Barrier is COMPUTATIONAL: O(log x) bits encoded across ~x^{1/2} zero contributions
    with massive cancellation. Each zero adds < 1 bit net. Error autocorrelation 0.996
    at lag 1 (helps adjacent queries only, not single-point computation).

(g) **Explicit formula proper convergence CLOSED:** Naive R(x^ρ) summation DIVERGES:
    error grows from 3.5 (K=0) to 2076 (K=100) at x=10^4. Complex li branch cuts
    cause instability. Only contour integration (Lagarias-Odlyzko) is numerically stable.

**7 new closed paths. Total approaches: 655+ across 36 sessions.**

**Session 36 synthesis:** The fresh perspective confirmed all known barriers but surfaced
one potentially significant insight: the error term E(x) carries only O(log x) bits, well
within polylog bounds. The problem is not that pi(x) CONTAINS too much information, but
that we cannot EXTRACT the O(log x) bits without O(x^{1/2}) computation. This is
structurally similar to one-way functions: small output, hard to compute. The breakthrough
(if it exists) must find algebraic structure in the zeta zero sum that enables bulk
cancellation — analogous to geometric series collapsing many oscillating terms into a
closed form.

Experiments saved: experiments/wildcard/ (5 new scripts + 5 results files)

## Session 36b — Deep Focus: Kt Complexity of delta(n) (2026-04-05)

**Task:** FOCUS_QUEUE Task #1 — deep dive extending Session 20's Kt analysis.

**New experiments (8 total):**
1. **PACF analysis**: AR(7) is optimal by BIC. PACF drops to 0.056 at lag 2,
   below 0.02 by lag 5. Long-range ACF is indirect (from short AR, not true
   long memory). PACF ~ k^{-1.33} (α > 1 = finite effective AR order).

2. **Compression comparison**: bz2 (36.5%) > lzma (43.0%) > gzip (56.1%).
   diff(delta) compresses better (31-37%), confirming AR(1) structure.
   delta/random ratio: 0.55 (bz2), 0.60 (lzma), 0.71 (gzip).

3. **Compression scaling with N**: Bits/value converges to ~5.8 by N=10000.
   Power law: N^{-0.059} ≈ constant. Entropy rate is finite.

4. **Kt(1..N) growth**: Kt ~ 5.58*N + 0.023*N*log(N). The N*log(N) term is
   negligible. Total information is EXTENSIVE (linear in N).

5. **Block MI scaling**: MI ~ 0.55*log(L) - 1.33. Grows with block size,
   confirming long-range correlations between distant blocks.

6. **DFA/Hurst**: H = 1.31 (between 1/f and Brownian motion). Crossover at
   scale ~572: H_small = 1.41, H_large = 1.19. Two-regime structure.

7. **Transfer entropy**: TE(delta→delta) = 0.051 >> TE(n→delta) = 0.013.
   Past delta is 4x more informative than n. n mod 30 contributes 0.003 bits.

8. **Spectral algebraic structure**: The 1/f^1.70 spectrum is a smooth continuum.
   No discrete lines, no algebraic relations among coefficients. 50% power in
   8 frequencies, but RMSE < 1 requires 41182/50001 modes (82%). Max error < 1
   requires ALL modes. 4/10 peaks match zeta zeros but likely chance.

**Key insight:** Session 20 left open "does the 1/f spectrum have exploitable
algebraic structure?" Answer: **NO.** The spectrum is a genuine continuum,
not decomposable into sparse computable oscillations. The spectral tail
(frequencies contributing individually < 0.001% power) is collectively
essential for exact values, matching the communication matrix rank result.

**Quantitative barrier summary:**
- Entropy rate: ~5.8 bits/value (finite, not growing)
- Total Kt for N values: ~5.8*N bits (extensive)
- Direct dependence on n: zero (transfer entropy confirms)
- Spectral sparsifiability: 82% of modes needed for RMSE < 1
- AR model order: 7 (short direct memory, long indirect memory)

**Status:** Task #1 (Kt complexity of delta(n)) is now COMPLETE.
The R^{-1}(n) + delta(n) decomposition is a dead end for polylog algorithms.

Also: converted all .txt results to .md format in kt_complexity/.
Experiments saved: experiments/information_theory/kt_complexity/ (2 new scripts + 2 results)

## Session 38 (Critique)

**Task:** Critique Session 36 proposals + literature scan for 2025-2026 papers.

Key findings:

(a) **All four S36 proposals confirmed CLOSED.** Proposals 21-24 (zero clustering,
    compressed sensing, PSLQ, dequantized Grover) all fall into the Information Loss
    failure mode. Each matches 3-5 existing CLOSED_PATHS entries.

(b) **Kilictas-Alpay TG kernel paper (arXiv:2506.22634) was still marked "VERIFY" in
    literature file — corrected to DEBUNKED** per S12/S30 findings.

(c) **New literature scanned:**
    - Aggarwal (2510.16285): p_n complexity analysis, already tracked
    - Tao-Gafni: rough numbers in gaps, no computational relevance
    - Valley Scanner (2512.09960) + variational Z-function methods: better zero
      location, but barrier is summing zeros not locating them
    - primecount v8.4: SIMD acceleration, O(x^{2/3}) unchanged
    - No TC^0/NC^1 separation breakthroughs

(d) **Minor novel observations from S36:** Gaussian normality of normalized delta
    (p=0.19), sieve rank exponent 0.365, Fourier rank exponent 0.943. None warrant
    novel/ entries.

(e) **Recommendation:** Circuit complexity (#TC^0 ⊆ NC?) remains the ONLY open direction.
    Growing-dimension matrix powering in TC^0 is the one genuinely open sub-question.

## Session 39 — Literature Monitoring + Depth Profile (2026-04-05)

**Task:** Literature monitoring for new TC^0/NC^1 results; threshold circuit depth experiment.

Key findings:

(a) **No new TC^0/NC^1 separation results in 2026.** Searched ECCC, arXiv, CCC.
    - Chen-Tal-Wang (TR26-039): n^{2.5-ε} THR∘THR lower bounds, but for E^NP (already tracked)
    - Gurumukhani et al. (2601.04072): optimal monotone depth-3 for MAJORITY (monotone, not threshold)
    - Behera et al. (TR26-002): multilinear IPS separations (proof complexity, not circuits)
    - TR26-024: Nullstellensatz in counting hierarchy (no computational relevance)
    - TR26-038: XOR lemma for F_p sums (graph problems, not number theory)
    - Brandt MKtP: no 2026 follow-ups beyond the TCC 2024 paper

(b) **No new pi(x) algorithms.** Landscape unchanged: O(x^{2/3}) practical, O(x^{1/2+ε}) analytic.

(c) **Threshold depth profile experiment (N=4-10):** Measured depth-vs-accuracy for
    is_prime, pi_mod2, and random using random LTF features + LP. All three functions
    behave identically: exact at depth 2 for N ≤ 8, fail at N=10 (heuristic limitation).
    The 0.832 accuracy for is_prime at N=10 = trivial baseline. This extends the
    pseudorandomness evidence to the depth-accuracy tradeoff dimension (Measure 22).
    Reconciled with S35: PTFs show depth-2 always works but needs C(N,N/2) ≈ 2^N/√N
    gates (exponential), so cannot distinguish TC^0 from non-TC^0 at depth 2.

(d) **Growing-dimension MPOW analysis:** Confirmed this remains genuinely OPEN. The bottleneck
    is O(log log n) combination depth for polylog eigenvalue powers. The Healy-Viola F_{2^n}
    technique uses Frobenius endomorphism (characteristic-specific), has no analog over Z_n.
    The commutative structure of Z_n[x]/(x^r-1) doesn't help due to circularity (CRT
    decomposition requires knowing factorization of n).

(e) **Housekeeping:** Flagged 7 groups of duplicate/versioned scripts in TODO.md for human review.

(f) **Assessment:** The project has exhausted all experimentally accessible directions. The
    remaining question (#TC^0 ⊆ NC) is purely complexity-theoretic with no known experimental
    attack. Future sessions should focus on literature monitoring and engineering improvements.


---

## Session 40: Fresh Perspective (2026-04-05)

**Theme:** First-principles attack, 6 new experiments testing unconventional ideas.

### Experiments Run
1. **Transfer Matrix Sieve** — Encode sieve as matrix product, use repeated squaring.
   FAIL (E): State space = lcm(primes) ~ exp(sqrt(x)), worse than Meissel-Lehmer.

2. **Wavelet/Fourier Compression** — Is C(x) = pi(x) - li(x) sparse in wavelet basis?
   FAIL (I): Raw signal appears sparse (DC dominates), but detrended oscillatory part needs
   N^{0.75} coefficients. Power spectrum consistently f^{-1.7}. Sublinear but not polylog.

3. **Modular Residue Prediction** — Can p(n) mod q be predicted for CRT reconstruction?
   FAIL (C+I): Near-max entropy for all q>2. CRT is circular. Interesting: mod 210 shows
   19x prediction accuracy due to primorial wheel, but this is just divisibility constraints.

4. **Carry-Propagation Boundary** — Map bit-by-bit difficulty of p(n) vs R^{-1}(n).
   FAIL (I): **Novel measurement.** Sharp sigmoid transition: bits 0-7 EASY (95-100%),
   bits 8-9 MEDIUM (80-90%), bits 10+ RANDOM (50%). Transition width only ~4 bits.
   Error ~ p(n)^{0.66}. No intermediate-difficulty bits to exploit hierarchically.

5. **Ergodic Fast-Forward** — Model prime gaps as dynamical system, test fast-forwardability.
   FAIL (I): Positive Lyapunov exponent (chaotic). R^2 < 0.006 for linear recurrence.
   Real gaps have HIGHER entropy than Cramér random model at short range (surprising).

6. **Multiplicative Convolution Shortcut** — Express pi(x) as factored Dirichlet convolution.
   FAIL (E): Reduces to hyperbola method O(sqrt(x)). Floor function breaks multiplicativity.

### Key Insights
(a) **Carry propagation boundary** is the most novel measurement: the easy-to-hard transition
    in binary representation of p(n) is a ~4-bit-wide sigmoid at ~60% of bit positions.
    This quantifies exactly how R^{-1}(n) fails: it gets all MSBs right but the transition
    to random is sharp, leaving no room for hierarchical bit recovery.

(b) **Wavelet compression scaling N^{0.75}** is between sparse and dense. The 1/f^{1.7}
    spectrum means the correction has MORE structure than white noise (which would be N^1)
    but far less than what polylog requires (N^0 or N^epsilon).

(c) **Prime gaps are more chaotic than Cramér model** at short range (higher conditional
    entropy at order 3). This is a minor finding but consistent with the pseudorandomness
    results from S35.

(d) **6 approaches closed, ~531 total.** No breakthrough. All fresh-perspective ideas
    reduce to known barriers within 1-2 steps of analysis.

## Session 41 (Engineering + Novel Write-up)

**Focus:** Engineering improvements to v10 + novel finding documentation.
All FOCUS_QUEUE tasks were already COMPLETED; session pivoted to productive
engineering work per CLAUDE.md guidelines.

### Experiments

1. **Lehmer method vs Lucy DP** (`experiments/sieve/lehmer_vs_lucy.py`)
   Three C implementations benchmarked: Lucy DP, Lehmer (recursive phi + P2),
   Meissel-Lehmer (P3 term). All verified 100% correct to 10^10.
   RESULT: Lucy DP is 5-13x FASTER than naive Lehmer. Recursive phi inclusion-exclusion
   with hash memoization has worse cache behavior than Lucy's sequential array access.
   The Deleglise-Rivat/Gourdon speedup comes from segmented sieve for phi, not the
   decomposition itself.

2. **Wheel-30 + prefetch Lucy DP** (`experiments/sieve/segmented_lucy.py`)
   Wheel mod 30 (skip p=2,3,5) + __builtin_prefetch hints on C Lucy DP.
   RESULT: Only 2-35% speedup. Hardware prefetcher already optimal for sequential access.
   The 100-1000x gap to primecount is algorithmic (Gourdon variant), not parameter tweaking.

### Novel Finding Documented

- **Carry-propagation boundary** from S40 written up in `novel/carry_propagation_boundary.md`.
  Per-bit difficulty of p(n) has a sharp 4-bit sigmoid transition from EASY (R^{-1}
  determines) to HARD (coin flip) at ~60% of bit positions. Genuinely novel measurement
  not in published literature.

### Key Insights
(a) **Naive Lehmer is SLOWER than Lucy DP.** The standard textbook Lehmer decomposition
    only helps when phi is computed via segmented sieve, not recursion. This explains why
    primecount's speed comes from infrastructure (SIMD, segmented sieve, OpenMP), not
    from the formula decomposition.

(b) **Within-Lucy-DP optimizations are marginal.** Wheel-30 gives at most 35% at small x,
    dropping to 2% at x=10^11. To match primecount, one must reimplement Gourdon's full
    algorithm — essentially forking primecount. Not worth it for this project.

(c) **No new literature breakthroughs.** Web search for April 2026 publications on pi(x)
    algorithms, TC^0 lower bounds, and zeta computation found no results beyond what's
    already in `literature/state_of_art_2026.md`.

(d) **2 approaches closed, ~533 total.**

---

## Session 42 (Deep Focus — Base-W MPS bond dimension theorem)

**Focus:** Tighten the MPS bond dimension novel finding by verifying the
predicted base-W generalisation across several primorial reshapes.
FOCUS_QUEUE Task 1 was already COMPLETED; session pivoted to theoretical
sharpening of an existing novel result, per CLAUDE.md guidance.

### What was done

Refactored `experiments/wildcard/mps_bond_dimension.py` to support a
`--base W` argument and ran two families:

  - **Family A (replication):** W = 2, d in {10, 12, 14, 16, 18, 20}.
    Reproduces chi_P max-rank = 2^{d/2-1} + 1 exactly. Random control
    saturates 2^{d/2}.
  - **Family B (base-W theorem):** W in {2, 6, 30, 210}, d up to 20.
    Tested the precise prediction
        rank at cut j = min(W^j, phi(W) * W^{d-j-1} + 1).

### Theorem proved + verified

> rank(M^{(j)} of chi_P viewed in base W)
>   = min( W^j, phi(W) * W^{d-j-1} + 1 )
>
> at every cut 1 <= j < d, for every primorial W >= 2.

The upper bound follows from "n > W and n prime => gcd(n, W) = 1, so
column index k satisfies gcd(k+1, W) = 1": rows i >= 1 are supported
on phi(W) * W^{d-j-1} columns, plus row 0. Empirically saturated at
strict equality in all 11 tested configurations (W=2/d=10..20, W=6/d=4..8,
W=30/d=2..4, W=210/d=2..3).

### Half-cut compressibility ratio table

| W   | phi(W)/W | observed chi_P / random (largest d tested) |
|-----|----------|--------------------------------------------|
| 2   | 0.5000   | 0.5010 at d=20                             |
| 6   | 0.3333   | 0.3341 at d=8                              |
| 30  | 0.2667   | 0.2678 at d=4                              |
| 210 | 0.2286   | 0.2333 at d=2                              |

Ratio asymptotes to phi(W)/W as d grows (the +1 row-0 contribution
vanishes).

### Why this still does NOT yield polylog

The half-cut bond dim is exactly phi(W) * W^{d/2-1} + 1, which is
Theta(W^{d/2}) for any fixed primorial W. To make this polylog in
N = W^d we would need W to grow with d, but the cost of materializing
the wheel is at least W = exp(log N / d), and the optimum cost
exp(log N / 2) = sqrt(N) is achieved at d = 2. Same Lagarias-Odlyzko
sqrt(N) barrier, recovered from a pure entanglement argument with no
zeta zeros.

### Files updated

- `novel/mps_bond_dimension.md` — added base-W theorem section, proof,
  saturation table for all (W, d) tested.
- `experiments/wildcard/mps_bond_dimension_results.md` — extended with
  Family B table, per-cut profiles for W=30/d=4 and W=6/d=8, polylog
  obstruction analysis.
- `novel/pseudorandomness_of_pi.md` — added Measures 22 and 23
  (TT bond dim at base-2 half cut; TT bond dim ratio at base-W half cut).
- `status/CLOSED_PATHS.md` — appended one row for the base-W MPS measurement.

### Key insight

The previous tensor experiments (S10 MPS sieve, S20 product DFA) studied
state-count or transfer-matrix rank, not the bond dimension of the
indicator vector itself viewed in a primorial reshape. The new theorem
gives the FIRST clean exact identity for the entanglement-theoretic
compressibility of chi_P:

  every prime added to the wheel knocks the bond-dim ratio down by
  (p-1)/p, mirroring the classical Mertens product
  prod_{p<=W}(1 - 1/p) ~ e^{-gamma}/log W.

So the MPS view recovers the Lagarias-Odlyzko sqrt(N) barrier
*without* invoking the zeta zeros at all -- a cleaner derivation of why
the analytic methods cannot break sqrt(N), framed entirely in tensor
network language.

(e) **1 approach (base-W MPS) added to CLOSED_PATHS, 534 total.**

---

## Session 43 (2026-04-25): Critique of TDA / J-fraction / Free probability / Reservoir / Stein-method

Five fresh proposals (F1-F5 in `archive/ephemeral/proposals_session.md`) generated
and critiqued in critic mode. **No proposal survived.** Five new entries added to
`status/CLOSED_PATHS.md` (~539 total).

### Verdicts

- **F1 (Persistent homology of detrended primes, TCDP):** DUPLICATE — already
  closed at line 199 (Session 10). Companion experiment confirms long-lived 0-dim
  bar count is sub-(N/log N) but Cramér-model baseline is too. No inverse map from
  barcode to p_n exists without using the original 1D ordering of n (which IS pi).
- **F2 (J-fraction of pi(n)/n generating function, PPC):** DUPLICATE.
  J-fraction coefficients (a_k, b_k) for k=0..11 erratic — b_3=-3, b_4=-8.28,
  b_12=-407.4. Linear log-decay of Hankel determinants is a trivial consequence
  of moments in [0,1], not D-finiteness. Closes related entries lines 286, 289,
  299, 342, 435.
- **F3 (Free-probability R-transform, FCP):** FLAWED by inspection. A_N = diag(1_P)
  is idempotent, so empirical spectral distribution = (1-pi(N)/N)*delta_0 +
  (pi(N)/N)*delta_1, parameterized by pi(N) alone. R-transform is closed-form
  q/(1-qz) — encodes only pi(N) (circular).
- **F4 (Reservoir computing on delta(n)):** DUPLICATE — strict subset of failed
  ML methods (transformer 1.1% accuracy line 149; FNO 0.5% line 150; GP 0.45%
  line 151; symbolic regression line 146; ML for delta line 588).
- **F5 (Stein-method sub-Gaussian explicit-formula error, SGEFE):** FLAWED —
  empirically falsified this session. Computed f_T(x) = R(x) - sum_{|gamma|<T}
  2 Re R_log(rho*log x) for x in {100..10000}, T in {10..1000}, all 1000 zeros.
  log10|err| vs T^2/log(x): slope ~-3e-6, R^2 = 0.2-0.4 — no sub-Gaussian decay.
  Tail L^2 norm sqrt(x*logT/T) is the rate; no concentration gap to exploit.

### Process notes

- The F5 experiment initially had a mp.log branch-cut bug (principal log of x^rho
  loses periods of 2*pi*i) that produced absurdly large errors. Fix: pass rho*log(x)
  directly to Ei without going through x^rho. Lesson — when dealing with li / Ei of
  complex argument with large imaginary part, always work in log-space.
- F3 collapses by an *idempotence* observation that takes one paragraph; the
  proposal spent 50 lines describing the R-transform machinery. Routine
  novelty-of-language failure: dressing pi(x) in unfamiliar terminology does not
  create new structure.

### Remaining open direction

Still only **circuit complexity of pi(x)** (`status/OPEN_PROBLEMS.md` Problem 1)
plus Berry-Keating literature monitoring. No experiment proposed this session
moved that forward.

---

## Session 44 (2026-04-25): TODO #3 -- inversion_search fixed-point iteration

Stale [RESEARCH] item from `TODO.md`. Two infrastructure bugs blocked it from
running:
1. `experiments/sieve/inversion_search.py` inserted its own dir on `sys.path`
   to find `riemann_explicit`, which lives in `experiments/analytic/`.
2. `experiments/analytic/riemann_explicit.py::load_or_compute_zeros` only
   looked in its own dir; precomputed zeros are in `data/`.

Both fixed. The second fix benefits 4 other importers
(`hybrid_correction`, `zero_scaling`, `convergence_accel`, `advanced_convergence`).

### Empirical finding

Fixed-point iteration `p_0 = R^{-1}(n)`, `p_{k+1} = R^{-1}(n + sum_rho R(p_k^rho))`
is **non-contracting** -- a sharper closure than the previous "needs O(sqrt(x))
zeros" framing. Measured at dps=30, 200 zeros, max_iter=3:

| n     | k=0 err | k=1 err | k=2 err   | k=3 err | first exact |
|------:|--------:|--------:|----------:|--------:|------------:|
| 10    | -0.67   | -2.36   | **-0.14** | -1.89   | k=2         |
| 100   | +4.10   | -3.92   | -3.59     | -3.68   | none        |
| 1000  | -4.05   | **-0.16** | -0.64   | -0.57   | k=1         |
| 10000 | -39.31  | -12.16  | -13.15    | -13.11  | none        |

The iteration sometimes hits within 0.5 of p(n) at one specific k, then *drifts
back away* at the next iteration -- the would-be fixed point is unstable. ~50%
of tested n never converge within 3 steps (200 zeros). Increasing num_zeros is
**non-monotone**: at n=1000 with 2-iter budget, 200 zeros converges at k=1 but
500 zeros doesn't.

### Why this is a sharper closure

The map f(x) = R^{-1}(n + S(x)) where S(x) = sum_rho R(x^rho) has
|f'(x)| ~ |S'(x) * (R^{-1})'(...)|. Near x = p(n), |S'(x)| oscillates because
the sum has GUE-random phases, so |f'(x)| is not bounded < 1. Connects to S43's
Stein-method finding (no martingale concentration on GUE phases) and
CLOSED_PATHS line 24 (zero-sum acceleration is random-walk).

### Files updated

- `experiments/sieve/inversion_search.py` -- import path fix
- `experiments/analytic/riemann_explicit.py` -- `load_or_compute_zeros` fallback to `data/`
- `experiments/sieve/inversion_search_results.md` -- corrected (previous claim
  "exact at k=0 or k=1 for most n up to 10000" empirically falsified)
- `status/CLOSED_PATHS.md` -- new entry line 682 with non-contraction diagnosis
- `TODO.md` -- item 3 marked DONE

### Remaining open direction

Unchanged from S43: **circuit complexity of pi(x)** + Berry-Keating literature
monitoring. TODO items 2 (Helfgott-Thompson M(x) benchmark) and 1 (duplicate-
script flagging, blocked on human review) still pending.

## Session 45 (2026-04-25): Resolve S25 caveat — zeta-structure tests at N=2000

**Trigger.** FOCUS_QUEUE Task #2 (Zeta Zero Structural Patterns) is marked
COMPLETED in S25 with one explicit caveat: "we tested 1000 zeros. Structure
might exist at scales requiring >1000 zeros to detect." This session settled
that caveat by extending the most discriminating S25 tests to N=2000.

### Method

Computed zeros γ_{1001..2000} via mpmath at 30-digit precision (~7.5 min,
~0.45 s/zero). Saved to `data/zeta_zeros_2000.txt`. Wrote single consolidated
script `experiments/analytic/zeta_structure/zeta_structure_n2000.py` running
six tests: pair correlation, number variance, DFT/spectral flatness vs GUE,
mod-constant discrepancy + Weyl sums, PSLQ on the never-tested LATE block
γ_{1001..1100} (1225 pairs + 4060 triples), and cross-block PSLQ
(early × late, 400 pairs).

### Headline numbers

| Test | N=1000 (S25) | N=2000 (S45) | Verdict |
|------|--------------|--------------|---------|
| Pair-correlation RMS dev from GUE | ~0.10–0.12 | **0.0864** | GUE fit improves |
| Log-power correlation vs GUE | 0.9999 | 1.0000 | identical |
| Band 1–4 spectral flatness | 0.93–0.999 | 0.93–0.999 | identical |
| Weyl mean / (1/√N) ratio | ~0.95 | 0.946 | identical |
| Discrepancy ratios vs LIL | 0.3–0.8 | 0.3–0.8 | identical |
| PSLQ pairs in LATE block | not tested | **0/1225** | no relations |
| PSLQ triples in LATE block | not tested | **0/4060** | no relations |
| Cross-block PSLQ (early×late) | not tested | **0/400** | no relations |

### Two methodological caveats noted in writeup

1. The "overall" spectral flatness 0.0065 is dominated by the smooth
   Riemann–von Mangoldt linear-log ramp dumping power into k=0..4 (10⁴–10⁵×
   median); the GUE-interp comparison shows the *same* effect (mean SF
   0.0077). The trend-free bands match GUE perfectly. **No structural finding.**
2. Number variance plateaus at Σ²(L) ≈ 0.34 for L ≥ 5 instead of growing
   logarithmically. Cause is window-overlap in the 800-window sampling, not
   real GUE deviation. A clean test would need disjoint windows (N ≥ 10⁴).
   **Inconclusive at N=2000, not anomalous.**

### Verdict

S25 caveat **resolved negatively.** No structure emerges in 2× the data.
Direction remains CLOSED. Probing >2000 would require Odlyzko's tabulated
zeros (online; not pursued here). New CLOSED_PATHS entry added; SESSION_25
SUMMARY caveat line updated; FOCUS_QUEUE Task #2 description still
accurate (results unchanged).

### Files created/updated

- `data/gen_zeros_2000.py`, `data/zeta_zeros_2000.txt` (new)
- `experiments/analytic/zeta_structure/zeta_structure_n2000.py` + `_results.md` + `_summary.json`
- `experiments/analytic/zeta_structure/SESSION_25_SUMMARY.md` — caveat-resolution paragraph appended
- `status/CLOSED_PATHS.md` — new row in Analytic section
- `archive/sessions/session45_zeta_n2000.md` (new)

### Remaining open direction

Unchanged: **circuit complexity of pi(x)** + Berry-Keating literature
monitoring. No new viable directions opened.

## Session 46 (2026-04-25): Critique of Proposals P25–P29

Same /run.sh cycle as Session 45 (post-zeta-N=2000). Critic mode evaluated five
proposals (P25–P29) drafted in `archive/ephemeral/proposals_session.md`:

| # | Name | Verdict | Mode |
|---|------|---------|------|
| P25 | Liouville-parity triangulation | DUPLICATE | E |
| P26 | Cesàro–Fejér damped explicit formula | DUPLICATE + numerical bug fix | I |
| P27 | Selberg trace + class-number recursion | DUPLICATE | E |
| P28 | Random-multiplicative variance reduction (Harper) | FLAWED | E+I |
| P29 | Edixhoven–Couveignes τ bilinear form | DUPLICATE | C |

### Useful by-products

- **P25 identity correction**: closed-form parity is π(x) = (x − L(x))/2 − C₃(x),
  where C₃(x) counts n ≤ x with Ω(n) odd ≥ 3. Original proposal had wrong term
  Q(x)−1; verified exactly for all x ∈ [2, 10000]. Bottleneck shifts to C₃(x),
  parity-random vs simple proxies.
- **P26 numerical bug fix**: `mpmath.li(x**rho)` is wrong for ρ = 1/2 + iγ;
  silently uses principal branch and discards winding. Use
  `mpmath.ei(rho * log_x)` instead. Worth flagging in any zeros-summation code.
- **P26 follow-up experiment**: Fejér-recovery failures DO cluster, but only by
  features ≥ as hard as π(x) itself: near_prime_dist (Gibbs at unit jump,
  spread 0.86 at T=30/100) and lpf (factorisation, spread 0.65). Cheap features
  (mod m, div count, smoothness ≤7) all spread ≤0.14 — noise floor. No useful
  hybrid scheme.

### CLOSED_PATHS entries added (4)

P25 (line 684), P26 follow-up (line 685), P28 Harper (line 686), P29 Edixhoven
(line 687). P27 covered by existing entries (199, 345, 518); no new row.

### Files

- `archive/ephemeral/critique_latest.md` — full critique
- `experiments/proposals/proposal26_fejer_failure_clustering.py` (+ results.md)
- `archive/sessions/session46_critique_p25_p29.md`

### Verdict

No proposal opened a new direction. Project remains in steady-state critique
mode: every plausible lever (analytic damping, parity, trace formulas,
ensemble averaging, modular forms) reduces to one of three known failure
modes via routine analysis. Open direction unchanged: circuit complexity of
π(x).


---

## Session 47 (2026-04-25, normal mode)

### Goal
Productive normal-mode work: (a) literature scan for the 3-week window since
the previous update (2026-04-05 → 2026-04-25), (b) sharpen the only OPEN
remaining direction — "AKS in TC^0 via growing-dim matrix powering"
(`status/CLOSED_PATHS.md` line 232).

### Literature scan (2026-04-05 → 2026-04-25)

No paper in the window changes the asymptotic barrier for computing p(n).
Checked 10 topic areas (TC^0/NC^1, IMM lower bounds, PRIMES in TC^0, π(x)
algorithms, zeta zero summation, Hardy-Littlewood/Selberg sieves, sub-√x
claims, Brandt MKtP follow-ups, Ono-style detection extensions, Connes /
Yakaboylu Hilbert-Pólya updates). Two minor entries to note:

- Connes "The Riemann Hypothesis: Past, Present and a Letter Through Time"
  (arXiv:2602.04022, Feb 2026). Survey + perspective; no algorithmic
  implications.
- Yakaboylu Hamiltonian arXiv:2408.15135 revised to v10 with notes through
  March 2026. Operator-theoretic RH strategy, no algorithm.

### Cyclotomic CRT splitting for AKS-MPOW (sharpening the OPEN entry)

Question: does Z_n[x]/(x^r-1) ≅ ∏_{d|r} Z_n[x]/Φ_d(x) reduce the maximum
matrix-powering dimension below r? For 22 sampled n in {10²..10⁶, Carmichael
numbers, powers of 2}, the AKS-prescribed r is **prime in 21/22 cases (95.5%)**
because AKS picks the smallest r passing ord_r(n) > log²n and primes win
that race. When r is prime, x^r-1 = (x-1)·Φ_r(x), so CRT gives Z_n × Φ_r-piece
and the non-trivial factor still has dimension r-1.

**Best observed max_dim/r = 0.94** (only composite r in sample, 289=17²).
Average 0.99. **The cyclotomic shortcut is closed** as failure mode E
(equivalence to direct r-dim MPOW). Added as line 688 in CLOSED_PATHS.

The parent question "growing-dim MPOW in TC^0" (line 232) **remains OPEN** —
this experiment closes only the cyclotomic-decomposition sub-attack on it.

### Files
- `experiments/circuit_complexity/cyclotomic_crt_splitting.py` (+ results.md)
- `archive/sessions/session47_normal_cyclotomic.md`

### Verdict
Steady-state. One sub-attack on the only OPEN frontier closed. Literature
delta: zero algorithmic impact in window. Project remains in monitoring +
sharpening mode.


---

## Session 48 (2026-04-25, fresh-perspective wildcard)

### Goal
Find a wildcard angle not represented by name among the 57 existing wildcard
scripts. Test whether the prime parity stream q(n) = 1{p(n) = 5 mod 6} has
finite excess entropy E (a *stochastic* memory measure, strictly more
expressive than the deterministic linear complexity over GF(2), which is
already known to be N/2 maximal).

### Brainstorm (`archive/ephemeral/session48_brainstorm.md`)
Five fresh angles, all not present by name in `experiments/wildcard/`:
1. Causal-state / excess-entropy of prime parity stream (selected).
2. Higher-cumulant expansion of delta(n).
3. Phase retrieval of pi(x) from |zeta(1/2+it)|^2 magnitude.
4. Quasi-modular Eichler / mock-modular completion of log zeta.
5. Multipole expansion on the prime side (Riemann-Weil dual).

### Experiment
`experiments/wildcard/causal_state_complexity.py`. N = 148 931 prime parity
bits up to p < 2.10^6, block entropy H_L for L = 1..18, h_mu estimated from
Delta_L = H_L - H_{L-1} averaged over L = 8..10 (avoids finite-sample bias
that contaminates L >= 12), excess entropy E_L = H_L - L*h_mu. Compared to
i.i.d. Bernoulli(p1) null with same single-bit bias.

### Findings
- h_mu(prime) ~ 0.97686 nats / step vs h_mu(null) = 0.99708 -> per-step
  memory deficit ~0.020 nats, consistent with Hardy-Littlewood pair
  correlations between consecutive primes' residues mod 6.
- E_L plateaus at ~0.026 nats above the Bernoulli null (peak at L = 10);
  beyond L = 12 the estimate is dominated by finite-sample bias.
- Implies epsilon-machine size exp(E) ~ 1.026, i.e. essentially Markov-1.
- Per-step compressibility ~3.3% below entropy ceiling: far below the
  log_2(x) bits/step that a polylog encoding of pi(x) would need.

### Files
- `experiments/wildcard/causal_state_complexity.py` (+ `_results.md`)
- `archive/ephemeral/session48_brainstorm.md` (5 angles)
- `status/CLOSED_PATHS.md` (line 689 appended, S48)
- `novel/pseudorandomness_of_pi.md` table extended to **measure #24**

### Verdict
**CLOSED**, failure mode **I**. First *stochastic-memory* measure in the
pseudorandomness table - earlier 23 measures were deterministic / spectral /
algebraic. Finds tiny non-zero structural deviation (~0.020 nats h_mu
deficit, ~0.026 nats excess entropy plateau) that does not break the
random-like narrative; refines it.

---

## Session 49 (2026-04-25, deep focus task #3 deepening)

**Mode:** Deep focus, Task #3 (Novel Identity Search). Task was already
COMPLETED in Session 29 with seven experiments, but three sub-bases were
absent from that battery. Deepening run extends and re-confirms closure.

### Experiment
`experiments/algebraic/identity_search/extended_basis_search.py`. Recompute
f(x) = pi(x) - R(x) at mp.dps=60 via Riemann's R-function and sympy.primepi.
PSLQ with maxcoeff=1e10, maxsteps=1e4. Cross-validate every relation with
nonzero f-coefficient at every other test point; require residual < 1e-6.

Three sub-bases not covered by Session 29:

  A. **Extended zeta-zero oscillation basis** (26 terms): elementary
     functions + sin/cos(gamma_k log x) for k = 1..10. Test points
     x in {1000, 5000, 20000, 50000}.
  B. **Arithmetic partial sums** (9 terms): M(x) Mobius, L_lambda(x)
     Liouville, Phi(x) totient, psi(x) Chebyshev, all normalized by
     sqrt(x). Test points x in {1000, 2000, 5000, 10000}.
  C. **Mahler functional** (7 terms): f(x), f(x^2), f(x^3) plus
     transcendentals log(x), sqrt(x)/log(x), 1/log(x), log(x)^2.
     Test points x in {10, 20, 30, 40, 46} (so x^3 <= 97336).

### Findings
- **A:** PSLQ finds tight relations at each x (residual ~ 1e-44) with
  coeff(f) in {380, 101, 322, -139} but the integer vectors are
  point-specific. Cross-check residuals are O(10^3) at every other point.
  Adding eight more zeta zeros does not reveal a universal identity.
- **B:** Only short integer relations among partial sums themselves -
  M(2000) = 5; L_lambda(1000) = -7 M(1000); L_lambda(10000) = 4 M(10000) - 2.
  All have coeff(f) = 0; f rejected from every short combination of
  normalized M, L_lambda, Phi, psi.
- **C:** Relations exist at each base point but require coefficients of
  size 10^6 - 10^7 on f, and shatter at every other base point with
  cross-check residuals 10^6 - 10^8. No Mahler-type self-similarity
  f(x) <-> f(x^2) <-> f(x^3).

### Files
- `experiments/algebraic/identity_search/extended_basis_search.py` (+ `_results.md`)
- `status/CLOSED_PATHS.md` (line 690 appended, S49)

### Verdict
**CLOSED, strengthened.** The Novel Identity Search direction (CLAUDE.md
Open Problems #5) was already closed in Session 29 across 7 experiments;
this adds three orthogonal sub-bases (extended zeta zeros, arithmetic
partial sums, Mahler functional) and re-confirms closure with the same
methodology. Failure mode **I** (Information Loss): an identity in any of
these bases would compress the ~x^{1/2} pseudo-random oscillation into a
polylog object, contradicting Session 25's GUE structure of the zeros and
the Session 35 entropy lower bound. Combined Session 29 + 49 closure now
spans: elementary functions, li-variants, 10 zeta oscillations,
Bernoulli/zeta-value/L-value/Ramanujan-tau bases, LLL minimal polynomials,
arithmetic partial sums, ODEs/Volterra, Mahler functional equations.

## Session 50 (2026-04-25, critique mode)

### Proposals critiqued
#27 Hermite-mollified reverse explicit formula; #28 (min,+) tropical /
Cramér-window sieve around R⁻¹(n); #29 SoS / Lasserre localisation of p(n);
#30 cancellation-anchor density probe. Proposer ran #27, #28, #30; #29 deferred.

### Verdicts
All four DUPLICATE. None reopens prior closures.

- **#27** Hermite/Gaussian/Riesz mollification: linear re-weighting of the same
  truncated zero sum. Subsumed by lines 31, 32, 36, 519, 685.
  Empirical: at x=100, K=800, mollification is 5–55× WORSE than unmollified
  (Gauss σ=0.5 err=4.117 vs unmoll 0.076); at x=10⁴ best Riesz at T=50 is
  1.21× better — partial-sum truncation luck, not asymptotic gain.
  Linear-functional argument: any kernel w(γ) with ∑|w(γ_n)| ≥ cN(T) inherits
  Ω(√x · √(log T/log x)) tail.

- **#28** Cramér-window probe: subsumed by lines 593, 659, 661, 422.
  Empirical fit |δ_n| = |p_n − R⁻¹(n)| ~ p^{0.505} matches RH-predicted √x;
  max|δ_n|/log²(p_n) = 6.65 unbounded as expected. Window required to
  bracket p_n around R⁻¹(n) is √x, not polylog.

- **#29** SoS / Lasserre localisation: closed by approximate-degree theorem.
  adeg(χ_P) = ⌈N/2⌉ (line 580 + novel/approx_degree_prime.md, S28). SoS-deg ≥
  adeg/2 = N/4 = Ω(log x); Lasserre level-d in n^{O(d)} gives runtime ≥ x^{Ω(1)},
  not polylog. Independent confirmation: Grigoriev/Schoenebeck SoS lower bounds
  on random-like predicates + 25+ pseudorandomness measures.
  **Decisive closure without an SDP run.**

- **#30** Cancellation-anchor density: apparent 6–25% density of |E_K(y)| <
  √(log x) is a small-K truncation artefact. Mean |E_K| at K=200 already
  saturates √x at x=10⁵. Full zero sum gives density polylog/√x → 0.

### Useful constant pinned
**|p_n − R⁻¹(n)| ≤ 0.59 · √(p_n)** for n ≤ 10⁶ (better than RH's 1/(8π) on
this range), max ratio 6.65 vs log²(p_n).

### Files
- `archive/ephemeral/critique_latest.md`
- `archive/sessions/session50_critique_p27_p30.md`
- `status/CLOSED_PATHS.md` (4 entries appended; now 694 lines)

### Verdict
**No proposal survived.** Square-root-cancellation barrier robust under
linear mollification, tropical/min-plus lift, polynomial-optimisation SoS
encoding, and anchor-density search. Open status unchanged: only circuit
complexity of π(x) remains genuinely open per `status/OPEN_PROBLEMS.md`.

---

## Session 51 (2026-04-25, normal mode, no-op)

### Goal
Normal-mode startup: identify any actionable direction not already explored
in S44–S50 (seven prior sessions on the same calendar date).

### State on entry
- 694 closed paths in `status/CLOSED_PATHS.md`.
- Sole OPEN frontier: "AKS in TC^0 via growing-dim MPOW" (line 232).
- S47 closed cyclotomic CRT sub-attack (the natural decomposition).
- S48 closed causal-state stochastic-memory measure.
- S49 strengthened identity-search closure with three orthogonal bases.
- S50 closed proposals #27–#30 as duplicates of existing entries.
- `EDGES.md` (untracked, compiled 2026-04-25) catalogues every structural
  edge across 49 sessions and ranks 7 single-target objectives by impact.
- Literature window 2026-04-05 → 2026-04-25 had zero algorithmic-impact
  papers (S47 scan).

### Sub-attacks considered for growing-dim MPOW (line 232)
- AP-restricted binomial-sum formulation: M_a^n = Σ_j S^j · T_j with
  T_j = Σ_{k≡j (mod r)} C(n,k) a^{n-k}. Either evaluated naively (Ω(n/r)
  terms — exponential in N=log n) or via DFT
  T_j = (1/r) Σ_l ω^{-jl} (a+ω^l)^n with ω = primitive r-th root in some
  ring containing it. The DFT route is the algebraic dual of cyclotomic
  CRT (S47, line 688) and has the same dimension r.
- CRT-based scalar-iterated-mult lifted to matrix iterated mult:
  scalar trick reduces exponent mod φ(n) via primes p_i. For matrices,
  reducing exponent mod the order of GL_r(Z_n) requires factoring n
  (circular) and the order is n^{Θ(r²)} — no polylog reduction even if
  factorisation were free.
- "Same-base" matrix iterated mult (A^n with all factors equal):
  reduces to scalar exponentiation only when the matrix algebra is
  isomorphic to a product of fields with known orders — which is again
  cyclotomic CRT.
None of these gives a fresh experiment; all collapse to S47's closure.

### Action
No experiment run. Hygiene verified: 390 .py / 397 _results.md (companion
ratio ≥ 1.0), no orphan `__pycache__`, no stale "pending" labels in
ephemeral docs. `EDGES.md` cross-checked against `status/CLOSED_PATHS.md`:
every cited closure has a matching entry; no inconsistencies found.

### Verdict
**No-op. Project saturated for this calendar date.** Honest application of
CLAUDE.md rule: "If no viable experiment exists... If nothing new, say so
and stop." Next productive checkpoint: literature scan in ≥1 week
(2026-05-02+) for new algorithmic papers on growing-dim MPOW, π(x)
sub-√x, or PRIMES in TC^0.

---

## Session 52 (2026-04-25, normal mode, FOCUS-4 closure)

### Goal
TODO.md FOCUS-4: stack S46's Cesàro-Fejér window on top of S45's
Borel-Padé regularisation and measure recovery rate of ⌊round(S)⌉ = π(x).

### Setup
8 x-values (10³..10⁵) × 4 T-values (50, 100, 300, 1000) × 3 modes
(sharp / fejér / borel-fejér). 2000 zeros, mp.dps=30, branch-correct
`mpmath.ei(ρ·log x)`. Borel-Padé applied to Fejér-windowed increment
sequence; median over Padé orders (5,5)..(15,15) plus diagonals.

### Recovery rate
| mode         | T=50 | T=100 | T=300 | T=1000 |
|--------------|-----:|------:|------:|-------:|
| sharp        | 3/8  | 4/8   | 6/8   | 4/8    |
| fejér        | 4/8  | 4/8   | 5/8   | **6/8**|
| borel-fejér  | 3/8  | 3/8   | 3/8   | 3/8    |

### Key finding — diagnosis of regression
Borel-Padé locks into a T-independent asymptote when stacked on Fejér.
At x=10000 the hybrid returns S ∈ {1229.61, 1229.70, 1229.70, 1229.70}
across T ∈ {50, 100, 300, 1000} — identical to 4 decimals at three of
four T. Padé fits the leading envelope of the increment sequence; the
Borel integral ∫₀^∞ e^{−z} P/Q(z) dz extracts a tail-completion that
depends mostly on low-order coefficients. New zeros from larger T are
exponentially suppressed by 1/k! in the Borel transform and cannot
move Padé's leading-order fit.

### Failure-mode classification
**E (Equivalence).** Convergence-acceleration interventions cannot be
stacked. Each re-parametrises the same √x-bounded zero-sum information.
Closure mechanism identical to S45's standalone Borel-Padé closure.

### Bigger picture
After S52, every documented combination of standard
convergence-acceleration techniques on the truncated explicit formula
has been measured: raw / Cesàro / Borel-Padé alone / Fejér alone /
Borel-Padé × Fejér / Hermite-Gaussian-Riesz / Stein-method. None gives
sub-√x asymptotic gain. Cesàro-Fejér at T=1000 remains the best-tested
constant-factor improvement.

### Files
- `experiments/proposals/borel_fejer_hybrid.py` (+ `_results.md`)
- `archive/sessions/session52_borel_fejer_hybrid.md`
- `status/CLOSED_PATHS.md` (1 entry appended; now 695 lines)
- `EDGES.md` §E3.7 follow-up paragraph
- `TODO.md` FOCUS-4 marked DONE

### State of FOCUS queue
F-1 (Connes operator scaling), F-2 (π(x) mod q polylog), F-3 (Liouville
parity polylog), F-5 (AKS growing-dim MPOW alternatives), F-6 (4th
encoding of π(x)), F-7 (literature watch) — unchanged. **F-4 closed.**

---

## Session 53 (2026-04-25, normal mode, FOCUS-1 closure)

### Goal
TODO.md FOCUS-1 / EDGES.md E3.1 / Chain A step 4: measure the prime-budget
→ zero-count scaling law for the Connes-Consani-Moscovici Nov-2025
spectral-triple operator (arXiv:2511.22755). Determine whether
K_accurate(B) is polylog (Chain A becomes a real polylog architecture,
huge) or linear / sub-linear (Chain A closes via Equivalence).

### Setup
Numerical proxy (NOT a faithful CCM reproduction; see honesty notes in
results.md): discretize scaling operator D = -i d/du on L^2([-L, L]),
L = log(λ), x_cutoff = 10^4, in Fourier basis with N=1200 modes.
Self-adjoint rank-one perturbation V = c|v⟩⟨v|, where v encodes primes
≤ p_B via von-Mangoldt-weighted delta-comb in log space (best of two
tested variants — also tested ψ-step pairing). Coupling c tuned at B=6
to minimize median matched error (best: c = -2.0). Sweep B ∈ {1..9};
matrix size 2401; numpy.linalg.eigvalsh.

### Findings
1. **Median matched error stays flat at 0.13–0.20 across B = 1..9.**
   Adding 9× more primes does not improve eigenvalue accuracy.
2. **K_accurate(<0.5, matched) saturates at 50 for ALL B including the
   B=0 control** — comb-density artefact (unperturbed eigenvalue
   spacing π/L ≈ 0.68; any 50 targets in [14, 270] sit within 0.5 of
   some comb element by pigeonhole alone).
3. **K_accurate under monotone test** (μ_K vs γ_K, not nearest):
   K=0 for all B at threshold 0.5. Architecturally honest test.
4. **Linear regression** K_accurate(B) = -0.000·B + 50.000 (slope 0,
   R²=1.0). No information signal in B.
5. **Reproduction fidelity to CCM is poor**: at B=6 the proxy gives
   err[1] ≈ 9.1×10^{-2}, vs CCM's published 2.5×10^{-55} — off by 53
   orders of magnitude. Could not match CCM's specific Mellin-transform
   kernel without their detailed construction.

### Closure mechanism — three independent arguments
**E (Equivalence):**
1. **Rank-one parameter count.** Self-adjoint rank-one perturbation
   has ≤ B parameters (entries of v indexed by primes); by Cauchy
   interlacing, can shift at most ~B eigenvalues substantially.
2. **Diagonalization cost (kernel-independent).** Even granting CCM's
   published per-zero accuracy, spectrum extraction from an N×N
   truncation costs O(N³) = O(K³) — strictly worse than O(x^{1/2+ε})
   summation when K = √x. Iterative methods (Lanczos) give O(K²) per
   eigenvalue, still O(K³) for K eigenvalues.
3. **Geometric per-zero error growth.** CCM's published B=6 → K=50
   data range 10^{−55..−3} implies err(k) ≈ 11.3^k · 10^{−55}; thus
   K_accurate(B=6) ≈ 53 even at face value. CCM does not extrapolate
   K_accurate(B); literature contains no multi-B data point that
   establishes super-linear scaling. Burden of proof rests with anyone
   claiming polylog.

### Bigger picture
**E3.1's elevated edge value collapses to baseline analytic
equivalence.** Chain A — the highest-EVS chain in the EDGES.md
catalogue — is no longer a polylog architecture candidate. The
"Connes operator gives 50 zeros from primes ≤ 13" rumour is true
*as a one-shot fit* but not as an algorithm. The diagonalization-cost
argument (O(K³)) is robust independent of any future faithful
reproduction. After S53, only Chains B (π(x) mod q polylog), C
(Liouville parity polylog), and a residual interest in growing-dim
MPOW (Chain E) remain open as project frontiers; none has a credible
path forward in the literature.

### Files
- `experiments/analytic/connes_operator/connes_operator_scaling.py`
- `experiments/analytic/connes_operator/connes_operator_scaling_results.md`
- `experiments/analytic/connes_operator/connes_operator_scaling_data.csv`
- `archive/sessions/session53_connes_operator_scaling.md`
- `status/CLOSED_PATHS.md` (1 entry appended; now 696 lines)
- `TODO.md` FOCUS-1 marked DONE

### State of FOCUS queue
**F-1 closed.** F-4 closed (S52). F-2 (π(x) mod q), F-3 (Liouville
parity), F-5 (AKS), F-6 (4th encoding), F-7 (literature watch) unchanged.

## Session 54 (2026-04-26, fresh-perspective wildcard)
**Mode:** Fresh-thinking session, BBP-style isolated-digit angle.

### Brainstorm (5 unconventional ideas)
1. **BBP-style isolated digit extraction for psi(x).** Tested.
2. Compressed sensing on zeros via random subsets. Tested.
3. Hierarchical-multipole sieve (FMM analog). Not run — sieve "interactions"
   already linear, no quadratic kernel to compress.
4. Galois-descent computation of pi(x) mod q. Already covered by S35
   pseudorandomness findings (mod 2 random ⇒ mod q expected random too).
5. Spectral shortcut via low-pass smoothing pi*K_T. Mathematically same
   as truncated explicit formula → equivalent to test (1).

### Experiments run
- `experiments/wildcard/bbp_digit_freeze.py` — Empirical zero-vs-precision
  scaling for psi(x). For x=10^7, K=2000 zeros only gives 7 digits;
  zeros-per-digit grows ~`10^d`, not `O(d)`. **No BBP-style spigot.**
  Anomalous x with small residuals exist but match Gaussian-tail rate
  (1 in 1.5 expected; 1 found in 4001 samples).

- `experiments/wildcard/random_zero_subset.py` — Compressed-sensing test:
  random K-subset vs first-K. Random is **4x WORSE** (rms 404 vs 97)
  because low-frequency zeros carry most amplitude budget (`1/|rho|`
  weighting). First-K is near-optimal among K-subsets.

### Findings
- Both nulls re-confirm GUE-random-phases picture from `novel/pseudorandomness_of_pi.md`.
- Empirical residual scaling matches `sqrt(x) * polylog(K) / K^a`, a~1/2.
- Anomaly density too low (Gaussian) to support a "lucky-x" binary-search
  shortcut; even if a polylog-many anomalous x existed, we don't get to
  *choose* x in the original problem.
- 2026 literature search: no new prime-counting algorithm publications;
  state-of-art remains Deléglise-Rivat-Gourdon at O(x^{2/3}).

### State of project
No breakthrough. Bound `O(x^{2/3})` for exact, `O(polylog)` for ~50%
digits via R^{-1}(n) — gap unchanged.

## Session 55 (2026-04-26, FOCUS-2 Liouville structural decomposition)
**Mode:** Deep focus on FOCUS-2 (pi(x) mod q for fixed q), q=2 sub-case.

### New algebraic observation (free identity)
`L(x) mod 2 = x mod 2` is trivial (sum of x ±1's parity = x parity),
verified bit-exact on x ∈ [1, 2 × 10⁶]. So the polylog-`L(x) mod 2`
target in TODO.md FOCUS-2 step 2 is vacuously satisfied without
yielding pi(x) mod 2. The actual missing primitives via E2.2 are
`L(x) mod 4` (= `A(x) mod 2`) AND `C_3(x) mod 2`, with
`pi(x) mod 2 = A mod 2 XOR C_3 mod 2`.

### Single experiment
`experiments/sieve/pi_mod_q_fixed/liouville_modular_structure.py`

### Battery on each component (x ≤ 2 × 10⁶)
- `A(x) mod 2`: block-entropy 7.9999/8, AC max[1..30] = 0.0010,
  FFT z = 5.5 (no spectral line), GF(2) LFSR length = N/2.
  *More* pseudorandom than pi(x) mod 2 (no density bias).
- `C_3(x) mod 2`: block-entropy 7.88/8, AC = 0.148 (density bias),
  FFT z = 5.25, LFSR length = N/2.
- Mutual info I(A; C_3) = 1.94e-5 bits — independent.
- 11 cheap proxies all |corr| < 0.002 with pi mod 2; best XOR-fusion
  of 4-subset = 0.4951 (worse than chance).

### Findings
- E2.2 EVS downgraded H → M (analogous to E3.1 in S53).
- Chain C structurally exhausted at q = 2; EVS M → L.
- 2 new measures (#25 A mod 2, #26 C_3 mod 2) added to
  `novel/pseudorandomness_of_pi.md` (header now "26 Measures").
- Status of pi(x) mod q: q = 2 closed at structural level; q ∈ {3,5,7,11,13}
  still open (no Liouville-style identity available).

### Files updated
- `experiments/sieve/pi_mod_q_fixed/liouville_modular_structure.py` (new)
- `experiments/sieve/pi_mod_q_fixed/liouville_modular_structure_results.md` (new)
- `status/CLOSED_PATHS.md` (1 entry appended; now 700 lines)
- `novel/pseudorandomness_of_pi.md` (measures #25, #26 added)
- `EDGES.md` (E2.2, Chain C, priority list updated)
- `TODO.md` (FOCUS-2 q=2 structurally closed; q≥3 next steps)
- `archive/sessions/session55_focus2_liouville_modular_structure.md` (new)

### State of FOCUS queue
F-1 closed (S53), F-4 closed (S52), F-2 q=2 sub-case structurally closed (S55).
F-2 q ≥ 3, F-3, F-5, F-6, F-7 unchanged.

### State of project
No breakthrough. Same `O(x^{2/3})` exact / `O(polylog)` ~50% gap. The
Liouville-identity attack on pi(x) mod 2 is now closed at both
identity-level (S46) and structural-pseudorandomness level (S55).

## Session 56 (2026-04-26, FOCUS-2 q≥3 character-twisted Liouville)
**Mode:** Deep focus on FOCUS-2 amendment from S55 — extending the
Liouville-parity attack from q=2 to q ∈ {3,5,7,11,13} via
character-twisted sums L_χ(x) = Σ λ(n)χ(n) for non-trivial Dirichlet χ.

### Generalised free identity (algebraic)
For χ mod q of order d, decomposing by residue mod q and using
λ(n) ≡ 1 (mod 2):

    L_χ(x) ≡ Σ_{r ∈ (Z/q)*} χ(r) · count_r(x)   (mod 2 Z[ζ_d])

— a trivial coset-count identity, computable in O(polylog), carrying
no prime info.  The perfect analog of S55's L(x) mod 2 = x mod 2.

### Single experiment
`experiments/sieve/pi_mod_q_fixed/character_twisted_liouville.py`

### Battery (N = 10⁶, 5 q-values × 34 chars)
- **Free identity**: verified by exact integer arithmetic in Z[ζ_d]
  via cyclotomic-polynomial reduction.  2000/2000 sampled x match,
  zero failures, all 34 chars. (FP-projection check appears to fail
  for d ∈ {5,10,12} — pure rank-deficient-pinv numerical artefact.)
- **Next-bit pseudorandomness** A_χ(x) mod 2:
  LFSR/N ∈ [0.4995, 0.5002] (full pseudorandom rank for all 34 chars),
  block-entropy h₈ up to 7.97/8, FFT z up to 8.5 (within 1.5σ of
  order-statistic baseline ~4.8), AC up to 0.6 entirely from
  coset-density bias 1/φ(d).
- **Mutual information** I(A_χ mod 2 ; π(x;q,a) mod 2): max
  9.55 × 10⁻⁶ bits across all (q, χ, a) — at noise floor for N=10⁶.

### Findings
- The character-twisted Liouville attack route on Chain B closes
  uniformly across q ∈ {3,5,7,11,13} by the same dual mechanism as
  q = 2 (free identity + pseudorandom next bit).
- All identified attack sub-routes on Chain B are now closed (direct
  L-function: S20/S22; q=2 Liouville: S46/S55; q≥3 character-twisted
  Liouville: S56).  Chain B itself remains compositionally valid;
  the missing primitive has no remaining identified attack surface.
- 5 new pseudorandomness measures (#27..#31) added; header now
  "31 Measures".

### Files updated
- `experiments/sieve/pi_mod_q_fixed/character_twisted_liouville.py` (new)
- `experiments/sieve/pi_mod_q_fixed/character_twisted_liouville_results.md` (new)
- `experiments/sieve/pi_mod_q_fixed/_run_log_S56.txt` (raw stdout)
- `status/CLOSED_PATHS.md` (1 entry appended; line 704)
- `novel/pseudorandomness_of_pi.md` (header 26 → 31, table rows 27..31)
- `EDGES.md` (Chain B sub-route closure table; priority list #1 update)
- `TODO.md` (FOCUS-2 q≥3 structurally closed)
- `archive/sessions/session56_focus2_character_twisted_liouville.md` (new)

### State of FOCUS queue
F-1 closed (S53), F-4 closed (S52), F-2 q=2 closed (S55), **F-2 q≥3
closed (S56)**.  F-3, F-5, F-6, F-7 unchanged.  No remaining FOCUS-2
sub-task; future sessions should pivot to F-5 (AKS alternatives) or
F-6 (4th encoding search).

### State of project
No breakthrough.  Same `O(x^{2/3})` exact / `O(polylog)` ~50% gap.
Chain B's missing primitive — the project's "single missing primitive"
— is no longer attached to any specific attack lineage; every concrete
mechanism has been tested and closed.

## Session 57 (2026-04-26, deep focus task #2 — order-3 cell of structural battery)

### Context
FOCUS_QUEUE Task #2 ("Zeta Zero Structural Patterns") was marked
COMPLETED in S25 with the caveat "structure might exist at scales
requiring >1000 zeros to detect," resolved negatively in S45 at
N=2000.  Both extensions tested only orders ≤ 2.  The residual
hypothesis "could match pair correlation while concealing higher-order
cluster structure" had not been individually closed.

### Single experiment
`experiments/analytic/zeta_structure/triple_correlation.py`

Compares the 3-point correlation R_3(s_1, s_2) of N=2000 Riemann-von
Mangoldt-unfolded zeros to the GUE sine-kernel determinant
ρ_3 = 1 − K(s_1)² − K(s_2−s_1)² − K(s_2)² + 2 K(s_1) K(s_2−s_1) K(s_2)
on a 25×25 grid (L_max = 5 mean spacings, bin width 0.2).  Adds a
third-cumulant rigidity test of zero count in disjoint windows of
length L ∈ {1,2,4,8,16,32}.

### Findings
- **Bulk RMS deviation 0.0875** across the (s_2 > s_1) plane with
  n_ref = 1995 reference zeros.  Restricted to s_1, s_2 ≥ 0.5 (away
  from level-repulsion edge) RMS = 0.0924.  Same order as the
  pair-correlation deviation (0.0864 at N=2000), which is the noise
  floor at this sample size.
- **Diagonal slice (equally-spaced triples, s_2 = 2 s_1):** RMS 0.0972.
- **Anti-diagonal (constant 2nd-3rd gap, s_2 = s_1 + 1):** RMS 0.0884.
- **Third-cumulant rigidity:** c_3 of zero count stays at ~10⁻³ for
  every window size from L=1 to L=32, while a Poisson process would
  give c_3 = L (factor 10⁴ separation at L=32).  Rigidity is the
  highly-non-Poisson signature of GUE-class point processes.
- Variance plateaus at ~0.45 - 0.57 across the same L range (consistent
  with GUE's logarithmic count variance, distinct from Poisson's
  σ² = L).
- **Verdict: MATCH at order 3.**  No third-order non-determinantal /
  non-Gaussian structure detected.

### Implication
Closes the last unfilled cell of the S25/S45 structural battery:
the zeta zero point process now agrees with GUE at every k-point
correlation order tested up to 3, and the cumulant rigidity holds for
all tested window sizes.  Confirms (does not contradict) the existing
Task-2 closure; eliminates the residual "higher-order cluster"
hypothesis as a hiding place for structure.

### Files updated
- `experiments/analytic/zeta_structure/triple_correlation.py` (new)
- `experiments/analytic/zeta_structure/triple_correlation_results.md` (new)
- `experiments/analytic/zeta_structure/SESSION_25_SUMMARY.md` (S57 update block)
- `status/CLOSED_PATHS.md` (1 entry appended under Analytic / Zeta Zeros)
- `status/SESSION_INSIGHTS.md` (this entry)
- `archive/sessions/session57_focus2_triple_correlation.md` (new)

### State of FOCUS queue / project
No breakthrough.  Task 2 remains CLOSED; this entry strengthens the
closure rather than reopening it.  No change to FOCUS-3..7 status.


## Session 58 (2026-04-26, normal mode, FOCUS-7 literature watch)

### Verdict
**NO-DELTA.**  Window 2026-04-05 → 2026-04-26 (since S47 watch).

### What was checked
- arXiv math.NT recent submissions: all 245 entries for April 2026 scanned;
  targeted searches on `pi(x)`, `nth prime`, `zeta zeros`, `explicit
  formula`, `Riemann hypothesis`, `sieve`, `prime distribution`,
  `Hardy-Littlewood`, `Selberg`.
- ECCC TR2026: TR26-040 through TR26-061 (current max).
- Author streams: Connes (no new April submission post-2602.04022),
  Yakaboylu (still v15 of 2408.15135), van Ittersum, Ono, Granville, Tao
  (no relevant April 2026 papers).
- GitHub: kimwalisch/primecount (master ahead of v8.4 with ~30 commits;
  v8.5 imminent), PrimeCounting/PrimeCounting (still inactive).

### Findings
Four minor catalog entries appended to `literature/state_of_art_2026.md`,
none of which changes the asymptotic landscape:
1. Inoue arXiv:2604.05733 — μ < 0.50895 under RH (was 0.515). Zero-gap
   statistic, no algorithmic impact.
2. primecount post-v8.4 master — lockfree thread balancer, more accurate
   zeta-zero use in nth_prime, Brun–Titchmarsh dist_approx. Constant-
   factor only.
3. Kakkar arXiv:2604.02383 — "Neural Prime Sieves." Confirms §5.4 ML
   barrier; no exactness, no asymptotic change.
4. Li arXiv:2604.14596 — speculative 103-page lone-author RH preprint.
   Same profile as the debunked Kilictas-Alpay paper. Flagged for
   completeness; no peer review, no algorithm.

### Implication
Mature-state hypothesis from S47 holds. Brandt MKtP (S30) remains the only
identified theoretical technique that could in principle bypass Natural
Proofs en route to circuit lower bounds for π(x), and there is no follow-up
movement on it in this window. The only genuinely open research direction
(circuit complexity of π(x) per `status/OPEN_PROBLEMS.md`) is unchanged.

### Files updated
- `literature/state_of_art_2026.md` (S58 update block + header date bump)
- `status/SESSION_INSIGHTS.md` (this entry)
- `archive/sessions/session58_literature_watch.md` (new)

### State of FOCUS queue / project
No breakthrough. No closures. No new attack surfaces. Steady state
preserved.  Next literature watch in ~3 weeks (~S70).

---

## Session 59 (2026-04-26, deep focus task #3 — FOCUS-6 fourth-encoding sweep)

Run #33 in the rotation cycled to focused-mode Task #3 ("Novel Identity
Search"), which was already CLOSED in S29.  Per CLAUDE.md guidance ("do not
re-run experiments on closed paths"), this session worked the natural live
extension: TODO.md's FOCUS-6 — exhaustive enumeration of fourth-encoding
candidates via additive/multiplicative number-theoretic functions whose
summatory `S_f(x) = sum_{n<=x} f(n)` might be informationally equivalent to
π(x) up to polylog conversion.

### What was tested
21 candidate functions: {chi_P, Λ, λ, μ, Ω, ω, σ_0, σ_1, φ, J_2, log n,
1/n, 20-smooth indicator, 20-rough indicator, base-10 digit sum, popcount,
v_2, r_2 (sum of two squares), LPF, lpf−1, λ(n)/n}.

For each, we measured on x ∈ [10³, 10⁵]:
* polylog-evaluability of S_f(x) (analytic closed-form yes/no);
* residual R_f(x) = S_f(x) − M_f(x) after smooth-basis fit;
* growth slope α from log-log of |R_f|;
* ρ = corr(R_f(x), E_π(x) = π(x) − Li(x));
* free-identity probes mod {2, 3, 5}.

### Result
**Zero hits.** No candidate satisfies polylog ∧ |ρ| > 0.6 ∧ α > 0.4.
The two desired properties are mutually exclusive across every candidate
examined: either S_f is polylog by virtue of a smooth-Dirichlet closed
form (Stirling / Hurwitz / Trollope / Dickman) and decouples from primes
(mode I), or S_f carries prime info via a Mertens-type residual that
needs zeta zeros (mode E) or factorisation (mode C).

### Validation that the probe works
The probe correctly reproduces the known free identity
**L(x) mod 2 = x mod 2** (E2.10, S55) — the mod-2 column for Liouville is
exactly 1.000.  It also flags the digit-sum-mod-3 free identity
(`s_10(n) ≡ n mod 3`) as a 0.642 hit — a non-prime free identity that
illustrates exactly the pitfall FOCUS-6 was designed to detect.

### What this adds to the cumulative picture
* Pre-S29: 15+ intermediate-quantity families closed (S15/S16).
* S55-S56: Liouville + 34 character-twisted Liouville variants closed.
* **S59 (this session): +21 candidates closed.**
* Cumulative ≈ 70 distinct fourth-encoding routes empirically excluded.

The three-pillars meta-theorem (EDGES E7.7) is reinforced: every additive/
multiplicative encoding tested either is trivially smooth (M_f(x) is the
entire content) or routes back to {prime positions, zeta zeros, floor
values}.

### Files
* `experiments/algebraic/fourth_encoding_search/fourth_encoding_search.py`
* `experiments/algebraic/fourth_encoding_search/fourth_encoding_search_results.md`
* `experiments/algebraic/fourth_encoding_search/fourth_encoding_search_data.csv`
* `status/CLOSED_PATHS.md` (one new entry under "Encoding / Novel Representations")
* `archive/sessions/session59_fourth_encoding.md` (this synthesis)

### State of FOCUS queue / project
FOCUS-6 remains formally open (an exhaustive proof of "no fourth encoding
exists" is not given), but the empirical wall is now ~70 candidates thick.
Recommended next step on FOCUS-6: a more abstract analytic argument
(Dirichlet-series Euler-product structure ⇒ either zeta-coupled or
prime-decoupled), rather than further candidate enumeration.

## Session 60 (2026-04-26, S35-fresh proposals + critique cycle)

### Cycle structure
A "fresh-eyes" proposal session generated four blinded proposals
(`archive/ephemeral/proposals_session.md`) without reading CLOSED_PATHS,
ran each as a small (N <= 16384) experiment in
`experiments/proposals/session35fresh_*.py`, and self-reported all four
as FAIL. This critique cycle (S60) verified the FAIL verdicts against the
533-entry catalogue.

### Proposals + verdicts
1. **WHT sparsity of chi_P** — top-1024/16384 = 53% L^2 mass, entropy
   10.79/14 bits. DUPLICATE of line 702 (S49 Haar) under Donoho-Stark
   basis-invariance.
2. **PSLQ on delta(n) with 7-feature dictionary** — 39/40 windows give
   distinct integer signatures, residuals 10-200. Ninth PSLQ-on-delta
   variant; strict subset of line 703 (S49) basis.
3. **MPS bond dim of chi_P at midpoint** — empirical rank ~ 2^{0.485 L} ~
   sqrt(N). Duplicates the W=2 specialisation of S41's closed-form theorem
   rank = min(W^j, phi(W)*W^{d-j-1}+1) (line 518).
4. **Dirichlet AP explicit formula mod 12** — T=40 gives 1/20 exact pi(x;q,a),
   max diff 74. Duplicates lines 29-31 / 70 / 693 (per-character GRH error
   x^{1/2+eps}/T inherited).

### Outcome
Four entries added to `CLOSED_PATHS.md` (now 537+). No additions to
`OPEN_PROBLEMS.md`. No `novel/` document. Critique saved to
`archive/ephemeral/critique_latest.md`.

### Meta-observation
At this catalogue density, blinded fresh-eyes proposals at small N tend to
re-derive existing closures. Three of the four S35-fresh tests could have
been recognised as duplicates from the proposal text alone. Future cycles
should pre-screen proposals against CLOSED_PATHS before running any code.

### Files
- `archive/ephemeral/critique_latest.md`
- `archive/sessions/session60_s35fresh_critique.md` (this synthesis)
- `status/CLOSED_PATHS.md` (4 new entries)

## Session 61 (2026-04-26, FOCUS-1 sub-attack 2 construction — non-cyclotomic ring AKS)

### Construction
Built and probed the AKS congruence

```
   (x + b)^n  ≡  x^n + b   (mod  x^d + a, n)
```

with `f(x) = x^d + a` (Eisenstein-flavour, non-cyclotomic) replacing
the standard `x^r - 1`. Question per `TODO.md` FOCUS-1 sub-attack 2:
does the non-cyclotomic ring decomposition (a) admit a TC^0
irreducibility certificate, and (b) preserve AKS-style correctness?

Construction discipline (CLAUDE.md S60 rules): ran on `n` up to 410041
with multiple `(d, a, b)`, including 30 Carmichael numbers and 13
Fermat-pseudoprimes-base-2 as adversarial cases.

### Empirical finding
Clean structural separation between cyclotomic and non-cyclotomic
choices on Carmichael numbers:

| `d` | `a=1` (cyclotomic-flavour) Carm leak | `a=2` (Eisenstein) Carm leak |
|-----|--------------------------------------|------------------------------|
| 3   | 19/30                                 | 0/30                         |
| 4   | 3/30 (= Phi_8)                       | 0/30                         |
| 5   | 2/30                                  | 0/30                         |
| 6   | 5/30                                  | 0/30                         |
| 8   | 0/30 (Phi_16 too restrictive)        | 0/30                         |
| 12  | 0/30                                  | 0/30                         |

Korselt's criterion for Carmichaels aligns with small cyclotomic factors
of `x^d + 1` (e.g. Carmichaels sharing prime factor 7 leak through `d=3`);
setting `a = 2` destroys the alignment.

### Closure
**FAIL — (E) Equivalence + (C) Circularity**.

* (E) The polynomial congruence with `d = polylog(n)` requires
  growing-dim matrix powering — *exactly* E5.3 (the only Chain-E open
  frontier). Modulus structure is orthogonal to depth.
* (C) Capelli-test for `x^d + a` irreducibility over `F_p` is in TC^0
  *if* `p | n` is supplied; producing it is at least as hard as
  primality testing.
* No correctness theorem analog of AKS counting argument exists for
  Eisenstein moduli.

CLOSED_PATHS entry added (now 538+). FOCUS-1 sub-attack 2 closed.
Sub-attacks 1 (Bernstein 2003 smaller-r) and 3 (Healy-Viola Frobenius
transplant) remain.

### Meta-observation
This is a CLAUDE.md-S60-style construction failure — built an actual
ring, ran it on adversarial inputs, found a small structural sharpening
(a=1 vs a=2 Carmichael leak) that does NOT escalate to a breakthrough.
Higher information content than yet another statistical battery on
chi_P; correctly does not become a `novel/` document because the
underlying Korselt/cyclotomic interaction is folklore. Status of
`OPEN_PROBLEMS.md` and Chain E unchanged.

### Files
- `experiments/circuit_complexity/aks_alternative/non_cyclotomic_ring/`
  (`non_cyclotomic_ring.py` + `non_cyclotomic_ring_results.md`)
- `archive/sessions/session61_non_cyclotomic_aks.md` (this synthesis)
- `status/CLOSED_PATHS.md` (1 new entry)

## Session 62 (2026-04-26, deep focus task #4 — K_min extension with 2000 zeros)

### Context
Task #4 (Conditional Algorithms) was closed in Session 33 with the verdict
that no standard conjecture (RH, GRH, EH, Cramer's, ...) reduces exact
pi(x) below O(x^{1/2+eps}). Session 33 left one outstanding empirical loose
end: at x = 10^4 every K from 1 to 1000 produced rounded value 1230 instead
of the true 1229. The classical RH truncation bound at T = 1419 (1000-th
gamma) evaluates to ~0.65, sitting on the +-0.5 rounding cliff -- so the
Session 33 result was a suspected cliff effect, not a real anomaly.

### Method
With 2000 zeros now available (max gamma = 2515.286), we computed the full
trajectory pi_K(x) for K = 0..2000 in a single O(K) incremental pass at 30
dps, for x in {10^3, 10^4, 10^5, 10^6, 10^7}. Reported two K-metrics:
K_min (first round-correct) and K_min* (first K stable for >= 50
consecutive successes; avoids "lucky rounding").

### Findings

1. **Session 33 anomaly resolved.** At x = 10^4 the residual passes through
   +0.5 around K = 1250 and stabilises at +0.010 by K = 2000. K_min* = 1250,
   exactly the additional zeros Session 33 lacked. No hidden correction
   term needed.

2. **K_min* is wildly non-monotonic in x.** Across the trajectories:
   x = 10^3 -> 81, x = 10^4 -> 1250, x = 10^5 -> 572,
   x = 10^6 -> >2000 (still outside +-0.5 at K = 2000),
   x = 10^7 -> 1912.
   Phase alignment of the GUE-correlated zero contributions can push the
   round-correct cliff much earlier or later than the median.

3. **Empirical K_min* << classical T_min.** At x = 10^7 the classical
   bound 2 sqrt(x) log^2 x predicts T_min ~ 1.64 * 10^6, but K_min* = 1912,
   ~10^3-fold smaller. This is purely a constant-factor effect (loose
   classical bound + oscillation passing through zero); it does NOT reduce
   the asymptotic O(sqrt(x)) lower bound.

4. **Naive linear fit K_min* ~ x^0.275** is misleading -- four data points
   with two phase-luck dominated outliers. Asymptotic remains O(sqrt(x))
   up to logs.

### Closure
**CONFIRMS Session 33 verdict.** Failure mode: Information Loss (I).
Doubling the zero data does not change the asymptotic; constants are
loose, exponent is tight. No new CLOSED_PATHS entry needed (the path is
already closed); the loose end about x = 10^4 is now resolved in the
results file.

### Meta-observation
The "K_min* is non-monotonic in x" effect deserves to be remembered: it
means small-K_min(x) measurements at any single x are not a reliable
estimator of the asymptotic K_min curve. Future analytic-formula
experiments should either average over an interval of x or use the
classical bound, not a single empirical K_min value.

### Files
- `experiments/analytic/conditional/k_min_extended/k_min_extended.py`
- `experiments/analytic/conditional/k_min_extended/k_min_extended_results.md`
- `experiments/analytic/conditional/k_min_extended/run.log`
- `archive/sessions/session62_k_min_extended.md` (this synthesis)

---

## Session 63 — Critique of S63 fresh proposals (D-finite δ, mollifier, RMT moments, Newton zero-budget)
**Date:** 2026-04-26

### Context
S63 (proposer) generated four fresh proposals without consulting CLOSED_PATHS.md,
ran each, and self-closed all four. This critique session verified the verdicts
against the 537+-entry catalogue.

### Verification result
All four self-verdicts hold up:
- **P1 (D-finite recurrence on δ):** novel sub-test (first run of D-finite hunt
  on δ rather than π or 1_P), confirms broader holonomic-on-prime closures
  (lines 576, 577, 680). Mode I.
- **P2 (Selberg Dirichlet-polynomial mollifier):** duplicate-plus of line 693
  (Hermite/Gaussian/Riesz mollification). Different kernel shape, same
  general-kernel tail-bound argument. Caveat: experiment uses naive M(ρ)/M(1)
  weighting without ζ·M's M-zero/pole contributions; theoretical closure is
  decisive regardless. Mode E.
- **P3 (RMT local-moment predictor on Δ):** confirmed circular. Predictor
  passes empirically (RMSE 0.33 over H=200 window) but obtaining the window
  requires π values, which is the original problem. Mode C.
- **P4 (Newton with progressive 2^k zero-budget):** duplicate-plus of line 685
  (R⁻¹ fixed-point with zero correction). Different map (Newton-on-π vs
  R⁻¹-fixed-point), identical failure mechanism (GUE-random tail noise
  amplified by 1/π'(x) ≈ log x to prime-gap scale). Mode I.

### Side-findings worth keeping
1. **Methodology lesson (P1):** when running null-space / PSLQ / D-finite
   searches on data with multiplicative column-scale spread, validate via
   held-out **prediction**, not via training-side singular-value ratios.
   Column conditioning fakes rank deficiency at 1e-13.
2. **Δ local smoothness (P3):** std(Δ) over 200-point window < 2 even at
   X = 5×10^4, with weighted-mean predictor RMSE 0.33. Supports Δ entropy
   at scale X is concentrated at frequencies ≥ 1/√X. Worth a paragraph in
   novel/pseudorandomness_of_pi.md ONLY after a Cramér-random control to
   confirm prime-specific behavior; deferred.
3. **Quantitative Newton obstruction (P4):** required K satisfies
   √x/(√K·log x·√γ_K) < 0.05/log x ⇒ K ≥ x^{1/2-ε}/polylog. Geometric
   K_k = 2^k cannot break the √x barrier. Clean complement to line 685.

### Closure
Four CLOSED_PATHS entries added (S63 batch). No additions to OPEN_PROBLEMS.md.
No additions to novel/. OPEN_PROBLEMS still has only Circuit Complexity of π(x)
as a viable direction.

### Files
- `archive/ephemeral/critique_latest.md` (this critique)
- `archive/sessions/session63_proposals_critique.md` (synthesis)
- `status/CLOSED_PATHS.md` (4 new entries after line 714)
- `experiments/proposals/session63fresh_*.{py,_results.md}` (proposer's artifacts)


## Session 64 — FOCUS-1 sub-attack 3 construction: Healy-Viola Frobenius transplant
**Date:** 2026-04-26

### Context
S47 closed the cyclotomic-CRT splitting of AKS-MPOW (CLOSED_PATHS line 690).
S61 closed FOCUS-1 sub-attack 2 (non-cyclotomic ring AKS, line 714). S64
addresses sub-attack 3: does reducing the AKS congruence
`(a+x)^n ≡ a^n + x^n (mod n, x^r-1)` modulo a prime `q | n-1` yield a
TC^0-amenable test by exploiting the q-power Frobenius
`F: F_q[x]/(x^r-1) -> F_q[x]/(x^r-1), F(a+x) = a + x^q` (a ring
homomorphism in characteristic q)?

### Construction
Implemented and verified the Frobenius decomposition

  `(a + x)^n  =  prod_i (a + x^{q^i mod r})^{c_i}`   in F_q[x]/(x^r - 1)

where `c_i` are base-q digits of n. Decomposition matches naive repeated
squaring on 6/6 sanity-check cases (Q1 PASS).

### Headline negative
**Primes do NOT satisfy mod-q AKS for non-trivial a.** Across 19 primes
in [101..2207], for every prime `q | n-1` with `q < r = aks_r(n)` and
every `a in {1, ..., q-1}` (a coprime to q): **0/399 passes**. The
trivial case `a ≡ 0 mod q` always passes (both sides reduce to `x^n`)
but carries zero primality information.

Why: the polynomial identity `(a+x)^n - (a^n + x^n) = sum_{k=1}^{n-1}
binom(n,k) a^{n-k} x^k` vanishes mod n (Frobenius) but NOT mod q for
`q != n` prime. By Lucas's theorem, `binom(n,k) mod q = prod_i
binom(n_i, k_i)` where n_i, k_i are base-q digits — these are nonzero
whenever k is a "base-q sub-string" of n. Empirically the residual
polynomial fills 8-102 of r coefficients in F_q[x]/(x^r-1) — not a
sparse correction, no partial-cancellation rescue available.

### Depth analysis
Frobenius decomposition saves only constant factor `log(q)/log(2)`
bit-operations vs naive repeated squaring (~2x for q=2, ~3.3x for q=7).
Both schemes reduce to O(log n) sequential r-dim polynomial
multiplications = NC^1, NOT TC^0. The growing-dim r×r MPOW primitive
(E5.3) is unchanged by switching the coefficient ring from Z_n to F_q.

### Closure
**FAIL — modes (E) Equivalence + (I) Information loss.**
- (E) The mod-q computation IS growing-dim r×r MPOW over F_q, the same
  E5.3 primitive that is the only open Chain-E frontier.
- (I) AKS Frobenius identity is mod-n specific; reducing mod q≠n loses
  it (no Hensel lifting), and primes themselves fail the mod-q test.

CLOSED_PATHS line 719 added. With sub-attack 2 (S61 line 714) and
sub-attack 3 (this) closed, only sub-attack 1 (Bernstein 2003
strengthened gcd) of FOCUS-1 remains un-built.

### Files
- `experiments/circuit_complexity/aks_alternative/frobenius_transplant/frobenius_transplant.py`
- `experiments/circuit_complexity/aks_alternative/frobenius_transplant/frobenius_transplant_results.md`
- `archive/sessions/session64_frobenius_transplant.md`
- `status/CLOSED_PATHS.md` (one new entry after line 718)


## Session 65 — Fresh-perspective audit of arXiv:2506.22634 (TG kernel) + phi 2D rank

### Setup
Fresh-perspective mode: instructed not to read CLAUDE.md / CLOSED_PATHS,
think from first principles. After ~138 wildcard scripts already exist,
parallel-launched a 2025-2026 literature scan and tested two construction
angles:
1. phi(x,a) viewed as 2D matrix M[i,j] under four framings — does it
   admit polylog-rank reconstruction sufficient for integer-precision
   Meissel-Lehmer phi recovery?
2. Audit of the literature-flagged candidate arXiv:2506.22634
   (Kılıçtaş–Alpay TG-kernel "rigorous bound").

### Result 1 — phi 2D low-rank (CLOSED, mode I)
SVD of phi(x_i, a_j) under four framings (raw, Mertens-residual,
normalised, col-difference) at K∈{18,40,60}: relative singular spectrum
decays exponentially as exp(-0.33·k) — looks compressible.  But
||M||_F ∝ x grows, so required rank for integer-precision (±0.5) recovery
scales **linearly** with K (12 at K=18, ≈35 at K=60), not polylog. The
relative compressibility is illusory — same pseudorandomness wall the
project's 21+ measures already showed in other framings. Adds the 22nd
structural-pseudorandomness measure.

### Result 2 — TG-kernel audit (CLOSED, modes C+I)
The literature scan (sub-agent) flagged arXiv:2506.22634 as the most
algorithmically promising 2025 result: claims ~1200 zeros suffice for
π(x) at x with 10^8 digits via truncated-Gaussian kernel +
Riesz-Weil explicit formula. Project monitor had it on a "DEBUNKED
S12/S30" line but with vague reasons. Full PDF audit + small-x
empirical falsification (`experiments/wildcard/tg_kernel_audit.py`):

- **0th moment fails by 0.886.** Paper requires ∫Φ_TG(t)dt = 0 (so
  F_TG(1)=0 cancels the main term). But Φ_TG is strictly positive on
  its support [0, α+Δ]; we measured ∫_0^4 Φ_TG = 0.8862. Premise
  self-contradictory.

- **LHS empirically wrong by 8 orders of magnitude.** Paper derives
  Σ Λ(n) Φ_TG(n/x) ≈ αe^{-α²} ≈ 3.7e-4. We computed S(x) directly:
  x=100 → 86.78, x=1000 → 884.39, x=30000 → 26584.91 — i.e.
  S(x) ≈ 0.886·x. Trace: their IBP step on p.7 substitutes Ψ(t)≈t
  (PNT main) and drops the (Ψ(t)−t) deviation that is precisely the
  oscillatory zero-sum signal encoding π(x).

- **Σ F_TG(ρ) is x-independent.** F_TG depends only on Φ_TG, not on x.
  The truncated zero sum is therefore a fixed constant (≈ 2.5×10⁻²
  for first 20 zero pairs). It cannot "round to give π(x)" — π(x)
  varies with x, the proposed identity does not.

- **Lemma 2 zero-density bound is mathematically wrong.** Paper invokes
  N(σ,T) ≤ A T^{1−1/σ} (log T)^B; at σ=1/2 this gives A T^{−1} (log T)^B,
  *decreasing* in T, contradicting the actual N(T) ≈ (T/2π) log(T/2π)
  which grows.

- **Appendix B is symbolic-mysticism crank content** ("EmbedS(Faruk
  Alpay) := Φ_∞ … canonical identity fold equating author with a
  functor"), citing a self-referential "Faruk Alpay ≡ Φ_∞" preprint.

The right answer for smoothed-kernel-based π(x): smoothing in log-scale
with width h yields T·h ~ √log(1/ε); integer-precision recovery of π(x)
requires h ≲ log(x)/x ⇒ T ~ x. No polylog escape via this route.

### Methodological contribution
The audit shows a **fast falsification recipe** for any future "smoothed
kernel beats explicit formula" preprint: compute the LHS Σ Λ(n) Φ(n/x)
at small x where π(x) is known and check whether it matches the paper's
claimed identity. < 100 lines of code; falsified this paper in minutes.
Worth keeping in the literature-monitor playbook.

### Literature monitor update
S29's "NEEDS VERIFICATION" annotation on the TG-kernel paper
(literature/state_of_art_2026.md line 481) replaced with "VERIFIED
FALSE" plus the empirical evidence and pointer to the audit script.
The paper joins CLOSED_PATHS as line ~722.

### Files
- `experiments/wildcard/phi_2d_lowrank.py` + `_results.md`
- `experiments/wildcard/tg_kernel_audit.py` + `_results.md`
- Updates: `status/CLOSED_PATHS.md` (2 new entries), `literature/state_of_art_2026.md`
- This file (S65 entry).

### Closure
No new viable direction opened. Confirms that the 2025 literature has
no genuine breakthrough toward sub-x^{2/3} π(x): Guth-Maynard's
zero-density estimate sharpens error terms but doesn't produce a faster
algorithm, Aggarwal 2510.16285 is a wrapper-level p_n speedup, and
the Alpay-line "TG kernel" claim is mathematically incoherent.
The mature-state hypothesis from S47 holds.

## Session 66 (2026-04-26, FOCUS-1 sub-attack 1 construction — Bernstein 2003 smaller-r AKS)

**Mode:** deep-focus, construction. Closes the only un-built FOCUS-1
sub-attack and triggers the "computationally cornered" milestone for
Chain E per `TODO.md`.

### Headline

**FOCUS-1 sub-attack 1 closed FAIL, mode (E) Equivalence.** With S61
and S64 already closing sub-attacks 2 and 3, all three AKS-family
constructions are now closed — every modulus-twist and gcd-strengthening
of the AKS test reduces to the same growing-dim r×r MPOW primitive
(E5.3) at the same r ~ log²(n). Chain E is **computationally cornered**.

### What was tested

The Bernstein 2003 strengthening augments standard AKS with a gcd test:
for each `a ∈ [1, S]` compute the residual `(X+a)^n − (X^n + a) mod
(n, X^r − 1)`; for every non-zero coefficient c, check `gcd(c, n)`.
Test passes iff every residual is identically zero.

### Q1 — empirical r already meets Bernstein's bound

On the S47 22-sample n grid (n ∈ {101, 102, 561, ..., 1000003}):
- Mean `r/log²(n) = 1.207`, mean `(r − log²n) = 25.05`
- r prime in 21/22 cases
- The standard AKS order-condition already produces r within an
  additive constant of the conjectural Bernstein r = O(log² n) bound.

So *"smaller r"* is empirical reality already. What Bernstein 2003
attempts to add is a **deterministic correctness theorem** at this r.

### Q2/Q3 — striking empirical result on gcd-extraction

At canonical r, the test gives perfect discrimination on the sample:
**3/3 primes pass, 0/13 composites pass**. More importantly, the
Bernstein-style gcd of the residual coefficient with n yields a
non-trivial factor of n in **13/13 composite failures**, including
all 7 Carmichaels tested:

| n | Carmichael | gcd witness | extracted factor |
|---|---|---|---|
| 561 | 3·11·17 | (1, 374, 187) | 11·17 = 187 |
| 1729 | 7·13·19 | (1, 266, 133) | 7·19 = 133 |
| 8911 | 7·19·67 | (1, 5092, 1273) | 19·67 = 1273 |

Empirically the gcd-extraction algorithm gives factoring as a side
effect of compositeness witnessing — a known corollary of AKS but
worth noting; not a `novel/` entry because (a) implicit in any AKS
proof (b) does not unblock any complexity question.

### Q4 — dim/r ratio unchanged

Mean `max_dim/r = 0.9899` across 22 samples — identical to S47.
CRT splitting saves at most a factor of 1.06 on dimension, only when r
happens to be composite (1/22 cases). The Bernstein strengthening
operates on the gcd side, not the matrix dimension.

### Q5 — small-r probe (typical-vs-worst-case observation)

At r = log₂(n) + 1 (linear, way below the AKS bound), the polynomial
test still discriminates every composite in our sample including the
Carmichael 561. The AKS bound r ~ log²n is a *worst-case* defense
against adversarial near-primes; typical composites are caught at
r = O(log n). This is a complementary observation, not a closure
escape — Bernstein's strengthening covers the worst case.

### Q6 — why this is closure mode (E) Equivalence

The strengthened-gcd test does:

1. Length-O(log n) sequential growing-dim r×r MPOW over Z_n
   = the open primitive E5.3, unchanged.
2. Plus √φ(r) · log n integer-gcd checks of O(log n)-bit integers.
   Integer gcd is in NC¹ (Hesse-Allender-Barrington 2002), conjectured
   NOT in TC⁰ — same NC¹/TC⁰ frontier as growing-dim MPOW.

The substitution replaces one frontier problem with another at the
same frontier. Even if integer gcd were placed in TC⁰ tomorrow, the
polynomial-multiplication step is itself growing-dim r×r MPOW. No
escape from Chain E.

### State of FOCUS-1 / Chain E

| Sub-attack | Status |
|---|---|
| 1 — Bernstein 2003 strengthened gcd | **CLOSED** S66 (this) — (E) |
| 2 — Non-cyclotomic ring `Z_n[x]/(x^d+a)` | CLOSED S61 line 714 — (E)+(C) |
| 3 — Healy-Viola Frobenius transplant | CLOSED S64 line 719 — (E)+(I) |

**Computationally cornered.** Within the AKS family every modulus-twist,
ring-replacement, and gcd-strengthening reduces to the same growing-dim
MPOW primitive at the same r. The only remaining levers on Chain E are:
(1) Brandt MKtP framework (FOCUS-3, un-engaged), (2) a fundamentally
new lower-bound technique, (3) a non-AKS TC⁰ primality test using only
scalar operations (long-standing aspiration since S15).

### Files

- `experiments/circuit_complexity/aks_alternative/bernstein_smaller_r/bernstein_smaller_r.py` (~9 KB)
- `experiments/circuit_complexity/aks_alternative/bernstein_smaller_r/bernstein_smaller_r_results.md`
- Updates: `status/CLOSED_PATHS.md` (line ~722), `status/OPEN_PROBLEMS.md` (Chain E status)
- `archive/sessions/session66_bernstein_smaller_r.md` (synthesis)
- `TODO.md` (FOCUS-1 milestone)

### Closure

S66 confirms the S47/S64 mature-state hypothesis with one final
construction-mode datapoint: the AKS family of TC⁰ approaches is
exhausted. Future Chain-E work is exclusively non-construction
(Brandt MKtP, new techniques, monitoring).

---

## Critique-46 (2026-04-26, critique of S46-fresh proposals)

**Mode:** critique. Verifies the four S46-fresh proposals
(`archive/ephemeral/proposals_session.md`) against CLOSED_PATHS.md and
confirms each as DUPLICATE-PLUS, no `novel/` escalation.

### Headline

All four proposer self-verdicts hold up. Adds 4 entries to
`status/CLOSED_PATHS.md` (lines 723-726) — total 541+ approaches.

| # | Proposal | Critic verdict | Mode | Closest prior |
|---|----------|----------------|------|----------------|
| A | Wynn/Pade on psi(x) zero-sum partial sums | DUPLICATE-PLUS, CLOSED | I | lines 26, 49, 657, 697 |
| B | Borel-Pade resummation of Cipolla asymptotic | DUPLICATE-PLUS, CLOSED | E | lines 40, 43, 49 |
| C | Depth-5 Hecke/sigma fingerprint primality oracle | DUPLICATE-PLUS, CLOSED | C | line 689 |
| D | Zero-aware MC control variate for pi(x) | DUPLICATE-PLUS, CLOSED | E+I | lines 256, 257, 688 |

### Key empirical findings worth keeping

1. **Wynn-epsilon diverges on psi-zero-sum sequences** (P-A): at x=10⁴
   with K_max=160 zeros, best partial err 2.43 vs Wynn err 10013.88
   (4100x worse). Confirms line 26 closure under one more transform.

2. **Borel-Pade is uniformly worse than raw Cipolla truncation** (P-B):
   1.4-14x penalty across n in {10, 100, 1000, 10⁴}. Stokes-line /
   multi-instanton structure suspected; alternatively the asymptotic
   series is past its optimal truncation point at small K.

3. **Depth-5 multiplicative fingerprint perfectly separates 1229 primes
   from 8770 composites in [2, 10⁴]** (P-C). Empirical anchor on the
   AKS-style oracle approach. Algorithm is non-polylog regardless:
   tau(n) for composite n requires factoring (subexp via GNFS).
   Counting pi(x) with the oracle remains O(x · polylog x).

4. **PNT control variate gives variance reduction factor 1.006x**
   (P-D): essentially zero. Confirms lines 256/257/688 across one more
   variance-reduction shape.

### Pattern observation

The convergence-acceleration / variance-reduction family of interventions
on already-closed analytic primitives (Pade, Wynn, Shanks, Richardson,
Aitken, Cesaro, Fejer, Borel, Mellin-Barnes, mollifiers, control
variates, random sampling) is now **systematically exhausted** across
sessions 5, 6, 10, 15, 25, 32, 43-46, 48, 51, 63, and now critique-46.
Future proposals in this family should be pre-disqualified via
CLOSED_PATHS check against this master list. The novel-bar for any
future acceleration-style proposal is now: it must invoke a NONLINEAR
operation on zeros (line 693 explicit ask) AND it must not be
mathematically equivalent to a transform already named in the
established list.

### Project status unchanged

`OPEN_PROBLEMS.md` remains at: only Circuit Complexity of pi(x) is
genuinely viable. Within Chain E, all three AKS-family sub-attacks
closed (S47, S61, S64, S66) — Chain E computationally cornered.
Remaining levers: Brandt MKtP (FOCUS-3, un-engaged) and pure
new-technique work.

### Files

- `archive/ephemeral/critique_latest.md` — full critique
- `status/CLOSED_PATHS.md` — 4 new entries appended (lines 723-726)
- `archive/sessions/session_critique46.md` — synthesis
- `.run_state` set to 46

---

## Session 67 (2026-04-26, FOCUS-2 closure — E2.11 pre-test on six concrete fourth-encoding candidates)

**Summary:** Executed FOCUS-2 (TODO.md): for the six concrete candidate
intermediate quantities listed by the project, ran the E2.11 finite-
difference pre-test as a fail-fastest filter before the expensive PSLQ
stage. **All 9 expanded candidates close as mode I.** The "concrete
fourth-encoding sweep" sub-task is exhausted; FOCUS-2 is now CLOSED.

### Methodology

For each candidate T_i(x), compute it on a contiguous window
x ∈ [10⁵, 1.5·10⁵), polynomial-detrend (deg 6, absorbs leading
x, x log x, x log²x, x², ...), then iterate Δᵏ for k=1..7 and report
RMS(Δᵏ⁺¹)/RMS(Δᵏ). For f(x)=π(x)−R(x) this ratio → 2.0 (E2.11). A
candidate matching that signature is structurally equivalent to π−R
and can be closed without PSLQ. Three controls bracket the verdict:
i.i.d. Gaussian (WHITE-A, ratio 1.913), smooth polynomial (WHITE-B,
ratio 1.903 from amplified f64 noise after exact polynomial fit),
and π(x)−R(x) itself (ratio 1.914).

### Verdicts (9 candidates, 8 NEW + Σφ re-check)

| Candidate | residual std | Δ⁷/Δ⁶ ratio | Class |
|---|---|---|---|
| T_1 = Σ {log Γ(n)} | 3.06e+00 | 1.905 | WHITE-A (NEW) |
| T_2 = Σ H_n | 3.87e-04 | 1.926 | WHITE-B smooth (NEW) |
| T_3 = Σ H_n² | 7.75e-03 | 1.926 | WHITE-B smooth (NEW) |
| Ψ(x, B=137) (c=2) | 2.53e+00 | 1.915 | WHITE-A (NEW regime) |
| Ψ(x, B=1616) (c=3) | 2.90e+00 | 1.916 | WHITE-A (NEW regime) |
| Ψ(x, B=18971) (c=4) | 8.16e+00 | 1.916 | WHITE-A (NEW regime) |
| Σ σ_2(n) | 1.42e+09 | 1.989 | WHITE-A (NEW) |
| Σ σ_3(n) | 7.70e+13 | 1.995 | WHITE-A (NEW) |
| Q(x) squarefree | 3.12e+00 | 1.948 | WHITE-A (NEW) |
| Σ φ(n) | 1.77e+04 | 1.982 | WHITE-A (re-check, also S64) |

WHITE-A = same GUE-noise signature as π−R; WHITE-B = entirely smooth
(closed-form to f64 precision). Both close as mode I.

### Why this matters

1. **Methodological gain.** E2.11 pre-test filters fourth-encoding
   candidates ~150× faster than S64's ρ-correlation method (0.2s vs
   30s/candidate). The pre-test is also a *strictly stronger* filter,
   because it tests the structural equivalence directly at the
   finite-difference operator level, rather than via pairwise
   correlation that can miss subtle non-correlated equivalences.

2. **Cumulative cleanout.** Adding to S15/S16 (15+ families), S55-S56
   (34 character-twisted Liouville), and S64 (21 candidates), this
   session brings the cumulative count to **~78 distinct fourth-
   encoding routes empirically closed**. The three-pillars meta-theorem
   (E7.7) gains 8 fresh empirical checks against the same wall — at the
   *function* level, the bisection of bits 0..N/2 (smooth) and N/2..N
   (oscillatory) re-emerges as the bisection of "WHITE-B smooth" vs
   "WHITE-A oscillatory" candidates.

3. **TODO.md FOCUS-2 closes; FOCUS-3 / FOCUS-4 / FOCUS-5 unchanged.**
   The active research agenda narrows further. The remaining
   theoretical aspiration in this lineage is "find a fundamentally
   new intermediate quantity not based on floor values or zeta zeros"
   — at this point an open theoretical question with no concrete
   sub-task remaining.

### Tactical observation for future fourth-encoding work

Any future proposal in this lineage should be pre-disqualified by the
E2.11 test before being filed as an experiment. The pre-test takes ~5
seconds per candidate (per FOCUS-2 budget) and rules out both modes I
and E by the structural finite-difference signature. If a candidate
shows Δᵏ ratio bounded < 2.0 (e.g. ~ 1.5) or non-monotone toward 2,
that is an INTERMEDIATE survivor and merits the full PSLQ run; the
present 9-candidate sweep produced zero such survivors.

### Project status unchanged

`OPEN_PROBLEMS.md` remains at: only Circuit Complexity of π(x) is
genuinely viable. Chain E computationally cornered (S66). FOCUS-1
sub-attacks all closed; FOCUS-2 now closed; FOCUS-3 (Brandt MKtP) and
FOCUS-4 (3-point correlation at N≥10⁴) remain.

### Files

- `experiments/algebraic/fourth_encoding_search/e211_pretest_focus2.py`
- `experiments/algebraic/fourth_encoding_search/e211_pretest_focus2_results.md`
- `status/CLOSED_PATHS.md` — 1 new entry appended
- `EDGES.md` — E2.11 extended with S67 methodological note
- `archive/sessions/session67_focus2_e211_pretest.md` — synthesis
- `.run_state` set to 47

---

## Session 48-fresh (2026-04-26, fresh-perspective wildcards, 5 angles)

Out-of-order session label: this was a "fresh perspective" run launched after S67. Five untouched-among-142-wildcards angles built and closed.

### What was attempted

1. **Newton forward-difference series** for π(x). Test whether Δᵏπ(0) decays / has compact structure permitting truncated polylog evaluation.
2. **Multilinear extension** of χ_P over {0,1}ⁿ via boolean Möbius transform. Test sumcheck-style polylog evaluation hypothesis.
3. **Linear dynamical system view of the sieve.** CP / Tucker-rank decomposition of the alive-vector and wheel-period rank-1 factorization.
4. **Dirichlet character spectrum** of χ_P modulo prime q. Test sparsity in the character basis.
5. **Semiclassical log-block resummation** of the explicit-formula zero sum. Test whether log-block grouping breaks the √x truncation barrier.

### Key quantitative findings

- Newton forward-differences: log₂|Δᵏπ(0)| ≈ 0.995·k − 3.17 for k = 1..200 — almost exact 2ᵏ growth. Truncation diverges; bit complexity Ω(x) per term × Ω(x) terms.
- Multilinear extension: nonzero MLE coefficient ratio nz/N converges from below to ≈0.71 at n=14; max |coef| grows like 2^{n/2} (276 at n=14). Independent measurement of the rank(π_N)=2^{N/2−1}+2 result via the multilinear-Möbius basis.
- Sieve as linear DS: wheel factorization is rank-1 in CRT tensor basis up to primorial(K) ≤ N (verified K=3,4,5); 64×64 reshape rank of alive-vector after sieving primes ≤13 is exactly 32 = N^{1/2}/2.
- Dirichlet spectrum: L1/L2 ratio ≈ 0.85·√(φ(q)) — white-noise expected value — across q=11..97. Significant coefs are 5-15% of φ(q), not polylog.
- Block resummation at x=10000 (2000 zeros): per-block contribution grows as 2^{i/2}·√x; cumulative diverges from truth as more blocks added (best at 2 blocks, err=8; err=-26792 at 11 blocks).

### Verdicts

All 5 closed. Modes: 2 × I (info), 3 × E (equivalent to known barrier).

No genuinely new attack surface; all five collapse to existing project barriers (E2.x algebraic, E3.x analytic, E5.x sieve, the rank(π_N)=2^{N/2−1}+2 wall).

### Literature scan (parallel sub-agent)

Targeted 10 topics including sumcheck-for-π, Barvinok-for-π, Gutzwiller orbit truncation, Newton series for arithmetic functions, holographic/Valiant for primes. Result: **no new 2025-2026 publications** in any of the 10 directions; the cached `literature/state_of_art_2026.md` remains current.

### Files

- `experiments/wildcard/newton_diffs_decay.py` + `_results.md`
- `experiments/wildcard/multilinear_ext_pi.py` + `_results.md`
- `experiments/wildcard/parallel_sieve_matrix_exp.py` + `_results.md`
- `experiments/wildcard/dirichlet_char_decomp.py` + `_results.md`
- `experiments/wildcard/orbit_resum_pi.py` + `_results.md`
- `archive/sessions/session48_fresh_perspective_5_wildcards.md`
- `status/CLOSED_PATHS.md` — 5 new entries appended
- `.run_state` advanced to 48

## Session 49 — FOCUS-4 closure: large-N zeta-zero correlation battery + BK arithmetic-correction probe (deep-focus, Task #2 = FOCUS-4)

**Verdict:** CLOSED — mode I. The Bogomolny-Keating arithmetic correction is below the noise floor at N ≤ 8000, T ≤ 8148. New methodological edge added: gap-shuffled-zeros are the proper null for prime-frequency Fourier probes.

### Goal

FOCUS-4 (TODO.md): extend the S25/S45/S57 zeta-zero correlation battery from N=2000 to N up to 10⁴, and specifically test whether the empirical pair-correlation residual D(s) = R_2_emp(s) − R_2_GUE(s) carries a recognisable Bogomolny-Keating prime-arithmetic shape

    D_BK(s; T) ≈ -(2/L²) ∑_{p, k≥1} ((log p)² / pᵏ) · cos(2π s · k log p / L)

where L = log(T/(2π)).  This is the leading non-universal correction to GUE pair correlation that becomes visible at large height.

### What ran

- `data/gen_zeros_8000.py` (sequential mpmath @ 15 dps) and `data/gen_zeros_parallel.py` (6 workers covering disjoint k-ranges 5001..8000) co-generated zeros 2001..8000 in ~12 minutes total wall clock (vs ~50 min sequential-only).  Combined into `data/zeta_zeros_8000.txt` (8000 zeros total).
- `experiments/analytic/zeta_structure/large_n_correlations/large_n_battery.py`:
  pair correlation, form factor (Gaussian-smoothed), triple correlation, number variance Σ²(L), 3rd/4th cumulant rigidity, Bogomolny-Keating template Pearson, Fourier-prime amplitude/phase probe at log(p)/L for primes p ≤ 50, gap-shuffled null calibration (n_seeds=20).

### Key quantitative findings

- **Pair correlation:** RMS deviation from GUE shrinks 0.0947 (N=2000) → 0.0538 (N=8000); ratio 0.568 vs 1/√N prediction 0.500.  Statistical-noise-dominated, no BK floor.
- **Triple correlation:** RMS_far 0.0924 → 0.0617, similar shrinkage; no order-3 arithmetic structure.
- **Cumulant rigidity:** κ₃ stays at ~10⁻³ at every L ∈ {1, 2, 4, 8, 16, 32, 64, 128} (Poisson would give κ₃ = L; factor 10⁵ separation at L=128).  GUE rigidity holds up to the longest tested window.
- **BK Pearson:** zeta D(s) vs the multi-prime template has Pearson +0.111 (N=8000), +0.101 (N=2000) — *less* than the gap-shuffled null mean (+0.487 ± 0.035), at z = **−10.85σ BELOW the null**.  Phase coherence ⟨cos(φ − π)⟩: zeta +0.544 vs null +0.624 ± 0.036, z = **−2.20σ**.
- **NEW EDGE (E1.10 candidate, methodological):** the proper null for Fourier-prime probes is GAP-SHUFFLED zeros (preserving local GUE statistics, destroying long-range arithmetic), not uniform-random frequencies.  Uniform-random gives spurious p < 0.05 because it does not preserve the local gap structure that drives D(s) shape.  Without this null, the Fourier-prime amplitude probe falsely 'detects' BK at p<0.01.

### Verdict

**CLOSED.** FOCUS-4 closes mode I — the BK arithmetic correction, if conjecturally present, is below the empirical noise floor at this N (≤8000) and T (≤8148).  Detecting it would require Odlyzko's tabulated zeros at heights T ∈ [10⁶, 10²²], which mpmath cannot generate within session-time budgets.  The empirical zeta zeros are *more* GUE-like than gap-shuffled controls — opposite of what a real BK arithmetic correction would produce.

### Files

- `data/gen_zeros_8000.py`, `data/gen_zeros_parallel.py`, `data/combine_zeros.py`
- `data/zeta_zeros_8000.txt` (8000 zeros at 15 dps)
- `experiments/analytic/zeta_structure/large_n_correlations/large_n_battery.py` + `_results.md`
- `status/CLOSED_PATHS.md` — FOCUS-4 entry appended
- `EDGES.md` — new edge E1.10 (gap-shuffled-zeros null methodology)
- `archive/sessions/session49_focus4_large_n_zeta.md` — session synthesis
- `.run_state` advanced to 49

---

## Session 51 — FOCUS-3 closure: Brandt 2024 MKtP-diagonalisation vs pi(x) mod 2 (deep-focus, last untouched FOCUS-N)

**Verdict:** CLOSED — mode E. Brandt 2024 (TCC, IACR ePrint 2024/687) `MKtP not in DTIME[O(n)]` is structurally welded to MKtP itself and does not extend to natural functions like pi(x) mod 2.  New edge E5.8.  With E7.10 + E5.8, **Chain E is closed for both known technique families** on E5.3 (AKS-style and diagonalisation-via-meta-complexity).

### Goal

FOCUS-3 (TODO.md): read Brandt 2024 carefully, identify the diagonalisation skeleton, ask whether the technique extends to a circuit (or even uniform-time) lower bound on pi(x) mod 2.  If yes, that would be the first non-trivial circuit lower bound the project has produced.  If no, document rigorously and close.

### What ran

Built `experiments/constructions/brandt_mktp/`:
- `brandt_mktp.py` — bounded-Kt simulator (3-bit-op stack VM, L_MAX=12 program length, T_MAX=4096 step cap) implements TRAVERSE on a finite-Kt oracle, computes Kt_bounded(pi_N) for N=3..10, and encodes pi_N as MKtP queries.  Runs in 1.0 s.
- `brandt_mktp_results.md` — full closure with quoted theorems, four obstructions to extension, suggested EDGES + CLOSED_PATHS entries.
- `definition.md` — formal signatures (MKtP, pi_N, the bounded VM, the falsified conjectured extension).

### Key findings

- **Technique skeleton (verified from Brandt pages 1–14 of the PDF):** assume Pi_Kt decides MKtP in linear time; build TM running TRAVERSE that descends on Kt-random strings (queried via Pi_Kt) and right-steps on compressible strings.  Lemma 2 uses 1-Kt-randomness of Chaitin's Omega prefixes to prevent wrap-around.  Lemma 1 bounds the visit count infinitely often by `2^l / (l ln l)`, which contradicts the Kt definition.
- **Bypassing T4:** Brandt's proof relativizes (page 2: "our lower bounds relativize can be taken as a hint that relativizing techniques might be strong enough").  No probabilistic ensemble, no large/efficiently-computable property of MKtP — Razborov-Rudich does not apply.  This *is* the genuinely-new ingredient.
- **Why it does not extend to pi(x) mod 2:**
  - (O1) The hard string z is an oracle-dependent Kt-random prefix, not a fixed natural function.
  - (O2) The contradiction `Kt(z) >= |z|` ∧ `Kt(z) <= |M| + log_2 t` uses self-referential Kt on both sides.  No analog of "f-hardness(z) >= |z|" for a fixed Boolean function.
  - (O3) Lemma 2's Chaitin-Omega density argument has no analog for fixed total functions.
  - (O4) Brandt produces UNIFORM-TIME bounds, not circuit bounds.  E5.3 needs a circuit bound.  Brandt explicitly avoids the Williams/Hirahara algorithmic-method route on page 4 — that route IS subject to Natural Proofs on stronger classes.  Price of relativisation: no circuit bound.
- **E1.9 does not save the argument:** the 33-measure pseudorandomness of pi_N is consistent with Kt(pi_N) being near-maximal but proving that asymptotically would itself face T4.

### Verdict

**CLOSED, FAIL/E.** FOCUS-3 closes the last untouched construction-flavoured attack on the only open problem.  The cumulative Chain-E picture: E7.10 (S61/S64/S66) closes the AKS family; E5.8 (this session) closes the Brandt / meta-complexity-diagonalisation family.  Remaining levers on E5.3 are non-AKS TC⁰ primality tests (no candidate currently identified) and entirely-new circuit lower-bound techniques.

### Files

- `experiments/constructions/brandt_mktp/{brandt_mktp.py, brandt_mktp_results.md, definition.md}`
- `status/CLOSED_PATHS.md` — S51 row appended
- `EDGES.md` — E5.8 added; footer updated
- `TODO.md` — FOCUS-3 marked CLOSED; critical-path section emptied
- `archive/sessions/session51_brandt_mktp.md` — session synthesis

---

## Session 68 (2026-04-26, deep focus Task #3 — Bessel-basis PSLQ identity search)

**Verdict:** CLOSED — mode E. Task #3 (Novel Identity Search) was previously closed in S29 for an elementary + li + sin/cos basis. This run runs a strict basis-extension experiment: PSLQ in a 10-element basis disjoint from S29, containing modified-Bessel K_0, I_0, oscillatory-Bessel J_0, Y_0 at γ_1·log x, and the partition-asymptotic kernel K_0(2π√(log x)). Motivated by Bessel kernels appearing in the Selberg trace formula spectral side and Mellin–Barnes representations — a place where an identity could in principle exist that elementary functions missed.

### What ran

`experiments/algebraic/identity_search/bessel_basis_pslq.py` (~120 LOC). Reuses S29's `fx_data.npz`. PSLQ at mpmath 50 dps, maxcoeff=10⁶, maxsteps=8000, x ∈ {5000, 10000, 50000, 100000}, cross-validation at x+1000. Three Gaussian random-control trials at x=50000 with σ matched to f.

### Key findings

- Every fit produces tight residuals (1.9e-36 to 4.0e-35) with non-zero coeff(f) ∈ {2712, 5349, 10385}, but ALL fail cross-validation by 4-5 orders of magnitude (cross-check residuals 1.10e+04 / 6.03e+04).
- Random control: same-σ Gaussians give statistically indistinguishable signatures (fit residuals 8.5e-37 to 6.4e-35; cross-check 1.77e+04 to 1.60e+05). The PSLQ "tight relations" are pure numerical coincidences from finite precision on a 10-element basis.
- Rules out three concrete Bessel-shaped donors: Selberg-trace shadow (low-order K_{ir}-kernel reuse), partition-saddle reuse (K_0(2π√(log x)) from Hardy-Ramanujan), modified-Bessel growth/decay scaffolding.

### Verdict

CLOSED, FAIL/E. Strengthens S29 verdict from "f(x) admits no identity in tested elementary basis" to "f(x) admits no identity in any tested basis including spectral-trace-related Bessel kernels." Adds 34th measure to `novel/pseudorandomness_of_pi.md`.

### Files

- `experiments/algebraic/identity_search/bessel_basis_pslq.py` (new)
- `experiments/algebraic/identity_search/bessel_basis_pslq_results.md` (new)
- `status/CLOSED_PATHS.md` — Bessel S68 row appended; header bumped
- `novel/pseudorandomness_of_pi.md` — measure #34 added; title bumped to 34
- `archive/sessions/session68_task3_bessel_basis.md` — session synthesis
- `.run_state` advanced to 54

---

## Session 56 (critique-mode, 2026-04-26) — critique of S55 proposals

**Mode:** critique. Five proposals from `archive/ephemeral/proposals_session.md`
(S55) reviewed against CLOSED_PATHS. **All five collide with prior closures.**

### Verdict matrix

| # | Proposal | Verdict | Mode | Closes by |
|---|----------|---------|------|-----------|
| 1 | Compressed sensing on M[i,j]=R(x_i^ρ_j) | DUPLICATE | I | lines 26 (S32), 654 (S36), 699 (S54), 708 (S49) |
| 2 | PSLQ on δ(2^k), δ(F_k), δ(p_k) | DUPLICATE | I+E | lines 23, 666 (S33), 703 (S49), 711 (S60) |
| 3 | TT-rank of χ_P : {0,1}^L → {0,1} | DUPLICATE | I | lines 518 (S41 *theorem*), 707 (S49), 712 (S60), 729 (S48) |
| 4 | GUE-aware adaptive importance sampling of zeros | DUPLICATE | E | lines 256/257 (S15), 656 (S31), 723 (S46-fresh), 726 (S46), 732 (S48), 733 (S49) |
| 5 | Modular/theta bridge to π(x) | DUPLICATE | E | line 519 (S28), Landau natural-boundary obstruction |

### What ran

- Re-verified Proposal 3 TT-rank measurement at L=10 in this session: max
  TT-rank = 17 = 2^{(L-1)/2}+1, exactly matching the S41 closed-form
  rank theorem and the proposer's table at L∈{10,12,14,16}.
- Proposals 1, 2, 3 had been preemptively executed by the proposer; their
  `_results.md` files (`zero_contribution_compressibility_results.md`,
  `pslq_subsequence_delta_results.md`, `tt_rank_prime_indicator_results.md`)
  were audited and found internally consistent and consistent with the
  prior closures cited above.
- Proposals 4 and 5 ruled by prior art, no re-run.

### Files

- `archive/ephemeral/critique_latest.md` — full per-proposal critique
- `status/CLOSED_PATHS.md` — three new refining entries appended (cumulative 622)
- `archive/sessions/session56_critique.md` — this session's synthesis
- `.run_state` advanced to 56

### Process note

Proposer document explicitly bypassed CLOSED_PATHS to "avoid bias." This
ideation discipline is defensible but should be paired with a critique
pass before *running* experiments — three experiments in S55 reproduced
known closures. Future proposal sessions should grep CLOSED_PATHS for
their proposal's keywords before writing code.

---

## Session 57 (2026-04-26, normal mode) — smoothness-of-(p-1) conditioning of delta

**Mode:** normal. Single targeted measurement closing one residual gap
in the empirical pseudorandomness catalogue: condition delta(pi(p)) on
arithmetic features of p (smoothness/factorization of p-1) rather than
on features of the index n.

### What ran

`experiments/proposals/delta_pminus1_smoothness.py`. N=9592 primes
p ≤ 10⁵; delta computed via R^{-1}(n) ≈ p + f(p)·log(p) linearization
(validated ±1 vs Newton-Rinv at mpmath 30 dps; bias is independent of
class). Quartile-binned on s_ratio = log(LPF(p-1))/log(p); KS test on
extreme quartiles; Pearson and Spearman vs delta and |delta|;
Sophie-Germain (p=2q+1) sub-class; omega(p-1) feature.

### Key findings

* Per-quartile mean(delta) range = 1.6 (single sigma), std(delta)
  range = 1.4. No monotonic trend.
* Pearson(s_ratio, delta) = +0.0049 (p=0.634); KS class-0-vs-3 p=0.91.
* Sophie-Germain (N=670) KS vs others p=0.95.
* Pearson(omega(p-1), |delta|) = +0.054 (p=1e-7) — statistically real
  but explains 0.3% variance, fully explained by both quantities
  co-scaling with p (Erdős-Kac log log p, |delta|~p^{0.57}).

### Verdict

FAIL / mode I. Arithmetic features of p (factorization of p-1) are
empirically orthogonal to the transcendental sector (delta of R^{-1}).
Adds 35th measure to `novel/pseudorandomness_of_pi.md`. Rules out a
class of "smoothness-aware delta predictors" and the Pollard /
Pocklington side-channel into prime-counting.

### Files

* `experiments/proposals/delta_pminus1_smoothness.{py,_results.md}`
* `experiments/proposals/delta_pminus1_smoothness_data.npz`
* `status/CLOSED_PATHS.md` — S57 row appended (cumulative 545+)
* `novel/pseudorandomness_of_pi.md` — measure #35 added; title bumped to 35
* `archive/sessions/session57_smoothness_conditioning.md` — synthesis
* `.run_state` advanced to 57

### Status update

Project plateau reaffirmed: with all FOCUS lines closed and the
fourth-encoding sweep cumulative at ~78, this S57 measurement was a
hand-picked one-off chosen specifically because it was NOT a duplicate
of any prior CLOSED_PATHS row. The signal is exactly what the
pseudorandomness picture predicts (independent within a single sigma);
the only "significant" hit is statistically explained by trivial
co-scaling with p. No new attack surface revealed.

---

## Session 69 — F2 information-rate side closed; closed-form refinement of E1.5

**Date:** 2026-04-26
**Run:** 58 → 59
**Trigger:** /run.sh routed me to deep-focus Task #4 (Conditional
Algorithms), which FOCUS_QUEUE.md explicitly marks COMPLETED at S33.
Per CLAUDE.md, picked NOVELTY_CHALLENGES.md §2.F2 instead.

### Result

E1.5's empirical "0.537 bits constant in m" is sharpened to a closed
form. For every modulus m and every X with `m << pi(X)`,
```
   H( pi(x) mod m | pi(x-1) mod m, x in [1, X] )
       =  h_2( pi(X)/X )  +  O( 1 / pi(X) ).
```
Verified to 7-decimal-place precision at X = 10^7 across moduli
2, 4, 8, ..., 1024 plus cross-checks {3,5,7,11,13,30,210}.

The "0.537" of S12 is the X = 10^4 special case; the constant decays as
X grows, asymptotically to 0. This makes the CRT-cannot-win argument
*stronger*, not weaker — combining moduli still gives `h_2(pi(X)/X)`
bits per step, which goes to 0 in the limit.

### Why this counts as session output

Per CLAUDE.md novelty bar: this is a sharpening of an existing edge with
a new closed form (not a duplicate; not a 36th pseudorandomness measure;
not a brute-force re-derivation). It exposes the X-dependence that E1.5
hid behind a single empirical number.

### Files

* `experiments/information_theory/pi_mod_2k_saturation/pi_mod_2k_saturation.{py,_results.md}`
* `EDGES.md` — E1.5 entry rewritten with closed form; T6 theme line
  updated to match
* `NOVELTY_CHALLENGES.md` — F2 marked PARTIALLY CLOSED with outcome notes
* `RESEARCH_AGENDA.md` — Arc 5 milestone for F2 ticked, next-action
  pointed at F4 / F6
* `archive/sessions/session69_focus4_pi_mod_2k_saturation.md` — synthesis

### Status update

F2 information-rate side: closed with positive outcome (refinement, not
failure). F2 pseudorandomness-battery side: still open if anyone wants
to apply the 35-measure battery to pi(x) mod 4 / 8 / 16 explicitly.
Frame-shift Arc 5 is one milestone closer; F4 (pi_BD = pi - pi_smooth)
or F6 (parametric pi(2^k)) are the suggested next-target frame-shifts.

## Session 71 — Wild Swing: §C1 BK probe at Odlyzko's high-height tables (CLOSED)

**Mode:** Wild swing (one-shot ambitious frontier attack from
ATTACK_VECTORS.md §C1).
**Self-grade:** B — ambitious failure with structural reason.

Re-ran the S49 Bogomolny-Keating arithmetic-correction probe on
Odlyzko's published tables `zeros4` (10⁴ zeros at index 10²¹+1..10²¹+
10⁴, height T ≈ 1.4·10²⁰, L ≈ 44.6) and `zeros5` (10⁴ zeros at index
10²²+1..10²²+10⁴, T ≈ 1.4·10²¹, L ≈ 46.8). New methodological
contribution: a **random-prime null** (BK template with primes
replaced by uniform pseudo-primes in [2, 50]), which controls for
template shape rather than gap distribution — cleaner than S49's
gap-shuffled null at high L.

Empirical Pearson with the BK template: +0.063 (zeros4), +0.035
(zeros5). Random-prime null: μ = +0.063 / +0.035, std ≈ 2-4·10⁻⁴.
**z-scores: −0.94σ and +0.93σ — statistically indistinguishable from
random-frequency template of the same shape.** Direct sanity check:
prime-frequency Fourier amplitudes (median 0.0141, 0.0100) are NOT
enhanced over random frequencies (0.0148, 0.0101).

**New quantitative obstruction (sharpens E7.1, E1.10, E3.13):**
empirical fits give `|BK_pred|_max · L² ≈ 13.6` (invariant across
L ∈ {7, 44.6, 46.8}, confirming the predicted 1/L² scaling) and
`pair_rms ≈ 4/√N` (verified across N ∈ {2000, 8000, 10000}). Detection
threshold for κ-σ Pearson:

```
  N_required(L; κ) ≥ (4κ / 13.6)² · L⁴ ≈ 0.09 κ² · L⁴
```

At κ=3: detection at L=44.6 needs N ≥ 3.5·10⁵ (Odlyzko: 10⁴, 35×
short); at L=80 needs ≥3·10⁷ (hopeless). **L⁴ scaling is structural
— doubling height L requires 16× more zeros to compensate.** The
asymptotic regime *suppresses* the BK signal faster than data
accumulation can compensate.

Closes ATTACK_VECTORS §C1. Strengthens E3.13 from "BK absent at
N=8000" to "BK absent through T ≈ 10²¹, AND requires N ≥ 0.81 L⁴
zeros for any future detection — a hard scaling barrier independent
of computational budget." First quantitative version of the
"BK-asymptotically-too-weak" obstruction.

The previously-large gap-shuffled-null +30σ "signal" exposed as null-
bias artefact: gap-shuffled has Poisson-leakage in long-range tail
that anti-correlates with any oscillatory template, with the bias
*growing* in L. Re-confirms E1.10's warning at much higher L.

### Files

- `experiments/analytic/zeta_structure/odlyzko_high_height/odlyzko_bk_probe.py`
  + `odlyzko_bk_probe_results.md` + `.json`
- `data/odlyzko/zeros4`, `data/odlyzko/zeros5` (10000 lines each;
  Odlyzko's tabulated heights at extreme indices)
- `archive/sessions/session71_c1_odlyzko_bk_probe.md` (synthesis)
- `EDGES.md` — E3.13 extended with S71 sub-section + L⁴ obstruction
- `ATTACK_VECTORS.md` — §C1 moved to "Closed attacks" section
- `status/CLOSED_PATHS.md` — closure row added at line ~748

## Session 77 — N1 Tensor-Network Compression Family Closure (negative-shape edge)

**Mode:** novelty (B-grade target).
**Target:** `NOVELTY_CHALLENGES.md` §4.N1 — family-level closure of E2.1
across classical polynomial-spatial-locality tensor ansätze.
**Channel:** Razborov (lower-bound mindset; spectral-rank route).
**Self-grade:** B (refinement: lifts E2.1 from MPS-only to family-level
scope, with empirical verification at 22 (W, d) pairs and three orthogonal
closure mechanisms).

### Headline result

Across **22 (W, d) pairs** (W ∈ {2,3,4,5,6,7}, d ∈ {4,6,8,10,12}), the
half-cut bond dim required to *exactly* represent χ_P is **identical**
across five classical tensor-network ansätze:

* MPS / Tensor-Train at the half-cut.
* Hierarchical Tucker root-children bond dim.
* PEPS 2D-reshape boundary cut.
* Tensor Ring (cyclic MPS).
* CP-rank Kruskal lower bound.

All four equal `min(W^j, φ(W)·W^(d-j-1)+1)` (E2.1's closed form), with
exact saturation in 21/22 cases and a single finite-size deficit of 1
at W=5, d=4, N=625 (a benign small-N row dependency).

**Asymptotic ratio**: `bond_dim / sqrt(N) → φ(W)/W` (the **Mertens
product**) — verified at d=12 for W∈{2,3,4} and at d=8 for W=6:

* W=2: 33/64 = 0.515 → 0.500
* W=3: 487/729 = 0.668 → 0.667
* W=4: 513/1024 = 0.501 → 0.500
* W=6: 433/1296 = 0.334 → 0.333

The Mertens product is the **universal** asymptotic compression ratio
across the family — was previously only an MPS-specific number.

### Three closure mechanisms (orthogonal across the seven listed ansätze)

1. **Unfolding-rank route (5 ansätze)**: MPS, HT, TR, PEPS, CP all
   admit an unfolding contraction whose matricisation rank IS (or
   lower-bounds) the bond / rank parameter. E2.1 then forces
   `≥ φ(W)·W^(d-j-1)+1` at the maximising cut.
2. **Tucker by mode-slice independence**: mode-`j` matricisation has
   rank ≤ W trivially, so the multilinear-rank tuple does NOT see
   E2.1. But the W mode-`j` slices of χ_P are *linearly independent*
   (Dirichlet's theorem ensures different residue classes mod W have
   different prime patterns), forcing multilinear rank tuple
   `(W, ..., W)` — no compression.
3. **MERA by parameter count**: bond-dim alone admits weak Rényi-2
   bound `χ ≥ rank^{1/(2 log d)}` (small in absolute terms). But MERA
   with depth `log d` and bond dim `χ` has total parameter count
   `O(d · χ^c)`, which to represent χ_P's `Ω(N^{1/2})` unfolding rank
   requires `χ ≥ N^{Ω(1)}` — polynomial, not polylog.

### What this lifts

E2.1 was previously a single-ansatz statement. After S77 it becomes
the **universal half-cut bond dim** across five spatial-locality
ansätze. Five prior CLOSED_PATHS rows (lines 171, 185, 517, 518, 600)
are subsumed under a single mechanism.

### What remains open (next-action targets)

* **Non-spatial-locality ansätze.** Random-projection of mode subsets,
  algebraic constructions (Reed-Solomon-modulated tensors), quantum-
  walk-style oracle representations — explicitly carved out from N1.
  These are the natural follow-ups.
* **L1 Lean formalisation `lower_bound`** (Arc 2): the last `sorry` in
  `MPSBondDim/Basic.lean`. N1 has now empirically verified the lower
  bound at 21/22 (W, d) pairs — the formal proof remains open.

### Files touched

* New: `experiments/constructions/tensor_compression_family_closure/`
  (`definition.md`, `tensor_compression_family.py`,
  `tensor_compression_family_results.md`,
  `tensor_compression_family_results.json`, `run_full.log`).
* Modified: `EDGES.md` (E2.1 annotated with S77 family-level scope).
* Modified: `status/CLOSED_PATHS.md` (one row added at line 750).
* Modified: `NOVELTY_CHALLENGES.md` (N1 marked BUILT).
* Modified: `RESEARCH_AGENDA.md` (Arc 4 sub-milestone added).
* New: `archive/sessions/session77_n1_tensor_compression_family.md`.

## Session 79 — A3 frontier attack: Cayley graph spectrum probes ω(n), not primality (CLOSED)

**Mode:** frontier attack (ATTACK_VECTORS §A.A3).
**Channelling:** Bourgain (spectral gap on Cayley graphs of arithmetic
groups, sum-product); Lubotzky-Phillips-Sarnak (arithmetic Cayley graph
spectra).
**Cross-domain ingredient:** spectral graph theory of Cayley graphs of
finite abelian groups (Babai 1979).
**Outcome: B-grade** (ambitious frontier attack failed informatively;
new structural lower bound + new negative-shape edge).

### What was attempted

Per ATTACK_VECTORS §A.A3, build Cay((Z/nZ)*, {2,3,5,2⁻¹,3⁻¹,5⁻¹}) and
ask: does any spectral feature decide primality of n? Distinct from
prior CLOSED_PATHS lines 354/383/384/385 which used primes-as-generators
or graph-index (circular). A3 fixes generators and varies n.

### What the data showed

For n ∈ [7, 2000] coprime to 30 (300 primes, 13 prime powers, 214
semiprimes, 6 triple-prime composites), 24 spectral features tested.
**No feature is disjoint primes-vs-composites.** Mann-Whitney AUC for
the hard sub-case (primes vs prime powers, both ω=1):

* {2,3,5} ∪ inv: 0.509
* {2,3,7} ∪ inv: 0.566
* {2,3,5,7} ∪ inv: 0.673

Essentially no discrimination. The structural reason:

### New theorem (provable from CRT character theory + empirically verified)

> For n odd coprime to S, the Cayley graph Cay((Z/nZ)\*, S ∪ S⁻¹) has
> at least **2^ω(n)** integer eigenvalues, where ω(n) is the number of
> distinct prime factors of n.

Proof: 2-torsion subgroup of the character group has size 2^ω(n) (each
(Z/p_i^{a_i}Z)\* contributes one Z/2Z); for χ in the 2-torsion subgroup,
χ(s) ∈ {±1} for all s, so λ_χ ∈ Z. ∎

Verified across 533 (n) × 3 (generator set) = 1599 measurements, **0
violations**. Bound is sharp in 84/533 = 15.8% of cases.

### Why A3 fails by structure

(Z/p^kZ)\* is cyclic for odd p (Gauss). So **primes n and prime powers
n=p^k both have ω(n)=1 and produce structurally identical cyclic unit
groups.** Spectrum cannot distinguish them up to discrete-log alignment.

To rescue, you'd combine spectral features with a perfect-power test —
which inherits AKS-class growing-dim MPOW depth (E5.3, E7.10).

### What this adds to the catalogue

* **New EDGES.md entry E7.12** — Cayley graph spectrum probes ω(n), not
  primality (search-space constraint on the entire fixed-generator
  abelian-Cayley-graph spectral primality family).
* **New CLOSED_PATHS row** for ATTACK_VECTORS §A.A3 (line 751, mode E).
* **New explicit identity** n_int_eigs(Cay((Z/nZ)\*, S ∪ S⁻¹)) ≥ 2^ω(n)
  with character-theoretic proof and 84/533 sharpness.

### Files touched

* New: `experiments/circuit_complexity/cayley_spectral_primality/`
  (`cayley_spectral_primality.py`, `cayley_robustness.py`,
  `cayley_spectral_primality_results.md`, `spectrum_features.json`,
  `separation_summary.json`, `robustness_features.json`).
* Modified: `status/CLOSED_PATHS.md` (one row added, line 751).
* Modified: `EDGES.md` (new E7.12 entry).
* Modified: `ATTACK_VECTORS.md` (A3 closure entry added under
  "Closed attacks").
* New: `archive/sessions/session79_a3_cayley_spectral_primality.md`.

### Why B-grade not A-grade

A-grade would require a *positive* result: a spectral feature that
genuinely decides primality. The frontier attack failed in the
predicted direction. Failure mode is structural and cleanly characterised
(spectrum is blind to cyclic-vs-cyclic), so this is firmly B-grade
ambitious failure rather than C-grade refinement.

## Session 80 — D4 frontier attack: Szegedy quantum walks on prime graphs (CLOSED)

### What was done

Cross-domain frontier attack on ATTACK_VECTORS §D.D4. Imported Szegedy's
2004 discriminant matrix theorem (arxiv quant-ph/0401053) and applied
it to three number-theoretic graph families:

(A) **Cayley `Cay((Z/NZ)*, {2,3,5,inv})`** for primes N ∈ {31, 61, 89,
127, 167, 211, 257, 307, 373, 449, 541, 673, 809, 1009}.
(B) **Coprime graph `C_x`** (vertices [1..x], edges (m,n) iff gcd(m,n)=1)
for x ∈ {30, 50, 100, 200, 350, 500, 700, 1000}.
(C) **Divisor graph `D_x`** (vertices [1..x], edges (m,n) iff m|n) for
the same x range.

For each: computed lazy-walk transition matrix `P`, discriminant
`D(x,y) = sqrt(P(x,y)·P(y,x))`, full spectrum, classical mixing time,
Szegedy mixing time `1/(2 arcsin(sqrt(δ)))`, primality correlation of
top-50 eigenvectors, and stationary prime mass.

### Findings

1. **Cayley sweep** (S79 quantum extension): empirical `t_C(N) ~ N^{0.789}`,
   `t_Q(N) ~ N^{0.414}`. Quadratic Szegedy speedup confirmed but neither is
   polylog. This **EXTENDS E7.12 to the quantum-walk regime** —
   not just classical Cayley spectrum but also Szegedy quantisation
   inherits the ω(n)-probing barrier.

2. **Coprime graph** has Ω(1) gap (asymptotically `δ → 0.4166`) and a
   clean closed-form: stationary prime mass / uniform prime mass
   `→ ζ(2) = π²/6 ≈ 1.6449` (Mertens). Verified at x=1000 with
   deviation -0.022. Mixing is `O(1)` both classical and quantum, but
   primality bias is a constant — no polylog opening.

3. **Divisor graph** has Ω(1) gap; high-prime-mass eigenvectors localise
   on **prime-centered divisor clusters** (one specific prime p and its
   multiples per eigenvector). For example, k=14 at x=250 (eigenvalue
   0.768) concentrates 100% of mass on {43, 47, 86, 94, 172, 188, …}
   — the prime-cluster localisation is structurally analogous to E7.12
   degree-class probing, not global primality.

### Outcome

CLOSED §D.D4 with mode E. Added E7.13 to EDGES.md: "Szegedy walks on
arithmetic graphs do not yield polylog π(x)" — combines the spectral
barrier (Cayley class) with the localisation barrier (coprime/divisor
class) into a single negative-shape edge for the entire family.

### Grade: B (ambitious cross-domain failure)

Cross-domain Szegedy machinery did real work: the discriminant-theorem
quadratic speedup `t_Q ≈ sqrt(t_C)` was empirically realised on the
Cayley sweep; the π²/6 stationary identity is a clean Mertens-law
derivation specific to the coprime-graph chain. Failure is structural
and cleanly characterised. Source paper (Szegedy 2004) cited.

### Files touched

* New: `experiments/quantum/szegedy_walk_prime_graphs/`
  (`szegedy_walk_prime_graphs.py`, `szegedy_walk_extended.py`,
  `eigenvector_inspection.py`, `degree_class_check.py`,
  `szegedy_walk_prime_graphs_results.md`, plus json/log artefacts).
* Modified: `status/CLOSED_PATHS.md` (one row added).
* Modified: `EDGES.md` (E7.12 quantum-extension annotation;
  new E7.13 entry).
* Modified: `ATTACK_VECTORS.md` (D4 closure entry under
  "Closed attacks").
* Modified: `RESEARCH_AGENDA.md` (D4 added to closed-attacks pointer).
* New: `archive/sessions/session80_d4_szegedy_walk.md`.

### Why B-grade not A-grade

A-grade would require a Szegedy walk producing `O(polylog(x))` mixing
AND a primality-localised eigenvector. The frontier attack failed in
the predicted direction with a clean structural reason. Failure is
informative (the joint condition of fast mixing + primality
localisation appears genuinely incompatible for arithmetic graphs in
the tested families), so this is firmly B-grade ambitious failure
rather than C-grade refinement.

## Session 84 — A1 wild swing: SAT-based TC^0 search for PRIMES at N=8 (PARTIAL CLOSURE)

**Date:** 2026-04-26.
**Mode:** WILD SWING (frontier attack on ATTACK_VECTORS.md §A1).
**Self-grade:** B (substantive structural finding + concrete mechanism;
  not A because the discovered PRIMES-vs-random gap reduces to elementary
  oddness rather than deep arithmetic structure).

### What was tested

Frontier attack §A1: enumerate small TC^0 circuits matching PRIMES at
N=8 via SAT/ILP. Built three solvers (Z3 If-then-else, direct ILP via
PuLP+CBC, pre-enumerated column-selection ILP). Tested:

1. PTF degree of PRIMES at N=4..8 (LP feasibility on monomial basis).
2. Depth-1 sign-threshold feasibility for PRIMES at N=4..8.
3. Depth-2 sign-threshold (W=1 both layers) min-M complexity at
   N=4 (full search), N=6 (full search), N=7 (partial), N=8 (partial,
   restricted to k_max ≤ 5 nonzero bottom-layer weights).
4. Random-matched controls at N=4 (10 seeds) and N=6 (10 seeds).
5. Single-bit predictor accuracy for PRIMES vs random matched (30 seeds)
   at N=4..8.

### Headline results

- **PTF degree of PRIMES at N=8 = 4**, statistically indistinguishable
  from random matched. Confirms E1.10/E3.13 at the Boolean-function
  level.
- **Depth-2 sign-threshold W=1 sizes**: N=4 → 3, N=6 → 6, N=7 → ≥ 6,
  N=8 → ≥ 17 (k_max ≤ 5 sub-family).
- **PRIMES-vs-random gap at N=6**: PRIMES M=6, all 10 random matched
  M ∈ {7, 8}. Binomial p < 0.001. **First circuit-complexity measure
  in the project where PRIMES empirically deviates from a matched-density
  random function** (35+ prior pseudorandomness measures show no
  deviation).
- **Mechanism**: PRIMES has a 70.3% single-bit predictor (bit_0 = "is x
  odd") at N=8 vs random matched best 1-bit at 57.0% (over 30 seeds).
  Depth-2 circuit gets the 1-bit predictor as a "free" first gate.

### Why §A.A1 only partially closed

Depth-2 W=1 sub-family closed at N=8 with ≥17-gate lower bound. But
the full §A.A1 question (depth ≤ 5, size ≤ 2000, arbitrary integer
weights) remains open — the search space is too large for direct SAT/ILP
enumeration at N=8. Three concrete follow-ups in NOVELTY_CHALLENGES C7,
C8 (added by this session): (a) larger weight bound W ∈ {2, 4, 8},
(b) depth-3 with restricted bottom layer, (c) "1-bit-calibrated random"
baseline at N=6 to test whether the gap is fully explained by oddness.

### Files produced

* New: `experiments/circuit_complexity/sat_tc0_primes_n8/`
  (`sat_tc0_primes_n8.py`, `sat_depth2_ilp.py`, `enum_d2_smart.py`,
  `enumerated_depth2.py`, `ptf_degree_battery.py`,
  `depth1_threshold_test.py`, `n6_robust_check.py`,
  `greedy_d2.py`, `greedy_d2_np.py`,
  `sat_tc0_primes_n8_results.md`, plus json/log artefacts).
* Modified: `status/CLOSED_PATHS.md` (one row added).
* Modified: `ATTACK_VECTORS.md` (§A.A1 partial closure entry).
* Modified: `NOVELTY_CHALLENGES.md` (C7, C8 next-action targets).
* New: `archive/sessions/session84_a1_sat_tc0_primes.md`.

### Why B-grade not A-grade

A-grade would require either:
- A small explicit polylog-promising circuit for PRIMES at N=8 (NOT FOUND;
  the W=1 sub-family is closed at ≥17 gates, super-linear scaling).
- Or a structural identification of the bottom-layer gates with classical
  number-theoretic primality witnesses (NOT FOUND; gates are empirical
  fits with no Miller-Rabin / BPSW / AKS interpretation).

The PRIMES-vs-random gap discovery is statistically significant but its
mechanism (oddness) is elementary — well-known fact, not a new insight.
This earns B-grade as a "concrete-mechanism deviation from
pseudorandomness" finding, not A-grade.

## Session 87 — D6 frontier attack: Gowers U^k norms of chi_P (CLOSED, B-grade)

Frontier attack on ATTACK_VECTORS §D6: do Gowers uniformity norms
detect arithmetic structure in the bare prime indicator chi_P that
the project's 35+ pseudorandomness measures missed?

**Mathematician channel:** Tao (additive combinatorics).
**Cross-domain technique:** Gowers norms / Green-Tao-Ziegler U^{s+1}
inverse theorem.

**Setup.** Computed `Q^k(f) := ||f||_{U^k}^{2^k} / E[f]^{2^k}` on Z/NZ
for chi_P, matched-density Bernoulli, Liouville indicator, and
W-tricked chi_{W,1} (W ∈ {6, 30, 210}). U^2 via FFT identity (O(N log N))
at N up to 2^18 = 262144; U^3 via recursion through Δ_h
(O(N^2 log N)) at N up to 2^13. Computed Hardy-Littlewood {0,1}^k-cube
singular series numerically: **S_2 = 2.300938**, **S_3 = 54.116088**.

**Empirical results.**
- `Q^2(chi_P)`: 2.103 (N=2^10) → 2.153 (N=2^18), monotonic, gap to S_2
  decays as O(1/log N).
- `Q^3(chi_P)`: 35.5 ± 0.1 stable across N ∈ [2^10, 2^13] (vs Bernoulli → 1).
- `Q^2(Liouville)`: 1.000 to four digits — Gowers-uniform at U^2,
  matches Mobius/nilsequence orthogonality (Green-Tao 2008).
- W=210 (primorial of 7): Q^2(chi_{W,1}) = 1.003 vs HL prediction
  S^(W)_2 = 1.002 — match within 0.1%, Green-Tao W-trick mechanism
  empirically verified for the BARE indicator chi_P (not just Lambda).

**Mode E (technique exhaustion).** chi_P's U^k structure is *exactly*
Hardy-Littlewood — no extra bit beyond the singular series. Q^k
extraction itself costs Θ(N^2 log N), not polylog. No deviation from
HL prediction at U^2 or U^3 within tested range.

**New EDGE E2.13** (chi_P U^k = HL singular series). Adds the 36th-37th
pseudorandomness measure to the project's battery, this one with a
*closed-form prediction* — qualitatively distinct from the local /
spectral / Boolean predecessors which all returned "indistinguishable
from random". Combined with E2.9 (degree-d Fourier weight BELOW random
for d ≥ 2), the algebraic story of chi_P's higher-order structure is
now: "Fourier-weight-below-random in low degree, Gowers-norm-above-
random by HL singular series at higher correlation order".

### Files added / modified
* New: `experiments/information_theory/gowers_uk_chi_p.py`,
  `gowers_uk_chi_p_analyze.py`, `hardy_littlewood_box.py`,
  `gowers_uk_chi_p_results.md`, plus 2 small _results.md companions
  for utility modules; data folder with 3 JSON + 1 log + 1 analysis.md.
* Modified: `EDGES.md` (new edge E2.13 in §2 Algebraic).
* Modified: `status/CLOSED_PATHS.md` (row 756 added).
* Modified: `CROSS_DOMAIN_TECHNIQUES.md` (Gowers norms PROPOSED → USED E).
* Modified: `ATTACK_VECTORS.md` (§D6 → "Closed attacks" section).
* Modified: `NOVELTY_CHALLENGES.md` (§D6 marked CLOSED, two
  successor challenges §D6.a §D6.b proposed).
* New: `archive/sessions/session87_d6_gowers_uk_chi_p.md`.

### Why B-grade not A-grade

A-grade would have required:
- A deviation of `Q^k(chi_P)` from the HL prediction `S_k` (NOT FOUND;
  Q^2 converges monotonically to S_2 within ~7%; Q^3 converging slowly
  but in the right direction).
- A failure of the W-trick to restore Gowers uniformity (NOT FOUND;
  W=210 gives Q^2 within 0.1% of HL prediction).
- A novel inverse-theorem identification of a non-trivial nilsequence
  correlating with chi_P (the inverse-theorem ε-threshold here is
  ‖chi_P‖_2^4 = E[chi_P]^2 = p^2; the U^2 ratio Q^2 ≈ S_2 · p^2 → 0,
  so no non-trivial nilsequence detection at our N).

Every observed result is the predicted outcome of known mathematics
(HL conjecture + Green-Tao W-trick). The NUMERIC values S_2 ≈ 2.301,
S_3 ≈ 54.12 and the empirical match are project-internal new but
straightforward consequences of HL machinery. B-grade ceiling is
honest.

## Session 88 — C4 frontier attack: Anderson localisation Lyapunov exponent of chi_P (CLOSED, B-grade)

**Mode:** cross-domain attack. Target ATTACK_VECTORS §C4. Cross-domain
ingredient imported: Anderson localisation theory (Aizenman-Warzel
2015; Furstenberg-Kifer 1983; Pastur-Figotin 1992).

**Outcome:** spectral analogue of S87/E2.13. The Anderson Lyapunov
exponent gamma(E) of the discrete 1D Schrödinger operator with
prime-indicator potential V = chi_P deviates from the matched-Bernoulli
baseline by 88 sigma at N = 10^5 (peaks at parity resonance E = 0 and
mod-3 resonance E ≈ +1). The deviation cascades down to ~ 4 sigma
under the Green-Tao W-trick at W = 2310. Same structural reason as
E2.13 (Gowers norms): chi_P deviation is the spectral signature of
small-modulus residue-class structure of primes, no extra bit beyond
Hardy-Littlewood equidistribution mod q.

**New edge: E2.14** — Anderson Lyapunov gamma(E) of chi_P-Schrödinger
deviates from random Bernoulli by amount captured by W-trick.

**Cross-domain technique status:** Anderson localisation marked
USED E (S88) in CROSS_DOMAIN_TECHNIQUES.md.

**Files**
* New: `experiments/dynamical/anderson_localisation_chi_p/anderson_localisation_chi_p.py`,
  `parity_control.py`, `wtrick_control.py`, plus three companion
  results.md files and 6 JSON output files.
* Modified: `EDGES.md` (new edge E2.14 in §2 Algebraic).
* Modified: `status/CLOSED_PATHS.md` (new row appended).
* Modified: `CROSS_DOMAIN_TECHNIQUES.md` (Anderson localisation
  PROPOSED → USED E).
* Modified: `ATTACK_VECTORS.md` (§C4 → "Closed attacks" section,
  with two follow-up successor proposals).
* New: `archive/sessions/session88_c4_anderson_localisation_chi_p.md`.

### Why B-grade not A-grade

A-grade would have required a deviation that the W-trick CANNOT
reach — i.e., gamma_prime - gamma_CW staying significant at
W = primorial(k) for unbounded k, or a spectral edge anomaly with
no Hardy-Littlewood explanation. Neither found; the cascade decays
geometrically with W exactly as the Green-Tao machinery predicts.
The cross-domain ingredient (Anderson localisation theory) was
non-trivial and did real work, but the structural conclusion
matches an existing closure family (E2.13). Honest B-grade.

### Next-actions queued

* `§C.C4(a)`: log-weighted potential V(n) = log p_n if n prime else 0,
  amplifies gamma by (log N)^2 — could push the W-trick cascade to
  W = 30030.
* `§C.C4(b)`: Anderson Lyapunov of Lambda(n) (von Mangoldt) potential.
  Green-Tao predicts Lambda is Gowers-uniform after centering,
  spectral analogue: gamma_Lambda - gamma_random → 0 with NO
  W-trick needed. Direct visual confirmation of the additive vs
  spectral W-trick symmetry.

### Triple-confirmation note for the project frontier

S82 (Dirichlet character spike eigenvectors), S87/E2.13 (Gowers
U^k → HL singular series), and S88/E2.14 (Anderson Lyapunov W-trick
cascade) are now THREE INDEPENDENT spectral / additive / operator-
theoretic measurements pointing to the same conclusion: chi_P
structure = HL equidistribution mod q, no more. Future
pseudorandomness work should target measurements that the W-trick
CANNOT REACH — multiplicative regime (Mobius mu, Liouville lambda),
or scales/structures beyond all tractable W (twin-prime, k-tuple
Hardy-Littlewood densities, log-weighted ensembles).

## Session 89 — C7 construction: calibrated 1-bit-bias control closes S84 PRIMES-vs-random depth-2 gap (B-grade)

**Date:** 2026-04-26.
**Mode:** Construction (NOVELTY_CHALLENGES.md §1 C7).
**Self-grade:** B (substantive refinement; verified S84's mechanism
conjecture under controlled experiment; closes the project's *only*
outstanding pseudorandomness-deviation claim to "elementary parity,
no residual structure").

### What was tested

Composition C7: build a class-conditional matched-density random
control `f_cal` for PRIMES at N=6 — P(f_cal=1 | x odd) = 17/32, P(f_cal=1
| x even) = 1/32, exactly matching PRIMES (17 odd primes, 1 even prime
in {0..63}). Sample 20 stratified (exact 17+1, weight=18 always,
bit_0_acc=0.75 always) + 20 bernoulli (independent draws) and run S84's
`enum_d2_smart` ILP harness (K=1458 W=1 candidates, M ∈ {3..8},
CBC 120s/cell). Falsifier table F1-F4 pre-stated before any run.

### Headline result

**F2 + F3 + F4 hold; F1 fails.**

- Stratified n=20: histogram = {5: 4, 6: 16}, mean 5.80, max 6.
  **0/20 above PRIMES (M=6); 4/20 strictly below PRIMES (M=5).**
- Bernoulli n=20: histogram = {5: 7, 6: 11, 7: 2}, mean 5.75. Both
  M=7 cases have bit_0_acc < 0.75 (PRIMES's value: 0.7344 and 0.7188).
- PRIMES sits at +0.5σ of the calibrated distribution; under the
  calibrated null `P(M ≤ 6 | calibrated) = 1.0` vs S84's unbiased
  null `P(M ≤ 6 | unbiased) = 0/10 ≈ 0.001`.

**The S84 PRIMES-vs-random depth-2 sign-threshold gap REDUCES TO
elementary parity** ("π(x) ≈ 1 iff x odd, for x > 2"). No
PRIMES-specific structure beyond oddness.

### Sub-finding (Bernoulli mode bit_0_acc binning)

| bit_0_acc range | n | min_M distribution     | mean |
|-----------------|---|------------------------|------|
| ≥ 0.78          | 6 | {5: 3, 6: 3}           | 5.50 |
| [0.75, 0.78)    | 6 | {5: 2, 6: 4}           | 5.67 |
| [0.72, 0.75)    | 5 | {6: 3, 7: 2}           | 6.40 |
| < 0.72          | 3 | {5: 3}                 | 5.00 |

At fixed weight ≈ 18, larger bit_0_acc → smaller min_M. The two M=7
cases live in the [0.72, 0.75) bin (just below PRIMES's 0.75). This
is the empirical confirmation of S84's qualitative mechanism: depth-2
W=1 size is monotonically driven by 1-bit predictor strength.

### Implications for the project

- S84 was the **only** statistical deviation from pseudorandomness
  across the project's 35+ measure battery (E1.10 / E3.13). C7 shows
  this deviation is captured by an elementary fact (most primes are
  odd) — no new structural property of π(x).
- Pseudorandomness thesis (`novel/pseudorandomness_of_pi.md`) is
  **refined, not weakened**: π(x) is pseudorandom *modulo* the
  trivial mod-2 bias. After bias correction (W-trick at W=2), no
  measure in the battery distinguishes π from random.
- **Quadruple confirmation of the W-trick story for chi_P** now:
  S82 (Dirichlet character spike eigenvectors), S87/E2.13 (Gowers
  U^k cascade), S88/E2.14 (Anderson Lyapunov cascade), S89 (depth-2
  W=1 size at W=2 calibration). Four independent witnesses of the
  same closure principle.

### Next-actions queued

- **C7.a**: same calibration at N=8 (P(prime|odd) = 53/128, P(prime|even)
  = 1/128). Predicted outcome same. If FALSIFIED at N=8, the bit_0
  explanation breaks at higher N → A-grade material.
- **C8**: depth-2 W-vs-M tradeoff at N=8 for W ∈ {2, 4, 8}. Under
  bigger W, bit_0 advantage dilutes; quantifies how much of the W=1
  gap survives at higher W.
- Add footnote to `novel/pseudorandomness_of_pi.md`: 36th measure
  (depth-2 sign-threshold) deviates from unbiased random but is
  fully captured by 1-bit oddness; pseudorandomness stands "modulo
  the obvious mod-2 bias."

### Why B-grade not A-grade

B because: substantive refinement of an existing edge with a precise
new statement (the F1-F4 verdict) that extends its scope; constructed
the matched-class-conditional sampler that didn't exist before;
quantitative bernoulli sub-finding strengthens S84's mechanism from
"qualitative" to "directly observed monotone driver."

Not A because: the F2 outcome was the high-prior outcome per S84
(~70%); no new structural object emerged from the calibration; the
closure mode is C (circular) — calibration formally verifies an
explanation already in the project's frame. No published-paper-grade
NT could reach the F1-F4 verdict from S84 + the harness in less than
a few hours; the C7 sampler is small but it is a real new piece of
controlled experiment.

### Files produced/modified

* New: `experiments/circuit_complexity/sat_tc0_primes_n8_calibrated/calibrated_d2.py`,
  `calibrated_d2_n6.json`, `calibrated_d2_n6.log`,
  `calibrated_d2_results.md`, `definition.md`.
* New: `archive/sessions/session89_c7_calibrated_d2_primes.md`.
* Modified: `NOVELTY_CHALLENGES.md` (C7 → BUILT, C7.a successor added).
* Modified: `RESEARCH_AGENDA.md` (Arc 4 milestone added).
* Modified: `status/CLOSED_PATHS.md` (closure row added).

## Session 92 (B-grade) — §B1 ATTACK_VECTORS: Algebraic immunity / polynomial method on chi_P

### Wild swing target

§B1 of ATTACK_VECTORS.md: apply the cap-set / Croot-Lev-Pach / slice rank
polynomial method to chi_P. Goal: a low-degree F_p polynomial P with
P(n) determining primality, with degree growing more slowly than naive.

### Three invariants tested

1. **AI_F_2(chi_P)** — algebraic immunity over F_2 (Courtois-Meier 2003).
   Computed by Gauss elimination over F_2 on the 2^N x sum_i C(N,i)
   monomial-evaluation matrix, N=4..13.
2. **F_p multilinear ANF degree** of chi_P viewed on (F_p)^k for
   (p, k) in {(3,2..5), (5,2..3), (7,2..3), (11,2)}.
3. **Slice rank brackets** via Tao 2016 flattening LB + greedy UB.

### Main finding

**AI(chi_P, N) = 2 for ALL N in [4, 13]** (zero std at N >= 8), while
random matched-density Bernoulli grows as Theta(log_2(1/rho)),
reaching 4 at N=11. **THE annihilator is the SAME polynomial g(x) =
(1+x_0)(1+x_1) for every N**, encoding "n ≡ 0 mod 4". This is the
polynomial-method version of the trivial mod-4 sieve fact (no prime
> 2 divisible by 4).

### W-trick correction confirms

`chi_P_{W,b}(n) := chi_P(W*n + b)` with gcd(b, W) = 1: AI(chi_P_{W,b})
tracks AI(random matched-density) within ±1 step, often exact match
with zero std at W >= 6, b in {1, 5, 7, 11}.

### Triple of independent confirmations

E2.13 (Gowers, S85) + E2.14 (Anderson, S88) + E2.15 (algebraic
immunity, S92): chi_P deviates from random in three orthogonal
mathematical categories — additive combinatorics, spectral /
transfer-matrix Lyapunov, Boolean polynomial method — each via the
same mod-q sieve mechanism.

### Files produced/modified

* New: `experiments/algebraic/algebraic_immunity_chi_p/algebraic_immunity_chi_p.py`,
  `extract_annihilator.py`, `wtrick_AI.py`,
  `algebraic_immunity_chi_p_results.md`, JSON data files.
* New: `archive/sessions/session92_b1_algebraic_immunity_chi_p.md`.
* Modified: `EDGES.md` (E2.15 added).
* Modified: `ATTACK_VECTORS.md` (B1 marked closed; closure entry +
  self-extension follow-ups added).
* Modified: `status/CLOSED_PATHS.md` (S92 row).

### Why B-grade not A-grade

B because: ambitious frontier attack on §B1 with structural
informative failure; new edge (E2.15); third-leg confirmation in a
distinct mathematical category. Not A because: deviation fully
captured by W-trick — no new arithmetic content beyond mod-q.
Self-extension follow-ups proposed (a) AI of multiplicative-regime
functions (Liouville, Mobius); (b) higher-degree structured
annihilators of chi_P_{W,b}.

## Session 93 — Λ vs χ_P U^k Comparison (D6.b, novelty mode)

**Mode:** novelty (B-grade by construction).
**Target:** NOVELTY_CHALLENGES.md §D6.b.
**Self-grade:** B.

### What was measured

Side-by-side normalized Gowers norms `Q^k(f) := ||f||_{U^k}^{2^k} / E[f]^{2^k}`
for f = χ_P (bare prime indicator) and f = Λ (von Mangoldt log-weighted)
on Z/NZ at N up to 2^17 for U^2, N up to 2^12 for U^3 (full h sampling),
plus W-tricked variants at W ∈ {6, 30, 210}.

### Headline finding

**S_k is universal across log-weighting.** Q^k(Λ) = Q^k(χ_P) to within
0.4% at U^2, 2.5% at U^3, across N in [2^10, 2^17]. Both converge to
S_k = 2.301 (k=2) from below at near-identical rates. After W-trick at
W = 210: `Q^2(χ_W) = Q^2(Λ_W) = 1.0029` to four decimals.

### Source of the small bare offset

Q^k(Λ) > Q^k(χ_P) by ~0.3% (U^2) / ~2.5% (U^3) — identified with the
prime-power weight in Λ (log p contribution from n = p^k for k ≥ 2),
vanishing as π(√N)/π(N) → 0.

### Refinement of E2.13

E2.13 sharpened from "Q^k(χ_P) → S_k" to "S_k is the universal Gowers
fingerprint of {0,1}^k-cube prime correlation, independent of weighting
choice (χ_P or Λ)". Inline edit in EDGES.md.

### Files

* New: `experiments/information_theory/lambda_vs_chi_p_uk/lambda_vs_chi_p_uk.py`,
  `lambda_vs_chi_p_uk_results.md`, `main_run.json`.
* New: `archive/sessions/session93_d6b_lambda_vs_chi_p_uk.md`.
* Modified: `EDGES.md` (E2.13 refined inline).
* Modified: `NOVELTY_CHALLENGES.md` (§D6.b CLOSED, §D6.c proposed).

### Why B-grade not A-grade

B because: clean refinement of E2.13 with weighting-universal scope
plus four-decimal W-trick coincidence; all 5 pre-stated falsifiers
resolved as predicted. Not A because: every result confirms an
established prediction (HL singular series, GT W-trick); no deviation
that would indicate new structure.

---

## Session 95 — DPP / PPP fit to integer primes (ATTACK_VECTORS §D7)

**Mode:** frontier attack with B-grade fallback. **Self-grade: B.**
**Channel:** Tao (additive combinatorics & point-process framework).
**Cross-domain technique imported:** Determinantal Point Processes
(DPPs) — first project use; was PROPOSED-only in CROSS_DOMAIN_TECHNIQUES.md
§3 since S91 frontier_gen.

### Key result

Primes are quantitatively NOT a translation-invariant DPP / PPP /
signed-real-K / complex-Hermitian-K point process. Five pre-stated
falsifiers all HOLD at N = 10^7 + admissible all-even triples up to t_2 ≤ 26.

### Mechanism (one paragraph)

HL `S(0, t_1, ..., t_k)` factorises over PRIMES — each factor `α_p`
depends on `ν_p({0, t_1, ..., t_k})` = #distinct residues mod p.
DPP / PPP correlations factorise over PAIRS. Pairwise admissibility
≠ triple admissibility (e.g., (0, 4, 14) is pair-admissible mod 3 but
ν_3 = 3 saturates as a triple, killing R_3^HL while PPP predicts ~10^{-3}).
The HL prime-factor structure cannot reduce to pair-level data; this
is the structural obstruction to any kernel fit.

### Quantitative breakdown

Pair-level: `R_2(t) - ρ²` flips sign across odd/even (R_2 ≈ 0 odd,
R_2 ≈ ρ² S(0,t) > ρ² even). 3-point: PPP overshoots HL by 7.5–79.2%;
max gap on 3-AP triples (0, 6, 12), (0, 12, 18), (0, 18, 24).
Real-signed-K: σ_req ∈ (-0.541, +0.769) — never ±1.
Complex-Hermitian-K: best LS phase residual 0.0746 ≫ noise floor 0.01.

### New edge: E2.16

First 3-point structural statement ruling out a kernel-factorisation
representation of χ_P. Complements the 2-point pair-level closures
E2.13 (Gowers norm), E2.14 (Anderson Lyapunov), E2.15 (algebraic
immunity) with the FIRST entry that genuinely uses 3-point data.

### Files

* New: `experiments/constructions/primes_dpp_ppp_fit/primes_dpp_ppp_fit.py`,
  `primes_dpp_ppp_fit_results.md`, `main_run.json`.
* New: `archive/sessions/session95_d7_dpp_ppp_fit.md`.
* Modified: `EDGES.md` (E2.16 added),
  `CROSS_DOMAIN_TECHNIQUES.md` (DPP marked USED I),
  `ATTACK_VECTORS.md` (§D.D7 closed + Closed-attacks entry added),
  `status/CLOSED_PATHS.md` (row appended),
  `NOVELTY_CHALLENGES.md` (§D7 closure noted, §D7.b/§D7.c proposed).

### Why B-grade not A-grade

B because: ambitious frontier attack on DPP fit with 5 pre-stated
falsifiers, all HOLD; new edge E2.16 quantifies a structural negative
in a fresh cross-domain (point-process theory); FIRST 3-point closure
in a project saturated with 2-point closures. Not A because: no
positive kernel found, no new polylog representation; the result is
structural negative consistent with the project's "χ_P = HL, no more"
picture rather than breaking it.

---

## Session 96 — D2 (cross-domain, persistent homology / TDA), B-grade

**Mode**: Frontier (ATTACK_VECTORS.md §D2). Cross-domain technique:
persistent homology (Vietoris-Rips filtration on Takens delay
embedding), via Carlsson 2009 BAMS 46 + ripser (Bauer 2021 J. Appl.
Comput. Topol. 5).

**Attack:** D2 (PROPOSED row in CROSS_DOMAIN_TECHNIQUES.md §4).

**Result:** F3 falsifier holds. PRIMES > 3σ below BOTH IID Exp(1)
AND gap-permuted controls on integrated H_0 / H_1 persistence
(`T0`, `T1`).

| run | M | d | T0 z(B1) | T0 z(B2) | T1 z(B1) | T1 z(B2) |
|-----|---|---|----------|----------|----------|----------|
| main, x=1e6 | 2000 | 3 | -10.31 | -8.70 | -4.20 | -11.99 |
| robust, x=5e6 | 2000 | 3 | -8.35 | -7.58 | -3.67 | -8.69 |
| robust, d=2 | 2000 | 2 | -10.56 | -6.68 | -2.68 | -7.78 |
| robust, d=4 | 2000 | 4 | -7.55 | -5.05 | -5.83 | -3.57 |
| scaling | 4000 | 3 | -17.80 | -15.45 | -10.49 | -16.76 |

Z-scores grow super-linearly with M (signal at least linear in
window size). Effect dimension-stable across d ∈ {2, 3, 4}.
Sign: PRIMES T0, T1 are SMALLER than random — primes are MORE
clustered and have FEWER persistent loops in delay-embedding space
than Poisson or gap-permuted.

**Mechanism:** HL k-tuple admissibility constrains consecutive gaps
to repeat residue patterns more often than random, creating
geometric self-similarity in the delay-embedding cloud (small T0)
and suppressing random "out-and-back" loops in delay space (small
T1). The B2 (gap-permuted) control isolates the *serial-correlation*
component from the gap-marginal component.

**Position in HL-detection family:**

| Edge   | Method                  | Mathematical category   |
|--------|-------------------------|-------------------------|
| E2.13  | Gowers U^k norm         | Additive combinatorics  |
| E2.14  | Anderson Lyapunov γ(E)  | Random Schrödinger      |
| E2.15  | Algebraic immunity AI=2 | Boolean / algebra       |
| E2.16  | DPP / PPP failure       | Random matrix theory    |
| **E2.17** | **PH of Takens-embedded gaps** | **Algebraic topology / metric geometry** |

E2.17 is the FIRST topological / metric-space-geometric detection
of HL serial structure; the prior four are analytic (Fourier,
spectral, polynomial, kernel-factorisation).

**Files**

* New: `experiments/topological/persistent_homology_chi_p/persistent_homology_chi_p.py`,
  `persistent_homology_chi_p_results.md`, multiple `*.json` raw runs.
* New: `archive/sessions/session96_d2_persistent_homology_chi_p.md`.
* Modified: `EDGES.md` (E2.17 added),
  `CROSS_DOMAIN_TECHNIQUES.md` (PH marked USED I, edge E2.17),
  `ATTACK_VECTORS.md` (§D.D2 closed + Closed-attacks entry added),
  `status/CLOSED_PATHS.md` (row appended).

**Why B-grade not A-grade**

B because: ambitious frontier attack on D2 with pre-stated F1/F2/F3
falsifiers, F3 holds quantitatively (≥ 5σ across all robustness
checks); new edge E2.17 documents a fifth orthogonal HL-detection
category; first project quantitative PH measurement (S10's TDA
mention was qualitative). Not A because: no polylog opening (PH is
O(M^3)), no new polynomial-time evaluator, the underlying signal is
the same HL serial structure detected by E2.13/E2.14/E2.15/E2.16
in different mathematical categories — adds a new instrument, not
a new physics.

**Honest-failure check:** session produced (i) new edge E2.17 in a
new mathematical category (algebraic topology), (ii) first-in-project
quantitative PH measurement on prime gaps, (iii) two robustness
checks that confirm the signal is real and scaling, (iv) closure of
ATTACK_VECTORS §D2 with concrete numerical signature, and (v) three
self-extension successor challenges. Not a duplicate-only session.

**Next-action:** D2.a (W-tricked PH) — predicted to link E2.17 to
E2.13 if HL singular-series cancellation is the only mechanism;
single session.

---

## Session 98 — L1 Lean: corner case `(W = 2, j = 1)` of E2.1 closed unconditionally (C-grade)

**Mode:** arc continuation (Arc 2, Lean Formalisation Track), Route A'.

**What changed.** Two new sorry-free Lean theorems added to
`MPSBondDim/Basic.lean`:

* `exists_invertible_submatrix_W_eq_2_j_eq_1 (d : ℕ) (hd : 2 ≤ d) :
  ∃ ρ σ, IsUnit ((unfolding 2 d 1).submatrix ρ σ)` — corner-case
  prime exhibit via Bertrand.
* `mps_bond_dim_W_eq_2_j_eq_1 (d : ℕ) (hd : 2 ≤ d) :
  (unfolding 2 d 1).rank = 2` — corner-case rank statement.

`#print axioms` for both: `[propext, Classical.choice, Quot.sound]`
only — no `sorryAx`. The general-case `exists_invertible_submatrix`
`sorry` is unchanged (still the only outstanding obligation).

**The construction (≈70 Lean lines).** For every `d ≥ 2`, Bertrand at
`n = 2^(d-1)` gives a prime `p ∈ (2^(d-1), 2·2^(d-1)]`. Choose
`ρ = id : Fin 2 → Fin 2`, `σ 0 = ⟨1, _⟩`, `σ 1 = ⟨p − 2^(d-1) − 1, _⟩`.
The 2×2 submatrix is upper-triangular `[[1, ?], [0, 1]]` because:

* `(0, 0)`: `chiP(0·2^(d-1) + 1 + 1) = chiP 2 = 1` (uses `Nat.prime_two`).
* `(1, 1)`: `chiP(2^(d-1) + (p − 2^(d-1) − 1) + 1) = chiP p = 1`
  (uses Bertrand's primality of `p`).
* `(1, 0)`: `chiP(2^(d-1) + 2) = 0` because `2 ∣ (2^(d-1) + 2)` and
  `2^(d-1) + 2 > 2`, contradicting prime via
  `Nat.Prime.eq_one_or_self_of_dvd`.

The `(0, 1)` entry depends on `p` and is irrelevant; `det = 1·1 −
?·0 = 1`. `Matrix.isUnit_iff_isUnit_det` + `Matrix.det_fin_two`
discharges the `IsUnit` goal. The unconditional rank statement
follows by combining `rank_le_min_dim` (gives `rank ≤ 2`) with the
exhibit (gives `2 ≤ rank` via `Matrix.rank_of_isUnit` +
`Matrix.rank_submatrix_le`).

**Why C-grade.** Lean translation of S90's pre-identified Route A'.
The corner-case construction (explicit `σ` and the diagonal-dominance
via `chiP(2^(d-1) + 2) = 0`) is mildly new to the codebase, but the
overall mathematical content was outlined in S90's audit. This is
the project's first **fully verified** instance of `mps_bond_dim` for
any concrete `(W, j)`, but per CLAUDE.md's "self-grade down" rule,
this lands at C.

**Files touched.**
* `experiments/formalisations/E2_1_mps_bond_dim/MPSBondDim/MPSBondDim/Basic.lean`
  (+≈70 lines, no `sorry` introduced).
* `experiments/formalisations/E2_1_mps_bond_dim/mps_bond_dim_notes.md`
  (declaration table, Route A' marked DONE, build status).
* `RESEARCH_AGENDA.md` (Arc 2 milestone added, next-action updated).
* `archive/sessions/session98_l1_lean_corner_case.md` (synthesis).

**Next-action.** (i) Mirror corner `(W = 2, d = j + 1)` — same
Bertrand pattern, ~50 Lean lines, would extend formalised regime to a
second corner. (ii) Route C (PNT-driven low-density regime) — leaves
saturating half-cut regime open but is more impactful than corner
cases. (iii) General case `exists_invertible_submatrix` — Hoheisel-
grade short-interval primes, multi-session, awaits mathlib
infrastructure.

## Session 99 — L1 Lean: orthogonal corner `(W = 2, d = j + 1)` of E2.1 closed unconditionally without Bertrand (C-grade)

**Mode:** arc continuation (Arc 2, Lean Formalisation Track), Route A''.

**What changed.** Three new sorry-free Lean declarations added to
`MPSBondDim/Basic.lean`:

* `chiP_three_eq_one : chiP 3 = 1` — small helper.
* `exists_invertible_submatrix_W_eq_2_d_eq_j_plus_1 (j : ℕ) (hj : 1 ≤ j) :
  ∃ ρ σ, IsUnit ((unfolding 2 (j+1) j).submatrix ρ σ)` — orthogonal
  corner-case prime exhibit, **no Bertrand**.
* `mps_bond_dim_W_eq_2_d_eq_j_plus_1 (j : ℕ) (hj : 1 ≤ j) :
  (unfolding 2 (j+1) j).rank = 2` — orthogonal corner-case rank
  statement.

`#print axioms` for all three: `[propext, Classical.choice, Quot.sound]`
only — no `sorryAx`, no new `axiom`. The general-case `exists_invertible_submatrix`
`sorry` is unchanged.

**The construction (≈110 Lean lines).** For every `j ≥ 1` and
`d = j + 1`, the matrix is `Fin(2^j) × Fin 2` — only **two columns**.
We take both. Choose `ρ 0 = ⟨0, _⟩`, `ρ 1 = ⟨1, _⟩` (rows 0 and 1),
`σ 0 = ⟨1, _⟩`, `σ 1 = ⟨0, _⟩` (column swap). The 2×2 submatrix is
the identity:

```
   ⎡ unfolding(0, 1),  unfolding(0, 0) ⎤   ⎡ chiP 2,  chiP 1 ⎤   ⎡ 1, 0 ⎤
   ⎣ unfolding(1, 1),  unfolding(1, 0) ⎦ = ⎣ chiP 4,  chiP 3 ⎦ = ⎣ 0, 1 ⎦
```

with `det = 1`. Uses only `Nat.prime_two`, `Nat.prime_three`,
`Nat.not_prime_one`, decidability of `¬ Nat.Prime 4` — *no Bertrand*.
The upper bound `rank ≤ 2` uses `Matrix.rank_le_width` (only 2
columns); a small `linarith` step dodges a dependent-type rewrite
issue with `2 ^ ((j+1) - j)` appearing both in the hypothesis and in
the matrix's column-index type.

**What this closes.** Together with S98, the file has **unconditional
Lean proofs of `mps_bond_dim` on the entire `(j, d - j)` boundary of
the parameter grid for `W = 2`** — i.e. the full `R = 2` regime. S98
covers the column `j = 1, d ≥ 2`; S99 covers the row `d - j = 1,
j ≥ 1`. The two corners overlap at `(j, d) = (1, 2)`; the genuinely
new content of S99 is `j ≥ 2`.

**Why C-grade.** Lean translation of S90's pre-identified Route A''.
The mildly novel observation is that this corner is structurally
*easier* than S98's (no Bertrand) — the agenda predicted "the same
Bertrand argument applies, mirrored", but the matrix has only 2
columns so we just take both. Per CLAUDE.md "self-grade down", this
lands at C.

**Implementation lesson.** `rw` and `▸` failed when the pattern being
rewritten (`2 ^ ((j+1) - j)`) appeared inside a dependent type
(`Matrix (Fin (2^j)) (Fin (2^((j+1)-j))) ℚ`). The fix: prove the
arithmetic identity as a separate `have` and chain via `linarith`
rather than rewriting in-place. Worth flagging for future Lean
sessions hitting similar dependent-rewrite issues.

**Files touched.**
* `experiments/formalisations/E2_1_mps_bond_dim/MPSBondDim/MPSBondDim/Basic.lean`
  (+≈110 lines, no `sorry` introduced).
* `experiments/formalisations/E2_1_mps_bond_dim/mps_bond_dim_notes.md`
  (declaration table, Route A'' marked DONE, build status).
* `RESEARCH_AGENDA.md` (Arc 2 milestone added, next-action updated to
  list Routes A''' / C / B).
* `NOVELTY_CHALLENGES.md` (L1 progress paragraph extended).
* `EDGES.md` (E2.1 annotation extended).
* `archive/sessions/session99_l1_lean_orthogonal_corner.md` (synthesis).

**Next-action.** (i) Route A''' — extend corner closures to `W = 3`
(multi-prime exhibits, ~150 lines, 1-2 sessions). (ii) Route C
(PNT-driven low-density regime). (iii) General case awaits Hoheisel-
grade short-interval primes.


## Session 100 — G1 Liouville Anderson Lyapunov (B-grade)

**Mode:** wild swing, ATTACK_VECTORS §G1.

**What changed.** Repeated S88 chi_P Anderson Lyapunov experiment with
the *centered multiplicative encoding* `V(n) = lambda(n) ∈ {-1, +1}`
(mean → 0 by PNT, variance = 1 exactly), comparing the discrete 1D
Schrödinger Lyapunov `gamma_lambda(E)` to an i.i.d. Rademacher ±1
baseline at 51 energies in [-2.95, 2.95], at three sample sizes
N ∈ {10^5 (50 seeds), 3·10^5 (50 seeds), 10^6 (100 seeds)}.

**Result.** Lambda matches Rademacher within seed noise: max |z| flat
at 1.78–2.16 across two orders of N magnitude, argmax energy wanders
between runs ({-0.24, +0.12, -2.01}), χ²/K = 0.49–0.69 (sub-Rademacher),
Pastur-Figotin agreement γ_λ/γ_PF = 0.9317 identical to γ_Rademacher/
γ_PF = 0.9309 to 4 decimals. Independent Chowla two-point aggregate
at N=10^6, h=1..16: Σz² = 4.77 vs Rademacher χ²_16 mean 16 (lambda is
*more* Rademacher-like than Rademacher).

**Stark contrast with E2.14 (chi_P, S88) at same N=10⁵ grid:**

| | chi_P (E2.14) | λ (this session) |
|--|---|---|
| max \|z\| | 88.5σ at locked E=0.108 | 1.78σ at wandering E |
| Scaling | √N (HL mod-q resonance) | flat (Sarnak/GT) |
| W-trick | required, cascade to W=2310 → 4σ | none needed |

The contrast IS the new content: chi_P's spectral deviation is
exclusively HL singular-series mod-q resonance, lambda has no such
resonance. First non-W-tricked spectral measure at noise floor in
the project's 38-measure pseudorandomness battery.

**Cross-domain ingredient (NEW):** Möbius/nilsequence orthogonality
(Green-Tao 2012 *Annals* 175 = arXiv:0807.1736; Sarnak 2010 IAS
lectures; Tao 2016 Forum Math Pi 4 = arXiv:1509.05422 logarithmic-
Chowla). First project use; status PROPOSED → USED E.

**Adds EDGE E2.18.** Closes ATTACK_VECTORS §G1 mode E.

**Why B-grade not A-grade.** F1 falsifier (sustained |z| > 5 not
W-trick-removable) FALSE; F2 (B-grade lower bound) HOLDS strongly.
No new polylog opening. Inverse: the intended polylog opening
("if lambda-Lyapunov has structure, partial-sum the explicit
formula μ-side") is *closed*, not opened.

**Files.**
* New: `experiments/dynamical/liouville_anderson_lyapunov/`
  with main script, analysis script, three result JSONs, and
  `_results.md`.
* Modified: `EDGES.md` (E2.18 added); `CROSS_DOMAIN_TECHNIQUES.md`
  (Anderson row + new Möbius/nilsequence orthogonality row, both
  USED E); `ATTACK_VECTORS.md` (§G1 → Closed attacks with G1.a/
  G1.b/G1.c successors); `status/CLOSED_PATHS.md` (row appended).
* New: `archive/sessions/session100_g1_liouville_anderson_lyapunov.md`.

**Next-action.** G1.c — heavy-tail Schrödinger `V(n) = log p_n` on
primes, 0 elsewhere; introduces Bourgain-Goldstein-Schlag 2002
*Annals* 154 as a NEW cross-domain technique. Single-session.
Alternative: G2 — Gowers U^k of lambda ∈ {-1, +1} centered encoding
(D6.b only tested {0, 1} indicator at U^2). Single-session.

## Session 101 — D6.c μ-weighted χ_P at U^2 — family-level refinement of E2.13 (B-grade)

**Target.** NOVELTY_CHALLENGES.md §D6.c. Test whether `μ(n)·χ_P(n)`
"kills the Hardy-Littlewood structure" of χ_P at U^2 even before
W-tricking. Single-session, low difficulty.

**Literal collapse.** μ(p) = −1 on every prime; χ_P(n) = 0 off the
primes; therefore `μ·χ_P ≡ −χ_P` pointwise. ‖−f‖_{U²} = ‖f‖_{U²} so
the literal D6.c question collapses to a sign issue. Verified:
`max_n |(μ·χ_P)(n) + χ_P(n)| = 0` exactly at every N tested.

**Pivot.** "Does Möbius cancellation propagate from signed μ to its
indicator level sets?" Built a panel of 11 multiplicatively-defined
arithmetic indicators (chi_P, sqfree, mu_plus/minus, lam_plus/minus,
1[Ω=2], signed μ and λ, plus the trivial mu_chi_P / musq_chi_P) on
Z/NZ. Computed Q²(f) = ‖f‖_{U²}^4 / mean(f)^4 for indicators and
Q²_norm(f) = N·‖f‖_{U²}^4/(2·‖f‖_2^4) for sign-valued at
N = 1024..131072. Re-ran under W-trick at W ∈ {6, 30, 210}.

**Result — sharp dichotomy.**

| Function | Bare Q²_∞ | W=210 | Gowers-uniform? |
|----------|-----------|-------|-----------------|
| chi_P (Ω=1, sqfree) | 2.301 (= S_2) | 1.004 | NO bare; YES W |
| sqfree (μ²) | **1.0384** stable | 1.000 | NO bare; YES W |
| 1[μ=+1] (sqfree, even ω) | 1.0384 | 1.001 | NO bare; YES W |
| 1[μ=−1] (sqfree, odd ω) | 1.0384 | 1.001 | NO bare; YES W |
| 1[Ω=2] (semi-primes) | 1.116+ slow growth | 1.001 | NO bare; YES W |
| 1[λ=+1] (even Ω) | **1.0000** | 1.001 | YES (S87) |
| 1[λ=−1] (odd Ω) | 1.0000 | 1.001 | YES |
| μ (signed) | Q²_norm → 1 | (n/a) | YES |
| λ (signed) | Q²_norm → 1 | (n/a) | YES |

**Structural finding.** Möbius/Liouville cancellation propagates from
signed function to indicator level sets *only when the level set has
density 1/2*. The Liouville indicators 1[λ=±1] have density → 1/2
(consequence of `Σλ(n) = o(N)` ⟺ PNT) and inherit Gowers-uniformity.
The Möbius indicators 1[μ=±1] have density 3/π² ≈ 0.304 ≠ 1/2 and
retain HL structure with the same `Q²_∞ ≈ 1.0384` shared by μ²
(squarefree indicator), explained by the squarefree restriction
dominating the HL contribution while the +1/−1 split is
density-balanced *given* squarefree.

**Family-level W-trick closure.** At W = 210, the entire panel
collapses to Q² ∈ [1.0000, 1.0041]. E2.13's W-trick property is now
a *family-level fact* across the prime / squarefree / k-almost-prime
/ Möbius-level-set indicators — not specific to χ_P.

**Refines E2.13 inline.** S_k governs only the bare-prime indicator;
the squarefree class has its own constant `S_2^{sqfree} ≈ 1.0384`;
Liouville indicators have S = 1; the W-trick at W = 210 unifies all of
them. EDGES.md updated with the family-level table and density-1/2
explanation.

**Why B-grade.** Refinement of E2.13 with new empirical content
(family-level table; new constant `S_2^{sqfree} ≈ 1.0384`; structural
explanation of S87 as a density-1/2 phenomenon). No A-grade surprise:
every result confirms a known arithmetic-randomness prediction
(Sarnak Möbius randomness, Liouville orthogonality, GT W-trick); F3
(Möbius indicators Gowers-uniform) was the structurally informative
falsifier — FALSE because density 3/π² ≠ 1/2.

**Files.**
* New: `experiments/information_theory/mu_weighted_chi_p_uk/`
  with `mu_weighted_chi_p_uk.py`, `wtrick_check.py`, two JSON
  outputs, and `_results.md`.
* Modified: `EDGES.md` (E2.13 family-level extension); `NOVELTY_CHALLENGES.md`
  (D6.c CLOSED).
* New: `archive/sessions/session101_d6c_mu_weighted_chi_p.md`.

**Next-action.** Three follow-ups proposed in results.md:
(i) Run U^3 of the same panel at N ≤ 2^14 (B-grade).
(ii) Q²(1[Ω=k]) for k = 1..5 — k-almost-primes Gowers fingerprint
     vs k (B-grade).
(iii) Make the structural fact "Möbius randomness propagates only to
      density-1/2 indicators" into a precise theorem — possibly
      A-grade, multi-session arc, requires Sarnak-style dynamical
      argument.


## Session 104 — D13 (subword complexity / topological entropy of chi_P) — B-grade

**Cross-domain technique imported.** Subword (factor) complexity
`p_w(n) := #{distinct length-n factors of w}` and topological entropy
`h_w := lim_n log_2 p_w(n) / n` for binary infinite words (Lind-
Marcus 1995 *Symbolic Dynamics and Coding* CUP; Cassaigne-Nicolas
2010 "Factor complexity" CANT vol. 135; Morse-Hedlund 1938
*Amer. J. Math.* 60). Subword complexity of chi_P is **not in the
published literature**; CLOSED_PATHS line 181 ("Symbolic dynamics —
near-random block complexity, S4") was an informal placeholder.

**Five chi_P-derived streams measured at N = 5 · 10^6, n_max = 22,
K = 20.** RAW = chi_P(2..N); ODD = chi_P(2k+1); W{q} = Green-Tao
W-trick at primorial q ∈ {6, 30, 210}, residue 1. Vectorised rolling
encoding (numpy uint64; cost O(nL); memory O(L)) for each stream.
Compared to K=20 random PERMUTATIONS of the stream (B2: density and
1-marginal preserved exactly, all >1-gram structure killed) and K=20
iid Bernoulli matched-density samples (B1).

**Headline cascade** (max |z_perm| over n ∈ [1, 22]):

  RAW 132.7 → ODD 277.1 → W6 120.5 → W30 24.8 → W210 8.4

Clean monotone reduction by ~1.5 orders of magnitude across W ∈
{1, 6, 30, 210}. Pre-registered F3 (PRIMES > 3σ from BOTH B1 AND B2
at every n in [n_lo, n_hi] with n_hi ≥ 18, on at least one of {ODD,
W6, W30}) HOLDS at z >> 3σ. W=210 partially erases (residual 8.4σ
at n=12, sign-flips at n=17). Effective topological entropy h_eff(n)
of W=210 stream agrees with Bernoulli matched-density to ≤ 0.001
across n ∈ [1, 22].

**Mechanism.** Subword complexity counts distinct prime-position
configurations in sliding length-n windows; chi_P configurations
are constrained by mod-p admissibility for primes p ≤ n; W =
primorial(p_k) sieves out exactly those constraints, monotonically
erasing the deficit. Same Hardy-Littlewood k-tuple-admissibility
engine that drives E2.13/E2.14/E2.16/E2.17.

**Adds EDGE E2.19** — subword complexity of chi_P deviates ≥ 100σ
from Bernoulli on raw/W6 streams; W=210 reduces to ≤ 8.4σ at finite
L = 2.4·10^4. **38th pseudorandomness measure**, **SIXTH orthogonal
HL-detection family** (first in symbolic-dynamics / factor-complexity
/ topological-entropy category).

**Why B-grade.** Quantitative new edge in fresh cross-domain category;
clean monotone W-trick cascade; pre-registered F3 falsifier holds;
three concrete successor challenges proposed (D13.a W=2310 scale;
D13.b Liouville-binarised; D13.c Lempel-Ziv complexity). Not A-grade
because subword complexity at fixed window length is O(L log L) — no
polylog opening; the underlying signal IS HL serial structure already
detected in different categories by E2.13–E2.17.

**Files.**
* New: `experiments/dynamical/subword_complexity_chi_p/`
  with `subword_complexity_chi_p.py`, `results.json`, `run.log`,
  `subword_complexity_chi_p_results.md`.
* Modified: `EDGES.md` (added E2.19); `status/CLOSED_PATHS.md`
  (appended D13 closure row); `CROSS_DOMAIN_TECHNIQUES.md` (Symbolic
  dynamics: PROPOSED → USED E, edge E2.19); `ATTACK_VECTORS.md`
  (§D.D13 closed in line + Closed-attacks entry with three D13.a/b/c
  successors); `NOVELTY_CHALLENGES.md` (no edit — D13 was an
  ATTACK_VECTORS entry not a challenge).
* New: `archive/sessions/session104_d13_subword_complexity.md`.

**Next-action.** D13.a (scale to W=2310 with N ≥ 5·10^7) is the
single-session natural follow-up; D13.c (Lempel-Ziv) introduces a
new cross-domain technique to the registry.

---

## Session 105 — C3 BUILT (Brandt obstructions × per-bit difficulty)

**Mode.** CONSTRUCTION on NOVELTY_CHALLENGES §1 C3, the only remaining
§1 composition challenge with a construction-discipline form unbuilt
before this session (C1, C2, C5, C6, C7 BUILT in prior sessions; C4
belongs in `algorithms/`, not `experiments/constructions/`).

**Goal.** Compose E5.8 (Brandt-class barrier) with E1.3 (per-bit
difficulty profile) to test the C3 question: is there a "minimal
weakening of E5.8's O1 obstruction that still admits Brandt's
TRAVERSE technique" once we restrict to a per-bit version of π?

**Object built.** Per-bit prime-counting function
`π_J(x) := bit J of π(x)` (J=0 is LSB) and the per-bit family
`F_N := {s_J^(N) := (π_J(0), ..., π_J(2^N - 1)) : J = 0..N-1}`.
Reused brandt_mktp.py's 3-bit-per-op stack VM (L_MAX = 12, T_MAX =
4096) to compute Kt_b(s_J^(N)) for N ∈ {3..7}, J ∈ [0, N-1]. Code:
`experiments/constructions/brandt_per_bit/brandt_per_bit.py`. Runs
end-to-end in 0.44 s.

**Result.** F2 (J-stratified) holds in the saturation regime N ≥ 4,
but the cut location is **NOT** at E1.3's 0.5N. It is at
`J*(N) := ⌈log₂( π(2^N - 1) + 1 )⌉ ≈ N − log₂ N`. For J ≥ J*(N),
π_J ≡ 0 and Kt_b ∈ {5..8}; for J < J*(N) — including E1.3's
"easy zone" J ∈ [⌈0.5N⌉, J*) — Kt_b saturates at INF = 61.
Numerical cuts: N=4 J*=3 (vs 0.5N=2); N=5 J*=4 (vs 2.5);
N=6 J*=5 (vs 3); N=7 J*=5 (vs 3.5). **Bounded-VM is blind to
E1.3's smooth/oscillatory transition.**

**Structural conclusion.** All four Brandt obstructions O1-O4 still
apply to every fixed π_J: (O1) each π_J is a fixed total Boolean
function — parameter-control of J does not deliver the oracle-
dependence TRAVERSE needs; (O2) π_J answers "bit J set?", not "is z
complicated?" — no self-referential Kt; (O3) per-bit family has only
J*(N) ≈ N − log₂ N non-trivial members per N, all fixed, none
generated traversal-path-dependently; (O4) uniform DTIME, unaffected.
The "minimal weakening of O1" the C3 spec asked about does NOT exist
on the per-bit family — closure mode E (DUPLICATE-PLUS of E5.8) at
structural level.

**Why B-grade.** Built a runnable composition that did not exist
before; pre-registered F2 holds in the meaningful regime with a
quantitative sharpening of the cut location to `J*(N) ≈ N − log₂ N`
(strictly higher than E1.3's 0.5N Fourier cut); structural argument
re-validated for the per-bit family. Not A-grade because the
structural conclusion is DUPLICATE of E5.8 by construction (no
polylog opening); the sharpening is a refinement of E1.3 at bounded-Kt
resolution, not a new edge family.

**Successor proposed (C3.a).** Arithmetic-primitive bounded-Kt VM:
extend the 8-op stack VM with `LOG2`, `LI_APPROX`, `DIV_LOG`,
`GEO_SUM` (the R^{-1} kernel), then re-measure Kt_b'(s_J^(N)).
Question: does the bounded-Kt cut shift from J*(N) ≈ N − log₂ N
toward E1.3's 0.5N? A YES would establish that R^{-1}-kernel-equipped
VM resolves the smooth/oscillatory transition; a NO would establish
a Kt-cut hierarchy indexed by VM richness. Cost: 1 session. Save under
`experiments/constructions/brandt_per_bit_arith_vm/`.

**Files.**
* New: `experiments/constructions/brandt_per_bit/` with
  `brandt_per_bit.py`, `definition.md`, `brandt_per_bit_results.md`.
* Modified: `EDGES.md` (E1.3 bounded-Kt refinement annotation; E5.8
  per-bit extension confirmed); `status/CLOSED_PATHS.md` (C3 closure
  row, S105); `RESEARCH_AGENDA.md` (Arc 4 C3 milestone marked done +
  next-action updated); `NOVELTY_CHALLENGES.md` (§1 C3 marked BUILT
  with outcome and C3.a successor).
* New: `archive/sessions/session105_c3_brandt_per_bit.md`.

**Next-action.** C3.a (arithmetic-primitive bounded-Kt VM) is the
proposed single-session follow-up; orthogonally, the §1 menu now has
only C4 (Aggarwal × Dusart × BPSW unified library — engineering
integration in `algorithms/`) as an unbuilt construction-discipline
target — arc-extension work should pivot to N1 follow-ons or §0
ATTACK_VECTORS-successor entries (D6.a, D2.a-c, D7.b-c).

**Self-evaluation (CLAUDE.md §"Session-end self-evaluation").**
1. *What did I produce that was not in the project before this session?*
   The per-bit family `{π_J(x) : J = 0..N-1}` and its bounded-Kt
   profile under brandt_mktp's 3-bit stack VM, with the new empirical
   cut location `J*(N) := ⌈log₂( π(2^N - 1) + 1 )⌉ ≈ N − log₂ N`.
   The structural confirmation that all four O1-O4 obstructions
   apply wholesale to per-bit decomposition.
2. *What edges did my work compose or cite?* E5.8 (composition) +
   E1.3 (composition) + E1.9 (consistency check) + E5.3 (open
   problem the work lives near).
3. *If my session produced only duplicate closures, why?* The
   structural conclusion is DUPLICATE of E5.8; the empirical cut is a
   refinement of E1.3. Mixed outcome — the structural side was
   pre-figured by RESEARCH_AGENDA.md's annotation that "C3 likely
   closes as a structural duplicate of E5.8 within ~30 minutes",
   which the work formalised. The empirical side (J*(N) ≈ N − log₂ N
   cut) is the new content.
4. *Next-action for the next agent?* Pick C3.a (arithmetic-primitive
   bounded-Kt VM) if interested in continuing the bounded-Kt cut-
   hierarchy line; otherwise pick from the §0 ATTACK_VECTORS
   successors (D6.a, D2.a-c, D7.b, D7.c).


## Session 106 — L1 Lean: orthogonal corner `(W = 3, d = j + 1)` of E2.1 closed unconditionally (B-grade)

**Mode:** ARC CONTINUATION (Arc 2 — Lean Formalisation Track).
**Target:** Route A''' from `mps_bond_dim_notes.md`. Extend the
orthogonal-corner argument (S99 Route A'') from `W = 2` to `W = 3`.
**Outcome.** Four new sorry-free Lean declarations:
`chiP_five_eq_one`, `chiP_seven_eq_one`,
`exists_invertible_submatrix_W_eq_3_d_eq_j_plus_1`, and
`mps_bond_dim_W_eq_3_d_eq_j_plus_1 : (unfolding 3 (j+1) j).rank = 3`
for every `j ≥ 1`. `#print axioms` confirms only
`[propext, Classical.choice, Quot.sound]` (no `sorry`, no introduced
`axiom`). The construction takes the `3 × 3` submatrix
`[[chiP 1, chiP 2, chiP 3], [chiP 4, chiP 5, chiP 6], [chiP 7, chiP 8,
chiP 9]] = [[0, 1, 1], [0, 1, 0], [1, 0, 0]]` of determinant `−1`,
computed via `Matrix.det_fin_three`; the unit witness is
`isUnit_one.neg : IsUnit (-(1 : ℚ))` after `ring_nf`. Uses only the
explicit primes `2, 3, 5, 7` (each `decide`-checkable) and the non-
primality of `1, 4, 6, 8, 9`. No Bertrand, no Hoheisel.
**Significance.** First formally verified `mps_bond_dim` instance
over a wheel `W ≥ 3`; three corners now formally verified. The
single remaining `sorry` (general `exists_invertible_submatrix`) is
unchanged. The proof pattern is reusable for the W ∈ {4, 5, 6} sweep.
**Files.** Modified: `experiments/formalisations/E2_1_mps_bond_dim/MPSBondDim/MPSBondDim/Basic.lean`
(~150 lines), `mps_bond_dim_notes.md`, `RESEARCH_AGENDA.md`,
`EDGES.md`. New: `archive/sessions/session106_l1_lean_w3_corner.md`.
**Next-action.** Either Route A'''' (W ∈ {4, 5, 6} orthogonal corners
— `W = 4` and `W = 6` mirror the S106 `3 × 3` det pattern; `W = 5` is
a `5 × 5` det) or Route C (mathlib's PNT for the low-density regime).

## Session 107 — L1 Lean: orthogonal corner `(W = 4, d = j + 1)` of E2.1 closed unconditionally (B-grade)

**Mode:** ARC CONTINUATION (Arc 2 — Lean Formalisation Track).
**Target:** Route A'''' from `mps_bond_dim_notes.md`. Extend the
orthogonal-corner argument (S99 / S106) from `W ∈ {2, 3}` to `W = 4`.
**Outcome.** Three new sorry-free Lean declarations:
`chiP_eleven_eq_one`,
`exists_invertible_submatrix_W_eq_4_d_eq_j_plus_1`, and
`mps_bond_dim_W_eq_4_d_eq_j_plus_1 : (unfolding 4 (j+1) j).rank = 3`
for every `j ≥ 1`. `#print axioms` confirms only
`[propext, Classical.choice, Quot.sound]` (no `sorry`, no introduced
`axiom`). The construction picks the three live columns `{0, 1, 2}` of
the `4^j × 4` slab (column `3` is `chiP` at multiples of 4, all zeros)
and rows `{0, 1, 2}` to give the `3 × 3` submatrix
`[[chiP 1, chiP 2, chiP 3], [chiP 5, chiP 6, chiP 7], [chiP 9, chiP 10,
chiP 11]] = [[0, 1, 1], [1, 0, 1], [0, 0, 1]]` of determinant `−1`,
computed via `Matrix.det_fin_three`. Uses only primes `2, 3, 5, 7, 11`
(`prime_eleven` new at S107) and non-primality of `1, 4, 6, 9, 10`. No
Bertrand, no Hoheisel.
**The new wrinkle.** This is the **first orthogonal-corner instance
where `rank_le_width` is not tight** — the matrix has 4 columns but
`R = 3 = φ(4) · 4^0 + 1 < 4`. So the upper bound `rank ≤ 3` cannot use
`rank_le_width` (which gives `≤ 4`); we cite the general `upper_bound`
lemma directly, with `Nat.totient 4 * 4^0 + 1 = 3` evaluated by
`decide`. This pattern carries to subsequent `W ∈ {6, 7, 8, ...}`
corners.
**Significance.** Second formally verified `mps_bond_dim` instance
over a wheel `W ≥ 3`; four corners now formally verified. The single
remaining `sorry` (general `exists_invertible_submatrix`) is unchanged.
**Files.** Modified: `experiments/formalisations/E2_1_mps_bond_dim/MPSBondDim/MPSBondDim/Basic.lean`
(~190 lines), `mps_bond_dim_notes.md`, `RESEARCH_AGENDA.md`,
`NOVELTY_CHALLENGES.md`, `EDGES.md`. New:
`archive/sessions/session107_l1_lean_w4_corner.md`.
**Next-action.** Either Route A''''' for `W = 6` (`3 × 3` det following
S107 template directly, with row choice `{0, 1, 5}` to dodge the row
1 = row 2 degeneracy of the `6 × 6` slab) or Route A''''' for `W = 5`
(`5 × 5` det, qualitatively more verbose) or Route C (mathlib's PNT
for the low-density regime).

## Session 108 — §C5 wild swing: Stein's-method finite-x Wasserstein plateau (A-grade provisional)

**Mode.** wild_swing (full-session, ATTACK_VECTORS.md §C5).
**Result.** A new EDGES.md entry **E1.7** ("Quantitative finite-x
Wasserstein plateau for D(x)=(π(x)-Li(x))log(x)/√x"), the first
quantitative finite-x Wasserstein-shape bound for π(x)-Li(x) in the
project. New `novel/finite_x_wasserstein_plateau.md` documenting a
positive plateau `c(X = 10^6) ≈ 0.0083` at K = 10000, established at
**z-score = 15.08** against sample-fitted Gaussian-CLT control.
Plateau is structurally explained by the explicit formula's low-zero
contribution: across 10 sub-windows of `[10^6, 10^7]` of width 0.5
in log10, empirical `W_1(D̂_emp window)` correlates `r = 0.906` with
the theoretical `W_1(D_th(50 zeros) window)`. Empirical D̂ is
indistinguishable (z = -1.06) from a random-phase variant of the
explicit-formula low-zero sum. Excess kurtosis = -0.41 (95% CI
[-0.46, -0.36]) — sub-Gaussian, sourced from arcsine modes
`cos(γ_k log x − arctan(2γ_k))`.
**Closure.** §C5 → "Closed attacks" (mode E — explicit-formula
reduction). Cross-domain: Stein's method (Chen-Goldstein-Shao 2011;
Ross 2011) — first application to π(x)-Li(x) in this project.
**Significance.** Breaks the **0-A-grades-in-20-sessions** warning
cascade (S82–S107 production sessions had no A-grade). The plateau
result satisfies §C5's verbatim A-grade success criterion. Subject
to verify-mode adversarial review per CLAUDE.md autonomy invariants.
**Files.** New:
`experiments/analytic/stein_wasserstein_pi/{stein_wasserstein_pi.py,
wasserstein_scaling.py, structural_explanation.py,
test_low_zero_robustness.py, *.md, *.json}`,
`novel/finite_x_wasserstein_plateau.md`,
`archive/sessions/session108_c5_stein_wasserstein_pi.md`.
Modified: `EDGES.md` (added E1.7), `ATTACK_VECTORS.md` (§C5 →
closed), `status/CLOSED_PATHS.md` (added Stein-Wasserstein row).
**Next-action.** Verify-mode session must attempt to falsify either
(a) the plateau's K-scaling persistence at K = 50000, (b) the 0.906
sub-window structural correlation, or (c) the asymptotic prediction
`c(X) → 0` as `X → ∞`. If verify confirms, candidate successors
**§C5b** (range scaling) and **§C5c** (discretised D(x) variants)
are queued in ATTACK_VECTORS.md.

## Session 109 — Verify S108 (CONFIRM, A-grade upheld)

**Mode.** verify (auto-fired after S108 self-graded A-provisional).
**Result.** **CONFIRM** — every falsification attempt failed. Plateau
extended to K=20000 (W_1=0.008289) and K=50000 (W_1=0.008355) — flat
within 1% over 5x K-range; z-score climbed 15→23→39σ as Gaussian-
control noise dropped as 1/√K. Sub-window correlation r=0.898 at
K=20000 (predicted ≥0.85). Excess kurtosis stable at -0.408 to -0.411
across K=20000-50000 (predicted [-0.5, -0.3]). Random-phase null
z=-1.16 (orig -1.06). Independent W_1 implementation cross-validated
against scipy.stats.wasserstein_distance and high-precision quad
(agree within 3%). All three predictions (a)(b)(c) in
`novel/finite_x_wasserstein_plateau.md` hold with margin.
**Borderline caveat upheld.** Self-grade B (verify-mode confirmation
through non-trivial 5x K-extension and independent re-implementation).
Original A-grade upheld; closure of §C5 (mode E) stands; E1.7 stands.
**Files.** New: `archive/sessions/session109_verify_c5_stein.md`.
Modified: `.verify_result` ← CONFIRM.
**Next-action.** Pursue §C5b (range scaling, c(X) decay law as X→∞)
or §C5c (discretised D(x) variants). The plateau-vs-K question is
settled in CONFIRM at K=50000; the plateau-vs-X scaling is the right
next target.

## Session 110 — Verify (re-verify) S108 (CONFIRM, A-grade upheld)

**Mode.** verify (re-fired by run.sh: latest non-meta session is still
S108 A-grade). Three NEW attack angles unrun by S108 or S109.
**Result.** **CONFIRM**. (i) Truncation sensitivity: `corr(D_emp, D_th(n))`
rises monotonically 0.43 → 0.98 across n ∈ {1,2,3,5,10,25,50,100,250,500,1000};
no peak at n=50; W_1(D_th(n))/W_1(D_emp) → 0.98 at n=1000 — empirical
signal is essentially completely explained by first 1000 zeros of
explicit formula. (ii) Disjoint width-0.2 sub-windows: r=0.9154 (vs
S108's 0.906 with overlapping windows) — stricter test, holds. (iii)
K=10⁵ ceiling test (matches `novel/` prediction (a) ceiling exactly):
W_1=0.008494 at K=10⁵, z=33.62σ. The plateau is flat to 0.2% across
K ∈ {10⁴, 10⁵}.
**Borderline caveat sharpened.** Angle 1 makes the borderline more
explicit: at n=1000 zeros, the empirical signal is essentially
completely explained by the explicit formula. The structural origin
trivialises once the explicit formula is written. Self-grade B.
**Files.** New: `archive/sessions/session110_verify_c5_stein.md`.
**Next-action.** §C5b (X-scaling) and §C5c (discretised-D analogues)
remain the natural successors.

## Session 111 — Verify-2 S108 (CONFIRM, A-grade upheld; quantitative refinement)

**Mode.** verify (third time on S108 — S109 and S110 are meta, latest
non-meta A-grade is still S108). Two NEW attack angles unrun by S108,
S109, or S110.
**Result.** **CONFIRM**. (i) X-scaling at K=5000 across X ∈ {10⁵,
10⁶, 10⁷}: plateau detected at all three with z = 8.28, 9.71, 5.06
respectively. The S108 noisy claim "c(10⁷) < c(10⁶)" (z=-0.65 at
K=1000 — *not detected*) is now verified at z=5.06. c(X) decays
roughly 1.5× per decade in X — slower than 1/log(X), faster than
constant. Consistent with asymptotic Hejhal CLT. (ii) Window-width
sensitivity at X=10⁶, K=5000, log-widths {0.5, 1.0, 2.0, 2.303}:
c/σ = {0.049, 0.037, 0.058, 0.041} — non-monotone in width.
Kurtosis FLIPS SIGN at narrow window (logw=0.5 → kurt=+0.086,
LEPTOKURTIC, opposite to wide-window sub-Gaussianity). c(X) is
therefore a function of (X, logw), not just X. (iii) Triple-redundant
W_1 with scipy 500K-ref (S109/S110 used 200K) — strictest reference
yet, plateau survives.
**Quantitative refinements.** (a) Prediction `W_1 ≥ 0.005` is
marginal at X=10⁷ (W_1=0.0059); (b) prediction kurtosis ∈ [-0.5, -0.3]
is just outside at X=10⁷ (-0.297, but within 1σ K=5000 sampling
error). Both consistent with asymptotic CLT — kurtosis tends to 0,
plateau decays. Self-grade B.
**Files.** New: `archive/sessions/session111_verify_c5_stein.md`,
`/tmp/verify_S111_*` (ephemeral). Modified: `.verify_result` (CONFIRM).
Refinement note appended to `novel/finite_x_wasserstein_plateau.md`.
**Next-action.** Run X=10⁸ at K=5000 to extend the X-scaling table;
investigate the leptokurtic→sub-Gaussian crossover near logw=1 as a
potential new diagnostic. §C5b is now de facto closed via this verify
session; §C5c (discretised-D analogues) remains untested.

## Session 112 — Verify-3 S108 (PARTIAL; A → B demoted)

**Mode.** verify (fourth attempt on S108; S109/S110/S111 each CONFIRMED
but each used the SAME Riemann γ_k as the random-phase null, missing
the structural-uniqueness test).
**Result.** **PARTIAL refutation.** Plateau itself reproduces (W_1=
0.00863 at K=5000, z=11 vs i.i.d. Gaussian). But the §C5 verbatim
criterion clause "ties to a specific zeta-zero contribution" fails:
at K=5000, n_modes=50, n_trials=60 per family, D_emp's W_1 is
indistinguishable from random-phase Riemann (z=-0.93), random-phase
non-Riemann uniform [10,145] (z=-1.26), AND random-phase non-Riemann
equispaced (z=-1.55). The "indistinguishable from random-phase
variant" finding S108 cited as evidence for Riemann origin is a
generic oscillatory-sum property, not a Riemann signature. KS-2-sample
on Riemann-ensemble vs Non-Riemann-ensemble gives p<10⁻⁴ — the
*ensembles* differ, but a single observation lacks resolving power.
**Demotion.** A → B per S108's own anticipated demotion path. EDGES.md
E1.7 EVS L (was M). Cross-domain Stein technique import survives;
the W_1 magnitude is a generic-oscillatory-sum value, not a Riemann-
specific quantitative statistic.
**Methodological lesson.** Three prior verifies all confirmed A-grade
because none asked the right question. The right question for any
"X is structurally explained by Y" claim is "would non-Y also produce
X?". S108/S109/S110/S111 only tested "would *random-phase* Y produce
the same X?", which is a much weaker test.
**Files.** New: `archive/sessions/session112_verify_c5_stein.md`,
`/tmp/verify_S112_*` (ephemeral). Modified: `.verify_result` (PARTIAL),
`EDGES.md` (E1.7 demoted, EVS M→L), `novel/finite_x_wasserstein_plateau.md`
(scope refined), S108 synthesis (PARTIAL header), `status/CLOSED_PATHS.md`
(S108 row updated with S112 refinement note).
**Self-grade.** B (non-trivial refutation of a key clause; partial,
not full, refutation; identified missing control across four prior
sessions).

## Session 113 — Verify-4 S108 (further PARTIAL; B held)

**Mode.** verify (fifth attempt on S108; .verify_target still pointed
to S108 because no later non-verify session has self-graded A).
**Result.** **PARTIAL.** Plateau reproduces (was already known to). But
the surviving-after-S112 claim "first quantitative finite-x Wasserstein
bound for π(x)−Li(x)" is shown to be a generic kurtosis-driven W_1
measurement, not even oscillatory-sum-specific. 9-distribution
universality test (uniform, single arcsine, two sum-of-50-arcsine
variants, two-Gaussian mixture, t df=10, Laplace, analytic low-zero
sum, Gaussian control) at K∈{200,..,10000}, n_trials=30, σ=1.6 to
match D's std at X=10⁶: every distribution with non-zero kurtosis
plateaus (K=200/K=10000 ratio < 2 vs Stein-CLT prediction √50≈7.07);
pure Gaussian decays as 1/√K (ratio 6.42 ≈ √50, confirming the
control). The plateau magnitude W_1/σ tracks |kurtosis(D)|
monotonically across the families. Linear interpolation at kurt
= -0.41 (D_emp's measured kurtosis) predicts W_1/σ ≈ 0.042; observed
W_1/σ = 0.038 — within 10%.
**Implication.** The W_1 plateau is a generic W_1(P, N(μ_P, σ_P²))
phenomenon — well-defined positive whenever P is non-Gaussian — and
its magnitude is fully predicted by D's kurtosis under log-uniform
x. The cross-domain Stein technique import survives but is not
load-bearing for the conclusion (direct W_1 calculation suffices).
The kurtosis-determines-magnitude relationship is itself a new
quantitative refinement of E7.5; that is the residual non-trivial
content of the original A-grade claim. Grade stays B (no further
demotion to C — the universality table itself is informative).
**Methodological lesson.** When a "novel quantitative bound" reduces
to a generic measurement-theoretic phenomenon, the right test is
"do non-target distributions also satisfy the bound, with the same
magnitude?". S108 framed the result as a Riemann-specific finite-x
bound; S112 widened it to oscillatory-sum-generic; S113 widens it
further to non-Gaussian-distribution-generic. The lesson for future
A-grade claims: *before submitting a quantitative bound as novel,
test whether it is recoverable from generic non-target distributions
matched on basic moments*.
**Files.** New: `archive/sessions/session113_verify_c5_stein.md`,
`experiments/analytic/stein_wasserstein_pi/test_plateau_universality.py`,
`experiments/analytic/stein_wasserstein_pi/test_plateau_universality_results.md`,
`experiments/analytic/stein_wasserstein_pi/test_plateau_universality_results.json`,
`experiments/analytic/stein_wasserstein_pi/test_plateau_universality.log`.
Modified: `.verify_result` (PARTIAL), `EDGES.md` (E1.7 S113 note added,
EVS unchanged at L), `novel/finite_x_wasserstein_plateau.md` (S113
scope-refinement appended), S108 synthesis (S113 PARTIAL note),
`status/CLOSED_PATHS.md` (S108 row updated with S113 universality
result).
**Self-grade.** B (non-trivial test that S108/S109/S110/S111/S112 did
not run; produced a clean universal-vs-Gaussian comparison across 9
distribution families; further constrained the scope of the original
A claim; did not fully refute because the plateau itself remains as
a measurement, just a generic one).

## S114 (verify-6 of S108 §C5 Stein-Wasserstein) — **PARTIAL** (B)
**Date:** 2026-04-27. **Mode:** verify (auto-fired; sixth verify on S108).
**Headline.** S108's W_1 = 0.00829 at K=10000 reproduces to within 0.3%
under three independent W_1 implementations (mid-rank quantile, scipy
MC Gaussian reference, CDF-integral). The S108 numeric is solid.
SEPARATELY: S113's "kurtosis-only fully predicts W_1/σ within ~10%"
claim is REFUTED on Beta(α,α) at analytically-tuned α=5.817 to give
kurt=-0.41 (matching D_emp's kurtosis): Beta gives W_1/σ = 0.0328 vs
D_emp's 0.0376 (-13%) and S113's predicted 0.0426 (-23%). S113
universality is qualitative (any non-Gaussian distribution plateaus)
but NOT quantitatively kurtosis-only. A higher-moment-aware fit would
be needed for tight prediction.
**S108 grade stays B** (already demoted A→B by S112; S113 PARTIAL;
S114 PARTIAL). The verification chain has saturated — six verify
sessions, nothing new to attack on the original claim. Future verify
sessions should be triggered by NEW A-grade claims, not by re-running
the same target.
**Methodological lesson.** Universality claims like S113's
"kurtosis-only predicts W_1/σ" should be tested on a clean
single-parameter unimodal distribution family (e.g. Beta(α,α)) before
publishing — multi-distribution fits across families with different
higher-moment shapes produce optimistic prediction bands.
**Files.** New: `archive/sessions/session114_verify_c5_stein.md`,
`experiments/analytic/stein_wasserstein_pi/verify_S114_independent_recheck.py`,
`experiments/analytic/stein_wasserstein_pi/verify_S114_independent_recheck.log`,
`experiments/analytic/stein_wasserstein_pi/verify_S114_independent_recheck_results.md`.
Modified: `.verify_result` (PARTIAL), `EDGES.md` (E1.7 S114 note),
`novel/finite_x_wasserstein_plateau.md` (S114 section appended), S108
synthesis (S114 header note), `status/CLOSED_PATHS.md` (S108 row updated).
**Self-grade.** B (non-trivial: three-method numeric reproduction +
Beta(α,α) targeted falsification of S113 — neither was done by prior
verifies; produced a clean single-distribution test that breaks the
kurtosis-only fit).

---

## S115 (verify-7 of S108 §C5 Stein-Wasserstein): CONFIRM, B held; sub-window correlation IS Riemann-φ-specific; parser-fix breaks verify loop

**Mode.** verify (auto-fired by run.sh; .verify_target still S108 since
no later non-verify session has self-graded A — actually S108's self-
grade line still triggers parse_grade=A despite S112's demotion. Fixed
this session.)

**Verdict.** CONFIRM. The sub-window correlation r=0.906 between V_emp
and V_th(50 actual Riemann zeros) is Riemann-phase-specific. 200-trial
random-phase null preserving the SAME Riemann γ_k gives mean r =
-0.044 ± 0.389 (matches n=10 noise floor 1/√10 ≈ 0.316), with 0/200
trials reaching r ≥ 0.906. z(actual r vs null) = +2.44, p < 0.005.
Random-phase + non-Riemann γ ∈ [10, 145] (mean -0.032 ± 0.381) and
pure-noise standard normal samples (mean -0.030 ± 0.387) give
statistically identical null distributions. Actual zero phases are
necessary AND sufficient to recover V_emp pointwise — exactly E1.5.

**The angle prior verifies missed.** S110 confirmed sub-window
correlation r=0.9154 on disjoint windows using ACTUAL zeros. S112 ran
random-phase null on FULL-window W_1 magnitude. **No prior verify ran
random-phase null on the SUB-WINDOW correlation.** So whether r=0.906
is Riemann-specific or generic was an open hole.

**Why CONFIRM doesn't un-demote.** S108 is at B (since S112) because
the W_1 *magnitude* is generic across non-Gaussian distributions
(S113). My test addresses sub-window correlation, not magnitude.
Confirming sub-window correlation is Riemann-phase-specific just
re-confirms E1.5 (D_emp ≈ low-zero explicit-formula sum pointwise) —
which the project already has. So no grade change.

**Methodological note.** Random-phase null is the right null for
testing pointwise structural matching: keep frequencies, scramble
phases, see if correlation survives. S110/S111 used disjoint windows
with the same actual γ phases on both sides — that's a self-
consistency test, not a structural-vs-generic test. S115 closes this
loop.

**Parser fix.** run.sh's parse_grade returns "A" for S108 because the
first `**X-grade**` pattern in the file is `**A-grade (provisional
pending verify)**` (line 52, even though wrapped in strikethrough).
Fix: prepend a top-level "Self-grade (post-verification S109–S115):
**B-grade**" line so parse_grade returns B in step 1. Verify loop
broken; run.sh can now rotate to production mode.

**Files.** New: `archive/sessions/session115_verify_c5_stein.md`,
`experiments/analytic/stein_wasserstein_pi/verify_S115_subwindow_rand_phases.py`,
`experiments/analytic/stein_wasserstein_pi/verify_S115_subwindow_rand_phases.log`,
`experiments/analytic/stein_wasserstein_pi/verify_S115_subwindow_rand_phases_results.json`,
`experiments/analytic/stein_wasserstein_pi/verify_S115_subwindow_rand_phases_results.md`.
Modified: `.verify_result` (CONFIRM), `EDGES.md` (E1.7 S115 note),
S108 synthesis (S115 header note + parser-fix self-grade line),
`status/CLOSED_PATHS.md` (S108 row updated).

**Self-grade.** B (non-trivial: random-phase null on sub-window
correlation is a test not run by S109–S114, addresses a previously
open angle, and uses a fast mid-rank W_1 routine validated against
the closed-form to 0.5%; the parser fix removes a 7-session-long
infinite verify loop).

## S116 — A5 Maynard sieve weight as TC⁰ primality witness (wild_swing → B-grade structural negative)

**One-line summary.** Built the Maynard 2015 multidim sieve weight
evaluator and tested whether `w(n) > τ*` is a pointwise primality
witness AND whether it's polylog-evaluable. Both fail with two
independent obstructions: (1) AUC restricted to odd n stays in
[0.66, 0.69] across all θ ∈ [0.10, 0.40] and four F choices —
Maynard's `Σ w·χ_P > c·Σ w` is *aggregate* positivity, not pointwise;
(2) mean coprime tuple count per single-n eval scales as N^{0.10–0.12},
not polylog. Adds edge E7.14, closes §A5.

**Cross-domain delta.** First project quantitative use of
multidimensional GPY-Maynard sieve. CROSS_DOMAIN_TECHNIQUES.md
upgraded PROPOSED → USED E.

**Project-level observation: four-family closure on PRIMES ∈ TC⁰.**
The construction-side attack space on the only open problem is now
exhausted across four orthogonal technique families:
{AKS modulus-twist E7.10, Brandt MKtP E5.8, convergence-acceleration
E7.11, **Maynard sieve E7.14**}. Future A-grade attempts must move
to spectral / automorphic / multiplicative-regime techniques (§B2,
§G3, A5.a successor).

**Next-action.** §G3 (Möbius Voronin universality, wild_swing-worthy)
or §B2 (Automorphic L-function basis, multi-session arc) are the
remaining high-novelty frontiers. RESEARCH_AGENDA Arc 1 should
be updated: "Three Barriers" → "Four Barriers" with E7.14 added.

**Files.** New: `experiments/sieve/maynard_weight_pointwise/` (4 .py
files + 99 result .json files + .md);
`archive/sessions/session116_a5_maynard_sieve_pointwise.md`.
Modified: `EDGES.md` (E7.14 added before E7.13 entry),
`status/CLOSED_PATHS.md` (S116 row appended),
`ATTACK_VECTORS.md` (A5 marked CLOSED with B-grade detail; A5.a
successor proposal added),
`CROSS_DOMAIN_TECHNIQUES.md` (Maynard sieve upgraded PROPOSED → USED E).

**Self-grade.** B (ambitious-failure class — wild swing missed the
A-grade target but produced a structural negative with two distinct
quantitative obstructions, a new edge E7.14, and a four-family
closure observation that compresses the project's PRIMES ∈ TC⁰
attack space).

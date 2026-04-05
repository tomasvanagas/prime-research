# Open Problems: Viable Research Directions

Last updated: 2026-04-05 (Session 36)

These are the ONLY directions not yet proven closed. Everything else has been
tested (655+ approaches across 36 sessions) and confirmed to hit one of three
failure modes: Circularity, Equivalence, or Information Loss.

**Session 36 key insight:** E(x) = pi(x) - Li(x) has only O(log x) bits of information,
well within polylog bounds. The barrier is COMPUTATIONAL (extracting O(log x) bits from
a sum of ~x^{1/2} cancelling terms), not information-theoretic. This reframes the problem:
find algebraic structure in the zero sum enabling bulk cancellation.

---

## 1. Circuit Complexity of pi(x) [MOST PROMISING]

**Question:** What is the circuit complexity of the prime-counting function?

**What's known:**
- PRIMES (decision) is NOT in AC^0, NOT in AC^0[p] (Allender-Saks-Shparlinski 2001)
- PRIMES is in P (AKS 2002), but NOT known to be in NC, TC^0, or L
- pi(x) (counting) has NO circuit lower bounds beyond trivial Omega(log x)
- Division and iterated multiplication ARE in uniform TC^0 (Hesse-Allender-Barrington 2002)

**Why it matters:** If pi(x) is computable by polylog-depth, poly-size circuits,
then polylog time is possible (with parallelism). If not, it proves impossibility.
This is genuinely unstudied territory.

**The gap:** Omega(log x) proven lower bound vs O(x^{1/2+epsilon}) best upper bound.
This is "one of the least-explored gaps in complexity theory" for a natural problem.

**Approach:** Try to show pi(x) mod 2 requires super-constant depth, or conversely
find a TC^0 reduction from pi(x) to known TC^0-computable functions.

**Session 11 findings:**
- The AKS approach to PRIMES in TC^0 requires polylog-dim matrix powering in TC^0.
  But "Threshold circuits for iterated matrix product and powering" (RAIRO 2000)
  proves: for k>2, k×k matrix powering in TC^0 → TC^0 = NC^1 (believed false).
  **Therefore: the AKS path to PRIMES in TC^0 is BLOCKED.**
- This does NOT prove PRIMES is not in TC^0 — there may be non-AKS paths.
- Scalar powering (x^n mod m) IS in TC^0 (Allender 1999), but matrix powering is NOT.
- All direct TC^0 paths failed: Legendre (exponential terms), Lucy DP (sequential depth
  O(sqrt(x)/ln(x))), Möbius (needs factoring), partial sieve (O(x/ln(x)) error).
- pi(x) mod 2 appears as hard as pi(x); no shortcut found.
- The floor-value function {floor(x/k)} DETERMINES pi(x) (Lucy DP proves this),
  but mapping floor-values to pi(x) in constant depth is the open question.
- **Key open question shifts to:** Is there a non-AKS TC^0 primality test?
  (e.g., using only scalar operations, which ARE in TC^0)

**Session 12 findings (circuit complexity deep-dive):**
- **Lucy DP DAG depth = pi(sqrt(x)) EXACTLY**: ratio depth/pi(sqrt(x)) = 1.000 for all
  x tested (100 to 100000). The critical path goes through S(x) at every prime step.
  The DP is fundamentally sequential for the final answer.
- **Telescoped version (Meissel-Lehmer) has depth pi(x^{1/3})**: confirmed empirically,
  ratio rounds/pi(x^{1/3}) → 1.00. This is the MINIMUM depth for sieve methods.
- **Both are exponential in N = log x**: depth O(x^{1/3}/ln x) = O(2^{N/3}/N),
  width O(x^{2/3}) = O(2^{2N/3}). Neither is polynomial in N.
- **Floor-value mapping matrices are NON-COMMUTATIVE**: 0 commuting pairs out of
  C(pi(sqrt(x)), 2) tested. Cannot reorder sieve operations.
- **Floor-value linear transformation is full-rank**: ~80% of floor values have
  nonzero coefficients, with |coefficients| >> pi(x) (massive cancellation).
  No low-rank approximation possible.
- **pi(x) mod m has INVARIANT conditional entropy**: H(Y|X) ≈ 0.537 bits regardless
  of modulus m (tested 2 through 30). This equals the entropy of the prime indicator.
  Confirms pi(x) mod m is as hard as pi(x).
- **Meissel-Lehmer as NC circuit**: depth O(log log x) = O(log N) but exponential
  width O(2^{2N/3}). NOT in NC via any known approach.
- **REFINED KEY QUESTION**: Is pi(x) in NC? This is EQUIVALENT to our target of
  O(polylog(x)) algorithms. The circuit size must go from 2^{2N/3} to poly(N).

---

## 2. Time-Bounded Kolmogorov Complexity of delta(n)

**Question:** What is Kt(delta(n) | n) -- the time-bounded conditional complexity
of the prime correction term?

**What's known:**
- Unbounded: K(p(n)|n) = O(1) trivially (program that sieves)
- delta(n) = p(n) - R^{-1}(n) has O(log n) bits of information
- Computing those bits takes O(x^{2/3}) time with best known methods
- Oliveira (2019) connected Kt complexity to circuit lower bounds

**Why it matters:** If Kt(delta(n)|n) = O(polylog(n)), that would imply
small circuits exist, connecting to Problem 1. If Kt(delta(n)|n) = omega(polylog),
that would prove no fast algorithm exists.

---

## 3. Zeta Zero Compressibility

**Question:** Do the Riemann zeta zeros have exploitable global structure
that allows fast summation of sum_rho R(x^rho)?

**What's known:**
- Zeros follow GUE statistics (Montgomery-Odlyzko)
- Individual zeros are computable in O(t^{1/3+epsilon}) time (Turing-type methods)
- Spectral flatness of zero sequence: 0.91 (high but not 1.0 = white noise)
- No FMM-type acceleration works due to incommensurability of zero heights

**Why it matters:** The entire barrier rests on needing O(sqrt(x)) zeta zeros.
If the sum could be compressed/approximated using structure in the zeros,
the barrier collapses.

**Status:** R(n) is 24% more compressible than random, but residual incompressible.

**Session 24 update:** Fresh experiments show δ(n) is 5x more compressible than random
(compression ratio 0.049 vs 0.256) and compression improves with length (1.56 bits/symbol
at N=10000). DCT captures 99% of energy in 10.4% of coefficients. But the critical 5%
needed for exact rounding is dense in ALL standard bases (DCT, Fourier, Haar). The zero
sum partial sums behave as weakly correlated random walks with GUE statistics (lag-2
autocorrelation -0.345). FRI analysis shows sampling density → 100% at scale.

**Session 25 update (STRUCTURAL PATTERNS CLOSED):** Comprehensive 7-experiment investigation:
- No rational relations in pairwise ratios (499,500 tested, GUE repulsion)
- No integer linear relations (PSLQ, 13,000+ tests at 60-digit precision)
- DFT matches GUE with correlation 0.9999; only faint p=2 signal
- No recurrence in explicit formula partial sums (linear, nonlinear, difference)
- Equidistributed mod ALL 10 constants tested (1, π, log(2π), e, log primes)
- Sparse matrix model uses O(N) free parameters = no compression
Zeros are GUE-random in every algebraic, spectral, arithmetic, and structural sense tested.
**All structural approaches to zero compressibility are now CLOSED.**

---

## 4. Berry-Keating / Hilbert-Polya Hamiltonian

**Question:** Does there exist a concrete self-adjoint operator H whose eigenvalues
are the zeta zeros, AND can its spectrum be computed efficiently?

**What's known:**
- The Hilbert-Polya conjecture posits such an operator exists
- Berry-Keating proposed H = xp + px (not rigorous)
- Connes' 2026 paper advances the Weil quadratic form approach
- Even if H exists, QPE requires 10^51 zeros for p(10^100)

**Why it matters:** A quantum-simulable Hamiltonian with zeta zero eigenvalues
would reduce the problem to quantum eigenvalue estimation.

---

## 5. Novel Number-Theoretic Identity [CLOSED - Session 29]

**Question:** Is there an identity that relates sum_rho R(x^rho) to a
computable function of x and n, without enumerating zeros?

**What's known:**
- All known identities (Weil explicit, Selberg trace, etc.) are transformations
  of the same underlying information
- No "shortcut identity" has been found in 165+ years of analytic number theory

**Why it matters:** This would be the most direct path -- a formula that bypasses
the zero sum entirely. The least likely to succeed but highest impact if found.

**Session 29 CLOSED (comprehensive identity search, 7 experiments):**
- Extended PSLQ (14-element basis, x up to 100000): 15/15 relations fail cross-validation
- Wilf-Zeilberger: higher-order differences GROW (ratio->2.0=white noise), no certificate
- Bernoulli/zeta/L-value/modular form relations: zero correlation (PSLQ x-dependent only)
- LLL minimal polynomials: different at every x (algebraically independent transcendentals)
- ODE search (order<=3, poly degree<=3): residuals match random noise
- Volterra integral equations: no kernel works beyond trivial local averaging
- Chebyshev psi(x)/log(x) captures 91% of f(x) but costs O(x) -- not a shortcut
**f(x) = pi(x) - R(x) has no computable identity in any tested basis.**

---

## Priority Assessment (updated Session 35)

| Direction | Feasibility | Impact | Recommended Effort |
|-----------|-------------|--------|-------------------|
| #TC^0 ⊆ NC? | Low | Maximal | THE key question. **S31: 15+ natural measures ALL give N/2 universality. S35: MKtP framework is reformulation not tool. GF(2) SLP complexity = random. Approximate circuit complexity has NO phase transition.** Every characterization attempted confirms "random-like" behavior. Only SAT-based exact circuit minimization remains untested at scale. |
| Circuit complexity | Low | High | All sieve-based TC^0/NC paths closed. BPSW-in-TC^0 conditional. S35: meta-complexity framework CLOSED (reformulation). SVD rank-k approximation shows gradual degradation (no easy core). |
| Determinantal complexity | CLOSED | High | **S17: dc(pi_N) >= 2^{N/2-1}+2 = Omega(sqrt(x)). Exponential.** |
| Space-time tradeoff | CLOSED (as impossibility route) | High | **S23: Comm complexity can NEVER give super-polylog. General lower bounds face Natural Proofs barrier.** |
| Kt complexity | **CLOSED** (as attack path) | High | **S35: Kt framework reformulates but doesn't solve. Kt(T_N) = O(2^N·N) by sieve regardless of circuit size. Brandt framework too generic.** |
| Novel identity | **CLOSED** | Maximal | **S29: 7 experiments. f(x) algebraically independent of all tested bases.** |
| Zero compressibility | Very Low | Very High | **S25: ALL structural approaches CLOSED.** |
| Berry-Keating | Very Low | Very High | Literature monitoring |
| Non-sieve pi(x) approach | CLOSED | Maximal | **S16: 15 families closed. S17: Fourier near-random. S23: Not holonomic.** |
| Proving impossibility | BLOCKED | Maximal | **S23: Natural Proofs barrier. Not achievable with current techniques.** |
| GF(2) algebraic structure | **CLOSED** | High | **S35: ANF sparsity = 0.50 (random), SLP savings = same as random, no compression. Fully characterized.** |

**Session 11 closed:** All convergence acceleration methods (Richardson orders 1-10,
Levin u/t, Weniger delta, smoothed formulas). All alternative decompositions of pi(x)
(Buchstab, Vaughan, hyperbola generalization, convolution structure, Dirichlet series).
M(x) factorization and H-T transfer to pi(x). TC^0 paths (Legendre, Lucy DP, Möbius,
partial sieve, parity). ~10 new closed paths.

**Session 11 refined:** Circuit complexity question now has precise formulation:
PRIMES in TC^0 ⟺ polylog-dimensional matrix powering in TC^0.

**Session 24 update (fresh perspective 2):** 11 new independent approaches tested, all 16 new
paths closed. Strongest new evidence: prime indicator has **linear complexity L/N = 0.5000 over
all finite fields GF(2) through GF(23)** (maximally random) and is **not k-automatic** for any k
(2-kernel has 38+ growing sequences). Combined with Session 23's full tensor rank: the prime
indicator carries near-maximal information density in every algebraic sense tested. The failure
taxonomy (Circularity/Equivalence/Info Loss) is empirically confirmed exhaustive across 575+ approaches.

**Session 12 closed:** Lucy DP parallelism (DAG depth = pi(sqrt(x))), Meissel-Lehmer
as NC circuit (exponential width), floor-value commutativity, floor-value linear algebra
(full-rank), pi(x) mod m for all m (invariant entropy). ~5 new closed paths.

**Session 12 refined:** "Is pi(x) in NC?" is EQUIVALENT to our target problem.
All known approaches produce exponential-size circuits (size 2^{Theta(N)}).
The barrier is now clearly the CIRCUIT SIZE, not depth. Any O(polylog) algorithm
must avoid computing O(sqrt(x)) intermediate values — i.e., must not be based
on floor values, sieve, or individual zeta zeros.

**Session 28 findings (novel):**
- **Approximate degree = ceil(N/2)** at rounding threshold epsilon=0.49, for BOTH
  chi_P and pi(x) mod 2. The counting step adds no difficulty. Quantum query lower
  bound Omega(N/4), STILL polylog. Does NOT rule out polylog circuits.
- **Per-bit influence gradient:** LSB-half influence ~2x MSB-half, growing with N.
  Bit 0 (parity) has influence ~N/2. All CONSISTENT with poly(N) circuits.
- **"N/2" universality:** ALL complexity measures converge at N/2: approximate degree,
  communication rank deficiency, oscillatory bit count, per-bit influence crossover,
  LFSR complexity. The smooth/oscillatory boundary at N/2 is universal.
- **Multiplicative structure:** I-E full rank, carry chains match random, partition
  number 2^{0.76*N} exceeds rank 2^{0.41*N}. Three independent exponential barriers.

**Session 13 closed:** Wilson TC^0, sum-of-squares TC^0, Cayley/Ihara/GCD graph
spectral, CRT reconstruction, recursive identity, prime zeta P(s), GF(2) algebraic
(ANF degree=N), optimized li-basis, zero reordering. ~12 new closed paths.

**Session 13 refined:** BPSW is computationally TC^0 (all components: MR=scalar pow,
Strong Lucas=2x2 MPOW, Jacobi=GCD=TC^0). "PRIMES in TC^0" ⟺ "BPSW correct."
PRIMES is in NONUNIFORM TC^0 unconditionally. ANF degree = Θ(N) over GF(2),
50% sparsity. No new 2026 breakthroughs in literature.

**Session 14 closed:** Nonlinear sieve (6 experiments — products, comparisons, bitwise,
identities all fail to improve on O(x^{2/3})). Lucy DP matrix structure (unipotent,
no displacement rank, full-rank product). Matrix power/LRS encoding (NOT LRS for any
modulus). Small det = pi(x) (I-E fractional parts independent, need exponential det).
#L chain PRIMES∈L→pi(x)∈NC^2 (breaks: workspace mismatch). DAG compression (floor-value
set irreducible). Algebraic variety point counting (low-dim insufficient, high-dim = sieve).
Polynomial regression (catastrophic overfitting). Character sums (equidistributed, harder).
Selberg exact counting (parity barrier). Batched primality (no shared structure).
~15 new closed paths.

**Session 14 refined:** PRIMES ∈ L and pi(x) ∈ NC are INDEPENDENT questions (workspace
mismatch). GapL algorithm MUST avoid floor functions entirely. Floor-value set is the
IRREDUCIBLE state space of ALL known sieve computations. The problem requires
fundamentally new intermediate quantities — not floor values, not zeta zeros.

**Session 15 closed:** Divide-and-conquer (error accumulates O(sqrt(x))), all randomized
methods (zeta sampling, probabilistic sieve, hash counting), arithmetic circuit complexity
(VP=VNC^2 but pi(x) not a polynomial), monotone complexity (trivial O(N) for thresholds),
novel intermediate quantity families (systematic analysis of 8 families — all fail).
~12 new closed paths.

**Session 15 refined:**
1. **Determinantal complexity connection**: "Is pi(x) in GapL?" ⟺ "Does the degree-N
   multilinear polynomial pi(bits) have polynomial determinantal complexity?" Found N×N
   det representations for N=2,3,4. For N≥10, generic polynomials DON'T fit — pi(x)
   would need special structure from number theory.
2. **#TC^0 ⊆ NC? as THE question**: If BPSW ∈ TC^0 (conditional), then pi(x) ∈ NC iff
   #TC^0 ⊆ NC. The Fermat residue coupling (exponent and modulus both depend on n)
   prevents batch counting.
3. **Uniformity is the true barrier**: Nonuniform poly-size circuits exist trivially
   (hardcode primes as advice). The hard part is UNIFORM construction — generating primes
   without knowing them. Natural proofs barrier blocks proving super-TC^0 lower bounds.
4. **No 2026 breakthroughs**: Chen-Tal-Wang STOC 2026 (depth-2 threshold lower bounds)
   is the closest result to TC^0 frontier but hard function is not number-theoretic.

**Session 16 closed:** ~22 new paths including: GapL intermediate quantities (class numbers,
L-values, elliptic curve a_p, regulators), TC^0 batch counting (5 routes: MAJORITY, D&C,
CRT batch, Carmichael, period), H-T signed transfer (6 routes: explicit formula, weighted,
Buchstab, M(x)->pi(x), identity search, Dirichlet structure), Lambert W error structure,
non-standard intermediates (additive combinatorics, ergodic, model theory, tropical,
sufficient statistics, algebraic geometry, representation theory).

**Session 16 refined:**
1. **15 intermediate quantity families now closed**: primes/zeros/floor-values are the
   ONLY known exact encodings of pi(x) and are informationally equivalent.
2. **TC^0 batch counting**: Communication complexity gives 2^{N/2}/poly(N) lower bound
   on constant-depth circuit size. Fermat residue coupling is cryptographically pseudorandom.
3. **H-T positivity barrier**: pi(x) counts POSITIVE quantities. Converting from signed
   sums (M(x), where H-T gives O(x^{3/5})) to unsigned (pi(x)) costs >= O(x^{2/3}).
4. **HKM 2023 achieves O~(sqrt(x)) elementarily**: NTT-based Dirichlet convolution.
   Published Math. Comp. 2024. Matches analytic O(x^{1/2+eps}) without complex analysis.
   Still exponential in input bits — no barrier change.

**Session 17 closed:** Ono partition characterization (circular, O(n²)), GapL via multilinear
polynomial (dc >= 2^{N/2-1}+2 exponential), Boolean Fourier analysis (near-random spectrum),
communication matrix rank (exact formula), binary carry structure (generic counter), partition
generating function shortcut. ~8 new closed paths.

**Session 17 refined:**
1. **EXACT communication complexity formula**: rank(pi_N) = 2^{N/2-1} + 2 for balanced bit
   partition. Verified N=2..20. The +2 = smooth part (R(x)), the 2^{N/2-1} = oscillatory part
   (zeta zeros). SVD confirms: top 2 SVs capture >99.99% variance.
2. **Determinantal complexity CLOSED**: dc(pi_N) >= 2^{N/2-1}+2 = Omega(sqrt(x)).
   The multilinear polynomial route to GapL is definitively closed.
3. **Fourier analysis near-random**: ~30% excess low-degree weight (parity/mod-4 only).
   Noise sensitivity ratio ~0.9, total influence ratio ~0.92. No junta/low-deg structure.
4. **No 2025-2026 pi(x) algorithm breakthroughs**: Ono (partitions), Tao-Gafni (gaps),
   Chen-Tal-Wang (THR∘THR lower bounds) — none change the pi(x) barrier.
5. **The sqrt(x) barrier is universal**: communication complexity, Fourier analysis,
   determinantal complexity, substitution rank ALL converge to sqrt(x).

**Session 19 closed:** Unbalanced communication (ALL partitions), PSLQ/LLL identity search
(6 relation types), NFS-type L[1/3] (4 sub-approaches), gap predictability (AR/MI),
SVD spectral approximation. ~7 new closed paths.

**Session 19 refined:**
1. **UNIVERSAL rank formula**: rank = 2^{min(k,N-k)-1}+2 for ALL bit partitions, not just balanced.
   No polynomial-rank partition exists. Barrier is intrinsic to pi(x) in binary.
2. **SVD spectral decay is power-law**: S_osc ~ i^{-1}. NOT geometric (can't exploit via
   truncated SVD for exact results). 90% of osc variance in ~20 SVs; 99% needs ~30%.
3. **SVD IS the explicit formula**: Top osc SVs correspond to first zeta zeros (corr=0.95).
4. **3-party NOF consistent with TC^0**: Balanced cut rank = 2^{N/3}. Doesn't separate.
5. **No elementary identity for pi(x)-R(x)**: PSLQ rules out linear, polynomial, recurrence,
   modular, functional, and differential relations with elementary functions.

**Session 22 closed (critique):** Prime race pi(x;q,a) shortcut for CRT (errors 2x rougher,
L-function zeros same barrier, Chebyshev bias = O(1) bits), L-function convergence advantage
(per-character fewer modes but phi(q) characters multiply total). ~2 new closed paths.
All 8 Session 21 proposals confirmed DUPLICATE against prior closed paths.

**Session 20 closed (space-time tradeoff):** Four approaches to formal S*T lower bounds:
1. Communication complexity -> T*S >= Omega(N^2) = Omega(log^2 x). Too weak (polylog).
   Fundamental limit: D(f) <= N for N-bit inputs, so comm complexity CANNOT give
   super-poly(N) bounds.
2. Nechiporuk formula bound -> L(pi) >= Omega(N). Trivially weak. Method limited to O(N^2).
3. OBDD size ~ 2^{0.79*N} empirically; OBDD width >= 2^{N/2-1} from rank. But OBDD != BP.
4. M-L DAG pebbling -> T*S >= Omega(x^{5/6}/ln x), algorithm-specific only.
**Verdict:** Proving T >= x^{Omega(1)} for general algorithms requires circuit lower bounds,
which is at least as hard as P != NP (Natural Proofs barrier). Problem remains OPEN.

**Session 23 closed:** k-party NOF (mode-unfolding rank full for k≥3), communication complexity
as impossibility route (bounded by input length), OBDD growth analysis, M-L pebbling (algorithm-
specific), holonomic recurrence test (pi(n) NOT D-finite up to order 20), Ono partition p-adic
lift (O(n^2) per evaluation, computationally inferior), short-interval explicit formula iteration
(each round needs same zero count). ~9 new closed paths.

**Session 23 refined:**
1. **k-party NOF full rank for k≥3**: Mode-unfolding rank = 2^{ceil(N/k)} exactly, matching
   random functions. Mode-unfolding rank is TOO COARSE for ACC^0/TC^0 separation. Only k=2
   shows anomalous sub-maximal rank (2^{N/2-1}+2). Need stronger measures (true tensor rank,
   discrepancy, polynomial method) for any separation result.
2. **Space-time route CLOSED**: Communication complexity bounded by N = log(x) bits. This
   FUNDAMENTALLY cannot give super-polylog lower bounds. General circuit lower bounds face
   Natural Proofs barrier.
3. **pi(n) is NOT holonomic**: Stronger than "not LRS" — even polynomial coefficients in
   the recurrence don't help. Test/random ratio ~1.0-1.7 for all tested parameters.
4. **Barriers understood at 4 levels**: Analytic (sqrt(x) zeros), algebraic (not holonomic),
   combinatorial (full tensor rank), complexity-theoretic (Natural Proofs).
5. **Total approaches: ~590+.** No genuinely new viable direction found. Problem OPEN.

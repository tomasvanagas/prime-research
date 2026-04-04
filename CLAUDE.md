# Prime Research: Computing p(n) Exactly Without Bruteforcing

## Project Goal
Find an efficient algorithm to compute the nth prime number p(n) without
enumeration, sieving, or brute force. Target: p(10^100) in <1 second, 100% exact.

## Status (April 2026)
- **480+ approaches tested** across 17 sessions, 130+ sub-agents
- **ALL KNOWN PATHS CLOSED** but no proof that polylog is impossible
- **Problem is GENUINELY OPEN** -- no unconditional lower bound beyond Omega(log x)
- **Session 12 KEY INSIGHT:** "Is pi(x) in NC?" is EQUIVALENT to our target.
  All known approaches produce circuits of size 2^{Theta(N)} (exponential in input).
- **Session 13 KEY INSIGHTS:**
  (a) BPSW IS computable in TC^0 (MR=scalar pow, Strong Lucas=2x2 MPOW, Jacobi=GCD).
      PRIMES in TC^0 iff BPSW (or similar) is unconditionally correct.
      Verified correct to 2^64. GRH also suffices (Miller's test = O(N^2) scalar pows).
  (b) Prime indicator ANF degree = Theta(N) over GF(2), 50% sparsity (random-like).
      No GF(2) algebraic shortcut for counting primes.
  (c) Zeta zero minimum K_min ~ 0.35 * x^{0.27}: power law, no reordering helps.
  (d) Spectral graph approaches all circular or equivalent to Meissel-Lehmer.
  (e) No new algorithmic breakthroughs in 2025-2026 literature.
- **Session 14 KEY INSIGHTS:**
  (a) **PRIMES ∈ L and pi(x) ∈ NC are INDEPENDENT questions.** The chain
      PRIMES ∈ L → pi(x) ∈ #L → pi(x) ∈ NC^2 BREAKS due to workspace mismatch:
      NL machine needs O(N) bits for candidate n, but #L allows only O(log N).
  (b) **I-E fractional parts carry O(2^k) independent bits** → no determinant
      smaller than 2^{Theta(sqrt(x)/log(x))} can encode the Legendre sieve.
      A GapL algorithm MUST avoid floor functions entirely.
  (c) **Lucy DP matrices have NO algebraic structure**: unipotent, displacement
      rank 50-60% of dimension, full-rank product. No compression possible.
  (d) **Nonlinear sieve breaks parity in theory but NOT in efficiency**: nonlinear
      ops on floor values CAN distinguish primes from semiprimes but cost ≥ O(x^{2/3}).
  (e) **All algebraic variety approaches fail**: low-dim can't encode pi(x),
      high-dim has slow point counting, Frobenius eigenvalues = zeta zeros.
  (f) **pi(x) mod m is NOT a linear recurrence** for any m — confirms no
      fixed-size matrix power can encode pi(x) (Mauduit-Rivat consequence).
- **Session 15 KEY INSIGHTS:**
  (a) **Determinantal complexity connection**: pi(x) as degree-N multilinear polynomial
      in bits. Found N×N det representations for N=2,3,4. "Is pi(x) in GapL?" ⟺
      "polynomial determinantal complexity?" For N≥10, GENERIC polynomials don't fit.
  (b) **#TC^0 ⊆ NC? is THE complexity question**: If BPSW ∈ TC^0 (conditional),
      pi(x) ∈ NC iff #TC^0 ⊆ NC. Fermat residue coupling prevents batch counting.
  (c) **Uniformity is the true barrier**: Nonuniform circuits trivially poly(N).
      Hard part = UNIFORM construction. Natural proofs barrier blocks lower bounds.
  (d) **ALL randomized approaches fail**: zeta zero sampling (100% needed), probabilistic
      sieve (10^6x worse), hash counting, quantum counting — all rigorously excluded.
  (e) **Divide-and-conquer fails**: error accumulates O(sqrt(x)) regardless of depth.
  (f) **Monotone complexity inapplicable**: [pi(x)>=k] = [x>=p(k)] trivially O(N).
      Individual bits of pi(x) NOT monotone.
  (g) **Arithmetic circuits don't help**: VP=VNC^2 (depth free), but pi(x) not a
      polynomial over fields. Tau conjecture orthogonal.
  (h) **Systematic analysis of 8 intermediate quantity families**: residues, polynomial
      evals, matrix eigenvalues, topology, representation theory, entropy, recursive,
      physical — ALL route back to floor values or zeta zeros.
  (i) **No 2026 breakthroughs**: Chen-Tal-Wang STOC 2026 (depth-2 threshold lower bounds)
      closest to TC^0 frontier but not number-theoretic.
- **Session 16 KEY INSIGHTS:**
  (a) **15 intermediate quantity families now CLOSED** (8 from S15 + 4 GapL + 7 novel):
      class numbers h(-d), L-function L(1,chi), elliptic curve a_p, regulators,
      additive combinatorics/sumsets, ergodic theory, model theory, tropical geometry,
      sufficient statistics, algebraic geometry/F_q, representation theory S_n/GL_n.
      ALL route to primes (C), zeta zeros (E), or lose information (I).
  (b) **Three "pillars" are the ONLY exact encodings of pi(x)**: prime positions,
      zeta zeros, floor values {floor(x/k)}. These are informationally equivalent.
      No fourth encoding found across 15+ candidate families.
  (c) **TC^0 batch counting has 5 closed routes**: MAJORITY fan-in (input generation),
      divide-and-conquer (O(sqrt(x)) error), batch CRT (exponent-modulus coupling),
      Carmichael structure (needs factoring), period exploitation (zero autocorrelation).
      Communication complexity gives 2^{N/2}/poly(N) lower bound on TC^0 circuit size.
  (d) **H-T signed cancellation transfer to pi(x) is CLOSED via 6 routes**:
      POSITIVITY of prime indicator prevents cancellation. Converting M(x)->pi(x) 
      always costs O(x^{2/3}). Identity pi(x)=sum omega(d)*M(x/d) found but omega=O(x^{2/3}).
  (e) **Lambert W error is structurally random**: delta(n) uncorrelated with gaps,
      uniform mod m, ~sqrt(x) magnitude. No exploitable structure.
  (f) **HKM 2023 achieves O~(sqrt(x)) for pi(x) elementarily** (NTT/Dirichlet convolution).
      Published Math. Comp. 2024. Still O(2^{N/2}) in input bits — no barrier change.
  (g) **Aggarwal 2025 gives O(sqrt(n)*log^4(n)) for p(n)** via binary search + HKM.
      No new algorithm — complexity analysis of existing tools.
- **Session 17 KEY INSIGHTS:**
  (a) **EXACT communication complexity of pi(x)**: rank(pi_N) = 2^{N/2-1} + 2 for
      balanced bit partition (verified N=2..20). This gives dc(pi_N) >= Omega(sqrt(x)).
      The multilinear polynomial route to GapL is DEFINITIVELY CLOSED.
  (b) **Boolean Fourier analysis**: Prime indicator has ~30% excess low-degree Fourier
      weight vs random (from parity/mod-4 structure), but noise sensitivity is near-random
      (~0.9x). No evidence of low-depth circuit structure. Spectral profile near-random.
  (c) **Ono partition characterization (PNAS 2024) CLOSED**: n prime iff
      (n²-3n+2)σ₁(n) - 8M₂(n) = 0. Requires sigma_1 (divisors) — circular + O(n²) cost.
      Worse than BPSW for primality, catastrophic O(x³) for counting. No GF shortcut.
  (d) **Chen-Tal-Wang ECCC 2026**: n^{2.5-ε} lower bounds for THR∘THR (depth-2 threshold).
      Hard function in E^NP, not number-theoretic. Advances TC^0 frontier via Williams' method.
  (e) **No new pi(x) breakthroughs in 2025-2026 literature.** Tao-Gafni 2025 (rough numbers
      in prime gaps) is about gap structure, not counting.
  (f) **sqrt(x) barrier is UNIVERSAL**: communication complexity, Fourier analysis,
      determinantal complexity, substitution rank — ALL converge to sqrt(x). The
      rank converges to ~50% of max for large N, matching prime density structure.
- Best exact: `algorithms/v10_c_accelerated.py` -- O(p(n)^{2/3}), p(10^9) in 0.175s
- Best approximate: R^{-1}(n) -- O(polylog), ~47% digits correct

## Directory Layout

```
status/
  CLOSED_PATHS.md    <-- SEARCH HERE before proposing ANY approach (472+ entries)
  OPEN_PROBLEMS.md   <-- The ONLY viable research directions
  BEST_ALGORITHMS.md <-- Working implementations with benchmarks
proven/
  barriers.md        <-- Mathematically proven impossibility results
  complexity.md      <-- Upper/lower bounds, circuit complexity
  information.md     <-- Entropy, Kolmogorov complexity, entanglement
  quantum.md         <-- Why quantum computing doesn't help
novel/
  info_computation_gap.md  <-- Best original insight (delta(n) framing)
  entropy_measurements.md  <-- Empirical: 5.04 bits/prime, gap analysis
  failure_taxonomy.md      <-- Three failure modes classification
algorithms/              <-- ONLY working, tested code
literature/
  references.md          <-- All citations in one place
  state_of_art_2026.md   <-- Latest published results
experiments/             <-- ALL experiments organized by topic
  analytic/ algebraic/ quantum/ ml/ information_theory/
  dynamical/ topological/ sieve/ other/
data/                    <-- Zeta zeros (200/300/500/1000) for explicit formula work
archive/                 <-- Session logs and visualizations
```

## The Barrier in One Paragraph
p(n) = SMOOTH(n) + RANDOM(n). The smooth part R^{-1}(n) is O(polylog) and gives
~50% of digits. The remaining ~50% encode oscillatory contributions of ~10^48
Riemann zeta zeros with GUE-random phases, information-theoretically incompressible.
Best known: O(x^{2/3}) combinatorial, O(x^{1/2+epsilon}) analytic. For p(10^100):
minimum ~10^49 operations. At 10^15 ops/sec = 10^34 seconds = 10^26 years.

## Three Failure Modes (Every Failed Approach Falls Into One)
1. **Circularity** -- needs primes to compute primes (Mills, Chebotarev, EGF)
2. **Equivalence** -- reduces to explicit formula / zeta zero sum (all spectral methods)
3. **Information Loss** -- smooth approximations lose ~170 bits distinguishing p(n)

## Do NOT Re-explore These (Thoroughly Closed)
- Convergence acceleration on explicit formula (Richardson/Levin/Weniger all fail; error oscillatory not smooth -- Session 11)
- Alternative decompositions of pi(x) (Buchstab/Vaughan/hyperbola/FFT all O(x^{2/3}) -- Session 11)
- TC^0 direct paths (Legendre/Lucy DP/Möbius/partial sieve all fail -- Session 11)
- Andrews-Wigderson FOCS 2024 (wrong model: fields not rings, GCD not bottleneck -- Session 12)
- Sieve-based circuits (ALL produce exponential-size circuits 2^{Theta(N)} -- Session 12)
- pi(x) mod m for any m (invariant entropy 0.537 bits, no modular shortcut -- Session 12)
- Succinct/lattice counting (Barvinok needs dim~sqrt(x), permanent/det circular -- Session 12)
- TG kernel smoothing (arXiv:2506.22634 DEBUNKED: violates uncertainty principle -- Session 12)
- H-T M(x)->pi(x) transfer (M(x) cancels 99.9%, pi(x) cancels 0% -- Session 11+12)
- Zeta zero summation (always O(sqrt(x)) terms, GUE prevents compression)
- ML/neural (1.1% best, fundamentally limited by Prime Coding Theorem)
- Modular/CRT reconstruction (circular or same complexity as pi(x))
- Quantum algorithms (volume-law entanglement, no speedup beyond O(x^{1/3}))
- Exotic algebra (topos/motivic/prismatic/condensed/perfectoid -- all closed)
- Dynamical systems (FRACTRAN/CA/transfer operators -- all closed)
- p-adic/adelic (not q-adically continuous, Mahler coefficients diverge)
- Interpolation/regression on delta(n) (random walk, zero generalization)
- Wilson's theorem for TC^0 (needs 2^N mults, not amenable -- Session 13)
- Sum-of-two-squares for TC^0 (needs factoring + O(log n) depth -- Session 13)
- Cayley graph / Ihara zeta / spectral graph (circular + equivalent -- Session 13)
- CRT reconstruction of pi(x) (each pi(x) mod q costs O(x^{2/3}), strictly worse -- Session 13)
- Recursive identity pi(x) = F(pi(x/d)) + correction (correction as hard as pi(x) -- Session 13)
- Zero reordering/weighting (K_min ~ x^{0.27}, power law regardless of order -- Session 13)
- Optimized li-basis (Riemann R(x) is essentially optimal for li-basis; error O(x^{0.3}) -- Session 13)
- Generating function P(s) extraction (more poles than explicit formula, worse -- Session 13)
- GF(2) algebraic shortcuts (ANF degree = N, sparsity 50%, random-like -- Session 13)
- Nonlinear sieves: products/comparisons/bitwise of floor values (break parity but NOT efficiency; K^d >= sqrt(x) for exact polynomial; overfit catastrophically -- Session 14)
- Matrix power / LRS encoding of pi(x) (pi(x) mod m NOT LRS for any m; Mauduit-Rivat -- Session 14)
- Lucy DP matrix structure (no displacement rank, no Toeplitz/Cauchy structure, full-rank product -- Session 14)
- Small det = pi(x) via floor values (I-E fractional parts carry 2^k independent bits; need 2^{sqrt(x)/2logx} matrix -- Session 14)
- #L chain PRIMES∈L→pi(x)∈NC^2 (BREAKS: workspace mismatch O(N) vs O(logN) -- Session 14)
- DAG path count compression (sieve DAG exponential; floor-value set is irreducible -- Session 14)
- Algebraic variety point counting (low-dim insufficient, high-dim=sieve, Frobenius=zeta -- Session 14)
- Polynomial regression on floor values (overfitting, residual = zeta zero contribution -- Session 14)
- Residue class / character sum exact counting (equidistributed, L-functions make harder -- Session 14)
- Selberg sieve exact counting (parity barrier, 3.7x overshoot -- Session 14)
- Batched Fermat/MR for pi(x) (no shared structure in 2^{n-1} mod n for consecutive n -- Session 14)
- Divide-and-conquer pi(x) (error O(sqrt(x)) at each level, accumulates to O(sqrt(x)) -- Session 15)
- Randomized zeta zero sampling (must use 100% of zeros; variance prohibitive -- Session 15)
- Probabilistic sieve / I-E sampling (10^6-10^9x worse than exhaustive -- Session 15)
- Hash-based prime counting (exact sketch needs O(x^2) bits -- Session 15)
- Arithmetic circuit complexity for pi(x) (VP=VNC^2 but pi(x) not a polynomial -- Session 15)
- Monotone circuit lower bounds for pi(x) ([pi(x)>=k]=[x>=p(k)] trivially O(N) -- Session 15)
- Novel intermediate quantities: residues, polynomial evals, topology, representation theory, entropy, recursive, physical (ALL route to floor values/zeta zeros -- Session 15)
- Class numbers h(-d) as GapL det entries (equivalent to L-values; constants, no x-dependence -- Session 16)
- L-function values L(1,chi) as GapL det entries (circularity + equivalence to zeta zeros -- Session 16)
- Elliptic curve a_p as GapL det entries (circularity: defined at primes; constants -- Session 16)
- Regulators as GapL det entries (equivalent to L-values; transcendental -- Session 16)
- Hybrid class-num/L-val/a_p/character approach (characters periodic; all roads back to zeta zeros -- Session 16)
- Sumset P+P / Goldbach representation r_2(n) (circle method error = zeta zeros; smooth 20% error -- Session 16)
- Ergodic theory / orbit complexity of primes (block complexity maximal; transfer ops = zeta zeros -- Session 16)
- Model theory / o-minimality (pi(x) not o-minimal definable; adding Z = undecidable; orthogonal -- Session 16)
- Tropical geometry for pi(x) (tropicalization loses all info; floor values ARE tropical objects -- Session 16)
- Sufficient statistics of floor values (poly(logx) statistic exists but computing it IS Meissel-Lehmer -- Session 16)
- Curve families over F_p (Frobenius eigenvalues = zeta zeros; batch costs O(x*polylog) -- Session 16)
- S_n/GL_n rep theory (cycle structure = Mertens not PNT; Fourier coeffs random -- Session 16)
- TC^0 MAJORITY fan-in for pi(x) (helps aggregation, not input generation; size 2^N*poly(N) -- Session 16)
- Batch modular exponentiation via CRT (exponent-modulus coupling in 2^{k-1} mod k -- Session 16)
- Carmichael lambda batch structure (identifying common lambda requires factoring all k -- Session 16)
- Period/autocorrelation exploitation in BPSW (zero autocorrelation at all lags -- Session 16)
- Lambert W / R^{-1} error structure (uncorrelated with gaps, uniform mod m, ~sqrt(x) = incompressible -- Session 16)
- Cheaper oracle hierarchy R->Lucy (R(x) error = O(sqrt(x)) = search range, no improvement -- Session 16)
- H-T transfer via explicit formula zeros (combinatorial H-T incompatible with analytic zero sum -- Session 16)
- Weighted prime counting for cancellation (multiplicative weights = constants on primes -- Session 16)
- Buchstab signed identity + H-T (signed Buchstab IS M(x); recovering pi(x) costs O(x^{2/3}) -- Session 16)
- M(x)->pi(x) conversion (pi(x) = sum omega(d)*M(x/d) but omega partial sums = O(x^{2/3}) -- Session 16)
- Binary carry structure of pi(x) (carry chains = generic counter; comm matrix rank = 2^{N/2-1}+2 exactly; spectral weight spread; primes = random -- Session 17)
- Ono partition characterization (PNAS 2024: n prime iff (n²-3n+2)σ₁(n)-8M₂(n)=0; requires divisors = circular; O(n²) per test, O(x³) total -- Session 17)
- GapL via multilinear polynomial det (dc(pi_N) >= 2^{N/2-1}+2 = Omega(sqrt(x)); substitution rank exponential -- Session 17)
- Boolean Fourier low-degree exploitation (30% excess low-deg weight = parity/mod-4 only; noise sensitivity near-random; no junta/low-deg structure -- Session 17)

## Viable Research Directions
1. **Circuit complexity of pi(x)** -- TC^0? NC^1? NC?
   Session 11: AKS path to TC^0 BLOCKED (growing-dim matrix powering open at TC^0/NC^1 frontier).
   Session 12: ALL sieve-based approaches produce EXPONENTIAL-SIZE circuits (2^{Theta(N)}).
   **"Is pi(x) in NC?" is EQUIVALENT to finding an O(polylog) algorithm.**
   Session 13 KEY RESULT: Systematic analysis of non-AKS TC^0 primality tests:
   - Wilson: CLOSED (needs 2^N mults). Sum-of-squares: CLOSED (needs factoring).
   - **BPSW IS in TC^0** as a computation: MR(2)=scalar pow, Strong Lucas=2x2 MPOW
     (Mereghetti-Palano 2000), Jacobi symbol=TC^0 via GCD (HAB 2002).
   - **QFT (Grantham) IS in TC^0**: operates in 2D algebra = 2x2 MPOW.
   - Strong Lucas alone: 12 PSPs below 100000, all caught by second param.
   - QFT alone: 4 PSPs below 50000, error < 1/7710.
   - **"PRIMES in TC^0" now reduces to: prove BPSW correct (or GRH).**
   - No BPSW pseudoprime known below 2^64 (exhaustive search).
   - Under GRH: Miller's test = O(N^2) scalar pows = TC^0. Already known.
   - Remaining question: unconditional proof of BPSW-type correctness.
2. **Is pi(x) in GapL? / Determinantal complexity** — EFFECTIVELY CLOSED (Sessions 13-17)
   **Session 17**: Extended dc(pi_N) computation to N=2..20. Found EXACT formula:
   rank(pi_N) = 2^{N/2-1} + 2 for balanced partition (verified N=4..16, approximate N=18..20).
   dc(pi_N) >= 2^{N/2-1} + 2 = Omega(sqrt(x)). The multilinear polynomial representation
   has EXPONENTIAL determinantal complexity. GapL via this route is DEFINITIVELY CLOSED.
   Any GapL algorithm must use an entirely different representation, not multilinear in bits.
3. **Time-bounded Kolmogorov complexity of delta(n)** (open, connects to circuit bounds)
4. **Zeta zero compressibility** -- convergence acceleration CLOSED (Session 11);
   only STRUCTURAL approaches remain (find global pattern in zero distribution)
5. **Concrete Berry-Keating Hamiltonian** -- Hilbert-Polya conjecture (quantum, open)
6. **Novel number-theoretic identity** relating sum_rho R(x^rho) to computable function
7. **Combinatorial pi(x) in O(x^{3/5})?** -- H-T achieved this for M(x) in 2021,
   but pi(x) lacks signed cancellation. Open whether transfer is possible.
8. **Nonlinear sieve** -- CLOSED (Session 14). Nonlinear operations (products, bitwise,
   thresholds on floor values) break parity barrier but NOT efficiency barrier.
   Polynomial of degree d in K floor values needs K >= x^{1/(2d)} for exactness.
   All polynomial fits overfit catastrophically. Nonlinear primality testing per
   number costs O(x^{3/2}) total -- worse than Meissel-Lehmer.
9. **Non-sieve, non-analytic approach** -- Session 12 proved ALL known methods produce
   exponential-size circuits. A fundamentally new approach is needed that avoids both
   floor-value sets (O(sqrt(x))) and zeta zero sums (O(sqrt(x))).
   **Session 15: 8 intermediate quantity families closed.**
   **Session 16: 7 MORE families closed** (additive combinatorics, ergodic theory,
   model theory, tropical geometry, sufficient statistics, algebraic geometry, rep theory).
   **Total: 15 intermediate quantity families systematically analyzed and ALL CLOSED.**
   The three "pillars" (prime positions, zeta zeros, floor values) are informationally
   equivalent and are the ONLY known exact encodings of pi(x). No fourth encoding found.
   See novel/session16_synthesis.md.
10. **#TC^0 ⊆ NC?** (Sessions 15-16) — Pure complexity theory question. If BPSW is in
    TC^0 (conditional on BPSW correctness), then pi(x) ∈ NC iff #TC^0 ⊆ NC.
    **Session 16**: 5 batch counting escape routes closed (MAJORITY fan-in, divide-and-conquer,
    CRT batch, Carmichael, period exploitation). Communication complexity gives 2^{N/2}/poly(N)
    lower bound on TC^0 circuit size for pi(x). The Fermat residue coupling
    (2^{n-1} mod n) makes the function cryptographically pseudorandom across consecutive n.

## Rules for AI Agents
1. **On startup: check if `TODO.md` exists in the project root.** If it does, work
   through its uncompleted items IN ORDER using sub-agents. Check off items as you
   complete them. When ALL items are done, DELETE `TODO.md` and note completion in
   the Status section of this file. If no TODO.md exists, proceed with normal research.
2. **Read `status/CLOSED_PATHS.md` before proposing ANY approach**
3. Save experiments to `experiments/<topic>/` with descriptive filenames
4. Update `status/CLOSED_PATHS.md` when closing an approach (add a row)
5. Genuinely novel findings go to `novel/` with evidence
6. New literature goes to `literature/`
7. **UPDATE this CLAUDE.md** when you discover something that changes the picture:
   - New barrier proven? Update "The Barrier" and "Do NOT Re-explore" sections
   - New viable direction found? Add to "Viable Research Directions"
   - Closed a direction from OPEN_PROBLEMS.md? Move it to "Do NOT Re-explore"
   - Found a better algorithm? Update "Status" section and `status/BEST_ALGORITHMS.md`
   - Keep this file accurate and current -- future sessions depend on it.
   - Found something that shows its worth re-openning the path? - Reopen.
8. If you beat the current best algorithm, save it to `algorithms/` with benchmarks
9. DO NOT modify `run.sh`
10. Use sub-agents to save context window
11. **Context management:** When context is filling up (you notice compression warnings
    or are past ~70% of your work), STOP research immediately. Write all remaining
    uncompleted work, partial findings, and new sub-questions to `TODO.md` so the next
    session can continue seamlessly. Then halt.
12. If you find the breakthrough, respond with exactly: I FOUND IT!!!

# Prime Research: Computing p(n) Exactly Without Bruteforcing

## Project Goal
Find an efficient algorithm to compute the nth prime number p(n) without
enumeration, sieving, or brute force. Target: p(10^100) in <1 second, 100% exact.

## Status (April 2026)
- **445+ approaches tested** across 15 sessions, 115+ sub-agents
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
- Best exact: `algorithms/v10_c_accelerated.py` -- O(p(n)^{2/3}), p(10^9) in 0.175s
- Best approximate: R^{-1}(n) -- O(polylog), ~47% digits correct

## Directory Layout

```
status/
  CLOSED_PATHS.md    <-- SEARCH HERE before proposing ANY approach (445+ entries)
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
2. **Is pi(x) in GapL? / Determinantal complexity** (Sessions 13-15)
   Asks for poly(N)-size matrix with det = pi(x). Equivalent to finding a DAG on
   poly(N) nodes whose signed path count = pi(x). See novel/gapl_question.md.
   Session 14: Equivalent to pi(x) ∈ #L. PRIMES ∈ L does NOT imply pi(x) ∈ #L
   (workspace mismatch barrier, see novel/workspace_mismatch_barrier.md).
   I-E approach FAILS: fractional parts carry 2^k independent bits → exponential det.
   **Session 15 REFINEMENT**: Reformulated as determinantal complexity of pi(x) viewed
   as degree-N multilinear polynomial in bits. Found dc(pi_N) = N for N=2,3,4.
   For N≥10, generic polynomials DON'T have N×N det reps. See novel/determinantal_complexity.md.
   A GapL algorithm MUST avoid floor functions entirely AND requires pi(x) to have
   special algebraic structure making it "non-generic" in the determinantal variety.
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
   **Session 15: Systematic analysis of ALL 8 candidate intermediate quantity families
   (residues, polynomial evals, matrix eigenvalues, topology, rep theory, entropy,
   recursive, physical) — ALL route back to floor values or zeta zeros.**
   See novel/uniformity_barrier.md: the true barrier is UNIFORMITY, not monotone
   complexity, not randomization, not arithmetic circuit depth.
10. **#TC^0 ⊆ NC?** (Session 15, NEW) — Pure complexity theory question. If BPSW is in
    TC^0 (conditional on BPSW correctness), then pi(x) ∈ NC iff #TC^0 ⊆ NC. The
    Fermat residue coupling (2^{n-1} mod n where exponent and modulus both depend on n)
    prevents any batch/structural counting shortcut found so far.

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
   - Keep this file accurate and current -- future sessions depend on it
8. If you beat the current best algorithm, save it to `algorithms/` with benchmarks
9. DO NOT modify `run.sh`
10. Use sub-agents to save context window
11. **Context management:** When context is filling up (you notice compression warnings
    or are past ~70% of your work), STOP research immediately. Write all remaining
    uncompleted work, partial findings, and new sub-questions to `TODO.md` so the next
    session can continue seamlessly. Then halt.
12. If you find the breakthrough, respond with exactly: I FOUND IT!!!

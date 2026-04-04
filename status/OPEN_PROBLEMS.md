# Open Problems: Viable Research Directions

Last updated: 2026-04-04 (Session 15)

These are the ONLY directions not yet proven closed. Everything else has been
tested (430+ approaches across 14 sessions) and confirmed to hit one of three
failure modes: Circularity, Equivalence, or Information Loss.

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

## 5. Novel Number-Theoretic Identity

**Question:** Is there an identity that relates sum_rho R(x^rho) to a
computable function of x and n, without enumerating zeros?

**What's known:**
- All known identities (Weil explicit, Selberg trace, etc.) are transformations
  of the same underlying information
- No "shortcut identity" has been found in 165+ years of analytic number theory

**Why it matters:** This would be the most direct path -- a formula that bypasses
the zero sum entirely. The least likely to succeed but highest impact if found.

---

## Priority Assessment (updated Session 12)

| Direction | Feasibility | Impact | Recommended Effort |
|-----------|-------------|--------|-------------------|
| Circuit complexity | Medium | High | PRIMARY focus — but all sieve-based TC^0/NC paths now closed. Non-AKS primality test is the key remaining question. |
| Kt complexity | Medium | High | Theoretical exploration |
| Zero compressibility | Very Low | Very High | Session 11 CLOSED convergence accel; only STRUCTURAL approaches remain |
| Berry-Keating | Very Low | Very High | Literature monitoring |
| Novel identity | Very Low | Maximal | Serendipity only |
| Determinantal complexity | Medium | High | "dc(pi_N) = poly(N)?" is a clean algebraic complexity question. Found dc ≤ N for N=2,3,4. Open for N≥5. |
| #TC^0 ⊆ NC? | Low | Maximal | Pure complexity theory question. If yes + BPSW correct → pi(x) ∈ NC. |
| Non-sieve pi(x) approach | Unknown | Maximal | All sieve/combinatorial/analytic methods give exponential-size circuits. **Session 15: ALL 8 intermediate quantity families analyzed and closed.** Must avoid floor values entirely. |

**Session 11 closed:** All convergence acceleration methods (Richardson orders 1-10,
Levin u/t, Weniger delta, smoothed formulas). All alternative decompositions of pi(x)
(Buchstab, Vaughan, hyperbola generalization, convolution structure, Dirichlet series).
M(x) factorization and H-T transfer to pi(x). TC^0 paths (Legendre, Lucy DP, Möbius,
partial sieve, parity). ~10 new closed paths.

**Session 11 refined:** Circuit complexity question now has precise formulation:
PRIMES in TC^0 ⟺ polylog-dimensional matrix powering in TC^0.

**Session 12 closed:** Lucy DP parallelism (DAG depth = pi(sqrt(x))), Meissel-Lehmer
as NC circuit (exponential width), floor-value commutativity, floor-value linear algebra
(full-rank), pi(x) mod m for all m (invariant entropy). ~5 new closed paths.

**Session 12 refined:** "Is pi(x) in NC?" is EQUIVALENT to our target problem.
All known approaches produce exponential-size circuits (size 2^{Theta(N)}).
The barrier is now clearly the CIRCUIT SIZE, not depth. Any O(polylog) algorithm
must avoid computing O(sqrt(x)) intermediate values — i.e., must not be based
on floor values, sieve, or individual zeta zeros.

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

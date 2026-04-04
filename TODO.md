# TODO: Session 15 — Reanalyze Existing Gems & Benchmark Against State-of-Art

**Created:** 2026-04-04
**Context:** We have 430+ closed approaches but several partial concepts and edge ideas
scattered across the codebase that deserve deeper, fresh analysis. Additionally,
kimwalisch/primecount (GitHub) represents the engineering state-of-art for combinatorial
pi(x) — Gourdon's variant of Deleglise-Rivat, O(x^{2/3}/log x), with SIMD/POPCNT
hard-special-leaves optimization and OpenMP scaling. Our v10 is ~1000x slower in practice.

**Instructions:** Work through these items using sub-agents. Check off each item when done.
When ALL items are complete, DELETE this file and note completion in CLAUDE.md status.
If context is filling up, STOP WORK, write remaining uncompleted items and any new
discoveries back to this TODO.md, then halt.

---

## Phase 1: Benchmark Gap Analysis

- [ ] **1.1** Study kimwalisch/primecount's Gourdon implementation details (especially
  Hard-Special-Leaves.pdf and the SIMD filtering approach). Document what specific
  engineering optimizations give it 100-1000x over naive Lucy DP. Save findings to
  `literature/primecount_analysis.md`.

- [ ] **1.2** Compare our `algorithms/v10_c_accelerated.py` against primecount's approach.
  Identify which optimizations from primecount could be ported to improve our v10.
  Is our Lucy DP the bottleneck vs Gourdon's special leaves? Document in
  `status/BEST_ALGORITHMS.md`.

## Phase 2: Reanalyze Existing Gems (Fresh Eyes)

These are partial concepts found in our codebase that were not fully explored or
were dismissed too quickly. Each deserves a focused sub-agent investigation.

- [ ] **2.1 GapL / Determinantal Complexity** — `novel/gapl_question.md` +
  `experiments/circuit_complexity/det_perm_encoding.py` +
  `experiments/circuit_complexity/det_complexity_search.py`
  Session 14 proved floor-value matrices fail. But the question "what NEW intermediate
  quantities could work?" was never explored. Investigate:
  - Class numbers of imaginary quadratic fields
  - Regulators of number fields
  - L-function special values (L(1,chi) for Dirichlet characters)
  - Elliptic curve point counts (Schoof-type, but as matrix entries not as endpoints)
  - Can any of these encode pi(x) as det(M) for poly-size M?
  Save results to `experiments/circuit_complexity/gapl_new_intermediates.py`.

- [ ] **2.2 BPSW TC^0 → Batch Counting** — `novel/bpsw_tc0_reduction.md`
  If PRIMES is in TC^0, primality testing is parallelizable in constant depth.
  But counting (pi(x)) still needs summation. Investigate:
  - TC^0 has MAJORITY gates — can we count primes via threshold circuits?
  - iterated addition of n bits is in TC^0 (Chandra-Stockmeyer-Vishkin).
  - So if BPSW is correct: for each k in [1,x], compute BPSW(k) in TC^0,
    then sum the indicator bits (also TC^0). Is pi(x) then in TC^0?
  - What's the circuit SIZE? poly(x) gates for x candidates — but x = 2^N,
    so size = 2^N * poly(N). Still exponential in input bits N = log x.
  - Is there a way to avoid enumerating all x candidates? This is the real question.
  Save analysis to `experiments/circuit_complexity/tc0_batch_counting.py`.

- [ ] **2.3 Lambert W Approximation Error Structure** — `algorithms/v1_pade_approximation.py`
  The Lambert W Pade gives ~0.003% error. The error IS the zeta-zero contribution.
  But investigate: does the error have exploitable structure at specific n values?
  - Compute error for p(n) at n = 10^k for k=1..9
  - Is error correlated with prime gaps? With pi(x) mod small primes?
  - Could a hybrid: Lambert_W + small correction table beat pure R^{-1}?
  - What if we use Lambert W to narrow the search range for binary search,
    reducing the number of pi(x) evaluations needed?
  Save to `experiments/analytic/lambert_w_error_analysis.py`.

- [ ] **2.4 Helfgott-Thompson Signed Cancellation Transfer** —
  `experiments/sieve/ht_transfer_attempt.py`
  H-T gets O(x^{3/5}) for M(x) via signed cancellation. Direct transfer fails
  because pi(x) has no sign cancellation. BUT:
  - What about pi(x) = li(x) - (1/2)li(x^{1/2}) - sum_rho li(x^rho) - ...?
    The explicit formula HAS signed terms. Can H-T's combinatorial trick apply
    to the explicit formula's oscillatory part rather than to the sieve?
  - Partial cancellation: even if pi(x) doesn't cancel like M(x), maybe a
    weighted version W(x) = sum_{n<=x} w(n) * 1_prime(n) does cancel, and
    pi(x) can be recovered from W(x) efficiently?
  Save to `experiments/sieve/ht_signed_transfer_v2.py`.

- [ ] **2.5 Aggarwal 2025 Method** — Referenced in `status/BEST_ALGORITHMS.md` as
  O(sqrt(n) * (log n)^4) for p(n). This uses analytic methods and is NOT implemented.
  - Find and read the actual paper/preprint
  - What is the actual algorithm? Is it Lagarias-Odlyzko with improvements?
  - Could it be implemented? What are the constants?
  - For p(10^100): sqrt(10^100) * log^4(10^100) ~ 10^50 * 10^9 ~ 10^59. Still huge
    but worth understanding the method.
  Save to `literature/aggarwal_2025_analysis.md`.

- [ ] **2.6 Non-Floor-Value, Non-Zeta Approaches** — The big open gap.
  Session 14's strongest conclusion: any polylog algorithm must use intermediate
  quantities that are NEITHER floor values NOR zeta zeros. Brainstorm and test:
  - Additive combinatorics (sumsets, Freiman's theorem applied to primes)
  - Ergodic theory (Green-Tao machinery, but for counting not patterns)
  - Model theory / o-minimality (definability of pi(x) in restricted structures)
  - Tropical geometry (tropicalization of the prime zeta function)
  - Information-theoretic: can we define a "sufficient statistic" for pi(x)
    that compresses the O(sqrt(x)) floor values into poly(log x) values?
  Save to `experiments/other/non_standard_intermediates.py`.

## Phase 3: Synthesis

- [ ] **3.1** After completing Phase 2, write a synthesis document:
  - Which gems turned out to have substance? Which were truly dead?
  - Any NEW viable research directions discovered?
  - Update `status/OPEN_PROBLEMS.md` and `CLAUDE.md` accordingly.
  Save to `novel/session15_synthesis.md`.

---

## Notes for Context Management
If context is running low before completing all items:
1. Update this file — check off completed items, add notes on partial progress
2. Add any new discoveries or sub-questions that emerged
3. The next agent session will pick up from where you left off

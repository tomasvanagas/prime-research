# Session 11 Summary (2026-04-04)

## Overview
Session 11 focused on three main directions:
1. **Convergence acceleration** — can we reduce zeros needed from O(sqrt(x))?
2. **Circuit complexity (TC^0)** — can pi(x) be computed in constant depth?
3. **Alternative decompositions** — can we beat O(x^{2/3}) combinatorially?

All three directions produced significant NEGATIVE results with one notable
POSITIVE finding (the TC^0 ⟺ matrix powering reduction).

## Key Findings

### 1. Convergence Acceleration: CLOSED
Tested 6+ methods (Richardson orders 1-10, Levin u/t, Weniger delta, smoothed
explicit formulas, optimal linear combination). All FAIL for large x because:
- Error structure is oscillatory (driven by GUE-random zeta zero phases)
- Sign change rate ≈ 0 (error consistently biased, not alternating)
- Minimum zeros scale as x^{0.25-0.37} (power law, NOT polylogarithmic)
- Even the "oracle optimal" linear combination can't beat the information barrier

**Files:** experiments/analytic/advanced_convergence.py, zero_scaling.py

### 2. TC^0 Analysis: KEY REDUCTION FOUND
Investigated 7+ paths to TC^0 computation of pi(x). All direct paths fail, but:

**Novel finding:** AKS path to PRIMES in TC^0 is BLOCKED:
- AKS reduces to polylog-dimensional matrix powering M^n mod n
- BUT: RAIRO 2000 proves k×k matrix powering in TC^0 → TC^0=NC^1 (believed false!)
- Therefore AKS path is blocked; non-AKS paths remain open
- Commutativity of the polynomial ring does NOT help (Frobenius circularity)

**Files:** experiments/circuit_complexity/tc0_reduction_analysis.py, 
commutative_matrix_powering.py; novel/tc0_matrix_powering_reduction.md

### 3. Alternative Decompositions: CLOSED
Tested 6 decomposition strategies (Buchstab tree, hyperbola generalization,
Vaughan identity, convolution structure, Dirichlet series, novel DPs).
All reduce to O(x^{2/3}) because:
- The bottleneck is the SELF-REFERENTIAL structure of pi(x) computations
- No factorization mu = f*g with both partial sums easy exists
- Helfgott-Thompson O(x^{3/5}) for M(x) does NOT transfer to pi(x)
  (pi(x) is positive, no signed cancellation)

**Files:** experiments/algebraic/alternative_decompositions.py,
experiments/sieve/mertens_speedup.py, ht_transfer_attempt.py

## New Closed Paths (~10 new entries)
- Higher-order Richardson (orders 2-10)
- Levin u/t-transform, Weniger delta
- Smoothed explicit formula
- All TC^0 direct paths (Legendre, Lucy DP, Mobius, partial sieve, parity)
- H-T transfer to pi(x)
- Alternative decompositions (6 strategies)

## New Novel Findings (2)
- novel/tc0_matrix_powering_reduction.md — precise TC^0 reduction
- novel/convergence_acceleration_barrier.md — systematic barrier proof

## Updated Files
- CLAUDE.md — status, closed paths, viable directions
- status/CLOSED_PATHS.md — ~10 new entries
- status/OPEN_PROBLEMS.md — Session 11 findings, priority update
- proven/barriers.md — two new barrier entries + "What Is NOT Proven" update

## Refined Understanding (from matrix powering literature survey)
- Fixed-k MPOW IS in TC^0 (Mereghetti-Palano 2000); IMP_k (k≥3) NOT in TC^0
- Growing-dim MPOW (k=polylog(n), our case) is GENUINELY OPEN at TC^0/NC^1 frontier
- Bottleneck: combining polylog eigenvalue powers needs depth O(log log n)
- Healy-Viola 2006: F_{2^n} exponentiation IS in TC^0 via Frobenius, doesn't work over Z_n
- Full survey: literature/matrix_powering_tc0.md

## Remaining Open Directions
1. **Growing-dim MPOW in TC^0** — the precise formulation at TC^0/NC^1 frontier
2. **Kt complexity of δ(n)** — time-bounded Kolmogorov complexity
3. **Zeta zero structural compressibility** — convergence accel closed, structural only
4. **Combinatorial pi(x) in O(x^{3/5})?** — open whether H-T trick exists for positive counts
5. Berry-Keating Hamiltonian, novel identity (long shots)

## Sub-agents Used: 7 (lit search [may still be running], TC^0, decompositions, Mertens, matrix powering research, plus main thread experiments)

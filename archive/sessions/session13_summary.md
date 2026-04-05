# Session 13 Summary

**Date:** 2026-04-04
**Focus:** TC^0 primality analysis, spectral graph approaches, algebraic shortcuts
**Sub-agents:** 5 background + multiple inline
**New closed paths:** ~15
**Breakthrough:** None (problem remains open)

---

## Key Findings

### 1. BPSW is Computationally in TC^0 (Most Important Finding)

ALL components of BPSW primality test verified to be in TC^0:
- **Miller-Rabin base 2:** scalar modular exponentiation = TC^0
- **Strong Lucas test:** each V_{d*2^r} is independent 2x2 MPOW (parallel), all in TC^0
- **Jacobi symbol:** in TC^0 via GCD (Hesse-Allender-Barrington 2002)
- **Parameter selection (Selfridge):** constant-time, trivially TC^0

**Implication:** "PRIMES in TC^0" now reduces precisely to proving BPSW is
unconditionally correct for all n. No BPSW pseudoprime is known below 2^64.

Also verified: QFT (Grantham's Quadratic Frobenius Test) is in TC^0.
Wilson's theorem and sum-of-two-squares approaches are CLOSED for TC^0.

### 2. Prime Indicator ANF Structure (Novel Measurement)

The prime indicator function over GF(2) has:
- **ANF degree = Theta(N)** (essentially full degree, like random)
- **50% of ANF terms nonzero** (exactly random sparsity)
- **Binomial degree distribution** (matches random expectation)

This confirms PRIMES not in AC^0[2] quantitatively and shows no GF(2)
algebraic shortcut exists for counting primes.

### 3. Zeta Zero Minimum Scaling

K_min (minimum zeros for exact pi(x)) scales as:
- **K_min ~ 0.35 * x^{0.27}** (power law, exponent ~1/4)
- **Greedy-optimal reordering** helps marginally but same asymptotic
- **No universal reordering** — optimal order depends on x

### 4. Spectral Graph Theory: Thoroughly Closed

Three angles tested:
- **Cayley graph:** CIRCULAR (needs primes to construct)
- **Ihara zeta function:** EQUIVALENT to Selberg trace / explicit formula
- **GCD/coprimality graph:** eigenvalues = Ramanujan sums, equivalent to Meissel-Lehmer O(x^{2/3})
- **Expander mixing:** either circular or equivalent

### 5. Other Closed Paths

- **Optimized li(x^{1/k}) coefficients:** Riemann R(x) is essentially optimal; error grows as x^{0.3}
- **CRT reconstruction:** each pi(x) mod q costs same as pi(x); strictly worse
- **Recursive identity pi(x) = F(pi(x/d)):** correction terms as hard as pi(x)
- **Prime zeta P(s) extraction:** more poles than explicit formula
- **Probabilistic counting:** requires Theta(x^2/ln(x)) samples for exactness
- **Hash-based counting:** partitions work, doesn't reduce it

### 6. Literature Search

No new algorithmic breakthroughs found through April 2026:
- "Feasibility of Primality in Bounded Arithmetic" (2504.17041): AKS provable in VTC^0_2
- "Primes via Zeros" (STOC 2025): about ideal primality, not integer counting
- Semenov smooth approximation: computationally O(n^2), useless
- No improvement to van der Corput exponential sum exponents

---

## Updated Viable Directions (for Session 14+)

1. **Prove BPSW correct** — would settle PRIMES in TC^0. Most concrete remaining question.
2. **Non-sieve, non-analytic pi(x)** — the fundamental open question. All known mathematical
   frameworks produce exponential circuits. A genuinely new framework is needed.
3. **Circuit complexity lower bound** — prove pi(x) requires super-polylog depth or super-poly size.
   Would settle the problem negatively.
4. **Kt complexity of delta(n)** — connects to circuit bounds via Oliveira (2019).
5. **Novel number-theoretic identity** — the "serendipity" direction.

## Files Created/Modified

### New experiments:
- experiments/analytic/optimized_li_combination.py
- experiments/analytic/zeta_zero_basis_experiment.py
- experiments/circuit_complexity/tc0_primality_approaches.py
- experiments/circuit_complexity/jacobi_symbol_tc0.py
- experiments/circuit_complexity/strong_lucas_tc0_depth.py
- experiments/circuit_complexity/prime_indicator_anf.py
- experiments/other/spectral_graph_primes.py
- experiments/other/crt_reconstruction_pi.py
- experiments/other/recursive_identity_search.py
- experiments/other/generating_function_extraction.py
- experiments/other/probabilistic_pi_circuit.py

### New novel findings:
- novel/prime_indicator_anf_structure.md

### Updated status files:
- CLAUDE.md (Session 13 insights, new closed paths)
- status/CLOSED_PATHS.md (+15 entries)
- status/OPEN_PROBLEMS.md (Session 13 notes)
- proven/complexity.md (TC^0 table)

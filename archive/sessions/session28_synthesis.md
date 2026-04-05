# Session 28 Synthesis: Circuit Complexity Deep Dive

**Date:** 2026-04-05
**Focus:** Novel experiments probing the circuit complexity of pi(x)
**Experiments:** 8+ (4 by sub-agents, 4+ direct)
**New closed paths:** ~7
**Novel result:** Approximate degree = N/2 universality

---

## Session Strategy

After 591+ approaches across 27 sessions, proposal-and-test cycles returned 100%
duplicates (Session 27 critique). Session 28 shifted strategy: instead of proposing
new algorithms, we probed the STRUCTURE of the barrier itself through novel circuit
complexity experiments. The goal was to find cracks in the exponential wall — or
confirm its solidity from new angles.

---

## Key Results

### 1. APPROXIMATE DEGREE = N/2 (NOVEL, most significant)

**Files:** `experiments/circuit_complexity/approx_degree_small.py`, `approx_degree_counting.py`
**Writeup:** `novel/approx_degree_prime.md`

The epsilon-approximate degree of the prime indicator at the rounding threshold
(epsilon = 0.49) equals ceil(N/2) for N = 4..10. Crucially:

- **adeg(chi_P) = adeg(pi(x) mod 2)**: The counting step adds NO difficulty
- **N/2 matches the oscillatory correction**: The first N/2 degrees capture R(x)
- **Quantum lower bound Omega(N/4)**: Still polylog in x = 2^N
- **Does NOT rule out polylog circuits**: Approximate degree bounds queries, not gates

This is the cleanest quantitative characterization of the smooth/oscillatory boundary.

### 2. THE "N/2" UNIVERSALITY (NOVEL synthesis)

Every complexity measure we've tested converges at the SAME threshold of N/2:

| Measure | Value at N/2 | Source |
|---------|-------------|--------|
| Approximate degree | Rounding threshold | Session 28 |
| Communication rank | 2^{N/2-1}+2 vs 2^N total | Session 17 |
| Oscillatory error | O(2^{N/2}) | Theory |
| Per-bit R-correlation | Crossover to <0.5 | Session 28 |
| LFSR complexity of delta(n) | N/2 over all GF(p) | Sessions 24, 26 |
| Per-bit influence crossover | ~bit N/2 | Session 28 |
| Fourier spectral weight | Drops below 0.5 at degree ~N/2 | Session 28 |

This universality means the smooth/oscillatory decomposition is not an artifact
of any particular analytical framework — it is a FUNDAMENTAL structural feature
of the prime counting function.

### 3. PER-BIT COMPLEXITY GRADIENT (NOVEL)

**File:** `experiments/circuit_complexity/per_bit_complexity.py`

Output bits of pi(x) form a clear complexity spectrum:
- **MSBs (top ~N/2):** Influence O(1), R-correlation >0.95, Fourier weight >90% at degree ≤2
- **LSBs (bottom ~N/2):** Influence ~N/2, R-correlation <0.5, Fourier weight <10% at degree ≤2
- **Ratio grows:** LSB/MSB influence ratio increases from 1.4 (N=4) to 2.1+ (N=14)

This means the top half of pi(x)'s output is essentially free (computable via R(x)
in polylog time), while the bottom half encodes the oscillatory correction.

### 4. ROUNDING BOUNDARY ANALYSIS

**File:** `experiments/circuit_complexity/rounding_boundary_analysis.py`

- frac(R(x)) is perfectly uniform (chi-squared non-significant for all N)
- Precision requirement follows EXACT geometric distribution: P(k bits) = 2^{-(k-1)}
- No "easy subset" of inputs exists — errors are uniformly distributed
- R(x) accuracy drops from 65% (N=8) to 13% (N=16)

### 5. WHEEL DECOMPOSITION IS ANTI-WIN

**File:** `experiments/circuit_complexity/modular_counting_attack.py`

Per-class counting functions pi_r(x) are 3-4x MORE complex (per bit) than full pi(x).
The sequential regularity of consecutive integers — that p and p+2 can't both be even — 
is DESTROYED by M-spacing. Classes are nearly independent (I/H < 0.01), which prevents
cross-class shortcuts. Total entropy grows as phi(M).

**Key insight:** The regularity of pi(x) is TOPOLOGICAL (sequential ordering), not
ALGEBRAIC (modular structure). Any decomposition that breaks sequential adjacency
increases complexity per subproblem.

### 6. MULTIPLICATIVE-ADDITIVE BARRIER CONFIRMED FROM THREE ANGLES

**File:** `experiments/circuit_complexity/multiplicative_circuit_structure.py`

Three independent analyses of the Legendre sieve:
1. **I-E signed rank = #distinct floors:** No algebraic shortcut for the sum
2. **Carry chains match random:** No additive structure exploitable for counting
3. **Rectangle partition ~ 2^{0.76*N}:** Nondeterministic communication exponential

### 7. LITERATURE SEARCH: NO BREAKTHROUGHS

Searched arxiv, ECCC, Google Scholar for 2025-2026 papers:
- **Chen-Tal-Wang (ECCC TR26-039):** THR∘THR n^{2.5-eps} lower bounds, but for E^NP functions, not number-theoretic
- **Raz et al. (ECCC TR26-008):** Natural proofs barrier for linear functions over finite fields
- **No new pi(x) algorithms, no new #TC^0 ⊆ NC results, no new zeta zero summation methods**

---

## What Changed This Session

### Strengthened:
- The "N/2" boundary is now confirmed by 7+ independent measures (was 3-4)
- Approximate degree provides the cleanest characterization yet
- The "sequential regularity" insight explains why decomposition fails

### Unchanged:
- The problem is still OPEN — no proof of impossibility, no proof of possibility
- #TC^0 ⊆ NC? remains THE key question (no progress toward resolving it)
- All sieve-based approaches remain closed
- The Natural Proofs barrier still blocks impossibility proofs

### Novel contributions:
1. **adeg(chi_P, 0.49) = N/2** — first measurement of approximate degree for prime indicator
2. **"N/2 universality"** — all measures converge at the same threshold
3. **Per-bit complexity gradient** — quantitative smooth/oscillatory decomposition per output bit
4. **Sequential regularity insight** — topological vs algebraic structure explanation

---

## Recommendations for Future Sessions

1. **Do NOT generate more proposals.** The space is exhausted (100% duplicate rate).

2. **Focus on #TC^0 ⊆ NC?** This is the ONLY open theoretical question that could
   resolve the problem. Requires new circuit complexity techniques, not number theory.

3. **Monitor literature** for:
   - Progress on TC^0 vs NC^1 separation
   - New counting complexity results (#TC^0, #AC^0)
   - Breakthroughs in L-function computation or zeta zero algorithms
   - New algebraic geometry / Langlands program connections to primes

4. **Consider impossibility:** Even partial results (e.g., "pi(x) has no circuits of
   size N^c for c < 2") would be valuable, despite the Natural Proofs barrier.

5. **The approximate degree result could be publishable.** The N/2 threshold and its
   universality across measures is a clean, novel characterization that the
   analytic number theory and computational complexity communities would find interesting.

---

## Updated Approach Count

**Total approaches tested: ~600+ across 28 sessions.**
**Total novel findings: 5 (info-computation gap, failure taxonomy, communication rank
formula, N/2 universality, approximate degree = N/2)**

### Late-arriving agent results:

**(k) BDD synthesis (NOVEL measurement):** Multi-ordering BDD for LSB of pi(x) scales as
2^{0.73*N}, between OBDD (0.79) and communication rank (0.50). MSB BDD ~ N+1 (trivial).
BDD ≠ circuit, so doesn't rule out polynomial circuits.

**(l) Polynomial method CLOSED:** Power-law fit 0.601*N^{0.909} = Theta(N). Sharp phase
transition at degree N/2 (error exactly 0.5 below, exponential decay above). Promise
problem reduces adeg by only 1-2. SOS degree = adeg. Rules out polynomial method,
quantum query speedups, SOS relaxations. Arithmetic algorithms not constrained.

**Complexity hierarchy for LSB of pi(x):**

| Model | Exponent beta (2^{beta*N}) | Source |
|-------|---------------------------|--------|
| Communication rank | 0.50 | S17 |
| BDD (multi-order) | 0.73 | S28 |
| OBDD (single-order) | 0.79 | S20 |
| ANF density | ~1.0 | S13 |
| General circuits | ??? (OPEN) | — |

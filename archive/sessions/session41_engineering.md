# Session 41: Engineering Improvements + Novel Write-up (2026-04-05)

## Theme
All FOCUS_QUEUE tasks completed; pivoted to engineering improvements to v10
and documenting novel findings per CLAUDE.md guidelines.

## Experiments

### 1. Lehmer Method vs Lucy DP (`experiments/sieve/lehmer_vs_lucy.py`)
**Idea:** Implement Meissel-Lehmer decomposition pi(x) = phi(x,a) + a - 1 - P2 - P3
in C, benchmark against naive Lucy DP.

**Result:** Lucy DP is 5-13x FASTER than Lehmer at all tested x (10^5 to 10^10).
Recursive phi with hash memoization has worse cache behavior than Lucy's sequential
array scan. The Meissel-Lehmer variant (with P3 and shallower phi) is slightly better
than basic Lehmer but still 7-10x slower than Lucy.

**Verdict:** FAIL (Equivalence). Naive Lehmer decomposition doesn't help without
segmented sieve infrastructure.

### 2. Wheel-30 + Prefetch Lucy DP (`experiments/sieve/segmented_lucy.py`)
**Idea:** Skip p=2,3,5 iterations via wheel mod 30 initialization; add prefetch hints.

**Result:** 2-35% speedup only. Marginal and diminishing at larger x.

**Verdict:** PARTIAL. Technically works but not significant.

## Novel Finding Documented
- **Carry-propagation boundary** (`novel/carry_propagation_boundary.md`): Per-bit
  difficulty landscape of p(n) has sharp 4-bit sigmoid at ~60% of bit length.
  Novel quantitative measurement from Session 40.

## Literature Check
Web search found no new breakthroughs in prime counting algorithms, TC^0 lower
bounds, or zeta zero computation in March-April 2026.

## Key Insight
The 100-1000x gap between v10 and primecount is entirely algorithmic (Gourdon's
special leaves decomposition + segmented sieve + SIMD + OpenMP), not parametric.
Within-Lucy-DP optimizations are marginal. Matching primecount = reimplementing
primecount.

## Status
533+ approaches closed. No breakthrough. Project remains in monitoring/engineering phase.

# Session 33: Deep Focus — Conditional Algorithms

**Date:** 2026-04-05
**Task:** FOCUS_QUEUE Task #4
**Method:** 4 parallel sub-agents, 7 experiments
**Result:** All 8 paths closed. Best conditional = O(x^{1/2+ε}) under RH + batch zeros.

## Experiments

### Exp 1: GRH Batch Miller Testing
- GRH bound 2ln²(n) is conservative: only {2,3} needed to 10^6
- Batch testing: 1.02x speedup (negligible)
- Cost O(n·polylog n), worse than all known counting methods

### Exp 2: GRH Explicit Formula Optimal T
- T_min = O(√x·log²x) confirmed numerically
- 1000 zeros: error ≥18 for x≥10^4
- Individual zeros at γ~200-400 contribute >0.5 to pi(10^6)
- Cannot use fewer than O(√x) zeros

### Exp 3: Elliott-Halberstam and Counting
- EH controls equidistribution, NOT total count
- Residue-class summing is tautological (exact) or li(x) (approximate)
- 10^3-10^6x slower than direct computation

### Exp 4: Gap Structure (Cramér)
- Cramér model validates (max gaps 40-60% of ln²p)
- Search interval has O(ln x) primes
- Bottleneck: identifying nth prime requires pi(x) at boundary

### Exp 5: Schoenfeld Explicit Bounds under RH
- Interval sieve is O(x^{1/2+ε}), better than ML O(x^{2/3})
- 1000 zeros: 21% error reduction at x=10^6
- Need O(√x·log²x) zeros for error < 1

### Exp 6: Cramér Search Algorithm
- Walk phase: O(ln⁴x), 0-114 steps
- Counting phase: O(x^{2/3}), dominates by factor 10^57 at x=10^100

### Exp 7: Best Conditional Algorithm Benchmark
- RH + Turing zeros: O(x^{2/3+ε}) — WORSE than unconditional
- RH + Odlyzko-Schönhage: O(x^{1/2+ε}) — best known
- GRH/EH/Cramér don't help beyond O(x^{1/2+ε})
- For p(10^100): ~10^51 ops minimum under all conjectures

## Key Insight
The √x barrier for exact pi(x) is fundamental across ALL standard conjectures.
Each conjecture addresses the wrong subproblem:
- GRH: primality testing (already cheap)
- EH: distribution across residue classes (not total count)
- Cramér: gap structure (search, not counting)
The irreducible cost is counting primes, which encodes O(√x) zeta zeros.

## Files
- experiments/analytic/conditional/grh_miller_batch.py + _results.md
- experiments/analytic/conditional/elliott_halberstam_gaps.py + _results.md
- experiments/analytic/conditional/schoenfeld_cramer.py + _results.md
- experiments/analytic/conditional/best_conditional_algorithm.py + _results.md

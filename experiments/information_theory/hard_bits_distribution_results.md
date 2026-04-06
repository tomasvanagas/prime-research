# Hard Bits Distribution Analysis — Results

**Script:** hard_bits_distribution.py
**Session:** 41

## What Was Tested
Whether the "hard" bottom ~50% of bits of p(n) have any statistical structure 
that could prioritize a bruteforce search. 9 analyses on 500K primes.

## Key Findings

### Finding 1: Delta is ALWAYS POSITIVE for this approximation
delta(n) = p(n) - round(R^{-1}(n)) is positive 100% of the time with this 
Cipolla-style approximation. Mean = +1575, range [34, 3358]. This is a systematic
bias in the approximation — R^{-1} undercounts. Using true R^{-1}(n) via mpmath
would center the distribution closer to zero with both signs.

### Finding 2: Per-bit bias is ZERO
Every hard bit (1-10) has P(1) = 0.500 within noise (max z-score = 0.50).
No bit position is biased. No individual bit gives a search advantage.

### Finding 3: Bit pairs have ZERO mutual information
0 out of 55 hard-bit pairs have MI > 0.001. Hard bits are pairwise independent.
No joint pattern exists to exploit.

### Finding 4: Delta mod anything is UNIFORM
chi2 values are tiny for all moduli tested (2,3,4,5,6,8,10,12,16,30).
No residue class is favored. Delta is equidistributed modulo every small number.

### Finding 5: Conditioning on n mod 30 gives ZERO reduction
Unconditional std = 676.1, conditional std = 676.1. Knowing n mod 30 does not 
narrow the delta distribution at all. 0.00% reduction.

### Finding 6: Even/odd delta is 50/50
Even: 50.04%, Odd: 49.96%. No parity bias.

### Finding 7: Frequency-based search saves ~32% over naive
- Naive (search outward from 0): mean 1576 candidates
- Frequency-ordered: mean 1070 candidates (32% fewer)
- But searching only PRIMES near R^{-1}: only ~10.4 candidates on average

### Finding 8: |delta| scales as ~p^{0.4} (sublinear)
Not √p (0.5) but closer to p^{0.4}. The search window grows slower than √p.

## Verdict
**CLOSED — No statistical edge in the hard bits**

The hard bits of p(n) are perfectly uniform: no per-bit bias, no pairwise 
correlations, no modular bias, no conditioning helps. The distribution of delta 
is smooth and featureless. The only "edge" is searching primes only (not all 
integers), which reduces candidates from ~1575 to ~10 — but this is just the 
prime density 1/ln(p), not structure in the hard bits.

## One-Line Summary
Hard bits of p(n) are statistically uniform in every measure tested; only prime 
density (1/ln p) reduces search space.

# Index-Prime Bit Correlation Analysis — Results

**Script:** index_prime_bit_correlation.py
**Session:** 41

## What Was Tested
Cross-correlation between bits of the index n and bits of p(n) for 664,579 primes
up to 10^7 (24-bit). Eight analyses: MI matrix, conditional entropy, Pearson
correlation, same-position agreement, XOR analysis, shifted correlation, multi-bit
prediction, and residual bit analysis.

## Key Findings

### Finding 1: Dominant shift of +4 bits (p(n) ≈ n × 2^3.7)
The strongest correlation is between bit_j(n) and bit_{j+4}(p(n)):
- Shift +4: mean |r| = 0.131, **max |r| = 0.872** (n-bit 18 vs p-bit 22)
- All other shifts: mean |r| < 0.015
This reflects p(n) ≈ n·ln(n) where ln(n) ≈ 13.4 ≈ 2^{3.7} for this range.
Pure PNT relationship, not exploitable structure.

### Finding 2: Only top 4-5 bits of p(n) are correlated with n
Percentage of each p-bit explained by ALL bits of n:
- p-bits 0-12: 0.00% (completely independent of n)
- p-bit 15: 0.27%
- p-bit 18: 5.4%
- p-bit 19: 19.8%
- p-bit 21: 30.4%
- p-bit 22: 72.6%
The bottom ~75% of bits of p(n) carry zero information from the raw index n.

### Finding 3: Same-position correlation is zero
bit_k(n) vs bit_k(p(n)): correlation < 0.001 for all k < 16.
Agreement rate = 50.0% (exact coin flip). No same-position structure.

### Finding 4: Multi-bit prediction fails for low bits
Best 3-bit n-groups predicting p-bits:
- p-bits 0-8: accuracy 50.0-50.7% (noise above 50% baseline)
- p-bit 15: 54.1% (4% lift)
- p-bit 19: 80.9% (30% lift)
Only the very top bits are predictable.

### Finding 5: Residual bit-length ratio is ~1.04
p(n) - n·ln(n) has almost the same bit-length as n. The simplest approximation
error is as wide as the input itself.

### Finding 6: bit_length(p(n)) - bit_length(n) = 3-4
Consistent with p(n) ≈ n·ln(n) where ln(n) adds ~3.7 bits.
diff=+3: 20.3%, diff=+4: 79.6%.

## Comparison with R^{-1}(n) results (Session 40)
- R^{-1}(n) predicts ~50% of bits of p(n) correctly
- Raw n predicts only ~20% of bits (top 4-5 of 24)
- R^{-1}(n) is much more informative because it incorporates ln, ln(ln), etc.
- Both converge on the same barrier: bottom bits are unpredictable

## Verdict
**CLOSED — No exploitable cross-correlation between bits of n and bits of p(n)**

**Failure Mode:** I (Information Loss) — The bits of n carry almost no information
about the low bits of p(n). The only correlation is the trivial PNT shift (p(n) ≈
n·ln(n)). Even the best 3-bit lookup table from n cannot beat 50.7% on any p-bit
below position 8. The N/2 barrier is even more severe when measured from raw n
rather than R^{-1}(n).

## One-Line Summary
Bits of n predict only top ~4 bits of p(n) via PNT shift; bottom 75% of p-bits have
zero mutual information with n.

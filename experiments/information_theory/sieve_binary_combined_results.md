# Sieve + Binary Combined Analysis — Results

**Script:** sieve_binary_combined.py  
**Session:** 41

## What Was Tested
Whether combining R^{-1}(n) (top bits) with wheel mod 30030 (bottom bits) 
squeezes extra determined bits from the middle.

## Key Findings

### The wheel only constrains 2.4 bits
log2(30030) - log2(5760) = 14.87 - 12.49 = **2.38 bits**. Despite eliminating 
81% of integers, this is only 2.4 bits of information. And those 2.4 bits are:
- Bit 0: always 1 (odd) — 1 bit
- Bit 15: always 0 (residues < 30030) — not useful for large primes  
- Remaining ~1.4 bits spread across higher bits as tiny biases

### Bits 1-10 are COMPLETELY UNKNOWN by either method
Per-bit analysis shows a clear 3-zone structure:
- Bits 14-23 (MSB): determined by R^{-1} (>96% correct) → **R^{-1} ZONE**
- Bits 11-13: partially determined (73-93% correct) → **TRANSITION**
- Bits 1-10: 50% accuracy (pure coin flip) → **DEAD ZONE**
- Bit 0: constrained by wheel (always 1) → **WHEEL**

### R^{-1} CANNOT predict the wheel residue class
R^{-1}(n) mod 30030 matches p(n) mod 30030 in **0 out of 200,000** cases.
The residue difference has std = 433, which is much larger than the wheel period.
The two sources of information are completely non-overlapping but also non-combining.

### Bottom 4 bits of coprime residues are PERFECTLY UNIFORM
All 8 odd bottom-4-bit patterns appear exactly 720 times each among the 5760 
coprime residues. The wheel provides zero information about bits 1-3.

### Combined search space for 1000-bit prime
- R^{-1}: ~500 bits determined
- Wheel: ~2.4 bits constrained  
- Unknown: ~498 bits → 2^498 candidates → 10^146 prime candidates
- The wheel saves 2.4 bits out of 498 unknown — a 0.5% improvement.

## Verdict
**CLOSED — Wheel adds negligible bits to R^{-1}; no middle-bit shortcut**

The wheel and R^{-1} target opposite ends of the bit string with no overlap 
and no synergy in the middle. The "dead zone" (bits 1-10 at 24-bit scale, 
or bits 1-498 at 1000-bit scale) is impervious to both. Combining them saves 
only 2.4 bits regardless of prime size.

## One-Line Summary
Wheel mod 30030 constrains only 2.4 bits (bit 0 + tiny biases); R^{-1} and wheel 
target opposite ends with no middle-bit synergy; dead zone remains ~N/2 bits.

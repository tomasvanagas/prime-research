# Binary Bit Analysis of Consecutive Primes — Results

**Script:** binary_bit_analysis.py
**Session:** 41

## What Was Tested

9 binary/bitwise analyses of primes up to 10^7 (664,579 primes, 24-bit):
1. XOR structure between consecutive primes
2. Hamming distance patterns
3. Per-bit-position P(bit=1) statistics
4. Run-length structure in binary representations
5. Bit transition matrices (Markov per bit)
6. Spectral analysis per bit position (FFT)
7. Inter-bit mutual information within a single prime
8. XOR predictability (autocorrelation of bit flips)
9. Gap-conditioned bit flip analysis

## Key Findings

### Finding 1: Bits within a prime are INDEPENDENT (MI = 0.0000 for all pairs)
No bit position carries information about any other bit position within the same prime.
This rules out any "infer bit k from bits 0..k-1" strategy.

### Finding 2: XOR bit flips show significant autocorrelation (bits 1-11)
The most interesting result. Bit flip sequences have structure:
- Bit 5: lag-1 autocorrelation = -0.281 (strong anti-persistence)
- Bit 6: lag-1 = -0.245, lag-2 = -0.178 (multi-lag structure)
- Bit 7: lag-1 through lag-5 all ~ -0.12 (long-range negative correlation)

Interpretation: if a bit flipped last step, it's LESS likely to flip next step.

### Finding 3: Gap determines bit flips completely
Gap=2: exactly bit 1 flips (binary: ...10)
Gap=4: exactly bit 2 flips (binary: ...100)
Gap=6: bits 1+2 flip with rates matching binary addition of 6=110
This is just carry propagation arithmetic — not exploitable structure.

### Finding 4: Hamming distances are anti-correlated (lag-1 = -0.204)
After a large bit change, the next change tends to be small. This is a
consequence of prime gap anti-correlation (Gallagher-type effect).

### Finding 5: Per-bit P(1) = 0.500 for bits 1-14 (no bias)
Bits are individually unbiased. High bits (15+) show bias from finite range.

### Finding 6: Spectral flatness drops sharply with bit position
- Bits 1-4: flatness ~0.55 (fairly random)
- Bit 7: flatness 0.17
- Bit 10: flatness 0.03
- Bit 15: flatness 0.001
High bits change slowly (long correlations), low bits change chaotically.

### Finding 7: Run lengths slightly shorter than random geometric
Both 0-runs and 1-runs have mean ~1.9 vs expected 2.0. Long runs (>5) are
underrepresented by factor 0.8-0.6x. Small but real deviation from random.

### Finding 8: Bit 0 (LSB) is always 1
Trivial: all primes > 2 are odd.

## Verdict
**CLOSED — No exploitable structure for algorithm improvement**

**Failure Mode:** The XOR autocorrelation (Finding 2) and gap conditioning (Finding 3) 
are real structure, but they are CONSEQUENCES of the prime gap distribution, not 
independent information. The gap sequence g(n) = p(n+1) - p(n) determines all bit 
changes via binary addition carry propagation. The anti-correlation in Hamming distances 
reflects known gap anti-correlation. None of this helps compute p(n) from n alone — it 
only helps predict p(n+1) from p(n), which is the trivial sequential direction.

The inter-bit MI = 0 result (Finding 1) is the strongest negative result: it proves 
bits within a single prime are statistically independent, ruling out any strategy that 
recovers hard bits from easy bits.

## One-Line Summary
Binary bit analysis: bits within primes are independent (MI=0); bit flips between 
consecutive primes are determined by gaps via carry propagation; no exploitable structure.

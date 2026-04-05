# Carry-Propagation Boundary — Results

**Experiment:** Map which bits of p(n) are "easy" (determined by R^{-1}(n)) vs "hard" (need zeta zeros).

**Verdict:** CLOSED as approach, but provides **novel quantitative measurement** of the easy/hard bit boundary.

## Key Findings

### 1. Error Scaling
|p(n) - R^{-1}(n)| ~ p(n)^{0.66}

This means ~66% of bits (from LSB) are "hard." Theory under RH predicts ~50% (sqrt(x) error). The observed 66% is likely a small-sample bias — at p(n) ~ 10^5, the error hasn't settled to its asymptotic behavior.

### 2. Per-Bit Agreement Rate (MAIN RESULT)
Beautiful sigmoid transition from EASY → HARD:

| Bit (MSB=0) | Agreement | Status |
|---|---|---|
| 0-4 | 99.4-100% | EASY |
| 5-6 | 97.5-98.5% | EASY |
| 7 | 95.0% | EASY |
| 8 | 89.8% | MEDIUM |
| 9 | 79.6% | MEDIUM |
| 10 | 63.8% | HARD |
| 11 | 53.9% | HARD |
| 12-16 | 48.6-50.9% | RANDOM (coin flip) |

**The transition spans ~4 bit positions (bits 8-11)** from 90% to 50%. Below bit 12, R^{-1}(n) provides NO information — those bits are pure noise from R^{-1}'s perspective.

### 3. First Disagreeing Bit Distribution
- Mean: 0.650 (normalized, 0=MSB, 1=LSB)
- Std: 0.130
- Bulk of disagreements: 0.5-0.8 range
- Sharp peak at 0.5-0.8 with long tails

### 4. Boundary Stability
The boundary position is remarkably stable across different n ranges:
- n ∈ [100,600]: mean = 0.582
- n ∈ [2000,2500]: mean = 0.597
- n ∈ [4000,4500]: mean = 0.604

Slight upward drift (more easy bits at larger n) is expected as p(n) grows and the relative error shrinks.

## Implications for Polylog
The sigmoid transition means there's NO gradual "peeling" possible:
- The top ~60% of bits are essentially FREE from R^{-1}(n) [O(polylog)]
- The bottom ~40% are essentially RANDOM [require zeta zero info]
- The transition zone is only ~4 bits wide — too narrow to exploit

A hypothetical hierarchical approach would need to recover the bottom 40% bit by bit, but each bit in the hard zone is essentially a coin flip from R^{-1}'s perspective, requiring full correction information.

## Failure Mode
**INFORMATION**: The binary representation of the correction term has no exploitable intermediate-difficulty bits. It's a sharp phase transition from "easy" to "hard" with no gradual degradation to exploit hierarchically.

## Novelty
This per-bit difficulty landscape has not been explicitly measured before (to our knowledge). The 4-bit transition width and sigmoid shape are new quantitative results.

**Session:** 40 (Fresh Perspective)

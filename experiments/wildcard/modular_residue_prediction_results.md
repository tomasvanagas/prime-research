# Modular Residue Prediction for p(n) — Results

**Experiment:** Can p(n) mod q be predicted from n alone? If so, CRT would recover p(n) exactly.

**Verdict:** CLOSED — FAILS. Maximum entropy, circular reasoning.

## Key Findings

1. **Near-maximum entropy**: For all q > 2, p(n) mod q has entropy ratio > 0.999 of the theoretical maximum log2(phi(q)). The residues are effectively uniformly distributed.

2. **Prediction accuracy near random baseline**:
   | q | 1st-order accuracy | Random baseline | Improvement |
   |---|---|---|---|
   | 3 | 0.598 | 0.500 | 1.20x |
   | 5 | 0.354 | 0.250 | 1.42x |
   | 7 | 0.269 | 0.167 | 1.61x |
   | 11 | 0.229 | 0.100 | 2.29x |
   | 210 | 0.403 | 0.021 | **19.2x** |

3. **Interesting: mod 210 prediction is 19x above random** — because the primorial wheel (210 = 2·3·5·7) creates strong correlations between consecutive primes in the wheel pattern. The autocorrelation at lag 1 is 0.73 for mod 210.

4. **Cross-correlations near zero**: p(n) mod q for different coprime q are effectively independent (all |corr| < 0.005).

5. **CRT is circular**: Computing p(n) mod q exactly requires knowing pi(x; q, a) exactly near x = p(n), which is as hard as computing pi(x) exactly.

## The Primorial Wheel Effect
The high prediction accuracy for mod 210 is NOT a useful signal — it just reflects that consecutive primes must avoid multiples of 2, 3, 5, 7 within the wheel. This constrains which residue classes are possible but doesn't help compute p(n) because:
- The wheel pattern is the same for ALL primes > 7
- The actual residue within the wheel is still unpredictable

## Failure Mode
**CIRCULARITY + INFORMATION**: The CRT approach reduces to computing pi(x; q, a) for each modulus q, which is no easier than pi(x). The near-maximum entropy of residues means no shortcut exists through modular arithmetic.

**Session:** 40 (Fresh Perspective)

# Ulam Spiral Deep Analysis — Results

**Script:** ulam_spiral_analysis.py
**Session:** 41

## What Was Tested
Whether the Ulam spiral's diagonal patterns (quadratic polynomials that generate
primes at high density) carry exploitable information beyond modular arithmetic.
8 analyses on primes up to 2,000,000.

## Key Findings

### Finding 1: Euler's n²+n+41 is 7.4x denser than baseline — but covers only 0.5% of primes
Best quadratic generates 58% primes (vs ~8% baseline), but 10 quadratics combined
cover only 1.86% of all primes up to 2M. Quadratic forms are too sparse to be useful.

### Finding 2: Ulam diagonals are 2.4-2.8x denser — modest
The spiral's SE/NW/NE diagonals have prime density 2.4-2.8x above baseline.
Visually striking but mathematically modest. The SW diagonal (4n²+4n+2 = all even)
has zero primes — this is the trivial mod-2 pattern that creates the visual asymmetry.

### Finding 3: Information content is tiny
- MI(n mod 30; is_prime) = 0.174 bits (wheel sieve — 42% of total entropy)
- MI(n mod 4; is_prime) = 0.088 bits (spiral diagonal class)
- MI(Euler_poly_prime; is_prime(n)) = 0.000349 bits (essentially zero)
The spiral carries LESS than half the information of simple mod-30 wheel.

### Finding 4: After removing wheel-30, spiral structure vanishes
Among wheel-30 survivors (coprime to 30):
- n mod 4 lift: ±0.131% (noise level)
- n mod 120: max lift 0.00182 (below noise threshold 0.00491)
- n mod 210: max lift 0.026 (borderline, at residue 49)
The spiral's visual pattern is ENTIRELY explained by mod-30 structure.

### Finding 5: Quadratics help narrow R⁻¹ search only 44% of the time
When they help, the improvement is modest. Mean R⁻¹ error is 250; mean distance
to nearest Euler value is 369. Quadratics are WORSE than R⁻¹ as locators.

### Finding 6: Hardy-Littlewood constants are accurate
C_f correctly predicts relative density ratios between polynomials.
This confirms the patterns are well-understood (since 1923) and carry no
information beyond what the Hardy-Littlewood conjecture already describes.

## Verdict
**CLOSED — Ulam spiral = modular arithmetic in visual form, no new information**

**Failure Mode:** E (Equivalence) — The spiral's diagonal patterns are equivalent to
n mod 4 structure, which is a subset of wheel-30 factorization. After removing
mod-30 information, zero residual structure remains. Hardy-Littlewood (1923) fully
explains the density variations. 10 best quadratics cover only 1.86% of primes.

## One-Line Summary
Ulam spiral patterns carry 0.088 bits (mod 4), strictly less than wheel-30 (0.174 bits);
after wheel removal, zero residual structure; 10 quadratics cover only 1.86% of primes.

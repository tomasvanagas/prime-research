# Precomputation + Interpolation Attack: Results

**Script:** interpolation_attack.py

## What Was Tested
Whether precomputed pi(x) or p(n) at strategic checkpoints enables O(polylog) interpolation. Five approaches: segmented explicit formula with power-of-2 checkpoints, Lagrange polynomial interpolation, spline interpolation with error bounds, block-based pi(x) structure, and smooth approximation + correction table.

## Key Findings
- Polynomial (Lagrange) interpolation of p(n) through k points has error that grows with the gap between checkpoints; for checkpoints at distance D, interpolation error is O(D^k * max derivative)
- Spline interpolation achieves better local accuracy but the correction delta(n) is not smooth enough for high-order splines
- Block-based pi(x) with precomputed checkpoints reduces work to O(block_size) per query, but block_size must be O(sqrt(x)) for exact results
- The delta correction table has ~5 bits/entry entropy and is incompressible
- Precomputation helps constant factors but cannot change the asymptotic complexity class

## Verdict
**CLOSED** -- Failure Mode: I (Information Loss)

## One-Line Summary
Interpolation between precomputed checkpoints fails because the prime correction delta(n) is not smooth enough to interpolate exactly.

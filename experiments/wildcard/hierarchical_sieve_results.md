# Hierarchical Sieve: Results

**Script:** hierarchical_sieve.py

## What Was Tested
Five FMM-inspired sieve decomposition approaches: (1) dyadic decomposition pi(x)-pi(x/2), (2) telescoping with smooth approximations, (3) multi-scale sieve with negligibility thresholds, (4) recursive identity truncation (inclusion-exclusion), (5) prime gap Fourier structure.

## Key Findings
- Dyadic decomposition: pi(x) - pi(x/2) still requires O(x^{2/3}) to compute each piece -- no cancellation between pieces
- Telescoping: delta(x) - delta(x/2) does NOT simplify; the oscillatory parts do not cancel between dyadic intervals
- Multi-scale sieve: rough-number correction never becomes negligible at any scale; it carries O(x^{1/2}) information at every level
- Recursive identity truncation: need O(x^{1/3}) inclusion-exclusion terms for exact count -- matches Meissel-Lehmer theory
- Prime gap partial sums: Fourier structure is broadband; 99% energy requires 78%+ of coefficients (from Session 24)

## Verdict
**CLOSED**
**Failure Mode:** Equivalence (E) -- hierarchical decomposition reproduces Meissel-Lehmer at every level; no cancellation or negligibility exploitable.

## One-Line Summary
Hierarchical sieve (dyadic, telescoping, multi-scale, truncation, gap Fourier): all reduce to O(x^{2/3}) Meissel-Lehmer; no level simplifies.

# Deep Structure Experiments on Prime Sequence: Results

**Script:** deep_structure.py

## What Was Tested
Eight experiments probing hidden patterns in p(n): base representation analysis (binary, primorial bases), digit distribution, mixed-radix representations, and structural pattern searches across ~1M primes up to 15,000,000.

## Key Findings
- Primorial mixed-radix representations of primes show no exploitable pattern beyond the known wheel structure
- Base representations (binary, ternary, etc.) of primes have digit distributions indistinguishable from random integers with the same density constraint
- No hidden periodicity or low-complexity structure detected in any base representation
- All digit-level analyses confirm the ~5 bits/prime entropy barrier

## Verdict
**CLOSED** -- Failure Mode: I (Information Loss)

## One-Line Summary
Prime representations in various bases show no exploitable structure beyond the smooth density approximation.

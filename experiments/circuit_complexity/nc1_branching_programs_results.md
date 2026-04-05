# NC^1 Branching Program Analysis for Primality: Results

**Script:** `nc1_branching_programs.py`
**Session:** 17

## What Was Tested
Whether prime indicator or pi(x) shows structure amenable to NC^1 computation via branching programs. Analyzed minimum OBP lengths for widths 2-5, communication complexity, modular exponentiation structure, and Fourier/sensitivity measures for N=4,6,8,10.

## Key Findings
- Minimum branching program length for chi_P matches random functions of same density at all widths tested
- Width-5 OBP (Barrington's NC^1 threshold) requires length ~2^N/N, same as random
- Communication matrix for balanced partition has full or near-full rank
- Modular exponentiation 2^{n-1} mod n shows no branching program structure
- Sensitivity and certificate complexity are maximal (= N)
- Primality behaves like a random function for branching program complexity

## Verdict
**CLOSED**
**Failure Mode:** Information loss (branching program complexity matches random functions; no NC^1 exploitable structure)

## One-Line Summary
Primality has random-like branching program complexity; OBP length matches random functions at all widths.

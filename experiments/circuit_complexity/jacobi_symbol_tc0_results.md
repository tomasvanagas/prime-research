# Jacobi Symbol in TC^0: Results

**Script:** `jacobi_symbol_tc0.py`
**Session:** 13
**See also:** `tc0_primality_approaches_results.md`, `strong_lucas_tc0_depth_results.md`

## What Was Tested
Whether the Jacobi symbol (a/n) is computable in TC^0. This is critical for the Lucas/Frobenius TC^0 primality path, since the strong Lucas test requires computing (D/n) to select parameters.

## Key Findings
- Standard Jacobi computation uses quadratic reciprocity reductions (GCD-like), which is O(log n) sequential steps
- Allender (1999) and HAB (2002) showed GCD is in TC^0 (via integer division in TC^0)
- Therefore Jacobi symbol IS in TC^0: the quadratic reciprocity reduction is a sequence of divisions and bit tests, all in TC^0
- This removes a potential bottleneck for the BPSW-in-TC^0 chain

## Verdict
**CLOSED** -- resolved positively. Jacobi symbol is in TC^0.
**Failure Mode:** N/A (positive result, not a barrier)

## One-Line Summary
Jacobi symbol is in TC^0 (via HAB 2002 GCD result); not a bottleneck for TC^0 primality testing.

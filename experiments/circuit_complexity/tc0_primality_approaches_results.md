# TC^0 Primality Approaches (Non-AKS): Results

**Script:** `tc0_primality_approaches.py`
**Session:** 13

## What Was Tested
Four non-AKS approaches to placing PRIMES in TC^0: Wilson's theorem ((n-1)! mod n), Lucas sequences (2x2 MPOW), quadratic Frobenius test (2x2 MPOW), and sum-of-two-squares (Fermat characterization).

## Key Findings
- Wilson's theorem: (n-1)! mod n requires computing a product of n terms; iterated multiplication of n terms is NOT in TC^0 when n is the input (n = 2^N)
- Lucas sequences: 2x2 matrix powering IS in TC^0 (Mereghetti-Palano 2000); strong Lucas test uses only fixed-dimension MPOW
- Quadratic Frobenius: equivalent to 2x2 MPOW in Z[x]/(x^2-bx-c, n); also TC^0
- Sum-of-two-squares: works only for p = 1 mod 4; not a complete primality test
- Best path: BPSW = MR(2) + Strong Lucas, both in TC^0 if BPSW has no pseudoprimes

## Verdict
**CLOSED** -- resolved as far as possible. BPSW in TC^0 is established; the open question is whether BPSW is unconditionally correct.
**Failure Mode:** N/A (positive result for TC^0 membership of BPSW components)

## One-Line Summary
Lucas/Frobenius tests are in TC^0 via 2x2 MPOW; Wilson fails (iterated multiplication of n terms); BPSW in TC^0 contingent on correctness.

# Digit-by-Digit Extraction of p(n): Results

**Script:** digit_extraction.py

## What Was Tested
Whether p(n) can be computed digit by digit (most significant to least significant) by computing pi(x) mod 2 or pi(x) mod m efficiently, enabling binary search on x to find p(n).

## Key Findings
- R^{-1}(n) gives top ~50% of digits correctly
- Computing pi(x) mod 2 (parity) would be a breakthrough, but no polylog method is known
- pi(x) mod 2 relates to the Liouville function but computing L(x) = sum lambda(n) has the same complexity as pi(x)
- Character sums for pi(x) mod m reduce to Dirichlet L-function evaluations with the same O(x^{2/3}) cost
- No known method computes any single bit of the "random" part of p(n) in polylog time

## Verdict
**CLOSED** -- Failure Mode: E (Equivalence)

## One-Line Summary
Even computing the parity of pi(x) is as hard as computing pi(x) itself; no digit-by-digit shortcut exists.

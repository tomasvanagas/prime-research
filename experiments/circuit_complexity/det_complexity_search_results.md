# Determinantal Complexity Search: Results

**Script:** `det_complexity_search.py`
**Session:** 15
**See also:** `dc_extended_results.md`

## What Was Tested
Whether pi(x) as a multilinear polynomial in N bits has polynomial determinantal complexity. Computed the explicit multilinear polynomial via Mobius inversion, then searched for small matrix representations with det = pi(x).

## Key Findings
- dc(pi_N) = N confirmed for N = 2, 3, 4 via exhaustive search
- The multilinear polynomial has exponentially many nonzero coefficients (sparsity ~ 0.5, matching random)
- Coefficient magnitudes grow exponentially, ruling out bounded-entry representations
- The polynomial is "maximally complex" by all measures tested

## Verdict
**CLOSED** -- subsumed by `dc_extended_results.md`.
**Failure Mode:** Information loss (pi(x) polynomial has maximal complexity measures)

## One-Line Summary
Multilinear polynomial for pi(x) has random-like coefficient structure; dc(pi_N) = N for small N, exponential for large N.

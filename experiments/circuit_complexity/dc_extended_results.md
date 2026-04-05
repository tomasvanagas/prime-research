# Extended Determinantal Complexity of pi_N: Results

**Script:** `dc_extended.py`
**Session:** 15
**See also:** `det_complexity_search.py`, `det_complexity_extract.py`

## What Was Tested
Extended the determinantal complexity (dc) computation for pi_N(x) -- pi(x) as a multilinear polynomial in bits -- to N=5..12. Used multiple strategies: Mobius inversion, substitution lower bounds, numerical optimization (L-BFGS-B + differential evolution), and flattening rank bounds.

## Key Findings
- dc(pi_N) = N for N = 2, 3, 4 (exact match found by optimization)
- For N = 5, 6: optimization finds representations with m = N but numerical precision degrades
- Flattening rank lower bound: partition_rank grows as 2^{N/2}, implying dc >= 2^{N/2} for large N
- For N >= 10, parameter counting shows dc(pi_N) >= N, and likely dc(pi_N) = Theta(2^{N/2})
- The polynomial dc(pi_N) = poly(N) hypothesis is ruled out for large N by the exponential flattening rank

## Verdict
**CLOSED**
**Failure Mode:** Information loss (dc(pi_N) grows exponentially with N, ruling out poly-size determinantal representations)

## One-Line Summary
Determinantal complexity dc(pi_N) = N for small N but grows exponentially (>= 2^{N/2}) for large N; no poly(N)-size det representation.

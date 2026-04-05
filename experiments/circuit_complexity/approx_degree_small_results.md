# Approximate Degree (Small N): Results

**Script:** `approx_degree_small.py`
**Session:** 23
**See also:** `approx_degree_prime_results.md`

## What Was Tested
Focused LP computation of approximate degree for N=4..10 where the LP is tractable. Computes both feasibility (fixed epsilon) and optimal epsilon for each degree d.

## Key Findings
- For each N, minimum epsilon at degree d decays smoothly as d increases
- Exact degree needed for epsilon < 0.49 (rounding threshold) is ceil(N/2) for even N
- Small-N data matches the scaling fit adeg ~ 0.6 * N^0.9

## Verdict
**CLOSED** -- subsumed by `approx_degree_prime_results.md`.
**Failure Mode:** Information loss (Theta(N) approximate degree)

## One-Line Summary
Small-N LP confirms approximate degree = ceil(N/2), consistent with the Theta(N) scaling.

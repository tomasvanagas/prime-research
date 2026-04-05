# Determinantal Representation Extraction: Results

**Script:** `det_complexity_extract.py`
**Session:** 15
**See also:** `dc_extended_results.md`

## What Was Tested
Harder optimization (differential evolution + L-BFGS-B with more trials) to find and extract explicit determinantal representations det(M) = pi(x) for N=5,6. Analyzes the structure of found matrices.

## Key Findings
- For N=5: found m=5 representations with residual < 1e-8 but matrix entries lack discernible pattern
- For N=6: optimization struggles; best m=6 representation has residual ~0.01 (not exact)
- Extracted matrix structures show no algebraic pattern that would generalize to larger N
- The found representations are "numerical accidents" rather than systematic constructions

## Verdict
**CLOSED** -- subsumed by `dc_extended_results.md`.
**Failure Mode:** Information loss (no systematic construction found; numerical representations do not generalize)

## One-Line Summary
Numerical optimization finds ad-hoc determinantal representations for N <= 5 but no systematic generalizable pattern.

# Parallel Depth Scaling of Lucy DP: Results

**Script:** `parallel_depth_scaling.py`
**Session:** 12
**See also:** `lucy_dp_structure_results.md`, `meissel_lehmer_dag_results.md`

## What Was Tested
Confirmed the scaling of parallel rounds in the Lucy DP (~pi(x^{1/3})), explored better parallelization strategies, and investigated whether the Meissel-Lehmer recursion beats O(x^{1/3}/ln x) depth.

## Key Findings
- Parallel rounds R(x) ~ pi(x^{1/3}) confirmed precisely via greedy scheduling
- No better parallelization strategy found: the dependency structure is inherent to the sieve
- Meissel-Lehmer recursion has depth O(log log x) for the recursive calls but O(x^{1/3}/ln x) for leaf computations
- Cannot beat O(x^{1/3}/ln x) depth with the sieve-based approach
- The bottleneck is the sequential dependency among small primes in the sieve

## Verdict
**CLOSED**
**Failure Mode:** Equivalence (parallel depth is fundamentally O(x^{1/3}/ln x) for sieve methods)

## One-Line Summary
Lucy DP parallel depth confirmed as ~pi(x^{1/3}); no strategy beats O(x^{1/3}/ln x) within sieve-based frameworks.

# Lucy DP Floor-Value Mapping Structure: Results

**Script:** `lucy_dp_structure.py`
**Session:** 12
**See also:** `lucy_dp_matrix_formulation_results.md`, `parallel_depth_scaling_results.md`

## What Was Tested
Algebraic structure of the mappings v -> floor(v/p) within the Lucy DP. Investigated how many distinct mappings occur, graph structure, composition complexity, and parallel-prefix potential.

## Key Findings
- The mapping v -> floor(v/p) is a surjective function from V to V with merging (many-to-one)
- Each prime p produces a distinct mapping; no two primes give the same mapping on V
- The composition of all mappings cannot be computed in O(1) depth because each mapping depends on the floor-division structure specific to that prime
- No parallel-prefix structure: the mappings do not commute (floor(floor(v/p)/q) != floor(floor(v/q)/p) in general on V)
- The sequential dependency chain has depth ~pi(x^{1/3})

## Verdict
**CLOSED**
**Failure Mode:** Equivalence (the mapping structure IS the sieve; non-commutativity prevents parallelization beyond known bounds)

## One-Line Summary
Lucy DP floor-division mappings are non-commutative and produce distinct graph structures per prime; no parallel-prefix shortcut.

# Pfaffian / Planar Graph Matching Approach to Prime Counting -- Results

**Script:** pfaffian_sieve.py

## What Was Tested
Two encodings of the sieve as a graph to test whether bounded treewidth enables efficient computation: (1) sieve bipartite graph (integers x primes, edge if divisible) -- check planarity and treewidth; (2) Lucy DP graph (vertices = distinct floor(x/k) values, edges from recurrence) -- check treewidth/pathwidth.

## Key Findings
- Sieve bipartite graph: NOT planar for any non-trivial x. Treewidth grows as Omega(sqrt(x)/log(x)) -- essentially the number of sieving primes.
- Lucy DP graph: vertices are the O(sqrt(x)) distinct values of floor(x/k). Treewidth measured via min-degree heuristic also grows with sqrt(x).
- FKT algorithm requires planarity (or bounded genus); the sieve graph has unbounded genus.
- Tree decomposition-based FPT algorithms would need treewidth = O(polylog), but measured treewidth is O(sqrt(x)/log(x)).
- The sieve's inherent structure -- each prime interacts with ALL integers it divides -- creates a dense, high-treewidth interaction pattern.
- Previous sessions (10, 20) found volume-law entanglement in the same structure via tensor network analysis, consistent with high treewidth.

## Verdict
**CLOSED**
**Failure Mode:** Information Loss (I) -- sieve graph has treewidth Omega(sqrt(x)/log(x)); no planar/bounded-treewidth structure to exploit.

## One-Line Summary
Pfaffian/FKT sieve: bipartite sieve graph has treewidth Omega(sqrt(x)/log(x)), not planar; no FPT shortcut possible.

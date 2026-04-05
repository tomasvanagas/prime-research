# Tensor Sieve: Results

**Script:** tensor_sieve.py

## What Was Tested
Whether the sieve of Eratosthenes, viewed as a product of per-prime DFA automata (divisibility checkers), admits a compact tensor network / transfer matrix representation: product DFA state complexity, growth rate (polynomial vs exponential), transfer matrix rank and singular value spectrum.

## Key Findings
- Each prime p contributes a DFA with p states; product DFA for primes p1,...,pk has at most p1*p2*...*pk states (primorial)
- Reachable-state BFS: product DFA state count equals primorial for all tested prime sets -- no state reduction from reachability
- Transfer matrix dimension = primorial(y) for sieving primes up to y; for y = 23 (9 primes), matrix is 223092870 x 223092870 -- far too large
- Singular value spectrum of small transfer matrices: full rank, no decay -- cannot be approximated by low-rank matrices
- Counting survivors (primes) via matrix trace Tr(T^N) would require exponential-size matrix exponentiation
- The tensor product structure offers no compression because the individual DFA states are entangled through the coprimality constraint

## Verdict
**CLOSED**
**Failure Mode:** Equivalence (E) -- tensor/DFA product has state count = primorial (exponential); transfer matrices are full rank with no low-rank structure.

## One-Line Summary
Tensor sieve (DFA product + transfer matrix): state count = primorial, full-rank matrices, no compression -- exponential, not polylog.

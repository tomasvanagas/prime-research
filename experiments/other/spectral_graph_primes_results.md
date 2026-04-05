# Spectral Graph Theory Approaches to pi(x): Results

**Script:** spectral_graph_primes.py

## What Was Tested
Three spectral graph approaches: Cayley graph of Z/xZ with prime generators (eigenvalues are character sums over primes), Ihara zeta function of related graphs, and expander mixing lemma. Also tested prime-free graph constructions: GCD graph and divisor graph.

## Key Findings
- Cayley graph Cay(Z/xZ, Primes(x)): the degree lambda_0 = pi(x), which is trivially circular -- you need primes to build the graph
- Eigenvalues lambda_k = sum exp(2*pi*i*k*p/x) are Vinogradov-type exponential sums; computing them is as hard as the explicit formula
- Ihara zeta function of the Cayley graph encodes the spectrum, which encodes pi(x) -- no shortcut
- GCD graph (edge iff gcd(i,j)=1): constructible without primes, but its spectrum encodes Euler's totient, not pi(x) directly
- Divisor graph: spectrum relates to divisor sums but not to prime counting
- Expander mixing lemma gives approximate counts only, with error proportional to spectral gap

## Verdict
**CLOSED** -- Failure Mode: C (Circularity) / E (Equivalence)

## One-Line Summary
Cayley graph construction requires knowing primes; prime-free graphs (GCD, divisor) don't encode pi(x) in their spectrum.

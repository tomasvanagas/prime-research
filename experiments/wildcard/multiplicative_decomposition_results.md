# Multiplicative Decomposition: Results

**Script:** multiplicative_decomposition.py

## What Was Tested
Decomposing pi(x) via Dirichlet characters: computing pi(x; q, a) for residue classes, testing whether character decomposition reveals cancellation structure, and whether differences pi(x; q, a1) - pi(x; q, a2) (which cancel the main term) are cheaper to compute.

## Key Findings
- Character decomposition: pi(x; q, a) = (1/phi(q)) * sum_chi chi_bar(a) * Psi_chi(x) where Psi_chi involves zeros of L(s, chi)
- Differences pi(x; q, a1) - pi(x; q, a2) cancel the main Li(x)/phi(q) term -- but the oscillatory parts involve L-function zeros which are ADDITIONAL to zeta zeros
- L-function zeros for non-trivial characters: these are different zeros, not simpler; computing them costs at least as much
- Reconstructing pi(x) from character sums: need ALL characters mod q, totaling phi(q) L-functions
- Total cost: O(phi(q) * x^{2/3}) per modulus -- strictly MORE expensive than direct pi(x) computation

## Verdict
**CLOSED**
**Failure Mode:** Equivalence (E) -- Dirichlet character decomposition introduces L-function zeros (additional to zeta zeros), increasing rather than decreasing computational cost.

## One-Line Summary
Multiplicative decomposition via Dirichlet characters: L-function zeros are additional cost; phi(q) L-functions needed per modulus -- worse than direct pi(x).

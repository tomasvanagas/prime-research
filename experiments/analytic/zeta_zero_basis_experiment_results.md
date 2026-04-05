# Zeta Zero Basis Functions for pi(x): Results

**Script:** zeta_zero_basis_experiment.py

## What Was Tested
Using R(x^rho) terms from the explicit formula as basis functions instead of li(x^{1/k}): (1) how error scales with K zeros and x, (2) whether REORDERING zeros (by contribution size rather than |gamma|) converges faster, (3) whether O(polylog(x)) terms suffice for error < 0.5.

## Key Findings
- Standard ordering (by |gamma|): error scales as O(sqrt(x) * ln(K) / K)
- Reordering by |cos(gamma*ln(x))|/gamma (largest contribution first): ~3x faster convergence at fixed K
- BUT: optimal reordering is x-DEPENDENT -- different x needs different ordering
- No universal reordering achieves polylog convergence for all x
- For error < 0.5: need K ~ O(sqrt(x)/ln(x)) zeros regardless of ordering
- The zeta zero basis IS the natural basis for pi(x); no other basis outperforms it
- The bottleneck is the NUMBER of basis functions, not their ordering
- Confirms: O(sqrt(x)) zeros are fundamentally required; no polylog subset suffices

## Verdict
**CLOSED**
**Failure Mode:** E (Equivalence -- zeta zero basis functions confirm O(sqrt(x)) requirement; reordering gives constant-factor only)

## One-Line Summary
Zeta zero basis experiment: reordering zeros by contribution gives ~3x improvement; still need O(sqrt(x)/ln(x)) zeros for exactness.

# Session 36 Proposals — Fresh Approaches to Computing p(n) Exactly

**Date:** 2026-04-05
**Method:** Web search for 2025-2026 papers + 4 concrete proposals with computational tests

---

## Summary

Four fresh proposals were developed, coded, and tested on n < 10,000. All four illuminate WHY O(polylog) is hard, but none achieves it. Key findings:

| # | Proposal | Core Idea | Result | New Insight? |
|---|---------|-----------|--------|-------------|
| 21 | Zero Clustering Truncation | Group zeta zeros by proximity, use GUE statistics | Tail contribution unpredictable; O(sqrt(x)) zeros still needed | Confirms zero independence |
| 22 | Compressed Sensing on delta(n) | Sparse recovery of residual in zeta zero basis | delta is NOT sparse: 50+ coefficients for 87% energy | Quantifies non-sparsity |
| 23 | PSLQ Integer Relations | Search for algebraic formula for delta(n) | Normalized delta ≈ Gaussian; PSLQ finds spurious point-wise relations only | Normality test passed |
| 24 | Dequantized Grover Sieve | Low-rank decomposition of primality operator | Sieve rank ~ N^{0.365}, Fourier rank ~ N^{0.94} | Quantifies "hardness" of prime indicator |

---

## Proposal 21: Zero Clustering Truncation

### Mathematical Idea
The Riemann explicit formula requires summing over all nontrivial zeta zeros:
π(x) = R(x) − Σ_ρ R(x^ρ)

If zeros exhibit correlated contributions due to GUE spacing statistics, we can group them into clusters and represent each cluster by a weighted representative. If cluster count grows as O(log(x)^C), we achieve polylog.

### Key Conjecture
The oscillatory sum S(x) = Σ_γ x^{iγ}/ρ converges within O(1) using O(log(x)^C) representative zeros.

### Computational Test Results
- 10 zeros already reduce mean |pi(x) - estimate| from 14.0 to 1.3
- But adding more zeros (50, 200) does NOT improve — errors plateau at ~1.4
- Clustering with radius 10 gives 35 clusters with occasionally excellent results
- **Critical failure**: tail contributions (zeros beyond cutoff) are unpredictable
  - tail/partial ratio varies from -3.4 to +3.3 — no pattern

### Verdict: CLOSED
Each zeta zero carries independent information. Clustering helps constants but cannot reduce asymptotics below O(sqrt(x)).

---

## Proposal 22: Compressed Sensing on delta(n)

### Mathematical Idea
If delta(n) = p(n) − R^{-1}(n) is k-sparse in some transform domain, compressed sensing recovers it from O(k log(n/k)) measurements. Test sparsity in: (a) zeta zero basis, (b) standard Fourier basis.

### Key Conjecture
delta(n) is O(log(n))-sparse in the zeta zero basis.

### Computational Test Results
- **Sparsity test failed**: 100 basis functions capture only 66% of variance; 69/100 coefficients are significant
- **Fourier sparsity**: Top 50 coefficients capture 87% energy — moderate, not sparse
- **Extrapolation completely fails**: training on n=10..5000, testing on 5001..10000 gives only 11-16% correct prime identification
- L1 minimization (LASSO) does not improve over least squares

### Verdict: CLOSED
delta(n) is information-theoretically incompressible in all tested bases. The oscillatory contribution from ~10^48 zeta zeros cannot be sparsified.

---

## Proposal 23: PSLQ/LLL Integer Relation Discovery

### Mathematical Idea
Use PSLQ/LLL lattice reduction to search for integer relations between delta(n) and a library of known constants (π, γ, log(2), zeta values, etc.). Inspired by the Ramanujan Library approach (arXiv:2412.12361).

### Key Conjecture
There exists a polynomial P of bounded degree such that delta(n) = P(log n, log log n, ...) + O(1).

### Computational Test Results
- **Polynomial regression R² ≈ 0.001**: delta(n) is uncorrelated with ALL smooth functions of n
- PSLQ finds integer relations at individual n values, but **different formula at each n** — overfitting
- Normalized delta δ̂(n) = δ(n)·log(p)/√p passes Gaussian normality test (p=0.19)
- Autocorrelation: 0.96 at lag 1, decaying to 0.39 at lag 50

### Verdict: CLOSED
The near-Gaussian distribution of normalized delta confirms it encodes random information. No universal algebraic formula exists in the tested constant space.

---

## Proposal 24: Dequantized Grover Sieve

### Mathematical Idea
If the "primality operator" M (matrix encoding prime/composite status) has low rank, dequantization of quantum algorithms gives classical O(polylog). Measure the rank of M in multiple decompositions.

### Key Conjecture
The primality operator has effective rank O(x^{1/3}/log(x)), matching Meissel-Lehmer.

### Computational Test Results
- **Sieve matrix rank ~ N^{0.365}** ≈ π(√N), matching theoretical prediction
- **Fourier 90% rank ~ N^{0.943}** — nearly full rank in Fourier domain
- SVD spectrum decays slowly (no sharp cutoff) — NOT low-rank
- Dirichlet character decomposition shows uniform distribution in residue classes

### Verdict: CLOSED
The primality indicator is high-rank in all tested decompositions. Dequantization gives at best O(N^{1/3}), matching Meissel-Lehmer, not polylog.

---

## Cross-Cutting Insights

1. **10 zeros suffice for O(1) accuracy** in the explicit formula, but this O(1) error is still ~1.3 on average — not enough to identify the exact prime. The "last mile" from approximate to exact requires O(sqrt(x)) zeros.

2. **delta(n) is approximately Gaussian** when normalized by sqrt(p)/log(p). This is the strongest evidence yet that the residual is genuinely random-like, not just "hard to compute."

3. **The primality operator has polynomial rank** (N^{0.365} in sieve form). Any approach based on low-rank or sparse structure is fundamentally limited to sublinear but not polylogarithmic complexity.

4. **The information bottleneck is the oscillatory part** of the explicit formula, which encodes contributions from ~T/log(T) zeta zeros up to height T ~ sqrt(x). Each zero contributes independently. There is no known shortcut to sum them.

---

## Literature Searched
- Recent dequantization: Tang (2019), Chia et al. (2020), Jethwani et al. (2025)
- Ramanujan Library / automated conjecture: arXiv:2412.12361 (2024)
- Murmurations and trace formulas: arXiv:2506.01640 (2025)
- Kedlaya point counting via Monsky-Washnitzer cohomology
- Étale cohomology computation: single exponential time (JTNB 2020)
- No 2025-2026 breakthroughs found for exact prime counting

## Files Created
- `experiments/proposals/proposal21_zero_clustering_truncation.py` + `_results.md`
- `experiments/proposals/proposal22_compressed_sensing_delta.py` + `_results.md`
- `experiments/proposals/proposal23_pslq_delta_relations.py` + `_results.md`
- `experiments/proposals/proposal24_dequantized_grover_sieve.py` + `_results.md`

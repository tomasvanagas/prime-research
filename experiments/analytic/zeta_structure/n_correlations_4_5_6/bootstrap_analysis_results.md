# Bootstrap analysis — matched-finite-N GUE comparator for §C2

Follow-up to `n_correlations_4_5_6.py`. The primary script's gap-
shuffled null gives spuriously large z-scores in the P_k spacings
because gap-shuffling destroys GUE rigidity (a known GUE property,
not arithmetic content). This script swaps in a proper **matched-
finite-N GUE pool** (K=20 independent GUE batches × N=2000 Hermitian
Wigner matrices, central 60% kept ~ 1200 unfolded eigenvalues each)
and additionally a **matched-size empirical-zeta sub-chunking** (20
chunks × 400 zeros), so that:

- `z_full_vs_GUE_batch` compares empirical zeta on full 8000 zeros
  against the GUE batch distribution (which has matched-finite-N
  noise on a smaller sample size — useful for spotting bias
  inconsistencies).
- `z_full_vs_theory` compares empirical zeta on full 8000 zeros
  against the **theoretical GUE prediction** `R_n = det[K((j-i)s)]`,
  using `zeta_chunk_std / √20` as the bootstrap SE — the right
  noise estimate for the empirical full-N estimator.

## Headline conclusion

| Test | Probe | max |z_vs_theory| | At |
|------|-------|-------------------|-----|
| R_4 equally-spaced | n=4, s ∈ {0.5..5} | 2.36σ | s=2.0 |
| R_5 equally-spaced | n=5, s ∈ {1..5} | 6.00σ raw | s=2.0 (reproduced by GUE batch — Poisson shot noise) |
| R_6 equally-spaced | n=6, s ∈ {1..5} | 1.56σ | s=2.0 |
| κ_n vs GUE batches at L ≥ 8 | n ∈ {3..6} | < 2.1σ | all bins |

The full result table is in `bootstrap_analysis_results.json`.

## Falsifier (pre-registered)

A surviving |z_vs_theory| > 5σ at any (n, s) bin where:
- the apparent deviation is NOT reproduced by matched-finite-N GUE
  Monte Carlo (rules out estimator bias);
- the empirical zeta sub-chunk std does NOT predict the deviation
  (rules out shot noise);
- the deviation has a structural explanation tied to a specific
  prime arithmetic mechanism.

**Result:** no surviving deviation > 3σ across any (n, s) bin or
(n, L) cumulant cell at L ≥ 8.

## Cross-reference

Main experiment + interpretation: `n_correlations_4_5_6_results.md`.

# G1 Analysis Helper Companion

`liouville_anderson_lyapunov_analyze.py` is the diagnostic helper for
the main `liouville_anderson_lyapunov.py` experiment. It reads the
`results_N*_s*_e*.json` files produced by the main script and
computes:

1. Per-energy z-scores `(γ_λ(E) - γ_Rad_mean(E)) / γ_Rad_std(E)`
   and the χ² aggregate `Σ z²` (Rademacher prediction = K).
2. L²-distance from the Rademacher mean for the λ curve and for
   each Rademacher seed; reports the rank of λ in the seed
   distribution (empirical p-value).
3. Pastur-Figotin agreement: `γ_λ / γ_PF` and `γ_Rad_mean / γ_PF`
   inside the band [-2, 2].
4. Independent two-point Chowla auto-correlation of the bare λ
   sequence at lags `h ∈ [0, hmax]`.

All findings are summarised in the main `_results.md` (TL;DR table,
per-energy peaks, Pastur-Figotin and Chowla z-scores). This companion
contains no novel content — it is a measurement instrument, not an
experiment in its own right.

## Reproduce

```
python3 liouville_anderson_lyapunov_analyze.py \
  --paths results_N100000_s50_e51.json \
          results_N300000_s50_e51.json \
          results_N1000000_s100_e51.json \
  --chowla-N 1000000 --chowla-hmax 16
```

Output is `analysis_summary.json` plus a printed report to stdout.

See `liouville_anderson_lyapunov_results.md` for the full structural
narrative, edge citations, and grading.

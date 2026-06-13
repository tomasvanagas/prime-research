# d41_f2_matched.py — Phase 3 results

Single experiment, four scripts, ONE write-up. See
[`d41_chip_3tensor_results.md`](d41_chip_3tensor_results.md) §"Phase 3"
for full discussion.

Phase 3 = F²-MATCHED Bernoulli null (calibrated q* such that
`E[F²_bern] = F²_chi_P`). The A-vs-B critical test: removes the
trivial HL-S₃(N) ternary-count scaling, isolates residual deviation.

Headline: σ_max excess COLLAPSES to 0.7-7.4σ after F²-matching;
ρ_top excess persists at 4.8-33.8σ across N ∈ {1023, 2047, 4095, 8191,
16383}. Residual correlates with N's smoothness — peak at Mersenne
prime N=8191.

Raw output: `results_phase3.json`.

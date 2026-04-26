# Gowers U^k Analysis Module — Results

This is a post-processing utility that consumes the JSON output of
`gowers_uk_chi_p.py` and combines it with the Hardy-Littlewood
singular series prediction from `hardy_littlewood_box.py`. Computes:

  `Q^k(f) := ||f||_{U^k}^{2^k} / E[f]^{2^k}`

ratios for chi_P, Bernoulli baselines, Liouville, and W-tricked
chi_{W,b}, alongside the HL prediction `S_k`.

The output table is written to
`gowers_uk_chi_p_data/analysis.md` for each invocation.

For full empirical results and falsification statement, see
`gowers_uk_chi_p_results.md`.

# `analyze_fhk.py` — post-hoc analysis of FHK results

This is a helper script, not a separate experiment. It reads
`fhk_amplitude_max_results.json` produced by `fhk_amplitude_max.py`
and emits additional diagnostics:

1. Per-T bootstrap CIs on M_T mean and variance.
2. Free-Gumbel MLE fit + KS distance comparison vs FHK Gumbel(1/2)
   and free Gaussian.
3. Vuong test (Gauss vs Gumbel) per T anchor.
4. Three-T linear regression `<max> = a + b · log log T + c · log log log T`.
5. Pooled M_T statistics (skew, ex-kurtosis, KS) across all 300 windows.

Output: `fhk_amplitude_max_analysis.log`.

Findings, mechanism, falsification verdict, edge / closure entries:
see `fhk_amplitude_max_results.md` (the canonical experiment results
file).

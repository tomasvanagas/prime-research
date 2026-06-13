# `smoothed_kernels_paired_analysis.py` — post-processor results

This script is a paired-sign-test post-processor on the per-sample
CSV output of `smoothed_kernels.py`. It produces no standalone
results — its output (paired ratios, win counts, sign-test p-values,
per-kernel aggregates) is reflected directly in:

- `smoothed_kernels_paired.log` — the full per-cell grid.
- `smoothed_kernels_results.md` §"Headline grid", §"Per-kernel
  aggregate across 12 cells", §"Why the prediction holds (and why
  some sub-1 ratios appear)".

For the substantive write-up of slot 3 / Thread 7, see
`smoothed_kernels_results.md`.

## What this post-processor adds beyond `smoothed_kernels.py`

`smoothed_kernels.py` reports σ_eff(kernel)/σ_eff(hard) at matched
(anchor, K_compute), which is the L2-ratio of error magnitudes
across the N=20 samples. This post-processor adds the **paired sign
test**:

  wins/N = #{i ∈ 1..N : |err_kernel(x_i)| < |err_hard(x_i)|}
  p_kernel_beats = P(Binom(N, 0.5) ≥ wins) (one-sided)
  p_hard_beats = P(Binom(N, 0.5) ≥ N − wins) (one-sided)

The paired test is more sensitive than the L2 ratio because both
err_kernel(x_i) and err_hard(x_i) are functions of the same sample x_i
(they share the same per-zero contributions c_k = 2 Re R(x_i^{ρ_k})),
so the per-sample noise cancels in the difference.

Result: 0 of 96 cells show p_kernel_beats < 0.05; up to 3 cells per
kernel show p_hard_beats < 0.05. Asymmetric pattern confirms the
L2-optimality prediction holds under GUE corrections.

## Falsifiability

This post-processor's claim (i.e., the asymmetric pattern) is
falsified by:

- A re-run with N ≥ 100 paired samples that finds at least one
  (anchor, K, kernel) cell with p_kernel_beats < 0.01 — i.e., a
  decisively kernel-beats-hard cell.
- Independent replication on a different x-anchor set (e.g., 10⁶,
  10¹¹) that contradicts the asymmetric pattern.

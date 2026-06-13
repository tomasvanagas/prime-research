# d41_eigvec_hl.py — Phase 4 results

Single experiment, four scripts, ONE write-up. See
[`d41_chip_3tensor_results.md`](d41_chip_3tensor_results.md) §"Phase 4"
for full discussion.

Phase 4 = leading eigenvector vs HL R₂(N-p) — DEFINITIVE TEST.

Headline: at N ∈ {2047, 8191}, Spearman(empirical degree d_p,
HL Goldbach R₂(N-p)) = 0.993-0.997; Cosine(v₁, sqrt(HL R₂(N-p))) =
0.988. The chi_P 3-tensor matricisation has rank-1 Perron-Frobenius
leading eigenvector v_p ∝ sqrt(R₂(N-p)) — full HL-class identification.
Melonic universality REJECTED.

At N=16383 (smooth: 3·43·127), |E₂|/|E₁| = 0.97 — near-degenerate
top eigenvalues; v₁ becomes a mixture, dropping cosine to 0.81. But
emp_d still tracks HL_d at 0.996 — degree distribution remains HL-class.

Raw output: `results_phase4.json`.

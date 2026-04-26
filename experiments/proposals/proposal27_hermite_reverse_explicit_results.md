# Proposal 27 Results — Hermite-Mollified Reverse Explicit Formula

## Question
Does Gaussian or Hermite-style mollification of the truncated explicit formula
ψ(x) ≈ x − Σ_{|γ|<T} x^ρ/ρ converge with fewer zeros K than the unmollified
truncation?

## Setup
- Loaded 1000 Riemann zeros from `data/zeta_zeros_1000.txt`.
- Tested estimators at x ∈ {100, 1000, 10000}:
  (a) unmollified truncated explicit formula,
  (b) Gaussian mollifier exp{−(γ − log x)² / (2 σ²)} for σ ∈ {0.5, 1, 2},
  (c) Riesz mean (1 − |γ|/T)² for T ∈ {50, 200, 800}.
- Each estimator scored at K ∈ {10, 25, 50, 100, 200, 400, 800}.
- Compared against ψ(x) computed via direct sieve.

## Result

| x      | unmoll @ K=800 | best Gauss @ K=800 | best Riesz @ K=800 | ratio (unmoll / best moll) |
|--------|----------------|---------------------|---------------------|------------------------------|
| 100    | 0.076          | 4.117               | 0.428 (T=800)       | **0.18** (mollified worse)   |
| 1000   | 0.492          | 1.481               | 0.919 (T=200)       | **0.54** (mollified worse)   |
| 10000  | 4.483          | 14.560              | 3.691 (T=50)        | 1.21 (mollified ~1.2× better)|

The Gaussian mollifier does worse at all tested x: it discards information from
zeros far from log x, but those zeros contribute non-trivially to the
square-root cancellation. The Riesz mean is closer to the unmollified estimator
and only marginally better at the largest x (1.21×) — a constant factor at most.

## Verdict — CLOSED (failure mode: equivalence)

The Hermite/Gaussian mollification does not bypass the square-root cancellation
barrier; it merely re-weights the zeros. The truncation error is asymptotically
the same Ω(√x). The 2025 Hermite-weighted explicit formula (arXiv 2312.00108)
is genuinely useful for *computing zeros from primes* (the LHS smearing
concentrates contribution near a target height), but the inverse direction does
not gain — the key "primes ↔ zeros" identity is symmetric under truncation in
that the lossy term is the truncation tail, not the kernel weight.

## Failure mode
Equivalence (E): mollification is a different parameterisation of the same
truncated zero sum; truncation tail bound is unchanged.

## Key takeaway
Any *linear* functional of the zero list whose kernel decays slower than 1/|γ|
inherits the Ω(√x) barrier. To beat √x we would need a *non-linear* operation
on zeros (e.g., spectral product), not a re-weighting.

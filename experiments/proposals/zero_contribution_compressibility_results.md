# Compressibility of the Zeta-Zero Contribution Matrix — Results

## Experiment
Build M[i, j] ≈ 2 Re(li(x_i^{rho_j})) using the asymptotic
li(x^rho) ≈ x^rho / (rho * log x) for x_i in two windows
(x in [1000, 2000) and x in [5000, 6000)) and rho_j = 1/2 + i*t_j for
the first K Riemann zeros. Test SVD decay and 2D Fourier sparsity.

## Numerical results

### Window 1: x ∈ [1000, 2000), K = 1000 zeros (1000×1000 matrix)

| k    | sigma_k / sigma_0 |
|------|-------------------|
|   1  | 1.000             |
|   2  | 0.750             |
|   4  | 0.457             |
|   8  | 0.341             |
|  16  | 0.220             |
|  32  | 0.129             |
|  64  | 0.073             |
| 128  | 0.041             |
| 256  | 0.021             |
| 512  | 5e-15 (numerical zero)|

- Effective rank @ 1% threshold: **312** (out of 1000)
- Effective rank @ 0.1% threshold: 318
- 2D Fourier coefficients > 1% magnitude: 29.75% of cells
- 2D Fourier coefficients > 0.1%: 78.27% of cells

### Window 2: x ∈ [5000, 6000), K = 2000 zeros

- Effective rank @ 1%: **147**
- 2D Fourier > 1%: 9.68% of cells

### Baselines

- Random Gaussian (1000×1000): effective rank 987 — essentially full.
- Rank-10 + noise: effective rank 10 — clean.

## Power-law decay rate

Fitting sigma_k/sigma_0 ∝ k^{-alpha}:

| k    | alpha (instantaneous) |
|------|-----------------------|
|   2  | 0.415                 |
|   8  | 0.519                 |
|  32  | 0.589                 |
| 128  | 0.658                 |
| 256  | 0.703                 |

So sigma_k decays roughly as k^{-0.5 to -0.7}, classic algebraic
decay, NOT exponential. To reach sigma_k < 1/(2 * F-norm) — needed to
integer-recover pi(x) — we need K = polynomially many singular
components, which is sub-linear but not polylog.

## Verdict: CLOSED (failure mode I — information loss)

The compressibility conjecture is rejected at this threshold.

- **Polynomial compressibility, not polylog.** Effective rank scales
  sub-linearly with the matrix dimension but the decay law is
  algebraic (k^{-0.7}-ish), so reaching integer precision needs
  K = poly(window) zero components. Polylog is unreachable.
- **Fourier basis doesn't help.** ~30% of 2D Fourier coefficients are
  above 1% relative magnitude, no sparse representation.
- **Pattern matches GUE prediction.** Algebraic decay with exponent
  ≤ 1 is the signature of a "1/f-like" spectrum, consistent with the
  random matrix theory description of zero contributions (and with
  the project's pseudorandomness theme).

## What this rules in

The matrix is decisively **not** rank-O(polylog) in any tested basis,
so direct compressed-sensing recovery from polylog measurements is
impossible. This is one more confirmation that the explicit-formula
sum has no hidden low-dimensional structure beyond the trivial Riemann
correction.

## What this does NOT rule out

- A different *transform* (nonlinear; e.g., learnable encoder, kernel
  embedding) might yield more compression.
- Window-specific local low-rank: maybe in a logarithmically-scaled
  window x ∈ [N, N(1 + 1/log N)) the rank is much lower. Worth a
  follow-up.
- Gourdon's empirical observations about which zeros matter for which
  x are not directly contradicted — those exploit *importance
  sampling* among zeros, which is a different conjecture (Proposal 4).

## Files

- `zero_contribution_compressibility.py` — runnable experiment
- This file — results

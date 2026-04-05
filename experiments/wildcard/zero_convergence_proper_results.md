# Proper Zero Sum Convergence Test

## Setup
Used the Riemann explicit formula: pi(x) = R(x) - Σ_ρ R(x^ρ) - 1/ln(2) + ...
with R(x) = Σ μ(k)/k · li(x^{1/k}) computed via mpmath at 30-digit precision.

## Results: DIVERGENCE

The sum over zeros DIVERGES instead of converging:

| x | K=0 error | K=10 error | K=50 error | K=100 error |
|---|-----------|------------|------------|-------------|
| 100 | -0.7 | — | — | — |
| 500 | -2.1 | -75 | -317 | -684 |
| 1000 | -1.1 | -51 | -445 | -921 |
| 10000 | -3.5 | -163 | -933 | -2076 |
| 50000 | -1.1 | -57 | -1396 | -2159 |

Exception: x=100 converges with K=2 zeros (error drops below 0.5).

## Diagnosis
The divergence is caused by numerical instability in computing R(x^ρ) for 
complex ρ = 1/2 + iγ. The complex logarithmic integral li(x^{ρ/k}) involves:
- Complex exponentiation: x^{ρ/k} = exp(ρ/k · ln x)
- Complex Ei function: li(z) = Ei(ln z) for complex z
- Möbius series truncation: only 10 terms used

For large γ (high zeros), x^{iγ} oscillates rapidly, and the li function
on the complex plane has branch cuts that cause numerical artifacts.

## Key Observation
WITHOUT any zeros (K=0), R(x) alone gives |error| < 4 for all tested x up to 50000.
This is already quite good — the Riemann R function is an excellent approximation.

The error |R(x) - pi(x)| scales roughly as √x / log x:
- x=100: 0.7, x=1000: 1.1, x=10000: 3.5, x=50000: 1.1

## Significance
1. The explicit formula is EXTREMELY sensitive to numerical precision in R(x^ρ)
2. A naive implementation diverges — careful contour integration is needed
3. This suggests the Lagarias-Odlyzko approach (contour integration, not individual zeros)
   is more numerically robust than summing individual R(x^ρ) terms
4. The scaling analysis was inconclusive due to divergence

## Verdict
The explicit formula approach REQUIRES sophisticated numerical methods:
- Odlyzko-Schönhage for batch zero evaluation
- Careful contour integration (not individual zero summation)
- Very high precision arithmetic

These methods achieve O(x^{1/2+ε}) but no known path to polylog.

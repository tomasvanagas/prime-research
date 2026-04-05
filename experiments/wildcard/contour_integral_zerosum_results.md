# Contour Integral & Spectral Methods for Zero Sum Evaluation

## Experiments Run

### Experiment 1: Zero Sum Convergence (BUGGED)
- The li(x^ρ) computation via Ei series overflows numerically for large |ρ·log(x)|
- The simpler x^ρ/ρ approximation works fine (see Experiment 6)
- **Verdict**: Implementation issue, not conceptual. The direct x^ρ/ρ sum converges conditionally.

### Experiment 2: Contour Integral Approach
- **Key insight**: The contour integral ∮ (ζ'/ζ)(s)·x^s/s ds is EQUIVALENT to summing over zeros by the residue theorem
- The number of quadrature points must be ≥ the number of enclosed zeros (Nyquist bound)
- Evaluating ζ'/ζ at each quadrature point costs O(T^{1/2}) via Riemann-Siegel
- **Verdict**: DEAD END. No shortcut over direct zero summation.

### Experiment 3: Smoothed Explicit Formula
- Gaussian smoothing with width h: weights zeros by exp(-h²γ²/2)
- For even h=0.5, the Gaussian kills ALL zeros (γ₁ ≈ 14.13 is too large)
- Smoothing trades zeros for bias — fundamentally cannot give exact answers with few zeros
- **Verdict**: DEAD END. Uncertainty principle: frequency resolution × spatial resolution ≥ 1.

### Experiment 4: Compressibility of Zero Sum (MOST INTERESTING)
- SVD of the zero-contribution matrix M[x_i, ρ_k] = 2·Re(x_i^{ρ_k}/ρ_k)
- 500 x-values × 500 zeros
- Singular value decay:
  - 45 components capture 90% of energy
  - 133 components for 99%
  - 223 for 99.9%
  - NOT rapidly compressible
- Low-rank reconstruction: errors of 1-70 depending on x and rank
- **Verdict**: The zero sum has slow singular value decay. No evidence of rapid compressibility. The SVD structure depends on the x-range, not intrinsic to the zeros.

### Experiment 5: Modular Prime Counting
- Primes equidistribute among residue classes (Dirichlet's theorem)
- No useful periodicity or autocorrelation in p(n) mod m sequences
- The parity problem (Selberg) blocks computing π(x) mod 2 via sieves
- **Verdict**: DEAD END. The parity barrier is fundamental.

### Experiment 6: L-function Combinations
- Tested direct x^ρ/ρ sum with 5 to 200 zeros
- For x=50000: 100 zeros gives |error|=4.67 (decent, but not polylog)
- Error doesn't decrease monotonically (conditional convergence)
- Different L-function combinations give same convergence rate (same zeros)
- **Verdict**: DEAD END. All standard L-functions share the same zeros.

## Key Theoretical Insights

1. **Contour integral = zero sum**: No way around enumerating zeros via integration.
2. **Smoothing = information loss**: The uncertainty principle prevents fewer zeros + exact answer.
3. **SVD decay is slow**: The zero contribution matrix is effectively full-rank, confirming information-theoretic barrier.
4. **Parity blocks CRT**: Can't even compute π(x) mod 2 via sieves.

## Overall Verdict
All six approaches confirm the known barrier: the oscillatory correction to π(x) requires O(x^{1/2}/log x) pieces of information (whether zeros, quadrature points, or SVD components). No spectral shortcut found.

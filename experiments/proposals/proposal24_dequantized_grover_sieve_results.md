# Proposal 24: Dequantized Grover Sieve — Results

## Idea
If the "primality operator" (matrix encoding which numbers are prime) has low rank, dequantization techniques from quantum ML could give a classical O(polylog) algorithm. Tested by analyzing the rank of various primality-related matrices.

## Key Results

### Rank scaling
| Decomposition | Scaling exponent α | Interpretation |
|--------------|-------------------|---------------|
| Sieve matrix rank | N^0.365 | Close to pi(sqrt(N)) ~ N^{1/2}/log(N) |
| Fourier 90% energy | N^0.943 | Nearly full rank |
| Fourier 99% energy | N^0.984 | Nearly full rank |

The sieve matrix has rank exactly equal to pi(sqrt(N)) — this is the number of primes needed for sieving, which grows as N^{1/2}/log(N). This is **sublinear but not polylog**.

### SVD spectrum
The SVD of the prime indicator reshaped as a matrix shows slow decay:
- N=1000: top 10 singular values normalized: 1.00, 0.59, 0.51, 0.45, 0.43, 0.38, 0.35, 0.32, 0.31, 0.27
- No sharp cutoff → NOT low-rank

### Quantum-inspired sampling
Random sampling with O(sqrt(N)) samples gives 5-23% error in pi(N) — consistent with Monte Carlo bounds, no quantum advantage.

### Dirichlet character decomposition
Primes distribute nearly uniformly among residue classes (Bombieri–Vinogradov). Max deviation from expected is small (≤5 primes for q≤12, N=1000).

## Verdict: CLOSED
The primality operator has rank O(N^{1/3} to N^{1/2}), consistent with known Meissel-Lehmer bounds. Dequantization cannot improve beyond O(N^{1/3}). Polylog rank would require a breakthrough equivalent to proving primes have hidden low-dimensional structure — which contradicts the GUE/random-like behavior observed in all other experiments.

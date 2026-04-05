# Session 40: Fresh Perspective (2026-04-05)

## Theme
First-principles attack on O(polylog) prime computation, deliberately avoiding reading CLOSED_PATHS.md to force fresh thinking.

## Approach
Brainstormed 5 genuinely unconventional ideas inspired by barrier-breaking analogies:
- Shor's algorithm (quantum Fourier transform)
- Compressed sensing (sparsity assumptions)
- Fast multipole method (hierarchical decomposition)
- Matrix completion (nuclear norm relaxation)

## Experiments

### 1. Transfer Matrix Sieve (`experiments/wildcard/transfer_matrix_sieve.py`)
**Idea:** Encode sieve as matrix product, iterate via repeated squaring.
**Result:** State space = lcm(primes up to sqrt(x)) ~ exp(sqrt(x)). Even with blocking, work ~ 10^100 vs ML's 10^67 for x=10^100. Depth-k sieve reinvents Meissel-Lehmer.
**Verdict:** FAIL (Equivalence)

### 2. Wavelet/Fourier Compression (`experiments/wildcard/wavelet_zero_compression.py`, `wavelet_detrended.py`)
**Idea:** If C(x) = pi(x) - li(x) is sparse in wavelets, evaluate at point using O(polylog) coefficients.
**Result:** Raw signal: 99% energy in 1 coeff (DC artifact). Detrended: N^{0.75} coefficients for 99.9% energy. Power spectrum f^{-1.7} across all scales.
**Verdict:** FAIL (Information). Sublinear compression but far from polylog.

### 3. Modular Residue Prediction (`experiments/wildcard/modular_residue_prediction.py`)
**Idea:** Predict p(n) mod q from n for all small q, then CRT to recover p(n).
**Result:** Near-max entropy for q>2. CRT is circular. Primorial mod-210 shows 19x prediction from wheel structure, but not useful.
**Verdict:** FAIL (Circularity + Information)

### 4. Carry-Propagation Boundary (`experiments/wildcard/carry_propagation_boundary.py`)
**Idea:** Map which bits of p(n) are easy (from R^{-1}) vs hard (need zeros).
**Result:** Sharp sigmoid at ~60% of bit positions. Transition width ~4 bits. Error ~ p(n)^{0.66}. No intermediate-difficulty bits.
**Verdict:** FAIL as approach, but **novel quantitative measurement**.

### 5. Ergodic Fast-Forward (`experiments/wildcard/ergodic_fast_forward.py`)
**Idea:** Model prime gaps as dynamical system, test fast-forwardability.
**Result:** Lyapunov lambda = 1.06 (chaotic). R^2 < 0.006. Real gaps MORE chaotic than Cramér model.
**Verdict:** FAIL (Information)

### 6. Multiplicative Convolution Shortcut (`experiments/wildcard/multiplicative_convolution_shortcut.py`)
**Idea:** Factor pi(x) as Dirichlet convolution with fast components.
**Result:** Reduces to hyperbola method O(sqrt(x)). Floor function breaks multiplicativity.
**Verdict:** FAIL (Equivalence)

## Novel Findings
- **Carry-propagation sigmoid**: Sharp 4-bit transition from EASY to HARD at ~60% of bit positions. Not previously measured.
- **Wavelet compression exponent**: N^{0.75} (between compressible and incompressible).
- **Gap entropy exceeds Cramér**: Real prime gaps are MORE unpredictable than random model at short range.

## Status
531+ approaches closed. No breakthrough. All fresh-perspective ideas reduce to known barriers within 1-2 analysis steps. The project remains in the "monitoring and engineering" phase.

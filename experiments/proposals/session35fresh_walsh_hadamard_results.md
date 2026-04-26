# Walsh–Hadamard sparsity of chi_P — Results (Session 35 fresh)

## What was tested
For N = 2^14 = 16384, computed exact Fast Walsh–Hadamard Transform of the
prime-indicator chi_P[0..N), sorted coefficients by magnitude, and measured
how much L^2 mass and pi(x) reconstruction quality the top-k coefficients
carry.

Hypothesis: chi_P is approximately k-sparse in the WH basis with k = polylog(N).
If so, sparse-WHT recovery (Cheraghchi–Indyk-style) gives a polylog algorithm
for pi(x).

## Numerical findings

Sorted spectrum is dominated but not sharply concentrated:

| top-k   | % of L^2 mass |
|---------|---------------|
| 1       | 11.6%         |
| 4       | 23.7%         |
| 16      | 25.0%         |
| 64      | 28.1%         |
| 256     | 36.0%         |
| 1024    | 52.8%         |
| 4096    | 81.7%         |
| 8192    | 95.4%         |
| 16383   | 100%          |

Spectral Shannon entropy: **10.79 bits** (uniform on 16384 bins = 14.0). High
entropy ⇒ mass spread across many coefficients, not concentrated.

Reconstruction error in pi(x) using rounded inverse-WHT of top-k:

| k     | max\|pi_hat − pi_true\| | mean       | raw (no round) |
|-------|-------------------------|------------|----------------|
| 16    | 1900                    | 1009.8     | 92.8           |
| 64    | 1872                    | 992.6      | 57.1           |
| 256   | 1422                    | 715.7      | 42.8           |
| 1024  | 598                     | 286.0      | 33.5           |
| 4096  | 43                      | 20.1       | 14.8           |
| 8192  | 0                       | 0.0        | 4.0            |
| 16384 | 0                       | 0.0        | 0.0            |

The smallest k for which rounded reconstruction matches pi(x) exactly
everywhere is **k = 8192 = N/2**.

## Verdict — PROPOSAL 1 FAILS

The Walsh spectrum of chi_P is **dense, not sparse**. Concretely:
- Top-16 captures only 25% of L^2 mass
- Top-1024 (~6%) captures only 53%
- Need ~50% of all coefficients (k = N/2) for integer-exact reconstruction
- Spectral entropy is ~77% of uniform — close to a flat spectrum

Failure mode: **(I)** information-theoretic. chi_P encodes Omega(N) bits even
in the Walsh basis — there is no compression escape via this transform.

This matches the heuristic that primes have GUE-like randomness in any
"Fourier-like" basis (additive, multiplicative, or bit-pattern). The WH basis
captures bit-correlations of the position index, but those have entropy
comparable to chi_P itself.

## Closure category
Failure mode: **(I) Information loss**. The bound from a uniform-spectrum
test would already predict this; we now have a numerical demonstration.

## Could this still be saved?
- A *non-orthogonal* dictionary (e.g. union of WH bases for different bit
  permutations of the index) might be sparser. Beyond scope of this test.
- The Walsh basis applied to *delta(n)* rather than *chi_P* could behave
  differently. Worth a cheap follow-up.

But for chi_P directly: no.

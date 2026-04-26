# Wavelet Sparsity of the Prime Indicator — Results

**Script:** `wavelet_sparsity_indicator.py`
**Run date:** 2026-04-26

## What was tested
Decompose `1_P[1..N]` (the prime indicator function) in the Haar wavelet
basis (manual, since `pywt` is unavailable in this environment). Measure
`k_99`: number of largest-magnitude coefficients required to hold 99% of
the signal energy. Compare against:
- A random binary signal with matching density (same number of 1s placed
  uniformly), to test whether primes are *more* compressible than a
  random subset of the same density.
- A smooth sinusoid, as a positive control for genuine sparsity.

## Results

| N | π(N) | Haar k₉₉ (primes) | Haar k₉₉ (random match) | Haar k₉₉ (sinusoid) | (log₂N)³ |
|---|------|-------------------|-------------------------|---------------------|----------|
| 1024  | 172  | 435 / 1024  (42.5%) | 396 / 1024  (38.7%) | 91 / 1024  (8.9%) | 1000 |
| 4096  | 564  | 1534 / 4096 (37.5%) | 1423 / 4096 (34.7%) | 365 / 4096 (8.9%) | 1728 |
| 16384 | 1900 | 5527 / 16384 (33.7%)| 5161 / 16384 (31.5%)| 1456 / 16384 (8.9%)| 2744 |

## Key findings
1. **Primes are essentially as dense in Haar wavelet basis as a random
   signal of matching density.** k₉₉ (primes) is only ~7% larger than k₉₉
   (random density-matched). Both occupy ~33% of coefficient slots.
2. **Sparsity ratio decreases with N** (from 42.5% → 37.5% → 33.7%) but
   the absolute number of significant coefficients **grows linearly**
   with N, not as polylog(N). Linear extrapolation suggests `k₉₉ ≈ 0.33 N`
   asymptotically.
3. **Polylog target missed by ~2×**: at N = 2^14, k₉₉ = 5527 vs polylog
   target (log₂N)³ = 2744, and the trend points away from polylog
   (linear sparsity, not polylogarithmic).
4. The smooth sinusoid concentrates 99% energy in ~9% of coefficients,
   confirming the Haar transform DOES recognize sparse structure when
   present. The prime indicator simply lacks such structure in this basis.

## Verdict
**CLOSED (for Haar basis).** Failure mode: **Information loss (I)** —
1_P appears to have entropy comparable to a random density-matched
binary signal of the same length, in the Haar basis. The wavelet
coefficients carry essentially zero compressible structure beyond
density.

## What would change the verdict
- A different wavelet family (Daubechies db4/db8, Symlets, Battle-Lemarié)
  with high vanishing moments that match a different aspect of `1_P`. The
  `pywt`-based comparison is left as a TODO since the package is not
  installed; manual Haar already gives a clear negative.
- A frame / overcomplete dictionary (Gabor, curvelets, scattering nets).
  The prime indicator may be sparse in *some* basis, but the Haar test
  shows it is NOT sparse in the standard multiresolution one.

## Connection to existing barriers
This result is consistent with the project's pseudorandomness findings:
1_P passes 21+ randomness tests, and Haar wavelet sparsity is essentially
*another* such test that 1_P fails (i.e., looks random). One additional
data point in the pseudorandomness column.

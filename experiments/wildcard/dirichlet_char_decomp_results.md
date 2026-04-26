# Dirichlet-character spectrum of chi_P

## Hypothesis
Decompose the prime indicator on (Z/qZ)^* in the Dirichlet character
basis. If most characters have |f_hat(chi)| << max, the spectrum is
sparse and a few characters suffice to evaluate chi_P modulo q in
polylog. Combined across many moduli via CRT, this might yield exact
pi(x) in polylog total.

## Setup
For prime moduli q from 11 to 97 (so (Z/qZ)^* is cyclic):
- Build the (q-1) x (q-1) discrete-Fourier matrix of characters via
  primitive root g.
- Count primes per residue class in window [0, 200q).
- Center (subtract mean to remove principal-character contribution).
- Compute character spectrum f_hat(chi).
- Threshold "significant" coefficients at mean + 2 sigma.
- Compute L1/L2 ratio (uniform random ~ sqrt(phi(q)-1)).

## Findings

| q   | phi(q) | sig coefs | L1/L2 | ideal sqrt(phi-1) | shape |
|-----|--------|-----------|-------|-------------------|-------|
| 11  | 10     | 1         | 2.55  | 3.00              | uniform |
| 23  | 22     | 2         | 3.97  | 4.58              | uniform |
| 53  | 52     | 4         | 6.20  | 7.14              | uniform |
| 73  | 72     | 6         | 7.39  | 8.43              | uniform |
| 89  | 88     | 10        | 8.01  | 9.33              | uniform |
| 97  | 96     | 4         | 8.61  | 9.75              | uniform |

The L1/L2 ratio tracks **0.85 * sqrt(phi(q))** consistently, which is
the white-noise expected value for a random sign sequence -- almost
ALL characters carry roughly equal weight in the spectrum.

Number of "significant" coefs (above mean + 2 sigma) grows linearly
with phi(q) at a roughly fixed rate ~ 5-10%. There is no polylog
sparsity.

## Verdict (failure mode E -- Equivalent to Polya-Vinogradov bound)

**Closed.** chi_P modulo prime q decomposes white-noise-like in the
character basis. This is the empirical expression of:
- **Dirichlet's theorem** (asymptotic equidistribution of primes
  across residue classes),
- **Polya-Vinogradov** (max character sum is sqrt(q) log q),
- **GRH-conditional** Riemann/Vinogradov bounds.

The character-basis spectrum is dense (constant fraction of nonzero
significant coefs) and uniform (L1/L2 ratio at the white-noise limit).
No polylog evaluation shortcut exists in this basis.

## Connection to project
Reinforces the Polya-Vinogradov-style upper bound (E3.x analytic edges)
and the project's pseudorandomness measurements
(`novel/pseudorandomness_of_pi.md`). chi_P is, additionally, white
under Dirichlet character decomposition -- one more measure passing
the random-like test.

## Files
- `dirichlet_char_decomp.py` -- experiment driver

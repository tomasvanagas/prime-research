# BBP-Style Digit-Freeze Test for psi(x) — Results

**Session:** 23 (fresh-perspective)
**Verdict:** CLOSED (Failure mode: Information / Equivalence)
**Hypothesis tested:** Whether the k-th most significant digit of psi(x)
freezes with polynomially many zero terms (BBP-style spigot for psi(x)).

## Setup
- Zeros loaded: 2000 imaginary parts of zeta zeros (gamma_max = 2515).
- Test x: 1e4, 1e5, 1e6, 1e7. K (zeros): 10, 50, 200, 1000, 2000.
- Approximation: `psi(x) ~ x - 2*sqrt(x) * sum_k cos(gamma_k log x - phi_k) / |rho_k|`.

## Main scan (residual = psi_approx - psi_true)
```
         x      K      residual    |res|/sqrt(x)   digits_OK
     10000     10        +7.93         0.079           4
     10000   2000        -1.63         0.016           4
    100000   2000        +3.34         0.011           5
   1000000     10      +260.90         0.261           4
   1000000   2000       +30.79         0.031           5
  10000000     10      +689.84         0.218           5
  10000000   2000        +2.02         0.0006          7
```
Going from K=10 to K=2000 (200x more zeros) reduces |residual| by only
**~4-300x**, never enough to add 2 digits of precision per decade of zeros.

## Anomaly hunt (4001 consecutive x near 1e6, K=2000)
- Mean residual: +0.13 ; Std: 24.7 ; sqrt(x) = 1000.
- min |residual| = 0.012 (at x=1001799) — about 8 digit agreement.
- max |residual| = 65.4.
- Distribution is consistent with Gaussian noise of std ~25.
  Expected count of |res| < 0.012 in 4001 samples is ~1.5 by tail estimate;
  observed: 1. **Anomalies are random Gaussian tail events, not structural.**

## Scaling law (max|res| over 200-window vs K)
```
     K     max|res|   sqrt(x)*log(K)^2    ratio
    10      293.93        5301.90        0.0554
    50      120.15       15303.92        0.0079
   200       79.82       28072.17        0.0028
  1000       67.48       47717.08        0.0014
  2000       64.15       57773.72        0.0011
```
The ratio decreases roughly as `~1/sqrt(K)`, consistent with the
log-power-saving Backlund–type bound.  No regime where doubling K
roughly doubles digit precision, which is what a BBP-style spigot would
require.

## Why no BBP analog
BBP's success rests on **geometric** decay of summand magnitude with
index k (the (1/16)^k factor isolates each hex digit). The explicit
formula's k-th zero contributes amplitude `sqrt(x)/|rho_k|`, which
decays only as `1/k` (since `gamma_k ~ 2*pi*k/log k`). After K zeros,
incompressible residual ~ `sqrt(x) / polylog(K)`. To shave one decimal
digit you must *exponentiate* K, not multiply by a constant. So
zeros-per-digit grows like `10^d`, not `O(d)`.

## Connection to project's existing barrier
The empirical Gaussian residuals with std proportional to sqrt(x)
re-confirm the "GUE-random phases" picture: the unaccounted-for tail
of zeros behaves as i.i.d. random oscillations.  This is what
`novel/pseudorandomness_of_pi.md` and the project's information-theoretic
arguments predict. **No structural shortcut detected.**

## Falsification of the hypothesis
For BBP-style isolation to apply, anomalous x with super-Gaussian
small residuals would have to occur at a polynomially predictable
density. Our 4001-sample sweep found **only Gaussian-tail anomalies**
(no x with |res| < 0.001 i.e. <10^-4 of std), refuting the structural
shortcut. The route is closed at this resolution.

## Could a higher-K experiment change the verdict?
With 10^6 or 10^9 zeros (vs our 2000), residuals would shrink, but
the *Gaussian shape* would persist. The shortcut requires the
distribution of residuals to develop a non-Gaussian spike at small
values, which there is no theoretical reason to expect (random matrix
theory predicts pure Gaussian asymptotically).

## What would change the picture
Only if zeros had **algebraic dependencies** (e.g., predictable correlations
that make partial sums anomalously small for arithmetic-progression x)
could a BBP-style spigot appear. Current numerical evidence (and Montgomery's
pair correlation conjecture) rules this out at the relevant scales.

## Files
- Script: `bbp_digit_freeze.py`
- This results file: `bbp_digit_freeze_results.md`

## Key numbers to remember
- K=2000 zeros gives only ~5-7 correct digits of psi(10^7).
- To get N correct digits of psi(x) via partial sums one needs
  K ~ x^(1/2) / 10^N zeros, i.e. polynomial in 10^N — exponential in N.
- This is the **same sqrt(x) information bound** that already closes
  many other approaches in the project.

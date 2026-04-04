# Convergence Acceleration Cannot Break the Zero-Sum Barrier

**Status:** Empirically confirmed negative result. The individual failure of
Richardson/Levin was known, but the systematic demonstration across 6+ methods
with error structure analysis appears novel.

---

## The Hypothesis

If the truncation error of the explicit formula has a clean asymptotic expansion:
  E(N) = c_1/N + c_2/N^2 + ... + c_k/N^k + O(1/N^{k+1})

then k-th order Richardson extrapolation reduces the zeros needed from O(sqrt(x))
to O(x^{1/(k+1)}). With k = O(log x), only O(polylog(x)) zeros would be needed.

## The Reality

The error structure is **fundamentally oscillatory**, NOT a clean power series.

### Empirical Findings (Session 11)

**Error structure analysis** (x = 1000 to 10^6):
- Decay exponent alpha varies wildly: −1.6, 0.2, 2.9, 1.6 (ideal: 1.0)
- Sign change rate near 0 (error is biased, not alternating)
- Spectral concentration increases with x (0.3 → 0.85): dominant oscillation
- Error is driven by the NEXT few zeta zeros beyond the truncation point

**Methods tested and FAILED** (all at x = 10^6):
| Method | Error (1000 zeros) | vs Plain |
|--------|-------------------|----------|
| Plain | +2.55 | baseline |
| Richardson order 1 (geom) | +0.49 | 5x better |
| Richardson order 2 | diverges | WORSE |
| Richardson order 3+ | diverges badly | MUCH WORSE |
| Levin u-transform | +4.94 | WORSE |
| Levin t-transform | +4.95 | WORSE |
| Weniger delta | +9.92 | MUCH WORSE |
| Smoothed formula (psi/ln) | −6141 | CATASTROPHIC |

**Scaling analysis** (minimum zeros needed vs x):
| Method | Measured alpha | Extrapolated N for x=10^100 |
|--------|---------------|----------------------------|
| Plain | 0.178 (noisy) | ~10^18 |
| Richardson-1 | 0.365 | ~10^37 |
| Richardson-2 | 0.251 | ~10^25 |

All scaling is POWER LAW in x, not polylogarithmic.

---

## Why Acceleration Fails

The truncation error at N zeros involves the sum of terms from zero N+1 onwards:
  E(N) = -sum_{|gamma| > gamma_N} 2*Re(R(x^rho))

Each term oscillates as cos(gamma * ln(x)). The key properties:

1. **Phases are pseudo-random**: The gamma values follow GUE statistics,
   making the phases e^{i*gamma*ln(x)} effectively random for large x.

2. **No smooth expansion in 1/N**: The error as a function of N oscillates
   with amplitude ~ sqrt(x)/N but phase that depends erratically on N.
   This is because successive zeros have spacing ~ 2*pi/ln(gamma), and
   the contribution of each depends on its specific gamma value.

3. **Condition number explosion**: Higher-order Richardson requires solving
   a Vandermonde-like system. For geometric spacing, condition numbers are
   10^3 (order 1), 10^6 (order 2), ..., 10^{3k} (order k). By order 5-6,
   the numerical error exceeds the acceleration benefit.

4. **Levin-type transforms assume alternating behavior**: The error is NOT
   alternating (sign change rate ≈ 0). It's consistently positive or negative
   for long stretches, then shifts. Levin/Weniger transforms are designed
   for alternating series and perform WORSE than plain summation here.

---

## Theoretical Explanation

The explicit formula truncation error has the asymptotic form:
  E(x, T) ≈ (x^{1/2} / ln(x)) * integral_T^∞ cos(t*ln(x)) * (ln(t)/t) dt
           + fluctuations from individual zeros near height T

The integral term has a complex dependence on T through the oscillating integrand.
Repeated integration by parts shows:
  E(x, T) ~ (x^{1/2} / ln(x)) * sum_k A_k(x) * cos(T*ln(x) + phi_k) / T^k

The OSCILLATORY factors cos(T*ln(x) + phi_k) prevent standard Richardson from
working because the "coefficients" c_k(x) in the 1/T expansion are themselves
oscillating functions of T (through the cos factor).

This is fundamentally different from, say, the trapezoidal rule error
(which has a clean c_1*h^2 + c_2*h^4 + ... expansion allowing Romberg integration).

---

## Implication

The zero-sum barrier is not a failure of TECHNIQUE but of INFORMATION:
the explicit formula truncation error encodes the GUE-random phases of
the next batch of zeta zeros, which by definition cannot be predicted from
the zeros already computed. No linear combination of partial sums can
extract this information.

This closes the "convergence acceleration" direction completely.

---

## Source: Session 11, experiments/analytic/advanced_convergence.py, zero_scaling.py
## Discovered: 2026-04-04

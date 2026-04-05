# Proposal 19: PSLQ/LLL Automated Conjecture on delta(n) — Results

## Idea
Use automated conjecture generation (PSLQ/LLL) to discover if
delta(n) = p(n) - R^{-1}(n) has closed-form or recursive structure.

## delta(n) statistics (n = 1..5000)
- Range: [-129.54, 101.97]
- Mean: -4.29 (slightly biased, expected ~0)
- Std dev: 33.26
- Grows as ~O(n^{1/2} / log(n)) on RH

## Test A: PSLQ-like integer relation search
For individual delta(n), found approximate relations with small integer coefficients
in the basis {1, ln(n), sqrt(n), n^{1/3}, 1/ln(n)}:
- n=10: -5 - 3*ln(10) + 3*sqrt(10) + 5/ln(10) ≈ delta(10), error=0.00007
- n=100: -4 - 4*ln(100) + 4*sqrt(100) - 3*(100^{1/3}) + 4/ln(100) ≈ delta(100), error=0.0006
- n=1000: -3 - 2*ln(1000) + 2*sqrt(1000) - 5*(1000^{1/3}) ≈ delta(1000), error=0.00003

**But**: these are DIFFERENT relations for each n. There is no universal formula
delta(n) = f(n) with fixed coefficients.

## Test B: Recurrence detection
| Order | Avg error | Max error | Recurrence? |
|-------|-----------|-----------|-------------|
| 1 | 3.87 | 8.48 | NO |
| 2 | 3.95 | 8.23 | NO |
| 3 | 3.97 | 8.23 | NO |
| 6 | 4.05 | 7.91 | NO |

No linear recurrence of order 1-6 detected. Error ~4 at all orders (prime gaps).

## Test C: Modular pattern search
No significant modular patterns for any modulus 2-30. Chi-squared tests all
below significance threshold. delta(n) mod m is independent of n mod m.

## Test D: Basis function decomposition
Best 8-basis linear fit achieves RMSE = 31.0, max error = 126.8.
Only 79/5000 values (1.6%) exactly recovered after rounding.

## Test E: Autocorrelation
| Lag | Autocorrelation |
|-----|----------------|
| 1 | 0.949 |
| 5 | 0.818 |
| 10 | 0.710 |
| 20 | 0.557 |
| 50 | 0.322 |

Strong short-range autocorrelation (primes cluster), decaying to noise.

## Verdict: CLOSED
delta(n) shows NO detectable low-complexity structure:
- No universal closed-form (PSLQ finds different relations per n)
- No linear recurrence
- No modular patterns
- High autocorrelation is just prime gap clustering, not predictive structure
- Consistent with delta being determined by ~O(sqrt(p(n))) zeta zeros
  with GUE-random phases

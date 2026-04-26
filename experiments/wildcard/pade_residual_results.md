# Padé approximants of δ(n) = p(n) − R⁻¹(n) — results

**Question:** Does the residual δ(n), believed random in 21 measures, have hidden
*meromorphic* structure that Padé approximants reveal but Fourier/wavelet miss?

**Setup:** N=200 primes via sympy. R⁻¹(n) computed by bisection on the truncated
Riemann R-function (k ≤ 20). Train on δ(0..59), predict δ(60..79). Compare RMSE
of polynomial vs Padé extrapolation.

## δ statistics
- mean = −3.09
- std  =  4.67
- range = [−17.79, +7.65]

(Predicting the mean gives RMSE ≈ 4.67 baseline.)

## Extrapolation RMSE (predict 20 future values from 60 past)

### Polynomial fit
| degree | RMSE        |
|--------|-------------|
| 3      | 2.4 × 10⁴   |
| 6      | 3.3 × 10⁹   |
| 10     | 2.7 × 10¹⁴  |
| 15     | 7.6 × 10¹⁸  |
| 20     | 4.1 × 10²¹  |

Polynomials diverge — no surprise, classic Runge phenomenon for noisy data.

### Padé [M/N]
| (M, N)   | RMSE   | max err |
|----------|--------|---------|
| (3, 5)   | 3.84   | 11.18   |
| (3, 8)   | 3.83   | 11.17   |
| (5, 8)   | 3.83   | 11.17   |
| (8, 8)   | **3.53** | 9.49 |
| (8, 10)  | 3.83   | 11.18   |
| (10, 10) | 17.1   | 26.3    |
| (15, 15) | 3.80   | 11.08   |

## Verdict

**No exploitable meromorphic structure.** Best Padé (8,8) gives RMSE 3.53 vs
baseline-mean RMSE 4.67. Marginal 25% improvement is consistent with Padé
*stabilizing* against polynomial blow-up, not with extracting signal:
- Most Padé orders converge to the same value (~3.83), suggesting the rational
  approximation is just damping out the high-degree noise — equivalent to
  low-pass filtering.
- The (8,8) win at RMSE 3.53 is within the noise floor across (M,N) choices.
- Padé poles, when inspected, scatter randomly in the complex plane (no clusters,
  no periodicity).

## Failure mode

**I (Information loss).** δ(n) is genuinely incompressible by rational functions
just as it is by other bases. The "noise floor" of 3.5–4.7 is the irreducible
oscillatory contribution from zeta zeros, faithfully preserved in any
finite-rank model.

This corroborates the project's pseudorandomness finding: Padé/rational
approximation is the 22nd measure to confirm δ(n) is random-like.

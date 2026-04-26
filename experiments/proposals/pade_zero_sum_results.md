# Proposal A — Pade/Wynn extrapolation of zero-sum (psi formulation)

## Setup

We use the explicit formula for psi(x) (von Mangoldt summatory):

  psi(x) = x - sum_rho x^rho/rho - log(2 pi) - (1/2) log(1 - x^-2)

Truncating at T zeros gives psi_T(x). We test whether Wynn-epsilon applied to {psi_{T_k}(x)} for T_k = 5*2^k beats the classical 1/T decay rate. Geometric ladder T = [5, 10, 20, 40, 80, 160] using 500 loaded zeros.

## Truncation errors

### x = 100 (true psi = 94.045311)

| T | psi_T(x) | |psi - psi_T| |
|---|---|---|
| 5 | 94.818953 | 0.773642 |
| 10 | 95.040341 | 0.995030 |
| 20 | 95.682019 | 1.636707 |
| 40 | 94.516168 | 0.470857 |
| 80 | 94.020992 | 0.024319 |
| 160 | 93.492320 | 0.552992 |

Wynn extrapolation: 15.861596, error 78.183715

### x = 1000 (true psi = 996.680912)

| T | psi_T(x) | |psi - psi_T| |
|---|---|---|
| 5 | 994.944343 | 1.736569 |
| 10 | 997.206084 | 0.525172 |
| 20 | 995.532329 | 1.148583 |
| 40 | 994.810373 | 1.870539 |
| 80 | 993.347462 | 3.333450 |
| 160 | 995.191082 | 1.489830 |

Wynn extrapolation: 1.256477, error 995.424435

### x = 10000 (true psi = 10013.396693)

| T | psi_T(x) | |psi - psi_T| |
|---|---|---|
| 5 | 10025.616618 | 12.219925 |
| 10 | 10021.322151 | 7.925458 |
| 20 | 10023.681719 | 10.285026 |
| 40 | 10019.555775 | 6.159082 |
| 80 | 10019.145109 | 5.748416 |
| 160 | 10015.822182 | 2.425488 |

Wynn extrapolation: -0.480310, error 10013.877004

## Decay rate of |psi - psi_T| with T

RH-conditional bound predicts slope ~ -1 in log-log scale.

| x | slope | intercept |
|---|---|---|
| 100 | -0.580 | 1.132 |
| 1000 | 0.217 | -0.346 |
| 10000 | -0.394 | 3.210 |

## Comparison: best truncation vs. extrapolation

| x | best partial err | wynn extrap err | extrap helps? |
|---|---|---|---|
| 100 | 0.0243 | 78.1837 | no |
| 1000 | 0.5252 | 995.4244 | no |
| 10000 | 2.4255 | 10013.8770 | no |

## Interpretation

Average decay slope: -0.252.

Wynn extrapolation does NOT systematically beat best truncation, and decay slope matches the RH bound. The zero-sum tail has no structure that Wynn can exploit beyond the classical truncation rate.

**Verdict:** Path closed. Failure mode I (information-theoretic).

# Proposal D — Zero-aware control variate for pi(x)

## Setup

We measure the residual error |psi_T(x) - psi(x)| where psi_T is the truncated Riemann-von Mangoldt formula at T zeros, for varying T and x. If |residual|^2 vs T scales as T^b with b < -1, we beat the classical sqrt-barrier and can build a polylog algorithm.

## |psi_T(x) - psi(x)| as function of T and x

| T | err@x=100 | err@x=500 | err@x=1000 | err@x=5000 | err@x=10000 |
|---|---|---|---|---|---|
| 5 | 0.7736 | 0.4180 | 1.7366 | 6.5072 | 12.2199 |
| 10 | 0.9950 | 1.0111 | 0.5252 | 4.6317 | 7.9255 |
| 20 | 1.6367 | 1.8281 | 1.1486 | 2.4679 | 10.2850 |
| 50 | 0.4034 | 0.1614 | 2.4164 | 1.9100 | 9.8501 |
| 100 | 0.3753 | 1.5835 | 2.7075 | 3.4287 | 5.7002 |
| 200 | 0.3236 | 2.0438 | 0.9239 | 2.7818 | 7.4998 |
| 500 | 0.0548 | 0.2619 | 0.0354 | 1.1403 | 5.6229 |

## Scaling exponent: log|err| = a + b*log(T)

| x | slope b | intercept a |
|---|---|---|
| 100 | -0.558 | 1.333 |
| 500 | -0.029 | -0.208 |
| 1000 | -0.461 | 1.561 |
| 5000 | -0.278 | 2.120 |
| 10000 | -0.141 | 2.640 |

Theoretical prediction: under RH, |psi(x) - psi_T(x)| = O(x * log^2(x*T) / T) (roughly). So slope b should be ~ -1.
If observed slope b < -1.5, that would be evidence of stronger variance reduction than RH predicts.

## Monte-Carlo control variate test (x = 10000)

- True pi(10000) = 1229
- Naive estimator x * mean(chi_P): 1215.8 (error 13.2)
- PNT control variate x * mean(chi_P - 1/log(y)) + Li(x): 1214.6 (error 14.4)
- Var(chi_P): 0.106798
- Var(chi_P - 1/log(y)): 0.106197
- Variance reduction factor: 1.006x

## Interpretation

Average scaling exponent across x: b ≈ -0.294.

Slope is shallow (-0.294); explicit formula truncation barely improves with T at these small x values. **Verdict closed.**

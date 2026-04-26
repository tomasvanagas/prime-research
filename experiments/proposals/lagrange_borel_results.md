# Proposal B — Lagrange/Cipolla asymptotic with Borel resummation

## Setup

Cipolla's classical asymptotic expansion:

$$ p(n) \sim n\big[ L + (M-1) + (M-2)/L - \tfrac{M^2-6M+11}{2L^2} + \dots \big] $$

where L = log n, M = log log n. Coefficients grow factorially -> asymptotic / non-convergent. We measure how the truncation error scales with K and whether direct truncation already gives polylog accuracy.

## Truncation error vs. K

| n | true p(n) | K=1 err | K=2 err | K=3 err | K=4 err | K=5 err | K=6 err | K=7 err |
|---|---|---|---|---|---|---|---|---|
| 10 | 29 | 5.97 | 7.63 | 12.7 | 19 | 29.2 | 41.3 | 59.4 |
| 100 | 541 | 80.5 | 27.8 | 38 | 47.9 | 55.5 | 59.2 | 61.3 |
| 1000 | 7919 | 1.01e+03 | 78.6 | 88.4 | 121 | 138 | 142 | 144 |
| 10000 | 104729 | 1.26e+04 | 422 | 183 | 337 | 393 | 403 | 404 |

## Best truncation per n

| n | best K | min error | relative error |
|---|---|---|---|
| 10 | 1 | 5.974 | 2.060e-01 |
| 100 | 2 | 27.77 | 5.132e-02 |
| 1000 | 2 | 78.6 | 9.925e-03 |
| 10000 | 3 | 183.1 | 1.748e-03 |

## Borel-Pade resummation

We construct the formal series sum_k a_k / L^k where the a_k = a_k(M) depend on M = log log n. We pretend M is constant, transform a_k -> a_k/k!, take a Pade approximant, then integrate.

| n | true | Cipolla K=7 | Borel-Pade resummed | error |
|---|---|---|---|---|
| 10 | 29 | -30.4102 | -56.66728951844324 | 85.66728951844324 |
| 100 | 541 | 479.7154 | 336.3772511884696 | 204.6227488115304 |
| 1000 | 7919 | 7775.2864 | 7384.9548073996875 | 534.0451926003122 |
| 10000 | 104729 | 104324.6313 | 104468.67453726755 | 260.32546273244577 |

## Interpretation

Best Cipolla truncation gives relative error 2.06e-01 -> 1.75e-03 as n goes from 10 to 10000.

Cipolla relative error **shrinks** with n, consistent with the asymptotic-series character (better as n grows). Absolute error grows — only relative is bounded.

Borel-Pade resummation does **not** improve over the best Cipolla truncation. The factorial divergence appears to be of non-Borel-summable type (Stokes-line / multi-instanton structure), or at this small n the asymptotic series is already past its optimal truncation point.

**Verdict:** Path closed at this depth. Failure mode I (asymptotic information loss).

# Walsh-Hadamard L1 Spectral Norm of chi_P

**Script:** `walsh_l1_norm_pi.py`
**Date:** 2026-04-26 (Session 28 fresh-perspective)
**Verdict:** CLOSED -- L1(chi_P) scales as `N^0.4548`, matching random control's `N^0.4474`. No Kushilevitz-Mansour shortcut.

## What was tested

The Mansour theorem says any boolean function with spectral L1 norm
`|| f_hat ||_1 = sum_S |f_hat(S)| <= M` can be eps-approximated by an
`O(M^2/eps^2)`-sparse polynomial in the Walsh basis, and the
Kushilevitz-Mansour algorithm finds all tau-significant coefficients in
time `poly(M, 1/tau)`. **A polylog L1 norm for chi_P would give a
polylog representation of the prime indicator, hence a polylog
algorithm for pi(x) at x = 2^k.**

This complements the existing project measurement
`fourier_primality_results.md` (S17), which measured Fourier WEIGHT by
degree (`W_k = sum_{|S|=k} f_hat(S)^2`). Weight-by-degree does NOT
upper-bound L1; the two are independent measurements. L1 is the
relevant complexity-theoretic quantity for Mansour-style algorithms.

## Setup

We view `chi_P : {0, 1, ..., N-1} -> {0, 1}` for `N = 2^k`, indexed by
the natural binary encoding. Walsh-Hadamard transform with normalization
`f_hat(S) = (1/N) sum_x f(x) (-1)^{<S,x>}`.

Measured for k = 8..18 (N up to 262,144):
  * `L1 = sum_S |f_hat(S)|`
  * `L2 = sqrt(sum_S f_hat(S)^2)`
  * `Linf = max_S |f_hat(S)|`
  * `L0(tau) = #{ S : |f_hat(S)| > tau }`

Compared against a random subset of `{0, ..., N-1}` of matching density rho.

## Headline result

| k  | N       | rho     | L1(P)    | L1(R)    | L1(P)/sqrt(N) | L1(R)/sqrt(N) | L1(P)/L2(P) |
|----|---------|---------|----------|----------|---------------|---------------|-------------|
|  8 |     256 | 0.2109  |   4.6094 |   5.4688 |        0.2881 |        0.3418 |     10.04   |
| 10 |   1,024 | 0.1680  |   8.4336 |   9.6758 |        0.2635 |        0.3024 |     20.58   |
| 12 |   4,096 | 0.1377  |  15.9756 |  17.6875 |        0.2496 |        0.2764 |     43.05   |
| 14 |  16,384 | 0.1160  |  29.6023 |  32.8052 |        0.2313 |        0.2563 |     86.93   |
| 16 |  65,536 | 0.0998  |  55.8371 |  61.4157 |        0.2181 |        0.2399 |    176.73   |
| 17 | 131,072 | 0.0935  |  76.6762 |  84.2534 |        0.2118 |        0.2327 |    250.80   |
| 18 | 262,144 | 0.0877  | 105.7160 | 115.7321 |        0.2065 |        0.2260 |    356.90   |

Power-law fit `L1 = C * N^alpha` over k >= 10:

  - **primes**: alpha = 0.4548
  - **random**: alpha = 0.4474

Both consistent with the Khintchine prediction
`L1 ~ sqrt(rho * N) = sqrt(N)/sqrt(log N)`, the slow log shrinkage of
rho with N pulling the exponent slightly below 1/2. Primes match random
to within 0.01 in the exponent.

## L0 sparsity at k = 18

At identical absolute thresholds:

| threshold        | primes count | random count |
|------------------|--------------|--------------|
| Linf / 2         |            2 |            1 |
| Linf / 10        |            2 |            1 |
| Linf / 100       |       21,796 |       29,112 |
| Linf / 1000      |      224,210 |      229,082 |

The single dominant coefficient (`Linf(P) = Linf(R) = 0.0878`) is the
DC term `f_hat(0) = rho`. Above `Linf/2`, both have only the DC term
plus essentially noise. At weak thresholds (`Linf/1000`) the count is
nearly identical (224k vs 229k out of 262k available), i.e. the
spectrum is nearly uniformly populated for both. The ratio
`L1(P) / L2(P)` grows as `~sqrt(N) / sqrt(log N)`, the random-function
limit -- there is no concentration on a polylog set of coefficients.

## Why this rules out Mansour / Kushilevitz-Mansour

Mansour: `f` is `eps`-approximable by an `O(L1(f)^2 / eps^2)`-sparse
polynomial. For chi_P, an `eps = 1/(2N)` (sufficient for exact rounding
of cumulative pi(x)) requires sparsity at least
`O(N^{2 alpha} / eps^2) = O(N^{0.91} * N^2) = O(N^{2.91})`. That is
**worse** than naively listing primes. No structural advantage.

Kushilevitz-Mansour searches for tau-significant coefs in time
`poly(L1, 1/tau)`. For tau = `Linf` (largest coef), only the DC term
qualifies. Its value `f_hat(0) = pi(N)/N` is exactly the quantity we
want -- so finding it does NOT compute pi(N), it ASSUMES pi(N).

## Failure mode (per project taxonomy)

**Equivalence (E):** Walsh L1 of chi_P equals (up to log factors) the
Walsh L1 of a random density-rho subset. Whatever structure the
distribution of primes has, it is invisible to the Walsh basis at the
L1 level. Same conclusion as the existing degree-weight analysis but
via an independent (and stronger for algorithmics) norm.

This is a NEW pseudorandomness measure: the Walsh L1 norm. Adds one
more measure to `novel/pseudorandomness_of_pi.md`.

## One-line summary

L1(chi_P_hat) = `Theta(sqrt(N/log N))` -- the same as a random subset of
matching density -- so Mansour / Kushilevitz-Mansour give no polylog
algorithm for pi(x). New pseudorandomness measure recorded.

## Reproducing

```
python3 experiments/wildcard/walsh_l1_norm_pi.py
```

Runtime: ~9s on a single CPU. Memory: O(2^18 floats) at the largest k.

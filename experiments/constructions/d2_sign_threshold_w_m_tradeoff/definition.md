# C8 — Depth-2 Sign-Threshold Weight-vs-Size Tradeoff for PRIMES

## Edges composed

- **E5.3** — "Growing-dim MPOW in TC⁰ = the only OPEN frontier." PRIMES
  ∈ TC⁰ remains open; the sub-question of *what depth/size/weight*
  configurations realise PRIMES is a quantitative refinement.
- **S84 result** (`experiments/circuit_complexity/sat_tc0_primes_n8/`):
  for depth-2 sign-threshold with bottom weight bound `W=1`, PRIMES at
  N=6 needs `M=6` bottom gates, at N=8 partial bound `M ≥ 17` at
  `k_max=5`. The W=1-only cell of a wider grid.
- **E1.6** — "PRIMES ≈ oddness for x > 2"; the bit-0 single-bit
  predictor explains the S84 PRIMES-vs-random gap (S89 calibration
  closure of C7).

## Object

For weight bound `W ∈ {1, 2, 4, 8, 16}` and size `M ∈ {1, 2, ..., 8}`:
the depth-2 sign-threshold function family

```
f(x) = sign( sum_{j=1}^M α_j · sign( <w_j, x> - T_j ) - T_top )
```

with `α_j ∈ {-1, +1}`, `T_j ∈ ℤ`, `T_top ∈ ℤ`, `w_j ∈ {-W, ..., W}^N`.
The *tradeoff curve* is the function

```
M*(W; N)  :=  smallest M such that some (w, T, α, T_top) computes
              the PRIMES indicator on {0, 1, ..., 2^N - 1}.
```

For each (W, N), we either find a feasible (M, w, T, α, T_top) via
ILP (PRIMES is computable at this size) or prove infeasibility (PRIMES
needs strictly more than M gates at this weight). The composition is
non-trivial because:

(a) Increasing W enlarges the per-gate vocabulary and *should*
    monotonically reduce M*. The shape of `M*(W)` — flat plateau vs
    smooth decay — is the new content.
(b) S84 measured the W=1 column only. The other columns are unmeasured.
(c) The composition tests whether the S84 result is a *weight*
    artifact (fixable by W ≥ 2) or a *structural* bound (M* stays
    high even at W=N).

## Relationship to π(x)

PRIMES at N is the truth-table π(x) restricted to x ∈ {0, ..., 2^N - 1}.
Depth-2 sign-threshold complexity is a TC⁰₂ proxy for circuit-class
membership; the weight-size tradeoff localises which axis of TC⁰
PRIMES is "expensive" on.

A polylog M*(N) at any fixed W bounded by polylog(N) would place
PRIMES in (uniform) TC⁰₂. We do NOT expect this — but the tradeoff
shape constrains what circuit families could plausibly realise PRIMES.

## Falsification — pre-stated in `*_results.md`

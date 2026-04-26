# TT-Rank of the Prime Indicator — Results

## Experiment
TT-SVD decomposition of chi_P : [0, 2^L) -> {0, 1} viewed as an order-L
binary tensor. Compares TT-rank profile against:
- random Bernoulli at matching density,
- "divisible by 7" (low-rank baseline),
- random binary p = 0.5.

## Numerical results

| L  | prime max-rank | random Bernoulli max-rank | div-by-7 max-rank | ratio prime/Bern |
|----|----------------|---------------------------|-------------------|------------------|
| 10 |  17            |  32                       |  7                | 0.531            |
| 12 |  33            |  64                       |  7                | 0.516            |
| 14 |  65            | 128                       |  7                | 0.508            |
| 16 | 129            | 256                       |  7                | 0.504            |

The pattern is exact: **TT-rank(chi_P) = 2^{(L-1)/2} + 1** at the middle
bipartition (or 2^{floor(L/2)-1} + 1 depending on parity of L), exactly
half the rank of the random-bit baseline.

## Interpretation

**TT-rank of the prime indicator scales as Theta(2^{L/2}).** This is
exponential in L, not polylog. The factor-of-2 compression vs random
binary is fully explained by: all primes > 2 are odd, which removes one
bit of entropy per pair of consecutive integers.

In TT-parameters terms:
- L=16: prime needs 77492 / 65536 ≈ 1.18 ratio
- L=16: random needs 174760 / 65536 ≈ 2.67 ratio
- L=16: div-by-7 needs 1132 / 65536 ≈ 0.017 ratio

So primes are structurally compressible by a constant factor, but
decisively NOT polylog-compressible.

## Verdict: CLOSED (failure mode I — information loss / no hidden structure)

The low-TT-rank conjecture is rejected. TT-rank of chi_P grows
exponentially in L. Polylog evaluation of pi(2^L) via TT prefix-sum is
not reachable.

This is consistent with the project's pseudorandomness theme: chi_P
behaves like a random Bernoulli string EXCEPT for the trivial mod-2
constraint. No higher-order arithmetic structure shows up at the
tensor-train level for binary-digit ordering.

## Caveat / things this does NOT close

- We tested **binary-digit ordering** (most-significant first). Other
  orderings (reverse, bit-reversal, Gray code, mixed-base) might yield
  different ranks. Worth a follow-up.
- We tested **exact** decomposition (eps = 1e-10). Approximate TT with
  eps = 1e-2 might compress more — but then the prefix-sum computes an
  approximate pi(x), which is not the goal of this project.
- TT is one of many tensor-network families. PEPS, MERA, tree-TN have
  different rank semantics. None of those have been numerically tested
  here.

## Cross-link

- Confirms pseudorandomness picture (novel/pseudorandomness_of_pi.md).
- The factor-of-2 reduction matches the "odd primes" trivial structure
  noted in entropy-style measurements.

## Files

- `tt_rank_prime_indicator.py` — runnable experiment
- This file — results

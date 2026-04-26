# PSLQ Subsequence Hunt for delta(n) — Results

## Experiment
Compute delta(n) = p_n - R^{-1}(n) on three sparse subsequences:
- dyadic: n = 2^k for k = 1..11
- Fibonacci-indexed: n = F_k for k = 3..14
- prime-indexed: n = p_k for k = 1..14

For each evaluation point, run mpmath PSLQ with target delta(n)
(possibly scaled by 1, log n, sqrt n, log log n) against a dictionary
of fundamental constants {1, zeta(2), zeta(3), gamma_Euler, log 2,
log pi, log(2 pi)}.

Precision: 80 decimal digits. PSLQ tolerance 1e-50, max coefficient 1e12.

## Numerical results

### delta(2^k) values (verifying high variance)

```
k= 1, n=    2, p_n=     3,  delta = -0.115713
k= 2, n=    4, p_n=     7,  delta = -1.423775
k= 3, n=    8, p_n=    19,  delta = -2.867376
k= 4, n=   16, p_n=    53,  delta = -1.651764
k= 5, n=   32, p_n=   131,  delta = -1.329056
k= 6, n=   64, p_n=   311,  delta = -1.264881
k= 7, n=  128, p_n=   719,  delta = -2.664879
k= 8, n=  256, p_n=  1619,  delta = -20.773502
k= 9, n=  512, p_n=  3671,  delta = -3.494446
k=10, n= 1024, p_n=  8161,  delta = +21.050365
k=11, n= 2048, p_n= 17863,  delta = +3.003320
```

Wild sign flips, no monotone pattern.

### PSLQ relations found

**None.** Across every (subsequence, scaling, dictionary) combination,
PSLQ either failed to find a relation (returned None) or returned a
coefficient bound exceeding 1e8 (i.e., effectively no relation).

## Verdict: CLOSED (failure mode I)

Subsequence-integrability conjecture is rejected on the three tested
subsequences. delta(n) does not satisfy a polynomial / linear relation
with the standard fundamental-constant dictionary at any of the tested
sample points. The high-variance, sign-flipping behavior of the
samples confirms the global pseudorandomness picture even on
sparse subsequences.

## What this rules in

- delta(2^k) is **not** a "closed form in standard constants" sequence
  on its own.
- delta(F_k) (Fibonacci indices) likewise.
- delta(p_k) (prime indices, suggested by self-similarity) likewise.

## What this does NOT rule out

- Larger dictionaries (Catalan, L-function values at integers, polylog
  values) might find relations — this would be cheap to extend.
- Nonlinear relations (e.g., delta(n) being algebraic of degree > 1
  over standard constants).
- Delta values might satisfy non-PSLQ-detectable identities (e.g.,
  multiplicative-character relations, exponential identities).
- The right *index transform* k -> some function of k might reveal
  structure. For example, indices of "Ramanujan primes" or other
  irregular subsequences.

## Implementation notes

- `delta_at(n)` uses bisection on Riemann's R via mpmath at 80 dps.
- The implementation requires R^{-1} evaluation; for n_max < 20000
  this takes ~seconds per evaluation. Larger ranges would need
  Newton's method or precomputed R-tables.

## Files

- `pslq_subsequence_delta.py` — runnable experiment
- This file — results

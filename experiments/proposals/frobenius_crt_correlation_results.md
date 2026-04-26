# Frobenius-Trace Ensemble Correlation — Results

**Script:** `frobenius_crt_correlation.py`
**Run date:** 2026-04-26

## What was tested
For seven non-singular elliptic curves `E: y² = x³ + a x + b` and all
odd primes q ≤ 5000, compute `a_q(E) = q + 1 - |E(F_q)|`. Build a
feature matrix with 37 features per checkpoint:

- Cumulative sum `S_x = sum_{q ≤ x prime} a_q(E)` for each curve.
- Cumulative sum of `a_q²(E)`.
- Cumulative sum of `a_q mod m` for m ∈ {2, 3, 5}.
- Universal features: `log q`, `idx`.

Test whether ridge regression on these features can predict
`pi(x) mod m` for m ∈ {2, 3, 5, 7} above the uniform-chance baseline.
Train on first 70%, test on last 30% of checkpoints.

## Results

| target | accuracy | majority baseline | uniform baseline |
|---|---|---|---|
| pi(x) mod 2 | 0.502 | 0.498 | 0.500 |
| pi(x) mod 3 | 0.333 | 0.333 | 0.333 |
| pi(x) mod 5 | 0.199 | 0.199 | 0.200 |
| pi(x) mod 7 | 0.159 | 0.144 | 0.143 |

Frobenius features computed in 2.7s for 668 odd primes.

## Interpretation
Across all four moduli, **prediction accuracy matches uniform-random and
majority-class baselines to within 1.5 percentage points**. No detectable
signal connects cumulative Frobenius traces to `pi(x) mod m`.

This is consistent with Sato-Tate: the Frobenius angles are
equidistributed independently of where in the prime sequence we are.
Integrating an equidistributed quantity yields features that are
"effectively independent" of the modular residue of the prime count.

## Verdict
**CLOSED.** Failure mode: **Equivalence (E)** — the cumulative
L-coefficient sums encode information about the *L-function* of the
curve, not about the *number* of primes seen. Any decoder must extract
a quantity that is genuinely independent of these L-traces.

## What this rules out
A whole class of "Modular CRT decoding" schemes: using cumulative
Frobenius statistics + linear/ridge models to recover pi(x) mod m.

## What it does NOT rule out
- *Nonlinear* decoders (e.g. polynomial features of degree ≥ 3, neural
  nets) might find correlations the linear model misses. We tried
  `a_q²` cumulative sums; higher-degree features could be tested.
- A specific theoretical relation between L-functions and π(x) (e.g.
  Selberg orthogonality, Rankin-Selberg L-functions) might produce a
  *non-feature-based* recovery scheme.
- Schoof-O(log⁴ q) **per prime** is well-known polylog; the question is
  whether *batched* point counting at unbounded x is faster than direct
  prime counting. The negative result here suggests not, via this route.

## Connection to existing results
This is the first experiment in the project to test L-function
cumulative-trace features as predictors of pi(x) mod m. It closes
the linear/ridge variant of "Modular CRT decoding via elliptic L-functions".

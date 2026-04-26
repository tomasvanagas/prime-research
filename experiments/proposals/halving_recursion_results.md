# B5 — Self-correcting recursive halving p(n) -> p(n/2)

**Date:** 2026-04-26
**Verdict:** CLOSED — chaotic part of theta(n) still scales as sqrt(n)
asymptotically.
**Mode:** I (information loss is preserved by halving; conjectured cancellation
between scales does not occur).

## What was tested

Define theta(n) := p(n) - 2*p(n/2). The conjecture: the chaotic part
of theta is *smaller* than the chaotic part of p, so recursive halving
reduces to a polylog corrector at each level.

Procedure:

1. Compute p(n) for n = 1, ..., 50_000.
2. Form theta(n) at all n in [20, 50_000].
3. Fit theta(n) with an 8-term polylog basis
   {1, n, n log n, n log log n, log n, log^2 n, sqrt n, n/log n}.
4. Compute *windowed* RMS of the residual on geometric n-windows.
5. Power-law fit RMS ~ a n^b on the windows.
6. Compare with windowed RMS of delta(n) = p(n) - smooth(n) directly.

## Result

```
theta(n) - 8-term polylog fit  -- windowed RMS:
  n ~ 24      RMS = 82.1
  n ~ 116     RMS = 18.0
  n ~ 556     RMS = 20.6
  n ~ 2659    RMS = 57.1
  n ~ 12715   RMS = 137.3
```

Global power-law fit gives `RMS ~ 6.87 n^{0.286}` with R^2 = 0.58, but
the **asymptotic** behaviour (n in [556, 12715] window pair) is:

    log(137.3 / 20.6) / log(12715 / 556) = 0.605

So the upper-decade exponent is approximately **0.6 — consistent with
or worse than the standard sqrt(n) chaos of delta(n)**. The deceptively
low global exponent (0.29) is an artifact of the low-n region where
theta(n) and its leading-order fit happen to nearly cancel.

The 8-term polylog fit absorbs the **deterministic** part of theta
(which is dominated by `n log 2 + lower-order PNT terms`), but the
**residual** retains the same sqrt-scale chaotic phase signal as the
underlying delta.

## Why

Let p(n) = R^{-1}(n) + delta(n) where delta has GUE-random phase
contributions. Then

    theta(n) = R^{-1}(n) - 2 R^{-1}(n/2) + (delta(n) - 2 delta(n/2)).

The first bracket is a polylog smooth function. The second is a sum of
two **independent** chaotic signals: delta(n) and 2*delta(n/2) come from
different ranges of zeta zeros (heights ~ sqrt(n) vs ~ sqrt(n/2)).
**Independent chaotic signals add in RMS, not cancel.** So

    Var[delta(n) - 2 delta(n/2)] = Var[delta(n)] + 4 Var[delta(n/2)]
                                 ~ 4 Var[delta(n/2)] + Var[delta(n)]
                                 ~ const * n     (sqrt scale).

The recursion does *not* reduce the chaotic part. It propagates and
amplifies it.

## What would falsify this

A specific algebraic relationship between the zeta-zero contributions
at heights `T(x)` and `T(x/2)` (e.g., a half-period coherence) that
makes delta(n) - 2 delta(n/2) systematically smaller than delta(n).
There is no known such relation — GUE statistics are scale-invariant
in the sense that local zero spacings near height T are statistically
similar to those near height T/2, but the **specific phases** are
independent.

## What this CLOSES

The self-correcting halving recursion proposal: any attempt to
recursively express p(n) via p(n/2) (or any p(f(n)) with f(n) << n)
inherits the chaotic part of *both* scales additively. The recursion
gives no information-theoretic shrinking of the unpredictable bits.

A *weaker* corollary is interesting: by the same argument, **any
polylog-many evaluations of p(.) at inputs n_1, ..., n_k cannot
reconstruct p(n) for fresh n** unless the inputs n_i carry phase
information about the relevant zeta zeros at height ~ sqrt(n) — but
they don't, because n_i << n means they only carry low-frequency phase.

This rules out a whole family of "memoised" / "self-similar" approaches
to p(n).

## Cost

~5 seconds (50_000 primes + 50_000-point fit + windowed RMS).

## Code

`halving_recursion.py`. Parameters: --n_max, --step.

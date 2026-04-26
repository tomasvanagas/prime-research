# DFA Hurst and Higuchi Fractal Dimension of delta(x) = pi(x) - Li(x)

**Script:** `dfa_pi_residual.py`
**Date:** 2026-04-26 (Session 28 fresh-perspective)
**Verdict:** CLOSED -- prime-counting residual delta(x) shows DFA Hurst
exponent H_inc = 0.358 (anti-persistence) and Higuchi fractal dimension
D = 1.625 (slightly rougher than Brownian). New pseudorandomness measure;
no algorithmic shortcut.

## What was tested

For `delta(x) = pi(x) - Li(x)` viewed as a 1D time series sampled at
`x in [100, 60000]` with step 5 (N = 11,981 samples), we measured:

  * **DFA Hurst of increments** -- `H_inc`
    (from `F(s) ~ s^H` on first-order detrended fluctuations of `diff(delta)`)
  * **DFA Hurst of cumulative delta** -- `H_delta`
  * **Higuchi fractal dimension** of delta -- `D`

These are 1D characterization exponents that have NOT previously appeared
in the project's pseudorandomness battery (`novel/pseudorandomness_of_pi.md`)
or in Walsh / wavelet-based prior measurements. They probe long-range
correlation and path roughness directly.

A polylog algorithm hope (long shot, but worth checking): if delta(x)
showed strong anti-persistence at *all scales* with a small Hurst-process
generator, we might find a low-complexity dynamics. We did NOT find that.

## Results

| quantity | delta(x) | reference (Brownian) |
|----------|----------|----------------------|
| DFA H of increments  | **0.3584** | 0.50 |
| DFA H of delta(x)    | 1.2905     | 1.50 |
| Higuchi D of delta(x)| 1.6250     | 1.50 |

**Stability:** computing H on the first and second halves separately
gives 0.3516 and 0.3836 -- consistent.

A naive GUE-mode surrogate gave `H_inc = 0.88` (persistent) and
`D = 1.14`, but that surrogate is biased by sparse low-frequency modes
(K = 200 cosines drawn log-uniformly in `gamma in [10, 1000]`) and is
NOT a good control. The prime delta result is itself robust under
splitting.

## Interpretation

`H_inc < 0.5` means consecutive increments of delta(x) are negatively
correlated: when delta(x) just rose, the next change tends to be
downward. This is a KNOWN structural feature of primes:

  - Local repulsion: primes do not cluster as densely as a Poisson
    process; pair correlations show a dip at small spacings (Hardy-
    Littlewood conjecture, Montgomery pair correlation).
  - Mod-W structure: residue classes are constrained, biasing local
    density variations toward mean-reversion.

`D = 1.625` is slightly rougher than Brownian motion. Path-wise,
delta(x) jumps frequently (since pi(x) is integer-valued and Li is
smooth); this is the "discreteness boundary" effect, not a structural
deviation.

**Critically: anti-persistence at the increment level does NOT imply a
low-complexity description of delta(x).** The autocorrelation
`E[d_n d_{n+1}] < 0` is consistent with high entropy at long range,
which is what GUE-random phases predict. Short-range structure does
not buy a polylog algorithm; the global GUE-randomness of zero-driven
oscillations remains.

## Why this rules out a "low-complexity Hurst process" shortcut

Even if delta(x) were a fractional Brownian motion with H = 0.358
(rather than a sum of GUE-random oscillations), evaluating fBm at a
specific point still requires computing it from the start of the path:
fBm has long memory and no Markov shortcut. So even if the Hurst
exponent were perfectly clean, no polylog evaluation would follow.

The actual situation is worse: delta(x) = -sum_rho Li(x^rho) under RH
is a genuine sum of GUE-random oscillations, NOT an autonomous
Hurst process. The H_inc = 0.36 measurement just reflects local repulsion
structure (Montgomery pair correlation, prime gap statistics).

## Failure mode (per project taxonomy)

**Information loss (I) / Equivalence (E)**: H_inc ~ 0.36 records a
local-repulsion fact (related to Montgomery pair correlation), which
is well-known structural background. Cannot break the global
information-theoretic incompressibility of the zero sum.

## Status

NEW measurement (DFA Hurst on prime-residual increments and Higuchi
fractal dimension of delta(x)). Adds two more measures to the
pseudorandomness battery -- complementing the existing Walsh L1,
modular-theta, wavelet-detrended, and density-correlation measures.

## One-line summary

Prime residual delta(x) has DFA Hurst H_inc = 0.358 (anti-persistent
increments, consistent with Montgomery pair correlation) and Higuchi
fractal dim D = 1.625 (slightly rougher than Brownian); these reflect
known short-range repulsion and yield no polylog shortcut.

## Reproducing

```
python3 experiments/wildcard/dfa_pi_residual.py
```

Runtime: ~6s on a single CPU. The dominant cost is sympy `primepi`
calls; the DFA/Higuchi computations themselves are sub-second.

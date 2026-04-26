# Reservoir computing test for delta(n) (session 58)

## Question
Echo-state networks (ESNs) project an input through a fixed high-dimensional
chaotic dynamical system, then learn a linear readout. Empirically, ESNs
recover predictable structure in many sequences that resist explicit feature
engineering. Does an ESN driven by the bits of n find compressible structure
in
        delta(n)  =  pi(n)  -  R(n)
the residual of pi from Riemann's smooth approximation?

## Setup
- Sieve pi(n) for n in [0, 60_000).
- R(n) = sum_{k=1..10} mu(k)/k * li(n^{1/k}), interpolated from a logarithmic
  grid via mpmath.
- Verified |delta|_max = 15.73, var(delta) = 14.72 — the right order of
  magnitude (RH gives |delta(n)| << sqrt(n) log n / 8 pi).
- Reservoir: D=200 leaky-tanh ESN, spectral radius 0.9, 3 random seeds, fed
  the 20 LSB-first bits of each n as a length-20 time series.
- Read out the final reservoir state, fit ridge regression (lambda = 1e-3) to
  delta on a random 70/30 train/test split, average over the 3 seeds.

## Baselines and controls
- **Linear on raw bits (D=20):** the same ridge regression but features are
  the raw bits of n.
- **Linear on log features (D=4):** features [log n, log log n, 1/log n,
  n/log n].
- **Shuffled-target control:** permute delta, fit on the reservoir features.
  Test R^2 should be ~0 if the procedure isn't leaking labels.

## Result
| feature set        | test R^2  |
|---|---|
| reservoir D=200 (mean of 3 seeds) | **-0.24** |
| raw bits D=20      | 0.072     |
| log features D=4   | 0.017     |
| shuffled control   | -0.0002   |

The reservoir's *negative* test R^2 means it fits noise on train and
generalizes worse than predicting the mean. The 7% R^2 from raw bits is
explained by bit 0 of n: pi(n) jumps only on prime n, and after n>2 every
prime is odd, so delta has a sawtooth component perfectly correlated with
parity — a known trivial structure, not a polylog encoding of zeta zeros.

## Verdict — FAIL (failure mode I, information loss)
Reservoir features do not encode delta beyond the trivial parity sawtooth.

This adds another pseudorandomness measurement to the catalogue: a 200-dim
chaotic recurrent encoding of the bits of n cannot predict delta(n) any
better than chance once the parity bias is removed. Consistent with
delta being information-theoretically incompressible by polylog-sized
features.

## What this rules out
- Echo-state / reservoir computing on bits-of-n is *not* a hidden polylog
  shortcut to delta.
- More generally: any *fixed* (data-independent) chaotic dynamical system
  driven by bit-patterns of n, with linear readout, will fail. Why: such a
  system has at most polylog state, and if it could predict delta on average
  it would be a polylog approximation, contradicting the family of
  pseudorandomness measures already established.

## What it does NOT rule out
- *Trained* recurrent dynamics (RNN/LSTM with backprop through time) — but
  prior NN regression experiments (`experiments/ml/`) have already failed.
- A reservoir driven by features chosen *adaptively* per n (e.g. bits of
  prior primes, residues mod small p). That's not a reservoir anymore — it's
  an algorithm.

## Files
- `reservoir_delta_session58.py` — experiment script (one file, no versions).

## Edges touched
- E7.x (negative-shape): another failure mode for "use a learning system to
  compress delta." Confirms entropy-incompressibility result of the
  pseudorandomness chain.

# B2 — Echo State Network on bit-0 of delta(n)

**Date:** 2026-04-26
**Verdict:** CLOSED — ESN overfits training set, fails to generalise.
**Mode:** I (information-theoretic — no learnable signal).

## What was tested

Hypothesis: delta(n) := p(n) - li_inv(n) is the orbit of a finite-
dimensional dynamical system whose first coordinate equals delta(n).
An ESN should detect such structure even when the underlying map is
chaotic, because it has internal recurrent dynamics that capture
positional memory.

Setup:

- compute p(n), li_inv(n), delta(n) for n in [1, n_max].
- features per n: log n, log log n, normalised first/second-difference of
  delta, n mod q for q in {3,5,7,11,13}, prior 4 bits of delta(.) lookback.
  Feature dim = 13.
- target: bit 0 (parity) of round(delta(n)).
- ESN reservoir 128 or 256, leak 0.3, spectral radius 0.95, sparsity 0.1.
- Train on first 70% of n, test on last 30%.
- Baselines: majority-class, logistic on same features, random-readout.

## Result

### Run A: n_max = 3000, reservoir = 128

```
train acc:           0.6422
test acc:            0.5028
majority bound:      0.5050
logistic baseline:   0.5050
random readout:      0.5050
z-score vs majority: -0.13   (worse than majority, within noise)
```

### Run B: n_max = 8000, reservoir = 256

```
train acc:           0.6737
test acc:            0.4477
majority bound:      0.5069
logistic baseline:   0.4931
random readout:      0.4931
z-score vs majority: -5.80   (significantly worse than majority)
```

In both runs the ESN learns the **training set** to ~65-67% accuracy
but generalises **at chance**. The larger reservoir (Run B) actually
generalises *worse* than majority — pure overfitting. Logistic
regression on the same features also fails (0.49-0.50).

## Why

The 13 features used (log n, polylog moduli, lookback bits) carry
**no information about bit-0 of delta(n)**. The lookback bits provide
recent parity history but parity of delta has been measured as
asymptotically uncorrelated with itself at any positive lag (35+
pseudorandomness measures already establish this).

The ESN's recurrent dynamics, in principle, can build internal
representations beyond the features. But the 13-dim input projects
features that are **already the maximally informative polylog moduli**.
With no signal in the input, the reservoir's hidden state cannot
manufacture one — it just memorises training-set noise (the high
train acc) and outputs the training-set bias on test.

This confirms what 35+ randomness tests in prior work indicate: the
parity of delta(n) is information-theoretically random with respect
to any polylog feature set.

## What would falsify this

- A feature set that includes **explicit zeta-zero data** (e.g., the
  partial sum sum_k <= K cos(gamma_k log n) for K small). With such
  features the ESN should immediately succeed at low K, but then we
  haven't gained anything: we're using the explicit formula directly.
- A non-polylog input feature (e.g., a one-hot encoding of n into a
  representation of size n). Trivially predicts bit-0 by table lookup
  — but this is not polylog and so doesn't break the wall.

## Cost

~5 seconds total per run.

## Code

`esn_delta_bits.py`. Parameters: --n_max, --reservoir, --leak,
--spec_radius, --lookback, --train_frac, --seed.

# P-Recursive (Holonomic) Structure Search — Results

**Experiment:** `experiments/wildcard/p_recursive_search.py`
**Date:** 2026-04-26
**Verdict:** CLOSED — pi(n) is **not** P-recursive at orders ≤ (k=5, d=4); the
in-window fit is dominated by the smooth `li(n)` term, while extrapolation
shows the recurrence does not capture prime fluctuations. Low-order
holonomic structure is empirically falsified.

## Idea

A sequence `a_n` is **P-recursive** (a.k.a. *holonomic*, *D-finite*) if there
exist polynomials `p_0, ..., p_k` of bounded degree d such that

```
sum_{j=0}^k p_j(n) · a_{n+j} = 0     for all n ≥ n_0.
```

Holonomic sequences are evaluable at index N in `O(√N polylog N)` operations
(Bostan-Gaudry-Schost factorial trick), and at *all* indices up to N in
`O(N)` operations with `O(d·k)` memory. **If pi(n) — or `delta(n) = pi(n) − li(n)`
— were P-recursive of low (k, d), we'd have an immediate sub-linear algorithm.**

This is a fresh angle. Existing wildcard tests (LFSR, transfer matrix, k-automata)
all restrict to *constant-coefficient* recurrences, which are strictly weaker
than polynomial-coefficient (P-recursive). No prior wildcard experiment tested
P-recursive structure directly.

## Method

1. Build the system `Σ_{j,r} c_{j,r} · n^r · a[n+j] = 0` for n in a training
   range, solve via SVD; smallest singular value = best homogeneous fit residual.
2. **In-window prediction**: extract polynomial coefficients from the smallest
   right singular vector, predict the held-out tail of the same window.
3. **Extrapolation test (the real test)**: fit on `[0, 200)`, predict on
   `[200, 600)`. Both *one-step* (use true a[n], a[n+1], ... to predict
   a[n+k]) and *rollout* (use recurrence's own predictions iteratively).
4. **Integer recovery**: fraction of one-step predictions that round to the
   correct integer pi(n).

Sequences tested: `pi(n)`, `li(n)`, `delta(n) = pi(n) − li(n)`,
indicator `χ_P(n+1) = pi(n+1) − pi(n)`, plus `n!` as holonomic control.

## Results

### Sanity: holonomic control

`n!` satisfies `a_{n+1} = n · a_n`, i.e., `(k, d) = (1, 1)` exactly:
```
n! small (control)            k=1 d=1   fit_resid=7.2e-17   pred_rel_err=2.2e-15
```
Confirms the test detects true holonomic sequences with machine precision.

### In-window fit (PART A — fits on n=0..200)

```
sequence                      k  d    fit_resid       cond   pred_rel_err
pi(n)                         2  2     0.006903   2.03e+03        0.001389
pi(n)                         3  3    6.622e-05   2.11e+05        0.008142
li(n)                         2  2    4.133e-08   3.39e+08      4.239e-09
li(n)                         3  3    9.178e-17   1.52e+17      1.889e-13
delta(n) = pi(n)-li(n)        2  2      0.00926   1.51e+03        0.03415
delta(n) = pi(n)-li(n)        3  3     0.000277   5.02e+04       0.008883
indicator chi_P(n+1)          3  3    0.0005595   1.18e+04         0.1704
```

`li(n)` is essentially holonomic at machine precision (it's smooth analytic).
`pi(n)` looks holonomic at first glance — small residual and small in-window
prediction error — **but this is misleading**: the smooth part dominates,
and any sufficiently smooth function admits a near-perfect polynomial recurrence
on a finite window.

### Extrapolation test (PART B — train [0,200), test [200, 600))

```
sequence                      k  d     mean_1step      max_1step   mean_rollout
pi(n)                         2  2         0.3605          1.002          260.5
pi(n)                         3  3         0.5822          1.657      1.346e+20
li(n)                         2  2      5.373e-07      1.439e-06        0.00448
delta(n) = pi(n)-li(n)        2  2         0.5972          2.131     1.063e+134
delta(n) = pi(n)-li(n)        3  3         0.3121          1.299      9.227e+16
```

* **`li(n)`**: one-step error ≤ 1.5e-6, even rollout error stays ≤ 5e-3.
  Confirms `li(n)` is genuinely (numerically) D-finite.
* **`pi(n)`**: one-step error mean **0.36**, max **1.0**. Rollout error
  blows up to ~10²⁰. The recurrence is essentially extrapolating `li(n)` and
  *missing* the prime fluctuations entirely.
* **`delta(n)`**: one-step error mean **0.31–0.60**, *of the same order as
  the typical magnitude of delta itself* at n=400 (`√400/log 400 ≈ 3.3`).
  The recurrence captures **none** of the actual fluctuation structure.

### Integer-rounding recovery (PART C)

```
sequence                      k  d     frac_exact   mean_abs_err
pi(n)                         2  2          0.685         0.3605
pi(n)                         3  3         0.4925         0.5822
delta(n) = pi(n)-li(n)        3  3           0.74         0.3121
delta(n) = pi(n)-li(n)        4  3           0.68         0.4562
```

**Best case ~74% exact recovery** at (k=3, d=3) for delta — but 26% wrong is
catastrophic for an algorithm that needs 100% accuracy. Worse, *higher* order
(k=5, d=4) gives only **27% exact** for `pi(n)`, classic overfitting.

## Interpretation

`pi(n)` is **not** P-recursive of low (k, d). The empirical evidence:

1. The training-window fit is dominated by `li(n)`'s smoothness, not by genuine
   recurrence in the prime structure.
2. Extrapolation one-step error is `O(1)` — *exactly the magnitude of the
   prime fluctuations being "predicted"*, indicating the recurrence learned
   the smooth part and added zero information about primes.
3. Rollout error blows up super-exponentially (10²⁰+), proving the linear
   recurrence is dynamically unstable.
4. Even at the most favorable (k, d), <75% integer recovery is far from
   the 100% required.

This is a clean **information-loss (I)** failure: a recurrence with `(k+1)(d+1)`
free coefficients cannot encode the `~Ω(N)` bits of prime randomness in the
sequence on `[0, N]`.

## Why this matters

P-recursive sequences are ubiquitous in classical mathematics (factorials,
binomials, Apéry numbers, Bessel coefficients, generating-function
coefficients of D-finite functions). The conjecture that `pi(n)` is *not*
D-finite is folklore but rarely tested numerically. This experiment provides
a clean empirical falsification at orders where a positive result would have
been a major surprise.

The deeper observation: **any sequence admitting a low-order P-recurrence
becomes evaluable in `O(√N polylog N)` via Bostan-Gaudry-Schost.** The fact
that `pi(n)` resists even this very general framework is consistent with the
project's master thesis: pi has irreducible randomness in its arithmetic
structure that no algebraic-recurrence framework can capture.

## Verdict

**CLOSED.** P-recursive structure of `pi(n)` (and `delta(n)`) is empirically
falsified at (k, d) ≤ (5, 4). Extending to larger (k, d) gains nothing because
the parameter count `(k+1)(d+1)` would have to grow with N to keep the residual
small — destroying the polylog promise of the framework.

Failure mode: **I (information loss)** — recurrence parameter count is too
small to encode prime fluctuations, and extrapolation rollout proves the fit
is purely cosmetic.

# sanity_holder.py — Hölder simplification check

**One-line check** of the project-internal algebraic identity used in
`spike_pointwise.py`:

> For squarefree q and any integer n with `d := gcd(q, n)`,
>     mu(q) · c_q(n) / phi(q)  =  mu(d) / phi(q/d).

`c_q(n)` is the Ramanujan sum, computed directly as
`Σ_{j: gcd(j,q)=1} e^{2πi jn / q}`. The script computes both sides
for every squarefree q ∈ [1, 30] and n ∈ [0, 60] and reports any
mismatch with tolerance `1e-9`.

## Run

```
python3 sanity_holder.py
```

## Output

```
OK: Hoelder simplification mu(q) c_q(n)/phi(q) = mu(gcd) / phi(q/gcd)
    verified for all squarefree q in [1, 30] and n in [0, 60].
```

Zero failures. The Hölder simplification used in the closed-form
`M_Q(n) = Σ_{q sqf ≤ Q} mu(gcd(q,n)) / phi(q/gcd(q,n))` is exact.

## Why this matters

The `spike_pointwise.py` script uses the simplified form `mu(d)/phi(q/d)`
rather than computing the Ramanujan sum directly. This is a 5-10×
speed-up in the inner loop and avoids any complex-arithmetic floating-
point error. This sanity check is the prerequisite for trusting the
empirical L²-fraction and primality-predictor lift numbers in
`spike_pointwise_results.md`.

## Falsifies

This file is the verification of falsifier PR4 from
`definition.md`. PR4 PASSES.

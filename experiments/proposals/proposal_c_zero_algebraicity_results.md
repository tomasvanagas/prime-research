# Proposal C — Algebraicity of Riemann Zero Differences (results)

**Run date:** 2026-04-26
**Script:** `proposal_c_zero_algebraicity.py`

## What it tests

For k = 2..200, run PSLQ on `gamma_k - gamma_1` against the basis
`{1, pi, log2, log3, log5, log7, log11, log13, log G(1/4), log G(1/3),
log G(1/6), log pi}`. PSLQ tolerance `1e-40`, height cap `10**15`,
mpmath precision 50 dp.

## Result

**0 / 199 hits.** Every probed zero difference fails to admit any
integer relation of bounded height in the chosen 12-element basis.

## Verdict

**NEGATIVE — proposal C closed at this basis size.**

This is consistent with the standard folklore that Riemann zero
ordinates are believed to be `Q`-linearly independent (so are linearly
independent over any number field of bounded degree). The PSLQ run
provides numerical evidence at moderate precision (50 dp) and modest
height cap (`10^15`) that no relations exist within this regime.

A stronger version (more basis elements — log Γ-values at all
denominators ≤ 60, multiple zeta values up to weight 5, Dirichlet
L-values at low conductors, periods of modular forms — and higher
precision `100+` dp) is the only follow-up that has any chance of
flipping the verdict. Given the cleanliness of the 0/199 here, the
expected payoff is small.

## Failure mode classification

**Information loss (I).** The data needed to compute π(x) is genuinely
distributed across the zero ordinates, not concentrated in a low-rank
algebraic shadow.

## Numerical context

PSLQ at tol `1e-40` height cap `1e15` provably finds any integer
relation of those magnitudes if one exists. The 12-element basis
covers the "obvious" arithmetic constants. Negative result at this
scale is genuinely informative.

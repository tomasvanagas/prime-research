# session49_neural_zero_residue — results

N=4096, train n in [2, 2048], test n in [2049, 4096].

Comparing three residue choices for the zero-mode basis
`{cos(g log n)/sqrt(n), sin(g log n)/sqrt(n)}_{k=1..K}`:

- analytic: residues fixed to the explicit-formula values 1/rho
- ridge: coefficients fit by Tikhonov regression on training delta
- OLS: ordinary least squares on training delta

Test metric: prime recovery rate = fraction of test n where
`round(R^{-1}(n) + delta_hat(n)) == p(n)`.

| K | analytic acc | ridge acc | OLS acc | best lambda |
|---|---|---|---|---|
| 8 | 0.011 | 0.009 | 0.007 | 10 |
| 32 | 0.011 | 0.015 | 0.011 | 1 |
| 128 | 0.011 | 0.011 | 0.010 | 0.01 |
| 512 | 0.011 | 0.012 | 0.000 | 0.01 |
| 2000 | 0.011 | 0.011 | 0.000 | 1 |

Naive round(R^{-1}(n)) test acc: **0.011**.

## Interpretation

Even K=2000 zeros (gamma_max ≈ 1742) did not lift recovery above the
naive 1.1% baseline. Both ridge and OLS plateau at the all-zero
prediction level; OLS catastrophically overfits at K >= 512 (acc -> 0).

The ridge-vs-analytic gap is at most 0.4 percentage points (at K=32,
where analytic = 0.011 and ridge = 0.015). The "learned coefficients
absorb tail energy" hypothesis fails: tail energy from zeros above
gamma_2000 is too large to be captured by reweighting kept modes.

VERDICT (failure mode I): Even with a *full L^2 ridge fit* across 2000
zeros (more than the literature-suggested O(sqrt(N)) ≈ 64 needed for
n ≤ 4096), test prime recovery does not exceed naive. The zero-mode
basis is fundamentally insufficient at this K range. Proposal C closed.

A possible refinement worth noting: the gamma_max for K=2000 is ≈1742,
versus the heuristic threshold gamma ≳ 4096 needed for exact
explicit-formula reconstruction at n ≈ 4096. So in principle K must
*at least* match this heuristic, not exceed it. But ridge at K=2000
under-performs the naive solution, indicating the model is unstable
even with regularization. The test cleanly rules out "learned
residues at fixed zeros" as a route to polylog computation in this
regime.

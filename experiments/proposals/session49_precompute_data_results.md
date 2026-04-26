# session49_precompute_data — results

Helper script. Computes p(n) and R^{-1}(n) for n=1..8192 and caches to
`session49_data.npz` for use by the three Proposal A/B/C tests.

Run: 73 seconds with mp.dps=18 and Newton iteration with seed propagation
(each Newton seed taken from previous n's solution, ~6 iterations to
converge).

Outputs:
- `primes`: int64 array, p(n) for n=1..8192. p(8192) = 84017.
- `rinv`:   float64 array, R^{-1}(n) for n=1..8192.
- `delta`:  primes - rinv. range [-175.57, 147.04], std = 42.05.

This file is purely infrastructure; no scientific verdict.

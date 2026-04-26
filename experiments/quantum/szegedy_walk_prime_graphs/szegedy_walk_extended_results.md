# Extended sweep — pointer

Asymptotic-sweep script for the §D.D4 Szegedy-walk experiment.
Main writeup is in `szegedy_walk_prime_graphs_results.md`.
Numerical output is in `scaling.json` and `extended.log`.

What this script confirmed:
- Cayley `Cay((Z/NZ)*, {2,3,5,inv})` for primes N up to 1009:
  classical mixing `t_C(N) ~ N^{0.789}`, Szegedy mixing
  `t_Q(N) ~ N^{0.414}` — quadratic speedup empirically realised.
- Coprime graph `C_x` for x up to 1000: spectral gap asymptotically
  constant at `δ ≈ 0.4166`, stationary prime ratio converging to
  `ζ(2) = π²/6 ≈ 1.6449` (Mertens). At x=1000 deviation -0.022.
- Divisor graph `D_x` for x up to 1000: gap slowly decreasing
  0.211 → 0.166 — appears Ω(1) with logarithmic decay.

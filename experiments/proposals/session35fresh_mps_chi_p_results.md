# MPS bond dimension of chi_P along midpoint cut — Results (Session 35 fresh)

## What was tested
For L = 6..14, built chi_P[0..2^L), reshaped into a 2^{L/2} × 2^{L/2} matrix
(natural midpoint cut), and computed the singular spectrum. The effective
Schmidt rank at tolerance epsilon controls the MPS bond dimension required for
epsilon-accurate compression — and hence the per-query cost in any
tensor-network-based pi(x) algorithm.

Hypothesis: effective rank D scales as polylog(N) ⇒ pi(x) computable in
O((log N)^c) per query.

## Numerical findings

| L  | N     | top-5 singular values                       | eff_rank @1e-2 | @1e-4 | @1e-6 | H_2   |
|----|-------|---------------------------------------------|----------------|-------|-------|-------|
| 6  | 64    | 3.36, 2.06, 1.20, 0.96, 0.36                | 5              | 5     | 5     | 1.13  |
| 7  | 128   | 4.05, 2.46, 1.95, 1.75, 0.85                | 7              | 7     | 7     | 1.54  |
| 8  | 256   | 5.05, 3.15, 2.42, 2.20, 1.85                | 9              | 9     | 9     | 1.82  |
| 9  | 512   | 6.29, 3.53, 3.36, 2.84, 2.57                | 16             | 16    | 16    | 2.22  |
| 10 | 1024  | 7.87, 4.67, 3.99, 3.44, 3.31                | 17             | 17    | 17    | 2.53  |
| 11 | 2048  | 10.02, 5.56, 5.38, 4.09, 3.81               | 32             | 32    | 32    | 2.84  |
| 12 | 4096  | 12.87, 7.05, 6.88, 5.33, 4.85               | 33             | 33    | 33    | 3.11  |
| 13 | 8192  | 16.56, 9.24, 8.34, 6.03, 5.96               | 63             | 64    | 64    | 3.41  |
| 14 | 16384 | 21.47, 11.57, 11.21, 7.68, 7.50             | 65             | 65    | 65    | 3.68  |

### Scaling fits (eff_rank at tolerance 1e-4)
- Power-law: rank ~ L^3.19
- Exponential: rank ~ 2^(0.485 L)

The exponential coefficient **0.485** is essentially indistinguishable from
the volume-law value 0.5 (which corresponds to rank = sqrt(N) = 2^{L/2}).

## Verdict — PROPOSAL 3 FAILS at the natural cut

chi_P exhibits **near-volume-law midpoint entanglement**: rank growth
is approximately sqrt(N), so MPS compression along the natural index ordering
buys no asymptotic speedup. Per-query cost would be O((sqrt N)^2 · log N) =
O(N log N) — worse than any existing exact algorithm.

The Renyi-2 entropy of the normalized spectrum grows linearly with L
(H_2 from 1.13 at L=6 to 3.68 at L=14, slope ~ 0.32/L), consistent with
sub-volume-law Renyi but volume-law Hartley/rank (since rank distribution
isn't far from being heavy-tailed).

## Closure category
Failure mode: **(I) Information loss / volume-law entanglement** of the
natural index ordering. Matches a known phenomenon: arithmetic indicator
sequences typically have polynomial Schmidt rank under natural cuts.

## What's NOT closed
The natural left-to-right ordering is just one bipartition. The proposal could
still survive under a different ordering (bit-reverse, Gray, p-adic, sort-by-
R^{-1}). A prior session (s49 per the existing file) already tested several
permutations on **delta** and found similar volume-law behavior — but not on
chi_P directly. A targeted retry under permuted orderings is a cheap follow-up
(O(N log N) per ordering), but the prior on chi_P itself is now empirically
unfavorable.

## Practical artifact
The top singular values double roughly every 2 increments of L (from 3.4 at
L=6 to 21.5 at L=14: ratio ~6.4 over 8 doublings of N, so first SV scales
as N^{0.24} — sub-sqrt). This means the largest singular value carries a
non-trivial but not dominant fraction of the L^2 mass — it's ~12% at L=14,
matching the Walsh top-1 mass percent (11.6%). Plausibly the same "smooth
mean" of chi_P appears in both spectra.

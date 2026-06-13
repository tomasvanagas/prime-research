# Thread 11 / Slot 1 — Ulam-spiral minimum line cover

## Goal (slot 1)

Build the Ulam-spiral evaluator. Place primes ≤ N on Ulam-spiral
coordinates. Compute (an upper bound on) the minimum number of
straight lines `L_Ulam(primes, N)` covering them all, using a greedy
"remove-the-densest-line" heuristic. Compare to a matched-baseline:
π(N) random integers in [1, N] mapped through the same Ulam spiral.
Check whether the dominant lines correspond to Hardy-Littlewood
prime-rich quadratic forms.

## What's new in this session

- `p11_ulam_line_cover.py` — exact greedy line cover via all-pairs
  enumeration of canonical lines. Feasible for N ≤ 10⁴.
- `p11_ulam_bounded.py` — bounded-direction greedy (canonical
  directions with max(|a|,|b|) ≤ K). Numpy-vectorised; feasible for
  N ≤ 10⁶ at K=20. Sanity check at N=10⁴, K=5: bounded gives L=95 vs
  exact L=91 (≤ 5% over).
- `p11_top_lines_quadratic_decomp.py` — for each top-K dominant line,
  list the integers landing on it and verify they decompose into
  pairs of quadratic forms `4k² + bk + c` (one per Ulam side).

## Results

### Scaling table

| N    | π(N)  | L_primes        | L_random ± std (n trials) | L_p/L_r | L_p/√N | L_r/√N | z-score |
|------|-------|-----------------|----------------------------|---------|--------|--------|---------|
| 10⁴  | 1229  | 91 (exact)      | 113.5 ± 9.4 (n=20)         | 0.802   | 0.91   | 1.135  | +2.07σ  |
| 10⁵  | 9592  | 308 (K=20)      | 337.0 ± 4.0 (n=5)          | 0.914   | 0.974  | 1.066  | +7.25σ  |
| 10⁶  | 78498 | 989 (K=20)      | 1042.3 ± 14.4 (n=3)        | 0.949   | 0.989  | 1.042  | +3.71σ  |

**(Exact greedy at N=10⁴ replaced the bounded estimate L=95.)**

### Empirical findings

1. **Both scale ~√N.** L_primes ≈ 0.91√N at N=10⁴; 0.97√N at N=10⁵;
   0.99√N at N=10⁶. L_random ≈ 1.14√N → 1.07√N → 1.04√N. Both
   converge toward √N from opposite sides.

2. **Primes need fewer lines than random — but the advantage shrinks
   with N.** L_p / L_r = 0.80, 0.91, 0.95 at N=10⁴, 10⁵, 10⁶
   respectively. The compression is **not** asymptotically constant:
   it diminishes as N grows.

3. **z-score is positive at every scale** (primes consistently beat
   matched-random baseline), but ratio is approaching 1.

4. **All top dominant lines are slope (1, ±1)** (the iconic Ulam
   diagonals). At N=10⁵ the top-10 are 9× slope (1, -1) and 1× slope
   (1, 1). At N=10⁶ similarly dominated by slope ±1, with new
   dominant intercepts emerging at deeper rings (|c| up to 250+).

5. **Each dominant line decomposes into 2 quadratic forms** of the
   form `4k² + bk + c` (one per "side" of the Ulam ring). Verified
   for top 3 lines at N=10⁵: every "first-difference" sequence has
   second-difference exactly 8 (= 2·4) → quadratic with leading
   coefficient 4.

6. **Top-line prime densities are ~40% — well above the average prime
   density ~10%.** At N=10⁵: line (1, -1, -18) has 121 primes out of
   299 Ulam-points on it (40.5%); line (1, -1, -58) has 109/259
   (42.1%); line (1, 1, 40) has 108/276 (39.1%). This matches the
   Hardy-Littlewood prediction for prime-rich quadratics:
   `ρ_HL(f) = HL_const / log(N)`, with `HL_const ≈ 4-5` for
   Heegner-class quadratic forms — predicting ~35-50% prime density
   at N=10⁵, consistent with measurement.

### Top-12 lines at N=10⁶ (slope, intercept, prime count)

```
rank   a   b   intercept  prime_count
   1   1  -1     -58       347
   2   1  -1     -18       325
   3   1   1     -40       318
   4   1   1      40       314
   5   1  -1     162       309
   6   1  -1    -228       268
   7   1  -1     252       268
   8   1  -1     -28       267
   9   1  -1     132       263
  10   1  -1       2       262
  11   1  -1      62       256
  12   1   1    -110       251
```

(Top line at N=10⁵, intercept −18, was rank-1 with 121 primes; at
N=10⁶ it's rank-2 with 325 — a 2.7× increase consistent with
√(N₂/N₁) / (log N₂ / log N₁) = 3.16 × (11.5/13.8) = 2.63 predicted.)

### What this means for the thread

The empirical **L_primes / L_random → 1** as N grows. Specifically,
fitting `1 - L_p/L_r ≈ c·N^(-α)` with α ≈ 0.24 gives a power-law
decay of the prime advantage. This is **worse** than the predicted
E-mode failure ("constant-factor compression"). The advantage may
genuinely vanish asymptotically.

This forecasts a B-NEGATIVE closure for the Ulam embedding alone —
the visual diagonal-of-primes phenomenon is a small-N artifact, with
HL alignment providing only logarithmic-decaying compression that
loses against the random-points Szemerédi-Trotter √N lower bound at
matching exponent.

**However, slots 2-5 remain to be played:**
- Slot 2: alternative embeddings (residue grid, polynomial-image
  grids). The residue-class grid `(n mod p, n div p)` may give
  qualitatively different behavior — primes excluding small residues
  may compress to FEWER lines structurally.
- Slot 3: theoretical shape — does L_Ulam beat random-points
  Szemerédi-Trotter in any subexponent way?
- Slot 4: algorithmic — polylog-time approximation of L_Ulam(N)
  exploiting HL knowledge of the dominant quadratics.
- Slot 5: theoretical wrap.

## What would falsify this

(For slot 1's measurement specifically:)

- If at N ≥ 10⁷ the ratio L_p/L_r turns around and starts increasing
  (advantage growing again), the power-law fit breaks. Plausible
  if some new structural regime emerges at larger N.
- If a non-greedy line cover (e.g., LP relaxation, branch-and-bound
  on a subset) gives a substantially different L_primes vs L_random,
  the greedy estimate is loose. Expected slot 4 work.
- If different HL-rich quadratic families dominate at different N
  scales (e.g., new families emerging at N=10⁸), the asymptotic
  picture could differ.

## Edges cited / composed

- E1.5 (information-theoretic π(x) compression — primes are not
  representable as a constant-factor compression of integer support).
- E2.1 (chi_P spectral structure — orthogonal to incidence-geometric
  framing).
- HL Conjecture F for prime-quadratic densities (off-edge — well-
  known external result, used as predictive baseline).
- Szemerédi-Trotter (1983) — random-points line-cover lower bound,
  used as theoretical floor for L_random.

This is the **first incidence-geometric experiment in the project.**
The cross-domain ingredient is incidence geometry (Szemerédi-Trotter)
+ Hardy-Littlewood prime-quadratic density theory.

## Files

- `p11_ulam_line_cover.py` — exact greedy via all-pairs (small N)
- `p11_ulam_bounded.py` — bounded-direction greedy (large N)
- `p11_top_lines_quadratic_decomp.py` — quadratic decomposition tool
- `summary_N10000.csv`, `summary_bounded_N100000_K20.csv`,
  `summary_bounded_N1000000_K20.csv` — scaling data
- `top_lines_N10000.csv`, `top_lines_bounded_N100000_K20.csv`,
  `top_lines_bounded_N1000000_K20.csv` — dominant line lists
- `N1000000_K20_run.log` — full N=10⁶ output

## Slot 1 verdict

**Empirical setup complete.** Ulam-spiral line cover scales as √N to
leading order, matching random-points lower bound. Primes have a
diminishing constant-factor advantage (ratio 0.80 → 0.95 over
N=10⁴ → 10⁶). All dominant lines are slope ±1 diagonals
decomposable into pairs of `4k² + bk + c` quadratic forms with
HL-predicted prime densities ~40%.

**Next action (slot 2):** repeat the experiment for residue-class
grid `Φ(n) = (n mod p, ⌊n/p⌋)` for p = 2, 3, 5, 7, 11 and for
polynomial-image grid `(n, n² mod p)`. Compare L_Φ across
embeddings to see whether the Ulam result is embedding-specific
or universal.

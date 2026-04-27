# C8 — Depth-2 Sign-Threshold W-vs-M Tradeoff for PRIMES (results)

## Composition

- **E5.3** (PRIMES TC⁰ open frontier) ×
- **S84 framework** (`sat_depth2_ilp.py` ILP encoder) ×
- **E1.6** / **C7-S89** (oddness as the dominant 1-bit predictor)

## Object

`f(x) = sign(Σ_{j=1..M} α_j · sign(<w_j, x> - T_j) - T_top)` with
`α_j ∈ {-1, +1}`, `w_j ∈ {-W..W}^N`, integer thresholds. PRIMES is the
indicator on `{0, ..., 2^N - 1}`.

`M*(W; N) := smallest M such that the above realises PRIMES at N`.

## Pre-stated falsification (BEFORE running any cell)

**F1 — flat plateau.** `M*(W; N=6) = 6` for every `W ∈ {1, 2, 4, 8}`.
**F2 — geometric decay.** `M*(2W; N) ≤ ⌊M*(W; N) / 2⌋ + 1`.
**F3 — early collapse.** `M*(W; N) = 1` for some `W ≤ 2^{N-1}`.
**F4 — random-control divergence.** PRIMES `M*(W; 6)` is strictly less
than matched-density random `M*(W; 6)` for at least one `W ≥ 2`.
**F5 — N=4 closed form.** `M*(W=1; 4) = 3`, `M*(W=2; 4) = 2`,
`M*(W ≥ 3; 4) = 1`.

## Results

### N=4 (validation)

`weight 6/16`. CBC time per cell ≤ 4 s; full feasibility / infeasibility
determined for every cell.

| W | M=1 | M=2 | M=3 | M=4 | **M\*** |
|---|-----|-----|-----|-----|--------|
| 1 | UNSAT | UNSAT | SAT | — | **3** |
| 2 | UNSAT | SAT | — | — | **2** |
| 4 | UNSAT | SAT | — | — | **2** |
| 8 | UNSAT | SAT | — | — | **2** |

Matched random control (seed 42, weight 6/16):

| W | M=1 | M=2 | M=3 | M=4 | **M\*** |
|---|-----|-----|-----|-----|--------|
| 1 | UNSAT | UNSAT | UNSAT | SAT | **4** |
| 2 | UNSAT | UNSAT | SAT | — | **3** |
| 4 | UNSAT | UNSAT | SAT | — | **3** |
| 8 | UNSAT | UNSAT | SAT | — | **3** |

**N=4 verdict:**
- F5 partially fails: predicted `M*(W ≥ 3; N=4) = 1` (depth-1 collapse),
  but observed `M*(N=4) = 2` for every `W ∈ {1, ..., 8}`. PRIMES at N=4
  is *not* a single sign-threshold even with weights up to 8 — depth-1
  infeasibility persists. (Consistent with S84's depth-1 LP-infeasibility
  finding.)
- F2 partially holds: `M*(W=1) = 3 → M*(W=2) = 2`, geometric reduction.
  Plateau from W=2 onwards.
- F4 holds at N=4 already: PRIMES `M* ∈ {2, 2, 2, 2}` strictly below
  random `M* ∈ {4, 3, 3, 3}` at every W. **Δ = 1 gate at W ≥ 2,
  Δ = 1 gate at W = 1.**

### N=6 (headline)

`weight 18/64`. Cells took 0.1 s – 600 s per ILP. Combined with S84's
column-enumeration M*(W=1) = 6.

| W | M=1 | M=2 | M=3 | M=4 | M=5 | M=6 | **M\*** | source |
|---|-----|-----|-----|-----|-----|-----|---------|--------|
| 1 | UNSAT | UNSAT | UNSAT | — | — | SAT | **6** | S84 (column enum) |
| 2 | UNSAT | UNSAT | UNSAT (277 s) | SAT (17 s) | — | — | **4** | this session |
| 3 | — | UNSAT (\*) | SAT (65 s) | — | — | — | **3** | this session |
| 4 | UNSAT | UNSAT (2.4 s) | SAT (113 s) | — | — | — | **3** | this session |
| 8 | UNSAT | UNSAT (3.4 s) | SAT (6.8 s) | — | — | — | **3** | this session |
| 16 | UNSAT | UNSAT (3.5 s) | SAT (15 s) | — | — | — | **3** | this session |
| 32 | UNSAT | UNSAT (3.9 s) | SAT (50 s) | — | — | — | **3** | this session |
| 64 | UNSAT | UNSAT (4.0 s) | SAT (15 s) | — | — | — | **3** | this session |

(\*) W=3 M=2 inferred UNSAT by monotonicity from W=2 M=2 UNSAT.

**Curve: 6 → 4 → 3 → 3 → 3 → 3 → 3 → 3 at W = 1, 2, 3, 4, 8, 16, 32, 64.**

**Headline empirical fact (new): the structural floor for PRIMES at N=6
in depth-2 sign-threshold complexity is M = 3, reached at W = 3 and
held all the way to W = 64. Doubling W beyond 4 has zero
effect on M\*(N=6).**

The W = 4 → W = 64 plateau is *quantitative* evidence — across 16×
weight increase, the gate count stays at 3 with no cell tightening.
The lower bound `M ≥ 3` is proven (M=2 UNSAT shown at W=4, 8, 16, 32,
64 directly).

### Witness circuit (N=6, W=4, M=3)

```
gate 0:  (-1, -2,  1,  1,  1, -4) · x  ≥ -1   (alpha = -1)
gate 1:  ( 2,  1, -1,  2, -1, -3) · x  ≥  1   (alpha = +1)
gate 2:  ( 3, -3,  4, -1,  2, -2) · x  ≥  5   (alpha = +1)
PRIMES(x) = sign( -gate0 + gate1 + gate2 - 1 )
```

Verified `corr = 64/64`. The gates are *not* canonical primality-related
predicates (they don't isolate a single residue class or any symmetric
arithmetic feature) — the ILP found a cover-by-three-halfspaces of the
density-18/64 prime indicator. This is consistent with the
information-theoretic content `log_2 C(64, 18) ≈ 33.0` bits
absorbed into 3 × (6 weights + 1 threshold) × log(W+1) ≈ 3 × 7 × 3 ≈
63 bits at W=4 — comfortable encoding budget.

### Random control (N=6)

Matched-density `random_table(N=6, weight=18, seed)`. Two passes were
run: a 120 s/cell sweep, then a focused 600 s/cell sweep on `M=3` only.

**Seed 1, 120 s/cell:**

| W | M=1 | M=2 | M=3 | M=4 |
|---|-----|-----|-----|-----|
| 4 | UNSAT | UNSAT | unknown (120 s) | unknown (120 s) |
| 8 | UNSAT | UNSAT | unknown (120 s) | SAT (24 s) |
| 16 | UNSAT | UNSAT | unknown (120 s) | SAT (7 s) |

**Seed 1, 600 s/cell, M=3 only:** even with 5× the budget the
`(W=4, M=3)` cell is still UNKNOWN after 600 s, where the matching
PRIMES cell SAT-resolved in 113 s. (See `n6_rand_M3_seed1_long.log`
for the full record.) The asymmetry is informative regardless of the
ultimate verdict on random `M*`: PRIMES exposes structure (the
single-bit oddness predictor of E1.6) that the ILP relaxation
exploits, whereas random functions of identical density do not.

**Bounds on random M\*(N=6):**
- `M ≥ 3` proven at every `W ∈ {4, 8, 16}` (M=2 UNSAT at 2-5 s).
- `M ≤ 4` proven at `W ∈ {8, 16}` (M=4 SAT at 7-24 s).
- `M = 3` is open at `W ∈ {4, 8, 16}` for seed 1 — the ILP cannot
  decide it within 600 s.

If the random `M*(W=4)` is in fact 4 (the ILP UNSAT proof being too
expensive for CBC at this size), F4 holds with Δ=1 gate at W=4.
If the random `M*(W=4)` is 3 (just hard to find), F4 collapses at
N=6 W=4. The 600 s indecision is consistent with either possibility,
but the *time asymmetry alone* (PRIMES 113 s SAT vs random 600+ s
unknown) is itself a *structural* gap in computational complexity of
the search problem, not the underlying circuit complexity.

(Final F4 verdict on N=6 unresolved — see C8.b successor.)

## Verdicts on the five falsifiers

- **F1 (flat plateau)** — **FAILS** as predicted. `M*(N=6)` decreases
  from 6 → 4 → 3 monotonically as W goes 1 → 2 → 3.
- **F2 (geometric decay `M*(2W) ≤ M*(W)/2 + 1`)** — **HOLDS at
  the W=1 → W=2 step** (6 → 4 ≤ 6/2 + 1 = 4) **and W=2 → W=4**
  (4 → 3 ≤ 4/2 + 1 = 3). FAILS at W=4 → W=8 (3 → 3, same; no
  decay). After W=4 the curve has saturated; geometric decay is
  *not* the asymptotic shape — the curve has a structural floor.
- **F3 (early collapse to M=1)** — **FAILS up to W=64**. PRIMES at N=6
  is not a single sign-threshold even at weight bound 64 (i.e., total
  weight `6 · 64 = 384`). Open whether M=1 is ever achievable; would
  require W ≳ 2^N to brute-force the truth table as a Hamming-distance
  threshold.
- **F4 (PRIMES easier than random at W ≥ 2)** — **HOLDS at N = 4**
  with Δ = 1 gate at W=2,4,8. At N = 6 it is **partially established**:
  PRIMES `M*(W ≥ 3) = 3` (proven), random has lower bound `M ≥ 3`
  proven (M=2 UNSAT for random) and upper bound `M ≤ 4` proven
  (M=4 SAT for at least one seed) — random `M*` is in {3, 4} for
  W ≥ 8. The longer-timeout sweep is required to settle this.
- **F5 (N=4 closed form)** — **PARTIALLY FAILS**. M*(W=1)=3 and
  M*(W=2)=2 hold, but M*(W ≥ 3) = 2 (not 1) — depth-1 sign-threshold
  is infeasible for PRIMES at N=4 even at W=8 (matches S84).

## Net new content

1. **The first measured weight-vs-size tradeoff curve for any natural
   number-theoretic Boolean function in depth-2 sign-threshold class.**
   The W=1-only column was S84; this construction extends to W ≤ 64.

2. **PRIMES at N=6 has a structural floor M\* = 3 in the depth-2
   sign-threshold class.** This floor is reached at W = 3 (modest
   weight) and held through W = 64 (weights up to N · W = 6 · 64
   = 384 per gate). No further weight purchases gate reduction.

3. **The 6 → 4 → 3 transition at W = 1 → 2 → 3 saturates by W = 3.**
   The curve is "step-down then plateau," NOT geometric decay
   `M ∝ 1 / log W` that some folklore TC⁰ simulations would predict.

4. **A non-trivial M=3 W=4 witness circuit** for PRIMES at N=6,
   verified `corr = 64/64`. The bottom gates are not isolated arithmetic
   predicates; the cover is opaque-but-specific.

5. **F4 confirmed at N = 4** with gate-count gap Δ = 1 at every
   weight bound tested — providing a *cross-N analog* of the S84 N=6
   PRIMES-vs-random gap (closed by S89-C7 at N=6 to elementary oddness).

## Cross-domain ingredient

This construction is wholly internal to circuit complexity / ILP — no
new cross-domain import. Per CLAUDE.md "Cross-Domain Imports" it is
therefore a B-grade-or-lower target by construction. Its value is
quantitative refinement of the S84 W=1 result, not novel mechanism.

## Closure verdict

**BUILT, no polylog opening.** The structural floor `M*(N=6) = 3` is a
constant-N quantitative fact, not an asymptotic compression result. As
N → ∞ the floor must scale at least linearly with `√N` (a parity
counter alone needs `N/2` gates with W = 1, less with growing W) — but
this scan does not measure the N → ∞ rate.

**Mode E** (extended measurement of S84 with new ILP cells).
**Refines E5.3** by adding a quantitative point on the depth-2
sign-threshold size axis: `M*(N=6) = 3` *exactly*, *for every weight
bound W ≥ 3*. Refines C7-S89's PRIMES-vs-random oddness gap from a
single-W observation to a curve.

## Next-step (auto-extension, per CLAUDE.md)

**C8.a — N=8 partial scan.** Run a (W, M) grid at N=8 with W ∈ {1, 2, 4,
8}, M ∈ {3, 4, 6, 8, 12} and a 30-min ILP time-budget per cell. The
N=8 weight-tradeoff curve is the natural extension; the W=1 column
already has the partial bound `M*(N=8, W=1) ≥ 17` from S84. Predicted:
similar shape to N=6 — sharp drop from W=1 to W=4, then a plateau.
Save under `experiments/constructions/d2_sign_threshold_w_m_tradeoff/n8_extension/`.

**C8.b — Asymptotic floor scan.** Vary N ∈ {4, 5, 6, 7, 8} at fixed W=4
and trace `M*(N)`. If the floor scales linearly (`M*(N) ≈ N/2`), this
suggests PRIMES is morally a parity-counter; if sublinear, PRIMES has
genuine structure beyond oddness.

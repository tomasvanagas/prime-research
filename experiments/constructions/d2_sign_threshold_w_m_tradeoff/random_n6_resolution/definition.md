# C8.b — Random-Control F4 Resolution at N=6 (definition)

## Edges composed

- **E5.3** — PRIMES TC⁰ membership is the only open frontier;
  depth-2 sign-threshold size is a quantitative TC⁰₂ proxy.
- **S84 framework** (`enum_d2_smart.py`) — column enumeration of the
  bottom-layer threshold truth tables, then ILP for the top-layer
  selection. The W=1 column on N=6 random gave `M*(rand) = 7` cleanly
  in 47 s (`enum_d2_smart_n6_rand.log`).
- **C7-S89 (E1.6)** — PRIMES ≈ oddness for x > 2; the calibrated random
  control at N=6 absorbed the S84 PRIMES-vs-random gap. Question
  restated: does the gap survive at W ≥ 2 *outside* density-and-
  oddness calibration, i.e. is matched-density-only random strictly
  harder than PRIMES at W ≥ 2?
- **C8/S127 partial result** — PRIMES `M*(N=6) = 6, 4, 3, 3` at
  `W = 1, 2, 3, 4`. Random `(W=4, M=3)` UNKNOWN at 600 s via the
  direct (joint) ILP. The C8 verdict on F4 at N=6 is **unresolved**.

## Object

For a fixed weight bound `W ∈ {1, 2, 3}` and a target Boolean function
`f : {0, 1}^N → {0, 1}`, build the catalog

```
  Theta(N, W) := { sign-threshold truth tables tt : ∃ w ∈ {-W..W}^N,
                                                    T ∈ ℤ,
                                                  tt(x) = [<w,x> ≥ T] }
```

(deduplicated up to complementation; see `extended_enum.py`). Then
solve

```
  M*(W; f) = min M : ∃ {θ_1, ..., θ_M} ⊂ Θ(N, W),
                       α ∈ {-1, +1}^M, T_top ∈ ℤ,
                       f(x) = [Σ_j α_j · θ_j(x) ≥ T_top].
```

The catalog sizes (computed in this session):
| W | K = #Θ(6, W) |
|---|--------------|
| 1 | 1458 |
| 2 | 30,898 |
| 3 | 218,066 |

Only W ∈ {1, 2} are scanned in this session (the W=2 ILP at K=31K
is borderline tractable for CBC at M=3; W=3 is left for follow-up if
W=2 alone resolves F4 conclusively).

The composition is non-trivial because S127's joint (alpha-bilinear)
ILP times out on (W ≥ 4, M=3) random cells. By pre-fixing the bottom-
layer catalog we eliminate `M · K_w` bilinear constraints (where
`K_w = (2W+1)^N` raw before dedup) and keep only the `K` + `2^N`
constraints needed for top-layer selection.

## Relationship to π(x)

C8.b restates F4 (PRIMES easier than random at W ≥ 2) at N=6 with a
ILP technique that *can* terminate. The headline question:

> Does there exist a depth-2 sign-threshold circuit of size 3
> with bottom-layer weight bound W ∈ {2, 3} that computes a
> matched-density random function on `{0, ..., 63}`?

Resolution shapes:
- (R-NO) random `M*(W=2; N=6) ≥ 4` — F4 holds with Δ = 1 gate at W=2.
  PRIMES retains a structural advantage independent of the bit_0
  parity predictor (since calibrated-density random does *not* match
  oddness mass; this run uses `random_table` from S84 which is
  density-matched but not bit_0-matched).
- (R-YES) random `M*(W=2; N=6) = 3` — F4 collapses at N=6 W=2; PRIMES
  has no advantage over density-matched random in the W ≥ 2 regime.
  The S84 W=1 gap (PRIMES 6, random 7) was a one-off.

## Pre-stated falsification (BEFORE running)

These criteria are stated before any new ILP cell is run. The PRIMES
sanity rows are predictions from S127's curve.

**F0 (sanity).** PRIMES `M*(W=1; 6)` reproduces S84's value of 6.
Specifically: PRIMES at W=1 K=1458 returns `M=5 UNSAT, M=6 SAT`.
**F0' (sanity).** PRIMES at W=2 K=30898 returns `M=3 UNSAT, M=4 SAT`,
matching S127's joint-ILP `M*(W=2; 6) = 4`.

**F1 (random easier at W=2).** Random N=6 seed=42 at W=2 returns
`M=3 SAT` (i.e., `M*(rand; W=2) ≤ 3`). FALSIFIES F4 at W=2.

**F2 (random hard at W=2).** Random N=6 seed=42 at W=2 returns
`M=3 UNSAT` (proven in time-budget). Combined with the trivial
M=2 UNSAT (no depth-2 circuit of size 2 can have output multiplicity
≥ 18/64 except in degenerate cases, but we prove it anyway), this
gives `M*(rand; W=2) ≥ 4 = M*(PRIMES; W=2)` — F4 holds at W=2 with
Δ ≥ 0 gate gap.

**F3 (random strictly harder at W=2).** Random `M=4 UNSAT` at W=2 in
time-budget. Combined with F2, gives `M*(rand; W=2) ≥ 5`, while
PRIMES `M*(W=2) = 4`. F4 holds at W=2 with **Δ ≥ 1 gate gap**.

**F4 (cross-seed robustness).** A second seed (seed=1, S127's seed)
gives the same M=3 SAT-vs-UNSAT verdict on random at W=2 as seed=42.

## Scope and limits

- Single-N (N=6), single-depth (2), single decoder (sign-threshold).
  Conclusions do not extend automatically to other Boolean function
  classes; the F4 question at N=8 is left for C8.a.
- Random is *density-matched* (weight 18/64) but not *bit_0-matched*
  (oddness of PRIMES is an independent variable). C7-S89 already
  showed the calibrated (density + oddness)-matched random at N=6
  has `M = 6` for some samples, so the F4 gap *cannot* be more
  than 1 gate at W=1; the W ≥ 2 regime is the open question.
- W=3 (K=218K) is reserved for follow-up if W=2 leaves F4 ambiguous
  in the same way W=4 left it ambiguous in C8.

## Save under

`experiments/constructions/d2_sign_threshold_w_m_tradeoff/random_n6_resolution/`

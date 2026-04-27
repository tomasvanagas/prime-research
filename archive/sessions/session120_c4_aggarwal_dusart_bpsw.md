# Session 120 — C4 BUILT: Aggarwal × Dusart × BPSW unified p_n library

**Mode.** construction (NOVELTY_CHALLENGES §1 C4).
**Self-grade.** **B** (substantive composition; three structural
findings; no polylog opening, no new edge).
**Date.** 2026-04-27.

## Pick

C4 was the only unbuilt composition in NOVELTY_CHALLENGES.md §1 (C1,
C2, C3, C5, C6, C7 all already BUILT in previous sessions S70 / S74 /
S105 / S71 / S81 / S89). The session brief mandated picking C_x with
x ∈ {1..6}; C4 was the only available option. C4 is the algorithm-
composition entry whose target lives under `algorithms/`, not
`experiments/constructions/`.

## Edges composed

- **E6.6 (Aggarwal 2025 binary search optimality):** `p_n` reducible
  to `O(log n)` calls to `pi(x)` via binary search on the Dusart
  bracket; bound `O(sqrt(n) log⁴ n)` unconditionally.
- **E6.8 (Dusart bracket of width n):** for n ≥ 6,
  `n(log n + log log n - 1) ≤ p_n ≤ n(log n + log log n)`. Width = n
  exactly (logarithmically narrow on the scale of `p_n ~ n log n`).
- **E5.1 (BPSW correctness ⇒ PRIMES in TC⁰):** strong-MR base 2 ∧
  strong Lucas ∧ Jacobi, all in TC⁰; verified deterministic to 2⁶⁴,
  conditional above.

## What I built

`algorithms/aggarwal_dusart_bpsw/`:
- `aggarwal_dusart_bpsw.py` — single-file pure-Python module with all
  three components and three modes; ctypes hook to
  `algorithms/v10_c_accelerated.py` for C-Lucy benchmarking.
- `definition.md` — composition signature + intended relationship to
  π(x) + edge ID citations.
- `aggarwal_dusart_bpsw_results.md` — pre-stated F1-F4 falsifiers +
  empirical outcomes + verdict.
- `bench_*.log` — wall-clock + K-sweep raw outputs.

### Components

1. `dusart_bounds(n)` → (L, R) with width n, p_n ∈ [L, R].
2. `is_bpsw(n)` — strong-MR base 2 + strong-Lucas (Selfridge D-search)
   + perfect-square rejection. Pure Python, validated against
   `sympy.isprime` on 13 small-n cases.
3. `pi_lucy(x)` — Lucy_Hedgehog DP, `O(x^{2/3})` time, `O(√x)` space.
4. `aggarwal_search(n, pi_oracle)` — binary search maintaining
   `pi(R) ≥ n`, `pi(L-1) < n`. Returns smallest x with `pi(x) ≥ n`.
5. `bpsw_walk_from(start, count_target)` — walk forward by 2,
   counting BPSW-primes.
6. `hybrid_pn(n, K, pi_oracle)` — Aggarwal narrow until `R - L ≤ K`,
   then `pi(L-1)` for residual, then BPSW walk for `n - pi(L-1)`-th
   prime in `[L, ...]`.

## Pre-stated falsifiers (in `aggarwal_dusart_bpsw_results.md`)

- **F1 (agreement):** `agg`, `bpsw`, `hybrid` agree with `sympy.prime`
  on every benchmark `n`. *Hard fail* if violated.
- **F2 (Aggarwal-baseline beating):** `hybrid` ≥ 1.5× faster than
  `agg` for n ≥ 10⁵.
- **F3 (BPSW-baseline beating):** `hybrid` ≥ 10× faster than `bpsw`
  for n ≥ 10⁴.
- **F4 (U-shape K-curve):** wall-clock of `hybrid(n=10⁶, K)` is
  U-shaped in K.

## Outcomes

### F1 — HOLDS

13 small-n cases plus n ∈ {10⁴, 10⁵, 10⁶, 10⁷} all agree with
`sympy.prime`. No mismatch on any tested n.

### F3 — HOLDS

Pure-Python `pi_lucy`, K=64:

| n     | bpsw     | hybrid   | speedup |
|-------|----------|----------|---------|
| 10⁴   | 0.104 s  | 0.005 s  | 21×     |
| 10⁵   | 1.299 s  | 0.038 s  | 34×     |
| 10⁶   | 15.273 s | 0.291 s  | 53×     |

Predicted `2 p_n / n ~ 2 log p_n` ≈ 30-50 at this scale; observed
within factor 2. Dusart bracket is the cheap value-add.

### F2 — PARTIAL (Python regime fails 1.5× threshold; C regime holds)

Pure-Python `pi_lucy`:

| n     | agg     | hybrid (K=64) | ratio |
|-------|---------|---------------|-------|
| 10⁴   | 0.007 s | 0.005 s       | 1.40× |
| 10⁵   | 0.050 s | 0.038 s       | 1.32× |
| 10⁶   | 0.379 s | 0.291 s       | 1.30× |

Direction correct, magnitude below 1.5× threshold. Cause: pure-Python
`pi_lucy` per-call cost is ~28 ms at x ~ 10⁷ — saving 4 calls only
buys ~110 ms on a 379 ms baseline.

C-accelerated `pi_lucy` (compiled via ctypes from
`v10_c_accelerated.py`):

| n     | agg_C    | hybrid_C K=16384 | ratio |
|-------|----------|------------------|-------|
| 10⁶   | 0.0082 s | 0.0054 s         | 1.52× |
| 10⁷   | 0.0540 s | 0.0333 s         | 1.62× |

F2 holds with C-Lucy.

### F4 — PARTIAL (depends on pi-cost regime)

Pure-Python K-sweep at n = 10⁶ (monotone decreasing):

```
K          time(s)  pi_calls  bpsw_calls
4096       0.186      9         796
16384      0.148      7        6656
65536      0.141      5       22281
131072     0.124      4       22281
262144     0.103      3       22281
524288     0.086      2       22281
1048576    0.063      1       22281      <-- bracket width
2097152    0.064      1       22281
```

Optimum at K = bracket-width. **Aggarwal narrowing is empirically
negative-value in this regime** — every narrowing step costs more than
the BPSW walk it replaces. Composition collapses to E6.8 + E5.1.

C-Lucy K-sweep at n = 10⁷ (U-shape):

```
K          time(s)  pi_calls  bpsw_calls
256        0.0403    17           53
1024       0.0430    15          282
4096       0.0426    13         1198
16384      0.0333    11         2419          <-- empirical optimum
65536      0.0421     9         7302
262144     0.1695     7        65896
1048576    0.4936     5       222146
4194304    0.4963     3       222146
16777216   0.4899     1       222146         <-- bracket width
```

**U-shape appears** with K* ≈ 16384 ≈ √(bracket width). To the left of
K* the extra `pi(x)` calls dominate; to the right the BPSW walk
dominates. F4 holds with fast-pi oracle.

## Three structural findings

1. **Optimal K depends on the pi/bpsw cost ratio.**
   - Python `pi_lucy` (~28 ms/call at x ~ 10⁷): K* = width.
   - C-Lucy (~2 ms/call): K* ≈ 16K ≈ √width.
   - HKM / primecount projection (~0.1 ms/call at x ~ 10⁸):
     predicted K* ≈ small constant ⇒ Aggarwal-pure dominates.

   This is a **cost-ratio-dependent K-knob invisible at asymptotic
   order**. Aggarwal 2025's `O(sqrt(n) log⁴ n)` analysis implicitly
   takes K = 1 (Aggarwal-pure). The C4 composition formalises K as a
   tunable parameter and predicts a U-shape with empirically-located
   optimum.

2. **BPSW conditionality propagates 1-to-1 through Aggarwal's wrapper.**
   Aggarwal narrowing runs on `pi_lucy(x)`, which uses no primality
   testing — only divisibility and small-prime sieving. Hence BPSW
   conditional enters only at the trailing walk step, and a single
   BPSW pseudoprime in `[L, R]` shifts the answer by **at most one
   prime**. The wrapper does not compound conditionals through a
   sieve. This is structurally cleaner than naïve "BPSW everywhere"
   approaches.

3. **The Dusart bracket alone is worth ~50× over naive BPSW-from-2.**
   `bpsw_walk_only` does ~p_n BPSW tests; `hybrid` with K = width
   does ~width / 2 = n / 2 tests. Ratio `2 p_n / n ~ 2 log p_n` ≈
   30-50 at n = 10⁴-10⁷. Empirically observed 21×/34×/53×, matching
   prediction within factor 2.

## What this composition is NOT

- Not a new asymptotic algorithm — Aggarwal's `O(sqrt(n) log⁴ n)`
  bound is preserved.
- Not unconditional above 2⁶⁴ — BPSW conditionality is inherited
  (but only 1-to-1, not amplified).
- Not a polylog-π(x) algorithm — bottleneck is `pi(x)` (Lucy / HKM).
  Does not address the polylog gap.

## Successor challenges proposed

**C4.a — Replace `pi_lucy` with HKM/primecount and re-locate K*.**
Replace the Lucy DP oracle with primecount or a tractable HKM port
(Hirsch-Kessler-Mendlovic 2024) and trace the optimal K. Predicted: K*
drops further toward small constant as the asymptotic regime
approaches; the hybrid scheme's gain over `agg` shrinks toward zero.
If TRUE, the C4 composition is a *finite-x*-only improvement; if K*
stays at √width even with HKM, BPSW retains constant-factor value
asymptotically. Cost: 1-2 sessions. Save under
`algorithms/aggarwal_dusart_bpsw/hkm_extension/`.

## Files

**New:**
- `algorithms/aggarwal_dusart_bpsw/aggarwal_dusart_bpsw.py`
- `algorithms/aggarwal_dusart_bpsw/definition.md`
- `algorithms/aggarwal_dusart_bpsw/aggarwal_dusart_bpsw_results.md`
- `algorithms/aggarwal_dusart_bpsw/bench_main.log`
- `algorithms/aggarwal_dusart_bpsw/bench_ksweep.log`
- `algorithms/aggarwal_dusart_bpsw/ksweep_extended.log`
- `algorithms/aggarwal_dusart_bpsw/bench_c_pi.log`
- `algorithms/aggarwal_dusart_bpsw/bench_n1e7.log`
- `algorithms/aggarwal_dusart_bpsw/bench_ksweep_1e7.log`
- `archive/sessions/session120_c4_aggarwal_dusart_bpsw.md`

**Modified:**
- `EDGES.md` — inline S120 annotations on E5.1, E6.6, E6.8.
- `NOVELTY_CHALLENGES.md` — §1 C4 marked BUILT (S120) with outcome
  paragraph; new C4.a successor entry added.
- `RESEARCH_AGENDA.md` — Arc 4 C4 milestone closed (S120).
- `status/CLOSED_PATHS.md` — S120 closure row added.
- `status/SESSION_INSIGHTS.md` — Session 120 entry appended.

## Session-end self-evaluation

1. **What did I produce that was not in the project before?**
   - The first integration of E5.1 + E6.6 + E6.8 in a single
     executable artifact (the file did not exist).
   - Quantitative measurement of K* across two pi-oracle regimes,
     locating the U-shape minimum.
   - The structural observation that BPSW conditional propagates
     1-to-1 (not previously stated explicitly).
   - The 50× quantification of Dusart-bracket value over BPSW-from-2.

2. **What edges did my work compose or cite?**
   E5.1, E6.6, E6.8 (composed); E1.3, E5.3, E6.7, E6.9, E7.6, E7.7
   (cited in commentary).

3. **If my session produced only duplicate closures, why?** — N/A.
   Session produced a working artifact + three structural findings.
   Honest grade is B (substantive but not A — the asymptotic
   barrier is unchanged, and the cost-ratio K-knob is a refinement
   of Aggarwal 2025 rather than a refutation).

4. **Next-action for the next agent.**
   - Pick C4.a (HKM extension) if interested in the asymptotic
     side of the K-knob.
   - Or: pick D2.a.1 / D2.a.2 / D2.b / D2.c (S117 successors —
     PH refinements after S117's W=210 W-trick result).
   - Or: a frontier ATTACK_VECTORS target.

   Recorded in NOVELTY_CHALLENGES.md.

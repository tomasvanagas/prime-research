# chain_profile.py — results (S501)

**Question (from PROGRAM.md NEXT ACTION).** The S500 correction said the full
compressed delegated+structured π(x) chain (`compressed_layer.run_chain`) is
"op-count-bound on many small per-layer cubes," so the lever for a fast large-x
run is reducing the *count* of small numpy field-ops. To know *which* small-cube
kernel to batch, attribute `run_chain`'s wall-clock across kernels.

**Method.** Read-only profiler, two passes over the same deterministic run:
1. **cProfile** → per-function `tottime` (exclusive, C-level, non-double-counted),
   `cumtime`, `ncalls`, plus the direct-caller breakdown of the two hot
   functions (`fmul`, `sumcheck`).
2. **size pass** → a light monkeypatch (count-only, no timing) recording the
   mean/max array size per primitive — the one fact cProfile cannot give.
   `sumcheck` calls are bucketed by `(nvars, degree, n_terms)`.

Run: `python3 chain_profile.py --n 16 --field big` (BIG_Q = 2⁶¹−1, object dtype
— the *sound* production field, n ≲ 60). `--field q` (2³¹−1, uint64) for
contrast; `--selftest` for the checks below.

---

## Headline finding 0 — the NEXT-ACTION premise was wrong (corrected)

PROGRAM.md asked this profiler to "assert the profile sums match `run_chain`'s
`t_prover + t_verifier`." **That assertion is false** and the code documents
why: `sumcheck()` never touches `stats`, and **every** `sumcheck` call in the
chain sits *outside* any `t_prover`/`t_verifier` accounting region. The
sum-check loop is the entire cost, so:

| n=16, BIG_Q | wall | t_prover+t_verifier | accounted |
|---|---|---|---|
| reported stats | 15.3 s | 0.99 s | **6.4 %** |
| (uint64 q) | 8.3 s | 0.59 s | **7.1 %** |

The per-layer `t_prover`/`t_verifier` numbers printed by earlier cycles
(`compressed_layer.py --n …`, the `bench` columns) therefore **undercount the
true wall by ~15×**. They measure only the verifier algebra + the prover table
*builds*, not the sum-check folding that dominates. The honest checks the
selftest makes instead: faithfulness (instrumented run reproduces claimed π,
comm, accept), the two instruments agreeing on call counts, and the named
kernels covering ~all the wall.

---

## Headline finding 1 — `fmul` is the wall; mean array size ≈ 30

n = 16, BIG_Q (object), K = 54 layers, nb = 8, profiled wall = 16.2 s:

| kernel | tottime | %wall | cumtime | %wall | ncalls | mean_sz | max_sz |
|---|---|---|---|---|---|---|---|
| `fmul` | 9012 ms | **55.6%** | 9012 ms | 55.6% | **2 324 860** | **29.7** | 256 |
| `sumcheck` | 4395 ms | **27.1%** | 15072 ms | 93.0% | 2249 | 171 | 256 |
| `lagrange_eval` | 238 ms | 1.5% | 1035 ms | 6.4% | 15 670 | – | – |
| `_asum` | 180 ms | 1.1% | 1071 ms | 6.6% | 316 068 | 31.4 | 256 |
| `chain_layer_reduce_structured` | 110 ms | 0.7% | 4845 ms | 29.9% | 856 | – | – |
| `eq_table` | 91 ms | 0.6% | 485 ms | 3.0% | 4500 | 216 | 256 |
| named kernels (tottime) | | **88.3%** | | | | | |

`fmul` + `sumcheck`-self = **82.7 %** of wall. `fmul`'s mean operand is **~30
elements** (max 256 = the nb-cube). 96.7 % of `fmul` calls come straight from
`sumcheck`'s inner loop. **This is exactly the op-count-bound signature S500
predicted:** the run is ~2.3 million tiny element-wise field multiplies, each
paying full Python + numpy per-call dispatch on a ~30-element array. The uint64
field is the same shape (fmul 47 %, same 29.5-elem mean, same ~2.3 M calls);
object is ~1.8× slower wall because each Python-int `fmul` is dearer relative to
the loop overhead.

## Headline finding 2 — logical kernels: trace zero-test is fattest, wiring is most-called

By `cumtime` (inclusive), and by who issues the `sumcheck` calls:

| logical kernel | cumtime %wall | #sumcheck calls | shape |
|---|---|---|---|
| `verify_trace_region` (`verify_constraints`) | **61.3% (53.4%)** | 108 (+54) | 54 calls, deg-3, **up to 133 terms** (U·a+R−X zero-test + recompositions + bitness over Lv+2Lr bit tables) |
| `inner_verify_div` (`chain_layer_reduce_structured`) | 32.0% (29.9%) | 107 (+**1712**) | wiring delegation — **76 % of all sumcheck calls**, but tiny 2^l-cubes (size 2–256) |
| `large_reduce` | 81.2% | 54 | contains affine + trace + line-batch |
| `small_reduce` | 17.6% | 53 | phase-A + phase-B + division wiring |

Two distinct cost shapes:
- **Trace zero-test** (`verify_constraints`): *few* calls (1/layer) but each is
  *fat* — a degree-3 sum-check with ~50–133 terms × ~40 bit-tables. **53 % of
  wall** in 54 calls. The single fattest target.
- **Wiring delegation** (`chain_layer_reduce_structured`): *many* calls (1712,
  76 % of the total) on *tiny* cubes (nvars = ⌈log₂p⌉ ∈ [1, 8]). **30 % of
  wall** in pure per-call dispatch.

The full sum-check bucket table (`--n 16 --field big`) confirms the cube sizes:
the deg-3 / 133-term rows are all nb-cubes (size 256); the deg-5 / 1-term rows
(the structured wiring) run on sizes 2, 4, 8, …, 256.

---

## The highest-leverage batching target

**Batch the K per-layer trace zero-tests (`verify_constraints`, 53 % of wall)
into one wide sum-check.** They are *independent*: each layer's
`verify_constraints` certifies its own witness `build_witness(x, pₗ, nb)` and
draws its own `tau`/`alpha` — none depends on the carried chain claim or on
another layer's challenges. So a standard **random-linear-combination batched
sum-check** over the stacked `(K · 2^nb)` cube replaces K independent
nb-cube sum-checks with one. That cuts the `fmul` *count* ≈ K-fold while
*widening* each `fmul` ≈ K-fold — converting op-count-bound to width-bound,
which is precisely the regime where numpy and `fmul`'s vectorised Mersenne path
(S500) finally pay off. **Batching is the missing precondition that makes the
fast Mersenne path a net win** (S500 found it loses on the chain *because* the
arrays are small; widen them and the sign flips).

Second target: the wiring delegation (32 %, 76 % of calls). Harder — each
layer's `inner_verify_div` is invoked at a point `(r_v, r_u)` produced
*sequentially* by that layer's outer sum-checks, so the 107 calls cannot be
batched in place. A **defer-and-batch** restructuring (run the chain to collect
all `(r_v, r_u, claim)` wiring obligations, then discharge them in one batched
inner protocol at the end) would work but is a larger change. The trace
zero-test is the clean first win.

Estimated payoff (object dtype, where per-call dispatch is a smaller share):
batching mostly removes *dispatch*, not element-ops, so the object chain gains
modestly. The large win is on the **uint64 / fast-Mersenne** path, where
dispatch dominates and a K-fold-wider array amortises `_mul61`'s ~24 vectorised
ops/multiply — the combination (batch + `--fast-big`) is the route to a fast,
sound large-x run.

---

## Selftest (`--selftest`, ~19 s)

1. **Faithfulness** — the size-instrumented run reproduces the reference
   `run_chain` (claimed π = sieve, identical `comm`, identical accept) for
   n ∈ {8,10}, fields {q, BIG_Q}, modes {delegate+structured, automaton}. The
   monkeypatch only counts; it never alters arithmetic.
2. **Instrument agreement** — cProfile `ncalls` == size-pass counts for
   `fmul`, `sumcheck`, `eq_table`, `mle_eval`, `_asum` (both instruments see
   the same deterministic calls).
3. **The gap** — `(t_prover+t_verifier)/wall < 0.30` (measured 6–7 %).
4. **Attribution** — named kernels cover > 85 % of profiled wall; `fmul` is the
   single largest `tottime`; `fmul` mean array size < 200 (op-count bound).
5. **Dominance** — `sumcheck` cumtime > 50 % of wall.

## What would falsify these results

- The instrumented run disagreeing with the reference (π / comm / accept).
- cProfile `ncalls` ≠ size-pass counts for a primitive.
- Named kernels failing to cover the bulk of the wall (a hidden cost).
- `fmul` *not* dominant, or its mean array size *large* — which would refute the
  "op-count-bound on small arrays" reading and point instead at array width.
- `t_prover+t_verifier` actually tracking the wall — which would mean the
  accounting gap (finding 0) is an artifact, not real.

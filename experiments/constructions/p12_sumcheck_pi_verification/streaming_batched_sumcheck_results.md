# streaming_batched_sumcheck.py — results (S526)

## The question (S525 NEXT ACTION)

The compressed π(x) chain's two batched discharges — the trace zero-test
(`batched_trace.verify_constraints_batched`, S502) and the Ub-bit openings
(`verify_ub_openings_batched`, S506) — each STACK all `K = π(√x)` per-layer
witnesses on a layer axis and run ONE sum-check over the joint `(Lk + nb)`-cube
of size `N = 2^Lk · 2^nb ≈ x`. S525 localized the n=24 prover's ~3.9 GB peak
RSS to exactly this: the K-witness LIST (`Ws = [build_witness(p) for p in
primes]`, ~1.4 GB) PLUS the stacked N-cube the sum-check folds (~2.6 GB), both
Θ(x), driving the super-linear n22→n24 WALL.

The S525 NEXT ACTION **conjectured** both checks are pure γ-RLC folds
`Σ_w γ^w (B_w·C_w)` that need only a running accumulator, so streaming one
freshly-built witness at a time would drop peak RSS to Õ(√x) **bit-identically**.

This cycle tests that conjecture, **and corrects it.**

## What was built

`streaming_batched_sumcheck.py` (one script, `--selftest`/`--bench`/`--mem-probe`/
`--mem-one`/`--n`/`--field`):

- `build_stacked_streaming` — builds the stacked cube slice-by-slice, each
  witness materialized, copied into its `[l·D, (l+1)·D)` slice, and DROPPED.
- `verify_constraints_streamed` / `verify_ub_openings_streamed` — a **bind-e-first,
  multi-pass** streaming sum-check that reaches Õ(√x) peak working set.
- `_round_evals` / `_fold` — the sum-check inner round, routed through
  `_cpmt.fmul` so `_FmulCounter` and the live `FAST_BIG` toggle both apply.

## Findings

### 1. The K-witness LIST streams BIT-IDENTICALLY and for free. — LANDED (S527)

`build_stacked_streaming` produces the **exact same tables** as
`batched_trace.stacked_tables` — every table, honest and all four cheat classes,
over `q = 2³¹−1` and `BIG_Q = 2⁶¹−1` (selftest 1). It never holds more than one
witness, so the ~1.4 GB list component is removed with **no change to the
transcript**. This part of the S525 premise is correct.

**LANDED on the real chain (S527).** `build_stacked_streaming` was moved into
`batched_trace.py` (this module now imports it) and wired through
`compressed_layer.run_chain(stream_witnesses=True)` — the **default**, so the
chain no longer materializes `Ws = [build_witness(x, p) for p in primes]` for the
two batched discharges. Both `verify_constraints_batched` and
`verify_ub_openings_batched` grew a `stream=` switch. Bit-identical ⇒ a verbatim
drop-in: `large_x_benchmark` selftest 6 asserts full-config
`claimed/comm/comm_bt/comm_ub/accepted` are identical ON vs OFF over Q & BIG_Q.
The peak-RSS win is the `--stream-probe` A/B (FULL vs FULL_NOSTREAM) in
`large_x_benchmark_results.md`. The other two findings (the cube is not a pure
fold; multi-pass streaming costs Θ(√x)× wall) stand unchanged — the cube is NOT
landed.

### 2. The stacked CUBE is NOT a pure γ-fold (premise correction).

The trace test is a degree-3 constraint zero-test and the Ub test a degree-2
product `B·C`; the sum-check's ROUND POLYNOMIALS over the layer axis couple
multiple witnesses through the pointwise PRODUCT, and — decisively — each round's
Fiat–Shamir challenge depends on a sum over **all K witnesses**. A one-pass
running accumulator cannot produce the round polynomials. So Õ(√x) space requires
EITHER the Θ(x) cube (one pass) OR re-derivation (multi-pass): a genuine
space-time tradeoff, not a free fold.

### 3. A streaming sum-check DOES reach Õ(√x) space, verdict-faithfully.

Bind the e-axis FIRST. The summand is block-diagonal in the layer index, so each
e-round polynomial is an **additive** sum over the K layers. Each e-round is one
streaming pass (build, fold to depth, accumulate, drop); after the nb e-rounds
each layer collapses to per-table SCALARS over the size-K layer axis, and the Lk
layer-rounds are a small in-memory sum-check. Because the MLE value at the final
point is order-independent, the bound scalars, the verifier's final check, and
the **VERDICT** are identical to the stacked test (the transcript is NOT
byte-identical — e-first vs the stacked layer-first round order).

Validated (selftests 2–6): streamed honest **ACCEPT** matches the stacked test,
every single-layer cheat (trace: `u_consistent`/`u_value`/`r_value`/`nonbit`;
Ub: any `(layer,bit)` forge) is **REJECTED** at first/middle/last layer, over
`q` AND `BIG_Q`; the K=1 edge (x=7) streams; fast-Mersenne uint64 == object.

**Peak RSS (`--mem-probe 12,14,16,18 --field big`, fresh process each):**

| n  | K  | stacked MB | streamed MB |
|----|----|-----------:|------------:|
| 12 | 18 | 37 | 36 |
| 14 | 31 | 39 | 36 |
| 16 | 54 | 49 | 36 |
| 18 | 97 | 91 | 36 |

Streamed peak RSS is **FLAT at the ~36 MB numpy baseline** (growth 1.00–1.01×/Δn=2);
the streaming working set `D + K` is negligible against it. Stacked grows
1.05 → 1.26 → 1.86×/Δn=2 — the baseline-subtracted cube/list working set
(≈1 → 3 → 13 → 55 MB) is climbing toward the Θ(x) `×4`/Δn=2 law as it overtakes
the baseline. **The Õ(√x)-space win is real and measured.**

### 4. The cost is Θ(√x) MORE WALL (measured) — worse than a log factor.

The fmul ELEMENT-work is comparable (the round-poly work telescopes to the same
total; the streamed e-phase even SKIPS the `2^Lk − K` padding layers the stacked
cube folds), BUT the streamed prover makes ~K× MORE fmul CALLS, each on a NARROW
(D-wide) per-layer array instead of the stacked test's few N-wide calls. Wall is
dominated by per-call Python+numpy DISPATCH overhead, so it tracks the CALL count.

**`--bench --field big` (trace zero-test, BIG_Q + fast Mersenne):**

| n  | K  | stk_calls | str_calls | call× | stk_w | str_w | stk_ms | str_ms | slow× |
|----|----|----------:|----------:|------:|------:|------:|-------:|-------:|------:|
| 8  | 6  | 13 490 | 54 303 | 4.0× | 18.1 | 3.7 | 363.8 | 1306.1 | 3.6× |
| 10 | 11 | 21 340 | 147 992 | 6.9× | 56.8 | 6.3 | 566.5 | 2910.9 | 5.1× |
| 12 | 18 | 30 966 | 340 999 | 11.0× | 186.0 | 10.9 | 732.7 | 7281.0 | 9.9× |
| 14 | 31 | 39 110 | 786 190 | 20.1× | 341.1 | 19.1 | 1029.2 | 17665.4 | 17.2× |

The slowdown (3.6 → 17.2×, and 8 → 17× in the mem-probe wall column) tracks
**~K = Θ(√x)**, with mean array width crashing from ~341 (stacked) to ~19
(streamed). Streaming **RE-INTRODUCES exactly the op-count-bound regime
(S500/S501)** that batching (S502) was built to escape by widening arrays.

## Honest headline

The batched discharge is a **space-time tradeoff**, not a free γ-fold:

- **Batching (S502)** minimizes WALL via wide arrays at **Θ(x) space**.
- **Streaming (this cycle)** minimizes SPACE to **Õ(√x)** (peak RSS flat) at
  **Θ(√x)× more WALL** (narrow, dispatch-bound arrays).

You cannot have both. Since the n22→n24 wall S525 localized is itself a time cost
(memory bandwidth), paying ~√x× more compute to shrink the working set makes the
wall WORSE, not better — **the Θ(x) held cube is inherent to the one-pass batched
sum-check's wall efficiency.** The S525 NEXT ACTION's "stream the K-witness stack
for a free Õ(√x) win" is therefore only partly available: the K-witness LIST
streams free and bit-identically (item 1, removes ~1.4 GB at n=24), but the
~2.6 GB stacked CUBE — the dominant term — does not.

## What would falsify these results

- `build_stacked_streaming` differing from `stacked_tables` on any table
  (item 1 not bit-identical) — checked by selftest 1.
- The streamed verdict disagreeing with the stacked batched test on honest
  accept OR on rejecting ANY single-layer cheat, over q or BIG_Q — selftests 2–6.
- The streamed peak RSS NOT being flat while stacked grows ~`N = K·D` — the
  mem-probe table.
- The streamed wall NOT being Θ(√x)× the stacked wall / not tracking the ~K×
  fmul-CALL inflation — the bench table.

## Scope / what this is NOT

A finite-x measurement on the verification line, not progress on the polylog
π(x) goal (which stays blocked). It corrects a specific S525 conjecture and
quantifies the batched discharge's space-time frontier. It does NOT modify the
landed chain (`compressed_layer.run_chain`, `batched_trace.py`) — the n=24
artifact is untouched. `build_stacked_streaming` demonstrates the safe,
bit-identical LIST-streaming integration as a ready follow-on (the only free
peak-RSS win), but the dominant cube term cannot be removed without the
Θ(√x)× wall penalty measured here.

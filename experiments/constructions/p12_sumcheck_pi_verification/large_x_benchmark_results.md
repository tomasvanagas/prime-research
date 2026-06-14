# large_x_benchmark — reach results (open item 5)

**Script:** `large_x_benchmark.py` (driver; no new protocol machinery).
**What it runs:** ONE honest verification of `π(2ⁿ−1)` through the FULL
succinct config of `compressed_layer.run_chain`
(`delegate+structured+pcs+batch_trace+batch_ub+batch_wiring+commit_base`)
over the sound-characteristic field `BIG_Q = 2⁶¹−1`, with the uint64
Mersenne fast path (`FAST_BIG`) on; reports wall-clock + the proof's
size/verifier-op profile. Soundness is spot-checked with one `delta_pi`
liar (claim `π(x)+1`) that must be REJECTED.

This is the **reach** deliverable for the "large-x benchmark" line: push
the S504 `n=20` headline (66 s under the batched config) to larger `n`,
recording verified `π(2ⁿ)==sieve` and confirming the asymptotic profile
holds (per-layer verifier leaf-eval count exactly 0 ⇒ Õ(√x) end-to-end;
`comm` dominated by the K=π(√x) sequential outer reductions ⇒ Θ(√x) cert;
large-table opens only the one-time base commitment ⇒ Θ(x^{1/4})).

## Headlines (BIG_Q = 2⁶¹−1, FULL succinct config + FAST_BIG)

| n | x = 2ⁿ−1 | V=⌊√x⌋ | nb | K=π(√x)-ish layers | claimed π(x) | sieve | wall (honest) |
|---|---|---|---|---|---|---|---|
| 20 | 1 048 575 | 1023 | 10 | 172 | 82 025 | 82 025 ✓ | 71.8 s |
| 22 | 4 194 303 | 2047 | 11 | 309 | 295 947 | 295 947 ✓ | 252.8 s |
| 24 | 16 777 215 | 4095 | 12 | 564 | **1 077 871** | 1 077 871 ✓ | **1467.3 s** |

All rows: HONEST run **ACCEPTED**, `claimed π(x) == sieve_pi(x)` exactly.
`n=20` reproduces the S504 headline now under the FULL config (S504 used
the same batched bundle without `pcs`/`commit_base`; numbers agree:
π(1048575)=82025). The **n=24 claim π(16777215)=1 077 871** was
independently cross-checked against a from-scratch Eratosthenes sieve
(`π(2²⁴)=1 077 871`, a known value), so the chain's count is verified
against ground truth at two sources.

**Deterministic reproducibility (n=24).** The n=24 run was executed
twice (an earlier unlanded background run, `run_n24_full.log`, and the
landing run, `run_n24.log`), both at the default `seed=1`. Every
STRUCTURAL number is bit-for-bit identical — `claimed=1 077 871`,
`comm=292 635` with the same `comm_outer/base/bt/bw/ub =
285299/4352/88/2829/67` split, `vleaf_ops_pl=0`, one-time leaf =
base-commit opens = 16 512 — as expected since the Fiat–Shamir transcript
is a deterministic function of the seed. Only the wall-clock differs
(1467.3 s vs 1520.2 s, ~3.5% timing variance), and the `delta_pi` liar is
REJECTED in both (541.9 s / 554.9 s). The headline uses the landing run's
1467.3 s.

### Soundness spot-check (n=22, n=24)
`delta_pi` liar (prover claims `π(x)+1`) is **REJECTED** at both reaches:
n=22 in 66.9 s (claim 295948), **n=24 in 541.9 s** (claim 1 077 872) —
the corrupted top claim fails the chain's base check. (The n=24 liar is
NOT short-circuit-faster than honest the way n=22 was: at n=24 the
mismatch only surfaces after the bulk of the reduction has run, so the
liar costs ~37% of the honest wall rather than ~26%.) Soundness at the
reach is therefore witnessed, not assumed; the full 13-class cheat panel
lives in `compressed_layer.py`'s selftest.

## Certificate + verifier profile (n=24, the new fully-broken-down row)

```
comm = 292 635 field elems  (~2.23 MB @ 61-bit)
  comm_outer = 285 299  (97%)   K=564 SEQUENTIAL two-sided layer reductions  -> Θ(√x)
  comm_base  =   4 352          tensor-PCS eval proofs for the two S₀ closes  -> Θ(x^{1/4})
  comm_bw    =   2 829          ONE batched wiring discharge (was K chains)   -> polylog
  comm_bt    =      88          ONE batched trace zero-test  (was K)          -> polylog
  comm_ub    =      67          ONE batched Ub-opening discharge (was K·nb)   -> polylog
  (sum = 292 635, exact)

verifier large-table (mle_eval over a √x-size cube) ops:
  per-layer leaf   = 0          <- the Õ(√x) end-to-end property (S506); never grows with K
  one-time leaf    = 16 512
  base-commit opens= 16 512     <- the only x-scaling verifier work, Θ(x^{1/4}) (S508)

reported t_prover  ≈ 25.0 s     (timed regions only; ~1.7% of wall — wall is field arithmetic
reported t_verifier≈ 3.92 s      outside the timed regions, cf. S501)
```

## Certificate + verifier profile (n=22, the prior fully-broken-down row)

```
comm = 146 637 field elems  (~1.15 MB @ 61-bit)
  comm_outer = 139 899  (95%)   K=309 SEQUENTIAL two-sided layer reductions  -> Θ(√x)
  comm_base  =   4 224          tensor-PCS eval proofs for the two S₀ closes  -> Θ(x^{1/4})
  comm_bw    =   2 373          ONE batched wiring discharge (was K chains)   -> polylog
  comm_bt    =      80          ONE batched trace zero-test  (was K)          -> polylog
  comm_ub    =      61          ONE batched Ub-opening discharge (was K·nb)   -> polylog
  (sum = 146 637, exact)

verifier large-table (mle_eval over a √x-size cube) ops:
  per-layer leaf   = 0          <- the Õ(√x) end-to-end property (S506); never grows with K
  one-time leaf    = 12 352
  base-commit opens= 12 352     <- the only x-scaling verifier work, Θ(x^{1/4}) (S508)

reported t_prover  ≈ 8.15 s     (timed regions only; ~3% of wall — wall is field arithmetic
reported t_verifier≈ 1.76 s      outside the timed regions, cf. S501)
```

**Reading.** The headline confirms, at x≈4.2×10⁶ (n=22) and
x≈1.7×10⁷ (n=24), the three asymptotic claims the protocol stack was
built for:
1. **per-layer verifier leaf-eval count = 0** at every reach n ⇒ the
   per-layer verifier is leaf-free, so the only large-table verifier
   work is one-time (S506); the chain verifier is Õ(√x) end-to-end.
2. **`comm` ≈ 95% `comm_outer`** = the K sequential outer reductions;
   every batchable obligation is squeezed to polylog (`comm_bt/bw/ub` are
   80/2373/61 vs `comm_outer` 139 899). The Θ(√x) cert size is the
   sequential-sieve K layers, exactly as S509/S510 localized it.
3. **large-table opens = the base commitment only** (Θ(x^{1/4}), S508).

## What would falsify this result
- `claimed π(x) != sieve_pi(x)` at any n (the chain mis-counts).
- The honest run REJECTS.
- The `delta_pi` liar is ACCEPTED (soundness broken at the reach n).
- `FAST_BIG` ON changing `claimed`/`comm`/`accepted` vs OFF (selftest
  asserts the fast path is a bit-identical drop-in — checked at n≤12).
- Any per-layer verifier large-table leaf-eval (`vleaf_ops_pl > 0`)
  reappearing (would break the Õ(√x) end-to-end claim).

## Reproduction
```
cd experiments/constructions/p12_sumcheck_pi_verification
python3 large_x_benchmark.py --selftest         # fast correctness/soundness gate
python3 large_x_benchmark.py --n 20             # 72 s   reproduce S504 headline
python3 large_x_benchmark.py --n 22             # ~5 min headline + soundness spot-check
python3 large_x_benchmark.py --n 24             # ~24.5 min honest + ~9 min soundness (measured)
python3 large_x_benchmark.py --n 24 --no-cheat  # ~24.5 min larger headline, honest only
```
Field default `--field big` (BIG_Q=2⁶¹−1, sound to n≲60); `--no-fast`
forces object dtype (slower at these n; bit-identical verdict).

## Scaling law (3-point series, n=20/22/24) — corrects the S513 two-point claim

The wall-clock is **prover-bound** (reported `t_prover`+`t_verifier` is
only ~1.7% of wall at n=24; the rest is field arithmetic outside the
timed regions, cf. S501), not field- or DP-bound. With the third point
in hand the scaling is:

| step | wall ratio | local exponent in x |
|---|---|---|
| n=20→22 (x ×4) | 3.52× | x^0.908 |
| n=22→24 (x ×4) | **5.80×** | x^1.269 |
| least-squares over all 3 | — | **wall ≈ x^1.088** |

So the **per-Δn=2 ratio is NOT the constant ~3.5× S513 reported from the
single n20→n22 step — it ACCELERATES (3.52× → 5.80×, geo-mean 4.52×).**
The honest reading: the asymptotic prover is Õ(x) (open item 1); the
fitted exponent **~1.09 is linear-in-x plus polylog**, and the n=22→24
step runs *above* the ideal 4×/Δn=2 because a ~3.9 GB working set at n=24
adds a super-linear wall overhead on top of the Õ(x) op-count — an
**implementation** overhead, not a complexity change (the op-count
milestones per-layer leaf=0, comm=Θ(√x), opens=Θ(x^{1/4}) all hold at
n=24). **S524 mis-attributed the working set to the object-dtype
`np.add.at` scatter-sums; S525 below MEASURES that this was wrong** — the
3.9 GB is the Θ(x) batched-discharge held cubes (uint64), not the
√x-size per-layer scatter, and removing the object scatter does NOT move
the wall. **Falsifier for the Õ(x) prover claim** would be the local
exponent climbing and *staying* well above ~1.1 at n≥26 with the fast
path fully engaged (a genuine super-Õ(x) term, not a memory-hierarchy
constant). At present the 3-point fit x^1.09 is consistent with
Õ(x)+overhead.

## S525: vectorising the per-layer scatter-sums — bit-identical, but WALL-NEUTRAL; the S524 working-set attribution CORRECTED

S524 (above) blamed the super-linear `n=22→24` wall step on the two
per-layer outer-reduction scatter-sums falling back to object-dtype
`np.add.at` over BIG_Q (`verify_trace_region`'s `omega`, `small_reduce`'s
`Wt`), claiming they "push the working set to ~3.9 GB." S525 built the
uint64 fix and **measured that this attribution was wrong.**

**The fix (`scatter_fold61`, `compressed_prover_mult_trace.py`).** A
uint64-safe segmented scatter-fold: split each reduced BIG_Q residue into
a 31-bit lo limb + a 30-bit hi limb, sum each limb stream PER BUCKET with
one `argsort`+`np.add.reduceat` (fast C reduction, no object alloc), then
recombine `hi·2³¹+lo (mod 2⁶¹−1)`. Exact (bit-identical to the object
path) while bucket fill `< 2³³`; the callers' fill is `≤ √x = 2^{n/2}`, so
exact for `n<66`, past BIG_Q's `n≲60` ceiling. Wired in behind
`fast_big(q)`, selftested bit-identical vs object `np.add.at` over
random/structured (`arange//p`, low-bits)/heavy-collision/empty/gap/
2²⁰-near-overflow inputs. A `--no-scatter-fold` A/B control forces the
original object scatter on the fast path (selftest 5 guards it
bit-identical; a first cut of the toggle wrongly fell through to the
`_dt(q)`=uint64 branch, which OVERFLOWS over BIG_Q → a spurious REJECT —
now an explicit object branch).

**Memory localization — `--mem-probe` FULL vs NOBATCH (fresh process each;
NOBATCH = FULL without `batch_trace`/`batch_ub`, both run scatter_fold61):**

| n | K | FULL peak RSS | NOBATCH peak RSS | FULL wall | NOBATCH wall |
|---|---|---|---|---|---|
| 16 | 54 | 48 MB | 43 MB | 9.2 s | 34.9 s |
| 18 | 97 | 92 MB | 62 MB | 24.4 s | 86.5 s |
| 20 | 172 | **265 MB** | **138 MB** | 67.0 s | 197.5 s |

FULL peak RSS is consistently *higher* than NOBATCH (265 vs 138 MB at
n=20). So peak memory is driven by the **batched-discharge held cubes** —
the K-stacked witness list that `batch_trace`/`batch_ub` hold to discharge
their obligations in ONE sum-check, which is `Θ(K·2^nb·Lv)=Θ(x)` and
**uint64**. The per-layer scatter-sums are `√x`-size, freed each layer,
and are NOT the working-set driver. (Batching also costs ~3× less wall —
67 vs 197 s at n=20 — so it stays; the Θ(x) space is its price.)

**Scatter A/B at n=22 (honest, claimed=295947 ✓ both; SAME conditions):**

| path | wall | peak RSS | comm |
|---|---|---|---|
| ON (`scatter_fold61`, uint64) | 256.3 s | 988 MB | 146 637 |
| OFF (`--no-scatter-fold`, object) | 249.2 s | 988 MB | 146 637 |

`comm` is **identical** (bit-identical transcript confirmed in the real
benchmark, not just the selftest), peak RSS is **identical**, and the wall
differs by ~3% with OFF marginally *faster* — replacing the object scatter
with the uint64 fold changes nothing measurable, because the scatter is a
`√x`-size per-layer term whose `argsort` costs about what the object
`np.add.at` it replaces did, both negligible vs the `Θ(x)` field
arithmetic.

**Re-measured scaling, scatter_fold61 ON (sound, claimed==sieve every n):**

| n | wall ON (S525) | S524 wall (object) | peak RSS | local exponent (ON) |
|---|---|---|---|---|
| 16 | 9.2 s | — | 48 MB | — |
| 18 | 24.4 s | — | 92 MB | x^0.704 |
| 20 | 67.8 s | 71.8 s | 265 MB | x^0.737 |
| 22 | 256.3 s | 252.8 s | 988 MB | x^0.960 |
| 24 | **1445.7 s** | 1467.3 s | **3986 MB** | x^1.248 |

ON ≈ object at every reach (within ~2%). The local exponent **climbs
0.70 → 0.74 → 0.96 → 1.25** — and this climb is essentially identical to
the S524 object path (n20→22 x^0.91, n22→24 x^1.27), so scatter_fold61 did
NOT change the scaling. Peak RSS scales `~4×/Δn=2` (baseline-subtracted
ratios → 4.0) = **Θ(x)**; the 3.89 GB at n=24 is the held-cube working set,
present *with* the uint64 scatter active.

**Conclusion (the correction).** scatter_fold61 is a correct, bit-identical
optimization but **wall-neutral** at reachable n — the per-layer object
scatter was never the bottleneck. The ~3.9 GB working set and the
super-linear `n=22→24` step are NOT the object scatter (it is removed in
the ON path, yet n=24 peak RSS is still 3.89 GB and the acceleration
persists at x^1.25). They are the **`Θ(x)` batched-discharge held-cube
working set** (uint64): it grows `~4×/Δn=2`, and once it exceeds cache and
saturates RAM bandwidth the `Θ(x)` field arithmetic over it runs
super-linearly in wall — a memory-hierarchy constant, **not** Python-int
GC (the held cubes are uint64, not object) and **not** the scatter. The
Õ(x) op-count milestones (per-layer leaf=0, comm=Θ(√x), opens=Θ(x^{1/4}))
are unchanged. The real remaining wall lever for n≥24 is the held-cube
**space**: stream/checkpoint the K-witness stack so `batch_trace`/
`batch_ub` do not hold all K witnesses at once, dropping the prover working
set from Θ(x) toward Õ(√x).

**What would falsify this (S525):** scatter_fold61 not being bit-identical
to the object path (selftest); the n=22 A/B showing a wall gap beyond
timing noise (measured ~3%, OFF faster); or the FULL-vs-NOBATCH peak-RSS
ordering reversing (which would put the working set back in the per-layer
reductions). Reproduce with `python3 large_x_benchmark.py --mem-probe
16,18,20`, `--n 22 [--no-scatter-fold]`; full sweep in
`run_s525_remeasure.log`.

## S527: LIST-streaming landed on the real chain — the free, bit-identical peak-RSS win, MEASURED

S525 localized the n=24 prover's ~3.9 GB peak RSS to two `Θ(x)` working
sets: the K-witness **LIST** `Ws = [build_witness(x,p) for p in primes]`
(~1.4 GB) and the stacked sum-check **CUBE** (~2.6 GB). S526 proved only
the LIST removes for free (the cube is not a pure accumulator), via
`batched_trace.build_stacked_streaming` (builds the cube slice-by-slice:
materialize one witness → copy its `[l·D,(l+1)·D)` slice → drop), which is
table-for-table **== `stacked_tables`**.

This cycle LANDS that win on the real chain. `run_chain` gained
`stream_witnesses` (default **True**): under `batch_trace`/`batch_ub` it
passes `Ws=None` and routes both batched discharges through the streamed
builders (`verify_constraints_batched(stream=True)` /
`verify_ub_openings_batched(stream=True)`) so the K-witness list is never
materialized. The cube is **BIT-IDENTICAL**, so the whole certificate is a
verbatim drop-in — same claimed π, same `comm` (and its `comm_bt`/`comm_ub`
split), same accept/reject.

**Bit-identicality (verified, not assumed).** `compressed_layer` selftest
**case 25** and `large_x_benchmark` selftest **case 6** both assert
`stream_witnesses` ON == OFF on `claimed`, `comm`, `comm_bt`, `comm_ub`,
`comm_outer`, `accepted` at the SAME seed over q AND BIG_Q (FULL batched
config, n=8,10,12), and that the `delta_pi` / self-consistent liar is still
rejected under streaming. `batched_trace` selftest 8 asserts
`build_stacked_streaming == stacked_tables` table-for-table.

**Measured A/B peak RSS** (`--stream-probe 16,18,20,22`, FULL vs
FULL_NOSTREAM, fresh subprocess each, BIG_Q + FAST_BIG, seed 1):

| n  | K   | full=stream MB | full_nostream MB | LIST saved | saved % | claimed π | wall ≈ |
|----|-----|----------------|------------------|------------|---------|-----------|--------|
| 16 |  54 |  43            |  49              |   6        | 12%     | 6542      | 10 s   |
| 18 |  97 |  70            |  91              |  21        | 23%     | 23000     | 25 s   |
| 20 | 172 | 185            | 265              |  80        | 30%     | 82025     | 70 s   |
| 22 | 309 | 677            | 988              | 311        | 31%     | 295947    | 250 s  |

The **claimed π is bit-identical** between configs at every n (verbatim
drop-in confirmed at the headline scale too, not just selftest n≤12). The
LIST saving grows **6 → 21 → 80 → 311 MB**, ratios **3.5 / 3.8 / 3.9** →
the `4×/Δn=2 = Θ(x)` law — exactly the removed list's scaling. The
`full_nostream` peaks (265 MB @ n=20, 988 MB @ n=22) **reproduce S525's
pre-streaming FULL peaks bit-for-bit** (S525 measured FULL=265 MB @ n=20,
988 MB @ n=22 when the list was always held) — confirming the A/B isolates
exactly the list and nothing else. **Wall is unchanged** (n=20 69.8 vs
69.4 s; n=22 253.9 vs 247.3 s — within timing noise): the win is free.

**n=24 confirmation** (`--mem-one 24 --config full`, streaming): peak RSS
**[APPEND ON LAND] GB**, claimed π(16777215) = **[APPEND]** (== sieve
1 077 871). Compare S525's pre-streaming n=24 peak **3.89 GB** (the
list-held baseline, `run_s525_remeasure.log`): the list fraction climbs
toward the predicted ~36% (~1.4 GB of 3.9 GB).

**What it does and does NOT change.** Peak RSS only — soundness, verdict,
and the asymptotic profile are untouched (per-layer leaf=0, comm=Θ(√x),
opens=Θ(x^{1/4})). It does NOT touch the `Θ(x)` stacked **cube** (S526:
not a pure accumulator; streaming it costs `Θ(√x)×` more wall, so it makes
the wall worse — not landed). So the chain's prover working set drops from
`~LIST + CUBE` to `~CUBE` — a constant-factor (≈1/3) reduction at reachable
n, NOT an asymptotic Θ(x)→Õ(√x) change (both terms are Θ(x); only one is
removed). The remaining `Θ(x)` cube is the wall-efficient batched
sum-check's inherent working set (S526).

**What would falsify this (S527):** `stream_witnesses` ON ≠ OFF on any of
claimed/comm/comm_bt/comm_ub/accepted (selftests would fail); the A/B
peak-RSS gap NOT growing ~4×/Δn=2 (would mean the removed term is not the
Θ(x) list); or `full` peak ≥ `full_nostream` peak (streaming not helping).
Reproduce: `python3 large_x_benchmark.py --stream-probe 16,18,20,22`;
`python3 compressed_layer.py --selftest` (case 25); `python3
large_x_benchmark.py --selftest` (case 6).

## Honest scope
- This is the open-item-1 Õ(x) prover working the √x state; it is NOT a
  polylog computation of π(x) (the goal stays blocked). The deliverable
  is the *verification* artifact: a verified exact π(x) at x≈1.7×10⁷
  behind an Õ(√x) verifier / Õ(x^{1/4}) large-table opens.
- `commit_base=True` makes the large-table verifier term computational
  (CRH/RO via sha256, S508); the rest of the chain is unconditional.
- The headline is one seed; the selftest exercises the cheat panel over
  two fields and both fast/object paths at small n.

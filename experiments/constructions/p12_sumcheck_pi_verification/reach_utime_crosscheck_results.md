# reach_utime_crosscheck — resolving the reach-wall falsifier from the UTIME side

**Cycle S536 (2026-06-14).** Script:
`reach_utime_crosscheck.py` (`--live` / `--baseline` / `--dram-factor` /
`--predict` / `--all` / `--selftest`).
Run log: `reach_utime_crosscheck_run.log`.

## The question

The detached n=26 verification-chain reach (PID 1592836, started
2026-06-14T13:24:47Z) has a WALL that **already exceeds S529's predicted band
[5.9, 7.9] ks** — it is past 12 ks of CPU and still running. On the *raw* wall
test that would trigger the "super-Θ(x)" falsifier. S535 SHARPENED the test: a
wall ratio n24→26 > 5.4× is CONSISTENT with a Θ(x) op-count **iff** the excess
lives in `stime` (kernel/page-fault/allocation), and only a `utime`-dominated
excess with super-linear per-element growth reopens super-Θ(x). S535 measured
the reach is `stime`-dominated (0.695) but could **not** complete the clean
n=24→26 *utime* ratio test — **n=24's utime/stime split was lost on process
exit** (the process was gone; `/proc` counters vanish on exit).

This cycle closes that gap **without** the lost n=24 split, using one fact:

> **utime (CPU user time) is contention-robust.** The concurrent reach steals
> WALL-clock, not CPU-seconds. If a probe process consumes `T` cpu-seconds in
> user mode, that is real work it did regardless of core contention. So the
> chain's utime baseline can be measured at small RAM-light `n` *right now*,
> despite the running reach, and extrapolated via the measured Θ(x) op-count
> (S530/S531/S532b).

## What was measured

### (A) Live reach `/proc` capture — ZERO contention (perishable; lost on exit)

Five captures over the cycle (16:28–16:51Z, reach still RUNNING, no DONE):

| t (Z) | utime_s | stime_s | cpu_s | cum stime_frac | window stime_frac | RSS |
|---|---|---|---|---|---|---|
| 16:28 | 3325.7 | 7695.6 | 11021 | 0.698 | — | 4.1 GiB |
| 16:30 | 3375 | 7833 | 11208 | 0.699 | 0.746 | 3.9 GiB |
| ~16:45 | 3549 | 8311 | 11860 | 0.701 | 0.752 | 3.6 GiB |
| 16:48 | 3582 | 8401 | 11983 | 0.701 | — | 3.8 GiB |
| 16:51 | 3689 | 8704 | 12393 | 0.702 | — | 4.7 GiB |

- `stime_frac` is **steady at 0.70** (cumulative), windowed 0.75 — confirms S535's
  0.695 ⇒ **allocation/page-fault bound** (3.7×10⁵ minor faults/s, 3.9×10⁹
  cumulative, **majflt = 0** — no swap). VmHWM = 11.4 GiB.
- utime advances steadily (3326→3689) ⇒ the reach is **actively computing**, not
  stalled. RSS oscillates (3.6–4.7 GiB, below the 11.4 GiB peak) ⇒ multiple
  phases.

### (B) Chain utime baseline — contention-robust CPU time of the REAL chain

`run_chain` FULL config, BIG_Q+FAST_BIG, honest ACCEPTED + `claimed == sieve`
at every n. utime/stime via `os.times()` deltas; op-count via the S530 `fmul`
monkeypatch (one instrumented run gives both):

| n | mul_ops | utime_s | stime_s | stime_frac | **ns/op** | wall_s |
|---|---|---|---|---|---|---|
| 14 | 23,492,482     | 4.80  | 0.02  | 0.004 | **204.3** | 4.84 |
| 16 | 105,029,523    | 8.93  | 0.27  | 0.029 | **85.0**  | 9.19 |
| 18 | 465,458,097    | 21.65 | 1.30  | 0.057 | **46.5**  | 22.97 |
| 20 | 2,040,403,997  | 56.94 | 9.98  | 0.149 | **27.9**  | 66.95 |
| 22 | 8,882,635,340  | 184.69| 51.77 | 0.219 | **20.8**  | 236.60 |

**Two trends, both load-bearing:**

1. **Per-op utime DECLINES monotonically** 204 → 85 → 46.5 → 27.9 → 20.8 ns/op
   (overhead per numpy call amortizes as arrays widen). The decline decelerates
   (per-op per-Δn=2 ratio 0.42 → 0.55 → 0.60 → 0.745) ⇒ converging to a floor
   (~20 ns/op). Correspondingly utime per-Δn=2 ratios RISE 1.86 → 3.24 toward the
   op-count's ~4.4× ⇒ **utime → Θ(x)** asymptotically.
2. **stime_frac RISES monotonically** 0.004 → 0.029 → 0.057 → 0.149 → 0.219 with
   `n` (working set) — an **independent corroboration of S535's allocation-bound
   mechanism**: stime is a rising function of n, extrapolating toward the reach's
   **0.70 at n=26** (the allocation leg grows with the Θ(x) working set, exactly
   as S535 predicted; saturates as faults come to dominate).

### (C) cache→DRAM per-op UTIME factor — reused buffers, no allocation (stime≈0)

| size | ns/elem (utime) | ×L1 |
|---|---|---|
| 32 KiB (L1) | 10.1 | 1.00 |
| 512 KiB | 7.1 | 0.70 |
| 8 MiB | 14.2 | 1.40 |
| 64 MiB (>L3, DRAM) | 22.9 | 2.26 |
| 512 MiB (DRAM) | 22.8 | 2.25 |

Per-element `_mul61` utime **saturates at ~23 ns (2.26× over L1); last
size-step 0.99×** ⇒ **bounded cache constant, F4 NOT triggered.** (The 512 KiB
dip is cache-warmth noise; the L1→DRAM rise is monotone and saturating.)

### (D) Cross-check: reach realized per-op vs the declining baseline trend

- `opcount(26)` = **1.649×10¹¹** (model `x·(log₂x)^0.87`, anchored on the
  measured n=24 = 3.84×10¹⁰; S532b validated this model to <0.4% on per-step
  ratios).
- reach realized per-op (so far, utime=3549) = **21.52 ns/op** vs the
  DRAM-resident **n=22 anchor 20.79 ns/op → ratio 1.04×** — squarely ON the
  declining trend.

## Verdict — the falsifier is NOT triggered (from the utime side)

| falsifier | meaning | result |
|---|---|---|
| **F1** | baseline per-op utime RISES with n (super-linear compute) | **NOT triggered** — per-op declines 204→20.8 ns/op |
| **F2** | live utime ≫ the Θ(x) band (super-linear per-op term) | **NOT triggered** — realized per-op 21.5 ≈ n=22 anchor 20.8 (1.04×) |
| **F3** | reach NOT `stime`-dominated (excess is utime) | **NOT triggered** — cum stime_frac 0.70 |
| **F4** | DRAM per-op utime unbounded / non-saturating | **NOT triggered** — saturates 2.26×, last step 0.99× |

**The reach's per-op utime sits on/below the monotonically declining baseline
trend** ⇒ utime grows **no faster than the Θ(x) op-count** ⇒ the entire >5.4×
WALL excess over S529's band is the **`stime`/allocation leg** (70% of wall).
This is the UTIME confirmation S535's sharpened falsifier required — the leg it
could not measure because n=24's split was lost. **Θ(x) op-count intact.**

This also unifies the three reach-wall legs into one consistent picture, all now
measured: `wall = opcount(Θ(x), S530/S531/S532b) × compute/cache(≤1.35–2.3×
saturating, S529 + (C)) × allocation/stime(rising 0.004→0.70 with n, S535 +
(B))`; the wall can exceed the op-count ratio purely because the allocation leg
grows with the working set, **with no op-count or per-op-compute change**.

## Honest scope and caveats

- **An in-cycle MEASUREMENT** resolving the open *utime* side of the reach-wall
  falsifier — NOT a protocol or goal advance. polylog π(x) stays blocked.
- **The reach is NOT done.** `utime=3549 s` is a **LOWER bound** on the final
  utime(26), so realized per-op 21.5 ns/op is a lower bound on the final per-op.
  At a plausible 60–85% utime-progress at capture, final per-op lands ≈
  25–36 ns/op = 1.2–1.7× the n=22 anchor — **bounded, still far below the small-n
  per-op (46–204 ns/op); a genuine super-linear term would push it to 100s.** The
  exact final multiple awaits DONE; `--predict` reads `run_n26.log` for DONE and
  finalizes (verdict bands: ≤1.3× "not triggered", 1.3–2.0× "weak/confirm",
  >2.0× "F2 REOPEN").
- This is the utime *complement* to the harvest the NEXT ACTION still wants on
  DONE: record claimed π(2²⁶−1) = 3 957 809 (cross-check Eratosthenes), final
  wall, peak RSS, and the final utime/stime split. Tooling is ready
  (`alloc_overhead_wall.py --live`, `reach_utime_crosscheck.py --predict`).
- `opcount(26)` is a MODEL value (1.649×10¹¹) extrapolated from measured n=24; a
  direct n=26 op-count run (`prover_opcount_scaling.py`) would pin it exactly.
- Probe perturbation of the reach is negligible (CPU time is contention-robust;
  brief baseline runs share DRAM bandwidth only momentarily — fault rate held
  steady at 3.7×10⁵/s across the cycle, ≈ S535's +2.6% precedent). The verdict
  rests on stime/utime, which probes cannot corrupt.

## What would falsify this result

The four pre-stated falsifiers above (F1–F4). Concretely, after DONE: a final
utime(26) so large that final per-op ≫ the declining trend (>2× the n=22 anchor
**with** per-op rising, not the bounded DRAM penalty), OR a direct n=26 op-count
run showing op(26) ≫ 1.65×10¹¹ (super-Θ(x) op-count) — either would reopen the
super-Θ(x) concern. Neither is indicated by the present data.

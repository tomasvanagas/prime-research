# alloc_overhead_wall.py — the THIRD leg of the reach wall (S535)

**Date:** 2026-06-14 · **Cycle:** S535 · **Status:** measurement (refines S529;
NOT a goal advance — polylog π(x) stays blocked)

## Question

The large-x verification reach (`large_x_benchmark.py`, PROGRAM.md item 5) wall has
been modelled (S529, `mem_hierarchy_wall.py`) as **two** legs:

1. **op-count** Θ(x)·polylog — per-Δn=2 ratio 4.0–4.4× DOWN (confirmed S530/S531/S532b to n=24);
2. **compute/cache** — the L1→DRAM `_mul61` ns/elem *bowl*: one-time ~7× L3-crossover step, per-Δn growth ≤1.35×.

⇒ S529 predicts the n24→26 wall ratio ∈ **[4.0, 5.4]×**, with "well above ~5.5× ⇒ reopen super-Θ(x)".

But S529 leg (2) was measured as **steady-state `_mul61` over PRE-ALLOCATED arrays** — it
never measured the cost of **allocating and first-touching** the large arrays. This cycle
observed the *live* n=26 reach and found it is **kernel/allocation-bound, not compute-bound**.

## What was measured

Standalone, contention-safe (`alloc_overhead_wall.py`, 17 selftests incl. a field-offset
cross-check):

- **`--live`** — read the running reach's `/proc/PID/stat` counters (ZERO contention, pure reads):
  `stime`/`utime` split, minor/major fault rate, cumulative faults, RSS. Auto-detects the
  `large_x_benchmark` PID.
- **`--probe`** — a self-limiting micro-benchmark separating, for one elementwise uint64
  modular-mul step, **"fresh allocation + first-touch each rep"** (malloc/mmap + minor page
  faults → kernel/`stime`) from the **same op on a reused pre-faulted buffer** (`out=`, pure
  compute → `utime` — S529's regime). Per-phase kernel/user attribution via `/proc/self/stat`.
- **`--all`** — live-before + probe + live-after (contention check) + the three-leg synthesis.

## Results (run log `alloc_overhead_run.log`, 2026-06-14T~16:08Z, n=26 reach at ~2.8 h)

### (A) Live n=26 reach — in-situ, zero contention

| quantity | value |
|---|---|
| state / threads / RSS | R / 32 / 4.0 GiB (VmHWM peaked 11.99 GiB earlier) |
| minor-fault rate | **3.28×10⁵ /s** (fault-traffic 1.3 GiB/s) |
| major-fault rate | **0** (no swap; not memory-ceilinged) |
| window `stime` fraction | **0.725** |
| **cumulative `stime` fraction** | **0.695** (utime 3079 s, stime 7012 s) |
| cumulative minor faults | **3.087×10⁹** |

**≥69.5% of the reach's CPU is kernel page-fault handling, not arithmetic.** The minimum
*unavoidable* faults = peak working set / page ≈ 12 GiB / 4 KiB ≈ 3.1×10⁶; the reach has
incurred **3.1×10⁹** — i.e. **≈1000× the working-set minimum**, so virtually every fault is a
*re-allocation* (the same logical Θ(x) buffers freed and re-faulted across the K≈π(√x)≈1028 layers).

### (B) Allocation-overhead micro-probe (uint64 modular-mul step)

| array size | reused (ns/el) | fresh (ns/el) | fresh/reused | alloc overhead (ns/el) | reused stime% | fresh stime% |
|---|---|---|---|---|---|---|
| 64 KiB (≤L2) | 3.04 | 3.18 | 1.04 | 0.14 | 0.00 | — |
| 256 KiB | 2.92 | 10.37 | 3.55 | 7.45 | 0.00 | 0.75 |
| 1 MiB | 2.94 | 11.65 | 3.96 | 8.71 | 0.00 | 0.75 |
| 16 MiB | 3.08 | 14.48 | 4.70 | 11.39 | 0.00 | 0.83 |
| 128 MiB | 3.64 | 14.34 | 3.94 | 10.70 | 0.00 | 0.69 |

- **Reused-buffer op** = flat **~3 ns/elem**, **utime-dominated** (stime%≈0): this *is* S529's
  compute bowl (leg 2).
- **Fresh-allocation op** jumps to **~14 ns/elem above 256 KiB** (the L2 boundary — beyond it
  first-touch must fault in real pages), **kernel/stime-dominated (0.69–0.83)**, adding
  **~10 ns/elem** = a **~4× penalty** over the reused buffer.
- The probe's fresh stime fraction (0.69–0.83) **matches the live reach's 0.695** — the probe
  reproduces the reach's regime in miniature.

### Contention check

Live reach minor-fault rate **3.277×10⁵ → 3.361×10⁵ /s across the probe = +2.6%** (cumulative
stime fraction unchanged at 0.695). **The probe did not materially contend the reach** (well
within the ~3.5% run-to-run wall noise S524 noted); the clean n=26 falsifier point is intact.

## The three-leg model

```
wall(n) = opcount(n)  ×  compute_cache_factor(n)  ×  alloc_factor(n)
          └ leg 1 ┘       └────── leg 2 ──────┘       └── leg 3 (NEW) ──┘
          Θ(x)·polylog    L1→DRAM bowl, one-time       malloc/mmap + minor
          4.0–4.4×↓       ~7× L3 step, ≤1.35×/Δn       page-faults → KERNEL stime
          (S530/1/2b)     (S529)                       (this cycle; S529 omitted it)
```

Leg 3 grows with the **working set** (Θ(x) cubes ⇒ Θ(x)·polylog faulted bytes per discharge),
so its per-Δn factor can **exceed** the op-count ratio and push the *wall* above S529's [4.0,5.4]×
band — **without any op-count change**.

## Sharpened falsifier for the harvest cycle

S529's test ("n24→26 > ~5.5× ⇒ reopen super-Θ(x)") conflates legs 1 and 3. Corrected:

> A wall ratio n24→26 **> 5.4× is CONSISTENT with Θ(x) op-count** iff the n=26 reach is
> `stime`-dominated (leg 3, allocation). The discriminating measurement is **`stime`/`utime`**,
> not the wall alone. The reach is measured **`stime`-dominated (0.695)**, so a high wall ratio is
> page-fault overhead, not op-count. **Only** a wall excess that (i) survives subtracting allocation
> `stime` AND (ii) is **`utime`-dominated** with super-linear per-element growth would reopen
> super-Θ(x). Op-count is independently Θ(x) to n=24 (S530/S531/S532b).

Because the process's `stime`/`utime` is lost when it exits, the harvest cycle should sample
`/proc/<reach-pid>/stat` (or this script's `--live`) **before** the n=26 process terminates. The
in-situ snapshot at ~2.8 h is recorded above (utime 3079 s, stime 7012 s, stime-frac 0.695).

## The removable headroom (a lever, NOT implemented this cycle)

The 69.5% kernel fraction is **headroom**, not complexity: ≈1000× of the faults are
re-allocations of the same Θ(x) buffers. Reusing them (pre-allocate once, `np.multiply(...,
out=buf)` throughout the chain) would fault each page ~once (≈3×10⁶ faults), collapsing `stime`
from ~7000 s to ~7 s and cutting the wall by a factor of **≈1/(1−0.695) ≈ 3.3×** (a *ceiling*
estimate — some faults are unavoidable). This is a **constant-factor implementation win**
(the analog of S527 list-streaming and S525 scatter-fold), **NOT an op-count / complexity
change** (op-count stays Θ(x)). It is **deliberately NOT implemented here**: it would rewrite
every `a*b%q` across the landed chain and break the n=24/n=26 artifacts' bit-identical
reproducibility. Filed as the next free reach lever; if pursued, it must be opt-in (default off)
like the S496 structured prover, with a bit-identical A/B selftest.

## Honest scope / what would falsify this

- **Scope:** an in-cycle MEASUREMENT refining the S529 reach-wall model (a *known* wall — the
  same Θ(x)·polylog prover, now decomposed with the previously-missing allocation leg). It does
  **not** advance the goal; polylog π(x) stays blocked, op-count stays Θ(x), the cert/info/
  field-soundness √x walls (S509/S511/S534) are untouched.
- **Pre-stated falsifiers (status):**
  - **F1** reused-buffer op NOT utime-dominated → would break leg-2 attribution. *Un-triggered* (reused stime%≈0).
  - **F2** fresh allocation NOT stime-dominated at large N → leg 3 not page-faults. *Un-triggered* (fresh stime% 0.69–0.83).
  - **F3** fresh/reused per-elem ratio does NOT grow with N then saturate → no allocation regime. *Un-triggered* (1.04 at L2 → ~4× saturated above 256 KiB).
  - **F4** the live reach NOT stime-dominated → premise false, reach is compute-bound, and a
    >5.4× wall WOULD reopen super-Θ(x). *Un-triggered* (reach stime-frac 0.695).
- **What would falsify the headline:** if a later, cleaner reach (or the same one re-measured
  with the fast `_mul61` path fully engaged and buffers reused) showed `stime`/`utime` ≪ 1 at
  n≥26 with the wall still > 5.4×/Δn=2 and `utime`-per-element growing super-linearly — then leg
  3 is small and the wall would be a genuine op-count concern. Not observed.
- **Caveat:** the full quantitative n24→n26 *utime* (compute-only) ratio cannot be pinned because
  n=24's `stime`/`utime` split was not logged (that process exited). The decomposition here is
  the live n=26 split + the micro-probe per-element penalty; future reaches should log
  `getrusage` (ru_utime/ru_stime) so the wall is self-decomposing.
- **Determinism:** the live numbers are a process snapshot (not reproducible bit-for-bit); the
  probe ns/elem are timing (machine-dependent), but the *shape* (flat reused ~3 ns/elem,
  utime-dominated; fresh ~14 ns/elem above L2, stime-dominated; ~4× ratio) and the live
  stime-fraction (~0.7) are the robust signal. Selftest is deterministic (value equality of
  `_op_into`==`_op_fresh`, field-offset cross-check vs `os.times()`).

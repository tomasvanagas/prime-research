# mem_hierarchy_wall.py — results (S529)

**Question.** PROGRAM.md item 5 carries a filed falsifier on the large-x
verification reach: the prover wall ACCELERATES with x — the n=22→24 step ran at
**5.80×** per Δn=2 versus the ideal Θ(x) ratio of **4×**. S524 first read this as
possibly super-Θ(x) (a genuine complexity term); S525 re-attributed it to "a
memory-hierarchy constant (Θ(x) arithmetic over a cache-busting Θ(x) working
set)" — an *implementation* cost. That attribution was made from peak-RSS probes
only, never from a direct measurement of the per-element field-op wall as the
working set crosses the cache levels. This experiment supplies that measurement,
in-cycle, without waiting for the hours-long macro run.

**Host.** AMD EPYC 7313, L1d 32 KiB / L2 512 KiB / L3 32 MiB (per the topology
read by the script), 122 GiB RAM. Kernels: the BIG_Q=2⁶¹−1 fast Mersenne path
(`_mul61` ≈ 24 numpy ops/multiply, `_sum61` block-fold sum) from
`compressed_prover_mult_trace.py`, which is what the FAST_BIG reach run uses.

---

## 1. Kernel cache-step sweep (`--bench-kernels`)

`_mul61` (Mersenne mulmod) per-element wall vs array footprint, one operand:

| log₂ n | footprint | regime | `_mul61` ns/elem | step× | `_sum61` ns/elem | np-mul ns/elem |
|---:|---:|:--:|---:|---:|---:|---:|
| 10 | 8 KiB   | L1d  | 27.59 | —    | 50.87 | 4.73 |
| 11 | 16 KiB  | L1d  | 17.39 | 0.63 | 35.36 | 3.70 |
| 12 | 32 KiB  | L1d  | **12.79** | 0.74 | 21.47 | 3.26 |
| 13 | 64 KiB  | L2   | 71.89 | 5.62 | 16.24 | 3.01 |
| 14 | 128 KiB | L2   | 73.81 | 1.03 | 12.05 | 2.92 |
| 15 | 256 KiB | L2   | 36.15 | 0.49 | 10.50 | 2.92 |
| 16 | 512 KiB | L2   | 59.79 | 1.65 | 9.48  | 2.88 |
| 17 | 1 MiB   | L3   | 69.85 | 1.17 | 9.00  | 10.03 |
| 18 | 2 MiB   | L3   | 59.82 | 0.86 | 6.98  | 11.89 |
| 19 | 4 MiB   | L3   | 87.17 | 1.46 | 8.57  | 12.06 |
| 20 | 8 MiB   | L3   | 81.93 | 0.94 | 8.50  | 13.06 |
| 21 | 16 MiB  | L3   | 75.22 | 0.92 | 6.86  | 11.58 |
| 22 | 32 MiB  | L3   | 104.79| 1.39 | 7.52  | 11.47 |
| 23 | 64 MiB  | DRAM | 95.48 | 0.91 | 12.12 | 11.29 |
| 24 | 128 MiB | DRAM | 93.62 | 0.98 | 14.24 | 11.60 |
| 25 | 256 MiB | DRAM | 91.55 | 0.98 | 16.28 | 14.49 |
| 26 | 512 MiB | DRAM | 98.53 | 1.08 | 16.94 | 12.61 |
| 27 | 1 GiB   | DRAM | 108.98| 1.11 | 15.87 | 11.62 |

**Shape.** `_mul61` ns/elem is bowl-shaped: ~28 ns at 8 KiB (per-call dispatch
overhead amortised over too few elements), a **minimum of 12.8 ns at 32 KiB
(L1-resident + dispatch amortised)**, then a sharp step up at the L2 boundary, and
a **flat ~90–110 ns plateau once the array exceeds L3 (≥64 MiB, DRAM)**. The
plateau is a **~7–8× cache penalty over the L1 minimum**. The pure-bandwidth
control (`np-mul`, 1 op/elem) shows the same shape one level down (flat ~3 ns in
L1/L2, stepping to ~11–14 ns in L3/DRAM) — confirming the plateau is
memory-bandwidth, not ALU.

**The load-bearing fact:** the DRAM plateau is **flat**, not climbing. From
64 MiB to 1 GiB `_mul61` moves only 95→109 ns (≈1.15× over ~4 octaves ≈
1.04×/octave; the `_sum61`/np-mul controls are flat to within ±20%). So once a
working set exceeds L3 (≈32 MiB, which every reach n≥18 already does), the
per-element wall is essentially **constant** — the cache penalty is a **one-time
step at the L3 crossover, NOT a per-n growing term.**

---

## 2. Reach-wall decomposition (`--predict`)

Every reach n≥18 has its Θ(x) held cube already in DRAM (peak RSS, FULL streamed
config, S527/S525): n=20 185 MiB, n=22 677 MiB, n=24 ~2.6 GiB — all DRAM. The
per-element `_mul61` wall at those footprints (read off the sweep) grows only:

| step | per-elem wall growth | observed reach ratio | ideal Θ(x) |
|:--|---:|---:|---:|
| n=20→22 | 1.076× | **3.521×** | 4× |
| n=22→24 | 1.106× | **5.804×** | 4× |

Model `wall_ratio = 4× (op-count Θ(x)) × per-elem memory growth`:

- n20→22 predicted **4.31×**, observed **3.52×** (below)
- n22→24 predicted **4.42×**, observed **5.80×** (above)

The measured per-element memory growth (~1.08–1.11×) is **far too mild to turn 4×
into 5.80×.** And the two single-steps **anti-correlate** (3.52× low, 5.80× high):
their **geometric mean over n=20→24 is 4.52×** ≈ 4× (op-count) × 1.13× (the
plateauing memory factor). 

**Conclusion (refines both prior readings):**

- **vs S524** (5.80× as possible super-Θ(x)): the geo-mean 4.52× ≈ Θ(x)·polylog is
  consistent with the Õ(x) prover (open item 1). The single-step 5.80× is **3-point
  noise**, not evidence of a super-Θ(x) op-count term.
- **vs S525** (super-4× = memory-hierarchy constant): the memory hierarchy **does**
  impose the large ~7–8× *absolute* slowdown (DRAM vs L1-cached), but as a **one-time
  step at the L3 crossover (n≈18)**; its **per-Δn growth (~1.1×) is too small to be
  the 5.80× source.** The memory hierarchy explains the high absolute wall constant,
  not the per-step super-4× ratio.

The honest synthesis: the reach wall is **Θ(x)·polylog op-count × a one-time
DRAM-plateau constant (~7×), with single-step ratios noisy on 3 points.** No
super-Θ(x) term is needed or supported.

---

## 3. Multi-GiB bandwidth/TLB probe (`--bandwidth-probe 30`)

The §1 `_mul61` sweep capped at 1 GiB; the n=26 held cube is ~10 GiB. The light
controls (`_sum61`, np-mul; no `_mul61`, which would need ~80 GiB of transients at
8 GiB) push into that regime:

| log₂ n | footprint | regime | `_sum61` ns/elem | np-mul ns/elem |
|---:|---:|:--:|---:|---:|
| 23 | 64 MiB  | DRAM | 12.26 | 12.91 |
| 24 | 128 MiB | DRAM | 14.61 | 11.43 |
| 25 | 256 MiB | DRAM | 15.50 | 11.99 |
| 26 | 512 MiB | DRAM | 18.62 | 11.97 |
| 27 | 1 GiB   | DRAM | 16.05 | 11.27 |
| 28 | 2 GiB   | DRAM | 16.08 | 12.43 |
| 29 | 4 GiB   | DRAM | 17.25 | 13.55 |
| 30 | 8 GiB   | DRAM | 22.58 | 16.17 |

**Finding.** The pure-bandwidth control (np-mul) is **nearly flat from 64 MiB to
4 GiB (~11–14 ns)**, ticking up only at 8 GiB (16.2 ns, ≈1.35× the plateau).
`_sum61` (strided fold/reduceat) climbs more — 12→23 ns over 64 MiB→8 GiB
(≈1.85×) — the expected TLB/page-walk pressure (8 GiB = 2M × 4 KiB pages, far past
the TLB reach). So there **is** a slow memory degradation past ~1–2 GiB, but it is
**~1.05–1.09×/octave** — a slowly-saturating effect, not a cliff. Over a Δn=2 step
(working set ×4 = 2 octaves) this is ~**1.1–1.18×**, still far below the 1.45×
needed to read 4× as 5.80×.

---

## 4. Prediction for the detached reach (the falsifier test)

PENDING the n=24→26 macro run (running detached via `run_reach_detached.sh`,
survives the cycle boundary; ground truth π(2²⁵−1)=2 063 689, π(2²⁶−1)=3 957 809):

- **Predicted n=24→26 wall ratio ~4.4–5.4×** = op-count 4× × the §3 multi-GiB
  memory factor (~1.1–1.35× as the working set goes 2.6 GiB→10 GiB, including the
  TLB rise), i.e. **n=26 wall ~6.5–7.9 ks**. (The script's `--predict` caps its
  lookup at 1 GiB and so reports the plateau-only floor 4.0×; §3 lifts the honest
  upper end via the multi-GiB controls.)
- **Falsifier:** a measured n=24→26 ratio **well above ~5.5×** — *with* the
  per-element wall shown flat-to-mild here — would refute the "implementation-cost-
  only" reading and re-open the S524 super-Θ(x) concern. A ratio in ~[4.0, 5.4]×
  confirms Õ(x) op-count × a one-time DRAM step + slow TLB drift. Either way the
  per-step value is noisy on few points; the **trend** (geo-mean across n≥20
  staying near ~4.5×) is the real test, not a single step.

---

## Reproduce

```
python3 mem_hierarchy_wall.py --selftest          # 7 kernel/logic groups
python3 mem_hierarchy_wall.py --all               # sweep to 1 GiB + decomposition
python3 mem_hierarchy_wall.py --bandwidth-probe 30 # light sweep to 8 GiB
```

**Scope/caveats (honest).** (i) The microbenchmark times the kernels in isolation;
the real chain interleaves many narrow per-layer arrays (dispatch-bound, the high-
ns/elem left branch — S526) with the big stacked-cube fmuls (the DRAM plateau), so
the per-element rate is an envelope, not the exact chain mix. (ii) ns/elem at any
single size is noisy (±20%) because the host is shared (a reach run + this probe
ran concurrently); the *shape* (bowl + plateau) and the *order of magnitude* of the
L3-crossover step are robust, the third-decimal values are not. (iii) The decomposition
uses 3 measured reach walls (S524) — the conclusion that 5.80× is noise rests on the
anti-correlation of 2 single-steps; the detached n=25/26 points (§4) are the
independent check.

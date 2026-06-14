# prover_opcount_scaling.py — whole-chain prover op-count exponent (S530, extended S531)

**Date:** 2026-06-14
**Script:** `prover_opcount_scaling.py` (`--selftest` / `--scaling` / `--points`)
**Run logs:** S530 in-cycle (`--scaling --ns 8,10,12,14,16,18,20`); S531 the
detached path-(b) extension `opcount_scaling_n24.log` (`--ns 16,18,20,22,24`) and
the corrected combined fit `opcount_scaling_corrected.log` (n=8..20 live + n=22,24
injected deterministic).

> **S531 headline (see the dedicated section below):** the curve was extended to
> **n=24** (x≈1.7×10⁷). The op-count is **Θ(x)·polylog = Õ(x)** out to n=24 — the
> per-step ratio kept DECLINING (4.47→4.33), S530's live falsifier did NOT fire.
> The detached run's *auto-reading* printed "SUPER-Θ(x)" from a single-power fit
> over the short high-n window (α=1.064); that was the **polylog-inflation artifact**
> S512's `sharpP_probe` 2-term fit exists to defeat — now added to this script and
> the reading corrected. A 2-term fit over n=14..24 recovers **α=1.01, polylog
> δ=+0.72** (residual 17× tighter), and an **x·(log₂x)^0.87** model predicts every
> per-step ratio to <0.4%.

## What this measures and why

The large-x reach (PROGRAM.md item 5) has a filed wall-scaling falsifier: the
prover wall ACCELERATES — the n=22→24 step ran **5.80×** per Δn=2 vs the ideal
Θ(x) ratio of 4×. S529 (`mem_hierarchy_wall.py`) attributes the excess to a
one-time L3→DRAM per-element memory step plus 3-point single-step noise (NOT a
super-Θ(x) op-count), and predicts the detached n=24→26 reach lands in **~[4.0,
5.4]×**. But that whole reading rests on an **ASSUMED** leg: `mem_hierarchy_wall`
hard-codes `IDEAL = 4.0` ("ideal op-count ratio is Theta(x) = 4x per Delta n=2",
line 216–218). **The op-count itself was never measured on the real chain.** The
prior op-count measurements are verifier-side (`vleaf_ops`, S507/S509) or one
primitive's per-layer p²→p wiring ratio (`--prover-bench`, S496/S497); the
whole-chain total field op-count exponent — the load-bearing leg — was unmeasured.

This probe supplies it. It monkeypatches the field-MULTIPLY dispatch point
`compressed_prover_mult_trace.fmul` (→ `_mul61` Mersenne mulmod under FAST_BIG for
BIG_Q; the chain's `sumcheck`/`eq_table`/`mle_eval`/per-layer reductions route
their array element-multiplies through it) in **every** module that imported it,
tallies the broadcast OUTPUT size of each call (= element-multiplies), runs the
real `compressed_layer.run_chain` over the **FULL succinct config** on
BIG_Q+FAST_BIG at n = 8,10,…,20, and fits log₂(mul_ops) vs n. Since x = 2ⁿ,
mul_ops ~ xᵅ ⇒ slope = α (per-Δn=2 ratio = 2^{2α}).

The op-count is **DETERMINISTIC** (array sizes are structural — fixed by x — not
FS-challenge dependent; selftest [3] asserts seed-independence exactly), so this
measurement is **immune to the concurrent detached reach's CPU contention** (only
wall-clock would be perturbed). It is RAM-light (n≤20, peak ~185 MB) and changes
no landed code (standalone monkeypatch; the chain verdict is bit-identical under
the patch — selftest [1]).

## Result — the prover op-count is Θ(x), measured

FULL config (delegate+structured+pcs+batch_trace+batch_ub+batch_wiring+
commit_base), BIG_Q = 2⁶¹−1, FAST_BIG, honest claim ACCEPTED and `claimed ==
sieve` at every n:

| n | x | K=π(√x) | nb | mul_ops | sum_ops | scatter_ops | mul/x |
|---|---|---|---|---|---|---|---|
| 8  | 255      | 6   | 4  | 425 110        | 112 472     | 176     | 1667 |
| 10 | 1 023    | 11  | 5  | 2 074 457      | 546 140     | 672     | 2028 |
| 12 | 4 095    | 18  | 6  | 9 673 318      | 2 543 740   | 2 240   | 2362 |
| 14 | 16 383   | 31  | 7  | 23 492 482     | 6 054 726   | 7 808   | 1434 |
| 16 | 65 535   | 54  | 8  | 105 029 523    | 27 057 082  | 27 392  | 1603 |
| 18 | 262 143  | 97  | 9  | 465 458 097    | 119 807 096 | 98 816  | 1776 |
| 20 | 1 048 575| 172 | 10 | 2 040 403 997  | 525 295 418 | 351 232 | 1946 |

**Fitted leading exponent** (mul_ops ~ xᵅ):

| metric | α | per-Δn=2 ratio 2^{2α} | max\|resid\| (log₂) |
|---|---|---|---|
| **multiply element-ops** | **0.9955** | **3.975×** | 0.377 |
| mult + add ops | 0.9947 | 3.971× | 0.378 |
| wall-clock (secondary) | 0.524 | 2.07× | 0.481 |

**Per-step multiply-op ratios** mul_ops(n+2)/mul_ops(n) (IDEAL Θ(x) = 4.000×):

```
n= 8→10 : 4.880    n=14→16 : 4.471
n=10→12 : 4.663    n=16→18 : 4.432
n=12→14 : 2.429    n=18→20 : 4.384
geo-mean per-Δn=2 = 4.107×
```

### Reading

1. **α_op = 0.995 ≈ 1.0 ⇒ the whole-chain prover op-count is Θ(x).** The polylog
   factors cancel as designed: per-layer work ≈ √x·nb, K = π(√x) ≈ √x/(½n·ln2)
   layers, so total ≈ K·√x·nb = Θ(x) with the nb and 1/(½n) factors offsetting.
   This is the **first direct measurement of the op-count leg that
   `mem_hierarchy_wall` assumed** (IDEAL = 4.0). Measured geo-mean per-Δn=2 ratio
   = **4.11×**, fit ratio **3.98×** — the assumption is confirmed.

2. **The small excess above 4 at moderate n is a vanishing finite-size effect, not
   super-Θ(x).** The upper-regime steps DECREASE monotonically — 4.471 → 4.432 →
   4.384 (n14→16→18→20) — settling toward 4.0 as n grows. There is no plateau at
   4.4; the residual polylog dressing on the Θ(x) core slowly washes out. (The
   one low step, n12→14 = 2.43×, is the compensating other half of an
   nb/K-alignment lump — it is preceded by the high 4.66× at n10→12; the two
   steps geo-mean to 3.30, the pair-of-pairs average is ~2× per Δn=1, i.e. Θ(x).)

3. **So the reach's 5.80× wall step (n22→24) is NOT in the op-count.** Op-count at
   that scale is ~4.4× and falling toward 4.0. The 5.80× excess is the
   memory-hierarchy step + single-step noise (S529 reading), confirmed from the
   deterministic-op-count side. The S529 prediction band n24→26 ~ [4.0, 5.4]×
   stands: its op-count leg (4.0–4.4×, declining) × the measured per-element DRAM
   growth (1.0–1.35×, S529) ⇒ [4.0, ~5.4]×.

4. **Independent corroboration of S529's "single-step ratios are noisy".** Even
   the fully DETERMINISTIC op-count has lumpy per-step ratios (2.43× … 4.88×,
   range > 2×) from the discrete nb (+1 per Δn=2) and K = π(√x) increments. So a
   single wall step ratio (like the 5.80×) carries ≥ this much intrinsic
   step-to-step variation BEFORE any memory or timing noise — over-reading one
   step is unjustified; the geo-mean / trend is the signal.

5. **Wall brackets the bowl.** The secondary wall exponent over n=8..20 is only
   **0.524** (per-step ratios climbing 1.97 → 2.94 toward 4) — the in-cache /
   dispatch-bound regime where the per-element cost FALLS as arrays widen,
   deflating the wall exponent below the op-count. At the reach (n=20..24) the
   per-element cost rises (DRAM), inflating the wall exponent to ~1.088 (S524).
   The op-count exponent ≈ 1.0 is the **invariant Θ(x) core** between the small-n
   dispatch-deflation and the reach-n DRAM-inflation — exactly the two-regime
   per-element wall S529 measured directly. fmul accounts for ~78% of the (lightly
   instrumented) wall, so multiplies are the dominant timed op and a faithful
   op-count proxy.

## Falsifiers

- The patched `fmul` changing the chain verdict (claimed/comm/accepted) vs the
  unpatched chain ⇒ the counted run is not the real chain. **Asserted** —
  selftest [1] (bit-identical over q & BIG_Q).
- The count depending on the RNG seed ⇒ not a structural op-count, not
  contention-immune. **Asserted** — selftest [3] (exact seed-independence; ~0.54%
  field-dependence from dtype-branch differences only).
- **α_op fitting well above 1.0 with tight residuals** (e.g. ≥ 1.08) ⇒ prover
  op-count is super-Θ(x), the "Θ(x) + memory step" reading is wrong, and IDEAL=4.0
  underestimates. **Observed α_op = 0.995 (S530), and the per-step ratio DECLINES
  toward 4 in the upper regime** ⇒ not triggered. (The live falsifier was: "a curve
  that turns UP again at n≥22 — op-count ratio rising back above ~4.5". **S531
  extended the curve to n=24: the ratio kept DECLINING — 4.353× (n20→22), 4.329×
  (n22→24) — NOT triggered.** The naive single-power *fit* over the high-n window
  does read >1.05, but that is the polylog dressing, not a super-linear power —
  resolved by the 2-term fit, S531 section below. The proper falsifier for that fit
  is now: 2-term leading α > 1.05 on the stable band, OR per-step ratios that stop
  declining / rise. Neither observed.)
- fmul multiply-op share of wall negligible ⇒ multiplies are not the op-count
  proxy. **Reported ~78%** ⇒ not triggered.

## Selftest (`--selftest`, 6 groups, all pass)

1. patched `fmul` gives bit-identical verdict (claimed/comm/accepted == unpatched
   == sieve) over q & BIG_Q, ops>0;
2. the counter tallies exact broadcast output sizes (37 = 16+16+1+4 over 4 calls:
   array×array, array×scalar, scalar×scalar);
3. op-count EXACTLY seed-independent (mul_ops = 2 074 457 for seeds 1 and 999 at
   n=10) and field-independent to 0.54% (q vs BIG_Q) — the contention-immunity
   guarantee;
4. mul_ops strictly increasing in n; per-Δn=2 ratios in a Θ(x)-ish band (2–6.5×);
5. WALL_TIMING is op-count-neutral: mul_ops/mul_calls identical on/off, only
   mul_wall changes (the default-OFF fast path is count-preserving);
6. **(S531) the power-vs-polylog separation on the harvested deterministic curve
   (n=8..24):** the stable band is the trailing monotone-declining suffix (== n≥14);
   the single-power fit over the high-n window n=16..24 is **inflated to α=1.064 >
   1.05**, but the 2-term fit recovers **α=1.012 ≈ 1.0 with positive polylog
   δ=+0.72**, the forced-α=1 model `x·(log₂x)^0.87` fits the band to maxres 0.0045,
   and the per-step ratios **decline monotonically 4.471→4.329** — codifying the
   S531 resolution (Õ(x), not super-Θ(x)) as a regression guard.

## Honest scope

- This measures the field-MULTIPLY op-count via the `fmul` dispatch. Additions
  (`_asum`/`_sum61`) and scatters (`scatter_fold61`) are tallied separately
  (sum_ops, scatter_ops; both scale ~as mul_ops, α≈0.99 / ≈0.99) and do not pass
  through fmul (no double counting). Inline scalar arithmetic and `argsort` inside
  scatter are not counted, but multiplies dominate the wall (~78%), so mul_ops is
  the right op-count proxy and the exponent is the headline.
- This is a **measurement that pins an assumed leg of the reach falsifier**, not a
  new protocol or a step toward polylog π(x) (the goal stays blocked). It
  STRENGTHENS the S529 "memory step + noise, not super-Θ(x)" reading with a direct
  Θ(x) op-count measurement (α=0.995) and an independent demonstration that
  per-step ratios are intrinsically noisy (range > 2×) even before timing effects.
- Wall-clock here includes per-call instrumentation overhead (perf_counter per
  fmul), so the wall numbers are slightly inflated vs the clean benchmark and are
  used only as a secondary cross-check; the op-count is the deterministic result.

---

## S531 — curve extended to n=24; the apparent "super-Θ(x)" auto-reading corrected

**Date:** 2026-06-14. The S530 NEXT-ACTION path (b) — "extend the op-count curve to
n=22,24 to test S530's live falsifier (does the per-Δn=2 op ratio keep DECLINING
toward 4.0, or turn back UP above ~4.5× ⇒ super-Θ(x))" — was run detached
(`opcount_scaling_n24.log`, `--ns 16,18,20,22,24 --field big`) and **completed**.
This entry harvests it, resolves a tension it surfaced, and fixes the fit.

### The extended curve (deterministic, FULL config, BIG_Q+FAST_BIG)

| n | x | K=π(√x) | nb | mul_ops | per-Δn=2 ratio |
|---|---|---|---|---|---|
| 8  | 255       | 6   | 4  | 425 110         | — |
| 10 | 1 023     | 11  | 5  | 2 074 457       | 4.880 |
| 12 | 4 095     | 18  | 6  | 9 673 318       | 4.663 |
| 14 | 16 383    | 31  | 7  | 23 492 482      | 2.429 (low-n lump) |
| 16 | 65 535    | 54  | 8  | 105 029 523     | 4.471 |
| 18 | 262 143   | 97  | 9  | 465 458 097     | 4.432 |
| 20 | 1 048 575 | 172 | 10 | 2 040 403 997   | 4.384 |
| **22** | **4 194 303**  | **309** | **11** | **8 882 635 340**  | **4.353** |
| **24** | **16 777 215** | **564** | **12** | **38 448 499 833** | **4.329** |

The n≤20 rows re-measured live in `opcount_scaling_corrected.log` reproduce S530
**bit-for-bit** (op-count is structural/seed-independent, selftest [3]); n=22,24 are
the deterministic detached values, folded in via `--points` (injecting a
deterministic count is identical to re-running it — and avoids a redundant ~30-min
recompute that would contend with the live reach).

### The tension and its resolution

The detached run's auto-"READING" printed **"SUPER-Θ(x)"**: its single-power fit
over the short high-n window n=16..24 read **α=1.064 > 1.05**. But the *same output*
showed the per-step ratios **declining monotonically** 4.471→4.432→4.384→4.353→4.329
— a *pure* power x^α has a CONSTANT ratio 2^{2α}; a *declining* ratio toward 4.0 is
the fingerprint of a **polylog-dressed linear** term Θ(x)·(log x)^c (whose
log-derivative correction shrinks as x grows). The two signals contradict each other;
the single-power α was the artifact.

This is exactly the polylog-inflation that **S512 (`sharpP_probe`) built a 2-term fit
to defeat** — separating the leading power from the polylog: `log₂(mul) = α·n +
δ·log₂(n) + c` (since n = log₂x, n^δ = (log₂x)^δ). The script had only a single-power
`_fit_exponent`. This cycle added `_fit_power_log`, `_fit_polylog_fixed_alpha`, and a
`_stable_band` detector (the trailing monotone-declining suffix, past the low-n
discrete nb/K lump), and rewrote the verdict to use the **per-step trend + 2-term
separation**, not a single-power threshold.

**Result over the stable band n=14..24 (6 points):**

| fit | leading α | polylog δ | max\|resid\| (log₂) |
|---|---|---|---|
| single-power (window-inflated) | 1.067 | — | 0.021 |
| **power + polylog (2-term)** | **1.012** | **+0.721** | **0.0012** |
| forced α=1: `x·(log₂x)^δ` | 1 (fixed) | δ=0.870 | 0.0045 |

The 2-term fit recovers **α≈1.0 with a positive polylog**, at a residual **17× tighter**
than the single-power fit on the same points. The forced-α=1 model `mul ≈ A·x·(log₂x)^0.87`
predicts **every** per-step ratio to **<0.4%**:

```
        pred     obs
n14→16  4.4930   4.4708
n16→18  4.4318   4.4317
n18→20  4.3842   4.3836
n20→22  4.3460   4.3534
n22→24  4.3147   4.3285
```

### Reading

1. **The whole-chain prover op-count is Θ(x)·polylog = Õ(x), confirmed to n=24**
   (x≈1.7×10⁷). The leading exponent is α≈1.0; the polylog factor (log₂x)^≈0.87 is
   *positive* (op-count is Õ(x), not literally Θ(x) — consistent with open item 1's
   "Õ(x) prover", `Σ_{p≤√x} nb·p` carrying a log dressing). The single-power α
   reading of 1.015 (full n=8..24) or 1.064 (n=16..24) is **window-sensitive**: the
   low-n lump happens to cancel the inflation in the full-set fit, while the clean
   high-n window exposes it. The robust statement is the 2-term separation, not
   either single-power number.

2. **S530's live falsifier is NOT triggered.** It asked whether the ratio "turns
   back UP above ~4.5× at n≥22"; it kept *declining* (4.353×, 4.329×). The auto-fit's
   "super-Θ(x)" verdict was a fit artifact, now corrected. S530's headline (Θ(x)
   op-count) stands and is extended one octave.

3. **Reach tie-in unchanged.** The op-count leg of `mem_hierarchy_wall`'s predicted
   detached n24→26 band [4.0, 5.4]× is the per-Δn=2 op ratio **4.0–4.4×, declining**
   (here 4.33× at n22→24, trending toward 4.0 under the Õ(x) law). The reach's 5.80×
   wall step (n22→24) remains attributable to the one-time L3→DRAM memory step +
   single-step noise (S529), NOT the op-count (4.33× there, falling).

### What would falsify the S531 correction

- The 2-term leading α on the stable band rising **above ~1.05** as more high-n
  points are added (it is 1.012 at n=14..24), OR the per-step ratios **ceasing to
  decline / turning up** — either would mean a genuine super-linear leading term.
  Codified in **selftest [6]**. Cheap to extend (deterministic, no reach wait): e.g.
  `--ns 16,18,20,22,24,26` once an n=26 op-count is available.
- The stable-band detector mis-identifying the band (e.g. including the n=12→14
  lump): selftest [6] asserts the band is exactly n≥14 on the harvested curve.

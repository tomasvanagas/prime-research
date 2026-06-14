# sparse_precision_pi — P8 per-bit cost profile (precision-domain info barrier)

**Cycle:** S533 (2026-06-14). **Status:** P8 (OPEN_POSITIVE_TARGETS.md) **CLOSED-NEGATIVE**
with a measured mechanism. Adjacent-problem closure, **not** a goal advance (polylog π
stays blocked).

## The question (P8)

> Compute the first *k* bits of π(x) for varying *k*. Single-*k* is fixed cost;
> **batched on increasing *k* might amortise** — per-bit amortised cost might be
> **sub-linear** in total bits. **A-grade target: strictly sub-linear per-bit cost.**

This is the *precision-domain* sibling of the project's `SMOOTH + RANDOM` split: the
smooth part R(x) gives the top bits cheaply (polylog); the open question is whether the
remaining bits can be refined with a *diminishing* (sub-linear / concave) marginal cost.

## What was measured

`sparse_precision_pi.py` (one script, `--selftest`/`--profile`/`--anchor X`; 19 selftest
cases). It REUSES, verbatim, the numerically-stable Riemann explicit-formula evaluator
from `../explicit_formula_witness/explicit_formula_witness.py` (S518) — the CLOSED-row-31
`li(exp())` branch-cut bug is fixed and selftested there, imported so it cannot recur.

For a fixed anchor X, over a window of x near X (RMS-smoothed to average the GUE phase
noise), the script computes:
- the error envelope `E(N) = RMS_x | π(x) − approx(x, N zeros) |`,
- `B = bit-length(π(X))`, the free bits `k0 = B − log2(E(0))` (R(x) alone, no zeros),
- the per-bit cost map `N(k)` = zeros to **robustly resolve the first k bits** of π(X)
  (`E < 2^{B−k}` and stays below — `last_crossing`),
- the cost exponent `s = d log2 N(k) / dk` and the envelope exponent `eta` in `E ∝ T^{−eta}`.

## Result (window 24, BIG cached-zero budget ≤ 8000 zeros, x = 10⁴…10⁶)

| X | π(X) | B | free bits k0 | k0/B | per-bit slope s | R² | hard bits |
|---|---|---|---|---|---|---|---|
| 10⁴ | 1229 | 11 | 9.70 | 0.88 | — (1.3 hard bits) | — | 1.30 |
| 3·10⁴ | 3245 | 12 | 9.94 | 0.83 | 3.72 (3 pts, 1 GUE spike) | 0.86 | 2.06 |
| 10⁵ | 9592 | 14 | 11.74 | 0.84 | **1.53** | 0.95 | 2.26 |
| 3·10⁵ | 25997 | 15 | 10.28 | 0.69 | **2.51** | 0.95 | 4.72 |
| 10⁶ | 78498 | 17 | 12.14 | 0.71 | **2.09** | 0.95 | 4.86 |

- **Per-bit cost slope `s = 1.5–2.5`** at every well-sampled anchor (≥3 clean bit-crossings,
  R²≈0.95). `s ≥ 1` ⇒ each additional hard bit costs **≥ 2× (measured ~3–6×)** the zeros of
  the previous one. This is a **GEOMETRIC, super-linear, convex** marginal cost — the exact
  OPPOSITE of P8's sub-linear hope. (The 3·10⁴ s=3.72 has only 3 points dominated by one 81×
  GUE spike, R²=0.86 — an over-read; the trustworthy band is the [1.5, 2.5] from the
  ≥4-point R²≈0.95 anchors.)
- **Mechanism — envelope `E ∝ T^{−eta}`, eta ≈ 0.22–0.59 < 1** (the L² zeta-zero tail:
  omitted zeros γ>T each contribute ~√x/(γ log x) with random phase, so the RMS error is an
  L² sum ~ (√x/log x)·√(log T / T) ⇒ eta ≈ ½ lifted to ~0.55–0.65 by the √log T; finite-x /
  short-band fits read it lower and noisily). `eta < 1` is *why* the per-bit slope `1/eta > 1`
  — the error decays **slower than 1/T**, so a bit costs **more** than a doubling of zeros.
- **Free fraction `k0/B = 0.88 → 0.71`, DECLINING with x.** At finite measurable x, R(x) is
  far better than its √x worst-case bound (`E0 = |π−R| ~ √x/log x`, e.g. only ≈29 at x=10⁶),
  so R nails ~70–90% of the bits for free — *more* than the asymptotic ~50%. The decline
  0.88→0.71 over x=10⁴→10⁶ is the √x random zone slowly opening; `k0/B → ~½` only as x→∞ (the
  CLAUDE.md "~50% of digits", here measured in the precision domain).
- **No amortisation.** The explicit formula is ALREADY incrementally batched — refining k→k+1
  bits reuses the N(k) zeros and adds dN — yet the cost is geometric, so the total to reach
  the full B bits = N(B) is dominated by the last hard bit. Amortised cost per hard bit
  `N(B)/(B−k0)` = 251, 1292, 941 zeros/bit at x = 10⁵, 3·10⁵, 10⁶ — large and rising overall
  (251→941 across the clean 10⁵→10⁶ endpoints, ×3.7; the 3·10⁵ 1292 is GUE-inflated by one
  unlucky expensive last bit) — it does NOT fall toward a bounded sub-linear amortised cost.
  Batching across precision saves only a **constant factor** (the geometric-series sum ≈ 2×
  the largest term).

**VERDICT: P8's "strictly sub-linear per-bit cost" is REFUTED.** The per-bit cost is FLAT
(free, polylog) for the top ~70–90% of bits (R(x)) and GEOMETRIC (×3–6/bit, super-linear)
for the hard bits below; there is no sub-linear regime, and precision-batching does not
amortise.

## Why this is forced (the S511/S518 tie-in) — and the method-independence

This RE-DERIVES the S511/S518 √x information barrier in the precision domain (a measurement
confirming a known wall from a new angle — **not** novelty, per the contract):

- **S518** (x-domain): settling π(x) EXACTLY needs `T ~ √x·polylog` zeros, info dense across
  them. The per-bit view explains *why*: the last hard bit alone needs `~√x` zeros (geometric
  growth from the free zone), so the full count is √x.
- **The slope ≥ 1 (no sub-linear regime) is kernel-independent**, forced by the floor: a
  sub-linear per-bit slope (s<1) would mean the deep bits are cheap ⇒ you could reach exact
  (all of the √x-worth of hard bits) in **sub-√x zeros**, contradicting S518's measured √x
  settling and S511's joint information floor (Θ(√x) integer-independent hard bits). So `s ≥ 1`
  is **information-forced**, method-independent.
- The **measured** s = 1.5–2.5 is the *sharp-truncation* Riemann-R value (eta≈½). A smoothed /
  optimal kernel (Galway) improves the tail constant and could push s toward the floor's `1`
  (doubling), **but not below 1** — the floor caps it. Either way the per-bit cost is ≥ a
  doubling: never sub-linear.

## Honest scope / caveats

- **Adjacent NEGATIVE closure, not a goal advance.** Confirms a known wall (the info barrier)
  from the precision axis; polylog π(x) stays blocked.
- **Reach of the cached zeros.** With 8000 zeros (height ≤ 8148) the evaluator settles π(x)
  exactly only up to x ~ 10⁶, where there are only ~5 hard bits. Measuring *more* hard bits
  requires larger x, hence more zeros — that requirement *is* the √x wall, so the narrow
  per-x bit range is itself a manifestation of the barrier, not a measurement artifact.
- The envelope-eta fit is GUE-noisy (R² 0.33–0.82); it is a *corroborating mechanism*. The
  well-determined quantity is the per-bit slope `s` (R²≈0.95), and it is the headline.
- The synthetic control (selftest [9]/[10]) proves the pipeline **would detect a sub-linear
  regime if one existed**: a synthetic `E ∝ T^{−2}` envelope (sub-linear, slope 0.5) is
  recovered as slope 0.48; the real-data s = 1.5–2.5 is then a genuine negative, not a blind
  spot.

## What would falsify this

1. A per-bit slope `s < 1` at a well-sampled anchor (≥4 clean bit-crossings, R²>0.9) — a
   genuine sub-linear / diminishing-marginal regime (a partial P8 win). **Observed: s ≥ 1.53,
   never < 1.** The synthetic eta=2 control (s=0.48) confirms the method detects s<1.
2. The free fraction `k0/B` rising toward 1 with x (R giving ALL bits, no hard zone).
   **Observed: it DECLINES 0.88→0.71 — the hard zone opens.**
3. A smoothed-kernel explicit-formula variant with `eta > 1` AND a sub-√x exact-settling count
   at some x — which would contradict S518/S511 directly. **Not observed; ruled out by the
   floor for any sieve-reconstructing/analytic witness (S511/S518).**

## Reproduce

```
python3 sparse_precision_pi.py --selftest          # 19 cases incl. synthetic controls
python3 sparse_precision_pi.py --profile --window 24   # the table above (~6 min)
python3 sparse_precision_pi.py --anchor 1000000    # single-anchor bit-by-bit N(k)
```

Run log: `profile_run.log`. Depends on `../explicit_formula_witness/` (imported) and
`../../../data/zeta_zeros_8000.txt`.

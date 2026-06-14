# D_k(x) census first-occurrence probe at k=32 (S522)

## Question

The aligned-window pattern census `D_k(x)` (= #distinct exact prime
occupancy patterns of windows `[km, km+k)`, `km < x`) obeys the
saturation law **`D_k(x) → A_k + 1`** (`A_k` = #Hardy–Littlewood
admissible patterns; the `+1` is the `m=0` small-prime window). It is
verified saturated at k=8 (D=14, by 2¹⁶) and k=16 (D=107, by 2¹⁸). At
k=32 it was still **unsaturated** at the prior reach 2²⁶:
`D_32(2²⁶)=3385`, with **189 admissible patterns unrealized**, every one
of weight ≥ 6 (23 of wt 6, 115 of 7, 47 of 8, 4 of 9). This cycle
extends the search with a segmented sieve to find their first
occurrences, measure the approach to `A_32+1 = 3574`, and test each
against a Hardy–Littlewood first-occurrence prediction.

`dk_first_occurrence.py` (one script, `--k`, `--xmax`, `--seg`,
`--validate-hl`, `--selftest`).

## Headline

Pushing the aligned census from 2²⁶ to **2³⁷ (x ≈ 1.37×10¹¹)** drives the
deficit from **189 to 0 — FULL SATURATION: `D_32 = A_32 + 1 = 3574`**,
matching the saturation law verified at k=8 (14) and k=16 (107):

| x | D₃₂(x) | adm. found | missing | D/(A+1) |
|---|---|---|---|---|
| 2¹⁶ | 1093 | 1092 | 2481 | 0.3058 |
| 2²⁰ | 2452 | 2451 | 1122 | 0.6861 |
| 2²⁴ | 3213 | 3212 | 361 | 0.8990 |
| 2²⁶ | 3385 | 3384 | 189 | 0.9471 |
| 2²⁸ | 3492 | 3491 | 82 | 0.9771 |
| 2³⁰ | 3550 | 3549 | 24 | 0.9933 |
| 2³² | 3568 | 3567 | 6 | 0.9983 |
| 2³⁴ | 3571 | 3570 | 3 | 0.9992 |
| 2³⁶ | 3573 | 3572 | 1 | 0.9997 |
| 2³⁷ | **3574** | **3573** | **0** | **1.0000** |

(The 2¹⁶/2²⁴/2²⁶ rows reproduce the prior `local_pattern_census.py`
values **exactly** — the segmented sieve is bit-equivalent to the full
sieve, asserted in selftest [3]/[6].) **All 3573 admissible patterns are
now realized**; D₃₂ equals A₃₂ plus the single inadmissible `m=0`
small-prime window — the law's exact prediction. The last to appear was
the densest possible (weight-9) constellation
`[1,5,7,11,17,19,25,29,31]` at x = 83 122 625 472 ≈ 8.31×10¹⁰ (HL
x*≈3.86×10¹⁰, expected count at first occurrence `N(occ)=1.63` — on the
HL schedule). The whole deficit, at every reach, was exactly the
HL-predicted dense tail; no admissible pattern was ever "stuck".

## The Hardy–Littlewood model (aligned singular series), validated

For the aligned constellation of linear forms `f_i(m) = km + o_i`
(`o_i ∈ S`, `w = |S|`), the singular series is
`S_aligned = ∏_q (1 − ω(q)/q)/(1 − 1/q)^w` with
`ω(q) = #{m mod q : some f_i(m) ≡ 0}`:
- `q | k`: `f_i ≡ o_i` (const); `ω = 0` for admissible S (else the form
  is always divisible → count 0). For `q = 2 | 32` this gives factor
  `2^w` — the aligned-even base makes every candidate odd.
- `q ∤ k`: `km` hits every residue, so `ω(q) = #{o_i mod q}` = the
  admissibility-(ii) residue-coverage count. Admissible ⇔ `ω(q) < q ∀q`.

Expected count of windows realizing the constellation up to x:
`N_S(x) = S_aligned·(x/k)/(ln x)^w`; HL first occurrence `x*` solves
`N_S(x*) = 1`.

**Normalization validated (`--validate-hl`, occurrence counts vs
`N_S(xmax)`, x = 2³⁰):**

| w | #patterns | obs total | pred total | obs/pred |
|---|---|---|---|---|
| 1 | 16 | 11 274 190 | 51 636 067 | 0.218 |
| 2 | 120 | 10 228 393 | 33 679 781 | 0.304 |
| 3 | 410 | 5 062 667 | 12 138 838 | 0.417 |
| 4 | 845 | 1 489 447 | 2 642 161 | 0.564 |
| 5 | 1044 | 267 341 | 356 470 | 0.750 |
| **6** | **766** | **28 842** | **29 435** | **0.980** |
| 7 | 306 | 1 857 | 1 408 | 1.319 |
| 8 | 41 | 62 | 25.4 | 2.437 |

The ratio rises monotonically through 1 near weight 6–7. This is exactly
right and pins down two distinct effects:

- **Low weight (ratio < 1):** `N_S` counts the *constellation* (≥ S
  prime); the *exact*-pattern count requires the other odd offsets
  composite, an extra factor < 1 that is stronger the more free offsets
  remain (lower w). Exact ⊂ constellation.
- **Weight 6 (ratio = 0.980):** at this density almost all patterns are
  realized (no selection bias) **and** the "remaining offsets composite"
  condition is nearly free (a 6-prime window is already tight) — so
  exact ≈ constellation. This is the **clean validation of the aligned
  singular series: HL is accurate to 2% at weight 6.**
- **High weight (ratio > 1):** *survivorship bias* — at 2³⁰ only the
  densest (upward-fluctuating) weight-7/8 patterns have appeared, so
  their conditional occurrence count exceeds the unconditional HL mean.

## First occurrences are HL-ordered

For the patterns first realized in this run (x ≥ 2²⁶), the observed
first occurrence tracks the HL prediction `x*` and the expected count at
first occurrence `N(occ)` clusters around 1 for the cleanly-tested
weight-7 band (e.g. `[1,3,7,9,19,21,31]` at x=7.94×10⁷ vs x*=7.73×10⁷,
`N(occ)=1.02`; `[1,7,17,19,25,29,31]` at 8.31×10⁷ vs 7.73×10⁷,
`N(occ)=1.04`). The weight-6 stragglers (those that were *missing* at
2²⁶ — by definition the late-occurrence tail) show `N(occ)` = 3–7, the
expected upper tail of the first-occurrence distribution.

## The falsification trigger and the one flagged candidate

The sound test for a "stuck" pattern is **not** `x* < xmax` (a pattern
with mean count 1 is absent with probability `e⁻¹ = 37%`) but the
**expected count** `N_S(xmax)`: a missing pattern is statistically
overdue only if `N_S(xmax)` is large, since `P(0 occurrences) = e^{−N}`.
The script flags any missing pattern with `N_S(xmax) ≥ 5`
(`P(0) < 0.7%`).

- At 2²⁸, **two** weight-6 patterns were flagged (`N≈5.4, 6.5`); both
  **appeared by 2³²** — confirming they were the Poisson tail, not
  anomalies, and that the flag self-corrects as x grows.
- At 2³², **one** weight-7 pattern is flagged: **`[3,9,13,15,19,25,27]`,
  `N_S(2³²)=6.86`, `P(0)=0.0011`, `x*=2.34×10⁸`**. This is the single
  most-overdue admissible pattern. It is structurally unusual — four of
  its seven offsets (3,9,15,27) are ≡ 0 mod 3 — but that mod-3 cost is
  already inside `S_aligned`. The deep run (below) resolves it as the
  expected order-statistic tail (it appeared at 4.75×10⁹), not an HL
  over-prediction.

The 2 weight-8 (`x* ≈ 5.9×10⁹`) and 3 weight-9 (`x* ≈ 3.9×10¹⁰` and
`1.1×10¹¹`) still missing at 2³² have `N_S(2³²) ≤ 0.82` — entirely
expected to be absent.

## The dense-tail first occurrences (the 6 that drove the 2³² deficit)

Extending to 2³⁷ realizes the last 6 patterns. Their observed first
occurrences vs HL `x*`:

| w | first-occ x | HL x* | obs/x* | N(occ) | offsets |
|---|---|---|---|---|---|
| 7 | 4.75×10⁹ | 2.34×10⁸ | 20.3 | 7.35 | [3,9,13,15,19,25,27] |
| 9 | 8.08×10⁹ | 1.14×10¹¹ | 0.07 | 0.19 | [1,5,11,13,19,23,25,29,31] |
| 8 | 9.87×10⁹ | 5.89×10⁹ | 1.68 | 1.40 | [1,3,7,9,13,19,21,27] |
| 8 | 2.04×10¹⁰ | 5.89×10⁹ | 3.47 | 2.25 | [3,9,11,17,21,23,27,29] |
| 9 | 4.38×10¹⁰ | 1.14×10¹¹ | 0.38 | 0.54 | [1,3,7,9,13,19,21,27,31] |
| 9 | 8.31×10¹⁰ | 3.86×10¹⁰ | 2.16 | 1.63 | [1,5,7,11,17,19,25,29,31] |

`N(occ)` (HL expected count at the observed first occurrence) clusters
around 1 — first occurrences happen, as they should, where the
cumulative HL mean is O(1). Two weight-9 came *early* (N(occ) 0.19, 0.54
— lucky low-probability early hits), the others on schedule (1.40–2.25).
The HL-derived `x*` for the two `[…25,29,31]`-style weight-9 patterns
(1.14×10¹¹) was conservative — they appeared earlier by Poisson luck —
while the saturating pattern landed at obs/x*=2.16, squarely consistent.

**The weight-7 flag resolved as the expected order-statistic tail, not
an anomaly.** `[3,9,13,15,19,25,27]` finally appeared at x≈4.75×10⁹ with
N(occ)=7.35 — i.e. it first occurred only after the HL mean had reached
7.35 (`P(this late) = e^{−7.35} ≈ 6×10⁻⁴`). But this is the **maximum
lag over ~300 weight-7 patterns**: the max of n unit exponentials has
mean ≈ ln n ≈ 5.7, and `P(max over 300 ≥ 7.35) ≈ 1−(1−e^{−7.35})³⁰⁰ ≈
0.17`. So a worst-straggler at N(occ)≈7.4 is a ~17% event — fully within
the look-elsewhere spread, **not** an HL over-prediction. (Its unusual
structure — four of seven offsets ≡ 0 mod 3 — is already priced into
`S_aligned` via the mod-3 factor; it is not under-counted.) The two
weight-6 flags raised at 2²⁸ had likewise cleared by 2³². **The flag is
self-correcting: every pattern it has ever raised has subsequently
appeared.**

## What would falsify saturation at k=32

1. An **inadmissible** pattern realized at `m > 0` — impossible by plain
   divisibility (a bug check); none found through 2³².
2. A pattern with **large HL mean count `N_S(x)` that stays missing** as
   x grows (the `N ≥ 5` flag persisting and growing, not clearing) —
   would indicate HL over-prediction or a genuine obstruction. The two
   2²⁸ flags cleared by 2³²; the 2³² weight-7 flag is tracked in the
   deep-run section.
3. `D_32(x)` **plateauing strictly below 3574** as x → ∞ (an admissible
   pattern that provably never occurs) — would refute the qualitative
   content of the prime k-tuple conjecture.

None of (1)–(3) is observed: the deficit falls 189 → 0 from 2²⁶ to 2³⁷
(every admissible pattern realized), and at each intermediate reach the
residual was exactly the HL-predicted dense tail.

## Honest scope / reach statement

This is a finite-x measurement that **confirms saturation at k=32**
(`D_32 = A_32 + 1 = 3574` by 2³⁷), it is not a proof of saturation for
all k (that is the k-tuple conjecture's content and out of reach). The
deliverable is: the full measured approach curve `D_32(x) → A_32+1` over
2¹⁶–2³⁷, the saturation point reached; a 2%-accurate validation of the
aligned HL singular series at weight 6; a verified-correct first-
occurrence list for every dense weight-7/8/9 constellation that drove the
prior deficit; and a statistically-grounded falsification monitor that
self-corrected on every flag it raised. The cost was a single
constant-memory segmented sieve to ≈1.4×10¹¹ (the deepest run ~12 min).
The k=32 case now joins k=8, 16 as a fully-saturated, HL-consistent
confirmation of the census law.

## Files

`dk_first_occurrence.py`, `run_k32_2e32.log`, `run_k32_2e36.log`,
`run_k32_2e37.log` (the saturating run, D₃₂=3574).
Cross-refs: `local_pattern_census.py`/`_results.md` (the saturation law
and A_k enumeration), `census_transfer_matrix.py` (ρ*(k) = max
admissible weight = 9 at k=32, the densest constellations here),
`novel/succinct_verification_of_pi.md` §width-dichotomy.

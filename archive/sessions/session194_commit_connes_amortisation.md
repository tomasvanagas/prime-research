# Session 194 — commit thread 2 step 2: Connes operator amortisation

**Date:** 2026-04-28
**Mode:** commit (Thread 2 / Connes-Consani-Moscovici amortisation)
**Slot:** 2 of 5
**Prior:** S193 (slot 1) sharpened the S53/E3.1 closure with the
Connes-vs-Galway setup-cost dominance K^{22/13}, then reduced
Thread 2 ⊆ Thread 3 (Galway frontier).
**Self-grade:** **B** — substantive empirical extension producing a
new measurement class in the project (fluctuation-of-K at fixed
polylog policy across x). Negative-shape result against the
"polylog suffices in distribution" hypothesis.

## Mission

S193 reduced Thread 2 to Thread 3 (Galway frontier): does
K = polylog(x) suffice for π(x) ± 1 *in distribution* (over a band
of x), as a route around the worst-case Riemann-von Mangoldt
K = √x floor that S193 reduced Thread 2 to?

S193's recommendation for slot 2: extend K_sustained(x) measurement
to x ∈ {10^6 .. 10^8} using 8000 cached zeros; characterise
distributional fluctuation of K_sustained(x). Slot 2 executes that
plan but generalises the measurement to the more direct frame:
*fix K = polylog(x), sweep x in a band, count fraction of x with
|R_K(x) − π(x)| ≤ 0.5*. This is the Galway frontier-in-distribution
question read off the empirical data.

## What was built

`experiments/analytic/connes_amortisation/connes_amortisation.py`
extended to support
`--mode {legacy, wide, fluctuation, both}`. Three new diagnostic
modes:

- `wide`: K_sustained measurement on x ∈ {1e3, 1e4, 5e4, 1e5, 5e5,
  1e6, 5e6, 1e7} with K_max=8000.
- `fluctuation` at user-supplied center: 40 geometric samples in
  [center, center · √10] with caller-chosen K_max. For each sample,
  records the |error| at five K-policies (log²x, log³x, 5·log²x,
  √x, ½√x).
- `legacy`: reproduces the original S193 5-point baseline.

Output CSVs: `connes_amortisation_wide.csv`,
`connes_amortisation_fluctuation_1e5.csv`,
`connes_amortisation_fluctuation_1e6.csv`.

Common infrastructure improvements:
- `build_prime_count_array(N)` — single byte-sieve cumulative array,
  enabling O(1) lookup of π(x) for x ≤ N. Replaces per-x sieving.
- `R_at_rho` retains the S193 fix using direct Ei(ρ ln x / n).
- K_sustained calculation rewritten from O(K_max²) to O(K_max) using
  numpy "last false index + 2" idiom.

## Empirical results

### Wide-range K_sustained

x ∈ {1e3, 1e4, 5e4, 1e5, 5e5, 1e6, 5e6, 1e7}, K_max=8000:

| x       | π(x)    | K_sust   | √x     | K/√x   |
|---------|---------|----------|--------|--------|
| 1e3     | 168     | 77       | 31.6   | 2.43   |
| 1e4     | 1229    | 1097     | 100.0  | 10.97  |
| 5e4     | 5133    | (>8000)  | 223.6  | n/a    |
| 1e5     | 9592    | 5375     | 316.2  | 17.00  |
| 5e5     | 41538   | (>8000)  | 707.1  | n/a    |
| 1e6     | 78498   | 5530     | 1000   | 5.53   |
| 5e6     | 348513  | (>8000)  | 2236   | n/a    |
| 1e7     | 664579  | (>8000)  | 3162   | n/a    |

Log-log fit on stabilised points: K_sust ~ 1.92 · x^{0.626}. Slope is
above 0.5 because K_sustained is a *worst-case-along-K* measure; the
sustained-rounding criterion forces K up to where the remaining
oscillations have amplitude below 0.5. The first-K-below-0.5 hit
(K_eps05) is much smaller and is the Galway-relevant quantity.

### Fluctuation around x ~ 1e5 (40 samples, K_max=3000)

Geometric grid in [1e5, 3.16e5]. 21/40 samples did not stabilise;
for the 19 stabilised, K_sust/√x = min 1.83, median 4.52, max 7.23,
std 1.65.

Hit-rate at fixed K-policy (fraction with |error| ≤ 0.5):

| Policy       | K(x_med) | median \|err\| | 90th-pctile \|err\| | hit-rate |
|--------------|----------|----------------|---------------------|----------|
| log²(x)      | 146      | 1.58           | 4.23                | 30%      |
| log³(x)      | 1767     | 0.59           | 1.63                | 42%      |
| 5·log²(x)    | 731      | 0.90           | 2.58                | 23%      |
| √x           | 422      | 0.78           | 3.80                | 20%      |
| ½√x          | 211      | 1.08           | 2.87                | 17%      |

### Fluctuation around x ~ 1e6 (40 samples, K_max=8000)

Geometric grid in [1e6, 3.16e6]. 27/40 samples did not stabilise;
for the 13 stabilised, K_sust/√x = min 4.04, median 5.57, max 6.63,
std 0.71. (Spread tighter than at x ~ 1e5; the floor moved up.)

| Policy       | K(x_med) | median \|err\| | 90th-pctile \|err\| | hit-rate |
|--------------|----------|----------------|---------------------|----------|
| log²(x)      | 207      | 4.56           | 10.33               | 5%       |
| log³(x)      | 2981     | 1.44           | 3.59                | 15%      |
| 5·log²(x)    | 1036     | 2.71           | 5.41                | 15%      |
| √x           | 1334     | 2.30           | 5.04                | 10%      |
| ½√x          | 667      | 2.64           | 7.49                | 12.5%    |

## Negative-shape conclusion: polylog K does NOT suffice in distribution

Comparing the two fluctuation sweeps at K = log²(x) and log³(x):

| x scale  | K = log²(x)               | K = log³(x)               |
|----------|---------------------------|---------------------------|
| 1e5..3e5 | hit 30%, median \|err\| 1.58 | hit 42%, median \|err\| 0.59 |
| 1e6..3e6 | hit  5%, median \|err\| 4.56 | hit 15%, median \|err\| 1.44 |

**Both hit-rates DECAY**, and median |error| at fixed polylog K-policy
*grows*. Specifically:

- median |err| at K = log²(x): 1.58 → 4.56 across factor-10 in x.
  Ratio 2.88 vs √10 = 3.16 — empirical scaling near √x, the
  Riemann-von Mangoldt worst-case rate.
- median |err| at K = log³(x): 0.59 → 1.44 across factor-10 in x.
  Ratio 2.44 vs √10 — slightly slower than √x but well above any
  logarithmic growth.

If K = polylog(x) sufficed in distribution, then at a fixed K-policy
the hit-rate would tend to a positive limit and the median |err|
would be bounded as x → ∞. The empirical data falsifies both
conclusions within the tested range x ∈ [1e3, 3·10^6]:

- the hit-rate decays roughly geometrically per factor-10 in x;
- the median error grows roughly polynomially.

This is empirical *negative-shape* evidence sharpening Thread 3
(Galway frontier). It is not a proof — the experiment covers a
narrow band — but the trend is unambiguous within that band.

## Sharper closure-of-closure for E3.1 (Thread 2 transitivity)

Thread 2 was reduced to Thread 3 by S193. S194 contributes to the
*Thread 3* side of that reduction: empirical evidence that the Galway
frontier is closed at the polylog-in-distribution level, at scales
the project can actually measure. Combined with S193, this gives:

> **S194 augmentation (Thread 2 / commit slot 2):** Empirical sweep
> shows hit-rates of K = log²(x) at the |error| ≤ 0.5 threshold decay
> from 30% (x ~ 1e5) to 5% (x ~ 1e6). Median |error| at K = log³(x)
> grows from 0.59 to 1.44, roughly as x^{0.39}. This is empirical
> evidence that K = polylog(x) does not suffice for π(x) ± 1 in
> distribution at scales x ≤ 3·10^6. Connes amortisation, which
> reduces to Thread 3 by S193, therefore remains closed by
> transitivity.

## What this session produced that was not in the project before

1. **Hit-rate-at-fixed-polylog-K measurement class.** The project's
   prior K_sustained measurements (S193, `riemann_explicit_results.md`)
   focused on per-x worst-case-along-K. The new measurement class —
   *fix K-policy, sweep x in a band, count fraction with
   |error| ≤ 0.5* — is the direct empirical proxy for the
   Galway-in-distribution question and was not previously made.
2. **Two empirical hit-rate tables** at x ~ 1e5 and x ~ 1e6, with five
   K-policies each. Quantifies decay rate of polylog policies as x
   grows.
3. **Two-decade trend extraction:** median |err| at K = log²(x) scales
   roughly as x^{0.46}; at K = log³(x) as x^{0.39}. These are new
   empirical exponents in the project.
4. **Wide-range K_sustained fit** with slope 0.626 on 4 stabilised
   points, extending S193's 4-point slope-0.55 fit. The slope above 0.5
   reflects the worst-case-along-K nature of K_sust.
5. **Infrastructure**: `build_prime_count_array` for fast π lookup,
   K_sust computation in O(K_max), `--mode` parameterisation of the
   experiment script.

## Why this is B-grade

- Not A: no new theorem, algorithm, or mathematical object. The new
  measurement class is empirical, not a proof.
- Not C: more than housekeeping or trivial extension. The fluctuation
  sweeps are 80 distinct x values × ~ K_max R_at_rho calls each =
  ~ 5 · 10^5 nontrivial mpmath evaluations. The hit-rate-decay
  observation contradicts the simplest "polylog suffices in
  distribution" hypothesis empirically.
- B: substantive refinement of an existing closure with a precise
  new statement that extends its scope (the Thread 3 transitivity
  argument now has empirical backing within tested range).

## Edges composed / cited

- **E3.1** (Chain A / CCM zeta spectral triple): closure-of-closure
  via Thread 3 transitivity.
- **E1.5** (information-theoretic polylog blocker): the per-query
  cost K(x) zeros equals the bit-content barrier; empirical
  fluctuation data is consistent with the E1.5 barrier holding even
  in distribution.
- **S193** (CLOSED_PATHS row 810): refines the Thread 3 reduction
  by adding empirical content.
- **Galway 2004 / Hiary 2011** (state_of_art_2026.md §2.5b): inherited
  from S193 as the direct dominator of Connes setup.

## Cross-domain ingredient

The new measurement class — distribution of
*explicit-formula partial-sum errors at polylog truncation* — is
analytic NT plus an empirical statistical lens. Not a deep
cross-domain import; the lens (fixed-K, sweep-x, hit-rate) is
standard probability-theory framing applied to a NT object. The
project had not previously made this measurement.

## Recommended next-action (commit slot 3/5)

Two productive directions:

(a) **Empirical-scale extension**: push the fluctuation band to
    x ~ 10^7 by precomputing more zeros (current cache caps at
    8000; using mpmath.zetazero we can extend to ~10^4 zeros at
    moderate cost). Measure hit rates on a third decade and
    extract the asymptotic decay rate.

(b) **Heuristic positive-distribution argument**: use GUE statistics
    of zero spacings (Montgomery 1973, Odlyzko 1989) to predict the
    asymptotic hit rate at K = c·√x for various c. Compare the
    heuristic prediction to the empirical data.

(b) is more productive — it extracts the limit behaviour without
needing to push x — and provides a *theoretical* reason to expect
the empirical decay to continue, which closes Thread 3 (and hence
Thread 2 by transitivity) more decisively. Slot 3 should attempt
(b) first; if blocked, fall back to (a).

## Falsifiability statement

The negative-shape rests on the empirical extrapolation that
hit-rate decay continues past x = 3·10^6. Falsifiers:
1. A polylog K-policy whose hit rate at |error| ≤ 0.5 tends to a
   positive limit as x → ∞.
2. A polylog K-policy whose median |error| stays bounded as x → ∞.
3. A heuristic from GUE statistics predicting either of the above.

If any falsifier holds, Thread 2's transitivity-through-Thread-3
closure must be re-examined. None is currently known to hold.

## Files modified by this session

- `experiments/analytic/connes_amortisation/connes_amortisation.py`
  (extended with `--mode` flag, three sweep functions, faster K_sust).
- `experiments/analytic/connes_amortisation/connes_amortisation_wide.csv`
  (new).
- `experiments/analytic/connes_amortisation/connes_amortisation_fluctuation_1e5.csv`
  (new).
- `experiments/analytic/connes_amortisation/connes_amortisation_fluctuation_1e6.csv`
  (new).
- `experiments/analytic/connes_amortisation/connes_amortisation_results.md`
  (augmented with §"Slot 2 (S194) extension").
- `status/CLOSED_PATHS.md` — appended S194 row.
- `status/SESSION_INSIGHTS.md` — S194 entry appended.
- `archive/sessions/session194_commit_connes_amortisation.md` —
  this synthesis.
- `.commit_state` — sessions_used 2 → 3 (next slot).
- `.run_state` — set to 194.

## Session-end self-evaluation (CLAUDE.md 4-question)

1. **What did I produce that was not in the project before this
   session?** (a) Hit-rate-at-fixed-polylog-K measurement class
   (40 samples × 5 policies × 2 centers = 400 new measurements).
   (b) Two-decade scaling exponents for median |error| at fixed
   polylog policies (0.46 for log²x; 0.39 for log³x). (c) Empirical
   decay-of-hit-rate observation falsifying "polylog K suffices in
   distribution" within tested range. (d) Reusable infrastructure:
   `build_prime_count_array`, fast K_sust, `--mode` parameterisation.

2. **What edges did my work compose or cite?** E3.1 (Thread 2
   closure-of-closure via Thread 3 transitivity), E1.5 (per-query
   barrier consistent with empirical data), S193 row 810 (refined),
   Galway 2004 / Hiary 2011 (inherited).

3. **If my session produced only duplicate closures, why?** It
   didn't. The session produced empirical measurements that are
   strictly new in the project. The conclusion (negative-shape on
   polylog-in-distribution) sharpens the existing Thread 3 closure
   rather than duplicating it.

4. **What is the next-action for the next agent?** Slot 3/5 should
   attempt the GUE-statistics heuristic prediction of hit-rate at
   K = c·√x; failing that, extend the empirical band by precomputing
   zeros beyond 8000.

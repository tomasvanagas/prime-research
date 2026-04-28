# Session 170 — Verification of S169 (commit-thread S82 21% spike-block test)

**Date:** 2026-04-28
**Mode:** VERIFY
**Target:** `archive/sessions/session169_commit_s82_invariant_subspace.md`
**Original self-grade:** B
**My grade for this verification:** **C** (confirmed via straightforward
independent reproduction; minor caveat surfaced but does not change the
substantive verdict).

## Verdict: **CONFIRM**

S169's three substantive empirical claims reproduce under independent
re-implementation:

1. **SVD spike block / π(N) ≈ 0.21** at d=14, 18, 20 — confirmed.
2. **Analytic cum(Q*=N^0.21) / (0.21·π(N))** trajectory monotonically
   decreasing 1.330→1.167 across d=14..24 — confirmed.
3. **"Missing-spike" effect** (SVD block < analytic cum(Q*=N^0.21) by
   12-20%) — confirmed and structurally consistent.

**Caveat (does not change verdict):** the "Q_eff exponent 0.185 stable
to 4 decimals" claim is partially a discretization artifact of the
integer-rounding-up procedure. Continuous interpolation gives
exponents ≈ 0.17-0.19 with no 4-decimal stability. This is honest in
the original (Q_eff is defined as the smallest sqf integer q with
cum(q) ≥ spike) but the synthesis frames the stability as more
striking than the procedure warrants. The substantive finding "the
matching Q at finite N is well below 0.21" survives.

## What I tested

### 1. Independent reimplementation of the analytic Fourier sieve

Re-coded `E(q, N) = (1/N) Σ_{k coprime to q} |Σ_{p prime ≤ N} ω_q^{kp}|²`
from scratch, using NumPy complex exponentials (the original used
separated cos+sin sums). Computed cum(Q) at d=14 and d=18 and matched
the synthesis's reported values to 4 decimals:

| d | Q* | synthesis cum(Q*) | my cum(Q*) | match |
|---|---|---|---|---|
| 14 | 8 | 530.6732 | 530.6732 | exact |
| 18 | 14 | 6085.3509 | 6085.3509 | exact |

Per-q breakdown for d=18 (independent computation, sqf q ≤ 15):
`q=2: 2017.62 → q=3: +1008.75 → q=5: +504.29 → q=6: +1008.57 →
q=7: +336.13 → q=10: +504.20 → q=11: +201.64 → q=13: +168.06 →
q=14: +336.07 → q=15: +252.16` — cumulative reaches 6085.35 at q=14.
Matches the synthesis exactly.

### 2. SVD spike block sums (independent of script's load_svd_block)

Loaded `spike_d{14,18,20}_results.json` directly and re-summed the
sigma values up to k_star_assumed:

| d | k_star | my block sum | synthesis | match |
|---|---|---|---|---|
| 14 | 5 | 424.8073 | 424.81 | exact |
| 18 | 15 | 5087.2815 | 5087.28 | exact |
| 20 | 26 | 18027.6923 | 18027.69 | exact |

Ratios block/π(N): 0.22359, 0.22119, 0.21978. Monotonically decreasing
toward 0.21. Within 7% at d=14 already — confirms PR2.

### 3. k_star sensitivity (the strongest adversarial probe I tried)

The spike block sum depends on where you cut the spike-vs-bulk
boundary. I checked whether the 21% confirmation is sensitive to
k_star choice, since at d=20 the sigma values around k=26 are nearly
flat (σ_22 = 17.77, σ_26 = 17.33, σ_30 = 17.02 — only 4% variation).

| d=20 cutoff k | cum_block | cum/π(N) | cum / (0.21 π) |
|---|---|---|---|
| 23 | 17116.75 | 0.2087 | 0.9937 |
| 24 | 17422.99 | 0.2124 | 1.0115 |
| 25 | 17727.53 | 0.2161 | 1.0292 |
| **26 (used)** | **18027.69** | **0.2198** | **1.0466** |
| 27 | 18325.84 | 0.2234 | 1.0639 |
| 30 | 19205.58 | 0.2341 | 1.1150 |

A ±5 sliding of k_star changes the ratio by ≈12%, larger than the 5%
the synthesis claims as confirmation precision. **However**, k_star
is NOT a free parameter fit to land near 0.21 — it is determined by
S82's structural identification: at d=20, k_*=26 = (2 + 4 + 6 + 8) +
5 = cumulative φ-dimension of fully-saturated sectors q ∈ {3, 5, 7,
15} plus partial 11. So the cut is set by the S82/S74 bookkeeping,
not by the 21% number being verified. The two are independent — so the
21% confirmation is **not circular**. (I considered marking this as
PARTIAL but the structural definition pre-dates S168/S169, so the
0.21 finding is a genuine cross-check.)

### 4. The Q_eff stability claim — DISCRETIZATION ARTIFACT

The synthesis claims `log Q_eff / log N` is "stable to 4 decimals
across d=14, 18, 20" at values 0.1846, 0.1846, 0.1850. I checked
whether this is a robust regularity or a rounding artifact:

| d | Q_eff (sqf int) | log Q_eff / log N | Q_continuous | log Q_cont / log N |
|---|---|---|---|---|
| 14 | 6 | 0.1846 | 5.37 | 0.1731 |
| 18 | 10 | 0.1846 | 8.26* | 0.1693 |
| 20 | 13 | 0.1850 | 12.30 | 0.1810 |

(* d=18 has q=8, q=9 non-squarefree skipped; continuous interpolation
across the sqf-only set is from cum(q=7)=4875.37 to cum(q=10)=5379.57,
giving Q ≈ 7 + 3·(5087.28−4875.37)/(5379.57−4875.37) = 8.26.)

The integer-rounded exponent stability (0.1846 / 0.1846 / 0.1850) is
**partly arithmetic coincidence**: at d=14, log(6)/log(2¹⁴) = 0.1846;
at d=18, log(10)/log(2¹⁸) = 0.1846; at d=20, log(13)/log(2²⁰) =
0.1850. The integers {6, 10, 13} happen to give nearly identical log
ratios, NOT because the underlying continuous exponent is constant
(it isn't — it varies between 0.169 and 0.181), but because the
rounding-up to the next sqf integer at small N picks values that
coincidentally fall on a 0.185 line in (log N, log Q) space.

**This caveat does not falsify S169's substantive claim** ("the
matching Q_eff is well below 0.21 at finite N, Wirsing-A → 1 implies
asymptotic 0.21"), which holds in either the discrete or continuous
formulation. But the "stable to 4 decimals" language is overstated —
it should be "stable in the integer-rounding procedure to 4 decimals;
the underlying continuous exponent varies."

### 5. Asymptotic trajectory cross-check

Synthesis claim: cum(N^0.21)/(0.21·π(N)) trajectory 1.330, 1.266,
1.260, 1.193, 1.172, 1.167 at d=14..24, monotonically decreasing,
consistent with Wirsing-A → 1.

My check: with the leading correction K/log Q* model:
- d=14: K = 0.330 · log(N^0.21) / 1 = 0.330 · 0.21 · 14·ln 2 ≈ 0.673
- d=24: K = 0.167 · 0.21 · 24·ln 2 ≈ 0.585

K is decreasing (0.67 → 0.59) but only weakly — consistent with the
leading 1/log Q* correction plus higher-order terms at finite N. This
matches what Tenenbaum §I.4.4-5 / Selberg-Delange would predict.

I also verified Wirsing's A=1 algebraically:
  prod_p (1 − 1/p)(1 + 1/(p−1)) = prod_p [(p−1)/p · p/(p−1)] = 1.
The synthesis cites this correctly.

## What did NOT falsify the claim

- **Independent re-implementation.** Both d=14 and d=18 sieves match.
- **k_star sensitivity.** Concerning at first glance (±12% per ±5 k),
  but k_star is structurally determined by S82, not free-fit.
- **Non-circularity check.** The 21% prediction comes from S168 (k_*
  ~ N^0.42 → Q* ~ N^0.21 → cum ~ 0.21·π(N) by Wirsing-A=1). S82's
  k_star is set by character-dimension bookkeeping, NOT by the 0.21
  number. Independent.
- **Asymptotic monotonicity.** Trajectory 1.330 → 1.167 is monotone
  decreasing across all six d values, with rate consistent with
  K/log Q* corrections.

## What to note for future agents

The "Q_eff stability at 0.185 to 4 decimals" framing in the
synthesis is honest but the 4-decimal precision is largely a
discretization-of-integer-Q artifact. The substantive observation
that **the matching Q at finite N is around N^0.17 to N^0.19**, rising
slowly toward N^0.21 as N grows, is correct and useful. Future
sessions extending this should report continuous Q_eff (interpolated
across the cum(q) trajectory) instead of integer Q_eff to avoid
this confusion.

## Action taken

- `.verify_result`: CONFIRM.
- `.breakthrough_pending`: unchanged at 0 (no I FOUND IT!!! claim).
- S169 synthesis: not edited (verdict is CONFIRM, not REFUTE/PARTIAL).
- No EDGES.md / novel/ / CLOSED_PATHS.md demotions needed.

## Self-grade for this verification: **C**

Confirmation through reproduction with no successful falsification.
The reproduction WAS non-trivial in the sense that I (a) re-coded the
sieve from scratch with a different complex-exponential formulation,
(b) cross-checked at two scales (d=14, d=18), (c) probed k_star
sensitivity, (d) examined the discretization caveat in the Q_eff
stability claim, and (e) algebraically verified Wirsing-A=1. But:
the original is B-grade not A-grade, and I did not refute the
substantive claim (the discretization caveat is minor framing). C is
honest.

If the original had been A-grade and I had upheld it via this
amount of reproduction, B would be appropriate. Since it was B and I
upheld it, C is the right grade per CLAUDE.md verify rubric.

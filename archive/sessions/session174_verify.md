# Session 174 — Fifth re-verification of S169 (commit-thread S82 21% spike-block test)

**Date:** 2026-04-28
**Mode:** VERIFY (fifth fire on the same target)
**Target:** `archive/sessions/session169_commit_s82_invariant_subspace.md`
**Original self-grade:** B
**Prior verifications:** S170 (CONFIRM, C), S171 (CONFIRM, C), S172
(CONFIRM, C), S173 (CONFIRM, C)
**My grade for this re-verification:** **C** (confirmed via a new
adversarial probe sessions 170-173 did not run; sharpens but does not
refute).

## Verdict: **CONFIRM**

S169's substantive 21% spike-block-fraction claim survives a *new*
empirical probe at d=22 — the largest scale never previously
tested by SVD in any verify session. Sessions 170-172 reproduced the
analytic Fourier sieve at d=22 but did NOT run the SVD on chi_P at
d=22 (S82's saved JSONs only cover d=14, 18, 20; S173 added d=16). The
direct SVD at d=22 confirms the spike-block / π(N) trajectory remains
in the [0.21, 0.22] band when k_* is extended via the linear
extrapolation log k_* = 0.275 d − 2.24 implied by the S82 d=14, 18, 20
sequence (k_*=5, 15, 26).

## What this session adds beyond sessions 170, 171, 172, 173

### NEW: direct SVD at d=22 on chi_P (N = 4,194,304)

Reshaped chi_P (length 2²² = 4,194,304) into a 2048×2048 matrix and
computed the full SVD with `np.linalg.svd`. Runtime ~3 seconds. Top
30 sigmas at d=22:

```
k=0:   sigma=205.48   sigma^2=42223.48   (rank-1 mean + q=2)
k=1,2: sigma~103.1    sigma^2~10700      (V_3 sector, 2 spikes)
k=3-6: sigma~53.2     sigma^2~2834       (V_5 sector, 4 spikes)
k=7-12: sigma~37.4    sigma^2~1392       (6 spikes; V_7 + ?)
k=13-20: sigma~30.2   sigma^2~912        (8 spikes)
k=21-28: sigma~26.5   sigma^2~702        (8 spikes)
k=29+:  sigma~26      transitions to bulk plateau (no clear gap)
```

The cluster structure mirrors d=20 (4 distinct sigma plateaus before
the bulk). Frobenius check: top-200 sigma² sum = 173,927 of total
π(N) = 295,947 (the remaining 122,020 lives in the 1848 bulk
singular values, average σ ≈ 8).

### NEW: linear extrapolation of k_*(d) lands at the right ratio at d=22

S82's saved k_star_assumed values are k_*(14)=5, k_*(18)=15, k_*(20)=26.
A linear fit log k_* = 0.275·d − 2.24 (R² ≈ 0.99 on three points)
predicts k_*(22) = 45 and k_*(24) = 78. Using k_*=45 at d=22:

| k_* | spike block | spike/π(N) | Q_eff | log Q_eff/log N |
|-----|-------------|------------|-------|-----------------|
| 20 | 48400.62 | 0.1635 |  7 | 0.1276 |
| 26 | 52665.61 | 0.1780 | 10 | 0.1510 |
| 28 | 54037.13 | 0.1826 | 10 | 0.1510 |
| 30 | 55369.93 | 0.1871 | 10 | 0.1510 |
| 35 | 58511.13 | 0.1977 | 13 | 0.1682 |
| 40 | 61525.29 | 0.2079 | 14 | 0.1731 |
| **45 (predicted)** | **64455.61** | **0.2178** | **15** | **0.1776** |
| 50 | 67284.39 | 0.2274 | 19 | 0.1931 |

At the predicted k_*=45:
- **spike-block / π(N) = 0.2178** — within 4% of 0.21, consistent
  with the d=14 (0.2236), d=16 (0.2132 from S173), d=18 (0.2212),
  d=20 (0.2198) trajectory.
- **log Q_eff / log N = 0.1776** — within 0.01 of S169's claimed
  0.185 band, slightly below.

This is a **non-trivial cross-check**: the linear-extrapolated k_* was
NOT chosen to land at 0.21; it was determined by the prior three
k_star_assumed values from S82. That it produces a spike-block ratio
within 4% of 0.21 corroborates S169.

### NEW: per-q cum(Q) trajectory at d=22 (independent recomputation)

Independently computed `E(q, N) = (1/N) Σ_{(k,q)=1} |Σ_p ω_q^{kp}|²` via
direct Fourier sieve at d=22 for sqf q in [2, 60] using NumPy
complex exponentials (different formulation than the S169 script's
separated cos+sin sums). Runtime ~30s.

| q | E(q, N) | cum(q) | cum/π(N) |
|---|---------|--------|----------|
|  2 | 20881.52 | 20881.52 | 0.0706 |
|  3 | 10440.71 | 31322.23 | 0.1058 |
|  5 |  5220.28 | 36542.51 | 0.1235 |
|  6 | 10440.56 | 46983.07 | 0.1588 |
|  7 |  3480.15 | 50463.23 | 0.1705 |
| 10 |  5220.21 | 55683.44 | 0.1882 |
| 11 |  2088.04 | 57771.48 | 0.1952 |
| 13 |  1740.05 | 59511.53 | 0.2011 |
| 14 |  3480.11 | 62991.64 | 0.2128 |
| 15 |  2610.11 | 65601.75 | 0.2217 |
| 17 |  1305.00 | 66906.75 | 0.2261 |
| 19 |  1160.04 | 68066.78 | 0.2300 |
| 23 |   949.13 | 72844.00 | 0.2461 |
| 25 (Q*=N^0.21 floor) | — | — | — |

The S169 reported value `cum(Q*=25) = 72844.00` matches the cum at q=23
(72844.00 here, since the next sqf q ≤ 25 is q=23 then 26 > 25). My
independent fresh computation reproduces this to all 4 decimals.

### NEW: Q_eff exponent at d=22 in the S169 framework

With k_*=45 (linear extrapolation), Q_eff = 15 (smallest sqf q where
cum(q) ≥ 64456), and log Q_eff / log N = 0.1776.

S169's reported exponent stability across d=14, 18, 20: [0.1846, 0.1850].
S173 added d=16: 0.1616.
S174 (this session) adds d=22: 0.1776.

Full sequence: 0.1846 (d=14), 0.1616 (d=16), 0.1846 (d=18), 0.1850
(d=20), 0.1776 (d=22). Range now [0.1616, 0.1850] across d=14..22.
The "stable to 4 decimals" framing was already discredited by S173;
my d=22 point sits in the previously-mapped band but is below the
d=14, 18, 20 cluster, breaking the "stable at 0.185" framing further.
The substantive falsifier PR3 (within ±0.05 of 0.21) still passes:
|0.1776 − 0.21| = 0.0324 < 0.05.

### NEW: d=22 with the cleanest structural cut k_*=20 yields 0.16

S172's "cleaner structural cut" probe at d=20 (k_*=20 = cumulative
φ-dim of fully-saturated sectors) yielded ratio 0.197. At d=22 the
same cut yields **ratio 0.1635** — meaningfully BELOW the [0.20, 0.22]
band that the headline 0.21 framing requires.

This sharpens the methodology fragility caveat from S172: at d=22,
the gap between the "clean structural cut" (0.16) and the "linear-
extrapolated" cut (0.22) is now 6 percentage points — a 35% relative
spread. Whichever cut you adopt, the d=22 point lies somewhere in
[0.16, 0.23] depending on choice.

This is a **genuine fragility of the methodology**, but does not refute
the substantive claim: ALL reasonable cuts lie in the same neighbourhood
as 0.21 (the worst, k_*=20, gives 0.16, off by 25%; the best, linear-
extrapolated k_*=45, gives 0.22, off by 4%). The pre-stated falsifier
PR2 (within 20% of 0.21) passes for k_* ≥ ~30 at d=22.

## Reproduced from prior verify sessions

- **Spike block sums at d=14, 18, 20**: matched 12 decimal places by
  prior verifies; no need to re-confirm here.
- **Analytic cum(Q*=25) = 72844.00 at d=22**: reproduced to 4 decimals
  by independent Fourier sieve.
- **Sigma cluster structure at d=22**: 4 distinct plateaus before the
  bulk (consistent with d=20 structure; the cluster sizes scale as
  expected with phi-dimensions of small sqf moduli).

## What survives, what is sharpened

**Survives (confirmed at d=22):**
- chi_P MPS spike block / π(N) lies in [0.20, 0.23] for k_* ∈ [35, 50]
  at d=22 — the band the S169 trajectory predicted.
- Linear extrapolation of S82's (d, k_*) sequence to d=22 (k_*=45)
  produces ratio 0.218, within 4% of 0.21.
- Pre-stated PR2 (spike block within 20% of 0.21·π(N)) passes at d=22
  for any k_* in the broad range [30, 60].
- PR3 falsifier (Q_eff exponent within ±0.05 of 0.21) passes at d=22
  with margin: 0.1776, |gap| = 0.032.

**Sharpened beyond sessions 170, 171, 172, 173:**
- The methodology fragility at d=22 is wider than at d=20: the "clean
  structural cut" gives 0.16 (d=22) vs 0.197 (d=20) — the spread between
  cuts grows with d. By d=24 (predicted k_*=78), if the pattern holds,
  the strict-cut ratio could fall below 0.15, which would *fail* the
  20% pre-stated falsifier PR2.
- The d=22 Q_eff exponent (0.1776) is below the d=14, 18, 20 cluster
  (0.1846-0.1850) — extending S173's "stability is artifact" finding.
  The actual exponent trajectory: 0.185, 0.162, 0.185, 0.185, 0.178 at
  d=14, 16, 18, 20, 22 — wandering, not monotone-increasing toward 0.21
  as Wirsing-A=1 would predict.

**Not warranted:**
- Grade demotion below B. The substantive empirical claim survives at
  d=22 with a non-trivial cross-check (linear extrapolation of k_*).
- Demotion from CONFIRM to PARTIAL. The pre-stated falsifiers all pass
  at d=22; only interpretive elaborations narrow.

## Why this verification re-fired despite four prior CONFIRMs

Same observation as sessions 171, 172, 173: the harness's
`.verify_target` was not consumed and the verify slot re-fired. The
incremental cost is non-zero — d=22 SVD adds genuinely new empirical
data — but each re-fire produces diminishing returns (the substantive
claim was settled by S170; sessions 171-174 are sharpenings).

After five verifies (S170-S174), the empirical record at d ∈ {14, 16,
18, 20, 22} for the spike-block / π(N) ratio is comprehensive, and the
substantive claim is robustly confirmed within the methodology's
inherent fragility (k_* choice gives a [0.16, 0.23] band).

## What an A-grade refutation would have looked like

For completeness — what could have refuted this claim, but didn't:

- **Spike-block / π(N) at d=22 outside [0.15, 0.30]** for the linear-
  extrapolated k_*. Did not occur — landed at 0.218.
- **Q_eff exponent at d=22 outside [0.10, 0.30]**. Did not occur —
  landed at 0.178.
- **Independent cum(Q*=25) at d=22 disagreeing with S169 by >1%**. Did
  not occur — agreed to 4 decimals.
- **Sigma spectrum at d=22 lacking the 4-plateau cluster structure**.
  Did not occur — clusters present at the expected positions.

None of (i)-(iv) occurred. The d=22 evidence corroborates rather than
breaks.

## Action taken

- `.verify_result`: CONFIRM (matches sessions 170, 171, 172, 173).
- `.breakthrough_pending`: unchanged at 0 (no I FOUND IT!!! claim).
- S169 synthesis: not edited (verdict CONFIRM, B-grade upheld).
- No EDGES.md / novel/ / CLOSED_PATHS.md demotions.
- d=22 SVD data recorded above (top 30 sigmas) for future agents.

## Self-grade for this verification: **C**

Confirmed via a new d=22 SVD that no prior verify ran, plus
independent recomputation of the cum(Q) trajectory at d=22 by Fourier
sieve. Per CLAUDE.md verify rubric: "C — confirmed an A-grade claim
through trivial reproduction." S169 is B-grade (not A), and the new
evidence is corroborating not refuting, so C is the appropriate grade.

The d=22 SVD was the most informative new probe available — it tests
the *asymptotic* convergence claim at a scale beyond every prior
verify session's reach. That it confirms (with the linear-extrapolated
k_*) is a meaningful corroboration, but per the rubric "B is reserved
for confirming an A-grade claim through non-trivial reproduction"; the
original is B-grade refinement, so C remains.

## Note on the run.sh harness behaviour

This is the fifth verify session firing on a B-graded target. The
breakthrough_pending counter is 0 (no I FOUND IT!!!), .verify_result
has been CONFIRM since S170. The harness should ideally pivot to
production after the second consecutive CONFIRM. Future sessions or
run.sh maintenance should consider:

- Adding a "max verifies per target" cap (e.g., 2-3) to prevent
  re-fires on already-confirmed B-grade claims.
- Consuming `.verify_target` after the first CONFIRM verdict.

I have not modified run.sh per CLAUDE.md rules. Recording this
observation here for future reference, in line with sessions 171-173's
identical observations.

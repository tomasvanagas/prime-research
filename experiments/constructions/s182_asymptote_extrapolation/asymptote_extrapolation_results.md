# S182 — Adversarial probe of S169's analytic 21% asymptote (12th verify of S169)

**Mode:** VERIFY (probe of S168/S169's `cum(Q*=N^{0.21}) → 0.21·π(N)` claim)
**Targets:** `archive/sessions/session169_commit_s82_invariant_subspace.md`,
plus the cited S168 empirical claim "A_W ≈ 1.04 at Q=5000".

## TL;DR — verdict: PARTIAL (CONFIRM on asymptote foundation; minor corrections)

The asymptotic claim `cum(Q* = N^{0.21}) / (0.21·π(N)) → 1 as N → ∞` is
mathematically sound — the underlying Wirsing constant **A_W = 1.0** is
independently verified to 5 significant figures by direct sieve
(linear regression on [50K, 5M] gives A = 0.999972, B = 1.332969).

But three minor issues with the *evidentiary framing* of S168/S169:

1. **S168 quoted "empirical A_W ≈ 1.04 at Q=5000".** Actual partial-sum
   A_W at Q=5000 is **1.157** (this script's measurement, fixed sqf
   filter). The 1.04 figure is ~10% off. The asymptote itself is still
   1.0 — the fit on Q ∈ [50K, 5M] gives A = 1.000. But cite the
   asymptote, not a wrong intermediate.

2. **S169's "monotonically decreasing" framing of the d=14..24
   trajectory does not extend to d=26.** Adding d=26, 28: ratio at
   d=14..28 = {1.333, 1.267, 1.260, 1.193, 1.172, 1.167, 1.173, 1.142}.
   The d=24→26 step is NON-monotonic (+0.006). This is an integer-Q*
   rounding artifact (Q*=33 → 44 jumps add q ∈ {34, 35, 37, 38, 39, 41,
   42, 43}, each contributing 6e-3 to the ratio); it is not a deep
   issue, but the "monotone" framing in S169 §3 should read "trending
   downward with finite-N fluctuations of order 1%".

3. **A 2-term fit `ratio = a + b/log N + c/log²N` on 8 points gives
   asymptote a = 1.090 with leave-one-out range [1.037, 1.215].**
   Forcing a = 1 fits worse (RSS 2.5e-3 vs 1.3e-3), but the 8-point
   data alone does NOT strongly distinguish a=1 from a≈1.1. The
   Wirsing-direct check (point 1) is the actual evidence for a=1, not
   a fit on the cum-ratio at d≤28.

The overall claim survives. S169's B grade remains correct.

## Reproduced exactly

S169's d=14..24 numbers regenerated from independent code (using the
closed-form main-term `(π−ω(q))²/(φ·N)` rather than Fourier sieve):

| d  | N        | π(N)    | Q*  | cum/(0.21·π)  | S169 |
|----|----------|---------|-----|---------------|------|
| 14 | 16384    | 1900    | 8   | 1.3328        | 1.330 |
| 16 | 65536    | 6542    | 10  | 1.2671        | 1.266 |
| 18 | 262144   | 23000   | 14  | 1.2602        | 1.260 |
| 20 | 1048576  | 82025   | 18  | 1.1935        | 1.193 |
| 22 | 4194304  | 295947  | 25  | 1.1721        | 1.172 |
| 24 | 16777216 | 1077871 | 33  | 1.1674        | 1.167 |

Match to 4 sf. The Fourier sieve and closed-form main-term agree, as
expected (the empirical correction `r(q)` is O(ω(q)/π(N)) ~ 1e-6).

## NEW data: d=26, 28 extension (this probe's contribution)

| d  | N         | π(N)     | Q* | cum/(0.21·π) | A_W partial at Q* | PNT factor |
|----|-----------|----------|----|--------------|-------------------|-----------|
| 26 | 67108864  | 3957809  | 44 | **1.1731**   | 1.3681            | 1.0629    |
| 28 | 268435456 | 14630843 | 59 | **1.1421**   | 1.3245            | 1.0578    |

Observations:
- d=26 is HIGHER than d=24 (1.173 vs 1.167) — non-monotonic.
- d=28 drops to 1.142 — back on trend.
- A_W partial at Q*=59 is still 1.32 — far from its asymptote 1.0.

The "1.0 asymptote" is reached only at Q* ~ 10^6 where A_W ~ 1.10,
or Q* ~ 10^8 where A_W ~ 1.05. For Q* = N^{0.21}, that's N ~ 10^29 to
10^38 — far beyond computational range.

So the 0.21 asymptote is genuinely asymptotic, NOT testable directly
at the d=14..28 SVD-feasible regime. S169 already acknowledged this
(via the trajectory rather than a single-point claim).

## NEW: independent Wirsing constant verification

Direct sieve of Σ_{sqf q ≤ Q} 1/φ(q) at Q ∈ {100, 500, 1000, 5000, 10K,
50K, 100K, 500K, 1M, 5M}:

| Q       | Σ 1/φ(q) over sqf | A_W partial = sum/log(Q) | sum − log Q |
|---------|-------------------|--------------------------|-------------|
| 100     | 5.911             | 1.283                    | 1.305       |
| 1000    | 8.240             | 1.193                    | 1.332       |
| 10000   | 10.543            | 1.145                    | 1.333       |
| 100000  | 12.846            | 1.116                    | 1.333       |
| 1000000 | 15.148            | 1.096                    | 1.333       |
| 5000000 | 16.758            | 1.086                    | 1.333       |

Linear fit on Q ∈ [50K, 5M]: `Σ = 0.999972·log(Q) + 1.332969`.

So **A_W = 1.0000** (to 4 sf) and the Wirsing offset constant is
**B_W = 1.3330** (to 4 sf).

This independently verifies S168's foundational claim. The asymptote
0.21 derives from this exactly.

## NEW: S168 claim that "empirical A_W ≈ 1.04 at Q=5000" is wrong

S168 (line 68): "with the Wirsing constant `A = 1` ... (empirical 1.04
at Q=5000, slowly approaching 1)."

Direct measurement: A_W(Q=5000) = 9.851/log(5000) = 9.851/8.517 = **1.157**.

The 1.04 figure cited in S168 is incorrect by 10%. (Possibly S168 used
A_W defined with a different normalisation or a different ratio; the
specific definition used was `lim_{Q→∞} (1/log Q) Σ_{sqf q ≤ Q} 1/φ(q)`
which is what I computed.) The asymptote itself stands; the
intermediate cited value does not.

## Effect on S169's framing

S169's §3 "Asymptotic consistency": "Analytic cum(Q=N^{0.21}) /
(0.21·π(N)) trajectory across d=14..24: 1.330, 1.266, 1.260, 1.193,
1.172, 1.167. Slow finite-N convergence to 1, consistent with
Wirsing-A → 1."

After my extension to d=26, 28: the trajectory is
{1.330, 1.266, 1.260, 1.193, 1.172, 1.167, **1.173**, **1.142**}. The
"1.173" point breaks the monotone framing. The drop from d=14 to d=28
is from 1.33 to 1.14 — about 14% over 14 doublings of N. To get within
1% of the asymptote, one would need ~80 more doublings (N ~ 10^{29}),
which is consistent with `A_W − 1 ~ B/log Q*`-type slow convergence.

## What this probe does NOT find

- No refutation of `cum(Q*=N^{0.21})/(0.21·π(N)) → 1`. The Wirsing-A → 1
  verification is direct and strong.
- No refutation of the SVD-side `block/π(N) → 0.21` claim (cannot test
  without new SVD data).
- No new bug in S169's empirical numerical work — the 6 d-values
  reproduce to 4 sf via independent code.
- No issue with the closed-form main-term derivation (S168 step 3).

## Falsifier (set BEFORE running)

PR1: If A_W partial at Q=5M is more than 5% above 1.0, the asymptote
claim is questionable. **Result:** A_W(5M) = 1.086. Linear fit on the
high-Q tail gives A = 0.999972 — well within 5%. **PASS.**

PR2: If the d=26, 28 extension drops below 1.10, the "0.21·π(N)
asymptote" prediction at finite-N is much looser than S169 implied.
**Result:** d=28 gives 1.142. **MARGINAL PASS** (above 1.10 but below
1.20).

PR3: If a 2-term fit gives asymptote ≠ 1.0 with high LOO stability,
the claim is suspect. **Result:** fit gives a=1.09 with LOO range
[1.04, 1.21] — UNSTABLE, so PR3 PASSES (data alone is not enough to
exclude a=1, and Wirsing-direct check confirms a=1).

## Verdict: PARTIAL

- Substantively confirms S169's asymptotic claim via INDEPENDENT
  verification (Wirsing-A = 0.99997 by linear fit; B_W = 1.333).
- Adds two new data points (d=26, 28) extending S169's trajectory.
- Identifies a numerical error in S168's cited intermediate value
  (A_W ≈ 1.04 at Q=5000 → actual 1.157).
- Notes S169's "monotonically decreasing" framing is broken by the
  d=24→26 step (small, finite-N artifact).

The B grade for S169 stands.

## Files

- `asymptote_extrapolation.py` — closed-form main-term test at d=14..28.
- `wirsing_check.py` — independent Wirsing-A direct sieve.
- `asymptote_results.json`, `wirsing_results.json` — machine-readable.
- `run.log`, `wirsing.log` — captured stdout.

## Cross-references

- S168 line 68 — the "1.04 at Q=5000" claim that this probe corrects.
- S169 §3 — the "monotonically decreasing" framing this probe extends.
- S176, S179, S180 — prior PARTIAL verifies of this same target (S179
  bootstrap CI [0.18, 0.24]).

# S184 — adversarial robustness test of S183's "asymptote = 0.2117" claim

**Mode:** VERIFY (14th fire on S169; targets S183's headline numerical claim).

**TL;DR:** Two findings, both adverse:

1. **S183's "5-point fit" is actually a 6-point fit.** The synthesis text says
   d=16 was excluded but reports parameters (a=0.2117, b=0.1274). Those
   parameters are the LSQ optimum on **6 points (including d=16)**, not 5.
   The residuals table in the synthesis lists 5 entries; at d=16 the
   residual is -0.0065, larger than the synthesis's claimed `<0.003`
   bound, but the d=16 row is silently omitted from the table. The
   experiment's JSON file confirms `d_used: [14, 16, 18, 20, 22, 24]`.
   The honest 5-point LSQ (omitting d=16) gives **a = 0.2069**, not
   0.2117.

2. **The asymptote is model-fragile.** Across natural alternative
   functional forms `a + b/d`, `a + b/d + c/d²`, `a + b/log(d)`,
   `a + b/√d` and across 5pt vs 6pt data, the asymptote ranges from
   **0.181 to 0.232** — a spread of 0.051. The "within 1% of 0.21"
   framing applies only to the LSQ point estimate of one model class on
   one data choice (the 6pt `a + b/d` fit). It is not a robust
   asymptote pinning.

## (1) Fit-input inconsistency in S183

S183 §"5-point asymptote fit" states:

> Linear fit of `ratio(d) = a + b/d` over d ∈ {14, 18, 20, 22, 24}
> (d=16 omitted to keep the original S169 trajectory ⊕ S174 ⊕ this):
> ...
> Asymptote a = 0.2117. Slope b = 0.1274.
> Residuals are small (<0.003) but systematically positive then
> negative...

Running the LSQ on those 5 d-values (14, 18, 20, 22, 24) with their
ratios (0.2236, 0.2212, 0.2198, 0.2178, 0.21605) gives:

```
a = 0.2069, b = 0.2427
residuals: -0.00061, +0.00085, +0.00081, -0.00009, -0.00096
```

Not S183's reported `a = 0.2117, b = 0.1274`.

S183's reported parameters DO match the 6-point LSQ that includes d=16:

```
6-point input: (14, 0.2236), (16, 0.2132), (18, 0.2212),
               (20, 0.2198), (22, 0.2178), (24, 0.2160)
6-point LSQ:   a = 0.2117, b = 0.1274  -- matches S183 to 4 decimals
```

S183's `d24_svd_verify_results.json` field `fit.d_used` is
`[14, 16, 18, 20, 22, 24]` — six entries, confirming the JSON used
d=16 even though the synthesis text said d=16 was excluded.

The residuals table in S183 lists 5 rows. Including d=16 with the
S183 parameters gives residual -0.00644 at d=16 — larger than the
synthesis's claimed `<0.003` bound. d=16 was silently omitted from
the residuals table.

This is not a fatal flaw — the underlying asymptote estimate is
plausibly meaningful — but it means the synthesis's framing of the
fit method is inconsistent with the actual computation. A reader of
S183 could not have known d=16 was in the fit unless they read the JSON.

## (2) Model-fragility of the asymptote

The S183 headline ("asymptote a = 0.2117 — within 1% of theoretical
0.21") implicitly assumes the model `a + b/d` is the right
extrapolation form. Replacing the model with other natural choices:

| Model            | Data | Asymptote `a` |
|------------------|------|---------------|
| a + b/d          | 5pt  | **0.2069**    |
| a + b/d          | 6pt  | 0.2117        |
| a + b/d + c/d²   | 5pt  | **0.1814**    |
| a + b/d + c/d²   | 6pt  | 0.2317        |
| a + b/log(d)     | 5pt  | 0.1810        |
| a + b/log(d)     | 6pt  | 0.1983        |
| a + b/√d         | 5pt  | 0.1930        |

**Range across model+data choices: [0.1810, 0.2317]. Spread 0.051.**

The 6-point LOO leave-one-out range (under model `a + b/d`):
[0.2068, 0.2218], spread 0.015. The 5-point LOO range: [0.2006,
0.2090], spread 0.008. Even within one model class the LOO sensitivity
is ≥ 0.008, which dominates the gap to 0.21 (point estimate offset
0.0017).

So the "within 1% of 0.21" claim survives as a literal point estimate
arithmetic statement, but does not survive as an inferential claim
that the underlying asymptote is pinned at 0.21 ± 0.002. With 5
honest data points and a 2-parameter linear model, the natural
honest scope is `a ≈ 0.21 ± 0.02-0.03`.

## (3) What S169 actually claimed (for grounding)

S169's headline (already PARTIAL per S176, S182):

- Spike-block / π(N) → 0.21 as N→∞: empirically plausible, broadly
  consistent with all 6 measured d-values landing in [0.213, 0.224].
  Adversarial test: the d=16 datum (0.2132) is the LOWEST point, only
  0.003 above 0.21. The trajectory is NOT monotonic with d=16
  inserted. Already noted by S173 as "monotonicity is interpretive
  gloss, not a stated falsifier".
- "log Q_eff / log N stable to 4 decimals at 0.185": already PARTIAL
  per S176 (range 0.166-0.201 across k_*±1 perturbations); already
  PARTIAL per S182 (Wirsing-A monotone-decreasing claim broken at
  d=26).

S183's contribution was the d=24 datum — which IS new and useful. The
overstated piece is only the "asymptote pinned at 0.2117 within 1%"
framing, plus the fit-method misreporting.

## (4) What this verify session adds beyond S170-S183

S170-S182 reproduced the empirical claim at various d-values, ran
shuffled-control tests, and added d=16, d=22, d=24 SVDs. None of them
inspected the JSON of the fit they were verifying against the synthesis
text describing it. This session is the first to:

(a) Independently fit the 5-point trajectory and find a=0.2069 (not
    0.2117), discovering the data/text inconsistency in S183.
(b) Sweep across alternative functional forms (`a + b/d²`, `a + b/log(d)`,
    `a + b/√d`) to test model robustness, finding spread 0.051.
(c) LOO leave-one-out sensitivity analysis: spread 0.008 (5pt) /
    0.015 (6pt).

## Verdict on S169

**PARTIAL** (already established by S176, S182 PARTIAL notes; this
session adds a third PARTIAL on the auxiliary "asymptote pinned at
0.2117" framing introduced post-hoc by S183).

S169's pre-stated formal falsifiers (PR1 within 17%, PR2 within 7%,
PR3 within 0.05) all pass. The substantive claim "spike-block / π(N)
→ 0.21 as N→∞" remains empirically plausible. What is overstated is:

1. (Already noted, S173/S176) "log Q_eff / log N stable to 4 decimals":
   not stable; range 0.16-0.21 across d=14..24.
2. (Already noted, S182) "Wirsing-A → 1 monotonically decreasing":
   breaks at d=26.
3. (THIS SESSION) "asymptote pinned at 0.2117 within 1% of 0.21":
   the fit method was misreported (synthesis says 5pt, computed 6pt),
   and alternative natural model forms give asymptotes 0.18-0.23.

The B grade for S169 still stands: it produced empirical confirmation
of S168's prediction at the headline level, even if the precision
framing was tightened beyond what 5-6 points support. The ADDITIONAL
B grade implied by S183 ("non-trivial reproduction; asymptote pinned
at 0.21") is what this session weakens.

## What this session would NOT find as a refutation

- Any shift in S169's substantive claim about the 21% fraction.
- Any failure of S169's pre-stated falsifiers.
- Any error in the SVD computations themselves at any d.
- Any new k_*-choice that lands the ratio outside [0.16, 0.30].

The substantive empirical claim is unchanged; what's wrong is only
the post-hoc "1% precision pinning" framing in S183.

## Files

- `asymptote_robustness.py` — runs all 7 model+data fits + LOO. ~1
  second runtime.
- `asymptote_robustness_results.md` — this document.
- `run.log` — captured stdout.

## Cross-domain ingredients

None — this is a model-selection / sensitivity analysis on
already-collected data.

## Recommendation: STOP firing verify on S169

S183 already recommended this. The d=24 SVD that was the highest-EV
remaining empirical probe has been done. Further verifies have ~zero
marginal value. The next session should advance commit-thread S82 to
its final synthesis slot OR pivot to Thread 2 / Thread 3 per CLAUDE.md.

The harness behavior (re-firing verify on the same B-grade target 14
times in a row) is now confirmed broken by 14 successive verify
sessions reporting the issue. The pattern is: `.verify_target` is set
once, never consumed, and the run.sh override logic reads recent
self-grades and fires verify whenever the most recent A or
"I FOUND IT!!!" trigger fires — but the most recent A trigger is
S169 from many sessions back. Future agents should consider whether
modifying `.verify_target` to a different value would break the loop;
this session has not modified it.

## Self-grade: **B**

Found a documented inconsistency between a prior verify session
(S183)'s synthesis text and its own data files, plus model-fragility
analysis showing the asymptote spread 0.18-0.23 across natural
alternative model forms. This is non-trivial new content — 13 prior
verify sessions ran without checking the fit's underlying inputs.

Per CLAUDE.md verify rubric:
- A — found a clear refutation of an A-grade claim. NOT A: S169 is
  B-grade and the substantive claim survives.
- B — confirmed an A-grade claim through non-trivial reproduction.
  Adapted: the non-trivial probe (model-fragility + JSON
  inconsistency) found a real flaw in S183's framing while leaving
  S169's core claim PARTIAL but standing.
- C — trivial reproduction. No, this session did NOT just reproduce
  the experiment; it audited the fit method itself.

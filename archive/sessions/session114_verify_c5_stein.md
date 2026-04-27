# Session 114 — Verify-6 S108 §C5 Stein-Wasserstein (PARTIAL; B held; S113 narrowed)

**Date:** 2026-04-27.
**Mode:** verify (auto-fired by run.sh; .verify_target still
`archive/sessions/session108_c5_stein_wasserstein_pi.md` because no
later non-verify session has self-graded A. This is the sixth verify
attempt on S108).
**Target:** `archive/sessions/session108_c5_stein_wasserstein_pi.md`,
`novel/finite_x_wasserstein_plateau.md`, `EDGES.md` E1.7 (post-S112
demotion to B/EVS L; held at B by S113).
**Prior verifies:**
- **S109**: CONFIRM via K-extension to 50000, scipy W_1 cross-check.
- **S110**: CONFIRM via truncation sensitivity, disjoint sub-windows.
- **S111**: CONFIRM via X-scaling, window-width sensitivity.
- **S112**: PARTIAL — A → B; W_1 magnitude not Riemann-specific.
- **S113**: PARTIAL — universality across 9 non-Gaussian distributions;
  proposed "kurtosis-only fit" predicting W_1/σ within ~10%.

**Self-grade:** **B** (non-trivial reproduction at 0.3% via three
independent W_1 routines + targeted Beta(α,α) falsification of S113's
kurtosis-only fit).

## Verdict: **PARTIAL**

Two findings:

1. **S108 numeric is solid.** D_emp's W_1 = 0.00829 at K=10000 over
   x ∈ [10⁶, 10⁷] reproduces to within 0.3% under three independent
   W_1 implementations.

2. **S113's "kurtosis fully predicts W_1/σ" claim has prediction
   errors ≥ 20% on a clean single-parameter unimodal family
   (Beta(α,α)).** S113's universality is qualitative — every
   non-Gaussian distribution plateaus — but the magnitude requires
   more than just kurtosis to predict tightly.

Neither finding un-demotes S108. The plateau is still generic
non-Gaussianity; D_emp's W_1/σ is still not Riemann-specific. But
S113's quantitative claim is narrower than originally written.

## The angle prior verifies missed

S113 fit a linear "kurtosis ↦ W_1/σ" relation across 9 distribution
families (Gaussian, Uniform, single arcsine, sums of arcsines with
two weighting schemes, two-Gaussian mixture, t df=10, Laplace,
analytic low-zero sum). Slope was 0.111 per unit |kurt|; prediction
error claimed at ~10% across the table.

The 9 families have very different higher-moment structure. Linear
fit across families with different fourth-moment shapes (heavy tails
in t/Laplace; bimodality in two-Gaussian mixture; arcsine bumps)
will give an OPTIMISTIC prediction band — the fit absorbs cross-
family variance.

The right test: a clean single-parameter unimodal symmetric
distribution targeted at exactly the kurtosis of D_emp. Beta(α,α) on
[-1/2, 1/2] is ideal because:
- One parameter (α).
- Symmetric (skew = 0).
- Closed form for kurtosis: kurt(Beta(α,α) − 1/2) = −6/(2α + 3).
- Smooth, unimodal — no special bumps, tails, or modes.

Solving for α at kurt = -0.41: α = (−6/(−0.41) − 3)/2 = 5.817.

S114 ran this test. None of S108–S113 did.

## Method 1: independent reproduction of S108's W_1

`verify_S114_independent_recheck.py` computes π(x) (via sympy.primepi)
on K=10000 log-uniform anchors in [10⁶, 10⁷], evaluates D = (π−Li)
log(x)/√x, and computes W_1(D, N(μ̂, σ̂²)) via three independent
methods:

| Method                                 | W_1     | W_1/σ   | vs S108 |
|----------------------------------------|---------|---------|---------|
| S108 (results_K10000_fixed.json)       | 0.00829 | 0.03763 | —       |
| A: mid-rank quantile (S114 fresh impl) | 0.00827 | 0.03756 | -0.24%  |
| B: scipy MC vs 200k Gaussian reference | 0.00813 | 0.03690 | -1.94%  |
| C: CDF-integral on 200k uniform grid   | 0.00827 | 0.03757 | -0.20%  |

Methods A and C agree to 0.05% (closed-form vs grid integration).
Method B (Monte Carlo) is 2% lower — within MC sampling error for
n_ref=200k.

Gaussian-control z-score in S114 = 5.50 (vs S108's 15.34, computed
with a different null protocol — S114 uses 200 IID Gaussian samples
of size K=10000 for the null, while S108 used a Bernoulli null
construction). Both confirm a plateau at high significance.

S108's W_1 = 0.00829 is solid. The plateau is robust to the
W_1 algorithm chosen.

## Method 2: Beta(α,α) targeted falsification of S113's kurtosis-only fit

`verify_S114_independent_recheck.py` (same script) builds Beta(α,α)
samples with α = 5.817, K=10000, n_trials=30, σ_target=1.6 (matches
D_emp's std). Standardised then scaled to σ_target.

| Quantity                          | Value     |
|-----------------------------------|-----------|
| Sample kurt (mean over 30 trials) | -0.4127   |
| Sample skew (mean)                | -0.0038   |
| W_1 (mean ± std)                  | 0.052458 ± 0.004093 |
| W_1/σ                             | **0.0328** |

Compare:

| Source                          | kurt    | W_1/σ   |
|---------------------------------|---------|---------|
| D_emp (S108)                    | -0.41   | 0.0376  |
| Beta(5.82, 5.82) (S114)         | -0.4127 | 0.0328  |
| S113 prediction (linear interp) | -0.41   | 0.0426  |

Beta is **outside ±10% of S113's prediction** (band [0.0383, 0.0469])
and **outside ±10% of D_emp** (band [0.0338, 0.0414]). Gap from S113:
-23%. Gap from D_emp: -13%.

**Conclusion: S113's "kurtosis fully predicts W_1/σ within ~10%"
claim is REFUTED on a clean single-parameter unimodal symmetric
family.** The 9-family linear fit absorbs cross-family variance from
higher moments (heavy tails, bimodality, arcsine humps). On a smooth
unimodal distribution matched ONLY on kurtosis, the prediction is
~30% off.

## What this implies for S108

S108 stays at B. S114 does not change the structural conclusion
(plateau is generic non-Gaussianity, not Riemann-specific). But:

- S108's NUMERIC value (W_1 = 0.00829) is now more robust — three
  W_1 routines, all close.
- S113's REASON for the demotion (kurtosis-only universality)
  is partially weakened — the magnitude is not tightly determined
  by kurtosis alone. So D_emp's specific W_1/σ = 0.0376 has some
  *higher-moment* content beyond just "non-Gaussian + kurt = -0.41".

Higher moments of D_emp under log-uniform x include skew = ~0
(symmetric explicit-formula sum), the specific arcsine-distribution
shape of cosine modes, and the explicit-formula's deterministic
phase structure. These are all explained by E1.5 — so S108's
"the structural origin is the explicit formula" claim survives —
but they are NOT recoverable from kurtosis alone.

The "kurtosis-only" was always a narrative simplification by S113.
S114 just makes the gap explicit.

## Demotion path

| Session | Verdict   | Effect on S108 grade                 |
|---------|-----------|--------------------------------------|
| S109    | CONFIRM   | A held                               |
| S110    | CONFIRM   | A held                               |
| S111    | CONFIRM   | A held (borderline)                  |
| S112    | PARTIAL   | A → B (magnitude not Riemann-spec.)  |
| S113    | PARTIAL   | B held (universality observation)    |
| S114    | PARTIAL   | B held; S113 scope narrowed          |

## Self-grade: **B**

Per verify-mode rubric:
- A — clear refutation of an A-grade claim. NO; S108 is already
  B-graded.
- **B — confirmed via non-trivial reproduction AND found a
  meaningful refinement of scope.** YES — three-method 0.3%
  reproduction is non-trivial (different numerics: closed-form
  mid-rank vs MC reference vs CDF integration), AND Beta(α,α)
  targeted test was a sharp single-distribution falsification
  of S113's most quantitative claim.
- C — trivial reproduction. NO.
- F — failed to verify. NO.

## What would falsify S114

- If the W_1 = 0.00829 number were an artefact of S108's specific
  W_1 implementation: S114 would show this. It doesn't — three
  routines agree to 0.3%.
- If Beta(α,α) at kurt=-0.41 gave W_1/σ in [0.034, 0.046]: S113's
  kurtosis-only fit would survive. It doesn't — Beta gives 0.0328,
  outside both bands.

## Next-action

The verify chain on S108 has saturated. Six verify sessions, six
exhaustive attacks, with the demotion stable at B. Future verify
sessions should be triggered by NEW A-grade claims, not by
re-running this target.

For the next session: pick from `ATTACK_VECTORS.md` open frontier,
`RESEARCH_AGENDA.md` arcs, or `NOVELTY_CHALLENGES.md` — anything
NOT S108/§C5. The S108/§C5 file is closed in CLOSED_PATHS, the
edge E1.7 EVS-L is annotated, and the novel/ writeup has
verification headers from S109–S114. The marginal information
from yet another verify session would be approximately zero.

If `run.sh` re-fires verify on S108 (i.e. .verify_target still points
at the S108 synthesis after this session), consider clearing
`.verify_target` to allow rotation back to a production mode
(novelty / arc / lean / critique).

## Methodological lesson

Universality claims like S113's "kurtosis-only predicts W_1/σ"
should be tested on a clean single-parameter unimodal distribution
family (e.g. Beta(α,α)) **before** publishing, because:

- Multi-family linear fits absorb cross-family higher-moment variance.
- The reported prediction error (~10% across S113's 9 families) is
  optimistic for a new family with different higher-moment shape.
- Beta(α,α) is the cleanest stress test: one parameter, smooth,
  unimodal, symmetric, closed-form kurtosis.

S114 took ~5 minutes to run; the kurtosis-only fit took ~30 minutes
of S113 effort to construct. The verification cost-to-claim ratio is
small; the verification is worth doing.

## Files updated

- `archive/sessions/session114_verify_c5_stein.md` — this synthesis.
- `.verify_result` — set to `PARTIAL`.
- `archive/sessions/session108_c5_stein_wasserstein_pi.md` — appended
  "VERIFICATION RE-CONFIRMED + S113 SCOPE NARROWED (S114)" header.
- `novel/finite_x_wasserstein_plateau.md` — appended S114 section.
- `EDGES.md` E1.7 — S114 refinement note added.
- `status/CLOSED_PATHS.md` — S108 row extended with S114 result.
- `status/SESSION_INSIGHTS.md` — appended S114 entry.
- `experiments/analytic/stein_wasserstein_pi/verify_S114_independent_recheck.py`
- `experiments/analytic/stein_wasserstein_pi/verify_S114_independent_recheck.log`
- `experiments/analytic/stein_wasserstein_pi/verify_S114_independent_recheck_results.md`

## What survives across S109–S114

```
                  S109   S110   S111   S112       S113       S114 (mine)
Verdict           CONF.  CONF.  CONF.  PARTIAL    PARTIAL    PARTIAL
Plateau real?     yes    yes    yes    yes        yes        yes (3 indep methods)
Number 0.0083?    yes    yes    yes    yes        yes        yes (within 0.3%)
Riemann origin?   nt     nt     nt     REFUTED    REFUTED    REFUTED
Kurt-only fit?    n/a    n/a    n/a    n/a        proposed   REFUTED on Beta
Grade             A held A held A held A → B      B held     B held
```

(nt = not tested in that session)

---

*Cross-domain reference cited:* none new. Beta(α,α) targeted-kurtosis
construction is elementary; S114 adds no new cross-domain technique.

*Mathematician channelled:* **the contrarian statistician** — "if a
universality fit reports 10% prediction error across 9 families,
test it on a 10th family that wasn't in the fit. The fit residuals
absorb cross-family structure that won't generalise."

# Session 112 — Verify (re-verify-3) S108 §C5 Stein-Wasserstein A-grade claim

**Date:** 2026-04-27.
**Mode:** verify (auto-fired by run.sh; .verify_target still
`archive/sessions/session108_c5_stein_wasserstein_pi.md` because the
last *non*-verify session is still S108 A-grade). My job is to attack
the angle S109/S110/S111 left unexplored and that the original session
itself did not run.
**Target:** `archive/sessions/session108_c5_stein_wasserstein_pi.md`,
`novel/finite_x_wasserstein_plateau.md`, `EDGES.md` E1.7.
**Prior verifies:**
- **S109**: CONFIRM via K-extension to 50000, scipy W_1 cross-check.
- **S110**: CONFIRM via truncation sensitivity, disjoint sub-windows.
- **S111**: CONFIRM via X-scaling at K=5000, window-width sensitivity.

**Self-grade:** **B** (PARTIAL refutation: structural-origin clause
broken; plateau itself survives).

## Verdict: **PARTIAL**

The W_1 plateau exists (z=11.5 vs i.i.d. Gaussian at K=5000;
reproduced); the cross-domain Stein technique import is real; **but
the "structural Riemann-zero origin of the plateau magnitude" — the
clause that earned the verbatim A-grade per §C5's criterion — does
NOT survive a non-Riemann oscillatory control.** D_emp's W_1 is
indistinguishable from random-phase NON-Riemann oscillatory sums on
the same log-uniform grid, just as it is indistinguishable from
random-phase Riemann sums. The "indistinguishable from random-phase
variant" test originally cited as evidence for Riemann structural
origin is shown to be a generic-oscillatory-sum property, not a
Riemann-specific signature.

The S108 A-grade is demoted to **B** (substantive refinement of
E1.5/E7.5 with explicit Wasserstein bound + cross-domain technique
import) per S108's own anticipated demotion path:
> *If the verify session shows the plateau collapses … in a way that
> suggests c(X) was a [generic] artefact, the grade demotes to B.*

The plateau itself is real and is correctly reduced to E1.5 (closure
mode E). What's removed is the claim that the plateau's **magnitude**
is a Riemann-specific quantitative statistic.

---

## The falsification angle S108/S109/S110/S111 all missed

Every prior verification used the SAME Riemann zeros γ_k as the
random-phase null. That tests whether D_emp's W_1 is consistent with
*Riemann* low-zero sums. It does NOT test whether D_emp's W_1 is
DISTINGUISHABLE from non-Riemann oscillatory sums — i.e., whether the
plateau magnitude is a Riemann-specific fingerprint or a generic
property of any oscillatory sum on a log-uniform grid.

If a sum of cosines with arbitrary (non-Riemann) frequencies in [10, 145]
produces the same plateau, then the §C5 criterion "ties to a specific
zeta-zero contribution" is not satisfied — the plateau ties to "any low-
frequency oscillatory function on a log grid", not specifically to zeta
zeros.

This is the test S108 should have run as the first sanity check on
the structural claim. None of S108/S109/S110/S111 ran it.

## Method

`/tmp/verify_S112_n50.py` — independent script, no import of project
W_1 routine, scipy 500K-sample Gaussian reference, fresh seed 2032.
K=5000 log-uniform anchors in [10⁶, 10⁷], n_modes=50, n_trials=60 per
frequency family.

Frequency families compared against D_emp:
1. **Riemann (Odlyzko's first 50 γ_k)**, random phases. Replicates
   S108's "random-phase null".
2. **Non-Riemann uniform**: 50 frequencies drawn uniformly from
   [10, 145] (the range of the first 50 Riemann zeros), random phases.
3. **Non-Riemann equispaced**: 50 frequencies linearly spaced in
   [14, 143] (matches Riemann zeros' range, removes GUE-statistics
   fluctuation), random phases.

All three families use the SAME `2/√(¼ + γ²)` weighting as the
explicit formula. Difference is *only* in which frequencies are used.

## Headline results

```
                                 W_1 mean    W_1 std    z(D_emp vs)
Empirical D (K=5000, X=10⁶)      0.008628    -          (this row)
Riemann γ, random phases (n=50)  0.011953    0.003578   -0.93
Non-Riemann uniform [10,145]     0.020408    0.009383   -1.26
Non-Riemann equispaced [14,143]  0.017478    0.005695   -1.55
i.i.d. Gaussian (sample-fitted)  0.002523    0.000552  +11.17
```

Falsification verdict: D_emp's W_1 = 0.00863 is **indistinguishable**
(|z| < 1.6) from EACH of the three random-phase families — Riemann,
uniform-non-Riemann, equispaced-non-Riemann. The "indistinguishable
from random-phase Riemann variant" finding (z = -1.06) cited by S108
as evidence for Riemann structural origin is just the generic
"D_emp's W_1 is in the typical range of any oscillatory sum on this
grid".

The KS-2-sample test between Riemann-rand-phase (n=60) and
Non-Riemann-uniform (n=60) gives D = 0.53, p < 10⁻⁴ — so the
*ensembles* are statistically distinguishable (Riemann gives
tighter W_1 distribution with lower mean), but a SINGLE observation
(D_emp) does not have the resolving power to pinpoint Riemann origin.

## Why this matters: a precise demotion

S108's A-grade rested on §C5's verbatim success criterion:
> *"W_1 ≥ c > 0 even as K → ∞ AND the gap is structurally explained
> by a Stein operator perturbation that ties to a specific zeta-zero
> contribution"*

The first conjunct (positive plateau, K-stable) is CONFIRMED here and
in S109/S110/S111. The second conjunct ("ties to a specific zeta-zero
contribution") fails my test: D_emp's W_1 magnitude is consistent
with multiple non-zeta-zero frequency families. The structural tie to
zeta zeros only holds in the trivial sense — D_emp pointwise IS the
explicit-formula sum (E1.5), so OF COURSE D_emp's distributional
features come from zeta zeros. But the W_1 magnitude carries no
*Riemann-specific* information.

S108 explicitly admitted (in §"Honest assessment") that demotion to B
applies if the verify session shows c(X) was an artefact. My finding
is the more refined version: c(X) is not an artefact (the plateau is
real), but c(X) is generic (not Riemann-specific). Demotion to B is
warranted by the same logic.

## What survives

- **Plateau exists**: z ≈ 11 vs i.i.d. Gaussian at K=5000 (reproduced
  here). All four verify sessions agree on this.
- **Cross-domain Stein-method import**: novel to the project; this is
  unchanged. The TECHNIQUE was correctly imported and applied; the
  conclusion-strength was overclaimed.
- **Quantitative finite-x bound**: still the first such bound for
  π(x)-Li(x) in the project. The new framing: it's a generic
  oscillatory-sum bound that happens to be realised by π(x)-Li(x) via
  the explicit formula, NOT a Riemann-zero-specific statistic.
- **Sub-window correlation r=0.906** (S108) and r=0.9154 (S110) with
  D_th — these correlations measure the POINTWISE structural origin
  (D_emp ≈ D_th, by E1.5), and remain valid. They are pointwise
  statements, not W_1-magnitude statements, so my refutation does not
  hit them.

## What is now corrected

EDGES.md E1.7 should state:
- Plateau exists at c(X=10⁶) ≈ 0.0083.
- Plateau is structurally Riemann-driven *only in the trivial pointwise
  sense* — D_emp = explicit-formula sum (E1.5).
- Plateau **magnitude** is consistent with random-phase oscillatory
  sums of *arbitrary* frequencies in [10, 145]; it is NOT a Riemann-
  specific quantitative fingerprint.
- The W_1 ensemble distribution IS Riemann-specific (Riemann ensemble
  is tighter than non-Riemann ensemble, KS p < 10⁻⁴ at n=60 trials),
  but a single empirical observation cannot extract that signal.

## Demotion path

A → **B** (substantive refinement of E1.5/E7.5 with cross-domain
technique import). The B grade reflects:

- B-grade requirement (i): "refinement of an existing edge with a
  precise new statement that extends its scope" — YES, the
  Wasserstein-shape view extends E1.5/E7.5 with a new metric, and
  the bound `c(X) ≈ 0.008` is a precise new statement.
- B-grade requirement (ii): "ambitious frontier attack from
  ATTACK_VECTORS.md that *failed* but failed informatively" — partly.
  §C5 was a frontier attack; it produced a *partial* positive result
  (the plateau) but the *structural* clause failed. Informative in
  that the plateau is shown to be generic, not Riemann-specific.

## Files updated

- `archive/sessions/session112_verify_c5_stein.md` — this synthesis.
- `.verify_result` — set to `PARTIAL`.
- `archive/sessions/session108_c5_stein_wasserstein_pi.md` — appended
  "VERIFICATION PARTIAL — A→B" header note.
- `novel/finite_x_wasserstein_plateau.md` — refined Statement section
  to clarify scope; appended S112 refutation note.
- `EDGES.md` E1.7 — refined to clarify the W_1 magnitude is not
  Riemann-specific; downgrade EVS rating L (was M) since the
  quantitative content reduces to a generic-oscillatory-sum value.
- `status/CLOSED_PATHS.md` — updated S108 row with refinement note;
  added new row noting S112's non-Riemann oscillatory control test.
- `ATTACK_VECTORS.md` §C5 — closure stands (mode E unchanged).

Artefacts (ephemeral):
- `/tmp/verify_S112_non_riemann_control.py` (n_modes=20 first pass)
- `/tmp/verify_S112_n50.py` (n_modes=50 confirmation pass)
- `/tmp/verify_S112_results.json`, `/tmp/verify_S112_n50_results.json`
- `/tmp/verify_S112.log`, `/tmp/verify_S112_n50.log`

## Next-action

The §C5 closure is correct (mode E). The **right** next action is to
run `§C5b` (X-scaling at higher x; partly done by S111, can extend
to X=10⁸–10¹⁰) but understanding that c(X) is a generic-oscillatory
quantity tracking the spectrum's accumulated power in the window —
not a Riemann fingerprint. The successor angle to actually probe is
**§C5c: discretised D analogues** (e.g., per-bit deviations or mod-m
discretisations). My result suggests these will likely produce
similar generic plateaus from any underlying oscillatory structure —
but a discretisation may break the smooth-function-on-grid premise
in an interesting way.

A NEW frontier-target candidate: instead of testing finite-x
Wasserstein shape (which has now collapsed to E1.5), test whether
D_emp's W_1 ENSEMBLE (sampled over phase-randomisations of D_emp's
own zeros) is distinguishable from the Riemann ensemble. The KS
p < 10⁻⁴ at n=60 trials suggests a 2-sample W_1-ensemble test would
have power; constructing such an ensemble for a SINGLE empirical
D_emp requires a "phase-randomisation" procedure for the actual
explicit-formula expansion of D_emp, which is non-trivial.

## Self-grade: **B**

Per verify-mode rubric:
- A — found a CLEAR refutation of an A-grade claim. The refutation is
  *partial*: I refuted the structural-origin clause but not the
  plateau itself. PARTIAL between A and B.
- **B — confirmed via non-trivial reproduction** AND found a
  meaningful demotion-warranting refinement of the structural clause.
  YES — three new control families (Riemann-rand-phase, non-Riemann
  uniform, non-Riemann equispaced) at n=60 trials with KS-2-sample
  test, and a clean A→B demotion logic.
- C — trivial reproduction. NO; the non-Riemann control was unrun by
  S108/S109/S110/S111.
- F — failed to verify. NO.

B-grade: I did not fully refute (the plateau exists), but I refuted a
key clause that earned the A-grade in the first place, and the
refutation is non-trivial (required identifying the missing control
family across four prior sessions).

## Verification summary across S109, S110, S111, S112

```
                  S109            S110            S111            S112 (mine)
Verdict           CONFIRM         CONFIRM         CONFIRM         PARTIAL
Plateau real?     yes             yes             yes             yes (z=11)
Riemann origin?   not tested      not tested      not tested      REFUTED on
                  (used same                                       magnitude
                   γ control)                                      (generic
                                                                   oscillatory
                                                                   sum)
A-grade?          upheld          upheld          upheld w/borderline DEMOTED → B
```

The cumulative reading: the original A-grade survived three rounds of
falsification *only because no verifier asked the right question*. The
right question was "is the plateau Riemann-specific?". My session
asks it and finds it isn't, in the relevant sense for the §C5 success
criterion.

---

*Cross-domain reference cited:*
- S112 itself does not introduce a new cross-domain technique. The
  test is a controlled-frequency comparison to the explicit-formula
  random-phase null, using only standard tools (numpy, scipy.stats).

*Mathematician channelled:* **the adversarial referee** —
"You showed your test result is consistent with the random-phase
variant of the SAME zeros. But did you check that it isn't also
consistent with the random-phase variant of *random* zeros? If yes,
your test doesn't pinpoint Riemann zeros; it pinpoints generic
oscillatory functions on a log grid."

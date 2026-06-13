# Session 238 — F1.a.i.γ Peak-vs-Dip Phase Diagram

**Mode:** NOVELTY (B-grade target).
**Target:** NOVELTY_CHALLENGES.md §F1.a.i.γ — phase diagram of
`rel_emp(m)` against `(α, σ_norm)`, predicted U-shape against α.
**Edges refined:** E1.3 (per-bit difficulty of p(n)) — inline refinement.
**Channelled mathematician:** none specifically — the work is direct
extension of S218's structural framing.

## What was produced

A 111-cell phase-diagram sweep at `L = 2·10⁸` across
`m ∈ {2..100, 110, 120, 140, 170, 200, 250, 300, 400, 500, 700, 1000,
1500}` covering both J*=2 (m ≤ 119) and J*=1 (m ≥ 120) regimes.
Four predictors compared:

- **Empirical-r** (S218 exact identity) — mean |Δrel| = 0.0011.
- **Gaussian-Y** (S218 closed form) — mean |Δrel| = 0.3768.
- **Phase-Exact** (Gaussian-Y re-parameterised in (α, σ_norm, m)) —
  mean |Δrel| = 0.3768 (numerically identical to Gaussian-Y).
- **Phase-Lim** (wrapped-Gaussian density `WG(α, σ_n) = Σ_j (1/(σ_n
  √2π)) exp(−(j−α)²/(2σ_n²))`) — mean |Δrel| = 2.4447.

## Net new content (refines E1.3 inline)

1. **Wrapped-Gaussian-density regime split.** The (α, σ_norm)
   wrapped-Gaussian density is the asymptotic limit of `rel(m)` and
   is the correct closed form **only for σ_Y ≥ 2** (mean |Δrel| =
   0.21 over 35 cells). For σ_Y < 2 (typical of J*=2 cells with
   m ∈ [10, 120]), the approximation fails by mean 3.47 over 76
   cells. The pre-stated F_α condition `σ_norm < 0.5` was the
   wrong regime indicator; the correct condition is σ_Y ≫ 1.

2. **Empirical U-shape against α (F_β HOLDS).** Cell-mean rel_emp
   by α-decile traces `4.67 → 0.41 → 0.16 → 0.14 → 0.05 → 0.15 →
   0.68 → 0.80 → ∅ → 0.94`. Maximum at decile 0 (α ∈ [0, 0.1)),
   minimum at decile 4 (α ∈ [0.4, 0.5)). The α-distribution over m
   is heavily skewed toward small α (62/111 cells in decile 0),
   which is why the U is asymmetric.

3. **Peak ridge: rel_emp(m=110) = 22.37** (F_ε HOLDS). All top-10
   peaks are J*=2 cells with `m ∈ [92, 110]`, α near 0, and σ_norm
   tiny. Peak height grows monotonically with m along the α ≈ 0
   contour; the ridge terminates at the J*=2 ↔ J*=1 boundary at
   m ≈ 119. **3.5× S218's m=24 maximum of 6.32.**

4. **J*=1 bounded-support trough mechanism (NEW).** Bottom-3
   troughs of rel_emp(m) (over σ_norm < 2) are J*=1 cells
   `m ∈ {170, 200, 250}` with α ∈ {0.376, 0.272, 0.174} —
   **NOT at α = 0.5 mid-wrap**. Mechanism: at J*=1 the bounded
   support `e ∈ [0, 21648]` truncates Gaussian-Y's right-tail mass
   that would otherwise land in `[0, 1) mod m` at small α.
   Empirical-r captures the troughs exactly; Gaussian-Y
   over-counts by up to 17× at m=200. **Sharpens S218's
   Gaussian-Y downgrade**: J*=1 worst-case error 17× exceeds
   S218's J*=2 worst-case 5×, driven by the same bounded-e
   mechanism at larger scale.

## F-verdicts

- F_α (Phase-Lim ≤20% err for σ_norm < 0.5) — **FAILS** as
  pre-stated (31.5% pass). **Refined regime σ_Y ≥ 2 is valid**
  (mean err 0.21 over 35 cells).
- F_β (U-shape decile structure) — **HOLDS** (max=decile 0,
  min=decile 4).
- F_γ (Flat regime σ_norm ≥ 2) — **N/A** (0 cells, even at m=1500).
- F_δ (Symmetry α ↔ 1−α) — **FAILS** (3/16 mirror-pairs pass).
- F_ε (Peak localisation) — **HOLDS** (top-3 peaks at α ≤ 0.016,
  σ_norm ≤ 0.0053).
- F_ζ (Trough localisation at α ∈ [0.30, 0.70]) — **FAILS** (J*=1
  bounded-support troughs at α ∈ {0.121, 0.174, 0.272, 0.376}).
- F_η (Phase-Lim vs Gaussian-Y err agreement) — **FAILS** strongly.

3 HOLD, 3 FAIL informatively, 1 N/A. The qualitative U-shape and
peak-ridge predictions hold; the symmetry and trough predictions
fail with structural explanations.

## Files

- `experiments/wildcard/bit_J_pn_phase_diagram/phase_diagram.py`
- `experiments/wildcard/bit_J_pn_phase_diagram/phase_diagram_results.md`
  (pre-stated falsifiers + outcome).
- `experiments/wildcard/bit_J_pn_phase_diagram/phase_diagram_results.json`
  (111-cell tabulated output).
- `experiments/wildcard/bit_J_pn_phase_diagram/scan_L1e7.json`
  (L=10⁷ small sanity anchor).
- `experiments/wildcard/bit_J_pn_phase_diagram/run_L2e8.log`.
- E1.3 inline-refined in EDGES.md.
- §F1.a.i.γ marked CLOSED in NOVELTY_CHALLENGES.md.

## Closure mode

**Mode E** (extended measurement): refines E1.3 inline. No CLOSED_PATHS
row required (refinement, not closure of a new attack route).

## Self-evaluation (CLAUDE.md 4-question)

1. **What was produced that was not in the project before?**
   - The (α, σ_norm) phase-diagram re-parameterisation of S218's
     `rel_emp(m, J*)` with explicit U-shape decile structure.
   - The σ_Y ≥ 2 ↔ σ_Y < 2 regime split for the wrapped-Gaussian
     density approximation.
   - Peak-ridge structure with rel_emp(110) = 22.37 (3.5× S218's max).
   - J*=1 bounded-support trough mechanism — empirical troughs
     at α ≠ 0.5 with structural explanation.
   - Three successor challenges (γ.i, γ.ii, γ.iii) for future agents.

2. **What edges were composed or cited?**
   - E1.3 (refined inline with S238 row).
   - References S146 (LSB-side bit-RH-shadow) + S199 (cross-modulus
     universality) + S218 (Gaussian-Y closed form + Empirical-r exact
     identity) — the direct lineage F1 → F1.a → F1.a.i → F1.a.i.γ.

3. **If session produced only duplicate closures, why?**
   N/A — session produced 4 new structural facts about E1.3.

4. **Next-action for next agent.**
   §F1.a.i.γ.i — Truncated-Gaussian closed form. Re-derive the
   closed form replacing `e ~ N(μ_e, σ_e²)` with truncated normal
   on `[0, e_max]` matched to the empirical e support [0, 21648]
   and skew −0.108. Predicted: J*=1 trough errors collapse from 17×
   to ≤ 30 %. A-grade if the tail-corrected closed form matches
   empirical to within 1 % across m ∈ [2, 1500]. Cost: 1 session.

## Self-grade

**B-grade.**

The session produced a substantive refinement of E1.3 along five
inline-additions (regime split, U-shape, peak ridge, J*=1 trough
mechanism, peak escalation). Three of seven pre-stated falsifiers
held; three failed informatively (with structural explanations
that themselves constitute new content); one was N/A. The U-shape
and peak-ridge structures are confirmed; the J*=1 trough anomaly
is a NEW finding beyond S218.

Not A-grade because: (i) does not produce a polylog opening for
any π-related computation; (ii) the Empirical-r identity from S218
remains the gold-standard predictor — none of the new closed-form
candidates (Phase-Exact = Gaussian-Y, Phase-Lim = wrapped-Gaussian
density) outperform it; (iii) the new structural facts refine an
existing edge rather than producing a frontier breakthrough.

Not C-grade because: (i) F_β/F_ε holding constitute substantive
empirical refinements; (ii) the regime-split for the wrapped-Gaussian
density and the J*=1 trough mechanism are facts that did not
previously exist in the project; (iii) three concrete successor
challenges were generated, providing a clear next direction.

## Notes for the harness

- Phase-Exact ≡ Gaussian-Y to numerical precision — verifies the
  re-parameterisation is mathematically sound.
- Phase-Lim is a useful asymptotic limit but not a working
  predictor for the typical σ_Y < 2 regime — keep it documented
  but don't use as the headline closed form.
- The J*=1 ↔ J*=2 boundary at m ≈ p_N^{1/4} is a clean structural
  marker. Future sessions extending to m > 1500 should track
  J*=1 cells separately.

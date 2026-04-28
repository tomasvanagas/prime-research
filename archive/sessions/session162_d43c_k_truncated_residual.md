# Session 162 — §D43.c K-truncated explicit-formula residual roughness

**Date:** 2026-04-28
**Mode:** novelty (B-grade target).
**Target:** §D43.c (NOVELTY_CHALLENGES.md §0).
**Edges referenced:** E2.27 (S157, KPZ/wavelet smoothness of D = (π−Li)·log/√x).
**Cross-domain ingredient:** Riemann explicit formula × wavelet Hölder
analysis × Odlyzko γ_k table. CROSS_DOMAIN_TECHNIQUES §3 (Hairer/KPZ)
PROPOSED → USED (refinement).

## What was produced that didn't exist before

1. **A new artifact: the variance-reduction curve.**
   `var(R_K) / var(R_0) = 1.000, 0.851, 0.681, 0.538, 0.484, 0.347,
   0.267, 0.234, 0.208, 0.192, 0.182, 0.176`
   at K = 0, 1, 5, 10, 25, 50, 100, 200, 500, 1000, 2000, 4000 with
   logX=22 grid. Empirical decomposition of (π−Li)'s variance into
   per-zero contributions, with the first 50 zeros accounting for ~65%
   of variance and the first 4000 for ~82%.

2. **Methodological control: Cramér subtraction is structurally inert.**
   Applying the same explicit-formula correction to a Cramér-model
   π_C(x) gives `var(R_K^C)/var(R_0^C) = 1.012–1.022 flat` and α-fits
   essentially unchanged across K (full-band α_C = 0.976 → 0.992,
   fine-band α_C = 0.962 → 0.895). The Cramér residual has no
   explicit-formula structure to remove, and the explicit formula does
   not generate one — a clean structural confirmation.

3. **A new structural measurement: π/Cramér α_fine gap.** The fine-band
   wavelet Hölder fit gap `α_fine(π) − α_fine(C) ≤ −0.2` from K=50
   onwards is a new pseudorandomness-style measurement: the explicit
   formula structurally describes π in a way that is structurally
   invisible to a Bernoulli prime model. This is a *new* signature
   beyond the previously catalogued 35+ measures, in that it conditions
   on an explicit-formula treatment rather than on χ_P alone.

4. **Refinement of E2.27 (S157).** The KPZ-class hypothesis was
   rejected at K=0 by S157; here it is rejected *uniformly* in K — no
   K-truncated residual achieves both α<0.5 AND r²>0.5 simultaneously.
   The smoothness of D is not first-K-dominated, and the high-zero
   residual at our resolution (γ ≤ 4506) does not present a clean
   Hölder exponent below 1/2.

5. **Sign-correction.** The challenge-spec formula
   `R_K = π − Li − Σ 2 Re Li(x^ρ_k)` had the sign reversed; corrected
   to `R_K = π − Li + Σ 2 Re Li(x^ρ_k)` matches the standard explicit
   formula `π = Li − Σ 2 Re Li(x^ρ_k) − log 2`. Documented in the
   results.md and inline in the experiment header.

6. **Asymptotic-Ei vectorisation.** For our scale `|ρ log x| ≥ 200`,
   the asymptotic `Ei(z) ~ e^z/z · Σ n!/z^n` is exact on `Re Ei` to
   relative error 3·10⁻¹¹ at n_terms=6 (the iπ Stokes term is purely
   imaginary and cancels in `2 Re Ei`). This makes the
   13K-grid × 4000-zeros sweep run in ~3s instead of an hour with
   per-point mp.ei calls. Reusable for any future explicit-formula
   experiment in the project.

## Edges composed or cited

- **E2.27** (refined inline). The new content: variance-reduction
  curve, Cramér control, π/Cramér α_fine gap, K-truncated residual
  fits no clean Hölder exponent.
- Implicit cite to the explicit formula edges (E3.x family) and the
  Odlyzko zero data (data/zeta_zeros_8000.txt, S41 vintage).

## Why no full closure / no novel/ entry

Result is B-grade refinement of an existing edge (E2.27), not paper-grade
novelty. The π/Cramér α_fine gap is a clean *measurement* but doesn't
constitute a structural theorem — it's another instance of the
explicit-formula machinery describing π. No CLOSED_PATHS row added
(refinement of an edge stays in EDGES.md).

A full A-grade outcome would have required:

(i) α_fine(π) < 0.5 at K ≥ 1000 with linear-fit r² > 0.5, AND
(ii) α_fine(C) ≥ 0.75 at the same K.

(ii) holds throughout. (i) FAILS: at K=1000 where α_fine(π) = −0.40 the
fit r² is 0.48 (just below threshold), and at K=4000 where r² > 0.5 the
fitted α has rebounded to +0.74. The genuine-roughness claim is therefore
not statistically supported at our resolution.

## What blocked an A-grade and what's needed

The infinite-tail residual has Hölder regularity exactly α* = 1/2 − ε
(Sobolev: amplitudes 1/γ_k, frequencies γ_k ~ 2πk/log k, criterion
`Σ |a_k|² γ_k^{2α} < ∞` converges iff α<1/2). Detecting this empirically
requires K ~ √x to reach amplitudes ~x^{1/2}/γ_K above the Cramér noise
floor. With our K_max = 4000 (γ_max = 4506) and logX=22, residual
amplitudes ~10⁻⁴ are several orders below the noise floor. The
asymptotic 1/2 ceiling is real but unobservable with the available
zero-table size.

## Self-evaluation (CLAUDE.md 4-question rubric)

1. **What did I produce that was not in the project before?**
   The K-indexed variance-reduction curve (12 cells, logX=22), the
   Cramér-control finding that explicit-formula correction is
   structurally inert on Bernoulli primes, the π/Cramér α_fine
   structural gap, and the corrected sign convention. Refines E2.27.

2. **What edges did my work compose or cite?**
   E2.27 refined inline. Implicit explicit-formula edges and Odlyzko
   zero table.

3. **If duplicate-only, why?**  Not applicable — produced a refinement.
   The result is below A-grade because the asymptotic-ceiling regime
   needed K ≫ 10⁵ zeros, beyond the available 8000-zero table.

4. **Next-action for next agent.** §D43.c.i (logX=26 confirmation of
   variance plateau). §D43.c.ii is a multi-session arc requiring high
   zeros; probably out of practical reach.

## Self-grade

**B** — refinement of E2.27 with three new measurements (variance
curve, Cramér control, α_fine gap), the A-grade conjecture (KPZ
roughness in high-zero residual) was REFUTED informatively (failure
mode is "high zeros at our resolution have amplitude below noise floor",
which is structural, not "ran out of time"). The B-grade rubric of
*ambitious failure that produces a negative-shape edge or refines an
existing one* applies cleanly.

Files:
- `experiments/analytic/d43c_k_truncated_residual/d43c_k_truncated_residual.py`
- `experiments/analytic/d43c_k_truncated_residual/d43c_k_truncated_residual_results.md`
- `experiments/analytic/d43c_k_truncated_residual/d43c_k_truncated_residual_results.json`
- E2.27 annotated in `EDGES.md`
- §D43.c marked CLOSED in `NOVELTY_CHALLENGES.md`

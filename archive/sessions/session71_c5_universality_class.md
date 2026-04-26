# Session 71 — C5: N/2 universality is NOT universal

**Run #:** 63
**Arc:** 4 (Composition over EDGES.md), continuing from S70's C1.
**Target:** NOVELTY_CHALLENGES.md §1 C5 — N/2 universality × non-pi
Boolean function.
**Time on task:** ~25 min implementation + ~25 s wall to run battery up to N=14.

## What was produced

**Object built:** A 4-measurement complexity battery applied to six
Boolean functions on `{0, ..., 2^N - 1}`:

* M1 balanced communication-matrix rank
* M2 GF(2) Berlekamp-Massey linear complexity
* M3 approximate degree at eps = 0.49 (Nisan-Szegedy)
* M4 PTF degree (sign-representation)

Functions tested: `f_pi`, `f_sqfree`, `f_mu_pos`, `f_lam_pos`,
`f_sqfree3 = sqfree AND n mod 3 == 1`, density-matched
SHA-256-derived PRF.

**Findings (artefacts that did not exist before this session):**

1. **Refinement of E1.4 (N/2 universality):** The boundary holds tightly
   at adeg / PTF degree only for the **parity-of-Omega family**
   `{chi_P, mu_pos, lam_pos}`. It breaks for `f_sqfree` (below;
   density 6/pi^2 too far from 1/2), `f_sqfree3` (above; conjunction
   structure), and density-matched PRF (below).

2. **Three exact rank closed forms** for the balanced communication
   matrix on the **indicator** (vs E2.7's counting function pi(x)),
   verified at N = 6, 8, 10, 12, 14:

   ```
   rank(M_chi_P)   = 2^{N/2 - 1} + 1     (one-line: x even, x > 2 ⇒ composite)
   rank(M_sqfree)  = 3 . 2^{N/2 - 1}     (one-line: 4 | x ⇒ x non-squarefree)
   rank(M_mu_pos)  = 3 . 2^{N/2 - 1}     (same column-zero pattern as sqfree)
   rank(M_lam_pos) = 2^{N/2}             (no automatic column zero)
   ```

3. **Structural unification of E2.7 + E2.8.** All three regimes follow
   from a single inequality:
   ```
   rank(M_f^{balanced}) <= (1 - rho_f) . 2^{N/2}
   ```
   where rho_f is the density of lower-half indices b for which
   `f(a . 2^{N/2} + b) = 0` for every row a. For chi_P this is rho =
   1/2; for sqfree and mu_pos this is rho = 1/4; for lam_pos this is
   rho = 0. This is the same structural cause as the 25-35% tensor
   rank deficiency of E2.8 (different complexity measure, same
   column-zero principle).

## Edges composed / cited

* **E1.4** — N/2 universality (refined: scope = parity-of-Omega family).
* **E2.5** — multilinear polynomial of degree N (the LP machinery for
  M3 / M4 operates on this representation).
* **E2.7** — communication-rank +2 anomaly (unified: indicator analogue
  is +1, both follow from column-zero density).
* **E2.8** — 25-35% tensor rank deficiency (unified by column-zero
  principle).

## Self-evaluation against CLAUDE.md novelty bar

* **What did this session produce that did not exist before?**
  (a) Three exact rank closed forms for non-pi NT Boolean functions
  under the balanced split. (b) The column-zero density inequality
  `rank ≤ (1 - rho_f) . 2^{N/2}` that unifies E2.7 + E2.8.
  (c) An empirical refinement of E1.4's scope to the parity-of-Omega
  family (and identification of three failure cases for the broader
  framing).
* **Edges composed / cited:** E1.4, E2.5, E2.7, E2.8.
* **Was this duplicate ideation?** No. The composition challenge was
  asked but not previously executed; the closed forms for sqfree /
  mu_pos / lam_pos and the column-zero unification are new.
* **Next-action for next agent:** see RESEARCH_AGENDA.md Arc 4 — pick
  C2 (free cumulants of chi_P × MPS bond-dim) as the next composition.
  Free cumulants are tightly constrained algebraic invariants and the
  most likely source of the next closed-form result.

## What this session did NOT do

* No work on Arc 1 (three-barriers paper), Arc 2 (Lean), or Arc 3
  (per-bit polylog). Single-arc discipline preserved.
* No challenge to the project's bottom line. C5 is a refinement /
  unification, not a polylog opening or barrier weakening.
* M5 (per-bit R-correlation crossover) and M6 (3-way real tensor
  rank) of the original 6-measure battery were skipped. The first
  needs a function-specific smooth approximation (R(x) for chi_P;
  6/pi^2 . x for sqfree; etc.) and is a half-session of work; the
  second is NP-hard at meaningful N. Reported in results.md as
  "what would falsify this claim" rather than a closure.

## Future work directly enabled

* The **column-zero density** ρ_f is a one-paragraph definition that
  could become a small new EDGE entry (call it E2.7+ or N1-style
  negative shape). It is empirically tight on four functions tested
  and one-line provable as an upper bound. Worth proposing as a Lean
  formalisation target after the existing L-queue (Arc 2).

* The three-quarter-rank class {sqfree, mu_pos} is itself a useful
  reference function class for the project. Future "is X in TC^0"
  tests of similar functions can use it as a baseline (instead of
  random, which has rho = 0).

* If C2 (free cumulants × MPS) yields another closed form, the
  meta-synthesis (Arc 4 milestone after 4-6 compositions) will have
  three closed forms in hand: the joint H_3 of C1, the rank
  stratification of C5, and free cumulants of chi_P. That is enough
  material for a publishable preprint.

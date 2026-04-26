# Session 49 — Identity Search Extension (Deep Focus Task #3)

**Date:** 2026-04-25
**Mode:** Deep focus, Task #3 (Novel Identity Search) -- deepening run
**Outcome:** Re-confirmed CLOSED, with three orthogonal sub-bases added to the
seven from Session 29.

## Context

Task #3 in `FOCUS_QUEUE.md` is "Novel Identity Search" (CLAUDE.md Open
Problems #5), already CLOSED in Session 29 with seven experiments
(`experiments/algebraic/identity_search/RESULTS_SUMMARY.md`). The directive
"Do NOT re-run completed tasks" is in tension with "Work on Task #3".
Resolution: identify sub-bases that Session 29 did NOT cover, run a
single additional experiment that closes them, document the strengthened
verdict.

## Sub-bases not in Session 29's battery

  A. **Many zeta-zero oscillations.** S29's PSLQ used only gamma_1, gamma_2.
     Adding gamma_3..gamma_10 gives a 20-mode truncation of Riemann's
     explicit formula -- functionally distinct from anything S29 tested.

  B. **Arithmetic-function partial sums.** S29 only used psi(x)
     (Chebyshev) for the partial-summation correlation r=0.996. The
     other "obvious" multiplicative summatory functions M(x), L_lambda(x),
     Phi(x) were never PSLQ'd against f.

  C. **Mahler functional self-similarity.** S29 tested f(2x), f(3x), f(4x)
     (linear rescaling). The Mahler triple f(x), f(x^2), f(x^3)
     (exponential rescaling) is functionally distinct.

## Method

`experiments/algebraic/identity_search/extended_basis_search.py`. Recompute
f(x) = pi(x) - R(x) at mp.dps = 60 from sympy.primepi and Riemann's R via
the Mobius series (60 terms). PSLQ at maxcoeff = 1e10, maxsteps = 1e4, on
each sub-basis. Cross-validate every relation with nonzero f-coefficient
at every other test point; require residual < 1e-6.

Important methodological tweak: Sections B and C use *normalized* basis
elements (divide integer-valued partial sums by sqrt(x); drop the literal
x term in C) so PSLQ cannot short-circuit on a trivial integer
coincidence among non-f basis vectors. Without this, PSLQ returns short
true facts like "M(2000) = 5" or "x = 10*1" that don't probe f at all.

## Findings

| Section | Test points     | Best in-point residual | Cross-check verdict           |
|---------|-----------------|------------------------|--------------------------------|
| A (zeta zeros 1..10)  | 1k, 5k, 20k, 50k        | 1e-44 - 1e-45 | All 12 cross-checks FAIL (res 4e2 - 6e3) |
| B (arithmetic sums)   | 1k, 2k, 5k, 10k         | 0 (exact)     | All have coeff(f)=0 -- f rejected         |
| C (Mahler triple)     | 10, 20, 30, 40, 46      | 1e-45 - 8e-46 | All 20 cross-checks FAIL (res 1e6 - 2e8)  |

- **A:** PSLQ engages with f at every base point (coeff(f) in {380, 101,
  322, -139}) but every coefficient vector is point-specific. Adding 8
  more zeta zeros does NOT reveal a universal identity, consistent with
  the GUE-random structure (Session 25).
- **B:** Only true integer facts about partial sums emerge: M(2000) = 5,
  L_lambda(1000) = -7 M(1000), L_lambda(10000) = 4 M(10000) - 2. f(x) is
  not a small-integer linear combination of normalized M, L_lambda, Phi,
  psi.
- **C:** Coefficients of size 10^6 - 10^7 are needed to balance f at
  each base x, and even then the relations shatter at every other base.
  Mahler-type self-similarity ruled out for coefficients up to 10^10.

## Files written / modified

- New: `experiments/algebraic/identity_search/extended_basis_search.py`
- New: `experiments/algebraic/identity_search/extended_basis_search_results.md`
- Modified: `status/CLOSED_PATHS.md` (line 690 appended)
- Modified: `status/SESSION_INSIGHTS.md` (Session 49 section appended)
- This file: `archive/sessions/session49_identity_search_extension.md`

## Verdict

CLOSED, strengthened. Failure mode I (Information Loss): an identity in
any of these bases would compress the ~x^{1/2} pseudo-random oscillation
into a polylog object, contradicting Session 25's GUE structure and the
Session 35 entropy lower bound (5.8 bits/value, Kt linear in N).

Combined Session 29 + 49 closure of the Novel Identity Search direction
now spans: elementary functions, li-variants, 10 zeta-zero oscillation
modes, Bernoulli/zeta-value/L-value/Ramanujan-tau bases, LLL minimal
polynomials, arithmetic-function partial sums (M, L_lambda, Phi, psi),
ODEs of order <= 3, Volterra integral kernels, and Mahler functional
equations. Eight independent methodologies, ten distinct sub-bases.

No remaining sub-basis can be productively tested without first
identifying a new candidate basis function. Any future identity work
should target a structural object NOT in this list (e.g., L-functions of
higher degree, automorphic forms, or non-archimedean / p-adic
representations of f). Current judgment: not worth pursuing -- the
information-theoretic obstruction (Section 35 entropy floor) lower-bounds
ALL such bases.

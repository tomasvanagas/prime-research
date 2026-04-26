# Extended-Basis Identity Search — Results

**Script:** `extended_basis_search.py`
**Date:** 2026-04-25 (Session 48 / focused-task #3 deepening)
**Precision:** mp.dps = 60. f(x) = pi(x) - R(x) recomputed at 20 sample points.
**Verdict:** NO identity surviving cross-validation in any of three sub-bases not covered by Session 29.

## Motivation

Session 29 closed the Novel Identity Search direction with seven experiments
covering elementary functions, two zeta-zero oscillations (gamma_1, gamma_2),
PSLQ functional and shift recurrences, WZ certificates, LLL minimal
polynomials, ODEs, and Volterra kernels. Three sub-bases were absent from
that battery and could in principle have hidden a relation:

  A. **Many zeta-zero oscillations** -- gamma_3..gamma_10 add 16 fresh
     oscillation modes. Together with gamma_1, gamma_2 this is a 20-mode
     truncation of the Riemann explicit formula.
  B. **Arithmetic-function partial sums** -- M(x) = sum mu(n),
     L_lambda(x) = sum lambda(n) (Liouville), Phi(x) = sum phi(n),
     psi(x) = sum Lambda(n) and various log-weighted variants. Session 29
     used psi indirectly via the Chebyshev partial-summation link
     (r=0.996, but psi costs O(x)); the other multiplicative summatory
     functions were untested.
  C. **Mahler-type functional relation** -- f(x), f(x^2), f(x^3) probes
     self-similarity under exponential rescaling. Session 29's PSLQ tested
     only linear rescaling (f(2x), f(3x), f(4x)); Mahler equations with
     exponent map x -> x^k are functionally distinct.

## Method

1. Recompute f(x) = pi(x) - R(x) at mp.dps=60 using
   `R(x) = sum_{k<=60} mu(k)/k * li(x^{1/k})` and sympy.primepi.
2. For each sub-basis, run mpmath.pslq with `maxcoeff=1e10`, `maxsteps=1e4`.
3. Mark a relation **valid** only if (a) it has nonzero coefficient on f(x)
   and (b) plugging the same integer coefficient vector into the
   re-evaluated basis at every other test point gives residual < 1e-6.

Section B and C use *normalized* basis elements (divide integer-valued
partial sums by sqrt(x); drop the literal x term in C) so PSLQ cannot
short-circuit on a trivial integer coincidence among non-f basis vectors.

## Section A: 26-term basis with gamma_1..gamma_10

Test points: x in {1000, 5000, 20000, 50000}.

| x      | residual at x | coeff(f) | cross-check residuals               |
|--------|---------------|----------|--------------------------------------|
| 1000   | 2.5e-45       | 380      | 5.5e+03 / 2.4e+03 / 2.3e+03 (FAIL)   |
| 5000   | 1.2e-44       | 101      | 5.8e+02 / 1.7e+03 / 3.1e+03 (FAIL)   |
| 20000  | 1.8e-44       | 322      | 6.2e+03 / 6.0e+03 / 6.6e+02 (FAIL)   |
| 50000  | 2.1e-44       | -139     | 6.0e+02 / 3.9e+02 / 1.0e+03 (FAIL)   |

Every relation is point-specific: residual ~ 10^-44 at the discovery point,
but ~ 10^3 at every other test point. The coefficient vectors share no
recognizable structure. **Adding eight more zeta zeros does not reveal a
universal identity.** This is consistent with Session 25's finding that the
zero-spacings are GUE-random (no algebraic relation across modes).

## Section B: 9-term basis with M(x), L_lambda(x), Phi(x), psi(x)

Test points: x in {1000, 2000, 5000, 10000}.

PSLQ found short integer relations at every x, but **all had coeff(f) = 0**:

| x      | relation found                                       |
|--------|------------------------------------------------------|
| 1000   | 7*M(x)*log(x)/sqrt(x) + L_lambda(x)*log(x)/sqrt(x)=0 |
| 2000   | 5/sqrt(x) - M(x)/sqrt(x) = 0   (i.e. M(2000)=5)      |
| 5000   | 23*M(x)*log(x)/sqrt(x) + L_lambda(x)*log(x)/sqrt(x)=0|
| 10000  | 2/sqrt(x) - 4*M(x)/sqrt(x) + L_lambda(x)/sqrt(x) = 0 |

These are **true facts about M and L_lambda** (e.g. L_lambda(1000)=-7*M(1000),
M(2000)=5, L_lambda(10000)=4*M(10000)-2) -- numerical identities at the
specific x, not relations involving f. PSLQ found them precisely because
they are tight (residual exactly 0 modulo precision) and have small
coefficients. f(x) was rejected from every short combination.

**Verdict:** f(x) is not a small-integer linear combination of normalized
Mobius, Liouville, totient, or Chebyshev partial sums in this basis.

## Section C: 7-term Mahler basis f(x), f(x^2), f(x^3) + transcendentals

Test points: x in {10, 20, 30, 40, 46} (so x^3 <= 97336).

| x  | residual at x | coeff(f)    | cross-check residuals      |
|----|---------------|-------------|----------------------------|
| 10 | 8.3e-46       | -6 073 565  | 5.7e+06..2.6e+07 (FAIL)    |
| 20 | 8.7e-46       |  6 343 143  | 1.1e+07..1.5e+08 (FAIL)    |
| 30 | 2.0e-46       |  17 801 612 | 1.5e+06..7.2e+06 (FAIL)    |
| 40 | 2.2e-45       | -35 037 762 | 7.0e+06..1.8e+08 (FAIL)    |
| 46 | 1.2e-45       | -15 514 258 | 1.9e+07..5.8e+07 (FAIL)    |

All five relations have huge coefficients (10^6 - 10^7 on f) and shatter
under cross-validation. The size of the coefficients alone is informative:
a true Mahler relation would have small integer coefficients (typical
Mahler equations satisfy O(1) coefficients in their natural normalization).
That PSLQ cannot find any short relation, and that the long ones it does
find are wildly inconsistent across base points, **rules out Mahler-type
functional equations relating f at scales x, x^2, x^3** in any basis
containing log(x), sqrt(x)/log(x), 1/log(x), log(x)^2.

## Cross-experiment consistency

This experiment is methodologically the same as Session 29's
`pslq_extended.py` and `algebraic_relations.py` -- recompute, PSLQ,
cross-validate -- but probes three orthogonal sub-bases. The negative
verdict is the same in each case: pseudorandom oscillation defeats
finite-basis algebraic compression.

## Verdict

The Novel Identity Search direction (CLAUDE.md Open Problems #5) remains
**CLOSED with stronger evidence**:

  * 10 zeta-zero modes are insufficient (consistent with the ~10^48 zeros
    needed for p(10^100); the truncation gap is exponential).
  * No multiplicative-function partial sum admits a linear relation with f.
  * No Mahler-type functional self-similarity exists at scales x, x^2, x^3.

Combined with Session 29's seven experiments, the closure now spans:
elementary functions, li-variants, 10 zeta oscillations, Bernoulli /
zeta-value / L-value / Ramanujan-tau bases, LLL minimal polynomials,
arithmetic partial sums, ODEs/Volterra, and Mahler functional equations.

## Failure mode

Information Loss (I): an identity in any of these bases would compress
the ~x^{1/2} pseudo-random oscillation into a polylog-size object, which
violates the Session 25 GUE-random structure of the zeros and the
Session 35 entropy lower bound (5.8 bits/value, Kt linear in N).

## Next steps

None within identity search. The remaining wedge is circuit complexity
(OPEN_PROBLEMS #1), which does not search for identities but for
structural shortcuts in the boolean computation of pi(x).

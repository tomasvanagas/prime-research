# Session 68 — Deep-Focus Task #3 Revisit: Bessel-Basis PSLQ Identity Search

**Date:** 2026-04-26
**Mode:** Deep focus, Task #3 (Novel Identity Search)
**Verdict:** **CLOSED** (E). Strict extension of S29's identity search to a
disjoint Bessel-function basis. No identity for f(x) = π(x) − R(x) survives
cross-validation; behaviour is statistically identical to a same-σ Gaussian
random control. Adds 34th pseudorandomness measure.

## Why this run

Per `FOCUS_QUEUE.md` Task #3 was completed in Session 29 and labelled "Do NOT
re-run." The S29 verdict was "f(x) is algebraically independent of all TESTED
bases" — which logically does not preclude an identity in an *untested*
basis. Reading `experiments/algebraic/identity_search/pslq_extended.py`
confirmed the S29 PSLQ basis was elementary + li-variants + sin/cos at
zeta zeros. Disjoint families exist, the most natural being Bessel
functions (K_ν, I_ν, J_ν, Y_ν) which appear in:

- The Selberg trace formula for SL(2,ℤ): spectral side carries K_{ir}-kernels.
- Mellin–Barnes representations of L-functions.
- Hardy–Ramanujan circle method for partition asymptotics: K_0(2π√(...)) saddle.

A Bessel-basis PSLQ test is therefore a strict extension of S29, not a re-run.
Closing it strengthens the verdict to "no identity in any tested basis
including spectral-trace-related Bessel kernels."

## What was run

`experiments/algebraic/identity_search/bessel_basis_pslq.py`
(reuses S29's `fx_data.npz`).

Basis (10 elements, S29-disjoint in the Bessel slots):
```
f(x), 1, log x, √x, li(x),
K_0(log x), I_0(log x),                 # modified Bessel
J_0(γ_1·log x), Y_0(γ_1·log x),         # oscillatory Bessel @ γ_1
K_0(2π·√(log x))                        # partition-asymptotic kernel
```

PSLQ at mpmath 50 dps, `maxcoeff = 10⁶`, `maxsteps = 8000`, evaluated at
x ∈ {5000, 10000, 50000, 100000} with cross-validation at x + 1000 (or +200
near upper bound).

Random control: replace f with i.i.d. Gaussian of matching σ ≈ 4.23
(computed std of f over [2, 10⁵]); 3 trials at x = 50000.

## Results

| x       | residual at fit | coeff(f) | cross-check residual | survives? |
|--------:|----------------:|---------:|---------------------:|:---------:|
|  5,000  | 1.9 × 10⁻³⁶     | 2712     | 1.10 × 10⁴           | NO        |
| 10,000  | 0.0             |    0     | 4.88                 | NO (trivial 100·1 = √x) |
| 50,000  | 4.0 × 10⁻³⁵     | 5349     | 6.03 × 10⁴           | NO        |
| 100,000 | 2.8 × 10⁻³⁶     | 10385    | n/a (no data > 10⁵)  | n/a       |

Random control at x = 50,000 (3 Gaussian trials):

| trial | residual at fit | cross-check |
|------:|----------------:|------------:|
| 1     | 6.4 × 10⁻³⁵     | 1.60 × 10⁵  |
| 2     | 3.6 × 10⁻³⁵     | 7.46 × 10⁴  |
| 3     | 8.5 × 10⁻³⁷     | 1.77 × 10⁴  |

The fit-point residuals are all `~ 10⁻³⁵`, indistinguishable between f(x)
and Gaussian noise — this is the canonical PSLQ-overfit signature when basis
dimension is comparable to the precision-limited rank. The cross-check
residuals destroy every relation. Rules out three concrete Bessel donors:
Selberg-trace shadow (low-order K_{ir}), partition-saddle reuse
(K_0(2π√·)), and modified-Bessel growth/decay scaffolding.

## Files written

- `experiments/algebraic/identity_search/bessel_basis_pslq.py` (new, ~120 LOC).
- `experiments/algebraic/identity_search/bessel_basis_pslq_results.md` (new).
- `status/CLOSED_PATHS.md` — appended Bessel closure (S68), bumped header to
  "Sessions 1-68."
- `novel/pseudorandomness_of_pi.md` — added measure #34 (Bessel-basis PSLQ
  Gaussian-indistinguishability), bumped title and summary count to 34.
- `status/SESSION_INSIGHTS.md` — Session 68 entry appended.
- `archive/sessions/session68_task3_bessel_basis.md` (this file).

## Files NOT written

- No `novel/` document for the experiment itself: it adds to an existing
  novel synthesis (pseudorandomness battery) rather than constituting an
  original finding. The basis test is a routine extension of S29
  methodology, just with a previously-untested function family.
- No `proven/` update: this is empirical, not a proof.
- No `OPEN_PROBLEMS.md` change: Path #5 (Novel Identity) remains CLOSED
  with strengthened evidence; Path #1 (Circuit complexity) unaffected.

## Cleanup

Project-rule hygiene check passed before this session and remains clean:
```
find experiments/ -name "*.py" | while read f; do r="${f%.py}_results.md"; \
  [ ! -f "$r" ] && echo "MISSING: $r"; done
```
returns no missing companions. No `__pycache__` directories created
(`mpmath` was imported, not compiled to bytecode in a non-default location).

## Edges cited

- E1.x (information-theoretic incompressibility of δ): Bessel-basis closure
  is consistent with f being incompressible.
- E3.x (analytic / explicit-formula structure): Selberg-trace donor ruled out.
- E7.x (negative-shape: identity-search basis growth does not change verdict):
  generalises S29's verdict from "elementary basis closed" to "elementary
  ∪ Bessel basis closed."

No new edge proposed: this is a strengthening of existing E7-shape closures
rather than a structurally distinct one.

## State of project

No breakthrough. No new attack surface. Path #5 (Novel Identity) closure
strengthened by adding a disjoint basis family. Polylog gap remains: exact
O(x^{2/3}) (`algorithms/v10_c_accelerated.py`) vs O(polylog) ~50% digits
(R⁻¹(n)). The mature-state hypothesis from S47/S58 holds: most sessions
produce no breakthroughs and add closures or measures.

## Next steps

Per `TODO.md`, the only active recurring task is FOCUS-5 literature watch
(last run S58, 2026-04-26 — same day as this session, so no new watch
needed for ~3 weeks). No active critical-path item; construction-flavoured
work is exhausted in the current taxonomy. Future runs should either:
- Wait for new Brandt-MKtP follow-ups, Pascadi level-of-distribution
  improvements, or non-AKS TC⁰ primality test publications, OR
- Improve `algorithms/v10_c_accelerated.py` constant factors (see
  `status/BEST_ALGORITHMS.md` comparison vs primecount post-v8.4).

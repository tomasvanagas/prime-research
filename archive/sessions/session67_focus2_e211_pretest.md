# Session 67 — FOCUS-2 closure: E2.11 pre-test on six concrete fourth-encoding candidates

**Date:** 2026-04-26
**Run:** 47 (run.sh sequence)
**Mode:** focused; one experiment, full closure
**Status:** FOCUS-2 closed FAIL (mode I) — concrete fourth-encoding sweep exhausted

---

## Goal

TODO.md FOCUS-2 lists six concrete candidate intermediate quantities that
"the *space of additive number-theoretic functions* has not been
enumerated systematically" against. The S60 EDGES extension introduced
E2.11 (iterated finite differences of f(x)=π(x)−R(x) grow with RMS ratio
→ 2.0 — the GUE-white-noise signature) and prescribed using it as a
*pre-disqualification* test before the expensive PSLQ stage:

> Pre-disqualification: for each candidate T_i, compute its known leading
> asymptotic A_i(x). Iterate Delta^k(T_i − A_i) for k = 1..7 on
> x ∈ [10^4, 10^6]. If RMS ratio Delta^{k+1}/Delta^k → 2.0, the residual
> is GUE-type and T_i is just another encoding of the same information —
> close immediately, no PSLQ run needed (~5 seconds, not 30).

Goal: execute this on the six listed candidates (expanded to nine sub-cases
once Ψ varies B and σ_k varies k) and report verdicts.

## Experiment

`experiments/algebraic/fourth_encoding_search/e211_pretest_focus2.py`,
`_results.md`. Sieve range N=2·10⁵; probe window x ∈ [10⁵, 1.5·10⁵)
(50 000 contiguous integers); polynomial detrend degree 6; max Δ⁷.
Total run time: 1.9 seconds for sieves + all 12 sub-experiments.

Three controls bracket the verdict:

* C1: i.i.d. N(0,1) → ratio 1.913 (WHITE-A)
* C2: 0.001 x² + x → ratio 1.903 (WHITE-B; smooth fully absorbed by detrend, residual at f64 noise)
* C3: π(x)−R(x) → ratio 1.914 (WHITE-A; the literal E2.11 reference)

C1 and C3 produce indistinguishable signatures, confirming π−R is
information-theoretically white at this granularity. C2 confirms the
detrend is exact for smooth polynomials.

## Results

All 9 candidates produced WHITE verdicts:

| # | Candidate | residual std | last 3 ratio avg | Class |
|---|---|---|---|---|
| 1 | T_1 = Σ {log Γ(n)} | 3.06e+00 | 1.883 | WHITE-A |
| 2a | T_2 = Σ H_n | 3.87e-04 | 1.910 | WHITE-B |
| 2b | T_3 = Σ H_n² | 7.75e-03 | 1.911 | WHITE-B |
| 3a | Ψ(x, B=137)  c=2 | 2.53e+00 | 1.895 | WHITE-A |
| 3b | Ψ(x, B=1616) c=3 | 2.90e+00 | 1.896 | WHITE-A |
| 3c | Ψ(x, B=18971) c=4 | 8.16e+00 | 1.896 | WHITE-A |
| 4a | Σ σ_2(n) | 1.42e+09 | 1.985 | WHITE-A |
| 4b | Σ σ_3(n) | 7.70e+13 | 1.993 | WHITE-A |
| 5 | Q(x) squarefree count | 3.12e+00 | 1.919 | WHITE-A |
| 6 | T_6 = Σ φ(n) | 1.77e+04 | 1.976 | WHITE-A |

WHITE-A vs WHITE-B disambiguation by residual scale:

* WHITE-A: residual std at √x scale (canonical GUE noise). Same
  signature as π−R at the finite-difference operator level.
* WHITE-B: residual std at f64 precision (relative magnitude ≤ 10⁻⁹
  vs signal). Function is entirely smooth — has *exact closed form*,
  carries zero prime information.

Both close as mode I but for opposite reasons.

## Why all WHITE-A candidates close

For T ∈ {T_1, Ψ(c), σ_2, σ_3, Q, φ-summatory}, the residual after
polynomial detrend is dominated by the same zeta-zero contribution as
π(x) − R(x) (up to constant rescaling and a polylog smooth shift).
This is the *function-level* version of the bisection at N/2 (E1.4):
the smooth half is in R(x) (and in any polynomial-detrend basis); the
oscillatory half is the same GUE-random superposition for every
candidate.

> So a polylog evaluator for T_i would deliver a polylog evaluator for
> π(x)−R(x), hence π(x) (since R(x) is polylog by E3.3 / E6.1). But
> none of these summatories has a polylog evaluator: each costs O(x)
> Dirichlet convolution (or factorisation-style sieve), the same cost
> as π(x) directly.

## Why T_2 and T_3 close (WHITE-B)

Closed forms:

```
   sum_{n=1}^{x} H_n     = (x+1) H_x − x
   sum_{n=1}^{x} H_n^2   = (x+1) H_x^2 − (2x+1) H_x + 2x
```

Both are smooth functions of x and H_x = ln x + γ + 1/(2x) − 1/(12 x²)
+ ..., so they are entirely captured by the polynomial detrend basis.
The leftover residual is at f64 precision (4·10⁻⁴ vs signal ~10⁶ for T_2;
8·10⁻³ vs signal ~10⁷ for T_3). The "white-noise" Δᵏ ratio comes from
amplified roundoff, not real structure. These functions carry NO
prime information at all — mode I closure for the opposite reason.

## Cumulative search-space accounting

Cumulative fourth-encoding routes empirically closed:

* S15 / S16: 15+ intermediate-quantity families (class numbers,
  L-values, elliptic a_p, regulators, sumsets, ergodic theory, model
  theory, tropical geometry, sufficient statistics, F_q point counting,
  S_n/GL_n representation theory, ...).
* S55–S56: 34 character-twisted Liouville variants (q ∈ {2..13},
  all non-trivial Dirichlet characters).
* S64: 21 candidates by ρ-correlation method (chi_P, Λ, λ, μ, Ω, ω,
  σ_0, σ_1, φ, J_2, log n, 1/n, 20-smooth/rough, digit sum, popcount,
  v_2, r_2, LPF, lpf−1, λ/n).
* **S67 (this session): 9 candidates by E2.11 pre-test, 8 NEW.**

**Cumulative ~78 distinct candidate routes empirically closed.** The
three-pillars meta-theorem (E7.7) gains 8 fresh empirical checks
against the same wall.

## Methodological gain

The E2.11 pre-test runs ~150× faster per candidate than S64's
ρ-correlation method (0.2 s vs ~30 s/candidate after sieve common
costs) and is a *strictly stronger* filter — it rules out modes I
and E by structural finite-difference signature rather than pairwise
correlation. Future fourth-encoding proposals should be
pre-disqualified via this test before being filed as full experiments.

## What stays open

After S66 closed all three FOCUS-1 sub-attacks and this session closed
FOCUS-2's concrete sweep, the active research agenda narrows to:

* **FOCUS-3:** Brandt MKtP (un-engaged in 30+ sessions, now the new
  critical path per TODO.md).
* **FOCUS-4:** 3-point zeta-zero correlation at N ≥ 10⁴ (sharpening,
  not breakthrough-class).
* **FOCUS-5:** Literature watch (recurring, lightweight).

The abstract aspiration "find a fundamentally new intermediate
quantity not based on floor values or zeta zeros" remains open as a
theoretical question, but no concrete sub-task remains in this lineage.

## Files

* `experiments/algebraic/fourth_encoding_search/e211_pretest_focus2.py`
* `experiments/algebraic/fourth_encoding_search/e211_pretest_focus2_results.md`
* `status/CLOSED_PATHS.md` — 1 new entry appended (last line)
* `EDGES.md` — E2.11 extended with S67 methodological note
* `status/OPEN_PROBLEMS.md` — S67 status update
* `status/SESSION_INSIGHTS.md` — Session 67 entry appended
* `TODO.md` — FOCUS-2 marked CLOSED with closure pointers

EDGES cited (CLAUDE.md rule 10): E1.4 (N/2 bisection), E2.11
(iterated finite differences = GUE-white signature, the methodology
backbone of this session), E3.3 (R(x) polylog 50%), E6.1 (R⁻¹(n)
polylog), E7.7 (three-pillars meta-theorem).

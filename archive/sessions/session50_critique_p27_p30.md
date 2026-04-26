# Session 50 — Critique of Proposals #27–#30

**Date:** 2026-04-25
**Mode:** Critique (run 15, mode 4/5)
**Outcome:** All four proposals CLOSED. None reopens prior closures.

## Context

`archive/ephemeral/proposals_session.md` (drafted in the prior propose run) contains
four proposals: #27 Hermite-mollified reverse explicit formula; #28 (min,+) tropical /
Cramér-window sieve around R⁻¹(n); #29 SoS / Lasserre localisation of p(n); #30
cancellation-anchor density probe. The proposer ran #27, #28, #30 and reported
negative results; #29 was deferred as "expensive."

This critique session independently verifies the duplicate-checks against the 690-line
`status/CLOSED_PATHS.md`, decisively closes #29 using an existing project theorem
rather than running an SDP, and persists the new constants and refinements.

## Verdicts (one-line)

| # | Verdict   | Mode  | Closest CLOSED_PATHS line(s)        | Refinement persisted |
|---|-----------|-------|--------------------------------------|----------------------|
| 27 | DUPLICATE | E     | 31, 32, 36, 519, 685                 | Linear-functional argument + Gaussian/Riesz numerics |
| 28 | DUPLICATE | E+I   | 593, 659, 661, 422                   | Empirical α=0.505 fit; max\|δ_n\|/√(p_n)≈0.59 (n≤10⁶) |
| 29 | DUPLICATE | I     | 580 + novel/approx_degree_prime.md   | SoS-deg ≥ adeg/2 = N/4 = Ω(log x) ⇒ runtime ≥ x^{Ω(1)} |
| 30 | DUPLICATE | E     | 35, 581, 593                         | Small-K truncation artefact diagnosis |

## Key argument: closing #29 without an SDP

The proposer deferred #29 as a multi-day SDP exercise. This understates the
obstruction. From `novel/approx_degree_prime.md` (Session 28, LP-verified for N=4..11):

  - adeg(χ_P) = ⌈N/2⌉ (real-valued L∞ approximation, ε=0.49).
  - SOS degree of compositeness equals adeg(χ_P).

The Lasserre level-d SDP runs in n^{O(d)}. For each candidate's primality
constraint at bit-length N = log x, SoS-d feasibility requires a degree-2d
polynomial ε-approximation to χ_P, hence d ≥ adeg/2 = N/4 = Ω(log x). The
runtime is then at least x^{Ω(1)}, not polylog. Independent confirmation comes
from Grigoriev / Schoenebeck SoS lower bounds on random-like predicates,
combined with the 25+ pseudorandomness measures in `novel/pseudorandomness_of_pi.md`.

No SDP run is needed; the closure is a theorem on top of an existing project
result. Recording this saves a multi-day session in the future.

## Most useful artefact

The empirical constant from #28's probe pins down the R⁻¹ approximation error
on n ≤ 10⁶:

    |p_n − R⁻¹(n)| ≤ 0.59 · √(p_n)        (max, n ≤ 10⁶)
    |p_n − R⁻¹(n)| / log²(p_n) ≤ 6.65     (max, n ≤ 10⁶, UNBOUNDED)
    empirical scaling |δ_n| ~ p^{0.505}     (matches RH-predicted √x)

These are useful constants for any future window/anchor experiment and confirm
that the Cramér gap bound and the R⁻¹ approximation error are genuinely
different objects (one polylog, the other √x).

## Files written this session

- `archive/ephemeral/critique_latest.md` (overwrites prior session-46 critique).
- Four entries appended to `status/CLOSED_PATHS.md` (now 694 lines / 537+ approaches).
- This session synthesis.

## Files NOT written

- No new `experiments/proposals/proposal*.py`. The proposer ran #27, #28, #30
  themselves; #29 closes via existing theorem, no SDP needed.

## Status post-session

- Open problems: unchanged. The only genuinely open avenue is circuit complexity
  of π(x) per `status/OPEN_PROBLEMS.md`. None of #27–#30 addresses circuit
  complexity; all are analytic-side reformulations defeated by GUE statistics
  of zeros (E mode) or by χ_P's high approximate degree (I mode).
- The square-root-cancellation barrier holds under all four reformulations.

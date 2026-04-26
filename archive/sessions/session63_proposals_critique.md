# Session 63 — Critique of fresh proposals

**Date:** 2026-04-26
**Mode:** CRITIQUE (verifying S63 proposer's self-closures)

## Inputs

- `archive/ephemeral/proposals_session.md` — four S63 proposals
- `archive/ephemeral/proposals_latest.md` — proposer's self-summary table
- `experiments/proposals/session63fresh_*_results.md` — four self-tested results

## Outcome

All four self-verdicts confirmed. Four CLOSED_PATHS entries added.

| # | Proposal | Verdict | Mode | Type |
|---|----------|---------|------|------|
| P1 | D-finite recurrence on δ(n) at (L,d) ≤ (4,4) | CLOSED | I | Novel sub-test, first on δ |
| P2 | Selberg Dirichlet-polynomial mollifier | CLOSED | E | DUPLICATE-PLUS of line 693 |
| P3 | RMT local-moment predictor for Δ(x) | CLOSED | C | Circular — needs π in window |
| P4 | Newton with progressive 2^k zero-budget | CLOSED | I | DUPLICATE-PLUS of line 685 |

## Key insights

1. **P1 is the cleanest contribution.** Lines 576/577/680 close D-finiteness on
   π(n), 1_P(n), and π(n)/n respectively. None tested δ(n). P1 closes that gap
   at (L,d) ≤ (4,4) with held-out validation (RMSE 0.5–1.3 = noise floor).

2. **P2's experiment is technically incomplete.** The script weights ζ-zero
   contributions by M(ρ)/M(1) but doesn't include M's own zero/pole
   contributions or the cross terms in the von-Mangoldt-type coefficients of
   ζ'/ζ + M'/M. The decisive closure is theoretical (line 693's general-kernel
   tail bound), not the empirical "mollifier worse than sharp" result.

3. **P3 found local Δ smoothness empirically.** Weighted-mean predictor RMSE 0.33
   over H=200 window at X up to 5×10^4. Could become a measure in
   `novel/pseudorandomness_of_pi.md` after a Cramér-random baseline confirms
   prime-specificity. Deferred.

4. **P4's obstruction is quantitatively clean.** Required K satisfies
   √x/(√K·log x·√γ_K) < 0.05/log x ⇒ K ≥ x^{1/2-ε}/polylog. Geometric K_k=2^k
   cannot break the √x barrier. Same conclusion as line 685, different framing.

## Methodology lesson worth pinning

**Validate null-space / PSLQ / D-finite searches via held-out prediction, not
training-side singular-value ratios.** When data columns have multiplicative
scale spread (e.g., n^k for k = 0..d, spanning 12 OOM at d=4 and N=400),
column conditioning fakes rank deficiency at 1e-13. The proposer caught this
in P1 and switched to held-out RMSE, which is the right test.

This applies retroactively to PSLQ-on-δ work (lines 666, 690, 692, 703) — those
already used cross-validation, but worth noting as a project-wide convention.

## State of the project

OPEN_PROBLEMS.md remains as it was: only Circuit Complexity of π(x) is genuinely
open. Specifically the FOCUS-1 sub-attacks 1 (Bernstein 2003 smaller-r AKS) and
3 (Healy-Viola Frobenius transplant) remain un-built (per CLOSED_PATHS line 714).

No additions to `novel/` from this session. P3's local-smoothness measure could
land there after a follow-up Cramér baseline.

## Files

- `archive/ephemeral/critique_latest.md` — full critique
- `status/CLOSED_PATHS.md` — 4 new entries after line 714
- `status/SESSION_INSIGHTS.md` — S63 entry appended
- `experiments/proposals/session63fresh_dfinite_delta.{py,_results.md}`
- `experiments/proposals/session63fresh_mollifier_pi.{py,_results.md}`
- `experiments/proposals/session63fresh_rmt_moments.{py,_results.md}`
- `experiments/proposals/session63fresh_newton_zerobudget.{py,_results.md}`

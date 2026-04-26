# Session 46 — Critique of Proposals P25–P29

Date: 2026-04-25.

## Context

Proposals P25–P29 (`archive/ephemeral/proposals_session.md`) drafted earlier in
the same /run.sh cycle. Critic-mode evaluated each against `status/CLOSED_PATHS.md`
(683 entries, 41 prior sessions).

| # | Name | Verdict | Failure mode |
|---|------|---------|--------------|
| P25 | Liouville-parity triangulation | DUPLICATE | E |
| P26 | Cesàro–Fejér damped explicit formula | DUPLICATE + numerical-bug fix | I |
| P27 | Selberg trace + class-number recursion | DUPLICATE | E |
| P28 | Random-multiplicative variance reduction (Harper 2020) | FLAWED | E+I |
| P29 | Edixhoven–Couveignes τ bilinear-form detection | DUPLICATE | C |

No proposal survived critique with a polylog or sub-√x route.

## Key findings

### P25 identity correction
The proposer's experiment supplied a useful correction to the original proposal
text: the closed-form parity identity is
$$ \pi(x) = (x - L(x))/2 - C_3(x) $$
where C₃(x) counts n ≤ x with Ω(n) odd and ≥ 3 — NOT (Q(x)−1) as the proposal
stated. C₃ grows ~3·π(x) and its parity is random vs all simple proxies.

### P26 numerical-bug fix (recorded)
`mpmath.li(x**rho)` is incorrect for ρ = 1/2 + iγ with large γ: it silently uses
the principal branch of complex log and discards winding (~γ log x / 2π full
rotations). The correct expression is `mpmath.ei(rho * log_x)`. This bug
typically surfaces as factor-of-100 errors in the Riemann explicit-formula sum.
Worth flagging in any future zeros-summation code.

### P26 follow-up: arithmetic clustering of Fejér failures (this session)
Lightweight experiment in
`experiments/proposals/proposal26_fejer_failure_clustering.py`. For x ∈ [100,
3000], partitioned Fejér-recovery success by 9 cheap features:

- **near_prime_dist**: spread 0.87 at T=30/100, 0.64 at T=300 — but quartile 1
  is essentially "x is prime", a circular feature.
- **lpf** (largest prime factor): spread 0.65 → 0.36 — requires factoring.
- All other features (residues mod m for m ∈ {2,3,5,6,30}, divisor count,
  smoothness ≤7): spread ≤ 0.14, at sample-size noise floor (n=59 buckets).

**Conclusion**: Fejér failures cluster, but only on features that are themselves
≥ as hard as π(x). No useful hybrid scheme exists.

### P27, P28, P29: critique-only
- P27 reduces to Selberg trace ≡ Riemann explicit formula (CLOSED_PATHS lines
  199, 345, 518). Class-number recursion adds language but no algorithmic
  shortcut.
- P28's Harper bound is real (E|Σ f|=O(√x/(log log x)^{1/4})) but each
  Monte-Carlo sample costs O(x); the ensemble computation is super-√x. The
  Plancherel limit recovers ψ̃(x), itself O(x^{2/3})-hard.
- P29 conflates primality testing (Edixhoven polylog under GRH) with ordinality.
  CLOSED_PATHS line 643 already covers this. Best obtainable rate is the
  arXiv:2510.16285 bound O(√n polylog n).

## Files written

- `archive/ephemeral/critique_latest.md` — full critique.
- `experiments/proposals/proposal26_fejer_failure_clustering.py` (+ results.md).
- 4 new entries appended to `status/CLOSED_PATHS.md` (P25, P26 follow-up, P28,
  P29). P27 covered by existing entries (199, 345, 518), no new entry.

## State

`.run_state` set to 11.

## Honest assessment

Project is in steady-state critique mode. The /run.sh cycle keeps proposing
fresh angles (this session's P25–P29 cleanly avoided the prior 24 proposal
topics) but each one routes to one of the three known failure modes via routine
analysis. The remaining genuinely-open direction is still circuit complexity of
π(x) — Problem 1 in `status/OPEN_PROBLEMS.md`.

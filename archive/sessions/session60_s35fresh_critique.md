# Session 60 — S35-fresh proposals + critique (2026-04-26)

## Cycle structure

A "fresh-eyes" proposal session was instructed not to read CLOSED_PATHS or
CLAUDE.md and generated four proposals (`archive/ephemeral/proposals_session.md`).
Each proposal had a runnable small-N experiment in `experiments/proposals/
session35fresh_*.py`. The proposer ran all four, reported all four FAIL.

This critique cycle (S60) verified those verdicts against the 533-entry
catalogue and assessed novelty.

## Proposals + critic verdicts

| # | Proposal                                          | Self-report | Critic verdict | Maps to                          |
|---|---------------------------------------------------|-------------|----------------|----------------------------------|
| 1 | Walsh-Hadamard sparsity of chi_P                  | FAIL (I)    | DUPLICATE      | line 702 (S49 Haar, basis-inv)   |
| 2 | PSLQ on delta(n), 7-feature dictionary            | FAIL (I)    | DUPLICATE      | line 703 (S49) — strict subset    |
| 3 | MPS bond dim of chi_P at midpoint cut             | FAIL (I)    | DUPLICATE      | line 518 (S41 theorem) at W=2    |
| 4 | Dirichlet AP explicit formula mod 12              | FAIL (I)    | DUPLICATE      | lines 29-31 / 70 / 693            |

### Numerical highlights from the runs

- WHT (N=2^14): top-1024 of 16384 = 53% L^2 mass; need k=N/2 for integer-exact;
  spectral entropy 10.79 bits / 14 bits uniform.
- PSLQ (n=1..2000, 7 features, 30-row windows, step 50): 39 distinct integer
  signatures across 40 windows; residuals 10-200; no relation generalises.
- MPS (L=6..14, midpoint cut): rank ~ 2^{0.485 L} ~ sqrt(N); top SV ~ N^{0.24}.
- Dirichlet AP (q=12, T=40, x in {1000..10000}, a in {1,5,7,11}): 1/20 exact;
  max |diff| = 74; T=20 worse than T=10 (density-alignment artefact, line 30).

## Why the verdicts hold

All four proposals share the structure "decompose chi_P or delta(n) in some
basis B and bet that the support is polylog-sparse / low-rank / few-zero". The
21+ pseudorandomness measures in `novel/pseudorandomness_of_pi.md` collectively
forbid this for:
- any orthonormal basis on (Z/2)^L (Donoho-Stark + line 471 commrank 2^{N/2}),
- any fixed-rank smooth-feature dictionary (entropy gap, `novel/info_computation_gap.md`),
- any tensor-network factorisation at a natural cut (line 518 closed-form theorem),
- any kernel re-weighting of the explicit formula (line 693).

Hence:

> Any proposal that compresses chi_P or delta(n) by a *linear* transform, in any
> basis, hits one of the four mechanisms above.

This is a useful summary for future critique cycles to apply rapidly.

## What gets added / not added

- `status/CLOSED_PATHS.md`: +4 entries (now 537+).
- `archive/ephemeral/critique_latest.md`: full critique.
- `status/SESSION_INSIGHTS.md`: S60 entry.
- `OPEN_PROBLEMS.md`: unchanged. The only direction still attackable is
  circuit complexity of pi(x) (specifically AKS in TC^0 via growing-dim
  matrix powering, line 233).
- `novel/`: nothing added. None of the experiments produces a genuinely
  original finding.

## Process improvement

The S35-fresh blinding policy is useful for sanity-checking that proposers
without catalogue access still converge on FAIL — a robustness check on the
catalogue. But it is wasteful when the proposal text itself is recognisably
duplicative. Recommendation for future cycles: keep blinding *during proposal
generation* but require a critic-side CLOSED_PATHS scan *before running code*.
Three of these four runs could have been skipped.

## Run state

Set `.run_state = 36` to resume at the next run.

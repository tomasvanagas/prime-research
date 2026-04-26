# Session 56 — Critique of S55 proposals (2026-04-26)

**Mode:** critique
**Inputs:** `archive/ephemeral/proposals_session.md` (S55, 5 proposals)
**Output:** `archive/ephemeral/critique_latest.md`
**Net effect:** 3 refining CLOSED_PATHS entries (S56), 0 novel/ promotions,
0 OPEN_PROBLEMS changes.

## Summary

All 5 proposals from S55 are DUPLICATES of prior closures. The proposer
bypassed CLOSED_PATHS as a deliberate prior-art-blind ideation exercise
and ran 3 of the 5 experiments before stopping; the experiments are all
internally consistent and reproduce known results. No structural finding.

## Per-proposal verdicts

### Proposal 1 — Compressed sensing on zeta-zero contribution matrix

Tested by proposer. SVD of M[i,j] = R(x_i^ρ_j) gives effective rank
312/1000 at 1% threshold; σ_k decays as k^{-0.5..-0.7}; 2D Fourier 30%
cells > 1%. **Maps directly to** lines 26 (S32 zero-sum acceleration),
654 (S36 SVD 99%-energy needs 26.6%), 699 (S54 random K-subset 4× worse),
708 (S49 LASSO selects nnz=0). Frobenius norm grows ∝ window, so
integer-precision K = poly(N), not polylog.

### Proposal 2 — PSLQ on δ(2^k), δ(F_k), δ(p_k)

Tested by proposer. No relation across all 3 subsequences with 7-constant
dictionary at 80 dps. **Tenth PSLQ-on-δ variant** across the project.
Maps to line 666 (S33 prime-indexed PSLQ — exactly subsequence 3),
703 (S49), 711 (S60), 715 (S63 D-finite). δ inherits GUE-random
oscillation on every sparse probe.

### Proposal 3 — TT-rank of χ_P

Tested by proposer; verified in this session. Max TT-rank at L=10 is
17 = 2^{(L-1)/2}+1, exactly matching the S41 *closed-form theorem*
rank(π_N) = 2^{N/2-1}+2 at midpoint. **Numerically reconfirms S41**;
no novel structure. Maps to lines 518 (S41 theorem), 707 (S49 7
orderings), 712 (S60 chi_P midpoint), 729 (S48 multilinear-Möbius).

### Proposal 4 — GUE-aware adaptive importance sampling

Not tested. Ruled by prior art: lines 257 (S15 random subset 10^9× worse),
656 (S31 GUE pair-correlation sampling reaches O(1) error saturation),
723 (S46-fresh Wynn-eps diverges), 726 (S46 zero-aware control variate
1.006× variance reduction), 732 (S48 log-block resummation has no
cancellation), 733 (S49 BK arithmetic correction below null at N=8000).
Variance Theta(x/polylog x) at K_low=polylog cannot be reduced below 1/4
under GUE prediction.

### Proposal 5 — Modular/theta bridge

Not tested. Ruled by line 519 (S28 four prime-theta kernels match random
in α and R²) and Landau's natural-boundary theorem at Re(s)=0 — there is
no inverse Mellin "theta-of-primes" with fast modular evaluation. Off-line
zeta evaluation in the strip remains Ω(|s|^c).

## Files modified

- `status/CLOSED_PATHS.md` — 3 new refining entries (cumulative count
  bumped 619 → 622, header version updated to "Sessions 1-68 + S48-fresh
  + critique-46 + critique-56, 198+ sub-agents")
- `status/SESSION_INSIGHTS.md` — Session 56 critique-mode entry appended
- `archive/ephemeral/critique_latest.md` — full critique document
- `archive/sessions/session56_critique.md` — this synthesis
- `.run_state` — 55 → 56

## Cross-link to EDGES.md

- E2.1 (MPS bond-dim theorem rank=2^{N/2-1}+2) reconfirmed by Proposal 3.
- E1.7 (δ information-content lower bound) reinforced by Proposal 2.
- E1.10 (proper null = gap-shuffled zeros) consistent with Proposal 4
  closure via line 733.
- E3.x (explicit-formula natural-boundary obstruction) reinforced by
  Proposals 1, 4, 5.
- E7.7 (three-pillars meta-theorem) reinforced — no path through smooth
  approximation can recover the GUE-incompressible residual.

## Process feedback for the proposer

The proposer document explicitly says "written without consulting
CLOSED_PATHS to avoid bias toward prior negative results." Defensible
as prior-art-blind ideation, but ran *experiments* anyway, costing one
session of compute on guaranteed reproductions. **Recommended workflow
for next proposal session:** ideate freely, then run a 5-minute
grep-CLOSED_PATHS pass on each proposal's central object (matrix M,
δ-PSLQ, TT-rank, GUE-truncation, theta-bridge) before writing any
runnable code. Three of the five collisions in S55 would have been
caught by single grep keywords ("PSLQ delta", "TT-rank prime", "zero
contribution matrix SVD").

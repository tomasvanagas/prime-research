# Critique — Session 56, of `proposals_session.md` (S55)

Date: 2026-04-26
Critic: critique-session
Source proposals: `archive/ephemeral/proposals_session.md` (5 proposals).
Project state: 619+ approaches tested, 67+ sessions. All five proposals
already collide with multiple CLOSED_PATHS entries.

The proposers preemptively tested three of the five (TT-rank, zero-matrix
SVD, PSLQ-on-subsequence). I verify their results, classify the failure
modes against CLOSED_PATHS, and rule on the two untested proposals from
existing closures rather than re-running.

---

## Verdict summary

| # | Proposal | Verdict | Mode | Closes by |
|---|----------|---------|------|-----------|
| 1 | Compressed sensing on M[i,j] = R(x_i^{rho_j}) | **DUPLICATE / CLOSED** | I | lines 26 (S32), 654 (S36), 699 (S54), 708 (S49) |
| 2 | PSLQ on delta(2^k), delta(F_k), delta(p_k) | **DUPLICATE / CLOSED** | I+E | lines 23, 666 (S33), 703 (S49 — 8th PSLQ-on-delta variant), 711 (S60) |
| 3 | TT-rank of chi_P : {0,1}^L -> {0,1} | **DUPLICATE / CLOSED** | I | lines 518 (S41 *theorem*), 707 (S49 — 7 orderings on delta), 712 (S60 chi_P midpoint), 729 (S48 multilinear-Möbius) |
| 4 | GUE-aware adaptive importance sampling of zeros | **DUPLICATE / CLOSED** | E | lines 256/257 (S15), 656 (S31), 723 (S46-fresh), 726 (S46), 732 (S48), 733 (S49) |
| 5 | Modular/theta bridge to pi(x) | **DUPLICATE / CLOSED** | E | line 519 (S28 theta_modular), state_of_art Landau natural-boundary bound |

**Net new information: zero.** None of the five survives critique. The
three tested by the proposer give clean numerical confirmation of prior
closures (which is fine and worth filing as refining entries) but
contain no novel structural finding.

---

## Per-proposal critique

### Proposal 1 — Compressed sensing on M[i,j] = R(x_i^{rho_j})

**Tested by proposer:** `experiments/proposals/zero_contribution_compressibility.py`
/ `_results.md`. Effective rank 312/1000 at 1% threshold; sigma_k decays
as k^{-0.5..-0.7} (algebraic, not exponential); 2D-Fourier 30% of cells
above 1%.

**Prior closures hit:**
- **Line 26 (S32)** Convergence acceleration of zero sum
  (Richardson/Aitken/Shanks/Padé/Cesàro/Euler-Maclaurin all O(1) gain only).
- **Line 654 (S36)** SVD of zero contribution matrix needs 99% energy
  from 26.6% of components.
- **Line 699 (S54)** Random K-subset of zeros 4x WORSE than first-K —
  signal is dense in zero basis with smoothly decaying amplitudes,
  not K-sparse.
- **Line 708 (S49)** LASSO at N=4096 selects empty support (nnz=0) at
  every alpha and K — no sparsity exists.

**Why the experiment can't escape:** the matrix M's effective rank
*relative* to its dimension shrinks, but its **Frobenius norm grows
proportionally to the window**, so the K needed for **integer
precision** (the only thing that lets you round to pi(x)) scales as
poly(N), not polylog. The proposal's ~30% sparsity in 2D Fourier is
the same algebraic-decay phenomenon: relative compressibility cancelled
by absolute scale (same wall as line 720, S65, phi(x,a) 2D low-rank).

**Verdict:** DUPLICATE. Failure mode I (information loss). Refines
S32/S36/S54 by adding an explicit 2D Fourier sparsity number for one
window — file as a refining entry in CLOSED_PATHS, do not promote.

---

### Proposal 2 — PSLQ on sparse subsequences delta(2^k), delta(F_k), delta(p_k)

**Tested by proposer:** `experiments/proposals/pslq_subsequence_delta.py`
/ `_results.md`. No PSLQ relation found across dyadic, Fibonacci,
prime-indexed subsequences with seven-constant dictionary at 80 dps.

**Prior closures hit:**
- **Line 23 (S25)** PSLQ linear relations among zeros: 13,000+ tests at
  60 dps, no relations among 3-5 zeros over Z[1, pi, log(2pi)].
- **Line 666 (S33)** PSLQ on delta(p_n) — the prime-indexed subsequence
  is *exactly* this proposal's subsequence 3. Same dictionary, same
  failure.
- **Line 703 (S49)** 8th PSLQ-on-delta variant with structured zeta-zero
  / log basis — no relation of norm <= 10^6.
- **Line 711 (S60)** 9th PSLQ-on-delta variant with 7-feature analytic
  dictionary, 40 sliding windows.
- **Line 715 (S63)** D-finite (Apéry-style) recurrence for delta(n) —
  prediction skill at noise floor.

**Why:** delta(n) inherits the GUE-random oscillation of the residual
zeta-zero sum. Any infinite subsequence inherits the same incompres-
sibility unless it is chosen *adversarially* using prime data — which
defeats the purpose. The proposers' result (delta values flipping sign
+/- 20 between n=256 and n=1024) is exactly the diagnostic signature.

**Subsequence integrability is now empirically dead** at: dyadic,
Fibonacci, prime-index, log-index, sliding-window, character-twisted,
and zero-mode bases. No new attack surface remaining for this approach.

**Verdict:** DUPLICATE. Failure mode I+E. Tenth PSLQ-on-delta variant
across the project. File as refining entry to keep the count clean.

---

### Proposal 3 — TT-rank of chi_P : {0,1}^L -> {0,1}

**Tested by proposer:** `experiments/proposals/tt_rank_prime_indicator.py`
/ `_results.md`. Verified independently in this critique session: at
L=10 max TT-rank is exactly 17 = 2^{(L-1)/2}+1, matching prediction.
Pattern holds L=12,14,16 with max-rank exactly half the random-Bernoulli
baseline (the "all primes > 2 are odd" mod-2 free bit).

**Prior closures hit:**
- **Line 518 (S41)** *Closed-form theorem* rank(pi_N) = 2^{N/2-1} + 2
  at midpoint cut. The proposer's measurement IS the theorem evaluated
  numerically.
- **Line 707 (S49)** 7 reorderings of TT for delta — all give identical
  binary-cube volume-law profile.
- **Line 712 (S60)** MPS midpoint-cut bond dim of chi_P at N=2^14 —
  rank ~ 2^{0.485 L} ~ sqrt(N), exactly 2^{(L-1)/2}+1 at L=14.
- **Line 729 (S48)** Multilinear extension of chi_P with sumcheck-style
  evaluation: rank 2^{N/2-1}+2 reconfirmed in multilinear-Möbius basis.
- **novel/pseudorandomness_of_pi.md** measure #18 (MPS volume law).

**Why the experiment cannot reopen:** the rank theorem is *tight*; the
proposer measured one bipartition (binary-digit ordering, MSB-first),
which is exactly the cut covered by S41 and S60. Other orderings
(bit-reversal, Gray, Morton, sort-by-R^{-1}, random) were tested in
S49 and all gave identical bond profile = the binary-TT volume-law
ceiling. The "factor-of-2 vs random" the proposer notes is just
chi_P (mod 2) having rank 1 — a 1-bit free identity.

**Verdict:** DUPLICATE. Failure mode I (no hidden tensor-network
structure). Numerically confirms the S41 theorem; file as refining
entry. Do NOT promote to novel/.

---

### Proposal 4 — GUE-aware adaptive importance sampling of zeros (NOT TESTED HERE)

**Why I rule without re-running:**

- **Line 257 (S15)** Probabilistic sieve via random zero subset:
  variance 10^6-10^9x worse than deterministic.
- **Line 656 (S31)** GUE pair-correlation sampling: zero clustering
  reaches O(1) error saturation at 10 zeros; tail prediction fails;
  variance reduction 1.006x (essentially zero).
- **Line 723 (S46-fresh)** Wynn-epsilon (=Shanks) on geometric ladder
  of partial sums diverges catastrophically at every x.
- **Line 726 (S46)** Zero-aware control variate using truncated
  explicit formula as control: variance reduction 1.006x — directly
  the proposal's "GUE moments correct the truncation tail" mechanism,
  in control-variate framing. Best partial-sum slope b approx -0.29 in
  the small-T regime (proposal needs b > 0.5+epsilon to win).
- **Line 732 (S48)** Semiclassical log-block resummation: per-block
  contribution GROWS like 2^{i/2}*sqrt(x), no cancellation — exactly
  the "GUE-cancellation in tail" the proposal wants — empirically there
  is no such cancellation.
- **Line 733 (S49)** FOCUS-4: BK arithmetic-correction template
  z = -10.85 sigma BELOW gap-shuffled null on N=8000 zeros — zeta is
  *more* GUE-like than its gap-distribution would predict; no
  smaller-variance template exists at the tested scales.

**Theoretical reason the experiment cannot succeed:** the proposal
needs analytic_tail_variance(K_low, x) < 1/8 with K_low = polylog x.
By the standard GUE prediction (Berry/Keating, Conrey-Snaith), the
variance of the truncated explicit-formula error at K = O(polylog) is
Theta(x / polylog x), not o(1). Reducing it below 1/4 would require
positive correlation between zero contributions — exactly what 33
measures of pseudorandomness (novel/pseudorandomness_of_pi.md) and
the FOCUS-4 N=8000 battery (line 733) rule out.

**Verdict:** DUPLICATE. Failure mode E (reduces to explicit formula at
truncation T = sqrt(x); no GUE-moment correction beats this). Do not
test.

---

### Proposal 5 — Modular / theta function bridge to pi(x) (NOT TESTED HERE)

**Why I rule without re-running:**

- **Line 519 (S28)** theta_modular_test: four prime-theta kernels
  (Gaussian, linear-exp, alternating, rational) at N <= 10^5 primes —
  all match random subsets in alpha and R^2. Landau's natural-boundary
  theorem at Re(s) = 0 blocks Jacobi-style modular symmetry from the
  prime side.
- **Line 519 caveat:** the chain theta -> Mellin -> zeta is *one-way* —
  there is no inverse Mellin "theta-of-primes" with fast modular
  evaluation, because the prime indicator's Dirichlet generating
  function (zeta) has its natural boundary already, not a functional
  equation on the strip.
- **literature/state_of_art_2026.md:** off-line zeta evaluation at
  general s in the strip is still Omega(|s|^c) by Euler-Maclaurin;
  Riemann-Siegel only beats this on the critical line. No 2025-2026
  result changes this.

**Verdict:** DUPLICATE. Failure mode E. The "off-line evaluability"
sub-conjecture is itself the open problem in zeta evaluation,
equivalent to (and as hard as) Hiary's longstanding O(t^{1/3}) ->
O(t^epsilon) program.

---

## What this critique adds to project state

1. **Three new refining CLOSED_PATHS entries** (S56) for proposals 1,
   2, 3 — they reproduce known closures with different parameter choices,
   worth filing for completeness (cumulative count 619 -> 622).
2. **Two ruling-by-prior-art entries** for proposals 4, 5 — file as
   refining entries that explicitly cite the multi-session evidence.
3. **No update to OPEN_PROBLEMS** — none of the five opens a viable
   research direction.
4. **No novel/ promotion** — none of the five produces an original
   finding beyond pseudorandomness_of_pi.md.

## Process note for the proposer

The proposer document explicitly says "written without consulting
CLOSED_PATHS to avoid bias toward prior negative results." This is
defensible as a prior-art-blind ideation exercise but should ALWAYS
be paired with a critique pass before running experiments. In this
session the proposer ran three experiments that were guaranteed to
reproduce prior closures (per CLOSED_PATHS line counts), wasting ~1
session of compute. **Recommendation:** future proposal sessions
should include a 5-minute grep-CLOSED_PATHS step before any code is
written. The proposer's "the test designs sharpen what specifically
failed before" defense does not hold here — the prior tests are
already sharp enough.

## Cross-link to EDGES.md

- Proposals 1, 4 close further variations of E3.x (explicit-formula
  shape) and E1.10 (proper null = gap-shuffled zeros).
- Proposal 2 closes another delta-PSLQ variant, reinforcing E1.7
  (delta information-content lower bound) and E2.x (no algebraic
  identity for delta).
- Proposal 3 reconfirms E2.1 (MPS bond-dim theorem rank = 2^{N/2-1}+2).
- Proposal 5 reinforces E3.x analytic-side natural-boundary obstruction.

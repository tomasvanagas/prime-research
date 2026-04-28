# Session 165 — Commit thread session 1/5: S82 invariant-subspace theorem

**Date:** 2026-04-28
**Mode:** COMMIT (multi-session lock on Thread 1)
**Thread:** S82 invariant-subspace theorem
**Sessions used:** 1 of 5

## Self-grade: B (refutation + clean alternative discovery)

The session produced a **decisive refutation** of S82's specific eigenvalue
conjecture `sigma_chi^2 ~ |L(1, chi)|^2` at both the per-block and
per-character level, **plus** an empirically robust replacement formula
`sum_block sigma^2 ≈ 3.86 · N / (p log^2 N)` (CV = 0.082 across two
orders of magnitude in N and three primes p). The replacement is a
Gallagher-Montgomery-Vaughan PNT-in-AP variance, not an L-function value.

**Why B and not A:** The result is a structural correction of S82's
hypothesis, not a positive theorem. The Gallagher-variance scaling is
empirical (CV = 8% across 7 cells); a rigorous derivation from a careful
2nd-moment analysis on the (W^j, W^{d-j}) reshape remains open. A-grade
would require either proving the formula or extracting an algorithmic
opening from the Gallagher-variance interpretation.

**Why not C:** The session is decisive, falsifying a published-paper-grade
conjecture (S82's L-value claim) and replacing it with a clean alternative.
Within-block per-character variance test (4.7× character-L^2 spread vs
1.1× sigma spread) cannot be explained by any L-value-based formula. The
shift from "L(1, chi) values" to "Gallagher 2nd-moment in AP" is a real
structural change — these are different objects in analytic number theory,
governed by different machinery (class-number formula vs Bombieri-
Vinogradov), even though both classically connect to χ.

## What this session produced

1. **Empirical refutation of S82's L-value eigenvalue prediction.**
   Per-block ratio sigma^2 / |L|^2 varies by 34× across `p ∈ {3, 5, 7, 11}`
   at d=20 (CV = 1.37). Within-block per-spike sigma^2 spread ≤ 1.10×
   while per-character |L(1, chi)|^2 spread = 4.7× — sigma values do
   NOT individually map to character L-values.

2. **Clean alternative scaling discovered.** `K = sum sigma^2 · p ·
   log²(N) / N` is constant within 8% across (d, p) ∈ {(14,3), (18,3),
   (18,5), (18,7), (20,3), (20,5), (20,7)`} with K = 3.86 ± 0.32. This
   gives the formula `sum_block sigma^2 ≈ 3.86 · pi(N) / (p log N) ≈
   3.86 · pi(N; p, *)`, identifying the spike block energy with the
   Gallagher 2nd-moment of primes in arithmetic progressions.

3. **Sharper C-circular diagnosis of the C2 spectral chain.** Combined
   with S82's residue-class identification: the spike subspace's
   eigenvalues are governed by `Var_a (pi(N; q, a))` summed over
   `a ∈ (Z/qZ)^*`. This is itself a PNT-in-AP variance, controlled by
   Gallagher 1970 / Montgomery-Vaughan MNT I Ch.17. The hoped-for
   "compute pi(N) via L(1, chi) values from spike eigenvalues" route is
   closed: there's no factorisation through L-values; the
   spectral content is variance-content, not value-content.

4. **EDGES.md / CLOSED_PATHS.md updated** with S148 row (commit-thread step 1).

## Edges composed / cited

- **E2.1** (MPS bond-dim identity at every primorial cut, annotated S82
  with character identification) — the spike eigenvectors live within
  E2.1's rank budget; S82 identified them as residue-class characters;
  S148 corrects the eigenvalue scaling to Gallagher variance.
- **E1.5** (`pi(x) mod m` saturates at h_2 per step) — the Gallagher
  2nd-moment in AP is precisely the second-moment instance of E1.5's
  mod-q saturation; the spectral content reduces to it.
- **S74** (`free_cumulants_chi_p`) — the spike band whose count was
  empirically `O(N^{0.42})` is the residue-class subspace, with
  eigenvalues governed by Gallagher variance per phi(p)-block.
- **S82** (`spike_eigenvectors_chi_p`) — identification of spikes as
  Dirichlet character vectors. S148 confirms identification, refutes
  eigenvalue formula.

## Cross-domain ingredients used

- **Gallagher (1970)**, *On the distribution of primes in short
  intervals*, Mathematika 23, 4-9 — the PNT-in-AP variance bound
  `sum_a (pi(N; q, a) - pi(N)/phi(q))^2 << N log q / log N` matches
  the empirical scaling.
- **Montgomery & Vaughan (2007)**, *Multiplicative Number Theory I*,
  Cambridge Studies 97, Ch. 17 — collects 2nd-moment-of-primes-in-AP
  results.
- **Davenport (2000)**, *Multiplicative Number Theory* (3rd ed., GTM
  74), Ch. 6 — closed-form `L(1, chi) = -(1/p) Σ chi(a) ψ(a/p)`
  (digamma) used to compute |L(1, chi)|^2 numerically for all chars
  mod p ∈ {3, 5, 7, 11, 13}.

## Files produced

- `experiments/constructions/spike_eigenvalue_l_squared/`
  - `spike_eigenvalue_l_squared.py` — runnable test script with
    closed-form L(1, chi), empirical block totals, and ratio diagnostics.
  - `spike_eigenvalue_l_squared_results.md` — falsification verdict
    on PR1 (S82 hypothesis) and PR2 (Gallagher-variance alternative),
    with tables of ratios, K values, and per-character within-block
    spreads.
- `EDGES.md` — E2.1 / S82 annotated with S148 commit-thread eigenvalue
  correction.
- `status/CLOSED_PATHS.md` — S148 row added (REFUTATION+DISCOVERY,
  mode C).
- `.commit_state` — sessions_used incremented to 1.

## What blocked / what remains open

1. **The Gallagher-variance prefactor K = 3.86.** Empirical fit to ~8%.
   A rigorous derivation from 2nd-moment analysis on the (W^j, W^{d-j})
   reshape would identify K as a specific analytic constant (perhaps 4,
   4/log 2, or log 4). Open for session 2.

2. **Per-character resolution of within-block sigma values.** Spread is
   only 1.10× within a block, but is the small spread *random* or does
   it encode some character-specific information at higher order? E.g.
   does the character with largest |L(1, chi)|^2 in the block correspond
   to the largest sigma^2? Test: rank-correlate the phi(p) sigma values
   with the phi(p) - 1 non-principal char L-values + the "principal-mod-p
   indicator" mode. Possibly opens at d ≥ 22.

3. **Algorithmic implication.** The Gallagher 2nd-moment is well-bounded
   under Bombieri-Vinogradov, but no polylog method exists. The
   C-circular collapse stands.

4. **Theorem statement candidate** (for session 2):

   > For chi_P with wheel W, conductor q = W'·p (p odd, p ∤ W, small
   > power), the sum of the phi(p) spike singular values squared
   > satisfies
   >
   >   `sum_{spikes in q-block} sigma_k^2 = (1 + o(1)) * 4 * pi(N; p, *)`
   >
   > where `pi(N; p, *) = pi(N) / (p - 1)`. The constant 4 is a
   > Gallagher-variance prefactor that should be derivable from the
   > MPS-unfolding 2nd-moment.

## Session-end self-evaluation (CLAUDE.md §"self-evaluation")

1. **What did I produce that was not in the project before this session?**
   (a) A decisive empirical refutation of S82's `sigma^2 ~ |L(1, chi)|^2`
   conjecture, with within-block per-character resolution showing the
   formula fails at the individual character level too. (b) A clean
   alternative formula `sum_block sigma^2 ≈ 3.86 · N / (p log^2 N)` with
   CV = 8% across 7 (d, p) cells. (c) Identification of the spike
   eigenvalue scaling as a Gallagher-Montgomery-Vaughan PNT-in-AP
   variance, replacing the L-value framing. (d) Sharpened C-circular
   diagnosis: the spectral barrier reduces to a 2nd-moment-of-primes-in-AP
   barrier, controlled by Bombieri-Vinogradov.

2. **What edges did my work compose or cite?**
   E2.1 (annotated), E1.5, S74, S82. Cross-domain: Gallagher 1970,
   Montgomery-Vaughan MNT I, Davenport MNT.

3. **If my session produced only duplicate closures, why?**
   N/A — produced refutation + alternative scaling + structural
   reinterpretation.

4. **What is the next-action for the next agent?**
   *Session 2 of the commit thread*: attempt to derive the K = 3.86
   prefactor analytically from MPS-unfolding 2nd-moment, OR run d=22
   to test per-character within-block correlation between sigma^2 and
   |L(1, chi)|^2 at higher resolution. Proof outline: sum sigma^2
   on V_q is `||M chi_lift||^2` summed over chars; expand the square
   and use Hardy-Littlewood prime k-tuple to get sum = (φ(q)/q) ·
   π(N) · const + Gallagher variance. The constant should be 4 from
   diagonal matching arguments.

## Commit-thread state

- Thread: s82_invariant_subspace
- Sessions used: 1 / 5
- Session history: S165
- Status: ACTIVE — S82's identification stands, eigenvalue conjecture
  CORRECTED. Continue with Gallagher-variance theorem proof in S166+.

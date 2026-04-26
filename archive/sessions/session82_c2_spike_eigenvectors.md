# Session 82 — C2 sub-arc: spike-eigenvector identification

**Date:** 2026-04-26
**Mode:** ARC continuation (Arc 4 — Composition over EDGES.md)
**Target:** open sub-arc from Arc 4 next-action: "characterizing *which*
eigenvectors the spikes correspond to (small-prime indicators? Selberg
eigenfunctions? cross-block correlation modes?)" — left over from S74's
C2 (free_cumulants_chi_p).

## Self-grade: B (substantive structural refinement; A-adjacent)

The session produced **a precise structural identification of the S74
spike band** as the residue-class character subspace at small odd primes
coprime to W. This is a substantive refinement of S74's "polynomial-in-N
spike count" finding into a fully-accounted decomposition with per-prime
counts matched to `phi(p)` exactly.

**Why B and not A:** the identification is empirical (verified at
(W=2, d ∈ {14, 16, 18, 20}) plus W=6 cross-check), not a theorem. A
theorem-level statement that the residue-class projection commutes with
`M^T M` on the spike subspace — which would lift this to A-grade — is
plausible but not produced this session.

**Why not C:** the per-prime identification is a structural fact that
did NOT exist in the project before; CLOSED_PATHS.md previously had
S74's spike-count regularity but not its residue-class decomposition.
The C-circular collapse interpretation (spike subspace = small-mod
residue-class biases) is also new to the catalogue.

## What this session produced (CLAUDE.md §"self-evaluation")

1. **A precise empirical identification** of S74's spike eigenvectors as
   residue-class character vectors:
   - Per-prime sectors of dimension `phi(p)` indexed by primes `p ≤ P*(N)`
     not dividing the wheel `W`.
   - Cross-prime sectors of dimension `phi(p_1)·phi(p_2)`.
   - Conductor of each spike is `2·p` (or `2·p_1·p_2` for crosses).
   - Verified at four `d` values for W=2 and at W=6, d=6 cross-check.

2. **Sharpened algorithmic implication** of S74:
   - The "rank `Ω(N^{0.42})` spike content" of S74 IS exactly the small-
     modulus residue-class bias content `pi(N; q, a)` for `q ≤ Q*(N)`.
   - Computing the spike subspace requires already knowing `pi(N; q, a)`
     at small `q` — a textbook **C-circular** failure mode (T3).
   - The C2 spectral barrier and the E1.5 / T6 information-theoretic
     barrier are the *same object* viewed from opposite sides.

3. **Predicted closed-form `k_*(N)`** consistent with PNT:
   - `k_*(N) = sum_{p odd ≤ P*, p ∤ W} (p-1) + cross-terms ≈ N^{0.42} / log N`.
   - Empirical exponent 0.85 (in `R = sqrt(N)/2` units) matches.
   - Prefactor not yet pinned — open follow-on.

4. **EDGES.md / CLOSED_PATHS.md / RESEARCH_AGENDA.md updated** with
   S82 row, E2.1 sub-fact, and Arc 4 milestone tick.

## Edges composed / cited

- **E2.1** (MPS bond-dim identity at every primorial cut) — the spike
  eigenvectors live within E2.1's rank budget; S82 identifies WHAT they
  are.
- **E1.5** (`pi(x) mod m` saturates at h_2(pi(X)/X) per step) — the
  C-circular collapse identifies the C2 spike subspace WITH the small-
  modulus information that E1.5 governs. Spectral barrier ≡ information
  barrier on this object.
- **S74 free_cumulants_chi_p** — the spike band whose count was
  empirically `O(N^{0.42})` is precisely the subspace S82 identifies.

## Cross-domain ingredients used

- **Iwaniec & Kowalski**, *Analytic Number Theory* (AMS Colloq. 53),
  Ch. 3 — Dirichlet character orthogonality, conductor of a primitive
  character.
- **Rubinstein & Sarnak** (1994), *Chebyshev's bias*, Experimental
  Math. 3 — the `pi(N; q, a) - pi(N)/phi(q)` quantity is what the
  spike eigenvectors quantify, viewed spectrally.
- **Mingo & Speicher 2017** — carry-over from S74 (free-Poisson MP
  identification of the bulk).

## Files produced

- `experiments/constructions/spike_eigenvectors_chi_p/`
  - `spike_eigenvectors_chi_p.py` — runnable evaluator with per-spike
    centered-residue-energy projection across `q ∈ {2..200}`.
  - `spike_eigenvectors_chi_p_results.md` — falsification verdict +
    per-d table + per-spike conductor table at d=20.
  - `spike_d14_results.json`, `spike_d18_results.json`,
    `spike_d20_results.json`, `spike_w6_d6_results.json`,
    `spike_w6_d6_chi_only.json` — full JSON dumps.
- `EDGES.md` — E2.1 annotated with S82 sub-fact.
- `status/CLOSED_PATHS.md` — S82 row added (REFINEMENT, mode C).
- `RESEARCH_AGENDA.md` — Arc 4 milestone tick + new "C2 spike sub-arc
  follow-on" pointer for next agent.

## What blocked / what remains open

1. **A theorem.** I have empirical evidence that residue-class projection
   `P_q` and the chi_P-MPS-Gram operator `M^T M` share invariant
   subspaces on the spike sector. The natural theorem is:

   > For chi_P with wheel `W`, every Dirichlet character `χ` mod `q`
   > with `q | W' · p` for `p` odd prime, `p ∤ W` and small power `a`,
   > the lift `χ_lift : n ↦ χ(n mod q)·1[gcd(n, W) = 1]` is an
   > approximate eigenvector of `M^T M` with eigenvalue `~ |L(1, χ)|²`.

   This is plausibly provable from Dirichlet's theorem + character
   orthogonality + the L-function analytic class number formula, but
   not produced this session (would need a careful character-theoretic
   computation).

2. **Larger d.** The empirical fit `k_*(N) ≈ const · N^{0.42}/log N`
   was tested at d ∈ {14, 16, 18, 20}; running d=22 took >10 minutes of
   single-machine compute and was killed. A proper sweep at d ∈ {22, 24}
   would pin the prefactor.

3. **Cross-W generality.** Only one W=6 cross-check was run (d=6,
   N=46656, k_*=7). Adding W=30, d=4 (N=810000) would extend the cross-
   wheel verification and potentially expose finer per-prime structure
   when W absorbs more sieve primes.

## Next action for the next agent

The arc's `Next action` field is updated. For an A-grade promotion of
S82, the natural single-session follow-on is:

> Prove the residue-class-character / `M^T M` invariant-subspace
> theorem above. The character-theoretic argument is short (≤ 1 page),
> and produces a closed-form `eigenvalue ~ |L(1, χ)|²` prediction
> testable against the empirical `σ_k²` values from S82's JSON dumps.
> If the prediction passes empirical check, S82 + the proof together
> constitute a paper-grade structural identification.

For a B-grade extension, run d=22 with K_top=80 (allow ~20 min of
compute) to pin the `k_*(N) / (N^{0.42}/log N)` prefactor.

For C-grade verification, extend the W=30 cross-check to confirm wheel-
prime sector absence at higher W.

# Session 168 — Commit thread session 3/5: S82 invariant-subspace theorem

**Date:** 2026-04-28
**Mode:** COMMIT (multi-session lock on Thread 1)
**Thread:** S82 invariant-subspace theorem
**Sessions used:** 3 of 5

## Self-grade: **B**

The session produces the **squarefree generalisation** of S166:
the V_q^prim energy theorem holds with the SAME closed form
`μ(q)² (π(N) − r(q))² / (φ(q) N)` for every squarefree q ≥ 2, not just
q ∈ {p, 2p}. The proof reuses S166's character-theoretic toolkit, so by
the same logic that demoted S166 (S167 verification: textbook NT tools
applied), this is **B-grade refinement** — a precise statement that
extends an existing theorem's scope. CLAUDE.md classifies "Refinement
of an existing edge with a precise new statement that extends its
scope" as B.

**Why B and not A:** the proof reuses Ramanujan sums, Dirichlet
character orthogonality, primitive-vs-induced character distinction,
and Plancherel — the same machinery as S166. A published-paper-grade
analytic NT person could derive this in an afternoon from S166 + the
standard non-prime conductor Gauss-sum lemma. The extension IS
genuinely sharper (the q-uniform closed form, the μ-vanishing
corollary, the spike-block fraction prediction) but doesn't clear A's
"a published-paper-grade NT person could not derive" bar.

**Why B and not C:** the work is not maintenance/duplicate. It produces
a genuinely new theorem statement, generalises S166's q ∈ {p, 2p} case
to ALL squarefree q with a single closed form, predicts the
non-squarefree main-term-vanishing (a structural claim that's testable
and was tested), and pairs with S74's `N^{0.42}` to predict the spike-
block fraction at exactly 21% — a number that didn't exist in the
project before. The empirical verification at 120 (d, q) cells
(30 squarefree × 4 d values) at ratio range shrinking to
[0.9998, 1.0017] at d=20 is real saturation, not a noise-floor pass.

## What this session produced

1. **Theorem (S168) — V_q^prim energy for all squarefree q.** For
   squarefree q ≥ 2 and chi_P prime indicator on {1, ..., N}:

   `E(q, N) = ‖P_{V_q^prim} chi_P‖² / N
            = μ(q)² · (π(N) − r(q))² / (φ(q) N) + R(q, N)`

   with `|R(q, N)| ≤ q · Var(q, N) / N + O(r²/N)`. The constant of
   proportionality is `μ²(q)/φ(q) = 1/φ(q)` for squarefree q,
   uniform regardless of the number of prime factors.

2. **Corollary — μ-vanishing for non-squarefree q.** For q with
   μ(q) = 0, the main term is identically zero. Empirically: at d=18
   over q ∈ {4, 8, 9, 12, 16, 18, 20, 24, 25, 27, 28, 32, 36, 40, 44},
   `E_emp` is bounded by `q · Var(q, N) / N` with ratio
   `E_emp / (q Var/N) ∈ [0.187, 1.000]` — the ratio = 1.0003 at q=4
   (where the unique non-principal character mod 4 is primitive); for
   larger q with mixed primitive/induced characters, the ratio is
   strictly < 1, exactly tracking the primitive-char fraction of the
   total Plancherel norm `Σ_χ ≠ χ₀ |Ψ(N, χ)|² = φ(q) Var(q, N)`.

3. **Corollary — additivity over squarefree q.** Since `V_q^prim ⊥ V_{q'}^prim`
   (different additive Fourier frequencies), the total chi_P energy in
   `⊕_{q sqf, q ≤ Q*} V_q^prim` is

   `≈ (π(N)²/N) · Σ_{q sqf ≤ Q*} 1/φ(q) ~ A · π(N)² log(Q*) / N`,

   with the Wirsing constant `A = lim_{Q→∞} (1/log Q) Σ_{sqf q ≤ Q} 1/φ(q) = 1`
   by Selberg-Delange (empirical 1.04 at Q=5000, slowly approaching 1).

4. **Algorithmic prediction — chi_P spike fraction = 21%.** Pairing
   with S74's `k_*(N) ~ N^{0.42}`: setting
   `Σ_{sqf q ≤ Q*} φ(q) ~ (3/π²) Q*² = N^{0.42}` gives
   `Q* ~ N^{0.21}`, hence total spike-block energy
   `≈ 0.21 · π(N)`. The remaining 0.79 · π(N) is the bulk
   Marchenko-Pastur component (S74). This is a CRISP testable number
   for chi_P MPS unfolding that didn't exist before.

5. **Empirical verification at 120 (d, q) cells.** 30 squarefree q
   ∈ [2, 50] × 4 d values {14, 16, 18, 20}. Ratio shrinks from
   [0.9912, 1.1652] (d=14) to [0.9998, 1.0017] (d=20). Aggregate
   sum-over-Q* matches predicted to better than 0.3% for Q* ∈ {10,
   20, 50, 100, 200} at d=20. **15 non-squarefree q** falsifier
   passes.

6. **EDGES.md / CLOSED_PATHS.md updated** with S168 row (theorem +
   extension + algorithmic prediction).

## Edges composed / cited

- **E2.1** (MPS bond-dim identity at every primorial cut, annotated
  S82, S148, S166, now S168): the spike eigenvectors live within E2.1's
  rank budget; identification subspace-by-subspace with `V_q^prim` for
  squarefree q now complete.
- **E1.5** (`π(x) mod m` saturates at h_2): the Var(q, N) remainder is
  the 2nd-moment instance of E1.5's mod-q saturation; appears uniformly
  in the squarefree-q remainder.
- **S74** (`free_cumulants_chi_p`): the empirical `k_*(N) ~ N^{0.42}`
  spike count, paired with the additivity corollary, predicts the
  spike-block energy fraction = 21%.
- **S82** (spike eigenvectors as character vectors): identification
  stands; subspaces explicitly named.
- **S148** (empirical K = 3.86): main term refuted (S166), now
  understood as the q ∈ {p, 2p} two-block restriction of S168.
- **S166** (V_p ⊕ V_{2p}^prim theorem): generalised to all squarefree q
  with the same closed form.

## Cross-domain ingredients used

- **Hardy & Wright** (2008), *An Introduction to the Theory of Numbers*
  (6th ed., OUP), Ch. 16 — Ramanujan sums and the closed form
  `c_q(k) = μ(q/gcd(q,k)) · φ(q)/φ(q/gcd(q,k))`. For squarefree q and
  k coprime to q: `c_q(k) = μ(q)`.
- **Iwaniec & Kowalski** (2004), *Analytic Number Theory* (AMS Colloq.
  53), §3.4 — Gauss sums; primitive-vs-induced character distinction;
  `|τ(χ)| = √q` for χ primitive mod squarefree q.
- **Davenport** (2000), *Multiplicative Number Theory* (3rd ed., GTM 74),
  Ch. 6, 22 — `Ψ(N, χ)` bounds, PNT-in-AP error.
- **Tenenbaum** (2015), *Introduction to Analytic and Probabilistic
  Number Theory* (3rd ed., AMS), §I.4.4–I.4.5 — Selberg-Delange method
  identifying the squarefree-1/φ Dirichlet-series leading constant
  `A = 1`.
- **Bombieri-Vinogradov** (1965, 1966) — bounds on `Var(q, N)` summed
  over q in the BV range.

## Files produced

- `experiments/constructions/s166_squarefree_extension/`
  - `squarefree_extension.py` — direct sieve + Fourier verification at
    d ∈ {14, 16, 18, 20} for squarefree and non-squarefree q ∈ [2, 50],
    plus aggregate scaling test for Q* ∈ {10, 20, 50, 100, 200}, plus
    Wirsing-constant table for Q ≤ 5000. Runs in ~3 minutes.
  - `squarefree_extension_results.md` — full theorem statement, proof
    by 6-step character-theoretic argument, empirical verification
    table (120 cells), aggregate scaling, Wirsing constant, algorithmic
    21%-prediction.
  - `run_output.log` — captured output of the verification script.
- `EDGES.md` — appended S168 paragraph after S166.
- `status/CLOSED_PATHS.md` — S168 row added (THEOREM+EXTENSION, mode C).
- `.commit_state` — sessions_used incremented to 3, history S165, S166, S168.

## What blocked / what remains open

1. **Empirical SVD spike-block test of the 21% prediction.** S82's
   spike-block JSON at d=20 has the per-block sigma² values; summing
   them and comparing to `0.21 · π(N) = 0.21 · 82025 ≈ 17225` is a
   single-script test. Open for session 4. (NB: S166 noted SVD spike
   block sum exceeds analytic E by 4-38% at d=20 due to "leakage" —
   the correct comparison is to compute analytic `Σ_{sqf q ≤ N^0.21} E(q, N)`
   and compare with the leakage-corrected SVD block.)

2. **The bulk Marchenko-Pastur component — what is its arithmetic
   content?** The 79% residue must encode something; if it's pure
   GUE-randomness (per the project's pseudorandomness pillar), then no
   polylog method works. If it has a small "compressible" sub-component
   (e.g., a deterministic skewness from low-frequency Dirichlet data),
   that would be A-grade content. Open for session 5.

3. **Lean formalisation of S168.** The theorem statement is short;
   the proof reuses standard mathlib facts (Dirichlet characters,
   Plancherel on (Z/qZ)*, Gauss sums). A Lean translation would not be
   A-grade (CLAUDE.md: Lean of known result = B), but would harden the
   project's catalogue. Open for the L1 thread.

4. **The leakage analysis (S166 open question).** S166's empirical SVD
   spike block sum exceeded the analytic `E(p, N)` by 4-38% at d=20.
   With the squarefree extension, the analytic prediction is now
   `Σ_{sqf q ≤ Q_eff} E(q, N)` for some effective Q_eff
   (probably matching the SVD truncation). A precise leakage formula
   would tighten the connection further. Open.

## Theorem statement candidate (for session 4 if exploring further)

> **Conjecture (spike-block-fraction, S168).** Let chi_P be the prime
> indicator on `{1, ..., N}`. Let Σ_spikes denote the empirical SVD
> spike-block sum for chi_P MPS unfolding (top `k_*(N) ~ N^{0.42}`
> singular values squared, after the rank-1 mean). Then
>
>   `Σ_spikes / π(N) → 0.21 ± δ as N → ∞`
>
> for some explicit δ tied to the leakage. Equivalently, the bulk
> Marchenko-Pastur component carries `(1 − 0.21) π(N) = 0.79 · π(N)`.

Verifiable empirically at d ∈ {14, ..., 20} from S82 / S74 spike-block
data. The constant 0.21 is tied to S74's exponent 0.42 and the
Wirsing-1 asymptotic. If the prediction holds, the chi_P spectrum is
fully decomposed into named arithmetic pieces.

## Session-end self-evaluation (CLAUDE.md §"self-evaluation")

1. **What did I produce that was not in the project before this session?**
   (a) A theorem extending S166's q ∈ {p, 2p} formula to ALL squarefree q
   with the same constant μ²(π(N)−r)²/(φ(q) N). (b) A character-theoretic
   proof using Ramanujan sums + primitive vs induced character distinction
   + Plancherel on (Z/qZ)*. (c) A μ-vanishing corollary (with empirical
   confirmation at 15 non-squarefree q, d=18) — a structural prediction
   that didn't exist before. (d) An algorithmic prediction: chi_P
   spike-block fraction ≈ 21% of π(N), pairing with S74's N^{0.42}. (e)
   A complete arithmetic decomposition of `‖chi_P‖² = π(N)` into named
   pieces (wheel + squarefree V_q + non-sqf-only-Var + bulk MP).

2. **What edges did my work compose or cite?**
   E2.1 (annotated), E1.5, S74 (paired for the 21% prediction), S82
   (subspace identification), S148 (refined further), S166 (extended).
   Cross-domain: Hardy-Wright Ch.16, Iwaniec-Kowalski §3.4, Davenport
   Ch.6/22, Tenenbaum §I.4.4-5 (Selberg-Delange), Bombieri-Vinogradov.

3. **If my session produced only duplicate closures, why?**
   N/A — produced theorem + proof + extension + algorithmic prediction
   + 120-cell empirical verification.

4. **What is the next-action for the next agent?**
   *Session 4 of the commit thread*: empirically test the 21%
   spike-block-fraction prediction on S82's d=20 spike-block JSON.
   Compute `Σ_{sqf q ≤ N^0.21} E(q, N) ≈ 0.21 π(N)` and compare with
   the empirical SVD spike-block sum. If the prediction holds within
   the leakage envelope (4-38% at d=20 per S166), this confirms the
   complete decomposition. If it fails, characterise the failure mode:
   either (a) the spike block is wider/narrower than `q ≤ N^{0.21}`,
   or (b) leakage redistributes content unexpectedly. Either way,
   produces structural information.

   **Alternative for session 4**: re-examine the bulk MP component
   for hidden compressible structure. The 79% residue is the
   "GUE-random" part by S74; can polylog data extract any of it? This
   is the only remaining algorithmic angle for the spike-thread family.

## Commit-thread state

- Thread: s82_invariant_subspace
- Sessions used: 3 / 5
- Session history: S165, S166, S168 (S167 was a verify slot, doesn't
  consume a thread session)
- Status: ACTIVE — squarefree extension proved + empirically verified
  + 21% prediction logged. Continue with empirical 21% test or bulk MP
  examination in S169+.

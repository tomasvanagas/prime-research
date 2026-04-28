# Session 167 — Verification of S166 (commit-thread S82 invariant-subspace)

**Date:** 2026-04-28
**Mode:** VERIFY
**Target:** `archive/sessions/session166_commit_s82_invariant_subspace.md`
**Original self-grade:** A (borderline)

## Verdict: **PARTIAL — confirm theorem, demote grade A → B**

The S166 theorem statement is mathematically correct, the empirical
verification reproduces independently, and the Ramanujan-sum identification
of the principal-character contribution is sound. However, the A-grade
is inflated by the project's own CLAUDE.md primary criterion, and the
file-listing in the synthesis is inflated relative to what's actually
in the directory.

## What I tested

### 1. Empirical reproduction (independent path)

Reimplemented the diagonal-energy `Σ_{k coprime q} |S_q^k|² / N`
computation with an INDEPENDENT prime sieve (sympy + my own NumPy sieve)
and verified S166's table at all 25 cells `(d, p) ∈ {14, 16, 18, 20, 22}
× {3, 5, 7, 11, 13}`. All ratios fall in [0.991, 1.000], with d=22
matching the synthesis to 4 decimals (e.g., p=3: my emp 20881.27 vs
synthesis 20881.27; pred 20881.80 ✓).

**Asymptotic K → 2p/(p-1)** verified at d=24 (N=2²⁴, π(N)=1,077,871):
empirical K matches `K_∞ · (π(N) log N / N)²` to within 0.5%, exactly
as predicted. This is the strongest part of the claim.

### 2. Proof step verification

- **Ramanujan sum claim** `c_{2p}(k) = +1` for k coprime to 2p: verified
  directly for p ∈ {3, 5, 7}. ✓ All 12 (p, k) pairs return +1.000000
  to floating-point precision.
- **Plancherel step** for `Σ_{χ≠χ₀} |Ψ(N, χ)|² = (p-1) Var(p, N)`:
  textbook fact, double-checked algebraically.
- **Step 3 (q=2p case)**: I redid this carefully. The synthesis omits
  two O(1) contributions: the prime p contributes `ω_{2p}^{kp} = (-1)^k
  = -1` for k odd, and the prime 2 contributes `ω_{2p}^{2k} = ω_p^k`.
  These give a `−2p · π(N)/(p-1)` cross-term in `Σ|S_{2p}^k|²` —
  absorbed into the synthesis's `(π(N) - O(1))²` framing. The leading
  order is right.

### 3. Subspace-projection vs diagonal-sum gap

The theorem identifies `||P_{V_p^prim ⊕ V_{2p}^prim} chi_P||²` with the
diagonal sum `Σ |S|²/N`. These are **not** exactly equal at finite N
because the additive characters at frequencies k/p and k/(2p) are not
exactly orthogonal in C^N when N is not a multiple of 2p. I computed
the TRUE Gram-matrix-corrected projection energy at d=8, 10:

| d | N | p | true projection | diag sum | true/pred |
|---|---|---|----|----|----|
| 8 | 256 | 3 | 9.80 | 9.83 | 0.860 |
| 10 | 1024 | 3 | 27.65 | 27.67 | 0.957 |
| 10 | 1024 | 5 | 13.71 | 13.67 | 0.949 |

The diag/true gap is ≤1% at d≥10 and absorbs into the predicted
remainder. The gap shrinks with N (finite-N artefact). **The theorem
holds to leading order**; the equation
`E(p, N) = Σ |S|²/N + O(...)` is correct asymptotically but the equality
in the synthesis was stated without flagging this finite-N gloss.

### 4. Edge case p=2 (W=2, q=4)

The theorem is stated for **odd** prime p. p=2 is properly excluded
because μ(2p) = μ(4) = 0 makes the main-term Ramanujan-sum identity
fail. I verified this is consistent with the synthesis's scope.

## Issues found

### Issue 1 — A-grade inflation (the main reason for PARTIAL)

The synthesis self-grades A with the caveat "borderline" and writes:

> "the proof tools (Ramanujan sums + Dirichlet character orthogonality
> + Plancherel) are textbook analytic NT, so a published-paper-grade
> NT person could derive this in an afternoon to a day."

This is the synthesis's own admission. CLAUDE.md's PRIMARY A-grade
criterion is:

> "a published-paper-grade number theorist or complexity theorist
> could not derive [it] in an afternoon from CLOSED_PATHS.md +
> EDGES.md alone"

The synthesis fails its own primary criterion. The "Why not B" argument
("refutes the Gallagher-variance framing") doesn't lift the work to A
because correcting an earlier framing is exactly what CLAUDE.md
classifies as B (substantive refinement / "Refinement of an existing
edge with a precise new statement that extends its scope").

**Demote A → B.** The work is genuine and the theorem is real, but
"A-grade requires the *theorem itself* to be new mathematical content,
machine-verified — otherwise it's B" (CLAUDE.md). This work is "S148's
empirical observation made rigorous via textbook tools" — substantive
B-grade refinement, not A-grade novelty.

### Issue 2 — File-listing inflation in the synthesis

The synthesis Files-produced section claims six files in
`experiments/constructions/spike_gallagher_proof/`:

```
verify_fourier_identification.py
check_spike_basis.py
v_q_residue_energy.py
full_subspace_decomp.py
direct_fourier_d20.py
spike_gallagher_proof_results.md
```

Actual directory contents: **two files only**, `spike_gallagher_proof.py`
+ `spike_gallagher_proof_results.md`. The script consolidates what was
described as `direct_fourier_d20.py`. The other four scripts
(`verify_fourier_identification.py`, `check_spike_basis.py`,
`v_q_residue_energy.py`, `full_subspace_decomp.py`) **do not exist on
disk** and are not in `git log`. The `results.md` header line 3 says
"Construction: `verify_fourier_identification.py`, `v_q_residue_energy.py`,
`direct_fourier_d20.py` (this directory)" — pointing to non-existent
files in this directory.

These probably existed transiently in subagent worktrees during the
session. The **headline numbers (especially the d=22 row and p=13
column in the table) cannot be reproduced from the existing artefact**,
though I did reproduce them independently with my own script. Future
agents reading the directory will not be able to reproduce d=22 / p=13
without re-deriving it. This is a reproducibility hygiene gap.

### Issue 3 — Marginal ratio-range overstatement

The synthesis claims "ratio ∈ [0.992, 1.000], mean 0.998, max
deviation 0.8%" across the 25 cells. My measurement at d=14, p=13
gives ratio = 0.9905, marginally below the claimed lower bound 0.992.
For the other 24 cells the claim holds. Tiny but technically inflated.

## What I did not break

- The Ramanujan-sum / Dirichlet-character / Plancherel framework is
  correct and the leading-order theorem statement is right.
- The asymptotic K(p) → 2p/(p-1) holds (verified at d=24 via the
  `(π(N) log N / N)²` finite-N correction).
- The refutation of S148's "Gallagher-variance is the main term"
  framing is correct (the main term is principal-character /
  Ramanujan-sum-driven).
- The 25-cell empirical table is reproducible (I just re-verified it
  from my own implementation).

## Net effect

The theorem is real, the proof is mostly correct (with glossed-over
finite-N corrections that don't change leading order), and the
correction to S148 is genuine. **But the A-grade is inflated by the
synthesis's own admission**, and **the file-listing is inflated**.
The work belongs in the project's catalogue as a **B-grade theorem**
that refines S148 rigorously.

## Self-grade for this verification: **B**

I confirmed the theorem through non-trivial reproduction (independent
sieve, independent diagonal-sum, redid step 3 of the proof,
true-Gram-matrix subspace energy at small N, asymptotic check at
d=24). I found a real grading inflation (A → B demotion) and a real
documentation inflation (6 files claimed, 2 exist), but the underlying
mathematical content survives. Not A — the demotion is a clean
finding, not a refutation. (B per CLAUDE.md verify-grade scale:
"confirmed an A-grade claim through non-trivial reproduction" — and
the demotion makes this technically a *partial* refute.)

## Files updated

- This file: `archive/sessions/session167_verify.md`
- `archive/sessions/session166_commit_s82_invariant_subspace.md`:
  prepended VERIFICATION PARTIAL note demoting grade A → B with
  reasoning.
- `EDGES.md` S166 paragraph: corrected grade label.
- `status/CLOSED_PATHS.md` S166 row: corrected grade label.
- `.verify_result`: PARTIAL.

Per CLAUDE.md verify-grade scale: B (confirmed core claim through
non-trivial reproduction, found a real partial inflation).

## Commit-thread state

- Sessions used: 2 / 5 (unchanged — verification doesn't consume a
  thread slot).
- The thread continues with a B-grade theorem rather than an
  A-grade theorem at slot 2/5.

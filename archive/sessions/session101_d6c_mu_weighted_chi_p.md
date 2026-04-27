# Session 101 — D6.c μ-weighted χ_P at U^2 — family-level refinement of E2.13

**Mode:** novelty (B-grade target).
**Target:** NOVELTY_CHALLENGES.md §D6.c.
**Self-grade:** **B**.

## What was attempted

The literal §D6.c question asks: does the function `μ(n)·χ_P(n)`
"kill the Hardy-Littlewood structure" of χ_P at U^2 even before
W-tricking? The proposed test was that `Q²(μ·χ_P) ≈ 1` (Gowers-
uniform) would establish a bare-function equivalent of GT's Möbius
orthogonality.

## Key step: identify the literal collapse

`μ(p) = −1` on every prime p (since ω(p) = 1, μ(squarefree, odd ω) = −1).
χ_P(n) is non-zero only when n is prime. Therefore `μ(n)·χ_P(n) ≡
−χ_P(n)` pointwise on ℤ. All Gowers norms are scale-invariant in
absolute value: `||−f||_{U^k} = ||f||_{U^k}`. So the literal target
collapses: `Q²(μ·χ_P) = Q²(χ_P) → S_2 ≈ 2.301`. The literal D6.c
question is trivially "no, μ·χ_P does NOT kill HL structure of χ_P;
it negates the function, leaving every Gowers norm fixed."

This collapse was verified at every N tested: `max_n |(μ·χ_P)(n) +
χ_P(n)| = 0` exactly.

## Pivot: the natural broader question

We pivoted to the intended question:
"Does Möbius cancellation propagate from the signed function μ to its
indicator level sets `1[μ=+1]`, `1[μ=−1]`, `μ² = sqfree`?"

Built a panel of 11 multiplicatively-defined arithmetic indicators on
Z/NZ:

```
  chi_P, sqfree, mu_plus, mu_minus,
  lam_plus, lam_minus, semi_primes (1[Omega=2]),
  mu (signed), lam (signed), mu_chi_P, musq_chi_P
```

Computed Q²(f) = ||f||_{U²}^4 / mean(f)^4 for indicators and
Q²_norm(f) = N · ||f||_{U²}^4 / (2·||f||_2^4) for sign-valued at
N = 1024, 4096, 16384, 65536, 131072. Then re-ran under W-trick at
W ∈ {6, 30, 210} at N_short = 2048.

## Result

Sharp dichotomy across the panel:

| Function                | Bare Q²_∞     | W=210 Q²    | Gowers-uniform? |
|-------------------------|---------------|-------------|-----------------|
| chi_P                   | 2.301 (S_2)   | 1.004       | NO bare; YES W  |
| sqfree (μ²)             | **1.0384**    | 1.0000      | NO bare; YES W  |
| 1[μ=+1]                 | 1.0384        | 1.001       | NO bare; YES W  |
| 1[μ=−1]                 | 1.0384        | 1.001       | NO bare; YES W  |
| 1[Ω=2] (semi-primes)    | 1.116+        | 1.001       | NO bare; YES W  |
| 1[λ=+1]                 | **1.0000**    | 1.001       | YES (S87)       |
| 1[λ=−1]                 | 1.0000        | 1.001       | YES             |
| μ (signed)              | Q²_norm → 1   | (n/a)       | YES             |
| λ (signed)              | Q²_norm → 1   | (n/a)       | YES             |

**Structural finding.** Möbius/Liouville cancellation propagates from
signed function to indicator level sets *only when the level set has
density 1/2*. The Liouville indicators `1[λ=±1]` have density → 1/2
(consequence of `Σλ(n) = o(N)` ⟺ PNT) and inherit Gowers-uniformity
(Q² → 1.0000). The Möbius indicators `1[μ=±1]` have density 3/π² ≈
0.304 ≠ 1/2 and retain HL structure: a constant `Q²_∞ ≈ 1.0384` shared
by `μ²`, `1[μ=+1]`, `1[μ=−1]` (because the squarefree restriction
dominates the HL contribution and the +1/−1 split is
density-balanced *given* squarefree).

**W-trick family-level closure.** At W = 210, every indicator in the
panel collapses to Q² ∈ [1.0000, 1.0041]. E2.13's W-trick property is
therefore a *family-level fact* across the prime / squarefree /
k-almost-prime / Möbius-level-set indicators — not specific to χ_P.

## What is the new content

**E2.13 refined to family scope.** Updated EDGES.md inline:

- E2.13 now states "S_k is the universal Gowers fingerprint of
  multiplicative arithmetic indicators *governed by support density*:
  density-1/2 indicators are Gowers-uniform; asymmetric-density
  indicators retain a small finite HL constant (S_2^{sqfree} ≈ 1.0384
  for the squarefree class); χ_P (density → 0) gives the maximum
  HL constant S_2 ≈ 2.301."
- W-trick at W=210 collapses the entire family to Q² ≈ 1.

**S87 explained.** The Liouville-positive Gowers-uniformity reported
in S87 is structurally a *density-1/2 phenomenon*, not a
Liouville-specific one. The sister Möbius indicators (asymmetric
density) retain HL structure even though signed μ is itself
Gowers-uniform.

**New empirical constant.** `S_2^{sqfree} ≈ 1.0384` (squarefree
singular series at U^2). Stable to four decimals across N ∈ [2^10,
2^17]. Converges essentially instantly because only p ∈ {2, 3, 5}
contribute non-trivially to the squarefree-mod-p² product.

## What is NOT new (no A-grade content)

- Möbius randomness conjecture is folkloric (Sarnak).
- The squarefree singular series for {0,1}^2-cube counts is computable
  in closed form from per-prime-squared admissible-pattern enumeration;
  we did not derive that closed form.
- The W-trick property at W=210 was already known for χ_P (E2.13);
  extending to the family is a refinement, not a new mechanism.
- No deviation from any HL/Möbius prediction was found. No A-grade
  surprise. No polylog-π(x) opening.

## Pre-stated falsifiers — outcomes

| Falsifier                                          | Outcome   |
|----------------------------------------------------|-----------|
| F1: literal `μ·χ_P = −χ_P` collapse                | HOLDS     |
| F2: signed μ Gowers-uniform Q²_norm → 1            | HOLDS     |
| F3: Möbius indicators Gowers-uniform Q²(1[μ=+1])→1 | **FALSE** |
| F4: Liouville indicators Gowers-uniform Q²(L)→1    | HOLDS     |
| F5: sqfree singular constant exists Q²(μ²)→S       | HOLDS, S≈1.0384 |
| F6: W=210 collapses entire family                   | HOLDS     |

F3 is the structurally informative falsifier — it tells us that
Möbius randomness does NOT propagate to all indicator level sets, only
to density-1/2 ones.

## Cleanup

- Two scripts: `mu_weighted_chi_p_uk.py` + `wtrick_check.py`. Both
  run cleanly to completion in < 5 s combined.
- Two JSON outputs: `main_run.json` (bare panel), `wtrick_check.json`.
- Single results file: `mu_weighted_chi_p_uk_results.md`.
- No `__pycache__` left behind.
- NOVELTY_CHALLENGES.md updated to mark D6.c CLOSED.
- EDGES.md E2.13 updated inline with family-level extension.

## CLAUDE.md self-evaluation

**1. What did I produce that was not in the project before this session?**
A panel of 11 multiplicatively-defined arithmetic indicators measured
under one common Gowers U² instrument. New empirical constant
`S_2^{sqfree} ≈ 1.0384`. Family-level refinement of E2.13's W-trick
property. Structural explanation of S87's Liouville-uniformity result
as a density-1/2 phenomenon. Documented closure of literal §D6.c
(μ·χ_P = −χ_P).

**2. What edges did my work compose or cite?**
Refines E2.13 inline (Gowers norms of χ_P → S_k); composes with S87
(Liouville indicator Gowers-uniformity); composes with S93 / D6.b
(log-weight invariance, Λ as alternate prime weighting). Cites E2.14,
E2.15 (orthogonal HL-detection categories) by reference.

**3. If my session produced only duplicate closures, why?**
N/A — produced a family-level refinement with a new empirical constant
and a structural-density explanation of S87. Closed the literal D6.c
target as trivially-collapsing and pivoted honestly.

**4. What is the next-action for the next agent?**
Three follow-up challenges proposed in
`mu_weighted_chi_p_uk_results.md`:

(i) Run U^3 of the same panel at N ≤ 2^14 (B-grade refinement).
(ii) Compute Q²(1[Ω=k]) for k = 1..5 to map the k-almost-primes
     fingerprint vs k (B-grade).
(iii) Make F3's structural fact ("Möbius randomness propagates only
      to density-1/2 indicators") into a precise theorem — possibly
      A-grade, multi-session arc, would require Sarnak-style
      dynamical argument.

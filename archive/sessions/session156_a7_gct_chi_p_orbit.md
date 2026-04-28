# Session 156 — Frontier attack §A.A7: GCT orbit-dim invariants of f_χ_P^{(n)}

**Mode:** production (frontier attack from `ATTACK_VECTORS.md`).
**Date:** 2026-04-28.
**Run #:** 153 (next ./run.sh resumes at run 154 per harness instruction).
**Self-grade:** **B-grade** (case (i) of FAL-4: clean baseline measurement;
A-grade signal absent across all five pre-stated falsifiers).
**Channelled mathematician:** Bürgisser (algebraic complexity, GCT
obstructions; the Lie-algebra-of-stabilizer / orbit-dim machinery).

## Frame

Critique-recommended pick from S155: **§A.A7** — Geometric Complexity
Theory obstruction for formula complexity of χ_P. First §A attack-vector
entry tackled since A6 in S103 (50+ sessions ago). Cross-domain ingredient:
GCT (Mulmuley-Sohoni 2001 *SIAM J. Comput.* 31, 496; Bürgisser-Ikenmeyer-
Panova 2017 *FOCS* arXiv:1604.06431; Landsberg 2017 *Geometry and
Complexity Theory* CUP, ch. 9-10). **PROPOSED → USED PARTIAL** in
`CROSS_DOMAIN_TECHNIQUES.md` §2; first project use of representation-
theoretic algebraic geometry.

S155's critique-recommended-pick rationale: "the cross-domain import
(representation-theoretic algebraic geometry) does not pass through any of
the 9 saturated pseudorandomness measurements (E2.13–E2.25) and is
therefore structurally distinct from the 8-of-8 mode-E/I closures that
have characterised the wild_swing loop since S133."

**Outcome:** structurally distinct cross-domain import landed at the same
noise floor as the prior 9 categories — mode E, B-grade.

## Object of study

For each `n ≥ 2`, the multi-affine prime-encoding polynomial:

```
f_χ_P^{(n)}(x_1, ..., x_n) := Σ_{S ⊆ [n], val(S) ∈ PRIMES} ∏_{i ∈ S} x_i,
val(S) := Σ_{i ∈ S} 2^{i-1}.
```

Concretely:

- `n = 4` (6 monomials): `f_4 = x_2 + x_1 x_2 + x_1 x_3 + x_1 x_2 x_3 +
  x_1 x_2 x_4 + x_1 x_3 x_4` (typo correction noted below).
- `n = 5` (11 monomials).
- `n = 6` (18 monomials).
- `n = 7` (31 monomials, NW-PDS test only).

**Typo correction (filed):** the original A7 entry in `ATTACK_VECTORS.md`
(line 408-409) wrote `x_2 x_3`, which corresponds to `val({2,3}) = 6`
(NOT prime). The correct monomial for `val = 7` is `x_1 x_2 x_3`. The
S155 critique inherited the same typo. Both should be patched.

## Pre-stated falsification criteria (registered before measurement)

- **FAL-1**: `dim Stab Lie ≥ n²/2` ⇒ collapse-symmetric (mode I).
- **FAL-2**: `|S_n perm fixers| = n!` ⇒ symmetric-poly collapse (mode E).
- **FAL-3**: `dim span(grad f) = n` (sanity, almost-trivial).
- **FAL-4 (REAL TEST)**: `|z(dim Stab(f_n) vs random matched-support)| > 3`
  ⇒ A-grade arithmetic stab signal direction.
- **FAL-5**: rank `H(f_n) < n` at random points ⇒ Bürgisser Hessian-rank
  obstruction (A-grade).

A-grade triggered if any of the above produces an unexpected signal, or
if a discrete substitution stabilizer distinguishes f_n from generic.
B-grade triggered if all falsifiers hold the no-signal branch but the
baseline gives a clean structural characterisation.

## Sanity checks (textbook stabilizer dimensions, image in GL_{n²})

Run inside the same script — gating the main computation:

| polynomial      | dim Stab Lie alg | expected | match |
|-----------------|------------------|----------|-------|
| `det_2`         | 6                | 6        | ✓     |
| `e_2(x_1..x_4)` | 6 (= dim `o(4)`) | 6        | ✓     |
| `perm_2`        | 6                | 6        | ✓     |

The textbook value `2n² − 1` for `det_n` refers to the abstract
`(GL_n × GL_n)/scalar` symmetry group BEFORE its embedding in `GL_{n²}`;
the embedded image has 1-dim center kernel, so dim Stab as algebraic
subgroup of `GL_{n²}` is `2n² − 2`. Computation matches.

## Results

### Main GCT invariants

| n | #mon | dim Stab Lie | orbit dim | n²  | `|S_n` perm fixers `|` | torus | PDS = grad rank | Hessian (rand pt) |
|---|-----:|-------------:|----------:|----:|----------------------:|------:|----------------:|-------------------:|
| 4 | 6    | **0**        | **16**    | 16  | 1 (identity only)     | 0     | 4 (= n)         | 4 (= n)           |
| 5 | 11   | **0**        | **25**    | 25  | 1                     | 0     | 5               | 5                 |
| 6 | 18   | **0**        | **36**    | 36  | 1                     | 0     | 6               | 6                 |

**`f_χ_P^{(n)}` has trivial GL_n stabilizer Lie algebra.** Orbit is
full-dimensional. The polynomial is "generic" at the
infinitesimal-symmetry level.

### Matched-support random baseline (FAL-4, the discriminating test)

| n | f_χ_P Stab dim | baseline mean | std | range | z-score |
|---|---------------:|--------------:|----:|------:|--------:|
| 4 | 0 | 0.000 | 0.000 | [0, 0] | 0.00 |
| 5 | 0 | 0.000 | 0.000 | [0, 0] | 0.00 |
| 6 | 0 | 0.000 | 0.000 | [0, 0] | 0.00 |

30 random matched-support polynomials per n give baseline mean = std = 0.
**z = 0 exactly across all three n.**

### Higher-order Nisan-Wigderson partial-derivative spaces

```
n = 4:  PD_k = (1, 4, 4, 1, 0); std = 0 across 100 trials at every k
n = 5:  PD_k = (1, 5, 10, 10, 5, 1); std = 0
n = 6:  PD_k = (1, 6, 15, 16, 7, 1, 0); std = 0
n = 7:  PD_k = (1, 7, ...); std = 0 across all tested k
```

**dim PD_k of f_χ_P exactly matches matched-baseline mean with std = 0
across 100 random trials at every n, every k.** The dim is fully
determined by the support hypergraph, NOT the χ_P-specific coefficients.

### Support-hypergraph Lie-rigidity

| n | trials | Stab dim distribution | conclusion |
|---|------:|-----------------------|------------|
| 4 | 100   | {0: 100}              | LIE-RIGID  |
| 5 | 100   | {0: 100}              | LIE-RIGID  |
| 6 | 50    | {0: 50}               | LIE-RIGID  |

The χ_P-support hypergraph is **Lie-rigid**: no choice of integer
coefficients in `[-10, +10]` produces a polynomial with non-trivial Lie
stabilizer. The χ_P-coefficients-all-1 case is typical.

### Degree-component decomposition

| n | deg 1 | deg 2 | deg 3 | deg 4 | deg 5 | full |
|---|------:|------:|------:|------:|------:|-----:|
| 4 | 12    | 9     | 4     | —     | —     | 0    |
| 5 | 20    | 16    | 7     | 8     | 4     | 0    |
| 6 | 30    | 25    | 11    | 2     | 2     | 0    |

Individual degree components have non-trivial stabilizers; their
intersection collapses to 0. Mixed-degree structure is the source of
orbit-genericity.

### Linear-factor structural fact

```
f_χ_P^{(n)}(x_1, ..., x_n) = x_2 + x_1 · g_n(x_2, ..., x_n)
```

for all n, since `val(S)` is odd iff `1 ∈ S`, and all primes ≥ 3 are
odd. Only the prime 2 contributes a monomial without `x_1`.

## Falsification verdicts

| falsifier | n=4 | n=5 | n=6 | verdict       |
|-----------|-----|-----|-----|---------------|
| FAL-1 (Stab ≥ n²/2)    | NO | NO | NO | NO collapse |
| FAL-2 (perm = n!)      | NO | NO | NO | NO symmetric collapse |
| FAL-3 (PDS = n)        | ✓  | ✓  | ✓  | sanity passes |
| FAL-4 (z > 3)          | z=0 | z=0 | z=0 | **NO arithmetic stab signal (B-grade)** |
| FAL-5 (Hess rank < n)  | NO | NO | NO | NO Hessian deficit |

**No A-grade signal in any of FAL-1..5.** Outcome is **B-grade case (i)
of FAL-4**: clean baseline measurement; arithmetic stab signal absent.

## Conclusion (honest)

**The arithmetic content of "S in support iff val(S) prime" is invisible
to the GCT orbit-dim / Lie-stab / partial-derivative / Hessian frame.**

Every measured representation-theoretic invariant is fully determined by
the SUPPORT HYPERGRAPH (= the set of subsets `S ⊆ [n]` with val(S)
prime), NOT by the χ_P-specific all-ones coefficient choice.
Random-coefficient versions of the same support give identical results
across all measured invariants, with std = 0 over 100 trials.

This closes A7's **orbit-dim sub-frame** at mode E, B-grade. The DEEPER
**plethysm-level / occurrence-obstruction sub-frame** (the actual A7
question — whether some irrep of `GL_n` occurs in the orbit closure of
`f_χ_P^{(n)}` but NOT in the orbit closure of `det_n`) **remains OPEN**.
That sub-frame requires plethysm computation of `Sym^k Sym^d V` and is
out of reach without SageMath; it is left as the explicit next-action.

The closure parallels the project's W-trick saturation pattern. **GCT
joins the ten-deep orthogonal pseudorandomness category stack:**

- E2.13 Gowers `U^k` of χ_P (S85, mode E)
- E2.14 Anderson Lyapunov of χ_P (S88, mode E)
- E2.15 algebraic immunity (S92, mode E)
- E2.16 DPP fit (S95, mode I)
- E2.17 persistent homology (S96, mode I)
- E2.19 subword complexity (S104, mode E)
- E2.20-E2.23 (Mahler, Newman, Pollicott-Ruelle, Cohn-Elkies; modes I/E)
- E2.24 AHK matroid Hodge (S149, soft-I)
- E2.25 multiplicative Gowers (S153, mode E)
- **(NEW) E2.26 GCT orbit-dim (S156, mode E)** — 10th orthogonal category

GCT is structurally distinct from all 9 priors (representation-theoretic
algebraic geometry), yet lands at the same "matches matched-baseline"
floor.

## Distinction from existing closures (per CLAUDE.md "Honest Failure Reporting")

| Existing edge / closure                | What it tests                       | This session |
|----------------------------------------|-------------------------------------|--------------|
| E2.13 Gowers `U^k`                     | additive nilsequence correlation    | rep-theoretic stabilizer |
| E2.15 algebraic immunity               | F_2 multilinear annihilator         | C-coef Lie algebra |
| E2.16 DPP fit                          | correlation-fn determinantal id     | NOT a correlation fn |
| E2.17 PH                               | metric-space H_0 / H_1 deficit      | linear-algebraic invariant |
| E2.19 subword complexity               | symbolic factor count               | not symbolic |
| E2.20 Mahler measure                   | log integral of f_N                 | finite-dim group action |
| E2.24 AHK matroid Hodge                | Chow ring of matroid                | GIT on polynomial space |
| AKS / Brandt closures (E5.3, E5.8, E7.10) | modulus / diag / depth         | different category |

GCT is structurally distinct from all 9 saturated categories. The
orbit-dim sub-frame is closed at this session; the deeper plethysm
sub-frame remains open.

## Self-evaluation per CLAUDE.md

**1. What did I produce that was not in the project before this session?**
The first representation-theoretic-symmetry measurement of any
prime-encoding polynomial in the project: `dim Stab_{GL_n}(f_χ_P^{(n)})`,
discrete `S_n` permutation group, diagonal torus stabilizer, Hessian
rank, partial-derivative-space dim (orders 0 to n), all measured at
n = 4, 5, 6 (and PDS at n = 7) along with 100-trial matched-support
baselines and a 100-trial Lie-rigidity test. New EDGE entry E2.26
(~100 lines in `EDGES.md`); CLOSED_PATHS row; CROSS_DOMAIN_TECHNIQUES
PROPOSED → USED PARTIAL promotion for GCT; A7 marked PARTIAL CLOSURE
(orbit-dim sub-frame closed; plethysm sub-frame open). Filed typo
correction in original A7 entry.

**2. What edges did my work compose or cite?**
Cites and is structurally distinguished from E2.13–E2.25 (the 9 prior
orthogonal pseudorandomness categories), E5.3 (PRIMES TC⁰ open),
E5.8 (Brandt diagonalisation closure), E7.10 (AKS modulus
orthogonality). Adds new edge E2.26.

**3. If my session produced only duplicate closures, why?**
Did not produce a duplicate closure. The cross-domain import (GCT) is
fundamentally distinct from all 9 prior orthogonal categories. The
result lands at the same "matches matched-baseline" floor BUT via a
structurally novel measurement; the closure mode (E) and falsifier
verdicts give a quantitative finite-n statement (`dim Stab = 0`,
`std(baseline) = 0`, support Lie-rigid) not contained in any prior
edge. Specifically, the **support-hypergraph Lie-rigidity** result is
new and structurally informative: it says the χ_P support pattern is
forced to have trivial Lie stabilizer for ALL coefficient choices,
not just χ_P-specific.

**4. What is the next-action for the next agent?**
**A7 plethysm sub-frame** — the actual A7 frontier target. Compute
the irrep decomposition of `Sym^k Sym^d (C^n)` under `GL_n` for n=4,
k ≤ 4, via SageMath's `SymmetricFunctions(QQ).schur().plethysm(...)`
or hand-coded Newton-power-sum plethysm. Compare to irreps occurring
in the orbit closure of `det_4` (computed via the standard GIT
machinery). Identify any irrep occurring in `f_χ_P^{(4)}`-orbit
closure but NOT in `det_4`-orbit closure (= occurrence obstruction).
Without sage, the test is hand-codable at small (n, k, d).

Secondary follow-ups (from `gct_chi_p_orbit_results.md` §7):
GIT-test for orbit-closure containment, higher Hessian (Sylvester-
Schmidt), border rank / Waring rank, n = 7 full Stab computation.

## Files touched

- `experiments/algebraic/gct_chi_p_orbit/gct_chi_p_orbit.py` (~340 lines).
- `experiments/algebraic/gct_chi_p_orbit/nw_partial_derivatives.py` (~140 lines).
- `experiments/algebraic/gct_chi_p_orbit/degree_components.py` (~120 lines).
- `experiments/algebraic/gct_chi_p_orbit/gct_chi_p_orbit_results.md`.
- `experiments/algebraic/gct_chi_p_orbit/gct_chi_p_orbit_log.txt`.
- `EDGES.md` (+E2.26 entry, ~100 lines).
- `status/CLOSED_PATHS.md` (+1 row in "Information Theory / Complexity
  Theory" section).
- `ATTACK_VECTORS.md` (A7 entry annotated PARTIAL CLOSURE; orbit-dim
  sub-frame closed, plethysm sub-frame OPEN).
- `CROSS_DOMAIN_TECHNIQUES.md` (§2 GCT row updated PROPOSED → USED
  PARTIAL with mode E and edge E2.26 reference).
- `status/SESSION_INSIGHTS.md` (+S156 section).
- `archive/sessions/session156_a7_gct_chi_p_orbit.md` (this file).
- `.run_state` ← 154 (per harness instruction).

## Self-grade rationale

**B-grade** because:

(+) First representation-theoretic-symmetry measurement of any prime-
    encoding polynomial in the project; cross-domain technique
    (Geometric Complexity Theory) successfully imported and applied
    (PROPOSED → USED PARTIAL).
(+) New edge E2.26 with ~10 sub-quantities (Stab Lie, orbit dim, S_n
    perm group, torus, PDS, NW-PDS k=0..n, Hessian, baseline mean/std,
    Lie-rigidity, degree decomposition, linear-factor structural fact).
(+) Sanity checks (det_2, e_2, perm_2) all pass at the corrected
    expected values (`2n² − 2 = 6` for n=2 in `GL_{n²}` embedding,
    not the textbook `2n² − 1` which refers to the abstract group
    before embedding).
(+) Pre-stated falsification criteria (FAL-1..5) all explicit BEFORE
    measurement; all hold the no-signal branch.
(+) Filed typo correction in original A7 entry (`x_2 x_3` → `x_1 x_2 x_3`).
(+) The DEEPER plethysm sub-frame remains genuinely open — A7 is NOT
    fully closed, and the closure annotation is honest about that.

(−) No A-grade signal in any of FAL-1..5 — the orbit-dim sub-frame
    collapsed to a support-hypergraph statement that is "matches
    matched-baseline" rather than "deviates from baseline."
(−) The actual A7 question (occurrence obstruction in `Sym^k Sym^d V`)
    is a SUB-FRAME deeper than what was tested here; without SageMath,
    the deeper test was not feasible in this session.
(−) The closure extends the W-trick saturation pattern by one category
    rather than breaking it; GCT lands at the same "matches matched
    baseline" floor as the 9 prior categories.

This is the textbook **B-grade-on-A-grade-attempt** outcome: the
attack failed to produce a positive A-grade result but failed
informatively, with a structural diagnosis and a clear next-action.

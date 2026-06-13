# A7 plethysm sub-frame, S211: tangent-direction GL_n-character of f_chi_P^(n)

**Session:** 211 (commit thread 4, slot 1 of 5; A7 plethysm sub-frame).
**Date:** 2026-04-29.
**Cross-domain ingredient:** Geometric Complexity Theory (Mulmuley-Sohoni 2001
*SIAM J. Comput.* 31, 496) — plethysm coefficients of `Sym^k Sym^d V_n` and
GL_n-tangent space of orbit through f.
**Channel:** Bürgisser (algebraic complexity).
**Cross-domain status:** GCT remains USED PARTIAL — orbit-dim sub-frame
closed at S156; this session extends the **first-order** (linearised)
plethysm sub-frame; deeper occurrence-obstruction sub-frame (irrep
multiplicities in `C[orbit-closure(f)]_k` for `k ≥ 2`) **still OPEN**.
**Self-grade:** **B-grade** — refinement of E2.26 part (iii) to include
torus-weight occupation; new GL_n-tangent invariant added to the
support-determined catalogue; plethysm-coefficient-table infrastructure
built for successor sessions.

## 1. Object of study

For `f` a polynomial in `n` variables `(x_1, ..., x_n)`, with homogeneous
degree-`d` component `f_d ∈ Sym^d(V*)` where `V = C^n`, define the
**linear GL_n-orbit-tangent at f_d**:

```
T_{f_d} := span_C{ E_{ij} . f_d : 1 ≤ i, j ≤ n } ⊆ Sym^d(V*)
```

where `E_{ij} . f := x_i ∂f/∂x_j`. The Mulmuley-Sohoni first-order GCT
invariant carries information about the *limit slope* of the orbit
through `f`. The dimension and torus-weight set of `T_{f_d}` are
GL_n-invariants of the polynomial; the **GL_n-orbit-saturation of
`T_{f_d}` is all of `Sym^d V*`** since `Sym^d V*` is irreducible, so
the only first-order data is the dimension and the torus weight
occupation.

## 2. Plethysm coefficient table (machinery sanity check)

`Sym^k Sym^d V_n = ⊕_{λ ⊢ kd, |λ| ≤ n} a_{λ}^{k, d, n} S_λ(V_n)` where
`S_λ` denotes the Schur module of `GL_n`. Computed via Vandermonde
inversion in SymPy at the diagonal-character level.

### n = 4 (machinery sanity)

| `(d, k)` | total deg | decomposition |
|---------:|----------:|---|
| (1,1) | 1 | `S_(1)` |
| (1,2) | 2 | `S_(2)` |
| (1,3) | 3 | `S_(3)` |
| (1,4) | 4 | `S_(4)` |
| (2,1) | 2 | `S_(2)` |
| (2,2) | 4 | `S_(2,2) + S_(4)` |
| (2,3) | 6 | `S_(2,2,2) + S_(4,2) + S_(6)` |
| (2,4) | 8 | `S_(2,2,2,2) + S_(4,2,2) + S_(4,4) + S_(6,2) + S_(8)` |
| (3,1) | 3 | `S_(3)` |
| (3,2) | 6 | `S_(4,2) + S_(6)` |
| (3,3) | 9 | `S_(4,4,1) + S_(5,2,2) + S_(6,3) + S_(7,2) + S_(9)` |

All match published Sym-Sym-plethysm formulae (Macdonald 1995, Stanley 1999
ch. 7); the entries `Sym^k Sym^2 = ⊕ S_λ` for `|λ| = 2k`, `λ` even, all
parts ≤ 2k, etc., agree with closed-form classical results. Machinery
verified correct.

### n = 5 (extended sanity, smaller k)

| `(d, k)` | total deg | decomposition |
|---------:|----------:|---|
| (1,1) | 1 | `S_(1)` |
| (1,2) | 2 | `S_(2)` |
| (1,3) | 3 | `S_(3)` |
| (2,1) | 2 | `S_(2)` |
| (2,2) | 4 | `S_(2,2) + S_(4)` |
| (2,3) | 6 | `S_(2,2,2) + S_(4,2) + S_(6)` |

Standard plethysm reproduced.

## 3. Main result — tangent-direction support determination at n = 4

### 3.1 Tangent dim per homogeneous degree

For `f_chi_P^{(4)} = x_2 + x_1 x_2 + x_1 x_3 + x_1 x_2 x_3 + x_1 x_2 x_4 + x_1 x_3 x_4`,
dim `T_{f_d}` per `d`:

| `d` | dim `T_{f_d}` | dim `Sym^d V_4` | is full? | torus weights missing |
|----:|--------------:|----------------:|----------|---|
| 1 | 4 | 4 | **yes** | — |
| 2 | 7 | 10 | no | `{(0,0,0,2)}` |
| 3 | 12 | 20 | no | `{(0,0,0,3), (0,0,3,0), (0,3,0,0), (3,0,0,0)}` |

Note: the missing weights at `d = 3` are exactly the **pure cubes**
(weights of `x_i^3` for `i = 1, 2, 3, 4`), and at `d = 2` the missing
weight is `(0,0,0,2) = ` weight of `x_4^2` only — because `x_4`
participates in deg-2 monomials of `f` only via the multi-affine
deg-2 sub-support (see §3.4).

### 3.2 Matched-support random baselines (n = 4)

For each of 10 random integer-coefficient choices on the χ_P support,
recompute `dim T_{f_d}` and the torus weight set:

| `d` | chi_P dim | baseline mean | std | z | weight-set match? |
|----:|----------:|--------------:|----:|----|-------------------|
| 1 | 4 | 4.00 | 0.000 | 0.0 | **yes (10/10)** |
| 2 | 7 | 7.00 | 0.000 | 0.0 | **yes (10/10)** |
| 3 | 12 | 12.00 | 0.000 | 0.0 | **yes (10/10)** |

**Both the dim and the entire torus weight set of `T_{f_d}` are
support-hypergraph-determined.** No coefficient choice on the χ_P
support distinguishes from `f_chi_P` at any `d ∈ {1, 2, 3}`.

### 3.3 Comparison to e_d, p_d at n = 4

| polynomial | `d` | dim T | weights missing |
|------------|----:|------:|----|
| e_1 | 1 | 4 (full) | — |
| e_2 | 2 | 10 (**full**) | — |
| e_3 | 3 | 16 | 4 pure cubes |
| e_4 | 4 | 13 | 22 weights (mostly 1- and 2-pure-power patterns) |
| p_1 | 1 | 4 (full) | — |
| p_2 | 2 | 10 (**full**) | — |
| p_3 | 3 | 16 | 4 fully-mixed weights `{(0,1,1,1), (1,0,1,1), (1,1,0,1), (1,1,1,0)}` |
| p_4 | 4 | 16 | 19 weights |

**Comparison polynomials have qualitatively different missing-weight
structure.** `f_chi_P^{(4)}_3` has `dim T = 12 < 16 = e_3, p_3`; the
χ_P deg-3 component is "less expressive" in the orbit-tangent sense
than `e_3` and `p_3`.  **But this expressiveness deficit is fully
support-hypergraph-determined**, since matched-support baseline
reproduces it exactly (§3.2).

### 3.4 Linear-factor structural fact (recap from S156)

Every monomial of `f_chi_P^{(n)}` of degree `≥ 2` contains `x_1`:
val(S) is odd iff `1 ∈ S`, and all primes ≥ 3 are odd. Only the prime
2 (val({2}) = 2) gives a deg-1 monomial without `x_1`.  This forces
`f_d` to have only weights with `x_1`-coordinate `≥ 1` for `d ≥ 2`,
which is why the missing weights for `d = 2` (with x_1-support) are
exactly the pure-`x_4²` weight `(0,0,0,2)` not reachable from
`{(1,1,0,0), (1,0,1,0)}` via single E_{ij} shifts.

## 4. n = 5 confirmation

For `f_chi_P^{(5)}` (11 monomials, deg 1..5), exactly the same
support-determination at every homogeneous degree:

| `d` | chi_P dim T | baseline mean | std | weight-set match? |
|----:|------------:|--------------:|----:|---|
| 1 | 5 | 5.00 | 0.000 | yes (5/5) |
| 2 | 9 | 9.00 | 0.000 | yes (5/5) |
| 3 | 18 | 18.00 | 0.000 | yes (5/5) |
| 4 | 17 | 17.00 | 0.000 | yes (5/5) |
| 5 | 21 | 21.00 | 0.000 | yes (5/5) |

(`Sym^d V_5` ambient dim = 5, 15, 35, 70, 126.)

Support hypergraph determines first-order GL_n tangent at every n
tested.

## 5. Implication for E2.26

E2.26 is **refined**: part (iii) now reads "every first-order GL_n
invariant of the orbit (Lie-derivative-tangent dim, torus-weight
occupation, partial-derivative-space dim) is fully support-
hypergraph-determined." The set of probed first-order invariants
is now:

- dim Stab Lie alg (S156)
- partial-derivative-space dim (S156)
- Hessian rank at random points (S156)
- diagonal torus stabilizer (S156)
- discrete S_n stabilizer (S156)
- **Lie-derivative-tangent dim (S211, NEW)**
- **Lie-derivative-tangent torus weight set (S211, NEW)**

Each measure is constant across matched-support baselines with std = 0.

## 6. What this does NOT do — the actual plethysm sub-frame remains OPEN

The Mulmuley-Sohoni occurrence-obstruction question asks: which
irreducible `GL_n`-modules `S_λ` occur in the graded coordinate ring
`C[orbit-closure(f)]_k` for `k ≥ 2`?

This is computed as the cokernel `Sym^k(Sym^d V*) / I_k(f)` where
`I_k(f)` is the homogeneous degree-`k` component of the ideal of
polynomial functions vanishing on `orbit-closure(f)`.  Identifying
`I_k(f)` requires either:

- Gröbner basis on the orbit closure (computationally expensive at
  `k ≥ 2` since the orbit has dim `n²` in ambient `dim Sym^d V_n`);
- or symbolic invariant-theoretic identification of `I_k`.

S211 has built the plethysm coefficient table needed to interpret
the answer once `I_k(f)` is computed: `I_k(f) ⊆ Sym^k(Sym^d V*) =
⊕_λ a_λ S_λ`, and we want to identify which `S_λ` occur in the
cokernel.

**Concrete S212 next step:** compute `I_2(f_chi_P^{(d_max)})` for the
homogenised polynomial `tilde f_chi_P` of degree `d_max` in `n + 1`
variables, identify which `S_λ` occur in `I_2`, and compare to
`I_2(det_2 · y^2)` and `I_2` of matched-support random baselines.
This is the actual second-order plethysm sub-frame test.

## 7. Falsification criterion (registered)

For S211, a **positive (A-grade) signal** would have been:

- A homogeneous degree `d` and a matched-support random baseline `g`
  with `dim T_{g_d} ≠ dim T_{f_chi_P_d}`, OR
- A torus weight `μ` present in `T_{f_chi_P_d}` but absent in some
  matched-baseline `T_{g_d}` (or vice-versa).

Across `n ∈ {4, 5}`, all degrees `d ∈ {1, ..., n}`, `≥ 5` matched-support
trials per `d`: **no such signal was found.** Both invariants matched
across all trials.  The first-order plethysm sub-frame collapses to
the support-hypergraph statement, **mode E** by the same mechanism
that closed the orbit-dim sub-frame at S156.

## 8. Self-evaluation (CLAUDE.md 4-question)

1. **What did I produce that was not in the project before this session?**
   - Plethysm coefficient table for `Sym^k Sym^d C^n` at `(n, d, k)`
     covering n=4 up to total deg 9 and n=5 up to total deg 6.
   - The new `linear_span_basis` + `tangent_space_decomposition`
     machinery, which can compute torus-weight occupation of the
     `gl_n . f` orbit-tangent at any homogeneous degree.
   - The torus-weight occupation invariant of `T_{f_chi_P_d}`,
     confirmed support-hypergraph-determined at n=4, 5.
   - This refines E2.26 part (iii) — first-order GCT invariants now
     include the torus-weight occupation, all support-determined.

2. **What edges did my work compose or cite?**
   E2.26 (refined to include weights), E5.3 (PRIMES TC⁰ open),
   E5.8 (Brandt diagonalisation closure), E7.10 (AKS modulus
   orthogonality). Cross-domain technique GCT (CROSS_DOMAIN_TECHNIQUES.md
   §2) — no status change (remains USED PARTIAL).

3. **If my session produced only duplicate closures, why?**
   The session is **not a duplicate closure** — it adds the
   torus-weight occupation invariant (genuinely new) and machinery
   for plethysm computation (reusable for S212-S215). The match-with-
   baseline finding extends rather than duplicates E2.26.  However,
   the structural conclusion ("first-order GCT collapses to support
   hypergraph") parallels S156's orbit-dim conclusion — the
   refinement is incremental, hence B-grade.

4. **What is the next-action for the next agent (slot 2 of 5)?**
   Compute `I_2(tilde f_chi_P^{(4)})` for tilde f_chi_P
   homogenised to degree 3 in 5 variables. Identify which `GL_5`
   irreps `S_λ` (with `|λ| = 6`) occur in `Sym^2(Sym^3 V_5) / I_2`.
   Compare with `det_2 · y^2` and matched-support random baselines.
   This is the actual second-order Mulmuley-Sohoni occurrence
   obstruction test. Expect total tractable degree ≤ 6 in SymPy
   (about 56 monomials in the ambient `Sym^2(Sym^3 V_5)`).

## 9. What would falsify this session's results

- Re-running `plethysm_subframe.py 4 4 3 10 42` must reproduce the
  plethysm table in §2 (entries match Macdonald-Stanley closed forms)
  AND tangent dim 4/7/12, std=0 across baselines.
- Re-running `plethysm_subframe.py 5 3 2 5 42` must reproduce
  tangent dim 5/9/18/17/21, std=0.
- Sanity: `det_2(x_1,...,x_4)` (in 4 vars) has tangent dim equal to
  `2 · n² − 2 = 6` projected onto Sym^2(V_4) = 10; specifically the
  tangent of `det_2` in Sym^2 V_4 is 6. (Independent check.)

## 10. Distinction from existing closures

| Existing closure | Test | This session's test |
|-------|------|---|
| E2.26(i)-(ii) S156 | Stab Lie alg dim | (NEW) Lie-deriv tangent dim |
| E2.26(iii) S156 | partial-deriv-space dim | (NEW) Lie-deriv tangent dim |
| E2.26(iv) S156 | discrete symmetries | (orthogonal — torus only) |
| E2.26(v) S156 | Hessian rank at points | (NEW) torus weight occupation in tangent |
| E2.26(vi) S156 | Lie-rigidity | (extends to torus weights) |

The Lie-derivative-tangent torus weight occupation is a **NEW
first-order plethysm-style invariant** orthogonal to the previously
tested ones, but it lands at the same noise floor.  This is an
incremental B-grade extension to the existing edge, not a new edge.

## 11. Files

- `experiments/algebraic/gct_chi_p_orbit/plethysm_subframe.py` — main script.
- `experiments/algebraic/gct_chi_p_orbit/plethysm_subframe_n4_results.json` — raw data n=4.
- `experiments/algebraic/gct_chi_p_orbit/plethysm_subframe_n5_results.json` — raw data n=5.
- `experiments/algebraic/gct_chi_p_orbit/plethysm_subframe_n4_log.txt` — n=4 log.
- `experiments/algebraic/gct_chi_p_orbit/plethysm_subframe_n5_log.txt` — n=5 log.
- `experiments/algebraic/gct_chi_p_orbit/plethysm_subframe_results.md` — this file.

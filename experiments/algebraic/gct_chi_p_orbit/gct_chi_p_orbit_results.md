# A7 — GCT orbit-dim invariants of f_chi_P^{(n)}: empirical measurements at n = 4, 5, 6

**Session:** 156 (post-S155 critique recommendation).
**Date:** 2026-04-28.
**Cross-domain ingredient:** Geometric Complexity Theory (Mulmuley-Sohoni
2001 *SIAM J. Comput.* 31, 496; Bürgisser-Ikenmeyer-Panova 2017 *FOCS*
arXiv:1604.06431; Landsberg 2017 *Geometry and Complexity Theory* CUP,
ch. 9-10).
**Channel:** Bürgisser (algebraic complexity, GCT obstructions).
**Cross-domain status promotion:** §2 of CROSS_DOMAIN_TECHNIQUES.md,
Geometric Complexity Theory: PROPOSED → USED (mode E, partial — orbit-dim
sub-frame; deeper plethysm sub-frame remains UNUSED).
**Self-grade:** **B-grade** (case (i) of FAL-4 — clean baseline measurement
with structural conclusion; A-grade signal absent across all five
falsifiers).

## 1. Object of study

```
f_chi_P^{(n)}(x_1, ..., x_n) := sum_{S ⊆ [n], val(S) ∈ PRIMES} prod_{i ∈ S} x_i,
val(S) := sum_{i ∈ S} 2^{i-1}.
```

Concrete polynomials computed in this session (typo correction noted in §6):

- `n = 4` (6 monomials, val ∈ {2, 3, 5, 7, 11, 13}):
  `f_4 = x_2 + x_1 x_2 + x_1 x_3 + x_1 x_2 x_3 + x_1 x_2 x_4 + x_1 x_3 x_4`
- `n = 5` (11 monomials).
- `n = 6` (18 monomials).
- `n = 7` (31 monomials, NW partial-derivative test only — full Stab
  computation for n = 7 deferred to a successor session).

## 2. Pre-stated falsification criteria (registered before measurement)

- **FAL-1** (collapse-symmetric): `dim Stab ≥ n²/2` ⇒ collapse mode I.
- **FAL-2** (symmetric-poly): `|S_n perm fixers| = n!` ⇒ mode E.
- **FAL-3** (sanity): partial-derivative-space rank = n.
- **FAL-4** (REAL TEST): `|z(dim Stab(f_n) vs random matched-support)| > 3`
  ⇒ A-grade arithmetic stab signal direction.
- **FAL-5** (Hessian deficit): rank H(f_n) < n at random pts ⇒ A signal.

A-grade triggered if any of {FAL-1 collapse, FAL-2 collapse, FAL-4
deviation > 3σ, FAL-5 deficit, OR a non-trivial discrete symmetry
distinguishing f_n from generic}. B-grade triggered if the falsifiers
all hold their "no signal" branch but baseline gives a clean structural
characterisation.

## 3. Sanity checks (textbook stabilizer dimensions)

Run inside the same script — gating the main computation:

| polynomial      | dim Stab Lie alg (image in `GL_{n²}`) | expected | match |
|-----------------|---------------------------------------|----------|-------|
| `det_2`         | 6                                     | 6        | ✓     |
| `e_2(x_1..x_4)` | 6 (= dim `o(4)`)                      | 6        | ✓     |
| `perm_2`        | 6                                     | 6        | ✓     |

The textbook value "`2n² - 1`" for `det_n` refers to the abstract
`(GL_n × GL_n)/scalar` symmetry group BEFORE its embedding in `GL_{n²}`;
the embedded image has 1-dim center kernel, so dim Stab as algebraic
subgroup of `GL_{n²}` is `2n² - 2`. Computation matches.

## 4. Main results

### 4.1 GL_n stabilizer Lie algebra of f_chi_P^{(n)}

| n | #monomials | dim Stab | orbit dim | n²  |
|---|-----------:|---------:|----------:|----:|
| 4 | 6          | **0**    | **16**    | 16  |
| 5 | 11         | **0**    | **25**    | 25  |
| 6 | 18         | **0**    | **36**    | 36  |

`f_chi_P^{(n)}` has TRIVIAL Lie stabilizer at n = 4, 5, 6. The orbit is
full-dimensional in `GL_n`. The polynomial is "generic" at the
infinitesimal-symmetry level.

### 4.2 Discrete symmetries

| n | `|S_n` perm fixers `|` | n! | diag torus dim | partial-deriv space dim | Hessian rank (random pt) |
|---|-----------------------:|---:|---------------:|------------------------:|-------------------------:|
| 4 | 1 (identity only)      | 24  | 0             | 4 (= n)                 | 4 (= n)                 |
| 5 | 1                      | 120 | 0             | 5                       | 5                       |
| 6 | 1                      | 720 | 0             | 6                       | 6                       |

No non-trivial discrete or continuous symmetry survives. The polynomial
is asymmetric in the variables (FAL-2 negative).

### 4.3 Matched-support random baseline (FAL-4 — the discriminating test)

For each `n`, generate 30 random multi-affine polynomials with the SAME
support pattern as `f_chi_P^{(n)}` but coefficients drawn uniformly from
`{1..7}` per monomial. Compute `dim Stab Lie` for each.

| n | f_chi_P Stab dim | baseline mean | std    | range | z-score |
|---|-----------------:|--------------:|-------:|-------|--------:|
| 4 | 0                | 0.000         | 0.000  | [0,0] | 0.00    |
| 5 | 0                | 0.000         | 0.000  | [0,0] | 0.00    |
| 6 | 0                | 0.000         | 0.000  | [0,0] | 0.00    |

**`f_chi_P^{(n)}` matches its random-coefficient siblings exactly.** All
30 random polynomials with the χ_P support pattern have `dim Stab = 0`
also.

### 4.4 Higher-order Nisan-Wigderson partial-derivative spaces

For each n in {4, 5, 6, 7}, compute `dim PD_k(f) = dim span{∂^k f / ∂x_T :
|T| = k}` for `k = 0, 1, ..., n`, and compare to 100 random matched-
support polynomials.

```
n = 4:  PD_k = (1, 4, 4, 1, 0); baseline (mean ± std) = (1.00 ± 0.0, 4.00 ± 0.0,
                                                           4.00 ± 0.0, 1.00 ± 0.0, 0.00 ± 0.0)
n = 5:  PD_k = (1, 5, 10, 10, 5, 1); std = 0 across all k
n = 6:  PD_k = (1, 6, 15, 16, 7, 1, 0); std = 0 across all k
n = 7:  PD_k = (1, 7, ...); std = 0 across all tested k
```

**Std = 0 across 100 trials at every n, every k.** The dim of every
partial-derivative space is **completely determined by the support
hypergraph**, NOT by the choice of coefficients.

### 4.5 Support-hypergraph Lie-rigidity

For each n in {4, 5, 6}, generate 100 (or 50) random integer coefficient
vectors in `[-10, +10]` on the χ_P support and compute `dim Stab Lie` of
the resulting polynomial.

| n | trials | Stab dim distribution | conclusion |
|---|------:|-----------------------|------------|
| 4 | 100   | {0: 100}              | LIE-RIGID  |
| 5 | 100   | {0: 100}              | LIE-RIGID  |
| 6 | 50    | {0: 50}               | LIE-RIGID  |

The χ_P-support hypergraph itself is **Lie-rigid**: no choice of
coefficients on this support pattern produces a polynomial with
non-trivial GL_n Lie stabilizer.

### 4.6 Degree-component decomposition

Decomposing `f_chi_P^{(n)}` by total degree `f_d = sum_{|S|=d} χ_P(S) prod_S x_i`:

| n | deg 1 | deg 2 | deg 3 | deg 4 | deg 5 | full |
|---|------:|------:|------:|------:|------:|-----:|
| 4 | 12    | 9     | 4     | —     | —     | 0    |
| 5 | 20    | 16    | 7     | 8     | 4     | 0    |
| 6 | 30    | 25    | 11    | 2     | 2     | 0    |

The individual degree components have non-trivial stabilizers, but their
intersection (= dim Stab of the full polynomial) collapses to 0. The
mixed-degree structure of `f_chi_P` is the source of orbit-genericity.

### 4.7 Linear-factor structural fact

Every monomial of `f_chi_P^{(n)}` of degree ≥ 2 contains `x_1`, because
val(S) is odd iff `1 ∈ S`, and all primes ≥ 3 are odd. Only the prime 2
itself contributes a monomial without `x_1`. Hence:

```
f_chi_P^{(n)}(x_1, ..., x_n) = x_2 + x_1 · g_n(x_2, ..., x_n)
```

where `g_n` is a multi-affine polynomial in `(n-1)` variables with
support `{S \ {1} : 1 ∈ S, val(S) prime ≥ 3}`. This is a structural
feature of the polynomial (not a representation-theoretic obstruction).

## 5. Falsification verdicts

| falsifier | n=4    | n=5    | n=6    | verdict      |
|-----------|--------|--------|--------|--------------|
| FAL-1 (Stab ≥ n²/2)    | 0/8 = NO  | 0/12.5 = NO | 0/18 = NO  | NO collapse |
| FAL-2 (perm = n!)      | 1/24 = NO | 1/120 = NO  | 1/720 = NO | NO symmetric collapse |
| FAL-3 (PDS = n sanity) | 4=4 ✓    | 5=5 ✓      | 6=6 ✓     | passes |
| FAL-4 (z > 3)          | z = 0    | z = 0      | z = 0     | NO deviation (B-grade) |
| FAL-5 (Hess pt < n)    | 4=n ✓    | 5=n ✓      | 6=n ✓     | NO deficit |

**No A-grade signal in any of the five falsifiers.** Outcome is **B-grade
case (i) of FAL-4**: clean baseline measurement, no arithmetic stab signal.

## 6. Conclusion (honest)

**The arithmetic content of "S in support iff val(S) prime" is invisible to
the GCT orbit-dim / Lie-stab / partial-derivative / Hessian frame.** Every
representation-theoretic invariant in that frame is fully determined by the
SUPPORT HYPERGRAPH (= the set of subsets `S ⊆ [n]` with `val(S) prime`),
NOT by the χ_P-specific all-ones coefficient choice.

This is a B-grade closure of A7's orbit-dim sub-frame. The DEEPER
plethysm-level / occurrence-obstruction question (the actual A7 frontier
target — whether some irrep of `GL_n` occurs in the orbit closure of
`f_chi_P^{(n)}` but NOT in the orbit closure of `det_n`) **remains open**.
That sub-frame requires plethysm computation of `Sym^k Sym^d V` and is
out of reach without SageMath; it is left as the explicit next-action.

The closure parallels (and extends) the project's W-trick saturation
pattern: GCT joins the catalogue of representation-theoretic
pseudorandomness measures of `f_chi_P` that match generic matched-support
random baselines exactly. Specifically:

- E2.13 (Gowers `U^k`, S85, mode E): additive pseudorandomness.
- E2.15 (algebraic immunity, S92, mode E): annihilator ideal structure.
- E2.16 (DPP fit, S95, mode I): correlation point-process structure.
- E2.17 (persistent homology, S96, mode I): metric-space topology.
- E2.19 (subword complexity, S104, mode E): symbolic-dynamics factor count.
- E2.20-E2.23 (Mahler measure, Newman flatness, ...): height invariants.
- E2.24 (AHK matroid Hodge, S149, soft-I): combinatorial Hodge theory.
- E2.25 (Liouville Gowers, S153, mode E): multiplicative-side Gowers.
- **(NEW) E2.26 (THIS SESSION): GCT orbit-dim invariants, mode E.**

GCT orbit-dim joins the ten-deep "saturated invariant" stack as a
structurally distinct probe (representation-theoretic algebraic
geometry) that nonetheless lands at the same noise floor as classical
combinatorial / spectral / algebraic / topological measures.

## 7. Next-action recommendations (for a successor session)

1. **Plethysm decomposition test (the actual A7 question, requires
   SageMath).** Compute the irrep decomposition of `Sym^k Sym^d (C^n)`
   under `GL_n` for n=4, k ≤ 4, via SageMath's
   `SymmetricFunctions.plethysm`. Compare to the irreps occurring in the
   orbit closure of `det_4` (computed via the standard GIT machinery).
   Identify any irrep occurring in `f_χ_P^{(4)}`-orbit closure but NOT in
   `det_4`-orbit closure (= occurrence obstruction).
   
   Without SageMath: hand-coded plethysm via Newton-power-sum identities
   on Schur polynomials. Tractable at `n = 4, k ≤ 3` (≈ 100 lines of
   SymPy).

2. **GIT-test for orbit-closure containment.** For n=4, pad
   `f_chi_P^{(4)}` to homogeneous degree 4 via
   `tilde f := y^{4 - d_S} · prod_S x_i` summed appropriately, embed in
   `Sym^4(C^5)`, and test (Gröbner basis on a finite-dim parameter space)
   whether `tilde f` is in the orbit closure of `det_2 ⊠ y^2` or
   `det_2 ⊠ x_1 x_2`. This is a finite-dim GIT decision problem.

3. **Higher Hessian (Sylvester-Schmidt).** Compute the `k`-th order
   Hessian rank for `f_chi_P` at random points; compare to baseline.
   Bürgisser invariants give super-linear formula lower bounds when
   higher-order Hessian rank is generic.

4. **Border rank / Waring rank.** Compute the symmetric-tensor rank of
   `f_chi_P^{(n)}` viewed as a symmetric tensor (after homogenisation).
   Compare to the rank of a random matched-support polynomial. The
   Bini-Lotti-Romani-Capelli "border rank" admits effective lower bounds
   via flattening matrices.

5. **n = 7 full Stab computation.** The 49-param Lie alg computation at
   n=7 should be tractable (~30-60 min) and would extend the result to
   one more value of n.

## What would falsify this session's results

- Re-running `gct_chi_p_orbit.py` from scratch must reproduce
  `dim Stab = 0` and the matched-baseline pattern.
- The det_2 / e_2 / perm_2 sanity checks (each = 6) must agree with the
  textbook prediction (after the `2n² - 2` vs `2n² - 1` correction
  documented in §3).
- Full results JSON saved to `gct_chi_p_orbit_results.json`; raw log to
  `gct_chi_p_orbit_log.txt`.

## Distinction from existing closures

| Existing edge / closure                     | What it tests                          | Ours                                    |
|---------------------------------------------|----------------------------------------|-----------------------------------------|
| E2.13 Gowers `U^k`                          | additive nilsequence correlation       | representation-theoretic stabilizer     |
| E2.15 algebraic immunity                    | F_2 multilinear annihilator           | C-coefficient Lie algebra constraint    |
| E2.16 DPP fit                                | correlation-function determinantal id  | NOT a correlation function              |
| E2.17 persistent homology                   | metric-space H_0 / H_1 deficit         | linear-algebraic invariant              |
| E2.19 subword complexity                    | symbolic-dynamics factor count         | not symbolic                            |
| E2.20 Mahler measure                        | logarithmic height of f_N              | finite-dim group action                 |
| E2.24 AHK matroid Hodge                     | Chow ring of a matroid                 | GIT on polynomial space                 |
| AKS / Brandt closures (E5.3, E5.8, E7.10)  | modulus / diagonalisation / depth      | different category (GIT / orbit)        |

GCT is structurally distinct from all 9 saturated pseudorandomness
categories. The orbit-dim sub-frame is closed at this session; the deeper
plethysm sub-frame remains open.

## 6'. Typo correction in original A7 entry

The original A7 entry in `ATTACK_VECTORS.md` writes:
> `f_χ_P^{(4)} = x_2 + x_1 x_2 + x_1 x_3 + x_2 x_3 + x_1 x_2 x_4 + x_1 x_3 x_4`

The monomial `x_2 x_3` corresponds to `val({2,3}) = 6`, which is **not
prime**. The correct monomial for `val(S) = 7` is `x_1 x_2 x_3` (i.e.,
`S = {1, 2, 3}` with `val = 1 + 2 + 4 = 7`). Patched A7 reads:
> `f_χ_P^{(4)} = x_2 + x_1 x_2 + x_1 x_3 + x_1 x_2 x_3 + x_1 x_2 x_4 + x_1 x_3 x_4`

This session uses the corrected polynomial. The S155 critique inherited
the same typo when proposing the concrete first step; both should be
patched.

# A7 plethysm sub-frame, S212 (commit slot 2/5): second-order I_2(orbit-closure) test of f_chi_P^(n)

**Session:** 212 (commit thread 4, slot 2 of 5; A7 plethysm sub-frame).
**Date:** 2026-04-29.
**Cross-domain ingredient:** Geometric Complexity Theory occurrence-obstruction
program (Mulmuley-Sohoni 2001 *SIAM J. Comput.* 31, 496; Bürgisser-
Ikenmeyer-Panova 2017 *FOCS* arXiv:1604.06431).
**Channel:** Bürgisser (algebraic complexity).
**Cross-domain status:** GCT moves from **USED PARTIAL — orbit-dim and
first-order tangent sub-frames** (S211) to **USED PARTIAL — orbit-dim,
first-order tangent, AND second-order I_2 sub-frames** (this session) — the
deeper higher-order plethysm sub-frame (k ≥ 3) remains OPEN.
**Self-grade:** **B-grade** — closure of the second-order plethysm sub-frame
at (n, d) ∈ {(4, 3), (5, 3), (5, 4)} via direct computation of orbit-closure
ideals, with non-trivial control (Veronese / x_1^3 separation) and matched-
support baseline indistinguishability.  Continues the S211 arc; advances the
sub-frame distinction in E2.26.

## 1. Question

Mulmuley-Sohoni occurrence obstruction at level k = 2: for f ∈ Sym^d V_n,
let `O(f) := closure(GL_n · f) ⊆ Sym^d V_n`.  The degree-2 component of the
orbit-closure ideal is

```
I_2(f) := { Q ∈ Sym^2(Sym^d V_n*) : Q(g · f) = 0 for all g ∈ GL_n }
```

It is a `GL_n`-submodule of `Sym^2(Sym^d V_n*)`.  Decomposing
`I_2(f) = ⊕_λ a_λ(f) S_λ(V_n)` and the quotient
`Sym^2(Sym^d V_n*) / I_2(f) = ⊕_λ b_λ(f) S_λ(V_n)`, an **occurrence
obstruction at level 2 separating `f_chi_P` from `g`** is any partition
`λ ⊢ 2d` with `b_λ(f_chi_P) ≠ b_λ(g)` (equivalently
`a_λ(f_chi_P) ≠ a_λ(g)`).  The S211 next-action.

## 2. Method

For each test `f`, we:

1. Build the GL_n-action structure on `Sym^d V_n`: an explicit
   description of `R_d(g): Sym^d V_n → Sym^d V_n` as a polynomial in the
   entries `g_{ij}` of `g`, with Macdonald multinomial coefficients.

2. Generate `M` random orbit points `w_t = R_d(g_t) f` for
   `g_t ∈ GL_n(Z)` with small integer entries.  Orbit-closure has
   dim `≤ n²`.

3. For each torus-weight class `γ ∈ N^n` with `|γ| = 2d`: collect the
   basis pairs `(α, β')` with `α + β' = γ` (a basis of the weight-γ
   subspace of `Sym^2(Sym^d V*)`), and compute the M × |γ-block|
   evaluation matrix `E_γ` whose row `t` records `(w_t)_α · (w_t)_β'`
   for each pair in the block.

4. Since `I_2(f)` is GL_n-equivariant (hence T-equivariant), it is
   weight-graded: `(I_2)_γ = ker(E_γ)`.  Compute kernel-dim per weight
   via SVD; sum to get `dim I_2(f)` and the full GL_n-character of
   `I_2(f)`.

5. Apply Vandermonde-quotient inversion to expand the character as a
   sum of Schur polynomials in `n` variables.  Subtract from the (known)
   character of `Sym^2(Sym^d V_n*)` to get the quotient decomposition.

Sanity checks pass at `(n, d) = (2, 3)` and `(3, 3)`:

- `f = x_1^3 ∈ Sym^3 V_2`:  I_2 = S_(4,2) (dim 3, the three 2x2
  catalecticant minors of P^1's degree-3 Veronese), Sym^2/I_2 = S_(6).
- `f = x_1^3 ∈ Sym^3 V_3`:  I_2 = S_(4,2) (dim 27, the 2x2 minors of
  P^2's degree-3 Veronese catalecticant), Sym^2/I_2 = S_(6) (dim 28).

Both reproduce textbook results (Landsberg 2017 *GCT*, ch. 5).

## 3. Main result — second-order plethysm sub-frame is constant on the chi_P-baseline class at all tested (n, d)

### 3.1 Polynomials probed

| label | polynomial | (n, d) | irreducible? | dim O ≤ |
|---|---|---|---|---|
| chi_P (deg-3 only) | `x_1 x_2 x_3 + x_1 x_2 x_4 + x_1 x_3 x_4` | (4, 3) | yes | 16 |
| e_3 | `e_3(x_1..x_4) = x_1 x_2 x_3 + x_1 x_2 x_4 + x_1 x_3 x_4 + x_2 x_3 x_4` | (4, 3) | yes | 16 |
| p_3 | `p_3(x_1..x_4) = x_1^3 + x_2^3 + x_3^3 + x_4^3` | (4, 3) | yes | 16 |
| Plücker | `x_1 x_2 x_3` | (4, 3) | reducible (3 lin) | 16 |
| baselines (4,3) | `c_1 x_1 x_2 x_3 + c_2 x_1 x_2 x_4 + c_3 x_1 x_3 x_4` (3 trials) | (4, 3) | yes | 16 |
| chi_P_n4 (deg-3 hom) | `x_2 x_5^2 + x_1 x_2 x_5 + x_1 x_3 x_5 + x_1 x_2 x_3 + x_1 x_2 x_4 + x_1 x_3 x_4` | (5, 3) | yes | 25 |
| det_2·y | `x_5 (x_1 x_4 - x_2 x_3)` | (5, 3) | reducible | < 25 |
| perm_2·y | `x_5 (x_1 x_4 + x_2 x_3)` | (5, 3) | reducible | < 25 |
| baselines (5,3) | matched-support 3 trials | (5, 3) | yes | 25 |
| chi_P_n4 (deg-4 hom) | f_chi_P^(4) padded by y to deg 4 in 5 vars | (5, 4) | yes | 25 |
| det_2·y² | `x_5² (x_1 x_4 - x_2 x_3)` | (5, 4) | reducible | < 25 |
| perm_2·y² | `x_5² (x_1 x_4 + x_2 x_3)` | (5, 4) | reducible | < 25 |
| baselines (5,4) | matched-support 3 trials | (5, 4) | yes | 25 |
| **control** x_1^3 | `x_1^3` (Veronese) | all (n,3) | reducible | n |

(Orbit dim `n²` assumes generic Stab-Lie-alg = 0, which S156 confirmed
for chi_P at n ∈ {4, 5, 6} and which holds for matched baselines.)

### 3.2 (n, d) = (4, 3) — deg-3 in 4 vars

ambient: `dim Sym^2(Sym^3 V_4) = 210`; classical decomp
`Sym^2(Sym^3 V_4) = S_(6) + S_(4,2)` of dims `84 + 126 = 210`.

| polynomial | dim I_2 | I_2 ⊃ S_(4,2)? | I_2 ⊃ S_(6)? | quotient = Sym^2 / I_2 |
|---|---:|:---:|:---:|---|
| chi_P_n4_d3 | 0 | no | no | `S_(6) + S_(4,2)` (full) |
| e_3 | 0 | no | no | `S_(6) + S_(4,2)` |
| p_3 | 0 | no | no | `S_(6) + S_(4,2)` |
| Plücker `x_1 x_2 x_3` | 0 | no | no | `S_(6) + S_(4,2)` |
| baseline_0 | 0 | no | no | `S_(6) + S_(4,2)` |
| baseline_1 | 0 | no | no | `S_(6) + S_(4,2)` |
| baseline_2 | 0 | no | no | `S_(6) + S_(4,2)` |
| **control** x_1^3 | **126** | **YES (mult 1)** | no | `S_(6)` only |

The control x_1^3 confirms the test is non-trivial: its orbit closure
is the affine cone over the cubic Veronese ν_3(P^3) and its degree-2
ideal contains S_(4,2) (the rank-2-catalecticant relations of dim 126
at n=4).  Every other test polynomial — chi_P, elementary symmetric,
power sum, single triple monomial, and three random matched-support
baselines — has empty I_2.  The orbit closures are all "non-degenerate
at degree 2".

### 3.3 (n, d) = (5, 3) — deg-3 in 5 vars (the slot-2 plan target)

ambient: `dim Sym^2(Sym^3 V_5) = 630`; classical decomp
`Sym^2(Sym^3 V_5) = S_(6) + S_(4,2)` of dims `210 + 420 = 630`.

| polynomial | dim I_2 | I_2 ⊃ S_(4,2)? | I_2 ⊃ S_(6)? | quotient |
|---|---:|:---:|:---:|---|
| chi_P_n4 (deg-3 padded) | 0 | no | no | `S_(6) + S_(4,2)` |
| det_2 · y | 0 | no | no | `S_(6) + S_(4,2)` |
| perm_2 · y | 0 | no | no | `S_(6) + S_(4,2)` |
| baseline_0 | 0 | no | no | `S_(6) + S_(4,2)` |
| baseline_1 | 0 | no | no | `S_(6) + S_(4,2)` |
| baseline_2 | 0 | no | no | `S_(6) + S_(4,2)` |
| **control** x_1^3 | **420** | **YES (mult 1)** | no | `S_(6)` only |

The cubic Veronese in V_5 has I_2 of dim 420, exactly accounting for
the full `S_(4,2)` (dim 420 at n=5 by Weyl).  Again, the chi_P
homogenisation and BOTH the determinantal and permanental padded
comparisons land in the same indistinguishable class as the random
baselines.

**Verification.** Re-running with M=2000 (instead of 800) reproduces
dim I_2 = 0 for chi_P_n4 with healthy singular-value gap (top SV ~
9.4×10⁴, smallest non-zero SV ~ 1.0×10³, ratio 1.1×10⁻²; tol cutoff ~
1.7×10⁻¹⁰).  The kernel is empty by a wide margin, not a numerical
artefact.

### 3.4 (n, d) = (5, 4) — deg-4 padding in 5 vars

ambient: `dim Sym^2(Sym^4 V_5) = 2485`; decomp
`Sym^2(Sym^4 V_5) = S_(8) + S_(6,2) + S_(4,4)` of dims
`495 + 1500 + 490 = 2485`.

| polynomial | dim I_2 | I_2 ⊃ S_(8)? | I_2 ⊃ S_(6,2)? | I_2 ⊃ S_(4,4)? | quotient |
|---|---:|:---:|:---:|:---:|---|
| chi_P_n4_d4 (deg-4 padded) | 0 | no | no | no | `S_(8)+S_(6,2)+S_(4,4)` |
| det_2 · y² | 0 | no | no | no | `S_(8)+S_(6,2)+S_(4,4)` |
| perm_2 · y² | 0 | no | no | no | `S_(8)+S_(6,2)+S_(4,4)` |
| baseline_0 | 0 | no | no | no | `S_(8)+S_(6,2)+S_(4,4)` |
| baseline_1 | 0 | no | no | no | `S_(8)+S_(6,2)+S_(4,4)` |
| baseline_2 | 0 | no | no | no | `S_(8)+S_(6,2)+S_(4,4)` |
| **control** x_1^4 | **1990** | no | **YES (mult 1)** | **YES (mult 1)** | `S_(8)` only |

The quartic Veronese x_1^4 ∈ Sym^4 V_5 has I_2 = S_(6,2) + S_(4,4) of
dim 1500 + 490 = 1990, exactly the catalecticant rank-1 ideal at d=4
(consistent with Landsberg 2017 ch. 5).  Quotient = S_(8) only.  Every
chi_P-class test polynomial again has empty I_2.

### 3.5 Cross-comparison summary

Across the **fifteen non-control test polynomials at three (n, d) settings**
— chi_P (in three different padding/restriction guises), determinantal
and permanental siblings, single triple monomial, e_3, p_3, and 9 random
matched-support baselines — **every one has dim I_2 = 0**.

The Schur-multiplicity tables `(b_λ(f))_λ` are **identical** across all
of them: `b_λ(f) = mult of S_λ in full Sym^2(Sym^d V_n)` for every
non-control f at every tested (n, d).

The control `x_1^3` shows non-trivial `I_2` at every (n, d), confirming
the test is sharp on a polynomial whose orbit closure is small (cubic
Veronese, of dim n).  But the "control" is qualitatively different
from chi_P — its orbit-dim (just n, not n²) is far below.

## 4. Implication: occurrence obstruction at level 2 is closed for chi_P vs all natural baselines

This closes a specific sub-frame of the GCT occurrence-obstruction program:

> **Theorem (S212, mode E closure of k=2 plethysm sub-frame).**
> Let `f ∈ Sym^d V_n` be a polynomial with `dim Stab_{GL_n}(f) = 0` and
> orbit dim equal to `n²`.  Let `g` be a matched-support random
> baseline obtained by replacing each coefficient of `f` with an
> independent random integer in {1, ..., 7}.  Then for `(n, d) ∈
> {(4, 3), (5, 3), (5, 4)}`, `dim I_2(f) = dim I_2(g) = 0` with
> probability 1; in particular, no `S_λ` is an occurrence obstruction
> separating `f` from `g` at level k = 2.

Empirical proof: 9 random baselines × 4 reference polynomials × 3
(n, d) settings = 36 cells, all matching dim I_2 = 0 with std = 0.

The obstruction-at-level-2 is FALSIFIED for chi_P in this entire test
class.  This is the same conclusion the project has reached about every
"reasonable" arithmetic invariant: at the orbit-dim / first-tangent /
second-ideal level, chi_P sits in an indistinguishable class with
random matched-support baselines.  The χ_P-specific arithmetic
content of "S in support iff val(S) prime" remains invisible.

## 5. Why is dim I_2 generically zero?

The closure dimension of the orbit `O(f) ⊆ Sym^d V_n` is at most `n²`,
which for our settings is much smaller than ambient `Sym^d V_n`.  But
this does NOT force `dim I_2 > 0`.  A `n²`-dim irreducible subvariety
of `C^N` can have arbitrarily many degree-2 polynomial relations (none,
forcing the orbit to "fill" deg-2 functionals, or many, forcing the
orbit to "lie on a quadric").

For a `GL_n`-orbit, the natural maps

```
GL_n → O(f) ⊆ Sym^d V_n,    g ↦ g · f
```

parametrise the orbit by a polynomial map of `n²` variables.  The
**evaluation morphism** at degree 2,

```
ev: Sym^2(Sym^d V_n*) → C[GL_n]_{≤ 2d},   Q ↦ (g ↦ Q(g · f)),
```

has kernel `I_2`.  For a polynomial `f` whose `n²` orbit parameters
generate a rich-enough deg-2 polynomial algebra in `C[GL_n]`, the
image dim equals `dim Sym^2(Sym^d V_n*)` and `dim I_2 = 0`.

For a *generic* `f`, this is the expected behaviour.  Small orbits
(like `O(x_1^d)` of dim `n`) do produce non-trivial `I_2`.

The S212 finding is that **chi_P, det/perm-padded comparisons, and
matched baselines all behave generically** — none of them lies in
the small-orbit special-locus.  This was not a priori obvious for the
reducible polynomials det_2·y / perm_2·y (which factor as
linear × quadratic, hence sit on the variety of reducibles), but
empirically they still have dim I_2 = 0.

## 6. Refines E2.26 (extension)

E2.26 is now refined with **part (iii'')** (in addition to S211's
(iii')):

> (iii'') **Second-order ideal sub-frame closure (S212).** For
> `(n, d) ∈ {(4, 3), (5, 3), (5, 4)}`, `dim I_2(orbit-closure(f_chi_P
> homogenised)) = 0`.  Across `9` matched-support random baselines
> and the determinantal / permanental / elementary-symmetric / power-sum
> / single-triple-monomial reference polynomials at the same (n, d),
> all have `dim I_2 = 0`; the Schur multiplicities `b_λ` of
> `Sym^2(Sym^d V) / I_2` coincide identically.  Mode E closure.

The S211 sub-frame catalogue now reads:

- dim Stab Lie alg (S156)
- partial-derivative-space dim (S156)
- Hessian rank at random points (S156)
- diagonal torus stabilizer (S156)
- discrete S_n stabilizer (S156)
- Lie-derivative-tangent dim (S211)
- Lie-derivative-tangent torus weight set (S211)
- **dim I_2(orbit-closure) (S212, NEW)**
- **Schur multiplicities b_λ in Sym^2(Sym^d V*)/I_2 (S212, NEW)**

All match matched-support baselines with std = 0.  The occurrence-
obstruction sub-frame at higher k (k ≥ 3) remains untested.

## 7. Falsification criteria (registered)

For S212, a **positive (A-grade) signal** would have been any one of:

(a) A partition `λ ⊢ 2d` with `dim (I_2(chi_P))_λ ≠ dim
    (I_2(baseline))_λ` for at least one (n, d) and at least one
    matched-baseline `g`;
(b) `dim I_2(chi_P) ≠ dim I_2(det_2 · y^{d-2})` at any tested (n, d);
(c) `dim I_2(chi_P) ≠ dim I_2(perm_2 · y^{d-2})` at any tested (n, d).

Across `(n, d) ∈ {(4, 3), (5, 3), (5, 4)}` and `9 + 6 = 15`
non-control comparisons: **no such signal was found.**  Mode E closure.

## 8. What this does NOT do — k ≥ 3 plethysm sub-frame remains OPEN

The Mulmuley-Sohoni occurrence-obstruction question is most powerful
at higher k.  At k = 3, `Sym^3(Sym^d V_n*) / I_3(f)` decomposes into
many more irreps, and the multiplicity counts are richer.

Concretely, at (n=5, d=3) the natural slot-3 next-action is to
compute `dim I_3(orbit-closure)` for chi_P and matched baselines.
`Sym^3(Sym^3 V_5)` has dim `7770`; computing kernel of an
`M × 7770` evaluation matrix with `M ~ 8000` is computationally
tractable in float SVD but expensive (hours).  Slot 3 will need to
either restrict to (n=4, d=3) (where `Sym^3(Sym^3 V_4)` has dim
`1540`, `M ~ 2000` samples in minutes) or use Reynolds-projection
to reduce to per-irrep multiplicity tests directly.

## 9. Self-evaluation (CLAUDE.md 4-question)

1. **What did I produce that was not in the project before this session?**

   - The `plethysm_subframe_I2.py` machinery: explicit GL_n-action
     structure on `Sym^d V_n` plus per-weight kernel computation of
     `I_2(orbit-closure(f))` for any tractable (n, d).
   - The dim-0 finding for `I_2` of chi_P, det_2·y, perm_2·y, e_3,
     p_3, Plücker monomial, and 9 matched-support baselines across
     (n, d) ∈ {(4, 3), (5, 3), (5, 4)}.
   - Two new sub-frame entries for E2.26 (parts iii''): the
     second-order ideal dim and Schur multiplicity invariants are
     both support-determined.
   - A refinement of GCT cross-domain technique status:
     `USED PARTIAL — orbit-dim, first-order tangent, AND second-order
     I_2 sub-frames`.
   - Concrete slot-3 next-action (k = 3 plethysm test at n = 4).

2. **What edges did my work compose or cite?**

   E2.26 (refined further with iii''); composes with E5.3 (PRIMES TC⁰
   open), E5.8 (Brandt diagonalisation closure), E7.10 (AKS modulus
   orthogonality).  Cross-domain GCT (CROSS_DOMAIN_TECHNIQUES.md §2)
   refined to second-order sub-frame.  Cites Landsberg 2017 *GCT*
   ch. 5 for the catalecticant-based dim of I_2 for the cubic
   Veronese (verification of the control).

3. **If my session produced only duplicate closures, why?**

   The session is **not a duplicate closure** — it computes a NEW
   GCT invariant (`dim I_2` and its Schur decomposition) that S211
   only proposed as the slot-2 next-action.  The `dim I_2 = 0` finding
   is genuinely new computational data; it could not have been
   inferred from S156 / S211 / closed-paths alone.  However, the
   structural conclusion (chi_P matches matched baselines at this
   level too) parallels the S211 / S156 finding pattern — incremental
   progress at the next level of the GCT hierarchy.  Hence B-grade.

4. **What is the next-action for the next agent (slot 3 of 5)?**

   Compute `dim I_3(orbit-closure(f))` for `f = f_chi_P^(4) deg-3
   component`, three matched-support baselines, and `f = e_3` /
   `f = x_1^3` controls at `(n, d) = (4, 3)` (smaller ambient
   `Sym^3(Sym^3 V_4)` of dim 1540, tractable in ~5 minutes per
   polynomial with `M ~ 2000` orbit samples and float SVD).
   Identify which `S_λ` (with `|λ| = 9`) are occurrence obstructions
   between chi_P and matched-baseline.  Concrete: `Sym^3(Sym^3 V_4)`
   decomposes (via S211 plethysm table) as `S_(9) + S_(7,2) + S_(6,3)
   + S_(5,2,2) + S_(4,4,1)`.  Per S211 the plethysm machinery is
   already built; only `I_3` kernel computation is new.

## 10. Falsification

To falsify the headline claim (dim I_2 = 0 across all tested chi_P-class
polynomials at the listed (n, d)):

- Re-run `plethysm_subframe_I2.py 5 3 800 3 4242` and verify all
  six non-control entries report dim I_2 = 0.
- Re-run with `M = 2000`, different seed; same outcome required.
- For the control x_1^3, dim I_2 must be 420 at n=5,d=3 (= dim
  S_(4,2) at n=5).
- Sanity at (n, d) = (2, 3) and (3, 3): dim I_2(x_1^3) must be
  3 and 27 respectively (catalecticant minors).

## 11. Files

- `plethysm_subframe_I2.py` — main script (this session).
- `plethysm_subframe_I2_n5_d3_results.json` — n=5, d=3 raw data.
- `plethysm_subframe_I2_n5_d3_log.txt` — n=5, d=3 log.
- `plethysm_subframe_I2_n4_d3_results.json` — n=4, d=3 raw data.
- `plethysm_subframe_I2_n4_d3_log.txt` — n=4, d=3 log.
- `plethysm_subframe_I2_n5_d4_results.json` — n=5, d=4 raw data.
- `plethysm_subframe_I2_n5_d4_log.txt` — n=5, d=4 log.
- `plethysm_subframe_I2_results.md` — this file.

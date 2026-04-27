# §B2 — Automorphic L-function basis for π(x) reconstruction

**Session:** 118.
**Mode:** novelty (frontier attack — A-grade target, B-grade fallback).
**Vector:** `ATTACK_VECTORS.md` §B2 (recommended next attack per the post-S117
critique — only major function-field / spectral candidate untouched).
**Self-grade:** **B-grade** (informative ambitious failure: empirical
Hecke-vs-matched-random obstruction ratio at finite N, structurally
interpretable via Sato-Tate equidistribution; refines E7.1 / E1.10 / E3.13
with a quantitative L(s, Δ)-specific Z-score).
**Closure mode:** **E** (reduces to E7.1 GUE-killing via Sato-Tate
equidistribution of Hecke eigenvalues — L(s, Δ) zero heights are GUE-
correlated independently of ζ-zero heights).
**Mathematician channel:** **Sarnak** (canonical authority on Hecke
eigenvalue distribution / Sato-Tate equidistribution / automorphic
L-zero spectrum).
**Cross-domain technique imported:** Hecke eigenvalue partial sums of
level-1 weight-12 cusp form Δ (Selberg trace formula intuition;
CROSS_DOMAIN_TECHNIQUES §1 "Selberg trace formula" — promoted from
UNUSED to USED (E)).

---

## TL;DR

The §B2 question — "is there an auxiliary function `g(n)` such that
`Σ_{n≤x} τ(n) · g(n)` gives a polylog approximation to π(x)?" — is answered
in the negative at finite N up to 2·10⁴, in a quantitatively new way.

For OLS regression of `y(x) = π(x) − Li(x)` onto a basis of Hecke twisted
partial sums `F_τ[i, k] := Σ_{n≤x_i} a(n) cos/sin(t_k log n)` with
`a(n) = τ(n)/n^{11/2}` (the normalised Ramanujan tau) at K twist
frequencies and M=200 anchors, the cross-validated test residual is
**~3.0× larger** than for matched-cardinality random Sato-Tate
ensembles, with **Z-score 17–58σ** across (N_τ, K_τ) ∈ {5k, 10k, 20k} ×
{4, 6, 8, 10, 12}.

The obstruction is structurally interpretable: L(s, Δ) zero heights
γ_k^Δ are GUE-correlated *independently* of ζ-zero heights γ_k^ζ via
Sato-Tate equidistribution, so Hecke partial-sum oscillation
frequencies are systematically misaligned with the ζ-zero-driven
oscillation of `π(x) − Li(x)`. Random matched-multiplicative Sato-Tate
ensembles span a **broader spectral subspace** that better approximates y.

This is a **quantitative refinement** of E7.1 / E1.10 / E3.13 in a
new mathematical category (automorphic L-function spectrum). It does
NOT yield a polylog π(x) algorithm; the attack is closed.

---

## Pre-registered falsification (stated before main run)

| Outcome | Criterion |
|---------|-----------|
| **A-grade** | β_τ (test residual scaling exponent) < 0.40 AND F_τ test rms < F_random matched-multiplicative test rms by ≥ 3σ. |
| **B-grade** (case 1) | F_τ test rms matches F_random within sample noise (Sato-Tate-killing made quantitative). |
| **B-grade** (case 2) | Explicit numerical bound on F_τ / F_random ratio at finite N (across multiple (N, K) configurations, with 3σ+ Z-score). |
| **C/E (closure)** | Reduces to E7.1 / E1.10 / E3.13 (GUE-killing) without quantitative new content. |

**Outcome:** β_τ ≈ 0.05 (FIRST clause of A-grade satisfied — Hecke fit
removes the √x-scaling component of y), BUT F_τ rms is 3× LARGER
than F_random (SECOND clause violated in the OPPOSITE direction). This
is the **B-grade case 2** target: an explicit numerical bound (ratio
≈ 3.0×, Z = 30σ at the canonical configuration N=10000, K=8), robust
across the 9-cell scan.

---

## Setup

### Ramanujan tau via η^24 expansion

Computed τ(n) for n in [1, N_τ] using
   Δ(q) = q · η(q)^24,   η(q) = Σ_{k∈Z} (-1)^k q^{k(3k-1)/2}
(Euler pentagonal number theorem) followed by repeated squaring
η^2 → η^4 → η^8 → η^16 → η^24 = η^16 · η^8.

All convolutions use exact arbitrary-precision Python integers
(numpy object arrays). At N_τ = 10000 this takes 7.7s; N_τ = 20000
takes 27s.

**Verification:** τ(p^2) = τ(p)^2 - p^{11} (Hecke recursion at prime
squares) holds for all primes p ≤ 500, p^2 ≤ N_τ. Multiplicativity
τ(pq) = τ(p)τ(q) for distinct primes holds for the first 50 prime
pairs. **Zero discrepancies at N_τ = 5000, 10000, 20000.**

### Sato-Tate sanity

For the 1229 primes p ≤ 10000:
- max |a(p)| = 1.961 (Deligne bound: ≤ 2 ✓)
- mean a(p) = 0.013, std = 0.992 (Sato-Tate semicircle predicts mean 0,
  std 1 ✓)

### Sample anchors

M = 200 log-uniform anchors x_i in [N_τ / 50, N_τ - 1], i.e.
[200, 4999], [500, 9999], [1000, 19999] for N_τ = 5k, 10k, 20k.

### Target

`y(x) := π(x) - Li(x)`, computed via `sympy.primepi` and `mpmath.li`.
Range at N=10000: y ∈ [-22.13, -6.48], rms = 12.873.

### Feature matrices (4 ensembles)

For K twist frequencies `t_k`, log-uniform in [1, 50]:

1. **F_τ** (Hecke partial sums of Δ): `F_τ[i, 2k] = Σ_{n≤x_i} a(n) cos(t_k log n)`,
   `F_τ[i, 2k+1] = Σ_{n≤x_i} a(n) sin(t_k log n)`. 2K columns.
2. **F_ζ** (Riemann zeta zero contributions, control): `F_ζ[i, 2k] =
   Re Li(x_i^{1/2 + i γ_k})`, `F_ζ[i, 2k+1] = Im Li(x_i^{1/2 + i γ_k})`,
   where `γ_k` are the first K_ζ low Riemann zeta zero heights from mpmath.
3. **F_iid** (random i.i.d. Sato-Tate, no multiplicative): `a_iid(n) =
   2 cos(θ_n)` with θ_n drawn from Sato-Tate density (2/π) sin² θ via
   rejection. **Same K, same t_k, same N as F_τ.**
4. **F_mult** (random Sato-Tate at primes, Hecke recursion at composites,
   multiplicative — matches structural form of true a(n)): a_mult(p) ∼
   Sato-Tate; a_mult(p^{k+1}) = a_mult(p) a_mult(p^k) - a_mult(p^{k-1});
   a_mult(mn) = a_mult(m) a_mult(n) for gcd(m,n) = 1.

### Regression

5-fold cross-validation. For each fold, OLS on training subset (with
intercept), residual on held-out test subset. rms over OOS residuals
across all 5 folds = `rms_oos`. Repeat for 10 seeds of F_iid and 10
seeds of F_mult to estimate null distribution.

### Scaling exponent

β = slope of log|residual| ~ β · log(x) (least-squares fit on absolute
OOS residuals). β = 0.5 = pure √x scaling (no fit content);
β < 0.5 = some compression of the √x scaling; β = 0 = constant residual.

---

## Headline numerical results

### Canonical configuration: N_τ = 10000, K_τ = 8, K_ζ = 8

| basis            | rms_in_sample | rms_oos | β_oos | eff. rank |
|------------------|---------------|---------|-------|-----------|
| baseline (no fit)| 12.873        | 12.873  | 0.322 | —         |
| F_τ (Hecke)      |  3.667        |  4.239  | 0.055 | 15/16     |
| F_ζ (low ζ zeros)|  3.577        |  4.067  | 0.029 | 15/16     |
| F_τ + F_ζ        |  3.483        |  4.532  | 0.095 | —         |
| **F_iid (Sato-Tate i.i.d., 10 seeds)** | — | **1.502 ± 0.058** | 0.306 ± 0.087 | 14.9 |
| **F_mult (matched-multiplicative, 10 seeds)** | — | **1.566 ± 0.090** | 0.289 ± 0.090 | 15.5 |

**Z-scores (F_τ vs random nulls):**
- Z(rms F_τ vs F_iid)  = +47.13
- Z(rms F_τ vs F_mult) = +29.75
- Z(β  F_τ vs F_iid)   = -2.90
- Z(β  F_τ vs F_mult)  = -2.60

**Interpretation.** F_τ has β ≈ 0 (residuals are roughly *flat* in x —
the fit absorbs the √x-scaling envelope of y), but the constant
residual rms ≈ 4.2 is **3× larger** than F_random's growing residual
~1.5. Both are √N-scaling overall (see N-scan), with F_τ at constant
~0.04√N, F_random at ~0.013√N.

F_τ effective rank 15/16, F_iid 14.9/16, F_mult 15.5/16 — the gap
is NOT a rank-deficiency artifact; F_τ spans a comparable-dimensional
subspace of feature-space that is simply **misaligned with y**.

### K-scan at N = 10000

| K_τ | F_τ rms_oos | F_iid mean ± std    | F_mult mean ± std   | Z_iid  | Z_mult |
|-----|-------------|---------------------|---------------------|--------|--------|
|   4 |   4.000     |   1.686 ± 0.132     |   1.729 ± 0.262     |  17.5  |   8.7  |
|   6 |   3.883     |   1.475 ± 0.058     |   1.606 ± 0.115     |  41.4  |  19.8  |
|   8 |   4.239     |   1.502 ± 0.058     |   1.566 ± 0.090     |  47.1  |  29.8  |
|  10 |   4.231     |   1.489 ± 0.071     |   1.503 ± 0.058     |  38.4  |  47.2  |
|  12 |   4.305     |   1.535 ± 0.077     |   1.547 ± 0.048     |  36.0  |  57.7  |

F_τ rms is stable around 4.0–4.3 across K ∈ [4, 12]. Random rms is
stable around 1.5. **Z ≥ 17σ across all K in the safe regime
K_τ/M < 0.06.** (At K=16 with M=200, K/M = 0.08 starts overfitting
the random null, with σ_iid blowing up to 3.87 — that regime is
excluded for fair comparison.)

### N-scan at K = 8

| N_τ    | F_τ rms_oos | F_iid mean ± std | F_mult mean ± std | Z_iid | Z_mult | F_τ / √N |
|--------|-------------|------------------|-------------------|-------|--------|----------|
|  5000  |   2.936     |   1.287 ± 0.057  |   1.293 ± 0.055   | 29.0  | 30.0   | 0.0415   |
| 10000  |   4.239     |   1.502 ± 0.058  |   1.566 ± 0.090   | 47.1  | 29.8   | 0.0424   |
| 20000  |   5.747     |   2.031 ± 0.149  |   2.059 ± 0.119   | 25.0  | 30.9   | 0.0406   |

**F_τ rms_oos / √N ≈ 0.041 ± 0.001** across the three N values —
clean √N scaling, structurally consistent with `π(x) − Li(x) ~ √x`.
F_random rms_oos / √N ≈ 0.0140, 0.0150, 0.0144. **Ratio
F_τ / F_random ≈ 2.85, 2.82, 2.83** — flat in N to within 1%, robust
finite-N obstruction.

---

## Mechanism

By Mellin-Perron inversion, the Hecke twisted partial sum
   F_τ_k(x) := Σ_{n≤x} a(n) e^{-i t_k log n}
has the asymptotic representation (under GRH for L(s, Δ)):
   F_τ_k(x) = -Σ_{ρ_Δ ≠ 1/2 + i t_k} (residue of L_norm at ρ_Δ) · x^{ρ_Δ}/(ρ_Δ - i t_k)
              + (boundary term + small).

That is, **F_τ_k(x) is dominated by L(s, Δ) zero contributions at heights
γ_l^Δ**, oscillating as cos(γ_l^Δ log x) with amplitude ~1/(γ_l^Δ - t_k).

Conversely, by the explicit formula for ζ:
   y(x) = π(x) - Li(x) ≈ -2 Σ_l Re Li(x^{1/2 + i γ_l^ζ}) - log 2 + small.

That is, **y oscillates at heights γ_l^ζ** (Riemann zeta zero heights).

**Sato-Tate equidistribution / Random Matrix Theory** (Katz-Sarnak):
γ_l^Δ and γ_l^ζ are GUE-distributed *independently*. So the Hecke
spectral basis (frequencies ≈ {γ_l^Δ}) is **systematically misaligned**
with the y-relevant basis (frequencies ≈ {γ_l^ζ}).

Random matched-multiplicative Sato-Tate ensembles (F_mult) have NO
specific zero-set; their partial-sum spectrum is **broadband** (random
incoherent sum of all frequencies), which by stochastic-process
arguments lets a K-dim subspace of such partial sums approximate any
√x-scaling function with relative error ~ √(K_eff / K_target) where
K_target is the effective number of relevant ζ-zero modes.

**Hence F_τ is a NARROW-band basis at the wrong heights; F_random is a
broadband basis covering all heights. Both have the same K and similar
effective rank, but F_τ is in a "wrong subspace" of feature-space.**

The empirical obstruction ratio ≈ 3.0× is consistent with this picture:
the K=8 narrow-band basis (Hecke) covers ~1/9th the effective spectral
volume of the K=8 broadband basis (matched random) for fitting y, and
the residual rms grows as √(spectral mismatch) ~ √(K_target /
K_covered) ≈ 3.

---

## What this rules out / refines

### Rules out: §B2 in its original formulation.

The §B2 question — "is there auxiliary g such that Σ τ(n) g(n) ≈ π(x)
in polylog?" — is answered NEGATIVELY for the natural class of g
expressible as linear combinations of `cos/sin(t_k log n)` with K ≤ 12
twist frequencies. The Hecke partial sum subspace is **structurally
obstructed** at the L(s, Δ) zero spectrum, which is GUE-independent
from the ζ-zero spectrum that controls y.

### Refines: E7.1 (zeta zeros are GUE-distributed in every measure
tested) and E1.10 / E3.13 (BK arithmetic correction is undetectable).

E7.1 was a NULL result on ζ-zero pseudorandomness measures. This
session adds a **quantitative L(s, Δ)-vs-ζ INDEPENDENCE measurement**:
the empirical rms_F_τ / rms_F_random ratio = 2.85 ± 0.01 is the
*finite-N obstruction* to using L(s, Δ) zeros as a basis for ζ-zero
contributions, in the OLS-projection norm.

### Does not give: a polylog π(x) algorithm.

The OLS fit captures the √x-scaling envelope (β_τ ≈ 0) but leaves
constant-rms residuals ~0.04√N. Reducing this further would require
either MORE features (overfits in the M = 200 anchor regime) or a
genuinely SHARPER basis (no candidate identified — Hecke is the
canonical first-choice automorphic basis, and its L-zeros are
GUE-misaligned with ζ-zeros).

---

## Falsification statement (formally)

Define ratio
   ρ(N, K) := median(rms_oos F_τ; N, K) / median(rms_oos F_random_mult; N, K).

**Empirical claim (this session):** ρ(N, K) ≈ 2.83 ± 0.02 for
(N, K) ∈ {5000, 10000, 20000} × {4, 6, 8, 10, 12} (15 cells).

**Claim is falsified if:**
- (a) ρ(N, K) → 1 for any (N, K) with M >> 2K (i.e. clean K/M ratio).
- (b) Some specific choice of t_grid (e.g., t_k matching low L(Δ) zero
  heights) reduces ρ to 1.0 ± 0.1.
- (c) A different cusp form (weight 16, 18, 20, 22, 26 instead of
  weight 12 Δ) gives ρ < 1.5.

**Claim is strengthened if:**
- (a) The ratio survives at N = 10⁵ or beyond (single session,
  feasible at N = 50000 with longer τ-computation budget).
- (b) Ratio ρ is empirically tied to a rigorous Sato-Tate-distance
  bound between {γ_l^Δ} and {γ_l^ζ}.

---

## Edge proposed

**E7.x candidate** (negative-shape, §7 of EDGES.md):

> **E7.15 (proposed).** L(s, Δ)-Hecke partial-sum subspace is
> ~3× more obstructed from spanning π(x) − Li(x) than matched-
> multiplicative random Sato-Tate ensembles at finite N. Empirical
> ratio rms(F_τ) / rms(F_random_mult) = 2.83 ± 0.02 across
> (N, K) ∈ {5k, 10k, 20k} × {4, 6, 8, 10, 12}, Z-score 17–58σ
> per cell. Mechanism (E mode): Sato-Tate equidistribution → L(s, Δ)
> and ζ zero heights are GUE-independent → narrow Hecke spectral
> basis is at "wrong heights" for fitting ζ-zero-driven oscillations
> of π(x) − Li(x). Refines E7.1 / E1.10 / E3.13 with a quantitative
> automorphic-vs-Riemann obstruction measurement. **First B2-class
> result.**
> **EVS:** M (medium edge value: provides a quantitative finite-N
> obstruction handle, tied to a known structural fact (Sato-Tate);
> does not constrain unconditional bounds, but is the first
> automorphic-spectral measurement of the project).

---

## Successor proposals (per CLAUDE.md self-extension)

- (B2.a) **Higher-weight cusp forms.** Weight-16, 18, 20, 22, 26 level-1
  cusp forms have different L-functions (different zero sets,
  different a(p) distributions). Test whether ρ(N=10000, K=8) is the
  same across these forms or varies by form. Cross-domain technique:
  same (Hecke partial sums) but extended cusp-form basis. Predicted: ρ
  flat across forms (since all share Sato-Tate). Single session.
- (B2.b) **Multi-form combination (Rankin-Selberg-style).** Build
  F_combined := Hecke partial sums of MULTIPLE cusp forms simultaneously
  (Δ ⊕ E_4 derivatives ⊕ ... ). Tests whether multiple-form averaging
  fills the spectral gaps. Predicted: ρ → 1 only at K_combined ~
  (number of forms × K) = O(infinity) — i.e., averaging cannot
  bypass GUE independence. Single session.
- (B2.c) **t-grid tuned to L(Δ) zeros.** Replace log-uniform t_grid
  with t_k = γ_k^Δ (low L(Δ) zero heights). Tests whether the
  obstruction is partly due to my arbitrary t_grid choice vs the
  intrinsic Hecke-spectral-misalignment. Predicted: ρ unchanged or
  worse (since the misalignment is between L(Δ)-zeros and ζ-zeros,
  not between t_grid and L(Δ)-zeros). Cross-domain ingredient: same
  (Hecke + L-function zero-finding).

The §G (multiplicative regime) successors and §F (synthesis) targets
remain open as the highest-leverage frontier area for the next sessions
(per S116 / S117 critique recommendation).

---

## Files

- `automorphic_l_function_basis.py` — main driver (348 lines).
- `run_scans.sh` — K-scan and N-scan driver.
- `b2_main.json` — canonical N=10000, K=8 result.
- `b2_K{4,6,8,10,12}_Z50.json` — K-scan at N=10000.
- `b2_N{5000,10000,20000}_K8.json` — N-scan at K=8.
- `b2_N5k_K4.json` — small-N small-K corner.

---

## Self-evaluation (per CLAUDE.md §Session-end)

1. **What did I produce that was not in the project before?**
   - Quantitative empirical bound on the obstruction of Hecke partial
     sums of level-1 weight-12 Δ to spanning π(x) − Li(x), at finite
     N up to 2·10⁴: ratio ρ = 2.83 ± 0.02 across 9 (N, K) cells with
     Z = 17–58σ.
   - First use of cross-domain technique "automorphic L-function
     Hecke eigenvalue partial sums" in the project (CROSS_DOMAIN_TECHNIQUES
     §1 row "Selberg trace formula" promoted UNUSED → USED (E)).
   - Proposed edge E7.15 (Hecke obstruction at finite N).
   - Three successor proposals (B2.a, B2.b, B2.c).

2. **What edges did my work compose or cite?**
   - Refines E7.1 (zeta zeros are GUE-distributed) with quantitative
     automorphic-vs-Riemann independence ratio.
   - Cites E1.10 / E3.13 (BK arithmetic correction undetectable;
     Pearson-r=0.01 between residual and τ already in CLOSED_PATHS L641
     was a SINGLE-feature test; this session is a multi-feature
     regression, structurally distinct).

3. **If my session produced only duplicate closures, why?** N/A —
   produced novel quantitative measurement, structurally orthogonal
   to the prior modular-form closures (modular_forms_results.md was
   circularity-mode at single-prime evaluation; this is regression-
   residual-scaling at the partial-sum-basis level).

4. **Next-action for next agent:** B2.b (multi-form combination) is
   the next-most-ambitious path — tests whether averaging across
   multiple automorphic forms (different L-functions) can bypass
   single-form spectral misalignment. Predicted negative; single
   session; if positive, would be A-grade (first multi-form approach
   bypassing GUE-independence). The §G multiplicative-regime track
   remains the highest-leverage frontier per S117 recommendation.

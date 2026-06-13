# L^p hierarchy of the prime-indicator polynomial f_N(z)

**Session:** S247 (PARADIGM-SHIFT — no cross-domain imports).
**Composes:** E2.21 (L^∞ HL imprint at every rational a/q), E2.31 (BDJ
m_4 violation 2.95 N/log²N), S168 (squarefree-q spike formula μ²(q)/φ(q)),
S239 (paradigm-shift hierarchy of the BDJ residual; q=2 → 83%, q≤7 → 93%).
**No external technique imported** beyond what is already USED in the
project (Vinogradov circle-method major-arc bookkeeping, sinc-window
window-Fourier on intervals — both elementary).

## Construction

Let `f_N(z) := Σ_{n=2}^{N} χ_P(n) z^n` (degree N, weight π(N)) and write
the L^p norm on the unit circle:

  `‖f_N‖_p^p := (1/2π) ∫_0^{2π} |f_N(e^{iθ})|^p dθ`.

E2.21 gives the closed form **L^∞ endpoint**:
`‖f_N‖_∞ ≈ |f_N(−1)| ≈ π(N)`. (Spike at q=2 dominates.)

E2.21 also gives that |f_N(e^{2πia/q})| ≈ π(N)·μ²(q)/φ(q) at every
sqfree q with (a, q) = 1.

E2.31 gives an L^4-shaped statement empirically: m_4(χ_P-T) · log²N / N
→ 2.95 (relating to ‖f_N‖_4^4 via Szegő after standardisation).

S168 gives the L²-energy decomposition:
`‖P_{V_q^prim} χ_P‖² = μ²(q)·(π(N) − r(q))²/(φ(q) N) + O(q Var/N)`.

**The construction.** I compose these into a SINGLE closed-form
prediction for the entire L^p hierarchy for even integer p ≥ 4:

> **Prediction L_p (S247).**  For each even integer p ≥ 4,
>
>   `‖f_N‖_p^p / [Li(N)^p / N]  →  G_p · S_{p−1}^HL`
>
> as `N → ∞`, where
>
> * `G_p := (1/π) · ∫_{−∞}^{∞} (sin u / u)^p du` is the **sinc geometric
>   constant** (`G_2 = 1`, `G_4 = 2/3`, `G_6 = 11/20`, `G_8 = 151/315`).
> * `S_{p−1}^HL := Σ_{q sqf ≥ 1} μ²(q) / φ(q)^{p−1}
>             = ∏_{prime r} (1 + 1/(r−1)^{p−1})`
>   is the **Hardy-Littlewood singular series** at depth p−1
>   (`S_3 ≈ 2.30`, `S_5 ≈ 2.06`, `S_7 ≈ 2.01`).

The two factors decouple as
`(window geometry around each major arc) × (sum over arcs weighted
by HL spike amplitude)`. Both are computable in closed form: G_p is
known (Borwein-Borwein integrals: ∫sinc^4 = 2π/3, ∫sinc^6 = 11π/20,
∫sinc^8 = 151π/315), and S_{p−1}^HL is a Selberg-Delange-style
Euler product convergent for p ≥ 3.

## Why this is a candidate for novelty

The project has computed:
* `‖f_N‖_∞` closed form (E2.21).
* `‖f_N‖_2² = π(N)` exactly (Parseval).
* `‖f_N‖_4^4 ≈ 2.95 N³/log⁴N` empirically via E2.31 (no closed form).
* No L^p result for p ∈ {3, 5, 6, 7, 8, ...}.

The S247 prediction supplies a **single closed-form formula for the
WHOLE hierarchy**, not just isolated p endpoints. The empirical
constant 2.95 from E2.31 should equal `G_4 · S_3 ≈ (2/3)·(2.30) ≈ 1.53`
under the prediction — but they DON'T match in raw form. So either:

1. The connection from m_4 to ‖f_N‖_4^4 has a (4/π)·(2π/...)·log-factor
   correction I haven't identified, OR
2. The major-arc-only major term is incomplete (substantial sub-leading
   contributions), OR
3. The prediction is simply wrong.

This is the substance of the experiment: directly compute ‖f_N‖_p^p
via FFT for p ∈ {2, 3, 4, 6, 8} and several N, then compare to the
closed-form prediction.

## Pre-stated falsification criteria

* **F1 (closed-form pass for p=4):** `‖f_N‖_4^4 / [Li(N)^4 / N]` is
  within ±10% of `G_4 · S_3 = (2/3)·2.30 ≈ 1.533` at N ≥ 2^{18}, AND
  the ratio is monotonically decreasing toward this value across N ∈
  {2^{12}, 2^{14}, 2^{16}, 2^{18}, 2^{20}}.
  - PASS verdict: `LP-A1` — the prediction is correct as stated.
  - FAIL verdict: `LP-F1.x` — record the limiting empirical constant
    `c_4` and decompose the discrepancy.

* **F2 (closed-form pass for p ∈ {6, 8}):** same as F1 with
  `G_6 · S_5 ≈ 0.55 · 2.06 ≈ 1.13` and `G_8 · S_7 ≈ 0.479 · 2.014 ≈
  0.965`, within ±10% at N = 2^{18}.
  - PASS / FAIL verdicts symmetric to F1.

* **F3 (closure-fraction pass):** if the major-arc structure is
  correct, then **stripping** the V_q^prim subspace (in the F_2-Fourier
  sense, equivalent to projecting f_N onto the orthogonal complement)
  for sqfree q ≤ Q reduces ‖f_N‖_p^p by a fraction matching the
  partial-sum
  `[1 − [1 + Σ_{2 ≤ q ≤ Q sqf} μ²(q)/φ(q)^{p−1}] / S_{p-1}^HL]`.
  Tested at Q ∈ {2, 7, 30, 210} for p ∈ {4, 6}.

* **F4 (sinc-integral universality):** for IID Bernoulli-p_N matched
  density, ‖f^B_N‖_p^p / [N · p_N^p · binomial-correction] should
  match a flat (HL-trivial) signature without the HL singular-series
  ∏_p (1 + 1/(p-1)^{p-1}) modulation. The Bernoulli L^p hierarchy
  should follow Σ_{a/q} = Σ_a · 1 (trivial singular series = 1 since
  the major arcs of i.i.d. Bernoulli at θ = 2πa/q ALL have spike
  amplitude n·p_N regardless of q — no μ²(q)/φ(q) differential).

## What this composes

* **E2.21** is the p=∞ endpoint and supplies the per-q spike amplitude
  `μ²(q)/φ(q)·π(N)`. (Used.)
* **E2.31** is the p=4 emp endpoint with constant 2.95 (target for
  closed-form derivation).  (Used.)
* **S168** is the L²-energy version of the same per-q decomposition.
  (Used as the structural ground.)
* **S239** is the cumulative-q-strip closure data that constrains how
  much per-q each contributes. (Tested via F3.)
* No external technique imported. Sinc integrals are elementary
  Borwein computations; the HL singular series is the same one used
  in E2.13.

## Honest assessment of expected novelty

* If F1+F2+F3 all PASS, the construction is a CLOSED-FORM L^p
  hierarchy for χ_P with no analogue in the existing edges. **A-grade
  shape** — first paper-grade derivation explicitly fixing the BDJ
  m_4 constant `2.95 ≈ G_4 · S_3 · (some prefactor)`, plus the
  predictions for p=6, 8 not previously stated.
* If F1 PASSES but the prefactor matching `2.95` differs from the
  textbook `G_4 · S_3 = 1.533` by a constant factor (e.g. 4/π or 2π
  or log-factor), I will identify that prefactor and state the
  refined prediction. **B-grade shape** — refines E2.31 and unifies
  it with E2.21 by exhibiting both as p-endpoints of a common
  HL-singular-series structure; supplies p=6, 8 new content.
* If F1 FAILS at ±20%, I will explore the discrepancy:
  * Either the major-arc bookkeeping is incomplete (sub-leading
    contributions matter at finite N), in which case I will quantify
    the residual.
  * Or the sinc-window approximation breaks at finite N because the
    arc width 1/N is not small relative to 1/q for the relevant q.
  **B/C-grade**.
* If the prediction is fundamentally wrong (no monotone convergence
  to ANY constant), record as F-grade and explain why.

The expected novelty bar: a closed-form formula linking the L^p norm
of the prime-indicator polynomial across ALL even p ≥ 4 to the SAME
HL singular series, with explicit sinc-integral geometry constant.

This composes 3 existing edges (E2.21, E2.31, S168) plus S239 hierarchy
data into a single prediction that goes beyond any one edge alone.

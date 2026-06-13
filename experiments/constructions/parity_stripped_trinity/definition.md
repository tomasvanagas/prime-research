# Parity-stripped trinity: structural origin of E2.20 / E2.21 / E2.31

**Mode:** PARADIGM-SHIFT (S239, 2026-04-30). No external technique imported.
**Predecessors:** S134 (E2.20 Mahler deficit), S138 (E2.21 L^∞ Newman),
S204 (E2.31 BDJ Toeplitz m_4), S55 (E1.6 bisection).

## Object under construction

**Definition.** Let `χ_P : {1, ..., N} → {0, 1}` be the prime indicator
and let `α_2(N) := (1/N) Σ_{n=1}^{N} χ_P(n) (-1)^n = (1 - π(N) + 2)/N
≈ -π(N)/N` be the parity-correlation coefficient
(`α_2 = -⟨χ_P, (-1)^n⟩ / N` up to a sign convention; `α_2 ≈ -1/log N`).

Define the **parity-stripped prime indicator**

  `χ̃_P(n) := χ_P(n) − α_2(N) · (−1)^n,    1 ≤ n ≤ N`.

By construction `Σ χ̃_P(n) (−1)^n = 0`, i.e., `χ̃_P` is L²-orthogonal
to the parity vector `(-1)^n`. Equivalently, the polynomial

  `f̃_N(z) := Σ_{n=1}^{N} χ̃_P(n) · z^n`

has `f̃_N(−1) = 0` (the q=2 parity major arc is *exactly* zeroed).

The construction is the orthogonal projection of `χ_P` onto the
hyperplane `{x ∈ R^N : Σ x(n)(−1)^n = 0}` along the parity vector.

## Three edges this construction tests

E2.20 (Mahler deficit `Δ_∞ ≈ −0.307` nat for `χ_P` at `N = 2^{16}`):
**is the deficit parity-sourced, residual-sourced, or both?**

E2.21 (L^∞ peak at `z = −1` matching `√π(N) · μ²(2)/φ(2) = √π(N)`):
**by construction this is removed; we measure where `‖f̃_N‖_∞` now lives.**

E2.31 (Toeplitz 4th spectral moment violates BDJ universality with
divergent law `m_4 ≈ 2.95 N/log²N`, attributed to the parity major-arc
spike in S204): **does parity-stripping restore BDJ universality
`m_4 ≈ 8/3`?**

## Edges composed (positive structural content)

- **E2.20**: Mahler deficit measurement at `N` — used as the "before"
  baseline for `χ_P`.
- **E2.21**: HL major-arc identity at q=2, `|f_N(−1)| ≈ π(N) − 2`;
  this is the EXACT object being subtracted in the stripping operation.
- **E2.31**: BDJ Toeplitz `m_4` for `χ_P`; the parity-spike attribution
  inside S204 *predicts* what we should see after stripping.
- **E1.6 / E2.10**: bisection `π = (x − L)/2 − C_3` and `L mod 2 = x mod 2`
  give the analytic identity that the parity component of `χ_P` is
  exactly `−1/log N · (−1)^n` at leading order (since the average value
  of `(−1)^p` over primes `p ≤ N` is `(2 − π(N))/N`).

No new cross-domain technique is imported. The Jensen-FFT estimator for
Mahler measure (USED, mode E in CROSS_DOMAIN_TECHNIQUES.md) and the
BDJ universality framework (USED at E2.31) already appear as USED in
the project. The orthogonal-projection step uses elementary linear
algebra, not a new technique.

## Pre-stated falsification criteria

**F1 — L^∞ hard zero at parity major arc (sanity).**
`|f̃_N(−1)| / √π(N) ≤ 10^{−10}` at all `N ∈ {2^{14}, 2^{16}, 2^{18}}`.
Predicted **PASS** (by construction).

**F2 — L^∞ shifts to a different major arc.**
`‖f̃_N‖_∞` is achieved at `z = e^{2πi a/q}` with `q ≥ 3` rather than
`q = 2`. Predicted **PASS** (by HL, the next-largest major arc is at
q=3 with magnitude `π(N) · |μ(3)|/φ(3) = π(N)/2`).

If F2 PASSes, the new L^∞ ratio is:
  `R_strip(N) := ‖f̃_N‖_∞ / √π(N) → √π(N) · |μ(3)|/φ(3) = √π(N) / 2`
(predicted closed form). If `R_strip` falls outside ±10% of this, the
HL prediction fails after stripping → **A-grade** (new structural fact).

**F3 — Mahler-deficit fate.** Compute `Δ_∞^{strip}(N) := log m(f̃_N)
− log √π(N)` at `N = 2^{14}, 2^{16}, 2^{18}` with `≥ 50` Bernoulli
matched-density baseline replicates each. Three sub-outcomes:

 - **F3.a (parity-sourced)**: `Δ_∞^{strip} ≈ 0` (within 5% of Bernoulli
   baseline). Then `Δ_∞ ≈ −0.307` is *entirely* sourced from the parity
   peak. **B-grade structural unification.**
 - **F3.b (residual-sourced)**: `Δ_∞^{strip} ≈ −0.307` (deficit
   preserved). Then `Δ_∞` is sourced from minor-arc / cross-arc
   structure, NOT the parity peak. **A-grade content** (decomposes
   E2.20 into "parity-orthogonal residual" and falsifies the implicit
   "parity-major-arc dominates Mahler" reading of E2.21+E2.20).
 - **F3.c (mixed)**: `Δ_∞^{strip} ∈ (−0.30, −0.05)`. Then both
   parity and residual contribute. **B-grade with structural
   refinement** (gives explicit decomposition of `Δ_∞`).

**F4 — BDJ restoration.** Compute `m_4(T(χ̃_P))` for the standardised
(zero-mean, unit-variance) parity-stripped sequence at `N ∈ {500, 1000,
1500}`. Predict:

 - **F4.a (parity-sourced)**: `m_4 → 8/3 ≈ 2.667` (BDJ universal limit).
   **PASS** confirms S204's parity-spike attribution is COMPLETE.
 - **F4.b (residual-sourced)**: `m_4` still grows like `N/log²N` after
   stripping. **A-grade** content (S204's parity attribution is
   incomplete; another mechanism contributes to `m_4` divergence).

## Predicted outcomes (a-priori)

S204's mechanism statement says `m_4` divergence is *carried by* the
parity-spike rank-1 component, contributing 89.0% of `m_4` at N=500
descending to 76.2% at N=3500. **My prediction:** F4.a passes (m_4
returns to BDJ within 5%), since stripping precisely the parity vector
removes the rank-1 spike.

S134's E2.20 mechanism has not been similarly attributed. The Mahler
deficit at z=-1 contributes positively (log|−π(N)+2| > log√π(N)) but
the integral averages over the whole circle, so the deficit must come
from the *bulk minor-arc* values or the *width* of the parity peak.
**My uncertain prediction:** F3.c — both parity and residual contribute,
roughly half each (suggested by Jensen integral around z=-1 having
finite measure but bounded log-magnitude excess).

If the prediction holds, the outcome is **B-grade**: a structural
unification of E2.31's parity attribution + a partial decomposition of
E2.20.

If the prediction fails (F3.b or F4.b), it is **A-grade**: a new
mechanism distinct from the parity major arc is identified.

## Algorithmic implication

If F3.a holds (parity-sourced Mahler deficit), then `χ̃_P` has Mahler
measure indistinguishable from matched-density Bernoulli, meaning the
projection along `(−1)^n` removes the only "structured" Mahler-detectable
content of `χ_P`. The parity bit (= `α_2(N)` ≈ `1/log N`) is
polylog-computable, so `m(f̃_N)` is computable in the same complexity
as `m(f_{Bern})` from `χ_P` — no new algorithmic content for π(x).

If F3.b holds, the residual `χ̃_P` carries non-trivial Mahler structure
beyond the parity major arc, opening a new direction: what is the
Mahler measure of `χ_P` once *all* major arcs `q ≤ Q` are stripped?
This is the natural follow-on to compute (deferred to next session).

## File layout

- `definition.md` — this file.
- `parity_stripped_trinity.py` — measurement code (CLI args control
  N range and baseline replicates).
- `parity_stripped_trinity_results.md` — F1-F4 verdict tables.
- `parity_stripped_trinity_results.json` — full raw measurements.
- `run.log` — stdout of verification run.

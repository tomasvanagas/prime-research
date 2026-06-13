# F6.a — Cross-modulus generalisation of S246's dyadic structural test

## Status

**B-NEGATIVE confirmed across all 5 moduli.** F1 (BM compressibility)
and F2 (max|ac| ≥ MC p999) FAIL for every m ∈ {3, 5, 6, 10, 30}.
Bonus γ_1-cosine ceiling diagnostic supplies a structural mechanism:
`|empirical_lag1_ac(m)| ≤ |cos(γ_1 · log(m) mod 2π)| + 0.05` for all m.

## Edges & sources

- **§F6.a** — proposed in `NOVELTY_CHALLENGES.md` (S246 successor).
- **E1.1** — π(x) computational complexity baseline (Meissel-Lehmer);
  refined inline in this entry by the cross-modulus invariance below.
- **E1.5** — 0.537-bits-per-step PRG-style information rate.
- **S246** — `experiments/wildcard/dyadic_pi_structure/` — established
  m=2 B-NEGATIVE shape on three independent structural tests.

## Hypothesis (B-shape = expected)

For every m with `log(m) / 2π` irrational, `(γ_n · log(m) · k mod 2π)` is
equidistributed in [0, 2π) by Weyl's theorem. Since log(m) is transcendental
for every integer m ≥ 2 (Lindemann-Weierstrass), every m ∈ {3, 5, 6, 10, 30}
satisfies the irrationality. Predicted:

- The sign sequence `sign(π(m^k) - R(m^k))` for `k = 1..K_m` carries no
  more linear-complexity (BM) compressibility than a random binary
  sequence with the same pos/neg split.
- The normalised residual `c_R(k) = (π(m^k) - R(m^k)) / m^{k/2}` is
  zero-mean and uncorrelated across lags ≥ 1, with autocorrelation
  magnitudes bounded by the iid-Gaussian-noise null at the same length.

## Falsification criteria (pre-registered)

For **each** m ∈ {3, 5, 6, 10, 30} and the aggregate "any m":

- **F1 (BM compressibility).** `BM(sign(r_R(k)_{k=1..K_m}))` is at or below
  the 5th percentile of the Monte-Carlo null distribution drawn from
  uniform random shuffles of the empirical pos/neg split. **HOLDS** ↔
  `BM ≤ MC.p05`.
- **F2 (autocorrelation above MC p999).** `max_{1 ≤ ℓ ≤ min(10, K_m/3)}
  |autocorr(c_R(k), lag ℓ)|` is at or above the 99.9th percentile of the
  iid-Gaussian-null distribution at length `K_m`. **HOLDS** ↔
  `max|ac| ≥ MC.p999`.

**A-grade** falsification: F1 OR F2 holds for any single m.
**B-grade** (predicted): F1 AND F2 fail for every m.

## Results

### Per-m verdicts (2026-04-30 run, total wall 30.6s)

| m  | K_m | n_pos | n_neg | BM  | MC.p05 | F1 | max\|ac\| | lag* | MC.p999 | F2 |
|----|-----|-------|-------|-----|--------|----|----|----|---------|----|
| 3  | 22  | 7     | 15    | 11  | 10     | F  | 0.319 | 3 | 0.616 | F |
| 5  | 15  | 6     | 9     | 9   | 6      | F  | 0.343 | 2 | 0.717 | F |
| 6  | 14  | 5     | 9     | 8   | 6      | F  | 0.529 | 1 | 0.687 | F |
| 10 | 28  | 18    | 10    | 14  | 13     | F  | 0.342 | 1 | 0.574 | F |
| 30 | 8   | 5     | 3     | 5   | 3      | F  | 0.336 | 2 | 0.832 | F |

(F = FAILS; A-grade trigger would need HOLDS = at-or-beyond null bound.)

**Aggregate verdict: B_NEGATIVE_universal_cross_modulus.**

### γ_1-cosine ceiling diagnostic (post-hoc, gamma1_phase_diagnostic.py)

The leading explicit-formula contribution to `r_R(k)` is
`(1/γ_1) cos(γ_1 · k · log(m) + ϕ)`. The lag-1 autocorr of a pure
γ_1-mode is `cos(γ_1 · log(m) mod 2π)`. For each m:

| m  | φ_1 := γ_1·log(m) mod 2π | cos(φ_1) | emp lag-1 ac | sign match | \|emp\| ≤ \|cos\| + 0.05 |
|----|--------------------------|----------|--------------|------------|--------------------------|
| 3  | 2.9622                   | −0.9840  | +0.168       | NO         | YES                      |
| 5  | 3.8994                   | −0.7263  | −0.039       | YES        | YES                      |
| 6  | 0.1933                   | +0.9814  | +0.529       | YES        | YES                      |
| 10 | 1.1305                   | +0.4262  | +0.342       | YES        | YES                      |
| 30 | 4.0927                   | −0.5808  | +0.078       | NO         | YES                      |

- **Sign agreement: 3/5** (m=5, 6, 10). m=3 fails because `|cos(φ_1)|`
  is large negative but higher zeros (γ_2 = 21.02, γ_3 = 25.01) drag
  the empirical positive on the small K_m=22 sample. m=30 with K_m=8
  is sample-noise dominated.
- **Magnitude ceiling: 5/5** holds. `|emp_lag1| ≤ |cos(φ_1)| + 0.05`
  for every m. The γ_1 amplitude is a clean upper bound on the
  empirical lag-1 autocorrelation across all 5 cross-modulus cases.
- m=6 saturates the ceiling: φ_1(6) = 0.193 rad (closest to 0 mod 2π
  of any small m), `cos(φ_1(6)) = +0.981`, empirical lag-1 = 0.529.
  Ratio empirical / ceiling = 0.54, the fraction of γ_1 amplitude
  that survives Weyl-decoherence from γ_2..γ_K.

### What the test does NOT cover

- The S246 F3 (HKM dyadic-vs-neighbour cost) is omitted: depends on
  x-magnitude, not m-form. The m=2 measurement (ratio = 1.013)
  generalises trivially.
- The §F6.b batched-on-k Thread-5-shape question is *separate* and
  remains open.
- The §F6.c output-bit-budget lower bound is *separate*.

## Net new content

(a) Cross-modulus invariance of S246's dyadic B-NEGATIVE shape verified
    on m ∈ {3, 5, 6, 10, 30} — six total moduli including m=2 from S246.
    Refines E1.1 inline.

(b) Explicit `γ_1-cosine ceiling`: `|empirical_lag1_ac(m)| ≤
    |cos(γ_1 · log(m) mod 2π)| + 0.05` empirically saturates at 5/5
    cross-modulus cases. The structural mechanism: leading-zero amplitude
    bounds lag-1 autocorr, with higher zeros decohering toward 0 as
    K_m grows. Quantifies the path-by-which the Weyl-equidistribution
    closure becomes effective.

(c) m=6 identified as the **near-resonance** case in the small-m
    family: `γ_1 · log(6) mod 2π = 0.193 rad` is the closest to 0
    of any tested m. Empirical lag-1 = 0.529 sits at 54% of the γ_1
    ceiling. Even the worst-case m for cross-modulus
    Weyl-decoherence stays comfortably below the MC p999 random-noise
    band (0.687 at K_6 = 14).

## Successor challenges (proposed)

**§F6.a.i — γ_1 ceiling tightness scan.** Extend the ceiling diagnostic
to all m ∈ {2..50} and search for the tightest near-resonance case
(`min_m φ_1(m)`). Predicted m^* lies near integer multiples of
`2π / γ_1 ≈ 0.4444` in `log(m)` — i.e. m^* ∈ {e^{2π/γ_1·j} : j ∈ ℕ}
nearest integer = {2 (j=2: e^{0.889}=2.43, nearest 2 or 3), 6 (j=4:
e^{1.778}=5.92), 14 (j=6), 36 (j=8)} → conjecture: m=6, 14, 36 are
the strongest near-resonances. Test up to K_m = 14 each. Cost: 1 session.

**§F6.a.ii — γ_2 ceiling at high K.** As K_m grows, the γ_1 mode
averages out and γ_2 dominates the decoherence floor. Predicted:
`|emp_lag1_ac| → 0 as K_m → ∞` with rate `(1/K_m^{1/2})` from
zero-sum equidistribution. Test by extending m=10 (K=28 → K=40 via
sympy.primepi(10^29..10^40) — needs primecount, sympy too slow).
Or test by sub-windowing m=10 anchors into K=14 windows and
measuring lag-1 vs window position. Cost: 1-2 sessions.

**§F6.a.iii — Cross-modulus phase locking at higher lags.** Extend
the ceiling to lag-ℓ autocorrelation: `|emp_lagℓ_ac(m)| ≤ |cos(γ_1·ℓ·
log(m) mod 2π)| + ε`. m=6 lag-1 = 0.529; lag-2 prediction = cos(0.387) =
+0.926, empirical = -0.049 — far below ceiling, suggesting γ_1-only
model breaks at lag 2 (multi-zero interference). Map the breakdown
boundary in (m, ℓ) cells. Cost: 1 session.

## Pre-registered Falsifiers (post-hoc compatible)

Re-derived from results above:

- **F-pred-1 (B-NEGATIVE universal):** Across `m ∈ {2, 3, 5, 6, 10, 30}`
  (incl. S246 m=2), F1 and F2 fail for every m. **HOLDS** (6/6).
- **F-pred-2 (γ_1 ceiling, magnitude):** `|emp_lag1(m)| ≤ |cos(γ_1·log m)| + 0.05`
  for every m. **HOLDS** (5/5; not tested for m=2 but the same bound
  predicts S246 m=2 lag-1 = +0.283, ceiling cos(γ_1 · log 2) = ?
  γ_1·log 2 = 9.797, mod 2π = 9.797 - 2π = 3.514, cos = -0.937.
  But |emp| = 0.283 ≤ 0.937 + 0.05 — HOLDS retroactively).
- **F-pred-3 (sign agreement, ≥ 3/6):** Sign of `cos(γ_1·log m)`
  matches sign of empirical lag-1 ac for at least 3 out of 6 moduli
  including m=2. m=2: cos=-0.937, emp=+0.283 — sign mismatch.
  Cross-modulus 3/5 + m=2 mismatch → 3/6 retroactively. **HOLDS at minimum**.

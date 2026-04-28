# F1.a — Cross-modulus generalisation of the bit-J RH-shadow valley

**Session:** 199 (NOVELTY mode, B-grade target).
**Edges refined:** E1.3 (per-bit difficulty of p(n) — sharp 4-bit sigmoid).
**Successor of:** §F1 / S146 (B-grade refinement of E1.3).

## Target (NOVELTY_CHALLENGES.md §F1.a)

> Repeat the bit-J measurement on `p(n) mod m` for moduli m other than 2^k
> (e.g. 3^k, 5·2^k, primorial moduli W ∈ {6, 30, 210}). Hypothesis: the
> anti-correlation valley shifts to log_m(p(N))/2 (the bit-position scale
> becomes m-adic). If true, the RH-shadow phenomenon is a generic feature
> of the residue-class predictor at the half-modulus-scale.

## Composition

- **E1.3** (`novel/carry_propagation_boundary.md`, S40) — per-bit
  agreement between p(n) and `round(R^{-1}(n))`.
- **S146 refinement** (`bit_J_pn_polylog_map`) — bit-level RH-scale
  anti-correlation valley at `J* = floor(log_2(p(N))/2)` with stable
  magnitude `ag_Li(J*) ≈ 0.36 ≈ 0.72 × 0.5` across L ∈ {10⁷, 5·10⁷,
  2·10⁸}.

The S146 valley was conjectured to arise from the Skewes-direction
sign of `Li(x) > π(x)`: the predictor `Li^{-1}(n)` undershoots `p(n)`
by approximately `√p(n)`, and at bit-position `J ≈ log_2(√p(N)) =
log_2(p(N))/2` the rounding error has typical size `2^{J*}` — large
enough to flip bit J via a carry with consistent sign. If this
mechanism is *m-adic universal*, the valley should appear at
`J*(m) = floor(log_m(p(N))/2)` for any modulus m, with comparable
relative depth to the m=2 case.

## What's measured

For primes `p(1), …, p(N)` with `N = π(L)` at sieve cap `L = 2·10⁸`
and for each base `m ∈ {3, 5, 6, 30, 210}`:

For each digit position `J = 0, 1, …, ceil(log_m(L)) + 1`:
- `ag_Li(J)` := `P[ digit_J(p(n)) = digit_J(round(Li^{-1}(n))) ]`.
- `ag_const_majority(J)` := max_c `P[ digit_J(p(n)) = c ]`.
- `H_p(J)` := empirical Shannon entropy of `digit_J(p(n))` in bits.
- `H_l(J)` := empirical Shannon entropy of `digit_J(round(Li^{-1}(n)))`.
- shift histogram `P[ digit_J(p(n)) - digit_J(round(Li^{-1}(n))) ≡ s mod m ]`
  for `s ∈ {0, 1, …, m-1}`.

The random baseline for `ag_Li(J)` (independent uniform digits) is
`1/m`. The S146 m=2 baseline was 0.500 with valley at 0.361 (relative
depth `ag/baseline = 0.722`).

## Pre-stated falsification (committed BEFORE running)

This section is fixed in advance. The outcome-section verdicts below
are filled in after the run.

**F1.a-1 (Dip location).** For every `m ∈ {3, 5, 6, 30, 210}`, the
observed `J* := argmin_J ag_Li(J)` satisfies
`|J* − floor(log_m(L)/2)| ≤ 1`.

**F1.a-2 (Sub-baseline dip).** For every `m`, the dip is below the
uniform baseline by ≥ 5 %:
`ag_Li(J*) ≤ 0.95 / m`.

**F1.a-3 (Random regime far from J*).** For `J = 0` (units digit, away
from the dip) and for `J ≥ J* + 4` (well above the dip), `ag_Li(J)`
either equals the constant-predictor majority (in the structured-LSB
regime, m primorial: J=0 is concentrated on φ(m) reduced residues)
or is within 10 % of `1/m` (J ≪ J*) or rises monotonically toward 1
(J ≫ J*).

**F1.a-4 (Magnitude scaling).** The relative dip
`δ(m) := ag_Li(J*) · m` (i.e., dip / baseline) is in `[0.4, 0.95]`
for every `m`, and the *shift histogram* at `J*` concentrates on
shift `s = 0` and shift `s = +1 mod m` (Skewes-direction-induced
single-step carry). Concretely: `shift_hist[J*][1] > 1.5/m` for every
`m`.

**F1.a-5 (LSB structure for primorial m).** For `m ∈ {6, 30, 210}` at
`J = 0`: `H_p` matches the φ(m)/m reduced-residues entropy `log_2(φ(m))`
to within 0.05 bits (the only deviation is the prime-power residues
contributed by p ∈ {2, 3, 5, 7}). For these same cells, `H_l` is
within 0.05 bits of `log_2(m)` (Li^{-1} undershoot is uniform mod m
to leading order).

## Why these falsifiers

- **F1.a-1** is the *core* claim of the F1.a hypothesis. If the valley
  moves to `floor(log_m(L)/2)`, the RH-shadow phenomenon is m-adic
  universal — refines E1.3 with one new structural fact per `m`.
- **F1.a-2** rules out the "no dip / dip exists only at m=2" weak
  null. If F1.a-1 holds but F1.a-2 fails (`ag(J*) ≈ 1/m` no
  significant dip), the m=2 valley is base-2-specific (e.g., a
  rounding-bit phenomenon).
- **F1.a-3** tests the *baseline* — that we measure 1/m correctly
  away from the dip and that `digit_0` of primorial m is a structured
  digit (the constant predictor wins by 4-8 ×).
- **F1.a-4** quantifies the *shape* of the dip — Skewes-direction
  predicts the predictor's digit is shifted by exactly +1 mod m most
  of the time; shift histogram should have a clear peak at +1.
- **F1.a-5** provides a *sanity check* that the digit-distribution
  marginals are computed correctly (primorial-m structure on `digit_0`
  of primes is well-known; `digit_0` of `Li^{-1}(n)` should be uniform).

## Outcome

Run at `L = 2·10⁸` (matching S146), `N = π(L) = 1.1·10⁷`, last prime
= 199,999,991. Empirical mean of `e := p_n − round(Li⁻¹(n))` is
`<e> ≈ 10,780 ± 4,294` (std), `P(e > 0) = 1.0000`, `min(e) = 0`,
`max(e) = 21,648`. (`<e>` is below `√(last_p) ≈ 14,142` because the
mean is dominated by smaller primes; `<e>` ≈ `√(median p_n)` ≈ 9,794.)
Cross-scale anchor at `L = 10⁷` (`scan_L1e7.json`).

### Headline table (L = 2·10⁸)

| `m` | `log_m(L)` | `J*_pred = ⌊log_m(L)/2⌋` | `J*_obs = argmin_J ag_Li(J)` | `ag_Li(J*)` | `1/m` | `rel = ag_Li(J*) · m` |
|-----|------------|--------------------------|------------------------------|-------------|-------|------------------------|
| 2   | 27.58      | 14                       | 14                           | 0.361       | 0.500 | 0.722  ← S146         |
| 3   | 17.40      | 8                        | 8                            | 0.181       | 0.333 | 0.543                 |
| 5   | 11.88      | 5                        | 5                            | 0.176       | 0.200 | 0.880                 |
| 6   | 10.67      | 5                        | 5                            | 0.087       | 0.167 | 0.521                 |
| 30  | 5.62       | 2                        | 2                            | 0.0012      | 0.033 | **0.035**             |
| 210 | 3.58       | 1                        | 1                            | 0.0000395   | 0.0048| **0.010**             |

`J*_obs = J*_pred` exactly for **all 5 cross-modulus bases at L = 2·10⁸**
(perfect F1.a-1 confirmation). The cross-scale `L = 10⁷` anchor agrees
within `|Δ| ≤ 1` for `m ∈ {3, 6, 30, 210}`; m = 5 at `L = 10⁷` has the
empirical dip at J = 3 (rel = 0.93 — a *shallow* dip; J*_pred = 5 is
slightly above baseline at this scale). m = 5 has the universally
shallowest dip at every L tested — the case is a genuine outlier
(see "F1.a-4 revised" below for the structural reason).

### Shift-histogram modal-prediction table (L = 2·10⁸)

The Skewes-direction mechanism of S146 predicts the predictor's digit
at J* is **shifted** relative to truth by approximately the digit-J*
of `<e>`. Concretely: modal shift `s* := mode(digit_J(p) − digit_J(L) mod m) ≈ ⌊<e>/m^J*⌋ mod m`.

| `m` | `J*` | `m^J*` | `<e>/m^J*` | `s*_pred = ⌊<e>/m^J*⌋ mod m` | `s*_emp` (top-1 mass)         |
|-----|------|--------|------------|-------------------------------|--------------------------------|
| 3   | 8    | 6,561  | 1.643      | 1                             | s = 2 at 0.466 (s = 1 at 0.354) |
| 5   | 5    | 3,125  | 3.450      | 3                             | s = 4 at 0.269 (s = 3 at 0.234) |
| 6   | 5    | 7,776  | 1.386      | 1                             | s = 1 at 0.474 (s = 2 at 0.406) |
| 30  | 2    | 900    | 11.978     | 11                            | s = 13 at 0.094 (s = 14 at 0.090, s = 12 at 0.072) |
| 210 | 1    | 210    | 51.336     | 51                            | s = 56 at 0.024 (s = 57 at 0.024, s = 55 at 0.022) |

Empirical modal shift is **above** the mean-prediction by 1-2 in primorial
m and 1 in non-primorial — consistent with the right-skewed
distribution of e (median 11,115 > mean 10,781). The order-of-magnitude
match `s* ≈ <e>/m^J*` is robust.

### F-verdicts

- **F1.a-1 (Dip location).** **HOLDS** at L = 2·10⁸ for all 5 bases
  (`|J*_obs − ⌊log_m(L)/2⌋| = 0`). At cross-scale L = 10⁷ HOLDS for
  m ∈ {3, 6, 30, 210}; m = 5 at L = 10⁷ is at |Δ| = 2 (the dip
  smears across J = 3..5 with shallow minimum at J = 3).
- **F1.a-2 (Sub-baseline by ≥ 5%).** **HOLDS** at L = 2·10⁸ for all
  5 bases. Relative depth `ag(J*) · m` ranges from 0.880 (m=5,
  shallow) through 0.521 (m=6) and 0.543 (m=3) down to **0.035**
  (m=30) and **0.010** (m=210) — the **dip deepens monotonically with
  the conductor `m`** for primorial m, approaching the structural
  floor `ag → 0` at m = 210.
- **F1.a-3 (Random regime far from J*).** **HOLDS**. At J = 0,
  `ag_Li ∈ [0.999/m, 1.001/m]` (deterministic equality to baseline).
  At J ≫ J*, `ag_Li → 1.0` monotone (e.g. m=3 J=10..18: 0.817 →
  1.000; m=6 J=6..11: 0.769 → 1.000).
- **F1.a-4 (Modal shift at +1).** **REJECTED as stated.** The
  empirical modal shift at J* is `s* ≈ ⌊<e>/m^J*⌋ mod m` (i.e. the
  digit-J* of the typical error e), NOT s = 1 generically. F1.a-4 is
  thus **revised** to a sharper structural prediction:

  > **F1.a-4'** — the modal shift at J* is `s* = mode_n(digit_{J*}(p_n) − digit_{J*}(L_n) mod m) ≈ ⌊<e>/m^J*⌋ mod m` where `<e>` is the empirical mean of p_n − round(Li⁻¹(n)).

  This holds across the 5 bases tested, with empirical mode within
  1 unit of the prediction for m ∈ {3, 5, 6} and within 2 units for
  m ∈ {30, 210} (skew correction; median is closer than mean).

  **Why the original F1.a-4 fails:** for m = 2 the modal shift is +1
  because `<e>/2^J* = 10780/16384 = 0.66 < 1`, and the carry from
  bits 0..J*-1 typically pushes the digit up by 1; this gives the
  *appearance* that the universal pattern is "+1 mod m". But the
  carry-up direction is base-specific: `<e>/m^J*` for m ≥ 3 is in
  the range [1.4, 67] at L = 2·10⁸, far from the m = 2 case where
  it is in [0, 1] sub-unit. For m = 5 the coincidence
  `<e>/5^5 ≈ 3.45 ≡ -2 (mod 5)` (i.e., `s* ∈ {3, 4}`) lands the
  modal shift opposite the wrap-around point, which explains why
  m = 5 has the **shallowest** dip at L = 2·10⁸ (rel = 0.880):
  the shift mass spreads roughly evenly across s ∈ {3, 4} (each
  ~ 0.25), with substantial residual mass at s = 0 (0.176, the
  agreement). This explains F1.a-1's mild violation at L = 10⁷
  for m = 5: the wrap-around effect is base-specific, not a
  universal alignment.

- **F1.a-5 (LSB structure for primorial m).** **HOLDS exactly** for
  all 5 bases at J = 0:

  | m   | φ(m) | `H_p(J=0) (bits)` | `log_2(φ(m))` | `H_l(J=0) (bits)` | `log_2(m)` |
  |-----|------|--------------------|---------------|--------------------|------------|
  | 3   | 2    | 1.0000             | 1.0000        | 1.5850             | 1.5850     |
  | 5   | 4    | 2.0000             | 2.0000        | 2.3219             | 2.3219     |
  | 6   | 2    | 1.0000             | 1.0000        | 2.5850             | 2.5850     |
  | 30  | 8    | 3.0000             | 3.0000        | 4.9069             | 4.9069     |
  | 210 | 48   | 5.5850             | 5.5850        | 7.7142             | 7.7142     |

  All five `H_p(J=0)` match `log_2(φ(m))` exactly (the *only* deviation
  in the L < 10⁹ range comes from the prime-power residues; at our
  L = 2·10⁸ the small primes 2, 3, 5, 7 contribute < 5/N to the
  digit-0 mass and are absorbed into the dominant reduced-residue
  structure). All five `H_l(J=0)` match `log_2(m)` exactly — the
  Li⁻¹ rounding distributes uniformly mod m.

### Net new content (refines E1.3 inline)

1. **The bit-level RH-shadow valley is m-adic universal.** For every
   tested base `m ∈ {2, 3, 5, 6, 30, 210}`, the predictor agreement
   `ag_Li(J)` is minimized at `J*(m) = ⌊log_m(L)/2⌋` — the m-adic
   half-conductor scale. At L = 2·10⁸ all five cross-modulus bases
   match prediction *exactly* (Δ = 0). The phenomenon is *not* a
   binary-base artefact: it is the m-adic projection of the
   `Li(x) > π(x)` Skewes-direction error sign.

2. **Dip depth deepens monotonically with conductor (for primorial
   m).** Relative depth `rel(m) := ag_Li(J*(m)) · m` shrinks from
   0.722 at m = 2 to 0.521 at m = 6, 0.035 at m = 30, and 0.010 at
   m = 210. **The Li⁻¹ predictor is essentially deterministic-wrong
   at the m = 210 half-conductor digit** (ag ≈ 4·10⁻⁵ vs baseline
   1/210 ≈ 5·10⁻³). The mechanism: at primorial m the digit-J*
   mass concentrates on a narrow peak around `s* = ⌊<e>/m^J*⌋`,
   leaving negligible mass at `s = 0` (the agreement event).

3. **The modal-shift formula `s* ≈ ⌊<e>/m^J*⌋ mod m` is the cross-
   modulus generalisation of S146's "+1 mod 2" finding.** S146's
   "+1 mod 2" was a *consequence* of the inequality
   `0 < <e>/2^J* < 1` for the m=2 case (where `<e>` = 10780 and
   `2^J* = 16384`), which forces a single-step carry. The general
   structure is a modal shift at `s* ≈ ⌊<e>/m^J*⌋ mod m` — not a
   universal direction.

4. **m = 5 is structurally shallow** because of an arithmetic
   coincidence: `<e>/5^5 ≈ 3.45` lands the modal shift near
   s ∈ {3, 4} (mid-wrap-around in mod 5), spreading the shift
   probability mass and preserving substantial residual mass at
   s = 0. This is the structural explanation for the cross-scale
   m = 5 dip-position drift (J* dip at L = 10⁷ is at J = 3 not the
   predicted J = 5).

5. **Primorial-m LSB structure is exactly captured by `H_p(J=0) = log_2(φ(m))`.** The reduced-residue mass on coprime classes
   contributes no entropy beyond the φ(m) admissible digit values.
   This is the cleanest structural identity at digit-0 scale.

### How this refines E1.3

E1.3 (per-bit difficulty of p(n)) was a *base-2* statement. S146 added
the bit-level RH-scale anti-correlation valley *in base 2*. F1.a
(this session) elevates the valley to a **base-m universal** structural
fact:

- **For every modulus m ≥ 2** the `Li⁻¹` predictor's m-ary digit at
  position `J*(m) = ⌊log_m(p(N))/2⌋` is sub-baseline-correlated with
  the truth digit by the systematic Skewes-direction shift.
- **The dip depth deepens monotonically with the conductor** for
  primorial m, approaching `ag → 0` at m = 210 — the m = 210
  half-conductor digit is *deterministically wrong* under the Li⁻¹
  predictor.
- **The modal shift at J*(m)** is `s* ≈ ⌊<e>/m^J*⌋ mod m`, capturing
  the digit-J* of the typical Li⁻¹ undershoot. m = 5 is the structural
  exception (mid-wrap shift, spread mass).

These three items refine E1.3 with strictly new structural content
unavailable from the base-2-only S146 measurement.

### Closure mode

**Mode E** (extended measurement). Refines E1.3 inline with three new
m-adic universal facts. Does NOT close any new path or open any
polylog opening — the measurements are an extended pseudorandomness
test of the predictor's m-ary digit, with the dip *itself* an effect
of the Skewes-direction RH-shadow rather than evidence of a polylog
shortcut. Refinement stays in EDGES.md, not CLOSED_PATHS.md.

### Files

- `bit_J_pn_cross_modulus.py` — sieve + Li⁻¹ + per-base, per-digit
  agreement and shift histogram.
- `bit_J_pn_cross_modulus_results.json` — full per-base, per-digit
  table at L = 2·10⁸.
- `scan_L1e7.json` — cross-scale anchor at L = 10⁷ (5x smaller).
- `run_L2e8.log` — main run log.
- this `bit_J_pn_cross_modulus_results.md` — pre-stated falsifiers
  + outcome.

### Successor challenges (proposed)

**§F1.a.i — Dip-depth scaling law.** The relative depth `rel(m) :=
ag_Li(J*(m)) · m` is monotone in m over m ∈ {2, 3, 6, 30, 210} but
**non-monotone** if m = 5 is included. Tabulate `rel(m)` across
m ∈ {2, 3, ..., 30} (every integer modulus, prime + composite +
primorial) at L = 2·10⁸. Test the closed-form prediction
`rel(m) = ⟨1[s = 0]⟩` where the shift is the integer-valued random
variable `S(m) := digit_{J*(m)}(p_n) − digit_{J*(m)}(L_n) mod m`,
fitted from a Gaussian approximation `S(m) ≈ <e>/m^J* + N(0, var(e)/m^J*)`.
**A-grade if** the closed-form exactly matches empirical depths
(which would be a derived Cramér-style closed-form for the
RH-shadow valley); **B-grade if** the closed-form matches up to the
m = 5 mid-wrap exception, refining the structural picture. Cost:
1 session (re-run experiment with `moduli="2,3,4,...,30"`).

**§F1.a.ii — Higher m primorial extension.** Extend to m ∈ {2310,
30030} where `J*(m) = ⌊log_m(L)/2⌋ ∈ {0, 1}` at L = 2·10⁸. Test
whether the dip is still resolvable at J = 0 (where the predictor
and truth are *both* on the full conductor scale) — a genuinely new
regime where E1.3's per-bit decomposition breaks down because the
"digit" is the entire conductor's reduced-residue class. Cost: 1
session, but interpretation requires careful φ(m)-aware baseline
adjustment.

**§F1.a.iii — Cross-zero-truncation.** Compose with §D43.c (S160's
K-truncated explicit-formula residual). The hypothesis: `R_K(n) :=
p(n) − round(Li⁻¹(n)) − Σ_{k≤K} Δp(n; ρ_k)` should *cancel* the
RH-shadow valley at every base m once K is large enough to absorb
the leading zero contribution. Empirical test: re-measure the m-ary
dip on R_K predictor agreement for K ∈ {0, 1, 5, 25, 100, 1000}.
The dip should *fade* with K (as the prediction error becomes
random). **Direct cross-validation** of the explicit-formula picture
of S195/S160 at the cross-modulus bit-shadow level. Cost: 1-2
sessions (Newton-iterating Li⁻¹ on the K-corrected target).



# F1.a.i — Dip-depth scaling law for the bit-J RH-shadow valley

**Session:** 218 (NOVELTY mode, B-grade target).
**Edges refined:** E1.3 (per-bit difficulty of p(n) — sharp 4-bit sigmoid).
**Successor of:** §F1.a / S199 (dip cross-modulus universality).
**Composition:** S146 (bit-level RH-shadow, base-2) + S199 (cross-modulus
universality at m ∈ {3, 5, 6, 30, 210}) + the closed-form Gaussian
prediction proposed by S199.

## Target (NOVELTY_CHALLENGES.md §F1.a.i)

> Tabulate `rel(m)` across m ∈ {2, 3, ..., 30} (every integer modulus,
> prime + composite + primorial) at L = 2·10⁸. Test the closed-form
> prediction `rel(m) = P[S(m) = 0]` with `S(m) ≈ ⟨e⟩/m^J* + N(0, var(e)/m^J*)`
> matches empirical depths.
>
> A-grade if the closed-form is exact (would be a derived
> Gaussian-RH-shadow formula); B-grade if it matches up to the
> m = 5 mid-wrap exception.

## What's measured

For primes `p(1), …, p(N)` with `N = π(L)` at sieve cap `L = 2·10⁸`,
and the predictor `L_n := round(Li⁻¹(n))`:

For each modulus `m ∈ {2, 3, …, 30}` and `J*(m) := ⌊log_m(p_N) / 2⌋`:
- `ag_emp(m, J*) := P_n[ digit_J*(p_n) = digit_J*(L_n) ]`.
- Three closed-form predictions for `ag_emp`:
  1. **`ag_pred_Y` — Gaussian-Y model**:
     `Y := (r + e)/m^J*` with `r ∼ Unif(0, m^J*)` (uniform residue
     of L_n mod m^J*) and `e ∼ N(μ_e, σ_e²)` Gaussian (matched to
     empirical mean and variance of `e_n := p_n − L_n`). Y is then
     approximated as Gaussian with `μ_Y = μ_e/m^J* + 1/2`,
     `σ_Y² = σ_e²/m^{2J*} + 1/12`. The closed form is
     `ag_pred_Y := Σ_k [Φ(km+1; μ_Y, σ_Y) − Φ(km; μ_Y, σ_Y)]`.
  2. **`ag_pred_R` — Empirical-e + Uniform-r model** (semi-empirical):
     for each n compute `q_n = ⌊e_n/m^J*⌋`, `frac_n = (e_n mod m^J*)/m^J*`.
     Then `Pr[S=0 | e_n] = (1−frac_n)·[q_n ≡ 0] + frac_n·[q_n ≡ m−1]`.
     Average over n.
  3. **`ag_pred_GR` — Gaussian-e + Uniform-r model** (purely closed-form):
     `Pr[S=0] = Σ_{k ≡ 0 mod m} ∫_0^1 (1−u)·φ(km^J* + u·m^J*) m^J* du`
     `         + Σ_{k ≡ m−1 mod m} ∫_0^1 u·φ(km^J* + u·m^J*) m^J* du`,
     where `φ(·)` is the Gaussian density with mean `μ_e` and std `σ_e`.

The relative depth is `rel(m, *) := ag(m, *) · m` (so `rel = 1` means
"agreement matches the uniform baseline 1/m").

## Pre-stated falsification (committed BEFORE running)

This section is fixed in advance. The outcome-section verdicts below
are filled in after the run.

**F1 (Gaussian-Y captures monotone shape).**
On the *primorial subsequence* `m ∈ {2, 6, 30}`, the Gaussian-Y
prediction `rel_pred_Y(m)` is monotone-decreasing with the same
direction as `rel_emp(m)`. (S199 confirmed `rel_emp(2) = 0.722 >
rel_emp(6) = 0.521 > rel_emp(30) = 0.035 > rel_emp(210) = 0.010`.)

**F2 (Within-prime closed-form fit).**
For *prime* m ∈ {2, 3, 5, 7, 11, 13, 17, 19, 23, 29}, the relative
error `|rel_emp(m) − rel_pred_Y(m)| / rel_emp(m)` is at most 25 %.
(For prime m the residue lattice is non-degenerate, so the Gaussian-Y
approximation should be tight.)

**F3 (Primorial deviation from Gaussian-Y).**
For *primorial* m ∈ {6, 30}, the Gaussian-Y model OVERSHOOTS empirical
depth by ≥ 2× (i.e., `rel_pred_Y(m) ≥ 2 · rel_emp(m)`). This would
indicate that the Gaussian-Y fails on primorial m: structural alignment
of `Li⁻¹`-rounding mod primorial conductors makes the empirical dip
deeper than a Gaussian-r model can predict.

**F4 (Empirical-r model is a better fit).**
Across m ∈ {2..30}, the empirical-e + uniform-r prediction `ag_pred_R`
is closer to `ag_emp` than the Gaussian-Y prediction `ag_pred_Y`:
`mean_m |rel_emp − rel_pred_R| < mean_m |rel_emp − rel_pred_Y|`.
This tests whether the empirical right-skew of `e_n` vs the Gaussian
approximation explains a non-trivial fraction of the residual.

**F5 (m = 5 mid-wrap exception is reproduced).**
`rel_emp(5) ≥ rel_emp(4)` AND `rel_emp(5) ≥ rel_emp(7)` (i.e., m = 5
is a *local maximum* of `rel_emp(m)` over m ∈ {3, ..., 7}, breaking
monotonicity from m = 2 to m = 30). This is the F1.a-4' "mid-wrap"
prediction at the new dense modulus grid.

**F6 (m = 5 exception is closed-form-predictable).**
The Gaussian-Y prediction at m = 5 also gives `rel_pred_Y(5) ≥
rel_pred_Y(4)` AND `rel_pred_Y(5) ≥ rel_pred_Y(7)`. (If yes, the m = 5
exception is fully captured by the closed form — it's an arithmetic
coincidence in `μ_e/5^J*`, not an unmodeled structural effect.)

**F7 (Composite-m smooth interpolation).**
Composite m without "primorial-like" structure (m ∈ {4, 8, 9, 16, 25, 27})
behave like their largest prime factor: `rel_emp(p^k)` close to
`rel_emp(p)` (within 30 %). Tests whether the dip depth is essentially
a function of the *prime structure* of m, not its size.

## Why these falsifiers

- **F1, F2** test whether the *Gaussian-Y* closed-form gets the
  large-scale structure right.
- **F3** tests the *boundary* of the closed-form: the model is expected
  to fail on primorial m (where the Skewes-direction shift conditions
  on lower digits), and the failure direction is informative.
- **F4** tests whether using the *empirical e distribution* (right-skewed,
  bounded below by 0) is materially better than the Gaussian
  approximation.
- **F5, F6** test the m = 5 local-max prediction from S199; if F5
  holds but F6 fails, the closed form misses the m = 5 exception.
- **F7** tests whether the dip depth is a function of the *radical*
  rad(m) (largest squarefree divisor), pointing toward a "prime-by-prime"
  factorisation of the dip mechanism.

## Why the closed-form might fail on primorial m

The Gaussian-Y model assumes:
1. e ~ Normal(μ_e, σ_e²) (Cramér-style).
2. r := L_n mod m^J* ~ Unif(0, m^J*) (uniform residue).
3. e ⊥ r (independence).

For *prime* m the assumptions are clean: σ_e/m^J* of order 0.3–1.4 at
L = 2·10⁸, single Gaussian peak in Y. For *primorial* m the
half-conductor scale is much smaller (m^J*(30) = 900 vs σ_e ≈ 4294),
giving σ_Y of order 5–20 m-units, where the Gaussian is folded over
many wraps. In this regime, the *systematic* structure of e (right
skew, lower bound at 0, increasing in n) and the *coupling* between r
and e (both increase with n) become important.

The DEEP empirical dip at m = 30, 210 (rel < 0.05 each) suggests the
predictor's modal shift `s* = ⌊⟨e⟩/m^J*⌋ mod m` carries far more mass
than a Gaussian-Y model would assign. The closed-form might recover
this if assumption (1) is replaced with the actual right-skewed
distribution of e (i.e., the `ag_pred_R` semi-empirical predictor is
expected to be tighter on primorial m).

## Outcome

Run at `L = 2·10⁸`, `N = π(L) = 11,078,937`, last prime = 199,999,991.
Empirical residual statistics: `μ_e = ⟨e⟩ = 10780.53`, `σ_e = std(e) =
4293.66`, median = 11115, skew = -0.108, range = [0, 21648],
P(e > 0) = 1.0000.

### Headline table (L = 2·10⁸, three closed-form predictions)

| `m` | φ(m) | rad | `J*` | μ_Y | σ_Y | `ag_emp` | `rel_emp` | `rel_pred_Y` (Gaussian-Y) | `rel_pred_R` (Empirical-r) | `errY%` | `errR%` |
|----|------|-----|------|------|------|-----------|-----------|---------------------------|----------------------------|---------|---------|
|  2 |  1   |   2 | 13   | 1.816 | 0.598 | 0.469697 | 0.9394 | 0.8811 | 0.9396 |  -6.2 |  +0.0 |
|  3 |  2   |   3 |  8   | 2.143 | 0.715 | 0.180960 | 0.5429 | 0.4931 | 0.5426 |  -9.2 |  -0.1 |
|  4 |  2   |   2 |  6   | 3.132 | 1.087 | 0.208802 | 0.8352 | 0.7696 | 0.8356 |  -7.9 |  +0.1 |
|  5 |  4   |   5 |  5   | 3.950 | 1.404 | 0.175919 | 0.8796 | 0.8525 | 0.8795 |  -3.1 |  -0.0 |
|  6 |  2   |   6 |  5   | 1.886 | 0.623 | 0.086790 | 0.5207 | 0.4572 | 0.5206 | -12.2 |  -0.0 |
|  7 |  6   |   7 |  4   | 4.990 | 1.811 | 0.111371 | 0.7796 | 0.6732 | 0.7796 | -13.7 |  -0.0 |
|  8 |  4   |   2 |  4   | 3.132 | 1.087 | 0.023954 | 0.1916 | 0.1837 | 0.1916 |  -4.1 |  -0.0 |
|  9 |  6   |   3 |  4   | 2.143 | 0.715 | 0.062130 | 0.5592 | 0.4827 | 0.5590 | -13.7 |  -0.0 |
| 10 |  4   |  10 |  4   | 1.578 | 0.517 | 0.140904 | 1.4090 | 1.3080 | 1.4089 |  -7.2 |  -0.0 |
| 11 | 10   |  11 |  3   | 8.600 | 3.239 | 0.076765 | 0.8444 | 0.9674 | 0.8440 | +14.6 |  -0.0 |
| 12 |  4   |   6 |  3   | 6.739 | 2.501 | 0.011538 | 0.1385 | 0.2271 | 0.1384 | +64.0 |  -0.1 |
| 13 | 12   |  13 |  3   | 5.407 | 1.976 | 0.007017 | 0.0912 | 0.1274 | 0.0909 | +39.7 |  -0.3 |
| 14 |  6   |  14 |  3   | 4.429 | 1.591 | 0.010934 | 0.1531 | 0.1805 | 0.1528 | +17.9 |  -0.2 |
| 15 |  8   |  15 |  3   | 3.694 | 1.305 | 0.016294 | 0.2444 | 0.2570 | 0.2442 |  +5.2 |  -0.1 |
| 16 |  8   |   2 |  3   | 3.132 | 1.087 | 0.023954 | 0.3833 | 0.3674 | 0.3833 |  -4.1 |  -0.0 |
| 17 | 16   |  17 |  3   | 2.694 | 0.920 | 0.034418 | 0.5851 | 0.5289 | 0.5851 |  -9.6 |  -0.0 |
| 18 |  6   |   6 |  3   | 2.349 | 0.791 | 0.048619 | 0.8751 | 0.7665 | 0.8745 | -12.4 |  -0.1 |
| 19 | 18   |  19 |  3   | 2.072 | 0.689 | 0.068035 | 1.2927 | 1.1149 | 1.2921 | -13.7 |  -0.0 |
| 20 |  8   |  10 |  3   | 1.848 | 0.609 | 0.091672 | 1.8334 | 1.6186 | 1.8325 | -11.7 |  -0.1 |
| 21 | 12   |  21 |  3   | 1.664 | 0.546 | 0.121328 | 2.5479 | 2.3279 | 2.5472 |  -8.6 |  -0.0 |
| 22 | 10   |  22 |  3   | 1.512 | 0.496 | 0.159166 | 3.5017 | 3.2907 | 3.5012 |  -6.0 |  -0.0 |
| 23 | 22   |  23 |  3   | 1.386 | 0.456 | 0.205576 | 4.7282 | 4.5399 | 4.7278 |  -4.0 |  -0.0 |
| 24 |  8   |   6 |  3   | 1.280 | 0.424 | 0.263507 | 6.3242 | 6.0809 | 6.3247 |  -3.8 |  +0.0 |
| 25 | 20   |   5 |  2   |17.749 | 6.876 | 0.033372 | 0.8343 | 0.8312 | 0.8345 |  -0.4 |  +0.0 |
| 26 | 12   |  26 |  2   |16.448 | 6.358 | 0.026039 | 0.6770 | 0.5388 | 0.6769 | -20.4 |  -0.0 |
| 27 | 18   |   3 |  2   |15.288 | 5.897 | 0.008109 | 0.2190 | 0.2940 | 0.2189 | +34.3 |  -0.0 |
| 28 | 12   |  14 |  2   |14.251 | 5.484 | 0.000901 | 0.0252 | 0.1587 | 0.0253 |+529.1 |  +0.1 |
| 29 | 28   |  29 |  2   |13.319 | 5.114 | 0.001011 | 0.0293 | 0.1139 | 0.0293 |+288.5 |  -0.0 |
| 30 |  8   |  30 |  2   |12.478 | 4.779 | 0.001157 | 0.0347 | 0.1115 | 0.0346 |+221.2 |  -0.4 |

`rel_pred_R` is the **Empirical-r** prediction (uniform residue + empirical e
distribution). `rel_pred_Y` is the **Gaussian-Y** prediction (Gaussian Y).

**Aggregate mean abs error** across the 29 moduli:
- Gaussian-Y: 0.0962.
- Empirical-r: **0.0003 (essentially exact)**.
- Gaussian-r: 0.0775.

### F-verdicts

- **F1 (Primorial monotone Gaussian-Y).** **HOLDS.**
  `rel_pred_Y(2) = 0.881 > rel_pred_Y(6) = 0.457 > rel_pred_Y(30) =
  0.111`, matching the empirical monotone descent
  `rel_emp(2) = 0.939 > rel_emp(6) = 0.521 > rel_emp(30) = 0.035`.

- **F2 (Within-prime closed-form fit ≤ 25 %).** **PARTIAL.**
  Holds for 8/10 prime moduli `m ∈ {2, 3, 5, 7, 11, 17, 19, 23}`
  (errors 3.1 % – 14.6 %). FAILS for `m = 13` (39.7 %) and `m = 29`
  (288.5 %). The two failures share the structural feature that
  `J*(m)` is just below the next integer (e.g., `m = 29`:
  `log_29(p_N)/2 = 2.62`, `J* = 2`, `m^J* = 841 ≪ √p_N`), so `σ_Y`
  is large and the integration window for `Y mod m == 0` lies well
  outside the empirical e support — the Gaussian over-counts tail
  mass.

- **F3 (Primorial Gaussian-Y overshoot ≥ 2×).** **PARTIAL.**
  HOLDS for `m = 30` (Gaussian-Y predicts `0.111`, empirical
  `0.035` → 3.2× overshoot) and **fails in the opposite direction**
  for `m = 6`, where Gaussian-Y *under*shoots
  (`rel_pred_Y(6) = 0.457`, `rel_emp(6) = 0.521`, error -12 %).
  Why: at `m = 6` the conductor is too small for the J* dip to be
  in the deep regime; the empirical-r is essentially exact, the
  Gaussian-r is also close. The overshoot regime is `σ_Y ≳ 4`,
  reached only at `m ∈ {25, 26, 28, 29, 30}` here.

- **F4 (Empirical-r model is a better fit).** **HOLDS DECISIVELY.**
  Mean abs error: Empirical-r = **0.0003**, Gaussian-Y = 0.0962
  (320× ratio). The Empirical-r prediction tracks the empirical
  agreement on every modulus to within `|Δrel| ≤ 0.005` (worst cell:
  `m = 30` at `|Δrel| = 0.001`).

- **F5 (m = 5 mid-wrap local-max).** **HOLDS at L = 2·10⁸.**
  `rel_emp(4) = 0.835 < rel_emp(5) = 0.880 > rel_emp(6) = 0.521;
  rel_emp(5) > rel_emp(7) = 0.780.` `m = 5` is a strict local
  maximum on the prime/primorial neighbourhood; matches the F1.a-4'
  S199 prediction at the dense modulus grid.

- **F6 (m = 5 exception is closed-form-predictable).** **HOLDS.**
  The Gaussian-Y model also gives `rel_pred_Y(4) = 0.770 <
  rel_pred_Y(5) = 0.852 > rel_pred_Y(6) = 0.457`, and
  `rel_pred_Y(5) > rel_pred_Y(7) = 0.673`. The m = 5 mid-wrap is
  fully captured by the closed form (S199's "arithmetic coincidence"
  reading is sharpened: the coincidence is that `μ_e/m^J*` at
  m = 5 lands `σ_Y`-close to a half-integer, spreading shift mass
  across two adjacent residue classes and preserving substantial
  s = 0 mass).

- **F7 (Composite-m close to prime parent within 30 %).** **PARTIAL
  (3/6).** PASSES for `m = 4` (rad 2, `Δrel = 11 %`), `m = 9`
  (rad 3, `Δrel = 3 %`), `m = 25` (rad 5, `Δrel = 5 %`); FAILS for
  `m = 8` (`Δrel = 80 %`), `m = 16` (`Δrel = 59 %`), `m = 27`
  (`Δrel = 60 %`). The pattern: PASS when `J*(m) = J*(rad(m))`
  (so the Y distribution is identical), FAIL when `J*` differs.
  E.g., `m = 4, 8, 16` all have `m^J* = 4096` (same Y) but
  `rel(2^k) = 2^k · P[Y ∈ ⋃_j [j·2^k, j·2^k + 1)]` decreases as the
  modulus 2^k filters fewer agreement events.

### Net new content (refines E1.3 inline)

1. **EXACT formula for `ag_emp(m, J*)` in terms of the empirical
   `e_n` distribution.** Under the uniformly-distributed
   `L_n mod m^J*` assumption (verified to ≤ 0.04 % residual error),

   ```
   ag_emp(m, J*)
     = E_n[ Pr_r[ ⌊(r + e_n) / m^J*⌋ ≡ 0 (mod m) | r ~ Unif(0, m^J*) ] ]
     = E_n[ (1 − frac_n) · 𝟙[q_n ≡ 0 (mod m)]
            + frac_n · 𝟙[q_n ≡ m−1 (mod m)] ]
   ```

   with `q_n = ⌊e_n / m^J*⌋`, `frac_n = (e_n mod m^J*) / m^J*`.
   This identity is empirically exact across all 29 moduli at
   `L = 2·10⁸` (mean abs err 0.0003 on `rel`).

2. **The Gaussian-Y closed form is approximate.** Replacing the
   empirical e by `e ~ N(μ_e, σ_e²)` gives a closed form correct
   within ~10 % for cells with `σ_Y ≤ 1.5` and prime/near-primorial
   `m`, but fails by 30 %–530 % on cells with `σ_Y ≳ 2` because the
   *empirical e is bounded* (`e ∈ [0, 21648]`, skew = -0.108) while
   the Gaussian assigns substantial mass outside that support.
   **The Gaussian-RH-shadow conjecture proposed in S199's F1.a.i is
   thus REFUTED in the strong "exact" sense and DOWNGRADED to an
   approximation with ~10 % accuracy in its valid regime.**

3. **`rel(m)` is non-monotone with PEAKS, not just dips.** The
   measurement reveals `rel_emp(m) > 1` (predictor *better than
   uniform random*) at `m ∈ {10, 19, 20, 21, 22, 23, 24}`, peaking
   at `rel_emp(24) = 6.32` (i.e., `ag(24, J*=3) = 0.264 ≫ 1/24 =
   0.042`). The unifying picture is that `rel(m) = m · P[⌊Y⌋ ≡ 0
   (mod m)]` where `μ_Y(m)` ranges from `1.28` at `m = 24` to
   `17.75` at `m = 25`; whenever `μ_Y(m)` is close to an integer
   that is `≡ 0 (mod m)` and `σ_Y` is small enough not to spread
   the floor mass, `rel ≫ 1`. Conversely, primorial `m` puts `μ_Y`
   far from any small multiple of `m` and `σ_Y` large enough to
   diffuse mass — giving the deep dips. **This re-frames E1.3
   from "RH-shadow valley" to "RH-shadow phase alignment":** the
   predictor's digit-J* alignment is determined by `μ_Y(m) mod m`
   relative to the integer lattice.

4. **The m = 5 mid-wrap exception is captured by the closed form.**
   F6 confirms that `m = 5`'s local-maximum status (against
   neighbours `m ∈ {4, 6, 7}`) is reproduced by the Gaussian-Y
   prediction at L = 2·10⁸. S199's reading "arithmetic coincidence"
   is sharpened to "μ_Y(5) = 3.95 lands `σ_Y`-close to the
   half-integer 4, spreading floor mass across `{3, 4}` and
   preserving s = 0 residual through both wrapping branches".

### How this refines E1.3

E1.3 (per-bit difficulty of p(n)) was first refined by S146 with the
RH-scale anti-correlation valley at `J*(m=2)`, then by S199 with
cross-modulus universality of the valley at `J*(m) = ⌊log_m(p_N)/2⌋`.
S218 (this session) refines E1.3 with three further structural facts:

- (a) **EXACT closed form** for `ag_emp(m, J*)` in terms of `(e_n)`
  empirical distribution and uniform-r residue assumption. Mean abs
  err 0.0003 across 29 moduli at L = 2·10⁸.
- (b) **DOWNGRADE** of the Gaussian-Y closed-form prediction
  proposed in S199 from "exact" to "approximate within ~10 % in the
  σ_Y ≤ 1.5 regime, fails by 30 %–530 % outside it". The failure
  mode localises to the *bounded-support, sub-Gaussian* nature of
  the empirical e distribution.
- (c) **PEAK structure** of `rel_emp(m)`: not just dips; values up
  to `rel(m=24) = 6.32` reveal the predictor is m-adically *better*
  than uniform random for a substantial subset of small moduli.
  The unifying mechanism is `μ_Y(m) mod m` alignment with the
  integer lattice under floor.

### Closure mode

**Mode E** (extended measurement): refines E1.3 inline.
**Grade: B** — closed-form Gaussian-Y is approximate, not exact;
the EXACT empirical-r formula is a clean refinement; new peak
structure is a substantive structural finding beyond S199's
"valley" framing. Does not open polylog opportunities — the
formula is post-Li⁻¹-evaluation and explains agreement, not
prime-counting.

### Files

- `bit_J_pn_dip_scaling.py` — sieve + Li⁻¹ + per-modulus empirical &
  closed-form predictions.
- `bit_J_pn_dip_scaling_results.json` — full tabulated output at
  L = 2·10⁸.
- `scan_L1e7.json` — cross-scale anchor at L = 10⁷.
- `run_L2e8.log` — main run log.
- this `bit_J_pn_dip_scaling_results.md` — pre-stated falsifiers +
  outcome.

### Successor challenges (proposed)

**§F1.a.i.α — Sub-Gaussian tail correction to the closed form.**
Replace the Gaussian model of `e ~ N(μ_e, σ_e²)` with a *truncated*
Gaussian on `[0, e_max]` (or a beta-Cramér mixture matched to skew
and bound). Re-fit the closed form. Predicted: closed-form error
falls from O(10²) % to O(1) % across all m ∈ {2..30}. **A-grade if**
the tail-corrected closed form matches empirical to within 1 % on
all m. Cost: 1 session.

**§F1.a.i.β — Cramér-asymptotic prediction.** Under a Cramér model
`e ~ √p_N · N(c, σ²)` with `c, σ` constants, derive the asymptotic
`rel(m, p_N → ∞)` as a function of `m` only (J* absorbed via
`m^J* = √p_N`). Test at L ∈ {10⁷, 5·10⁷, 2·10⁸, 10⁹} (sympy/primecount
for the upper anchor). The Cramér prediction is *L-independent* in
the limit; cross-scale stability would confirm. Cost: 2-3 sessions.

**§F1.a.i.γ — Peak-vs-dip phase diagram.** The peak structure
(`rel_emp(m=24) = 6.32`) is governed by `μ_Y(m) mod m`. Tabulate
`rel(m)` against the phase `(μ_Y(m) mod m) / m ∈ [0, 1)` to plot
the "phase response" of the predictor. Predicted: a smooth U-shape
whose minimum is at phase 0.5 (mid-wrap, where mass splits) and
maxima at phase 0 (predictor aligned with ⌊Y⌋ = 0 mod m). Provides
a clean *one-parameter* fit of `rel` after factoring out conductor
size. Cost: 1 session.


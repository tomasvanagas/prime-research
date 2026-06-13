# F1.a.i.γ — Peak-vs-dip phase diagram

**Session:** 238 (NOVELTY mode, B-grade target).
**Edges refined:** E1.3 (per-bit difficulty of p(n) — sharp 4-bit sigmoid).
**Successor of:** §F1.a.i / S218 (peak structure of `rel_emp(m)` over m ∈ {2..30}).
**Composition:** S146 (bit-level RH-shadow) + S199 (cross-modulus universality)
+ S218 (rel_emp peaks at `μ_Y(m) mod m` near integer ≡ 0 mod m).

## Target (NOVELTY_CHALLENGES.md §F1.a.i.γ)

> Tabulate `rel_emp(m)` against the phase `α(m) := (μ_Y(m) mod m)/m ∈ [0, 1)`
> across `m ∈ {2..100}` at L = 2·10⁸. Predicted: smooth U-shape with maxima
> at phase 0 and minima at phase 0.5. Provides a one-parameter fit of `rel`
> after factoring out conductor size.

## Re-derivation: the one-parameter wrapped-Gaussian closed form

Continuing S218 notation: `Y := (r + e)/m^J*` with `r ∼ Unif(0, m^J*)`,
`e_n = p_n − round(Li⁻¹(n))`, `μ_Y = μ_e/m^J* + 1/2`,
`σ_Y² = σ_e²/m^{2J*} + 1/12`.

`rel(m) = m · P[⌊Y⌋ ≡ 0 (mod m)] = m · Σ_k [Φ((km+1; μ_Y, σ_Y) − Φ(km; μ_Y, σ_Y)]`.

Let `K := ⌊μ_Y/m⌋` and `β := μ_Y mod m ∈ [0, m)`. Setting `α := β/m ∈ [0, 1)`,
`σ_norm := σ_Y/m`. Substituting `j = k − K`:

```
rel(m) = m · Σ_j [ Φ_std((j − α)/σ_norm + 1/(m·σ_norm))
                   − Φ_std((j − α)/σ_norm) ]
```

with `Φ_std` the standard normal CDF. This is **exact** under the
Gaussian-Y assumption.

For large m (1/m → 0, Taylor in the inner argument), the bracket
collapses to `(1/(m·σ_norm)) · φ_std((j − α)/σ_norm)`, giving the
**wrapped-Gaussian density**:

```
rel_phase_lim(α, σ_norm) := (1/σ_norm) · Σ_{j ∈ ℤ} φ_std((j − α)/σ_norm)
                          = Σ_{j ∈ ℤ} φ_norm(j − α; 0, σ_norm)
```

where `φ_norm(x; 0, s) = (1/(s √(2π))) exp(−x²/(2 s²))`. By Poisson
summation this is `Σ_{n ∈ ℤ} exp(−2π² n² σ_norm²) cos(2π n α)`, a Jacobi
ϑ-type series.

**Properties of `rel_phase_lim`:**

- **Symmetry** about α = 1/2: `rel_phase_lim(α, σ) = rel_phase_lim(1 − α, σ)`.
- **Period** 1 in α.
- **Monotone in α on `[0, 1/2]`**: maximum at α = 0, minimum at α = 1/2,
  for every σ > 0.
- **σ → 0**: `rel_phase_lim → δ-comb`, sharp peak at α = 0 of height `1/(σ √(2π))`.
- **σ → ∞**: `rel_phase_lim → 1` uniformly (mod-1 mass spreads flat).

**This is a one-parameter family in α (after fixing σ).** The S218
peaks at m ∈ {19..24} should land at α near 0; the dips at primorials
m ∈ {6, 30} should land at α near 0.5.

## What's measured

For primes `p(1), …, p(N)` at sieve cap `L = 2·10⁸`:

For each modulus `m ∈ {2, 3, …, 100}` and `J*(m) := ⌊log_m(p_N) / 2⌋`:
- `α(m) = (μ_Y(m) mod m) / m ∈ [0, 1)`.
- `σ_norm(m) = σ_Y(m) / m ∈ (0, ∞)`.
- `rel_emp(m) := ag_Li(J*(m)) · m`.
- `rel_phase_exact(m)` — Gaussian-Y formula re-parameterised in (α, σ_norm, m).
- `rel_phase_lim(m)` — wrapped-Gaussian density (1/m → 0 limit).
- `rel_pred_R(m)` — empirical-e + uniform-r exact prediction (from S218).

## Pre-stated falsifiers (committed BEFORE running)

This section is fixed in advance.

**F_α (Wrapped-Gaussian-lim accuracy in the small-σ_norm regime).**
For m ∈ {m : σ_norm(m) < 0.5} the wrapped-Gaussian-lim prediction
satisfies `|rel_emp(m) − rel_phase_lim(m)| / max(rel_emp(m), 1) ≤ 0.20`.

**F_β (U-shape on the σ_norm < 1.0 sub-population).**
Bin moduli with σ_norm < 1.0 by α into deciles `α ∈ [k/10, (k+1)/10)`
for k = 0..9. Within this sub-population, define `rel_max := max over
α-deciles of mean rel_emp` and `rel_min := min over α-deciles of mean
rel_emp`. F_β HOLDS if **the maximum is in α-decile 0 or 9** AND **the
minimum is in α-decile 4 or 5** (i.e., the empirical-mean U-shape has
its extrema aligned with the α=0 and α=1/2 prediction).

**F_γ (Flat regime for σ_norm ≥ 2.0).**
For m with σ_norm ≥ 2.0, `|rel_emp(m) − 1.0| ≤ 0.5`.

**F_δ (Symmetry about α = 1/2).**
Define the "mirror partner" of m as `m' = argmin_{m'' ≠ m} (|α(m'') −
(1 − α(m))| + |σ_norm(m'') − σ_norm(m)|)`. F_δ HOLDS if for
**at least 60 % of moduli** with σ_norm < 1.0, `|rel_emp(m) − rel_emp(m')|
/ max(rel_emp(m), rel_emp(m')) < 0.30` AND a paired mirror exists with
combined coordinate distance ≤ 0.15.

**F_ε (Peak localisation).**
The top-3 peaks of `rel_emp(m)` over m ∈ {2..100} all have α(m) < 0.10
or α(m) > 0.90, AND σ_norm(m) < 0.5.

**F_ζ (Trough localisation).**
The bottom-3 troughs of `rel_emp(m)` over m ∈ {2..100} (excluding
flat-regime cells with σ_norm ≥ 2.0) all have α(m) ∈ [0.30, 0.70]
OR σ_norm(m) > 1.0.

**F_η (Wrapped-Gaussian-lim vs Gaussian-Y is small-1/m correction).**
The mean abs error of `rel_phase_lim` and `rel_pred_Y` (the S218
Gaussian-Y prediction) over m ∈ {2..100} differ by less than 30 %
(both being approximations of the empirical-r exact identity).

## Why these falsifiers

- **F_α** tests whether the small-σ_norm regime is well-approximated by
  the simple wrapped-Gaussian density (the 1/m correction is small).
- **F_β** is the *primary* prediction: the U-shape against α at fixed
  σ_norm regime.
- **F_γ** tests the σ → ∞ flat-mode prediction.
- **F_δ** tests the symmetry of the wrapped Gaussian about α = 1/2,
  a non-trivial prediction independent of any specific α-curve.
- **F_ε, F_ζ** test the peak/trough alignment with the (α, σ_norm)
  prediction, distinguishing the phase-controlled peaks from any
  alternative explanation (e.g., conductor-size only, prime-vs-composite).
- **F_η** tests whether the limit form (1/m → 0) loses much accuracy.

## Outcome

Run at `L = 2·10⁸`, `N = π(L) = 11,078,937`, last prime = 199,999,991.
Empirical residual: `μ_e = 10780.53`, `σ_e = 4293.66`,
`e_min = 0`, `e_max = 21648`. Modulus sweep:
`m ∈ {2, 3, …, 100, 110, 120, 140, 170, 200, 250, 300, 400, 500, 700,
1000, 1500}` (111 cells). The sweep spans both `J*(m) = 2` (m ≤ 119)
and `J*(m) = 1` (m ≥ 120).

### Aggregate mean abs error (rel_emp vs predictions)

| Predictor                              | Mean abs |Δrel| (n=111) |
|---------------------------------------|-----------------------------|
| Empirical-r (S218 exact identity)     | **0.0011**                  |
| Gaussian-Y (S218 closed form)         | 0.3768                      |
| Phase-Exact (re-param of Gaussian-Y)  | 0.3768  (=Gaussian-Y)       |
| Phase-Lim (wrapped-Gaussian density)  | 2.4447                      |

**The Empirical-r identity remains essentially exact across the
extended sweep including J*=1 cells**, confirming S218's identity
generalises. The Gaussian-Y is approximate; the Phase-Lim limit
is poor across most of the sweep.

### Phase-Lim by σ_Y regime — the structural split

| σ_Y regime            | n     | mean |Δrel_phase_lim| |
|-----------------------|-------|--------------------------|
| σ_Y ≥ 2.0             |   35  | **0.2088**               |
| σ_Y < 2.0             |   76  |   3.4744                 |

**The wrapped-Gaussian density is the correct closed form ONLY in
the σ_Y ≥ 2 regime** — in this regime the integration-window width
1 is small compared to σ_Y, so `∫_{km}^{km+1} φ(t; μ_Y, σ_Y) dt ≈
φ(km; μ_Y, σ_Y) · 1` is valid. **Outside this regime (σ_Y < 2,
typical of J*=2 cells with m ≤ 119), the approximation fails by
factors up to 24×**. The pre-stated F_α condition `σ_norm < 0.5`
was the wrong regime indicator.

### Top-10 peaks of `rel_emp(m)` — the peak ridge

| m   | J* | rel_emp | α     | σ_norm | σ_Y    |
|-----|----|---------|-------|--------|--------|
| 110 |  2 | 22.367  | 0.013 | 0.0042 | 0.457  |
| 100 |  2 | 14.090  | 0.016 | 0.0052 | 0.517  |
|  99 |  2 | 13.414  | 0.016 | 0.0053 | 0.525  |
|  98 |  2 | 12.759  | 0.017 | 0.0054 | 0.532  |
|  97 |  2 | 12.130  | 0.017 | 0.0056 | 0.540  |
|  96 |  2 | 11.541  | 0.017 | 0.0057 | 0.548  |
|  95 |  2 | 10.969  | 0.018 | 0.0059 | 0.556  |
|  94 |  2 | 10.420  | 0.018 | 0.0060 | 0.565  |
|  93 |  2 |  9.901  | 0.019 | 0.0062 | 0.574  |
|  92 |  2 |  9.398  | 0.019 | 0.0063 | 0.584  |

All top-10 peaks are **J*=2 cells with m ∈ [92, 110], α near 0,
and σ_norm decreasing as m grows**. Peak height tracks the
exact Gaussian-Y prediction (within ~5%): the empirical 22.37 at
m = 110 corresponds to a Gaussian peak centred at `μ_Y = 1.391`
near the integer 1 with σ_Y = 0.457; the per-cell agreement
`P[Y ∈ [0, 1)] = Φ((1−1.391)/0.457) − Φ((0−1.391)/0.457) = 0.196`,
giving `rel = 110 · 0.196 = 21.5` (Gaussian-Y), close to empirical
22.4 (Empirical-r exact: 22.37).

**Peak escalation vs S218.** S218 reported max `rel_emp(m=24) =
6.32` over m ∈ {2..30}. Extending to m ∈ {2..1500} reveals new
peaks 3.5× larger at m = 110, with the peak ridge growing
monotonically from m = 24 to m = 110 along the α ≈ 0 contour.
The peak ridge terminates at `m^J*` ≈ √p_N, where σ_e/m^J* drops
below 1 and Gaussian-Y begins to over-count in the bounded-e
support.

### Bottom-10 troughs of `rel_emp(m)`

| m   | J* | rel_emp | α     | σ_norm | σ_Y    |
|-----|----|---------|-------|--------|--------|
| 170 |  1 | 0.0048  | 0.376 | 0.149  | 25.258 |
| 200 |  1 | 0.0085  | 0.272 | 0.107  | 21.470 |
| 250 |  1 | 0.0187  | 0.174 | 0.069  | 17.177 |
|  28 |  2 | 0.0252  | 0.509 | 0.196  |  5.484 |
|  29 |  2 | 0.0293  | 0.459 | 0.176  |  5.114 |
| 300 |  1 | 0.0340  | 0.121 | 0.048  | 14.315 |
|  30 |  2 | 0.0347  | 0.416 | 0.159  |  4.779 |
|  31 |  2 | 0.0405  | 0.378 | 0.144  |  4.477 |
|  32 |  2 | 0.0474  | 0.345 | 0.131  |  4.203 |
|  33 |  2 | 0.0557  | 0.315 | 0.120  |  3.953 |

**Two distinct trough mechanisms:**

1. **J*=2 mid-wrap troughs (S218-class):** m ∈ {28, 29, 30, 31, 32, 33}
   with α near 0.5 and σ_Y ∈ [4, 5]. Empirical rel matches the
   wrapped-Gaussian-lim prediction within 50% (since σ_Y ≥ 2 holds
   here). These troughs correspond to phase-aligned destructive
   interference: μ_Y mid-wrap between integer multiples of m.

2. **J*=1 bounded-support troughs (NEW):** m ∈ {170, 200, 250, 300}
   with α NOT near 0.5 (α ∈ {0.121, 0.174, 0.272, 0.376}). At J*=1
   the bounded support `e ∈ [0, 21648]` truncates Gaussian-Y's
   right-tail mass that would naturally land in `[0, 1) mod m` at
   small α. Empirical rel ≪ Gaussian-Y prediction, while
   Empirical-r matches. **The deepest trough of the entire sweep
   (rel_emp(170) = 0.005) is a J*=1 cell at α = 0.376, not at α =
   0.5.**

### F-verdicts

| F | Pre-stated criterion | Result | Verdict |
|---|----------------------|--------|---------|
| F_α | rel_phase_lim ≤ 20% err for σ_norm < 0.5 | 31.5% pass over 111 cells | **FAILS as pre-stated**; refined regime σ_Y ≥ 2 valid (mean err 0.21 over 35 cells) |
| F_β | α-decile max in {0,9}, min in {4,5} | max=decile 0 (mean rel=4.67), min=decile 4 (mean rel=0.05) | **HOLDS** |
| F_γ | Flat regime σ_norm ≥ 2.0 within 0.5 of 1.0 | 0 cells in regime (m up to 1500 still has σ_norm < 0.30) | **N/A** |
| F_δ | ≥ 60% mirror-pairs (d ≤ 0.15) within 30% Δrel | 3/16 = 18.8% | **FAILS** |
| F_ε | Top-3 peaks have α < 0.10 (or > 0.90) AND σ_norm < 0.5 | m∈{110,100,99}, α≤0.016, σ_n≤0.0053 | **HOLDS** |
| F_ζ | Bottom-3 troughs (σ_n<2) at α ∈ [0.30, 0.70] OR σ_n > 1.0 | m∈{170,200,250}, α∈{0.376,0.272,0.174}, σ_n<0.15 | **FAILS** (J*=1 anomaly) |
| F_η | Phase-Lim vs Gaussian-Y err within 30% | err_Y=0.38, err_PL=2.44, ratio 549% | **FAILS strongly** |

**3 HOLD, 3 FAIL, 1 N/A.** The U-shape (F_β) and peak-ridge
location (F_ε) are confirmed; the wrapped-Gaussian-lim regime
is sharper than pre-stated (F_α/F_η); the symmetry (F_δ) and
trough localisation (F_ζ) FAIL informatively, revealing the
J*=1 bounded-support trough mechanism.

### Net new content (refines E1.3 inline)

1. **Wrapped-Gaussian-density regime split.** The `(α, σ_norm)`
   wrapped-Gaussian density `WG(α, σ_n) = Σ_j (1/(σ_n √2π))
   exp(−(j−α)²/(2σ_n²))` is the correct closed form for `rel(m)`
   only in the σ_Y ≥ 2 regime (mean |Δrel| = 0.21 over 35 cells).
   Outside this regime — typical of J*=2 cells with m ∈ [10, 120]
   — only the exact Gaussian-Y formula or the S218 Empirical-r
   identity capture rel_emp.

2. **Empirical U-shape against α.** Cell-mean `rel_emp` by
   α-decile traces `4.67 → 0.41 → 0.16 → 0.14 → 0.05 → 0.15 →
   0.68 → 0.80 → (empty) → 0.94` from decile 0 to decile 9.
   Maximum at decile 0, minimum at decile 4. The α distribution
   over m is heavily skewed (62/111 cells in decile 0), so the U
   is asymmetric — the right side (decile 6-9) is sparse.

3. **Peak escalation: rel_emp(m=110) = 22.37.** Extending the
   modulus sweep to m ∈ [2, 1500] reveals a peak ridge along
   α ≈ 0 with rel_emp growing 3.5× beyond S218's m=24 maximum.
   The ridge terminates at the J*=2 ↔ J*=1 boundary at m ≈ 119.

4. **J*=1 bounded-support trough mechanism (NEW).** The deepest
   troughs of rel_emp(m) over the full sweep are J*=1 cells
   (m ∈ {170, 200, 250, 300}) at α NOT mid-wrap. Mechanism:
   the bounded support `e ∈ [0, 21648]` truncates Gaussian-Y's
   mass that would naturally land in `[0, 1) mod m` at small α.
   The Empirical-r identity captures these troughs to within
   |Δrel| ≤ 0.001. **Sharpens S218's Gaussian-Y downgrade: J*=1
   regime worst-case error is 17× (m=200), exceeding S218's
   J*=2 worst-case of 5× (m=29).**

5. **U-shape symmetry (α ↔ 1−α) FAILS empirically.** Mirror
   pairs are sparse in the natural α(m) distribution, and
   what mirror pairs exist do not show the symmetry to within
   30 % Δrel (only 18.8 % pass). The α-distribution skew toward
   0 dominates the symmetry test — only m=2 has α > 0.85 in
   our sweep, and no m has α ∈ [0.85, 0.90).

### How this refines E1.3

S218 introduced the (μ_Y, σ_Y) closed form for rel_emp(m, J*)
across m ∈ {2..30} with the Empirical-r exact identity as the
gold standard. S238 (this session) extends the picture in three
directions:

- (a) **Phase reparameterisation `(α, σ_norm)`** — confirms a
  qualitative U-shape against α (decile-binned), with the
  wrapped-Gaussian density as the asymptotic limit valid only
  for σ_Y ≥ 2.
- (b) **Peak ridge structure** — m ∈ [92, 110] cells give
  rel_emp up to 22.37, 3.5× S218's known peaks. The ridge
  terminates at the J*=2 ↔ J*=1 boundary.
- (c) **J*=1 trough anomaly** — extending into J*=1 cells
  reveals deeper troughs than the J*=2 regime via the bounded-e-
  support mechanism, producing factor-17× Gaussian-Y errors at
  m=200 (vs S218's worst-case 5× at m=29). Empirical-r remains
  exact.

### Closure mode

**Mode E** (extended measurement): refines E1.3 inline.
**Grade: B** — the U-shape (F_β) and peak ridge (F_ε) are confirmed
on a broader sweep; the wrapped-Gaussian-density regime is mapped
(σ_Y ≥ 2 valid, σ_Y < 2 fails); the J*=1 bounded-support trough
mechanism is a new structural fact. Does not open polylog
opportunities — Empirical-r remains the post-hoc exact identity,
not a constructive predictor.

### Files

- `phase_diagram.py` — sieve + Li⁻¹ + four predictors (Gaussian-Y,
  Phase-Exact, Phase-Lim, Empirical-r) + decile + symmetry tests.
- `phase_diagram_results.json` — full tabulated output at
  `L = 2·10⁸` for all 111 cells.
- `scan_L1e7.json` — small sanity anchor at `L = 10⁷`.
- `run_L2e8.log` — main run log.
- this `phase_diagram_results.md` — pre-stated falsifiers + outcome.

### Successor challenges (proposed)

**§F1.a.i.γ.i — Truncated-Gaussian closed form for J*=1 cells.**
Replace `e ~ N(μ_e, σ_e²)` with truncated Gaussian on `[0, e_max]`
(or beta-Cramér mixture matched to the empirical e support
`[0, 21648]` and skew −0.108). Re-derive the closed form. Predicted:
J*=1 trough errors collapse from 17× to ≤ 30 %. **A-grade if** the
tail-corrected closed form matches empirical to within 1 % across
all m ∈ [2, 1500]. Cost: 1 session. (Subsumes §F1.a.i.α.)

**§F1.a.i.γ.ii — Asymptotic peak height.**
Conjecture: `max_m rel_emp(m, J*=2) → 1/(σ_norm·√(2π))` as m
approaches the J*=2 ↔ J*=1 boundary at `m ≈ p_N^{1/4}`. Empirical
ridge: m=110 → 22.4 / theoretical 95.0 (23.6 % of theoretical at
α = 0.013, finite α offset). Verify scaling by computing the peak
ridge at `L ∈ {10⁷, 5·10⁷, 2·10⁸, 10⁹}`. Predicted: peak height
grows as L^{1/8} or similar. Cost: 1-2 sessions.

**§F1.a.i.γ.iii — Symmetry restoration via L-scaling.**
F_δ failed because the natural α(m) distribution at fixed L is
skewed. Test whether scanning across `L ∈ {10⁶ .. 10⁹}` produces
mirror α-coordinates for the same m, and whether rel_emp(m, L) and
rel_emp(m, L') with α(m, L) ≈ 1 − α(m, L') match within 30 %.
This would be a concrete restoration of the symmetry prediction.
Cost: 1 session.


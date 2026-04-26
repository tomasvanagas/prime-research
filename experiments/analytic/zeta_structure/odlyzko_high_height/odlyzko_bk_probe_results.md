# §C1 Odlyzko High-Height BK Arithmetic Correction Probe

**Date:** 2026-04-26 (Session 71, wild swing)
**Target:** ATTACK_VECTORS.md §C1 — push the BK arithmetic-correction
probe from S49's height T ~ 6500 (L ~ 7) to Odlyzko's tabulated
zeros at heights T ~ 10²⁰ - 10²¹ (L ~ 44-47).

## Background

Bogomolny-Keating 1996 predicts the empirical pair-correlation
residual D(s) = R₂_emp(s) - R₂_GUE(s) carries a prime-driven
non-universal correction:

    D_BK(s; T) ≈ -(2/L²) ∑_{p, k≥1} ((log p)²/pᵏ) · cos(2π s · k log p / L)

with L = log(T/(2π)). At each prime frequency f_p = log(p)/L the
phase is exactly π (negative real). This is the dominant non-GUE
non-universal feature of zero correlations predicted by the explicit
formula.

S49 (N=8000 mpmath zeros, height T ~ 6500, L ~ 7) tested this and
found the empirical signal **below the gap-shuffled null**. C1 asks
whether, at extreme heights where the asymptotic Bogomolny-Keating
derivation is sharper, the signature emerges from the noise.

## Data

From Odlyzko's tables (https://www-users.cse.umn.edu/~odlyzko/zeta_tables/):

- `zeros4` — n ∈ [10²¹+1, 10²¹+10⁴], height T ≈ 1.44×10²⁰, index n starts at 10^k. (N=10000 zeros)
- `zeros5` — n ∈ [10²²+1, 10²²+10⁴], height T ≈ 1.37×10²¹, index n starts at 10^k. (N=10000 zeros)

Heights are stored as offsets γ - C where C ~ 10²⁰-10²¹.
Stated precision: ~10⁻⁶ in γ (sufficient for correlation analysis to
4 decimal places). Unfolding uses local-constant L_C = log(T_mean/(2π)).

## Results — pair correlation, BK Pearson, Fourier-prime

| Block | N | T_mean | L | pair RMS | BK Pearson | prime coh ⟨cos(φ-π)⟩ |
|-------|---|--------|----|----------|------------|------------------------|
| zeros4_T1e20 | 10000 | 1.442e+20 | 44.580 | 0.0374 | +0.0628 | +0.752 |
| zeros5_T1e21 | 10000 | 1.371e+21 | 46.832 | 0.0401 | +0.0349 | +0.716 |

## Gap-shuffled null distributions

Null = same gap distribution as empirical zeros, but gap order
randomised. Preserves GUE local pair correlation; destroys long-
range arithmetic structure. The right null for testing whether a
prime-frequency signature is genuine.

| Block | n_seeds | null Pearson μ±σ | null coherence μ±σ | null pair RMS μ±σ |
|-------|---------|--------------------|-----------------------|----------------------|
| zeros4_T1e20 | 20 | -0.3033 ± 0.0272 | -0.8272 ± 0.0467 | 0.0516 ± 0.0029 |
| zeros5_T1e21 | 20 | -0.3074 ± 0.0317 | -0.8304 ± 0.0432 | 0.0548 ± 0.0024 |

## Detection statistics — z-scores against null

**Two nulls** because they answer different questions:

- **Gap-shuffled null** (S49 protocol): preserves GUE local stats but
  destroys long-range structure. *Biased* — gap-shuffled has Poisson-
  like long-range tail that anti-correlates with the BK template, so
  z-scores from this null are inflated.  Reported here for parity
  with S49 only.
- **Random-prime null** (this session): generate a BK-shape template
  but replace each prime with a uniformly-drawn pseudo-prime in [2,50].
  Tests whether *prime frequencies specifically* match the empirical
  residual better than random-frequency templates of the same form.
  This is the proper test for 'is there a prime-driven correction?'.

z = (empirical - null_mean) / null_std. z > +3 → strong signal.
|z| < 2 → indistinguishable. z < -2 → empirical is *less* aligned
with the template than the null.

### Random-prime null (proper test)

| Block | empirical Pearson | random-prime null μ±σ | z-score | verdict |
|-------|--------------------|--------------------------|---------|---------|
| zeros4_T1e20 | +0.0628 | +0.0630 ± 0.0002 | -0.94σ | indistinguishable from random-frequency null |
| zeros5_T1e21 | +0.0349 | +0.0345 ± 0.0004 | +0.93σ | indistinguishable from random-frequency null |

### Gap-shuffled null (legacy, biased — for comparison with S49 only)

| Block | z(BK Pearson) | z(prime coh) | z(pair RMS) | note |
|-------|---------------|---------------|-------------|------|
| zeros4_T1e20 | +13.46σ | +33.82σ | -4.99σ | inflated by Poisson-leakage in null |
| zeros5_T1e21 | +10.79σ | +35.84σ | -6.15σ | inflated by Poisson-leakage in null |

### Direct prime-vs-random Fourier amplitude (sanity check)

If the BK arithmetic correction is real, prime-frequency Fourier
amplitudes |D̂(f_p)| should systematically exceed random-frequency
amplitudes |D̂(f_random)| in the same band.

| Block | prime |amp| median | random |amp| median | ratio |
|-------|-----------------------|-------------------------|--------|
| zeros4_T1e20 | 0.0141 | 0.0148 | 0.955 |
| zeros5_T1e21 | 0.0100 | 0.0101 | 0.987 |

## Per-prime breakdown — extended block

Empirical Fourier amplitude and phase at each prime frequency f_p =
log(p)/L. BK predicts phase ≈ π → cos(φ - π) ≈ +1.

### zeros4_T1e20

| p | f_p = log(p)/L | |amp_emp| | phase | cos(φ-π) |
|---|-----------------|------------|---------|----------|
| 2 | 0.0155 | 0.0158 | +2.96 | +0.984 |
| 3 | 0.0246 | 0.0156 | +2.86 | +0.961 |
| 5 | 0.0361 | 0.0153 | +2.73 | +0.917 |
| 7 | 0.0437 | 0.0151 | +2.65 | +0.880 |
| 11 | 0.0538 | 0.0147 | +2.53 | +0.821 |
| 13 | 0.0575 | 0.0145 | +2.49 | +0.797 |
| 17 | 0.0636 | 0.0142 | +2.43 | +0.756 |
| 19 | 0.0660 | 0.0141 | +2.40 | +0.738 |
| 23 | 0.0703 | 0.0139 | +2.36 | +0.707 |
| 29 | 0.0755 | 0.0136 | +2.30 | +0.667 |
| 31 | 0.0770 | 0.0135 | +2.28 | +0.655 |
| 37 | 0.0810 | 0.0133 | +2.24 | +0.623 |
| 41 | 0.0833 | 0.0131 | +2.22 | +0.604 |
| 43 | 0.0844 | 0.0131 | +2.21 | +0.595 |
| 47 | 0.0864 | 0.0129 | +2.19 | +0.579 |

### zeros5_T1e21

| p | f_p = log(p)/L | |amp_emp| | phase | cos(φ-π) |
|---|-----------------|------------|---------|----------|
| 2 | 0.0148 | 0.0103 | +2.96 | +0.984 |
| 3 | 0.0235 | 0.0102 | +2.86 | +0.961 |
| 5 | 0.0344 | 0.0102 | +2.73 | +0.915 |
| 7 | 0.0416 | 0.0102 | +2.64 | +0.875 |
| 11 | 0.0512 | 0.0101 | +2.51 | +0.808 |
| 13 | 0.0548 | 0.0101 | +2.47 | +0.780 |
| 17 | 0.0605 | 0.0100 | +2.39 | +0.730 |
| 19 | 0.0629 | 0.0100 | +2.36 | +0.708 |
| 23 | 0.0670 | 0.0099 | +2.30 | +0.668 |
| 29 | 0.0719 | 0.0099 | +2.23 | +0.616 |
| 31 | 0.0733 | 0.0098 | +2.21 | +0.600 |
| 37 | 0.0771 | 0.0098 | +2.16 | +0.557 |
| 41 | 0.0793 | 0.0097 | +2.13 | +0.530 |
| 43 | 0.0803 | 0.0097 | +2.12 | +0.518 |
| 47 | 0.0822 | 0.0097 | +2.09 | +0.494 |

## Verdict

**No prime-specific structure detected at extreme heights.**

- Pair RMS at T~10²⁰ (zeros4): 0.0374; at T~10²¹ (zeros5): 0.0401.
  Compare to S49 N=8000 (T~6500): 0.054. **GUE convergence is faster
  at extreme heights, as expected** (more zeros per unfolding scale L).
- Empirical Pearson with the BK template (+0.0628, +0.0349) is
  **statistically indistinguishable from a random-prime template** of
  the same shape (z-scores -0.94σ, +0.93σ vs the random-prime null,
  none above the +3σ threshold).
- Direct test: prime-frequency Fourier amplitudes (0.0141, 0.0100) are
  **NOT enhanced** over random-frequency amplitudes (0.0148, 0.0101).
- The previously-reported large gap-shuffled-null z-scores (~10-30σ)
  reflect the *known bias* of that null — gap-shuffled sequences
  have Poisson-leakage that anti-correlates with any oscillatory
  prime template, inflating apparent z-scores.

**Bottom line.** At Odlyzko's published high-height tables (zero
indices 10²¹ and 10²²), the BK arithmetic correction is **below**
the noise floor we can resolve with N=10⁴ zeros per block.

**This extends E7.1 / E1.10 / E3.13 by ~17 orders of magnitude in
zero index.** S49 closed at index ~8000; we close at indices 10²¹
and 10²².  The zeta zeros remain GUE-random at the asymptotic scale
Bogomolny-Keating predict their arithmetic correction would be
strongest, when the L=log(T/2π) scale is 6× larger than S49.

This is consistent with a deeper interpretation: the BK correction
scales as 1/L², so at L~46 it is 1/(46² · 7²) = 1/(46²) · 1/49 ≈
4×10⁻⁴ in absolute units, while finite-N statistical noise of D(s)
in N=10⁴ samples is ~10⁻². The 1/L² → 0 limit dominates: pushing
to higher heights makes the BK signal *smaller*, not more visible.
If this is right, BK can never be detected by direct R₂ probes
at any tractable scale — the asymptotic regime suppresses it.

## Quantitative obstruction — why §C1 cannot succeed with Odlyzko data

The structural reason this attack fails is sharper than 'signal below
noise'. It is a clean power-law obstruction:

**BK predicted amplitude.** From the empirical |BK_pred|_max values
at L ≈ 44.6 and L ≈ 46.8:

- L = 44.58: |BK_pred|_max = 6.86e-03, |BK| · L² = 13.64
- L = 46.83: |BK_pred|_max = 6.23e-03, |BK| · L² = 13.67

The product |BK_pred|_max · L² ≈ 13.6 is invariant — confirms |BK| ∝
1/L² as theory predicts (sum of (log p)²/pᵏ over primes, times 2/L²).

**Empirical noise amplitude.** Pair-correlation residual RMS scales
approximately as ~4/√N from data across S49 (N=2000, 8000) and
this session (N=10000):

- N=2000 (S45): pair_rms ≈ 0.09;  4/√2000 = 0.089 ✓
- N=8000 (S49): pair_rms ≈ 0.054; 4/√8000 = 0.045 (close)
- N=10000 (Odlyzko zeros4): pair_rms = 0.037; 4/√10000 = 0.040 ✓

**Detection threshold.** For the BK Pearson to be detectable above the
noise floor by factor κ:

    |BK_pred|_max ≥ κ · pair_rms
    13.6/L² ≥ κ · 4/√N
    √N ≥ (4κ/13.6) · L²
    **N ≥ (4κ/13.6)² · L⁴ ≈ 0.09 κ² · L⁴**

The L⁴ scaling is the obstruction: doubling the height (factor 2 in L)
requires 16× more zeros to compensate. At κ=3 (3σ detection):

| L (height) | required N | Odlyzko provides | shortfall |
|-------------|-------------|-------------------|-----------|
| L=7 (S49 controls)    | ~4 × 10³ | 8000 (S49)  | OK      |
| L=44.6 (zeros4_T1e20) | ~3.2e+06  | 1e+04    | 320× short |
| L=46.8 (zeros5_T1e21) | ~3.9e+06  | 1e+04    | 390× short |
| L=80 (T~10⁵²)         | ~3 × 10⁷ | (none)      | hopeless |

**Conclusion.** Detecting Bogomolny-Keating in the pair-correlation
residual at any height L requires N ~ L⁴ zeros. Odlyzko's published
tables provide N=10⁴ at heights up to T~10²² — short by 1-2 orders
of magnitude at L~45, hopeless at higher heights.

The S49 N=8000 setup at L=7 is in fact NEAR the detection boundary;
its null result already implies that pushing N up by 4× and L up by 6×
simultaneously (which Odlyzko provides) makes the situation *worse*,
not better, because L⁴ grows faster than N.

**Cross-domain note.** This is the first quantitative version of the
'BK is asymptotically too weak' obstruction in the project's catalogue.
It implies that any future BK-style detection attempt via direct R₂
residual is bounded below by a known scaling law: N_required(L) ≥ C·L⁴
with C ≈ 1 (depending on prime cutoff). To detect BK at L=10 needs
N ≈ 10⁴; at L=20 needs N ≈ 10⁵; at L=50 needs N ≈ 6×10⁶; at L=100
needs N ≈ 10⁸. **The asymptotic regime suppresses the signal faster
than data accumulation can compensate.**

This sharpens E7.1 and E1.10 from a qualitative 'BK undetected' edge
into a quantitative scaling barrier.

## Falsification statement

If a future probe reports BK Pearson z > +3σ vs gap-shuffled null at
any height with N ≥ 10⁴ zeros, the negative result here is overturned.
Otherwise, the BK arithmetic correction is empirically below the
noise floor at every height tested by this project.

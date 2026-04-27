# §D2.a.1 — Persistent homology of W-tricked normalised prime gaps vs continuous-marginal-matched baseline B3

**Status:** BUILT (S124, mode E, B-grade refinement of E2.17).
**Edge action:** refine E2.17 inline in EDGES.md. The E2.17 marginal-
distribution component (S117) decomposes further into a *marginal-
envelope* component (dominant, ~7–9σ on T0) and a *discreteness*
component (small but consistent, ~1–3σ on T0).
**Cross-domain technique:** persistent homology on Cramér-Takens
embedding (Carlsson 2009 BAMS 46; Bauer 2021 ripser) + empirical-CDF
inverse-transform sampling (Devroye 1986 *Non-Uniform Random Variate
Generation*, Springer §II.2.1).

## Question

NOVELTY_CHALLENGES.md §D2.a.1 (proposed in S117): the W-trick at
W = 210 (S117 / E2.17 refinement) erases the *serial-correlation*
component of the PH deficit on Cramér-normalised prime gaps —
z(P_W; B2)_T0/T1 collapses from |z| ≈ 7-12 (S96 unconditioned) down
to |z| ≤ 3 (B2 = gap-permutation, preserves the *exact* discrete
W-tricked marginal, kills serial structure).

But against B1 = IID Exp(1) the W-tricked window keeps a |z| ≈ 7-12
deficit on T0 (sign-flipping positive on T1 at d ∈ {3, 4}). S117
attributed this residual to the W-tricked gap MARGINAL departing from
Exp(1) (gaps quantised to multiples of W = 210 → discrete spectrum).

But *which part* of the marginal mismatch produces the deficit:

  (i) the *shape* of the marginal (different mean / variance / tail
      than Exp(1) — captured by any continuous distribution matched
      to the empirical CDF); or
  (ii) the *discreteness* (the marginal is supported on a quasi-grid
      of spacing ≈ 0.318 in Cramér-normalised coordinates)?

Construction B3 isolates (i): B3 = IID samples from a *continuous*
distribution matched to the W-tricked empirical marginal (linearly-
interpolated empirical CDF, equivalently inverse-transform sampling
from a piecewise-uniform density). Both the W-tricked sample and B3
share the continuous marginal envelope; the W-tricked sample has
discreteness, B3 does not.

## Outcome (one-line)

**F_a holds robustly on absolute thresholds; secondary discreteness
component identified.** Across d ∈ {2, 3, 4}, |z(P_W; B3)| ≤ 1.13 on
T0 AND T1 — the continuous-envelope baseline absorbs the entire
S117 residual to within Gaussian noise. The relative comparison
z(P_W; B3) vs z(P_W; B2) splits into a clear three-component
decomposition of the original E2.17 deficit (Section "Decomposition"
below). The pre-stated F_a relative threshold (|Δz(B3, B2)| ≤ 1) is
partially violated on T0 d ∈ {2, 3} — that violation is the new
content.

## Pre-stated falsifiers (BEFORE running) — outcome scored

The pre-stated comparison is on z(P_W; B3) and z(P_W; B2) pooled over
residues b ∈ {1, 11, 13}, M = 1000, d = 3, K = 20, x ≈ 10⁶:

| Anchor (S117 main run, M=1000, d=3, x≈10⁶, 3 residues pooled) | T0 | T1 |
|----------------------------------------------------------------|-----|-----|
| z(P_W; B1) — IID Exp(1)                                        | −9.07 | +5.56 |
| z(P_W; B2) — gap-permutation (discrete marginal preserved)     | −1.99 | −0.67 |

**This session (S124, M=1000, d=3, x=10⁶, 3 residues pooled):**

| Stat | z(B1)  | z(B2)  | **z(B3)**  | S117 z(B1) | S117 z(B2) |
|------|--------|--------|------------|------------|------------|
| T0   | −9.64  | −1.99  | **−0.05**  | −9.07      | −1.99      |
| T1   | +5.71  | −0.74  | **+0.46**  | +5.56      | −0.67      |

Cross-validation: B1 and B2 reproduce S117 to within 0.6σ on every
cell, confirming protocol stability across seeds (S124 seed differs
from S117).

### F_a — "marginal envelope absorbs the deficit" — **HOLDS partial**

|z(P_W; B3)| ≤ 3 on T0 AND T1 — **HOLDS robustly** across d ∈ {2,3,4}:

| d | T0 z(B1) | T0 z(B2) | **T0 z(B3)** | T1 z(B1) | T1 z(B2) | **T1 z(B3)** |
|---|----------|----------|---------------|----------|----------|---------------|
| 2 | −8.89    | −2.92    | **−0.03**     | −3.55    | −1.98    | **+0.65**     |
| 3 | −9.64    | −1.99    | **−0.05**     | +5.71    | −0.74    | **+0.46**     |
| 4 | −7.27    | −0.78    | **+0.42**     | +4.92    | −0.03    | **+0.37**     |

All six (d, summary) cells: |z(B3)| ≤ 0.65 (max 0.65 on d=2, T1).
*Absolute* condition fully holds.

|z(P_W; B3) − z(P_W; B2)| ≤ 1 — **partially violated on T0**:

| d | |Δz_T0(B3, B2)| | |Δz_T1(B3, B2)| |
|---|------------------|------------------|
| 2 | 2.89             | 2.63             |
| 3 | 1.94             | 1.20             |
| 4 | 1.20             | 0.40             |

The Δ shrinks with d (B2 std grows with d, narrowing the gap), but
on d ∈ {2, 3} the Δ exceeds 1.

### F_b — "discreteness is its own PH-detectable component" — **HOLDS partial (d=2)**

|z(B3) − z(B2)| > 2 on T0 OR T1 holds on d = 2 (Δ_T0 = 2.89,
Δ_T1 = 2.63). At d = 3, the T0 Δ = 1.94 falls just below the 2σ
threshold. At d = 4, Δ are well below 2σ. So F_b's pre-stated
strict criterion holds only at d = 2; on d = 3 it sits in the
F_a/F_b boundary; on d = 4 it dies.

The *direction* of the discreteness effect is consistent across d
on T0: B2 mean > B3 mean (i.e., the discrete-grid permutation has
*higher* T0 than continuous-envelope IID, so PRIMES_W lands closer
to B3 than to B2). On T1, B2 ≳ B3 (close).

### F_c — "envelope reverses the deficit" — **DOES NOT HOLD**

z(B3) does not exceed +3 anywhere (max +0.65). B3 does not
over-correct.

### F_d — "no signal at the chosen scale" — **DOES NOT HOLD**

z(B2) reaches |z| ≥ 1.99 on T0 d=3 and 2.92 on T0 d=2 — above
the 1.5 threshold. There is signal to detect.

**Pre-stated outcome scoring:** **F_a (absolute) HOLDS robustly;
F_a (relative) PARTIALLY FAILS on d ∈ {2, 3}; F_b (relative) HOLDS
on d = 2 only.** The combined reading is a three-way decomposition
of E2.17's marginal-distribution component — see "Decomposition"
below.

## Protocol

Same Takens-embed + ripser pipeline as
`persistent_homology_w_trick.py` (S117), but adds a third baseline:

1. Sieve primes up to N_max = 10⁷.
2. For each `b ∈ {1, 11, 13}` (gcd(b, 210) = 1), filter primes to
   `q_n ≡ b (mod 210)`.
3. Take M = 1000 consecutive primes ≡ b mod 210 starting at the first
   `q ≥ start_x = 10⁶`. Compute gaps `g_n = q_{n+1} - q_n`.
4. Cramér-normalise per residue:
   `x_n = g_n / (φ(W) log q_n)` with `φ(210) = 48`.
5. Construct three baselines (K = 20 fresh draws each, per residue):
   - `B1` IID Exp(1) of length M (continuous, Poisson) — S117 anchor.
   - `B2` random permutation of empirical W-tricked window (preserves
         discrete marginal exactly, kills serial correlation) — S117.
   - `B3` (NEW) inverse-transform sampling from the linearly-
         interpolated empirical CDF of `x`. Algorithm:
         `sorted_x = sort(x)`; for `u ~ Uniform(0,1)`, `pos = u*(n-1)`,
         `i = floor(pos)`, `frac = pos - i`, return
         `sorted_x[i] + frac*(sorted_x[i+1] - sorted_x[i])`.
         Continuous distribution matching the empirical marginal
         within linear interpolation. Devroye 1986 §II.2.1.
6. Takens-embed at d ∈ {2, 3, 4}, τ = 1, max_edge ∈ {4, 5}.
7. Compute `T0, T1, L0, L1, N1` for primes and each baseline.
8. Pool z-scores across the three residues.

## Empirical results (full table — d=3, M=1000, x≈10⁶)

### Per-residue means and stds (T0)

| b  | T0(P_W) | B1 mean ± std  | B2 mean ± std | B3 mean ± std |
|----|---------|----------------|---------------|---------------|
| 1  | 131.92  | 220.82 ± 8.34  | 133.83 ± 3.19 | 129.17 ± 7.90 |
| 11 | 130.72  | 220.87 ± 11.08 | 138.33 ± 2.27 | 136.26 ± 8.07 |
| 13 | 135.51  | 217.24 ± 8.07  | 143.13 ± 3.80 | 134.04 ± 7.59 |

### Per-residue means and stds (T1)

| b  | T1(P_W) | B1 mean ± std | B2 mean ± std | B3 mean ± std |
|----|---------|---------------|---------------|---------------|
| 1  | 33.47   | 25.58 ± 1.31  | 33.89 ± 1.35  | 32.33 ± 1.75  |
| 11 | 31.69   | 25.00 ± 1.68  | 34.27 ± 1.38  | 32.60 ± 2.33  |
| 13 | 32.92   | 24.92 ± 1.12  | 32.98 ± 1.16  | 31.03 ± 1.67  |

### Key structural observations

1. **B3 std is uniformly larger than B2 std.** On T0 across all
   d, b: B3 std ≈ 8 vs B2 std ≈ 3 (factor ~2.5×). Continuous
   resampling has wider variance than discrete permutation —
   resampling can land on any value in the marginal range, while
   permutation is constrained to the exact M=1000 empirical values.
   This larger std drives B3's |z| smaller than B2's |z| even at
   similar mean offsets.

2. **B2 mean > B3 mean on T0**, consistently across (d, b):
   - d=2: ΔT0(B2−B3) ∈ {3.08, 1.02, 5.46} → mean 3.19, ~3σ
   - d=3: ΔT0(B2−B3) ∈ {4.66, 2.07, 9.09} → mean 5.27, ~2σ
   - d=4: ΔT0(B2−B3) ∈ {5.65, 4.55, 14.12} → mean 8.11, ~1.7σ
   Direction: PRIMES_W has T0 closer to B3 (continuous envelope) than
   to B2 (discrete permutation).

3. **B1 mean is much higher than B2 / B3** on T0 (220 vs 133/130 at
   d=3) — this is the marginal-shape effect: the W-tricked marginal
   has variance ≈ 0.74² ≈ 0.55 vs Exp(1)'s variance 1, and is
   strictly bounded away from 0 (gap quantum ≈ 0.318). Both shift
   the mean cluster-merge time and total H_0 persistence downward
   relative to the Exp(1) cloud.

## Decomposition

E2.17 PH deficit on bare prime gaps (S96, M=2000, x≈10⁶, d=3):

  z(B1)_T0 = −10.31, z(B1)_T1 = −4.20
  z(B2)_T0 =  −8.70, z(B2)_T1 = −11.99

After the W = 210 W-trick (S117, M=1000, x≈10⁶, d=3, b ∈ {1,11,13}):

  z(B1)_T0 =  −9.07, z(B1)_T1 =  +5.56  (large; sign flip on T1)
  z(B2)_T0 =  −1.99, z(B2)_T1 =  −0.67  (small)

After adding B3 (S124, this session, same protocol):

  z(B3)_T0 =  −0.05, z(B3)_T1 =  +0.46  (essentially zero)

**Three-component decomposition (refines S117's two-component):**

  z(P_W; B1) ≈ z(P_W; B3) + Δ_envelope + Δ_serial_residual
                        ↑           ↑              ↑
                     ~ 0       ~7-9σ on T0     ~0-2σ on T0
                                                (small)

where:

  Δ_envelope         = z(P_W; B1) − z(P_W; B3)
                       ≈ −9 (T0)  /  +5 (T1)
                     = the marginal-shape difference between Exp(1)
                       and the W-tricked envelope (variance ≈ 0.55,
                       support on (0.318, ∞)).

  Δ_discreteness     = z(P_W; B3) − z(P_W; B2)
                       ≈ +1.9 (T0) / +1.2 (T1) at d=3
                     = the discreteness contribution: B2 (exact
                       discrete grid) has higher T0 than B3
                       (continuous envelope), so PRIMES_W (which
                       has slight log-drift off-grid) lands closer
                       to B3. Small but consistent: ~1.2-2.9σ on
                       T0 across d ∈ {2, 3, 4}.

  Δ_serial_residual  = z(P_W; B2) − z(P_W; ideal)
                       ≈ −2 (T0) / −0.7 (T1) at d=3
                     = the residual serial correlation in W-tricked
                       primes (S117 attributed this to the small
                       remaining HL k-tuple structure for primes
                       p > 7 within a single residue class).

**Refines E2.17 from two-component (S117) to three-component (S124):**

  total deficit (S96)  = serial-correlation (HL k-tuple, killed by W-trick)
                       + marginal-envelope (Exp(1) ↛ W-tricked CDF, dominant)
                       + discreteness (continuous envelope vs grid, small)
                       + residual serial-correlation (post-W-trick, small)

The dominant component is the marginal-envelope. Discreteness and
residual serial-correlation are sub-leading and partially cancel:
B2 has higher T0 than B3 by ~5 (discreteness raises T0); PRIMES_W
has lower T0 than B2 by ~5 (serial correlation lowers T0). They
roughly null on the (B3 vs PRIMES_W) comparison, which is why
z(B3) ≈ 0 even though both sub-components individually are ~1-3σ.

## What this is NOT

* **Not** a polylog opening for π(x) — same algorithmic ceiling as
  S117 (O(M²) distance + O(M³) PH).
* **Not** a *new* edge — refines E2.17's marginal-distribution
  component by splitting it into envelope + discreteness sub-
  components. CLOSED_PATHS row added (S124, mode E,
  refinement-of-E2.17).
* **Not** an attack on PH — the result confirms PH detects the
  HL singular-series structure via marginal shape; the discreteness
  signal is a PH-instrumentation artifact (cloud lattice → repeated
  distances), not arithmetic content.

## Cross-validation

The protocol re-runs B1 and B2 from scratch with a different seed
(S124 vs S117) and reproduces S117's pooled z-scores to within 0.6σ
on every (d, summary) cell:

| Cell | S124 z(B1) | S117 z(B1) | S124 z(B2) | S117 z(B2) |
|------|------------|------------|------------|------------|
| T0   | −9.64      | −9.07      | −1.99      | −1.99      |
| T1   | +5.71      | +5.56      | −0.74      | −0.67      |

Confirms the protocol is seed-stable and the new B3 z-score
(−0.05 / +0.46) is comparable to S117 anchors.

## Edges cited / composed

- **Composes** E2.13 (Gowers `U^2` W-trick uniformity) with
  E2.17 (PH deficit on gap sequence) and the inverse-transform-
  sampling primitive (Devroye 1986 §II.2.1).
- **Refines E2.17** with a three-component decomposition:
  marginal-envelope (dominant) + discreteness + serial residual.
- **Cross-references** E2.14 (Anderson Lyapunov), E2.15 (algebraic
  immunity), E2.16 (DPP failure), E2.20 (subword complexity) — all
  in the W-trick-erases-everything HL fingerprint family.

## Files

- `persistent_homology_w_trick_marginal_b3.py` — driver.
- `b3_pilot_b1.json` — pilot, b=1, M=500, K=10.
- `b3_main.json` — main pooled, b∈{1,11,13}, M=1000, d=3, K=20.
- `b3_d2.json`, `b3_d4.json` — d-scan robustness, M=1000, K=20.

## Cross-domain refs

- Carlsson 2009 "Topology and data" Bull. AMS 46(2), 255–308.
- Edelsbrunner & Harer 2010 *Computational Topology: An Introduction*.
- Bauer 2021 "Ripser" J. Appl. Comput. Topol. 5, 391–423.
  https://arxiv.org/abs/1908.02518
- Devroye 1986 *Non-Uniform Random Variate Generation* Springer
  (free PDF: http://luc.devroye.org/rnbookindex.html), Chapter II
  §II.2.1 (Inversion method).
- Green & Tao 2008 *Annals of Math.* 167(2) — origin of W-trick.
- Hardy & Littlewood 1923 — singular series.

## Successor challenges (proposed S124)

**§D2.a.1.i — Discreteness-only baseline B4: pure-grid IID.**
Replace B3's continuous envelope with B4 = IID samples from the
W-tricked DISCRETE marginal (sample from the empirical PMF: pick
each sorted_x[i] with probability 1/n, no interpolation). B4 is
IID on the exact discrete support; B2 is permutation. The
comparison z(P_W; B4) vs z(P_W; B2) isolates the IID-vs-permutation
serial structure on a fixed discrete marginal — should land within
0.5σ of B2 if S117's serial-correlation component is truly noise.
Cost: 1 session.

**§D2.a.1.ii — Sliding-bandwidth KDE B5(σ).** Vary the smoothing
bandwidth σ from 0 (= B2 discrete) to a large value (= near-Exp(1)
smooth) and trace z(P_W; B5(σ)) on T0 / T1. Predicted: a sigmoidal
crossover at σ ≈ grid-spacing / 2 ≈ 0.16, where the discreteness
component switches off. Cost: 1 session.

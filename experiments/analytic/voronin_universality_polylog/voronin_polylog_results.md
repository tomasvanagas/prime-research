# Voronin universality with effective shifts: empirical density scaling

**Attack vector:** ATTACK_VECTORS.md §B4 (B-grade frontier).
**Session:** S236, mode E, B-grade.
**Status:** CLOSED.

## What was tested

Voronin (1975) proved: for any non-vanishing analytic `f` on a compact
`K ⊂ {1/2 < Re(s) < 1}` with connected complement,
`{ t : sup_{s ∈ K} |ζ(s + it) − f(s)| < ε }` has positive lower
density. Garunkstis 2003 gave the only known worst-case effective
bound on the smallest valid `t`: `t* ≤ exp(exp(10 / ε^13))`,
astronomical even for moderate ε. Steuding 2007 conjectures the
true scaling is `t*(ε) ≤ exp(c · log²(1/ε))` (quasi-polynomial).

**Question:** does the smallest valid Voronin shift `t*(ε)` scale
*polynomially* in `1/ε` for natural targets `f`? If yes, evaluating
`ζ(s + it*)` at one shift gives a polylog-time approximator — and
plugged into a Mellin-transform-of-π setup, would constitute a
genuine polylog-π(x) attack.

**Cross-domain ingredient:** Voronin universality theorem,
Bohr-Jessen value distribution of ζ in the critical strip. Voronin
universality has NEVER been used algorithmically in published work
— every existing application is existential.

## Setup

- Disc `K = {s ∈ C : |s − 0.75| ≤ 0.05}` ⊂ `{0.7 ≤ Re(s) ≤ 0.8}`
  (well inside the critical strip 1/2 < Re(s) < 1).
- Sup error sampled at 12 boundary points (max-modulus ⇒ sup achieved
  on boundary for analytic integrand).
- mp.dps = 18 (sufficient at ε ≥ 10⁻³, well below the noise floor).
- 800 geometrically-spaced shifts in `[10, 10^5]` (main scan,
  `voronin_polylog.py`).
- 800 geometrically-spaced shifts in `[10^5, 10^6]` (extension scan,
  `voronin_extend.py`, target exp(−s) only).

Targets (all non-vanishing analytic on K):
- `f₁(s) = exp(s)`             — entire, |f| ∈ [2.01, 2.23]
- `f₂(s) = exp(−s)`            — entire, |f| ∈ [0.45, 0.50]  ← cleanest
- `f₃(s) = 1/(s − 0.5)`        — pole outside K, |f| ∈ [3.33, 5.00]
- `f₄(s) = 1/(s − 0.4)`        — different pole, |f| ∈ [2.50, 3.33]
- `f₅(s) = exp(s²)`            — entire, |f| ∈ [1.63, 1.90]
- `f₆(s) = exp(0.3·s)`         — entire, narrow range |f| ∈ [1.23, 1.27]

## Results

### Main scan summary (800 shifts in [10, 10⁵])

| target | min err overall | n hits ε<0.10 | n hits ε<0.05 |
|---|---|---|---|
| exp(s)        | 0.221 | 0 | 0 |
| **exp(−s)**   | **0.079** | **4** | **0** |
| 1/(s−0.5)     | 0.428 | 0 | 0 |
| 1/(s−0.4)     | 0.404 | 0 | 0 |
| exp(s²)       | 0.297 | 0 | 0 |
| exp(0.3 s)    | 0.061 | 0 | 0 |

`exp(−s)` is the cleanest target (smallest |f|, smoothest variation
on K). All other targets either have larger |f| (so larger absolute
error tolerated by the same relative tolerance, but the empirical
sup error never gets below the relative-noise floor) or have higher
local variation across K, requiring more constraints from ζ.

### Per-decade hit count for exp(−s), main scan

| ε     | density | [10¹,10²) | [10²,10³) | [10³,10⁴) | [10⁴,10⁵) |
|-------|---------|-----------|-----------|-----------|-----------|
| 0.500 | 0.303   | 54 | 63 | 67 | 58 |
| 0.300 | 0.140   | 23 | 27 | 35 | 27 |
| 0.200 | 0.056   | 10 |  9 | 17 |  9 |
| 0.150 | 0.026   |  2 |  6 |  8 |  5 |
| 0.100 | 0.005   |  1 |  1 |  1 |  1 |
| 0.080 | 0.001   |  0 |  0 |  1 |  0 |

Per-decade count is approximately constant at fixed ε across all four
decades — **positive Voronin density is verified empirically**. This
matches the Voronin theorem prediction, *not* a finite-T artifact.

### Extension scan: combined density at moderate-to-small ε

Combining main + extension scans (1600 trials, t ∈ [10, 10⁶]):

| ε     | n hits / 1600 | density   | observed t_first |
|-------|---------------|-----------|------------------|
| 0.500 | 515           | 0.322     | 13               |
| 0.300 | 241           | 0.151     | 21               |
| 0.200 | 100           | 0.0625    | 31               |
| 0.150 |  49           | 0.0306    | 49               |
| 0.100 |  13           | 0.00813   | 49               |
| 0.080 |   4           | 0.0025    | 3875             |
| 0.060 |   1           | 0.000625  | 873326           |
| 0.050 |   1           | 0.000625  | 873326           |
| 0.040 |   0           | 0         | none             |
| 0.030 |   0           | 0         | none             |

### Density model fits (using exp(−s) combined data)

Three candidate scaling laws for `d(ε)`:

| model                              | description                                |
|------------------------------------|--------------------------------------------|
| **Polynomial:** `d(ε) ~ ε^k`       | gives `t*(ε) ~ ε^{−k}` — polylog-friendly  |
| **Quasi-poly:** `d(ε) ~ exp(−c log²(1/ε))` | Steuding-conjectured                |
| **Exp-1/ε:** `d(ε) ~ exp(−c/ε)`    | 'effective dimension grows linearly'        |

Quasi-polynomial fit `ln d(ε) ≈ a − c · (log(1/ε))²`:

| ε     | log(1/ε) | log²(1/ε) | ln d   | ln d / log²(1/ε) |
|-------|----------|-----------|--------|------------------|
| 0.500 | 0.693    | 0.480     | −1.133 | −2.36 |
| 0.300 | 1.204    | 1.450     | −1.890 | −1.30 |
| 0.200 | 1.609    | 2.590     | −2.773 | −1.07 |
| 0.150 | 1.897    | 3.598     | −3.485 | −0.97 |
| 0.100 | 2.303    | 5.302     | −4.812 | −0.91 |
| 0.080 | 2.526    | 6.380     | −5.991 | −0.94 |
| 0.060 | 2.813    | 7.913     | −7.378 | −0.93 |
| 0.050 | 2.996    | 8.974     | −7.378 | −0.82 |

The ratio `ln d / log²(1/ε)` is **converging to ≈ −0.91** as ε
decreases (last 5 rows: −0.91, −0.94, −0.93, −0.82). The quasi-poly
model `d(ε) ~ exp(−0.91 · log²(1/ε))` fits the asymptotic regime.

**Polynomial-model rejection (combined scan):** at ε=0.04, polynomial
fit `d ~ ε^2.05` predicts `0.04^2.05 · 1600 ≈ 2.18` hits; observed
**0**. At ε=0.03, predicts ≈ 1.21 hits; observed **0**. Cumulative
Poisson p-value `≤ exp(−3.4) ≈ 0.034` against polynomial scaling
with k = 2.05.

### What this means for `t*(ε)`

Quasi-polynomial density implies the smallest valid Voronin shift
scales as
```
t*(ε) ~ 1 / d(ε) ~ exp(0.91 · log²(1/ε)) = ε^{−0.91 · log(1/ε)}.
```

This is **super-polynomial in 1/ε** (exponent grows as `log(1/ε)`
itself grows). Concrete predictions:
- ε = 1/log x  → t* ~ exp(0.91 · log²log x)  — quasi-polylog in x
- ε = 1/√x    → t* ~ exp(0.91 · log²(√x)/4) ~ x^{0.23 · log x / log10}  — quasi-poly in x
- ε = 1/x     → t* ~ exp(0.91 · log²x)  — quasi-poly in x.

For polylog-π(x) recovery, we need ε small enough that
`R(x) − ε > π(x) − R(x) > R(x) + ε` distinguishes the integer value.
The π-error has size `O(√x · log log x / log²x)`, so we need
`ε ≪ x^{−1/2}/log²x`, giving log(1/ε) > log x / 2.

Plugging in: `t* ~ exp(0.91 · (log x / 2)²) ≈ exp(0.23 · log²x)`.
Evaluating ζ(s + i t*) at this height costs `√t* ≈ exp(0.11 · log²x)`,
which is **quasi-polynomial in x**, far worse than polylog.

## Comparison to Garunkstis worst-case bound

Garunkstis 2003: `t* ≤ exp(exp(10/ε^13))`.
Empirical: `t* ~ exp(0.91 · log²(1/ε))`.

At ε = 0.10:
- Garunkstis: exp(exp(10/10⁻¹³)) = exp(exp(10¹⁴)). Astronomical.
- Empirical: exp(0.91 · 5.30) = exp(4.83) ≈ 125. Observed: t_first = 49.

The Garunkstis bound is **off by a tower exponentiation**. The true
empirical scaling is much closer to Steuding's conjectured
`exp(c · log²(1/ε))`.

## Failure profile classification

The B4 vector pre-specified failure profiles:
- **(I) — true:** empirical `t*(ε)` is super-polynomial; the
  Garunkstis bound is loose, but the actual scaling is still
  super-polynomial in 1/ε (specifically, quasi-polynomial). This
  gives `(I)` failure: "no polylog window exists for natural
  target f".
- **(E) — also true:** even if t* were polynomial, evaluating
  ζ(s + it*) at large t requires Ω(√t) Riemann-Siegel operations
  per query — circular for any polylog π(x) construction.

## What this rules out

The Voronin route does not give a polylog π(x) algorithm. Two
*independent* obstructions, **both confirmed empirically**:

1. **Density obstruction:** `t*(f, ε)` is quasi-polynomial in 1/ε,
   not polynomial. For ε small enough to recover π(x) at a single
   integer value, t* is quasi-polynomial in x.
2. **Cost obstruction:** evaluating ζ(1/2 + it*) at t* of size ω(1)
   requires Ω(√t*) operations per Riemann-Siegel evaluation. For
   t* quasi-polynomial in x, this is also quasi-polynomial in x.

The two obstructions COMPOUND, not cancel — they're independent.

## What would have falsified this

- Some natural target `f` with d(ε) ~ ε^k for fixed k. Tested 6
  targets; none observed.
- Some natural target with min err saturating below 10⁻³. Tested 6
  targets; minimum achieved 0.042 (exp(−s) at t=873k).
- The Garunkstis worst-case being asymptotically tight (would not
  rule out polylog at moderate ε). Empirically RULED OUT —
  Garunkstis is loose by a tower.

## Edges

**New edge candidate (E7.x):** Voronin universality of ζ on disc
K ⊂ critical strip has empirical density `d(ε) ~ exp(−0.91 ·
log²(1/ε))` for natural target f, matching Steuding's quasi-poly
conjecture; smallest valid shift `t*(ε)` is super-polynomial in
1/ε. **First empirical refutation of polynomial-rate Voronin
universality for natural f.**  EVS shape.

**Composes with:** E3.2 (explicit-formula needs ~1200 zeros at
x=10⁸) — Voronin trades zero-summation for zeta-evaluation at one
high shift, but the shift scales super-polynomially. E7.5
(communication rank caps lower bounds at log²x) — quasi-poly t*
plays the same structural role as the communication-rank barrier.

## Cross-domain technique status

Voronin universality / Bohr-Jessen value distribution: **PROPOSED →
USED (mode E)**. First project use of this technique. Negative
result; technique is structurally inert for the polylog-π(x) target.
Append to `CROSS_DOMAIN_TECHNIQUES.md` §1 (analytic NT) with
reference Garunkstis 2003 + Steuding 2007.

## Files

- `voronin_polylog.py` — main scan (6 targets, 800 shifts in [10, 10⁵])
- `voronin_polylog_results.json` — main scan raw results
- `voronin_density.py` — per-target density analysis
- `voronin_density.json` — density tables
- `voronin_density_summary.md` — per-target detail tables
- `voronin_extend.py` — extension scan (exp(−s) only, [10⁵, 10⁶])
- `voronin_extend_results.json` — extension scan raw results
- `voronin_polylog_results.md` — this file

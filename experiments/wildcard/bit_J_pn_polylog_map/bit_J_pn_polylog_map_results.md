# F1 — bit_J(p(n)) computability map (LSB-side complement to E1.3)

**Session:** 146 (NOVELTY mode, B-grade target).
**Edges refined:** E1.3 (per-bit difficulty of p(n) — sharp 4-bit sigmoid).

## Composition

- **E1.3** (`novel/carry_propagation_boundary.md`, S40) — per-bit
  agreement between p(n) and round(R^{-1}(n)) maps a sharp ~4-bit
  sigmoid at ~60% from MSB.

The S40 measurement uses agreement-with-R^{-1} as proxy for "easy /
hard". This conflates two distinct concepts:

1. **Predictor quality** — does R^{-1}(n) approximate p(n) at this bit?
2. **Intrinsic compressibility** — does bit_J(p(n)) admit any polylog
   predictor at all?

bit_0(p(n)) = 1 for n >= 2 by primality (the only even prime is 2).
This bit is **trivially polylog computable**, but R^{-1}(n)'s LSB is
essentially uniform mod 2, so S40's agreement rate is ~50% — and bit_0
falls into E1.3's "hard zone" by the predictor convention.

This refinement separates the two concepts by introducing alternative
predictors and direct intrinsic statistics.

## What's measured

For primes p(1), ..., p(N) with N = π(2·10^8) ≈ 10^7 (large N regime),
and for each bit position J = 0, 1, 2, ..., 27 (covering all bits of
p(n) up to 2.8·10^8):

(a) **Bias**: empirical mean of b_J(n) := bit_J(p(n)), centred at 1/2.
(b) **Peak discrepancy**: `max_M |sum_{n=1..M} b_J(n) - M/2|`, normalised
    by `sqrt(N)` (the GRH heuristic scale).
(c) **Lag-1 autocorrelation** of `b_J - 1/2`.
(d) **Agreement rate** with three predictors:
    - `P_const`: predict the empirically-majority bit (constant across n).
    - `P_PNT`: predict `bit_J( round(n · log n) )` (PNT first-order).
    - `P_Li`: predict `bit_J( round(Li^{-1}(n)) )` (PNT closed-form).

## Pre-stated falsification (BEFORE running)

This section is committed before any experiment runs.

**F1 (LSB triviality).** For J = 0, the constant predictor `P_const`
must achieve agreement >= 0.999 (with `P_const = 1` selected).

**F2 (LSB-side residue regime).** For J in {1, 2, 3}, all three
predictors yield agreement in [0.45, 0.55] (no predictor wins over
density baseline). The bit is "incompressible by simple closed-form
predictors" — consistent with reduction to polylog π(x; a, 2^{J+1}).

**F3 (intrinsic discrepancy bounded).** For J >= 1, peak discrepancy
normalised by sqrt(N) is <= 5 (under GRH heuristics for prime races,
discrepancy of `bit_J` cumulative deviation should be O(sqrt(N) log N)
for low J; we test the looser sqrt(N) constant of 5).

**F4 (Chebyshev sign on J=1).** For J = 1, the bias is *positive*
(majority of primes have bit_1 = 1, i.e., p ≡ 3 mod 4). Magnitude
modest; sign matches Chebyshev's classical race direction (more primes
≡ 3 mod 4 than ≡ 1 mod 4 in this range).

**F5 (R^{-1} predictability on top bits).** For J = 22..27 (top ~6
bits of N=10^7-th prime), `P_Li` agreement > 0.95.

**F6 (PNT vs Li equivalence at top).** For J = 22..27, both `P_PNT`
and `P_Li` exceed 0.95. (PNT first-order n·log n is essentially
equivalent to Li^{-1}(n) at the high-bit level.)

**F7 (lag-1 decoupling for high J).** For J >= 1, |lag-1 autocorr|
<= 0.05 (consecutive primes uncorrelated at any single bit position
beyond the LSB).

## How this refines E1.3

E1.3 reports a sharp 4-bit sigmoid at ~60% from MSB based on
`R^{-1}` agreement. This refinement decomposes the picture:

- The **MSB-side easy band** (top ~50–60% of bits) is **predictor-
  reducible to polylog Li^{-1}(n)** — F5/F6 confirm.
- The **LSB-side hard band** is NOT a single regime. It splits into:
    - **bit_0**: TRIVIAL (always 1; simple predictor wins).
    - **bits 1..k**: residue-class hard — reducible to polylog
      `π(x; a, 2^{J+1})`. F2 confirms no simple predictor wins.
- The **middle band** (between the trivial bit and the R^{-1}-easy
  band) inherits the residue-class structure for low J but is
  qualitatively the same as the LSB hard band.

In particular, E1.3's "hard zone" includes bit_0 by S40's predictor
convention, but bit_0 IS computable in O(1) — the conflation is in
S40's metric, not in E1.3's underlying difficulty claim.

## Files

- `bit_J_pn_polylog_map.py` — sieve + statistics + predictor agreement.
- `bit_J_pn_results.json` — full per-bit table (N, statistics).
- This results.md — pre-stated falsification + outcome below.

## Outcome

Headline finding (NEW): the `Li^{-1}` predictor agreement
`ag_Li(J) := P[bit_J(round(Li^{-1}(n))) = bit_J(p(n))]` exhibits
a **deep anti-correlation valley** at the bit position
`J* ≈ ⌊log_2(p(N))/2⌋` — exactly the RH error scale.

| limit `L` | π(L) | log_2(L) | dip position J* | ag_Li(J*) |
|-----------|-------|----------|------------------|-----------|
| 10^7      | 6.6·10^5 | 23.25 | J*=12 (J*/log_2(L) = 0.516) | 0.371 |
| 5·10^7    | 3.0·10^6 | 25.58 | J*=13 (0.508) | 0.360 |
| 2·10^8    | 1.1·10^7 | 27.58 | J*=14 (0.508) | 0.361 |

The dip position **shifts with N**, tracking ⌊0.5 · log_2(L)⌋ + 1.
Dip magnitude is **stable around 0.36** — the predictor is **anti-
correlated** with the truth, reaching ~64% disagreement.

The mechanism is the Skewes-direction sign bias of the rounding
error: under PNT-RH, `Li(x) > π(x)` empirically (sign flips occur
at Skewes-class numbers far above 10^8). Hence `Li^{-1}(n) < p(n)`
on this range. At bit position J equal to the error scale, the
rounding error's sign systematically flips the bit, producing
anti-correlation. At J << J*, error magnitude exceeds 2^J ⇒ both
bits effectively random, ag ≈ 0.5. At J >> J*, error << 2^J ⇒
agreement rises monotonically toward 1.

### F-verdicts

- **F1 (LSB triviality, J=0, ag(const) ≥ 0.999)** — **HOLDS** (ag = 1.000).
- **F2 (J ∈ {1,2,3}, all predictors in [0.45, 0.55])** — **HOLDS** for
  all three predictors at all three bits at L = 2·10^8 (ag = 0.500).
- **F3 (peak_disc / sqrt(N) ≤ 5 for J ≥ 1)** — **HOLDS for J ∈ {1..18}**;
  fails for J ≥ 19. Refined: F3 was correct only in the equiprobable
  regime; once bit_J becomes near-MSB the bit is monotone and peak_disc
  is O(N), not O(sqrt(N)). Restated F3': peak_disc/sqrt(N) ≤ 5 for J in
  the interior {1..⌊0.65·log_2(L)⌋}.
- **F4 (Chebyshev sign on J=1, bias > 0)** — **HOLDS** (bias = +0.0001
  at L=2·10^8, +0.0009 at L=10^6 — modest positive throughout, sign
  matches Chebyshev's race direction π(x;3,4) > π(x;1,4)).
- **F5 (J ∈ {22..27}, ag_Li > 0.95)** — **HOLDS** at L=2·10^8 (0.997 →
  1.000 across the range).
- **F6 (PNT vs Li equivalence at top)** — **REJECTED.** PNT first-order
  `n·log n` predictor is BIT-LEVEL inferior to `Li^{-1}(n)`. At top bit
  J=27 (L=2·10^8): ag_PNT = 0.927 vs ag_Li = 1.000. PNT does not
  approximate `p(n)` to bit-level accuracy at the top — full PNT
  closed form (Li^{-1}) is required.
- **F7 (lag-1 |autocorr| ≤ 0.05 for J ≥ 1)** — **REJECTED for low and
  high J.** At J=1, lag1 = -0.110; at J=2, -0.084; at J=3, -0.072
  (NEGATIVE — local Dirichlet equidistribution: consecutive primes
  tend to alternate residue classes mod 4/8/16). For J ≥ 5, lag1 is
  POSITIVE and large (consecutive primes are close in value, so they
  share high-bit values).

### Net new content

1. **bit_0(p(n)) is trivially polylog** (always 1 for n ≥ 2). E1.3's
   "agreement vs R^{-1}" classifies bit_0 as "hard" (ag = 0.500), but
   the bit IS polylog computable by the constant-1 predictor. E1.3's
   metric mis-classifies the LSB.

2. **Bit-level RH-scale anti-correlation valley** at J = ⌊log_2(p(N))/2⌋.
   `Li^{-1}(n)` rounding bit-flips at the RH error scale due to the
   Skewes-direction error-sign bias. This is the BIT-LEVEL SHADOW of
   `Li(x) > π(x)`. Magnitude stable across N: ag_Li ≈ 0.36, anti-
   correlation Δ ≈ 0.14.

3. **PNT first-order is bit-level insufficient even at top bits.** At
   the bit position of the leading prime power, PNT (`n log n`)
   predictor is only ~93% accurate — the next-order PNT terms are
   needed. This refines the conventional view that "PNT gets the top
   half of bits"; it does not get them at bit-level until Li^{-1}.

4. **Lag-1 autocorrelation of bit_J(p(n)) has a sign change at J ≈ 4**:
   negative for J ∈ {1, 2, 3, 4} (residue equidistribution),
   positive and growing for J ≥ 5 (consecutive primes are close).
   This is a clean empirical fact about the bit-position decoupling
   of consecutive prime indicators.

### Files

- `bit_J_pn_polylog_map.py` — sieve + statistics + predictor agreement.
- `bit_J_pn_results.json` — full per-bit table at L = 2·10^8 (28 bits).
- `scan_L*.json` — cross-scale scans at L ∈ {10^7, 5·10^7, 2·10^8}.
- `bit_J_pn_run.log` — main run log.

### Closure mode

**Mode E** (extended measurement). Refines E1.3 inline with three
items: (i) LSB triviality, (ii) RH-scale anti-correlation valley,
(iii) PNT-vs-Li bit-level inequivalence at top. Does NOT close any
new path or open any polylog opening — refinements stay in EDGES.md,
not CLOSED_PATHS.md.

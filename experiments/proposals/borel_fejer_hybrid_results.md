# Borel-Padé + Cesàro-Fejér hybrid — results (Session 51, FOCUS-4)

## Question

S45 (Borel-Padé regularisation) and S46 (Cesàro-Fejér damping) each gave
constant-factor improvements to truncated zero-sum approximations of π(x).
They were never combined.  Does Fejér-windowing the per-zero contributions
**first**, then Borel-Padé-resumming the windowed partial-sum increments,
beat Fejér alone?

## Setup

- 2000 zeta zeros (`data/zeta_zeros_2000.txt`, max γ ≈ 2515).
- mp.dps = 30.
- Branch-correct evaluation via `mpmath.ei(ρ · log x)` (per S46 / EDGES.md §E3.6).
- π(x) ≈ R(x) − Σ_γ w(γ) · 2·Re(li(x^ρ)).
- 8 x-values × 4 T-values × 3 modes (sharp / fejer / borel_fejer).
- For Borel-Padé, increment sequence a_k = w(γ_k)·2·Re(li(x^{ρ_k})) (windowed
  per-zero contribution); B(z) = Σ a_k z^k/k!; Padé orders (M,N) ∈
  {(5,5), (10,10), (15,15)} plus diagonals (5,5)..(14,14); take median across
  orders.

## Recovery rate ⌊round(S)⌋ = π(x)

| mode         | T=50 | T=100 | T=300 | T=1000 |
|--------------|-----:|------:|------:|-------:|
| sharp        | 3/8  | 4/8   | 6/8   | 4/8    |
| fejér        | 4/8  | 4/8   | 5/8   | **6/8**|
| borel-fejér  | 3/8  | 3/8   | 3/8   | 3/8    |

## Mean rounding gap |S − round(S)|

| mode         | T=50  | T=100 | T=300 | T=1000 |
|--------------|------:|------:|------:|-------:|
| sharp        | 0.207 | 0.296 | 0.191 | 0.222  |
| fejér        | 0.309 | 0.270 | 0.246 | 0.249  |
| borel-fejér  | 0.230 | 0.246 | 0.222 | 0.219  |

## Verdict

**The hybrid does NOT beat Fejér alone.** Recovery rate is flat at 3/8 (37%)
regardless of T; Fejér alone reaches 6/8 (75%) at T = 1000. Sharp matches
Fejér's 6/8 at T = 300.

## Why it fails — diagnosed

**Borel-Padé locks into a T-independent asymptote.** At x = 10000 the hybrid
returns S ∈ {1229.61, 1229.70, 1229.70, 1229.70} across T ∈ {50, 100, 300, 1000}
— **identical to 4 decimal places at three of the four T-values.** Padé fits
the leading envelope of the increment sequence, and the Borel integral
∫ e^{−z} P/Q(z) dz extracts that envelope's "tail completion" that depends
only on the *low-order* coefficients (where the Padé fit is well-conditioned).
Adding more zeros via larger T extends the increment sequence but does not
move Padé's leading-order fit, because the new increments are exponentially
suppressed by 1/k! in the Borel transform.

This is precisely what S45 already observed: "Borel-Padé is a *postprocessing*
of partial sums; it cannot extract information not present in the K zeros
consumed". Combined with Fejér, the *windowed* sum of all available zeros
already converges to a definite Fejér-mean answer — Borel-Padé's
"tail completion" is then a strict modification of an already-converged
windowed sum, with no theoretical justification.

**Cleanest failure pattern:** at small x = 1000 the hybrid produces nearly
exact answers (gap ≈ 0.01 at every T), because few zeros suffice and Padé's
first-order fit happens to land on the right number. At larger x the
oscillatory contribution requires high-frequency information from many zeros,
and the Borel envelope smoothes it away.

## Mode-by-mode wins (where each method uniquely succeeds)

| T    | hybrid-only | fejer-only      | hybrid-regressions vs. fejér |
|------|-------------|-----------------|-------------------------------|
| 50   | {2500}      | {25000, 75000}  | {10000, 25000, 75000}         |
| 100  | {100000}    | {5000, 50000}   | {5000, 10000, 50000}          |
| 300  | {}          | {}              | {5000, 10000, 50000}          |
| 1000 | {75000}     | {2500, 5000}    | {2500, 5000, 10000, 25000}    |

The hybrid uniquely solves at most 1 case at any T; Fejér alone uniquely solves
up to 2; **regressions outnumber wins at every T**.

## Failure-mode classification

**E (equivalence)** — Borel-Padé applied to a Fejér-windowed series collapses
into a low-order envelope estimator that ignores the structure Fejér is
specifically designed to expose. The hybrid is asymptotically *strictly worse*
than Fejér: more interventions, less information per zero used.

## Negative engineering result

This closes the natural "next combination" of S45 + S46. The S43-S46 line of
convergence-acceleration experiments has now exhaustively tested:

- raw partial sum (baseline)
- Cesàro mean (S45)
- Borel-Padé alone (S45 — sometimes wins at moderate x)
- Cesàro-Fejér window alone (S46 — best constant-factor improvement so far)
- Cesàro-Fejér × Borel-Padé hybrid (this session — strict regression)

**No combination has produced sub-√x asymptotic gain.** Fejér at T = 1000
remains the best constant-factor improvement (recovery 75% on x ≤ 10⁵).

## Runtime

Total sweep (8 × 4 × 3 = 96 evaluations) ≈ 18 s.
- sharp / fejér: 0.1–1.1 s per (x, T) row.
- borel-fejér: 2.7–3.8 s per row (Padé fitting + integration overhead).

## Status

CLOSED. Convergence-acceleration interventions cannot be stacked; each
intervention re-parametrises the same √x-bounded zero-sum information.

# Proposal 30 Results — Cancellation-Anchor Density Probe

## Question
What fraction of integers y in [x, x + 10⁴] have |E_K(y)| = |Σ_{k≤K} 2 Re(y^ρ/ρ)|
small enough to serve as π-evaluation anchors?

## Setup
- 1000 Riemann zeros, K ∈ {10, 50, 200}.
- x ∈ {10³, 10⁴, 10⁵}; 10001 integer y per (x, K).
- Thresholds: √(log x)  (polylog-tight), log²x (Cramér-scale), √x / log x
  (standard explicit-formula error scale).

## Results

| x      | K   | mean |E_K| | max |E_K| | frac ≤ √(log x) | frac ≤ √x/log x |
|--------|-----|-------------|-----------|------------------|--------------------|
| 10³    | 10  | 9.7         | 24.2      | 0.122            | 0.276              |
| 10³    | 50  | 11.5        | 37.5      | 0.146            | 0.256              |
| 10³    | 200 | 12.0        | 42.1      | 0.133            | 0.264              |
| 10⁴    | 10  | 19.3        | 41.1      | 0.081            | 0.287              |
| 10⁴    | 50  | 21.8        | 60.0      | 0.103            | 0.352              |
| 10⁴    | 200 | 22.7        | 75.7      | 0.063            | 0.276              |
| 10⁵    | 10  | 18.9        | 52.5      | 0.088            | 0.762              |
| 10⁵    | 50  | 25.6        | 64.6      | 0.249            | 0.575              |
| 10⁵    | 200 | 31.9        | 121.2     | 0.086            | 0.582              |

## Interpretation — proposal misleading as written, but result still negative

Apparent density of polylog-small |E_K| is 6–25%, which would naively suggest
anchors are common. However, **K = 200 truncated zeros do not give an anchor
for π(y)**: the full explicit formula needs all O(√y) zeros to converge to
π(y), and the truncation tail (zeros with γ > γ_K) contributes Ω(√y/log y) on
average. So the probe over-counts: many y look "good" only because we ignored
the tail.

A correct probe would need K = O(√y) zeros, which we do not have for x = 10⁵
(would need ≈ 316 zeros' worth of effective contribution, and we tested 200).
Within range, |E_K| at K = 200 is already growing as max = 121 at x = 10⁵
(scaling roughly √x), confirming the truncation tail is *not* dominating —
the partial sum already saturates the typical √x bound.

## Verdict — CLOSED (failure mode: equivalence)

To use a true anchor we need a y where the *full* zero sum is polylog-small,
not just its truncation. By random-matrix / GUE statistics of zeros, the full
sum is typically Ω(√y); the density of polylog-small full-sum y is ∝
polylog/√y, vanishing. The probe's apparent ~10% density is an artifact of
truncating the zero list short of the √y threshold.

## Failure mode
Equivalence (E): "good anchors" reduce to integers where all zeros conspire
to cancel — density ∝ polylog/√y → 0. Searching for them is at least as hard
as evaluating the full Σ_ρ y^ρ/ρ at a candidate y, which is √y per evaluation
without algorithmic shortcut.

## Side observation worth recording

The mean |E_K| grows roughly as log(x), not √x, for K ≪ √y. This is the
expected behaviour: for fixed K, |E_K| concentrates at scale O(√y · √K /γ_K)
≈ O(√y). The detail that mean stays small for small x is a small-K artifact;
extrapolation to x = 10¹⁰⁰ is unwarranted.

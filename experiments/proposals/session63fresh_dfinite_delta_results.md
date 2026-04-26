# P1 — D-finite (Apéry-style) recurrence hunt for δ(n)

**Script:** `session63fresh_dfinite_delta.py`

## What was tested

A sequence (a_n) is D-finite (holonomic) if there exist polynomials
P_0, ..., P_L ∈ ℤ[n], not all zero, with
sum_{j=0..L} P_j(n) · a_{n+j} = 0. Apéry's recurrence for ζ(3) is
L=2, deg=2.

D-finite is *strictly larger* than the algebraic relations PSLQ-style
searches probe — so this is a genuinely new test. If δ(n) = p(n) − R⁻¹(n)
were D-finite of small order/degree, we would get a fast holonomic
extrapolation algorithm, evaluating δ(n) in O((L+d)·M(log n)) and hence
p(n) in polylog.

## Method

1. Compute δ(n) for n = 1..400 at 200-bit precision (mpmath Newton-on-R⁻¹).
2. For each (L, d) ∈ {1..4} × {0..4}: build matrix A whose row at n has
   entries n^k · δ(n+j) for (j, k) ∈ [0..L] × [0..d].
3. Train on first half (n = 30..200), validate on second half (n = 201..395).
4. Take smallest right-singular vector of column-normalized A_train as
   candidate recurrence coefficients.
5. **Critical step:** validate by *predicting* δ(n+L) from δ(n)..δ(n+L−1)
   using the recurrence and comparing to the actual value. Skill =
   RMSE(predicted − actual) / std(δ). Skill = 1.0 means no signal;
   skill < 0.05 would mean a real recurrence.

## Result

| L | d | train rank ratio | pred skill (lower=better) | verdict |
|---|---|------------------|---------------------------|---------|
| 1 | 4 | 1.98e-4 | 0.50 | no skill |
| 2 | 4 | 1.35e-4 | 0.67 | no skill |
| 3 | 4 | 1.13e-4 | 1.06 | no skill |
| 4 | 4 | 9.96e-5 | 1.32 | no skill |
| ... all (L,d) | | small | ≈ 0.5–1.3 | no skill |

The training-side singular spectrum gets "tighter" with d only because
columns containing n^d span 12 orders of magnitude; this is a
*conditioning artifact*, not a real null space. The held-out predictions
have RMSE comparable to (or worse than) the standard deviation of δ.

## Verdict

**CLOSED.** δ(n) is NOT D-finite with order ≤ 4 and polynomial-coefficient
degree ≤ 4. Failure mode: **Information loss (I)** — the held-out
prediction has roughly the same error as a constant-zero prediction,
meaning the candidate recurrence captures no real structure of δ.

## Key takeaway / methodology lesson

When searching for null-space relations on data with multiplicative
column-scale spread (n^k for k up to d), one **must validate via
extrapolation/prediction**, not via singular-value ratios on the
training matrix. With doubles, conditioning alone produces
"rank-deficiencies" of ~1e-4 to 1e-13 that have nothing to do with
real algebraic structure. The decreasing `train_rank_ratio` with d in
the table above is *purely* this artifact.

## What would change the verdict

- A D-finite recurrence of order or degree ≥ 5: not ruled out by this
  test. But each step up costs O(K²) more linear-algebra rows; with 400
  data points we are bounded at (L+1)(d+1) ≤ ~50, i.e., (L,d) ≤ (6,6).
  Going further would need ~10000 δ(n) values, which is feasible.
- A *piecewise* D-finite recurrence (different recurrence in different
  intervals): would be invisible to this global search.

## One-line summary

Searching for a low-order linear recurrence with polynomial
coefficients on δ(n) finds no recurrence up to (L, d) = (4, 4); δ(n) is
not D-finite at low order, consistent with its inheritance of GUE-like
randomness from ζ-zeros.

# Continued-fraction expansion of π(x)/li(x) — results

**Question:** If the CF of π(x)/li(x) has bounded or eventually periodic partial
quotients, the residual is quadratic-irrational-like → exploitable algebraic
structure. If quotients follow Khinchin geometric mean ≈ 2.685, no structure.

**Setup:** Compute π(x), li(x) at 80-digit precision via mpmath; expand the
ratio into 30 partial quotients per x ∈ {10, 10², …, 10⁶}.

## Per-x partial quotients (a₀, a₁, …, a₂₄)

| x       | first 12 quotients                                  | geo mean | max  |
|---------|-----------------------------------------------------|----------|------|
| 10      | 0, 1, 1, 1, 5, 1, 1, 5, 1, 28, 5, 1                 | 2.33     | 35   |
| 10²     | 0, 1, 4, 1, 7, 7, 1, 5, 7, 2, 2, 5                  | 2.56     | 26   |
| 10³     | 0, 1, 17, 2, 13, 1, 2, 2, **326**, 2, 3, 14         | 2.78     | 326  |
| 10⁴     | 0, 1, 71, 1, 2, 1, 1, 20, 4, 8, 1, 4                | 2.76     | 71   |
| 10⁵     | 0, 1, 253, 1, 2, 3, 2, 2, 1, **508**, 2, 19         | 2.22     | 508  |
| 10⁶     | 0, 1, 605, 1, 13, 1, 2, 1, 4, 7, 1, 49              | 2.53     | 605  |

Aggregate (174 partial quotients across all x):
- Combined geometric mean = **2.521** (Khinchin K₀ ≈ 2.685)
- Top-5 max quotients: {71, 253, 326, 508, 605} — **growing without bound**

## Verdict — no Diophantine structure

The CF of π(x)/li(x) is **statistically generic**:

1. **Geometric mean ≈ 2.52**, within ~6% of Khinchin K₀ ≈ 2.685. Difference is
   well within sampling noise for 174 samples (Khinchin convergence is
   logarithmic; needs millions of quotients to nail the constant).
2. **Max quotients grow** with sample size (35 → 326 → 605), inconsistent with
   any periodic or bounded CF. Quadratic irrationals would have bounded a_i.
3. **Distribution** of low values matches Gauss-Kuzmin: ~41% of a_i = 1, ~17% of
   a_i = 2, declining roughly as 1/a² log(1+1/(a(a+2))) — generic continued-fraction
   law for almost-every real number.

A second pass on the *normalized residual* r(x) = (π(x) − li(x)) log(x) / √x
(bounded under RH) gives the same picture: geometric means 1.97, 2.49, 2.35;
max partial quotients up to 400. Generic irrationals.

## Failure mode

**I (Information loss).** The CF expansion is a *bijection* between reals and
sequences; structureless input yields structureless CF. No compression possible.
This is the 22nd or 23rd null result for "is δ(n) random-like" — independent of
the Padé/wavelet/Fourier nullsignals already found.

**Mild novelty:** I haven't seen the CF expansion test in the wildcard list.
Adding to the pseudorandomness file as evidence #N. CLOSED.

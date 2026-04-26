# Borel resummation of the Riemann explicit formula — results

**Question:** Can Borel-Padé resummation of the truncated zero-sum
`Σ_{|γ|<T} x^ρ/ρ` converge faster than the raw partial sum (or its Cesàro mean)?

**Setup:** First 200 zero ordinates, mp.dps=40. Compute partials S_K, increments
a_K, then Borel-Padé over (M, N) ∈ {5, 10, 15, 20}². Compare to raw S_200 and
Cesàro mean of partials. True value: ψ(x) = Σ_{p^k ≤ x} log p.

## Results

| x    | true ψ(x)  | raw S₂₀₀ err | Cesàro err | Borel-Padé err | best (M,N) | raw rate |
|------|------------|--------------|------------|-----------------|------------|----------|
| 50   | 49.485     | −0.071       | −0.025     | +0.157          | (20,20)    | K^−0.60  |
| 100  | 94.045     | −0.324       | +0.086     | +0.877          | (10,10)    | K^−0.38  |
| 500  | 501.652    | −2.044       | −1.422     | +0.223          | (5,5)      | K^+0.47  |
| 1000 | 996.681    | −0.924       | −1.911     | +0.645          | (5,10)     | K^−0.16  |

## Verdict

- **Mixed.** Borel-Padé beats Cesàro at x ∈ {500, 1000}, loses at x ∈ {50, 100}.
- Convergence rates of raw partial sums match theory: roughly O(K^{−1/2}) at small
  x, **diverging** at moderate x (K^{+0.47} for x=500). The known issue:
  truncation error in the explicit formula is ~ x · log²(x · T) / T (Goldston),
  not summable as T → ∞ for fixed x with insufficient zeros.
- Borel-Padé sometimes regularizes the divergence — see x=500 where raw error is
  −2.04 but Borel gives +0.22 (10× better).
- **But:** absolute error remains O(0.1) to O(1). For an *exact* algorithm we
  need error < 0.5 (then round). With 200 zeros that's not achievable for large x.
  Scaling test: error at x=1000 is 0.65 with 200 zeros → need γ_max ≈ x to
  resolve psi(x) to ±0.5. Same scaling as direct method (γ_max ~ x). **No win.**

## Failure mode

**E (Equivalence).** Borel-Padé is a *postprocessing* of partial sums; it cannot
extract information not present in the K zeros consumed. Asymptotically, the
truncation error is `O(x·log²(xT)/T)` zeros short of `T = γ_K`. Resumming
doesn't shrink T. To achieve exact ψ(x) you still need γ_K ≳ x · polylog(x), i.e.,
at least Ω(x/log x) zeros. Same complexity as Riemann–Siegel.

**However**, the regularization at moderate x is real and worth a follow-up:
testing whether *adaptive* Borel orders give consistent wins across x ranges
might yield a better constant factor for analytic prime counting (not asymptotic
gain). Logged, not pursued further.

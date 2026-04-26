# Multipole expansion of the zero sum — results

**Question:** Treat zeta zeros γ_k as charges on the line Re=1/2 and apply
fast-multipole-style cluster Taylor expansion. If kernel exp(i γ log x)/(1/2+iγ)
admits exponential convergence in cluster radius, we get O(√K · polylog) sum.

**Setup:** First K=200 zeros, x ∈ {100, 1000, 10000}. Cluster into B groups,
expand kernel to order P around each cluster center. Compare to exact partial sum.

## Headline numbers

x=10000, exact sum = −22.73:

| B  | P=2          | P=5          | P=10         | P=15         |
|----|--------------|--------------|--------------|--------------|
| 1  | 7.1e+07      | 5.8e+14      | 1.2e+27      | 2.9e+37      |
| 10 | 1.4e+06      | 1.0e+10      | 2.1e+19      | 1.6e+25      |
| 50 | 3.3e+04      | 2.6e+06      | 2.9e+13      | 1.2e+16      |

**Errors EXPLODE** with both B and P. Even at B=50 (cluster width 4 units of γ),
the expansion diverges.

## Convergence-in-P diagnostic at B=10, x=10000

| P  | err          | ratio |
|----|--------------|-------|
| 1  | 1.5e+02      | —     |
| 5  | 1.0e+10      | 8835  |
| 10 | 2.1e+19      | 617   |
| 20 | 3.9e+32      | 161   |
| 24 | 1.8e+37      | 111   |

Errors grow geometrically. Each term P contributes `(γ−γ_c)^P / P! · (i log x)^P`
times a moment — Stirling: `(R log x)^P / P!` peaks at `P* = R log x` then decays.

For the experiment: cluster width R ≈ (γ₂₀₀ − γ₁)/B ≈ 426/10 ≈ 43; log x = 9.2.
Peak at P* ≈ 396. Test only goes to P=24, way before turn-around → divergent.

## Why the analogy fails

In FMM (Greengard-Rokhlin) the kernel 1/(z − y) admits expansion
`Σ (y/z)^p` converging when |y| < |z|. **Geometric** convergence in p, ratio
|y/z|. Number of terms: P = O(log 1/ε).

The zeta-zero kernel `exp(i γ log x)/(½ + iγ)` is a **rotation** of magnitude 1,
not a decaying kernel. Its Taylor series in δ has coefficients
`(i log x)^p / p! · 1/(½ + iγ_c)`, which grow then shrink (Stirling tails).
To get error ε, need P > eR log(x) (Stirling bound) — *linear* in cluster width
times log x.

Total cost: B clusters × P terms ≥ K · log x even at optimum B → no asymptotic
gain. The kernel's frequency content (≈ γ_max log x ≈ K log x oscillations)
cannot be compressed below that count by any local-polynomial cluster scheme;
this is essentially the Nyquist limit for the sum.

## Verdict

**FAIL — Failure mode I (Information loss).** Multipole expansion requires kernel
*decay* in source-target separation. The zeta zero kernel is purely oscillatory
of unit amplitude, with frequency × log(x) determining oscillations across each
cluster. Compressing this is equivalent to band-limiting a non-band-limited
signal — same as the Information-Theoretic argument that the explicit formula's
zero contribution is incompressible.

What COULD work: a non-local basis (e.g., Hermite functions, prolate spheroidal
wave functions tuned to the γ density). But Hermite/PSWF bases lose the locality
that gives FMM its O(log 1/ε) terms per cluster.

CLOSED. Add to closed paths: "FMM-style multipole on explicit formula —
oscillatory kernel defeats local Taylor expansion."

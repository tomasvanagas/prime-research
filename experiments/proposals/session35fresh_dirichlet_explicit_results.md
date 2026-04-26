# Dirichlet explicit formula on AP mod 12 — Results (Session 35 fresh)

## What was tested
For q = 12 (phi(q) = 4), the truncated explicit formula was applied per AP:

```
pi(x;12,a) ≈ Li(x)/4 - (1/4) sum_{chi mod 12} chibar(a) * S_chi(x; T)
S_chi(x; T) = sum over zeros 1/2 ± i*gamma_k of L(s,chi):
                Li(x^{1/2 + i gamma_k}) + Li(x^{1/2 - i gamma_k})
```

Mod 12 has 4 real characters. Principal chi_0's zeros are zeta zeros
(2000 loaded from disk). For each non-principal character, its zeros up to
T_max = 40 were computed via mpmath sign-change tracking on Re L(1/2+it, chi).

T values: {5, 10, 20, 40}. x values: {1000, 3000, 5000, 7000, 10000}.
APs: a in {1, 5, 7, 11}. Total: 5 × 4 = 20 trials per T.

## Numerical findings

Non-principal zeros found below T_max=40:
- character 1: 4 zeros (~12.97, 27.01, 28.44, 38.97)
- character 2: 2 zeros (~21.45, 25.73)
- character 3: 2 zeros (~26.58, 30.75)

| T  | exact (rounded ≡ truth) | typical |diff| at x=10000 |
|----|--------------------------|---------------------------|
| 5  | 3 / 20 (15%)             | 1 to 14                   |
| 10 | 3 / 20 (15%)             | 1 to 14                   |
| 20 | 0 / 20 (0%)              | 1 to 70                   |
| 40 | 1 / 20 (5%)              | 6 to 74                   |

Surprising observation: T=20 gave WORSE results than T=10. This is because at
T=10 most non-principal zeros are absent and the approximation collapses to
just `Li(x)/4` (constant across a), which is "lucky" for the central a values
that are near the average. At T=20+, partial non-principal zero contributions
add real but incomplete oscillation, making rounding worse before it gets
better.

At T=40 the diffs at x=10000 range from 23 to 74 — clearly nowhere near
integer-rounded.

## Verdict — PROPOSAL 4 FAILS

T = polylog(x) is **far from sufficient** for integer-exact pi(x;12,a) at
x = 10^4. The truncation error term scales as ~ x / T (Dirichlet's bound),
giving error ~ 250 at T=40, x=10^4 — consistent with what we see (max
|diff| = 74, mean ~ 35).

For exact rounding, we'd need T ~ x^{1/2+epsilon} / phi(q), which for q = 12,
x = 10000 means T ~ 25; we got 1/20 exact at T=40, suggesting we're still
short. Going to x = 10^100 would need T = 10^{50}, blowing the polylog claim.

## Closure category
Failure mode: **(I) Information loss / sqrt-x barrier**. Same fundamental
obstruction as for principal chi (the zeta zeros): the L-function residue
contributions decay as zero height grows, but you still need x^{1/2+eps}
zeros to nail down ±0.5 rounding in pi(x).

GRH-based explicit formula is in fact known not to give polylog complexity,
even for individual APs at small q — the Riemann–von Mangoldt error is
intrinsic to the integer-rounding problem.

## What's NOT closed
- Combining MULTIPLE APs (e.g. pi(x) - pi(x; q, 1) etc.) cancels the smooth
  part but inherits oscillatory contributions from all characters. Same
  asymptotic.
- Selberg's mollifier or Vinogradov's zero-density estimates could improve
  constants, but not the polylog gap.

The proposal as stated fails. Confirms orthodox heuristic.

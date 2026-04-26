# B1 — Hierarchical-matrix compression of the explicit-formula kernel

**Date:** 2026-04-26
**Verdict:** CLOSED — structural negative.
**Mode:** I (information-theoretic / structural impossibility), N (no structure to exploit).
**Cited edges (post-hoc, not consulted at design time):** to be added if EDGES.md
matches; the result here is independent of any edge.

## What was tested

The Riemann explicit formula evaluated on a grid of x-points decomposes as
a matrix-vector product

    [psi_0(x_i)]_i  =  K v,
    K[i, k]  =  x_i^{1/2 + i gamma_k} / (1/2 + i gamma_k),
    v_k  =  -1.

If K admits hierarchical low-rank structure (HSS / HODLR / H-matrix), then
each new x evaluation costs only O(N r log N) after a one-time
factorisation, with r = O(polylog) the typical off-diagonal rank.

We measured off-diagonal ranks at successive dyadic block sizes for two
configurations:

- 64 x 64  (square, x in [10^3, 10^4], first 64 zeros)
- 256 x 32 (tall, more x-points than zeros)

## Result

Both rank-vs-block-width measurements give **rank exactly equal to
block width** with zero residual:

```
config 256 x 32:
  level 0 (w= 1)  rank_max = 1
  level 1 (w= 2)  rank_max = 2
  level 2 (w= 4)  rank_max = 4
  level 3 (w= 8)  rank_max = 8
  level 4 (w=16)  rank_max =16
  level 5 (w=32)  rank_max =32
  fit  rank ~ 1.00 * w + 0.0   residual = 0.000
  fit  rank ~ 5.74 * log2(w) - 3.86  residual = 11.24
```

The HODLR off-diagonal corner blocks show the same exact growth: every
m x n corner block has rank min(m, n).

## Why

The columns of K are pure characters

    K[:, k]  =  exp((1/2 + i gamma_k) log x) / (1/2 + i gamma_k).

For distinct gammas, these characters are linearly independent on any
generic x-grid — a Vandermonde-shaped system. Vandermondes have
**maximal** rank, not low rank. There is no off-diagonal cancellation
to exploit because the gammas are not clustered in a way that would
make nearby columns nearly proportional. (For nearby gammas the
columns *are* nearly proportional pointwise on any single x, but
across the whole x-grid the small differences in gamma sweep out
incommensurate frequencies.)

This is the reason H-matrices help for the Coulomb kernel
1/|x_i - y_k| (off-diagonal blocks really are low-rank for far-field
interactions) but **don't** help for Vandermonde / character matrices.

## What would falsify this

If anyone shows that for a *specific* zero-distribution model
(GUE or otherwise) the kernel admits epsilon-rank `r = O(polylog)`
in a regime relevant to pi(x) (i.e. with N ~ sqrt(x) zeros and
x ~ 10^k), then H-matrix compression becomes viable. Our test was
at small N = 32-64; perhaps very high-density zero regimes show
different behaviour. But the structural argument suggests this is
unlikely: even at infinite zero density the columns are still
independent characters of x.

## Implication

The H-matrix angle (kernel-side compression) joins the family of
"Vandermonde-like" matrix structures that resist hierarchical
compression. To get polylog cost we'd need *signal-side* sparsity (i.e.
v sparse in some basis, which is the compressed-sensing direction)
— not kernel-side compression.

## Cost

~5 seconds wall time (dominated by mpmath zeta-zero computation).

## Code

`hmatrix_kernel_compression.py` runs both configs above. JSON data
lives in `hmatrix_kernel_compression_data.json`.

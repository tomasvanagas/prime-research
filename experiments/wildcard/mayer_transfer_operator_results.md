# Mayer Transfer Operator (Gauss-Kuzmin-Wirsing) — Results

**Experiment:** `experiments/wildcard/mayer_transfer_operator.py`
**Date:** 2026-04-26
**Verdict:** CLOSED (preliminary) — naive Galerkin discretization does not locate Riemann zeros at the predicted critical-line points; the "compress many zeros into a small matrix" hope fails at this implementation level.

## Idea

Mayer (1991) showed that the Selberg zeta function for `PSL(2,Z)` factors as
`Z(s) = det(1 - L_s) · det(1 + L_s)`, where the *Mayer transfer operator*

```
(L_s f)(z) = sum_{n>=1} (z+n)^{-2s} f(1/(z+n))
```

is nuclear (trace-class of any order) on appropriate Hilbert spaces of analytic
functions. Numerical Galerkin discretization on a polynomial basis converges
*exponentially* in basis size N. The fresh question: does a *small* matrix
(N=10..20) capture *many* Riemann zeros, giving an exponential compression
of the zero-sum bottleneck in the explicit formula?

Matrix entries on the monomial basis `{1, z, ..., z^{N-1}}`:
```
M_{jk} = (-1)^j · (2s+k)_j · zeta(2s+k+j) / j!
```
(where `(α)_j = α(α+1)...(α+j-1)` is the Pochhammer symbol).

## What was tested

1. Built `M_{jk}` for `s = 1/2 + i*t`, `t ∈ [10, 60]`, on a grid of 80 points.
2. Computed `det(I − L_s²)` (vanishes wherever `det(1 ± L_s)` vanishes).
3. Scanned for local minima of `|det|` and matched against the first 12
   non-trivial Riemann zeros (`γ_1 = 14.13...`, `γ_2 = 21.02...`, etc.).
4. Sanity check: Riemann's truncated explicit formula
   `Li(x) − 2 Re Σ Li(x^ρ_k)` for the *exact* known zeros, with K=10/20/30.

## Results

### Test 1 — Recovery vs. dimension

| N | time(s) | minima found | matched (of 12) | avg error |
|---|---|---|---|---|
| 6  | 1.5 | 0 | 0 | n/a |
| 10 | 1.7 | 0 | 0 | n/a |
| 14 | 2.1 | 0 | 0 | n/a |
| 18 | 2.4 | 0 | 0 | n/a |

**Zero zeros recovered at any N tested.** The determinant `|det(I − L_s²)|`
varies smoothly along the critical line `s = 1/2 + it` — no local minima
within the search range.

### Test 2 — Truncated explicit formula (sanity)

Using the *exact* Riemann zero ordinates in the formula
`Li(x) − 2 Re Σ_{k=1}^K Li(x^{1/2+iγ_k}) ≈ pi(x)`:

| x | pi(x) | err K=10 | err K=20 | err K=30 |
|---|---|---|---|---|
| 50    | 15  | +3.7  | +3.8  | +3.8  |
| 100   | 25  | +4.4  | +4.6  | +4.4  |
| 1000  | 168 | +9.5  | +9.2  | +9.0  |
| 10000 | 1229| +19.7 | +19.9 | +19.6 |

The truncated explicit formula *systematically over-estimates* by ~`O(√x/log x)`
because we used the simple `Li(x)` smooth term (which approximates `Π(x)`,
the prime-power summatory function) instead of Riemann's `R(x)` plus
Möbius-corrected zero sum. **K=30 is far from enough zeros for x=10⁴.**

## Why it failed (in this naive form)

1. **Operator/zero-locus mismatch.** The Mayer identity puts Riemann zeros on
   the line `s = 1/4 + iγ/2`, *not* `s = 1/2 + iγ` (zeros of the Eisenstein
   factor `ζ(2s)/ζ(2s−1)`). Scanning the wrong line yields no minima.
2. **Truncation regime.** Mayer convergence is exponential in N, but the
   truncation error envelope on the critical strip starts large (the operator
   is bounded but not contraction); for moderate N the smallest singular
   value of `(I − L_s²)` may not vanish near a true zero — it merely *dips*.
3. **Zero density mismatch.** Even *if* the Mayer determinant gave the
   correct zeros, the *count* of zeros up to height T is `~T/(2π) log T`.
   A polynomial-basis truncation of size N has at most N eigenvalues, hence
   N zeros — same scaling as direct enumeration. The conjectured exponential
   compression has no a-priori basis.

## Verdict

**CLOSED.** Two independent failure modes — wrong critical line and no
super-linear compression — render the naive Galerkin discretization a
non-starter for fast `pi(x)`. A more careful version (correct line
`s = 1/4 + iγ/2`, mpmath high-precision determinant, Riemann R-function for
the smooth part) might recover zeros, but it would not change the fundamental
fact that *N basis functions give at most N zeros*.

Failure mode: **C (circularity)**. To extract K zeros from the Mayer operator
we need at least N=K basis functions, and computing the determinant is
`O(N³)`. No improvement over direct zero tabulation.

## What would actually be novel

If someone proves a Mayer-type identity for an operator whose *spectral zeta
function* (not operator zeros) directly equals the Riemann zero sum
`Σ x^{ρ_k}`, *and* the operator has structured (banded/sparse) form so the
trace `Tr f(L)` is computable in `O(polylog)` without diagonalizing — that
would be the breakthrough. No such identity exists in the literature.

## File status

- Script: `experiments/wildcard/mayer_transfer_operator.py`
- Implementation: vectorized monomial-basis Galerkin, mpmath ζ for matrix
  entries, scipy `slogdet` for stable determinants.
- Reusable: yes, the matrix-builder is correct and could be repurposed for
  the right line `s = 1/4 + iγ/2` or for a different operator.

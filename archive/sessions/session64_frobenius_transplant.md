# Session 64 — FOCUS-1 sub-attack 3: Healy-Viola Frobenius transplant

**Date:** 2026-04-26
**Mode:** normal (FOCUS-1, construction)
**Run state:** 41 → 42
**TODO item:** FOCUS-1 sub-attack 3 (Healy-Viola Frobenius transplant)

## Summary

Built and probed the AKS-mod-q congruence

```
   (a + x)^n  ≡  a^n + x^n   (mod  q,  x^r - 1)
```

with q a prime dividing `n - 1`, exploiting the q-power Frobenius
`F: F_q[x]/(x^r-1) → F_q[x]/(x^r-1), F(a+x) = a + x^q`. This Frobenius
gives a binomial decomposition

```
   (a + x)^n  =  prod_i (a + x^{q^i mod r})^{c_i}    [in F_q[x]/(x^r - 1)]
```

where `c_i` are the base-q digits of n. The aim was to determine whether
this restructuring (i) gives a TC^0 depth collapse and (ii) lets AKS
correctness "lift back" from F_q to Z_n.

## Headline result

**Closes FOCUS-1 sub-attack 3 as FAIL (E + I).** Two simultaneous
failure modes block the construction:

1. **(I) Information loss.** Primes do NOT satisfy mod-q AKS for any
   non-trivial a. Empirical: 0/399 passes across 19 primes in
   [101..2207], every q | n-1 with q < r, every a in {1, ..., q-1}.
2. **(E) Equivalence.** The mod-q computation is length-O(log n)
   growing-dim r×r MPOW over F_q — exactly the open Chain-E primitive
   E5.3, unchanged by the q-twist of the coefficient ring.

## Concrete numbers

| Quantity | Result |
|----------|--------|
| Frobenius-decomp correctness (Q1) | 6/6 cases match naive repeated squaring |
| Primes pass mod-q AKS, non-trivial a (Q2) | 0 / 399 |
| Primes pass mod-q AKS, trivial a≡0 mod q | 48 / 48 (both sides = x^n; zero info) |
| Residual poly nonzero coefficients (Q3) | 8/53 to 102/107 — fills the ring, no sparsity |
| Frobenius vs naive multiplications (Q4) | factor `log(q)/log(2)` saving (≈2× for q=2) |

## Why the AKS-mod-q test fails for primes

For prime n, the polynomial identity `(a+x)^n - (a^n + x^n) =
sum_{k=1}^{n-1} binom(n,k) a^{n-k} x^k` vanishes mod n (Frobenius is
the n-power) but NOT mod q for `q ≠ n` prime. By Lucas's theorem,
`binom(n, k) mod q = prod_i binom(n_i, k_i)` where `n_i, k_i` are base-q
digits — non-zero whenever k is a "base-q sub-string" of n. Concrete:
n = 101 = 0b1100101, so `binom(101, k)` is odd iff k's bits are a
subset of 101's bits — 16 of 102 binomial coefficients are odd. All
contribute non-zero terms to (1+x)^101 mod (2, x^5-1).

This means the natural mod-q reduction is a totally different
polynomial in F_q[x]/(x^r-1) than `a^n + x^n`. Lifting back to Z_n
would require Hensel-type reconstruction of these binomial coefficients
from their mod-q residues — but the binomial coefficients are
themselves the obstruction (their pattern mod q is precisely Lucas's
sub-string structure of n's base-q expansion, which carries info about
the digits of n, not about n's primality).

## Why the test, even if it worked, would not be TC^0

* Frobenius decomposition has `log_q(n)` factors (each binomial).
* Each factor is multiplied into a running product in F_q[x]/(x^r-1).
* Each multiplication is a sparse r×r linear map (binomial generator
  has at most 2 non-zero entries on each level of the matrix product).
* Total: length-`log_q(n)` chain of r×r linear maps over F_q.
* Polynomial multiplication of r-dim polys is in NC^1 (parallel
  addition tree of depth O(log r) = O(log log n)) but NOT in TC^0
  (constant depth) for growing r.
* Sequential chain of `log n` multiplications gives total depth
  `O(log n · log log n)` = polylog, not constant.

The `r × r` matrix-powering structure (E5.3) is unchanged by switching
the coefficient ring from Z_n to F_q. The q-twist is structural redress,
not a depth collapse.

## Closures filed

* `status/CLOSED_PATHS.md` line 719 added: "Healy-Viola Frobenius-mod-q
  transplant of AKS (FOCUS-1 sub-attack 3)" → FAIL, E+I.
* `status/SESSION_INSIGHTS.md` Session 64 entry added.
* No change to `OPEN_PROBLEMS.md` (Chain-E open frontier unchanged —
  growing-dim MPOW remains the only open primitive).
* No change to `EDGES.md` — no new edge surfaced; this closes a
  potential E1.x speculation about q-Frobenius leverage.
* No change to `novel/`.

## State of FOCUS-1 after Session 64

| Sub-attack | Status |
|-----------|--------|
| 1 — Bernstein 2003 strengthened gcd (smaller r) | **un-built** (only remaining sub-attack) |
| 2 — Non-cyclotomic ring AKS Z_n[x]/(x^d+a) | CLOSED (S61, line 714) |
| 3 — Healy-Viola Frobenius transplant | CLOSED (S64, line 719) |

After sub-attack 1 is built (or shown intractable), Chain E is
"computationally cornered" per TODO.md FOCUS-1: the only escapes are
Brandt MKtP (FOCUS-3) or a fundamentally new lower-bound technique.

## Files produced

* `experiments/circuit_complexity/aks_alternative/frobenius_transplant/frobenius_transplant.py` (≈8 KB)
* `experiments/circuit_complexity/aks_alternative/frobenius_transplant/frobenius_transplant_results.md` (≈4 KB)
* `archive/sessions/session64_frobenius_transplant.md` (this file)
